#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

! get rid of autoindent   :setl noai nocin nosi inde=
#define STANDALONE 0

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "MARCHING_TETRA_F.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

#define AREAZERO (1.0D-12)

      module marching_tetra_module

      contains

      subroutine vinterp(valu,gridval,gridx,gridy,gridz, &
        IV0,IV1,XX,YY,ZZ,icomp)

      IMPLICIT NONE

      REAL_T, intent(in) :: valu
      REAL_T, intent(in) :: gridx(8),gridy(8),gridz(8),gridval(8)
      INTEGER_T, intent(in) :: IV0,IV1
      REAL_T, intent(out) :: XX(SDIM),YY(SDIM),ZZ(SDIM)
      INTEGER_T, intent(in) :: icomp
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
      end subroutine vinterp

! in 2d, NORMAL(3)=0
! fictitious node at x1,y1,1
      subroutine addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
      IMPLICIT NONE

      INTEGER_T, intent(inout) :: itri
      INTEGER_T, intent(in) :: imaxtri
      REAL_T, intent(inout) :: trianglelist(SDIM,imaxtri)
      REAL_T, intent(in) :: XX(SDIM),YY(SDIM),ZZ(SDIM),NORMAL(SDIM)
      REAL_T CPROD(SDIM),VEC1(3),VEC2(3),DOTPROD

#if (AMREX_SPACEDIM==3)
      VEC1(1)=XX(SDIM)-XX(1)
      VEC1(2)=YY(SDIM)-YY(1)
      VEC1(3)=ZZ(SDIM)-ZZ(1)
#elif (AMREX_SPACEDIM==2)
       ! assume the "triangle element" in 2D is a fictitious 3D element
       ! which is perfectly vertical.
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
      VEC2(3)=zero  ! since it is assumed that ZZ(2)=ZZ(1)
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

      trianglelist(1,itri+1)=XX(1)
      trianglelist(2,itri+1)=YY(1)
#if (AMREX_SPACEDIM==3)
      trianglelist(SDIM,itri+1)=ZZ(1)
#endif
#if (AMREX_SPACEDIM==3)
      if (DOTPROD.gt.0.0) then
       trianglelist(1,itri+2)=XX(2)
       trianglelist(2,itri+2)=YY(2)
       trianglelist(SDIM,itri+2)=ZZ(2)
       trianglelist(1,itri+SDIM)=XX(SDIM)
       trianglelist(2,itri+SDIM)=YY(SDIM)
       trianglelist(SDIM,itri+SDIM)=ZZ(SDIM)
      else if (DOTPROD.le.zero) then
       trianglelist(1,itri+2)=XX(SDIM)
       trianglelist(2,itri+2)=YY(SDIM)
       trianglelist(3,itri+2)=ZZ(SDIM)
       trianglelist(1,itri+SDIM)=XX(2)
       trianglelist(2,itri+SDIM)=YY(2)
       trianglelist(SDIM,itri+SDIM)=ZZ(2)
      else
       print *,"DOTPROD is NaN"
       stop
      endif
#elif (AMREX_SPACEDIM==2)
      trianglelist(1,itri+2)=XX(2)
      trianglelist(2,itri+2)=YY(2)
#else
      print *,"dimension bust"
      stop
#endif


      itri=itri+SDIM

      return
      end subroutine addtrianglelist

      subroutine polysegment(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,IV0,IV1,IV2,IV3,imaxtri)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(inout) :: itri
      INTEGER_T, intent(in) :: IV0,IV1,IV2,IV3
      INTEGER_T, intent(in) :: imaxtri
      REAL_T, intent(inout) :: trianglelist(SDIM,imaxtri)
      REAL_T, intent(in) :: valu
      REAL_T, intent(in) :: gridx(8),gridy(8),gridz(8),gridval(8)
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

        ! NORMAL is proportional to the gradient of the level
        ! set function, considering IV0,IV1,IV2
        ! (levelset function is constant in the z direction for
        ! 2D case)
       if (istat.ne.0) then
        NORMAL(1)=sourcex(1)
        NORMAL(2)=sourcex(2)
        if ((idxtri.eq.1).or.(idxtri.eq.6)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.2).or.(idxtri.eq.5)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.3).or.(idxtri.eq.4)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        endif
       endif ! istat.ne.zero
      endif ! idxtri.ne.7 or 0

      return
      end subroutine polysegment

      subroutine polytri(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,IV0,IV1,IV2,IV3,imaxtri)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(inout) :: itri
      INTEGER_T, intent(in) :: IV0,IV1,IV2,IV3
      INTEGER_T, intent(in) :: imaxtri
      REAL_T, intent(inout) :: trianglelist(SDIM,imaxtri)
      REAL_T, intent(in) :: valu
      REAL_T, intent(in) :: gridx(8),gridy(8),gridz(8),gridval(8)
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
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.2).or.(idxtri.eq.13)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,3)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.3).or.(idxtri.eq.12)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV3,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,3)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
 
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,3)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.4).or.(idxtri.eq.11)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV1,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.5).or.(idxtri.eq.10)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV3,XX,YY,ZZ,3)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
  
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.6).or.(idxtri.eq.9)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
 
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.7).or.(idxtri.eq.8)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV3,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV3,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV3,IV1,XX,YY,ZZ,3)
         call addtrianglelist(trianglelist,itri,imaxtri,XX,YY,ZZ,NORMAL)
        endif
       endif ! istat.ne.zero
      endif ! idxtri.ne.15 or 0

      return
      end subroutine polytri

      subroutine copy_levelset_node_to_gridval( &
           levelset_node,gridval,ISUM,kkhi,nodehi,valu)
      IMPLICIT NONE

      REAL_T, intent(in) :: levelset_node(0:1,0:1,0:1)
      REAL_T, intent(out) :: gridval(8)
      INTEGER_T :: N8(8)
      INTEGER_T, intent(out) :: ISUM
      INTEGER_T, intent(in) :: kkhi,nodehi
      REAL_T, intent(in) :: valu
      INTEGER_T :: itemp

      gridval(1)=levelset_node(0,0,0)
      gridval(2)=levelset_node(1,0,0)
      gridval(3)=levelset_node(1,1,0)
      gridval(4)=levelset_node(0,1,0)
      gridval(5)=levelset_node(0,0,kkhi)
      gridval(6)=levelset_node(1,0,kkhi)
      gridval(7)=levelset_node(1,1,kkhi)
      gridval(8)=levelset_node(0,1,kkhi)

      ISUM=0
      do itemp=1,nodehi
       if (gridval(itemp).ge.valu) then
        N8(itemp)=1
       else 
        N8(itemp)=0
       endif
       ISUM=ISUM+N8(itemp)
      enddo

      end subroutine copy_levelset_node_to_gridval

      subroutine add_to_triangle_list( &
           gridx,gridy,gridz,gridval,valu,trianglelist, &
           itri,imaxtri,xnode,kkhi,nodehi)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: kkhi
      INTEGER_T, intent(in) :: nodehi
      REAL_T, intent(inout) :: gridval(8)
      REAL_T, intent(inout) :: gridx(8)
      REAL_T, intent(inout) :: gridy(8)
      REAL_T, intent(inout) :: gridz(8)
      REAL_T, intent(in) :: valu
      REAL_T, intent(inout) :: trianglelist(SDIM,200)
      REAL_T, intent(in) :: xnode(0:1,0:1,0:1,SDIM)
      INTEGER_T, intent(inout) :: itri
      INTEGER_T, intent(inout) :: imaxtri
      INTEGER_T :: dir

      if (SDIM.eq.2) then
       if (kkhi.eq.0) then
        ! do nothing
       else
        print *,"kkhi invalid"
        stop
       endif
      else if (SDIM.eq.3) then
       ! do nothing
      else
       print *,"SDIM invalid"
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
       call polysegment(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,1,2,3,3,imaxtri)
       call polysegment(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,1,4,3,3,imaxtri)
      else if (SDIM.eq.3) then
       call polytri(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,1,3,4,8,imaxtri)
       call polytri(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,1,3,7,8,imaxtri)
       call polytri(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,1,5,7,8,imaxtri)
       call polytri(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,1,7,2,3,imaxtri)
       call polytri(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,1,7,2,5,imaxtri)
       call polytri(gridx,gridy,gridz,gridval,valu,trianglelist, &
        itri,6,7,2,5,imaxtri)
      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine add_to_triangle_list

#if (STANDALONE==0)

      subroutine fort_isogrid( &
       tid, &
       visual_tessellate_vfrac, &
       recon,DIMS(recon), &
       xlo,dx, &
       mask,DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level,gridno,nmat) &
      bind(c,name='fort_isogrid')

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: visual_tessellate_vfrac
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: level,gridno
      INTEGER_T :: im
       ! vof,refcen,order,slope,int
      REAL_T, intent(in), target :: recon(DIMV(recon),nmat*ngeom_recon) 
      REAL_T, pointer :: recon_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))
      REAL_T valu
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      character*2 matstr

      character*25 namestr25 !./temptecplot ..
      character*29 cennamestr29 !./temptecplot/ ...
      character*3 levstr
      character*5 gridstr
      character*33 filename33
      character*37 cenfilename37 ! ./temptecplot ...

      real(8), dimension(:,:), allocatable :: Node  ! dir,index
      integer(4), dimension(:,:), allocatable :: IntElem  ! 1 or 2, index
      integer(4) :: NumNodes
      integer(4) :: NumIntElems
      integer(4) :: CurNodes
      integer(4) :: CurIntElems
      INTEGER_T imaxtri,itri,curtri,ipass
      REAL_T trianglelist(SDIM,200)
      REAL_T lnode(0:1,0:1,0:1)
      REAL_T xnode(0:1,0:1,0:1,SDIM)
      INTEGER_T i,j,k,ii,jj,kk,dir,ISUM,itemp,jtemp
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

      REAL_T, dimension(D_DECL(:,:,:),:), target, allocatable :: reconlocal
      REAL_T, pointer :: reconlocal_ptr(D_DECL(:,:,:),:)

      INTEGER_T nmax
      INTEGER_T vofcomp
      REAL_T vfrac_sum_solid
      REAL_T mofdata(nmat*ngeom_recon)
      INTEGER_T mask_local
      INTEGER_T vcdeb,imdeb

      nhalf=3
      nmax=POLYGON_LIST_MAX ! in: fort_isogrid

      recon_ptr=>recon
      mask_ptr=>mask

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
      reconlocal_ptr=>reconlocal

      call checkbound_array(fablo,fabhi,reconlocal_ptr,1,-1,4111)
      call checkbound_array(fablo,fabhi,recon_ptr,1,-1,4112)
      call checkbound_array1(fablo,fabhi,mask_ptr,0,-1,4113)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
       do dir=1,nmat*ngeom_recon
        mofdata(dir)=recon(D_DECL(i,j,k),dir)
       enddo

       if ((visual_tessellate_vfrac.eq.1).or. &
           (visual_tessellate_vfrac.eq.3)) then
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
          print *,"is_rigid invalid MARCHING_TETRA_3D.F90"
          stop
         endif
        enddo ! im=1..nmat

        if ((vfrac_sum_solid.gt.VOFTOL).and. &
            (vfrac_sum_solid.le.1.1)) then
         ! before (mofdata): fluids tessellate
         ! after  (mofdata): fluids and solids tessellate
         call multi_get_volume_tessellate( &
          visual_tessellate_vfrac, &
          bfact, &
          dx,xsten,nhalf, &
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

#if (STANDALONE==0)
       write(namestr25,'(A14,A7,A2,A2)') &
               './temptecplot/','tempmat',matstr,'ls'
       write(cennamestr29,'(A14,A10,A2,A3)') &
               './temptecplot/','temprefcen',matstr,'pos'
#elif (STANDALONE==1)
       write(namestr25,'(A14,A7,A2,A2)') &
               './temptecplot_','tempmat',matstr,'ls'
       write(cennamestr29,'(A14,A10,A2,A3)') &
               './temptecplot_','temprefcen',matstr,'pos'
#else
       print *,"STANDALONE invalid"
       stop
#endif

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
        ! namestr25=./temptecplot ...
        ! filename33=./temptecplot ..
       write(filename33,'(A25,A3,A5)') namestr25,levstr,gridstr
       print *,"filename33 ",filename33
        !cennamestr29="./temptecplot/ ... "
        !cenfilename37="./temptecplot/ ... "
       write(cenfilename37,'(A29,A3,A5)') cennamestr29,levstr,gridstr
       print *,"cenfilename37 ",cenfilename37

       NumNodes=0
       NumIntElems=0

       open(unit=12,file=cenfilename37) !./temptecplot ...
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
             ! declared in GLOBALUTIL.F90
            call distfunc(bfact,dx,xsten,nhalf, &
             intercept,nn,xtarget,lnode(ii,jj,kk),SDIM)
           endif
          enddo
          enddo
          enddo ! ii,jj,kk

          call copy_levelset_node_to_gridval( &
            lnode,gridval,ISUM,kkhi,nodehi,valu)

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
        
          call add_to_triangle_list( &
           gridx,gridy,gridz,gridval,valu,trianglelist, &
           itri,imaxtri,xnode,kkhi,nodehi)

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
              Node(dir,CurNodes)=trianglelist(dir,curtri)
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

         ! filename33=./temptecplot ...
       open(unit=11,file=filename33)
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
      end subroutine fort_isogrid


      subroutine fort_combinetriangles( &
       grids_per_level,finest_level,nsteps,im, &
       arrdim,time,plotint) &
      bind(c,name='fort_combinetriangles')
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE


      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: plotint
      INTEGER_T :: strandid
      INTEGER_T, intent(in) :: arrdim,finest_level,nsteps,im
      INTEGER_T, intent(in) :: grids_per_level(arrdim)

      character*25 namestr25 !./temptecplot ...
      character*7 newnamestr7 ! mat??ls ...

      character*29 cennamestr29 !./temptecplot ...
      character*11 newcennamestr11 ! refcen ...

      character*3 levstr
      character*5 gridstr

      character*33 filename33
      character*37 cenfilename37 !./temptecplot ...

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
      
#if (STANDALONE==0)
      write(namestr25,'(A14,A7,A2,A2)') &
              './temptecplot/','tempmat',matstr,'ls'
      write(cennamestr29,'(A14,A10,A2,A3)') &
              './temptecplot/','temprefcen',matstr,'pos'
#elif (STANDALONE==1)
      write(namestr25,'(A14,A7,A2,A2)') &
              './temptecplot_','tempmat',matstr,'ls'
      write(cennamestr29,'(A14,A10,A2,A3)') &
              './temptecplot_','temprefcen',matstr,'pos'
#else
       print *,"STANDALONE invalid"
       stop
#endif

      write(newnamestr7,'(A3,A2,A2)') 'mat',matstr,'ls'
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
         !newnamestr7=mat??ls ...
         !newfilename17=mat??ls ...
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
!        write(11,'(A33,D25.16,A10,I10)')  &
         write(11,'(A33,E25.16,A10,I10)')  &
          'ZONETYPE=FETRIANGLE SOLUTIONTIME=',round_time(time), &
          " STRANDID=",strandid
        else if (SDIM.eq.2) then
!        write(11,'(A32,D25.16,A10,I10)')  &
         write(11,'(A32,E25.16,A10,I10)')  &
          'ZONETYPE=FELINESEG SOLUTIONTIME=',round_time(time), &
          " STRANDID=",strandid
        else
         print *,"dimension bust"
         stop
        endif

          ! newcennamestr='refcen ...'
          ! newcenfilename21='refcen ...'
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

!       write(12,'(A19,I14,A26,D25.16,A10,I10)') & 
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

          ! ./temptecplot_tempmat
         write(filename33,'(A25,A3,A5)') namestr25,levstr,gridstr
         print *,"filename33 ",filename33
         open(unit=4,file=filename33)

         read(4,*) PartNumNodes
         read(4,*) PartNumIntElems

          !cennamestr29='./temptecplot ..'
          !cenfilename37='./temptecplot ..'
         write(cenfilename37,'(A29,A3,A5)') cennamestr29,levstr,gridstr
         print *,"cenfilename37 ",cenfilename37
         open(unit=5,file=cenfilename37)

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

       rmcommand_mat='rm ./temptecplot_tempmat*'
       rmcommand_refcen='rm ./temptecplot_temprefcen*'
#if (STANDALONE==0)
       ! do nothing
#elif (STANDALONE==1)
       print *,"issuing command ",rmcommand_mat
       print *,"and issuing command ",rmcommand_refcen

       call execute_command_line(rmcommand_mat,exitstat=sysret)
       call execute_command_line(rmcommand_refcen,exitstat=sysret)
#else
       print *,"STANDALONE invalid"
       stop
#endif
      endif

      if (sysret.ne.0) then
       print *,"execute_command_line has sysret=",sysret
       stop
      endif
 
      return
      end subroutine fort_combinetriangles

#endif

       ! For finding closest point map, the levelset function
       ! will have to have ngrow ghost cells where ngrow=thickness of the
       ! narrow band.
      subroutine fort_isogridsingle( &
       levelset, &
       DIMS(levelset), &
       xlo,dx, &
       mask, &
       DIMS(mask), &
       lo,hi, &
       bfact, &
       level, &
       gridno) &
      bind(c,name='fort_isogridsingle')

      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: lo(SDIM), hi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: DIMDEC(levelset)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: level,gridno
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in), target :: levelset(DIMV(levelset))
      REAL_T, pointer :: levelset_ptr(D_DECL(:,:,:))
      REAL_T, intent(in), target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      REAL_T valu

      character*3 levstr
      character*5 gridstr
      character*32 filename32

      real(8), dimension(:,:), allocatable :: Node  ! dir,index
      integer(4), dimension(:,:), allocatable :: IntElem  ! 1 or 2, index
      integer(4) :: NumNodes
      integer(4) :: NumIntElems
      integer(4) :: CurNodes
      integer(4) :: CurIntElems
      INTEGER_T imaxtri,itri,curtri,ipass
      REAL_T trianglelist(SDIM,200)
      REAL_T lnode(0:1,0:1,0:1)
      REAL_T xnode(0:1,0:1,0:1,SDIM)
      INTEGER_T i,j,k,ii,jj,kk,dir,ISUM,itemp,jtemp
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

      nhalf=3
      if (bfact.lt.1) then
       print *,"bfact invalid151"
       stop
      endif

      levelset_ptr=>levelset
      mask_ptr=>mask

      imaxtri=200
      call checkbound_array1(lo,hi,levelset_ptr,1,-1,4114)
      call checkbound_array1(lo,hi,mask_ptr,0,-1,4115)

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
#if (STANDALONE==0)
      write(filename32,'(A14,A10,A3,A5)') &
            './temptecplot/','templssing',levstr,gridstr
#elif (STANDALONE==1)
      write(filename32,'(A14,A10,A3,A5)') &
            './temptecplot_','templssing',levstr,gridstr
#else
      print *,"STANDALONE invalid"
      stop
#endif
      print *,"filename32 ",filename32

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
           if (SDIM.eq.3) then
            ! do nothing
           else
            print *,"SDIM invalid"
            stop
           endif
           dir=SDIM 
           xtarget(dir)=xsten(1,dir)
          endif
           ! xnode: nodes of the hexahedron corresponding to cell (i,j,k)
           ! xnode(0,0,0,dir)=the dir'th component of the position
           ! of the node with the smallest coordinates with respect
           ! to the other 7 nodes.
           ! xnode(1,1,1,dir)=the dir'th component of the position
           ! of the node with the largest coordinates with respect
           ! to the other 7 nodes.
          do dir=1,SDIM
           xnode(ii,jj,kk,dir)=xtarget(dir)
          enddo

          lnode(ii,jj,kk)=( &
            levelset(D_DECL(i+ii-1,j+jj-1,k+kk-1))+ &
            levelset(D_DECL(i+ii,j+jj-1,k+kk-1))+ &
            levelset(D_DECL(i+ii-1,j+jj,k+kk-1))+ &
            levelset(D_DECL(i+ii,j+jj,k+kk-1))+ &
            levelset(D_DECL(i+ii-1,j+jj-1,k+kk))+ &
            levelset(D_DECL(i+ii,j+jj-1,k+kk))+ &
            levelset(D_DECL(i+ii-1,j+jj,k+kk))+ &
            levelset(D_DECL(i+ii,j+jj,k+kk)))/eight

         enddo
         enddo
         enddo

         call copy_levelset_node_to_gridval( &
           lnode,gridval,ISUM,kkhi,nodehi,valu)

         if ((ISUM.eq.0).or.(ISUM.eq.nodehi)) then
          itri=0
         else

          call add_to_triangle_list( &
           gridx,gridy,gridz,gridval,valu,trianglelist, &
           itri,imaxtri,xnode,kkhi,nodehi)

         endif

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
             Node(dir,CurNodes)=trianglelist(dir,curtri)
            enddo
            IntElem(jtemp,CurIntElems)=CurNodes
           enddo  ! jtemp=1..sdim
          enddo  ! itemp=1..itri/sdim
         else
          print *,"ipass invalid"
          stop
         endif
  
        endif ! mask=1
       enddo
       enddo
       enddo
      enddo ! ipass

      open(unit=11,file=filename32)
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
      end subroutine fort_isogridsingle


      subroutine fort_closest_point_map( &
       tid, & ! thread id
       sweep, & ! sweep=0 or 1
       ngrow, & 
       minLS, &
       maxLS, &
       ls_old, &
       DIMS(ls_old), &
       ls_new, &
       DIMS(ls_new), &
       ls_grad_new, &
       DIMS(ls_grad_new), &
       xlo,dx, &
       mask, &
       DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level) &
      bind(c,name='fort_closest_point_map')

      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: sweep ! sweep=0 or 1
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(ls_old) !ls_oldxlo,ls_oldxhi,...
      INTEGER_T, intent(in) :: DIMDEC(ls_new)
      INTEGER_T, intent(in) :: DIMDEC(ls_grad_new)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: level
      REAL_T, intent(in), target :: ls_old(DIMV(ls_old))
      REAL_T, intent(inout), target :: ls_new(DIMV(ls_new))
      REAL_T, intent(inout), target ::  &
        ls_grad_new(DIMV(ls_grad_new),SDIM)
      REAL_T, pointer :: ls_old_ptr(D_DECL(:,:,:))
      REAL_T, pointer :: ls_new_ptr(D_DECL(:,:,:))
      REAL_T, pointer :: ls_grad_new_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      REAL_T, intent(inout) :: minLS
      REAL_T, intent(inout) :: maxLS

      REAL_T valu

      INTEGER_T imaxtri,itri
      REAL_T trianglelist(SDIM,200)
      REAL_T lnode(0:1,0:1,0:1)
      REAL_T xnode(0:1,0:1,0:1,SDIM)
      INTEGER_T i,j,k
      INTEGER_T ilocal,jlocal,klocal
      INTEGER_T ii,jj,kk,ISUM
      INTEGER_T ibase
      INTEGER_T local_dir
      REAL_T gridval(8)
      REAL_T gridx(8)
      REAL_T gridy(8)
      REAL_T gridz(8)

      REAL_T xtarget(SDIM)
      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_local(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T kklo,kkhi,nodehi,ngrow_k
      REAL_T degenerate_face_tol
      REAL_T check_tol
      REAL_T dxmin
      REAL_T xcc(3)
      REAL_T save_LS
      REAL_T initial_LS
      REAL_T p_triangle(3,3) ! (ipoint,dir)
      INTEGER_T in_plane
      REAL_T xcp_0(3)
      REAL_T xcp_0_project(3)
      REAL_T dist_pij(3)
      REAL_T tan_vec(2,3)
      REAL_T tan1(3)
      REAL_T tan2(3)
      REAL_T n_triangle(3)
      REAL_T normal_closest(3)
      REAL_T n_magnitude
      REAL_T phi_x
      REAL_T save_x_cp(SDIM)
      INTEGER_T ipoint
      INTEGER_T ipoint_p1
      REAL_T xnode_seg(2,3)
      REAL_T nnode_seg(2,3)
      REAL_T xnode_point(3)
      REAL_T dist_point
      REAL_T unsigned_mindist
    
      if (tid.ge.0) then
       ! do nothing
      else 
       print *,"tid invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid151"
       stop
      endif

      call get_dxmin(dx,bfact,dxmin)
      degenerate_face_tol=1.0D-8*dxmin
 
      nhalf=3

      ls_old_ptr=>ls_old
      ls_new_ptr=>ls_new
      ls_grad_new_ptr=>ls_grad_new
      mask_ptr=>mask

      imaxtri=200

      if (ngrow.ge.1) then
       ! do nothing
      else
       print *,"ngrow invalid"
       stop
      endif

      call checkbound_array1(fablo,fabhi,ls_old_ptr,ngrow,-1,4116)
      call checkbound_array1(fablo,fabhi,ls_new_ptr,0,-1,4117)
      call checkbound_array(fablo,fabhi,ls_grad_new_ptr,0,-1,4118)
      call checkbound_array1(fablo,fabhi,mask_ptr,0,-1,4119)

        ! ok for 2D or 3D
        ! find closest point on the interface for each interior
        ! point x_{i,j,k}
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

         ! xcc is the coordinate at which a corresponding closest point
         ! is sought.
       do local_dir=1,SDIM
        xcc(local_dir)=xsten(0,local_dir)
       enddo
       if (SDIM.eq.3) then
        ! do nothing
       else if (SDIM.eq.2) then
        xcc(3)=zero
       else
        print *,"sdim invalid"
        stop
       endif

        ! sweep==0: find closest point for cells in which "ngrow" neighborhood
        !   contains a marching tetrahedra generated interface.
        ! sweep==1: set uninit distances to maxLS,minLS respectively.
       if (sweep.eq.0) then

        valu=zero

        if (SDIM.eq.3) then
         ngrow_k=ngrow-1
         kklo=0
         kkhi=1
         nodehi=8
        else if (SDIM.eq.2) then
         ngrow_k=0
         kklo=0
         kkhi=0
         nodehi=4
        else
         print *,"dimension bust"
         stop
        endif

        ! only process cell (i,j,k) if the cell is not covered 
        ! by a finer level.
        if (mask(D_DECL(i,j,k)).eq.one) then

         do local_dir=1,SDIM
          ls_grad_new(D_DECL(i,j,k),local_dir)=zero
          save_x_cp(local_dir)=xcc(local_dir)
         enddo

         save_LS=ls_old(D_DECL(i,j,k))
         if (save_LS.gt.zero) then
          initial_LS=1.0D+20
         else if (save_LS.lt.zero) then
          initial_LS=-1.0D+20 
         else if (save_LS.eq.zero) then
          initial_LS=zero
          ls_new(D_DECL(i,j,k))=zero
          do local_dir=1,SDIM
           ls_grad_new(D_DECL(i,j,k),local_dir)=zero
          enddo
         else
          print *,"save_LS is NaN"
          stop
         endif

         if (save_LS.eq.zero) then
          ! do nothing
         else if (save_LS.ne.zero) then

          do ilocal=i-ngrow+1,i+ngrow-1
          do jlocal=j-ngrow+1,j+ngrow-1
          do klocal=k-ngrow_k,k+ngrow_k

            ! klocal=0 in 2D
           if (SDIM.eq.2) then
            if (klocal.eq.0) then
             ! do nothing
            else
             print *,"klocal invalid"
             stop
            endif
           else if (SDIM.eq.3) then
            ! check nothing
           else
            print *,"SDIM invalid"
            stop
           endif

           call gridsten(xsten_local,xlo,ilocal,jlocal,klocal, &
             fablo,bfact,dx,nhalf)

           itri=0
           do ii=0,1
           do jj=0,1
           do kk=kklo,kkhi

            if (SDIM.eq.2) then
             if (kk.eq.0) then
              ! do nothing
             else
              print *,"kk invalid"
              stop
             endif
            else if (SDIM.eq.3) then
             ! check nothing
            else
             print *,"SDIM invalid"
             stop
            endif

             ! initialize xtarget to be the LLL corner of the box.
            do local_dir=1,SDIM
             xtarget(local_dir)=xsten_local(-1,local_dir)
            enddo
            if (ii.eq.1) then
             local_dir=1 
             xtarget(local_dir)=xsten_local(1,local_dir)
            endif
            if (jj.eq.1) then
             local_dir=2 
             xtarget(local_dir)=xsten_local(1,local_dir)
            endif

            if (kk.eq.1) then
             if (SDIM.eq.3) then
              ! do nothing
             else
              print *,"SDIM invalid"
              stop
             endif
             local_dir=SDIM 
             xtarget(local_dir)=xsten_local(1,local_dir)
            endif

            do local_dir=1,SDIM
             xnode(ii,jj,kk,local_dir)=xtarget(local_dir)
            enddo

              ! D_DECL(i,j,k)  => "i,j,k" in 3D
              !                   "i,j" in 2D
            lnode(ii,jj,kk)=( &
              ls_old(D_DECL(ilocal+ii-1,jlocal+jj-1,klocal+kk-1))+ &
              ls_old(D_DECL(ilocal+ii,jlocal+jj-1,klocal+kk-1))+ &
              ls_old(D_DECL(ilocal+ii-1,jlocal+jj,klocal+kk-1))+ &
              ls_old(D_DECL(ilocal+ii,jlocal+jj,klocal+kk-1))+ &
              ls_old(D_DECL(ilocal+ii-1,jlocal+jj-1,klocal+kk))+ &
              ls_old(D_DECL(ilocal+ii,jlocal+jj-1,klocal+kk))+ &
              ls_old(D_DECL(ilocal+ii-1,jlocal+jj,klocal+kk))+ &
              ls_old(D_DECL(ilocal+ii,jlocal+jj,klocal+kk)))/eight

           enddo
           enddo
           enddo

            !gridval,ISUM are outputs.
            !3D data assumed
           call copy_levelset_node_to_gridval( &
            lnode,gridval,ISUM,kkhi,nodehi,valu)

           if ((ISUM.eq.0).or.(ISUM.eq.nodehi)) then
            itri=0
           else
            call add_to_triangle_list( &
             gridx,gridy,gridz,gridval,valu,trianglelist, &
             itri,imaxtri,xnode,kkhi,nodehi)
           endif

            ! REAL_T :: trianglelist(SDIM,200)
            ! in 2D:
            !  first segment:
            ! trianglelist(1,1)
            ! trianglelist(2,1)
            ! trianglelist(1,2)
            ! trianglelist(2,2)
            !  second segment:
            ! trianglelist(1,3)
            ! trianglelist(2,3)
            ! trianglelist(1,4)
            ! trianglelist(2,4)
            !  ...
           ibase=0
           do while (ibase.lt.itri)

             ! line segment in 2D, triangle in 3D
             ! the ordering of the points does not matter here
             ! since we already know from "save_LS" what the sign
             ! of the level set function is.
            do ipoint=1,SDIM
             do local_dir=1,SDIM
              p_triangle(ipoint,local_dir)= &
                trianglelist(local_dir,ibase+ipoint)
             enddo
             if (SDIM.eq.2) then
              p_triangle(ipoint,3)=zero
             else if (SDIM.eq.3) then
              ! do nothing
             else
              print *,"sdim invalid"
              stop
             endif
            enddo ! ipoint=1..sdim
            if (SDIM.eq.2) then
             do local_dir=1,SDIM
              p_triangle(3,local_dir)= &
                trianglelist(local_dir,ibase+2)
             enddo
             p_triangle(3,3)=one
            else if (SDIM.eq.3) then
             ! do nothing
            else
             print *,"sdim invalid"
             stop
            endif
            do ipoint=1,3
             ipoint_p1=ipoint+1
             if (ipoint_p1.gt.3) then
              ipoint_p1=1
             endif 
             dist_pij(ipoint)=zero
             do local_dir=1,3
              dist_pij(ipoint)=dist_pij(ipoint)+ &
               (p_triangle(ipoint,local_dir)- &
                p_triangle(ipoint_p1,local_dir))**2
             enddo
             if (dist_pij(ipoint).ge.zero) then
              dist_pij(ipoint)=sqrt(dist_pij(ipoint))
             else
              print *,"dist_pij is NaN"
              stop
             endif
            enddo !ipoint=1,3

            in_plane=0
            check_tol=zero

            do local_dir = 1,3
             n_triangle(local_dir) = zero
             xcp_0(local_dir)=zero
             xcp_0_project(local_dir)=zero
            enddo
            unsigned_mindist=1.0D+20

            if ((dist_pij(1).ge.degenerate_face_tol).and. &
                (dist_pij(2).ge.degenerate_face_tol).and. &
                (dist_pij(3).ge.degenerate_face_tol)) then

             do ipoint=1,2
              do local_dir=1,3
               tan_vec(ipoint,local_dir)=p_triangle(ipoint,local_dir)- &
                      p_triangle(ipoint+1,local_dir)
              enddo
             enddo
             do local_dir=1,3
              tan1(local_dir)=tan_vec(1,local_dir)
              tan2(local_dir)=tan_vec(2,local_dir)
             enddo
             call crossprod(tan1,tan2,n_triangle)

             n_magnitude = &
                sqrt(n_triangle(1)**2+n_triangle(2)**2+n_triangle(3)**2)
             if (n_magnitude.gt.zero) then
              do local_dir = 1,3
               n_triangle(local_dir) = n_triangle(local_dir)/n_magnitude 
              enddo
             else
              print *,"n_magnitude invalid"
              stop
             endif
           
             phi_x=zero
             do local_dir = 1,SDIM
              phi_x=phi_x+ &
                n_triangle(local_dir)*(xcc(local_dir)- &
                   p_triangle(1,local_dir))
             enddo
            
             do local_dir = 1,SDIM
              xcp_0(local_dir) = xcc(local_dir)-phi_x*n_triangle(local_dir) 
             enddo
             if (SDIM.eq.2) then
              xcp_0(3)=zero
             else if (SDIM.eq.3) then
              ! do nothing
             else
              print *,"sdim invalid"
              stop
             endif
         
             do local_dir=1,3
              normal_closest(local_dir)=n_triangle(local_dir)
             enddo 
             call global_checkinplane(p_triangle,xcp_0,check_tol, &
                    xcp_0_project,in_plane)

             if (in_plane.eq.0) then
              ! do nothing
             else if (in_plane.eq.1) then
              unsigned_mindist=zero
              do local_dir=1,SDIM
               unsigned_mindist=unsigned_mindist+ &
                  (xcp_0_project(local_dir)-xcc(local_dir))**2
              enddo
              unsigned_mindist=sqrt(unsigned_mindist)
             else
              print *,"in_plane invalid"
              stop
             endif

             if (unsigned_mindist.ge.zero) then
              ! do nothing
             else
              print *,"unsigned_mindist invalid"
              stop
             endif

            else if ((dist_pij(1).le.degenerate_face_tol).or. &
                     (dist_pij(2).le.degenerate_face_tol).or. &
                     (dist_pij(3).le.degenerate_face_tol)) then
             ! do nothing
            else
             print *,"dist_pij NaN"
             stop
            endif
              
            do ipoint=1,3
             ipoint_p1=ipoint+1
             if (ipoint_p1.gt.3) then
              ipoint_p1=1
             endif 
             do local_dir=1,3
              xnode_seg(1,local_dir)=p_triangle(ipoint,local_dir)
              xnode_seg(2,local_dir)=p_triangle(ipoint_p1,local_dir)
              nnode_seg(1,local_dir)=n_triangle(local_dir)
              nnode_seg(2,local_dir)=n_triangle(local_dir)
             enddo

             call global_checkinline(nnode_seg,xnode_seg,check_tol,xcc, &
               in_plane,unsigned_mindist,xcp_0_project,normal_closest)

             do local_dir=1,3
              xnode_point(local_dir)=p_triangle(ipoint,local_dir)
             enddo

             call global_xdist(xnode_point,xcc,dist_point)

             if (dist_point.ge.zero) then
              ! do nothing
             else
              print *,"dist_point invalid"
              stop
             endif

             if ((dist_point.lt.unsigned_mindist).or. &
                 (in_plane.eq.0)) then
              in_plane=1
              unsigned_mindist=dist_point
              do local_dir=1,3
               xcp_0_project(local_dir)=xnode_point(local_dir)
              enddo
             else if ((dist_point.ge.unsigned_mindist).and. &
                      (in_plane.eq.1)) then
              ! do nothing
             else
              print *,"dist_point or in_plane invalid"
              stop
             endif
            enddo !ipoint=1,3

            ibase=ibase+SDIM 

            if (unsigned_mindist.lt.abs(initial_LS)) then
             if (save_LS.gt.zero) then
              initial_LS=unsigned_mindist
             else if (save_LS.lt.zero) then
              initial_LS=-unsigned_mindist
             else
              print *,"save_LS invalid"
              stop
             endif
             do local_dir=1,SDIM
              save_x_cp(local_dir)=xcp_0_project(local_dir)
             enddo
            else if (unsigned_mindist.ge.abs(initial_LS)) then
             ! do nothing
            else
             print *,"unsigned_mindist is NaN"
             stop
            endif
           enddo ! while (ibase.lt.itri)
          enddo !enddo klocal
          enddo !enddo jlocal
          enddo !enddo ilocal

          if (minLS.gt.initial_LS) then
           minLS=initial_LS
          endif
          if (maxLS.lt.initial_LS) then
           maxLS=initial_LS
          endif
          ls_new(D_DECL(i,j,k))=initial_LS
          if (abs(initial_LS).gt.zero) then
           do local_dir=1,SDIM
            ls_grad_new(D_DECL(i,j,k),local_dir)= &
              (xcc(local_dir)-save_x_cp(local_dir))/initial_LS
           enddo
          else if (abs(initial_LS).eq.zero) then
           do local_dir=1,SDIM
            ls_grad_new(D_DECL(i,j,k),local_dir)=zero
           enddo
          else
           print *,"initial_LS is NaN"
           stop
          endif
         else
          print *,"save_LS is NaN"
          stop
         endif

        else if (mask(D_DECL(i,j,k)).eq.zero) then
         ! do nothing
        else
         print *,"mask invalid"
         stop
        endif 

       else if (sweep.eq.1) then

        !set uninit distances to maxLS,minLS respectively
        if (mask(D_DECL(i,j,k)).eq.one) then
         if (ls_new(D_DECL(i,j,k)).lt.minLS) then
          ls_new(D_DECL(i,j,k))=minLS
         endif
         if (ls_new(D_DECL(i,j,k)).gt.maxLS) then
          ls_new(D_DECL(i,j,k))=maxLS
         endif
        else if (mask(D_DECL(i,j,k)).eq.zero) then
         ! do nothing
        else
         print *,"mask invalid"
         stop
        endif 
       else
        print *,"sweep invalid"
        stop
       endif

      enddo !k
      enddo !j
      enddo !i
 
      return
      end subroutine fort_closest_point_map

      subroutine fort_combinetrianglessingle( &
       grids_per_level,finest_level,nsteps,arrdim) &
      bind(c,name='fort_combinetrianglessingle')
      IMPLICIT NONE

      INTEGER_T, intent(in) :: arrdim,finest_level,nsteps
      INTEGER_T, intent(in) :: grids_per_level(arrdim)

      character*3 levstr
      character*5 gridstr
      character*32 filename32

      character*6 stepstr
      character*16 newfilename16

      character*80 rmcommand_LS

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
#if (STANDALONE==0)
         write(filename32,'(A14,A10,A3,A5)') &
            './temptecplot/','templssing',levstr,gridstr
#elif (STANDALONE==1)
         write(filename32,'(A14,A10,A3,A5)') &
            './temptecplot_','templssing',levstr,gridstr
#else
         print *,"STANDALONE invalid"
         stop
#endif

         print *,"filename32 ",filename32
         open(unit=4,file=filename32)

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

          rmcommand_LS='rm ./temptecplot_templssing*'
          sysret=0

#if (STANDALONE==0)
          ! do nothing
#elif (STANDALONE==1)
          print *,"issuing command ",rmcommand_LS

          call execute_command_line(rmcommand_LS,exitstat=sysret)
          if (sysret.ne.0) then
           print *,"execute_command_line has sysret=",sysret
           stop
          endif
#else
          print *,"STANDALONE invalid"
          stop
#endif

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
      end subroutine fort_combinetrianglessingle


      end module marching_tetra_module

#undef STANDALONE
