#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

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

! probtype==411 (see run2d/inputs.CAVTEST2D or run3d/inputs.CAVTEST3D)
module CAV3D_module
use amrex_fort_module, only : amrex_real

implicit none                   

contains

   ! do any initial preparation needed
 subroutine INIT_CAV3D_MODULE()
 IMPLICIT NONE

 return
 end subroutine INIT_CAV3D_MODULE

 subroutine OPEN_CAV3D_CASFILE(part_id,unit_id,file_format)
 use probcommon_module
 IMPLICIT NONE

 integer, INTENT(in) :: part_id
 integer, INTENT(in) :: unit_id
 integer, INTENT(out) :: file_format
 character(40) :: casname

 file_format=0

 if ((num_materials.eq.3).and.(probtype.eq.411)) then

  if (part_id.eq.1) then
   casname="nozzle.cas"
  else
   print *,"part_id invalid"
   stop
  endif

  print *,"opening ",casname
  OPEN(unit=unit_id,file=casname,access='sequential', &
      form="formatted",status='old')
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine OPEN_CAV3D_CASFILE

  ! i    j   k 
  ! u1  u2  u3
  ! v1  v2  v3
  ! curl(1)=u2 v3 - v2 u3
  ! curl(2)=u3 v1 - u1 v3
  ! curl(3)=u1 v2 - u2 v1
 subroutine CAV3D_ORDER_NODES(nodes,nodemap)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real), INTENT(in) :: nodes(3,3) ! dir,nodenum
 integer, INTENT(inout) :: nodemap(3)
 real(amrex_real) a(3),b(3),c(3)
 integer inode
 integer dir
 integer dircrit
 real(amrex_real) mag
 real(amrex_real) xcrit
 real(amrex_real) ncrit
 integer swap

 if ((num_materials.eq.3).and.(probtype.eq.411)) then

  do dir=1,3
   a(dir)=nodes(dir,2)-nodes(dir,1)
   b(dir)=nodes(dir,3)-nodes(dir,1)
  enddo
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=-(a(1)*b(3)-b(1)*a(3))
  c(3)=a(1)*b(2)-b(1)*a(2)
  mag=zero
  do dir=1,3
   mag=mag+c(dir)**2
  enddo
  mag=sqrt(mag)
  if (mag.gt.zero) then
   do dir=1,3
    c(dir)=c(dir)/mag
   enddo
   dircrit=1
   if (abs(c(2)).gt.abs(c(dircrit))) then
    dircrit=2
   endif 
   if (abs(c(3)).gt.abs(c(dircrit))) then
    dircrit=3
   endif 
   xcrit=zero
   do inode=1,3
    xcrit=xcrit+nodes(dircrit,inode)
   enddo
   xcrit=xcrit/three
   if (dircrit.eq.1) then
    if (abs(xcrit-5.0).le.VOFTOL) then
     ncrit=one
    else if (abs(xcrit+5.0).le.VOFTOL) then
     ncrit=-one
    else if (abs(xcrit-4.0).le.VOFTOL) then
     ncrit=-one
    else if (abs(xcrit+4.0).le.VOFTOL) then
     ncrit=one
    else if (abs(xcrit-2.06D0).le.VOFTOL) then
     ncrit=one
    else
     print *,"xcrit invalid(1) xcrit=",xcrit
     do inode=1,3
     do dir=1,3
      print *,"inode,dir,nodes(dir,inode) ",inode,dir,nodes(dir,inode)
     enddo
     enddo
     do inode=1,3
      print *,"inode,nodemap(inode) ",inode,nodemap(inode)
     enddo
     stop
    endif
   else if (dircrit.eq.2) then
    if (abs(xcrit-0.0).le.VOFTOL) then
     ncrit=one
    else if (abs(xcrit+16.0).le.VOFTOL) then
     ncrit=-one
    else if (abs(xcrit+8.0).le.VOFTOL) then
     ncrit=one
    else
     print *,"xcrit invalid(2) xcrit=",xcrit
     do inode=1,3
     do dir=1,3
      print *,"inode,dir,nodes(dir,inode) ",inode,dir,nodes(dir,inode)
     enddo
     enddo
     do inode=1,3
      print *,"inode,nodemap(inode) ",inode,nodemap(inode)
     enddo
     stop
    endif
   else if (dircrit.eq.3) then
    if (abs(xcrit-2.0).le.VOFTOL) then
     ncrit=one
    else if (abs(xcrit+2.0).le.VOFTOL) then
     ncrit=-one
    else if (abs(xcrit+0.97D0).le.VOFTOL) then
     ncrit=one
    else if (abs(xcrit-0.97D0).le.VOFTOL) then
     ncrit=-one
    else
     print *,"xcrit invalid(3) xcrit=",xcrit
     do inode=1,3
     do dir=1,3
      print *,"inode,dir,nodes(dir,inode) ",inode,dir,nodes(dir,inode)
     enddo
     enddo
     do inode=1,3
      print *,"inode,nodemap(inode) ",inode,nodemap(inode)
     enddo
     stop
    endif
   else
    print *,"dircrit invalid"
    stop
   endif
   if (ncrit*c(dircrit).gt.zero) then
    ! do nothing
   else if (ncrit*c(dircrit).lt.zero) then
    swap=nodemap(1)
    nodemap(1)=nodemap(3)
    nodemap(3)=swap
   else
    print *,"ncrit or c(dircrit) invalid"
    stop
   endif
  else
   print *,"mag invalid CAV3D_ORDER_NODES"
   stop
  endif

 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine CAV3D_ORDER_NODES

 subroutine CAV3D_SLICE(xmap3D,xslice3D,problo3D,probhi3D,dx_slice)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real), INTENT(in) :: dx_slice
 integer, INTENT(inout) :: xmap3D(3)
 real(amrex_real), INTENT(inout) :: xslice3D(3)
 real(amrex_real), INTENT(out) :: problo3D(3)
 real(amrex_real), INTENT(out) :: probhi3D(3)

 if ((num_materials.eq.3).and.(probtype.eq.411)) then

  if (dx_slice.gt.zero) then
   xmap3D(1)=1
   xmap3D(2)=2
   xmap3D(3)=0

   xslice3D(3)=zero
   problo3D(3)=-half*dx_slice
   probhi3D(3)=half*dx_slice
  else
   print *,"dx_slice invalid"
   stop
  endif

 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine CAV3D_SLICE

  ! fluids tessellate the domain, solids are immersed. 
 subroutine CAV3D_LS(x,t,LS,nmat)
 use probcommon_module
 IMPLICIT NONE

 integer, INTENT(in) :: nmat
 real(amrex_real), INTENT(in) :: x(SDIM)
 real(amrex_real), INTENT(in) :: t
 real(amrex_real), INTENT(out) :: LS(nmat)
 integer im

 if ((num_materials.eq.3).and.(probtype.eq.411)) then

  do im=1,num_materials
   if (im.eq.1) then !liquid
    LS(im)=99999.0
   else if (im.eq.2) then !ambient
    LS(im)=-99999.0
   else if (im.eq.num_materials) then ! geometry (placeholder)
    LS(im)=-99999.0
   else
    print *,"im invalid"
    stop
   endif
  enddo ! im=1..num_materials
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine CAV3D_LS

 subroutine CAV3D_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
 use probcommon_module
 IMPLICIT NONE

 integer, INTENT(in) :: nmat
 real(amrex_real), INTENT(in) :: x(SDIM)
 real(amrex_real), INTENT(in) :: t
 real(amrex_real), INTENT(in) :: dx(SDIM)
 real(amrex_real), INTENT(in) :: LS(nmat)
 real(amrex_real), INTENT(out) :: VEL(SDIM)
 integer, INTENT(in) :: velsolid_flag
 integer dir

 if ((velsolid_flag.eq.0).or. &
     (velsolid_flag.eq.1)) then
  ! do nothing
 else 
  print *,"velsolid_flag invalid"
  stop
 endif

  ! flow enters wide part (yhi) and exits narrow part (ylo)?
 if ((adv_dir.ge.1).and.(adv_dir.le.SDIM)) then
  do dir=1,SDIM
   VEL(dir)=zero
  enddo
  if ((abs(x(1)).lt.4.5).and.(velsolid_flag.eq.0)) then
   VEL(adv_dir)=-abs(adv_vel)
  else if ((abs(x(1)).ge.4.5).or.(velsolid_flag.eq.1)) then
   VEL(adv_dir)=zero
  else
   print *,"x(1) bust"
   stop
  endif
 else
  print *,"adv_dir invalid in CAV3D_VEL"
  stop
 endif

 return 
 end subroutine CAV3D_VEL

 subroutine CAV3D_PRES(x,t,LS,PRES,nmat)
 use probcommon_module
 IMPLICIT NONE

 integer, INTENT(in) :: nmat
 real(amrex_real), INTENT(in) :: x(SDIM)
 real(amrex_real), INTENT(in) :: t
 real(amrex_real), INTENT(in) :: LS(nmat)
 real(amrex_real), INTENT(out) :: PRES
 real(amrex_real) ymid

 ymid=half*(probloy+probhiy)
 if (x(2).le.ymid) then
  PRES=outflow_pressure
 else if (x(2).ge.ymid) then
  PRES=inflow_pressure
 else
  print *,"x(2) or ymid invalid"
  stop
 endif

 return 
 end subroutine CAV3D_PRES


 subroutine CAV3D_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
 use probcommon_module
 IMPLICIT NONE

 integer, INTENT(in) :: bcflag
 integer, INTENT(in) :: nmat
 integer, INTENT(in) :: nstate_mat
 real(amrex_real), INTENT(in) :: x(SDIM)
 real(amrex_real), INTENT(in) :: t
 real(amrex_real), INTENT(in) :: LS(nmat)
 real(amrex_real), INTENT(out) :: STATE(nmat*nstate_mat)
 integer im,ibase,n

 if ((num_materials.eq.3).and. &
     (num_state_material.ge.2).and. &
     (probtype.eq.411)) then
  do im=1,num_materials
   ibase=(im-1)*num_state_material
   STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im)
   if (t.eq.zero) then
    STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im)
   else if (t.gt.zero) then
    STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
   else
    print *,"t invalid"
    stop
   endif
   do n=1,num_species_var
    STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
   enddo
  enddo ! im=1..num_materials
 else
  print *,"num_materials,num_state_material, or probtype invalid"
  stop
 endif
  
 return
 end subroutine CAV3D_STATE

  ! dir=1..sdim  side=1..2
 subroutine CAV3D_LS_BC(xwall,xghost,t,LS, &
    LS_in,dir,side,dx,nmat)
 use probcommon_module
 IMPLICIT NONE

 integer, INTENT(in) :: nmat
 real(amrex_real), INTENT(in) :: xwall
 real(amrex_real), INTENT(in) :: xghost(SDIM)
 real(amrex_real), INTENT(in) :: t
 real(amrex_real), INTENT(inout) :: LS(nmat)
 real(amrex_real), INTENT(in) :: LS_in(nmat)
 integer, INTENT(in) :: dir,side
 real(amrex_real), INTENT(in) :: dx(SDIM)

 if ((dir.ge.1).and.(dir.le.SDIM).and. &
     (side.ge.1).and.(side.le.2)) then
  call CAV3D_LS(xghost,t,LS,nmat)
 else
  print *,"dir or side invalid"
  stop
 endif
 
 return
 end subroutine CAV3D_LS_BC


  ! dir=1..sdim  side=1..2 veldir=1..sdim
 subroutine CAV3D_VEL_BC(xwall,xghost,t,LS, &
    VEL,VEL_in,veldir,dir,side,dx,nmat)
 use probcommon_module
 IMPLICIT NONE

 integer, INTENT(in) :: nmat
 real(amrex_real), INTENT(in) :: xwall
 real(amrex_real), INTENT(in) :: xghost(SDIM)
 real(amrex_real), INTENT(in) :: t
 real(amrex_real), INTENT(in) :: LS(nmat)
 real(amrex_real), INTENT(inout) :: VEL
 real(amrex_real), INTENT(in) :: VEL_in
 integer, INTENT(in) :: veldir,dir,side
 real(amrex_real), INTENT(in) :: dx(SDIM)
 real(amrex_real) local_VEL(SDIM)
 integer velsolid_flag

 velsolid_flag=0
 if ((dir.ge.1).and.(dir.le.SDIM).and. &
     (side.ge.1).and.(side.le.2).and. &
     (veldir.ge.1).and.(veldir.le.SDIM)) then

  call CAV3D_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
  VEL=local_VEL(veldir)

 else
  print *,"dir,side, or veldir invalid"
  stop
 endif

 return
 end subroutine CAV3D_VEL_BC


  ! dir=1..sdim  side=1..2
 subroutine CAV3D_PRES_BC(xwall,xghost,t,LS, &
    PRES,PRES_in,dir,side,dx,nmat)
 use probcommon_module
 IMPLICIT NONE

 integer, INTENT(in) :: nmat
 real(amrex_real), INTENT(in) :: xwall
 real(amrex_real), INTENT(in) :: xghost(SDIM)
 real(amrex_real), INTENT(in) :: t
 real(amrex_real), INTENT(in) :: LS(nmat)
 real(amrex_real), INTENT(inout) :: PRES
 real(amrex_real), INTENT(in) :: PRES_in
 integer, INTENT(in) :: dir,side
 real(amrex_real), INTENT(in) :: dx(SDIM)

 if ((dir.ge.1).and.(dir.le.SDIM).and. &
     (side.ge.1).and.(side.le.2)) then

  call CAV3D_PRES(xghost,t,LS,PRES,nmat)

 else
  print *,"dir or side invalid"
  stop
 endif

 return
 end subroutine CAV3D_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine CAV3D_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real) :: local_STATE(nmat*num_state_material)
real(amrex_real), INTENT(inout) :: STATE
real(amrex_real), INTENT(inout) :: STATE_merge
real(amrex_real), INTENT(in) :: STATE_in
integer, INTENT(in) :: dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)
integer, INTENT(in) :: istate,im

integer ibase,im_crit
integer local_bcflag

 if (nmat.eq.num_materials) then
  ! do nothing
 else
  print *,"nmat invalid"
  stop
 endif
 local_bcflag=1

 if ((istate.ge.1).and. &
     (istate.le.num_state_material).and. &
     (im.ge.1).and. &
     (im.le.num_materials)) then
  call CAV3D_STATE(xghost,t,LS,local_STATE,local_bcflag, &
          nmat,num_state_material)
  ibase=(im-1)*num_state_material
  STATE=local_STATE(ibase+istate)
  call get_primary_material(LS,im_crit)
  ibase=(im_crit-1)*num_state_material
  STATE_merge=local_STATE(ibase+istate)
 else
  print *,"istate invalid"
  stop
 endif

 return
 end subroutine CAV3D_STATE_BC

 subroutine CAV3D_HEATSOURCE( &
      im,VFRAC,time, &
      x, &
      xsten, & ! xsten(-nhalf:nhalf,SDIM)
      nhalf, &
      temp, &
      heat_source,den,CV,dt, &
      nmat)
 use probcommon_module
 IMPLICIT NONE

 integer, INTENT(in) :: nmat
 integer, INTENT(in) :: im
 real(amrex_real), INTENT(in) :: VFRAC(nmat)
 real(amrex_real), INTENT(in) :: time
 integer, INTENT(in) :: nhalf
 real(amrex_real), INTENT(in) :: x(SDIM)
 real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
 real(amrex_real), INTENT(in) :: temp(nmat)
 real(amrex_real), INTENT(in) :: den(nmat)
 real(amrex_real), INTENT(in) :: CV(nmat)
 real(amrex_real), INTENT(in) :: dt
 real(amrex_real), INTENT(out) :: heat_source

 if ((num_materials.eq.3).and.(probtype.eq.411)) then
  heat_source=zero
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine CAV3D_HEATSOURCE

end module CAV3D_module
