#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

#include "AMReX_ArrayLim.H"


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

implicit none                   

contains

   ! do any initial preparation needed
 subroutine INIT_CAV3D_MODULE()
 IMPLICIT NONE

 return
 end subroutine INIT_CAV3D_MODULE

 subroutine OPEN_CAV3D_CASFILE(part_id,unit_id)
 IMPLICIT NONE

 INTEGER_T part_id,unit_id
 character(40) :: casname

 if (part_id.eq.1) then
  casname="nozzle.cas"
 else
  print *,"part_id invalid"
  stop
 endif

 print *,"opening ",casname
 OPEN(unit=unit_id,file=casname,access='sequential', &
      form="formatted",status='old')

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

 REAL_T nodes(3,3) ! dir,nodenum
 INTEGER_T nodemap(3)
 REAL_T a(3),b(3),c(3)
 INTEGER_T inode
 INTEGER_T dir
 INTEGER_T dircrit
 REAL_T mag
 REAL_T xcrit
 REAL_T ncrit
 INTEGER_T swap

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

 return
 end subroutine CAV3D_ORDER_NODES

 subroutine CAV3D_SLICE(xmap3D,xslice3D,problo3D,probhi3D,dx_slice)
 use probcommon_module
 IMPLICIT NONE

 REAL_T dx_slice
 INTEGER_T xmap3D(3)
 REAL_T xslice3D(3)
 REAL_T problo3D(3)
 REAL_T probhi3D(3)

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

 return
 end subroutine CAV3D_SLICE

  ! fluids tessellate the domain, solids are immersed. 
 subroutine CAV3D_LS(x,t,LS)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 INTEGER_T im
 REAL_T LS(num_materials)

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

 subroutine CAV3D_VEL(x,t,LS,VEL,velsolid_flag)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T VEL(SDIM)
 INTEGER_T dir
 INTEGER_T velsolid_flag

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

 subroutine CAV3D_PRES(x,t,LS,PRES)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T PRES
 REAL_T ymid

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


 subroutine CAV3D_STATE(x,t,LS,STATE)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T STATE(num_materials*num_state_material)
 INTEGER_T im,ibase,n

 if ((num_materials.eq.3).and. &
     (num_state_material.ge.2).and. &
     (probtype.eq.411)) then
  do im=1,num_materials
   ibase=(im-1)*num_state_material
   STATE(ibase+1)=fort_denconst(im)
   if (t.eq.zero) then
    STATE(ibase+2)=fort_initial_temperature(im)
   else if (t.gt.zero) then
    STATE(ibase+2)=fort_tempconst(im)
   else
    print *,"t invalid"
    stop
   endif
   do n=1,num_species_var
    STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
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
    LS_in,dir,side,dx)
 use probcommon_module
 IMPLICIT NONE

 REAL_T xwall
 REAL_T xghost(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T LS_in(num_materials)
 INTEGER_T dir,side
 REAL_T dx(SDIM)

 if ((dir.ge.1).and.(dir.le.SDIM).and. &
     (side.ge.1).and.(side.le.2)) then
  call CAV3D_LS(xghost,t,LS)
 else
  print *,"dir or side invalid"
  stop
 endif
 
 return
 end subroutine CAV3D_LS_BC


  ! dir=1..sdim  side=1..2 veldir=1..sdim
 subroutine CAV3D_VEL_BC(xwall,xghost,t,LS, &
    VEL,VEL_in,veldir,dir,side,dx)
 use probcommon_module
 IMPLICIT NONE

 REAL_T xwall
 REAL_T xghost(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T VEL
 REAL_T VEL_in
 INTEGER_T veldir,dir,side
 REAL_T dx(SDIM)
 REAL_T local_VEL(SDIM)
 INTEGER_T velsolid_flag

 velsolid_flag=0
 if ((dir.ge.1).and.(dir.le.SDIM).and. &
     (side.ge.1).and.(side.le.2).and. &
     (veldir.ge.1).and.(veldir.le.SDIM)) then

  call CAV3D_VEL(xghost,t,LS,local_VEL,velsolid_flag)
  VEL=local_VEL(veldir)

 else
  print *,"dir,side, or veldir invalid"
  stop
 endif

 return
 end subroutine CAV3D_VEL_BC


  ! dir=1..sdim  side=1..2
 subroutine CAV3D_PRES_BC(xwall,xghost,t,LS, &
    PRES,PRES_in,dir,side,dx)
 use probcommon_module
 IMPLICIT NONE

 REAL_T xwall
 REAL_T xghost(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T PRES
 REAL_T PRES_in
 INTEGER_T dir,side
 REAL_T dx(SDIM)

 if ((dir.ge.1).and.(dir.le.SDIM).and. &
     (side.ge.1).and.(side.le.2)) then

  call CAV3D_PRES(xghost,t,LS,PRES)

 else
  print *,"dir or side invalid"
  stop
 endif

 return
 end subroutine CAV3D_PRES_BC

  ! dir=1..sdim  side=1..2
 subroutine CAV3D_STATE_BC(xwall,xghost,t,LS, &
    STATE,STATE_merge,STATE_in,im,istate,dir,side,dx)
 use probcommon_module
 IMPLICIT NONE

 REAL_T xwall
 REAL_T xghost(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T local_STATE(num_materials*num_state_material)
 REAL_T STATE
 REAL_T STATE_merge
 REAL_T STATE_in
 INTEGER_T dir,side
 REAL_T dx(SDIM)
 INTEGER_T istate,im
 INTEGER_T ibase,im_crit,im_loop

 if ((istate.ge.1).and. &
     (istate.le.num_state_material).and. &
     (im.ge.1).and. &
     (im.le.num_materials)) then
  call CAV3D_STATE(xghost,t,LS,local_STATE)
  ibase=(im-1)*num_state_material
  STATE=local_STATE(ibase+istate)
  im_crit=1
  do im_loop=2,num_materials
   if (LS(im_loop).gt.LS(im_crit)) then
    im_crit=im_loop
   endif
  enddo
  ibase=(im_crit-1)*num_state_material
  STATE_merge=local_STATE(ibase+istate)
 else
  print *,"istate invalid"
  stop
 endif

 return
 end subroutine CAV3D_STATE_BC

 subroutine CAV3D_HEATSOURCE(im,VFRAC,time,x,temp, &
      heat_source,den,CV,dt)
 use probcommon_module
 IMPLICIT NONE

 INTEGER_T im
 REAL_T VFRAC(num_materials)
 REAL_T time
 REAL_T x(SDIM)
 REAL_T temp(num_materials)
 REAL_T den(num_materials)
 REAL_T CV(num_materials)
 REAL_T dt
 REAL_T heat_source

 if ((num_materials.eq.3).and.(probtype.eq.411)) then
  heat_source=zero
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine CAV3D_HEATSOURCE

end module CAV3D_module
