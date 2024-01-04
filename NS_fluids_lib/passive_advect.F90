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

! probtype==28,29,31
module passive_advect_module
use amrex_fort_module, only : amrex_real

implicit none                   

contains

  ! do any initial preparation needed
subroutine INIT_passive_advect_MODULE()
use probcommon_module
use global_utility_module
IMPLICIT NONE

 if ((probtype.eq.28).or. &
     (probtype.eq.29).or. &
     (probtype.eq.31)) then
  ! do nothing
 else
  print *,"probtype invalid"
  stop
 endif

return
end subroutine INIT_passive_advect_MODULE

 ! fluids tessellate the domain, solids are immersed. 
subroutine passive_advect_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(out) :: LS(nmat)
real(amrex_real) :: xstar,ystar,zstar
real(amrex_real) :: xprime,yprime
real(amrex_real) :: distline

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  xstar=x(1)
  ystar=x(2)
  zstar=x(SDIM)
 
  if (probtype.eq.28) then 

   if (adv_vel.ne.zero) then
    if (SDIM.eq.2) then
     if ((adv_dir.eq.1).or.(adv_dir.eq.SDIM+1)) then
      xstar=xstar-adv_vel*t
      do while (xstar.lt.problox)
       xstar=xstar+probhix-problox
      enddo
      do while (xstar.gt.probhix)
       xstar=xstar-probhix+problox
      enddo
     else if ((adv_dir.eq.2).or.(adv_dir.eq.SDIM+1)) then
      ystar=ystar-adv_vel*t
      do while (ystar.lt.probloy)
       ystar=ystar+probhiy-probloy
      enddo
      do while (ystar.gt.probhiy)
       ystar=ystar-probhiy+probloy
      enddo
     else
      print *,"adv_dir invalid probtype==28 (4)"
      stop
     endif
     zstar=ystar
    else if (SDIM.eq.3) then
     if ((adv_dir.eq.1).or.(adv_dir.eq.SDIM+1)) then
      xstar=xstar-adv_vel*t
      do while (xstar.lt.problox)
       xstar=xstar+probhix-problox
      enddo
      do while (xstar.gt.probhix)
       xstar=xstar-probhix+problox
      enddo
     else if ((adv_dir.eq.2).or.(adv_dir.eq.SDIM+1)) then
      ystar=ystar-adv_vel*t
      do while (ystar.lt.probloy)
       ystar=ystar+probhiy-probloy
      enddo
      do while (ystar.gt.probhiy)
       ystar=ystar-probhiy+probloy
      enddo
     else if ((adv_dir.eq.3).or.(adv_dir.eq.SDIM+1)) then
      zstar=zstar-adv_vel*t
      do while (zstar.lt.probloz)
       zstar=zstar+probhiz-probloz
      enddo
      do while (zstar.gt.probhiz)
       zstar=zstar-probhiz+probloz
      enddo
     else
      print *,"adv_dir invalid probtype==28 (5)"
      stop
     endif
    else
     print *,"dimension bust"
     stop
    endif
   else if (adv_vel.eq.zero) then
    ! do nothing
   else
    print *,"adv_vel is NaN"
    stop
   endif 

   if (axis_dir.eq.0) then ! dist<0 in the object
    call zalesakdist(LS(1),xstar,ystar)
    LS(2)=-LS(1)
   else if (axis_dir.eq.1) then
    xprime=(xstar-xblob)/10.0d0
    yprime=(ystar-yblob)/10.0d0+2.0d0
    call Adist(xprime,yprime,LS(1))
    LS(2)=-LS(1)
   else if (axis_dir.eq.2) then ! dist<0 in the circle
    if (SDIM.eq.2) then
     LS(1)=sqrt( (xstar-xblob)**2 + (ystar-yblob)**2 ) - radblob
    else if (SDIM.eq.3) then
     LS(1)=sqrt((xstar-xblob)**2+(ystar-yblob)**2+(zstar-zblob)**2) - radblob
    else
     print *,"dimension bust"
     stop
    endif
    LS(2)=-LS(1)
   else if (axis_dir.eq.3) then
    ! dist<0 inside the triangle.
    call triangledist(xstar,ystar,xblob,xblob2,yblob,yblob2,LS(1))
    LS(2)=-LS(1)
   else if (axis_dir.eq.4) then
    ! dist<0 inside the polygon
    call polygondist(xstar,ystar,xblob,xblob2,yblob,yblob2, &
      xblob3,yblob3,LS(1))
    LS(2)=-LS(1)
   else
    print *,"axis_dir invalid probtype=28"
    stop
   endif

  else if (probtype.eq.29) then ! single vortex

   if (SDIM.eq.2) then
    call deformdist(LS(1),xstar,ystar) ! dist<0 in the object
    LS(2)=-LS(1)

    if ((axis_dir.eq.3).or.(axis_dir.eq.4)) then

     if (denfact.eq.one) then
      ! do nothing - single vortex 2 materials
      if (num_materials.ne.2) then
       print *,"nmat invalid"
       stop
      endif
     else if (denfact.eq.-one) then ! split deforming circle in half
      if (num_materials.ne.3) then
       print *,"num_materials invalid"
       stop
      endif
      LS(3)=LS(1)  ! negative in the circle
      LS(1)=-LS(3) ! positive in the circle
      distline=half-x(1) ! positive left side, negative right side
      if (distline.lt.LS(1)) then
       LS(1)=distline ! LS(1) is negative right side of circle
      endif
      LS(2)=-LS(3) ! positive in the circle
      distline=x(1)-half !positive right side, negative left side
      if (distline.lt.LS(2)) then ! LS(2) is negative left side of circle
       LS(2)=distline             ! LS(2) is positive right side of circle
      endif
     else
      print *,"denfact invalid"
      stop
     endif
    endif

   else if (SDIM.eq.3) then

    LS(1)=sqrt((xstar-xblob)**2+(ystar-yblob)**2+(zstar-zblob)**2)-radblob
    LS(2)=-LS(1)

    if (denfact.eq.zero) then
     ! do nothing - single vortex 2 materials
     if (num_materials.ne.2) then
      print *,"nmat invalid"
      stop
     endif
    else if (denfact.eq.-one) then ! split deforming sphere in half
     if (num_materials.ne.3) then
      print *,"nmat invalid"
      stop
     endif
     LS(3)=LS(1)
     LS(1)=-LS(3)
     distline=xblob-x(1)
     if (distline.lt.LS(1)) then
      LS(1)=distline ! positive left side of circle
     endif
     LS(2)=-LS(3)
     distline=x(1)-xblob
     if (distline.lt.LS(2)) then
      LS(2)=distline  ! positive right side of circle
     endif

    else
     print *,"denfact invalid"
     stop
    endif

   else
    print *,"dimension bust"
    stop
   endif

  else if (probtype.eq.31) then ! translating circle
   if (SDIM.eq.2) then
    LS(1)=sqrt( (xstar-xblob)**2 + (ystar-yblob)**2 ) - radblob
   else if (SDIM.eq.3) then
    LS(1)=sqrt((xstar-xblob)**2+(ystar-yblob)**2+(zstar-zblob)**2)-radblob
   else
    print *,"dimension bust"
    stop
   endif
   LS(2)=-LS(1)

  else
   print *,"expecting probtype.eq. 28,29 or 31"
   stop
  endif

return
end subroutine passive_advect_LS

subroutine passive_advect_OVERRIDE_TAGFLAG( &
  i,j,k, &
  level,max_level, &
  snew_ptr,lsnew_ptr, &
  xsten,nhalf,time, &
  rflag,tagflag)
use amrex_fort_module, only : amrex_real
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: i,j,k
integer, INTENT(in) :: level,max_level
integer, INTENT(in) :: nhalf
real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) :: rflag
integer, INTENT(inout) :: tagflag
real(amrex_real), INTENT(in),pointer :: snew_ptr(D_DECL(:,:,:),:)
real(amrex_real), INTENT(in),pointer :: lsnew_ptr(D_DECL(:,:,:),:)
real(amrex_real), dimension(3) :: local_x
real(amrex_real), dimension(SDIM) :: local_delta
real(amrex_real) :: LS(D_DECL(-1:1,-1:1,-1:1))
integer :: dir
integer :: i1,j1,k1
real(amrex_real) :: curv
real(amrex_real) :: curv_cutoff

if (nhalf.lt.3) then
 print *,"nhalf invalid passive advect override tagflag"
 stop
endif
if ((level.ge.0).and.(level.lt.max_level)) then
 ! do nothing
else
 print *,"level and/or max_level invalid"
 print *,"level=",level
 print *,"max_level=",max_level
 stop
endif
do dir=1,SDIM
 local_x(dir)=xsten(0,dir)
enddo
local_x(3)=xsten(0,SDIM)
do dir=1,SDIM
 local_delta(dir)=xsten(1,dir)-xsten(-1,dir)
 if (local_delta(dir).gt.zero) then
  ! do nothing
 else
  print *,"local_delta invalid passive_advect_override_tagflag"
  stop
 endif
enddo !dir=1..sdim

if ((probtype.eq.28).and. & ! zalesak
    (axis_dir.eq.0)) then

 if (level.lt.max_level-1) then
  !do nothing
 else if (level.eq.max_level-1) then
  rflag=0.0d0
  tagflag=0
  if (abs(lsnew_ptr(D_DECL(i,j,k),1)).le.local_delta(1)) then
   do i1=-1,1
   do j1=-1,1
   do k1=-1,1
    LS(D_DECL(i1,j1,k1))=lsnew_ptr(D_DECL(i+i1,j+j1,k+k1),1)
   enddo
   enddo
   enddo
   curv_cutoff=one/(four*local_delta(1))
   call curverr(curv,LS,xsten,nhalf)
   if (abs(curv).lt.curv_cutoff) then
    !do nothing
   else if (abs(curv).ge.curv_cutoff) then
    rflag=1.0d0
    tagflag=1
   else
    print *,"curv=",curv
    print *,"curv_cutoff ",curv_cutoff
    stop
   endif
 
  endif
 else
  print *,"level invalid"
  stop
 endif

endif

end subroutine passive_advect_OVERRIDE_TAGFLAG


! initial velocity is some kind of shear flow
subroutine passive_advect_VEL(xvec,time,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xvec(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: VEL(SDIM)
integer dir
integer, INTENT(in) :: velsolid_flag
real(amrex_real) :: x,y,z

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

if ((velsolid_flag.eq.0).or. &
    (velsolid_flag.eq.1)) then
 ! do nothing
else 
 print *,"velsolid_flag invalid"
 stop
endif

x=xvec(1)
y=xvec(2)
z=xvec(SDIM)

do dir=1,SDIM
 if (dx(dir).gt.zero) then
  ! do nothing
 else
  print *,"dx invalid"
  stop
 endif
enddo

do dir=1,SDIM
 VEL(dir)=zero
enddo

do dir=1,SDIM

 if (SDIM.eq.2) then

  if (probtype.eq.28) then
   if (dir.eq.1) then
    call zalesakuu(VEL(dir),x,y,z,time,dx)
   else if (dir.eq.2) then
    call zalesakvv(VEL(dir),x,y,z,time,dx)
   else
    print *,"dir invalid"
    stop
   endif
  else if (probtype.eq.29) then
   if (dir.eq.1) then
    call deformuu(VEL(dir),x,y,time,dx)
   else if (dir.eq.2) then
    call deformvv(VEL(dir),x,y,time,dx)
   else
    print *,"dir invalid"
    stop
   endif
  else if (probtype.eq.31) then 
   if (dir.eq.1) then
    call circleuu(VEL(dir),x,y,y)
   else if (dir.eq.2) then
    call circlevv(VEL(dir),x,y,y)
   else
    print *,"dir invalid"
    stop
   endif
  else
   print *,"probtype invalid"
   stop
  endif
 else if (SDIM.eq.3) then

  if (probtype.eq.28) then
   if (dir.eq.1) then
    call zalesakuu(VEL(dir),x,y,z,time,dx)
   else if (dir.eq.2) then
    call zalesakvv(VEL(dir),x,y,z,time,dx)
   else if (dir.eq.3) then
    call zalesakww(VEL(dir),x,y,z,time,dx)
   else
    print *,"dir invalid"
    stop
   endif
  else if (probtype.eq.29) then
   if (dir.eq.1) then
    call deform3duu(VEL(dir),x,y,z,time,dx)
   else if (dir.eq.2) then
    call deform3dvv(VEL(dir),x,y,z,time,dx)
   else if (dir.eq.3) then
    call deform3dww(VEL(dir),x,y,z,time,dx)
   else
    print *,"dir invalid"
    stop
   endif
  else if (probtype.eq.31) then
   if (dir.eq.1) then
    call circleuu(VEL(dir),x,y,z)
   else if (dir.eq.2) then
    call circlevv(VEL(dir),x,y,z)
   else if (dir.eq.3) then
    call circleww(VEL(dir),x,y,z)
   else
    print *,"dir invalid"
    stop
   endif
  else
   print *,"probtype invalid"
   stop
  endif

 else
  print *,"dimension bust"
  stop
 endif

enddo ! dir=1..sdim

return 
end subroutine passive_advect_VEL


! this routine used as a default when
! pressure boundary conditions are prescribed.
! For the case when only top wall is 
! "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine passive_advect_PRES(x,t,LS,PRES,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: PRES

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif


 PRES=zero

return 
end subroutine passive_advect_PRES



subroutine passive_advect_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
integer, INTENT(in) :: nmat
integer, INTENT(in) :: nstate_mat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: STATE(nmat*nstate_mat)
integer im,ibase,n

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if (nstate_mat.eq.num_state_material) then
 ! do nothing
else
 print *,"nstate_mat invalid"
 stop
endif

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

   ! initial species in inputs?
  do n=1,num_species_var
   STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo

 enddo ! im=1..num_materials

return
end subroutine passive_advect_STATE

 ! dir=1..sdim  side=1..2
subroutine passive_advect_LS_BC(xwall,xghost,t,LS, &
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
real(amrex_real), INTENT(in) ::  dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call passive_advect_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine passive_advect_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine passive_advect_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
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

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call passive_advect_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine passive_advect_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine passive_advect_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
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

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

 call passive_advect_PRES(xghost,t,LS,PRES,nmat)

return
end subroutine passive_advect_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine passive_advect_STATE_BC(xwall,xghost,t,LS, &
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
 call passive_advect_STATE(xghost,t,LS,local_STATE,local_bcflag, &
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
end subroutine passive_advect_STATE_BC

subroutine passive_advect_clamped_LS(x,t,LS,vel,temperature,prescribed_flag,dx)
use probcommon_module
use global_utility_module
IMPLICIT NONE

  real(amrex_real), INTENT(in) :: x(SDIM)
  real(amrex_real), INTENT(in) :: dx(SDIM)
  real(amrex_real), INTENT(in) :: t
  real(amrex_real), INTENT(out) :: LS
  real(amrex_real), INTENT(out) :: vel(SDIM)
  real(amrex_real), INTENT(out) :: temperature
  integer, INTENT(out) :: prescribed_flag
  integer :: velsolid_flag
  real(amrex_real) :: LSarray(num_materials)

if ((probtype.eq.28).or. &
    (probtype.eq.29).or. &
    (probtype.eq.31)) then
 LS=CLAMPED_EVERYWHERE_LS
 velsolid_flag=0
 call passive_advect_VEL(x,t,LSarray,vel,velsolid_flag,dx,num_materials)
 temperature=fort_tempconst(1)
 prescribed_flag=1
else
 print *,"probtype invalid"
 stop
endif

return
end subroutine passive_advect_clamped_LS


subroutine passive_advect_CFL_HELPER(time,dir,uu,dx)
use probcommon_module
IMPLICIT NONE
integer, INTENT(in) :: dir
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) :: uu
real(amrex_real), INTENT(in) :: dx(SDIM)

if ((dir.lt.0).or.(dir.ge.SDIM)) then
 print *,"dir invalid"
 stop
endif

if (dir.eq.0) then
 ! do nothing
else if (dir.eq.1) then
 ! do nothing
else if ((dir.eq.2).and.(SDIM.eq.3)) then
 ! do nothing
else
 print *,"dir invalid passitve_advect_CFL_HELPER"
 stop
endif
if (probtype.eq.31) then
 uu=max(abs(uu),abs(adv_vel))
else if (probtype.eq.28) then
 ! do nothing
else if (probtype.eq.29) then
 ! do nothing
else
 print *,"probtype invalid"
 stop
endif

return
end subroutine passive_advect_CFL_HELPER


end module passive_advect_module
