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

! probtype==55 
module GENERAL_PHASE_CHANGE_module

implicit none                   

REAL_T :: DEF_VAPOR_GAMMA

contains

  ! do any initial preparation needed
subroutine INIT_GENERAL_PHASE_CHANGE_MODULE()
IMPLICIT NONE

  DEF_VAPOR_GAMMA =  1.666666667D0

return
end subroutine INIT_GENERAL_PHASE_CHANGE_MODULE


 ! this is velocity boundary condition at the top of the domain.  
subroutine acoustic_pulse_bc(time,vel_pulse,x_vec,for_dt)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: for_dt
REAL_T, intent(in) :: time
REAL_T, intent(out) :: vel_pulse
REAL_T, intent(in) :: x_vec(SDIM)

if ((time.ge.zero).and.(time.le.1.0e+20)) then
 ! do nothing
else
 print *,"time invalid"
 stop
endif

if (abs(x_vec(1))+abs(x_vec(2))+abs(x_vec(SDIM)).le.1.0e+20) then
 ! do nothing
else
 print *,"x,y, or z invalid"
 stop
endif
if ((for_dt.eq.0).or.(for_dt.eq.1)) then
 ! do nothing
else
 print *,"for_dt invalid"
 stop
endif

if (probtype.eq.55) then

 if (axis_dir.eq.7) then

  vel_pulse=zero
 
  if (1.eq.1) then 

   if (for_dt.eq.1) then
    vel_pulse=200.0
   else if (for_dt.eq.0) then 
    if ((time.ge.5.0e-8).and.(time.le.5.0e-7)) then
     vel_pulse=-200.0
    endif
   else
    print *,"for_dt invalid"
    stop
   endif

  endif

 else
  print *,"unexpected axis_dir"
  stop
 endif

else
 print *,"unexpected probtype"
 stop
endif

return
end subroutine acoustic_pulse_bc 

subroutine GENERAL_PHASE_CHANGE_CFL_HELPER(time,dir,uu,dx)
use probcommon_module
implicit none
INTEGER_T, intent(in) :: dir
REAL_T, intent(in) :: time
REAL_T, intent(inout) :: uu
REAL_T, intent(in) :: dx(SDIM)

INTEGER_T dir2
REAL_T utest
REAL_T xvec_dummy(SDIM)
INTEGER_T for_dt


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
 print *,"dir invalid GENERAL_PHASE_CHANGE_CFL_HELPER"
 stop
endif

if (probtype.eq.55) then
 if (axis_dir.eq.7) then
  do dir2=1,SDIM
   xvec_dummy(dir2)=zero
  enddo
  for_dt=1
  call acoustic_pulse_bc(time,utest,xvec_dummy,for_dt)
  uu=max(abs(uu),abs(utest))
 endif
else
 print *,"unexpected probtype"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_CFL_HELPER

! dist>0 in the fluid
subroutine GENERAL_soliddist(x,dist,im) 
use probcommon_module
use global_utility_module
implicit none
REAL_T, intent(in), dimension(SDIM) :: x !spatial coordinates
INTEGER_T, intent(in) :: im
REAL_T, intent(out) :: dist
INTEGER_T :: nmat

nmat=num_materials
if (nmat.lt.1) then
 print *,"nmat invalid in soliddist"
 stop
endif

if ((im.lt.1).or.(im.gt.nmat)) then
 print *,"im invalid11"
 stop
endif

if (FSI_flag(im).eq.1) then ! prescribed solid (EUL)
 ! do nothing
else
 print *,"FSI_flag(im) invalid"
 stop
endif

dist=99999.0

! GENERAL_soliddist: dist>0 in fluid 2d or 3d
if (probtype.eq.55) then 

 if ((axis_dir.eq.0).or. &
     (axis_dir.eq.5).or. &  ! freezing
     (axis_dir.eq.6).or. &  ! boiling (incomp)
     (axis_dir.eq.7).or. &  ! boiling (comp)
     (axis_dir.eq.1)) then

  if ((radblob5.gt.zero).and. &
      (axis_dir.eq.0)) then ! ellipse
 
! Professor Yongsheng Lian was here: 
   if (SDIM.eq.2) then
    dist=sqrt((x(1)-xblob2)*(x(1)-xblob2)/(radblob3**2)+ &
        (x(2)-yblob2)*(x(2)-yblob2)/(radblob4**2))-radblob5
   else
    dist=sqrt((x(1)-xblob2)*(x(1)-xblob2)/(radblob3**2)+ &
        (x(2)-yblob2)*(x(2)-yblob2)/(radblob4**2)+ &
        (x(SDIM)-zblob2)*(x(SDIM)-zblob2)/(radblob9**2))-radblob5
   endif

  else if ((axis_dir.eq.1).and.(nmat.eq.3)) then
   ! do nothing: solid replaced by ice.

   ! axis_dir=6: boiling sites
   ! axis_dir=1: falling drop on substrate and then freezing.
  else if ((radblob5.eq.zero).or. &
           (axis_dir.eq.5).or. &
           (axis_dir.eq.6).or. &
           ((axis_dir.eq.1).and. &
            (nmat.eq.4))) then

    ! dist>0 in the substrate
   call ice_substrate_distance(x(1),x(2),x(SDIM),dist)
    ! now make dist<0 in the substrate.
   dist=-dist

  endif

 else 
  print *,"axis_dir invalid"
  stop
 endif

else
 print *,"expecting probtype.eq.55"
 stop
endif

end subroutine GENERAL_soliddist

subroutine GENERAL_PHASE_CHANGE_check_vel_rigid(x,t,vel,dir)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: vel
INTEGER_T, intent(in) :: dir

if (t.ge.0.0d0) then
 ! do nothing
else
 print *,"t invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM)) then
 ! do nothing
else
 print *,"dir invalid"
 stop
endif

if (probtype.eq.55) then
 if (vel.eq.0.0d0) then
  ! do nothing
 else if (abs(vel).le.1.0D+20) then
  ! do nothing
 else
  print *,"GENERAL_PHASE_CHANGE_check_vel_rigid: vel not expected"
  print *,"x,y,z ",x(1),x(2),x(SDIM)
  print *,"t ",t
  print *,"dir,vel ",dir,vel
  stop
 endif
else
 print *,"probtype invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_check_vel_rigid


! ice + water thickness = radblob
! water thickness = radblob3  usually radblob3 << radlob (initial water
! is "seed" for starting the melting process)
!   ---------
!   | water |
!   ---------
!   | ice   |
!-------------------
!   substrate

 ! fluids tessellate the domain, solids are immersed. 
subroutine GENERAL_PHASE_CHANGE_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(out) :: LS(nmat)
REAL_T :: dist_gas,dist_liquid,dist_ice,dist_liq2,distsolid
INTEGER_T :: im
INTEGER_T :: im_solid_materialdist
REAL_T :: initial_time

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  im_solid_materialdist=im_solid_primary()
  initial_time=zero

  if (probtype.eq.55) then

   do im=1,nmat
    if (FSI_flag(im).eq.1) then
     call GENERAL_soliddist(x,LS(im),nmat)  ! returns LS<0 in solid
     LS(im)=-LS(im)   ! now LS>0 in solid
     distsolid=LS(im)
    endif
   enddo

   ! drop on slope problem (2d or 3d)
   ! radblob4 is used as a "switch" in order to specify static
   ! solution at t=0.
   ! radblob5 is switch for drop hitting ellipse problem.
   if (axis_dir.eq.3) then

    if (SDIM.eq.2) then
     LS(1)=x(SDIM)-(zblob+radblob*cos(two*Pi*x(1)/xblob))
    else if (SDIM.eq.3) then
     LS(1)=x(SDIM)-(zblob+radblob*cos(two*Pi*x(1)/xblob)*cos(two*Pi*x(2)/yblob))
    else
     print *,"dimension bust"
     stop
    endif
    LS(2)=-LS(1)

   else if (axis_dir.eq.5) then

     ! "drop_slope_dist" declared in GLOBALUTIL.F90
     ! in: GENERAL_PHASE_CHANGE_LS (initial angle=static angle)
     ! maxtall==two*radblob => no ice in this call.
    call drop_slope_dist(x(1),x(2),x(SDIM),initial_time,nmat, &
      two*radblob,dist_ice,dist_liquid)

    dist_gas=-dist_liquid
    LS(1)=dist_liquid
    LS(2)=dist_gas

    ! in materialdistbatch:
    ! nucleate boiling 2D or 3D
    ! Sato and Niceno problem
   else if ((axis_dir.eq.6).or. &
            (axis_dir.eq.7)) then 
    if (n_sites.gt.0) then
     if (nucleation_init_time.eq.zero) then
       ! in: GLOBALUTIL.F90; negative if x in a bubble.
      call nucleation_sites(x,dist_liquid,pos_sites)
     else if (nucleation_init_time.gt.zero) then
      dist_liquid=9999.0
     else
      print *,"nucleation_init_time invalid"
      stop
     endif
    else if (n_sites.eq.0) then
     ! do nothing
    else
     print *,"n_sites invalid [4]",n_sites,nucleation_init_time
     stop
    endif

    dist_gas=-dist_liquid
    LS(1)=dist_liquid
    LS(2)=dist_gas
     ! thickness of initial gas layer at outflow
    if (radblob10.gt.zero) then
     if ((nmat.eq.4).and.(im_solid_materialdist.eq.nmat)) then
      if (gravity_dir.eq.1) then
       if (radblob10.lt.problenx) then 
        LS(3)=x(gravity_dir)-(probhix-radblob10)
       else
        print *,"radblob10 invalid"
        stop
       endif
      else if (gravity_dir.eq.2) then
       if (radblob10.lt.probleny) then 
        LS(3)=x(gravity_dir)-(probhiy-radblob10)
       else
        print *,"radblob10 invalid"
        stop
       endif
      else if ((gravity_dir.eq.3).and.(SDIM.eq.3)) then
       if (radblob10.lt.problenz) then 
        LS(3)=x(gravity_dir)-(probhiz-radblob10)
       else
        print *,"radblob10 invalid"
        stop
       endif
      else
       print *,"gravity_dir invalid"
       stop
      endif
      if (LS(3).lt.zero) then
       LS(1)=min(LS(1),-LS(3))
      else if (LS(3).ge.zero) then
       if ((LS(2).lt.zero).and.(LS(1).gt.zero)) then
        LS(1)=-LS(3)
       else 
        print *,"LS(2) or LS(1) invalid"
        stop
       endif
      else
       print *,"LS(3) invalid"
       stop
      endif
     else
      print *,"nmat or im_solid_materialdist invalid"
      print *,"nmat=",nmat
      print *,"im_solid_materialdist=",im_solid_materialdist
      stop
     endif
    else if (radblob10.eq.zero) then
     ! do nothing
    else
     print *,"radblob10 invalid"
     stop
    endif

   else if ((axis_dir.eq.0).or. &
            (axis_dir.eq.1)) then

    if ((radblob6.gt.zero).and.(radblob7.gt.zero)) then
     print *,"cannot have both radblob6 and radblob7 positive"
     stop
    endif

    if (radblob3.gt.zero) then
      ! negative on the inside of the square
     call squaredist(x(1),x(2),xblob-radblob,xblob+radblob,yblob, &
      yblob+radblob3,dist_liquid)
     dist_liquid=-dist_liquid
    else if (radblob3.eq.zero) then
     if (SDIM.eq.2) then
      dist_liquid=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)
     else
      dist_liquid=radblob- &
          sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)
     endif
     if (radblob5.gt.zero) then
      if (SDIM.eq.2) then
       dist_liq2=radblob5-sqrt((x(1)-xblob5)**2+(x(2)-yblob5)**2)
      else
       dist_liq2=radblob5- &
         sqrt((x(1)-xblob5)**2+(x(2)-yblob5)**2+(x(SDIM)-zblob5)**2)
      endif

      if (dist_liq2.gt.dist_liquid) then
       dist_liquid=dist_liq2
      endif
     endif 
    else
     print *,"radblob3 invalid"
     stop
    endif

    if (radblob4.eq.zero) then
     ! do nothing
    else if (radblob4.eq.one) then
     if ((radblob6.ne.zero).or.(radblob7.ne.zero).or. &
         (radblob5.ne.zero)) then
      print *,"conflicting parameters probtype=",probtype
      stop
     endif
     ! in: materialdistbatch (initial angle=static angle)
     call drop_slope_dist(x(1),x(2),x(SDIM),initial_time,nmat, &
      two*radblob,dist_ice,dist_liquid)
    else
     print *,"radblob4 invalid radblob4=",radblob4
     stop
    endif

    if (radblob6.gt.zero) then
     if (SDIM.eq.2) then
      dist_liq2=radblob6-sqrt((x(1)-xblob6)**2+(x(2)-yblob6)**2)
     else if (SDIM.eq.3) then
      dist_liq2=radblob6- &
           sqrt((x(1)-xblob6)**2+(x(2)-yblob6)**2+(x(SDIM)-zblob6)**2)
     endif
     if (dist_liq2.gt.dist_liquid) then
      dist_liquid=dist_liq2
     endif
    endif

    dist_gas=-dist_liquid
    LS(1)=dist_liquid
    LS(2)=dist_gas
   else
    print *,"axis_dir invalid probtype=55 materialdistbatch"
    stop
   endif

   if (axis_dir.eq.0) then

    if (radblob7.gt.zero) then ! drop collision?
     if (nmat.lt.3) then
      print *,"nmat invalid"
      stop
     endif
     if (radblob6.gt.zero) then
      print *,"cannot have both radblob6 and radblob7 positive"
      stop
     endif

     LS(3)=radblob7-sqrt( (x(1)-xblob7)**2+(x(2)-yblob7)**2 )  ! pos. in drop
     dist_gas=-dist_liquid
     if (dist_gas.gt.-LS(3)) then
      dist_gas=-LS(3)
     endif
     LS(1)=dist_liquid
     LS(2)=dist_gas

    endif  ! radblob7>0
    !ICE MEHDI
   else if (axis_dir.eq.1) then  ! drop falling on substrate
    if (nmat.lt.3) then
     print *,"nmat invalid"
     stop
    endif
 
    ! positive in the "ice" or "substrate"
    ! materialdistbatch
    call ice_substrate_distance(x(1),x(2),x(SDIM),distsolid)

    if (nmat.eq.4) then
     if (im_solid_materialdist.ne.nmat) then
      print *,"expecting im_solid_materialdist=nmat"
      stop
     endif
    else if (nmat.eq.3) then
     if (im_solid_materialdist.ne.0) then
      print *,"expecting im_solid_materialdist=0"
      stop
     endif
    else
     print *,"nmat invalid"
     stop
    endif

    LS(nmat)=distsolid
    if (is_rigid(nmat,nmat).ne.1) then
     print *,"expecting last material to be rigid"
     stop
    endif
    dist_gas=-dist_liquid
    LS(1)=dist_liquid
    LS(2)=dist_gas

    ! ICE MEHDI: compare with freezing singularity paper
   else if (axis_dir.eq.5) then
    if (nmat.ne.4) then
     print *,"nmat invalid"
     stop
    endif
    if (im_solid_materialdist.ne.4) then
     print *,"expecting im_solid_materialdist=4"
     stop
    endif
    ! material 3: ice
    ! material 1: water
    ! in: materialdist_batch (initial angle=static angle)
    call drop_slope_dist(x(1),x(2),x(SDIM),initial_time,nmat,radblob3, &
      dist_ice,dist_liquid)
    if (is_rigid(nmat,3).ne.0) then
     print *,"expecting material 3 to be ice"
     stop
    endif
    LS(1)=dist_liquid
    LS(3)=dist_ice
    if (LS(2).gt.-distsolid) then
     LS(2)=-distsolid
    endif
    if ((LS(3).le.zero).and.(distsolid.ge.zero)) then
     LS(3)=distsolid
    endif

    if ((LS(1).lt.zero).and. &
        (LS(2).lt.zero).and. &
        (LS(3).lt.zero)) then
     print *,"fluids should tessellate"
     print *,"probtype,axis_dir ",probtype,axis_dir
     print *,"x,y,z ",x(1),x(2),x(SDIM)
     print *,"LS(1) ",LS(1)
     print *,"LS(2) ",LS(2)
     print *,"LS(3) ",LS(3)
     print *,"LS(4) ",LS(4)
     stop
    endif
   else if ((axis_dir.eq.6).or. &  ! incompressible boiling
            (axis_dir.eq.7)) then  ! compressible boiling
    ! do nothing (inputs.boiling) (boiling sites)
    !  (or Sato and Niceno problem)
    !  (or Tryggvason problem)
   else
    print *,"axis_dir invalid for drop falling on ice"
    stop
   endif
  else
   print *,"expecting probtype.eq.55"
   stop
  endif

return
end subroutine GENERAL_PHASE_CHANGE_LS

! initial velocity is zero
subroutine GENERAL_PHASE_CHANGE_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, intent(in) :: velsolid_flag
REAL_T :: temp
REAL_T :: xmid,zmid

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

! in "initvelocity":
! drop on slope 2d or 3d
! if advbot<>0, then we have liquid sphere falling onto ramp
! if axis_dir=1 then we have liquid sphere falling onto ice.
! liquid is material 1, gas is material 2, solid/ice is material 3
if (probtype.eq.55) then

 do dir=1,SDIM
  VEL(dir)=zero
 enddo

  ! in: GLOBALUTIL.F90
 call default_rampvel(t,VEL(1),VEL(2),VEL(SDIM))

 if (axis_dir.eq.5) then
  do dir=1,SDIM
   VEL(dir)=zero
  enddo
  ! axis_dir=6,7 is boiling (e.g. Sato and Niceno)
 else if ((axis_dir.eq.6).or. & ! incompressible
          (axis_dir.eq.7)) then ! compressible
  if (yblob10.gt.zero) then
   if((x(2).ge.yblob2).and.(x(2).le.yblob10)) then
    ! Distance from substrate
    temp = x(2)-yblob2
    VEL(1)=VEL(1)*(1.5d0*temp/yblob10 -0.5*(temp/yblob10)**3)
   end if
  else if (yblob10.eq.zero) then
   VEL(1)=zero
  else
   print *,"yblob10 invalid"
   stop
  end if
  VEL(2)=zero
  VEL(SDIM)=zero

 else if ((axis_dir.eq.0).or. &
          (axis_dir.eq.1)) then
  do dir=1,SDIM
   VEL(dir)=zero
  enddo

  if (advbot.ne.zero) then
   if ((radblob6.ne.zero).or.(radblob7.ne.zero)) then
    print *,"parameters conflict"
    stop
   endif
   if (LS(1).gt.-dx(1)) then
    VEL(SDIM)=-abs(advbot)
   endif
  endif ! advbot <> 0

  if (radblob5.gt.zero) then  ! impact droplet on ellipse
   if (adv_dir.eq.1) then
    VEL(1)=adv_vel
   else if (adv_dir.eq.2) then
    VEL(2)=adv_vel
   else if ((adv_dir.eq.3).and.(SDIM.eq.3)) then
    VEL(SDIM)=adv_vel
   else
    print *,"adv_vel invalid"
    stop
   endif 
  else if (radblob5.lt.zero) then
   print *,"radblob5 invalid"
   stop
  endif
  if ((radblob6.gt.zero).or.(radblob7.gt.zero)) then
   if ((radblob6.gt.zero).and.(radblob7.gt.zero)) then
    print *,"cannot have both radblob6 and radblob7 positive"
    stop
   endif
   if (radblob6.gt.zero) then
    if (LS(1).gt.-dx(1)) then
     if (levelrz.eq.1) then
      zmid=half*(yblob6+yblob)
      if (x(2).lt.zmid) then
       VEL(2)=abs(advbot)
      else
       VEL(2)=-abs(vinletgas)
      endif
     else if (levelrz.eq.0) then
      xmid=half*(xblob6+xblob)
      if (x(1).lt.xmid) then
       VEL(1)=abs(advbot)
      else
       VEL(1)=-abs(vinletgas)
      endif
     else
      print *,"levelrz invalid init velocity"
      stop
     endif
    endif  ! liquid
   else if (radblob7.gt.zero) then
    if (LS(1).gt.-dx(1)) then
     if (levelrz.eq.1) then
      VEL(2)=abs(advbot)
     else if (levelrz.eq.0) then
      VEL(1)=abs(advbot)
     else
      print *,"levelrz invalid init velocity 2"
      stop
     endif
    else if (LS(3).gt.-dx(1)) then
     if (levelrz.eq.1) then
      VEL(2)=-abs(vinletgas)
     else if (levelrz.eq.0) then
      VEL(1)=-abs(vinletgas)
     else
      print *,"levelrz invalid probtype 55"
      stop
     endif
    endif
   else
    print *,"bust"
    stop
   endif
  endif ! drop collision (radblob6 or radblob7 > 0)
 else
  print *,"axis_dir invalid"
  stop
 endif
else
 print *,"expecting probtype==55"
 stop
endif

return 
end subroutine GENERAL_PHASE_CHANGE_VEL


! this routine used as a default when
! pressure boundary conditions are prescribed.
! For the case when only top wall is 
! "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine GENERAL_PHASE_CHANGE_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: PRES

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
PRES=zero

return 
end subroutine GENERAL_PHASE_CHANGE_PRES



subroutine GENERAL_PHASE_CHANGE_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, intent(in) :: nmat
INTEGER_T, intent(in) :: nstate_mat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: STATE(nmat*nstate_mat)
INTEGER_T im,ibase,n
REAL_T water_temp

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
if (probtype.eq.55) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im) ! density prescribed in the inputs file.
  if (t.eq.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im) !initial temperature in inputs
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

   ! initial temperature for boiling or freezing
   ! nucleate boiling: Sato and Niceno or Tryggvason
  if ((axis_dir.eq.6).or. &  ! incompressible
      (axis_dir.eq.7)) then  ! compressible
   ! water phase
   if (im.eq.1) then
    ! bcflag=0 (calling from FORT_INITDATA)
    call outside_temperature(t,x(1),x(2),x(SDIM),water_temp,im,0)
    STATE(ibase+ENUM_TEMPERATUREVAR+1)=water_temp  
   endif ! im=1
  else if (axis_dir.eq.5) then ! freezing drop on substrate
   if (nmat.lt.4) then
    print *,"nmat too small for freezing drop on substrate"
    stop
   endif
   ! ice  or substrate (initial temperature)
   if ((im.eq.3).or.(im.eq.4)) then
    ! bcflag=0 (calling from FORT_INITDATA)
    call outside_temperature(t,x(1),x(2),x(SDIM),water_temp,im,0)
    STATE(ibase+ENUM_TEMPERATUREVAR+1)=water_temp  
   endif
  endif

 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine GENERAL_PHASE_CHANGE_STATE

 ! dir=1..sdim  side=1..2
subroutine GENERAL_PHASE_CHANGE_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(inout) :: LS(nmat)
REAL_T, intent(in) :: LS_in(nmat)
INTEGER_T, intent(in) :: dir,side
REAL_T, intent(in) ::  dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call GENERAL_PHASE_CHANGE_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine GENERAL_PHASE_CHANGE_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(inout) :: VEL
REAL_T, intent(in) :: VEL_in
INTEGER_T, intent(in) :: veldir,dir,side
REAL_T, intent(in) :: dx(SDIM)
REAL_T local_VEL(SDIM)
INTEGER_T velsolid_flag
REAL_T temp
INTEGER_T for_dt

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if (probtype.eq.55) then
 ! do nothing
else
 print *,"expecting probtype==55"
 stop
endif
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call GENERAL_PHASE_CHANGE_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)

 if (adv_dir.eq.veldir) then
  call rampvel(t,local_VEL(veldir))  ! default is adv_vel
 else
  local_VEL(veldir) = zero
 endif

 if ((veldir.eq.1).and.(dir.eq.1).and.(SDIM.eq.2)) then
  if ((yblob10.gt.zero).and.(axis_dir.eq.6)) then
   if((xghost(2).ge.yblob2).and.(xghost(2).le.yblob10)) then
    temp = xghost(2)-yblob2
    local_VEL(1)=local_VEL(1)*(1.5d0*temp/yblob10 - half*(temp/yblob10)**3)
   end if
  end if
 else if ((veldir.eq.2).and.(dir.eq.2).and.(side.eq.2).and.(SDIM.eq.2)) then
  if (axis_dir.eq.7) then ! compressible
   for_dt=0
   call acoustic_pulse_bc(t,local_VEL(veldir),xghost,for_dt)
  endif
 else if ((veldir.eq.3).and.(dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then
  if (axis_dir.eq.7) then ! compressible
   for_dt=0
   call acoustic_pulse_bc(t,local_VEL(veldir),xghost,for_dt)
  endif
 endif
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine GENERAL_PHASE_CHANGE_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(inout) :: PRES
REAL_T, intent(in) :: PRES_in
INTEGER_T, intent(in) :: dir,side
REAL_T, intent(in) :: dx(SDIM)
REAL_T base_pres
REAL_T gravity_dz

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 if (gravity_dir.eq.1) then
  gravity_dz=xghost(1)-probhix
 else if (gravity_dir.eq.2) then
  gravity_dz=xghost(2)-probhiy
 else if ((gravity_dir.eq.3).and.(SDIM.eq.3)) then
  gravity_dz=xghost(SDIM)-probhiz
 else
  print *,"gravity_dir invalid"
  stop
 endif
 call GENERAL_PHASE_CHANGE_PRES(xghost,t,LS,PRES,nmat)
 if (probtype.eq.55) then

  base_pres=zero

 else
  print *,"expecting probtype.eq.55"
  stop
 endif

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine GENERAL_PHASE_CHANGE_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T :: local_STATE(nmat*num_state_material)
REAL_T, intent(inout) :: STATE
REAL_T, intent(inout) :: STATE_merge
REAL_T, intent(in) :: STATE_in
INTEGER_T, intent(in) :: dir,side
REAL_T, intent(in) :: dx(SDIM)
INTEGER_T, intent(in) :: istate,im
INTEGER_T ibase,im_crit,im_loop
INTEGER_T local_bcflag

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
 call GENERAL_PHASE_CHANGE_STATE(xghost,t,LS,local_STATE,local_bcflag, &
         nmat,num_state_material)
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

 if (probtype.eq.55) then

  STATE=STATE_in
  STATE_merge=STATE

   ! xlo or xhi
  if ((dir.eq.1).and.(SDIM.eq.2)) then
   if (istate.eq.1) then ! density
    ! do nothing 
   else if (istate.eq.2) then ! temperature
    ! bcflag=1 (calling from denBC - boundary conditions
    ! for density, temperature and species variables)
    call outside_temperature(t,xghost(1),xghost(2),xghost(SDIM),STATE,im,1) 
    STATE_merge=STATE
   else
    print *,"istate invalid"
    stop
   endif

   ! ylo
  else if ((dir.eq.2).and.(side.eq.1).and.(SDIM.eq.2)) then

   ! prescribe_temperature_outflow=3 =>
   !  Dirichlet for inflow, outflow, wall; ylo states
   if ((prescribe_temperature_outflow.eq.3).and. &
       (axis_dir.eq.1)) then

    if (num_materials.lt.3) then
     print *,"num_materials invalid probtype=55"
     stop
    endif
   
     ! ylo
    if (istate.eq.1) then
     ! do nothing (density)
    else if (istate.eq.2) then
     STATE=fort_tempconst(3)  ! ice temperature for bottom of substrate.
     STATE_merge=STATE
    else
     print *,"istate invalid"
     stop
    endif

      ! freezing singularity or nucleate boiling problem: ylo
   else if ((prescribe_temperature_outflow.eq.3).and. &
            ((axis_dir.eq.5).or. &  ! freezing drop on substrate
             (axis_dir.eq.6).or. &  ! incompressible boiling
             (axis_dir.eq.7))) then ! compressible boiling

    if (istate.eq.1) then ! density
     ! do nothing 
    else if (istate.eq.2) then ! temperature
     ! bcflag=1 (calling from denBC)
     call outside_temperature(t,xghost(1),xghost(2),xghost(SDIM),STATE,im,1) 
     STATE_merge=STATE
    else
     print *,"istate invalid"
     stop
    endif

   endif

   ! yhi
  else if ((dir.eq.2).and.(side.eq.2).and.(SDIM.eq.2)) then

   if (istate.eq.1) then ! density
    ! do nothing
   else if (istate.eq.2) then ! temperature
    ! bcflag=1 (calling from denBC)
    call outside_temperature(t,xghost(1),xghost(2),xghost(SDIM),STATE,im,1) 
    STATE_merge=STATE
   else
    print *,"istate invalid"
    stop
   endif

   ! zlo
  else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then

   if ((prescribe_temperature_outflow.eq.3).and. &
       (axis_dir.eq.1)) then
     
    if (num_materials.lt.3) then
     print *,"num_materials invalid probtype=55"
     stop
    endif
  
    if (istate.eq.1) then ! density
     ! do nothing 
    else if (istate.eq.2) then ! temperature
     STATE=fort_tempconst(3)  ! ice temperature at zlo
     STATE_merge=STATE
    else
     print *,"istate invalid"
     stop
    endif

    ! freezing singularity or nucleate boiling problem: zlo
   else if ((prescribe_temperature_outflow.eq.3).and. &
            ((axis_dir.eq.5).or. &  ! freezing drop on substrate
             (axis_dir.eq.6).or. &
             (axis_dir.eq.7))) then ! compressible boiling

    if (istate.eq.1) then ! density
     ! do nothing 
    else if (istate.eq.2) then ! temperature
      ! bcflag=1 (calling from denBC)
     call outside_temperature(t,xghost(1),xghost(2),xghost(SDIM),STATE,im,1) 
     STATE_merge=STATE
    else
     print *,"istate invalid"
     stop
    endif
   endif
  endif
 else
  print *,"expecting probtype ==55"
  stop
 endif
else
 print *,"istate invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_STATE_BC

! MITSUHIRO: THIS ROUTINE MUST BE CUSTOMIZED.
! water=material 1
! gas=material 2
! ice=material 3
! substrate=material 4
! if (VFRAC(3)>1/2) then
!  set MITSUHIRO_CUSTOM_TEMPERATURE accordingly
!  note: xsten(0,1),xsten(0,2),xsten(0,SDIM) is the coordinate at which a
!  heat source might be prescribed.
!  x(1)=xsten(0,1)
!  x(2)=xsten(0,2)
!  x(SDIM)=xsten(0,SDIM)
subroutine GENERAL_PHASE_CHANGE_HEATSOURCE(im,VFRAC,time,x, &
     xsten,nhalf,temp, &
     heat_source,den,CV,dt,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
INTEGER_T, intent(in) :: im
REAL_T, intent(in) :: VFRAC(nmat)
REAL_T, intent(in) :: time
INTEGER_T, intent(in) :: nhalf
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
REAL_T, intent(in) :: temp(nmat)
REAL_T, intent(in) :: den(nmat)
REAL_T, intent(in) :: CV(nmat)
REAL_T, intent(in) :: dt
REAL_T, intent(out) :: heat_source

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
 
return
end subroutine GENERAL_PHASE_CHANGE_HEATSOURCE


 ! MEHDI VAHAB HEAT SOURCE
 ! called from: GODUNOV_3D.F90, 
 !  HEATSOURCE_FACE when center cell is a solid
 ! and the opposite cell is a fluid cell.
 ! heat_dir=1,2,3
 ! heat_side=1,2
subroutine GENERAL_PHASE_CHANGE_EB_heat_source(time,dt,xsten,nhalf, &
      heat_flux,heat_dir,heat_side)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nhalf
REAL_T, dimension(-nhalf:nhalf,SDIM), intent(in) :: xsten
REAL_T, intent(in) :: time
REAL_T, intent(in) :: dt
REAL_T, intent(out) :: heat_flux
INTEGER_T, intent(out) :: heat_dir
INTEGER_T, intent(out) :: heat_side

if (time.lt.zero) then
 print *,"time invalid"
 stop
endif
if (dt.le.zero) then
 print *,"dt invalid"
 stop
endif

if (probtype.eq.55) then

 heat_flux=zero
 heat_dir=0
 heat_side=0

 ! Sato and Niceno
 if ((axis_dir.eq.6).and. &
     (zblob3.lt.zero)) then
  heat_flux=abs(zblob3)
  heat_dir=SDIM
  heat_side=2
 endif

else
 print *,"expecting probtype==55"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_EB_heat_source

  ! only called at faces with an adjoining solid cell and
  ! an adjoining fluid cell.
subroutine GENERAL_PHASE_CHANGE_microcell_heat_coeff(heatcoeff,dx,veldir)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: dx(SDIM)
INTEGER_T, intent(in) :: veldir
REAL_T, intent(inout) :: heatcoeff

if (probtype.eq.55) then

 if ((veldir.ge.0).and.(veldir.lt.SDIM)) then
  ! boiling: Sato and Niceno  or Tryggvason and Lu
  if ((axis_dir.eq.6).or. &
      (axis_dir.eq.7)) then ! compressible boiling
   if (zblob3.eq.zero) then  ! Dirichlet at z=zlo
    ! do nothing
   else if (zblob3.gt.dx(veldir+1)) then ! TSAT dirichlet
    heatcoeff=heatcoeff*dx(veldir+1)/zblob3
   else if (zblob3.gt.zero) then
    ! do nothing
   else if (zblob3.lt.zero) then
    ! do nothing - heat source is specified at solid cells
    ! that adjoin fluid cells. (Sato and Niceno)
   else
    print *,"zblob3 invalid"
    stop
   endif
  endif
 else
  print *,"veldir invalid"
  stop
 endif

else
 print *,"expecting probtype==55"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_microcell_heat_coeff

subroutine GENERAL_PHASE_CHANGE_velfreestream(problen,local_buffer)
use probcommon_module
IMPLICIT NONE
REAL_T, intent(inout) :: local_buffer(2*SDIM)
REAL_T, intent(in)    :: problen(SDIM)
REAL_T :: buf
INTEGER_T :: ibuf
INTEGER_T :: dirbc,side

if (probtype.eq.55) then
 if (axis_dir.eq.0) then
  buf=problen(1)/32.0

  if (SDIM.eq.3) then
   dirbc=2
   side=1 
   buf=problen(1)/128.0
   ibuf=(side-1)*SDIM+dirbc
   if (local_buffer(ibuf).eq.zero) then
    local_buffer(ibuf)=buf
   endif
   side=2
   ibuf=(side-1)*SDIM+dirbc
   if (local_buffer(ibuf).eq.zero) then
    local_buffer(ibuf)=buf
   endif
  endif ! sdim==3
  dirbc=1
  side=1
  ibuf=(side-1)*SDIM+dirbc
  if (local_buffer(ibuf).eq.zero) then
   local_buffer(ibuf)=buf
  endif
  side=2
  ibuf=(side-1)*SDIM+dirbc
  if (local_buffer(ibuf).eq.zero) then
   local_buffer(ibuf)=buf
  endif
 endif ! axis_dir.eq.0
else
 print *,"expecting probtype==55"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_velfreestream


subroutine GENERAL_PHASE_CHANGE_nucleation(nucleate_in,xsten,nhalf,make_seed)
use probcommon_module_types
use probcommon_module
use global_utility_module
IMPLICIT NONE
INTEGER_T, intent(in) :: nhalf
REAL_T, dimension(-nhalf:nhalf,SDIM), intent(in) :: xsten
INTEGER_T, intent(inout) :: make_seed
type(nucleation_parm_type_input), intent(in) :: nucleate_in
REAL_T :: LL
REAL_T :: dist
REAL_T :: x_point(SDIM)
INTEGER_T :: dir

LL=nucleate_in%LL

if (probtype.eq.55) then

 if ((axis_dir.eq.6).or. &  ! incompressible
     (axis_dir.eq.7)) then ! compressible

  if ((nucleate_in%im_source.ne.1).or. &
      (nucleate_in%im_dest.ne.2)) then
   print *,"im_source or im_dest invalid"
   stop
  endif
  if (LL.gt.zero) then
   ! do nothing
  else
   print *,"expecting latent heat to be positive"
   stop
  endif
  do dir=1,SDIM
   x_point(dir)=xsten(0,dir)
  enddo
  if ((nucleate_in%do_the_nucleate.eq.1).and. &
      (n_sites.gt.0)) then
   call nucleation_sites( &
          x_point, &
          dist, &
          nucleate_in%nucleate_pos)
   if (dist.le.zero) then
    make_seed=1
   endif
  else if ((nucleate_in%do_the_nucleate.eq.0).or. &
           (n_sites.eq.0)) then
   ! do nothing
  else
   print *,"do_the_nucleate or n_sites invalid"
   stop
  endif

 endif ! axis_dir.eq.6 or 7

else
 print *,"expecting probtype==55"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_nucleation

subroutine GENERAL_PHASE_CHANGE_ICE_SUBSTRATE_DISTANCE( &
                xtarget,dist)
use probcommon_module
use global_utility_module
IMPLICIT NONE

REAL_T, intent(in) :: xtarget(SDIM)
REAL_T, intent(out) :: dist

 if (probtype.eq.55) then
  ! dist > 0 in the substrate
  call ice_substrate_distance( &
    xtarget(1),xtarget(2),xtarget(SDIM),dist)
 else
  print *,"expecting probtype.eq.55"
  stop
 endif

end subroutine GENERAL_PHASE_CHANGE_ICE_SUBSTRATE_DISTANCE

! This routine is called from fort_summass which is declared in
! NAVIERSTOKES_3D.F90
! MITSUHIRO: modify this routine to get the Nusselt number
! set ns.ncomp_sum_int_user1=4
! GRID_DATA_IN%cellten expected to have 0 ghost cells.
! GRID_DATA_IN%lsfab expected to have 2 ghost cells.
! GRID_DATA_IN%slopes expected to have 2 ghost cells.
! GRID_DATA_IN%den expected to have 1 ghost cells.
! GRID_DATA_IN%vel expected to have 1 ghost cells.
! it is recommended to use "vel" if viscous force 
! is needed in ghost cells.
subroutine GENERAL_PHASE_CHANGE_SUMINT(GRID_DATA_IN,increment_out1, &
                increment_out2,nsum1,nsum2,isweep)
use probcommon_module_types
use probcommon_module
use global_utility_module

INTEGER_T, intent(in) :: nsum1,nsum2,isweep
type(user_defined_sum_int_type), intent(in) :: GRID_DATA_IN
REAL_T, intent(inout) :: increment_out1(nsum1)
REAL_T, intent(inout) :: increment_out2(nsum2)
INTEGER_T :: i,j,k,dir,im
REAL_T :: xlocal(SDIM)
REAL_T :: cell_dim(SDIM)
INTEGER_T :: temperature_component
REAL_T temperature_plus
REAL_T temperature_minus
REAL_T zplus,zminus,heat_flux,area_face
INTEGER_T im_solid

call checkbound_array(GRID_DATA_IN%fablo,GRID_DATA_IN%fabhi, &
        GRID_DATA_IN%lsfab,2,-1,411) 
call checkbound_array(GRID_DATA_IN%fablo,GRID_DATA_IN%fabhi, &
        GRID_DATA_IN%slopes,2,-1,413) 
call checkbound_array(GRID_DATA_IN%fablo,GRID_DATA_IN%fabhi, &
        GRID_DATA_IN%den,1,-1,413) 
call checkbound_array(GRID_DATA_IN%fablo,GRID_DATA_IN%fabhi, &
        GRID_DATA_IN%vel,1,-1,413) 
call checkbound_array(GRID_DATA_IN%fablo,GRID_DATA_IN%fabhi, &
        GRID_DATA_IN%visco,1,-1,413) 

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid

if (nsum1.eq.4) then
 if (isweep.eq.0) then
  do dir=1,nsum1
   increment_out1(dir)=zero
  enddo
  do dir=1,SDIM
   xlocal(dir)=GRID_DATA_IN%xsten(0,dir)
   cell_dim(dir)=GRID_DATA_IN%xsten(1,dir)-GRID_DATA_IN%xsten(-1,dir)
   if (cell_dim(dir).gt.zero) then
    ! do nothing
   else
    print *,"cell_dim(dir) invalid"
    stop
   endif
  enddo
  if (num_materials.ge.3) then

   im_solid=im_solid_primary()

   if ((im_solid.ge.1).and.(im_solid.le.num_materials)) then

    ! find the heat flux next to a horizontal hot plate
    do im=1,num_materials
     temperature_component=(im-1)*num_state_material+2
     if ((im.eq.1).or.(im.eq.2)) then ! liquid or gas
      if ((GRID_DATA_IN%lsfab(D_DECL(i,j,k),im).gt.zero).and. &
          (GRID_DATA_IN%lsfab(D_DECL(i,j,k),im_solid).le.zero)) then
       if (GRID_DATA_IN%lsfab(D_DECL(i,j,k-1),im_solid).gt.zero) then
        if (SDIM.eq.3) then
         temperature_plus= &
          GRID_DATA_IN%den(D_DECL(i,j,k+1),temperature_component)
         temperature_minus= &
          GRID_DATA_IN%den(D_DECL(i,j,k-1),temperature_component)
          ! fort_solidheat_flag:0=diffuse in solid 1=dirichlet 2=neumann
         zplus=GRID_DATA_IN%xsten(2,SDIM)
         if ((fort_solidheat_flag.eq.0).or. &
             (fort_solidheat_flag.eq.2)) then
          zminus=GRID_DATA_IN%xsten(-2,SDIM)
         else if (fort_solidheat_flag.eq.1) then
          zminus=GRID_DATA_IN%xsten(-1,SDIM)
         else
          print *,"fort_solidheat_flag invalid"
          stop
         endif

         heat_flux=fort_heatviscconst(im)* &
            (temperature_plus-temperature_minus)/ &
            (zplus-zminus)
         area_face=cell_dim(1)*cell_dim(2)
         if (im.eq.1) then ! liquid
          increment_out1(1)=-area_face*heat_flux
          increment_out1(2)=area_face
          increment_out1(3)=-area_face*heat_flux
          increment_out1(4)=area_face
         else if (im.eq.2) then
          increment_out1(3)=-area_face*heat_flux
          increment_out1(4)=area_face
         else
          print *,"im invalid"
          stop
         endif
! in the "run.out" file (./amr2d ... inputs... >& run.out &
! ... TIME= ....  user_comp (1..ncomp_sum_int_user) 1  <total heat flux value>
! ... TIME= ....  user_comp (1..ncomp_sum_int_user) 2  <total area>
! ....
! in the inputs file (probtype==55), 
! ns.ncomp_sum_int_user=4
        endif
       endif
      endif
     else if ((im.ge.2).and.(im.le.num_materials)) then
      ! do nothing
     else
      print *,"im bust"
      stop
     endif
    enddo ! im=1..num_materials
   else
    print *,"im_solid invalid GENERAL_PHASE_CHANGE_SUMINT"
    stop
   endif
  else
   print *,"num_materials invalid GENERAL_PHASE_CHANGE_SUMINT"
   stop
  endif
 else if (isweep.eq.1) then
  ! do nothing
 else
  print *,"isweep invalid"
  stop
 endif
else if (nsum1.eq.0) then
 ! do nothing
else
 print *,"nsum1 invalid"
 print *,"nsum1 ",nsum1
 stop
endif

end subroutine GENERAL_PHASE_CHANGE_SUMINT

end module GENERAL_PHASE_CHANGE_module
