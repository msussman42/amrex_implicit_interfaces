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


  ! in 2D, "y" is ignored
  ! distance>0 in the drop
  ! if maxtall<radnew-vert then the distance
  ! to the ice part of the droplet is returned in dist
  ! and the distance to the remaining liquid part *above* the ice
  ! is returned in dist_truncate.
  ! this routine is called if probtype=55, axis_dir=0,1, or 5
subroutine drop_slope_dist(x,y,z,time,nmat, &
   maxtall,dist,dist_truncate)
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x,y,z,time
REAL_T, intent(out) :: dist,dist_truncate
REAL_T, intent(in) :: maxtall
INTEGER_T im,im_opp,im_3,iten_13,iten_23,imloop
INTEGER_T iten
REAL_T cos_angle,sin_angle
REAL_T term1,Vtarget,radnew,vert,test_angle
REAL_T xprime,yprime,zprime,rprime,rtop,rbot
REAL_T xcheck,ycheck,zcheck
INTEGER_T nten
REAL_T xvec(SDIM)
REAL_T marangoni_temp(nmat)
INTEGER_T im_solid_substrate
REAL_T, allocatable, dimension(:) :: user_tension

if (probtype.eq.55) then

 im_solid_substrate=im_solid_primary()

 xvec(1)=x
 xvec(2)=y
 if (SDIM.eq.3) then
  xvec(SDIM)=z
 endif

 if (nmat.ne.num_materials) then
  print *,"nmat invalid"
  stop
 endif
 nten=( (nmat-1)*(nmat-1)+nmat-1 )/2
 allocate(user_tension(nten))

 if (SDIM.eq.2) then
  if (abs(y-z).gt.VOFTOL) then
   print *,"y=z in 2d expected: drop_slope_dist"
   print *,"x,y,z= ",x,y,z
   stop
  endif
 endif
 if (maxtall.le.zero) then
  print *,"maxtall invalid"
  stop
 endif

  ! in: drop_slope_dist 
 if ((axis_dir.eq.0).or. &
     (axis_dir.eq.5)) then

  xcheck=xblob-xblob2
  ycheck=yblob-yblob2
  zcheck=zero
  if (SDIM.eq.3) then
   zcheck=zblob-zblob2
  endif

  if ((num_materials.ge.3).and. &
      (im_solid_substrate.ge.3).and. &
      (abs(xcheck).lt.1.0E-7).and. &
      (abs(ycheck).lt.1.0E-7).and. &
      (abs(zcheck).lt.1.0E-7)) then
   im=1
   im_opp=2
   im_3=im_solid_substrate
   call get_iten(im,im_opp,iten,num_materials)
   do imloop=1,nmat
    marangoni_temp(imloop)=293.0
   enddo
   call get_user_tension(xvec,time, &
     fort_tension,user_tension, &
     marangoni_temp, &
     nmat,nten,1)
     ! find angle between materials "im" and "im_3"
   call get_CL_iten(im,im_opp,im_3,iten_13,iten_23, &
    user_tension,nten,cos_angle,sin_angle)

    ! angles other than 0 or pi are supported:
    ! 0 < angle < pi
   if (abs(cos_angle).lt.one-1.0E-2) then 

    if (((SDIM.eq.3).and.(levelrz.eq.0)).or. &
        ((SDIM.eq.2).and.(levelrz.eq.1))) then
     term1=two/three-cos_angle+(cos_angle**3)/three
     if (term1.le.zero) then
      print *,"term1 invalid"
      stop
     endif
         
     Vtarget=half*(four/three)*Pi*(radblob**3)
     radnew=(Vtarget/(Pi*term1))**(one/three)
     vert=-radnew*cos_angle
    else if ((SDIM.eq.2).and.(levelrz.eq.0)) then
     test_angle=acos(abs(cos_angle))  ! 0<test_angle<=pi/2
     if (cos_angle.ge.zero) then
      term1=test_angle-half*sin(two*test_angle)
     else
      term1=Pi-test_angle+half*sin(two*test_angle)
     endif
     if (term1.le.zero) then
      print *,"term1 invalid"
      stop
     endif
     Vtarget=half*Pi*(radblob**2)
     radnew=sqrt(Vtarget/term1)
     vert=-radnew*cos_angle
    else
     print *,"dimension bust"
     stop
    endif
      ! rotate clockwise
      ! and shift "center" of inclined plane to origin
      ! need to modify if rotate about y-z plane.
    if (SDIM.eq.2) then
     xprime=(x-xblob2)*cos(radblob2)+(z-yblob2)*sin(radblob2)
     zprime=-(x-xblob2)*sin(radblob2)+(z-yblob2)*cos(radblob2)
     rprime=abs(xprime)
    else if (SDIM.eq.3) then
     xprime=(x-xblob2)*cos(radblob2)+(z-zblob2)*sin(radblob2)
     yprime=y-yblob2
     zprime=-(x-xblob2)*sin(radblob2)+(z-zblob2)*cos(radblob2)
     rprime=sqrt(xprime**2+yprime**2)
    else
     print *,"dimension bust"
     stop
    endif

     ! dist>0 in the liquid drop
    dist=radnew-sqrt(rprime**2+(zprime-vert)**2)
    dist_truncate=dist

     ! find distance to ice part of this droplet; also
     ! find distance to remaining liquid part above the ice.
    if (maxtall-vert.lt.radnew) then
     rtop=sqrt(radnew**2-(maxtall-vert)**2)
     rbot=sqrt(radnew**2-vert**2)

      ! outside drop, and above the ice.
     if ((dist.le.zero).and.(zprime.ge.maxtall)) then
      dist_truncate=dist

      ! inside the original drop.
     else if ((zprime.le.maxtall).and.(rprime.le.rtop)) then
      dist_truncate=zprime-maxtall

      ! outside the original drop, off to side of ice.
     else if ((zprime.le.maxtall).and.(rprime.ge.rtop)) then
      dist_truncate=-sqrt((rprime-rtop)**2+(zprime-maxtall)**2)
     else if (dist.ge.zero) then
      if (dist.lt.zprime-maxtall) then
       dist_truncate=dist
      else
       dist_truncate=zprime-maxtall
      endif 
     else
      print *,"dist invalid drop_slope_dist"
      stop
     endif

     if ((dist.lt.zero).and.(zprime.gt.vert+radnew)) then
      dist=maxtall-zprime
     else if ((dist.ge.zero).and.(zprime.ge.maxtall)) then
      dist=maxtall-zprime
     else if ((dist.ge.zero).and.(zprime.le.maxtall).and. &
              (zprime.ge.half*maxtall)) then
      if (dist.lt.maxtall-zprime) then
       ! do nothing
      else
       dist=maxtall-zprime
      endif
     else if ((dist.ge.zero).and.(zprime.le.half*maxtall).and. &
              (zprime.ge.zero)) then
      if (dist.lt.zprime) then
       ! do nothing
      else
       dist=zprime
      endif
     else if ((dist.ge.zero).and.(zprime.le.zero)) then
      dist=zprime
     else if ((dist.lt.zero).and.(zprime.lt.vert-radnew)) then
      dist=zprime
     else if ((dist.lt.zero).and.(zprime.ge.maxtall)) then
      dist=-sqrt((rprime-rtop)**2+(zprime-maxtall)**2)
     else if ((dist.lt.zero).and.(zprime.ge.zero)) then
      ! do nothing
     else if ((dist.lt.zero).and.(zprime.le.zero)) then
      dist=-sqrt((rprime-rbot)**2+zprime**2)
     else
      print *,"dist or zprime invalid"
      stop
     endif 
    endif !  maxtall-vert<radnew
     
   else
    print *,"contact angle too close to 0 or pi for drop on slope"
    print *,"probtype=",probtype
    print *,"radblob=",radblob
    print *,"radblob2=",radblob2
    print *,"radblob4=",radblob4
    print *,"radblob5=",radblob5
    print *,"radblob6=",radblob6
    print *,"radblob7=",radblob7
    stop
   endif

  else
   print *,"parameter conflict for probtype=55"
   stop
  endif

 else
  print *,"axis_dir incorrect   probtype,axis_dir=",probtype,axis_dir
  stop
 endif

 deallocate(user_tension)

else
 print *,"probtype invalid in drop_slope_dist"
 stop
endif

return
end subroutine drop_slope_dist


 ! this is velocity boundary condition at the top of the domain.  
subroutine acoustic_pulse_bc(time,vel_pulse,xsten,nhalf,for_dt)
IMPLICIT NONE

INTEGER_T, intent(in) :: for_dt
REAL_T, intent(in) :: time
REAL_T, intent(out) :: vel_pulse
INTEGER_T, intent(in) :: nhalf
REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
REAL_T x,y,z

x=xsten(0,1)
y=xsten(0,2)
z=xsten(0,SDIM)

if ((time.ge.zero).and.(time.le.1.0e+20)) then
 ! do nothing
else
 print *,"time invalid"
 stop
endif

if (abs(x)+abs(y)+abs(z).le.1.0e+20) then
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

 ! density_at_depth previously initialized by:
 ! init_density_at_depth() 
 ! called from: FORT_DENCOR, general_hydrostatic_pressure_density,
 ! boundary_hydrostatic, EOS_air_rho2, EOS_air_rho2_ADIABAT,
 ! SOUNDSQR_air_rho2, EOS_error_ind, presBDRYCOND, FORT_INITDATA 
subroutine GENERAL_PHASE_CHANGE_hydro_pressure_density( &
  xpos,rho,pres)
IMPLICIT NONE

REAL_T, intent(in) :: xpos(SDIM)
REAL_T, intent(out) :: rho
REAL_T, intent(out) :: pres
REAL_T denfree,zfree
REAL_T den_top,z_top
REAL_T z_at_depth

 if (probtype.eq.55) then
  if ((axis_dir.ge.0).and.(axis_dir.le.5)) then
   ! do nothing (compressible drop) ??
   FIX ME (boundary_hydrostatic)
  endif

  if (axis_dir.eq.7) then

   if (SDIM.eq.2) then
    z_at_depth=probloy
    z_top=probhiy
   else if (SDIM.eq.3) then
    z_at_depth=probloz
    z_top=probhiz
   else
    print *,"dimension bust GENERAL_PHASE_CHANGE_hydro_pressure_density"
    stop
   endif

   if (z_at_depth.ne.zero) then
    print *,"z_at_depth must be 0 for compressible boiling problem"
    stop
   endif

   den_top=fort_denconst(1)

   if (xpos(SDIM).gt.z_top) then
    rho=den_top ! atmos pressure at top of domain
   else
    rho= &
     ((density_at_depth-den_top)/ &
      (z_at_depth-z_top))*(xpos(SDIM)-z_top)+den_top
   endif
   call EOS_tait_ADIABATIC_rhohydro(rho,pres)


  else if (fort_material_type(1).eq.13) then

   denfree=fort_denconst(1)
   if (SDIM.eq.2) then
    zfree=probhiy
    z_at_depth=probloy
   else if (SDIM.eq.3) then
    zfree=probhiz
    z_at_depth=probloz
   else
    print *,"dimension bust"
    stop
   endif

   ! density_at_depth is found so that
   ! (p(density_at_depth)-p(rho_0))/(rho_0 (z_at_depth-zfree))=g
   !
   if (xpos(SDIM).gt.zfree) then
    rho=denfree
   else
    rho= &
     ((density_at_depth-denfree)/ &
      (z_at_depth-zfree))*(xpos(SDIM)-zfree)+denfree
   endif
   call EOS_tait_ADIABATIC_rhohydro(rho,pres)
  else
   print *,"axis_dir invalid GENERAL_PHASE_CHANGE_hydro_pressure_density"
   stop
  endif
 else
  print *,"probtype invalid GENERAL_PHASE_CHANGE_hydro_pressure_density"
  stop
 endif

return
end subroutine GENERAL_PHASE_CHANGE_hydro_pressure_density

subroutine GENERAL_PHASE_CHANGE_CFL_HELPER(time,dir,uu,dx)
INTEGER_T, intent(in) :: dir
REAL_T, intent(in) :: time
REAL_T, intent(inout) :: uu
REAL_T, intent(in) :: dx(SDIM)

INTEGER_T dir2
INTEGER_T side
REAL_T utest,uscale
REAL_T xsten_dummy(-1:1,SDIM)
INTEGER_T for_dt
INTEGER_T nhalf

nhalf=1

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
  do side=-nhalf,nhalf
   xsten_dummy(side,dir2)=dx(dir2)*half*side
  enddo
  enddo
  for_dt=1
  call acoustic_pulse_bc(time,utest,xsten_dummy,nhalf,for_dt)
  uu=max(abs(uu),abs(utest))
 endif
else
 print *,"unexpected probtype"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_CFL_HELPER


! Phi>0 in the solid
subroutine BOTTOM_substrateLS(x,Phi) 
use probcommon_module
implicit none
REAL_T, intent(in), dimension(SDIM) :: x !spatial coordinates
REAL_T, intent(out) :: Phi !LS dist, Phi>0 in the substrate

REAL_T substrate_height

if (abs(zblob2-yblob2).le.1.0D-14) then
 substrate_height=zblob2  ! substrate thickness
else
 print *,"zblob2 or yblob2 invalid"
 stop
endif

if (abs(x(SDIM)).le.1.0D+20) then
 Phi=substrate_height-x(SDIM)
else
 print *,"x(SDIM) invalid"
 stop
endif

end subroutine BOTTOM_substrateLS

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

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  if (probtype.eq.55) then

   ! drop on slope problem (2d or 3d)
   ! radblob4 is used as a "switch" in order to specify static
   ! solution at t=0.
   ! radblob5 is switch for drop hitting ellipse problem.
   if (axis_dir.eq.3) then

    if (SDIM.eq.2) then
     LS(1)=z-(zblob+radblob*cos(two*Pi*x/xblob))
    else if (SDIM.eq.3) then
     LS(1)=z-(zblob+radblob*cos(two*Pi*x/xblob)*cos(two*Pi*y/yblob))
    else
     print *,"dimension bust"
     stop
    endif
    LS(2)=-LS(1)

   else if (axis_dir.eq.5) then

     ! in: materialdistbatch (initial angle=static angle)
     ! maxtall==two*radblob => no ice in this call.
    call drop_slope_dist(x,y,z,initial_time,nmat, &
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
      call nucleation_sites(xsten,nhalf,dx,bfact,dist_liquid,pos_sites)
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
        LS(3)=x-(probhix-radblob10)
       else
        print *,"radblob10 invalid"
        stop
       endif
      else if (gravity_dir.eq.2) then
       if (radblob10.lt.probleny) then 
        LS(3)=y-(probhiy-radblob10)
       else
        print *,"radblob10 invalid"
        stop
       endif
      else if ((gravity_dir.eq.3).and.(SDIM.eq.3)) then
       if (radblob10.lt.problenz) then 
        LS(3)=z-(probhiz-radblob10)
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
     call squaredist(x,y,xblob-radblob,xblob+radblob,yblob, &
      yblob+radblob3,dist_liquid)
     dist_liquid=-dist_liquid
    else if (radblob3.eq.zero) then
     if (SDIM.eq.2) then
      dist_liquid=radblob-sqrt((x-xblob)**2+(y-yblob)**2)
     else
      dist_liquid=radblob-sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)
     endif
     if (radblob5.gt.zero) then
      if (SDIM.eq.2) then
       dist_liq2=radblob5-sqrt((x-xblob5)**2+(y-yblob5)**2)
      else
       dist_liq2=radblob5-sqrt((x-xblob5)**2+(y-yblob5)**2+(z-zblob5)**2)
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
     if ((abs(z-y).gt.VOFTOL).and.(SDIM.eq.2)) then
      print *,"z<>y prior to drop_slope_dist"
      stop
     endif
     ! in: materialdistbatch (initial angle=static angle)
     call drop_slope_dist(x,y,z,initial_time,nmat, &
      two*radblob,dist_ice,dist_liquid)
    else
     print *,"radblob4 invalid radblob4=",radblob4
     stop
    endif

    if (radblob6.gt.zero) then
     if (SDIM.eq.2) then
      dist_liq2=radblob6-sqrt((x-xblob6)**2+(y-yblob6)**2)
     else if (SDIM.eq.3) then
      dist_liq2=radblob6-sqrt((x-xblob6)**2+(y-yblob6)**2+(z-zblob6)**2)
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

     LS(3)=radblob7-sqrt( (x-xblob7)**2+(y-yblob7)**2 )  ! pos. in drop
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
    call ice_substrate_distance(x,y,z,distsolid)

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
    call drop_slope_dist(x,y,z,initial_time,nmat,radblob3, &
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
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, intent(in) :: velsolid_flag

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

return 
end subroutine GENERAL_PHASE_CHANGE_VEL


subroutine EOS_GENERAL_PHASE_CHANGE(rho,internal_energy,pressure, &
  imattype,im)
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: pressure

 call EOS_material_CORE(rho,internal_energy,pressure,imattype,im)

 return
end subroutine EOS_GENERAL_PHASE_CHANGE

subroutine SOUNDSQR_GENERAL_PHASE_CHANGE(rho,internal_energy,soundsqr, &
  imattype,im)
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: soundsqr
 REAL_T pressure

 call SOUNDSQR_material_CORE(rho,internal_energy,soundsqr, &
   imattype,im)

 return
end subroutine SOUNDSQR_GENERAL_PHASE_CHANGE


subroutine INTERNAL_GENERAL_PHASE_CHANGE(rho,temperature, &
  local_internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: temperature 
 REAL_T, intent(out) :: local_internal_energy

 call INTERNAL_material_CORE(rho,temperature,local_internal_energy, &
   imattype,im)

 return
end subroutine INTERNAL_GENERAL_PHASE_CHANGE

subroutine TEMPERATURE_GENERAL_PHASE_CHANGE(rho,temperature,internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(out) :: temperature 
 REAL_T, intent(in) :: internal_energy

 call TEMPERATURE_material_CORE(rho,temperature,internal_energy, &
     imattype,im)

 return
end subroutine TEMPERATURE_GENERAL_PHASE_CHANGE


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
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
if ((num_materials.eq.4).and. &
    (num_state_material.eq.2).and. & ! density, temperature
    (probtype.eq.2001)) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  STATE(ibase+1)=fort_denconst(im) ! density prescribed in the inputs file.
  if (t.eq.zero) then
   STATE(ibase+2)=fort_initial_temperature(im) !initial temperature in inputs
  else if (t.gt.zero) then
   STATE(ibase+2)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif
   ! always assume Dirichlet boundary condition at zlo for temperature.
  call outside_temperature(t,x(1),x(2),x(SDIM),STATE(ibase+2),im,bcflag)

   ! initial species in inputs?
  do n=1,num_species_var
   STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
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

 if ((veldir.eq.1).and.(dir.eq.1).and.(SDIM.eq.2)) then
  if ((yblob10.gt.zero).and.(axis_dir.eq.6)) then
   if((y.ge.yblob2).and.(y.le.yblob10)) then
    temp = y-yblob2
    local_VEL(1)=x_vel*(1.5d0*temp/yblob10 - half*(temp/yblob10)**3)
   end if
  end if
 else if ((veldir.eq.2).and.(dir.eq.2).and.(side.eq.2).and.(SDIM.eq.2)) then
  if (axis_dir.eq.7) then
   for_dt=0
   call acoustic_pulse_bc(time,local_VEL(veldir),xsten,nhalf,for_dt)
  endif
 else if ((veldir.eq.3).and.(dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then
  if (axis_dir.eq.7) then
   for_dt=0
   call acoustic_pulse_bc(time,local_VEL(veldir),xsten,nhalf,for_dt)
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

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call GENERAL_PHASE_CHANGE_PRES(xghost,t,LS,PRES,nmat)
 if (probtype.eq.55) then

  base_pres=zero
  if (fort_material_type(2).ne.0) then
   call general_hydrostatic_pressure(base_pres)
   PRES=base_pres
   if (fort_material_type(1).eq.13) then 
    call GENERAL_PHASE_CHANGE_hydro_pressure_density(xpos,rhohydro,PRES)
   else if (axis_dir.eq.6) then
    PRES=-fort_denconst(1)*abs(gravity)*gravity_dz
   endif
  endif

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

  if ((dir.eq.1).and.(SDIM.eq.2)) then
   if (istate.eq.1) then
    ! do nothing (density)
   else if (istate.eq.2) then
    ! bcflag=1 (calling from denBC - boundary conditions
    ! for density, temperature and species variables)
    call outside_temperature(time,x,y,z,STATE,im,1) 
    STATE_merge=STATE
   else
    print *,"istate invalid"
    stop
   endif
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
            ((axis_dir.eq.5).or. &
             (axis_dir.eq.6).or. &
             (axis_dir.eq.7))) then

    if (istate.eq.1) then
     ! do nothing (density)
    else if (istate.eq.2) then
     ! bcflag=1 (calling from denBC)
     call outside_temperature(time,x,y,z,STATE,im,1) 
     STATE_merge=STATE
    else
     print *,"istate invalid"
     stop
    endif

   endif
  else if ((dir.eq.2).and.(side.eq.2).and.(SDIM.eq.2)) then

   if (istate.eq.1) then
    ! do nothing
   else if (istate.eq.2) then
    ! bcflag=1 (calling from denBC)
    call outside_temperature(time,x,y,z,STATE,im,1) 
    STATE_merge=STATE
   else
    print *,"istate invalid"
    stop
   endif

  else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then

   if ((prescribe_temperature_outflow.eq.3).and. &
       (axis_dir.eq.1)) then
     
    if (num_materials.lt.3) then
     print *,"num_materials invalid probtype=55"
     stop
    endif
  
    if (istate.eq.1) then
     ! do nothing (density)
    else if (istate.eq.2) then
     STATE=fort_tempconst(3)  ! ice temperature at zlo
     STATE_merge=STATE
    else
     print *,"istate invalid"
     stop
    endif

    ! freezing singularity or nucleate boiling problem: zlo
   else if ((prescribe_temperature_outflow.eq.3).and. &
            ((axis_dir.eq.5).or. &
             (axis_dir.eq.6).or. &
             (axis_dir.eq.7))) then

    if (istate.eq.1) then
     ! do nothing (density)
    else if (istate.eq.2) then
      ! bcflag=1 (calling from denBC)
     call outside_temperature(time,x,y,z,STATE,im,1) 
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
REAL_T :: MITSUHIRO_CUSTOM_TEMPERATURE

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
 
 ! set a hot temperature in the liquid
if ((num_materials.eq.4).and.(probtype.eq.2001)) then

 heat_source=zero
 if (VFRAC(1).ge.half) then ! in the liquid
  MITSUHIRO_CUSTOM_TEMPERATURE=fort_tempconst(1)
  heat_source=(MITSUHIRO_CUSTOM_TEMPERATURE-temp(1))* &
               den(1)*CV(1)/dt
 else if (VFRAC(1).le.half) then
  ! do nothing
 else
  print *,"VFRAC(1) invalid"
  stop
 endif 

else
 print *,"num_materials or probtype invalid"
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
IMPLICIT NONE

REAL_T, intent(in) :: dx(SDIM)
INTEGER_T, intent(in) :: veldir
REAL_T, intent(inout) :: heatcoeff

if (probtype.eq.55) then

 if ((veldir.ge.0).and.(veldir.lt.SDIM)) then
  ! boiling: Sato and Niceno  or Tryggvason and Lu
  if ((axis_dir.eq.6).or. &
      (axis_dir.eq.7)) then
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
INTEGER_T, intent(in) :: nhalf
REAL_T, dimension(-nhalf:nhalf,SDIM), intent(in) :: xsten
INTEGER_T, intent(inout) :: make_seed
type(nucleation_parm_type_input), intent(in) :: nucleate_in
REAL_T :: LL
REAL_T :: dist

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
  if ((nucleate_in%do_the_nucleate.eq.1).and. &
      (n_sites.gt.0)) then
   call nucleation_sites(xsten,nhalf, &
          nucleate_in%dx, &
          nucleate_in%bfact, &
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


end module GENERAL_PHASE_CHANGE_module
