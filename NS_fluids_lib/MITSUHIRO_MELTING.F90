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

! probtype==414 (see run2d/inputs.MITSUHIRO_block_ice_melt)
module MITSUHIRO_MELTING_module

implicit none                   

contains

  ! do any initial preparation needed
subroutine INIT_MITSUHIRO_MELTING_MODULE()
IMPLICIT NONE

return
end subroutine INIT_MITSUHIRO_MELTING_MODULE

! Phi>0 in the solid
subroutine MITSUHIRO_substrateLS(x,Phi) 
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

end subroutine MITSUHIRO_substrateLS
!   ---------
!   | ice   |
!   ---------
!   |water  |
!-------------------
!   substrate

 ! fluids tessellate the domain, solids are immersed. 
subroutine MITSUHIRO_LS(x,t,LS)
use probcommon_module
IMPLICIT NONE

REAL_T x(SDIM)
REAL_T t
INTEGER_T im
REAL_T LS(num_materials)

if (abs(zblob2-yblob2).le.1.0D-14) then
 substrate_height=zblob2  ! substrate thickness
else
 print *,"zblob2 or yblob2 invalid"
 stop
endif

if ((num_materials.eq.4).and.(probtype.eq.414)) then

 if (SDIM.eq.2) then
  ice_vertical=yblob
 else if (SDIM.eq.3) then
  ice_vertical=zblob
 else
  print *,"dimension bust"
  stop
 endif
 if (abs(yblob2-(ice_vertical-half*radblob)).gt.VOFTOL) then
  print *,"bottom of original ice block must coincide w/substrate"
  stop
 endif

  !dist<0 inside the square
  !water below the ice
 if (SDIM.eq.2) then
  call squaredist(x,y,xblob-half*radblob,xblob+half*radblob, &
     -yblob2-radblob3,yblob2+radblob3,LS(1))
  LS(1)=-LS(1) ! water
  !ice (is above water)
  call squaredist(x,y,xblob-half*radblob,xblob+half*radblob, &
    yblob2+radblob3,yblob2+radblob,LS(3))
    LS(3)=-LS(3)

  !air; dist<0 inside the square
  call squaredist(x,y,xblob-half*radblob,xblob+half*radblob, &
    -yblob2-radblob,yblob2+radblob,LS(2))

 else if (SDIM.eq.3) then

  !dist<0 inside the square
  !water below the ice
  call cubedist(xblob-half*radblob,xblob+half*radblob, &
     yblob-half*radblob,yblob+half*radblob, &
     -zblob2-radblob3,zblob2+radblob3,x,y,z,LS(1))
  LS(1)=-LS(1)
  !ice
  call cubedist(xblob-half*radblob,xblob+half*radblob, &
     yblob-half*radblob,yblob+half*radblob, &
     zblob2+radblob3,zblob2+radblob,x,y,z,LS(3))

  !air; dist<0 inside the square
  !air everywhere not ice or water.
  ! important: fluids tessellate domain, solids (i.e. substrate)
  ! are embedded.
  call cubedist(xblob-half*radblob,xblob+half*radblob, &
     yblob-half*radblob,yblob+half*radblob, &
     -zblob2-radblob,zblob2+radblob,x,y,z,LS(2))  ! air

 else
  print *,"dimension bust"
  stop
 endif

 call MITSUHIRO_substrateLS(x,LS(4))
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine MITSUHIRO_LS

! initial velocity is zero
subroutine MITSUHIRO_LS_VEL(x,t,LS,VEL,velsolid_flag,dx)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(num_materials)
REAL_T, intent(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, intent(in) :: velsolid_flag

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
end subroutine MITSUHIRO_LS_VEL

! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine MITSUHIRO_PRES(x,t,LS,PRES)
use probcommon_module
IMPLICIT NONE

REAL_T x(SDIM)
REAL_T t
REAL_T LS(num_materials)
REAL_T PRES

PRES=zero

return 
end subroutine MITSUHIRO_PRES


subroutine MO_substrate_temperature(time,x,y,z,temperature)
use global_utility_module
IMPLICIT NONE

REAL_T, intent(in) :: time
REAL_T, intent(in) :: x,y,z
REAL_T, intent(out) :: temperature
REAL_T :: temp_slope
REAL_T substrate_height
INTEGER_T im_solid_temperature
REAL_T zprime

 if ((time.ge.zero).and.(time.le.1.0D+20)) then
  ! do nothing
 else if (time.ge.1.0D+20) then
  print *,"WARNING time.ge.1.0D+20 in MO_substrate_temperature"
 else if (time.lt.zero) then
  print *,"time invalid in MO_substrate_temperature"
  stop
 else
  print *,"time bust in MO_substrate_temperature"
  stop
 endif

 if (SDIM.eq.2) then
  if (abs(z-y).gt.VOFTOL) then
   print *,"expecting z=y in 2d routine MO_substrate_temperature"
   stop
  endif
 endif
 zprime=z
FIX ME 
 if (yblob2.gt.zero) then ! yblob2 is height of substrate 

  if (abs(yblob2-(yblob-half*radblob)).gt.VOFTOL) then
   print *,"initial bottom of ice block should coincide with substrate"
   stop
  endif
        if (radblob3.ge.radblob) then
         print *,"already melted portion exceeds original ice dimensions"
         stop
        endif

         ! called from initdata
        if (bcflag.eq.0) then
         if ((im.eq.1).or.(im.eq.2)) then ! water or air
          temperature=get_user_temperature(time,bcflag,im)
         else if (im.eq.3) then ! ice 
          temperature=get_user_temperature(time,bcflag,im)
         else if (im.eq.4) then ! substrate
          if (im_solid_temperature.ne.4) then
           print *,"im_solid_temperature invalid"
           stop
          endif
          if (num_materials.lt.4) then
           print *,"num_materials invalid"
           stop
          endif
           
          if (zprime.ge.yblob2) then ! above the substrate
           temperature=get_user_temperature(time,bcflag,3) ! ice
          else if ((zprime.ge.zero).and. &
                   (zprime.le.yblob2)) then
           temperature= &
            get_user_temperature(time,bcflag,im_solid_temperature)+ &
            (get_user_temperature(time,bcflag,3)- &
             get_user_temperature(time,bcflag,im_solid_temperature))* &
            zprime/yblob2 
          else if (zprime.le.zero) then
           temperature= &
            get_user_temperature(time,bcflag,im_solid_temperature) !substrate
          else
           print *,"zprime failure"
           stop
          endif
         else
          print *,"im invalid63"
          stop
         endif
        else if (bcflag.eq.1) then ! called from denBC
         if (im_solid_temperature.ne.4) then
          print *,"im_solid_temperature invalid"
          stop
         endif
          ! radblob3=thickness of underside of block already melted.
         if (zprime.ge.yblob2+radblob3) then
          temperature=get_user_temperature(time,bcflag,3) ! ice
         else if (zprime.ge.yblob2) then
          temperature=get_user_temperature(time,bcflag,3) ! water
         else if ((zprime.ge.zero).and. &
                  (zprime.le.yblob2)) then
          temperature= &
           get_user_temperature(time,bcflag,im_solid_temperature)+ &
           (get_user_temperature(time,bcflag,3)- &
            get_user_temperature(time,bcflag,im_solid_temperature))* &
            zprime/yblob2
         else if (zprime.le.zero) then
          ! substrate
          temperature=get_user_temperature(time,bcflag,im_solid_temperature) 
         else
          print *,"zprime failure"
          stop
         endif
        else
         print *,"bcflag invalid"
         stop
        endif
       endif ! yblob2>0

       ! in: outside_temperature
      else if (probtype.eq.55) then

       if (axis_dir.eq.5) then ! freezing: solid, ice, water, air
         ! substrate: 0<y<yblob2
        if (yblob2.gt.zero) then  
          ! called from initdata
         if (bcflag.eq.0) then
          if ((im.eq.1).or.(im.eq.2)) then ! water or air
           temperature=get_user_temperature(time,bcflag,im)
          else if (im.eq.3) then ! ice 
           temperature=get_user_temperature(time,bcflag,im)
          else if (im.eq.4) then ! substrate
           if (im_solid_temperature.ne.4) then
            print *,"im_solid_temperature invalid"
            stop
           endif
           if (num_materials.lt.4) then
            print *,"num_materials invalid"
            stop
           endif
           if (zprime.ge.yblob2) then
            temperature=get_user_temperature(time,bcflag,3) ! ice
           else if ((zprime.ge.zero).and. &
                    (zprime.le.yblob2)) then
            temperature= &
             get_user_temperature(time,bcflag,im_solid_temperature)+ &
             (get_user_temperature(time,bcflag,3)- &
              get_user_temperature(time,bcflag,im_solid_temperature))* &
             zprime/yblob2 
           else if (zprime.le.zero) then
            temperature= &
             get_user_temperature(time,bcflag,im_solid_temperature) !substrate
           else
            print *,"zprime failure"
            stop
           endif
          else
           print *,"im invalid64"
           stop
          endif
         else if (bcflag.eq.1) then ! called from denBC
          if (im_solid_temperature.ne.4) then
           print *,"im_solid_temperature invalid"
           stop
          endif
          if (zprime.ge.yblob2+radblob3) then
           temperature=get_user_temperature(time,bcflag,1) ! water
          else if (zprime.ge.yblob2) then
           temperature=get_user_temperature(time,bcflag,3) ! ice
          else if ((zprime.ge.zero).and. &
                   (zprime.le.yblob2)) then
           temperature= &
            get_user_temperature(time,bcflag,im_solid_temperature)+ &
            (get_user_temperature(time,bcflag,3)- &
             get_user_temperature(time,bcflag,im_solid_temperature))* &
             zprime/yblob2
          else if (zprime.le.zero) then
           ! substrate
           temperature=get_user_temperature(time,bcflag,im_solid_temperature) 
          else
           print *,"zprime failure"
           stop
          endif
         else
          print *,"bcflag invalid"
          stop
         endif
        endif ! yblob2>0

        ! boiling sites problem
        ! For Sato and Niceno problem:
        ! solidheat_flag=0 (diffuse in solid)
        ! zblob3<0.0 => heat source in solid cells that adjoin fluid cells.

       else if ((axis_dir.eq.6).or. &
                (axis_dir.eq.7)) then

         ! liquid, vapor, substrate,  or,
         ! liquid, vapor, gas, substrate
        if (im_solid_temperature.ne.num_materials) then 
         print *,"im_solid_temperature invalid"
         stop
        endif
          ! thermal layer thickness
        if (yblob3.le.zero) then
         print *,"yblob3 invalid"
         stop
        endif
        if (SDIM.eq.2) then
         substrate_height=yblob2
        else if (SDIM.eq.3) then 
         substrate_height=zblob2
        else
         print *,"dimension bust"
         stop
        endif

        if (radblob2.lt.zero) then
         print *,"radblob2 invalid"
        else if (radblob2.eq.zero) then
         ! do nothing
        else if (radblob2.gt.zero) then ! angle of incline
         if (levelrz.eq.1) then
          if ((SDIM.ne.2).or.(xblob2.ne.zero)) then
           print *,"parameters not supported"
           stop
          endif
          zprime=yblob2+(z-yblob2)*cos(radblob2)-x*sin(radblob2)
         else if (levelrz.eq.0) then
          if (SDIM.eq.2) then
           zprime=yblob2+(z-yblob2)*cos(radblob2)-(x-xblob2)*sin(radblob2)
          else if (SDIM.eq.3) then
           if (radblob3.ne.zero) then
            print *,"radblob3.ne.zero is not supported"
            stop
           endif
           zprime=zblob2+(z-zblob2)*cos(radblob2)-(x-xblob2)*sin(radblob2)
          else
           print *,"dimension bust"
           stop
          endif
         else
          print *,"levelrz not supported"
          stop
         endif
        else
         print *,"radblob2 invalid"
         stop
        endif

         ! (xblob2,yblob2,zblob2) is the "center" of the heated plate.
        if (zprime.ge.substrate_height+yblob3) then
         temperature=get_user_temperature(time,bcflag,1)
        else if (zprime.le.substrate_height) then
         temperature=get_user_temperature(time,bcflag,im_solid_temperature)
        else if ((substrate_height.le.zprime).and. &
                 (zprime.le.substrate_height+yblob3)) then
         temp_slope=(get_user_temperature(time,bcflag,im_solid_temperature)- &
                     get_user_temperature(time,bcflag,1))/yblob3
         temperature= &
          get_user_temperature(time,bcflag,im_solid_temperature)- &
          temp_slope*(zprime-substrate_height)
        else
         print *,"zprime or substrate_height invalid"
         stop
        endif
       endif

        ! in: outside_temperature
      else if (probtype.eq.710) then

       if (im_solid_temperature.ne.3) then
        print *,"im_solid_temperature invalid"
        stop
       endif
       temperature=get_user_temperature(time,bcflag,1)
          ! thermal layer thickness
       if (yblob3.le.zero) then
        print *,"yblob3 invalid"
        stop
       endif
       if (zprime.le.yblob3) then
        temperature=get_user_temperature(time,bcflag,im_solid_temperature)
       endif

      else
       print *,"expecting probtype=55, 59, or 710 in outside_temperature"
       stop
      endif
 
      return 
      end subroutine outside_temperature



subroutine MITSUHIRO_STATE(x,t,LS,STATE)
use probcommon_module
IMPLICIT NONE

REAL_T x(SDIM)
REAL_T t
REAL_T LS(num_materials)
REAL_T STATE(num_materials*num_state_material)
INTEGER_T im,ibase,n

if ((num_materials.eq.4).and. &
    (num_state_material.ge.3).and. & ! density, temperature, vapor
    (probtype.eq.414)) then
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
end subroutine MITSUHIRO_STATE

 ! dir=1..sdim  side=1..2
subroutine ZEYU_droplet_impact_LS_BC(xwall,xghost,t,LS, &
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
 call ZEYU_droplet_impact_LS(xghost,t,LS)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine ZEYU_droplet_impact_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(num_materials)
REAL_T, intent(out) :: VEL
REAL_T, intent(in) :: VEL_in
INTEGER_T, intent(in) :: veldir,dir,side
REAL_T, intent(in) :: dx(SDIM)
REAL_T local_VEL(SDIM)
INTEGER_T velsolid_flag

velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call ZEYU_droplet_impact_LS_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine ZEYU_droplet_impact_PRES_BC(xwall,xghost,t,LS, &
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

 call ZEYU_droplet_impact_PRES(xghost,t,LS,PRES)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine ZEYU_droplet_impact_STATE_BC(xwall,xghost,t,LS, &
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
 call ZEYU_droplet_impact_STATE(xghost,t,LS,local_STATE)
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
end subroutine ZEYU_droplet_impact_STATE_BC

subroutine ZEYU_droplet_impact_HEATSOURCE(im,VFRAC,time,x,temp, &
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

if ((num_materials.eq.2).and.(probtype.eq.412)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_HEATSOURCE

end module ZEYU_droplet_impact_module
