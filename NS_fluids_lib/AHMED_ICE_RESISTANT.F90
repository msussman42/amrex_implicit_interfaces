#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_FORT_INTEGER.H"
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

! probtype==425 (see run2d/inputs.AHMED_ICE_RESISTANT or
! run3d/inputs.AHMED_ICE_RESISTANT3d)
module AHMED_ICE_RESISTANT_module

implicit none                   

REAL_T :: DEF_VAPOR_GAMMA

contains

  ! do any initial preparation needed
subroutine INIT_AHMED_ICE_RESISTANT_MODULE()
IMPLICIT NONE

  DEF_VAPOR_GAMMA =  1.666666667D0

return
end subroutine INIT_AHMED_ICE_RESISTANT_MODULE

! Phi>0 in the solid
subroutine AHMED_substrateLS(x,Phi) 
use probcommon_module
use global_utility_module
implicit none
REAL_T, INTENT(in), dimension(SDIM) :: x !spatial coordinates
REAL_T, INTENT(out) :: Phi !LS dist, Phi>0 in the substrate
REAL_T :: yhalf,xshift
REAL_T :: local_time
REAL_T :: ptb_dist_low
REAL_T :: ptb_dist_high
INTEGER_T :: im_substrate
INTEGER_T :: expected_nmat

!drop falling vertically with given velocity onto oil lubricated surface
if (axis_dir.eq.0) then
 expected_nmat=4
!drop falling at an angle with given velocity onto a 
!  non-lubricated surface with horizontal air flow.
else if (axis_dir.eq.1) then
 expected_nmat=3
!drop falling at an angle with given velocity onto a 
! lubricated surface with horizontal air flow.
else if (axis_dir.eq.2) then
 expected_nmat=4
!ns.tension liq-air, liq-ice, liq-substrate, air-ice, air-substrate, 
! ice-substrate
!thermal conductivity: all materials
!specific heat: all materials
!density,viscosity:  "  "
!latent heat
!saturation temperature
!drop falling at an angle with given velocity onto a 
!  non-lubricated surface with horizontal air flow and 
!  ice-seed (sheet) on the surface.
!  im=1 drop, im=2 air, im=3 ice, im=4 substrate
else if (axis_dir.eq.3) then
 expected_nmat=4
else
 print *,"axis_dir invalid"
 stop
endif

if ((num_materials.eq.expected_nmat).and.(probtype.eq.425)) then
 local_time=zero
 im_substrate=num_materials
 if (SDIM.eq.2) then
  yhalf=0.2d0
  xshift=x(1)+0.2d0
 else if (SDIM.eq.3) then
  yhalf=x(2)
  xshift=x(1)
 else
  print *,"dimension bust"
  stop
 endif

  ! patterned_substrates is declared in GLOBALUTIL.F90.
  ! Phi<0 in the substrate, Phi>0 in the fluid
 ptb_dist_low=0.05d0
 ptb_dist_high=0.1d0
 if ((axis_dir.eq.0).or. &
     (axis_dir.eq.1).or. &
     (axis_dir.eq.2)) then
  ! do nothing
 else if (axis_dir.eq.3) then
  ptb_dist_low=ptb_dist_high
 else
  print *,"axis_dir invalid"
  stop
 endif

  ! patterned_substrates is declared in: GLOBALUTIL.F90
 call patterned_substrates(xshift,yhalf,x(SDIM),Phi,local_time,im_substrate, &
         ptb_dist_low,ptb_dist_high)
 Phi=-Phi
else
 print *,"num_materials or probtype invalid"
 stop
endif

end subroutine AHMED_substrateLS


subroutine AHMED_ICE_RESISTANT_check_vel_rigid(x,t,vel,dir)
use probcommon_module
IMPLICIT NONE

REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: vel
INTEGER_T, INTENT(in) :: dir
INTEGER_T :: expected_nmat

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

if (axis_dir.eq.0) then
 expected_nmat=4
else if (axis_dir.eq.1) then
 expected_nmat=3
else if (axis_dir.eq.2) then
 expected_nmat=4
else if (axis_dir.eq.3) then
 expected_nmat=4
else
 print *,"axis_dir invalid"
 stop
endif

if ((num_materials.eq.expected_nmat).and.(probtype.eq.425)) then
 if (vel.eq.0.0d0) then
  ! do nothing
 else
  print *,"AHMED_ICE_RESISTANT_check_vel_rigid: vel not expected"
  print *,"x,y,z ",x(1),x(2),x(SDIM)
  print *,"t ",t
  print *,"dir,vel ",dir,vel
  print *,"axis_dir ",axis_dir
  stop
 endif
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_check_vel_rigid


subroutine AHMED_ICE_RESISTANT_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(out) :: LS(nmat)
REAL_T :: ice_vertical
REAL_T :: substrate_height
INTEGER_T :: expected_nmat

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

if (axis_dir.eq.0) then
 expected_nmat=4
else if (axis_dir.eq.1) then
 expected_nmat=3
else if (axis_dir.eq.2) then
 expected_nmat=4
else if (axis_dir.eq.3) then
 expected_nmat=4
else
 print *,"axis_dir invalid"
 stop
endif

if ((num_materials.eq.expected_nmat).and.(probtype.eq.425)) then
  ! fluids tessellate the domain, substrate is embedded.

  ! water is material 1
 if (SDIM.eq.2) then
  LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)
 else if (SDIM.eq.3) then
  LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+ &
                     (x(SDIM)-zblob)**2)
 else
  print *,"dimension bust"
  stop
 endif

 if ((axis_dir.eq.0).or. &
     (axis_dir.eq.2)) then
  ! oil is material 3
  ! zblob2 is the altitude of the oil layer
  LS(3)=zblob2-x(SDIM)
  ! air
  LS(2)=-(max(LS(1),LS(3)))
 else if (axis_dir.eq.1) then
  LS(2)=-LS(1)

  ! water, air, ice, substrate
 else if (axis_dir.eq.3) then
  ! zblob2 is the altitude of the ice layer
  LS(3)=zblob2-x(SDIM)
  ! air
  LS(2)=-(max(LS(1),LS(3)))
 else
  print *,"axis_dir invalid"
  stop
 endif

 call AHMED_substrateLS(x,LS(num_materials))
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_LS

! initial velocity is zero
subroutine AHMED_ICE_RESISTANT_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: dx(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, INTENT(in) :: velsolid_flag
INTEGER_T :: expected_nmat
INTEGER_T :: assert_oil_velocity

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

 if (axis_dir.eq.0) then
  expected_nmat=4
 else if (axis_dir.eq.1) then
  expected_nmat=3
 else if (axis_dir.eq.2) then
  expected_nmat=4
 else if (axis_dir.eq.3) then
  expected_nmat=4
 else
  print *,"axis_dir invalid"
  stop
 endif

 if ((num_materials.eq.expected_nmat).and.(probtype.eq.425)) then

  do dir=1,SDIM
   VEL(dir)=zero
  enddo

  if (axis_dir.eq.1) then
   assert_oil_velocity=0 ! there is no oil or ice.
  else if ((axis_dir.eq.0).or. &
           (axis_dir.eq.2).or. &
           (axis_dir.eq.3)) then
   if (LS(3).ge.-dx(SDIM)) then
    assert_oil_velocity=1 ! oil/ice vel =0 at inflow wall.
   else if (LS(3).lt.-dx(SDIM)) then
    assert_oil_velocity=0
   else
    print *,"LS(3) invalid"
    stop
   endif
  else
   print *,"axis_dir invalid"
   stop
  endif

  if ((LS(nmat).ge.zero).or. &
      (velsolid_flag.eq.1).or. &
      (assert_oil_velocity.eq.1)) then
   ! the current value of "VEL" is used (which is VEL=0)
  else if ((LS(nmat).lt.zero).and. &
           (velsolid_flag.eq.0).and. &
           (assert_oil_velocity.eq.0)) then

   if (axis_dir.eq.0) then

    if (LS(1).ge.-dx(1)) then
     VEL(SDIM)=-abs(adv_vel)
    else if (LS(1).le.-dx(1)) then
     ! do nothing
    else
     print *,"LS(1) invalid"
     stop
    endif

   else if ((axis_dir.eq.1).or. &
            (axis_dir.eq.2).or. &
            (axis_dir.eq.3)) then
    VEL(1)=adv_vel

    if (LS(1).ge.-dx(1)) then ! near or inside the liquid drop.
     VEL(1)=xblob3
     if (SDIM.eq.2) then
      VEL(SDIM)=yblob3
     else if (SDIM.eq.3) then
      VEL(SDIM)=zblob3
     else
      print *,"SDIM invalid"
      stop
     endif
    else if (LS(1).le.-dx(1)) then
     ! do nothing
    else
     print *,"LS(1) invalid"
     stop
    endif

   else
    print *,"axis_dir invalid"
    stop
   endif
  else
   print *,"LS, velsolid, or assert_oil_velocity invalid"
   stop
  endif

 else
  print *,"num_materials or probtype invalid"
  stop
 endif

return 
end subroutine AHMED_ICE_RESISTANT_VEL


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine AHMED_ICE_RESISTANT_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: PRES

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
PRES=zero

return 
end subroutine AHMED_ICE_RESISTANT_PRES



subroutine AHMED_ICE_RESISTANT_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, INTENT(in) :: nmat
INTEGER_T, INTENT(in) :: nstate_mat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: STATE(nmat*nstate_mat)
INTEGER_T im,ibase,n
INTEGER_T :: expected_nmat

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

if (axis_dir.eq.0) then
 expected_nmat=4
else if (axis_dir.eq.1) then
 expected_nmat=3
else if (axis_dir.eq.2) then
 expected_nmat=4
else if (axis_dir.eq.3) then
 expected_nmat=4
else
 print *,"axis_dir invalid"
 stop
endif

if ((num_materials.eq.expected_nmat).and. &
    (num_state_material.ge.2).and. & ! density, temperature, vapor spec
    (probtype.eq.425)) then
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
   STATE(ibase+ENUM_SPECIESVAR+n)= &
       fort_speciesconst((n-1)*num_materials+im)
  enddo
 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine AHMED_ICE_RESISTANT_STATE

 ! dir=1..sdim  side=1..2
subroutine AHMED_ICE_RESISTANT_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(inout) :: LS(nmat)
REAL_T, INTENT(in) :: LS_in(nmat)
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) ::  dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call AHMED_ICE_RESISTANT_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine AHMED_ICE_RESISTANT_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(inout) :: VEL
REAL_T, INTENT(in) :: VEL_in
INTEGER_T, INTENT(in) :: veldir,dir,side
REAL_T, INTENT(in) :: dx(SDIM)
REAL_T local_VEL(SDIM)
INTEGER_T velsolid_flag

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

 call AHMED_ICE_RESISTANT_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine AHMED_ICE_RESISTANT_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(inout) :: PRES
REAL_T, INTENT(in) :: PRES_in
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call AHMED_ICE_RESISTANT_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine AHMED_ICE_RESISTANT_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T :: local_STATE(nmat*num_state_material)
REAL_T, INTENT(inout) :: STATE
REAL_T, INTENT(inout) :: STATE_merge
REAL_T, INTENT(in) :: STATE_in
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)
INTEGER_T, INTENT(in) :: istate,im
INTEGER_T ibase,im_crit
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
 call AHMED_ICE_RESISTANT_STATE(xghost,t,LS,local_STATE,local_bcflag, &
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
end subroutine AHMED_ICE_RESISTANT_STATE_BC

subroutine AHMED_ICE_RESISTANT_HEATSOURCE(im,VFRAC,time,x, &
     xsten,nhalf,temp, &
     heat_source,den,CV,dt,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
INTEGER_T, INTENT(in) :: im
REAL_T, INTENT(in) :: VFRAC(nmat)
REAL_T, INTENT(in) :: time
INTEGER_T, INTENT(in) :: nhalf
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
REAL_T, INTENT(in) :: temp(nmat)
REAL_T, INTENT(in) :: den(nmat)
REAL_T, INTENT(in) :: CV(nmat)
REAL_T, INTENT(in) :: dt
REAL_T, INTENT(out) :: heat_source
INTEGER_T :: expected_nmat

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if (axis_dir.eq.0) then
 expected_nmat=4
else if (axis_dir.eq.1) then
 expected_nmat=3
else if (axis_dir.eq.2) then
 expected_nmat=4
else if (axis_dir.eq.3) then
 expected_nmat=4
else
 print *,"axis_dir invalid"
 stop
endif

if ((num_materials.eq.expected_nmat).and. &
    (probtype.eq.425)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_HEATSOURCE

REAL_T function get_LS_seed(xpos,LS_normal)
use probcommon_module
IMPLICIT NONE

REAL_T, INTENT(in) :: xpos(SDIM)
REAL_T, INTENT(out) :: LS_normal(SDIM)
REAL_T :: mag
INTEGER_T :: dir

if (SDIM.eq.2) then
 mag=sqrt((xpos(1)-xblob4)**2+(xpos(2)-yblob4)**2)
else if (SDIM.eq.3) then
 mag=sqrt((xpos(1)-xblob4)**2+(xpos(2)-yblob4)**2+ &
               (xpos(SDIM)-zblob4)**2)
else
 print *,"dimension problem"
 stop
endif
get_LS_seed=radblob4-mag

do dir=1,SDIM
 LS_normal(dir)=zero
enddo

if (mag.gt.zero) then
 LS_normal(1)=(xblob4-xpos(1))/mag
 LS_normal(2)=(yblob4-xpos(2))/mag
 if (SDIM.eq.3) then
  LS_normal(SDIM)=(zblob4-xpos(SDIM))/mag
 endif
else if (mag.eq.zero) then
 ! do nothing
else
 print *,"mag invalid"
 stop
endif

end function get_LS_seed

subroutine AHMED_ICE_RESISTANT_ASSIMILATE( &
     assimilate_in,assimilate_out, &
     i,j,k,cell_flag)
use probcommon_module
use geometry_intersect_module
IMPLICIT NONE

type(assimilate_parm_type), INTENT(in) :: assimilate_in
type(assimilate_out_parm_type), INTENT(inout) :: assimilate_out
INTEGER_T, INTENT(in) :: i,j,k,cell_flag

INTEGER_T :: nstate,nstate_test
REAL_T :: xcrit(SDIM)
REAL_T :: LS_normal(SDIM)
REAL_T :: LS_test
REAL_T :: t_upper,t_lower
REAL_T :: LS_seed,LS_seed_buffer

INTEGER_T :: number_intervals
INTEGER_T :: dir
INTEGER_T :: im
REAL_T VEL_DROP(SDIM)
REAL_T ldata(D_DECL(3,3,3))
INTEGER_T :: i1,j1,k1,k1lo,k1hi
REAL_T :: xdata(SDIM)
REAL_T :: LS_normal_data(SDIM)
REAL_T :: vfrac_seed
REAL_T :: volcell
REAL_T :: facearea_seed
REAL_T :: centroid_seed(SDIM)
REAL_T :: cencell(SDIM)

nstate=assimilate_in%nstate

nstate_test=STATE_NCOMP
if (nstate.eq.nstate_test) then
 ! do nothing
else
 print *,"nstate invalid"
 print *,"nstate=",nstate
 print *,"nstate_test=",nstate_test
 stop
endif

if ((num_materials.ge.3).and. &
    (num_state_material.ge.2).and. & 
    (probtype.eq.425)) then
 do dir=1,SDIM
  xcrit(dir)=assimilate_in%xsten(0,dir)
 enddo
 if (assimilate_in%nhalf.ge.2) then
  ! do nothing
 else
  print *,"(assimilate_in%nhalf.ge.2) violated"
  stop
 endif
 t_upper=assimilate_in%cur_time  ! cur_time_slab
 t_lower=t_upper-assimilate_in%dt
 if (t_lower.lt.t_upper) then
  if (t_lower.ge.-VOFTOL) then
   ! do nothing
  else
   print *,"t_lower invalid"
   stop
  endif
  if (radblob5.gt.zero) then
   number_intervals=0
   do while (number_intervals*radblob5.le.t_lower)
    number_intervals=number_intervals+1
   enddo
   if (number_intervals*radblob5.gt.t_lower) then
    if (number_intervals*radblob5.le.t_upper) then

     LS_seed=get_LS_seed(xcrit,LS_normal)

     LS_seed_buffer=LS_seed+three*assimilate_in%dx(1)
     if (LS_seed_buffer.ge.zero) then

      VEL_DROP(1)=xblob5
      VEL_DROP(2)=yblob5
      if (SDIM.eq.3) then
       VEL_DROP(SDIM)=zblob5
      endif
      if (cell_flag.eq.0) then ! MAC GRID X
       assimilate_out%macx(D_DECL(i,j,k))=VEL_DROP(1)
      else if (cell_flag.eq.1) then ! MAC GRID Y
       assimilate_out%macy(D_DECL(i,j,k))=VEL_DROP(2)
      else if ((cell_flag.eq.2).and.(SDIM.eq.3)) then ! MAC GRID Z
       assimilate_out%macz(D_DECL(i,j,k))=VEL_DROP(SDIM)
      else if (cell_flag.eq.-1) then
       do dir=1,SDIM
        assimilate_out%state(D_DECL(i,j,k),STATECOMP_VEL+dir)=VEL_DROP(dir)
       enddo
       assimilate_out%LS_state(D_DECL(i,j,k),1)=LS_seed
       assimilate_out%LS_state(D_DECL(i,j,k),2)=-LS_seed
       do dir=1,SDIM
        assimilate_out%LS_state(D_DECL(i,j,k), &
                num_materials+dir)=LS_normal(dir)
        assimilate_out%LS_state(D_DECL(i,j,k), &
                num_materials+SDIM+dir)=-LS_normal(dir)
       enddo
       k1lo=0
       k1hi=0
       if (SDIM.eq.3) then
        k1lo=-1
        k1hi=1
       else if (SDIM.eq.2) then
        ! do nothing
       else
        print *,"dimension bust"
        stop
       endif

       do i1=-1,1
       do j1=-1,1
       do k1=k1lo,k1hi
        xdata(1)=assimilate_in%xsten(2*i1,1)
        xdata(2)=assimilate_in%xsten(2*j1,2)
        if (SDIM.eq.3) then
         xdata(SDIM)=assimilate_in%xsten(2*k1,SDIM)
        endif
        ldata(D_DECL(i1+2,j1+2,k1+2))=get_LS_seed(xdata,LS_normal_data)
       enddo
       enddo
       enddo
       call getvolume( &
         assimilate_in%bfact, &
         assimilate_in%dx, &
         assimilate_in%xsten, &
         assimilate_in%nhalf, &
         ldata, &
         vfrac_seed, &
         facearea_seed, &
         centroid_seed, &
         VOFTOL, &
         SDIM)
       call CISBOX( &
         assimilate_in%xsten, &
         assimilate_in%nhalf, &
         assimilate_in%xlo, &
         assimilate_in%dx, &
         i,j,k, &
         assimilate_in%bfact, &
         assimilate_in%level, &
         volcell,cencell,SDIM)

       assimilate_out%state(D_DECL(i,j,k),STATECOMP_MOF+1)=vfrac_seed
       assimilate_out%state(D_DECL(i,j,k),STATECOMP_MOF+ngeom_raw+1)= &
               one-vfrac_seed
       do dir=1,SDIM
        assimilate_out%state(D_DECL(i,j,k),STATECOMP_MOF+1+dir)= &
          centroid_seed(dir)-cencell(dir)
        assimilate_out%state(D_DECL(i,j,k),STATECOMP_MOF+ngeom_raw+1+dir)= &
          -(centroid_seed(dir)-cencell(dir))
       enddo

       do im=3,num_materials
        LS_test=assimilate_out%LS_state(D_DECL(i,j,k),im)
         ! oil and substrate should not be near the incoming droplets.
        if (LS_test.le.-three*assimilate_in%dx(1)) then
         ! do nothing
        else
         print *,"LS_test invalid LS_test=",LS_test
         print *,"LS_seed=",LS_seed
         stop
        endif
       enddo !do im=3,num_materials

      else 
       print *,"cell_flag invalid"
       stop
      endif
     else if (LS_seed_buffer.lt.zero) then
      ! do nothing
     else
      print *,"LS_seed_buffer is NaN"
      stop
     endif
    else if (number_intervals*radblob5.gt.t_upper) then
     ! do nothing
    else
     print *,"radblob5 or t_upper is NaN"
     stop
    endif 
   else
    print *,"statement failed:(number_intervals*radblob5.gt.t_lower)"
    stop
   endif 
  else if (radblob5.eq.zero) then
   !do nothing
  else
   print *,"radblob5 invalid"
   stop
  endif
 else if (t_lower.eq.t_upper) then
  ! do nothing
 else
  print *,"t_lower or t_upper invalid"
  stop
 endif
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_ASSIMILATE


end module AHMED_ICE_RESISTANT_module
