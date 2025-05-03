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

! probtype==423 (see run2d/inputs.CRYOGENIC_TANK_MK)
!
! IT IS ASSUMED THAT THE VERTICAL DIRECTION (i.e. the direction at 
! which gravity is applied) IS IN THE "Z" DIRECTION, BUT THE GEOMETRY
! VERTICAL DIRECTION IS DEFINED IN THE "Y" DIRECTION.  SO "Z" and "Y" must
! be swapped in 3D.
module CRYOGENIC_TANK_MK_module
use amrex_fort_module, only : amrex_real

implicit none

integer, PARAMETER :: TANK_MK_MATERIAL_TYPE=24

integer, PARAMETER :: TANK_MK_AUX_THICK_WALLS=0
integer, PARAMETER :: ZBOT_FLIGHT_ID=1
integer, PARAMETER :: TPCE_ID=0
integer, PARAMETER :: TANK_MK_GEOM_DESCRIPTOR=ZBOT_FLIGHT_ID

integer :: num_aux_expect

! x3D(1)=x(dir_x)=x(1)
! x3D(2)=x(dir_z)=x(SDIM)  vertical direction of geometry: "y"
! x3D(3)=0.0d0 if SDIM==2
! x3D(3)=x(dir_y)=x(2) if SDIM==3
integer, PARAMETER :: dir_x = 1
integer, PARAMETER :: dir_z = SDIM ! vertical direction of gravity force
! dir_y should not be used if SDIM==2 (dir_y=-1 as a placeholder)
! dir_y=2 if SDIM==3
integer :: dir_y  

!! MIDDLE OF THE TANK IS AT Z=0 
! Tank inner radius
real(amrex_real) :: TANK_MK_RADIUS
! Tank inner height (inner wall height)
real(amrex_real) :: TANK_MK_HEIGHT
! Location of liquid-gas interface in respect to
! middle of tank at z=0
real(amrex_real) :: TANK_MK_INTERFACE_LOCATION
! Tank speherical end radius
real(amrex_real) :: TANK_MK_END_RADIUS
! Tank spherical end curvature center (0,C_z)
real(amrex_real) :: TANK_MK_END_CENTER
! Heater flux
real(amrex_real) :: TANK_MK_HEATER_WATTS
! Heater location in dim=2 direction
real(amrex_real) :: TANK_MK_HEATER_WALL_MODEL
real(amrex_real) :: TANK_MK_HEATER_LOW
real(amrex_real) :: TANK_MK_HEATER_HIGH
real(amrex_real) :: TANK_MK_HEATER_R
real(amrex_real) :: TANK_MK_HEATER_R_LOW

real(amrex_real) :: TANK_MK_HEATER_THICK
real(amrex_real) :: TANK_MK_HEATER_TOTAL_ANGLE

real(amrex_real) :: TANK_MK_INSULATE_R
real(amrex_real) :: TANK_MK_INSULATE_R_HIGH
real(amrex_real) :: TANK_MK_INSULATE_THICK

real(amrex_real) :: TANK_MK_NOZZLE_RAD
real(amrex_real) :: TANK_MK_NOZZLE_HT
real(amrex_real) :: TANK_MK_NOZZLE_THICK_OUTLET
real(amrex_real) :: TANK_MK_NOZZLE_BASE

! Flat or spherical interface
real(amrex_real) :: TANK_MK_INTERFACE_RADIUS
real(amrex_real) :: TANK_MK_BUBBLE_X
real(amrex_real) :: TANK_MK_BUBBLE_Y
real(amrex_real) :: TANK_MK_BUBBLE_Z

! Initial mixture pressure at the hight point of the tank
real(amrex_real) :: TANK_MK_INITIAL_PRESSURE
! Universal gas constant [J/(mol K)]
real(amrex_real) :: TANK_MK_R_UNIV

real(amrex_real) :: TANK_MK_GAS_GAMMA
real(amrex_real) :: TANK_MK_GAS_CP
real(amrex_real) :: TANK_MK_GAS_CV

integer, parameter :: n_data=500
integer :: n_data_temp
integer :: n_data_xyz(3)
!data file (T deg K,z mm) -> internal data (z mm,T deg K)
real(amrex_real) :: parabolic_tinit(n_data,2) 
! (t seconds,ax (alpha g0))  g0=9.81 m/s^2
real(amrex_real) :: parabolic_xyz_accel(n_data,6) 

contains

! vertical direction is "y"
subroutine CRYOGENIC_TANK_MK_OPEN_CASFILE(part_id,unit_id,file_format)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: part_id
integer, INTENT(in) :: unit_id
integer, INTENT(out) :: file_format
integer :: stat

 print *,"CRYOGENIC_TANK_MK_OPEN_CASFILE should not be called"
 print *,"THIS ROUTINE WOULD BE USED IF LAGRANGIAN DATA"
 print *,"NEEDS TO BE DISTRIBUTED?"
 stop

 if (part_id.ne.1) then
  print *,"part_id invalid"
  stop
 endif

 file_format=1 ! vtk format
 if (axis_dir.eq.2) then
  open(unit=unit_id, file= 'tpce_geometry.vtk',status='old',iostat=stat)
  if (stat.ne.0) then
   print *,"tpce_geometry.vtk can not be opened"
   stop
  endif
 else
  print *,"expecting axis_dir.eq.2"
  stop
 endif

return
end subroutine CRYOGENIC_TANK_MK_OPEN_CASFILE


! vertical direction is "y"
subroutine CRYOGENIC_TANK_MK_GET_OUTSIDE_POINT( &
 exterior_BB, &
 xcell,time,x_outside,im_part,part_id)
use probcommon_module
use global_utility_module
real(amrex_real), INTENT(in) :: exterior_BB(3,2)
real(amrex_real), INTENT(in) :: xcell(3)
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(out) :: x_outside(3)
integer, INTENT(in) :: im_part
integer, INTENT(in) :: part_id
integer :: dir
real(amrex_real) :: BB_len

 do dir=1,3
  BB_len=exterior_BB(dir,2)-exterior_BB(dir,1)
  if (BB_len.gt.zero) then
   if (dir.eq.2) then
    x_outside(dir)=xcell(dir)
   else if ((dir.eq.1).or.(dir.eq.3)) then
    if (xcell(dir).gt.zero) then
     x_outside(dir)=exterior_BB(dir,2)+0.05d0*BB_len
    else if (xcell(dir).le.zero) then
     x_outside(dir)=exterior_BB(dir,1)-0.05d0*BB_len
    else
     print *,"xcell invalid"
     stop
    endif
   else
    print *,"dir invalid"
    stop
   endif
  else
   print *,"BB_len invalid"
   stop
  endif
 enddo ! dir=1..3
 
end subroutine CRYOGENIC_TANK_MK_GET_OUTSIDE_POINT

! vertical direction is "y"
subroutine CRYOGENIC_TANK_MK_OVERRIDE_FSI_SIGN_LS_VEL_TEMP( &
 exterior_BB, &
 interior_BB, &
 xcell,time,LS,VEL,TEMP,MASK,lev77,im_part,part_id)
use probcommon_module
use global_utility_module

real(amrex_real), INTENT(in) :: exterior_BB(3,2)
real(amrex_real), INTENT(in) :: interior_BB(3,2)
real(amrex_real), INTENT(in) :: xcell(3)
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(out) :: LS
real(amrex_real), INTENT(out) :: VEL(3)
real(amrex_real), INTENT(out) :: TEMP
integer, INTENT(out) :: MASK
integer, INTENT(in) :: lev77 !lev77=-1 for aux, >=0 otherwise.
integer, INTENT(in) :: im_part
integer, INTENT(in) :: part_id
integer :: dir

 if ((lev77.eq.-1).or. &
     (lev77.ge.1)) then
  ! do nothing
 else 
  print *,"lev77 invalid"
  stop
 endif

 MASK=FSI_NOTHING_VALID
 do dir=1,3
  VEL(dir)=0.0d0
 enddo
 TEMP=room_temperature  ! 293.0d0 (if double precision)

  ! the buffer size for the auxiliary mesh is 0.2*max_side_len on each
  ! side.
 if (lev77.eq.-1) then !lev77==-1 => aux grid used.

  if (TANK_MK_GEOM_DESCRIPTOR.eq.TPCE_ID) then

   if (TANK_MK_AUX_THICK_WALLS.eq.1) then

    if (part_id.eq.1) then ! heater_a (top heater)
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.2) then ! side heater
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.3) then ! source
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.4) then ! sink
     call check_outside_box(xcell,exterior_BB,LS,MASK)
     call check_inside_box(xcell,interior_BB,LS,MASK)
    else if (part_id.eq.5) then ! tank
     call check_outside_box(xcell,exterior_BB,LS,MASK)
     call check_inside_box(xcell,interior_BB,LS,MASK)
    else if (part_id.eq.6) then ! nozzle housing
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else
     print *,"part_id invalid"
     stop
    endif

   else if (TANK_MK_AUX_THICK_WALLS.eq.0) then

    if (part_id.eq.1) then ! heater_a (top heater)
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.2) then ! side heater
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.3) then ! source
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.4) then ! sink (2 masked off parts: TORUS shape)
     call check_outside_box(xcell,exterior_BB,LS,MASK)
     call check_inside_box(xcell,interior_BB,LS,MASK)
    else if (part_id.eq.5) then ! tank
     call check_outside_box(xcell,exterior_BB,LS,MASK)
     call check_inside_box(xcell,interior_BB,LS,MASK)
    else if (part_id.eq.6) then ! nozzle housing
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.7) then ! LAD housing
     call check_outside_box(xcell,exterior_BB,LS,MASK)
     call check_inside_box(xcell,interior_BB,LS,MASK)
    else
     print *,"part_id invalid"
     stop
    endif

   else 
    print *,"TANK_MK_AUX_THICK_WALLS invalid"
    stop
   endif

  else if (TANK_MK_GEOM_DESCRIPTOR.eq.ZBOT_FLIGHT_ID) then

   if (TANK_MK_AUX_THICK_WALLS.eq.1) then

    print *,"this option not supported"
    stop

   else if (TANK_MK_AUX_THICK_WALLS.eq.0) then

    if (part_id.eq.1) then ! side heater_a 
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.2) then ! side heater b
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.3) then ! source
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.4) then ! sink
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.5) then ! tank
     call check_outside_box(xcell,exterior_BB,LS,MASK)
     call check_inside_box(xcell,interior_BB,LS,MASK)
    else if (part_id.eq.6) then ! nozzle housing
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else if (part_id.eq.7) then ! LAD housing
     call check_outside_box(xcell,exterior_BB,LS,MASK)
    else
     print *,"part_id invalid"
     stop
    endif

   else 
    print *,"TANK_MK_AUX_THICK_WALLS invalid"
    stop
   endif

  else
   print *,"TANK_MK_GEOM_DESCRIPTOR invalid"
   stop
  endif

 else
  print *,"expecting all parts defined via the aux paradigm"
  stop
 endif
 
end subroutine CRYOGENIC_TANK_MK_OVERRIDE_FSI_SIGN_LS_VEL_TEMP


! vertical direction is "y"
subroutine CRYOGENIC_TANK_MK_BOUNDING_BOX_AUX(auxcomp, &
    minnode,maxnode,LS_FROM_SUBROUTINE,aux_ncells_max_side)
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: auxcomp
real(amrex_real), INTENT(inout) :: minnode(3)
real(amrex_real), INTENT(inout) :: maxnode(3)
integer, INTENT(out) :: LS_FROM_SUBROUTINE
integer, INTENT(out) :: aux_ncells_max_side

 if ((auxcomp.ge.1).and. &
     (auxcomp.le.num_aux_expect)) then
  if (axis_dir.eq.2) then
   if (num_materials.eq.3) then
    LS_FROM_SUBROUTINE=0

    if (TANK_MK_GEOM_DESCRIPTOR.eq.TPCE_ID) then

     if (TANK_MK_AUX_THICK_WALLS.eq.1) then
      if (auxcomp.eq.1) then ! heater a (top heater)
       aux_ncells_max_side=128
      else if (auxcomp.eq.2) then ! heater b
       aux_ncells_max_side=128
      else if (auxcomp.eq.3) then ! source
       aux_ncells_max_side=64
      else if (auxcomp.eq.4) then ! sink
       aux_ncells_max_side=128
      else if (auxcomp.eq.5) then ! tank
       aux_ncells_max_side=128
      else if (auxcomp.eq.6) then ! nozzle
       aux_ncells_max_side=64
      else
       print *,"auxcomp invalid"
       stop
      endif
     else if (TANK_MK_AUX_THICK_WALLS.eq.0) then

      if (auxcomp.eq.1) then ! heater a (top heater)
       aux_ncells_max_side=128
      else if (auxcomp.eq.2) then ! heater b
       aux_ncells_max_side=128
      else if (auxcomp.eq.3) then ! source
       aux_ncells_max_side=64
      else if (auxcomp.eq.4) then ! sink
       aux_ncells_max_side=128
      else if (auxcomp.eq.5) then ! tank
       aux_ncells_max_side=256
      else if (auxcomp.eq.6) then ! nozzle
       aux_ncells_max_side=128
      else if (auxcomp.eq.7) then ! LAD housing
       aux_ncells_max_side=128
      else
       print *,"auxcomp invalid"
       stop
      endif

     else 
      print *,"TANK_MK_AUX_THICK_WALLS invalid"
      stop
     endif

    else if (TANK_MK_GEOM_DESCRIPTOR.eq.ZBOT_FLIGHT_ID) then

     if (TANK_MK_AUX_THICK_WALLS.eq.1) then

      print *,"this option not supported"
      stop

     else if (TANK_MK_AUX_THICK_WALLS.eq.0) then

      if (auxcomp.eq.1) then ! side heater a 
       aux_ncells_max_side=64
      else if (auxcomp.eq.2) then ! side heater b
       aux_ncells_max_side=64
      else if (auxcomp.eq.3) then ! source
       aux_ncells_max_side=256
      else if (auxcomp.eq.4) then ! sink
       aux_ncells_max_side=256
      else if (auxcomp.eq.5) then ! tank
       aux_ncells_max_side=256  !256 for production, 64 for testing.
      else if (auxcomp.eq.6) then ! nozzle
       aux_ncells_max_side=256
      else if (auxcomp.eq.7) then ! LAD housing
       aux_ncells_max_side=256
      else
       print *,"auxcomp invalid"
       stop
      endif

     else 
      print *,"TANK_MK_AUX_THICK_WALLS invalid"
      stop
     endif

    else
     print *,"TANK_MK_GEOM_DESCRIPTOR invalid"
     stop
    endif

   else
    print *,"num_materials invalid in CRYOGENIC_TANK_MK_BOUNDING_BOX_AUX"
    stop
   endif
  else
   print *,"axis_dir invalid in CRYOGENIC_TANK_MK_BOUNDING_BOX_AUX"
   stop
  endif

 else
  print *,"auxcomp invalid CRYOGENIC_TANK_MK_BOUNDING_BOX_AUX"
  stop
 endif


end subroutine CRYOGENIC_TANK_MK_BOUNDING_BOX_AUX


! vertical direction is "y"
subroutine CRYOGENIC_TANK_MK_OPEN_AUXFILE(part_id,unit_id,file_format)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: part_id
integer, INTENT(in) :: unit_id
integer, INTENT(out) :: file_format
integer :: stat

 file_format=1 ! vtk format

 if (axis_dir.eq.2) then

  if (fort_num_local_aux_grids.eq.num_aux_expect) then

   if (TANK_MK_GEOM_DESCRIPTOR.eq.TPCE_ID) then

    if (TANK_MK_AUX_THICK_WALLS.eq.1) then

     if (part_id.eq.1) then ! top heater
      open(unit=unit_id,file= 'heatera_coarse.vtk',status='old',iostat=stat)
     else if (part_id.eq.2) then ! side heater
      open(unit=unit_id,file= 'heaterb_coarse.vtk',status='old',iostat=stat)
     else if (part_id.eq.3) then
      open(unit=unit_id,file= 'nozzlesource_coarse.vtk', &
        status='old',iostat=stat)
     else if (part_id.eq.4) then
      open(unit=unit_id,file= 'sink_coarse.vtk',status='old',iostat=stat)
     else if (part_id.eq.5) then
      open(unit=unit_id, file= 'tank_coarse.vtk',status='old', &
         iostat=stat)
     else if (part_id.eq.6) then
      open(unit=unit_id, file= 'nozzle_coarse.vtk',status='old', &
         iostat=stat)
     else
      print *,"part_id invalid"
      stop
     endif

    else if (TANK_MK_AUX_THICK_WALLS.eq.0) then

     if (part_id.eq.1) then ! top heater
      open(unit=unit_id,file= 'tpce_heatera.vtk',status='old',iostat=stat)
     else if (part_id.eq.2) then ! side heater
      open(unit=unit_id,file= 'tpce_heaterb.vtk',status='old',iostat=stat)
     else if (part_id.eq.3) then
      open(unit=unit_id,file= 'nozzle_source_15deg.vtk', &
              status='old',iostat=stat)
     else if (part_id.eq.4) then
      open(unit=unit_id,file= 'tpce_sink.vtk',status='old',iostat=stat)
     else if (part_id.eq.5) then
      open(unit=unit_id, file= 'tpce_shell.vtk',status='old', &
         iostat=stat)
     else if (part_id.eq.6) then
      if (1.eq.0) then
       open(unit=unit_id, file= 'nozzle_15deg.vtk',status='old', &
         iostat=stat)
      else
       open(unit=unit_id, file= 'nozzle_15deg_01thick.vtk',status='old', &
         iostat=stat)
      endif
     else if (part_id.eq.7) then
      open(unit=unit_id, file= 'tpce_ladhousing.vtk',status='old', &
         iostat=stat)
     else
      print *,"part_id invalid"
      stop
     endif
    else 
     print *,"TANK_MK_AUX_THICK_WALLS invalid"
     stop
    endif

   else if (TANK_MK_GEOM_DESCRIPTOR.eq.ZBOT_FLIGHT_ID) then

    if (TANK_MK_AUX_THICK_WALLS.eq.1) then

     print *,"this option not supported"
     stop

    else if (TANK_MK_AUX_THICK_WALLS.eq.0) then

     if (part_id.eq.1) then ! side heatera
      open(unit=unit_id,file= 'zbot_flight_heatera.vtk',status='old', &
              iostat=stat)
     else if (part_id.eq.2) then ! side heaterb
      open(unit=unit_id,file= 'zbot_flight_heaterb.vtk',status='old', &
              iostat=stat)
     else if (part_id.eq.3) then
      open(unit=unit_id,file= 'zbot_flight_inflow_thick.vtk',status='old', &
              iostat=stat)
     else if (part_id.eq.4) then
      open(unit=unit_id,file= 'zbot_flight_outflow_thick.vtk', &
              status='old',iostat=stat)
     else if (part_id.eq.5) then
      open(unit=unit_id,file= 'zbot_flight_tank_thicknozzles.vtk', &
              status='old',iostat=stat)
     else if (part_id.eq.6) then
      open(unit=unit_id,file= 'zbot_flight_inletnozzle_thick.vtk', &
              status='old',iostat=stat)
     else if (part_id.eq.7) then
      open(unit=unit_id,file= 'zbot_flight_outletnozzle_thick.vtk', &
              status='old',iostat=stat)
     else
      print *,"part_id invalid"
      stop
     endif
    else 
     print *,"TANK_MK_AUX_THICK_WALLS invalid"
     stop
    endif

   else 
    print *,"TANK_MK_GEOM_DESCRIPTOR invalid"
    stop
   endif

  else
   print *,"expecting fort_num_local_aux_grids=",num_aux_expect
   stop
  endif 

 else
  print *,"expecting axis_dir.eq.2"
  stop
 endif

return
end subroutine CRYOGENIC_TANK_MK_OPEN_AUXFILE

! z input is in MKS
! z data is mm
subroutine interp_parabolic_tinit(z,temperature)
IMPLICIT NONE

real(amrex_real), intent(in) :: z
real(amrex_real), intent(out) :: temperature
integer :: idata,index_1,index_2

index_1=1
index_2=2

if ((n_data_temp.ge.2).and.(n_data_temp.lt.n_data)) then
 !do nothing
else
 print *,"n_data_temp invalid: ",n_data_temp
 stop
endif

idata=1
do while ((parabolic_tinit(idata,index_1)*1.0D-3.lt.z).and. &
          (idata.lt.n_data_temp)) 
 idata=idata+1
 if ((parabolic_tinit(idata,index_1).ge. &
      parabolic_tinit(idata-1,index_1)).and. &
     (parabolic_tinit(idata-1,index_1).ge.zero)) then
  !do nothing
 else
  print *,"parabolic_tinit invalid"
  stop
 endif
enddo
temperature=parabolic_tinit(idata-1,index_2)+ &
     (parabolic_tinit(idata,index_2)-parabolic_tinit(idata-1,index_2))* &
     (z-parabolic_tinit(idata-1,index_1))/ &
     (parabolic_tinit(idata,index_1)-parabolic_tinit(idata-1,index_1))

if (temperature.ge.zero) then
 !do nothing
else
 print *,"temperature invalid: ",temperature
 stop
endif

return
end subroutine interp_parabolic_tinit

subroutine interp_parabolic_xyz(time,axyz)
IMPLICIT NONE

real(amrex_real), intent(in) :: time
real(amrex_real), intent(out) :: axyz(3)
integer :: idata,index_1,index_2,dir

do dir=1,3
 index_1=2*dir-1
 index_2=index_1+1

 if ((n_data_xyz(dir).ge.2).and.(n_data_xyz(dir).lt.n_data)) then
  !do nothing
 else
  print *,"n_data_xyz invalid: ",dir,n_data_xyz(dir)
  stop
 endif
 idata=1
 do while ((parabolic_xyz_accel(idata,index_1).lt.time).and. &
           (idata.lt.n_data_xyz(dir))) 
  idata=idata+1
  if ((parabolic_xyz_accel(idata,index_1).ge. &
       parabolic_xyz_accel(idata-1,index_1)).and. &
      (parabolic_xyz_accel(idata-1,index_1).ge.zero)) then
   !do nothing
  else
   print *,"parabolic_xyz_accel invalid"
   stop
  endif
 enddo
 axyz(dir)=parabolic_xyz_accel(idata-1,index_2)+ &
  (parabolic_xyz_accel(idata,index_2)-parabolic_xyz_accel(idata-1,index_2))* &
  (time-parabolic_xyz_accel(idata-1,index_1))/ &
  (parabolic_xyz_accel(idata,index_1)-parabolic_xyz_accel(idata-1,index_1))

enddo !dir=1,3

return
end subroutine interp_parabolic_xyz


 ! do any initial preparation needed
 subroutine INIT_CRYOGENIC_TANK_MK_MODULE()
  use probcommon_module
  implicit none

  integer :: dir
  integer :: idata,index_1,index_2

  num_aux_expect=0
 
  if (SDIM.eq.2) then
   dir_y=-1  ! dir_y is not used if SDIM==2
  else if (SDIM.eq.3) then
    !x3D(1)=x(1)   
    !x3D(2)=x(SDIM)   SDIM=gravity vert. dir.  2=vtk vertical dir.
    !x3D(3)=x(dir_y)=x(2)
   dir_y=2
  else
   print *,"dimension bust"
   stop
  endif

  TANK_MK_RADIUS             = xblob
  TANK_MK_HEIGHT             = yblob
  TANK_MK_INTERFACE_LOCATION = zblob

  TANK_MK_END_RADIUS       = xblob4
  TANK_MK_END_CENTER       = yblob4

  TANK_MK_INTERFACE_RADIUS = radblob2
  TANK_MK_BUBBLE_X         = xblob2
  TANK_MK_BUBBLE_Y         = yblob2
  TANK_MK_BUBBLE_Z         = zblob2

  TANK_MK_HEATER_WATTS      = xblob3

   ! see Barsi and Kassemi 2013, Journal of Thermal Science and Engineering
   ! Applications.
   ! ZBOT
  if (axis_dir.eq.0) then ! volume=pi(.09^2-0.08^2)*.02=pi(.01)(.17)(0.02)

   print *,"axis_dir=0, tank geometry prescribed internally (not from file)"

   number_of_source_regions=1 ! side heater

   TANK_MK_HEATER_WALL_MODEL = 0.1683d0
   TANK_MK_HEATER_LOW       = -0.1683d0
   TANK_MK_HEATER_HIGH      = TANK_MK_HEATER_LOW+0.0254d0
   TANK_MK_HEATER_R_LOW     = 0.1016d0
   TANK_MK_HEATER_R         = TANK_MK_HEATER_R_LOW+0.027d0

   TANK_MK_INSULATE_R = xblob+0.027d0 !xblob is the tank cavity radius
   TANK_MK_INSULATE_THICK = 0.0508
   TANK_MK_INSULATE_R_HIGH = yblob/2.0d0 !yblob is height of cylindrical part

   ! TPCE
  else if ((axis_dir.eq.1).or. &
           (axis_dir.eq.2)) then ! heater on top

   number_of_source_regions=3 ! heater A, inflow, outflow

   if (axis_dir.eq.1) then
    num_aux_expect=0
   else if (axis_dir.eq.2) then

    if (TANK_MK_AUX_THICK_WALLS.eq.0) then
     num_aux_expect=7
    else if (TANK_MK_AUX_THICK_WALLS.eq.1) then
     num_aux_expect=6
    else
     print *,"TANK_MK_AUX_THICK_WALLS invalid"
     stop
    endif
    print *,"TANK_MK_AUX_THICK_WALLS=",TANK_MK_AUX_THICK_WALLS
    print *,"num_aux_expect ",num_aux_expect
   else
    print *,"axis_dir invalid: ",axis_dir
    stop
   endif
  
   TANK_MK_HEATER_WALL_MODEL = 0.0
   if (radblob4.gt.0.0d0) then
!   TANK_MK_HEATER_THICK = 0.01d0
    TANK_MK_HEATER_THICK = radblob4
   else
    print *,"init radblob4 to be TANK_MK_HEATER_THICK "
    print *,"currently radblob4=",radblob4
    stop
   endif

   TANK_MK_HEATER_TOTAL_ANGLE = 30.0d0*Pi/180.0d0

   if (radblob3.gt.0.0d0) then
!   TANK_MK_NOZZLE_RAD=0.005D0  !dx coarsest =0.0015625, "1cm diameter."
    TANK_MK_NOZZLE_RAD=radblob3
   else
    print *,"init radblob3 to be TANK_MK_NOZZLE_RAD "
    print *,"currently radblob3=",radblob3
    stop
   endif

   TANK_MK_NOZZLE_HT=0.064D0
   if (radblob5.gt.0.0d0) then
!   TANK_MK_NOZZLE_THICK_OUTLET=0.005d0
    TANK_MK_NOZZLE_THICK_OUTLET=radblob5
   else
    print *,"init radblob5 to be TANK_MK_NOZZLE_THICK_OUTLET "
    print *,"currently radblob5=",radblob5
    stop
   endif

   TANK_MK_NOZZLE_BASE=-half*TANK_MK_HEIGHT-TANK_MK_END_RADIUS

   !Experimental characterization of non-isothermal sloshing in microgravity
   !Monteiro et al 2024
  else if (axis_dir.eq.3) then
   number_of_source_regions=0
   num_aux_expect=0

   print *,"reading parabolic_tinit"
   open(unit=2, file='parabolic_tinit')
   read(2,*) n_data_temp
   print *,"n_data_temp=",n_data_temp
   index_1=1
   index_2=2

    !temperature K
    !z mm
    !data file (T deg K,z mm) -> internal data (z mm,T deg K)
   do idata=1,n_data_temp
    read(2,*) parabolic_tinit(idata,index_2),parabolic_tinit(idata,index_1)
   enddo
   close(2)

   print *,"reading parabolic_xyz_accel"
   open(unit=2, file='parabolic_xyz_accel')
   read(2,*) n_data_xyz(1),n_data_xyz(2),n_data_xyz(3)

   do dir=1,3
    index_1=2*dir-1
    index_2=index_1+1

    print *,"dir, n_data_xyz=",dir,n_data_xyz(dir)
    !time second
    !accel = alpha * g0
    do idata=1,n_data_xyz(dir)
     read(2,*) parabolic_xyz_accel(idata,index_1), &
               parabolic_xyz_accel(idata,index_2)
    enddo

   enddo !dir=1,3

   close(2)

  else
   print *,"axis_dir invalid: ",axis_dir
   stop
  endif

  ! ASSUMING IDEAL GAS => The gas heat cpacities should satisfy this
  ! R_spc = C_{p,spc}-C_{v,spc}
  ! to have ideal mixture gas as well =>
  ! Only C_p or C_v can be picked from table and the other one
  ! calculated from equation above.
  ! Here we pick C_{v,spc} from input file.
  ! C_{p,spc} = C_{v,spc} + R_spc

  TANK_MK_R_UNIV = 8.31446261815324D0 !J/(mol Kelvin)
   !molar_mass units: kg/mol

  TANK_MK_GAS_CV = fort_stiffCV(2) ![J/(kg K)]
!  TANK_MK_GAS_CP = fort_stiffCP(2) ![J∕(kg·K)]
  TANK_MK_GAS_CP=TANK_MK_GAS_CV+TANK_MK_R_UNIV/fort_molar_mass(2) ![J∕(kg·K)]
  TANK_MK_GAS_GAMMA = TANK_MK_GAS_CP / TANK_MK_GAS_CV

  ! Initial pressure based on the given density and temeprature
  ! at the highest point of the tank
  if(fort_material_type(2).eq.0) then
   ! incompressible
   TANK_MK_INITIAL_PRESSURE = outflow_pressure
  elseif(fort_material_type(2).eq.TANK_MK_MATERIAL_TYPE) then
   ! compressible
   ! P = rho R_sp T = rho (gamma-1) U
   TANK_MK_INITIAL_PRESSURE = &
    fort_denconst(2)*(TANK_MK_GAS_GAMMA-one)*&
    fort_initial_temperature(2)*TANK_MK_GAS_CV
  else
   print *,"material type invalid for pressure setup!"
   stop
  endif
  
  return
 end subroutine INIT_CRYOGENIC_TANK_MK_MODULE

  ! LS>0 in the nozzle 
  ! x(2) is the vertical direction
 subroutine CRYOGENIC_TANK_MK_LS_NOZZLE(x,LS)
  use probcommon_module
  use global_utility_module
  IMPLICIT NONE

  real(amrex_real), INTENT(in) :: x(3)
  real(amrex_real), INTENT(out) :: LS
  real(amrex_real) :: xlo,xhi,ylo,yhi,xcen

  if ((TANK_MK_NOZZLE_RAD.gt.0.0d0).and. &
      (TANK_MK_NOZZLE_BASE.lt.0.0d0).and. &
      (TANK_MK_NOZZLE_HT.gt.0.0d0).and. &
      (TANK_MK_NOZZLE_THICK_OUTLET.gt.0.0d0)) then
   xlo=-TANK_MK_NOZZLE_RAD
   xhi=TANK_MK_NOZZLE_RAD
   xcen=zero
   ! "TANK_MK_HEIGHT" term insures nozzle is flush against the
   ! bottom of the tank.
   ylo=TANK_MK_NOZZLE_BASE-TANK_MK_HEIGHT
   yhi=TANK_MK_NOZZLE_BASE+TANK_MK_NOZZLE_HT
   if (SDIM.eq.2) then
    call squaredist(x(1),x(2),xlo,xhi,ylo,yhi,LS)
   else if (SDIM.eq.3) then
    call cylinderdist(x(1),x(3),x(2),xcen,xcen,xhi,ylo,yhi,LS)
   else
    print *,"sdim invalid"
    stop
   endif

   LS=-LS !now, LS>0 in the nozzle.
  else
   print *,"CRYOGENIC_TANK_MK_LS_NOZZLE parameter problem"
   stop
  endif

  return
 end subroutine CRYOGENIC_TANK_MK_LS_NOZZLE

 ! fluids tessellate the domain, solids are immersed. 
 ! fluid interfaces are extended into solids.
 ! material 1 is liquid
 ! material 2 is gas
 ! material 3 is solid

 ! MK = Validation of two-phase CFD models for propellant tank 
 !    self-pressurization: Crossing fluid types, scales, and gravity levels 
 !  Mohammad Kassemi , Olga Kartuzova, Sonya Hylton
 !
 ! xphys(2) is the vertical direction. 
 ! prior to this call (3D):
 ! x3D(1)=x(1)
 ! x3D(2)=x(3)
 ! x3D(3)=x(2)


subroutine rigid_displacement(xfoot,t,xphys,velphys)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real), INTENT(out) :: xfoot(3)
 real(amrex_real) :: xfoot_translate(3)
 real(amrex_real), INTENT(in) :: t
 real(amrex_real), INTENT(in) :: xphys(3)
 real(amrex_real), INTENT(out) :: velphys(3)

 real(amrex_real) :: xdisp_amplitude   ! xblob7
 real(amrex_real) :: xdisp_freq        ! yblob7
 real(amrex_real) :: xdisp_angV        ! zblob7
                             ! radblob7 prob type
                             ! radblob8 rot axis
 integer :: rotx,roty,rotz
 integer :: rot_dir
 integer :: dir
 real(amrex_real)    :: sign_term
 real(amrex_real)    :: Q(2,2)

 if (t.ge.0.0d0) then
  !do nothing
 else
  print *,"t invalid: ",t
  stop
 endif
 
 do dir=1,3                  ! bug fixed        
  xfoot(dir)=xphys(dir) 
  xfoot_translate(dir)=xphys(dir) 
  velphys(dir)=zero
 enddo

 if(NINT(radblob7).eq.0) then
  !do nothing
 elseif(NINT(radblob7).eq.1) then
  xdisp_amplitude=xblob7   
  xdisp_freq=yblob7

  xfoot(1)=xphys(1)-  &
     xdisp_amplitude*sin(xdisp_freq*t)

   !xphys=xfoot+A*sin(wt)
   !velphys=Aw cos(wt)
  velphys(1)=  &
     xdisp_amplitude*xdisp_freq*cos(xdisp_freq*t)
 elseif(NINT(radblob7).eq.2)then         
  !angular velocity zblob7 
  xdisp_angV=zblob7 

   ! xphys(2) is the vertical direction.
  rot_dir=NINT(radblob8)

  if(SDIM.eq.3)then
   if(rot_dir.eq.3)then
    rotx=1
    roty=2
    rotz=rot_dir
   elseif(rot_dir.eq.1)then
    rotx=2
    roty=3
    rotz=rot_dir
   elseif(rot_dir.eq.2)then !axis is vertical (y) direction.
    rotx=1
    roty=3
    rotz=rot_dir
   else
    print *,"rot_dir invalid"
    stop
   endif

    !xfoot=Q(t) xphys
    !xphys=Q(t)^{T} xfoot
    !velphys=Q'(t)^{T} xfoot = Q'(t)^{T} Q(t) xphys
    !
    !Q= cos(wt)              -sign_term sin(wt)
    !   sign_term*sin(wt)    cos(wt)
    !
    !Q'=w ( -sin(wt)          -sign_term cos(wt)
    !       sign_term*cos(wt) -sin(wt)  )
    !
    !Q'^{T}=w ( -sin(wt)           sign_term cos(wt)
    !           -sign_term*cos(wt) -sin(wt)  )
    !
    !Q'^{T}Q=w ( 0           sign_term 
    !            -sign_term   0 )

   sign_term=(-1.0d0)**(mod(rot_dir,2))
   Q(1,1)=cos(xdisp_angV*t)
   Q(1,2)=-sign_term*sin(xdisp_angV*t)
   Q(2,1)=-Q(1,2)
   Q(2,2)=Q(1,1)

   xfoot(rotx)=Q(1,1)*xphys(rotx)+Q(1,2)*xphys(roty)
   xfoot(roty)=Q(2,1)*xphys(rotx)+Q(2,2)*xphys(roty) 
   xfoot(rotz)=xphys(rotz)

   velphys(rotx)=sign_term*xdisp_angV*xphys(roty)
   velphys(roty)=-sign_term*xdisp_angV*xphys(rotx)
   velphys(rotz)=0.0d0
  elseif(SDIM.eq.2)then
   rotx=1
   roty=2 ! vertical direction

   Q(1,1)=cos(xdisp_angV*t)
   Q(1,2)=sin(xdisp_angV*t)
   Q(2,1)=-Q(1,2)
   Q(2,2)=Q(1,1)

   xfoot(rotx)= Q(1,1)*xphys(rotx)+Q(1,2)*xphys(roty)
   xfoot(roty)=Q(2,1)*xphys(rotx)+Q(2,2)*xphys(roty) 
   xfoot(rotz)=xphys(rotz)

   velphys(rotx)=-xdisp_angV*xphys(roty)
   velphys(roty)=xdisp_angV*xphys(rotx)
   velphys(rotz)=0.0d0
  else
   print *,"SDIM invalid"
   stop
  endif
 elseif(NINT(radblob7).eq.3) then
  !angular velocity zblob7 
  xdisp_angV=zblob7 
  xdisp_amplitude=xblob7   
  xdisp_freq=yblob7
  xfoot_translate(1)=xphys(1)-  &
     xdisp_amplitude*sin(xdisp_freq*t)

   ! xphys(2) is the vertical direction.
  rot_dir=NINT(radblob8)

  if (SDIM.eq.3) then
   !do nothing
  else
   print *,"expecting 3d"
   stop
  endif
  if (rot_dir.eq.2) then
   rotx=1
   roty=3
   rotz=rot_dir
  else
   print *,"expecting rot_dir==2"
   stop
  endif

  sign_term=(-1.0d0)**(mod(rot_dir,2))
  Q(1,1)=cos(xdisp_angV*t)
  Q(1,2)=-sign_term*sin(xdisp_angV*t)
  Q(2,1)=-Q(1,2)
  Q(2,2)=Q(1,1)

  xfoot(rotx)=Q(1,1)*xfoot_translate(rotx)+Q(1,2)*xfoot_translate(roty)
  xfoot(roty)=Q(2,1)*xfoot_translate(rotx)+Q(2,2)*xfoot_translate(roty) 
  xfoot(rotz)=xfoot_translate(rotz)

  velphys(rotx)=sign_term*xdisp_angV*xfoot_translate(roty)
  velphys(roty)=-sign_term*xdisp_angV*xfoot_translate(rotx)
  velphys(rotz)=0.0d0

   !xphys=xfoot+A*sin(wt)
   !velphys=Aw cos(wt)
  velphys(1)=velphys(1)+  &
     xdisp_amplitude*xdisp_freq*cos(xdisp_freq*t)

 else
  print *,"radblob7 invalid (rigid motion type invalid): ",radblob7
  stop
 endif

 end subroutine rigid_displacement

  ! swap y and z in 3D.
  ! no change in 2D.
  ! x(SDIM) is the vertical direction for gravity
  ! x3D(2) is the vertical direction for the geometry files.
 subroutine convert_to_x3D(x,x3D)
 IMPLICIT NONE

 real(amrex_real), INTENT(in) :: x(SDIM)
 real(amrex_real), INTENT(out) :: x3D(3)

 if ((dir_x.eq.1).and.(dir_z.eq.SDIM)) then
  x3D(1)=x(dir_x)
  x3D(2)=x(dir_z)
 else
  print *,"expecting dir_x=1 and dir_z=sdim: ",dir_x,dir_z
  stop
 endif

 if (SDIM.eq.2) then
  x3D(3)=0.0d0
 else if (SDIM.eq.3) then

  if (dir_y.eq.2) then
   x3D(3)=x(dir_y)
  else
   print *,"expecting dir_y=2: ",dir_y
   stop
  endif

 else
  print *,"dimension bust"
  stop
 endif

 end subroutine convert_to_x3D

  !x(SDIM) is the vertical direction of gravity.
 subroutine CRYOGENIC_TANK_MK_LS(x,t,LS,nmat)
  use probcommon_module
  use global_utility_module
  IMPLICIT NONE

  integer, INTENT(in) :: nmat
  real(amrex_real), INTENT(in) :: x(SDIM)
  real(amrex_real), INTENT(in) :: t
  real(amrex_real), INTENT(out) :: LS(nmat)
  real(amrex_real) :: nozzle_dist,LS_A
  integer :: called_from_heater_source

  real(amrex_real) :: x3D(3)
  real(amrex_real) :: xfoot3D(3)
  real(amrex_real) :: xvel(3)
  integer auxcomp
  integer num_cyl
  real(amrex_real) :: LS_heater_a
  real(amrex_real) :: LS_heater_b
  real(amrex_real) :: LS_tank
  real(amrex_real) :: LS_nozzle
  real(amrex_real) :: LS_LAD_housing

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif
  if (t.ge.0.0d0) then
   ! do nothing
  else
   print *,"t invalid"
   stop
  endif

   !x(SDIM) is the vertical direction of gravity.
   !x3D(2) is the vertical direction of geometry.
   !in 3D:
   ! x3D(1)=x(1)
   ! x3D(2)=x(3)
   ! x3D(3)=x(2)
   !in 2D:
   ! x3D(1)=x(1)
   ! x3D(2)=x(2)
   ! x3D(3)=0.0
  call convert_to_x3D(x,x3D)

  call rigid_displacement( &
   xfoot3D, & !intent(out)
   t, & !intent(in)
   x3D, & !intent(in)
   xvel) !intent(out)

   ! material 1= liquid  (e.g. Freon 113)
   ! material 2= vapor  
   ! material 3= tank geometry (e.g. acrylic)
  if ((num_materials.eq.3).and.(probtype.eq.423)) then
   ! liquid
   !TANK_MK_INTERFACE_RADIUS = radblob2
   if (TANK_MK_INTERFACE_RADIUS.eq.0.0d0) then
    LS(1)=TANK_MK_INTERFACE_LOCATION-xfoot3D(2)
   else if (TANK_MK_INTERFACE_RADIUS.gt.0.0d0) then
    if (SDIM.eq.2) then
     LS(1)=sqrt((xfoot3D(1)-TANK_MK_BUBBLE_X)**2+&
                (xfoot3D(2)-TANK_MK_BUBBLE_Y)**2)&
               -TANK_MK_INTERFACE_RADIUS
    else if (SDIM.eq.3) then 
     !TANK_MK_BUBBLE_X         = xblob2
     LS(1)=sqrt((xfoot3D(1)-TANK_MK_BUBBLE_X)**2+&
                (xfoot3D(3)-TANK_MK_BUBBLE_Y)**2+&
                (xfoot3D(2)-TANK_MK_BUBBLE_Z)**2)&
               -TANK_MK_INTERFACE_RADIUS
    else
     print *,"dimension bust"
     stop
    endif
   else
    print *,"radblob2 invalid: ",radblob2
    stop
   endif 
   LS(2)=-LS(1)

   if (1.eq.0) then
    print *,"x,t,LS(1),LS(2) ",x,t,LS(1),LS(2)
   endif

   ! Solid

   if ((axis_dir.eq.0).or. & !ZBOT
       (axis_dir.eq.1)) then !TPCE

    if (FSI_flag(3).eq.FSI_PRESCRIBED_PROBF90) then
     ! do nothing
    else
     print *,"expecting FSI_flag(3).eq.FSI_PRESCRIBED_PROBF90"
     print *,"(CRYOGENIC_TANK_MK)"
     print *,"FSI_flag(3)=",FSI_flag(3)
     stop
    endif

  !R=abs(sqrt(P(1)**2+P(3)**2))
  !Z=abs(P(2))
    LS(3)=SOLID_TOP_HALF_DIST(xfoot3D)

    if (axis_dir.eq.0) then !ZBOT
     ! do nothing
    else if (axis_dir.eq.1) then ! TPCE
     call CRYOGENIC_TANK_MK_LS_NOZZLE(xfoot3D,nozzle_dist)
     if (nozzle_dist.gt.LS(3)) then ! nozzle_dist>0 in the nozzle
      LS(3)=nozzle_dist
     endif
     called_from_heater_source=0
     !LS_A>0 in heater
     call CRYOGENIC_TANK_MK_LS_HEATER_A(xfoot3D,LS_A,called_from_heater_source) 
     if (LS_A.gt.LS(3)) then
      LS(3)=LS_A
     endif
    else
     print *,"axis_dir invalid: ",axis_dir
     stop
    endif
   else if (axis_dir.eq.2) then !TPCE e.g. tpce_geometry.vtk

    if (FSI_flag(3).eq.FSI_PRESCRIBED_PROBF90) then
     ! do nothing
    else
     print *,"expecting FSI_flag(3).eq.FSI_PRESCRIBED_PROBF90"
     print *,"(CRYOGENIC_TANK_MK)"
     print *,"FSI_flag(3)=",FSI_flag(3)
     stop
    endif
    auxcomp=1
    call interp_from_aux_grid(auxcomp,xfoot3D,LS_heater_a)
    LS(3)=LS_heater_a
    auxcomp=2
    call interp_from_aux_grid(auxcomp,xfoot3D,LS_heater_b)
    LS(3)=max(LS(3),LS_heater_b)
    auxcomp=5
    call interp_from_aux_grid(auxcomp,xfoot3D,LS_tank)
    LS(3)=max(LS(3),LS_tank)
    auxcomp=6
    call interp_from_aux_grid(auxcomp,xfoot3D,LS_nozzle)
    LS(3)=max(LS(3),LS_nozzle)
    if (TANK_MK_AUX_THICK_WALLS.eq.1) then
     ! do nothing
    else if (TANK_MK_AUX_THICK_WALLS.eq.0) then
     auxcomp=7
     call interp_from_aux_grid(auxcomp,xfoot3D,LS_LAD_housing)
     LS(3)=max(LS(3),LS_LAD_housing)
    else 
     print *,"TANK_MK_AUX_THICK_WALLS invalid"
     stop
    endif

    !Experimental characterization of non-isothermal sloshing in microgravity
    !Monteiro et al 2024
   else if (axis_dir.eq.3) then
    if (SDIM.eq.3) then
     num_cyl=4
    else if (SDIM.eq.2) then
     num_cyl=1
    else
     print *,"dimension problem"
     stop
    endif


   else
    print *,"axis_dir invalid: ",axis_dir
    stop
   endif

  else
   print *,"num_materials ", num_materials
   print *,"probtype ", probtype
   print *,"num_materials or probtype invalid"
   stop
  endif
  ! print*,"X= ",x," LS= ", LS
  return
 end subroutine CRYOGENIC_TANK_MK_LS

   !LS>0 in heater A region.
   !x(2) is the vertical direction
 subroutine CRYOGENIC_TANK_MK_LS_HEATER_A(x,LS,called_from_heater_source)
  use probcommon_module
  use global_utility_module
  IMPLICIT NONE

  real(amrex_real), INTENT(in) :: x(3)
  real(amrex_real), INTENT(out) :: LS
  integer, INTENT(in) :: called_from_heater_source
  real(amrex_real) zdiff
  real(amrex_real) angle_end_center
  real(amrex_real) shell_R
  real(amrex_real) shell_center
  real(amrex_real) r_crit
  real(amrex_real) z_crit
  real(amrex_real) r_cyl

  zdiff=x(2)-TANK_MK_END_CENTER

  if (SDIM.eq.2) then
   r_cyl=abs(x(1))
  else if (SDIM.eq.3) then
   r_cyl=sqrt(x(1)**2+x(3)**2)
  else
   print *,"sdim invalid"
   stop
  endif

  if (zdiff.le.0.0d0) then
   LS=-TANK_MK_END_RADIUS
  else if (zdiff.gt.0.0d0) then
   angle_end_center=atan(r_cyl/zdiff)
   if ((angle_end_center.ge.0.0d0).and. &
       (angle_end_center.lt.0.5d0*Pi)) then
    shell_R=sqrt(r_cyl**2+zdiff**2)
    shell_center=TANK_MK_END_RADIUS-0.5d0*TANK_MK_HEATER_THICK

    if (called_from_heater_source.eq.1) then
     LS=0.5d0*TANK_MK_HEATER_THICK-abs(shell_R-shell_center)
    else if (called_from_heater_source.eq.0) then
     LS=shell_R-TANK_MK_END_RADIUS+TANK_MK_HEATER_THICK
    else
     print *,"called_from_heater_source invalid"
     stop
    endif

    if ((TANK_MK_HEATER_TOTAL_ANGLE.gt.0.0d0).and. &
        (TANK_MK_HEATER_TOTAL_ANGLE.le.0.5d0*Pi)) then
     if (angle_end_center.le.TANK_MK_HEATER_TOTAL_ANGLE) then
      ! do nothing
     else if (angle_end_center.ge.TANK_MK_HEATER_TOTAL_ANGLE) then
      r_crit=sin(TANK_MK_HEATER_TOTAL_ANGLE)*(TANK_MK_END_RADIUS- &
            0.5d0*TANK_MK_HEATER_THICK)
      z_crit=cos(TANK_MK_HEATER_TOTAL_ANGLE)*(TANK_MK_END_RADIUS- &
            0.5d0*TANK_MK_HEATER_THICK)
      if ((r_crit.gt.0.0d0).and.(z_crit.gt.0.0d0).and. &
          (r_crit.le.TANK_MK_END_RADIUS).and. &
          (z_crit.le.TANK_MK_END_RADIUS)) then
       LS=-sqrt( (r_cyl-r_crit)**2+(zdiff-z_crit)**2 )
      else
       print *,"r_crit or z_crit invalid"
       stop
      endif
     else
      print *,"angle_end_center is NaN"
      stop
     endif 
    else
     print *,"TANK_MK_HEATER_TOTAL_ANGLE invalid"
     stop
    endif
   else
    print *,"angle_end_center invalid"
    stop
   endif
  else
   print *,"zdiff is NaN"
   stop
  endif

  return
 end subroutine CRYOGENIC_TANK_MK_LS_HEATER_A

 ! if SOLID VELOCITY requested everywhere (including outside of the solid),
 ! then velsolid==1
 ! This routine called from:
 ! fort_initvelocity (velsolid_flag=0, time=0)
 ! velsolid (velsolid_flag=1, VEL=0 by default)
 ! CRYOGENIC_TANK_MK_VEL_BC (velsolid_flag=0)
 !
 ! x(SDIM) is the vertical direction.
 subroutine CRYOGENIC_TANK_MK_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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

  real(amrex_real) :: x3D(3)
  real(amrex_real) :: xfoot3D(3)
  real(amrex_real) :: xvel(3)
  integer, parameter :: impulsive_tank_forcing=0

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  ! swap y and z in 3D.
  ! no change in 2D.
  ! x(SDIM) is the vertical direction for gravity
  ! x3D(2) is the vertical direction for the geometry files.
  call convert_to_x3D(x,x3D)

   ! if w_dot<>0 then
   !  w=t * w_dot
   !  Omega=w_dot * t^2/2
   !  Assume Satellite frame of reference,
   !  w=t*w_dot  0<t<t_cut
   !  w=t_cut*w_dt t>t_cut
   !  x_centrifugal=x_tank
   !  y_centrifugal=y_tank+lever_arm 
  call rigid_displacement(xfoot3D,t,x3D,xvel)

  if ((velsolid_flag.eq.0).or. &
      (velsolid_flag.eq.1)) then
   ! do nothing
  else 
   print *,"velsolid_flag invalid"
   stop
  endif

  if ((num_materials.eq.3).and.(probtype.eq.423)) then

   if ((axis_dir.eq.0).or. &  !ZBOT problem
       (axis_dir.eq.1).or. &  !TPCE problem
       (axis_dir.eq.2)) then  !read e.g. tpce_geometry.vtk (TPCE problem)

    do dir=1,SDIM
     VEL(dir)=0.0d0
    enddo

    ! if SOLID VELOCITY requested everywhere (including outside of the solid),
    ! then velsolid==1

    if (((t.eq.0.0d0).or. &           ! called from fort_initvelocity
         (velsolid_flag.eq.0)).and. & ! called from boundary condition routine
        (LS(3).lt.zero).and. &
        (impulsive_tank_forcing.eq.1)) then

     do dir=1,SDIM
      VEL(dir)=0.0d0
     enddo

    else if (((t.gt.0.0d0).and. &
              (velsolid_flag.eq.1)).or. &
             (LS(3).ge.zero).or. &
             (impulsive_tank_forcing.eq.0)) then

 !   VEL(dir_x)=xvel(dir_x)
     if (SDIM.eq.2) then
      VEL(1)=xvel(1)
      VEL(2)=xvel(2)
     else if (SDIM.eq.3) then
      VEL(1)=xvel(1)
      VEL(SDIM)=xvel(2) !SDIM=vertical dir of grav.  2=vertical of geom files
      VEL(2)=xvel(3)
     else
      print *,"dimension bust"
      stop
     endif

    else
     print *,"LS(3) or velsolid_flag bust"
     print *,"LS(3)=",LS(3)
     print *,"velsolid_flag=",velsolid_flag
     print *,"t=",t
     print *,"impulsive_tank_forcing =",impulsive_tank_forcing
     stop
    endif

   else
    print *,"axis_dir invalid: ",axis_dir
    stop
   endif

  else
   print *,"num_materials ", num_materials
   print *,"probtype ", probtype
   print *,"num_materials or probtype invalid"
   stop
  endif

  return 
 end subroutine CRYOGENIC_TANK_MK_VEL

!P(2) is the vertical direction
real(amrex_real) function SOLID_TOP_HALF_DIST(P)
 ! Returns the signed distance function to the
 ! cylindrical tank with spherical ends.
 ! The tank is symmetrical to x(SDIM)=0;
 ! The axis of cylinder is along dim=SDIM direction
 ! Inside the tank < 0
 ! Outside the tank > 0
 implicit none

 real(amrex_real), INTENT(in), dimension(3) :: P
 real(amrex_real) R,Z,FRZ,D1,D2
 
 if (SDIM.eq.2) then
  R=abs(P(1))
  Z=abs(P(2))
 elseif (SDIM.eq.3) then
  R=abs(sqrt(P(1)**2+P(3)**2))
  Z=abs(P(2))
 else
  print *,"Dimension bust at DIST_FINITE_CYLINDER"
  stop
 endif

  ! R=Z=0 is the dead center of the tank; computation domain goes from
  ! -ZTOTAL/2, ... ,ZTOTAL/2
  ! TANK_MK_END_RADIUS       = xblob4
 if (TANK_MK_END_RADIUS.eq.0.0d0) then ! rectangular tank

  !TANK_MK_RADIUS             = xblob
  !TANK_MK_HEIGHT             = yblob
  if (Z.le.TANK_MK_HEIGHT/2.0) then
   if (R.ge.TANK_MK_RADIUS) then
    SOLID_TOP_HALF_DIST=R-TANK_MK_RADIUS
   else
    D1 = TANK_MK_RADIUS - R
    D2 = TANK_MK_HEIGHT/2.0d0-Z
    if (D1.lt.D2) then
     SOLID_TOP_HALF_DIST=-D1
    else
     SOLID_TOP_HALF_DIST=-D2
    endif
   endif
  else 
   if (R.le.TANK_MK_RADIUS) then
    SOLID_TOP_HALF_DIST=Z-TANK_MK_HEIGHT/2.0d0
   else
    SOLID_TOP_HALF_DIST= &
     sqrt( (Z-TANK_MK_HEIGHT/2.0d0)**2+ &
           (R-TANK_MK_RADIUS)**2 )
   endif
  endif
    
  ! TANK_MK_END_RADIUS       = xblob4
 else if (TANK_MK_END_RADIUS.gt.0.0d0) then
 
  ! Equation of the line passing through
  ! (0,TANK_MK_END_CENTER)
  !          and
  ! (TANK_MK_RADIUS,TANK_MK_HEIGHT/2)
  ! TANK_MK_END_CENTER       = yblob4
  ! TANK_MK_RADIUS             = xblob
  FRZ=R*(TANK_MK_END_CENTER-TANK_MK_HEIGHT/2.0d0)/TANK_MK_RADIUS &
     + Z - TANK_MK_END_CENTER

  if (FRZ.le.0.0d0) then
   ! Below line
   if (R.le.TANK_MK_RADIUS) then
    ! In the tank
    D1 = TANK_MK_RADIUS - R
    D2 = sqrt((R-TANK_MK_RADIUS)**2 + (Z-TANK_MK_HEIGHT/2)**2)
    SOLID_TOP_HALF_DIST = -min(D1,D2)
   elseif (R.gt.TANK_MK_RADIUS) then
    ! Out of the tank
    if(Z.le.TANK_MK_HEIGHT/2.0d0) then
     ! Below the cap base line
     SOLID_TOP_HALF_DIST = R - TANK_MK_RADIUS
    elseif (Z.gt.TANK_MK_HEIGHT/2) then
     ! Above the cap baseline
     !TANK_MK_RADIUS             = xblob
     !TANK_MK_HEIGHT             = yblob
     SOLID_TOP_HALF_DIST = &
      sqrt((R-TANK_MK_RADIUS)**2 + (Z-TANK_MK_HEIGHT/2.0d0)**2)
    else
     print *,"Z invalid!"
     stop
    endif ! Z
   else
    print *,"Z invalid!"
    stop
   endif ! R
  elseif (FRZ.gt.0.0d0) then
   ! Above line
   ! Distance to curvature center
   D2 = sqrt(R**2 + (Z-TANK_MK_END_CENTER)**2)
   if (D2.le.TANK_MK_END_RADIUS) then
    ! Inside tank
    if (Z.le.TANK_MK_HEIGHT/2) then
     ! Below the cap baseline
     D1 = TANK_MK_RADIUS - R
    elseif (Z.gt.TANK_MK_HEIGHT/2) then
     ! Above the cap baseline
     D1 = sqrt((R-TANK_MK_RADIUS)**2 + (Z-TANK_MK_HEIGHT/2.0d0)**2)
    else
     print *,"Z invalid!"
     stop
    endif ! Z
    SOLID_TOP_HALF_DIST = -min(TANK_MK_END_RADIUS-D2,D1)
   elseif (D2.gt.TANK_MK_END_RADIUS) then
    SOLID_TOP_HALF_DIST = D2-TANK_MK_END_RADIUS
   else
    print *,"Distance to curvature center invalid!"
    stop
   end if ! D2
  
  else
   print *,"Line equation invalid!",FRZ
   stop
  end if ! FRZ
 else 
  print *,"TANK_MK_END_RADIUS invalid: ",TANK_MK_END_RADIUS
  stop
 endif

end function SOLID_TOP_HALF_DIST

!***********************************************
! compressible material functions for 
!   (ns.material_type = TANK_MK_MATERIAL_TYPE)
! C_spc => Specific heat capacity [J(kg K)]
! C_m   =>    Molar heat capacity [J(mol K)]
!
! U = C_{v,spc} T
! [U] = J/kg= J/(kg K)  K
! R_spc = C_{p,scp}-C_{v,scp}
! R_unv = C_{p,m} - C_{v,m}
! gamma = C_{p,scp}/C_{v,scp}
!
! p = rho R_spc T 
!   = rho R_spc x U/C_{v,spc}
!   = rho (C_{p,scp}-C_{v,scp})/C_{v,spc} U
!   = rhp (gamma-1) U
!
! a = sqrt(gamma R_sp T) = sqrt(gamma p/rho)
!
! EOS=Equation Of State P=P(rho,e,Y)
subroutine EOS_CRYOGENIC_TANK_MK(rho,massfrac_var, &
  internal_energy,pressure, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(in) :: internal_energy
 real(amrex_real), INTENT(out) :: pressure
 integer :: dummy_input

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then  ! material 2 is vapor
   if (imattype.eq.TANK_MK_MATERIAL_TYPE) then
    ! p = rho (gamma-1) U
    pressure=rho * (TANK_MK_GAS_GAMMA-one) * internal_energy
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid EOS_CRYOGENIC_TANK_MK"
    print *,"(breakpoint) break point and gdb: "
    print *,"(1) compile with the -g option"
    print *,"(2) break CRYOGENIC_TANK_MK.F90:350"
    print *,"By pressing <CTRL C> during this read statement, the"
    print *,"gdb debugger will produce a stacktrace."
    print *,"type 0 then <enter> to exit the program"
    read *,dummy_input
    error stop
   endif
  else
   call EOS_material_CORE(rho,massfrac_var, &
         internal_energy,pressure,imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine EOS_CRYOGENIC_TANK_MK

subroutine dVdT_CRYOGENIC_TANK_MK(dVdT,massfrac_var, &
   pressure,temperature, &
   imattype,im,num_species_var_in)
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: imattype,im,num_species_var_in
real(amrex_real), INTENT(in) :: pressure
real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
real(amrex_real), INTENT(in) :: temperature
real(amrex_real), INTENT(out) :: dVdT
integer :: dummy_input

 if (pressure.gt.0.0d0) then
  ! do nothing
 else
  print *,"pressure invalid dVdT_CRYOGENIC_TANK_MK:",pressure
  stop
 endif
 if (temperature.gt.0.0d0) then
  ! do nothing
 else
  print *,"temperature invalid dVdT_CRYOGENIC_TANK_MK:",temperature
  stop
 endif

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then  ! vapor material
   if (imattype.eq.TANK_MK_MATERIAL_TYPE) then
    ! p = rho (gamma-1) e  rho=1/V
    ! V = (gamma-1) e/p=(gamma-1)Cv T/p
    ! dVdT=(gamma-1)Cv/p
    dVdT=(TANK_MK_GAS_GAMMA-one) * TANK_MK_GAS_CV/pressure
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid dVdT_CRYOGENIC_TANK_MK"
    print *,"(breakpoint) break point and gdb: "
    print *,"(1) compile with the -g option"
    print *,"(2) break CRYOGENIC_TANK_MK.F90:350"
    print *,"By pressing <CTRL C> during this read statement, the"
    print *,"gdb debugger will produce a stacktrace."
    print *,"type 0 then <enter> to exit the program"
    read *,dummy_input
    error stop
   endif
  else
   call dVdT_material_CORE(dVdT,massfrac_var, &
         pressure,temperature, &
         imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine dVdT_CRYOGENIC_TANK_MK


subroutine SOUNDSQR_CRYOGENIC_TANK_MK(rho,massfrac_var, &
  internal_energy,soundsqr, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(in) :: internal_energy
 real(amrex_real), INTENT(out) :: soundsqr
 real(amrex_real) pressure

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then ! vapor
   if (imattype.eq.TANK_MK_MATERIAL_TYPE) then
     ! a = sqrt(gamma R_sp T) = sqrt(gamma p/rho)
    call EOS_CRYOGENIC_TANK_MK(rho,massfrac_var, &
     internal_energy,pressure,imattype,im,num_species_var_in)
    if (rho.gt.0.0d0) then
     soundsqr=TANK_MK_GAS_GAMMA*pressure/rho
    else
     print *,"rho invalid: ",rho
     stop
    endif
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid SOUNDSQR CRYOGENIC TANK_MK"
    stop
   endif
  else
   call SOUNDSQR_material_CORE(rho,massfrac_var, &
    internal_energy,soundsqr, &
    imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine SOUNDSQR_CRYOGENIC_TANK_MK

! this routine returns: e(T)
! CV for the temperature solver is defined as:
! (e(T+DT)-e(T))/DT
subroutine INTERNAL_CRYOGENIC_TANK_MK(rho,massfrac_var, &
  temperature,local_internal_energy, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(in) :: temperature 
 real(amrex_real), INTENT(out) :: local_internal_energy

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then ! vapor
   if ((imattype.eq.TANK_MK_MATERIAL_TYPE).or. &
       (imattype.eq.0)) then 
    ! U_mix = C_{v,spc} T
    local_internal_energy=TANK_MK_GAS_CV*temperature
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid INTERNAL CRYOGENIC TANK_MK"
    stop
   endif
  else
   call INTERNAL_material_CORE(rho,massfrac_var, &
    temperature,local_internal_energy, &
    imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine INTERNAL_CRYOGENIC_TANK_MK

subroutine TEMPERATURE_CRYOGENIC_TANK_MK(rho,massfrac_var, &
  temperature,internal_energy, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(out) :: temperature 
 real(amrex_real), INTENT(in) :: internal_energy

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then  ! vapor
   if ((imattype.eq.TANK_MK_MATERIAL_TYPE).or. &
       (imattype.eq.0)) then 
    ! T = U / C_{v,spc}
    if (TANK_MK_GAS_CV.gt.0.0d0) then
     temperature=internal_energy/TANK_MK_GAS_CV
    else
     print *,"TANK_MK_GAS_CV invalid 1"
     stop
    endif
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid TEMPERATURE_CRYOGENIC_TANK_MK"
    stop
   endif
  else
   call TEMPERATURE_material_CORE(rho,massfrac_var, &
     temperature,internal_energy, &
     imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine TEMPERATURE_CRYOGENIC_TANK_MK

!***********************************************
! called by the boundary condition routine
! might be called at initialization, so put a placeholder pressure here.
!
!x(SDIM) is the vertical direction of gravity.
subroutine CRYOGENIC_TANK_MK_PRES_UTIL(x,PRES,rho_hyd)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(out) :: PRES
real(amrex_real), INTENT(out) :: rho_hyd
integer simple_hyd_p
integer gravity_dir

simple_hyd_p=1

call fort_derive_gravity_dir(gravity_vector,gravity_dir)

!PRES=TANK_MK_INITIAL_GAS_PRESSURE
if(fort_material_type(2).eq.0) then

 if (simple_hyd_p.eq.0) then
  rho_hyd=fort_denconst(2)
  ! incompressible gas
  ! Flat open top x_2: TANK_MK_HEIGHT/two
  ! Known pressure(P_1) at top (outflow_pressure)
  ! P_2=P_1 + rho*g*(z_1-z_2)  [g>0]
  if (x(SDIM).ge.TANK_MK_INTERFACE_LOCATION) then
   PRES=TANK_MK_INITIAL_PRESSURE+&
       fort_denconst(2)*(TANK_MK_HEIGHT/two-x(SDIM))* &
       (abs(gravity_vector(gravity_dir))) 
  elseif (x(SDIM).lt.TANK_MK_INTERFACE_LOCATION) then
   PRES=TANK_MK_INITIAL_PRESSURE+&
       fort_denconst(2)*(TANK_MK_HEIGHT/two-TANK_MK_INTERFACE_LOCATION)* &
       (abs(gravity_vector(gravity_dir)))+ &
       fort_denconst(1)*(TANK_MK_INTERFACE_LOCATION-x(SDIM))* &
       (abs(gravity_vector(gravity_dir)))
  else
   print *,"x(SDIM) is invalid in CRYOGENIC_TANK_MK_PRES!"
   stop
  endif

 else if (simple_hyd_p.eq.1) then
  rho_hyd=fort_denconst(1)
  PRES=-abs(gravity_vector(gravity_dir))*rho_hyd*(x(SDIM)-probhiy-probhiy)
 else
  print *,"simple_hyd_p invalid"
  stop
 endif
elseif (fort_material_type(2).eq.TANK_MK_MATERIAL_TYPE) then
 ! compressible gas
 ! Known pressure(P_1) at top (based on given density and temperature)
 ! P_2=P_1 * exp(g*(z_1-z_2)/(R_sp*T_0))  [g>0]
 rho_hyd=fort_denconst(2)
 if (x(SDIM).ge.TANK_MK_INTERFACE_LOCATION) then
  PRES=TANK_MK_INITIAL_PRESSURE*&
       exp((TANK_MK_END_CENTER+TANK_MK_END_RADIUS-x(SDIM))* &
       abs(gravity_vector(gravity_dir))/&
           (TANK_MK_R_UNIV/fort_molar_mass(2)*fort_initial_temperature(2)))
 elseif (x(SDIM).lt.TANK_MK_INTERFACE_LOCATION) then
  PRES=TANK_MK_INITIAL_PRESSURE*&
       exp((TANK_MK_END_CENTER+TANK_MK_END_RADIUS-TANK_MK_INTERFACE_LOCATION)*&
            abs(gravity_vector(gravity_dir))/&
           (TANK_MK_R_UNIV/fort_molar_mass(2)*fort_initial_temperature(2)))+&
       fort_denconst(1)*(TANK_MK_INTERFACE_LOCATION-x(SDIM))* &
                        (abs(gravity_vector(gravity_dir)))
 else
  print *,"x(SDIM) is invalid in CRYOGENIC_TANK_MK_PRES!"
  stop
 endif
else
 print  *,"invalid material type in pressure setup!"
 stop
endif

return 
end subroutine CRYOGENIC_TANK_MK_PRES_UTIL

!x(SDIM) is the vertical direction for gravity.
subroutine CRYOGENIC_TANK_MK_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: PRES
real(amrex_real) :: rho_hyd

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

!x(SDIM) is the vertical direction for gravity.
call CRYOGENIC_TANK_MK_PRES_UTIL(x,PRES,rho_hyd)

return 
end subroutine CRYOGENIC_TANK_MK_PRES

!x(SDIM) is the vertical direction for gravity.
subroutine CRYOGENIC_TANK_MK_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
real(amrex_real) den,temperature,internal_energy,pressure,Pgamma
real(amrex_real) massfrac_parm(num_species_var+1)
real(amrex_real) LL

 ! num_state_material=2 (default)  density and temperature
 ! num_state_material>2 if scalar (species) variables added.

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

if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.423)) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  ! density
  if(im.eq.2) then
   if(fort_material_type(2).eq.0) then
    ! incompressible
    STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im)
   elseif(fort_material_type(2).eq.TANK_MK_MATERIAL_TYPE) then
    ! compressible
    ! rho =P/(R_sp T)
    call CRYOGENIC_TANK_MK_PRES(x,t,LS,pressure,nmat)
    STATE(ibase+ENUM_DENVAR+1) = pressure/&
     (TANK_MK_R_UNIV/fort_molar_mass(2)*fort_initial_temperature(2))
   else
    print *,"material type invalid for density setup!"
    stop
   endif ! material_type
  else if ((im.eq.1).or.(im.eq.3)) then ! liquid or tank walls
   STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im)
  else
   print *,"im invalid 758 ",im
   stop
  endif
  ! temperature
  if (t.eq.0.0d0) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im)
  else if (t.gt.0.0d0) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif
  ! species
  do n=1,num_species_var
   STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
  if (t.eq.0.0d0) then
   if (im.eq.2) then
    den=STATE(ibase+ENUM_DENVAR+1)
    temperature=STATE(ibase+ENUM_TEMPERATUREVAR+1)
    call init_massfrac_parm(den,massfrac_parm,im)
    do n=1,num_species_var
     massfrac_parm(n)=STATE(ibase+ENUM_SPECIESVAR+n)
    enddo
    call INTERNAL_CRYOGENIC_TANK_MK(den,massfrac_parm, &
     temperature,internal_energy,TANK_MK_MATERIAL_TYPE,im,num_species_var)
    call EOS_CRYOGENIC_TANK_MK(den,massfrac_parm, &
     internal_energy,pressure,TANK_MK_MATERIAL_TYPE,im,num_species_var)
    if (abs(TANK_MK_R_UNIV-fort_R_Palmore_Desjardins).le. &
            EPS2*TANK_MK_R_UNIV) then
     ! do nothing
    else
     print *,"mismatch between TANK_MK_R_UNIV and fort_R_Palmore_Desjardins"
     print *,"TANK_MK_R_UNIV=",TANK_MK_R_UNIV
     print *,"fort_R_Palmore_Desjardins=",fort_R_Palmore_Desjardins
     stop
    endif
    LL=get_user_latent_heat(1,room_temperature,1)

    if (LL.ne.zero) then
     call Pgamma_Clausius_Clapyron(Pgamma, &
            fort_reference_pressure(1), &
            temperature, &
            fort_saturation_temp(1), &
            LL, &
            TANK_MK_R_UNIV,fort_molar_mass(2))
     if (abs(Pgamma-pressure).le.EPS2*pressure) then
      ! do nothing
     else
      print *,"mismatch between Pgamma and Pgas"
      print *,"Pgamma=",Pgamma
      print *,"Pgas=",pressure
      print *,"reference pressure=",fort_reference_pressure(1)
      print *,"modify reference pressure to be: ", &
        fort_reference_pressure(1)*pressure/Pgamma
      print *,"im=",im
      print *,"fort_R_Palmore_Desjardins=",fort_R_Palmore_Desjardins
      print *,"den=",den
      print *,"temperature=",temperature
      print *,"fort_stiffCV(2)=",fort_stiffCV(2)
      print *,"TANK_MK_GAS_CV=",TANK_MK_GAS_CV
      print *,"internal_energy (e=cv T)=",internal_energy
      print *,"TANK_MK_GAS_CP=",TANK_MK_GAS_CP
      print *,"fort_stiffCP(2)=",fort_stiffCP(2)
      print *,"TANK_MK_GAS_GAMMA=",TANK_MK_GAS_GAMMA
      print *,"pressure=rho (gamma-1) e= ",pressure
      print *,"latent heat= ",LL
      print *,"fort_saturation_temp(1)=",fort_saturation_temp(1)
      print *,"fort_molar_mass(2)=",fort_molar_mass(2)
      print *,"P_Clausius=Pref*exp(-(L*W/R)*(1/Tgamma-1/TSAT))=",Pgamma
      stop
     endif
    else if (LL.eq.zero) then
     !do nothing
    else
     print *,"LL invalid: ",LL
     stop
    endif 

   else if ((im.eq.1).or.(im.eq.3)) then
    ! do nothing
   else
    print *,"im invalid 814 ",im
    stop
   endif
  else if (t.gt.0.0d0) then
   ! do nothing
  else
   print *,"t invalid"
   stop
  endif
 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine CRYOGENIC_TANK_MK_STATE

! dir=1..sdim  side=1..2
! xghost(SDIM)=vertical direction for gravity.
subroutine CRYOGENIC_TANK_MK_LS_BC(xwall,xghost,t,LS, &
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

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call CRYOGENIC_TANK_MK_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_LS_BC


! dir=1..sdim  side=1..2 veldir=1..sdim
! xghost(SDIM)=vertical direction for gravity.
subroutine CRYOGENIC_TANK_MK_VEL_BC(xwall,xghost,t,LS, &
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

 call CRYOGENIC_TANK_MK_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_VEL_BC

! xghost(SDIM)=vertical direction for gravity.
! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK_MK_PRES_BC(xwall,xghost,t,LS, &
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

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call CRYOGENIC_TANK_MK_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_PRES_BC

! xghost(SDIM)=vertical direction for gravity.
! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK_MK_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real) local_STATE(nmat*num_state_material)
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
 call CRYOGENIC_TANK_MK_STATE(xghost,t,LS,local_STATE, &
         local_bcflag,nmat,num_state_material)
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
end subroutine CRYOGENIC_TANK_MK_STATE_BC

! x(SDIM)=vertical direction for gravity.
subroutine CRYOGENIC_TANK_MK_HEATSOURCE( &
     im,VFRAC, &
     time, &
     x, &
     xsten, & ! xsten(-nhalf:nhalf,SDIM)
     nhalf, &
     temp, &
     heat_source, & ! unit of tildeQ not Q
     den,CV,dt, &
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

integer dir
real(amrex_real) local_dx(SDIM)
real(amrex_real) flux_magnitude

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

do dir=1,SDIM
 local_dx(dir)=xsten(1,dir)-xsten(-1,dir)
 if (local_dx(dir).gt.0.0d0) then
  ! do nothing
 else
  print *,"local_dx invalid"
  stop
 endif
enddo

flux_magnitude=0.0d0

if ((num_materials.eq.3).and.(probtype.eq.423)) then
 heat_source=0.0d0
else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_HEATSOURCE

! SDIM=vertical direction of gravity
subroutine CRYOGENIC_TANK_MK_SUMINT(GRID_DATA_IN,increment_out1, &
                increment_out2,nsum1,nsum2,isweep)
use probcommon_module_types
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nsum1,nsum2,isweep
type(user_defined_sum_int_type), INTENT(in) :: GRID_DATA_IN
real(amrex_real), INTENT(inout) :: increment_out1(nsum1)
real(amrex_real), INTENT(inout) :: increment_out2(nsum2)

real(amrex_real) massfrac_parm(num_species_var+1)
real(amrex_real) T1_probe(SDIM)  !ZBOT
real(amrex_real) T4_probe(SDIM)  !TPCE
integer im
integer dir
integer dencomp,local_ispec
real(amrex_real) den,temperature,internal_energy,pressure
real(amrex_real) support_r
real(amrex_real) dx_coarsest
real(amrex_real) charfn
real(amrex_real) volgrid
real(amrex_real) denom

integer :: level,finest_level

integer :: i,j,k
integer :: ilev

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid
level=GRID_DATA_IN%level
finest_level=GRID_DATA_IN%finest_level

if ((level.le.finest_level).and.(level.ge.0)) then
 ! do nothing
else
 print *,"level invalid"
 stop
endif

if ((num_materials.eq.3).and.(probtype.eq.423)) then

 if ((nsum1.eq.1).and.(nsum2.eq.2)) then
  ! integral of region surrounding T1 in Figure 3 of Barsi and Kassemi, 2013
  ! T1: r=0.0  Z=0.2921 relative to bottom of tank cylindrical section.
  ! (note the very bottom of the tank corresponds to z=-0.2032-0.0254)
  ! bottom of cylindrical section: -0.2032
  ! so T1_probe_z=-0.2032+0.2921=0.0889
  T1_probe(1)=0.0d0
  T1_probe(2)=0.0889d0
  T4_probe(1)=0.0d0
  T4_probe(2)=half*TANK_MK_HEIGHT+TANK_MK_END_RADIUS-0.025d0
 
  if (axis_dir.eq.2) then ! TPCE aux files
   if (TANK_MK_GEOM_DESCRIPTOR.eq.ZBOT_FLIGHT_ID) then
    T1_probe(1)=xblob2
    T1_probe(2)=zblob2
    T4_probe(1)=xblob2
    T4_probe(2)=zblob2
   endif
  endif

  if (SDIM.eq.2) then
   ! do nothing
  else if (SDIM.eq.3) then
   T4_probe(SDIM)=T4_probe(2)
   T4_probe(2)=T4_probe(1)

   T1_probe(SDIM)=T1_probe(2)
   T1_probe(2)=T1_probe(1)
  else
   print *,"dimension bust"
   stop
  endif

  im=2 ! vapor
  dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
  den=GRID_DATA_IN%den(D_DECL(i,j,k),dencomp)
  temperature=GRID_DATA_IN%den(D_DECL(i,j,k),dencomp+1)
  call init_massfrac_parm(den,massfrac_parm,im)
  do local_ispec=1,num_species_var
   massfrac_parm(local_ispec)= &
       GRID_DATA_IN%den(D_DECL(i,j,k),dencomp+1+local_ispec)
  enddo
  call INTERNAL_CRYOGENIC_TANK_MK(den,massfrac_parm, &
    temperature,internal_energy,TANK_MK_MATERIAL_TYPE,im,num_species_var)
  call EOS_CRYOGENIC_TANK_MK(den,massfrac_parm, &
     internal_energy,pressure,TANK_MK_MATERIAL_TYPE,im,num_species_var)
  support_r=0.0d0
  do dir=1,SDIM
   if (axis_dir.eq.0) then
    support_r=support_r+(GRID_DATA_IN%xsten(0,dir)-T1_probe(dir))**2
   else if (axis_dir.eq.1) then
    support_r=support_r+(GRID_DATA_IN%xsten(0,dir)-T4_probe(dir))**2
   else if (axis_dir.eq.2) then
    support_r=support_r+(GRID_DATA_IN%xsten(0,dir)-T4_probe(dir))**2
   else
    print *,"axis_dir invalid"
    stop
   endif
  enddo
  support_r=sqrt(support_r) 
  dx_coarsest=GRID_DATA_IN%dx(SDIM)
  do ilev=0,level-1
   dx_coarsest=2.0d0*dx_coarsest
  enddo

  if (support_r.le.2.0d0*dx_coarsest) then
   charfn=one
  else if (support_r.gt.2.0d0*dx_coarsest) then
   charfn=0.0d0
  else
   print *,"support_r invalid"
   stop
  endif
  volgrid=GRID_DATA_IN%volgrid
  if (isweep.eq.0) then
   increment_out1(1)=charfn*volgrid
   if (1.eq.0) then
    print *,"nsum1,nsum2 ",nsum1,nsum2
    print *,"charfn,volgrid,pressure,temperature ", &
       charfn,volgrid,pressure,temperature
    print *,"i,j,k ",i,j,k
   endif

  else if (isweep.eq.1) then
   denom=increment_out1(1)
   if (denom.gt.0.0d0) then
    increment_out2(1)=charfn*volgrid*pressure/denom
    increment_out2(2)=charfn*volgrid*temperature/denom
   else
    print *,"expecting denom>0.0:",denom
    print *,"nsum1,nsum2 ",nsum1,nsum2
    print *,"charfn,volgrid,pressure,temperature ", &
       charfn,volgrid,pressure,temperature
    print *,"i,j,k ",i,j,k
    stop
   endif
  else
   print *,"isweep invalid"
   stop
  endif
 else
  print *,"nsum1 or nsum2 invalid"
  stop
 endif

else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_SUMINT

subroutine CRYOGENIC_TANK_MK_INIT_REGIONS_LIST(constant_density_all_time, &
      num_materials_in,num_threads_in)
use probcommon_module

IMPLICIT NONE

integer, INTENT(in) :: num_materials_in
integer, INTENT(in) :: num_threads_in
integer, INTENT(in) :: constant_density_all_time(num_materials_in)
integer :: im,iregion,dir

 if (num_materials_in.eq.num_materials) then
  ! do nothing
 else
  print *,"num_materials_in invalid"
  stop
 endif
 if (num_threads_in.ge.1) then
  ! do nothing
 else
  print *,"num_threads_in invalid: ",num_threads_in
  stop
 endif
 do im=1,num_materials
  if ((constant_density_all_time(im).eq.0).or. &
      (constant_density_all_time(im).eq.1)) then
   ! do nothing
  else
   print *,"constant_density_all_time(im) invalid"
   stop
  endif
 enddo ! im=1..num_materials

 !number_of_source_regions is initialized in: INIT_CRYOGENIC_TANK_MK_MODULE()
 if (number_of_source_regions.ge.0) then
  ! do nothing
 else
  print *,"number_of_source_regions invalid"
  stop
 endif 

 number_of_threads_regions=num_threads_in
 allocate(regions_list(1:number_of_source_regions, &
                       0:number_of_threads_regions))

 do iregion=1,number_of_source_regions
  regions_list(iregion,0)%region_material_id=0
  regions_list(iregion,0)%region_dt=0.0d0  ! timestep
  regions_list(iregion,0)%region_mass_flux=0.0d0
  regions_list(iregion,0)%region_volume_flux=0.0d0
    ! default region_temperature_prescribe=0.0 => homogeneous
    ! flux condition.
  regions_list(iregion,0)%region_temperature_prescribe=0.0d0
  do dir=1,SDIM
   regions_list(iregion,0)%region_velocity_prescribe(dir)=0.0d0
  enddo
  regions_list(iregion,0)%region_energy_flux=0.0d0
  regions_list(iregion,0)%region_volume_raster=0.0d0 
  regions_list(iregion,0)%region_volume=0.0d0 
  regions_list(iregion,0)%region_mass=0.0d0 
  regions_list(iregion,0)%region_energy=0.0d0 
  regions_list(iregion,0)%region_energy_per_kelvin=0.0d0 
  regions_list(iregion,0)%region_volume_after=0.0d0 
  regions_list(iregion,0)%region_mass_after=0.0d0 
  regions_list(iregion,0)%region_energy_after=0.0d0 
 enddo ! iregion=1,number_of_source_regions

 if (axis_dir.eq.0) then 
  regions_list(1,0)%region_material_id=3 !heater
  regions_list(1,0)%region_energy_flux=TANK_MK_HEATER_WATTS ! Watts=J/s
 else if ((axis_dir.eq.1).or. &
          (axis_dir.eq.2)) then 
  regions_list(1,0)%region_material_id=3 ! heater
  regions_list(1,0)%region_energy_flux=TANK_MK_HEATER_WATTS ! Watts=J/s
   ! inflow
  regions_list(2,0)%region_material_id=1
  regions_list(2,0)%region_volume_flux=xblob5
  regions_list(2,0)%region_mass_flux=xblob5*fort_denconst(1)
   ! make xblob6 = 0 if homogeneous flux condition.
  regions_list(2,0)%region_temperature_prescribe=xblob6
  if (TANK_MK_NOZZLE_RAD.gt.0.0d0) then
   regions_list(2,0)%region_velocity_prescribe(SDIM)= &
      xblob5/(Pi*(TANK_MK_NOZZLE_RAD**2.0d0))
  else
   print *,"TANK_MK_NOZZLE_RAD invalid"
   stop
  endif
   ! outflow
   ! default region_temperature_prescribe=0.0 => homogeneous
   ! flux condition.
  regions_list(3,0)%region_material_id=1
  regions_list(3,0)%region_volume_flux=-xblob5
  regions_list(3,0)%region_mass_flux=-xblob5*fort_denconst(1)
 else
  print *,"axis_dir invalid"
  stop
 endif
 
end subroutine CRYOGENIC_TANK_MK_INIT_REGIONS_LIST

! x(SDIM)=vertical direction for gravity
subroutine CRYOGENIC_TANK_MK_CHARFN_REGION(region_id,x,cur_time,charfn_out)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: region_id
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(out) :: charfn_out
real(amrex_real) :: TANK_MK_R_WIDTH
real(amrex_real) :: shell_R,shell_center,LS_SHELL,LS_A,LS_nozzle,zdiff
integer :: called_from_heater_source
real(amrex_real) :: r_cyl
real(amrex_real) :: x3D(3)
real(amrex_real) :: xvel(3)
real(amrex_real) :: xfoot3D(3)
integer auxcomp
real(amrex_real) :: LS_tank(num_materials)

 call convert_to_x3D(x,x3D)
 call rigid_displacement(xfoot3D,cur_time,x3D,xvel)

 if (SDIM.eq.2) then
  r_cyl=abs(xfoot3D(1))
 else if (SDIM.eq.3) then
  r_cyl=sqrt(xfoot3D(1)**2+xfoot3D(3)**2)
 else
  print *,"sdim invalid"
  stop
 endif

if ((num_materials.eq.3).and.(probtype.eq.423)) then

 if (axis_dir.eq.0) then ! ZBOT
  TANK_MK_R_WIDTH=TANK_MK_HEATER_R-TANK_MK_HEATER_R_LOW
  if (TANK_MK_R_WIDTH.gt.0.0d0) then
   if (region_id.eq.1) then
    if ((r_cyl.le.TANK_MK_HEATER_R).and.&
        (r_cyl.ge.TANK_MK_HEATER_R_LOW).and.&
        (xfoot3D(2).ge.TANK_MK_HEATER_LOW).and.&
        (xfoot3D(2).le.TANK_MK_HEATER_HIGH)) then
     charfn_out=one
    else if ((r_cyl.gt.TANK_MK_HEATER_R).or. &
             (r_cyl.lt.TANK_MK_HEATER_R_LOW).or. &
             (xfoot3D(2).lt.TANK_MK_HEATER_LOW).or. &
             (xfoot3D(2).gt.TANK_MK_HEATER_HIGH)) then
     charfn_out=0.0d0
    else
     print *,"position bust"
     stop
    endif
   else
    print *,"region_id invalid"
    stop
   endif
  else
   print *,"TANK_MK_R_WIDTH invalid"
   stop
  endif

 else if (axis_dir.eq.1) then !TPCE

  if (region_id.eq.1) then
   called_from_heater_source=1
   !LS_A>0 in heater
   call CRYOGENIC_TANK_MK_LS_HEATER_A(xfoot3D,LS_A,called_from_heater_source) 
   if (LS_A.ge.0.0d0) then
    charfn_out=one
   else if (LS_A.le.0.0d0) then
    charfn_out=0.0d0
   else
    print *,"LS_A invalid"
    stop
   endif
  else if (region_id.eq.2) then ! inflow
   if ((TANK_MK_NOZZLE_RAD.gt.0.0d0).and. &
       (TANK_MK_NOZZLE_BASE.lt.0.0d0).and. &
       (TANK_MK_NOZZLE_HT.gt.0.0d0).and. &
       (TANK_MK_NOZZLE_THICK_OUTLET.gt.0.0d0)) then
    if ((r_cyl.le.TANK_MK_NOZZLE_RAD).and. &
        (xfoot3D(2).gt.TANK_MK_NOZZLE_BASE+TANK_MK_NOZZLE_HT).and. &
        (xfoot3D(2).le.TANK_MK_NOZZLE_BASE+TANK_MK_NOZZLE_HT+ &
                 TANK_MK_NOZZLE_THICK_OUTLET)) then
     charfn_out=one
    else
     charfn_out=0.0d0
    endif
   else
    print *,"region_id==2 parameter problem"
    stop
   endif
  else if (region_id.eq.3) then ! outflow
   call CRYOGENIC_TANK_MK_LS_NOZZLE(xfoot3D,LS_nozzle)
   if (LS_nozzle.ge.0.0d0) then
    charfn_out=0.0d0
   else if (LS_nozzle.le.0.0d0) then
    if ((TANK_MK_END_CENTER.gt.0.0d0).and. &
        (TANK_MK_HEATER_THICK.gt.0.0d0).and. &
        (TANK_MK_END_RADIUS.gt.0.0d0)) then
     zdiff=xfoot3D(2)+TANK_MK_END_CENTER
     if (zdiff.ge.0.0d0) then
      charfn_out=0.0d0
     else if (zdiff.le.-TANK_MK_END_RADIUS) then
      charfn_out=0.0d0
     else if ((zdiff.lt.0.0d0).and. &
              (zdiff.gt.-TANK_MK_END_RADIUS)) then
      shell_R=sqrt(r_cyl**2+zdiff**2)
      shell_center=TANK_MK_END_RADIUS-0.5d0*TANK_MK_HEATER_THICK
      LS_SHELL=0.5d0*TANK_MK_HEATER_THICK-abs(shell_R-shell_center)
      if (LS_SHELL.ge.0.0d0) then
       charfn_out=one
      else if (LS_SHELL.le.0.0d0) then
       charfn_out=0.0d0
      else
       print *,"LS_SHELL invalid"
       stop
      endif
     else
      print *,"zdiff is NaN"
      stop
     endif
    else
     print *,"parameter bust"
     stop
    endif
   else
    print *,"LS_nozzle is NaN"
    stop
   endif

  else
   print *,"region_id invalid"
   stop
  endif

 else if (axis_dir.eq.2) then

  if (region_id.eq.1) then ! heater A (top)
   auxcomp=1
   call interp_from_aux_grid(auxcomp,xfoot3D,LS_A)
   if (LS_A.ge.0.0d0) then
    charfn_out=one
   else if (LS_A.le.0.0d0) then
    charfn_out=0.0d0
   else
    print *,"LS_A invalid"
    stop
   endif
  else if (region_id.eq.2) then ! inflow
   auxcomp=3
   call interp_from_aux_grid(auxcomp,xfoot3D,LS_A)
   call CRYOGENIC_TANK_MK_LS(x,cur_time,LS_tank,num_materials)
   LS_A=min(LS_A,-LS_tank(3))
   if (LS_A.ge.0.0d0) then
    charfn_out=one
   else if (LS_A.le.0.0d0) then
    charfn_out=0.0d0
   else
    print *,"LS_A invalid"
    stop
   endif
  else if (region_id.eq.3) then ! outflow
   auxcomp=4
   call interp_from_aux_grid(auxcomp,xfoot3D,LS_A)
   call CRYOGENIC_TANK_MK_LS(x,cur_time,LS_tank,num_materials)
   LS_A=min(LS_A,-LS_tank(3))
   if (LS_A.ge.0.0d0) then
    charfn_out=one
   else if (LS_A.le.0.0d0) then
    charfn_out=0.0d0
   else
    print *,"LS_A invalid"
    stop
   endif
  else
   print *,"region_id invalid"
   stop
  endif

 else
  print *,"axis_dir invalid"
  stop
 endif

else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_CHARFN_REGION


! x(SDIM)=vertical direction for gravity
subroutine CRYOGENIC_TANK_MK_THERMAL_K(x,dx,cur_time, &
  density, &
  temperature, &
  thermal_k, &
  im, &
  near_interface, &
  im_solid, &
  temperature_wall, &
  temperature_wall_max, &
  temperature_probe, &
  nrm) ! nrm points from solid to fluid
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: im
integer, INTENT(in) :: im_solid
integer, INTENT(in) :: near_interface
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(in) :: density
real(amrex_real), INTENT(in) :: temperature
real(amrex_real), INTENT(in) :: temperature_wall
real(amrex_real), INTENT(in) :: temperature_wall_max
real(amrex_real), INTENT(in) :: temperature_probe
real(amrex_real), INTENT(in) :: nrm(SDIM) ! nrm points from solid to fluid
real(amrex_real), INTENT(inout) :: thermal_k
real(amrex_real) :: Ra,Gr,Pr,psi,alpha
real(amrex_real) :: mu_w,rho_w,nu,thermal_conductivity,Cp,xi,R,thermal_diffusivity
real(amrex_real) :: gravity_local
integer :: gravity_dir
real(amrex_real) :: expansion_coefficient
integer :: turb_flag

real(amrex_real) :: LS_A
integer :: called_from_heater_source
real(amrex_real) :: r_cyl
real(amrex_real) :: x3D(3)
real(amrex_real) :: xfoot3D(3)
real(amrex_real) :: xvel(3)
integer auxcomp


 call fort_derive_gravity_dir(gravity_vector,gravity_dir)

 call convert_to_x3D(x,x3D)
 call rigid_displacement(xfoot3D,cur_time,x3D,xvel)

 if (SDIM.eq.2) then
  r_cyl=abs(xfoot3D(1))
 else if (SDIM.eq.3) then
  r_cyl=sqrt(xfoot3D(1)**2+xfoot3D(3)**2)
 else
  print *,"sdim invalid"
  stop
 endif

if (probtype.eq.423) then
 ! do nothing
else
 print *,"probtype invalid"
 stop
endif

if ((im_solid.ge.0).and.(im_solid.le.num_materials)) then
 ! do nothing
else
 print *,"im_solid invalid"
 stop
endif

if ((im.ge.1).and.(im.le.num_materials)) then

 if (fort_thermal_microlayer_size(im).gt.zero) then
  !do nothing
 else
  print *,"thermal_microlayer_size(im) invalid"
  stop
 endif  

 if (axis_dir.eq.0) then ! ZBOT

  if (im.eq.2) then ! vapor
   ! do nothing
  else if ((im.eq.1).or.(im.eq.3)) then ! liquid or solid
   if ((r_cyl.le.TANK_MK_HEATER_R).and.&
       (r_cyl.ge.TANK_MK_HEATER_R_LOW-dx(1)).and.&
       (xfoot3D(2).ge.TANK_MK_HEATER_LOW).and.&
       (xfoot3D(2).le.TANK_MK_HEATER_HIGH)) then
    thermal_k=fort_heatviscconst(im)* &
      max(one,dx(1)/fort_thermal_microlayer_size(im))
   else if ((abs(xfoot3D(2)).ge. &
             TANK_MK_INSULATE_THICK+TANK_MK_HEIGHT/2.0d0).or. &
            (r_cyl.ge.TANK_MK_INSULATE_R_HIGH)) then
    if (im.eq.3) then
     thermal_k=0.0d0
    else if (im.eq.1) then
     ! do nothing
    else
     print *,"im invalid"
     stop
    endif
   else if ((abs(xfoot3D(2)).le.TANK_MK_HEIGHT/2.0d0).and. &
            (r_cyl.ge.TANK_MK_INSULATE_R)) then
    if (im.eq.3) then
     thermal_k=0.0d0
    else if (im.eq.1) then
     ! do nothing
    else
     print *,"im invalid"
     stop
    endif
   else if (im.eq.3) then
    ! do nothing
   else if (im.eq.1) then
    mu_w=fort_viscconst(im) 
    rho_w=fort_denconst(im)
    if (rho_w.gt.zero) then
     ! do nothing
    else
     print *,"rho_w invalid"
     stop
    endif
    nu=mu_w/rho_w

    thermal_conductivity=fort_heatviscconst(im)
    Cp=fort_stiffCP(im)
    if (thermal_conductivity.gt.zero) then
     ! do nothing
    else
     print *,"expecting thermal_conductivity to be positive"
     stop
    endif
    if (Cp.gt.zero) then
     ! do nothing
    else
     print *,"Cp invalid"
     stop
    endif
    xi=xfoot3D(2)-TANK_MK_HEATER_LOW
    R=r_cyl
    if ((xi.gt.0.0d0).and. &
        (xi.le.TANK_MK_HEATER_WALL_MODEL).and. &
        (nrm(1).eq.-one).and. &
        (near_interface.eq.1).and. &
        (temperature_wall_max.gt.temperature_probe)) then
     thermal_diffusivity=thermal_conductivity/(rho_w*Cp)
     gravity_local=abs(gravity_vector(gravity_dir))

     call SUB_UNITLESS_EXPANSION_FACTOR(im,temperature_wall_max, &
       temperature_probe,expansion_coefficient)
     expansion_coefficient=abs(expansion_coefficient)

     Gr=gravity_local*expansion_coefficient*(xi**3.0)/(nu*nu)
     Pr=nu/thermal_diffusivity
     Ra=Gr*Pr

     if(Ra.lt.1.0e+9)then
      turb_flag=0
     else if(Ra.ge.1.0e+9)then
      turb_flag=1
     else
      print *,"Ra number invalid"
      stop
     endif
     psi=(1.0d0+(0.492/Pr)**(9.0d0/16.0d0))**(-16.0d0/9.0d0)
     if (turb_flag.eq.0) then
      alpha=(thermal_conductivity/xi)* &
            (0.68d0+0.503d0*((Ra*psi)**0.25d0))
     else if (turb_flag.eq.1) then
      alpha=(thermal_conductivity/xi)* &
            (0.15d0*((Ra*psi)**(1.0d0/3.0d0)))
     else
      print *,"turb_flag invalid"
      stop
     endif 
     thermal_k=max(thermal_k,alpha*dx(1))
    else if ((xi.le.0.0d0).or. &
             (xi.ge.TANK_MK_HEATER_WALL_MODEL).or. &
             (nrm(1).ne.-one).or. &
             (near_interface.eq.0).or. &
             (temperature_wall_max.le.temperature_probe)) then
     ! do nothing
    else
     print *,"xi,nrm,temp_wall, or temp_probe invalid"
     stop
    endif
   else
    print *,"im must be 1 or 3"
    stop
   endif

  else
   print *,"im invalid"
   stop
  endif

 else if ((axis_dir.eq.1).or. &
          (axis_dir.eq.2)) then !TPCE
  if (im.eq.2) then ! vapor
   ! do nothing
  else if ((im.eq.1).or.(im.eq.3)) then ! liquid or solid
   called_from_heater_source=1
   if (axis_dir.eq.1) then
    call CRYOGENIC_TANK_MK_LS_HEATER_A(xfoot3D,LS_A,called_from_heater_source)
   else if (axis_dir.eq.2) then
    auxcomp=1 ! heater A (top)
    call interp_from_aux_grid(auxcomp,xfoot3D,LS_A)
   else
    print *,"axis_dir invalid"
    stop
   endif

   if (LS_A.gt.-dx(1)) then
    thermal_k=fort_heatviscconst(im)* &
      max(one,dx(1)/fort_thermal_microlayer_size(im))
   else if (LS_A.le.-dx(1)) then
    if (im.eq.3) then
     thermal_k=0.0d0
    else if (im.eq.1) then
     ! do nothing
    else
     print *,"im invalid"
     stop
    endif
   else
    print *,"LS_A is NaN"
    stop
   endif
  else
   print *,"im invalid"
   stop
  endif

 else
  print *,"axis_dir out of range"
  stop
 endif

else 
 print *,"im invalid in CRYOGENIC_TANK_MK_THERMAL_K"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_THERMAL_K

! Temperature Stratification in a Cryogenic Fuel Tank
! JOURNAL OF THERMOPHYSICS AND HEAT TRANSFER
! Vol. 27, No. 1, January–March 2013
! Daigle, Smelyanskiy, Boschee, Foygel
!   also
! Upper Stage Tank Thermodynamic Modeling Using SINDA/FLUINT
! 42nd AIAA/ASME/SAE/ASEE Joint Propulsion Conference & Exhibit
! 9 - 12 July 2006, Sacramento, California
! AIAA 2006-5051
! Paul Schallhorn, D. Michael Campbell, Sukhdeep Chase,
! Jorge Piquero, Cindy Fortenberry, Xiaoyi Li, Lisa Grob
!
! This routine only called when "law_of_the_wall==1" associated with 
! im_fluid.
! SDIM=vertical direction for gravity.
subroutine wallfunc_thermocorrelation( &
  dir, & ! =1,2,3
  data_dir, & ! =0,1,2
  dxmin, &
  x_projection_raster, &
  dx, &
  n_raster, & ! points to solid
  u, & !INTENT(in) uimage_raster_solid_frame(dir)
  uimage_tngt_mag, & !INTENT(in) 
  wall_model_velocity, & ! INTENT(in)
  dist_probe, & ! INTENT(in)
  dist_fluid, & ! INTENT(in)
  temperature_image, & !INTENT(in) 
  temperature_wall, & ! INTENT(in)      
  temperature_wall_max, & ! INTENT(in)      
  viscosity_molecular, & ! INTENT(in)      
  viscosity_eddy_wall, & ! INTENT(in)      
  y, & !INTENT(in) distance from image to wall
  ughost_tngt, & ! INTENT(out)
  im_fluid, &  ! INTENT(in)
  critical_length) ! INTENT(in) used for sanity check
use probcommon_module
use global_utility_module
implicit none
integer, INTENT(in) :: dir ! 1,2,3
integer, INTENT(in) :: data_dir ! 0,1,2
real(amrex_real), INTENT(in) :: dxmin
real(amrex_real), INTENT(in), pointer :: x_projection_raster(:)
real(amrex_real), INTENT(in), pointer :: dx(:)
real(amrex_real), INTENT(in), pointer :: n_raster(:) ! points to solid
integer, INTENT(in) :: im_fluid
real(amrex_real), INTENT(in) :: u !uimage_raster_solid_frame(dir)
real(amrex_real), INTENT(in) :: uimage_tngt_mag
real(amrex_real), INTENT(in) :: wall_model_velocity
real(amrex_real), INTENT(in) :: dist_probe
real(amrex_real), INTENT(in) :: dist_fluid
real(amrex_real), INTENT(in) :: temperature_image
real(amrex_real), INTENT(in) :: temperature_wall
real(amrex_real), INTENT(in) :: temperature_wall_max
real(amrex_real), INTENT(in) :: viscosity_molecular
real(amrex_real), INTENT(in) :: viscosity_eddy_wall
real(amrex_real), INTENT(in) :: y !delta_r
real(amrex_real), INTENT(in) :: critical_length
real(amrex_real), INTENT(out) :: ughost_tngt  ! dir direction

real(amrex_real) :: rho_w !wall density
real(amrex_real) :: mu_w  !mu_w: wall molecular viscosity
real(amrex_real) :: thermal_conductivity
real(amrex_real) :: thermal_diffusivity
real(amrex_real) :: gravity_local
integer :: gravity_dir
real(amrex_real) :: expansion_coefficient
real(amrex_real) :: Cp
real(amrex_real) :: Jtemp,Jtemp_no_area,dtemp,vtemp
real(amrex_real),parameter :: local_pi=4.0*atan(1.0d0)
integer :: turb_flag  ! 1 for turb  0 for laminar
real(amrex_real) :: Ra,Gr,Pr
real(amrex_real) :: nu
real(amrex_real) :: xi
real(amrex_real) :: R
real(amrex_real) :: macro_scale_thickness

call fort_derive_gravity_dir(gravity_vector,gravity_dir)

if ((im_fluid.lt.1).or.(im_fluid.gt.num_materials)) then
 print *,"im_fluid invalid in wallfunc_thermocorrelation"
 stop
endif

mu_w=fort_viscconst(im_fluid) 
rho_w=fort_denconst(im_fluid)

if (dx(1).gt.0.0d0) then
 ! do nothing
else
 print *,"dx(1) invalid"
 stop
endif

if (y.gt.zero) then
 ! do nothing
else
 print *,"y should be positive"
 stop
endif
if (dxmin.gt.zero) then
 ! do nothing
else
 print *,"dxmin should be positive"
 stop
endif

if (rho_w.gt.zero) then
 ! do nothing
else
 print *,"rho_w invalid"
 stop
endif

if (viscosity_eddy_wall.eq.fort_viscconst_eddy_wall(im_fluid)) then
 ! do nothing
else
 print *,"viscosity_eddy_wall.eq.fort_viscconst_eddy_wall(im_fluid) == false"
 stop
endif

if (mu_w.eq.viscosity_molecular) then
 ! do nothing
else
 print *,"mu_w.eq.viscosity_molecular == false"
 stop
endif

nu=mu_w/rho_w

thermal_conductivity=fort_heatviscconst(im_fluid)
Cp=fort_stiffCP(im_fluid)

if (thermal_conductivity.gt.zero) then
 ! do nothing
else
 print *,"expecting thermal_conductivity to be positive"
 stop
endif

if (Cp.gt.zero) then
 ! do nothing
else
 print *,"Cp invalid"
 stop
endif

xi=x_projection_raster(SDIM)-TANK_MK_HEATER_LOW
R=x_projection_raster(1)

if (1.eq.0) then
 print *,"xi=",xi
 print *,"n_raster(1)=",n_raster(1)
 print *,"temperature_wall ",temperature_wall
 print *,"temperature_wall_max ",temperature_wall_max
 print *,"temperature_image ",temperature_image
endif

if (dir.eq.data_dir+1) then
 if (abs(n_raster(dir)).eq.one) then
  ! do nothing
 else
  print *,"abs(n_raster(dir)) invalid"
  stop
 endif
else if (dir.ne.data_dir+1) then
 if (n_raster(dir).eq.zero) then
  ! do nothing
 else
  print *,"n_raster(dir) invalid"
  stop
 endif
else
 print *,"dir or data_dir corrupt"
 stop
endif

if ((xi.gt.0.0d0).and. &
    (n_raster(1).eq.one)) then

 thermal_diffusivity=thermal_conductivity/(rho_w*Cp)
 gravity_local=abs(gravity_vector(gravity_dir))
 Pr=nu/thermal_diffusivity

 call SUB_UNITLESS_EXPANSION_FACTOR(im_fluid,temperature_wall_max, &
   temperature_image,expansion_coefficient)
 expansion_coefficient=abs(expansion_coefficient)

 if (temperature_wall_max.gt.temperature_image) then

  Gr=gravity_local*expansion_coefficient*(xi**3.0)/(nu*nu)
  Ra=Gr*Pr

  if(Ra.lt.1.0e+9)then
   turb_flag=0
  else if(Ra.ge.1.0e+9)then
   turb_flag=1
  else
   print *,"Ra number invalid"
   stop
  endif

  vtemp=1.185d0*(nu/xi)*(Gr/(1.0d0+0.494*(Pr**(2.0d0/3.0d0))))**0.5d0

   ! Jtemp is mass flux through horizontal face of length dtemp.
   ! Jtemp has units of mass/sec=rho * velocity * area_face
   ! We assume that J=rho * u_tan * 2 pi R * dr
   ! u_tan=J/(2 pi R rho dr)
  if(turb_flag.eq.1)then
   dtemp=xi*0.565*((1.0d0+0.494d0*Pr**(2.0d0/3.0d0))/Gr)**(0.1d0)/ &
        (Pr**(8.0d0/15.0d0))
   Jtemp_no_area=0.1436d0*rho_w*vtemp
  else if(turb_flag.eq.0)then
   dtemp=xi*3.93*((0.952+Pr)/(Gr*(Pr**2.0d0)))**0.25d0
   Jtemp_no_area=0.0833d0*rho_w*vtemp
  else
   print *,"wrong turb_flag"
   stop
  endif
  Jtemp=2.0d0*local_pi*R*Jtemp_no_area*dtemp

  ! it is known that converged solutions can be obtained on a 128x512 grid.
  ! on a 16x64 grid, choose a thickness associated to the finer grid.
  macro_scale_thickness=dx(1)/16.0d0
  if ((macro_scale_thickness.lt.dtemp).or.(1.eq.1)) then
   macro_scale_thickness=dtemp
  endif
! ughost_tngt=(Jtemp_no_area/rho_w)*dtemp/macro_scale_thickness
  ughost_tngt=vtemp ! vtemp is the maximum tangential velocity in the
                    ! boundary layer region.

 else if (temperature_wall_max.le.temperature_image) then
  Gr=zero
  Ra=zero
  turb_flag=0
  vtemp=zero
  dtemp=zero
  Jtemp_no_area=zero
  Jtemp=zero
  macro_scale_thickness=zero
  ughost_tngt=zero
 else
  print *,"temperature_wall_max or temperature_image is NaN"
  stop
 endif

 if (wall_model_velocity.eq.zero) then
  ! do nothing
 else if (wall_model_velocity.ne.zero) then
  ughost_tngt=wall_model_velocity
 else
  print *,"wall_model_velocity is NaN"
  stop
 endif

  ! do not prescribe a wall velocity if near (or in) another fluid.
  !
 if ((dist_probe.lt.dx(SDIM)).or. &
     (dist_fluid.lt.dx(SDIM))) then
  ughost_tngt=0.0d0
 else if ((dist_probe.ge.dx(SDIM)).and. &
          (dist_fluid.ge.dx(SDIM))) then
  ! do nothing
 else
  print *,"dist_probe or dist_fluid is NaN"
  stop
 endif

 if (1.eq.0) then
  print *,"xi=",xi
  print *,"Gr,Pr,Ra,vtemp,dtemp ",Gr,Pr,Ra,vtemp,dtemp
  print *,"Jtemp=",Jtemp
  print *,"R=",R
  print *,"rho_w=",rho_w
  print *,"dx(1)=",dx(1)
  print *,"ughost_tngt=",ughost_tngt
 endif

else if ((xi.le.0.0d0).or. &
         (n_raster(1).ne.one)) then
 ughost_tngt=0.0d0
else
 print *,"xi,n_raster,twall or timage became corrupt"
 stop
endif

end subroutine wallfunc_thermocorrelation


! SDIM=vertical direction for gravity.
subroutine CRYOGENIC_TANK_MK_wallfunc( &
  dir, & ! =1,2,3
  data_dir, & ! =0,1,2
  dxmin, &
  x_projection_raster, &
  dx, &
  n_raster, & ! points to solid
  u, & !INTENT(in) uimage_raster_solid_frame(dir)
  uimage_tngt_mag, & !INTENT(in) 
  wall_model_velocity, & ! INTENT(in)
  dist_probe, & ! INTENT(in)
  dist_fluid, & ! INTENT(in)
  temperature_image, & !INTENT(in) 
  temperature_wall, & ! INTENT(in)      
  temperature_wall_max, & ! INTENT(in)      
  viscosity_molecular, & ! INTENT(in)      
  viscosity_eddy_wall, & ! INTENT(in)      
  y, & !INTENT(in) distance from image to wall
  ughost_tngt, & ! INTENT(out)
  im_fluid, &  ! INTENT(in)
  critical_length) ! INTENT(in) used for sanity check
use probcommon_module
use global_utility_module
implicit none
integer, INTENT(in) :: dir ! 1,2,3
integer, INTENT(in) :: data_dir ! 0,1,2
real(amrex_real), INTENT(in) :: dxmin
real(amrex_real), INTENT(in), pointer :: x_projection_raster(:)
real(amrex_real), INTENT(in), pointer :: dx(:)
real(amrex_real), INTENT(in), pointer :: n_raster(:) ! points to solid
integer, INTENT(in) :: im_fluid
real(amrex_real), INTENT(in) :: u !uimage_raster_solid_frame(dir)
real(amrex_real), INTENT(in) :: uimage_tngt_mag
real(amrex_real), INTENT(in) :: wall_model_velocity
real(amrex_real), INTENT(in) :: dist_probe
real(amrex_real), INTENT(in) :: dist_fluid
real(amrex_real), INTENT(in) :: temperature_image
real(amrex_real), INTENT(in) :: temperature_wall
real(amrex_real), INTENT(in) :: temperature_wall_max
real(amrex_real), INTENT(in) :: viscosity_molecular
real(amrex_real), INTENT(in) :: viscosity_eddy_wall
real(amrex_real), INTENT(in) :: y !delta_r
real(amrex_real), INTENT(in) :: critical_length
real(amrex_real), INTENT(out) :: ughost_tngt  ! dir direction

 if (1.eq.0) then
  ! remark: "subroutine wallfunc_newtonsmethod" is 
  ! declared in GLOBALUTIL.F90
  call wallfunc_newtonsmethod( &
   dir, & ! =1,2,3
   data_dir, & ! =0,1,2
   dxmin, &
   x_projection_raster, &
   dx, &
   n_raster, & ! points to solid
   u, & !INTENT(in) uimage_raster_solid_frame(dir)
   uimage_tngt_mag, & !INTENT(in) 
   wall_model_velocity, & ! INTENT(in)
   dist_probe, & ! INTENT(in)
   dist_fluid, & ! INTENT(in)
   temperature_image, & !INTENT(in) 
   temperature_wall, & ! INTENT(in)      
   temperature_wall_max, & ! INTENT(in)      
   viscosity_molecular, & ! INTENT(in)      
   viscosity_eddy_wall, & ! INTENT(in)      
   y, & !INTENT(in) distance from image to wall
   ughost_tngt, & ! INTENT(out)
   im_fluid, &  ! INTENT(in)
   critical_length) ! INTENT(in) used for sanity check
 else
  call wallfunc_thermocorrelation( &
   dir, & ! =1,2,3
   data_dir, & ! =0,1,2
   dxmin, &
   x_projection_raster, &
   dx, &
   n_raster, & ! points to solid
   u, & !INTENT(in) uimage_raster_solid_frame(dir)
   uimage_tngt_mag, & !INTENT(in) 
   wall_model_velocity, & ! INTENT(in)
   dist_probe, & ! INTENT(in)
   dist_fluid, & ! INTENT(in)
   temperature_image, & !INTENT(in) 
   temperature_wall, & ! INTENT(in)      
   temperature_wall_max, & ! INTENT(in)      
   viscosity_molecular, & ! INTENT(in)      
   viscosity_eddy_wall, & ! INTENT(in)      
   y, & !INTENT(in) distance from image to wall
   ughost_tngt, & ! INTENT(out)
   im_fluid, &  ! INTENT(in)
   critical_length) ! INTENT(in) used for sanity check
 endif

end subroutine CRYOGENIC_TANK_MK_wallfunc

subroutine CRYOGENIC_TANK_MK_K_EFFECTIVE( &
  interface_mass_transfer_model, &
  ireverse, &
  iten, &        
  molar_mass, & ! index: 1..nmat
  species_molar_mass, & ! index: 1..num_species_var
  k_model_predict, &
  k_model_correct, &
  k_physical_base, &
  T_probe_src, &
  T_probe_dst, &
  dxprobe_src, &
  dxprobe_dst, &
  LL, &
  num_materials_in, &
  num_species_var_in)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: interface_mass_transfer_model
integer, INTENT(in) :: num_materials_in
integer, INTENT(in) :: num_species_var_in
integer, INTENT(in) :: ireverse
integer, INTENT(in) :: iten
real(amrex_real), INTENT(in) :: molar_mass(num_materials_in)
real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var_in)
real(amrex_real), INTENT(in) :: k_model_predict(2) ! src,dst
real(amrex_real), INTENT(inout) :: k_model_correct(2) ! src,dst
real(amrex_real), INTENT(in) :: k_physical_base(2) ! src, dst
real(amrex_real), INTENT(in) :: T_probe_src
real(amrex_real), INTENT(in) :: T_probe_dst
real(amrex_real), INTENT(in) :: LL
real(amrex_real), INTENT(in) :: dxprobe_src
real(amrex_real), INTENT(in) :: dxprobe_dst

real(amrex_real) :: RA1,RA2

if (interface_mass_transfer_model.eq.1) then

  !December 14: replaced 302 with 307.
 RA1=(1.7069d+9)*abs(307.0d0-T_probe_src)
 RA2=(5.8395d+8)*abs(307.0d0-T_probe_dst)

 if (RA1.gt.1.0d+4.and.RA1.le.1.0d+7)then
  k_model_correct(1)=dxprobe_src*k_model_predict(1)/0.1016d0* &
           0.54d0*(RA1**(1.0d0/4.0d0))
 else if (RA1.ge.1.0d+7.and.RA1.le.1.0d+11)then
  k_model_correct(1)=dxprobe_src*k_model_predict(1)/0.1016d0* &
           0.15d0*(RA1**(1.0d0/3.0d0))
 else
   print *,"invalid Ra1 number"
   stop
 endif

 if (RA2.gt.1.0d+4.and.RA2.le.1.0d+7)then
  k_model_correct(2)=dxprobe_dst*k_model_predict(2)/0.1016d0* &
           0.54d0*(RA2**(1.0d0/4.0d0))
 else if (RA2.ge.1.0d+7.and.RA2.le.1.0d+11)then
  k_model_correct(2)=dxprobe_dst*k_model_predict(2)/0.1016d0* &
           0.15d0*(RA2**(1.0d0/3.0d0))
 else
   print *,"invalid Ra2 number"
   stop
 endif

else if (interface_mass_transfer_model.eq.0) then
 ! do nothing
else if (interface_mass_transfer_model.gt.0) then
 ! do nothing
else
 print *,"interface_mass_transfer_model invalid"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_K_EFFECTIVE


! returns (1/w) where w>>1 in "trouble" regions
! SDIM is vertical direction of gravity
subroutine CRYOGENIC_TANK_MK_MAPPING_WEIGHT_COEFF(dir,wt,phys_x)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: dir
real(amrex_real), INTENT(out) :: wt
real(amrex_real), INTENT(in) :: phys_x

if ((dir.ge.0).and.(dir.lt.SDIM)) then
 ! do nothing
else
 print *,"dir invalid"
 stop
endif
if ((phys_x.ge.zero).or.(phys_x.le.zero)) then
 ! do nothing
else
 print *,"phys_x is NaN"
 stop
endif

wt=one

if (fort_grid_stretching_parameter(1).gt.zero) then

 use_identity_mapping=0

 if (dir.eq.0) then
  if (phys_x.le.problenx/10.0d0) then
   wt=1.0d0/100.0d0
  endif
  if (abs(phys_x-xblob).le.problenx/10.0d0) then
   wt=1.0d0/100.0d0
  endif
 endif
endif

return
end subroutine CRYOGENIC_TANK_MK_MAPPING_WEIGHT_COEFF

subroutine CRYOGENIC_TANK_MK_angular_velocity(x,cur_time, &
   angular_velocity,angular_velocity_custom, &
   angular_velocity_dot,lever_arm)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(in) :: angular_velocity
real(amrex_real), INTENT(out) :: angular_velocity_custom
real(amrex_real), INTENT(out) :: angular_velocity_dot
real(amrex_real), INTENT(out) :: lever_arm

 if (cur_time.ge.0.0d0) then
  ! do nothing
 else
  print *,"cur_time invalid: ",cur_time
  stop
 endif

 if (angular_velocity.ge.zero) then
  !do nothing
 else
  print *,"angular_velocity invalid"
  stop
 endif

 if (xblob8.ge.zero) then
  !do nothing
 else
  print *,"expecting xblob8>=0"
  stop
 endif

 angular_velocity_custom=angular_velocity
 angular_velocity_dot=zero
 lever_arm=radblob8
 if (cur_time.ge.xblob8) then
  !do nothing
 else if ((cur_time.ge.zero).and.(cur_time.le.xblob8)) then
  angular_velocity_custom=angular_velocity*cur_time/xblob8
  angular_velocity_dot=angular_velocity/xblob8
 else
  print *,"cur_time invalid"
  stop
 endif

end subroutine CRYOGENIC_TANK_MK_angular_velocity

subroutine CRYOGENIC_TANK_MK_gravity_vector(x,cur_time, &
   gravity_vector_in, &
   gravity_vector_out)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(in) :: gravity_vector_in(SDIM)
real(amrex_real), INTENT(out) :: gravity_vector_out(SDIM)
integer :: dir

 if (cur_time.ge.0.0d0) then
  ! do nothing
 else
  print *,"cur_time invalid CRYOGENIC_TANK_MK_gravity_vector: ",cur_time
  stop
 endif

 do dir=1,SDIM
  gravity_vector_out(dir)=gravity_vector_in(dir)
 enddo

end subroutine CRYOGENIC_TANK_MK_gravity_vector


end module CRYOGENIC_TANK_MK_module
