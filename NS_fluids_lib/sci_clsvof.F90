#undef BL_LANG_CC
#define BL_LANG_FORT

#define STANDALONE 0

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#define element_buffer_tol 0.01

! 10 seconds for tail to do a full period
#define WHALE_LENGTH 13.0
#define PERIOD_TAIL 10.0
#define DT_DUFFY 0.01
#define STEPS_DUFFY 800
#define injG 1
#define MAX_PARTS 100
#define flags_per_element 3
#define CTMLoverflow (1.0D+20)

#ifdef BL_USE_MPI
#define mpi_activate 1
#else
#ifdef BL_USE_MPI3
#define mpi_activate 1
#else
#define mpi_activate 0
#endif
#endif

module CLSVOFCouplerIO
use probcommon_module

implicit none

type lag_type
 INTEGER_T :: n_nodes,n_elems
 REAL_T, pointer :: nd(:,:)    ! nd(dir,node_id) dir=1..3
 REAL_T, pointer :: ndvel(:,:) ! ndvel(dir,node_id) dir=1..3
 REAL_T, pointer :: ndmass(:) ! ndmass(node_id) 
 REAL_T, pointer :: ndforce(:,:) ! ndforce(dir,node_id)  dir=1..6
 REAL_T, pointer :: ndtemp(:)  ! ndtemp(node_id)
  ! number of nodes in element=elemdt(1,elemid)
  ! part number=elemdt(2,elemid)
  ! doubly wetted=elemdt(3,elemid)
 INTEGER_T, pointer :: elemdt(:,:) 
  ! node_id=intelemdt(1..3,elemid)
  ! original_elem_id=intelemdt(4,elemid)
 INTEGER_T, pointer :: intelemdt(:,:) 
end type lag_type

type mesh_type
 INTEGER_T :: PartID
 INTEGER_T :: flag_2D_to_3D
 INTEGER_T :: refine_factor
 INTEGER_T :: bounding_box_ngrow
 REAL_T :: max_side_len
 REAL_T :: min_side_len
 REAL_T :: max_side_len_refined
 REAL_T :: min_side_len_refined
 INTEGER_T :: IntElemDim,IntElemDimPaddle,IntElemDimPool
 INTEGER_T :: NumNodes,NumNodesPaddle,NumNodesPool
 INTEGER_T :: NumIntElems,NumIntElemsPaddle,NumIntElemsPool
 INTEGER_T, pointer :: ElemData(:,:)
 INTEGER_T, pointer :: IntElem(:,:)
 INTEGER_T, pointer :: Eul2IntNode(:)
 INTEGER_T :: NumNodesBIG
 INTEGER_T :: NumIntElemsBIG
 INTEGER_T, pointer :: ElemNodeCountBIG(:)
 REAL_T, pointer :: NodeBIG(:,:)
 REAL_T, pointer :: NodeVelBIG(:,:)
 REAL_T, pointer :: NodeForceBIG(:,:)  ! (6,node_id)
 REAL_T, pointer :: NodeMassBIG(:)
 REAL_T, pointer :: NodeTempBIG(:)  
 REAL_T, pointer :: NodeNormalBIG(:,:)
 REAL_T, pointer :: ElemDataXnotBIG(:,:)
 INTEGER_T, pointer :: ElemDataBIG(:,:)
 INTEGER_T, pointer :: IntElemBIG(:,:) ! IntElemBIG(inode,ielem)
 REAL_T, pointer :: Node(:,:)  ! Node(dir,inode)
 REAL_T, pointer :: Node_old(:,:)
 REAL_T, pointer :: Node_new(:,:)
 REAL_T, pointer :: Node_current(:,:)
 REAL_T, pointer :: NodeVel(:,:)
 REAL_T, pointer :: NodeVel_old(:,:)
 REAL_T, pointer :: NodeVel_new(:,:)
 REAL_T, pointer :: NodeForce(:,:)
 REAL_T, pointer :: NodeForce_old(:,:)
 REAL_T, pointer :: NodeForce_new(:,:)
 REAL_T, pointer :: NodeMass(:)
 REAL_T, pointer :: NodeTemp(:)
 REAL_T, pointer :: NodeTemp_old(:)
 REAL_T, pointer :: NodeTemp_new(:)
 REAL_T soliddrop_displacement
 REAL_T soliddrop_speed
 REAL_T solid_displ(3)
 REAL_T solid_speed(3)
 INTEGER_T deforming_part
end type mesh_type


INTEGER_T :: use_temp
INTEGER_T :: istepB,sci_sdim,sci_istop,sci_istep
REAL_T :: sci_curtime,sci_dt
REAL_T :: timeB,tstart,tfinish
REAL_T :: dtB

type(lag_type), dimension(:), allocatable :: multi_lag

type(mesh_type), dimension(MAX_PARTS) :: FSI

INTEGER_T :: normal_invert

REAL_T :: radradblob,denpaddle,dampingpaddle,TorquePos,TorqueVel,radradblobwall
REAL_T :: raddust,dendust,adheredust,floordust,tempdust
REAL_T, dimension(3) :: torquePosDust,torqueVelDust, &
 centerDust,centerVelDust,centerStartDust
REAL_T, dimension(3) :: xxblob,newxxblob,xxblobwall,newxxblobwall
character(35) :: whalein,whaleout
INTEGER_T :: whale_nodes
REAL_T, dimension(:), allocatable :: RR,SS,TT,UU,VV,WW,XX,YY,ZZ 
REAL_T, dimension(:), allocatable :: whale_spring
REAL_T, dimension(:), allocatable :: whale_X_init,whale_Y_init,whale_Z_init
! Nodes,20
INTEGER_T, DIMENSION(:,:), allocatable :: whale_list
REAL_T, Dimension(22) :: whale_angle, whale_timestep
INTEGER_T whale_cells,whale_counter
REAL_T whale_counter_real
REAL_T CLSVOF_whale_time
INTEGER_T exclusive_doubly_wetted

REAL_T problo_ref(AMREX_SPACEDIM)
REAL_T probhi_ref(AMREX_SPACEDIM)
REAL_T problen_ref(AMREX_SPACEDIM)

REAL_T problo_act(3)
REAL_T probhi_act(3)
REAL_T problen_act(3)

INTEGER_T FSI_NPARTS
INTEGER_T CTML_NPARTS
INTEGER_T TOTAL_NPARTS
INTEGER_T im_solid_mapF(MAX_PARTS)
INTEGER_T CTML_partid_map(MAX_PARTS)
INTEGER_T FSI_partid_map(MAX_PARTS)

INTEGER_T ctml_n_fib_bodies
INTEGER_T ctml_max_n_fib_nodes
INTEGER_T, dimension(:), allocatable :: ctml_n_fib_nodes
INTEGER_T, dimension(:), allocatable :: ctml_n_fib_active_nodes
REAL_T, dimension(:,:,:), allocatable :: ctml_fib_pst
REAL_T, dimension(:,:,:), allocatable :: ctml_fib_vel
REAL_T, dimension(:,:,:), allocatable :: ctml_fib_frc
REAL_T, dimension(:,:), allocatable :: ctml_fib_mass

contains


      subroutine checkbound3D(lo,hi, &
      DIMS3D(data), &
      ngrow,dir,id)
      IMPLICIT NONE

      INTEGER_T lo(3), hi(3)
      INTEGER_T DIMDEC3D(data)
      INTEGER_T ngrow,dir,id

      INTEGER_T ii(3)

      INTEGER_T hidata(3)
      INTEGER_T lodata(3)
      INTEGER_T dir2

      hidata(1)=ARG3D_H1(data)
      hidata(2)=ARG3D_H2(data)
      hidata(3)=ARG3D_H3(data)
      lodata(1)=ARG3D_L1(data)
      lodata(2)=ARG3D_L2(data)
      lodata(3)=ARG3D_L3(data)

      do dir2=1,3
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound 3d id=",id
        print *,"dir2=",dir2
        stop
       endif
       ii(dir2)=0
      enddo
      if ((dir.ge.0).and.(dir.lt.3)) then
       ii(dir+1)=1
      else if (dir.eq.-1) then
       ! do nothing
      else
       print *,"dir invalid checkbound3d"
       stop
      endif
 
      if (ngrow.lt.0) then
       print *,"ngrow invalid in checkbound 3d"
       stop
      endif
      if (id.lt.0) then
       print *,"id invalid in checkbound 3d"
       stop
      endif

      do dir2=1,3
       if (lo(dir2).lt.0) then
        print *,"lo invalid in checkbound 3d id=",id
        print *,"dir2=",dir2
        stop
       endif
       if (lodata(dir2).gt.lo(dir2)-ngrow) then
        print *,"lo mismatch id=",id
        print *,"dir2=",dir2
        stop
       endif
       if (hidata(dir2).lt.hi(dir2)+ngrow+ii(dir2)) then
        print *,"hi mismatch id=",id
        print *,"dir2=",dir2
        stop
       endif
      enddo ! dir2

      return
      end subroutine checkbound3D


subroutine init2_FSI(part_id)
IMPLICIT NONE

INTEGER_T part_id
INTEGER_T inode,dir

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif

 do inode=1,FSI(part_id)%NumNodes
  do dir=1,3
   FSI(part_id)%NodeVel_old(dir,inode)=0.0
  enddo
  do dir=1,6
   FSI(part_id)%NodeForce_old(dir,inode)=0.0
  enddo
  do dir=1,3
   FSI(part_id)%Node_current(dir,inode)=FSI(part_id)%Node_new(dir,inode)
   FSI(part_id)%NodeVel_new(dir,inode)=FSI(part_id)%NodeVel_old(dir,inode)
  enddo
  do dir=1,6
   FSI(part_id)%NodeForce_new(dir,inode)=FSI(part_id)%NodeForce_old(dir,inode)
  enddo
   ! in: init2_FSI
  FSI(part_id)%NodeMass(inode)=one
 enddo  ! inode=1,NumNodes

return
end subroutine init2_FSI


subroutine init3_FSI(part_id,ifirst,do_2nd_part,ioproc,isout)
IMPLICIT NONE

INTEGER_T part_id
INTEGER_T inode,dir,ifirst,do_2nd_part,it,ioproc,isout
REAL_T x,y,z,z0,z90,t,dt,t1,t2,inflowvel
REAL_T YK,ZK,lift0,lift90
REAL_T, dimension(3) :: displ1,displ2
REAL_T dilated_time

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id invalid"
   stop
  endif
  if (FSI(part_id)%partID.ne.part_id) then
   print *,"FSI(part_id)%partID.ne.part_id"
   stop
  endif

   ! in: init3_FSI
  if (ifirst.eq.1) then
   allocate(FSI(part_id)%Node(3,FSI(part_id)%NumNodes))
   allocate(FSI(part_id)%NodeVel(3,FSI(part_id)%NumNodes))
   allocate(FSI(part_id)%NodeForce(6,FSI(part_id)%NumNodes))
   allocate(FSI(part_id)%NodeTemp(FSI(part_id)%NumNodes))
  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  do inode=1,FSI(part_id)%NumNodes
   do dir=1,3
    FSI(part_id)%Node(dir,inode)=FSI(part_id)%Node_current(dir,inode)
    FSI(part_id)%NodeVel(dir,inode)=FSI(part_id)%NodeVel_new(dir,inode)
   enddo
   do dir=1,6
    FSI(part_id)%NodeForce(dir,inode)=FSI(part_id)%NodeForce_new(dir,inode)
   enddo
    ! in: init3_FSI
   FSI(part_id)%NodeTemp(inode)=FSI(part_id)%NodeTemp_new(inode)
  enddo ! inode=1,FSI(part_id)%NumNodes

  do dir=1,3
   FSI(part_id)%solid_displ(dir)=zero
   FSI(part_id)%solid_speed(dir)=zero
  enddo

  FSI(part_id)%soliddrop_displacement=zero
  FSI(part_id)%soliddrop_speed=zero

  if (probtype.eq.531) then

   if (timeB.le.0.03) then  ! actual time x 10^3 ?
    FSI(part_id)%soliddrop_displacement=-0.115*timeB
    FSI(part_id)%soliddrop_speed=-0.115
   else
    FSI(part_id)%soliddrop_displacement=-0.115*0.03 ! actual time x 10^3 ?
    FSI(part_id)%soliddrop_speed=0.0
   endif

! MARK:
! if probtype=538 and part_id=1 and RZ, then 
! dist(x,y,z)=dist(x+0.0015,y-0.0045,z)
! if probtype=538 and part_id=2, then
! dist(x,y,z)= dist(x,y,z+0.01)
! IN FUTURE, DO NOT CALL "UNITE", INSTEAD
! USE BOTH FSI and FSI_NEEDLE when returning LS or VEL.
! ALSO IN FUTURE, DO NOT REPEATEDLY CALL GENERATE_NEW_TRIANGLES
! IN CLSVOF_ReadNodes.
!  
  else if ((probtype.eq.538).or.(probtype.eq.541)) then

   if (FSI(part_id)%PartID.eq.2) then

    if (probtype.eq.538) then
     dilated_time=timeB
     dilated_time = dilated_time*51.
    else if (probtype.eq.541) then
     if (AMREX_SPACEDIM.eq.3) then
      dilated_time=timeB-1.460e-3 ! skip first 200 mus
     else
      dilated_time=timeB ! skip first 200 mus
     endif
    endif

    ! load trajectory of needle tip

    OPEN(unit=15,file="needle_motion.dat",access='sequential', &
     form="formatted",status='old')
    dt = 0.
    do it=1,251

     if (probtype.eq.538) then
      READ(15,*) t,z0,z90,x,y
     else if (probtype.eq.541) then
      READ(15,*) t,lift0,lift90,YK,ZK
     endif

     if (probtype.eq.538) then
      t = t*1e-6 ! mus to s
      x = (x+24.32)*1e-4 ! mum to cm (IJFM 2013)
      y = -(y+33.56)*1e-4
      z0 = (z0-468.0)*1e-4 ! at t = 910 mus
      z90 = (z90-470.75)*1e-4 
     else if (probtype.eq.541) then
!     t = t*1e-6 ! mus to s
!     x = (YK+24.32)*1e-4 ! mum to cm ! settings for IJFM 2013
!     y = -(ZK+33.56)*1e-4
!     z0 = (lift0-468.0)*1e-4 ! at t = 910 mus
!     z90 = (lift90-470.75)*1e-4
      t = t*1e-6 ! mus to s
        !rotate by 9 deg and scale from mum to cm
      y = -(0.987688*YK+0.156434*ZK)*1e-4 
      x =  (0.987688*ZK-0.156434*YK)*1e-4
      z0 = -lift0*1e-4
      z90= -lift90*1e-4
     endif

     if (t.le.dilated_time) then
      t1 = t
      displ1(1) = x
      displ1(2) = y
      displ1(3) = 0.5*(z0+z90)
     else if (t.gt.dilated_time) then
      t2 = t
      displ2(1) = x
      displ2(2) = y
      displ2(3) = 0.5*(z0+z90)
      dt = t2-t1
      goto 10
     endif
    enddo ! it
    ! interpolate and differentiate
10  continue
    CLOSE(15)
    if (dt.le.0.) then
     if ((ioproc.eq.1).and.(isout.eq.1)) then
      print*," bad needle trajectory data dilated time ",dilated_time
     endif
     dt=1.0
    endif
    do dir=1,3
     if (injG.eq.0) then
      FSI(part_id)%solid_displ(dir)= &
       displ1(dir)+(dilated_time-t1)/dt*(displ2(dir)-displ1(dir))
      FSI(part_id)%solid_speed(dir)= &
       (displ2(dir)-displ1(dir))/dt
     else if (injG.eq.1) then
      FSI(part_id)%solid_displ(dir)=zero
      FSI(part_id)%solid_speed(dir)=zero
     else
      print *,"injG invalid"
      stop
     endif


     if ((AMREX_SPACEDIM.eq.2).and.(probtype.eq.541)) then
      FSI(part_id)%solid_displ(dir)=0.0
      FSI(part_id)%solid_speed(dir)=0.0
     endif
    enddo ! dir

    if (injG.eq.0) then
     ! do nothing
    else if (injG.eq.1) then
     FSI(part_id)%solid_displ(3)=-0.005  ! fully open sprayG
    else
     print *,"injG invalid"
     stop
    endif

    if (levelrz.eq.0) then
     ! do nothing
    else if (levelrz.eq.1) then
     FSI(part_id)%solid_speed(1)=zero
     FSI(part_id)%solid_speed(2)=zero
    else
     print *,"levelrz invalid init3 fsi"
     stop
    endif 

    if ((ioproc.eq.1).and.(isout.eq.1)) then
     print*,"sci_clsvof: interpolation at main time ",timeB
     print*,"sci_clsvof: interpolation at dilated time ",dilated_time, &
       " between ",t1," and ",t2
     print*,"solid displacement: ",FSI(part_id)%solid_displ
     print*,"solid velocity: ",FSI(part_id)%solid_speed
    endif
   else if (FSI(part_id)%PartID.eq.1) then

    do dir=1,3
     FSI(part_id)%solid_displ(dir)=0.
     FSI(part_id)%solid_speed(dir)=0.
    enddo

   else
    print *,"part id invalid"
    stop
   endif
    
  endif  ! probtype.eq.538 or 541

  if (do_2nd_part.eq.1) then

    ! solid body translation is default

   do inode=1,FSI(part_id)%NumNodes
    x=FSI(part_id)%Node(1,inode)
    y=FSI(part_id)%Node(2,inode)
    z=FSI(part_id)%Node(3,inode)

     ! viorel's sphere problem, cannot specify velocity at outer shell
     ! nodes because they are too far apart.
    if (probtype.eq.5601) then  ! viorel's sphere problem
     inflowvel=0.0
     FSI(part_id)%NodeVel_old(1,inode)=inflowvel
     FSI(part_id)%NodeVel_new(1,inode)=inflowvel
     FSI(part_id)%NodeVel(1,inode)=inflowvel
    else if (probtype.eq.5602) then  ! internal inflow problem
     inflowvel=0.0
     if ((z.le.-19.0).and.(z.ge.-21.0)) then
      if ((x.ge.14.0).and.(x.le.27.0).and.(y.ge.9.0).and.(y.le.16.0)) then
       inflowvel=1.0
      endif
      if ((x.ge.14.0).and.(x.le.21.0).and.(y.ge.49.0).and.(y.le.56.0)) then
       inflowvel=-1.0
      endif 
     endif
     FSI(part_id)%NodeVel_old(3,inode)=inflowvel
     FSI(part_id)%NodeVel_new(3,inode)=inflowvel
     FSI(part_id)%NodeVel(3,inode)=inflowvel
    else if (1.eq.0) then
     FSI(part_id)%NodeVel_old(3,inode)=FSI(part_id)%soliddrop_speed
     FSI(part_id)%NodeVel_new(3,inode)=FSI(part_id)%soliddrop_speed
     FSI(part_id)%NodeVel(3,inode)=FSI(part_id)%soliddrop_speed
     FSI(part_id)%Node(3,inode)=FSI(part_id)%Node_old(3,inode)+ &
       FSI(part_id)%soliddrop_displacement
     FSI(part_id)%Node_new(3,inode)=FSI(part_id)%Node(3,inode)
    else
     do dir=1,3
      FSI(part_id)%NodeVel_old(dir,inode)=FSI(part_id)%solid_speed(dir)
      FSI(part_id)%NodeVel_new(dir,inode)=FSI(part_id)%solid_speed(dir)
      FSI(part_id)%NodeVel(dir,inode)=FSI(part_id)%solid_speed(dir)

! MARK: displacement taken into account elsewhere.
! do not want to change the node positions since it will
! effect what happens in generate_new_triangles.

      if (probtype.eq.538) then
       FSI(part_id)%Node(dir,inode)=FSI(part_id)%Node_old(dir,inode)
      else
       FSI(part_id)%Node(dir,inode)=FSI(part_id)%Node_old(dir,inode)+ &
         FSI(part_id)%solid_displ(dir)
      endif
      FSI(part_id)%Node_new(dir,inode)=FSI(part_id)%Node(dir,inode)
     enddo ! dir

    endif

   enddo  ! inode
  else if (do_2nd_part.eq.0) then
   ! do nothing
  else
   print *,"do_2nd_part invalid"
   stop
  endif

return
end subroutine init3_FSI

subroutine init_FSI(part_id,allocate_intelem)
IMPLICIT NONE

INTEGER_T part_id
INTEGER_T inode,dir,allocate_intelem

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif

 allocate(FSI(part_id)%Eul2IntNode(FSI(part_id)%NumNodes))
 !(1,iface)=nodes per element (2,iface)=part num (3,iface)=doubly wet flag
 allocate(FSI(part_id)%ElemData(flags_per_element,FSI(part_id)%NumIntElems))
 if (allocate_intelem.eq.1) then
  allocate(FSI(part_id)%IntElem(FSI(part_id)%IntElemDim, &
           FSI(part_id)%NumIntElems))
 endif
 allocate(FSI(part_id)%Node_old(3,FSI(part_id)%NumNodes))
 allocate(FSI(part_id)%Node_new(3,FSI(part_id)%NumNodes))
 allocate(FSI(part_id)%Node_current(3,FSI(part_id)%NumNodes))
 allocate(FSI(part_id)%NodeVel_old(3,FSI(part_id)%NumNodes))
 allocate(FSI(part_id)%NodeVel_new(3,FSI(part_id)%NumNodes))
 allocate(FSI(part_id)%NodeForce_old(6,FSI(part_id)%NumNodes))
 allocate(FSI(part_id)%NodeForce_new(6,FSI(part_id)%NumNodes))
 allocate(FSI(part_id)%NodeTemp_old(FSI(part_id)%NumNodes))
 allocate(FSI(part_id)%NodeTemp_new(FSI(part_id)%NumNodes))

  ! in: init_FSI
 allocate(FSI(part_id)%NodeMass(FSI(part_id)%NumNodes))

 do inode=1,FSI(part_id)%NumNodes
  FSI(part_id)%Eul2IntNode(inode)=inode

   ! in: init_FSI
  FSI(part_id)%NodeMass(inode)=one

  FSI(part_id)%NodeTemp_old(inode)=0.0
  FSI(part_id)%NodeTemp_new(inode)=0.0
  do dir=1,3
   FSI(part_id)%Node_old(dir,inode)=0.0
   FSI(part_id)%Node_new(dir,inode)=0.0
   FSI(part_id)%Node_current(dir,inode)=0.0
   FSI(part_id)%NodeVel_old(dir,inode)=0.0
   FSI(part_id)%NodeVel_new(dir,inode)=0.0
  enddo
  do dir=1,6
   FSI(part_id)%NodeForce_old(dir,inode)=0.0
   FSI(part_id)%NodeForce_new(dir,inode)=0.0
  enddo
 enddo ! inode=1,FSI(part_id)%NumNodes

return
end subroutine init_FSI

subroutine xdist(x1,x2,dist)
IMPLICIT NONE
 
REAL_T, dimension(3),intent(in) :: x1,x2
REAL_T, intent(out) :: dist

 dist=sqrt( (x1(1)-x2(1))**2+ &
            (x1(2)-x2(2))**2+ &
            (x1(3)-x2(3))**2 )

return
end subroutine xdist


subroutine xdist_project(x1,x2,part_id,dist)
IMPLICIT NONE

INTEGER_T, intent(in) :: part_id
REAL_T, dimension(3),intent(in) :: x1,x2
REAL_T, intent(out) :: dist
INTEGER_T sdim_local
INTEGER_T dir

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%flag_2D_to_3D.eq.1) then
  sdim_local=2
 else if (FSI(part_id)%flag_2D_to_3D.eq.0) then
  sdim_local=3
 else
  print *,"FSI(part_id)%flag_2D_to_3D invalid"
  stop
 endif
 dist=zero
 do dir=1,sdim_local
  dist=dist+(x1(dir)-x2(dir))**2
 enddo
 if (dist.ge.zero) then
  dist=sqrt(dist)
 else
  print *,"dist bust"
  stop
 endif

return
end subroutine xdist_project


subroutine xdistmin(x1,x2,dist)
IMPLICIT NONE

REAL_T, dimension(3),intent(in) :: x1,x2
REAL_T, intent(out) :: dist
REAL_T, dimension(3) :: diff
INTEGER_T :: dir

 do dir=1,3
  diff(dir)=abs(x2(dir)-x1(dir))
 enddo
 dist=diff(1)
 if (dist.gt.diff(2)) then
  dist=diff(2)
 endif 
 if (dist.gt.diff(3)) then
  dist=diff(3)
 endif 

return
end subroutine xdistmin

! split triangles so that size is no bigger than "h_small"
subroutine generate_new_triangles(initflag,problo,probhi, &
  part_id,ioproc,isout,h_small)
IMPLICIT NONE

INTEGER_T :: part_id
REAL_T problo(3),probhi(3)
INTEGER_T :: initflag,ioproc,isout
INTEGER_T :: ielem,nodes_per_elem,inode,i,dir
INTEGER_T :: new_NumIntElems
REAL_T :: generate_time
REAL_T :: h_small
REAL_T, dimension(3) :: x1,x2,x3
REAL_T, dimension(3) :: vel1,vel2,vel3
REAL_T, dimension(6) :: force1,force2,force3,forcesplit
REAL_T, dimension(3) :: xsplit,velsplit
REAL_T :: biggest_h
REAL_T :: smallest_h
INTEGER_T :: first_measure
REAL_T :: temp_h
REAL_T :: mass1,mass2,mass3,mass_split
REAL_T :: d12,d23,d13
REAL_T :: temp1,temp2,temp3
REAL_T :: tempsplit
INTEGER_T :: iter
REAL_T, dimension(3) :: normal
INTEGER_T :: ilevel
INTEGER_T :: ilevel_current
INTEGER_T :: nsplit,esplit,isub,normal_cnt,base_ielem
INTEGER_T :: node1,node2,node3
INTEGER_T :: n_lag_levels
INTEGER_T :: save_n_elems,save_n_nodes
INTEGER_T :: local_refine_factor

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif
 if (h_small.le.zero) then
  print *,"h_small invalid"
  stop
 endif

 if ((ioproc.ne.1).and.(ioproc.ne.0)) then
  print *,"ioproc invalid"
  stop
 endif

 generate_time=zero

 local_refine_factor=FSI(part_id)%refine_factor

 if ((initflag.eq.1).and.(isout.eq.1)) then
  if (ioproc.eq.1) then
   print *,"in generate_new_triangles..."
  endif
 endif

! START OF LAGRANGIAN REFINEMENT SECTION ---------------------

 biggest_h=0.0
 smallest_h=0.0
 first_measure=0

 new_NumIntElems=0
 do ielem=1,FSI(part_id)%NumIntElems
  nodes_per_elem=FSI(part_id)%ElemData(1,ielem)
  do isub=0,nodes_per_elem-3
   new_NumIntElems=new_NumIntElems+1
   node1=FSI(part_id)%IntElem(1,ielem)
   node2=FSI(part_id)%IntElem(nodes_per_elem-isub-1,ielem)
   node3=FSI(part_id)%IntElem(nodes_per_elem-isub,ielem)
   if ((node1.gt.FSI(part_id)%NumNodes).or. &
       (node2.gt.FSI(part_id)%NumNodes).or. &
       (node3.gt.FSI(part_id)%NumNodes)) then
    print *,"node1,node2, or node3 invalid ",node1,node2,node3
    print *,"ielem=",ielem
    print *,"FSI(part_id)%NumNodes ",FSI(part_id)%NumNodes
    print *,"FSI(part_id)%NumIntElems ",FSI(part_id)%NumIntElems
    print *,"nodes_per_elem=",nodes_per_elem
    print *,"part_id=",part_id
    stop
   endif
    
   do dir=1,3
    x1(dir)=FSI(part_id)%Node(dir,node1)
    x2(dir)=FSI(part_id)%Node(dir,node2)
    x3(dir)=FSI(part_id)%Node(dir,node3)
    if ((abs(x1(dir)).lt.CTMLoverflow).and. &
        (abs(x2(dir)).lt.CTMLoverflow).and. &
        (abs(x3(dir)).lt.CTMLoverflow)) then
     ! do nothing
    else 
     print *,"x1,x2, or x3 overflow"
     print *,"dir,part_id,node1,node2,node3 ", &
      dir,part_id,node1,node2,node3
     print *,"x1,x2,x3 ",x1(dir),x2(dir),x3(dir)
     stop
    endif
   enddo ! dir=1..3

   call xdist_project(x1,x2,part_id,d12)
   call xdist_project(x2,x3,part_id,d23)
   call xdist_project(x3,x1,part_id,d13)
   if (first_measure.eq.0) then
    biggest_h=d12
    smallest_h=d12
   endif
   first_measure=1

   if (biggest_h.lt.d12) then
    biggest_h=d12
   endif
   if (biggest_h.lt.d23) then
    biggest_h=d23
   endif
   if (biggest_h.lt.d13) then
    biggest_h=d13
   endif

   if (smallest_h.gt.d12) then
    smallest_h=d12
   endif
   if (smallest_h.gt.d23) then
    smallest_h=d23
   endif
   if (smallest_h.gt.d13) then
    smallest_h=d13
   endif

  enddo ! isub=0..nodes_per_elem-3
 enddo ! ielem=1..FSI(part_id)%NumIntElems

 FSI(part_id)%max_side_len=biggest_h  
 FSI(part_id)%min_side_len=smallest_h  
 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"part_id,max_side_len,min_side_len ",part_id, &
   FSI(part_id)%max_side_len,FSI(part_id)%min_side_len
 endif

 if (local_refine_factor.eq.0) then
  n_lag_levels=2
 else if (local_refine_factor.ge.1) then
  n_lag_levels=1
  temp_h=biggest_h
  do while (temp_h.gt.(local_refine_factor-0.01)*h_small)
   temp_h=half*temp_h
   n_lag_levels=n_lag_levels+2  ! two sides might exceed limit
  enddo
 else
  print *,"local_refine_factor invalid"
  stop
 endif
 if (n_lag_levels.eq.1) then
  n_lag_levels=2
 else if (n_lag_levels.gt.1) then
  ! do nothing
 else
  print *,"n_lag_levels invalid"
  stop
 endif

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"original biggest_h ",biggest_h
  print *,"h_small ",h_small
  print *,"local_refine_factor ",local_refine_factor
  print *,"number of refinements needed for nodes/elements ",n_lag_levels
 endif

 allocate(multi_lag(n_lag_levels))
 multi_lag(1)%n_nodes=FSI(part_id)%NumNodes
 multi_lag(1)%n_elems=new_NumIntElems
 allocate(multi_lag(1)%nd(3,FSI(part_id)%NumNodes))
 allocate(multi_lag(1)%ndvel(3,FSI(part_id)%NumNodes))
 allocate(multi_lag(1)%ndforce(6,FSI(part_id)%NumNodes))
 allocate(multi_lag(1)%ndmass(FSI(part_id)%NumNodes))
 allocate(multi_lag(1)%ndtemp(FSI(part_id)%NumNodes))
 allocate(multi_lag(1)%elemdt(3,new_NumIntElems))
 allocate(multi_lag(1)%intelemdt(4,new_NumIntElems))
 do i=1,FSI(part_id)%NumNodes
  do dir=1,3
   multi_lag(1)%nd(dir,i)=FSI(part_id)%Node(dir,i)
   multi_lag(1)%ndvel(dir,i)=FSI(part_id)%NodeVel(dir,i)
  enddo
  multi_lag(1)%ndmass(i)=FSI(part_id)%NodeMass(i)
  do dir=1,6
   multi_lag(1)%ndforce(dir,i)=FSI(part_id)%NodeForce(dir,i)
  enddo
  multi_lag(1)%ndtemp(i)=FSI(part_id)%NodeTemp(i)
 enddo
  
 new_NumIntElems=0
 do ielem=1,FSI(part_id)%NumIntElems
  nodes_per_elem=FSI(part_id)%ElemData(1,ielem)
  do isub=0,nodes_per_elem-3
   new_NumIntElems=new_NumIntElems+1
   multi_lag(1)%intelemdt(1,new_NumIntElems)=FSI(part_id)%IntElem(1,ielem)
   multi_lag(1)%intelemdt(2,new_NumIntElems)= &
     FSI(part_id)%IntElem(nodes_per_elem-isub-1,ielem)
   multi_lag(1)%intelemdt(3,new_NumIntElems)= &
     FSI(part_id)%IntElem(nodes_per_elem-isub,ielem)
   multi_lag(1)%intelemdt(4,new_NumIntElems)=ielem
    ! FSI(part_id)%ElemData(3,ielem) is the doubly wetted flag
   do dir=1,3
    multi_lag(1)%elemdt(dir,new_NumIntElems)=FSI(part_id)%ElemData(dir,ielem)
   enddo
    ! always 3 nodes per element for refined surface
   multi_lag(1)%elemdt(1,new_NumIntElems)=3
  enddo ! isub=0..nodes_per_elem-3
 enddo ! ielem=1..FSI(part_id)%NumIntElems

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"numintelems,numnodes ",FSI(part_id)%NumIntElems,FSI(part_id)%NumNodes
  print *,"numintelems(3node),numnodes ",new_NumIntElems,FSI(part_id)%NumNodes
 endif

 do ilevel=1,n_lag_levels-1
  do iter=1,2
   if (iter.eq.1) then
    ilevel_current=ilevel
   else if (iter.eq.2) then
    ilevel_current=ilevel+1
    save_n_elems=multi_lag(ilevel+1)%n_elems
    save_n_nodes=multi_lag(ilevel+1)%n_nodes 
    allocate(multi_lag(ilevel+1)%nd(3,save_n_nodes))
    allocate(multi_lag(ilevel+1)%ndvel(3,save_n_nodes))
    allocate(multi_lag(ilevel+1)%ndforce(6,save_n_nodes))
    allocate(multi_lag(ilevel+1)%ndmass(save_n_nodes))
    allocate(multi_lag(ilevel+1)%ndtemp(save_n_nodes))
    allocate(multi_lag(ilevel+1)%elemdt(3,save_n_elems))
    allocate(multi_lag(ilevel+1)%intelemdt(4,save_n_elems))
    do i=1,multi_lag(ilevel)%n_nodes
     do dir=1,3
      multi_lag(ilevel+1)%nd(dir,i)=multi_lag(ilevel)%nd(dir,i)
      multi_lag(ilevel+1)%ndvel(dir,i)=multi_lag(ilevel)%ndvel(dir,i)
     enddo
     multi_lag(ilevel+1)%ndmass(i)=multi_lag(ilevel)%ndmass(i)
     do dir=1,6
      multi_lag(ilevel+1)%ndforce(dir,i)=multi_lag(ilevel)%ndforce(dir,i)
     enddo
     multi_lag(ilevel+1)%ndtemp(i)=multi_lag(ilevel)%ndtemp(i)
    enddo
   else
    print *,"iter invalid"
    stop
   endif  

   multi_lag(ilevel+1)%n_elems=0
   multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel)%n_nodes
 
   do ielem=1,multi_lag(ilevel)%n_elems
    node1=multi_lag(ilevel)%intelemdt(1,ielem)
    node2=multi_lag(ilevel)%intelemdt(2,ielem)
    node3=multi_lag(ilevel)%intelemdt(3,ielem)
    do dir=1,3
     x1(dir)=multi_lag(ilevel)%nd(dir,node1)
     x2(dir)=multi_lag(ilevel)%nd(dir,node2)
     x3(dir)=multi_lag(ilevel)%nd(dir,node3)
     vel1(dir)=multi_lag(ilevel)%ndvel(dir,node1)
     vel2(dir)=multi_lag(ilevel)%ndvel(dir,node2)
     vel3(dir)=multi_lag(ilevel)%ndvel(dir,node3)
    enddo
    mass1=multi_lag(ilevel_current)%ndmass(node1)
    mass2=multi_lag(ilevel_current)%ndmass(node2)
    mass3=multi_lag(ilevel_current)%ndmass(node3)
    do dir=1,6
     force1(dir)=multi_lag(ilevel)%ndforce(dir,node1)
     force2(dir)=multi_lag(ilevel)%ndforce(dir,node2)
     force3(dir)=multi_lag(ilevel)%ndforce(dir,node3)
    enddo
    temp1=multi_lag(ilevel)%ndtemp(node1)
    temp2=multi_lag(ilevel)%ndtemp(node2)
    temp3=multi_lag(ilevel)%ndtemp(node3)
    call xdist_project(x1,x2,part_id,d12)
    call xdist_project(x2,x3,part_id,d23)
    call xdist_project(x3,x1,part_id,d13)

    if ((local_refine_factor.gt.0).and. &
        (d12.ge.(local_refine_factor-0.01)*h_small).and. &
        (d12.ge.d23).and. &
        (d12.ge.d13)) then
     do dir=1,3
      xsplit(dir)=0.5*(x1(dir)+x2(dir))
      velsplit(dir)=0.5*(vel1(dir)+vel2(dir))
     enddo
     do dir=1,6
      forcesplit(dir)=0.5*(force1(dir)+force2(dir))
     enddo

     tempsplit=0.5*(temp1+temp2)
     mass_split=(mass1+mass2)/three  ! criterion: mass is conserved

     multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel+1)%n_nodes+1
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+2
     nsplit=multi_lag(ilevel+1)%n_nodes
     esplit=multi_lag(ilevel+1)%n_elems

     if (iter.eq.2) then
      do dir=1,3
       multi_lag(ilevel+1)%nd(dir,nsplit)=xsplit(dir)
       multi_lag(ilevel+1)%ndvel(dir,nsplit)=velsplit(dir)
      enddo 
      do dir=1,6
       multi_lag(ilevel+1)%ndforce(dir,nsplit)=forcesplit(dir)
      enddo
      multi_lag(ilevel+1)%ndmass(nsplit)=mass_split
      multi_lag(ilevel+1)%ndmass(node1)=two*mass1/three
      multi_lag(ilevel+1)%ndmass(node2)=two*mass2/three

      multi_lag(ilevel+1)%ndtemp(nsplit)=tempsplit

      base_ielem=multi_lag(ilevel)%intelemdt(4,ielem)
      multi_lag(ilevel+1)%intelemdt(1,esplit-1)=node1
      multi_lag(ilevel+1)%intelemdt(2,esplit-1)=nsplit
      multi_lag(ilevel+1)%intelemdt(3,esplit-1)=node3
      multi_lag(ilevel+1)%intelemdt(4,esplit-1)=base_ielem
       ! elemdt(3,ielem) is the doubly wetted flag
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit-1)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo

      multi_lag(ilevel+1)%intelemdt(1,esplit)=nsplit
      multi_lag(ilevel+1)%intelemdt(2,esplit)=node2
      multi_lag(ilevel+1)%intelemdt(3,esplit)=node3
      multi_lag(ilevel+1)%intelemdt(4,esplit)=base_ielem
       ! elemdt(3,ielem) is the doubly wetted flag
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo
     endif ! iter.eq.2
    else if ((local_refine_factor.gt.0).and. &
             (d23.ge.(local_refine_factor-0.01)*h_small).and. &
             (d23.ge.d12).and. &
             (d23.ge.d13)) then
     do dir=1,3
      xsplit(dir)=0.5*(x2(dir)+x3(dir))
      velsplit(dir)=0.5*(vel2(dir)+vel3(dir))
     enddo
     do dir=1,6
      forcesplit(dir)=0.5*(force2(dir)+force3(dir))
     enddo
 
     tempsplit=0.5*(temp2+temp3)
     mass_split=(mass2+mass3)/three ! criterion: mass is conserved

     multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel+1)%n_nodes+1
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+2
     nsplit=multi_lag(ilevel+1)%n_nodes
     esplit=multi_lag(ilevel+1)%n_elems

     if (iter.eq.2) then
      do dir=1,3
       multi_lag(ilevel+1)%nd(dir,nsplit)=xsplit(dir)
       multi_lag(ilevel+1)%ndvel(dir,nsplit)=velsplit(dir)
      enddo 
      do dir=1,6
       multi_lag(ilevel+1)%ndforce(dir,nsplit)=forcesplit(dir)
      enddo
      multi_lag(ilevel+1)%ndmass(nsplit)=mass_split
      multi_lag(ilevel+1)%ndmass(node2)=two*mass2/three
      multi_lag(ilevel+1)%ndmass(node3)=two*mass3/three

      multi_lag(ilevel+1)%ndtemp(nsplit)=tempsplit

      base_ielem=multi_lag(ilevel)%intelemdt(4,ielem)
      multi_lag(ilevel+1)%intelemdt(1,esplit-1)=node1
      multi_lag(ilevel+1)%intelemdt(2,esplit-1)=node2
      multi_lag(ilevel+1)%intelemdt(3,esplit-1)=nsplit
      multi_lag(ilevel+1)%intelemdt(4,esplit-1)=base_ielem
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit-1)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo

      multi_lag(ilevel+1)%intelemdt(1,esplit)=nsplit
      multi_lag(ilevel+1)%intelemdt(2,esplit)=node3
      multi_lag(ilevel+1)%intelemdt(3,esplit)=node1
      multi_lag(ilevel+1)%intelemdt(4,esplit)=base_ielem
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo
     endif ! iter.eq.2
    else if ((local_refine_factor.gt.0).and. &
             (d13.ge.(local_refine_factor-0.01)*h_small).and. &
             (d13.ge.d12).and. &
             (d13.ge.d23)) then
     do dir=1,3
      xsplit(dir)=0.5*(x1(dir)+x3(dir))
      velsplit(dir)=0.5*(vel1(dir)+vel3(dir))
     enddo
     do dir=1,6
      forcesplit(dir)=0.5*(force1(dir)+force3(dir))
     enddo

     tempsplit=0.5*(temp1+temp3)
     mass_split=(mass1+mass3)/three ! criterion: mass is conserved

     multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel+1)%n_nodes+1
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+2
     nsplit=multi_lag(ilevel+1)%n_nodes
     esplit=multi_lag(ilevel+1)%n_elems

     if (iter.eq.2) then
      do dir=1,3
       multi_lag(ilevel+1)%nd(dir,nsplit)=xsplit(dir)
       multi_lag(ilevel+1)%ndvel(dir,nsplit)=velsplit(dir)
      enddo 
      do dir=1,6
       multi_lag(ilevel+1)%ndforce(dir,nsplit)=forcesplit(dir)
      enddo
      multi_lag(ilevel+1)%ndmass(nsplit)=mass_split
      multi_lag(ilevel+1)%ndmass(node1)=two*mass1/three
      multi_lag(ilevel+1)%ndmass(node3)=two*mass3/three

      multi_lag(ilevel+1)%ndtemp(nsplit)=tempsplit

      base_ielem=multi_lag(ilevel)%intelemdt(4,ielem)
      multi_lag(ilevel+1)%intelemdt(1,esplit-1)=node1
      multi_lag(ilevel+1)%intelemdt(2,esplit-1)=node2
      multi_lag(ilevel+1)%intelemdt(3,esplit-1)=nsplit
      multi_lag(ilevel+1)%intelemdt(4,esplit-1)=base_ielem
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit-1)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo

      multi_lag(ilevel+1)%intelemdt(1,esplit)=node2
      multi_lag(ilevel+1)%intelemdt(2,esplit)=node3
      multi_lag(ilevel+1)%intelemdt(3,esplit)=nsplit
      multi_lag(ilevel+1)%intelemdt(4,esplit)=base_ielem
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo
     endif ! iter.eq.2
    else
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+1
     esplit=multi_lag(ilevel+1)%n_elems
     if (iter.eq.2) then
      base_ielem=multi_lag(ilevel)%intelemdt(4,ielem)
      multi_lag(ilevel+1)%intelemdt(1,esplit)=node1
      multi_lag(ilevel+1)%intelemdt(2,esplit)=node2
      multi_lag(ilevel+1)%intelemdt(3,esplit)=node3
      multi_lag(ilevel+1)%intelemdt(4,esplit)=base_ielem
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo
     endif  ! iter.eq.2
    endif
   enddo ! ielem
  enddo ! iter
  if ((ioproc.eq.1).and.(isout.eq.1)) then
   print *,"ilevel+1,numintelems,numnodes ",ilevel+1,  &
    multi_lag(ilevel+1)%n_elems,multi_lag(ilevel+1)%n_nodes
  endif
 enddo ! ilevel

 FSI(part_id)%NumNodesBIG=multi_lag(n_lag_levels)%n_nodes
 FSI(part_id)%NumIntElemsBIG=multi_lag(n_lag_levels)%n_elems

  ! in: generate_new_triangles
 if (initflag.eq.0) then 
  deallocate(FSI(part_id)%NodeBIG)
  deallocate(FSI(part_id)%NodeVelBIG)
  deallocate(FSI(part_id)%NodeForceBIG)
  deallocate(FSI(part_id)%NodeMassBIG)
  deallocate(FSI(part_id)%NodeTempBIG)
  deallocate(FSI(part_id)%NodeNormalBIG)
  deallocate(FSI(part_id)%ElemNodeCountBIG)

  deallocate(FSI(part_id)%ElemDataXnotBIG)
  deallocate(FSI(part_id)%ElemDataBIG)
  deallocate(FSI(part_id)%IntElemBIG)

 endif

  ! in: generate_new_triangles
 allocate(FSI(part_id)%NodeBIG(3,FSI(part_id)%NumNodesBIG))
 allocate(FSI(part_id)%NodeNormalBIG(3,FSI(part_id)%NumNodesBIG))
 allocate(FSI(part_id)%NodeVelBIG(3,FSI(part_id)%NumNodesBIG))
 allocate(FSI(part_id)%NodeForceBIG(6,FSI(part_id)%NumNodesBIG))
 allocate(FSI(part_id)%NodeMassBIG(FSI(part_id)%NumNodesBIG))
 allocate(FSI(part_id)%NodeTempBIG(FSI(part_id)%NumNodesBIG))
 allocate(FSI(part_id)%ElemNodeCountBIG(FSI(part_id)%NumNodesBIG))
 allocate(FSI(part_id)%ElemDataXnotBIG(3,FSI(part_id)%NumIntElemsBIG))
 allocate(FSI(part_id)%ElemDataBIG(3,FSI(part_id)%NumIntElemsBIG))
 allocate(FSI(part_id)%IntElemBIG(4,FSI(part_id)%NumIntElemsBIG))

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"NumNodes, NumIntElems ",FSI(part_id)%NumNodes, &
   FSI(part_id)%NumIntElems
  print *,"NumNodesBIG, NumIntElemsBIG ",FSI(part_id)%NumNodesBIG, &
   FSI(part_id)%NumIntElemsBIG
 endif

  ! in: generate_new_triangles
 do ielem=1,FSI(part_id)%NumIntElemsBIG
  do dir=1,4
   FSI(part_id)%IntElemBIG(dir,ielem)= &
     multi_lag(n_lag_levels)%intelemdt(dir,ielem)
  enddo
   ! ElemDataBIG(3,ielem) is the doubly wetted flag
  do dir=1,3
   FSI(part_id)%ElemDataBIG(dir,ielem)=multi_lag(n_lag_levels)%elemdt(dir,ielem)
  enddo
 enddo
  ! in: generate_new_triangles
 do inode=1,FSI(part_id)%NumNodesBIG
  do dir=1,3
   FSI(part_id)%NodeNormalBIG(dir,inode)=0.0
   FSI(part_id)%NodeBIG(dir,inode)=multi_lag(n_lag_levels)%nd(dir,inode)
   FSI(part_id)%NodeVelBIG(dir,inode)=multi_lag(n_lag_levels)%ndvel(dir,inode)
  enddo
  do dir=1,6
   FSI(part_id)%NodeForceBIG(dir,inode)= &
     multi_lag(n_lag_levels)%ndforce(dir,inode)
  enddo
  FSI(part_id)%NodeMassBIG(inode)=multi_lag(n_lag_levels)%ndmass(inode)
  FSI(part_id)%NodeTempBIG(inode)=multi_lag(n_lag_levels)%ndtemp(inode)
  FSI(part_id)%ElemNodeCountBIG(inode)=0
 enddo ! inode=1,FSI(part_id)%NumNodesBIG

 do i=1,n_lag_levels
  deallocate(multi_lag(i)%intelemdt)
  deallocate(multi_lag(i)%elemdt)
  deallocate(multi_lag(i)%nd)
  deallocate(multi_lag(i)%ndvel)
  deallocate(multi_lag(i)%ndforce)
  deallocate(multi_lag(i)%ndmass)
  deallocate(multi_lag(i)%ndtemp)
 enddo
 deallocate(multi_lag)
 
 if ((ioproc.eq.1).and.(isout.eq.1)) then 
  print *,"creating normals for the nodes and Xnot"
 endif

 do ielem=1,FSI(part_id)%NumIntElemsBIG
  nodes_per_elem=FSI(part_id)%ElemDataBIG(1,ielem)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem invalid"
   stop
  endif
  do dir=1,3
   FSI(part_id)%ElemDataXnotBIG(dir,ielem)=0.0
  enddo
  call scinormal(ielem,normal,part_id,generate_time)
  do i=1,nodes_per_elem
   inode=FSI(part_id)%IntElemBIG(i,ielem)
   do dir=1,3
    FSI(part_id)%NodeNormalBIG(dir,inode)= &
      FSI(part_id)%NodeNormalBIG(dir,inode)+normal(dir)
    FSI(part_id)%ElemDataXnotBIG(dir,ielem)= &
      FSI(part_id)%ElemDataXnotBIG(dir,ielem)+ &
      FSI(part_id)%NodeBIG(dir,inode)
   enddo
   FSI(part_id)%ElemNodeCountBIG(inode)= &
     FSI(part_id)%ElemNodeCountBIG(inode)+1
  enddo ! i=1,nodes_per_elem
  do dir=1,3
   FSI(part_id)%ElemDataXnotBIG(dir,ielem)= &
    FSI(part_id)%ElemDataXnotBIG(dir,ielem)/3.0
  enddo
 enddo

 do inode=1,FSI(part_id)%NumNodesBIG
  normal_cnt=FSI(part_id)%ElemNodeCountBIG(inode)
  if (normal_cnt.gt.0) then
   do dir=1,3
    FSI(part_id)%NodeNormalBIG(dir,inode)= &
     FSI(part_id)%NodeNormalBIG(dir,inode)/normal_cnt
   enddo
  endif
 enddo

 biggest_h=0.0
 smallest_h=0.0
 first_measure=0

 do ielem=1,FSI(part_id)%NumIntElemsBIG
  nodes_per_elem=FSI(part_id)%ElemDataBIG(1,ielem)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem invalid"
   stop
  endif
  node1=FSI(part_id)%IntElemBIG(1,ielem)
  node2=FSI(part_id)%IntElemBIG(2,ielem)
  node3=FSI(part_id)%IntElemBIG(3,ielem)
  do dir=1,3
   x1(dir)=FSI(part_id)%NodeBIG(dir,node1) 
   x2(dir)=FSI(part_id)%NodeBIG(dir,node2)
   x3(dir)=FSI(part_id)%NodeBIG(dir,node3)

   if ((abs(x1(dir)).lt.CTMLoverflow).and. &
       (abs(x2(dir)).lt.CTMLoverflow).and. &
       (abs(x3(dir)).lt.CTMLoverflow)) then
    ! do nothing
   else 
    print *,"x1,x2, or x3 overflow (BIG)"
    print *,"dir,part_id,node1,node2,node3 ", &
      dir,part_id,node1,node2,node3
    print *,"x1,x2,x3 ",x1(dir),x2(dir),x3(dir)
    stop
   endif

  enddo ! dir=1..3

  call xdist_project(x1,x2,part_id,d12)
  call xdist_project(x2,x3,part_id,d23)
  call xdist_project(x3,x1,part_id,d13)
  if (first_measure.eq.0) then
   biggest_h=d12
   smallest_h=d12
  endif
  first_measure=1

  if (biggest_h.lt.d12) then
   biggest_h=d12
  endif
  if (biggest_h.lt.d23) then
   biggest_h=d23
  endif
  if (biggest_h.lt.d13) then
   biggest_h=d13
  endif

  if (smallest_h.gt.d12) then
   smallest_h=d12
  endif
  if (smallest_h.gt.d23) then
   smallest_h=d23
  endif
  if (smallest_h.gt.d13) then
   smallest_h=d13
  endif
 enddo ! ielem=1..FSI(part_id)%NumIntElemsBIG

 FSI(part_id)%max_side_len_refined=biggest_h  
 FSI(part_id)%min_side_len_refined=smallest_h  

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"part_id,flag_2D_to_3D ",part_id,FSI(part_id)%flag_2D_to_3D
  print *,"part_id,max_side_len_refined,min_side_len_refined ",part_id, &
   FSI(part_id)%max_side_len_refined,FSI(part_id)%min_side_len_refined
  print *,"local_refine_factor ",local_refine_factor
  print *,"h_small ",h_small
 endif

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print*,"in allocate: part ID ",FSI(part_id)%PartID
 endif

! END OF LAGRANGIAN REFINEMENT SECTION -----------------------

return
end subroutine generate_new_triangles


subroutine scihandoffset(ofs,time)
IMPLICIT NONE

REAL_T :: ofs,time

 ofs=0.0
 if (time.le.0.4) then
  ofs=-time/8.0
 else if ((time.ge.0.4).and.(time.lt.0.8)) then
  ofs=-0.05
 else if ((time.ge.0.8).and.(time.lt.1.2)) then
  ofs=(time-0.9)/2.0
 else
  ofs=0.15
 endif

return
end subroutine scihandoffset

subroutine CTML_init_sci_node(ioproc,part_id,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T :: isout
INTEGER_T :: inode,inode_fiber,ioproc
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1,xval2
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
INTEGER_T :: dir

REAL_T, dimension(3) :: xxblob1,newxxblob1,xxblob2,newxxblob2
REAL_T, dimension(3) :: vel_local
REAL_T, dimension(6) :: stress_local
REAL_T :: mass_local
REAL_T :: radradblob1,radradblob2
INTEGER_T :: nmat
INTEGER_T :: ctml_part_id

  nmat=num_materials

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id out of range, part_id, TOTAL_NPARTS:",part_id,TOTAL_NPARTS
   stop
  endif
  if (FSI(part_id)%partID.ne.part_id) then
   print *,"FSI(part_id)%partID.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
   print *,"ioproc invalid"
   stop
  endif

  ctml_part_id=CTML_partid_map(part_id)

  if ((ctml_part_id.ge.1).and. &
      (ctml_part_id.le.CTML_NPARTS)) then

   xxblob1(1)=0.0
   xxblob1(2)=0.0
   xxblob1(3)=0.0
   newxxblob1(1)=0.0
   newxxblob1(2)=0.0
   newxxblob1(3)=0.0

   xxblob2(1)=0.0
   xxblob2(2)=0.0
   xxblob2(3)=0.0
   newxxblob2(1)=0.0
   newxxblob2(2)=0.0
   newxxblob2(3)=0.0
   radradblob1=1.0
   radradblob2=1.0

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   do inode=1,FSI(part_id)%NumNodes

    if (AMREX_SPACEDIM.eq.2) then
     if ((inode.ge.1).and. &
         (inode.le.ctml_n_fib_active_nodes(ctml_part_id))) then
      inode_fiber=inode
     else if ((inode.ge.ctml_n_fib_active_nodes(ctml_part_id)+1).and. &
              (inode.le.2*ctml_n_fib_active_nodes(ctml_part_id))) then
      inode_fiber=inode-ctml_n_fib_active_nodes(ctml_part_id)
     else
      print *,"node invalid"
      stop
     endif
    else if (AMREX_SPACEDIM.eq.3) then
     inode_fiber=inode
    else
     print *,"dimension bust"
     stop
    endif  
    do dir=1,AMREX_SPACEDIM
     xval(dir)=ctml_fib_pst(ctml_part_id,inode_fiber,dir)
     vel_local(dir)=ctml_fib_vel(ctml_part_id,inode_fiber,dir)
    enddo
    do dir=1,6
     stress_local(dir)=zero
    enddo
    do dir=1,2*AMREX_SPACEDIM
     stress_local(dir)=ctml_fib_frc(ctml_part_id,inode_fiber,dir)
    enddo
    mass_local=ctml_fib_mass(ctml_part_id,inode_fiber)
    if (mass_local.gt.zero) then
     ! do nothing
    else
     print *,"mass_local invalid, mass_local=",mass_local
     stop
    endif

    if (AMREX_SPACEDIM.eq.3) then
     ! do nothing
    else if (AMREX_SPACEDIM.eq.2) then
     vel_local(3)=zero
     if (inode.eq.inode_fiber) then
      xval(3)=problo_act(3)+problen_act(3)/ten
     else if (inode.eq.inode_fiber+ctml_n_fib_active_nodes(ctml_part_id)) then
      xval(3)=probhi_act(3)-problen_act(3)/ten
     else
      print *,"inode invalid"
      stop
     endif
    else
     print *,"dimension bust"
     stop
    endif
     
    do dir=1,3
     if ((minnodebefore(dir).gt.xval(dir)).or.(inode.eq.1)) then
      minnodebefore(dir)=xval(dir)
     endif
     if ((maxnodebefore(dir).lt.xval(dir)).or.(inode.eq.1)) then
      maxnodebefore(dir)=xval(dir)
     endif
    enddo

    do dir=1,3
     xval1(dir)=(xval(dir)-xxblob1(dir))/radradblob1 + newxxblob1(dir)
     xval2(dir)=(xval(dir)-xxblob2(dir))/radradblob2 + newxxblob2(dir)
    enddo

    do dir=1,3
     if ((minnode(dir).gt.xval1(dir)).or.(inode.eq.1)) then
      minnode(dir)=xval1(dir)
     endif
     if ((maxnode(dir).lt.xval1(dir)).or.(inode.eq.1)) then
      maxnode(dir)=xval1(dir)
     endif
     if ((minnode(dir).gt.xval2(dir)).or.(inode.eq.1)) then
      minnode(dir)=xval2(dir)
     endif
     if ((maxnode(dir).lt.xval2(dir)).or.(inode.eq.1)) then
      maxnode(dir)=xval2(dir)
     endif
    enddo ! dir
     
    do dir=1,3
     FSI(part_id)%Node_old(dir,inode)=xval1(dir)
     FSI(part_id)%Node_new(dir,inode)=xval1(dir)
     FSI(part_id)%Node_current(dir,inode)=xval1(dir)
     FSI(part_id)%Node(dir,inode)=xval1(dir)
     FSI(part_id)%NodeVel(dir,inode)=vel_local(dir)
     FSI(part_id)%NodeVel_old(dir,inode)=vel_local(dir)
     FSI(part_id)%NodeVel_new(dir,inode)=vel_local(dir)
    enddo
    do dir=1,6
     FSI(part_id)%NodeForce(dir,inode)=stress_local(dir)
     FSI(part_id)%NodeForce_old(dir,inode)=stress_local(dir)
     FSI(part_id)%NodeForce_new(dir,inode)=stress_local(dir)
    enddo

     ! in: CTML_init_sci_node
    FSI(part_id)%NodeMass(inode)=mass_local
    FSI(part_id)%NodeTemp(inode)=zero
    FSI(part_id)%NodeTemp_new(inode)=zero
       
   enddo  ! inode=1,NumNodes
     
   if ((ioproc.eq.1).and.(isout.eq.1)) then

    do dir=1,3 
     print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
    enddo
    do dir=1,3 
     print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
    enddo

   endif

  else
   print *,"ctml_part_id invalid"
   stop
  endif
  
return
end subroutine CTML_init_sci_node

subroutine CTML_init_sci(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  part_id,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T :: sdim,ifirst,isout
INTEGER_T :: iface,inode_base,ioproc
REAL_T :: curtime,dt
INTEGER_T :: dir,istep,istop
INTEGER_T :: nmat
INTEGER_T :: node_factor
INTEGER_T :: ctml_part_id
INTEGER_T :: inode_crit,inode

  nmat=num_materials

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id out of range, part_id, TOTAL_NPARTS:",part_id,TOTAL_NPARTS
   stop
  endif
  if (FSI(part_id)%partID.ne.part_id) then
   print *,"FSI(part_id)%partID.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
  endif

  ctml_part_id=CTML_partid_map(part_id)
  if ((ctml_part_id.ge.1).and. &
      (ctml_part_id.le.CTML_NPARTS)) then

   inode_crit=0  ! node index of first inactive node.
   inode=0
#ifdef MVAHABFSI
   call CTML_GET_POS_VEL_FORCE_WT( &
    ctml_fib_pst, &
    ctml_fib_vel, &
    ctml_fib_frc, &
    ctml_fib_mass, &
    ctml_n_fib_bodies, &
    ctml_max_n_fib_nodes, &
    ctml_part_id)

   do inode=1,ctml_n_fib_nodes(ctml_part_id)
    if (ctml_fib_mass(ctml_part_id,inode).gt.zero) then
     ! do nothing
    else if (ctml_fib_mass(ctml_part_id,inode).eq.zero) then 
     if (inode_crit.eq.0) then
      inode_crit=inode
     endif
    else
     print *,"ctml_fib_mass(ctml_part_id,inode) invalid"
     stop
    endif
   enddo ! inode=1,ctml_n_fib_nodes(ctml_part_id)
   if (inode_crit.eq.0) then
    ctml_n_fib_active_nodes(ctml_part_id)= &
      ctml_n_fib_nodes(ctml_part_id)
   else if (inode_crit.gt.1) then
    print *,"WARNING:inode_crit>1  inode_crit=",inode_crit
    ctml_n_fib_active_nodes(ctml_part_id)=inode_crit-1
    print *,"WARNING:ctml_part_id =",ctml_part_id
    print *,"WARNING:ctml_n_fib_active_nodes =", &
     ctml_n_fib_active_nodes(ctml_part_id)
   else
    print *,"inode_crit invalid"
    stop
   endif
#else
   print *,"define MVAHABFSI"
   stop
#endif

   timeB=curtime
   dtB=0.0
   TorquePos=0.0
   TorqueVel=0.0

   if ((ifirst.eq.0).or.(ifirst.eq.1)) then
    ! do nothing
   else
    print *,"ifirst invalid"
    stop
   endif

   if (AMREX_SPACEDIM.eq.3) then
    node_factor=1
    print *,"3D not supported yet"
   else if (AMREX_SPACEDIM.eq.2) then
    node_factor=2
   else
    print *,"dimension bust"
    stop
   endif
   FSI(part_id)%NumNodes=ctml_n_fib_active_nodes(ctml_part_id)*node_factor
    ! 2 triangular elements make up a rectangle.
   FSI(part_id)%NumIntElems= &
    (ctml_n_fib_active_nodes(ctml_part_id)-1)*node_factor
   FSI(part_id)%IntElemDim=3
 
    ! in: CTML_init_sci 
   if (ifirst.eq.1) then

     ! allocate_intelem=1
     ! allocates NodeMass
    call init_FSI(part_id,1)  

    allocate(FSI(part_id)%Node(3,FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%NodeVel(3,FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%NodeForce(6,FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%NodeTemp(FSI(part_id)%NumNodes))
   else
    print *,"something wrong, ifirst should be 1 here"
    stop
   endif  ! ifirst.eq.1

   call CTML_init_sci_node(ioproc,part_id,isout)

   do iface=1,FSI(part_id)%NumIntElems
    if (AMREX_SPACEDIM.eq.3) then
     print *,"3d not ready yet"
     stop
    else if (AMREX_SPACEDIM.eq.2) then
     if (iface.le.ctml_n_fib_active_nodes(ctml_part_id)-1) then
      inode_base=iface
      FSI(part_id)%IntElem(1,iface)=inode_base
      FSI(part_id)%IntElem(2,iface)= &
       inode_base+ctml_n_fib_active_nodes(ctml_part_id)
      FSI(part_id)%IntElem(3,iface)=inode_base+1
     else if ((iface.ge.ctml_n_fib_active_nodes(ctml_part_id)).and. &
              (iface.le.FSI(part_id)%NumIntElems)) then
      inode_base=iface+1-ctml_n_fib_active_nodes(ctml_part_id)
      FSI(part_id)%IntElem(1,iface)=inode_base+1
      FSI(part_id)%IntElem(2,iface)= &
       inode_base+ctml_n_fib_active_nodes(ctml_part_id)
      FSI(part_id)%IntElem(3,iface)= &
       inode_base+1+ctml_n_fib_active_nodes(ctml_part_id)
     endif
    else
     print *,"dimension bust"
     stop
    endif
    FSI(part_id)%ElemData(1,iface)=3   ! number of nodes in element
    FSI(part_id)%ElemData(2,iface)=1   ! part number
    FSI(part_id)%ElemData(3,iface)=1   ! doubly wetted (0=singly wetted)
   enddo ! iface=1,FSI(part_id)%NumIntElems

   do dir=1,3
    FSI(part_id)%solid_displ(dir)=0.
    FSI(part_id)%solid_speed(dir)=0.
   enddo
   FSI(part_id)%soliddrop_displacement=0.0
   FSI(part_id)%soliddrop_speed=0.0

   use_temp=0

  else
   print *,"ctml_part_id invalid"
   stop
  endif
  
return
end subroutine CTML_init_sci


subroutine initinjector(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  part_id,isout)
IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T :: sdim,ifirst,isout
INTEGER_T :: inode,iface,ioproc
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1,xval2
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
INTEGER_T :: dir,istep,istop
character(40) :: dwave

REAL_T, dimension(3) :: xxblob1,newxxblob1,xxblob2,newxxblob2
REAL_T :: radradblob1,radradblob2

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id invalid"
   stop
  endif
  if (FSI(part_id)%partID.ne.part_id) then
   print *,"FSI(part_id)%partID.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
  endif

  if ((probtype.eq.536).or.(probtype.eq.537).or. &
      (probtype.eq.538).or.(probtype.eq.539).or. &
      (probtype.eq.53).or.(probtype.eq.531).or. &
      (probtype.eq.5501).or. &
      (probtype.eq.541)) then

   timeB=curtime
   dtB=0.0
   TorquePos=0.0
   TorqueVel=0.0

   if ((ifirst.eq.0).or.(ifirst.eq.1)) then
    ! do nothing
   else
    print *,"ifirst invalid"
    stop
   endif

   if (ifirst.eq.1) then
    xxblob1(1)=0.0
    xxblob1(2)=0.0
    xxblob1(3)=0.0

    xxblob2(1)=0.0
    xxblob2(2)=0.0
    xxblob2(3)=0.0
    newxxblob2(1)=1.28
    newxxblob2(2)=1.28
    newxxblob2(3)=2.56  
    radradblob2=1.0/0.00000001

    if (probtype.eq.5501) then ! rough surface
     xxblob1(1)=0.0
     xxblob1(2)=0.0
     xxblob1(3)=0.02
     newxxblob1(1)=0.0
     newxxblob1(2)=0.0
     newxxblob1(3)=0.0
     radradblob1=1.0  
    else if (probtype.eq.53) then
     if (axis_dir.eq.100) then  ! fan injector
      xxblob1(1)=0.0
      xxblob1(2)=0.0
      xxblob1(3)=1.2
      newxxblob1(1)=0.0
      newxxblob1(2)=0.0
      newxxblob1(3)=0.0
      radradblob1=1.0  
     else
      newxxblob1(1)=1.28
      newxxblob1(2)=1.28
      newxxblob1(3)=0.15  
      radradblob1=1.0/0.1
     endif
    else if (probtype.eq.531) then
     newxxblob1(1)=0.0 
     newxxblob1(2)=0.0 
     newxxblob1(3)=0.055
     radradblob1=1.0/0.0115
    else if (probtype.eq.537) then ! 6 hole injector
     xxblob1(1)=0.0
     xxblob1(2)=0.0
     xxblob1(3)=0.0
     newxxblob1(1)=0.0
     newxxblob1(2)=0.0
     newxxblob1(3)=-0.125 ! this is so that tecplot shows the right thing
                          ! when a coarse mesh is used.
     radradblob1=1.0
    else if ((probtype.eq.536).or. &
             (probtype.eq.538).or. &
             (probtype.eq.539).or. &
             (probtype.eq.541)) then

       ! MARK:
       ! IF probtype=538 and partID=2, then add 0.01 in order to shift the
       ! needle all the way in the domain (newxxblob1(3)=0.01)
       ! this shift will be taken into account in the line that has:
       ! "xfoot(3)=xfoot(3)+0.01" (get_foot_from_target, get_target_from_foot)
     xxblob1(1)=0.0
     xxblob1(2)=0.0
     xxblob1(3)=0.0
     newxxblob1(1)=0.0
     newxxblob1(2)=0.0
     newxxblob1(3)=0.0

     if ((probtype.eq.538).or.(probtype.eq.541)) then
      if (FSI(part_id)%partID.eq.2) then !want the whole needle in the domain.
       if ((AMREX_SPACEDIM.eq.3).or.(probtype.eq.538)) then
        newxxblob1(3)=0.01
       else if ((AMREX_SPACEDIM.eq.2).and.(probtype.eq.541)) then
        newxxblob1(3)=0.0
       endif
      else if (FSI(part_id)%partID.eq.1) then
       newxxblob1(3)=0.0
      else
       print *,"part id invalid"
       stop
      endif
     endif

     radradblob1=1.0

      ! why do we duplicate this geometry???
     do dir=1,3
      xxblob2(dir)=xxblob1(dir)
      newxxblob2(dir)=newxxblob1(dir)
     enddo
     radradblob2=radradblob1

    else
     if (isout.eq.1) then
      print *,"probtype invalid in initinjector"
     endif
     stop
    endif 
   
    denpaddle=1.0
    dampingpaddle=0.0

     ! choice here regardless of axis_dir
    if (FSI(part_id)%partID.eq.1) then

      if (probtype.eq.5501) then
       dwave="rough.cas"
      else if ((probtype.eq.53).and.(axis_dir.eq.100)) then
       dwave="flat_fan_s.cas"
      else
       dwave="injectorgeom.dat"
      endif

    else if (FSI(part_id)%partID.eq.2) then
      dwave="injectorgeom_needle.dat"
    else
      print *,"part id invalid"
      stop
    endif

    if ((ioproc.eq.1).and.(isout.eq.1)) then
     print *,"opening ",dwave
    endif
    OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

    READ(14,*) FSI(part_id)%NumNodes,FSI(part_id)%NumIntElems
    FSI(part_id)%NumNodes=FSI(part_id)%NumNodes*2
    FSI(part_id)%NumIntElems=FSI(part_id)%NumIntElems*2
    FSI(part_id)%IntElemDim=3
  
    if (ifirst.eq.1) then
     ! allocates NodeMass
     ! allocate_intelem=1
     call init_FSI(part_id,1)  
    else
     print *,"something wrong, ifirst should be 1 here"
     stop
    endif  ! ifirst.eq.1

    do dir=1,3
     maxnode(dir)=0.0
     minnode(dir)=0.0
     maxnodebefore(dir)=-1.0e+10
     minnodebefore(dir)=1.0e+10
    enddo

    do inode=1,FSI(part_id)%NumNodes/2
     if ((probtype.ne.541).or.(AMREX_SPACEDIM.eq.3)) then
      READ(14,*) xval(1),xval(2),xval(3)
     else if ((probtype.eq.541).and.(AMREX_SPACEDIM.eq.2)) then
      READ(14,*) xval(1),xval(3),xval(2)
     endif

     do dir=1,3
      if ((minnodebefore(dir).gt.xval(dir)).or.(inode.eq.1)) then
       minnodebefore(dir)=xval(dir)
      endif
      if ((maxnodebefore(dir).lt.xval(dir)).or.(inode.eq.1)) then
       maxnodebefore(dir)=xval(dir)
      endif
     enddo

     do dir=1,3
      xval1(dir)=(xval(dir)-xxblob1(dir))/radradblob1 + newxxblob1(dir)
      xval2(dir)=(xval(dir)-xxblob2(dir))/radradblob2 + newxxblob2(dir)
     enddo

     do dir=1,3
      if ((minnode(dir).gt.xval1(dir)).or.(inode.eq.1)) then
       minnode(dir)=xval1(dir)
      endif
      if ((maxnode(dir).lt.xval1(dir)).or.(inode.eq.1)) then
       maxnode(dir)=xval1(dir)
      endif
      if ((minnode(dir).gt.xval2(dir)).or.(inode.eq.1)) then
       minnode(dir)=xval2(dir)
      endif
      if ((maxnode(dir).lt.xval2(dir)).or.(inode.eq.1)) then
       maxnode(dir)=xval2(dir)
      endif
     enddo ! dir
     
     do dir=1,3
      FSI(part_id)%Node_old(dir,inode)=xval1(dir)
      FSI(part_id)%Node_new(dir,inode)=xval1(dir)
      FSI(part_id)%Node_old(dir,inode+FSI(part_id)%NumNodes/2)=xval2(dir)
      FSI(part_id)%Node_new(dir,inode+FSI(part_id)%NumNodes/2)=xval2(dir)
     enddo
       
    enddo  ! inode=1,NumNodes
     
    do iface=1,FSI(part_id)%NumIntElems/2
     READ(14,*) FSI(part_id)%IntElem(1,iface), &
                FSI(part_id)%IntElem(2,iface), &
                FSI(part_id)%IntElem(3,iface)
     do dir=3,1,-1
      FSI(part_id)%IntElem(dir,iface+FSI(part_id)%NumIntElems/2)= &
        FSI(part_id)%IntElem(dir,iface)+FSI(part_id)%NumNodes/2 
     enddo
    enddo

    do iface=1,FSI(part_id)%NumIntElems
     FSI(part_id)%ElemData(1,iface)=3   ! number of nodes in element
     FSI(part_id)%ElemData(2,iface)=1   ! part number
     FSI(part_id)%ElemData(3,iface)=0   ! singly wetted (1=doubly wetted)
    enddo  ! iface, looping faces

    close(14)

    if ((ioproc.eq.1).and.(isout.eq.1)) then

     do dir=1,3 
      print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
     enddo
     do dir=1,3 
      print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
     enddo

    endif

    call init2_FSI(part_id) 
   else if (ifirst.eq.0) then
    ! do nothing
   else
    print *,"ifirst invalid"
    stop
   endif 

   use_temp=0

     ! this routine initializes FSI(part_id)%solid_displ and
     ! FSI(part_id)%solid_speed
     !do_2nd_part=1(initinjector)
   call init3_FSI(part_id,ifirst,1,ioproc,isout) 

  else
   print *,"probtype invalid in initinjector"
   stop
  endif
  
return
end subroutine initinjector


! called when probtype.eq.701
subroutine initflapping(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  part_id,isout)
IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T :: j1,sdim,ifirst,isout
INTEGER_T :: inode,iface,ioproc
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1
REAL_T, dimension(3) :: xvalbefore
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
INTEGER_T :: dir,istep,istop
character(40) :: dwave

REAL_T, dimension(3) :: xxblob1,newxxblob1,xxblob2,newxxblob2
REAL_T :: radradblob1,radradblob2
INTEGER_T tempelem1
INTEGER_T tempelem2
INTEGER_T tempelem3
INTEGER_T tempelem4
INTEGER_T numquads,quad_counter
INTEGER_T localElem(4)

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id invalid"
   stop
  endif
  if (FSI(part_id)%partID.ne.part_id) then
   print *,"FSI(part_id)%partID.ne.part_id"
   stop
  endif

  if (probtype.ne.701) then
   print *,"probtype invalid initflapping"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
  endif

  timeB=curtime
  dtB=0.0
  TorquePos=0.0
  TorqueVel=0.0

  if ((ifirst.eq.0).or.(ifirst.eq.1)) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  if (ifirst.eq.1) then
   xxblob1(1)=0.0
   xxblob1(2)=0.0
   xxblob1(3)=0.0

   xxblob2(1)=0.0
   xxblob2(2)=0.0
   xxblob2(3)=0.0

   newxxblob2(1)=0.0
   newxxblob2(2)=0.0
   newxxblob2(3)=0.0
   radradblob2=1.0

   newxxblob1(1)=0.0
   newxxblob1(2)=0.0
   newxxblob1(3)=0.0
   radradblob1=1.0  

   do dir=1,3
     xxblob2(dir)=xxblob1(dir)
     newxxblob2(dir)=newxxblob1(dir)
   enddo
   radradblob2=radradblob1

   denpaddle=one
   dampingpaddle=zero

   if (FSI(part_id)%partID.eq.1 ) then

    dwave="foreWing.sci"
    newxxblob1(1)=0.0  ! 1st wing at a "default" location.
    print *,"init_flapping partID=1"

   else if (FSI(part_id)%partID.eq.2) then

    if ((axis_dir.eq.0).or.(axis_dir.eq.2)) then
     print *,"part id invalid"
     stop
    else if (axis_dir.eq.1) then
!     dwave="naca0012.sci"
     dwave="hindWing.sci"
     ! position 2nd wing in "default" location too.
     ! kinematics will shift wing in the x direction.
     newxxblob1(1)=0.0  
     print *,"init_flapping partID=2"
    else
     print *,"axis_dir invalid"
     stop
    endif

   endif

   numquads=0

   print *,"opening to check to see how many quads to convert ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
   READ(14,*) FSI(part_id)%NumNodes
   print *,"NumNodes(check) ",FSI(part_id)%NumNodes
   READ(14,*) FSI(part_id)%NumIntElems
   print *,"NumIntElems(check) ",FSI(part_id)%NumIntElems
   do inode=1,FSI(part_id)%NumNodes
    READ(14,*) xvalbefore(1),xvalbefore(3),xvalbefore(2)
   enddo
   do iface=1,FSI(part_id)%NumIntElems
    numquads=numquads+1
    READ(14,*) tempelem1,tempelem2,tempelem3,tempelem4
   enddo
   close(14)
   print *,"numquads=",numquads

   print *,"opening ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

   READ(14,*) FSI(part_id)%NumNodes
   print *,"NumNodes ",FSI(part_id)%NumNodes
   READ(14,*) FSI(part_id)%NumIntElems
   print *,"NumIntElems ",FSI(part_id)%NumIntElems
   FSI(part_id)%IntElemDim=3
   FSI(part_id)%NumIntElems=FSI(part_id)%NumIntElems+numquads

   if (ifirst.ne.1) then
    print *,"ifirst bust"
    stop
   endif

   call init_FSI(part_id,1)  ! allocate_intelem=1

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   do inode=1,FSI(part_id)%NumNodes
    READ(14,*) xval(1),xval(3),xval(2)
    do dir=1,3
     if ((minnodebefore(dir).gt.xval(dir)).or.(inode.eq.1)) then
      minnodebefore(dir)=xval(dir)
     endif
     if ((maxnodebefore(dir).lt.xval(dir)).or.(inode.eq.1)) then
      maxnodebefore(dir)=xval(dir)
     endif
    enddo

    do dir=1,3
     xval1(dir)=(xval(dir)-xxblob1(dir))/radradblob1 + newxxblob1(dir)
    enddo

    do dir=1,3
     if ((minnode(dir).gt.xval1(dir)).or.(inode.eq.1)) then
      minnode(dir)=xval1(dir)
     endif
     if ((maxnode(dir).lt.xval1(dir)).or.(inode.eq.1)) then
      maxnode(dir)=xval1(dir)
     endif
    enddo ! dir
    
    do dir=1,3
     FSI(part_id)%Node_old(dir,inode)=xval1(dir)
     FSI(part_id)%Node_new(dir,inode)=xval1(dir)
    enddo
      
   enddo  ! inode=1,NumNodes
   
   quad_counter=0 
   do iface=1,FSI(part_id)%NumIntElems-numquads
    READ(14,*) localElem(4),localElem(3),localElem(2),localElem(1)
    do j1=0,1
     if (j1.eq.1) then
      quad_counter=quad_counter+1
     endif
     FSI(part_id)%ElemData(1,iface+quad_counter)=3 ! number of nodes in element
     FSI(part_id)%ElemData(2,iface+quad_counter)=1 ! part number
     FSI(part_id)%ElemData(3,iface+quad_counter)=0 ! singly wetted
     if (j1.eq.0) then
      do dir=1,3
       FSI(part_id)%IntElem(dir,iface+quad_counter)=localElem(dir)
      enddo
     else if (j1.eq.1) then
       FSI(part_id)%IntElem(1,iface+quad_counter)=localElem(3)
       FSI(part_id)%IntElem(2,iface+quad_counter)=localElem(4)
       FSI(part_id)%IntElem(3,iface+quad_counter)=localElem(1)
     else
       print *,"j1 invalid"
       stop
     endif
    enddo  ! j1=0..1
   enddo  ! iface, looping faces
   close(14)
   if (quad_counter.ne.numquads) then
     print *,"quad_counter.ne.numquads"
     stop
   endif

   do dir=1,3 
    print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
   enddo
   do dir=1,3 
    print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
   enddo

   call init2_FSI(part_id)
  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif 

  use_temp=0

   ! do_2nd_part=0  
   ! even if 2 wings, we want do_2nd_part=0 because
   ! the velocity is prescribed in the kinematic routines later on.
  call init3_FSI(part_id,ifirst,0,ioproc,isout) 


return
end subroutine initflapping

subroutine init_from_cas(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  part_id,isout)
#if (STANDALONE==0)
use CAV3D_module
#endif

IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T :: sdim,ifirst,isout
INTEGER_T :: inode,iface,ioproc
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
INTEGER_T :: dir,istep,istop

REAL_T, dimension(3) :: xxblob1,newxxblob1
REAL_T :: radradblob1
INTEGER_T localElem(3)
REAL_T :: local_nodes(3,3)  ! dir,node num

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id invalid"
   stop
  endif
  if (FSI(part_id)%partID.ne.part_id) then
   print *,"FSI(part_id)%partID.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
   print *,"ioproc invalid"
   stop
  endif

  timeB=curtime
  dtB=0.0
  TorquePos=0.0
  TorqueVel=0.0

  if ((ifirst.eq.0).or.(ifirst.eq.1)) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  if (ifirst.eq.1) then
   xxblob1(1)=0.0
   xxblob1(2)=0.0
   xxblob1(3)=0.0

   newxxblob1(1)=0.0
   newxxblob1(2)=0.0
   newxxblob1(3)=0.0
   radradblob1=1.0  

   denpaddle=one
   dampingpaddle=zero

   if (probtype.eq.411) then
#if (STANDALONE==0)
    call OPEN_CAV3D_CASFILE(part_id,14)
#else
    print *,"probtype not recognized in init_from_cas1"
    stop
#endif
   else
    print *,"probtype not recognized in init_from_cas2"
    stop
   endif

   READ(14,*) FSI(part_id)%NumNodes,FSI(part_id)%NumIntElems
   print *,"NumNodes ",FSI(part_id)%NumNodes
   print *,"NumIntElems ",FSI(part_id)%NumIntElems
   FSI(part_id)%IntElemDim=3
   FSI(part_id)%NumIntElems=FSI(part_id)%NumIntElems

   if (ifirst.ne.1) then
    print *,"ifirst bust"
    stop
   endif

   call init_FSI(part_id,1)  ! allocate_intelem=1

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   do inode=1,FSI(part_id)%NumNodes
    READ(14,*) xval(1),xval(2),xval(3)
    do dir=1,3
     if ((minnodebefore(dir).gt.xval(dir)).or.(inode.eq.1)) then
      minnodebefore(dir)=xval(dir)
     endif
     if ((maxnodebefore(dir).lt.xval(dir)).or.(inode.eq.1)) then
      maxnodebefore(dir)=xval(dir)
     endif
    enddo ! dir=1..3

    do dir=1,3
     xval1(dir)=(xval(dir)-xxblob1(dir))/radradblob1 + newxxblob1(dir)
    enddo

    do dir=1,3
     if ((minnode(dir).gt.xval1(dir)).or.(inode.eq.1)) then
      minnode(dir)=xval1(dir)
     endif
     if ((maxnode(dir).lt.xval1(dir)).or.(inode.eq.1)) then
      maxnode(dir)=xval1(dir)
     endif
    enddo ! dir
    
    do dir=1,3
     FSI(part_id)%Node_old(dir,inode)=xval1(dir)
     FSI(part_id)%Node_new(dir,inode)=xval1(dir)
    enddo
      
   enddo  ! inode=1,NumNodes
   
   do iface=1,FSI(part_id)%NumIntElems
    READ(14,*) localElem(1),localElem(2),localElem(3)
    do inode=1,3
     if ((localElem(inode).lt.1).or. &
         (localElem(inode).gt.FSI(part_id)%NumNodes)) then
      print *,"localElem(inode) out of range"
      stop
     endif
    enddo !inode=1..3
    if ((localElem(1).eq.localElem(2)).or. &
        (localElem(1).eq.localElem(3)).or. &
        (localElem(2).eq.localElem(3))) then
     print *,"duplicate nodes for triangle"
     stop
    endif
    do dir=1,3
     do inode=1,3
      local_nodes(dir,inode)=FSI(part_id)%Node_new(dir,localElem(inode))
     enddo
    enddo

#if (STANDALONE==0)
    call CAV3D_ORDER_NODES(local_nodes,localElem)
#else
    print *,"CAV3D not recognized in standalone"
    stop
#endif

    FSI(part_id)%ElemData(1,iface)=3 ! number of nodes in element
    FSI(part_id)%ElemData(2,iface)=1 ! part number
    FSI(part_id)%ElemData(3,iface)=0 ! singly wetted
    do dir=1,3
     FSI(part_id)%IntElem(dir,iface)=localElem(dir)
    enddo
   enddo  ! iface, looping faces
   close(14)

   do dir=1,3 
    print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
   enddo
   do dir=1,3 
    print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
   enddo

   call init2_FSI(part_id)

  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif 

  
  use_temp=0

   ! do_2nd_part=0  
  call init3_FSI(part_id,ifirst,0,ioproc,isout) 


return
end subroutine init_from_cas

subroutine convert_2D_to_3D_nodes_FSI(part_id,inode)
IMPLICIT NONE

INTEGER_T, intent(in) :: part_id
INTEGER_T, intent(in) :: inode
INTEGER_T local_nodes,orig_nodes,dir

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%flag_2D_to_3D.eq.1) then
  ! do nothing
 else
  print *,"FSI(part_id)%flag_2D_to_3D invalid"
  stop
 endif
 local_nodes=FSI(part_id)%NumNodes
 orig_nodes=local_nodes/2
 if (orig_nodes*2.eq.local_nodes) then
  if ((inode.ge.1).and.(inode.le.orig_nodes)) then
   FSI(part_id)%Node_old(3,inode)=-one
   do dir=1,3
    FSI(part_id)%Node_old(dir,inode+orig_nodes)= &
            FSI(part_id)%Node_old(dir,inode)
   enddo
   FSI(part_id)%Node_old(3,inode+orig_nodes)=one
   do dir=1,3
    FSI(part_id)%Node_new(dir,inode)= &
      FSI(part_id)%Node_old(dir,inode)
    FSI(part_id)%Node_new(dir,inode+orig_nodes)= &
      FSI(part_id)%Node_old(dir,inode+orig_nodes)
   enddo
  else
   print *,"inode out of range"
   stop
  endif
 else
  print *,"local_nodes not divisible by 2"
  stop
 endif

end subroutine convert_2D_to_3D_nodes_FSI


subroutine convert_2D_to_3D_elements_FSI(part_id,iface)
IMPLICIT NONE

INTEGER_T, intent(in) :: part_id
INTEGER_T, intent(in) :: iface
INTEGER_T local_nodes,orig_nodes
INTEGER_T iflag
INTEGER_T local_elements,orig_elements


 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%flag_2D_to_3D.eq.1) then
  ! do nothing
 else
  print *,"FSI(part_id)%flag_2D_to_3D invalid"
  stop
 endif

 local_nodes=FSI(part_id)%NumNodes
 orig_nodes=local_nodes/2
 local_elements=FSI(part_id)%NumIntElems
 orig_elements=local_elements/2

 if ((orig_nodes*2.eq.local_nodes).and. &
     (orig_elements*2.eq.local_elements)) then
  if ((iface.ge.1).and.(iface.le.orig_elements)) then
   do iflag=1,3
    FSI(part_id)%ElemData(iflag,iface+orig_elements)= &
      FSI(part_id)%ElemData(iflag,iface)
   enddo
   FSI(part_id)%IntElem(1,iface+orig_elements)= &
     FSI(part_id)%IntElem(3,orig_elements)
   FSI(part_id)%IntElem(2,iface+orig_elements)= &
     FSI(part_id)%IntElem(1,orig_elements)+orig_nodes
   FSI(part_id)%IntElem(3,iface+orig_elements)= &
     FSI(part_id)%IntElem(1,orig_elements)
  else
   print *,"iface out of range"
   stop
  endif
 else
  print *,"local_nodes or local_elements not divisible by 2"
  stop
 endif

end subroutine convert_2D_to_3D_elements_FSI


subroutine init_gingerbread2D(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  part_id,isout)

IMPLICIT NONE

INTEGER_T, intent(in) :: part_id
INTEGER_T :: sdim,ifirst,isout
INTEGER_T :: inode,iface,ioproc
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
INTEGER_T :: dir,istep,istop

REAL_T, dimension(3) :: xxblob1,newxxblob1
REAL_T :: radradblob1
INTEGER_T localElem(3)
character(40) :: dwave
INTEGER_T :: orig_nodes,local_nodes
INTEGER_T :: orig_elements,local_elements

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id invalid"
   stop
  endif
  if (FSI(part_id)%partID.ne.part_id) then
   print *,"FSI(part_id)%partID.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
   print *,"ioproc invalid"
   stop
  endif

  FSI(part_id)%flag_2D_to_3D=1

  timeB=curtime
  dtB=0.0
  TorquePos=0.0
  TorqueVel=0.0

  if ((ifirst.eq.0).or.(ifirst.eq.1)) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  if (ifirst.eq.1) then
   xxblob1(1)=0.0
   xxblob1(2)=0.0
   xxblob1(3)=0.0

   newxxblob1(1)=2.0/30.0
   newxxblob1(2)=2.0/30.0
   newxxblob1(3)=0.0
   radradblob1=1.0/30.0 

   denpaddle=one
   dampingpaddle=zero

   if (probtype.eq.400) then ! gingerbread
    dwave="gingeroutline"
    print *,"opening ",dwave
    OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
   else
    print *,"probtype not recognized in init_gingerbread2D"
    stop
   endif

   READ(14,*) FSI(part_id)%NumNodes,FSI(part_id)%NumIntElems
   print *,"doubling number of nodes and elements 2D -> 3D"

   orig_nodes=FSI(part_id)%NumNodes
   FSI(part_id)%NumNodes=orig_nodes*2
   local_nodes=FSI(part_id)%NumNodes

   orig_elements=FSI(part_id)%NumIntElems
   FSI(part_id)%NumIntElems=orig_elements*2
   local_elements=FSI(part_id)%NumIntElems

   print *,"NumNodes ",local_nodes
   print *,"NumIntElems ",local_elements
   FSI(part_id)%IntElemDim=3

   if (ifirst.ne.1) then
    print *,"ifirst bust"
    stop
   endif

   call init_FSI(part_id,1)  ! allocate_intelem=1

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   do inode=1,orig_nodes

    READ(14,*) xval(1),xval(2),xval(3)

    do dir=1,3
     if ((minnodebefore(dir).gt.xval(dir)).or.(inode.eq.1)) then
      minnodebefore(dir)=xval(dir)
     endif
     if ((maxnodebefore(dir).lt.xval(dir)).or.(inode.eq.1)) then
      maxnodebefore(dir)=xval(dir)
     endif
    enddo ! dir=1..3

    do dir=1,3
     xval1(dir)=(xval(dir)-xxblob1(dir))/radradblob1 + newxxblob1(dir)
    enddo

    do dir=1,3
     if ((minnode(dir).gt.xval1(dir)).or.(inode.eq.1)) then
      minnode(dir)=xval1(dir)
     endif
     if ((maxnode(dir).lt.xval1(dir)).or.(inode.eq.1)) then
      maxnode(dir)=xval1(dir)
     endif
    enddo ! dir
    
    do dir=1,3
     FSI(part_id)%Node_old(dir,inode)=xval1(dir)
     FSI(part_id)%Node_new(dir,inode)=xval1(dir)
    enddo
    
    call convert_2D_to_3D_nodes_FSI(part_id,inode)

   enddo  ! inode=1,orig_nodes
   
   do iface=1,orig_elements

!   READ(14,*) localElem(1),localElem(2),localElem(3)
    READ(14,*) localElem(1),localElem(2)

!   do inode=1,3
    do inode=1,2
     if ((localElem(inode).lt.1).or. &
         (localElem(inode).gt.orig_nodes)) then
      print *,"localElem(inode) out of range"
      stop
     endif
    enddo !inode=1..2
!   enddo !inode=1..3

    localElem(3)=localElem(2)+orig_nodes

    if ((localElem(1).eq.localElem(2)).or. &
        (localElem(1).eq.localElem(3)).or. &
        (localElem(2).eq.localElem(3))) then
     print *,"duplicate nodes for triangle"
     stop
    endif

    FSI(part_id)%ElemData(1,iface)=3 ! number of nodes in element
    FSI(part_id)%ElemData(2,iface)=1 ! part number
    FSI(part_id)%ElemData(3,iface)=0 ! singly wetted
    do dir=1,3
     FSI(part_id)%IntElem(dir,iface)=localElem(dir)
    enddo

    call convert_2D_to_3D_elements_FSI(part_id,iface)

   enddo  ! iface, looping faces

   close(14)

   do dir=1,3 
    print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
   enddo
   do dir=1,3 
    print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
   enddo

   call init2_FSI(part_id)

  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif 

  
  use_temp=0

   ! do_2nd_part=0  
  call init3_FSI(part_id,ifirst,0,ioproc,isout) 


return
end subroutine init_gingerbread2D


subroutine init_helix(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  part_id,isout)

IMPLICIT NONE

INTEGER_T, intent(in) :: part_id
INTEGER_T :: sdim,ifirst,isout
INTEGER_T :: inode,iface
INTEGER_T :: inode_read
INTEGER_T :: iface_read
INTEGER_T :: ioproc
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
INTEGER_T :: dir,istep,istop

REAL_T, dimension(3) :: xxblob1,newxxblob1
REAL_T :: radradblob1
INTEGER_T localElem(3)
character(40) :: dwave
INTEGER_T :: orig_nodes,local_nodes
INTEGER_T :: orig_elements,local_elements

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id invalid"
   stop
  endif
  if (FSI(part_id)%partID.ne.part_id) then
   print *,"FSI(part_id)%partID.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
   print *,"ioproc invalid"
   stop
  endif

  FSI(part_id)%flag_2D_to_3D=0

  timeB=curtime
  dtB=0.0
  TorquePos=0.0
  TorqueVel=0.0

  if ((ifirst.eq.0).or.(ifirst.eq.1)) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  if (ifirst.eq.1) then
   xxblob1(1)=0.0
   xxblob1(2)=0.0
   xxblob1(3)=0.0

   newxxblob1(1)=0.0
   newxxblob1(2)=0.0
   newxxblob1(3)=0.0
   radradblob1=1.0

   denpaddle=one
   dampingpaddle=zero

   if (probtype.eq.401) then ! helix
    dwave="helix.dat"
    print *,"opening ",dwave
    OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
   else
    print *,"probtype not recognized in init_helix"
    stop
   endif

   READ(14,*) FSI(part_id)%NumNodes,FSI(part_id)%NumIntElems

   orig_nodes=FSI(part_id)%NumNodes
   local_nodes=FSI(part_id)%NumNodes

   orig_elements=FSI(part_id)%NumIntElems
   local_elements=FSI(part_id)%NumIntElems

   print *,"NumNodes ",local_nodes
   print *,"NumIntElems ",local_elements
   FSI(part_id)%IntElemDim=3

   if (ifirst.ne.1) then
    print *,"ifirst bust"
    stop
   endif

   call init_FSI(part_id,1)  ! allocate_intelem=1

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   do inode=1,orig_nodes

    READ(14,*) inode_read,xval(1),xval(2),xval(3)

    if (inode_read.eq.inode) then
     ! do nothing
    else
     print *,"inode_read invalid"
     stop
    endif

    do dir=1,3
     if ((minnodebefore(dir).gt.xval(dir)).or.(inode.eq.1)) then
      minnodebefore(dir)=xval(dir)
     endif
     if ((maxnodebefore(dir).lt.xval(dir)).or.(inode.eq.1)) then
      maxnodebefore(dir)=xval(dir)
     endif
    enddo ! dir=1..3

    do dir=1,3
     xval1(dir)=(xval(dir)-xxblob1(dir))/radradblob1 + newxxblob1(dir)
    enddo

    do dir=1,3
     if ((minnode(dir).gt.xval1(dir)).or.(inode.eq.1)) then
      minnode(dir)=xval1(dir)
     endif
     if ((maxnode(dir).lt.xval1(dir)).or.(inode.eq.1)) then
      maxnode(dir)=xval1(dir)
     endif
    enddo ! dir
    
    do dir=1,3
     FSI(part_id)%Node_old(dir,inode)=xval1(dir)
     FSI(part_id)%Node_new(dir,inode)=xval1(dir)
    enddo
    
   enddo  ! inode=1,orig_nodes
   
   do iface=1,orig_elements

    READ(14,*) iface_read,localElem(1),localElem(2),localElem(3)

    if (iface_read.eq.iface) then
     ! do nothing
    else
     print *,"iface_read invalid"
     stop
    endif

    do inode=1,3
     if ((localElem(inode).lt.1).or. &
         (localElem(inode).gt.orig_nodes)) then
      print *,"localElem(inode) out of range"
      stop
     endif
    enddo !inode=1..3

    if ((localElem(1).eq.localElem(2)).or. &
        (localElem(1).eq.localElem(3)).or. &
        (localElem(2).eq.localElem(3))) then
     print *,"duplicate nodes for triangle"
     stop
    endif

    FSI(part_id)%ElemData(1,iface)=3 ! number of nodes in element
    FSI(part_id)%ElemData(2,iface)=1 ! part number
    FSI(part_id)%ElemData(3,iface)=0 ! singly wetted
    do dir=1,3
     FSI(part_id)%IntElem(dir,iface)=localElem(dir)
    enddo

   enddo  ! iface, looping faces

   close(14)

   do dir=1,3 
    print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
   enddo
   do dir=1,3 
    print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
   enddo

   call init2_FSI(part_id)

  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif 

  
  use_temp=0

   ! do_2nd_part=0  
  call init3_FSI(part_id,ifirst,0,ioproc,isout) 


return
end subroutine init_helix

subroutine initchannel(curtime,dt,ifirst,sdim,istop,istep)
IMPLICIT NONE

INTEGER_T :: j1,sdim,ifirst
INTEGER_T :: inode,iface
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1
REAL_T, dimension(3) :: maxnodebefore,minnodebefore,xvalbefore
INTEGER_T :: dir,istep,istop
INTEGER_T :: numquads,quad_counter
INTEGER_T :: tempelem1,tempelem2,tempelem3,tempelem4
INTEGER_T, dimension(4) :: localElem
character(40) :: dwave

REAL_T, dimension(3) :: xxblob1,newxxblob1
REAL_T :: radradblob1
INTEGER_T :: local_part_id

  local_part_id=1

  timeB=curtime
  dtB=0.0
  TorquePos=0.0
  TorqueVel=0.0

  if ((ifirst.eq.0).or.(ifirst.eq.1)) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  if (probtype.ne.5700) then
   print *,"probtype invalid initchannel"
   stop
  endif
  if (ifirst.eq.1) then
   if (axis_dir.eq.0) then
    xxblob1(1)=0.0
    xxblob1(2)=0.0
    xxblob1(3)=0.0
    newxxblob1(1)=0.0
    newxxblob1(2)=0.0
    newxxblob1(3)=0.0
    radradblob1=100.0
    dwave="trap2_device.sci"
   else if (axis_dir.eq.3) then
    xxblob1(1)=0.0
    xxblob1(2)=0.0
    xxblob1(3)=0.0
    newxxblob1(1)=0.0
    newxxblob1(2)=0.0
    newxxblob1(3)=0.0
    radradblob1=100.0
    dwave="trap2_device_simple.sci"
   else if (axis_dir.eq.4) then
    xxblob1(1)=0.0
    xxblob1(2)=0.0
    xxblob1(3)=0.0
    newxxblob1(1)=0.0
    newxxblob1(2)=0.0
    newxxblob1(3)=0.0
    radradblob1=100.0
    dwave="channel2013.sci"
   else if (axis_dir.eq.1) then
    xxblob1(1)=0.0
    xxblob1(2)=0.0
    xxblob1(3)=0.0
    newxxblob1(1)=0.0
    newxxblob1(2)=0.0
    newxxblob1(3)=0.0
    radradblob1=1.0
    dwave="square_T_junction.sci"
   else if (axis_dir.eq.2) then
    xxblob1(1)=0.0
    xxblob1(2)=0.0
    xxblob1(3)=0.0
    newxxblob1(1)=0.0
    newxxblob1(2)=0.0
    newxxblob1(3)=0.0
     ! modified from 1/100 on November 8, 2012
    radradblob1=1.0/10000.0  ! convert to cm  
! elements are ordered clockwise looking from the outside of the device
    dwave="squeeze_T_junction.sci"
   else
    print *,"axis_dir invalid"
    stop
   endif
  
   denpaddle=1.0
   dampingpaddle=0.0

   numquads=0

   print *,"opening to check to see how many quads to convert ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
   READ(14,*) FSI(1)%NumNodes
   print *,"NumNodes(check) ",FSI(1)%NumNodes
   READ(14,*) FSI(1)%NumIntElems
   print *,"NumIntElems(check) ",FSI(1)%NumIntElems
   do inode=1,FSI(1)%NumNodes
    READ(14,*) xvalbefore(1),xvalbefore(2),xvalbefore(3)
   enddo
   do iface=1,FSI(1)%NumIntElems
    numquads=numquads+1
    READ(14,*) tempelem1,tempelem2,tempelem3,tempelem4
   enddo
   close(14)
   print *,"numquads=",numquads

   print *,"opening ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

   READ(14,*) FSI(1)%NumNodes
   print *,"NumNodes ",FSI(1)%NumNodes
   READ(14,*) FSI(1)%NumIntElems
   print *,"NumIntElems ",FSI(1)%NumIntElems
   FSI(1)%IntElemDim=3
   FSI(1)%NumIntElems=FSI(1)%NumIntElems+numquads

   if (ifirst.eq.1) then
    call init_FSI(local_part_id,1)
   else
    print *,"ifirst bust"
    stop
   endif 

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   do inode=1,FSI(1)%NumNodes
    READ(14,*) xval(1),xval(2),xval(3)
    do dir=1,3
     if ((minnodebefore(dir).gt.xval(dir)).or.(inode.eq.1)) then
      minnodebefore(dir)=xval(dir)
     endif
     if ((maxnodebefore(dir).lt.xval(dir)).or.(inode.eq.1)) then
      maxnodebefore(dir)=xval(dir)
     endif
    enddo

    do dir=1,3
     xval1(dir)=(xval(dir)-xxblob1(dir))/radradblob1 + newxxblob1(dir)
    enddo
    if (axis_dir.eq.2) then
     if ((xblob2.le.0.0).or.(xblob2.gt.10.0)) then
      print *,"xblob2 invalid xblob2=",xblob2
      stop
     endif
     if (xval1(1).ge.10.0) then
      xval1(1)=xval1(1)-10.0+xblob2
     endif
    endif

    do dir=1,3
     if ((minnode(dir).gt.xval1(dir)).or.(inode.eq.1)) then
      minnode(dir)=xval1(dir)
     endif
     if ((maxnode(dir).lt.xval1(dir)).or.(inode.eq.1)) then
      maxnode(dir)=xval1(dir)
     endif
    enddo ! dir
    
    do dir=1,3
     FSI(1)%Node_old(dir,inode)=xval1(dir)
     FSI(1)%Node_new(dir,inode)=xval1(dir)
    enddo
      
   enddo  ! inode=1,NumNodes
   
   quad_counter=0 
   do iface=1,FSI(1)%NumIntElems-numquads
    READ(14,*) localElem(1),localElem(2),localElem(3),localElem(4)
    do j1=0,1
     if (j1.eq.1) then
      quad_counter=quad_counter+1
     endif
     FSI(1)%ElemData(1,iface+quad_counter)=3   ! number of nodes in element
     FSI(1)%ElemData(2,iface+quad_counter)=1   ! part number
     FSI(1)%ElemData(3,iface+quad_counter)=0   ! singly wetted
     if (j1.eq.0) then
      do dir=1,3
       FSI(1)%IntElem(dir,iface+quad_counter)=localElem(dir)
      enddo
     else if (j1.eq.1) then
       FSI(1)%IntElem(1,iface+quad_counter)=localElem(3)
       FSI(1)%IntElem(2,iface+quad_counter)=localElem(4)
       FSI(1)%IntElem(3,iface+quad_counter)=localElem(1)
     else
       print *,"j1 invalid"
       stop
     endif
    enddo  ! j1=0..1
   enddo  ! iface, looping faces
   close(14)
   if (quad_counter.ne.numquads) then
     print *,"quad_counter.ne.numquads"
     stop
   endif

   do dir=1,3 
    print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
   enddo
   do dir=1,3 
    print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
   enddo

   call init2_FSI(local_part_id)

  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop 
  endif 

  use_temp=0
  ! do_2nd_part=1  isout=1 (initchannel)
  call init3_FSI(local_part_id,ifirst,1,1,1)  

return
end subroutine initchannel




subroutine geominit(curtime,dt,ifirst,sdim,istop,istep)
IMPLICIT NONE

INTEGER_T :: i,j,k,i1,j1,m1,tempelem,itimecount
INTEGER_T :: ifirst ! parameter
INTEGER_T :: iread  ! local variable
INTEGER_T :: sdim
INTEGER_T :: numquads,quad_counter
INTEGER_T, dimension(4) :: localElem
REAL_T :: curtime,dt
INTEGER_T, dimension(0:1000) :: itimearr
REAL_T, dimension(0:1000) :: timearr
REAL_T, dimension(3) :: maxnode,minnode,xval,xtemp
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
REAL_T, dimension(3) :: xvalbefore
REAL_T, dimension(3) :: invertfactor
REAL_T :: tper,tcrit,theta,radradblobpool
INTEGER_T :: iper,icrit,ewave,fwave,dir,istep,istop
INTEGER_T :: gwave
INTEGER_T :: shift_from_zero_node
INTEGER_T :: shift_from_zero_face
INTEGER_T :: override_IntElemDim
character(100) :: dwave,poolname
REAL_T :: ofs
REAL_T :: plungerfreq,plungeramp
INTEGER_T :: local_part_id

  local_part_id=1

  timeB=curtime
  dtB=dt

  if ((ifirst.eq.0).or.(ifirst.eq.1)) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  if (probtype.eq.58) then
   iread=1
  else if (probtype.eq.55) then
   iread=ifirst
  else
 
   OPEN(unit=25,file="hcorresp.txt",access='sequential', &
     form="formatted",status='old')
   READ(25,115) itimecount
   print *,"itimecount= ",itimecount

   if (itimecount.gt.1000) then
    print *,"itimecoount too big"
    stop
   endif

   do i1=1,itimecount
    READ(25,*) itimearr(i1),timearr(i1)
    print *,"index,time ",itimearr(i1),timearr(i1)
   enddo
   close(25)

   timearr(0)=0.0
   itimearr(0)=itimearr(itimecount)

   tper=curtime/timearr(itimecount)
   iper=INT(tper)
   tcrit=curtime-iper*timearr(itimecount)
   if (tcrit.gt.timearr(itimecount)+1.0E-8) then
    print *,"tcrit invalid"
    stop
   endif

   icrit=0
   do i1=1,itimecount
    if ((tcrit.le.timearr(i1)+1.0E-8).and.(icrit.eq.0)) then
     icrit=i1
    endif
   enddo
   if (icrit.eq.0) then
    print *,"icrit invalid"
    stop
   endif

   if (ifirst.eq.1) then
    tstart=-1.0
    tfinish=-1.0
   endif

   if ((tcrit.ge.tstart).and.(tcrit.le.tfinish)) then
    iread=0
   else
    iread=1
    tstart=timearr(icrit-1)
    tfinish=timearr(icrit)
   endif
  endif

  if (iread.eq.1) then

   FSI(1)%NumNodesPool=0
   FSI(1)%NumIntElemsPool=0
   FSI(1)%IntElemDimPool=0

   do i1=1,2 
    if ((probtype.eq.55).or. &
        (probtype.eq.58)) then
     j=0
    else
     if (i1.eq.1) then
      j=itimearr(icrit-1)
     else
      j=itimearr(icrit)
     endif
    endif
    invertfactor(1)=1.0
    invertfactor(2)=1.0
    invertfactor(3)=1.0

    if (probtype.eq.55) then
     xxblob(1)=0.0
     xxblob(2)=0.0
     xxblob(3)=0.0
     newxxblob(1)=0.0
     newxxblob(2)=0.0
     newxxblob(3)=0.055
     radradblob=600.0
     dwave="Sphere.txt"
    else if (probtype.eq.58) then
     plungerfreq=14.0
     plungeramp=0.015  ! 0.015
     xxblob(1)=0.0
     xxblob(2)=0.0
     xxblob(3)=0.0
     newxxblob(1)=0.0
     newxxblob(2)=0.0
     newxxblob(3)=0.0
     radradblob=1.0
     newxxblob(1)=plungeramp*sin(2.0*3.14159*plungerfreq*curtime)
     dwave="square.txt"
    else if (probtype.eq.52) then
     xxblob(1)=260.0
     xxblob(2)=19.0
     xxblob(3)=-191.0
     newxxblob(3)=0.5
     newxxblob(1)=0.25
     newxxblob(2)=0.25
     radradblob=440.0

     if (j.lt.10) then
      dwave = "h/hand000"//char(48+j)//".txt"
     else if (j.lt.100) then
       ewave = j/10
       dwave = "h/hand00"//char(48+ewave)//char(48+j-ewave*10)//".txt"
     else
       ewave = j/100
       fwave = (j-ewave*100)/10
       dwave = "h/hand0"//char(48+ewave)//char(48+fwave)// &
              char(48+j-fwave*10-ewave*100)//".txt"
     endif
    else if (probtype.eq.57) then
!       maya generated heart
     if (axis_dir.eq.0) then 
      xxblob(1)=-375.0
      xxblob(2)=890.0
      xxblob(3)=-70.0
      radradblob=550.0
      print *,"maya heart"
     else if (axis_dir.eq.1) then
!       vue generated heart
      xxblob(1)=-15.0
      xxblob(2)=-50.0
      xxblob(3)=0.0
      radradblob=550.0
      print *,"vue heart"
     else if (axis_dir.eq.2) then

!       SCR heart
      xxblob(1)=98.0
      xxblob(2)=70.0
      xxblob(3)=82.0
      radradblob=200.0
      print *,"SCR heart"
     else 
      print *,"axis_dir invalid"
      stop
     endif

     newxxblob(1)=0.5
     newxxblob(2)=0.5
     newxxblob(3)=0.5

     if (j.lt.10) then
      dwave = "h/heartleft000"//char(48+j)//".txt"
     else if (j.lt.100) then
       ewave = j/10
       dwave = "h/heartleft00"//char(48+ewave)//char(48+j-ewave*10)//".txt"
     else
       ewave = j/100
       fwave = (j-ewave*100)/10
       dwave = "h/heartleft0"//char(48+ewave)//char(48+fwave)// &
              char(48+j-fwave*10-ewave*100)//".txt"
     endif
    else if (probtype.eq.562) then

! open mouth
     if (axis_dir.eq.3) then

! xval(dir)=invertfactor(dir)* &
!       (xvalbefore(dir)-xxblob(dir))/radradblob + newxxblob(dir)
!
! xval(1)=xtemp(2) xval(2)=xtemp(1) xval(3)=xtemp(3)

      xxblob(1)=0.0 ! will become y
      xxblob(2)=-1.0 ! will become -x
      invertfactor(2)=-1.0
      xxblob(3)=0.0  ! will become z
      newxxblob(1)=0.0
      newxxblob(2)=1.0
      newxxblob(3)=0.5
      radradblob=7.5     ! 12.5/1.67=7.5

! was 10 for whale file dated December, 2007  (whalenormal)
! xval(1)=xtemp(3) xval(2)=xtemp(1) xval(3)=xtemp(2)
     else if (axis_dir.eq.2) then
      xxblob(1)=0.0
      xxblob(2)=0.0
      xxblob(3)=0.0
      newxxblob(1)=0.0
      newxxblob(2)=0.5
      newxxblob(3)=0.9   ! was 1.0
      radradblob=10.87   ! was 12.5, now 10.87
! whalepregnant
     else if (axis_dir.eq.6) then
      xxblob(1)=0.0 ! will become y
      xxblob(3)=-1.0 ! will become -x
      invertfactor(3)=-1.0
      invertfactor(2)=-1.0
      xxblob(2)=0.0  ! will become z
      newxxblob(1)=0.0
      newxxblob(3)=1.0
      newxxblob(2)=0.5
      radradblob=10.87   ! was 12.5, now 10.87
! whaletailup 
     else if (axis_dir.eq.0) then
      xxblob(1)=0.0 ! will become y
      xxblob(2)=-1.0 ! will become -x
      invertfactor(2)=-1.0
      xxblob(3)=0.0  ! will become z
      newxxblob(1)=0.0
      newxxblob(2)=1.0
      newxxblob(3)=0.5
      radradblob=10.87   ! was 12.5, now 10.87
! whaletaildown
     else if (axis_dir.eq.5) then
      xxblob(1)=0.0 ! will become y
      xxblob(2)=-1.0 ! will become -x
      invertfactor(2)=-1.0
      xxblob(3)=0.0  ! will become z
      newxxblob(1)=0.0
      newxxblob(2)=1.0
      newxxblob(3)=0.5
      radradblob=10.87   ! was 12.5, now 10.87
     else
      print *,"axis_dir bad in geominit (for whale) probtype,axis_dir ", &
       probtype,axis_dir
      stop
     endif

     if (axis_dir.eq.0) then 
      if (1.eq.0) then
       if (j.lt.10) then
        dwave = "h/whale000"//char(48+j)//".txt"
       else if (j.lt.100) then
         ewave = j/10
         dwave = "h/whale00"//char(48+ewave)//char(48+j-ewave*10)//".txt"
       else
         ewave = j/100
         fwave = (j-ewave*100)/10
         dwave = "h/whale0"//char(48+ewave)//char(48+fwave)// &
               char(48+j-fwave*10-ewave*100)//".txt"
       endif
      else
       dwave = "h/whaletailup.txt"
       print *,"reading h/whaletailup.txt"
      endif 
     else if (axis_dir.eq.5) then
      dwave = "h/whaletaildown.txt"
      print *,"reading h/whaletaildown.txt"
     else if (axis_dir.eq.2) then
      dwave = "h/whalenormal.txt"
      print *,"reading h/whalenormal.txt"
     else if (axis_dir.eq.3) then
      dwave = "h/basic_whale_open_mouth08.txt"
      print *,"reading h/whale_open_mouth08.txt"
     else if (axis_dir.eq.6) then
      dwave = "h/whalepregnant.txt"
      print *,"reading h/whalepregnant.txt"
     else
      print *,"axis_dir bad in geominit (for whale) probtype axis_dir ", &
        probtype,axis_dir
      stop
     endif
! dog
! xval=invertfactor*(xvalbefore-xxblob)/radradblob + newxxblob
! x=z y=x z=y
! in the dog file: 
! length(z): -1.8 .. 2.7
! width(x) : .7 .. -.7
! height(y): 0 .. 3.6
! after transformation:
! length(x): 1..46
! width(y) : 17 .. 3
! height(z): 0..36
    else if (probtype.eq.5600) then
     xxblob(1)=0.0
     xxblob(2)=0.0
     xxblob(3)=0.0
     newxxblob(1)=10.0
     newxxblob(2)=0.0
     newxxblob(3)=19.0
     radradblob=0.1
     FSI(1)%NumNodesPool=0
     FSI(1)%NumIntElemsPool=0
     FSI(1)%IntElemDimPool=0

     if (j.lt.10) then
      dwave = "h/swimmer000"//char(48+j)//".txt"
     else if (j.lt.100) then
       ewave = j/10
       dwave = "h/swimmer00"//char(48+ewave)//char(48+j-ewave*10)//".txt"
     else if (j.lt.1000) then
       ewave = j/100
       fwave = (j-ewave*100)/10
       dwave = "h/swimmer0"//char(48+ewave)//char(48+fwave)// &
              char(48+j-fwave*10-ewave*100)//".txt"
     else
       ewave = j/1000
       fwave = (j-ewave*1000)/100
       gwave = (j-fwave*100-ewave*1000)/10
       dwave = "h/swimmer"//char(48+ewave)//char(48+fwave)// &
        char(48+gwave)//char(48+j-gwave*10-fwave*100-ewave*1000)//".txt"
     endif

    else if ((probtype.eq.56).or.(probtype.eq.561)) then
     xxblob(1)=0.2
     xxblob(2)=0.25
     xxblob(3)=0.5
     if (1.eq.0) then
      newxxblob(3)=0.65
      newxxblob(1)=0.375
      newxxblob(2)=0.375
      radradblob=40.0
      radradblobpool=40.0
     else
      newxxblob(3)=26
      newxxblob(1)=15   
      newxxblob(2)=15
      radradblob=1.0
      radradblobpool=1.0
     endif


! Viorel's test problem with a box.
     if (1.eq.0) then   
      xxblob(1)=0.0
      xxblob(2)=0.0
      xxblob(3)=0.0
      radradblob=3.0
     endif

! Viorel's ball and scythe problems 
     if (1.eq.0) then   
      xxblob(1)=2.5
      xxblob(2)=86.0
      xxblob(3)=-1.0
      newxxblob(3)=30.0
      newxxblob(1)=30.0
      newxxblob(2)=90.0
      radradblob=1.0
     endif
!Viorel's diving problem     
     if (1.eq.1) then   
      xxblob(1)=0.0
      xxblob(2)=36.0
      xxblob(3)=-0.5
      newxxblob(3)=30.0
      newxxblob(1)=30.0
      newxxblob(2)=18.0
      radradblob=1.0
     endif
!Viorel's paddle problem     
     if (1.eq.1) then   
      xxblob(1)=0.0
      xxblob(2)=0.0
      xxblob(3)=0.0
      newxxblob(3)=30.0
      newxxblob(1)=60.0
      newxxblob(2)=40.0
      radradblob=1.0
     endif

     if (probtype.eq.561) then
      poolname="h/Pool.txt"
      print *,"opening ",poolname
      OPEN(unit=141,file=poolname,access='sequential', &
        form="formatted",status='old')
      READ(141,*) FSI(1)%NumNodesPool
      print *,"NumNodesPool ",FSI(1)%NumNodesPool
      READ(141,*) FSI(1)%NumIntElemsPool
      print *,"NumIntElemsPool ",FSI(1)%NumIntElemsPool
      READ(141,*) FSI(1)%IntElemDimPool
      print *,"IntElemDimPool ",FSI(1)%IntElemDimPool
     else
      FSI(1)%NumNodesPool=0
      FSI(1)%NumIntElemsPool=0
      FSI(1)%IntElemDimPool=0
     endif
 
     if (j.lt.10) then
      dwave = "h/swimmer000"//char(48+j)//".txt"
     else if (j.lt.100) then
       ewave = j/10
       dwave = "h/swimmer00"//char(48+ewave)//char(48+j-ewave*10)//".txt"
     else if (j.lt.1000) then
       ewave = j/100
       fwave = (j-ewave*100)/10
       dwave = "h/swimmer0"//char(48+ewave)//char(48+fwave)// &
              char(48+j-fwave*10-ewave*100)//".txt"
     else
       ewave = j/1000
       fwave = (j-ewave*1000)/100
       gwave = (j-fwave*100-ewave*1000)/10
       dwave = "h/swimmer"//char(48+ewave)//char(48+fwave)// &
        char(48+gwave)//char(48+j-gwave*10-fwave*100-ewave*1000)//".txt"
     endif
    else
     print *,"probtype invalid in geominit"
     stop
    endif

    numquads=0
    if (FSI(1)%IntElemDimPool.eq.0) then
     print *,"opening to check to see how many quads to convert ",dwave
     OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
     READ(14,91) FSI(1)%NumNodes
     print *,"NumNodes(check) ",FSI(1)%NumNodes
     READ(14,91) FSI(1)%NumIntElems
     print *,"NumIntElems(check) ",FSI(1)%NumIntElems
     READ(14,91) FSI(1)%IntElemDim
     print *,"IntElemDim(check) ",FSI(1)%IntElemDim
     do i=1,FSI(1)%NumNodes
      READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
     enddo
     READ(14,*) j
     do i=1,FSI(1)%NumIntElems
      READ(14,*) j
      READ(14,*) k
      if (k.eq.FSI(1)%NumIntElems) then
       k=3
      endif
      if (k.eq.4) then
       numquads=numquads+1
      else if (k.ne.3) then
       print *,"elements must be triangles or quadrilaterals k=",k
       print *,"i,NumIntElems ",i,FSI(1)%NumIntElems
       stop
      endif
      do j1=1,k
       READ(14,*) tempelem
      enddo
     enddo

     close(14)
    endif
    print *,"numquads=",numquads

    print *,"opening ",dwave
    OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

    READ(14,91) FSI(1)%NumNodes
    print *,"NumNodes ",FSI(1)%NumNodes
    READ(14,91) FSI(1)%NumIntElems
    print *,"NumIntElems ",FSI(1)%NumIntElems
    READ(14,91) FSI(1)%IntElemDim
    print *,"IntElemDim ",FSI(1)%IntElemDim
    print *,"overriding IntElemDim to be 3"
    FSI(1)%IntElemDim=3
    FSI(1)%NumNodes=FSI(1)%NumNodes+FSI(1)%NumNodesPool
    FSI(1)%NumIntElems=FSI(1)%NumIntElems+FSI(1)%NumIntElemsPool

    FSI(1)%NumIntElems=FSI(1)%NumIntElems+numquads

    if (FSI(1)%IntElemDimPool.ne.0) then
     if (FSI(1)%IntElemDim.ne.FSI(1)%IntElemDimPool) then
      print *,"IntElemDim inconsistent"
      stop
     endif
    endif

    if ((ifirst.eq.1).and.(i1.eq.1)) then
      ! allocates NodeMass + other variables.
      ! it is assumed that the number of nodes and elements does not
      ! change from frame to frame.
     call init_FSI(local_part_id,1)
    endif

    do dir=1,3
     maxnode(dir)=0.0
     minnode(dir)=0.0
     maxnodebefore(dir)=-1.0e+10
     minnodebefore(dir)=1.0e+10
    enddo

    if (probtype.eq.561) then
     if (numquads.ne.0) then
      print *,"only triangles to be used here for now"
      stop
     endif
     do i=1,FSI(1)%NumNodesPool
      READ(141,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)

      do dir=1,3
       xval(dir)=invertfactor(dir)* &
        (xvalbefore(dir)-xxblob(dir))/radradblobpool + newxxblob(dir)
      enddo

      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
    
      if (i1.eq.1) then
       do dir=1,3
        FSI(1)%Node_old(dir,i)=xval(dir)
       enddo
      else if (i1.eq.2) then
       do dir=1,3
        FSI(1)%Node_new(dir,i)=xval(dir)
       enddo
      else
       print *,"i1 invalid"
       stop
      endif
      if (i.ne.j) then
       print *,"vertex mismatch reading pool file"
       stop
      endif
      
     enddo
     READ(141,*) j
     if (j.ne.FSI(1)%NumIntElemsPool) then
      print *,"face mismatch reading pool file"
      stop
     endif

     do i=1,FSI(1)%NumIntElemsPool
      READ(141,*) j
      if (i.ne.j) then
       print *,"face mismatch reading pool file"
       stop
      endif
      READ(141,*) k
      if (k.gt.FSI(1)%IntElemDim) then
       print *,"too many vertices k,IntElemDim ",k,FSI(1)%IntElemDim
       stop
      endif
      FSI(1)%ElemData(1,i)=k ! number of nodes in element
      FSI(1)%ElemData(2,i)=1 ! part number
      FSI(1)%ElemData(3,i)=2 ! singly wetted, but do not call "fill" for these
      do j1=1,k
       READ(141,*) FSI(1)%IntElem(j1,i) 
      enddo
     enddo  ! i, looping faces
     close(141)
    endif ! probtype.eq.561 


    shift_from_zero_node=0
    shift_from_zero_face=0
    override_IntElemDim=0
    do i=FSI(1)%NumNodesPool+1,FSI(1)%NumNodes
     READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
     if ((j.eq.0).and.(shift_from_zero_node.eq.0)) then
      shift_from_zero_node=1
      print *,"nodes in file start at 0; shifting to start at 1"
     endif
     if (shift_from_zero_node.eq.1) then
      j=j+1
     endif

     do dir=1,3
      xval(dir)=invertfactor(dir)* &
       (xvalbefore(dir)-xxblob(dir))/radradblob + newxxblob(dir)
     enddo

     if (probtype.eq.55) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(1)
      xval(2)=xtemp(2)
      xval(3)=xtemp(3)
     else if (probtype.eq.58) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(1)
      xval(2)=xtemp(2)
      xval(3)=xtemp(3)
     else if (probtype.eq.52) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
      if (i1.eq.1) then
       call scihandoffset(ofs,tstart)
      else if (i1.eq.2) then
       call scihandoffset(ofs,tfinish)
      else
       print *,"i1 invalid"
       stop
      endif
      xval(sdim)=xval(sdim)+ofs
     else if (probtype.eq.562) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
! open mouth whale
      if (axis_dir.eq.3) then

       xval(1)=xtemp(2)
       xval(2)=xtemp(1)
       xval(3)=xtemp(3)

! whale file used at beginning of summer 2007 (whalenormal)
      else if (axis_dir.eq.2) then
       xval(1)=xtemp(3)
       xval(2)=xtemp(1)
       xval(3)=xtemp(2)
! whaletaildown
      else if (axis_dir.eq.5) then
       xval(1)=xtemp(2)
       xval(2)=xtemp(1)
       xval(3)=xtemp(3)
! whaletailup 
      else if (axis_dir.eq.0) then
       xval(1)=xtemp(2)
       xval(2)=xtemp(1)
       xval(3)=xtemp(3)
! whale pregnant
      else if (axis_dir.eq.6) then
       xval(1)=xtemp(3)
       xval(2)=xtemp(1)
       xval(3)=xtemp(2)
      else
       print *,"bad axis_dir geominit (for whale) probtype,axis_dir ", &
         probtype,axis_dir
       stop
      endif
! dog
     else if (probtype.eq.5600) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
     else if ((probtype.eq.56).or.(probtype.eq.561)) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
! heart
     else if (probtype.eq.57) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
     else
      print *,"probtype invalid in geominit"
      stop
     endif
     do dir=1,3
      if ((minnode(dir).gt.xval(dir)).or.(i.eq.1)) then
       minnode(dir)=xval(dir)
      endif
      if ((maxnode(dir).lt.xval(dir)).or.(i.eq.1)) then
       maxnode(dir)=xval(dir)
      endif
      if ((minnodebefore(dir).gt.xvalbefore(dir)).or.(i.eq.1)) then
       minnodebefore(dir)=xvalbefore(dir)
      endif
      if ((maxnodebefore(dir).lt.xvalbefore(dir)).or.(i.eq.1)) then
       maxnodebefore(dir)=xvalbefore(dir)
      endif
     enddo
    
     if (i1.eq.1) then
      do dir=1,3
       FSI(1)%Node_old(dir,i)=xval(dir)
      enddo
     else if (i1.eq.2) then
      do dir=1,3
       FSI(1)%Node_new(dir,i)=xval(dir)
      enddo
     else
      print *,"i1 invalid"
      stop
     endif
     if (i-FSI(1)%NumNodesPool.ne.j) then
      print *,"vertex mismatch"
      print *,"NumNodesPool=",FSI(1)%NumNodesPool
      print *,"expected node index ",i-FSI(1)%NumNodesPool
      print *,"node index read from file ",j
      stop
     endif
      
    enddo
    READ(14,*) j
    if (j.ne.FSI(1)%NumIntElems-FSI(1)%NumIntElemsPool-numquads) then
     print *,"face mismatch"
     stop
    endif

    quad_counter=0
    do i=1+FSI(1)%NumIntElemsPool,FSI(1)%NumIntElems-numquads
     READ(14,*) j
     if ((j.eq.0).and.(shift_from_zero_face.eq.0)) then
      shift_from_zero_face=1
      print *,"faces in file start at 0; shifting to start at 1"
     endif
     if (shift_from_zero_face.eq.1) then
      j=j+1
     endif

     if (i-FSI(1)%NumIntElemsPool.ne.j) then
      print *,"face mismatch"
      print *,"NumIntElemsPool=",FSI(1)%NumIntElemsPool
      print *,"expected face index ",i-FSI(1)%NumIntElemsPool
      print *,"face index read from file ",j
      stop
     endif
     READ(14,*) k
     if (k.eq.FSI(1)%NumIntElems) then
      k=3
     endif
     if ((k.ne.3).and.(k.ne.4)) then
      print *,"k invalid k=",k
      stop
     endif

     do j1=1,k
      READ(14,*) localElem(j1)
     enddo

     do j1=0,k-3
      if (j1.eq.1) then
       quad_counter=quad_counter+1
      endif
      FSI(1)%ElemData(1,i+quad_counter)=3   ! number of nodes in element
      FSI(1)%ElemData(2,i+quad_counter)=1   ! part number
      FSI(1)%ElemData(3,i+quad_counter)=0   ! singly wetted
      if (probtype.eq.57) then  ! 57 heart    56 swimmer
       FSI(1)%ElemData(3,i+quad_counter)=1   ! doubly wetted
      endif
      if ((probtype.eq.56).and.(1.eq.1)) then  ! BOXSWIMMER
       FSI(1)%ElemData(3,i+quad_counter)=1   ! doubly wetted
      endif
      if (probtype.eq.562) then  ! 562 whale
       FSI(1)%ElemData(3,i+quad_counter)=0   ! 0=singly wetted 1=doubly wetted
       if (axis_dir.eq.6) then
        FSI(1)%ElemData(3,i+quad_counter)=0 
       endif
      endif
      if (probtype.eq.5600) then ! dog
       FSI(1)%ElemData(3,i+quad_counter)=0   ! 0=singly wetted 1=doubly wetted
      endif
      if (j1.eq.0) then
       do dir=1,3
        FSI(1)%IntElem(dir,i+quad_counter)=localElem(dir)
       enddo
      else if (j1.eq.1) then
       FSI(1)%IntElem(1,i+quad_counter)=localElem(3)
       FSI(1)%IntElem(2,i+quad_counter)=localElem(4)
       FSI(1)%IntElem(3,i+quad_counter)=localElem(1)
      else
       print *,"j1 invalid"
       stop
      endif
      do m1=1,3
       FSI(1)%IntElem(m1,i+quad_counter)=FSI(1)%IntElem(m1,i+quad_counter)+ &
         FSI(1)%NumNodesPool
       if (shift_from_zero_node.eq.1) then
        FSI(1)%IntElem(m1,i+quad_counter)=FSI(1)%IntElem(m1,i+quad_counter)+1
       endif
      enddo
     enddo  ! j1=0..k-3
    enddo  ! i, looping faces
    close(14)
    if (quad_counter.ne.numquads) then
     print *,"quad_counter.ne.numquads"
     stop
    endif

    print *,"i1=",i1
    do dir=1,3 
     print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
    enddo
    do dir=1,3 
     print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
    enddo
   enddo ! i1

   if ((probtype.eq.55).or. &
       (probtype.eq.58)) then
    theta=1.0
   else if (tcrit.le.tstart+1.0e-8) then
    theta=0.0
   else if (tcrit.ge.tfinish-1.0e-8) then
    theta=1.0
   else
    theta=(tcrit-tstart)/(tfinish-tstart)
   endif

   do i=1,FSI(1)%NumNodes
    if (probtype.eq.58) then 
     do dir=1,3
      FSI(1)%NodeVel_old(dir,i)=0.0
     enddo
     FSI(1)%NodeVel_old(1,i)=2.0*3.14159*plungeramp*plungerfreq*  &
       cos(2.0*3.14159*plungerfreq*curtime)
    else if (tfinish-tstart.gt.1.0e-8) then
     do dir=1,3
      FSI(1)%NodeVel_old(dir,i)= &
       (FSI(1)%Node_new(dir,i)-FSI(1)%Node_old(dir,i))/(tfinish-tstart)
     enddo
    else
     do dir=1,3
      FSI(1)%NodeVel_old(dir,i)=0.0
     enddo
    endif
    do dir=1,3
     FSI(1)%Node_current(dir,i)=theta*FSI(1)%Node_new(dir,i)+ &
      (1.0-theta)*FSI(1)%Node_old(dir,i)
     FSI(1)%NodeVel_new(dir,i)=FSI(1)%NodeVel_old(dir,i)
    enddo
   enddo  ! i=1,NumNodes

  else if (iread.eq.0) then
   ! do nothing
  else
   print *,"iread invalid"
   stop
  endif

  use_temp=0

   ! allocates Node,NodeVel,NodeForce,NodeTemp if ifirst==1
   ! It is assumed that the number of nodes and elements does not
   ! change from frame to frame.
   ! do_2nd_part=0 isout=1 (geominit)
  call init3_FSI(local_part_id,ifirst,0,1,1) 

  print *,"after geominit  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt,Theta_Dot,Theta ",timeB,dtB,TorqueVel,TorquePos
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
      xxblob(1),xxblob(3),newxxblob(1),newxxblob(3),radradblob
  print *," xnew=cos(theta)*(x-x_0)+sin(theta)*(z-z_0)+x_0"
  print *," znew=-sin(theta)*(x-x_0)+cos(theta)*(z-z_0)+z_0"

91     FORMAT(I7)
115    FORMAT(I4)

return
end subroutine geominit



! sphere centered almost at the origin with a radius of 5
! inflow at xlo, outflow at xhi
subroutine viorel_sphere_geominit(curtime,dt,ifirst,sdim,istop,istep)
IMPLICIT NONE

INTEGER_T :: ifirst
INTEGER_T :: i,j,k,i1,j1,m1,tempelem,sdim
INTEGER_T :: numquads,quad_counter
INTEGER_T, dimension(4) :: localElem
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
REAL_T, dimension(3) :: xvalbefore
REAL_T, dimension(3) :: invertfactor

INTEGER_T :: dir,istep,istop

INTEGER_T :: shift_from_zero_node
INTEGER_T :: shift_from_zero_face
INTEGER_T :: override_IntElemDim
character(100) :: dwave
INTEGER_T :: local_part_id

  local_part_id=1

  if (probtype.ne.5601) then
   print *,"probtype should be 5601"
   stop
  endif

  timeB=curtime
  dtB=dt

  if ((ifirst.eq.0).or.(ifirst.eq.1)) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  if (ifirst.eq.1) then
   dwave="viorel_sphere.txt"

   FSI(1)%NumNodesPool=0
   FSI(1)%NumIntElemsPool=0
   FSI(1)%IntElemDimPool=0

   invertfactor(1)=1.0
   invertfactor(2)=1.0
   invertfactor(3)=1.0

   xxblob(1)=0.0
   xxblob(2)=0.0
   xxblob(3)=0.0
   radradblob=1.0

   newxxblob(1)=0.0
   newxxblob(2)=0.0
   newxxblob(3)=0.0

   numquads=0
   print *,"opening to check to see how many quads to convert ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
   READ(14,91) FSI(1)%NumNodes
   print *,"NumNodes(check) ",FSI(1)%NumNodes
   READ(14,91) FSI(1)%NumIntElems
   print *,"NumIntElems(check) ",FSI(1)%NumIntElems
   READ(14,91) FSI(1)%IntElemDim
   print *,"IntElemDim(check) ",FSI(1)%IntElemDim
   do i=1,FSI(1)%NumNodes
    READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
   enddo
   READ(14,*) j
   do i=1,FSI(1)%NumIntElems
    READ(14,*) j
    READ(14,*) k
    if (k.eq.FSI(1)%NumIntElems) then ! use default 3 
     k=3
    endif
    if (k.eq.4) then
     numquads=numquads+1
    else if (k.ne.3) then
     print *,"elements must be triangles or quadrilaterals k=",k
     print *,"i,NumIntElems ",i,FSI(1)%NumIntElems
     stop
    endif
    do j1=1,k
     READ(14,*) tempelem
    enddo
   enddo ! sweeping elements
 
   close(14)

   print *,"numquads=",numquads

   print *,"opening ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

   READ(14,91) FSI(1)%NumNodes
   print *,"NumNodes ",FSI(1)%NumNodes
   READ(14,91) FSI(1)%NumIntElems
   print *,"NumIntElems ",FSI(1)%NumIntElems
   READ(14,91) FSI(1)%IntElemDim
   print *,"IntElemDim ",FSI(1)%IntElemDim
   print *,"overriding IntElemDim to be 3"
   FSI(1)%IntElemDim=3
 
   FSI(1)%NumIntElems=FSI(1)%NumIntElems+numquads

   if (ifirst.eq.1) then
    call init_FSI(local_part_id,1)
   else
    print *,"ifirst bust"
    stop
   endif

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   shift_from_zero_node=0
   shift_from_zero_face=0
   override_IntElemDim=0
   do i=1,FSI(1)%NumNodes
    READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
    if ((j.eq.0).and.(shift_from_zero_node.eq.0)) then
     shift_from_zero_node=1
     print *,"nodes in file start at 0; shifting to start at 1"
    endif
    if (shift_from_zero_node.eq.1) then
     j=j+1
    endif

    do dir=1,3
     xval(dir)=invertfactor(dir)* &
      (xvalbefore(dir)-xxblob(dir))/radradblob + newxxblob(dir)
    enddo

    do dir=1,3
     if ((minnode(dir).gt.xval(dir)).or.(i.eq.1)) then
      minnode(dir)=xval(dir)
     endif
     if ((maxnode(dir).lt.xval(dir)).or.(i.eq.1)) then
      maxnode(dir)=xval(dir)
     endif
     if ((minnodebefore(dir).gt.xvalbefore(dir)).or.(i.eq.1)) then
      minnodebefore(dir)=xvalbefore(dir)
     endif
     if ((maxnodebefore(dir).lt.xvalbefore(dir)).or.(i.eq.1)) then
      maxnodebefore(dir)=xvalbefore(dir)
     endif
    enddo
    
    do dir=1,3
     FSI(1)%Node_old(dir,i)=xval(dir)
    enddo
    do dir=1,3
     FSI(1)%Node_new(dir,i)=xval(dir)
    enddo
    if (i.ne.j) then
     print *,"vertex mismatch"
     print *,"NumNodesPool=",FSI(1)%NumNodesPool
     print *,"expected node index ",i-FSI(1)%NumNodesPool
     print *,"node index read from file ",j
     stop
    endif
   enddo !  reading nodes
      
   READ(14,*) j
   if (j.ne.FSI(1)%NumIntElems-numquads) then
    print *,"face mismatch"
    stop
   endif

   quad_counter=0
   do i=1,FSI(1)%NumIntElems-numquads
    READ(14,*) j
    if ((j.eq.0).and.(shift_from_zero_face.eq.0)) then
     shift_from_zero_face=1
     print *,"faces in file start at 0; shifting to start at 1"
    endif
    if (shift_from_zero_face.eq.1) then
     j=j+1
    endif

    if (i.ne.j) then
     print *,"face mismatch"
     print *,"NumIntElemsPool=",FSI(1)%NumIntElemsPool
     print *,"expected face index ",i-FSI(1)%NumIntElemsPool
     print *,"face index read from file ",j
     stop
    endif
    READ(14,*) k
    if (k.eq.FSI(1)%NumIntElems) then  ! use default 3
     k=3
    endif
    if ((k.ne.3).and.(k.ne.4)) then
     print *,"k invalid k=",k
     stop
    endif

    do j1=1,k
     READ(14,*) localElem(j1)
    enddo

    do j1=0,k-3
     if (j1.eq.1) then
      quad_counter=quad_counter+1
     endif
     FSI(1)%ElemData(1,i+quad_counter)=3   ! number of nodes in element
     FSI(1)%ElemData(2,i+quad_counter)=1   ! part number
     FSI(1)%ElemData(3,i+quad_counter)=0   ! singly wetted (set =1 doubly)
     if (j1.eq.0) then
      do dir=1,3
       FSI(1)%IntElem(dir,i+quad_counter)=localElem(dir)
      enddo
     else if (j1.eq.1) then
      FSI(1)%IntElem(1,i+quad_counter)=localElem(3)
      FSI(1)%IntElem(2,i+quad_counter)=localElem(4)
      FSI(1)%IntElem(3,i+quad_counter)=localElem(1)
     else
      print *,"j1 invalid"
      stop
     endif
     do m1=1,3
      if (shift_from_zero_node.eq.1) then
       FSI(1)%IntElem(m1,i+quad_counter)=FSI(1)%IntElem(m1,i+quad_counter)+1
      endif
     enddo
    enddo  ! j1=0..k-3
   enddo  ! i, looping faces
   close(14)
   if (quad_counter.ne.numquads) then
    print *,"quad_counter.ne.numquads"
    stop
   endif

   print *,"i1=",i1
   do dir=1,3 
    print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
   enddo
   do dir=1,3 
     print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
   enddo

   call init2_FSI(local_part_id)
  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif 

  use_temp=0 
  !do 2nd part=1 isout=1,viorel sphere geominit
  call init3_FSI(local_part_id,ifirst,1,1,1) 

  print *,"after viorel_sphere geominit  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt,Theta_Dot,Theta ",timeB,dtB,TorqueVel,TorquePos
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
      xxblob(1),xxblob(3),newxxblob(1),newxxblob(3),radradblob
  print *," xnew=cos(theta)*(x-x_0)+sin(theta)*(z-z_0)+x_0"
  print *," znew=-sin(theta)*(x-x_0)+cos(theta)*(z-z_0)+z_0"

91     FORMAT(I7)

return
end subroutine viorel_sphere_geominit




subroutine internal_inflow_geominit(curtime,dt,ifirst,sdim,istop,istep)
IMPLICIT NONE

INTEGER_T :: i,j,k,i1,j1,m1,tempelem,ifirst,sdim
INTEGER_T :: numquads,quad_counter
INTEGER_T, dimension(4) :: localElem
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
REAL_T, dimension(3) :: xvalbefore
REAL_T, dimension(3) :: invertfactor

INTEGER_T :: dir,istep,istop

INTEGER_T :: shift_from_zero_node
INTEGER_T :: shift_from_zero_face
INTEGER_T :: override_IntElemDim
character(100) :: dwave
INTEGER_T :: local_part_id

  local_part_id=1

  if (probtype.ne.5602) then
   print *,"probtype should be 5602"
   stop
  endif

  timeB=curtime
  dtB=dt

  if ((ifirst.eq.0).or.(ifirst.eq.1)) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  if (ifirst.eq.1) then
   dwave="internal_inflow.txt"

   FSI(1)%NumNodesPool=0
   FSI(1)%NumIntElemsPool=0
   FSI(1)%IntElemDimPool=0

   invertfactor(1)=1.0
   invertfactor(2)=1.0
   invertfactor(3)=1.0

   xxblob(1)=0.0
   xxblob(2)=0.0
   xxblob(3)=0.0
   radradblob=1.0

   newxxblob(1)=0.0
   newxxblob(2)=0.0
   newxxblob(3)=0.0

   numquads=0
   print *,"opening to check to see how many quads to convert ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
   READ(14,91) FSI(1)%NumNodes
   print *,"NumNodes(check) ",FSI(1)%NumNodes
   READ(14,91) FSI(1)%NumIntElems
   print *,"NumIntElems(check) ",FSI(1)%NumIntElems
   READ(14,91) FSI(1)%IntElemDim
   print *,"IntElemDim(check) ",FSI(1)%IntElemDim
   do i=1,FSI(1)%NumNodes
    READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
   enddo
   READ(14,*) j
   do i=1,FSI(1)%NumIntElems
    READ(14,*) j
    READ(14,*) k
    if (k.eq.FSI(1)%NumIntElems) then ! use default 3 
     k=3
    endif
    if (k.eq.4) then
     numquads=numquads+1
    else if (k.ne.3) then
     print *,"elements must be triangles or quadrilaterals k=",k
     print *,"i,NumIntElems ",i,FSI(1)%NumIntElems
     stop
    endif
    do j1=1,k
     READ(14,*) tempelem
    enddo
   enddo ! sweeping elements
 
   close(14)

   print *,"numquads=",numquads

   print *,"opening ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

   READ(14,91) FSI(1)%NumNodes
   print *,"NumNodes ",FSI(1)%NumNodes
   READ(14,91) FSI(1)%NumIntElems
   print *,"NumIntElems ",FSI(1)%NumIntElems
   READ(14,91) FSI(1)%IntElemDim
   print *,"IntElemDim ",FSI(1)%IntElemDim
   print *,"overriding IntElemDim to be 3"
   FSI(1)%IntElemDim=3
 
   FSI(1)%NumIntElems=FSI(1)%NumIntElems+numquads

   if (ifirst.eq.1) then
    call init_FSI(local_part_id,1)
   else
    print *,"ifirst bust"
    stop
   endif

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   shift_from_zero_node=0
   shift_from_zero_face=0
   override_IntElemDim=0
   do i=1,FSI(1)%NumNodes
    READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
    if ((j.eq.0).and.(shift_from_zero_node.eq.0)) then
     shift_from_zero_node=1
     print *,"nodes in file start at 0; shifting to start at 1"
    endif
    if (shift_from_zero_node.eq.1) then
     j=j+1
    endif

    do dir=1,3
     xval(dir)=invertfactor(dir)* &
      (xvalbefore(dir)-xxblob(dir))/radradblob + newxxblob(dir)
    enddo

    do dir=1,3
     if ((minnode(dir).gt.xval(dir)).or.(i.eq.1)) then
      minnode(dir)=xval(dir)
     endif
     if ((maxnode(dir).lt.xval(dir)).or.(i.eq.1)) then
      maxnode(dir)=xval(dir)
     endif
     if ((minnodebefore(dir).gt.xvalbefore(dir)).or.(i.eq.1)) then
      minnodebefore(dir)=xvalbefore(dir)
     endif
     if ((maxnodebefore(dir).lt.xvalbefore(dir)).or.(i.eq.1)) then
      maxnodebefore(dir)=xvalbefore(dir)
     endif
    enddo
    
    do dir=1,3
     FSI(1)%Node_old(dir,i)=xval(dir)
    enddo
    do dir=1,3
     FSI(1)%Node_new(dir,i)=xval(dir)
    enddo
    if (i.ne.j) then
     print *,"vertex mismatch"
     print *,"NumNodesPool=",FSI(1)%NumNodesPool
     print *,"expected node index ",i-FSI(1)%NumNodesPool
     print *,"node index read from file ",j
     stop
    endif
   enddo !  reading nodes
      
   READ(14,*) j
   if (j.ne.FSI(1)%NumIntElems-numquads) then
    print *,"face mismatch"
    stop
   endif

   quad_counter=0
   do i=1,FSI(1)%NumIntElems-numquads
    READ(14,*) j
    if ((j.eq.0).and.(shift_from_zero_face.eq.0)) then
     shift_from_zero_face=1
     print *,"faces in file start at 0; shifting to start at 1"
    endif
    if (shift_from_zero_face.eq.1) then
     j=j+1
    endif

    if (i.ne.j) then
     print *,"face mismatch"
     print *,"NumIntElemsPool=",FSI(1)%NumIntElemsPool
     print *,"expected face index ",i-FSI(1)%NumIntElemsPool
     print *,"face index read from file ",j
     stop
    endif
    READ(14,*) k
    if (k.eq.FSI(1)%NumIntElems) then  ! use default 3
     k=3
    endif
    if ((k.ne.3).and.(k.ne.4)) then
     print *,"k invalid k=",k
     stop
    endif

    do j1=1,k
     READ(14,*) localElem(j1)
    enddo

    do j1=0,k-3
     if (j1.eq.1) then
      quad_counter=quad_counter+1
     endif
     FSI(1)%ElemData(1,i+quad_counter)=3   ! number of nodes in element
     FSI(1)%ElemData(2,i+quad_counter)=1   ! part number
     FSI(1)%ElemData(3,i+quad_counter)=0   ! singly wetted (set =1 doubly)
     if (j1.eq.0) then
      do dir=1,3
       FSI(1)%IntElem(dir,i+quad_counter)=localElem(dir)
      enddo
     else if (j1.eq.1) then
      FSI(1)%IntElem(1,i+quad_counter)=localElem(3)
      FSI(1)%IntElem(2,i+quad_counter)=localElem(4)
      FSI(1)%IntElem(3,i+quad_counter)=localElem(1)
     else
      print *,"j1 invalid"
      stop
     endif
     do m1=1,3
      if (shift_from_zero_node.eq.1) then
       FSI(1)%IntElem(m1,i+quad_counter)=FSI(1)%IntElem(m1,i+quad_counter)+1
      endif
     enddo
    enddo  ! j1=0..k-3
   enddo  ! i, looping faces
   close(14)
   if (quad_counter.ne.numquads) then
    print *,"quad_counter.ne.numquads"
    stop
   endif

   print *,"i1=",i1
   do dir=1,3 
    print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
   enddo
   do dir=1,3 
     print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
   enddo

   call init2_FSI(local_part_id)
  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif 

  use_temp=0 
  call init3_FSI(local_part_id,ifirst,1,1,1)  ! do_2nd_part=1 isout=1

  print *,"after viorel_sphere geominit  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt,Theta_Dot,Theta ",timeB,dtB,TorqueVel,TorquePos
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
      xxblob(1),xxblob(3),newxxblob(1),newxxblob(3),radradblob
  print *," xnew=cos(theta)*(x-x_0)+sin(theta)*(z-z_0)+x_0"
  print *," znew=-sin(theta)*(x-x_0)+cos(theta)*(z-z_0)+z_0"

91     FORMAT(I7)

return
end subroutine internal_inflow_geominit


! overall_solid_advance and generate_new_triangles called every geom_interval
! time.
! This routine will be called just once at the very beginning or upon
! restart.
subroutine gearinit(curtime,dt,ifirst,sdim,istop,istep)
IMPLICIT NONE

INTEGER_T :: i,ifirst,sdim
INTEGER_T :: numquads
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
REAL_T, dimension(3) :: xvalbefore
REAL_T, dimension(3) :: invertfactor
INTEGER_T :: dir,istep,istop
character(100) :: dwave
INTEGER_T local_ifirst
INTEGER_T local_part_id

  local_part_id=1
  local_ifirst=1

  timeB=curtime
  dtB=dt

  invertfactor(1)=1.0
  invertfactor(2)=1.0
  invertfactor(3)=1.0

  xxblob(1)=0.0
  xxblob(2)=0.0
  xxblob(3)=0.0
  newxxblob(1)=0.0
  newxxblob(2)=0.0
  newxxblob(3)=0.0
  radradblob=0.01 ! units are expected in cm

  if (axis_dir.eq.0) then
   dwave="gear1_geom.cas"
  else if (axis_dir.eq.1) then
   dwave="cylinder_geom.cas"
  else if (axis_dir.eq.2) then
   dwave="gear1_geom_May27.cas"
  else
   print *,"axis_dir invalid"
   stop
  endif
 
  numquads=0

  print *,"opening ",dwave
  OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

  READ(14,*) FSI(1)%NumNodes,FSI(1)%NumIntElems
  FSI(1)%IntElemDim=3

  call init_FSI(local_part_id,1)

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
   maxnodebefore(dir)=-1.0e+10
   minnodebefore(dir)=1.0e+10
  enddo

  do i=1,FSI(1)%NumNodes
   READ(14,*) xvalbefore(1),xvalbefore(2),xvalbefore(3)
   do dir=1,3
    xval(dir)=invertfactor(dir)* &
      (xvalbefore(dir)-xxblob(dir))/radradblob + newxxblob(dir)
   enddo

   do dir=1,3
    if ((minnode(dir).gt.xval(dir)).or.(i.eq.1)) then
     minnode(dir)=xval(dir)
    endif
    if ((maxnode(dir).lt.xval(dir)).or.(i.eq.1)) then
     maxnode(dir)=xval(dir)
    endif
    if ((minnodebefore(dir).gt.xvalbefore(dir)).or.(i.eq.1)) then
     minnodebefore(dir)=xvalbefore(dir)
    endif
    if ((maxnodebefore(dir).lt.xvalbefore(dir)).or.(i.eq.1)) then
     maxnodebefore(dir)=xvalbefore(dir)
    endif
   enddo
    
   do dir=1,3
    FSI(1)%Node_old(dir,i)=xval(dir)
   enddo
   do dir=1,3
    FSI(1)%Node_new(dir,i)=xval(dir)
   enddo
      
  enddo

  do i=1,FSI(1)%NumIntElems
   READ(14,*) FSI(1)%IntElem(3,i),FSI(1)%IntElem(2,i),FSI(1)%IntElem(1,i)

   FSI(1)%ElemData(1,i)=3   ! number of nodes in element
   FSI(1)%ElemData(2,i)=1   ! part number
   FSI(1)%ElemData(3,i)=0   ! singly wetted
  enddo  ! i, looping faces
  close(14)

  do dir=1,3 
   print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
  enddo
  do dir=1,3 
   print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
  enddo

  call init2_FSI(local_part_id)

  use_temp=0 
  ! do_2nd_part=0  isout=1 (gearinit)
  call init3_FSI(local_part_id,local_ifirst,0,1,1)  

  do i=1,FSI(1)%NumNodes
   do dir=1,3
    FSI(1)%Node(dir,i)=FSI(1)%Node_new(dir,i)
   enddo
  enddo

return
end subroutine gearinit




! this subroutine generates a speed vs. time profile for 
! use with the first type of inflow 
subroutine timefluct(cur_time,value)
IMPLICIT NONE

  REAL_T, intent(in) :: cur_time
  REAL_T, intent(out) :: value
  REAL_T inittime, medtime, endtime, cur_timeW

  inittime = 0.0
  medtime = 0.2
  endtime = 0.43
!  inittime = 0.0
!  medtime = 0.15
!  endtime = 0.33
  cur_timeW = cur_time
  if (cur_time.ge.1.0) then
!  if (cur_time.ge.0.78) then
   cur_timeW = cur_time-floor(cur_time)
  endif
  if (cur_timeW.lt.inittime) then
   value = 0.05
  else if ((cur_timeW.ge.inittime).and.(cur_timeW.lt.medtime)) then
   value = (cur_timeW-inittime)/(medtime-inittime)
  else if ((cur_timeW.ge.medtime).and.(cur_timeW.lt.endtime)) then
   value = (endtime-cur_timeW)/(endtime-medtime)
  else
   value = 0.05
  endif

return
end subroutine timefluct

subroutine whale_geominit(curtime,dt,ifirst,sdim,istop,istep, &
  CLSVOF_curtime,CLSVOF_dt)
IMPLICIT NONE

INTEGER_T :: i,ifirst,sdim
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xtemp,vtemp
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
REAL_T, dimension(3) :: xvalbefore
REAL_T, dimension(3) :: invertfactor
INTEGER_T :: dir,istep,istop
INTEGER_T :: shift_from_zero_node
INTEGER_T :: shift_from_zero_face
INTEGER_T :: override_IntElemDim,whale_type
REAL_T :: CLSVOF_curtime,CLSVOF_dt
REAL_T :: STEPSPERIOD,LL_CLSVOF,UU_CLSVOF,TT_CLSVOF
REAL_T :: LL_DUFFY,TT_DUFFY,UU_DUFFY
INTEGER_T :: local_part_id

  local_part_id=1

  if (probtype.ne.562) then
   print *,"probtype must be 562 for animated whale"
   stop
  endif

  timeB=curtime
  dtB=dt

  FSI(1)%NumNodesPool=0
  FSI(1)%NumIntElemsPool=0
  FSI(1)%IntElemDimPool=0

  invertfactor(1)=1.0
  invertfactor(2)=1.0
  invertfactor(3)=1.0

  whale_type=2  ! whale normal

! open mouth
  if (whale_type.eq.3) then

! xval(dir)=invertfactor(dir)* &
!       (xvalbefore(dir)-xxblob(dir))/radradblob + newxxblob(dir)
!
! xval(1)=xtemp(2) xval(2)=xtemp(1) xval(3)=xtemp(3)

   xxblob(1)=0.0 ! will become y
   xxblob(2)=-1.0 ! will become -x
   invertfactor(2)=-1.0
   xxblob(3)=0.0  ! will become z
   newxxblob(1)=0.0
   newxxblob(2)=1.0
   newxxblob(3)=0.5
   radradblob=7.5     ! 12.5/1.67=7.5

! was 10 for whale file dated December, 2007  (whalenormal)
! xval(1)=xtemp(3) xval(2)=xtemp(1) xval(3)=xtemp(2)
  else if (whale_type.eq.2) then
   xxblob(1)=0.0
   xxblob(2)=0.0
   xxblob(3)=0.0
   newxxblob(1)=0.0
   newxxblob(2)=0.5
   newxxblob(3)=0.9   ! was 1.0
   radradblob=10.87   ! was 12.5, now 10.87
! whalepregnant
  else if (whale_type.eq.6) then
   xxblob(1)=0.0 ! will become y
   xxblob(3)=-1.0 ! will become -x
   invertfactor(3)=-1.0
   invertfactor(2)=-1.0
   xxblob(2)=0.0  ! will become z
   newxxblob(1)=0.0
   newxxblob(3)=1.0
   newxxblob(2)=0.5
   radradblob=10.87   ! was 12.5, now 10.87
! whaletailup 
  else if (whale_type.eq.0) then
   xxblob(1)=0.0 ! will become y
   xxblob(2)=-1.0 ! will become -x
   invertfactor(2)=-1.0
   xxblob(3)=0.0  ! will become z
   newxxblob(1)=0.0
   newxxblob(2)=1.0
   newxxblob(3)=0.5
   radradblob=10.87   ! was 12.5, now 10.87
! whaletaildown
  else if (whale_type.eq.5) then
   xxblob(1)=0.0 ! will become y
   xxblob(2)=-1.0 ! will become -x
   invertfactor(2)=-1.0
   xxblob(3)=0.0  ! will become z
   newxxblob(1)=0.0
   newxxblob(2)=1.0
   newxxblob(3)=0.5
   radradblob=10.87   ! was 12.5, now 10.87
  else
   print *,"whale_type bad in whale_geominit probtype,axis_dir ", &
       probtype,axis_dir
   stop
  endif

  FSI(1)%NumNodes=whale_nodes 
  FSI(1)%NumIntElems=whale_cells
  FSI(1)%IntElemDim=3

  if (ifirst.eq.1) then
   call init_FSI(local_part_id,0)
  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif 

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
   maxnodebefore(dir)=-1.0e+10
   minnodebefore(dir)=1.0e+10
  enddo

  shift_from_zero_node=0
  shift_from_zero_face=0
  override_IntElemDim=0

  do i=1,FSI(1)%NumNodes
   xvalbefore(1)=XX(i)
   xvalbefore(2)=YY(i)
   xvalbefore(3)=ZZ(i)

   do dir=1,3
    xval(dir)=invertfactor(dir)* &
      (xvalbefore(dir)-xxblob(dir))/radradblob + newxxblob(dir)
   enddo

   do dir=1,3
    xtemp(dir)=xval(dir)
   enddo

! open mouth whale
   if (whale_type.eq.3) then

    xval(1)=xtemp(2)
    xval(2)=xtemp(1)
    xval(3)=xtemp(3)

! whale file used at beginning of summer 2007 (whalenormal)
   else if (whale_type.eq.2) then
    xval(1)=xtemp(3)
    xval(2)=xtemp(1)
    xval(3)=xtemp(2)
! whaletaildown
   else if (whale_type.eq.5) then
    xval(1)=xtemp(2)
    xval(2)=xtemp(1)
    xval(3)=xtemp(3)
! whaletailup 
   else if (whale_type.eq.0) then
    xval(1)=xtemp(2)
    xval(2)=xtemp(1)
    xval(3)=xtemp(3)
! whale pregnant
   else if (whale_type.eq.6) then
    xval(1)=xtemp(3)
    xval(2)=xtemp(1)
    xval(3)=xtemp(2)
   else
    print *,"bad whale_type whale_geominit probtype,axis_dir ", &
      probtype,axis_dir
    stop
   endif

   do dir=1,3
    if ((minnode(dir).gt.xval(dir)).or.(i.eq.1)) then
     minnode(dir)=xval(dir)
    endif
    if ((maxnode(dir).lt.xval(dir)).or.(i.eq.1)) then
     maxnode(dir)=xval(dir)
    endif
    if ((minnodebefore(dir).gt.xvalbefore(dir)).or.(i.eq.1)) then
     minnodebefore(dir)=xvalbefore(dir)
    endif
    if ((maxnodebefore(dir).lt.xvalbefore(dir)).or.(i.eq.1)) then
     maxnodebefore(dir)=xvalbefore(dir)
    endif
   enddo ! dir
    
   do dir=1,3
    FSI(1)%Node_old(dir,i)=xval(dir)
    FSI(1)%Node_new(dir,i)=xval(dir)
   enddo
  enddo ! i=1,NumNodes

  do i=1,FSI(1)%NumIntElems
   FSI(1)%ElemData(1,i)=3   ! number of nodes in element
   FSI(1)%ElemData(2,i)=1   ! part number
   FSI(1)%ElemData(3,i)=0   ! singly wetted
      ! IntElem initialized in runonce
  enddo
 
  print *,"minnodebefore,maxnodebefore ", &
   minnodebefore(1),minnodebefore(2),minnodebefore(3), &
   maxnodebefore(1),maxnodebefore(2),maxnodebefore(3)
  print *,"minnode,maxnode ",minnode(1),minnode(2),minnode(3), &
     maxnode(1),maxnode(2),maxnode(3)

  STEPSPERIOD=4.0*STEPS_DUFFY
  LL_CLSVOF=WHALE_LENGTH/1.32
  UU_CLSVOF=1.0
  TT_CLSVOF=LL_CLSVOF/UU_CLSVOF
  LL_DUFFY=WHALE_LENGTH/(1.32*radradblob)
  TT_DUFFY=10.0/(STEPSPERIOD*DT_DUFFY)
  UU_DUFFY=LL_DUFFY/TT_DUFFY
  print *,"UU_DUFFY, UU_CLSVOF ",UU_DUFFY,UU_CLSVOF

  do i=1,FSI(1)%NumNodes
   vtemp(1)=UU(i)
   vtemp(2)=VV(i)
   vtemp(3)=WW(i)
   do dir=1,3
    vtemp(dir)=vtemp(dir)*invertfactor(dir)*UU_DUFFY/UU_CLSVOF
   enddo
   FSI(1)%NodeVel_old(1,i)=vtemp(3)
   FSI(1)%NodeVel_old(2,i)=vtemp(1)
   FSI(1)%NodeVel_old(3,i)=vtemp(2)
   do dir=1,3
    FSI(1)%Node_current(dir,i)=FSI(1)%Node_new(dir,i)
    FSI(1)%NodeVel_new(dir,i)=FSI(1)%NodeVel_old(dir,i)
   enddo
  enddo  ! i=1,NumNodes

  use_temp=0 
  ! do_2nd_part=0  isout=1 (whalegeominit)
  call init3_FSI(local_part_id,ifirst,0,1,1)  

  print *,"CLSVOF_curtime,CLSVOF_dt ",CLSVOF_curtime,CLSVOF_dt
  print *,"after whale_geominit  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt ",timeB,dtB
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
      xxblob(1),xxblob(3),newxxblob(1),newxxblob(3),radradblob


return
end subroutine whale_geominit


subroutine initpaddle(curtime,dt,sdim,istop,istep, &
  paddle_pos,paddle_vel)
IMPLICIT NONE

INTEGER_T :: i,j,k,j1,sdim
INTEGER_T :: inode,iface
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xtemp
INTEGER_T :: dir,istep,istop
character(20) :: dwave,dwave2
REAL_T :: xx,zz
REAL_T, intent(in) :: paddle_pos,paddle_vel
INTEGER_T :: local_ifirst
INTEGER_T :: local_part_id

  local_part_id=1
  local_ifirst=1

  timeB=curtime
  dtB=0.0
  TorquePos=paddle_pos
  TorqueVel=paddle_vel

  xxblob(1)=-1.8
  xxblob(3)=0.8
  xxblob(2)=4.3  ! vertical
  xxblobwall(1)=-1.8
  xxblobwall(3)=0.8
  xxblobwall(2)=4.3  ! vertical
  newxxblob(1)=1.45  
  newxxblob(3)=1.0
  newxxblob(2)=1.0  ! vertical
  newxxblobwall(1)=1.05  
  newxxblobwall(3)=1.0
  newxxblobwall(2)=1.0  ! vertical
  radradblob=60.0
  radradblobwall=60.0
  denpaddle=0.0035
  dampingpaddle=0.01

  dwave="h/paddle.txt"
  dwave2="h/wall.txt"

  print *,"opening ",dwave
  OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
  print *,"opening ",dwave2
  OPEN(unit=15,file=dwave2,access='sequential',form="formatted",status='old')

  READ(14,91) FSI(1)%NumNodesPaddle
  print *,"NumNodesPaddle ",FSI(1)%NumNodesPaddle
  READ(14,91) FSI(1)%NumIntElemsPaddle
  print *,"NumIntElemsPaddle ",FSI(1)%NumIntElemsPaddle
  READ(14,91) FSI(1)%IntElemDimPaddle
  print *,"IntElemDimPaddle ",FSI(1)%IntElemDimPaddle

  READ(15,91) FSI(1)%NumNodes
  print *,"NumNodes ",FSI(1)%NumNodes
  READ(15,91) FSI(1)%NumIntElems
  print *,"NumIntElems ",FSI(1)%NumIntElems
  READ(15,91) FSI(1)%IntElemDim
  print *,"IntElemDim ",FSI(1)%IntElemDim

  FSI(1)%NumNodes=FSI(1)%NumNodes+FSI(1)%NumNodesPaddle
  FSI(1)%NumIntElems=FSI(1)%NumIntElems+FSI(1)%NumIntElemsPaddle
  if (FSI(1)%IntElemDim.lt.FSI(1)%IntElemDimPaddle) then
   FSI(1)%IntElemDim=FSI(1)%IntElemDimPaddle
  endif

  call init_FSI(local_part_id,1)

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
  enddo

  do inode=1,FSI(1)%NumNodes
   if (inode.le.FSI(1)%NumNodesPaddle) then
    READ(14,93) j,xval(1),xval(2),xval(3)
    if (inode.ne.j) then
     print *,"inode,j mismatch"
     stop
    endif
    do dir=1,3
     xval(dir)=(xval(dir)-xxblob(dir))/radradblob + newxxblob(dir)
    enddo
   else
    READ(15,93) j,xval(1),xval(2),xval(3)
    if (inode.ne.j+FSI(1)%NumNodesPaddle) then
     print *,"inode,j mismatch"
     stop
    endif
    do dir=1,3
     xval(dir)=(xval(dir)-xxblobwall(dir))/radradblobwall + newxxblobwall(dir)
    enddo
   endif


   do dir=1,3
    xtemp(dir)=xval(dir)
   enddo
   xval(1)=xtemp(1)
   xval(2)=xtemp(3)
   xval(3)=xtemp(2)
   do dir=1,3
    if ((minnode(dir).gt.xval(dir)).or.(inode.eq.1)) then
     minnode(dir)=xval(dir)
    endif
    if ((maxnode(dir).lt.xval(dir)).or.(inode.eq.1)) then
     maxnode(dir)=xval(dir)
    endif
   enddo
    
   do dir=1,3
    FSI(1)%Node_old(dir,inode)=xval(dir)
    FSI(1)%Node_new(dir,inode)=xval(dir)
   enddo
      
  enddo  ! inode=1,NumNodes
 
  READ(14,91) j
  if (j.ne.FSI(1)%NumIntElemsPaddle) then
   print *,"face mismatch"
   stop
  endif

  READ(15,91) j
  if (j.ne.FSI(1)%NumIntElems-FSI(1)%NumIntElemsPaddle) then
   print *,"face mismatch"
   stop
  endif

  do iface=1,FSI(1)%NumIntElems
   if (iface.le.FSI(1)%NumIntElemsPaddle) then
    READ(14,91) j
    if (iface.ne.j) then
     print *,"face mismatch"
     stop
    endif
    READ(14,91) k
    if (k.gt.FSI(1)%IntElemDim) then
     print *,"14 too many vertices k,IntElemDim ",k,FSI(1)%IntElemDim
     stop
    endif
    do j1=1,k
     READ(14,91) FSI(1)%IntElem(j1,iface) 
    enddo
   else
    READ(15,91) j
    if (iface.ne.j+FSI(1)%NumIntElemsPaddle) then
     print *,"face mismatch"
     stop
    endif
    READ(15,91) k
    if (k.gt.FSI(1)%IntElemDim) then
     print *,"15 too many vertices k,IntElemDim ",k,FSI(1)%IntElemDim
     stop
    endif
    do j1=1,k
     READ(15,91) FSI(1)%IntElem(j1,iface) 
     FSI(1)%IntElem(j1,iface)=FSI(1)%IntElem(j1,iface)+FSI(1)%NumNodesPaddle
    enddo
   endif

   FSI(1)%ElemData(1,iface)=k   ! number of nodes in element
   FSI(1)%ElemData(2,iface)=1   ! part number
   FSI(1)%ElemData(3,iface)=0   ! singly wetted
   FSI(1)%ElemData(3,iface)=1  ! doubly wetted
   if (iface.gt.FSI(1)%NumIntElemsPaddle) then
    FSI(1)%ElemData(3,iface)=1  ! doubly wetted
   endif
  enddo  ! iface, looping faces

  close(14)
  close(15)
 
  print *,"minnode,maxnode ",minnode(1),minnode(2),minnode(3), &
      maxnode(1),maxnode(2),maxnode(3)

  call init2_FSI(local_part_id)

  do inode=1,FSI(1)%NumNodesPaddle
   xx=FSI(1)%Node_current(1,inode)-newxxblob(1)
   zz=FSI(1)%Node_current(3,inode)-newxxblob(3)
   FSI(1)%Node_new(1,inode)=xx*cos(TorquePos)+zz*sin(TorquePos)+newxxblob(1)
   FSI(1)%Node_new(3,inode)=-xx*sin(TorquePos)+zz*cos(TorquePos)+newxxblob(3)
   xx=FSI(1)%Node_new(1,inode)-newxxblob(1)
   zz=FSI(1)%Node_new(3,inode)-newxxblob(3)
   FSI(1)%NodeVel_new(1,inode)=TorqueVel*zz
   FSI(1)%NodeVel_new(3,inode)=-TorqueVel*xx
  enddo

  use_temp=0

  ! do_2nd_part=0 isout=1 (initpaddle)
  call init3_FSI(local_part_id,local_ifirst,0,1,1)  

  do i=1,FSI(1)%NumNodes
   do dir=1,3
    FSI(1)%Node(dir,i)=FSI(1)%Node_new(dir,i)
   enddo
  enddo

91     FORMAT(I7)
93     FORMAT(i7,f11.3,f11.3,f11.3)

return
end subroutine initpaddle


subroutine initship(curtime,dt,sdim,istop,istep, &
  paddle_pos,paddle_vel)
IMPLICIT NONE

INTEGER_T :: i,sdim
INTEGER_T :: inode,iface
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xtemp
INTEGER_T :: dir,istep,istop
INTEGER_T :: filler
character(40) :: dwave
REAL_T, intent(in) :: paddle_pos,paddle_vel
INTEGER_T :: local_ifirst
INTEGER_T :: local_part_id

  local_part_id=1
  local_ifirst=1

  timeB=curtime
  dtB=0.0
  TorquePos=paddle_pos
  TorqueVel=paddle_vel

  xxblob(1)=0.0
  xxblob(2)=0.0
  xxblob(3)=0.0
  newxxblob(1)=0.0
  newxxblob(2)=0.0
  newxxblob(3)=0.0  
  radradblob=1.0
  denpaddle=1.0
  dampingpaddle=0.01

  if (axis_dir.eq.1) then
   dwave="5415_froude4136.cas"
  else if (axis_dir.eq.2) then
   dwave="cavity.cas"
  else if (axis_dir.eq.3) then
   dwave="cutout.cas"
  else 
   print *,"bad axis_dir initship probtype,axis_dir ",probtype,axis_dir
   stop
  endif
   

  print *,"opening ",dwave
  OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

  READ(14,191) FSI(1)%NumNodes,FSI(1)%NumIntElems
  FSI(1)%IntElemDim=3

  call init_FSI(local_part_id,1)

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
  enddo

  do inode=1,FSI(1)%NumNodes
   READ(14,193) xval(1),xval(2),xval(3)
   do dir=1,3
     xval(dir)=(xval(dir)-xxblob(dir))/radradblob + newxxblob(dir)
   enddo

   do dir=1,3
    xtemp(dir)=xval(dir)
   enddo
   xval(1)=-xtemp(1)
   xval(2)=xtemp(2)
   xval(3)=xtemp(3)
   do dir=1,3
    if ((minnode(dir).gt.xval(dir)).or.(inode.eq.1)) then
     minnode(dir)=xval(dir)
    endif
    if ((maxnode(dir).lt.xval(dir)).or.(inode.eq.1)) then
     maxnode(dir)=xval(dir)
    endif
   enddo
    
   do dir=1,3
    FSI(1)%Node_old(dir,inode)=xval(dir)
    FSI(1)%Node_new(dir,inode)=xval(dir)
   enddo
      
  enddo  ! inode=1,NumNodes
 
  do iface=1,FSI(1)%NumIntElems
   READ(14,194) FSI(1)%IntElem(3,iface),FSI(1)%IntElem(2,iface), &
     FSI(1)%IntElem(1,iface),filler

   FSI(1)%ElemData(1,iface)=3   ! number of nodes in element
   FSI(1)%ElemData(2,iface)=1   ! part number
   if ((axis_dir.eq.1).or.(axis_dir.eq.3)) then
    FSI(1)%ElemData(3,iface)=0   ! singly wetted
   else if (axis_dir.eq.2) then
    FSI(1)%ElemData(3,iface)=1   ! doubly wetted
   else
    print *,"bad axis_dir initship2 probtype,axis_dir ",probtype,axis_dir
    stop
   endif
!   FSI(1)%ElemData(3,iface)=1   ! doubly wetted

  enddo  ! iface, looping faces

  close(14)
 
  print *,"minnode,maxnode ",minnode(1),minnode(2),minnode(3), &
      maxnode(1),maxnode(2),maxnode(3)

  call init2_FSI(local_part_id)

  use_temp=0

  ! do_2nd_part=0 isout=1 (initship)
  call init3_FSI(local_part_id,local_ifirst,0,1,1)  

  do i=1,FSI(1)%NumNodes
   do dir=1,3
    FSI(1)%Node(dir,i)=FSI(1)%Node_new(dir,i)
   enddo
  enddo

191    FORMAT(I12,I12)
193    FORMAT(E15.11,E15.11,E15.11)
194    FORMAT(I12,I12,I12,I12)

return
end subroutine initship

! if probtype==701,538,541 then generate_new_triangles is not needed.
subroutine overall_solid_advance(CLSVOF_curtime,CLSVOF_dt, &
  part_id,ioproc,isout)
IMPLICIT NONE

INTEGER_T :: part_id    
INTEGER_T :: ifirst,ioproc,isout
REAL_T :: CLSVOF_curtime,CLSVOF_dt
REAL_T :: STEPSPERIOD,LL_CLSVOF,UU_CLSVOF,TT_CLSVOF,whale_dt

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif

 ifirst=0
 sci_istep=sci_istep+1

  ! animated whale  
 if ((probtype.eq.562).and.(axis_dir.eq.4)) then
  STEPSPERIOD=4.0*STEPS_DUFFY
  LL_CLSVOF=WHALE_LENGTH/1.32
  UU_CLSVOF=1.0
  TT_CLSVOF=LL_CLSVOF/UU_CLSVOF
  whale_dt=PERIOD_TAIL/(TT_CLSVOF*STEPSPERIOD)

  do while (CLSVOF_whale_time.lt.CLSVOF_curtime)
   print *,"advancing whale geometry: CLSVOF_whale_time=",CLSVOF_whale_time
   print *,"advancing whale geometry: CLSVOF_curtime=",CLSVOF_curtime
   call new_geometry(RR,SS,TT,UU,VV,WW,XX,YY,ZZ,whale_list,whale_angle, &
    whale_timestep,whale_spring,whale_counter,whale_counter_real, &
    whale_nodes,whale_cells,whale_X_init,whale_Y_init,whale_Z_init)
   CLSVOF_whale_time=CLSVOF_whale_time+whale_dt
  enddo
   ! 562
  call whale_geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,CLSVOF_whale_time,CLSVOF_dt)
  ! animated heart
 else if (probtype.eq.57) then
  call geominit(CLSVOF_curtime,CLSVOF_dt,ifirst,sci_sdim,sci_istop,sci_istep)
 else if ((probtype.eq.52).or. &
     ((probtype.eq.56).and.  &
      ((axis_dir.eq.0).or.(axis_dir.eq.2).or.(axis_dir.eq.3))).or.  &
     (probtype.eq.561).or. &
     ((probtype.eq.55).and.(raddust.eq.zero)).or. &
     (probtype.eq.58).or. &
     (probtype.eq.562)) then
  call geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)
 else if (probtype.eq.5600) then ! dog
  call geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)
 else if (probtype.eq.5601) then ! viorel sphere
  call viorel_sphere_geominit(sci_curtime,sci_dt,ifirst,sci_sdim, &
     sci_istop,sci_istep)
 else if (probtype.eq.5602) then ! internal inflow
  call internal_inflow_geominit(sci_curtime,sci_dt,ifirst,sci_sdim, &
     sci_istop,sci_istep)
 else if (probtype.eq.5700) then ! microfluidics channel
  call initchannel(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep)
 else if (probtype.eq.53) then
  call initinjector(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout)
 else if ((probtype.eq.531).or. &
          (probtype.eq.5501)) then
  call initinjector(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout)

   ! ifirst=0 (in overall_solid_advance)
 else if ((probtype.eq.536).or.(probtype.eq.537).or. &
          (probtype.eq.538).or.(probtype.eq.539).or. &
          (probtype.eq.541)) then
  call initinjector(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout)    
   ! ifirst=0 (in overall_solid_advance)
 else if (probtype.eq.701) then
  call initflapping(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout)    
 else if ((probtype.eq.563).or. & ! gear
          (probtype.eq.50).or. & ! paddle
          (probtype.eq.9)) then ! ship
  call advance_solid(sci_sdim,sci_curtime,sci_dt,sci_istop,sci_istep,part_id)
 else if (probtype.eq.400) then
  call init_gingerbread2D(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout) 
 else if (probtype.eq.401) then
  call init_helix(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout) 
 else
   ! ifirst=0 (in overall_solid_advance)
  call init_from_cas(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout) 
 endif

return
end subroutine overall_solid_advance

subroutine overall_solid_init(CLSVOFtime,ioproc,part_id,isout)
use global_utility_module

IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T :: ifirst,ioproc,isout
REAL_T :: paddle_pos,paddle_vel,CLSVOFtime,CLSVOF_dt
REAL_T :: STEPSPERIOD,LL_CLSVOF,UU_CLSVOF,TT_CLSVOF,whale_dt
INTEGER_T :: nmat
INTEGER_T :: ctml_part_id
INTEGER_T :: fsi_part_id

 nmat=num_materials

 if ((ioproc.ne.1).and.(ioproc.ne.0)) then
  print *,"ioproc invalid"
  stop
 endif
 ctml_part_id=CTML_partid_map(part_id)
 fsi_part_id=FSI_partid_map(part_id)

 if (((ctml_part_id.ge.1).and. &
      (ctml_part_id.le.CTML_NPARTS)).or. &
     ((fsi_part_id.ge.1).and. &
      (fsi_part_id.le.FSI_NPARTS))) then

  normal_invert=0
  sci_sdim=3    
  exclusive_doubly_wetted=0
 
  if ((ctml_part_id.ge.1).and. &
      (ctml_part_id.le.CTML_NPARTS)) then

   exclusive_doubly_wetted=1
   sci_curtime=CLSVOFtime  ! modify this if restarting?
   sci_dt=0.0              ! modify this if restarting?
   sci_istop=0
   sci_istep=0
   ifirst=1
   call CTML_init_sci(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep, &
     ioproc,part_id,isout)

  else if (ctml_part_id.eq.0) then
   ! do nothing
  else
   print *,"ctml_part_id invalid"
   stop
  endif

  if ((fsi_part_id.ge.1).and. &
      (fsi_part_id.le.FSI_NPARTS)) then

   ! note: 5600 dog, 5601 viorel sphere, 5602 internal inflow
   ! heart
   if (probtype.eq.57) then
    exclusive_doubly_wetted=1
   endif

! sends header message to fluids code
! sends node message to fluids code 
! reads force message

   sci_curtime=CLSVOFtime
   sci_dt=0.0

! --------------- MODIFY THESE VARIABLES IF RESTARTING! ---------------
   paddle_pos=0.0
   paddle_vel=0.0
   sci_curtime=CLSVOFtime
   sci_dt=0.0
! --------------- END SECTION TO MODIFY IF RESTARTING! ---------------

   sci_istop=0 
   sci_istep=0
   ifirst=1 

   if ((probtype.eq.562).and.(axis_dir.eq.4)) then
    STEPSPERIOD=4.0*STEPS_DUFFY
    LL_CLSVOF=WHALE_LENGTH/1.32
    UU_CLSVOF=1.0
    TT_CLSVOF=LL_CLSVOF/UU_CLSVOF
    whale_dt=PERIOD_TAIL/(TT_CLSVOF*STEPSPERIOD)

    whalein="h/whalenormal.txt"
    call getinfo(whalein,whale_nodes,whaleout)
    allocate(RR(whale_nodes))
    allocate(SS(whale_nodes))
    allocate(TT(whale_nodes))
    allocate(UU(whale_nodes))
    allocate(VV(whale_nodes))
    allocate(WW(whale_nodes))
    allocate(XX(whale_nodes))
    allocate(YY(whale_nodes))
    allocate(ZZ(whale_nodes))
    allocate(whale_X_init(whale_nodes))
    allocate(whale_Y_init(whale_nodes))
    allocate(whale_Z_init(whale_nodes))
    allocate(whale_spring(whale_nodes))
    allocate(whale_list(whale_nodes,20))
   
    CLSVOF_whale_time=0.0 
    call runonce(RR,SS,TT,UU,VV,WW,XX,YY,ZZ,whale_list,whale_angle, &
      whale_timestep,whale_spring,whale_counter,whale_counter_real, &
      whale_nodes,whale_cells,whale_X_init,whale_Y_init, &
      whale_Z_init,whaleout) 
    do while (CLSVOF_whale_time.lt.CLSVOFtime)
     print *,"advancing whale geometry: CLSVOF_whale_time=",CLSVOF_whale_time
     print *,"advancing whale geometry: CLSVOFtime=",CLSVOFtime
     call new_geometry(RR,SS,TT,UU,VV,WW,XX,YY,ZZ,whale_list,whale_angle, &
       whale_timestep,whale_spring,whale_counter,whale_counter_real, &
       whale_nodes,whale_cells,whale_X_init,whale_Y_init,whale_Z_init)
     CLSVOF_whale_time=CLSVOF_whale_time+whale_dt
    enddo

    CLSVOF_dt=0.0
    call whale_geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
      sci_istep,CLSVOF_whale_time,CLSVOF_dt)
! heart
   else if (probtype.eq.57) then
    call geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)
   else if ((probtype.eq.52).or. &
       (probtype.eq.56).or.(probtype.eq.561).or. &
       (probtype.eq.58).or.(probtype.eq.562)) then
! pregnant whale
    if ((probtype.eq.562).and.(axis_dir.eq.6)) then
     normal_invert=0
    endif
    call geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)
   else if (probtype.eq.563) then
    call gearinit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)
   else if (probtype.eq.5600) then ! dog
    call geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)
   else if (probtype.eq.5601) then ! viorel sphere
    normal_invert=1
    call viorel_sphere_geominit(sci_curtime,sci_dt,ifirst, &
      sci_sdim,sci_istop,sci_istep)
   else if (probtype.eq.5602) then ! internal inflow
    normal_invert=1
    call internal_inflow_geominit(sci_curtime,sci_dt,ifirst, &
      sci_sdim,sci_istop,sci_istep)
   else if (probtype.eq.5700) then ! microfluidics channel
    call initchannel(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)    
   else if ((probtype.eq.53).or. &
            (probtype.eq.531).or. &
            (probtype.eq.5501)) then
    call initinjector(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep, &
      ioproc,part_id,isout) 

     ! ifirst=1 (overall_solid_init)
   else if ((probtype.eq.536).or.(probtype.eq.537).or. &
            (probtype.eq.538).or.(probtype.eq.539).or. &
            (probtype.eq.541)) then
    call initinjector(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep, &
      ioproc,part_id,isout) 
     ! ifirst=1 (overall_solid_init)
   else if (probtype.eq.701) then
    call initflapping(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep, &
      ioproc,part_id,isout) 
   else if (probtype.eq.50) then
    call initpaddle(sci_curtime,sci_dt,sci_sdim,sci_istop,sci_istep, &
      paddle_pos,paddle_vel)
   else if (probtype.eq.9) then
    call initship(sci_curtime,sci_dt,sci_sdim,sci_istop,sci_istep, &
      paddle_pos,paddle_vel)
   else if (probtype.eq.400) then
    call init_gingerbread2D(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop, &
     sci_istep,ioproc,part_id,isout) 
   else if (probtype.eq.401) then
    call init_helix(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop, &
     sci_istep,ioproc,part_id,isout) 
   else
    call init_from_cas(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep, &
      ioproc,part_id,isout) 
   endif

  else if (fsi_part_id.eq.0) then
   ! do nothing
  else
   print *,"fsi_part_id invalid"
   stop
  endif

 else if ((ctml_part_id.eq.0).and. &
          (fsi_part_id.eq.0)) then
  ! do nothing
 else
  print *,"ctml_part_id or fsi_part_id invalid"
  stop
 endif

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"after initialize solid dt=",sci_dt
 endif

return
end subroutine overall_solid_init


subroutine advance_solid(sdim,curtime,dt,istop,istep,part_id)
IMPLICIT NONE

INTEGER_T :: sdim
INTEGER_T :: part_id
REAL_T :: curtime,dt
INTEGER_T :: istop
INTEGER_T :: i,dir,istep

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif

 do i=1,FSI(part_id)%NumNodes
  do dir=1,3
   FSI(part_id)%Node_old(dir,i)=FSI(part_id)%Node_new(dir,i)
   FSI(part_id)%NodeVel_old(dir,i)=FSI(part_id)%NodeVel_new(dir,i)
  enddo
  do dir=1,6
   FSI(part_id)%NodeForce_old(dir,i)=FSI(part_id)%NodeForce_new(dir,i)
  enddo
 enddo

 do i=1,FSI(part_id)%NumNodes
  do dir=1,3
   FSI(part_id)%Node(dir,i)=FSI(part_id)%Node_new(dir,i)
   FSI(part_id)%NodeVel(dir,i)=FSI(part_id)%NodeVel_new(dir,i)
  enddo
  do dir=1,6
   FSI(part_id)%NodeForce(dir,i)=FSI(part_id)%NodeForce_new(dir,i)
  enddo

   ! in: advance_solid
  FSI(part_id)%NodeMass(i)=one
  FSI(part_id)%NodeTemp(i)=FSI(part_id)%NodeTemp_new(i)
 enddo

 dtB=0.0
 if (curtime.gt.timeB) then
  dtB=curtime-timeB
  timeB=curtime
 endif

  ! this first case will never happen since initinjector is called 
  ! instead of this routine. 
 if ((probtype.eq.538).or.(probtype.eq.541)) then  

  do i=1,FSI(part_id)%NumNodes
   do dir=1,3
    FSI(part_id)%Node_new(dir,i)=FSI(part_id)%Node_new(dir,i)+ &
     dtB*FSI(part_id)%NodeVel_new(dir,i)
   enddo
  enddo

  ! this second case will never happen since initflapping is called 
  ! instead of this routine. 
 else if (probtype.eq.701) then

  do i=1,FSI(part_id)%NumNodes
   do dir=1,3
    FSI(part_id)%Node_new(dir,i)=FSI(part_id)%Node_new(dir,i)+ &
     dtB*FSI(part_id)%NodeVel_new(dir,i)
   enddo
  enddo

 else

   ! default for ship, gear, ...
  do i=1,FSI(part_id)%NumNodes
   do dir=1,3
    FSI(part_id)%NodeVel_new(dir,i)=0.0
   enddo
   do dir=1,6
    FSI(part_id)%NodeForce_new(dir,i)=0.0
   enddo
  enddo

 endif

return
end subroutine advance_solid





! --------------------  SOLID ADVANCE STUFF ENDS HERE --------------


! calls CTML_DELTA or hsprime
subroutine check_force_weight( &
  xmap3D,inode,ielem, &
  xc,part_id,time,dx, &
  force_weight,force_vector)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif
IMPLICIT NONE

INTEGER_T, dimension(3), intent(in) :: xmap3D
INTEGER_T, intent(in) :: part_id
INTEGER_T, intent(in) :: inode,ielem
REAL_T, intent(in) :: time 
REAL_T, dimension(3), intent(out) :: force_vector
REAL_T, dimension(3), intent(in) :: dx
REAL_T, dimension(3), intent(in) :: xc
REAL_T, intent(out) :: force_weight
REAL_T, dimension(3) :: xtarget
REAL_T, dimension(3) :: xfoot
REAL_T, dimension(3) :: velparm
INTEGER_T nmat,nodes_per_elem,Node_repeat_count
INTEGER_T dir
INTEGER_T inode_raw
REAL_T dist_scale,df,support_size,line_mass

 nmat=num_materials

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif
 if (time.lt.zero) then
  print *,"time invalid"
  stop
 endif

 nodes_per_elem=FSI(part_id)%ElemDataBIG(1,ielem)

 if (nodes_per_elem.ne.3) then
  print *,"nodes_per_elem invalid"
  stop
 endif
 if ((inode.lt.1).or.(inode.gt.nodes_per_elem)) then
  print *,"inode invalid"
  stop
 endif

 inode_raw=FSI(part_id)%IntElemBIG(inode,ielem)
 if ((inode_raw.lt.1).or.(inode_raw.gt.FSI(part_id)%NumNodesBIG)) then
  print *,"inode_raw invalid"
  stop
 endif

 do dir=1,3
  if (dx(dir).le.zero) then
   print *,"dx invalid"
   stop
  endif
  xfoot(dir)=FSI(part_id)%NodeBIG(dir,inode_raw)
  velparm(dir)=zero
 enddo ! dir=1..3
 call get_target_from_foot(xfoot,xtarget, &
   velparm,time,part_id)

 Node_repeat_count=FSI(part_id)%ElemNodeCountBIG(inode_raw)

 if (Node_repeat_count.ge.1) then

  force_weight=one
  do dir=1,3

   if ((xmap3D(dir).eq.1).or. &
       (xmap3D(dir).eq.2).or. &
       (xmap3D(dir).eq.AMREX_SPACEDIM)) then
    dist_scale=abs(xc(dir)-xtarget(dir))/dx(dir)
    if (CTML_FSI_flagF(nmat).eq.1) then ! FSI_flag==4
#ifdef MVAHABFSI
     call CTML_DELTA(dir,dist_scale,df)
#else
     print *,"define MVAHABFSI"
     stop
#endif
    else if (CTML_FSI_flagF(nmat).eq.0) then
     support_size=two
     df=hsprime(dist_scale,support_size)
    else
     print *,"CTML_FSI_flagF(nmat) invalid"
     stop
    endif 
    if (df.ge.zero) then
     force_weight=force_weight*df/dx(dir)
    else
     print *,"df invalid"
     stop
    endif
   else if ((xmap3D(dir).eq.0).and.(AMREX_SPACEDIM.eq.2)) then
    ! do nothing
   else
    print *,"xmap3D(dir) invalid"
    stop
   endif
 
   force_vector(dir)=FSI(part_id)%NodeForceBIG(dir,inode_raw)
  enddo ! dir=1..3
 
  line_mass=FSI(part_id)%NodeMassBIG(inode_raw)
  force_weight=force_weight*line_mass/Node_repeat_count

 else
  print *,"Node_repeat_count invalid"
  stop
 endif
  
return
end subroutine check_force_weight



subroutine checkinpoint(xclosest,normal_closest, &
  inode,elemnum, &
  unsigned_mindist, &
  xc,inplane,part_id,time,dx)
IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T, intent(in) :: inode,elemnum
INTEGER_T, intent(inout) :: inplane
REAL_T, intent(in) :: time 
REAL_T, dimension(3), intent(inout) :: xclosest
REAL_T, dimension(3), intent(inout) :: normal_closest
REAL_T, dimension(3), intent(in) :: dx
REAL_T, dimension(3), intent(in) :: xc
REAL_T, intent(inout) :: unsigned_mindist
REAL_T :: curdist,dotprod,mag
INTEGER_T :: dir
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3) :: xfoot
REAL_T, dimension(3) :: xfoot_pert
REAL_T, dimension(3) :: xtarget
REAL_T, dimension(3) :: xtarget_pert
REAL_T, dimension(3) :: ntarget
REAL_T, dimension(3) :: velparm

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif
 if (time.lt.zero) then
  print *,"time invalid"
  stop
 endif

 nodes_per_elem=FSI(part_id)%ElemDataBIG(1,elemnum)

 if (nodes_per_elem.ne.3) then
  print *,"nodes_per_elem invalid"
  stop
 endif

 do dir=1,3
  if (dx(dir).le.zero) then
   print *,"dx invalid"
   stop
  endif
  xfoot(dir)=FSI(part_id)%NodeBIG(dir,FSI(part_id)%IntElemBIG(inode,elemnum))
  xfoot_pert(dir)=xfoot(dir)+0.1*dx(dir)* &
    FSI(part_id)%NodeNormalBIG(dir,FSI(part_id)%IntElemBIG(inode,elemnum))
  velparm(dir)=zero
 enddo
 call get_target_from_foot(xfoot,xtarget, &
   velparm,time,part_id)
 call get_target_from_foot(xfoot_pert,xtarget_pert, &
   velparm,time,part_id)
 mag=zero
 do dir=1,3
  ntarget(dir)=xtarget_pert(dir)-xtarget(dir)
  mag=mag+ntarget(dir)**2
 enddo
 mag=sqrt(mag)
 if (mag.le.zero) then
  print *,"mag invalid checkinpoint 0"
  stop
 endif
 do dir=1,3
  ntarget(dir)=ntarget(dir)/mag
 enddo

 dotprod=0.0
 do dir=1,3
  dotprod=dotprod+ntarget(dir)*(xc(dir)-xtarget(dir))
 enddo

 call xdist(xtarget,xc,curdist)
 if (curdist.lt.zero) then
  print *,"curdist invalid"
  stop
 endif

 if ((curdist.lt.unsigned_mindist).or. &
     (inplane.eq.0)) then
  inplane=1
  unsigned_mindist=curdist
  do dir=1,3
   xclosest(dir)=xtarget(dir)
   normal_closest(dir)=ntarget(dir)
  enddo
 endif

return
end subroutine checkinpoint



subroutine checkinline(xclosest,normal_closest, &
  inode,elemnum, &
  unsigned_mindist, &
  xc,inplane,part_id,time,dx)
IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T, intent(in) :: inode,elemnum
INTEGER_T, intent(inout) :: inplane
REAL_T, intent(in) :: time 
REAL_T, dimension(3), intent(inout) :: xclosest
REAL_T, dimension(3), intent(inout) :: normal_closest
REAL_T, dimension(3), intent(in) :: dx
REAL_T, dimension(3), intent(in) :: xc
REAL_T, dimension(2,3) :: xnode,nnode
REAL_T, intent(inout) :: unsigned_mindist
REAL_T, dimension(3) :: xnot,normal
REAL_T :: dottop,dotbot,t,curdist,dotprod,mag
INTEGER_T :: dir
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3) :: xfoot
REAL_T, dimension(3) :: xfoot_pert
REAL_T, dimension(3) :: xtarget
REAL_T, dimension(3) :: xtarget_pert
REAL_T, dimension(3) :: ntarget
REAL_T, dimension(3) :: velparm

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif
 if (time.lt.zero) then
  print *,"time invalid"
  stop
 endif

 if ((inode.ge.1).and.(inode.le.2)) then

  nodes_per_elem=FSI(part_id)%ElemDataBIG(1,elemnum)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem invalid"
   stop
  endif

  do dir=1,3
   if (dx(dir).le.zero) then
    print *,"dx invalid"
    stop
   endif
   xfoot(dir)=FSI(part_id)%NodeBIG(dir,FSI(part_id)%IntElemBIG(inode,elemnum))
   xfoot_pert(dir)=xfoot(dir)+0.1*dx(dir)* &
    FSI(part_id)%NodeNormalBIG(dir,FSI(part_id)%IntElemBIG(inode,elemnum))
   velparm(dir)=zero
  enddo
  call get_target_from_foot(xfoot,xtarget, &
    velparm,time,part_id)
  call get_target_from_foot(xfoot_pert,xtarget_pert, &
    velparm,time,part_id)
  mag=zero
  do dir=1,3
   ntarget(dir)=xtarget_pert(dir)-xtarget(dir)
   mag=mag+ntarget(dir)**2
  enddo
  mag=sqrt(mag)
  if (mag.le.zero) then
   print *,"mag invalid checkinline 0"
   stop
  endif
  do dir=1,3
   nnode(1,dir)=ntarget(dir)/mag
   xnode(1,dir)=xtarget(dir)
  enddo

  do dir=1,3
   xfoot(dir)=FSI(part_id)%NodeBIG(dir,FSI(part_id)%IntElemBIG(inode+1,elemnum))
   xfoot_pert(dir)=xfoot(dir)+0.1*dx(dir)* &
    FSI(part_id)%NodeNormalBIG(dir,FSI(part_id)%IntElemBIG(inode+1,elemnum))
   velparm(dir)=zero
  enddo
  call get_target_from_foot(xfoot,xtarget, &
    velparm,time,part_id)
  call get_target_from_foot(xfoot_pert,xtarget_pert, &
    velparm,time,part_id)
  mag=zero
  do dir=1,3
   ntarget(dir)=xtarget_pert(dir)-xtarget(dir)
   mag=mag+ntarget(dir)**2
  enddo
  mag=sqrt(mag)
  if (mag.le.zero) then
   print *,"mag invalid checkinline 1"
   stop
  endif
  do dir=1,3
   nnode(2,dir)=ntarget(dir)/mag
   xnode(2,dir)=xtarget(dir)
  enddo
 
  dottop=0.0
  dotbot=0.0
  do dir=1,3
   dottop=dottop+(xnode(2,dir)-xnode(1,dir))*(xnode(2,dir)-xc(dir))
   dotbot=dotbot+(xnode(2,dir)-xnode(1,dir))**2
  enddo
  dotbot=dotbot+1.0E-14
  t=dottop/dotbot
  if ((t.ge.-element_buffer_tol).and. &
      (t.le.1.0+element_buffer_tol)) then
   do dir=1,3
    xnot(dir)=t*xnode(1,dir)+(1.0-t)*xnode(2,dir)
    normal(dir)=t*nnode(1,dir)+(1.0-t)*nnode(2,dir)
   enddo

   dotprod=0.0
   do dir=1,3
    dotprod=dotprod+normal(dir)*(xc(dir)-xnot(dir))
   enddo

   call xdist(xnot,xc,curdist)
   if (curdist.lt.zero) then
    print *,"curdist invalid"
    stop
   endif

   if ((curdist.lt.unsigned_mindist).or. &
       (inplane.eq.0)) then
    inplane=1
    unsigned_mindist=curdist
    do dir=1,3
     xclosest(dir)=xnot(dir)
     normal_closest(dir)=normal(dir)
    enddo
   endif
  endif ! -tol<=t<=1+tol

 else if (inode.eq.3) then
  ! do nothing
 else
  print *,"inode invalid in checkinline"
  stop
 endif

return
end subroutine checkinline

subroutine checkinplane(xclosest,elemnum,inplane, &
  minnode,maxnode,element_scale,part_id,time)
IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T, intent(in) :: elemnum
REAL_T, intent(in) :: time
REAL_T, dimension(3), intent(in) :: minnode,maxnode,xclosest
REAL_T, intent(in) :: element_scale
INTEGER_T, intent(out) :: inplane
REAL_T, dimension(3,3) :: xnode,AA,AI
INTEGER_T :: dir,i,j,k
INTEGER_T :: nodes_per_elem
REAL_T :: det
REAL_T :: scaled_tol
REAL_T, dimension(3) :: tx
REAL_T, dimension(3) :: v1,v2,v1xv2
REAL_T, dimension(3) :: xfoot
REAL_T, dimension(3) :: xtarget
REAL_T, dimension(3) :: velparm

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif
 if (time.lt.zero) then
  print *,"time invalid"
  stop
 endif

 nodes_per_elem=FSI(part_id)%ElemDataBIG(1,elemnum)
 if (nodes_per_elem.ne.3) then
  print *,"nodes_per_elem invalid"
  stop
 endif

 scaled_tol=element_scale*element_buffer_tol
 if (scaled_tol.le.zero) then
  print *,"scaled_tol invalid"
  stop
 endif

 inplane=1
 do dir=1,3
  if ((xclosest(dir).lt.minnode(dir)-scaled_tol).or. &
      (xclosest(dir).gt.maxnode(dir)+scaled_tol)) then
   inplane=0
  endif
 enddo

 if (inplane.eq.1) then

  do i=1,3
   do dir=1,3
    xfoot(dir)=FSI(part_id)%NodeBIG(dir,FSI(part_id)%IntElemBIG(i,elemnum))
    velparm(dir)=zero
   enddo
   call get_target_from_foot(xfoot,xtarget, &
      velparm,time,part_id)
   do dir=1,3
    xnode(i,dir)=xtarget(dir)
   enddo
  enddo  ! i=1..3

    ! xnode(1)-xnode(1) is mapped to (0,0,0)
    ! xnode(2)-xnode(1)=v1 is mapped to (1,0,0)
    ! xnode(3)-xnode(1)=v2 is mapped to (0,1,0) 
    ! v1 x v2 is mapped to (0,0,1)
    ! Let A map from unit space to real space
    ! A (0,0,0) = (0,0,0)
    ! A (1,0,0) = v1 => first column of A is v1
    ! A (0,1,0) = v2 => second column of A is v2
    ! A (0,0,1) = v1 x v2 => third column of A is v1 x v2
    ! AI maps from real space back to unit space

  do dir=1,3
   v1(dir)=xnode(2,dir)-xnode(1,dir)
   v2(dir)=xnode(3,dir)-xnode(1,dir)
  enddo  
  v1xv2(1)=v1(2)*v2(3)-v1(3)*v2(2)  
  v1xv2(2)=v1(3)*v2(1)-v1(1)*v2(3)  
  v1xv2(3)=v1(1)*v2(2)-v1(2)*v2(1)  
  do dir=1,3
   AA(dir,1)=v1(dir)
   AA(dir,2)=v2(dir)
   AA(dir,3)=v1xv2(dir)
  enddo
  
  det=AA(1,1)*(AA(2,2)*AA(3,3)-AA(2,3)*AA(3,2))- &
      AA(1,2)*(AA(2,1)*AA(3,3)-AA(2,3)*AA(3,1))+ &
      AA(1,3)*(AA(2,1)*AA(3,2)-AA(2,2)*AA(3,1))
  if (abs(det).lt.1.0e-15) then
   inplane=0
  else
   det=1.0/det
   AI(1,1)=+(AA(2,2)*AA(3,3)-AA(2,3)*AA(3,2))
   AI(2,1)=-(AA(2,1)*AA(3,3)-AA(2,3)*AA(3,1))
   AI(3,1)=+(AA(2,1)*AA(3,2)-AA(2,2)*AA(3,1))
   AI(1,2)=-(AA(1,2)*AA(3,3)-AA(3,2)*AA(1,3))
   AI(2,2)=+(AA(1,1)*AA(3,3)-AA(1,3)*AA(3,1))
   AI(3,2)=-(AA(1,1)*AA(3,2)-AA(3,1)*AA(1,2))
   AI(1,3)=+(AA(1,2)*AA(2,3)-AA(2,2)*AA(1,3))
   AI(2,3)=-(AA(1,1)*AA(2,3)-AA(2,1)*AA(1,3))
   AI(3,3)=+(AA(1,1)*AA(2,2)-AA(1,2)*AA(2,1))
   do i=1,3
   do j=1,3
    AI(i,j)=AI(i,j)*det
   enddo
   enddo
  
   do i=1,3
    tx(i)=0.0
    do k=1,3
     tx(i)=tx(i)+AI(i,k)*(xclosest(k)-xnode(1,k))
    enddo
   enddo
  
   if ((tx(1).lt.-element_buffer_tol).or. &
       (tx(1).gt.1.0+element_buffer_tol).or. &
       (tx(2).lt.-element_buffer_tol).or. &
       (tx(2).gt.1.0+element_buffer_tol).or. &
       (tx(1)+tx(2).gt.1.0+element_buffer_tol)) then
    inplane=0
   endif
  endif ! det.ne.0
 
 endif  ! inplane.eq.1

return
end subroutine checkinplane



! normal points from solid to fluid
! solid nodes are ordered clockwise when viewed from the fluid.
! for 2d problems, it is assumed that the 3rd node is equal to the 2nd
! node, except that the 3rd node extends OUT of the paper. (positive z)
subroutine scinormal(elemnum,normal,part_id,time)
IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T, intent(in) :: elemnum
REAL_T, dimension(3), intent(out) :: normal
REAL_T, intent(in) :: time
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3,3) :: nodesave
REAL_T, dimension(3) :: nodeavg
REAL_T, dimension(3) :: xfoot
REAL_T, dimension(3) :: xtarget
REAL_T, dimension(3) :: velparm
REAL_T, dimension(2,3) :: vec
REAL_T :: dist
INTEGER_T :: i
INTEGER_T :: dir
INTEGER_T :: local_normal_invert

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif
 if (time.lt.zero) then
  print *,"time invalid"
  stop
 endif

 nodes_per_elem=FSI(part_id)%ElemDataBIG(1,elemnum)
 if (nodes_per_elem.gt.3) then
  print *,"nodes_per_elem>3 not supported"
  stop
 endif

 do dir=1,3
  nodeavg(dir)=0.0
 enddo

 do i=1,3
  do dir=1,3
   xfoot(dir)=FSI(part_id)%NodeBIG(dir,FSI(part_id)%IntElemBIG(i,elemnum))
   velparm(dir)=zero
  enddo
  call get_target_from_foot(xfoot,xtarget, &
      velparm,time,part_id)

  do dir=1,3
   nodesave(i,dir)=xtarget(dir)
   nodeavg(dir)=nodeavg(dir)+nodesave(i,dir)
  enddo
 enddo ! i=1..3

 do dir=1,3
  nodeavg(dir)=nodeavg(dir)/3.0
 enddo

 do i=1,2
  do dir=1,3
   vec(i,dir)=nodesave(i+1,dir)-nodesave(i,dir)
  enddo
 enddo

 normal(1)=vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2)
 normal(2)=vec(1,3)*vec(2,1)-vec(1,1)*vec(2,3)
 normal(3)=vec(1,1)*vec(2,2)-vec(2,1)*vec(1,2)

 dist=sqrt(normal(1)**2+normal(2)**2+normal(3)**2)
 if (dist.gt.1.0e-15) then
  do dir=1,3
   normal(dir)=normal(dir)/dist
  enddo
 endif

 local_normal_invert=normal_invert

 if (local_normal_invert.eq.1) then
  do dir=1,3
   normal(dir)=-normal(dir)
  enddo
 else if (local_normal_invert.ne.0) then
  print *,"local_normal_invert must be 0 or 1"
  stop
 endif

return
end subroutine scinormal


subroutine sciarea(elemnum,area,part_id)
IMPLICIT NONE

INTEGER_T :: part_id
INTEGER_T, intent(in) :: elemnum
REAL_T, intent(out) :: area
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3,3) :: nodesave
INTEGER_T :: i,j
REAL_T :: aa,bb,cc

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%partID.ne.part_id) then
  print *,"FSI(part_id)%partID.ne.part_id"
  stop
 endif

 nodes_per_elem=FSI(part_id)%ElemDataBIG(1,elemnum)
 if (nodes_per_elem.ne.3) then
  print *,"nodes in element?  elem,nodes_per_elem ",elemnum,nodes_per_elem   
  stop
 endif

 do i=1,3
  do j=1,3
   nodesave(i,j)=FSI(part_id)%NodeBIG(j,FSI(part_id)%IntElemBIG(i,elemnum))
  enddo
 enddo

 do i=1,2
 do j=1,3
   nodesave(i,j)=nodesave(i,j)-nodesave(3,j)
 enddo
 enddo
 aa=nodesave(1,2)*nodesave(2,3)-nodesave(2,2)*nodesave(1,3)
 bb=nodesave(1,1)*nodesave(2,3)-nodesave(2,1)*nodesave(1,3)
 cc=nodesave(1,1)*nodesave(2,2)-nodesave(2,1)*nodesave(1,2)
 area=0.5*sqrt(aa**2+bb**2+cc**2)

 if ((area.eq.0.0).and.(1.eq.0)) then
   print *,"area=0 x1,y1,x2,y2 ",nodesave(1,1),nodesave(1,2), &
    nodesave(2,1),nodesave(2,2)
   print *,"nodes_per_elem,elemnum ",nodes_per_elem,elemnum
   stop
 endif

return
end subroutine sciarea

! called from FORT_HEADERMSG when FSI_operation==4 and FSI_sub_operation==0
! isout==1 => verbose
subroutine CLSVOF_clear_lag_data(ioproc,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

INTEGER_T :: ioproc,isout
INTEGER_T :: nmat
INTEGER_T :: part_id
INTEGER_T :: ctml_part_id
INTEGER_T :: fsi_part_id
INTEGER_T :: dir,inode,num_nodes

 nmat=num_materials

 if (TOTAL_NPARTS.ge.1) then

  if (CTML_FSI_flagF(nmat).eq.1) then ! FSI_flag==4
#ifdef MVAHABFSI
   call CTML_RESET_ARRAYS(); ! vel_fib=zero  force_fib=zero
#else
   print *,"define MEHDI_VAHAB_FSI in GNUmakefile"
   stop
#endif
  else if (CTML_FSI_flagF(nmat).eq.0) then
   ! do nothing
  else 
   print *,"CTML_FSI_flagF(nmat) invalid"
   stop
  endif

  do part_id=1,TOTAL_NPARTS

   ctml_part_id=CTML_partid_map(part_id)
   fsi_part_id=FSI_partid_map(part_id)

   if ((ctml_part_id.gt.0).or. &
       (fsi_part_id.gt.0)) then

    if ((ctml_part_id.gt.0).and. &
        (ctml_part_id.le.CTML_NPARTS)) then
     do inode=1,ctml_n_fib_nodes(ctml_part_id)
      do dir=1,AMREX_SPACEDIM
       ctml_fib_vel(ctml_part_id,inode,dir)=zero
      enddo
      do dir=1,2*AMREX_SPACEDIM
       ctml_fib_frc(ctml_part_id,inode,dir)=zero
      enddo
      ctml_fib_mass(ctml_part_id,inode)=one
     enddo ! inode=1,ctml_n_fib_nodes(part_id)
    else if (ctml_part_id.eq.0) then
     ! do nothing
    else
     print *,"ctml_part_id invalid"
     stop
    endif

    num_nodes=FSI(part_id)%NumNodes
    if (num_nodes.le.0) then
     print *,"num_nodes invalid"
     stop
    endif
    do inode=1,num_nodes
     do dir=1,3
      FSI(part_id)%NodeVel(dir,inode)=zero
      FSI(part_id)%NodeVel_old(dir,inode)=zero
      FSI(part_id)%NodeVel_new(dir,inode)=zero
     enddo
     FSI(part_id)%NodeMass(inode)=one
     do dir=1,6
      FSI(part_id)%NodeForce(dir,inode)=zero
      FSI(part_id)%NodeForce_old(dir,inode)=zero
      FSI(part_id)%NodeForce_new(dir,inode)=zero
     enddo
    enddo ! inode=1..num_nodes
   else if ((ctml_part_id.eq.0).and. &
            (fsi_part_id.eq.0)) then
    ! do nothing
   else 
    print *,"ctml_part_id or fsi_part_id invalid"
    stop
   endif

  enddo ! part_id=1,TOTAL_NPARTS

 else
  print *,"TOTAL_NPARTS invalid"
  stop
 endif

return
end subroutine CLSVOF_clear_lag_data


! called from FORT_HEADERMSG when FSI_operation==4 and FSI_sub_operation==2
! isout==1 => verbose
subroutine CLSVOF_sync_lag_data(ioproc,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif
#if (mpi_activate==1)
use mpi
#endif

IMPLICIT NONE

INTEGER_T :: ioproc,isout
INTEGER_T :: ierr1
INTEGER_T :: nmat
double precision, dimension(:), allocatable :: sync_velocity
double precision, dimension(:), allocatable :: temp_velocity

INTEGER_T part_id
INTEGER_T ctml_part_id
INTEGER_T fsi_part_id
INTEGER_T num_nodes,sync_dim,inode,inode_fiber,dir

 nmat=num_materials

 if (TOTAL_NPARTS.ge.1) then

  do part_id=1,TOTAL_NPARTS

   ctml_part_id=CTML_partid_map(part_id)
   fsi_part_id=FSI_partid_map(part_id)

   if ((ctml_part_id.gt.0).or. &
       (fsi_part_id.gt.0)) then

    num_nodes=FSI(part_id)%NumNodes
    if (num_nodes.le.0) then
     print *,"num_nodes invalid"
     stop
    endif
    sync_dim=num_nodes*3
    allocate(sync_velocity(sync_dim))
    allocate(temp_velocity(sync_dim))
    do inode=1,num_nodes
    do dir=1,3
     sync_velocity(3*(inode-1)+dir)=FSI(part_id)%NodeVel(dir,inode)
     temp_velocity(3*(inode-1)+dir)=zero
    enddo
    enddo
    ierr1=0
#if (mpi_activate==1)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr1)
    call MPI_ALLREDUCE(sync_velocity,temp_velocity,sync_dim, &
     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr1)
    do inode=1,num_nodes
    do dir=1,3
     sync_velocity(3*(inode-1)+dir)=temp_velocity(3*(inode-1)+dir)
    enddo
    enddo
#elif (mpi_activate==0)
    ! do nothing
#else
    print *,"mpi_activate invalid"
    stop
#endif
    do inode=1,num_nodes
    do dir=1,3
     FSI(part_id)%NodeVel(dir,inode)=sync_velocity(3*(inode-1)+dir)
     FSI(part_id)%NodeVel_old(dir,inode)=sync_velocity(3*(inode-1)+dir)
     FSI(part_id)%NodeVel_new(dir,inode)=sync_velocity(3*(inode-1)+dir)
    enddo
    enddo

    deallocate(sync_velocity)
    deallocate(temp_velocity)

    if ((ctml_part_id.ge.1).and. &
        (ctml_part_id.le.CTML_NPARTS)) then
     do inode_fiber=1,ctml_n_fib_active_nodes(ctml_part_id)
      inode=inode_fiber
      do dir=1,AMREX_SPACEDIM
       ctml_fib_vel(ctml_part_id,inode_fiber,dir)= &
        FSI(part_id)%NodeVel(dir,inode)
      enddo
     enddo
    else if (ctml_part_id.eq.0) then
     ! do nothing
    else 
     print *,"ctml_part_id invalid"
     stop
    endif
   else if ((ctml_part_id.eq.0).and. &
            (fsi_part_id.eq.0)) then
    ! do nothing
   else
    print *,"ctml_part_id or fsi_part_id invalid"
    stop
   endif

  enddo !part_id=1,TOTAL_NPARTS

  if (CTML_FSI_flagF(nmat).eq.1) then ! FSI_flag==4
#ifdef MVAHABFSI
   call CTML_SET_VELOCITY(CTML_NPARTS, &
    ctml_max_n_fib_nodes,ctml_fib_vel) !vel_fib=ctml_fib_vel
#else
   print *,"define MEHDI_VAHAB_FSI in GNUmakefile"
#endif
  else if (CTML_FSI_flagF(nmat).eq.0) then
   ! do nothing
  else 
   print *,"CTML_FSI_flagF(nmat) invalid"
   stop
  endif

 else
  print *,"TOTAL_NPARTS invalid"
  stop
 endif

return
end subroutine CLSVOF_sync_lag_data



! called from FORT_HEADERMSG when FSI_operation==0
! isout==1 => verbose
subroutine CLSVOF_ReadHeader( &
  FSI_refine_factor, &
  FSI_bounding_box_ngrow, &
  nparts_in, &
  im_solid_map_in, &  ! im_part=im_solid_map_in(partid)+1
  h_small, &
  dx_max_level, &
  CTML_FSI_INIT, &
  CLSVOFtime,problo,probhi, &
  ioproc,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

INTEGER_T :: nparts_in
INTEGER_T im_solid_map_in(nparts_in)
INTEGER_T :: initflag,ioproc,isout
INTEGER_T :: CTML_FSI_INIT
REAL_T :: dx_max_level(AMREX_SPACEDIM)
REAL_T :: h_small
REAL_T :: CLSVOFtime
REAL_T problo(3),probhi(3)
INTEGER_T :: test_NPARTS
INTEGER_T :: part_id
INTEGER_T :: nmat
INTEGER_T :: dir
INTEGER_T :: im_part
INTEGER_T :: ctml_part_id,fsi_part_id
INTEGER_T FSI_refine_factor(num_materials)
INTEGER_T FSI_bounding_box_ngrow(num_materials)
INTEGER_T im_sanity_check

  nmat=num_materials

  do im_sanity_check=1,nmat
   if ((FSI_refine_factor(im_sanity_check).lt.0).or. &
       (FSI_refine_factor(im_sanity_check).gt.100)) then
    print *,"FSI_refine_factor(im_sanity_check) invalid"
    stop
   endif
   if (FSI_bounding_box_ngrow(im_sanity_check).ne.3) then
    print *,"FSI_bounding_box_ngrow(im_sanity_check) invalid"
    stop
   endif
  enddo

  if ((nparts_in.lt.1).or.(nparts_in.ge.nmat)) then
   print *,"nparts_in invalid"
   stop
  endif

  TOTAL_NPARTS=nparts_in
  ctml_part_id=0
  fsi_part_id=0

  do part_id=1,TOTAL_NPARTS

   im_solid_mapF(part_id)=im_solid_map_in(part_id)
   im_part=im_solid_mapF(part_id)+1
   if (CTML_FSI_mat(nmat,im_part).eq.1) then ! FSI_flag==4
    ctml_part_id=ctml_part_id+1
    CTML_partid_map(part_id)=ctml_part_id
   else if (CTML_FSI_mat(nmat,im_part).eq.0) then
    CTML_partid_map(part_id)=0
   else
    print *,"CTML_FSI_mat(nmat,im_part) invalid"
    stop
   endif
   if ((FSI_flag(im_part).eq.2).or. & ! prescribed solid (CAD)
       (FSI_flag(im_part).eq.6).or. & ! prescribed ice (CAD)
       (FSI_flag(im_part).eq.7)) then ! prescribed fluid (CAD)
    fsi_part_id=fsi_part_id+1
    FSI_partid_map(part_id)=fsi_part_id
   else if ((FSI_flag(im_part).eq.1).or. & ! prescribed solid (EUL)
            (FSI_flag(im_part).eq.4)) then ! CTML FSI
    FSI_partid_map(part_id)=0
   else
    print *,"FSI_flag(im_part) invalid"
    stop
   endif

   FSI(part_id)%flag_2D_to_3D=0

  enddo ! part_id=1,TOTAL_NPARTS

  CTML_NPARTS=ctml_part_id
  FSI_NPARTS=fsi_part_id
  if (FSI_NPARTS+CTML_NPARTS.gt.TOTAL_NPARTS) then
   print *,"FSI_NPARTS+CTML_NPARTS.gt.TOTAL_NPARTS"
   stop
  endif

  if (h_small.le.zero) then
   print *,"h_small invalid"
   stop
  endif
  do dir=1,3
   problo_act(dir)=problo(dir)
   probhi_act(dir)=probhi(dir)
   problen_act(dir)=probhi(dir)-problo(dir)
   if (problen_act(dir).le.zero) then
    print *,"problen_act(dir).le.zero"
    stop
   endif
  enddo ! dir=1..3
  do dir=1,AMREX_SPACEDIM
   problo_ref(dir)=problo(dir)
   probhi_ref(dir)=probhi(dir)
   problen_ref(dir)=probhi(dir)-problo(dir)
   if (problen_ref(dir).le.zero) then
    print *,"problen_ref(dir).le.zero"
    stop
   endif
   if (dx_max_level(dir).le.zero) then
    print *,"dx_max_level(dir).le.zero"
    stop
   endif
  enddo ! dir=1..AMREX_SPACEDIM

  use_temp=0

  if (CTML_FSI_flagF(nmat).eq.1) then ! FSI_flag==4
#ifdef MVAHABFSI
   if (CTML_FSI_INIT.eq.0) then
    call CTML_INIT_SOLID(dx_max_level, &
     problo_ref,probhi_ref,ioproc, &
     ctml_n_fib_bodies,ctml_max_n_fib_nodes)
    if (ctml_n_fib_bodies.eq.CTML_NPARTS) then
     allocate(ctml_n_fib_nodes(ctml_n_fib_bodies))
     allocate(ctml_n_fib_active_nodes(ctml_n_fib_bodies))
     call CTML_GET_FIB_NODE_COUNT( &
      ctml_n_fib_bodies, &
      ctml_n_fib_nodes)
     allocate(ctml_fib_pst(ctml_n_fib_bodies,ctml_max_n_fib_nodes,AMREX_SPACEDIM))
     allocate(ctml_fib_vel(ctml_n_fib_bodies,ctml_max_n_fib_nodes,AMREX_SPACEDIM))
     allocate(ctml_fib_frc(ctml_n_fib_bodies,ctml_max_n_fib_nodes, &
       2*AMREX_SPACEDIM))
     allocate(ctml_fib_mass(ctml_n_fib_bodies,ctml_max_n_fib_nodes))
    else
     print *,"ctml_n_fib_bodies out of range"
     stop
    endif
   else if (CTML_FSI_INIT.eq.1) then
    ! do nothing
   else
    print *,"CTML_FSI_INIT invalid"
    stop
   endif
   do part_id=1,TOTAL_NPARTS
    im_part=im_solid_mapF(part_id)+1
    if (CTML_FSI_mat(nmat,im_part).eq.1) then
     FSI(part_id)%deforming_part=1
    else if (CTML_FSI_mat(nmat,im_part).eq.0) then
     FSI(part_id)%deforming_part=0
    else
     print *,"CTML_FSI_mat(nmat,im_part) invalid"
     stop
    endif
   enddo ! part_id=1,TOTAL_NPARTS
   
#else
   print *,"define MVAHABFSI"
   stop
#endif
  else if (CTML_FSI_flagF(nmat).eq.0) then

   if((probtype.eq.538).or.(probtype.eq.541)) then ! needle, housing
    test_NPARTS=2
    FSI(1)%deforming_part=0
    FSI(2)%deforming_part=0
   else if (probtype.eq.701) then  ! flapping wing - ReadHeader
    if ((axis_dir.eq.0).or.(axis_dir.eq.2)) then
     test_NPARTS=1
     FSI(1)%deforming_part=0
    else if (axis_dir.eq.1) then
     test_NPARTS=2
     FSI(1)%deforming_part=0
     FSI(2)%deforming_part=0
    else
     print *,"axis_dir invalid"
     stop
    endif 
   else if ((probtype.eq.57).or. &
            (probtype.eq.562).or. & ! whale
            (probtype.eq.561).or. & 
            (probtype.eq.563).or. & 
            (probtype.eq.536).or. & 
            (probtype.eq.537).or. & 
            (probtype.eq.539).or. & 
            (probtype.eq.52).or. & 
            (probtype.eq.56).or. & 
            (probtype.eq.58).or. & 
            (probtype.eq.50).or. &  ! paddle
            (probtype.eq.5600).or. &  ! dog
            (probtype.eq.5601).or. &  ! viorel sphere
            (probtype.eq.5602).or. &  ! internal inflow
            (probtype.eq.5700).or. &  ! microfluidics
            (probtype.eq.53).or. &  
            (probtype.eq.531).or. &  
            (probtype.eq.5501)) then
    test_NPARTS=1
    FSI(part_id)%deforming_part=1
   else if (probtype.eq.9) then ! ship
    test_NPARTS=1
    FSI(part_id)%deforming_part=0
   else
    test_NPARTS=TOTAL_NPARTS
    do part_id=1,TOTAL_NPARTS
     FSI(part_id)%deforming_part=0
    enddo
   endif
   if (test_NPARTS.ne.TOTAL_NPARTS) then
    print *,"test_NPARTS.ne.TOTAL_NPARTS"
    stop
   endif
  else
   print *,"CTML_FSI_flagF invalid"
   stop
  endif 

  if ((TOTAL_NPARTS.ge.1).and. &
      (TOTAL_NPARTS.le.MAX_PARTS)) then

   do part_id=1,TOTAL_NPARTS

    im_part=im_solid_mapF(part_id)+1

    if ((FSI_partid_map(part_id).gt.0).or. &
        (CTML_partid_map(part_id).gt.0)) then
 
     FSI(part_id)%PartID=part_id

     if ((FSI_refine_factor(im_part).ge.0).and. &
         (FSI_refine_factor(im_part).le.100)) then
      FSI(part_id)%refine_factor=FSI_refine_factor(im_part)
     else
      print *,"FSI_refine_factor(im_part) invalid"
      stop
     endif 
     if (FSI_bounding_box_ngrow(im_part).eq.3) then
      FSI(part_id)%bounding_box_ngrow=FSI_bounding_box_ngrow(im_part)
     else
      print *,"FSI_bounding_box_ngrow(im_part) invalid"
      stop
     endif

     ! ifirst=1 (=>read from data files in 
     !           initinjector, initflapping, init_gingerbread2D,
     !           init_helix, or init_from_cas)
     call overall_solid_init(CLSVOFtime,ioproc,part_id,isout)  

     ! ReadHeader
     initflag=1
     call generate_new_triangles(initflag,problo,probhi, &
      part_id,ioproc,isout,h_small)
    else if ((FSI_partid_map(part_id).eq.0).and. &
             (CTML_partid_map(part_id).eq.0)) then
     ! do nothing
    else
     print *,"FSI_partid_map or CTML_partid_map invalid"
     stop
    endif

   enddo ! part_id=1..TOTAL_NPARTS
  else
   print *,"TOTAL_NPARTS invalid: ",TOTAL_NPARTS
   stop
  endif
 
return
end subroutine CLSVOF_ReadHeader

INTEGER_T function sign_valid(mask)
IMPLICIT NONE

INTEGER_T mask

if ((mask.eq.2).or. &  !singly wetted, sign and velocity init
    (mask.eq.3).or. &  !doubly wetted, sign and velocity init
    (mask.eq.103).or. &!doubly wetted, sign, value and velocity init
    (mask.eq.10).or. & !sign and velocity init from coarse level or prev step.
    (mask.eq.11)) then !sign from coarse level or prev step; velocity init.
 sign_valid=1
else if ((mask.eq.0).or.(mask.eq.1)) then
 sign_valid=0
else
 print *,"mask invalid"
 stop
endif

end function sign_valid


INTEGER_T function vel_valid(mask)
IMPLICIT NONE

INTEGER_T mask

if ((mask.eq.1).or. &  !velocity ok, sign not ok.
    (mask.eq.2).or. &  !singly wetted, sign and velocity init
    (mask.eq.3).or. &  !doubly wetted, sign and velocity init
    (mask.eq.103).or. &!doubly wetted, sign, value and velocity init
    (mask.eq.10).or. & !sign and velocity init from coarse level or prev step.
    (mask.eq.11)) then !sign from coarse level or prev step; velocity init.
 vel_valid=1
else if (mask.eq.0) then
 vel_valid=0
else
 print *,"mask invalid"
 stop
endif

end function vel_valid


subroutine check_overlap(part_id,ielem,time,minnode,maxnode, &
 tid,tilenum,dx3D,lev77,interior_flag,overlap,isweep)
IMPLICIT NONE

INTEGER_T part_id
INTEGER_T tid,tilenum
INTEGER_T ielem
INTEGER_T lev77
INTEGER_T interior_flag
INTEGER_T overlap
INTEGER_T isweep
INTEGER_T local_nelems
REAL_T time
REAL_T dx3D(3)
REAL_T minnode(3)
REAL_T maxnode(3)
INTEGER_T dir
INTEGER_T tilelo,tilehi
REAL_T xlo,xhi
INTEGER_T local_iband
REAL_T local_max_side
REAL_T local_buffer(3)

 if (lev77.lt.1) then
  print *,"lev77 invalid"
  stop
 endif
 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 local_iband=FSI(part_id)%bounding_box_ngrow
 local_max_side=FSI(part_id)%max_side_len_refined
 if ((local_iband.eq.3).and.(local_max_side.gt.zero)) then
  do dir=1,3
   local_buffer(dir)=local_iband*dx3D(dir)+local_max_side
  enddo
 else
  print *,"local_iband or local_max_side invalid"
  stop
 endif

 if (ielem.le.0) then
  print *,"ielem invalid"
  stop
 endif

 if ((tid.eq.0).or.(tilenum.eq.0)) then

  interior_flag=0
  overlap=0

 else if ((tid.ge.1).and.(tilenum.ge.1)) then

  if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.eq.0) then

   print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc=0"
   stop

  else if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.gt.0) then

   if (tilenum.gt.contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)) then
    print *,"tilenum.gt.contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)"
    stop
   endif

  else
   print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc bad"
   stop
  endif

  overlap=1
  interior_flag=1
  do dir=1,3
   xlo=contain_elem(lev77)%xlo3D(tid,tilenum,dir) 
   tilelo=contain_elem(lev77)%tilelo3D(tid,tilenum,dir) 
   tilehi=contain_elem(lev77)%tilehi3D(tid,tilenum,dir) 
   xhi=xlo+dx3D(dir)*(tilehi-tilelo+1)
   if (xhi.le.xlo) then
    print *,"xhi.le.xlo"
    stop
   endif
   if (maxnode(dir).lt.xlo-local_buffer(dir)) then
    overlap=0
   endif
   if (minnode(dir).gt.xhi+local_buffer(dir)) then
    overlap=0
   endif
   if ((minnode(dir).lt.xlo+local_buffer(dir)).or. &
       (maxnode(dir).gt.xhi-local_buffer(dir))) then
    interior_flag=0
   endif
   if (minnode(dir).gt.maxnode(dir)) then
    print *,"minnode(dir).gt.maxnode(dir)"
    stop
   endif
  enddo ! dir=1..3
  if (overlap.eq.1) then
   local_nelems=contain_elem(lev77)% &
                level_elem_data(tid,part_id,tilenum)%numElems
   local_nelems=local_nelems+1
   contain_elem(lev77)%level_elem_data(tid,part_id,tilenum)%numElems= &
      local_nelems
   if (isweep.eq.1) then
    ! do nothing
   else if (isweep.eq.2) then
    contain_elem(lev77)%level_elem_data(tid,part_id,tilenum)% &
      ElemData(local_nelems)=ielem
   else
    print *,"isweep invalid"
    stop
   endif
  else if (overlap.eq.0) then
   if (interior_flag.ne.0) then
    print *,"interior_flag.ne.0"
    stop
   endif
  else
   print *,"overlap invalid"
   stop
  endif

 else
  print *,"tid or tilenum invalid"
  stop
 endif 

return
end subroutine check_overlap

subroutine check_overlap_node(part_id,inode,time, &
 minnode, &
 tid,tilenum,dx3D,lev77,overlap)
IMPLICIT NONE

INTEGER_T part_id
INTEGER_T tid,tilenum
INTEGER_T inode
INTEGER_T lev77
INTEGER_T overlap
INTEGER_T local_nnodes
REAL_T time
REAL_T dx3D(3)
REAL_T minnode(3)
INTEGER_T dir
INTEGER_T tilelo,tilehi
REAL_T xlo,xhi

 if (lev77.lt.1) then
  print *,"lev77 invalid"
  stop
 endif
 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif

 if ((inode.ge.1).and. &
     (inode.le.FSI(part_id)%NumNodes)) then

  if ((tid.eq.0).or.(tilenum.eq.0)) then

   overlap=0

  else if ((tid.ge.1).and.(tilenum.ge.1)) then

   if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.eq.0) then

    print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc=0"
    stop

   else if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.gt.0) then

    if (tilenum.gt.contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)) then
     print *,"tilenum.gt.contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)"
     stop
    endif

   else
    print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc bad"
    stop
   endif

   overlap=1
   do dir=1,3
    xlo=contain_elem(lev77)%xlo3D(tid,tilenum,dir) 
    tilelo=contain_elem(lev77)%tilelo3D(tid,tilenum,dir) 
    tilehi=contain_elem(lev77)%tilehi3D(tid,tilenum,dir) 
    xhi=xlo+dx3D(dir)*(tilehi-tilelo+1)
    if (xhi.le.xlo) then
     print *,"xhi.le.xlo"
     stop
    endif
    if (minnode(dir).lt.xlo) then
     overlap=0
    endif
    if (minnode(dir).gt.xhi) then
     overlap=0
    endif
    if (1.eq.0) then
     print *,"dir,inode,xlo,xhi,tilelo,tilehi,xnode ", &
      dir,inode,xlo,xhi,tilelo,tilehi,minnode(dir)
    endif
   enddo ! dir=1..3
   if (overlap.eq.1) then
    local_nnodes=contain_elem(lev77)% &
                 level_node_data(tid,part_id,tilenum)%numNodes
    local_nnodes=local_nnodes+1
    contain_elem(lev77)%level_node_data(tid,part_id,tilenum)%numNodes= &
       local_nnodes
   else if (overlap.eq.0) then
    ! do nothing
   else
    print *,"overlap invalid"
    stop
   endif

  else
   print *,"tid or tilenum invalid"
   stop
  endif 

 else
  print *,"inode invalid"
  stop
 endif

return
end subroutine check_overlap_node


subroutine get_minmax_node(part_id,ielem,time,minnode,maxnode)
IMPLICIT NONE

INTEGER_T part_id
INTEGER_T ielem
REAL_T time
REAL_T minnode(3)
REAL_T maxnode(3)
INTEGER_T inode,dir,node_id
REAL_T nodetest
REAL_T xtarget(3)
REAL_T xfoot(3)
REAL_T velparm(3)

 if ((ielem.ge.1).and. &
     (ielem.le.FSI(part_id)%NumIntElemsBIG)) then
  do inode=1,3
   do dir=1,3
    node_id=FSI(part_id)%IntElemBIG(inode,ielem)
    if ((node_id.ge.1).and. &
        (node_id.le.FSI(part_id)%NumNodesBIG)) then
     xfoot(dir)=FSI(part_id)%NodeBIG(dir,node_id)
     if (abs(xfoot(dir)).le.CTMLoverflow) then
      ! do nothing
     else
      print *,"xfoot(dir) out of range:",xfoot(dir)
      print *,"dir=",dir
      print *,"node_id=",node_id
      print *,"inode=",inode 
      print *,"ielem=",ielem 
      stop
     endif 
    else
     print *,"node_id invalid"
     stop
    endif
    velparm(dir)=zero
   enddo ! dir=1..3
   call get_target_from_foot(xfoot,xtarget, &
    velparm,time,part_id)
   do dir=1,3
    nodetest=xtarget(dir)
    if (abs(nodetest).le.CTMLoverflow) then
     if (inode.eq.1) then
      minnode(dir)=nodetest
      maxnode(dir)=nodetest
     endif
     if (nodetest.lt.minnode(dir)) then
      minnode(dir)=nodetest
     endif
     if (nodetest.gt.maxnode(dir)) then
      maxnode(dir)=nodetest
     endif
    else
     print *,"nodetest out of range"
     stop
    endif
   enddo ! dir=1..3
  enddo  ! inode=1,3
 else
  print *,"ielem invalid"
  stop
 endif

return
end subroutine get_minmax_node


subroutine get_contained_node(part_id,inode,time,minnode)
IMPLICIT NONE

INTEGER_T part_id
INTEGER_T inode
REAL_T time
REAL_T minnode(3)
INTEGER_T dir
REAL_T xtarget(3)
REAL_T xfoot(3)
REAL_T velparm(3)


 if ((inode.ge.1).and. &
     (inode.le.FSI(part_id)%NumNodes)) then

  do dir=1,3
   xfoot(dir)=FSI(part_id)%Node(dir,inode)
   velparm(dir)=zero
  enddo
  call get_target_from_foot(xfoot,xtarget, &
   velparm,time,part_id)
  do dir=1,3
   minnode(dir)=xtarget(dir)
  enddo ! dir=1..3
 else
  print *,"inode invalid"
  stop
 endif

return
end subroutine get_contained_node

! elem_contain_type, node_contain_type,
! level_contain_type, contain_elem
! are declared in PROBCOMMON.F90
subroutine CLSVOF_FILLCONTAINER( &
 lev77, &
 sci_max_level, &
 nthread_parm, &
 dx3D, &
 part_id, &
 im_part, &
 nmat, &
 time, &
 dt)
use global_utility_module

IMPLICIT NONE

 INTEGER_T lev77
 INTEGER_T sci_max_level
 INTEGER_T nthread_parm
 REAL_T dx3D(3)
 INTEGER_T part_id
 INTEGER_T im_part
 INTEGER_T nmat
 REAL_T time
 REAL_T dt

 INTEGER_T interior_flag
 INTEGER_T overlap
 INTEGER_T cache_saved

 INTEGER_T ielem,inode,isweep
 INTEGER_T local_nelems,local_nnodes
 INTEGER_T num_elements,num_nodes
 INTEGER_T total_num_elements,total_num_nodes
 INTEGER_T total_num_elements_check,total_num_nodes_check

 INTEGER_T tid,tid_loop,tilenum,tilenum_loop
 INTEGER_T tid_predict,tilenum_predict

 REAL_T maxnode(3),minnode(3)

 INTEGER_T, dimension(:), allocatable :: tid_elem
 INTEGER_T, dimension(:), allocatable :: tilenum_elem
 INTEGER_T, dimension(:), allocatable :: interior_elem

 INTEGER_T, dimension(:), allocatable :: tid_node
 INTEGER_T, dimension(:), allocatable :: tilenum_node

 INTEGER_T ctml_part_id
 INTEGER_T fsi_part_id

 if ((lev77.lt.1).or.(lev77.gt.sci_max_level+1)) then
  print *,"lev77 invalid"
  stop
 endif
 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if ((part_id.lt.1).or.(part_id.ge.nmat)) then
  print *,"part_id invalid"
  stop
 endif
 if ((im_part.lt.1).or.(im_part.gt.nmat)) then
  print *,"im_part invalid"
  stop
 endif
 if (nmat.ne.num_materials) then
  print *,"nmat invalid"
  stop
 endif 
 if (nthread_parm.lt.1) then
  print *,"nthread_parm.lt.1"
  stop
 endif

 num_elements=FSI(part_id)%NumIntElemsBIG
 if (num_elements.le.0) then
  print *,"num_elements invalid"
  stop
 endif
 num_nodes=FSI(part_id)%NumNodes
 if (num_nodes.le.0) then
  print *,"num_nodes invalid"
  stop
 endif

 ctml_part_id=CTML_partid_map(part_id)
 fsi_part_id=FSI_partid_map(part_id)

 if (((ctml_part_id.ge.1).and. &
      (ctml_part_id.le.CTML_NPARTS)).or. &
     ((fsi_part_id.ge.1).and. &
      (fsi_part_id.le.FSI_NPARTS))) then

   ! store cache information for predicting the thread id and tilenum
   ! that the next element will lie in.
   ! If interior_elem==1 then the element can only belong to one 
   ! tile (tid,tilenum).
  allocate(tid_elem(num_elements))
  allocate(tilenum_elem(num_elements))
  allocate(interior_elem(num_elements))

  do ielem=1,num_elements
   tid_elem(ielem)=0
   tilenum_elem(ielem)=0
   interior_elem(ielem)=0
  enddo

   ! store cache information for predicting the thread id and tilenum
   ! that the next node will lie in.
  allocate(tid_node(num_nodes))
  allocate(tilenum_node(num_nodes))

  do inode=1,num_nodes
   tid_node(inode)=0
   tilenum_node(inode)=0
  enddo

  do tid=1,nthread_parm
   if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.eq.0) then
    ! do nothing
   else if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.gt.0) then
    if (contain_elem(lev77)%num_tiles_on_thread3D_proc(tid).gt. &
        contain_elem(lev77)%max_num_tiles_on_thread3D_proc) then
     print *,"contain_elem(lev77):ntiles>max_ntiles"
     stop
    endif
    do tilenum=1,contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)
     contain_elem(lev77)%level_elem_data(tid,part_id,tilenum)%numElems=0
    enddo
    do tilenum=1,contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)
     contain_elem(lev77)%level_node_data(tid,part_id,tilenum)%numNodes=0
    enddo
   else
    print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc bad"
    stop
   endif
  enddo ! tid=1,nthread_parm

  if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.eq.0) then
   ! do nothing
  else if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.gt.0) then

   do isweep=1,2

    if (isweep.eq.2) then

     total_num_nodes=0
     total_num_elements=0

     do tid=1,nthread_parm

      do tilenum=1,contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)

       local_nelems=contain_elem(lev77)% &
                   level_elem_data(tid,part_id,tilenum)%numElems

       total_num_elements=total_num_elements+local_nelems

       if (local_nelems.eq.0) then
        ! do nothing
       else if ((local_nelems.ge.1).and.(local_nelems.le.num_elements)) then
        allocate(contain_elem(lev77)%level_elem_data(tid,part_id,tilenum)% &
                ElemData(local_nelems))
       else
        print *,"local_nelems invalid"
        stop
       endif
       contain_elem(lev77)%level_elem_data(tid,part_id,tilenum)%numElems=0

       local_nnodes=contain_elem(lev77)% &
                   level_node_data(tid,part_id,tilenum)%numNodes

       total_num_nodes=total_num_nodes+local_nnodes

       if (local_nnodes.eq.0) then
        ! do nothing
       else if ((local_nnodes.ge.1).and.(local_nnodes.le.num_nodes)) then
        allocate(contain_elem(lev77)%level_node_data(tid,part_id,tilenum)% &
                NodeData(local_nnodes))
       else
        print *,"local_nnodes invalid"
        stop
       endif
       contain_elem(lev77)%level_node_data(tid,part_id,tilenum)%numNodes=0

      enddo ! tilenum

     enddo ! tid=1,nthread_parm

    else if (isweep.eq.1) then
     ! do nothing
    else
     print *,"isweep invalid"
     stop
    endif

    do ielem=1,num_elements

     if (isweep.eq.2) then

      tid=tid_elem(ielem)
      tilenum=tilenum_elem(ielem)
      interior_flag=interior_elem(ielem)
      if ((tid.eq.0).and.(tilenum.eq.0)) then
       ! do nothing (element outside of the present AMR level)
      else if ((tid.ge.1).and.(tid.le.nthread_parm).and. &
               (tilenum.ge.1)) then
       if (tilenum.gt.contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)) then
        print *,"tilenum out of range"
        stop
       endif
       local_nelems=contain_elem(lev77)% &
                    level_elem_data(tid,part_id,tilenum)%numElems
       local_nelems=local_nelems+1
       contain_elem(lev77)%level_elem_data(tid,part_id,tilenum)%numElems= &
          local_nelems
       contain_elem(lev77)%level_elem_data(tid,part_id,tilenum)% &
                           ElemData(local_nelems)=ielem
       if (interior_flag.eq.1) then
        ! do nothing
       else if (interior_flag.eq.0) then
        call get_minmax_node(part_id,ielem,time,minnode,maxnode) 
        do tid_loop=1,nthread_parm
        do tilenum_loop=1, &
               contain_elem(lev77)%num_tiles_on_thread3D_proc(tid_loop)
         if ((tid_loop.ne.tid).or.(tilenum_loop.ne.tilenum)) then
           ! isweep==2
          call check_overlap(part_id,ielem,time,minnode,maxnode, &
            tid_loop,tilenum_loop,dx3D,lev77,interior_flag,overlap,isweep)
         endif
        enddo
        enddo
       else
        print *,"interior_flag invalid"
        stop
       endif
      else
       print *,"tid or tilenum invalid"
       stop
      endif 

     else if (isweep.eq.1) then

      call get_minmax_node(part_id,ielem,time,minnode,maxnode)
     
      if (ielem.eq.1) then
       tid_predict=0
       tilenum_predict=0
      else if ((ielem.gt.1).and.(ielem.le.num_elements)) then
       tid_predict=tid_elem(ielem-1)
       tilenum_predict=tilenum_elem(ielem-1)
      else
       print *,"ielem invalid"
       stop
      endif

       ! isweep==1
       ! check_overlap increments
       !  contain_elem(lev77)%level_elem_data(tid,part_id,tilenum)%numElems
      call check_overlap(part_id,ielem,time,minnode,maxnode, &
        tid_predict,tilenum_predict, &
        dx3D,lev77,interior_flag,overlap,isweep)
      cache_saved=0
      if (overlap.eq.1) then
       cache_saved=1
       tid_elem(ielem)=tid_predict
       tilenum_elem(ielem)=tilenum_predict
       interior_elem(ielem)=interior_flag
      else if (overlap.eq.0) then
       ! do nothing
      else
       print *,"overlap invalid"
       stop
      endif

      if (interior_flag.eq.0) then

       do tid_loop=1,nthread_parm
       do tilenum_loop=1, &
               contain_elem(lev77)%num_tiles_on_thread3D_proc(tid_loop)
        if ((tid_loop.ne.tid_predict).or. &
            (tilenum_loop.ne.tilenum_predict)) then
          ! isweep==1
         call check_overlap(part_id,ielem,time,minnode,maxnode, &
          tid_loop,tilenum_loop,dx3D,lev77,interior_flag,overlap,isweep)
         if (overlap.eq.1) then
          if ((cache_saved.eq.0).or.(interior_flag.eq.1)) then
           tid_elem(ielem)=tid_loop
           tilenum_elem(ielem)=tilenum_loop
           interior_elem(ielem)=interior_flag
           cache_saved=1
          else if ((cache_saved.eq.1).and.(interior_flag.eq.0)) then
           ! do nothing
          else
           print *,"cache_saved or interior_flag invalid"
           stop
          endif
         else if (overlap.eq.0) then
          ! do nothing
         else
          print *,"overlap invalid"
          stop
         endif
        else if ((tid_loop.eq.tid_predict).and. &
                 (tilenum_loop.eq.tilenum_predict)) then
         ! do nothing
        else
         print *,"tid_loop or tilenum_loop invalid"
         stop
        endif
            
       enddo ! tilenum_loop
       enddo ! tid_loop

      else if (interior_flag.eq.1) then

       ! do nothing

      else
       print *,"interior_flag invalid"
       stop
      endif

     else
      print *,"isweep invalid"
      stop
     endif

    enddo ! ielem=1,num_elements

    do inode=1,num_nodes

     if (isweep.eq.2) then

      tid=tid_node(inode)
      tilenum=tilenum_node(inode)
      if ((tid.eq.0).and.(tilenum.eq.0)) then
       ! do nothing (node is not on the present AMR level)
      else if ((tid.ge.1).and.(tid.le.nthread_parm).and. &
               (tilenum.ge.1)) then
       if (tilenum.gt.contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)) then
        print *,"tilenum out of range"
        stop
       endif
       local_nnodes=contain_elem(lev77)% &
                    level_node_data(tid,part_id,tilenum)%numNodes
       local_nnodes=local_nnodes+1
       contain_elem(lev77)%level_node_data(tid,part_id,tilenum)%numNodes= &
          local_nnodes
       contain_elem(lev77)%level_node_data(tid,part_id,tilenum)% &
                           NodeData(local_nnodes)=inode
      else
       print *,"tid or tilenum invalid"
       stop
      endif 

     else if (isweep.eq.1) then

      if (inode.eq.1) then
       tid_predict=0
       tilenum_predict=0
      else if ((inode.gt.1).and.(inode.le.num_nodes)) then
       tid_predict=tid_node(inode-1)
       tilenum_predict=tilenum_node(inode-1)
      else
       print *,"inode invalid"
       stop
      endif

      call get_contained_node(part_id,inode,time,minnode)

       ! isweep==1
       ! check_overlap_node increments
       !  contain_elem(lev77)%level_node_data(tid,part_id,tilenum)%numNodes
      call check_overlap_node(part_id,inode,time, &
        minnode, &
        tid_predict,tilenum_predict, &
        dx3D,lev77,overlap)

      if (overlap.eq.1) then

       tid_node(inode)=tid_predict
       tilenum_node(inode)=tilenum_predict

      else if (overlap.eq.0) then

       do tid_loop=1,nthread_parm
       do tilenum_loop=1, &
               contain_elem(lev77)%num_tiles_on_thread3D_proc(tid_loop)

        if (overlap.eq.0) then
         if ((tid_loop.ne.tid_predict).or. &
             (tilenum_loop.ne.tilenum_predict)) then
           ! isweep==1
          call check_overlap_node(part_id,inode,time, &
           minnode, &
           tid_loop,tilenum_loop, &
           dx3D,lev77,overlap)
          if (overlap.eq.1) then
           tid_node(inode)=tid_loop
           tilenum_node(inode)=tilenum_loop
          else if (overlap.eq.0) then
           ! do nothing
          else
           print *,"overlap invalid"
           stop
          endif
         else if ((tid_loop.eq.tid_predict).and. &
                  (tilenum_loop.eq.tilenum_predict)) then
          ! do nothing
         else
          print *,"tid_loop or tilenum_loop invalid"
          stop
         endif
        else if (overlap.eq.1) then
         ! do nothing, node already added, do not want double counting,
         ! algorithm needs to be threead safe.
        else
         print *,"overlap invalid"
         stop
        endif
            
       enddo ! tilenum_loop
       enddo ! tid_loop

      else
       print *,"overlap invalid"
       stop
      endif

     else
      print *,"isweep invalid"
      stop
     endif

    enddo ! inode=1,num_nodes
  
    if (isweep.eq.1) then
     ! do nothing
    else if (isweep.eq.2) then

     total_num_nodes_check=0
     total_num_elements_check=0

     do tid=1,nthread_parm
      do tilenum=1,contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)

       local_nelems=contain_elem(lev77)% &
                    level_elem_data(tid,part_id,tilenum)%numElems

       total_num_elements_check=total_num_elements_check+local_nelems


       local_nnodes=contain_elem(lev77)% &
                    level_node_data(tid,part_id,tilenum)%numNodes

       total_num_nodes_check=total_num_nodes_check+local_nnodes

      enddo ! tilenum
     enddo ! tid

     if (total_num_nodes_check.ne.total_num_nodes) then
      print *,"total_num_nodes invalid"
      stop
     endif
     if (total_num_elements_check.ne.total_num_elements) then
      print *,"total_num_elements invalid"
      stop
     endif

     if (1.eq.1) then
      print *,"CLSVOF_FILLCONTAINER"
      print *,"lev77=",lev77
      print *,"part_id=",part_id
      print *,"total_num_elements ",total_num_elements 
      print *,"total_num_nodes ",total_num_nodes
     endif

    else 
     print *,"isweep invalid"
     stop
    endif
 
   enddo ! isweep=1,2

  else
   print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc bad"
   stop
  endif

  deallocate(tid_elem)
  deallocate(tilenum_elem)
  deallocate(interior_elem)

  deallocate(tid_node)
  deallocate(tilenum_node)

 else if ((ctml_part_id.eq.0).and. &
          (fsi_part_id.eq.0)) then
  ! do nothing
 else
  print *,"ctml_part_id or fsi_part_id invalid"
  stop
 endif

return
end subroutine CLSVOF_FILLCONTAINER

! nFSI==nparts*nFSI_sub
! nparts x (vel + LS + temperature + flag + stress)
! mask=0 prior to entry
! mask=1 velocity is init from fine lev, but sign is not
! mask=2 both sign and velocity are init on fine lev
! mask=3 both sign and velocity are init on fine lev (doubly wetted)
! mask=0 neither velocity, sign, or LS are valid.  
! mask=10 sign is init from coarse level
! mask=11 velocity is init from fine lev, sign is init from coarse level

! vel=temp=force=0 if no interfaces in cell's neighborhood
subroutine CLSVOF_InitBox(  &
  iter, &
  sdim_AMR, &
  lev77, &
  tid, &
  tilenum, &
  im_part, &
  nparts, &
  part_id, &
  ngrowFSI, &
  nmat, &
  nFSI, &
  nFSI_sub, &
  FSI_operation, &
  touch_flag, &
  h_small, &
  time,dt, &
  problo3D,probhi3D, &
  xmap3D, &
  xslice3D, &
  dx3D, &
  FSI_lo,FSI_hi, &
  FSI_growlo,FSI_growhi, &
  growlo3D,growhi3D, &
  xdata3D, &
  FSIdata3D, &
  masknbr3D, &
  DIMS3D(FSIdata3D), &
  CTML_force_model, &
  ioproc,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

  INTEGER_T :: iter
  INTEGER_T :: sdim_AMR
  INTEGER_T :: lev77
  INTEGER_T :: tid
  INTEGER_T :: tilenum
  INTEGER_T :: im_part
  INTEGER_T :: nparts
  INTEGER_T :: part_id
  INTEGER_T :: ngrowFSI
  INTEGER_T :: nmat
  INTEGER_T :: nFSI
  INTEGER_T :: nFSI_sub
  INTEGER_T :: FSI_operation
  INTEGER_T :: touch_flag
  INTEGER_T :: numtouch
  REAL_T :: h_small
  REAL_T :: time,dt
  REAL_T problo3D(3)
  REAL_T probhi3D(3)
  INTEGER_T xmap3D(3)
  REAL_T xslice3D(3)
  REAL_T dx3D(3)
  INTEGER_T FSI_lo(3),FSI_hi(3)
  INTEGER_T FSI_growlo(3),FSI_growhi(3)
  INTEGER_T growlo3D(3),growhi3D(3)
  INTEGER_T DIMDEC3D(FSIdata3D)
  REAL_T xdata3D(DIMV3D(FSIdata3D),3)
  REAL_T FSIdata3D(DIMV3D(FSIdata3D),nFSI)
  REAL_T masknbr3D(DIMV3D(FSIdata3D),2)
  REAL_T, dimension(:,:,:,:), allocatable :: old_FSIdata

  INTEGER_T CTML_force_model
  INTEGER_T ioproc,isout

  REAL_T dxBB(3) ! set in find_grid_bounding_box

  INTEGER_T :: ielem
  INTEGER_T :: ielem_container
  INTEGER_T :: nodes_per_elem,inode,nodeptr
  REAL_T, dimension(3) :: xc,xclosest
  REAL_T, dimension(3) :: normal
  REAL_T, dimension(3) :: normal_closest
  REAL_T, dimension(3) :: xnot
  REAL_T, dimension(3) :: xfoot
  REAL_T, dimension(3) :: xelem
  REAL_T, dimension(3) :: xtarget
  INTEGER_T, dimension(3) :: gridlo,gridhi,gridlen
  INTEGER_T :: i,j,k
  INTEGER_T :: i1,j1,k1
  INTEGER_T :: inplane
  REAL_T :: wallthick
  REAL_T, dimension(3) :: velparm
  REAL_T, dimension(6) :: forceparm
  REAL_T :: massparm
  INTEGER_T :: dir
  REAL_T :: dotprod
  REAL_T :: unsigned_mindist  ! unsigned
  REAL_T :: weighttotal,distwt,weight
    ! (vel + LS + temperature + flag + stress)
  REAL_T, dimension(nFSI_sub) :: weight_top 
  REAL_T :: weight_bot
  REAL_T, dimension(3) :: minnode,maxnode
  REAL_T, dimension(3) :: xx
  REAL_T, dimension(3) :: xcen
  REAL_T, dimension(3) :: xleft
  REAL_T, dimension(3) :: xright
  REAL_T, dimension(3) :: xside,xcrit
  INTEGER_T :: ii,jj,kk
  INTEGER_T :: hitflag
  REAL_T :: phiside,phicenter,testdist,hitsign,totaldist

  INTEGER_T modify_vel
  REAL_T mag_n,mag_n_test,n_dot_x,signtest
  REAL_T mag_x
  INTEGER_T in_sign_box
  INTEGER_T ibase
  REAL_T ls_local
  INTEGER_T mask_local,mask_node
  INTEGER_T new_mask_local
  REAL_T, dimension(3) :: vel_local
  REAL_T, dimension(6) :: stress_local
  REAL_T temp_local
  INTEGER_T nc
  INTEGER_T sign_defined
  INTEGER_T sign_defined_local
  INTEGER_T i_norm
  INTEGER_T ii_current
  INTEGER_T mminus,mplus
  REAL_T xminus,xplus
  REAL_T LSsign,LSMINUS,LSPLUS
  REAL_T data_minus(nFSI_sub)
  REAL_T data_plus(nFSI_sub)
  INTEGER_T sign_status_changed
  INTEGER_T num_elements
  INTEGER_T num_elements_container
  REAL_T LSMIN_debug
  REAL_T LSMAX_debug
  REAL_T element_scale
  REAL_T test_scale
  INTEGER_T mask_debug
  INTEGER_T debug_all
  INTEGER_T in_the_interior
  INTEGER_T used_for_trial
  INTEGER_T mask1,mask2
  INTEGER_T dirmax,sweepmax
  INTEGER_T dir_order,sweep,incr
  INTEGER_T dirmap(3)
  INTEGER_T index_first(3)
  INTEGER_T index_last(3)
  INTEGER_T map_index(3)
  INTEGER_T inverse_index(3)
  INTEGER_T imap,jmap,kmap
  INTEGER_T null_intersection
  REAL_T force_weight
  REAL_T force_vector(3)
  INTEGER_T ctml_part_id
  INTEGER_T fsi_part_id
  INTEGER_T local_iband

  call FLUSH(6)  ! 6=screen

  debug_all=0

  if ((part_id.lt.1).or.(part_id.gt.nparts)) then
   print *,"part_id invalid"
   stop
  endif
  if (nparts.ne.TOTAL_NPARTS) then
   print *,"nparts.ne.TOTAL_NPARTS"
   stop
  endif

  if ((lev77.lt.1).or.(tid.lt.0).or.(tilenum.lt.0)) then
   print *,"lev77 or tid or tilenum invalid"
   stop
  endif
  if (container_allocated.ne.1) then
   print *,"container_allocated.ne.1"
   stop
  endif
  if (level_container_allocated(lev77).ne.1) then
   print *,"level_container_allocated(lev77).ne.1"
   stop
  endif

  if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.eq.0) then

   print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc=0"
   stop

  else if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.gt.0) then
  
   if (tilenum+1.gt.contain_elem(lev77)% &
       num_tiles_on_thread3D_proc(tid+1)) then
    print *,"tilenum+1.gt.num_tiles_on_thread3D_proc(tid+1)"
    stop
   endif

  else
   print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc bad"
   stop
  endif

  local_iband=FSI(part_id)%bounding_box_ngrow
  if (local_iband.eq.3) then
   ! do nothing
  else
   print *,"local_iband invalid"
   stop
  endif

  do dir=1,3
   if (contain_elem(lev77)%tilelo3D(tid+1,tilenum+1,dir).ne. &
       FSI_lo(dir)) then
    print *,"tilelo3D(tid+1,tilenum+1,dir).ne.FSI_lo(dir)"
    stop
   endif 
   if (contain_elem(lev77)%tilehi3D(tid+1,tilenum+1,dir).ne. &
       FSI_hi(dir)) then
    print *,"tilehi3D(tid+1,tilenum+1,dir).ne.FSI_hi(dir)"
    stop
   endif 
   if (FSI_growlo(dir).ge.FSI_lo(dir)) then
    print *,"FSI_growlo(dir) invalid"
    stop
   endif
   if (FSI_growhi(dir).le.FSI_hi(dir)) then
    print *,"FSI_growhi(dir) invalid"
    stop
   endif
  enddo ! dir=1..3

  ctml_part_id=CTML_partid_map(part_id)
  fsi_part_id=FSI_partid_map(part_id)

  if (((ctml_part_id.ge.1).and. &
       (ctml_part_id.le.CTML_NPARTS)).or. &
      ((fsi_part_id.ge.1).and. &
       (fsi_part_id.le.FSI_NPARTS))) then

   num_elements=FSI(part_id)%NumIntElemsBIG

   if (num_elements.le.0) then
    print *,"num_elements invalid"
    stop
   endif
   if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
   endif
   if ((sdim_AMR.ne.2).and.(sdim_AMR.ne.3)) then
    print *,"sdim_AMR invalid"
    stop
   endif
   if (iter.lt.0) then
    print *,"iter invalid"
    stop
   endif
   if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
    print *,"part_id invalid"
    stop
   endif
   if ((nmat.lt.1).or.(nmat.gt.50)) then
    print *,"nmat out of range"
    stop
   endif
   if ((im_part.lt.1).or.(im_part.gt.nmat)) then
    print *,"im_part invalid"
    stop
   endif
   if ((CTML_force_model.ne.0).and.(CTML_force_model.ne.1)) then
    print *,"CTML_force_model invalid"
    stop
   endif

   if (isout.eq.0) then
    ! do nothing
   else if (isout.eq.1) then
    print *,"in (START): CLSVOF_InitBox"
    print *,"num_elements(NumIntElemsBIG)=",num_elements
    print *,"NumIntElems=",FSI(part_id)%NumIntElems
    print *,"sdim_AMR=",sdim_AMR 
    print *,"im_part=",im_part
    print *,"part_id=",part_id
    print *,"nmat=",nmat
    print *,"FSI_operation=",FSI_operation
    print *,"iter= ",iter
    print *,"touch_flag=",touch_flag
    print *,"h_small=",h_small
    print *,"time,dt ",time,dt
    print *,"ioproc=",ioproc
    print '(A9,3(f12.6))',"problo3D ",problo3D(1),problo3D(2),problo3D(3)
    print '(A9,3(f12.6))',"probhi3D ",probhi3D(1),probhi3D(2),probhi3D(3)
    print '(A7,3(I10))',"FSI_lo ",FSI_lo(1),FSI_lo(2),FSI_lo(3)
    print '(A7,3(I10))',"FSI_hi ",FSI_hi(1),FSI_hi(2),FSI_hi(3)
    print '(A11,3(I10))',"FSI_growlo ",FSI_growlo(1),FSI_growlo(2),FSI_growlo(3)
    print '(A11,3(I10))',"FSI_growhi ",FSI_growhi(1),FSI_growhi(2),FSI_growhi(3)
    print '(A11,3(I10))',"growlo3D ",growlo3D(1),growlo3D(2),growlo3D(3)
    print '(A11,3(I10))',"growhi3D ",growhi3D(1),growhi3D(2),growhi3D(3)
   else
    print *,"isout invalid1: ",isout
    stop
   endif

    ! nparts x (vel + LS + temperature + flag + stress)

   if (nFSI.ne.nparts*nFSI_sub) then
    print *,"nFSI invalid"
    stop
   endif
   if (nFSI_sub.ne.12) then
    print *,"nFSI_sub invalid"
    stop
   endif
   if ((nparts.lt.1).or.(nparts.ge.nmat)) then
    print *,"nparts invalid"
    stop
   endif

   ibase=(part_id-1)*nFSI_sub

   if ((touch_flag.ne.0).and.(touch_flag.ne.1)) then
    print *,"touch_flag invalid"
    stop
   endif
   if ((FSI_operation.ne.2).and.(FSI_operation.ne.3)) then
    print *,"FSI_operation invalid"
    stop
   endif
   if (ngrowFSI.ne.3) then
    print *,"ngrowFSI invalid"
    stop
   endif

   call checkbound3D(FSI_lo,FSI_hi, &
    DIMS3D(FSIdata3D), &
    ngrowFSI,-1,521)

   mask_debug=0
   LSMIN_debug=1.0D+10
   LSMAX_debug=-1.0D+10

    ! in NavierStokes::initData ():
    !  state solid velocity for FSI_flag=2,4,6,7 is init to 0.0
    !  state level set function for FSI_flag=2,4,6,7 is init to -99999.0
    !
    ! in NavierStokes::FSI_make_distance FSI_MF is initialized as
    ! follows:
    ! nparts x (vel + LS + temperature + flag + stress)
    ! velocity=0.0
    ! LS=-99999
    ! temperature=0.0
    ! mask=0.0
    ! stress=0.0
    !
    ! then in NavierStokes::ns_header_msg_level:
    ! 1. fill coarse patch if level>0 and copy into level state variable
    ! 2. do get state and overwrite FSI_MF (velocity, LS, temperature)
    ! 3. mask=10.0 if (level>0) or ((level==0)and(time>0.0))
    !
    ! the order of operations on startup:
    ! 1. initData is called and LS=-99999.0, vel=0.0
    ! 2. FSI_make_distance is called from initData
    ! 3. The rest of the (non FSI_flag=2,4,6,7) materials are initialized.
   if (FSI_operation.eq.2) then ! make distance in narrow band

    num_elements_container=contain_elem(lev77)% &
                           level_elem_data(tid+1,part_id,tilenum+1)% &
                           numElems

    if (isout.eq.0) then
     ! do nothing
    else if (isout.eq.1) then
     print *,"lev77=",lev77
     print *,"tid= ",tid
     print *,"part_id= ",part_id
     print *,"tilenum= ",tilenum
     print *,"num_elements_container= ",num_elements_container
    else
     print *,"isout invalid"
     stop
    endif

    if ((num_elements_container.lt.0).or. &
        (num_elements_container.gt.num_elements)) then
     print *,"num_elements_container invalid"
     stop
    endif

    do ielem_container=1,num_elements_container

     ielem=contain_elem(lev77)% &
           level_elem_data(tid+1,part_id,tilenum+1)% &
           ElemData(ielem_container)

     if ((ielem.lt.1).or. &
         (ielem.gt.num_elements)) then
      print *,"ielem invalid"
      stop
     endif

     do dir=1,3
      xelem(dir)=FSI(part_id)%ElemDataXnotBIG(dir,ielem)
      velparm(dir)=zero
     enddo 
     call get_target_from_foot(xelem,xnot, &
       velparm,time,part_id)

     nodes_per_elem=FSI(part_id)%ElemDataBIG(1,ielem)
     if (nodes_per_elem.lt.3) then
      print *,"elem,nodes_per_elem ",ielem,nodes_per_elem   
      stop
     endif
     ! normal points from solid to fluid
     ! phi=n dot (x-xnot)
     ! phi>0 in the fluid
     call scinormal(ielem,normal,part_id,time)

     if (debug_all.eq.1) then
      print *,"ielem=",ielem
      do dir=1,3
       print *,"dir,xelem ",dir,xelem(dir)
       print *,"dir,xnot  ",dir,xnot(dir)
       print *,"dir,normal ",dir,normal(dir)
      enddo
     endif

      ! for each node in the element, this routine calls:
      ! get_target_from_foot
     call get_minmax_node(part_id,ielem,time,minnode,maxnode)
      ! sanity check
     do dir=1,3
      test_scale=maxnode(dir)-minnode(dir)
      if (test_scale.lt.zero) then
       print *,"test_scale.lt.zero"
       print *,"dir,minnode,maxnode,test_scale ",dir,minnode(dir), &
        maxnode(dir),test_scale
       print *,"part_id,ielem,time ",part_id,ielem,time
       stop
      endif
     enddo ! dir=1..3

     element_scale=FSI(part_id)%min_side_len_refined

     if (element_scale.le.zero) then
      print *,"element_scale.le.zero"
      print *,"part_id= ",part_id
      print *,"ielem= ",ielem
      print *,"time= ",time
      do dir=1,3
       print *,"dir,minnode,maxnode ",dir,minnode(dir),maxnode(dir)
      enddo
      stop
     endif
    
     call find_grid_bounding_box( &
       part_id, &
       null_intersection, &
       minnode,maxnode, &
       FSI_lo,FSI_hi, &
       FSI_growlo,FSI_growhi, &
       growlo3D,growhi3D, &
       xdata3D, &
       DIMS3D(FSIdata3D), &
       gridlo,gridhi,dxBB) 

     do dir=1,3
      if (abs(dxBB(dir)-dx3D(dir)).gt.element_buffer_tol*dxBB(dir)) then
       print *,"abs(dxBB(dir)-dx3D(dir)).gt.element_buffer_tol*dxBB(dir)"
       stop
      endif
     enddo ! dir=1..3

     if (1.eq.0) then
      print *,"ielem=",ielem
      do dir=1,3
       print *,"dir,gridlo,gridhi ",dir,gridlo(dir),gridhi(dir)
      enddo
     endif

     if (null_intersection.eq.0) then

      in_the_interior=1
      do dir=1,3
       gridlen(dir)=gridhi(dir)-gridlo(dir)+1
       if ((gridlen(dir).ge.0).and. &
           (gridlen(dir).lt.2*local_iband)) then
        in_the_interior=0
       else if ((gridlen(dir).ge.2*local_iband).and. &
                (gridlen(dir).le.2048)) then
        ! do nothing
       else
        print *,"gridlen(dir) invalid"
        stop
       endif
      enddo ! dir=1..3

      used_for_trial=0

       ! this code is thread safe
       ! gridlo,gridhi restricted to growlo3D and growhi3D 
      do i=gridlo(1),gridhi(1)
      do j=gridlo(2),gridhi(2)
      do k=gridlo(3),gridhi(3)

        ! BEFORE: restrict (i,j,k) to growlo3D and growhi3D
   
       if ((i.ge.FSI_lo(1)).and.(i.le.FSI_hi(1)).and. &
           (j.ge.FSI_lo(2)).and.(j.le.FSI_hi(2)).and. &
           (k.ge.FSI_lo(3)).and.(k.le.FSI_hi(3))) then

        mask1=NINT(masknbr3D(i,j,k,1))
        mask2=NINT(masknbr3D(i,j,k,2))

         ! mask2==1 if (i,j,k) interior to tile
         ! mask1==0 if (i,j,k) is exterior to tile and is a
         !  coarse/fine ghost cell.
        if ((mask1.eq.0).or.(mask2.eq.1)) then

         do dir=1,3
          xx(dir)=xdata3D(i,j,k,dir)
         enddo

         if (CTML_force_model.eq.0) then
          do inode=1,nodes_per_elem
            ! calls either CTML_DELTA or hsprime
           call check_force_weight(xmap3D,inode,ielem, &
            xx,part_id,time,dxBB,force_weight,force_vector)
           ! nparts x (vel + LS + temperature + flag + stress)
           do dir=1,3
            FSIdata3D(i,j,k,ibase+6+dir)=  &
              FSIdata3D(i,j,k,ibase+6+dir)+ &
              dt*force_weight*force_vector(dir)
           enddo 
          enddo ! inode=1,nodes_per_elem
         else if (CTML_force_model.eq.1) then
          ! do nothing
         else
          print *,"CTML_force_model invalid"
          stop
         endif
        
         ! normal points from solid to fluid
         dotprod=0.0 
         do dir=1,3
          dotprod=dotprod+normal(dir)*(xx(dir)-xnot(dir))
         enddo
         do dir=1,3
          xclosest(dir)=xx(dir)-dotprod*normal(dir)
          normal_closest(dir)=normal(dir)
         enddo
         unsigned_mindist=abs(dotprod)
         phicenter=dotprod ! phicenter>0 in the fluid

         if (debug_all.eq.1) then
          print *,"ielem=",ielem
          do dir=1,3
           print *,"dir,xx ",dir,xx(dir)
           print *,"dir,xclosest  ",dir,xclosest(dir)
          enddo
          if (dotprod.ge.zero) then
           print *,"dotprodplus=",dotprod
          else
           print *,"dotprodminus=",dotprod
          endif
         endif

         call checkinplane(xclosest,ielem,inplane, &
          minnode,maxnode,element_scale,part_id,time)

! investigate using NodeNormalBIG (normal defined at nodes)
         do inode=1,nodes_per_elem
          ! check distance to the edges of a triangular element.
          call checkinline(xclosest,normal_closest, &
           inode,ielem, &
           unsigned_mindist, &
           xx,inplane,part_id,time,dxBB)
          ! check distance to the nodes of a triangular element.
          call checkinpoint(xclosest,normal_closest, &
           inode,ielem, &
           unsigned_mindist, &
           xx,inplane,part_id,time,dxBB)
         enddo ! inode=1,nodes_per_elem

         hitflag=0
         hitsign=zero

         do dir=1,3
          if (dir.eq.1) then
           xleft(dir)=xdata3D(i-1,j,k,dir)
           xright(dir)=xdata3D(i+1,j,k,dir)
          else if (dir.eq.2) then
           xleft(dir)=xdata3D(i,j-1,k,dir)
           xright(dir)=xdata3D(i,j+1,k,dir)
          else if (dir.eq.3) then
           xleft(dir)=xdata3D(i,j,k-1,dir)
           xright(dir)=xdata3D(i,j,k+1,dir)
          else
           print *,"dir invalid"
           stop
          endif
          if (xright(dir).le.xleft(dir)) then
           print *,"xright(dir).le.xleft(dir)"
           stop
          endif
         enddo ! dir=1..3

         ! normal points from solid to fluid
         ! phi>0 in the fluid
         ! xclosest=xx-phi n
         ! phi n = xx-xclosest
         ! phi = n dot (xx-xclosest)
         if (inplane.eq.1) then
          n_dot_x=zero
          mag_n=zero
          mag_n_test=zero
          mag_x=zero
          signtest=zero

          in_sign_box=1

          do dir=1,3
           if ((xclosest(dir).lt.xleft(dir)).or. &
               (xclosest(dir).gt.xright(dir))) then
            in_sign_box=0
           endif
           mag_n=mag_n+normal_closest(dir)**2
           mag_n_test=mag_n_test+normal(dir)**2
           !xx=cell center where LS needed
           mag_x=mag_x+(xx(dir)-xclosest(dir))**2 
           n_dot_x=n_dot_x+normal_closest(dir)*(xx(dir)-xclosest(dir))
           signtest=signtest+normal_closest(dir)*normal(dir)
          enddo ! dir=1..3
          mag_n_test=sqrt(mag_n_test)
          mag_n=sqrt(mag_n)
          mag_x=sqrt(mag_x)

          if (abs(mag_x).lt.unsigned_mindist) then
           unsigned_mindist=abs(mag_x)
          endif

          if (mag_n.gt.zero) then
           if (mag_n_test.gt.zero) then
            if (mag_x.gt.zero) then
             if (signtest.gt.zero) then
              if (signtest/(mag_n*mag_n_test).gt.one/sqrt(two)) then
               if (abs(n_dot_x)/(mag_n*mag_x).gt.one/sqrt(two)) then

                if (in_sign_box.eq.1) then

                 hitflag=1
                 used_for_trial=1
                 if (phicenter.eq.zero) then
                  hitsign=zero
                 else if (phicenter.gt.zero) then ! fluid
                  hitsign=one
                 else if (phicenter.lt.zero) then ! solid
                  hitsign=-one
                 else
                  print *,"phicenter bust"
                  stop
                 endif
    
                else if (in_sign_box.eq.0) then
                 ! do nothing
                else
                 print *,"in_sign_box invalid"
                 stop
                endif

               endif  ! abs(n_dot_x) big enough?
              endif ! signtest big enough?
             endif ! signtest>0 ?
            endif ! mag_x>0 ?
           endif ! mag_n_test>0?
          endif ! mag_n>0 ?

         else if (inplane.eq.0) then
          ! do nothing
         else
          print *,"inplane invalid"
          stop
         endif

! check crossing between cells

         do ii=-1,1
         do jj=-1,1
         do kk=-1,1
          if (abs(ii)+abs(jj)+abs(kk).gt.0) then

           do dir=1,3
            if (dir.eq.1) then
             ii_current=ii
            else if (dir.eq.2) then
             ii_current=jj
            else if (dir.eq.3) then
             ii_current=kk
            else
             print *,"dir invalid"
             stop
            endif
            if (ii_current.eq.-1) then
             xside(dir)=xleft(dir)
            else if (ii_current.eq.0) then
             xside(dir)=xx(dir)
            else if (ii_current.eq.1) then
             xside(dir)=xright(dir)
            else
             print *,"ii_current invalid"
             stop
            endif
           enddo ! dir=1..3

           ! normal points from solid to fluid
           ! phi=n dot (x-xnot)
           ! phi>0 in the fluid
           phiside=zero
           do dir=1,3
            phiside=phiside+normal(dir)*(xside(dir)-xnot(dir))
           enddo
           if (phiside*phicenter.le.zero) then
            do dir=1,3
             if (phiside.eq.zero) then
              xcrit(dir)=xside(dir)
             else if (phicenter.eq.zero) then
              xcrit(dir)=xx(dir)
             else
              xcrit(dir)=(abs(phiside)*xx(dir)+ &
                           abs(phicenter)*xside(dir))/  &
                          (abs(phicenter)+abs(phiside))
             endif
            enddo  ! dir

            call checkinplane(xcrit,ielem,inplane, &
             minnode,maxnode,element_scale,part_id,time)
! totaldist is the distance between xx and xside
! testdist  is the distance between xx and xcrit 
! xx=center point  xside=stencil point xcrit=crossing point
            if (inplane.eq.1) then  ! crossing point on the element?
             testdist=zero
             totaldist=zero
             do dir=1,3
              testdist=testdist+(xcrit(dir)-xx(dir))**2
              totaldist=totaldist+(xside(dir)-xx(dir))**2
             enddo
             testdist=sqrt(testdist)
             totaldist=sqrt(totaldist)
             if (testdist.le.unsigned_mindist) then
              unsigned_mindist=testdist
              hitflag=1
              used_for_trial=1
              if (phicenter.eq.zero) then
               hitsign=zero
              else if (phicenter.gt.zero) then
               hitsign=one
              else if (phicenter.lt.zero) then
               hitsign=-one
              else
               print *,"phicenter bust"
               stop
              endif
             endif

             if (totaldist.eq.zero) then
              print *,"totaldist invalid"
              print *,"totaldist= ",totaldist
              stop
             else if (totaldist.lt.testdist-1.0E-10) then
              print *,"cannot have totaldist<testdist"
              print *,"totaldist= ",totaldist
              print *,"testdist= ",testdist
              stop
             else
              testdist=testdist/totaldist
              if (testdist.gt.one) then
               testdist=one
              endif
             endif

            else if (inplane.eq.0) then
             ! do nothing
            else
             print *,"inplane invalid"
             stop
            endif  ! inplane
           endif ! phiside x phicenter <= 0
          endif ! abs(ii)+abs(jj)+abs(kk)>0
         enddo 
         enddo 
         enddo  ! ii,jj,kk

         modify_vel=0

         ! nmat x (vel + LS + temperature + flag + stress)
         ls_local=FSIdata3D(i,j,k,ibase+4)
         mask_local=NINT(FSIdata3D(i,j,k,ibase+6))
         do dir=1,3
          vel_local(dir)=FSIdata3D(i,j,k,ibase+dir)
         enddo
         do dir=1,6
          stress_local(dir)=FSIdata3D(i,j,k,ibase+6+dir)
         enddo
         temp_local=FSIdata3D(i,j,k,ibase+5)

         if ((mask_local.eq.0).or. &
             (mask_local.eq.10).or. &
             (unsigned_mindist.lt.abs(ls_local))) then

          modify_vel=1

          if ((ctml_part_id.ge.1).and. &
              (ctml_part_id.le.CTML_NPARTS)) then
           if (FSI(part_id)%ElemDataBIG(3,ielem).eq.1) then ! doubly wetted
            ! do nothing
           else
            print *,"expecting doubly wetted"
            stop
           endif
          else if (ctml_part_id.eq.0) then
           ! do nothing
          else
           print *,"ctml_part_id invalid"
           stop
          endif

          if (FSI(part_id)%ElemDataBIG(3,ielem).eq.1) then ! doubly wetted
           mask_local=103
           ls_local=-unsigned_mindist
          else if (FSI(part_id)%ElemDataBIG(3,ielem).eq.0) then !singly wetted
           if (mask_local.eq.0) then
            mask_local=1  ! vel init, sign not.
           else if (mask_local.eq.10) then
            mask_local=11 ! vel init, coarse sign init
           else if ((mask_local.eq.1).or. &
                    (mask_local.eq.11).or. &
                    (mask_local.eq.2).or. &
                    (mask_local.eq.3).or. &
                    (mask_local.eq.103)) then
            ! do nothing
           else
            print *,"mask_local invalid"
            stop
           endif
           if ((mask_local.eq.0).or. &
               (mask_local.eq.1).or. &
               (mask_local.eq.3).or. &
               (mask_local.eq.103)) then
            ls_local=unsigned_mindist
           else if ((mask_local.eq.10).or. &
                    (mask_local.eq.11).or. &
                    (mask_local.eq.2)) then
            ls_local=sign_funct(ls_local)*unsigned_mindist 
           else
            print *,"mask_local invalid"
            stop
           endif 
           if (hitflag.eq.1) then
            mask_local=2  ! sign init
            ls_local=hitsign*abs(ls_local)
           else if (hitflag.eq.0) then
            ! do nothing
           else
            print *,"hitflag invalid"
            stop
           endif 
          else
           print *,"wetted flag invalid"
           stop
          endif

          if (modify_vel.eq.1) then

           do dir=1,3
            vel_local(dir)=0.0d0
           enddo
           do dir=1,6
            stress_local(dir)=0.0d0
           enddo
           temp_local=0.0d0

           weighttotal=0.0d0
           do inode=1,nodes_per_elem
            nodeptr=FSI(part_id)%IntElemBIG(inode,ielem)
            do dir=1,3
             xfoot(dir)=FSI(part_id)%NodeBIG(dir,nodeptr)
             velparm(dir)=FSI(part_id)%NodeVelBIG(dir,nodeptr)
            enddo
            massparm=FSI(part_id)%NodeMassBIG(nodeptr)
            if (massparm.le.zero) then
             print *,"massparm invalid"
             stop
            endif
            do dir=1,6
             forceparm(dir)=FSI(part_id)%NodeForceBIG(dir,nodeptr)
            enddo
            call get_target_from_foot(xfoot,xtarget, &
             velparm,time,part_id)
    
            distwt=0.0d0
            do dir=1,3
             distwt=distwt+(xtarget(dir)-xx(dir))**2
            enddo
            weight=1.0d0/( (distwt+(1.0E-10)**2)**4 )
            weight=weight*massparm
            do dir=1,3
             vel_local(dir)=vel_local(dir)+weight*velparm(dir)
            enddo
            do dir=1,6
             stress_local(dir)=stress_local(dir)+ &
              weight*FSI(part_id)%NodeForceBIG(dir,nodeptr)
            enddo
            temp_local=temp_local+ & 
             weight*FSI(part_id)%NodeTempBIG(nodeptr)
            weighttotal=weighttotal+weight
           enddo ! inode=1..nodes_per_elem

           do dir=1,3
            vel_local(dir)=vel_local(dir)/weighttotal
           enddo
           do dir=1,6
            stress_local(dir)=stress_local(dir)/weighttotal
           enddo

           if ((probtype.eq.9).and.(axis_dir.gt.1)) then
            if ( xx(1) .gt. -0.4 .and. xx(1) .lt. -0.2 ) then
             if ( xx(3) .lt. 0.045 .and. xx(3) .gt. 0.03 ) then
              vel_local(3) = -1.0
             endif
            endif
           endif

           temp_local=temp_local/weighttotal

          else if (modify_vel.eq.0) then
           ! do nothing
          else
           print *,"modify_vel invalid"
           stop
          endif

         else if ( ((mask_local.eq.1).and. &
                    (unsigned_mindist.ge.abs(ls_local))).or. &
                   ((mask_local.eq.2).and. &
                    (unsigned_mindist.ge.abs(ls_local))).or. &
                   ((mask_local.eq.3).and. &
                    (unsigned_mindist.ge.abs(ls_local))).or. &
                   ((mask_local.eq.103).and. &
                    (unsigned_mindist.ge.abs(ls_local))).or. &
                   ((mask_local.eq.11).and. &
                    (unsigned_mindist.ge.abs(ls_local))) ) then
          ! do nothing
         else
          print *,"mask_local invalid"
          stop
         endif 

         ! nparts x (vel + LS + temperature + flag + stress)
         FSIdata3D(i,j,k,ibase+4)=ls_local
         FSIdata3D(i,j,k,ibase+6)=mask_local
         do dir=1,3
          FSIdata3D(i,j,k,ibase+dir)=vel_local(dir)
         enddo 
         if (CTML_force_model.eq.1) then
          do dir=1,6
           FSIdata3D(i,j,k,ibase+6+dir)=stress_local(dir)
          enddo 
         else if (CTML_force_model.eq.0) then
          ! do nothing
         else
          print *,"CTML_force_model invalid"
          stop
         endif
         FSIdata3D(i,j,k,ibase+5)=temp_local

        else if ((mask1.eq.1).and.(mask2.eq.0)) then
         ! do nothing
        else
         print *,"mask1 or mask2 invalid"
         stop
        endif
       else if ((i.ge.FSI_growlo(1)).and.(i.le.FSI_growhi(1)).and. &
                (j.ge.FSI_growlo(2)).and.(j.le.FSI_growhi(2)).and. &
                (k.ge.FSI_growlo(3)).and.(k.le.FSI_growhi(3))) then
        ! do nothing
       else
        print *,"i,j,k outside of FSI_growlo,FSI_growhi range"
        stop
       endif 
      enddo
      enddo
      enddo ! i,j,k=gridlo..gridhi

      if (1.eq.0) then
       if (in_the_interior.eq.1) then
        if (used_for_trial.eq.0) then
         print *,"used_for_trial.eq.0"
         print '(A8,3(f9.3))',"minnode ",minnode(1),minnode(2),minnode(3)
         print '(A8,3(f9.3))',"maxnode ",maxnode(1),maxnode(2),maxnode(3)
         print '(A7,3(I10))',"gridlo ",gridlo(1),gridlo(2),gridlo(3)
         print '(A7,3(I10))',"gridhi ",gridhi(1),gridhi(2),gridhi(3)
         print '(A6,3(f9.3))',"xelem ",xelem(1),xelem(2),xelem(3)
         print '(A5,3(f9.3))',"xnot ",xnot(1),xnot(2),xnot(3)
         print '(A7,3(f9.3))',"normal ",normal(1),normal(2),normal(3)
        else if (used_for_trial.eq.1) then
         ! do nothing
        else
         print *,"used_for_trial invalid"
         stop
        endif
       else if (in_the_interior.eq.0) then
        ! do nothing
       else
        print *,"in_the_interior invalid"
        stop
       endif
      endif

     else if (null_intersection.eq.1) then
      ! do nothing
     else
      print *,"null_intersection invalid"
      stop
     endif

    enddo ! ielem_container

    do dir=1,3
     if (dx3D(dir).le.zero) then
      print *,"dx3D(dir).le.zero"
      stop
     endif
    enddo ! dir=1..3

    do i=FSI_lo(1),FSI_hi(1)
    do j=FSI_lo(2),FSI_hi(2)
    do k=FSI_lo(3),FSI_hi(3)

     mask1=NINT(masknbr3D(i,j,k,1))
     mask2=NINT(masknbr3D(i,j,k,2))

     ! mask2==1 if (i,j,k) interior to tile
     ! mask1==0 if (i,j,k) is exterior to tile and is a
     !  coarse/fine ghost cell.
     if ((mask1.eq.0).or.(mask2.eq.1)) then

      ! nparts x (vel + LS + temperature + flag + stress)
      ls_local=FSIdata3D(i,j,k,ibase+4)
      mask_local=NINT(FSIdata3D(i,j,k,ibase+6))

      if (mask_local.eq.0) then
       ls_local=8.0*dx3D(1)
      else if (mask_local.eq.10) then
       ! do nothing, just use the value from fill coarse patch
      else if ((mask_local.eq.1).or. &
               (mask_local.eq.11).or. &
               (mask_local.eq.2).or. & ! valid sign (singly wetted)
               (mask_local.eq.3).or. & ! valid sign (doubly wetted)
               (mask_local.eq.103)) then  ! valid sign and d (doubly wetted)
       ! do nothing
      else
       print *,"mask_local invalid"
       stop
      endif

      if ((ctml_part_id.ge.1).and. &
          (ctml_part_id.le.CTML_NPARTS)) then
       if (exclusive_doubly_wetted.eq.1) then
        ! do nothing
       else
        print *,"expecting doubly wetted"
        stop
       endif
      else if (ctml_part_id.eq.0) then
       ! do nothing
      else
       print *,"ctml_part_id invalid"
       stop
      endif

      if (exclusive_doubly_wetted.eq.0) then
       ! do nothing
      else if (exclusive_doubly_wetted.eq.1) then
       if (mask_local.eq.3) then
        ! do nothing
       else if (mask_local.eq.103) then
        ! do nothing
       else if ((mask_local.eq.2).or. &
                (mask_local.eq.10).or. &
                (mask_local.eq.11).or. &
                (mask_local.eq.0).or. &
                (mask_local.eq.1)) then
        mask_local=3
       else
        print *,"mask_local invalid"
        stop
       endif
       ls_local=-abs(ls_local)
      else
       print *,"exclusive_doubly_wetted invalid"
       stop
      endif

      ! nparts x (vel + LS + temperature + flag + stress)
      FSIdata3D(i,j,k,ibase+4)=ls_local
      FSIdata3D(i,j,k,ibase+6)=mask_local

     else if ((mask1.eq.1).and.(mask2.eq.0)) then
      ! do nothing
     else
      print *,"mask1 or mask2 invalid"
      stop
     endif

    enddo   
    enddo   
    enddo   

! check "invert_solid_levelset"
! fix for doubly wetted 

    wallthick=1.0
    do i=FSI_lo(1),FSI_hi(1)
    do j=FSI_lo(2),FSI_hi(2)
    do k=FSI_lo(3),FSI_hi(3)

     mask1=NINT(masknbr3D(i,j,k,1))
     mask2=NINT(masknbr3D(i,j,k,2))

     ! mask2==1 if (i,j,k) interior to tile
     ! mask1==0 if (i,j,k) is exterior to tile and is a
     !  coarse/fine ghost cell.
     if ((mask1.eq.0).or.(mask2.eq.1)) then

      ! nparts x (vel + LS + temperature + flag + stress)
      ls_local=FSIdata3D(i,j,k,ibase+4)
      mask_local=NINT(FSIdata3D(i,j,k,ibase+6))
      if (mask_local.eq.103) then   !doubly wetted, d init
       ls_local=ls_local+wallthick*dx3D(1)
      else if (mask_local.eq.3) then !doubly wetted, d not init
       ! do nothing
      else if ((mask_local.eq.0).or. &
               (mask_local.eq.10).or. &
               (mask_local.eq.1).or. &
               (mask_local.eq.11).or. &
               (mask_local.eq.2)) then
       ! do nothing
      else
       print *,"mask_local invalid"
       stop
      endif

      if ((ctml_part_id.ge.1).and. &
          (ctml_part_id.le.CTML_NPARTS)) then
       if (invert_solid_levelset.eq.0) then
        ! do nothing
       else
        print *,"expecting invert_solid_levelset==0"
        stop
       endif
      else if (ctml_part_id.eq.0) then
       ! do nothing
      else
       print *,"ctml_part_id invalid"
       stop
      endif

       ! normals default to pointing from solid to fluid (LS>0 in the fluid)
      if (invert_solid_levelset.eq.0) then

       if (mask_local.eq.2) then ! sign valid for singly wetted
        ls_local=-ls_local
       else if (mask_local.eq.3) then ! sign valid for doubly wetted
        ! do nothing
       else if (mask_local.eq.103) then ! sign valid for doubly wetted
        ! do nothing
       else if ((mask_local.eq.0).or. &
                (mask_local.eq.10).or. &
                (mask_local.eq.1).or. &
                (mask_local.eq.11)) then
        ! do nothing
       else
        print *,"mask_local invalid"
        stop
       endif

      else if (invert_solid_levelset.eq.1) then

       if (mask_local.eq.2) then ! sign valid for singly wetted
        ! do nothing
       else if (mask_local.eq.3) then ! sign valid for doubly wetted
        ls_local=-ls_local
       else if (mask_local.eq.103) then ! sign valid for doubly wetted
        ls_local=-ls_local
       else if ((mask_local.eq.0).or. &
                (mask_local.eq.10).or. &
                (mask_local.eq.1).or. &
                (mask_local.eq.11)) then
        ! do nothing
       else
        print *,"mask_local invalid"
        stop
       endif

      else
       print *,"invert_solid_levelset invalid"
       stop
      endif
      FSIdata3D(i,j,k,ibase+4)=ls_local

      if (sign_valid(mask_local).eq.1) then
       if (ls_local.lt.LSMIN_debug) then
        LSMIN_debug=ls_local
        mask_debug=mask_local
       endif
       if (ls_local.gt.LSMAX_debug) then
        LSMAX_debug=ls_local
        mask_debug=mask_local
       endif
      else if (sign_valid(mask_local).eq.0) then
       ! do nothing
      else
       print *,"sign_valid(mask_local) invalid"
       stop
      endif

     else if ((mask1.eq.1).and.(mask2.eq.0)) then
      ! do nothing
     else
      print *,"mask1 or mask2 invalid"
      stop
     endif
  
    enddo
    enddo
    enddo

   else if (FSI_operation.eq.3) then

    allocate(old_FSIdata(DIMV3D(FSIdata3D),nFSI))

    do i=FSI_growlo(1),FSI_growhi(1)
    do j=FSI_growlo(2),FSI_growhi(2)
    do k=FSI_growlo(3),FSI_growhi(3)
     do nc=1,nFSI
      old_FSIdata(i,j,k,nc)=FSIdata3D(i,j,k,nc)
     enddo
    enddo
    enddo
    enddo

    if ((touch_flag.ne.0).and.(touch_flag.ne.1)) then
     print *,"touch_flag invalid"
     stop
    endif

    numtouch=0

    if ((iter.ge.0).and.(iter.le.ngrowFSI)) then
     dirmax=1
     sweepmax=1
    else if (iter.gt.ngrowFSI) then
     dirmax=6
     sweepmax=2
    else
     print *,"iter invalid"
     stop
    endif

    do dir_order=1,dirmax
    do sweep=1,sweepmax

     if (dir_order.eq.1) then
      dirmap(1)=1
      dirmap(2)=2
      dirmap(3)=3
     else if (dir_order.eq.2) then
      dirmap(1)=1
      dirmap(2)=3
      dirmap(3)=2
     else if (dir_order.eq.3) then
      dirmap(1)=2
      dirmap(2)=1
      dirmap(3)=3
     else if (dir_order.eq.4) then
      dirmap(1)=2
      dirmap(2)=3
      dirmap(3)=1
     else if (dir_order.eq.5) then
      dirmap(1)=3
      dirmap(2)=2
      dirmap(3)=1
     else if (dir_order.eq.6) then
      dirmap(1)=3
      dirmap(2)=1
      dirmap(3)=2
     else
      print *,"dir_order invalid"
      stop
     endif
    
     if (sweep.eq.1) then
      incr=1
     else if (sweep.eq.2) then
      incr=-1
     else
      print *,"sweep invalid"
      stop
     endif

      ! the interval must be growlo3D ... growhi3D
      ! if dirmax>1 sweepmax>1
     do dir=1,3
      if (sweep.eq.1) then
       index_first(dir)=growlo3D(dirmap(dir))
       index_last(dir)=growhi3D(dirmap(dir))
      else if (sweep.eq.2) then
       index_first(dir)=growhi3D(dirmap(dir))
       index_last(dir)=growlo3D(dirmap(dir))
      else
       print *,"sweep invalid"
       stop
      endif
     enddo ! dir=1..3

     if (((dir_order.eq.1).and.(sweep.eq.1)).or. &
         (numtouch.gt.0)) then

      numtouch=0

      do imap=index_first(1),index_last(1),incr 
      do jmap=index_first(2),index_last(2),incr 
      do kmap=index_first(3),index_last(3),incr 
       map_index(1)=imap
       map_index(2)=jmap
       map_index(3)=kmap
       do dir=1,3
        inverse_index(dirmap(dir))=map_index(dir)
       enddo     
       i=inverse_index(1)
       j=inverse_index(2)
       k=inverse_index(3)

       mask1=NINT(masknbr3D(i,j,k,1))
       mask2=NINT(masknbr3D(i,j,k,2))

       ! mask2==1 if (i,j,k) interior to tile
       ! mask1==0 if (i,j,k) is exterior to tile and is a
       !  coarse/fine ghost cell.
       if ((mask1.eq.0).or.(mask2.eq.1)) then

        ! nparts x (vel + LS + temperature + flag + stress)
        mask_local=NINT(old_FSIdata(i,j,k,ibase+6))
        ls_local=old_FSIdata(i,j,k,ibase+4)
        new_mask_local=mask_local

        do dir=1,3
         xcen(dir)=xdata3D(i,j,k,dir)
        enddo

        if (sign_valid(mask_local).eq.0) then

         sign_defined=0
         LSsign=zero
         sign_status_changed=0

         do dir=1,3

          ii=0
          jj=0
          kk=0
          if (dir.eq.1) then
           i_norm=i
           ii=1
          else if (dir.eq.2) then
           i_norm=j
           jj=1
          else if (dir.eq.3) then
           i_norm=k
           kk=1
          else
           print *,"dir invalid"
           stop
          endif

          if (i_norm.eq.FSI_growlo(dir)) then
           mminus=0
          else if ((i_norm.gt.FSI_growlo(dir)).and. &
                   (i_norm.le.FSI_growhi(dir))) then
           ! nparts x (vel + LS + temperature + flag + stress)
           do nc=1,nFSI_sub
            data_minus(nc)=old_FSIdata(i-ii,j-jj,k-kk,ibase+nc)
           enddo
           mminus=NINT(data_minus(6))
           LSMINUS=data_minus(4)
           xminus=xdata3D(i-ii,j-jj,k-kk,dir)
          else 
           print *,"i_norm invalid"
           stop
          endif
   
          if (i_norm.eq.FSI_growhi(dir)) then
           mplus=0
          else if ((i_norm.ge.FSI_growlo(dir)).and. &
                   (i_norm.lt.FSI_growhi(dir))) then
           ! nparts x (vel + LS + temperature + flag + stress)
           do nc=1,nFSI_sub
            data_plus(nc)=old_FSIdata(i+ii,j+jj,k+kk,ibase+nc)
           enddo
           mplus=NINT(data_plus(6))
           LSPLUS=data_plus(4)
           xplus=xdata3D(i+ii,j+jj,k+kk,dir)
          else 
           print *,"i_norm invalid"
           stop
          endif

          ! sign_valid==1 if mask=2,3,10,11   sign_valid==0 if mask=0,1
          if ((sign_valid(mplus).eq.1).and. &
              (sign_valid(mminus).eq.1)) then
           if (abs(LSPLUS).le.abs(LSMINUS)) then
            sign_defined_local=sign_funct(LSPLUS)
            if ((abs(LSPLUS).le.abs(LSsign)).or. &
                (sign_defined.eq.0)) then
             sign_defined=sign_defined_local
             LSsign=LSPLUS
             if (sign_valid(mask_local).eq.0) then
              ls_local=sign_defined*abs(ls_local)
              sign_status_changed=1
             else if (sign_valid(mask_local).eq.1) then
              ! do nothing
             else
              print *,"mask_local invalid"
              stop
             endif
            endif
           else if (abs(LSPLUS).ge.abs(LSMINUS)) then
            sign_defined_local=sign_funct(LSMINUS)
            if ((abs(LSMINUS).le.abs(LSsign)).or. &
                (sign_defined.eq.0)) then
             sign_defined=sign_defined_local
             LSsign=LSMINUS
             if (sign_valid(mask_local).eq.0) then
              ls_local=sign_defined*abs(ls_local)
              sign_status_changed=1
             else if (sign_valid(mask_local).eq.1) then
              ! do nothing
             else
              print *,"mask_local invalid"
              stop
             endif
            endif
           else
            print *,"LSPLUS or LSMINUS invalid"
            print *,"LSPLUS, LSMINUS ",LSPLUS,LSMINUS
            print *,"dir=",dir
            print *,"i,j,k ",i,j,k
            print *,"mask_local=",mask_local
            print *,"ls_local=",ls_local
            print *,"mplus=",mplus
            print *,"mminus=",mminus
            stop
           endif
          else if ((sign_valid(mplus).eq.1).and. &
                   (sign_valid(mminus).eq.0)) then
           sign_defined_local=sign_funct(LSPLUS)
           if ((abs(LSPLUS).le.abs(LSsign)).or. &
               (sign_defined.eq.0)) then
            sign_defined=sign_defined_local
            LSsign=LSPLUS
            if (sign_valid(mask_local).eq.0) then
             ls_local=sign_defined*abs(ls_local)
             sign_status_changed=1
            else if (sign_valid(mask_local).eq.1) then
             ! do nothing
            else
             print *,"mask_local invalid"
             stop
            endif
           endif
          else if ((sign_valid(mplus).eq.0).and. &
                   (sign_valid(mminus).eq.1)) then
           sign_defined_local=sign_funct(LSMINUS)
           if ((abs(LSMINUS).le.abs(LSsign)).or. &
               (sign_defined.eq.0)) then
            sign_defined=sign_defined_local
            LSsign=LSMINUS
            if (sign_valid(mask_local).eq.0) then
             ls_local=sign_defined*abs(ls_local)
             sign_status_changed=1
            else if (sign_valid(mask_local).eq.1) then
             ! do nothing
            else
             print *,"mask_local invalid"
             stop
            endif
           endif
          else if ((sign_valid(mplus).eq.0).and. &
                   (sign_valid(mminus).eq.0)) then
           ! do nothing
          else
           print *,"mplus or mminus invalid"
           stop
          endif

         enddo ! dir=1..sdim
        
         if (sign_status_changed.eq.1) then 
          new_mask_local=2
          touch_flag=1
          numtouch=numtouch+1
         else if (sign_status_changed.eq.0) then 
          ! do nothing
         else
          print *,"sign_status_changed invalid"
          stop
         endif

         if (vel_valid(mask_local).eq.0) then 

          if (nFSI_sub.ne.12) then
           print *,"nFSI_sub.ne.12"
           stop
          endif

          ! (vel + LS + temperature + flag + stress)
          do dir=1,nFSI_sub
           weight_top(dir)=zero
          enddo
          weight_bot=zero
          do i1=-3,3
          do j1=-3,3
          do k1=-3,3
           if ((i+i1.ge.FSI_growlo(1)).and.(i+i1.le.FSI_growhi(1)).and. &
               (j+j1.ge.FSI_growlo(2)).and.(j+j1.le.FSI_growhi(2)).and. &
               (k+k1.ge.FSI_growlo(3)).and.(k+k1.le.FSI_growhi(3))) then
            mask_node=NINT(old_FSIdata(i+i1,j+j1,k+k1,ibase+6))
            if (vel_valid(mask_node).eq.1) then
             weight=zero
             do dir=1,3
              xc(dir)=xdata3D(i+i1,j+j1,k+k1,dir)
              weight=weight+(xc(dir)-xcen(dir))**2
             enddo
             if (weight.le.zero) then
              print *,"weight invalid"
              stop
             endif
             weight=one/weight
             ! nparts x (vel + LS + temperature + flag + stress)
             do dir=1,3
              weight_top(dir)=weight_top(dir)+ &
               old_FSIdata(i+i1,j+j1,k+k1,ibase+dir)*weight
             enddo
              ! temperature
             weight_top(5)=weight_top(5)+ &
               old_FSIdata(i+i1,j+j1,k+k1,ibase+5)*weight
              ! stress
             do dir=1,6
              weight_top(6+dir)=weight_top(6+dir)+ &
               old_FSIdata(i+i1,j+j1,k+k1,ibase+6+dir)*weight
             enddo
     
             weight_bot=weight_bot+weight
            else if (vel_valid(mask_node).eq.0) then
             ! do nothing
            else
             print *,"vel_valid(mask_node) invalid"
             stop
            endif
           endif ! i1,j1,k1 in grid
          enddo
          enddo
          enddo ! i1,j1,k1
          if (weight_bot.gt.zero) then
           ! (vel + LS + temperature + flag + stress)
           do dir=1,3
            FSIdata3D(i,j,k,ibase+dir)=weight_top(dir)/weight_bot
           enddo
           FSIdata3D(i,j,k,ibase+5)=weight_top(5)/weight_bot
           if (CTML_force_model.eq.1) then
            do dir=1,6
             FSIdata3D(i,j,k,ibase+6+dir)=weight_top(6+dir)/weight_bot
            enddo
           else if (CTML_force_model.eq.0) then
            ! do nothing
           else
            print *,"CTML_force_model invalid"
            stop
           endif

           if (sweepmax.eq.1) then
            ! do nothing
           else if (sweepmax.eq.2) then
            ! (vel + LS + temperature + flag + stress)
            do dir=1,3
             old_FSIdata(i,j,k,ibase+dir)=weight_top(dir)/weight_bot
            enddo
            old_FSIdata(i,j,k,ibase+5)=weight_top(5)/weight_bot
            if (CTML_force_model.eq.1) then
             do dir=1,6
              old_FSIdata(i,j,k,ibase+6+dir)=weight_top(6+dir)/weight_bot
             enddo
            else if (CTML_force_model.eq.0) then
             ! do nothing
            else
             print *,"CTML_force_model invalid"
             stop
            endif
           else
            print *,"sweepmax invalid"
            stop
           endif
           
          endif 
         else if (vel_valid(mask_local).eq.1) then
          ! do nothing
         else
          print *,"vel_valid(mask_local) invalid"
          stop
         endif
        else if (sign_valid(mask_local).eq.1) then
         ! do nothing
        else
         print *,"mask_local invalid"
         stop
        endif

        if (sign_valid(new_mask_local).eq.1) then
         if (ls_local.lt.LSMIN_debug) then
          LSMIN_debug=ls_local
          mask_debug=new_mask_local
         endif
         if (ls_local.gt.LSMAX_debug) then
          LSMAX_debug=ls_local
          mask_debug=new_mask_local
         endif
        else if (sign_valid(new_mask_local).eq.0) then
         ! do nothing
        else
         print *,"new_mask_local invalid"
         stop
        endif

        ! (vel + LS + temperature + flag + stress)
        FSIdata3D(i,j,k,ibase+6)=new_mask_local 
        FSIdata3D(i,j,k,ibase+4)=ls_local 

        if (sweepmax.eq.1) then
         ! do nothing
        else if (sweepmax.eq.2) then
         old_FSIdata(i,j,k,ibase+6)=new_mask_local 
         old_FSIdata(i,j,k,ibase+4)=ls_local 
        else
         print *,"sweepmax invalid"
         stop
        endif

       else if ((mask1.eq.1).and.(mask2.eq.0)) then
        ! do nothing
       else
        print *,"mask1 or mask2 invalid"
        stop
       endif

      enddo 
      enddo 
      enddo  ! imap,jmap,kmap=index_first ... index_last,incr

     else if (((dir_order.gt.1).or.(sweep.eq.2)).and. &
              (numtouch.eq.0)) then
      ! do nothing
     else
      print *,"dir_order, sweep or numtouch invalid"
      stop
     endif

    enddo ! sweep
    enddo ! dir_order

    if (ioproc.eq.1) then
     print *,"numtouch ",numtouch
    endif

    deallocate(old_FSIdata)
   else
    print *,"FSI_operation invalid"
    stop
   endif

   if (isout.eq.0) then
    ! do nothing
   else if (isout.eq.1) then
    print *,"LSMIN_debug=",LSMIN_debug
    print *,"LSMAX_debug=",LSMAX_debug
    print *,"mask_debug=",mask_debug
    print *,"END: CLSVOF_InitBox"
    print *,"FSI_operation=",FSI_operation
    print *,"touch_flag=",touch_flag
   else
    print *,"isout invalid3: ",isout
    stop
   endif

   FLUSH(6)  ! 6=screen

  else if ((ctml_part_id.eq.0).and. &
           (fsi_part_id.eq.0)) then
   ! do nothing
  else 
   print *,"ctml_part_id or fsi_part_id invalid"
   stop
  endif


return
end subroutine CLSVOF_InitBox


      subroutine CLSVOF_Copy_To_LAG(  &
       sdim_AMR, &
       lev77, &
       tid, &
       tilenum, &
       im_part, &
       nparts, &
       part_id, &
       ngrowFSI, &
       nmat, &
       nFSI, &
       nFSI_sub, &
       FSI_operation, &
       time, &
       problo3D,probhi3D, &
       xmap3D, &
       xslice3D, &
       dx3D, &
       FSI_lo,FSI_hi, &
       FSI_growlo,FSI_growhi, &
       growlo3D,growhi3D, &
       xdata3D, &
       veldata3D, &
       masknbr3D, &
       maskfiner3D, &
       DIMS3D(FSIdata3D), &
       ioproc,isout)
       use global_utility_module
#ifdef MVAHABFSI
       use CTML_module
#endif

       IMPLICIT NONE

      INTEGER_T :: sdim_AMR
      INTEGER_T :: lev77
      INTEGER_T :: tid
      INTEGER_T :: tilenum
      INTEGER_T :: im_part
      INTEGER_T :: nparts
      INTEGER_T :: part_id
      INTEGER_T :: ngrowFSI
      INTEGER_T :: nmat
      INTEGER_T :: nFSI
      INTEGER_T :: nFSI_sub
      INTEGER_T :: FSI_operation
      REAL_T :: time
      REAL_T problo3D(3)
      REAL_T probhi3D(3)
      INTEGER_T xmap3D(3)
      REAL_T xslice3D(3)
      REAL_T dx3D(3)
      INTEGER_T FSI_lo(3),FSI_hi(3)
      INTEGER_T FSI_growlo(3),FSI_growhi(3)
      INTEGER_T growlo3D(3),growhi3D(3)
      INTEGER_T DIMDEC3D(FSIdata3D)
      REAL_T xdata3D(DIMV3D(FSIdata3D),3)
      REAL_T veldata3D(DIMV3D(FSIdata3D),3)
      REAL_T masknbr3D(DIMV3D(FSIdata3D),2)
      REAL_T maskfiner3D(DIMV3D(FSIdata3D))

      INTEGER_T ioproc,isout

      REAL_T dxBB(3) ! set in find_grid_bounding_box_node

      INTEGER_T :: inode
      INTEGER_T :: inode_container
      REAL_T, dimension(3) :: xnot
      REAL_T, dimension(3) :: xnode
      INTEGER_T, dimension(3) :: gridlo,gridhi
      INTEGER_T :: i,j,k
      INTEGER_T, dimension(3) :: idx
      REAL_T, dimension(3) :: velparm
      REAL_T, dimension(3) :: total_vel
      INTEGER_T :: dir
      REAL_T :: total_weight,wt,dist_scale,df,support_size

      INTEGER_T local_mask
      INTEGER_T num_nodes
      INTEGER_T num_nodes_container
      INTEGER_T debug_all
      INTEGER_T ctml_part_id
      INTEGER_T fsi_part_id

      debug_all=0

      if ((part_id.lt.1).or.(part_id.gt.nparts)) then
       print *,"part_id invalid"
       stop
      endif
      if (nparts.ne.TOTAL_NPARTS) then
       print *,"nparts.ne.TOTAL_NPARTS"
       stop
      endif

      if ((lev77.lt.1).or.(tid.lt.0).or.(tilenum.lt.0)) then
       print *,"lev77 or tid or tilenum invalid"
       stop
      endif
      if (container_allocated.ne.1) then
       print *,"container_allocated.ne.1"
       stop
      endif
      if (level_container_allocated(lev77).ne.1) then
       print *,"level_container_allocated(lev77).ne.1"
       stop
      endif

      if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.eq.0) then

       print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc=0"
       stop

      else if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.gt.0) then

       if (tilenum+1.gt.contain_elem(lev77)% &
           num_tiles_on_thread3D_proc(tid+1)) then
        print *,"tilenum+1.gt.num_tiles_on_thread3D_proc(tid+1)"
        stop
       endif

      else
       print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc bad"
       stop
      endif

      do dir=1,3
       if (contain_elem(lev77)%tilelo3D(tid+1,tilenum+1,dir).ne. &
           FSI_lo(dir)) then
        print *,"tilelo3D(tid+1,tilenum+1,dir).ne.FSI_lo(dir)"
        stop
       endif 
       if (contain_elem(lev77)%tilehi3D(tid+1,tilenum+1,dir).ne. &
           FSI_hi(dir)) then
        print *,"tilehi3D(tid+1,tilenum+1,dir).ne.FSI_hi(dir)"
        stop
       endif 
       if (FSI_growlo(dir).ge.FSI_lo(dir)) then
        print *,"FSI_growlo(dir) invalid"
        stop
       endif
       if (FSI_growhi(dir).le.FSI_hi(dir)) then
        print *,"FSI_growhi(dir) invalid"
        stop
       endif
      enddo ! dir=1..3

      ctml_part_id=CTML_partid_map(part_id)
      fsi_part_id=FSI_partid_map(part_id)

      if (((ctml_part_id.ge.1).and. &
           (ctml_part_id.le.CTML_NPARTS)).or. &
          ((fsi_part_id.ge.1).and. &
           (fsi_part_id.le.FSI_NPARTS))) then

       num_nodes=FSI(part_id)%NumNodes

       if (num_nodes.le.0) then
        print *,"num_nodes invalid"
        stop
       endif
       if ((ioproc.ne.1).and.(ioproc.ne.0)) then
        print *,"ioproc invalid"
        stop
       endif
       if ((sdim_AMR.ne.2).and.(sdim_AMR.ne.3)) then
        print *,"sdim_AMR invalid"
        stop
       endif
       if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
        print *,"part_id invalid"
        stop
       endif
       if ((nmat.lt.1).or.(nmat.gt.50)) then
        print *,"nmat out of range"
        stop
       endif
       if ((im_part.lt.1).or.(im_part.gt.nmat)) then
        print *,"im_part invalid"
        stop
       endif

       if (isout.eq.0) then
        ! do nothing
       else if (isout.eq.1) then
        print *,"in (START): CLSVOF_Copy_To_LAG"
        print *,"num_nodes(NumNodes)=",num_nodes
        print *,"sdim_AMR=",sdim_AMR 
        print *,"im_part=",im_part
        print *,"part_id=",part_id
        print *,"nmat=",nmat
        print *,"FSI_operation=",FSI_operation
        print *,"ioproc=",ioproc
        print '(A9,3(f12.6))',"problo3D ",problo3D(1),problo3D(2),problo3D(3)
        print '(A9,3(f12.6))',"probhi3D ",probhi3D(1),probhi3D(2),probhi3D(3)
        print '(A7,3(I10))',"FSI_lo ",FSI_lo(1),FSI_lo(2),FSI_lo(3)
        print '(A7,3(I10))',"FSI_hi ",FSI_hi(1),FSI_hi(2),FSI_hi(3)
        print '(A11,3(I10))',"FSI_growlo ", &
         FSI_growlo(1),FSI_growlo(2),FSI_growlo(3)
        print '(A11,3(I10))',"FSI_growhi ", &
         FSI_growhi(1),FSI_growhi(2),FSI_growhi(3)
        print '(A11,3(I10))',"growlo3D ",growlo3D(1),growlo3D(2),growlo3D(3)
        print '(A11,3(I10))',"growhi3D ",growhi3D(1),growhi3D(2),growhi3D(3)
       else
        print *,"isout invalid1: ",isout
        stop
       endif

        ! nparts x (vel + LS + temperature + flag + stress)

       if (nFSI.ne.nparts*nFSI_sub) then
        print *,"nFSI invalid"
        stop
       endif
       if (nFSI_sub.ne.12) then
        print *,"nFSI_sub invalid"
        stop
       endif
       if ((nparts.lt.1).or.(nparts.ge.nmat)) then
        print *,"nparts invalid"
        stop
       endif

       if (FSI_operation.ne.4) then
        print *,"FSI_operation invalid"
        stop
       endif
       if (ngrowFSI.ne.3) then
        print *,"ngrowFSI invalid"
        stop
       endif

       call checkbound3D(FSI_lo,FSI_hi, &
        DIMS3D(FSIdata3D), &
        ngrowFSI,-1,521)

       num_nodes_container=contain_elem(lev77)% &
                           level_node_data(tid+1,part_id,tilenum+1)% &
                           numNodes

       if (isout.eq.0) then
        ! do nothing
       else if (isout.eq.1) then
        print *,"lev77=",lev77
        print *,"tid= ",tid
        print *,"part_id= ",part_id
        print *,"tilenum= ",tilenum
        print *,"num_nodes_container= ",num_nodes_container
       else
        print *,"isout invalid"
        stop
       endif

       if ((num_nodes_container.lt.0).or. &
           (num_nodes_container.gt.num_nodes)) then
        print *,"num_nodes_container invalid"
        stop
       endif

       do inode_container=1,num_nodes_container

        inode=contain_elem(lev77)% &
              level_node_data(tid+1,part_id,tilenum+1)% &
              NodeData(inode_container)

        if ((inode.lt.1).or. &
            (inode.gt.num_nodes)) then
         print *,"inode invalid"
         stop
        endif

        do dir=1,3
         xnode(dir)=FSI(part_id)%Node(dir,inode)
         velparm(dir)=zero
        enddo 
        call get_target_from_foot(xnode,xnot, &
          velparm,time,part_id)

        if (debug_all.eq.1) then
         print *,"inode=",inode
         do dir=1,3
          print *,"dir,xnode ",dir,xnode(dir)
          print *,"dir,xnot  ",dir,xnot(dir)
         enddo
        endif

        call find_grid_bounding_box_node( &
         xnot, &
         FSI_lo,FSI_hi, &
         FSI_growlo,FSI_growhi, &
         growlo3D,growhi3D, &
         xdata3D, &
         DIMS3D(FSIdata3D), &
         gridlo,gridhi,dxBB) 

        do dir=1,3
         if (abs(dxBB(dir)-dx3D(dir)).gt.element_buffer_tol*dxBB(dir)) then
          print *,"abs(dxBB(dir)-dx3D(dir)).gt.element_buffer_tol*dxBB(dir)"
          stop
         endif
        enddo ! dir=1..3

        if (1.eq.0) then
         print *,"inode=",inode
         do dir=1,3
          print *,"dir,gridlo,gridhi ",dir,gridlo(dir),gridhi(dir)
         enddo
        endif

        do dir=1,3
         idx(dir)=(gridhi(dir)+gridlo(dir))/2
         if (gridhi(dir)-gridlo(dir).ne.4) then
          print *,"gridhi(dir)-gridlo(dir).ne.4"
          stop
         endif
        enddo

        local_mask=NINT(maskfiner3D(idx(1),idx(2),idx(3)))

        if (local_mask.eq.1) then
         total_weight=zero
         do dir=1,3
          total_vel(dir)=zero
         enddo
         do i=gridlo(1),gridhi(1)
         do j=gridlo(2),gridhi(2)
         do k=gridlo(3),gridhi(3)
          wt=one
          do dir=1,3

           if ((xmap3D(dir).eq.1).or. &
               (xmap3D(dir).eq.2).or. &
               (xmap3D(dir).eq.sdim_AMR)) then
            dist_scale=abs(xdata3D(i,j,k,dir)-xnot(dir))/dxBB(dir)
            if (CTML_FSI_flagF(nmat).eq.1) then ! FSI_flag==4
#ifdef MVAHABFSI
             call CTML_DELTA(dir,dist_scale,df)
#else
             print *,"define MVAHABFSI"
             stop
#endif
            else if (CTML_FSI_flagF(nmat).eq.0) then
             support_size=two
             df=hsprime(dist_scale,support_size)
            else
             print *,"CTML_FSI_flagF(nmat) invalid"
             stop
            endif 
            wt=wt*df
           else if ((xmap3D(dir).eq.0).and.(sdim_AMR.eq.2)) then
            ! do nothing
           else
            print *,"xmap3D(dir) invalid"
            stop
           endif

          enddo ! dir=1..3
          total_weight=total_weight+wt
          do dir=1,3
           total_vel(dir)=total_vel(dir)+wt*veldata3D(i,j,k,dir)
          enddo
         enddo
         enddo
         enddo
         if (total_weight.gt.zero) then
          do dir=1,3
           FSI(part_id)%NodeVel(dir,inode)=total_vel(dir)/total_weight
          enddo
         else
          print *,"total_weight invalid"
          stop
         endif
        else if (local_mask.eq.0) then
         ! do nothing
        else
         print *,"local_mask invalid"
         stop
        endif

       enddo ! inode_container=1,num_nodes_container

      else if ((ctml_part_id.eq.0).and. &
               (fsi_part_id.eq.0)) then
       ! do nothing
      else 
       print *,"ctml_part_id or fsi_part_id invalid"
       stop
      endif

      return
      end subroutine CLSVOF_Copy_To_LAG


      subroutine flappingKinematics(numMotion,motionPara,r,t)
      IMPLICIT NONE
      INTEGER_T numMotion
      REAL_T xPoint(3,numMotion),vTan(3,numMotion)
      REAL_T x0(3),v(3),motionPara(11,numMotion)
      REAL_T vNorm,theta,thetaMag,fTheta,phiTheta,theta0
      REAL_T t,hMag,fH,phiH,h0,ct,st
      REAL_T r(3,4),r1(3,4),r2(3,4)
      INTEGER_T motionType
      INTEGER_T rotateAlongALine,translateAlongALine
      INTEGER_T i,j,k,iMotion

!     pitching motion theta(t)=theta_mag*cos(2*pi*f_theta*t+phi_theta)+theta_0
!     plunging motion h(t)    =h_mag    *cos(2*pi*f_h    *t+phi_h    )+h_0
      rotateAlongALine   =0
      translateAlongALine=1 

      do iMotion=1,numMotion 

        if (1.eq.0) then
         do i=1,3
           x0(i)=xPoint(i,iMotion)
           v(i) =vTan(i,iMotion)
         enddo
         vNorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
         v(1)=v(1)/vNorm
         v(2)=v(2)/vNorm
         v(3)=v(3)/vNorm
        else
         do i=1,3
           x0(i)=zero
           v(i) =zero
         enddo
         vNorm=zero
        endif

        motionType=motionPara(1,iMotion)

        if (motionType.eq.rotateAlongALine) then

!         assign the pitching parameters

          thetaMag=motionPara(2,iMotion)
          fTheta  =motionPara(3,iMotion)
          phiTheta=motionPara(4,iMotion)
          theta0  =motionPara(5,iMotion)
          x0(1)   =motionPara(6,iMotion)
          x0(2)   =motionPara(7,iMotion)
          x0(3)   =motionPara(8,iMotion)
          v(1)    =motionPara(9,iMotion)
          v(2)    =motionPara(10,iMotion)
          v(3)    =motionPara(11,iMotion)

          theta=thetaMag*cos(2*Pi*fTheta*t+phiTheta)+theta0     

!         the rotationMatrix
          ct = cos(theta) 
          st=sin(theta)

!         Form the rotation matrix: (for a derivation, see the notes in RevolutionMapping.C)
!         R = v v^T + cos(theta) ( I -v v^T ) + sin(theta) ( v X )(  I -v v^T)
!    
!         x = R *( x(0) - x0 ) + x0
!           = R * x(0)  + (I-R)*x0

          r(1,1) = v(1)*v(1)*(1.-ct)+ct
          r(1,2) = v(1)*v(2)*(1.-ct)-st*v(3)
          r(1,3) = v(1)*v(3)*(1.-ct)+st*v(2);
          r(2,1) = v(1)*v(2)*(1.-ct)+st*v(3);  
          r(2,2) = v(2)*v(2)*(1.-ct)+ct
          r(2,3) = v(2)*v(3)*(1.-ct)-st*v(1);
          r(3,1) = v(1)*v(3)*(1.-ct)-st*v(2);  
          r(3,2) = v(3)*v(2)*(1.-ct)+st*v(1);  
          r(3,3) = v(3)*v(3)*(1.-ct)+ct     ;

          r(1,4) = (1.-r(1,1))*x0(1)    -r(1,2) *x0(2)    -r(1,3) *x0(3);
          r(2,4) =    -r(2,1) *x0(1)+(1.-r(2,2))*x0(2)    -r(2,3) *x0(3);
          r(3,4) =    -r(3,1) *x0(1)    -r(3,2) *x0(2)+(1.-r(3,3))*x0(3);

        elseif (motionType.eq.translateAlongALine) then

          hMag    =motionPara(2,iMotion)
          fH      =motionPara(3,iMotion)
          phiH    =motionPara(4,iMotion)
          h0      =motionPara(5,iMotion)
          x0(1)   =motionPara(6,iMotion)
          x0(2)   =motionPara(7,iMotion)
          x0(3)   =motionPara(8,iMotion)
          v(1)    =motionPara(9,iMotion)
          v(2)    =motionPara(10,iMotion)
          v(3)    =motionPara(11,iMotion)

          theta=hMag*cos(2*Pi*fH*t+phiH)+h0

          r(1,1) = 1.
          r(1,2) = 0.
          r(1,3) = 0.
          r(2,1) = 0.  
          r(2,2) = 1.  
          r(2,3) = 0.
          r(3,1) = 0.  
          r(3,2) = 0.  
          r(3,3) = 1.

          r(1,4) = x0(1) + v(1)*theta;
          r(2,4) = x0(2) + v(2)*theta;
          r(3,4) = x0(3) + v(3)*theta;

         !print*,theta,v(1),v(2),v(3)

         !pause

        else 
          print *, "motionType invalid.  wrong kinematics, stop here"
          stop
        endif

        if (iMotion.gt.1) then    ! compound motion
!        --- we compose the current motion with preMotion ---
!         Current:   x  = R1*x0 + g1
!         premotion: x0 = R2*x(0) + g2
!         Composed:
!         x = R1*( R2*x(0) + g2 ) + g1(t)
!           = R1*R2*x(0) + R1*g2+g1

!         store the current transformation matrix
         do i=1,3
         do j=1,4
           r1(i,j)=r(i,j)
         enddo
         enddo 

         do i=1,3
           do j=1,3
            r(i,j)=0.
            do k=1,3
              r(i,j)=r(i,j)+r1(i,k)*r2(k,j)
            enddo
            r(i,4)=r(i,4)+r1(i,j)*r2(j,4)    ! g=R1*g2+g1
           enddo ! j
         enddo ! i

       endif  ! iMotion>=1

!      store the temporary transformation matrix
       do i=1,3
       do j=1,4
          r2(i,j)=r(i,j)
       enddo
       enddo

      enddo  ! end the iMotion 

      return
      end subroutine flappingKinematics

      subroutine get_foot_from_target(xtarget,xfoot, &
       velparm,time,part_id)
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T part_id
      REAL_T velparm(3)
      REAL_T time
      INTEGER_T dir
      REAL_T xtarget(3)
      REAL_T xfoot(3)
      REAL_T xtargetsave(3)
      REAL_T xfootsave(3)
      REAL_T RPM,alpha,RR,radgear,theta

      INTEGER_T numMotion
      REAL_T r(3,4)
      REAL_T rinv(3,4)
      REAL_T rplus(3,4)
      REAL_T det
      REAL_T motionPara(11,2)
      REAL_T flapping_time
      REAL_T dt_flapping
      REAL_T flapping_time_plus


      if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
       print *,"part_id invalid"
       stop
      endif
      if (FSI(part_id)%partID.ne.part_id) then
       print *,"FSI(part_id)%partID.ne.part_id"
       stop
      endif

      do dir=1,3
       velparm(dir)=0.0
      enddo

       ! get_foot_from_target
      do dir=1,3
       xfoot(dir)=xtarget(dir)
       xfootsave(dir)=xtarget(dir)
       xtargetsave(dir)=xtarget(dir)
      enddo

      if (part_id.eq.1) then

       ! gear rotating clockwise
       ! theta(t)=theta_0-alpha t
       ! x=r cos(theta)
       ! y=r sin(theta)
       ! x'=alpha y
       ! y'=-alpha x
       ! xfoot=r cos(theta_0)=r cos(theta+alpha t)=
       !    r(costheta cos(alpha t)-sintheta sin(alpha t) )
       ! yfoot=r sin(theta_0)=r sin(theta+alpha t)=
       !    r(sintheta cos(alpha t)+costheta sin(alpha t))

       if (probtype.eq.563) then
        RPM=abs(vinletgas);
        alpha=two*Pi*RPM/60.0  ! radians/s

        ! at a distance "r" from the origin, a particle on the gear will
        ! travel 2 pi r distance in one revolution (2 pi radians)
        ! the time for one revolution is 2 pi/alpha
        ! so magnitude of velocity is 2 pi r /  ( 2pi /alpha )=alpha r

        velparm(1)=alpha*xtarget(3)
        velparm(2)=zero
        velparm(3)=-alpha*xtarget(1)

         ! get_foot_from_target
        RR=sqrt(xtarget(1)**2+xtarget(3)**2)
        if (axis_dir.eq.0) then  ! teeth
         radgear=4.0
        else if (axis_dir.eq.1) then  ! no teeth
         radgear=6.0
        else if (axis_dir.eq.2) then
         radgear=7.5
        else
         print *,"axis_dir invalid"
         stop
        endif

         ! get_foot_from_target
        if (RR.gt.radgear/10000.0) then
         call arctan2(xtarget(3),xtarget(1),theta)
         theta=theta+alpha*time  ! clockwise motion
         do while (theta.lt.zero)
          theta=theta+two*Pi
         enddo
         do while (theta.gt.two*Pi)
          theta=theta-two*Pi
         enddo
         xfoot(1)=RR*cos(theta)
         xfoot(3)=RR*sin(theta)
        else
         xfoot(1)=zero
         xfoot(3)=zero
        endif

       ! above 563=gear
       ! below 538,541=diesel injector
       else if ((probtype.eq.538).or.(probtype.eq.541)) then

        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then ! place inlet at center of domain
         xfoot(1)=xfoot(1)+0.0015
         xfoot(2)=xfoot(2)-0.0045
        else
         print *,"levelrz invalid get_foot_from_target"
         stop
        endif 
        
        ! flapping wing (get_foot_from_target)
       else if (probtype.eq.701) then 

        numMotion=2
  
        if ((axis_dir.eq.0).or.(axis_dir.eq.1)) then     
         motionPara(1,1)=0.0         ! 0: rotation; 1: translation
         motionPara(2,1)=30./180.*Pi ! amplitude (thetaMag)
         motionPara(6,1)=0.25        ! pivoting point
         motionPara(10,1)=0.
         motionPara(11,1)=1.
        else if (axis_dir.eq.2) then
          ! theta=hMag*cos(2*Pi*fH*t+phiH)+h0
          ! r=x0+v(3)*theta
         motionPara(1,1)=1.0
         motionPara(2,1)=1.0  ! hMag
         motionPara(6,1)=0.0  ! x0
         motionPara(10,1)=0.  ! v(2) 
         motionPara(11,1)=1.  ! v(3)
        else
         print *,"axis_dir invalid"
         stop
        endif

        motionPara(3,1)=1.          ! f=frequency 2*pi*f
        motionPara(4,1)=0.          ! offset alpha0
        motionPara(5,1)=0./180.*Pi  ! phase change (h0 if translate)
        motionPara(7,1)=0.
        motionPara(8,1)=0.
        motionPara(9,1)=0.  ! rotation direction

        motionPara(1,2)=1    ! motionType 
        motionPara(2,2)=0.0  ! hMag
        motionPara(3,2)=1.   ! hF (frequency)
        motionPara(4,2)=0.   ! hPhi
        motionPara(5,2)=0.   ! h0 
        motionPara(6,2)=0.0  ! x(1)
        motionPara(7,2)=0.   ! x(2)
        motionPara(8,2)=0.   ! x(3)
        motionPara(9,2)=0.   ! v(1)
        motionPara(10,2)=1.  ! v(2) 
        motionPara(11,2)=0.  ! v(3)

        dt_flapping=0.001

        if (time.ge.zero) then
         flapping_time=time-dt_flapping
        else
         print *,"time invalid"
         stop
        endif

          ! 0<t<1 
        call flappingKinematics(numMotion,motionPara,r,flapping_time)
        det=r(1,1)*r(2,2)-r(1,2)*r(2,1)
        if (det.eq.zero) then
         print *,"rotation matrix is singular"
         stop
        endif
         ! Ainv=( a22  -a12
         !        -a21 a11  )/det
         ! A   =( a11  a12
         !        a21  a22 )
        rinv(1,1)=r(2,2)/det
        rinv(1,2)=-r(1,2)/det
        rinv(2,1)=-r(2,1)/det
        rinv(2,2)=r(1,1)/det
         ! get_foot_from_target
        xfoot(1)=xtarget(1)-r(1,4)
        xfoot(3)=xtarget(3)-r(2,4)
        xfootsave(1)=xfoot(1)
        xfootsave(3)=xfoot(3)
        xfoot(1)=rinv(1,1)*xfootsave(1)+rinv(1,2)*xfootsave(3)
        xfoot(3)=rinv(2,1)*xfootsave(1)+rinv(2,2)*xfootsave(3)

        flapping_time_plus=flapping_time+dt_flapping
        call flappingKinematics(numMotion,motionPara,rplus,flapping_time_plus)
        velparm(1)=(rplus(1,1)-r(1,1))*xfoot(1)+ &
                   (rplus(1,2)-r(1,2))*xfoot(3)+ &
                   rplus(1,4)-r(1,4)
        velparm(2)=0.0
        velparm(3)=(rplus(2,1)-r(2,1))*xfoot(1)+ &
                   (rplus(2,2)-r(2,2))*xfoot(3)+ &
                   rplus(2,4)-r(2,4)
        do dir=1,3
         velparm(dir)=velparm(dir)/(2.0*dt_flapping)
        enddo
        
       endif   ! probtype.eq.701

      else if (part_id.eq.2) then

       ! get_foot_from_target

       ! needle for diesel injector (part 2)
       if ((probtype.eq.538).or.(probtype.eq.541)) then

         ! if injector is shifted, then needle should be shifted too.
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then ! place inlet at center of domain
         xfoot(1)=xfoot(1)+0.0015
         xfoot(2)=xfoot(2)-0.0045
        else
         print *,"levelrz invalid get_foot_from_target 2"
         stop
        endif 

        xfoot(3)=xfoot(3)+0.01
        do dir=1,3
         xfoot(dir)=xfoot(dir)-FSI(part_id)%solid_displ(dir)
        enddo

        do dir=1,3
         velparm(dir)=FSI(part_id)%solid_speed(dir)
        enddo

       else if (probtype.eq.701) then  ! 2nd flapping wing

        if ((axis_dir.eq.0).or.(axis_dir.eq.2)) then
         print *,"no part 2 should exist"
         stop
        else if (axis_dir.eq.1) then

         do dir=1,3
          xfoot(dir)=xtarget(dir)
          xfootsave(dir)=xtarget(dir)
          xtargetsave(dir)=xtarget(dir)
         enddo

         numMotion=2
             
         motionPara(1,1)=0.0
         motionPara(2,1)=0./180.*Pi ! was 30./180.*Pi
         motionPara(3,1)=1.
         motionPara(4,1)=0.
         motionPara(5,1)=0./180.*Pi
         motionPara(6,1)=0.0  ! was 0.25
         motionPara(7,1)=0.
         motionPara(8,1)=0.
         motionPara(9,1)=0.
         motionPara(10,1)=0.
         motionPara(11,1)=1.

         motionPara(1,2)=1   ! motionType 
         motionPara(2,2)=0.5 ! hMag (was 0.0)
         motionPara(3,2)=1.  ! hF
         motionPara(4,2)=0.  ! hPhi
         motionPara(5,2)=0.  ! h0 
         motionPara(6,2)=0.0 !x(1) the default translation point = (0,0,0)
         motionPara(7,2)=0.  !x(2) otherwise it shifts to the new location
         motionPara(8,2)=0.  !x(3) (x1,x2,x3)
         motionPara(9,2)=0.  !v(1)
         motionPara(10,2)=1. !v(2) 
         motionPara(11,2)=0. !v(3)

         ! NINT = round to nearest whole number 
         ! for 2nd wing motion, reverse direction of flapping motion.
         !   itime=NINT(time-0.49999999999999)
         !   flapping_time=time+half

         if (time.ge.zero) then
          flapping_time=time
         else
          print *,"time invalid"
          stop
         endif

           ! 0<t<1 
         call flappingKinematics(numMotion,motionPara,r,flapping_time)
         det=r(1,1)*r(2,2)-r(1,2)*r(2,1)
         if (det.eq.zero) then
          print *,"rotation matrix is singular"
          stop
         endif
         rinv(1,1)=r(2,2)/det
         rinv(1,2)=-r(1,2)/det
         rinv(2,1)=-r(2,1)/det
         rinv(2,2)=r(1,1)/det
          ! get_foot_from_target
         xfoot(1)=xtarget(1)-r(1,4)
         xfoot(3)=xtarget(3)-r(2,4)
         xfootsave(1)=xfoot(1)
         xfootsave(3)=xfoot(3)
         xfoot(1)=rinv(1,1)*xfootsave(1)+rinv(1,2)*xfootsave(3)
         xfoot(3)=rinv(2,1)*xfootsave(1)+rinv(2,2)*xfootsave(3)

         dt_flapping=0.01
         flapping_time_plus=flapping_time+dt_flapping
         call flappingKinematics(numMotion,motionPara,rplus,flapping_time_plus)
         velparm(1)=(rplus(1,1)-r(1,1))*xfoot(1)+ &
                    (rplus(1,2)-r(1,2))*xfoot(3)+ &
                    rplus(1,4)-r(1,4)
         velparm(2)=0.0
         velparm(3)=(rplus(2,1)-r(2,1))*xfoot(1)+ &
                    (rplus(2,2)-r(2,2))*xfoot(3)+ &
                     rplus(2,4)-r(2,4)
         do dir=1,3
          velparm(dir)=velparm(dir)/dt_flapping
         enddo

        else
         print *,"axis_dir invalid"
         stop
        endif

       else
        print *,"probtype invalid for part_id==2"
        stop
       endif  

      else if ((part_id.gt.2).and.(part_id.le.TOTAL_NPARTS)) then

       ! do nothing

      else
       print *,"part_id invalid"
       stop
      endif

      return
      end subroutine get_foot_from_target

       ! velparm is initialized with the foot velocity.
       ! if rigid body moves via prescribed mapping, then velparm
       ! is modified.
      subroutine get_target_from_foot(xfoot,xtarget, &
       velparm,time,part_id)
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T part_id
      REAL_T velparm(3)
      REAL_T time
      INTEGER_T dir
      REAL_T xtarget(3)
      REAL_T xfoot(3)
      REAL_T xtargetsave(3)
      REAL_T xfootsave(3)
      REAL_T RPM,alpha,RR,radgear,theta

      INTEGER_T numMotion
      REAL_T r(3,4)
      REAL_T rinv(3,4)
      REAL_T rplus(3,4)
      REAL_T det
      REAL_T motionPara(11,2)
      REAL_T flapping_time
      REAL_T dt_flapping
      REAL_T flapping_time_plus


      if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
       print *,"part_id invalid"
       stop
      endif
      if (FSI(part_id)%partID.ne.part_id) then
       print *,"FSI(part_id)%partID.ne.part_id"
       stop
      endif

       ! get_target_from_foot
      do dir=1,3
       xtarget(dir)=xfoot(dir)
       xfootsave(dir)=xtarget(dir)
       xtargetsave(dir)=xtarget(dir)
      enddo

      if (part_id.eq.1) then

       ! gear rotating clockwise
       ! theta(t)=theta_0-alpha t
       ! x=r cos(theta)
       ! y=r sin(theta)
       ! x'=alpha y
       ! y'=-alpha x
       ! xfoot=r cos(theta_0)=r cos(theta+alpha t)=
       !    r(costheta cos(alpha t)-sintheta sin(alpha t) )
       ! yfoot=r sin(theta_0)=r sin(theta+alpha t)=
       !    r(sintheta cos(alpha t)+costheta sin(alpha t))

       if (probtype.eq.563) then
        RPM=abs(vinletgas);
        alpha=two*Pi*RPM/60.0  ! radians/s

        ! at a distance "r" from the origin, a particle on the gear will
        ! travel 2 pi r distance in one revolution (2 pi radians)
        ! the time for one revolution is 2 pi/alpha
        ! so magnitude of velocity is 2 pi r /  ( 2pi /alpha )=alpha r
 
         ! get_target_from_foot
        RR=sqrt(xfoot(1)**2+xfoot(3)**2)
        if (axis_dir.eq.0) then  ! teeth
         radgear=4.0
        else if (axis_dir.eq.1) then  ! no teeth
         radgear=6.0
        else if (axis_dir.eq.2) then
         radgear=7.5
        else
         print *,"axis_dir invalid"
         stop
        endif

         ! get_target_from_foot
        if (RR.gt.radgear/10000.0) then
         call arctan2(xfoot(3),xfoot(1),theta)
         theta=theta-alpha*time  ! clockwise motion
         do while (theta.lt.zero)
          theta=theta+two*Pi
         enddo
         do while (theta.gt.two*Pi)
          theta=theta-two*Pi
         enddo
         xtarget(1)=RR*cos(theta)
         xtarget(3)=RR*sin(theta)
        else
         xtarget(1)=zero
         xtarget(3)=zero
        endif
        velparm(1)=alpha*xtarget(3)
        velparm(2)=zero
        velparm(3)=-alpha*xtarget(1)

       ! above 563=gear
       ! below 538,541=diesel injector
       else if ((probtype.eq.538).or. & ! inputs.injA
                (probtype.eq.541)) then

        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then ! place inlet at center of domain
         xtarget(1)=xtarget(1)-0.0015
         xtarget(2)=xtarget(2)+0.0045
        else
         print *,"levelrz invalid get_target_from_foot"
         stop
        endif 
        do dir=1,3
         velparm(dir)=zero
        enddo
        
        ! flapping wing (get_foot_from_target)
       else if (probtype.eq.701) then 

        numMotion=2
  
        if ((axis_dir.eq.0).or.(axis_dir.eq.1)) then     
         motionPara(1,1)=0.0         ! 0: rotation; 1: translation
         motionPara(2,1)=30./180.*Pi ! amplitude (thetaMag)
         motionPara(6,1)=0.25        ! pivoting point
         motionPara(10,1)=0.
         motionPara(11,1)=1.
        else if (axis_dir.eq.2) then
          ! theta=hMag*cos(2*Pi*fH*t+phiH)+h0
          ! r=x0+v(3)*theta
         motionPara(1,1)=1.0
         motionPara(2,1)=1.0  ! hMag
         motionPara(6,1)=0.0  ! x0
         motionPara(10,1)=0.  ! v(2) 
         motionPara(11,1)=1.  ! v(3)
        else
         print *,"axis_dir invalid"
         stop
        endif

        motionPara(3,1)=1.          ! f=frequency 2*pi*f
        motionPara(4,1)=0.          ! offset alpha0
        motionPara(5,1)=0./180.*Pi  ! phase change (h0 if translate)
        motionPara(7,1)=0.
        motionPara(8,1)=0.
        motionPara(9,1)=0.  ! rotation direction

        motionPara(1,2)=1    ! motionType 
        motionPara(2,2)=0.0  ! hMag
        motionPara(3,2)=1.   ! hF (frequency)
        motionPara(4,2)=0.   ! hPhi
        motionPara(5,2)=0.   ! h0 
        motionPara(6,2)=0.0  ! x(1)
        motionPara(7,2)=0.   ! x(2)
        motionPara(8,2)=0.   ! x(3)
        motionPara(9,2)=0.   ! v(1)
        motionPara(10,2)=1.  ! v(2) 
        motionPara(11,2)=0.  ! v(3)

        dt_flapping=0.001

        if (time.ge.zero) then
         flapping_time=time-dt_flapping
        else
         print *,"time invalid"
         stop
        endif

          ! 0<t<1 
        call flappingKinematics(numMotion,motionPara,r,flapping_time)
        det=r(1,1)*r(2,2)-r(1,2)*r(2,1)
        if (det.eq.zero) then
         print *,"rotation matrix is singular"
         stop
        endif
         ! Ainv=( a22  -a12
         !        -a21 a11  )/det
         ! A   =( a11  a12
         !        a21  a22 )
        rinv(1,1)=r(2,2)/det
        rinv(1,2)=-r(1,2)/det
        rinv(2,1)=-r(2,1)/det
        rinv(2,2)=r(1,1)/det
         ! get_target_from_foot
        xtargetsave(1)=xfoot(1)
        xtargetsave(3)=xfoot(3)
        xtarget(1)=r(1,1)*xtargetsave(1)+r(1,2)*xtargetsave(3)
        xtarget(3)=r(2,1)*xtargetsave(1)+r(2,2)*xtargetsave(3)
        xtarget(1)=xtarget(1)+r(1,4)
        xtarget(3)=xtarget(3)+r(2,4)

        flapping_time_plus=flapping_time+dt_flapping
        call flappingKinematics(numMotion,motionPara,rplus,flapping_time_plus)
        velparm(1)=(rplus(1,1)-r(1,1))*xfoot(1)+ &
                   (rplus(1,2)-r(1,2))*xfoot(3)+ &
                   rplus(1,4)-r(1,4)
        velparm(2)=0.0
        velparm(3)=(rplus(2,1)-r(2,1))*xfoot(1)+ &
                   (rplus(2,2)-r(2,2))*xfoot(3)+ &
                   rplus(2,4)-r(2,4)
        do dir=1,3
         velparm(dir)=velparm(dir)/(2.0*dt_flapping)
        enddo
       
       else ! above: probtype==701

        ! do nothing: velparm prescribed on input.
 
       endif   

      else if (part_id.eq.2) then

       ! get_target_from_foot

       ! needle for diesel injector (part 2)
       if ((probtype.eq.538).or. & ! inputs.injA
           (probtype.eq.541)) then

         ! if injector is shifted, then needle should be shifted too.
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then ! place inlet at center of domain
         xtarget(1)=xtarget(1)-0.0015
         xtarget(2)=xtarget(2)+0.0045
        else
         print *,"levelrz invalid get_target_from_foot 2"
         stop
        endif 

        xtarget(3)=xtarget(3)-0.01
        do dir=1,3
         xtarget(dir)=xtarget(dir)+FSI(part_id)%solid_displ(dir)
        enddo

        do dir=1,3
         velparm(dir)=FSI(part_id)%solid_speed(dir)
        enddo

       else if (probtype.eq.701) then  ! 2nd flapping wing

        if ((axis_dir.eq.0).or.(axis_dir.eq.2)) then
         print *,"no part 2 should exist"
         stop
        else if (axis_dir.eq.1) then

         do dir=1,3
          xfoot(dir)=xtarget(dir)
          xfootsave(dir)=xtarget(dir)
          xtargetsave(dir)=xtarget(dir)
         enddo

         numMotion=2
             
         motionPara(1,1)=0.0
         motionPara(2,1)=0./180.*Pi ! was 30./180.*Pi
         motionPara(3,1)=1.
         motionPara(4,1)=0.
         motionPara(5,1)=0./180.*Pi
         motionPara(6,1)=0.0  ! was 0.25
         motionPara(7,1)=0.
         motionPara(8,1)=0.
         motionPara(9,1)=0.
         motionPara(10,1)=0.
         motionPara(11,1)=1.

         motionPara(1,2)=1   ! motionType 
         motionPara(2,2)=0.5 ! hMag (was 0.0)
         motionPara(3,2)=1.  ! hF
         motionPara(4,2)=0.  ! hPhi
         motionPara(5,2)=0.  ! h0 
         motionPara(6,2)=0.0 !x(1) the default translation point = (0,0,0)
         motionPara(7,2)=0.  !x(2) otherwise it shifts to the new location
         motionPara(8,2)=0.  !x(3) (x1,x2,x3)
         motionPara(9,2)=0.  !v(1)
         motionPara(10,2)=1. !v(2) 
         motionPara(11,2)=0. !v(3)

         ! NINT = round to nearest whole number 
         ! for 2nd wing motion, reverse direction of flapping motion.
         !   itime=NINT(time-0.49999999999999)
         !   flapping_time=time+half

         if (time.ge.zero) then
          flapping_time=time
         else
          print *,"time invalid"
          stop
         endif

           ! 0<t<1 
         call flappingKinematics(numMotion,motionPara,r,flapping_time)
         det=r(1,1)*r(2,2)-r(1,2)*r(2,1)
         if (det.eq.zero) then
          print *,"rotation matrix is singular"
          stop
         endif
         rinv(1,1)=r(2,2)/det
         rinv(1,2)=-r(1,2)/det
         rinv(2,1)=-r(2,1)/det
         rinv(2,2)=r(1,1)/det
          ! get_target_from_foot
         xtargetsave(1)=xtarget(1)
         xtargetsave(3)=xtarget(3)
         xtarget(1)=r(1,1)*xtargetsave(1)+r(1,2)*xtargetsave(3)
         xtarget(3)=r(2,1)*xtargetsave(1)+r(2,2)*xtargetsave(3)
         xtarget(1)=xtarget(1)+r(1,4)
         xtarget(3)=xtarget(3)+r(2,4)

         dt_flapping=0.01
         flapping_time_plus=flapping_time+dt_flapping
         call flappingKinematics(numMotion,motionPara,rplus,flapping_time_plus)
         velparm(1)=(rplus(1,1)-r(1,1))*xfoot(1)+ &
                    (rplus(1,2)-r(1,2))*xfoot(3)+ &
                    rplus(1,4)-r(1,4)
         velparm(2)=0.0
         velparm(3)=(rplus(2,1)-r(2,1))*xfoot(1)+ &
                    (rplus(2,2)-r(2,2))*xfoot(3)+ &
                     rplus(2,4)-r(2,4)
         do dir=1,3
          velparm(dir)=velparm(dir)/dt_flapping
         enddo

        else
         print *,"axis_dir invalid"
         stop
        endif
 
       else
        ! do nothing: velparm prescribed on input.
       endif  

      else if ((part_id.gt.2).and.(part_id.le.TOTAL_NPARTS)) then
       ! do nothing
      else
       print *,"part_id invalid"
       stop
      endif

      return
      end subroutine get_target_from_foot


subroutine find_grid_bounding_box( &
 part_id, &
 null_intersection, &
 minnode,maxnode, &
 FSI_lo,FSI_hi, &
 FSI_growlo,FSI_growhi,  &
 growlo3D,growhi3D,  &
 xdata3D, &
 DIMS3D(xdata3D), &
 gridlo,gridhi,dxBB)
IMPLICIT NONE

 INTEGER_T part_id
 INTEGER_T null_intersection
 REAL_T minnode(3),maxnode(3)
 INTEGER_T FSI_lo(3),FSI_hi(3)
 INTEGER_T FSI_growlo(3),FSI_growhi(3)
 INTEGER_T growlo3D(3),growhi3D(3)
 INTEGER_T DIMDEC3D(xdata3D)
 REAL_T xdata3D(DIMV3D(xdata3D),3)
 INTEGER_T gridlo(3),gridhi(3)
 INTEGER_T dir
 INTEGER_T ii,jj,kk
 INTEGER_T i,j,k,incr,iter
 REAL_T xcontrol,xcost
 REAL_T dxBB(3)
 REAL_T xlo(3),xhi(3)
 INTEGER_T ngrow,ngrowtest
 INTEGER_T local_iband

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif

 local_iband=FSI(part_id)%bounding_box_ngrow
 if (local_iband.ne.3) then
  print *,"local_iband invalid"
  stop
 endif

 ngrow=0
 do dir=1,3
  ngrowtest=FSI_lo(dir)-FSI_growlo(dir)
  if ((ngrowtest.lt.0).or.(ngrowtest.gt.4)) then
   print *,"ngrowtest invalid1 ",ngrowtest
   stop
  endif
  if (ngrowtest.gt.ngrow) then
   ngrow=ngrowtest
  endif
  ngrowtest=FSI_growhi(dir)-FSI_hi(dir)
  if ((ngrowtest.lt.0).or.(ngrowtest.gt.4)) then
   print *,"ngrowtest invalid2 ",ngrowtest
   stop
  endif
  if (ngrowtest.gt.ngrow) then
   ngrow=ngrowtest
  endif
 enddo ! dir=1..3

 call checkbound3D(FSI_lo,FSI_hi, &
  DIMS3D(xdata3D), &
  ngrow,-1,123)

 do dir=1,3
  ii=0
  jj=0
  kk=0
  if (dir.eq.1) then
   ii=1
  else if (dir.eq.2) then
   jj=1
  else if (dir.eq.3) then
   kk=1
  else
   print *,"dir invalid"
   stop
  endif 
  i=FSI_lo(1)
  j=FSI_lo(2)
  k=FSI_lo(3)
  xlo(dir)=half*(xdata3D(i,j,k,dir)+xdata3D(i-ii,j-jj,k-kk,dir))
  i=FSI_hi(1)
  j=FSI_hi(2)
  k=FSI_hi(3)
  xhi(dir)=half*(xdata3D(i,j,k,dir)+xdata3D(i+ii,j+jj,k+kk,dir))
  dxBB(dir)=(xhi(dir)-xlo(dir))/(FSI_hi(dir)-FSI_lo(dir)+1)
 enddo ! dir=1..3

! NINT = round to nearest whole number
! x=(i-lo)dx+xlo+dx/2=(i-lo+1/2)dx+xlo
! i=(x-xlo)/dx+lo-1/2

 null_intersection=0

 do dir=1,3

  if (null_intersection.eq.0) then

   ii=0
   jj=0
   kk=0
   if (dir.eq.1) then
    ii=1
   else if (dir.eq.2) then
    jj=1
   else if (dir.eq.3) then
    kk=1
   else
    print *,"dir invalid"
    stop
   endif 

   i=FSI_lo(1)
   j=FSI_lo(2)
   k=FSI_lo(3)

   gridlo(dir)=NINT( (minnode(dir)-xlo(dir))/dxBB(dir)-half+FSI_lo(dir) )
   if (gridlo(dir).lt.growlo3D(dir)) then
    gridlo(dir)=growlo3D(dir)
   endif
   if (gridlo(dir).gt.growhi3D(dir)) then
    gridlo(dir)=growhi3D(dir)
    null_intersection=1
   endif

   if (null_intersection.eq.0) then

    incr=gridlo(dir)-FSI_lo(dir)
    xcontrol=xdata3D(i+ii*incr,j+jj*incr,k+kk*incr,dir)
    xcost=minnode(dir)-local_iband*dxBB(dir)
    iter=0
    do while ((gridlo(dir).gt.growlo3D(dir)).and. &
              (xcontrol.gt.xcost))
     gridlo(dir)=gridlo(dir)-1
     incr=gridlo(dir)-FSI_lo(dir)
     xcontrol=xdata3D(i+ii*incr,j+jj*incr,k+kk*incr,dir)
     iter=iter+1
     if (iter.gt.FSI_growhi(dir)-FSI_growlo(dir)+1) then
      print *,"failure to find xcontrol"
      stop
     endif
    enddo
  
    gridhi(dir)=NINT( (maxnode(dir)-xlo(dir))/dxBB(dir)-half+FSI_lo(dir) )
    if (gridhi(dir).gt.growhi3D(dir)) then
     gridhi(dir)=growhi3D(dir)
    endif
    if (gridhi(dir).lt.growlo3D(dir)) then
     gridhi(dir)=growlo3D(dir)
     null_intersection=1
    endif

    if (null_intersection.eq.0) then

     incr=gridhi(dir)-FSI_lo(dir)
     xcontrol=xdata3D(i+ii*incr,j+jj*incr,k+kk*incr,dir)
     xcost=maxnode(dir)+local_iband*dxBB(dir)
     iter=0
     do while ((gridhi(dir).lt.growhi3D(dir)).and. &
               (xcontrol.lt.xcost))
      gridhi(dir)=gridhi(dir)+1
      incr=gridhi(dir)-FSI_lo(dir)
      xcontrol=xdata3D(i+ii*incr,j+jj*incr,k+kk*incr,dir)
      iter=iter+1
      if (iter.gt.FSI_growhi(dir)-FSI_growlo(dir)+1) then
       print *,"failure to find xcontrol"
       stop
      endif
     enddo

    else if (null_intersection.eq.1) then
     ! do nothing
    else
     print *,"null_intersection invalid"
     stop
    endif
   else if (null_intersection.eq.1) then
    ! do nothing
   else
    print *,"null_intersection invalid"
    stop
   endif
  else if (null_intersection.eq.1) then
   ! do nothing
  else
   print *,"null_intersection invalid"
   stop
  endif

 enddo  ! dir=1..sdim

return
end subroutine find_grid_bounding_box


subroutine find_grid_bounding_box_node( &
 xnot, &
 FSI_lo,FSI_hi, &
 FSI_growlo,FSI_growhi,  &
 growlo3D,growhi3D,  &
 xdata3D, &
 DIMS3D(xdata3D), &
 gridlo,gridhi,dxBB)
IMPLICIT NONE

 REAL_T xnot(3)
 INTEGER_T FSI_lo(3),FSI_hi(3)
 INTEGER_T FSI_growlo(3),FSI_growhi(3)
 INTEGER_T growlo3D(3),growhi3D(3)
 INTEGER_T DIMDEC3D(xdata3D)
 REAL_T xdata3D(DIMV3D(xdata3D),3)
 INTEGER_T gridlo(3),gridhi(3)
 INTEGER_T dir,dirloc
 INTEGER_T idx(3),idxL(3),idxR(3)
 INTEGER_T iter,change
 REAL_T dist,distL,distR
 REAL_T dxBB(3)
 REAL_T xlo(3),xhi(3)
 INTEGER_T ngrow,ngrowtest
 INTEGER_T interp_support

 interp_support=2

 ngrow=0
 do dir=1,3
  ngrowtest=FSI_lo(dir)-FSI_growlo(dir)
  if ((ngrowtest.lt.0).or.(ngrowtest.gt.4)) then
   print *,"ngrowtest invalid1 ",ngrowtest
   stop
  endif
  if (ngrowtest.gt.ngrow) then
   ngrow=ngrowtest
  endif
  ngrowtest=FSI_growhi(dir)-FSI_hi(dir)
  if ((ngrowtest.lt.0).or.(ngrowtest.gt.4)) then
   print *,"ngrowtest invalid2 ",ngrowtest
   stop
  endif
  if (ngrowtest.gt.ngrow) then
   ngrow=ngrowtest
  endif
 enddo ! dir=1..3

 call checkbound3D(FSI_lo,FSI_hi, &
  DIMS3D(xdata3D), &
  ngrow,-1,123)

 do dir=1,3
  do dirloc=1,3
   idxL(dirloc)=FSI_lo(dirloc)
   idxR(dirloc)=FSI_lo(dirloc)
  enddo
  idxL(dir)=idxL(dir)-1
  xlo(dir)=half*(xdata3D(idxL(1),idxL(2),idxL(3),dir)+ &
                 xdata3D(idxR(1),idxR(2),idxR(3),dir))

  do dirloc=1,3
   idxL(dirloc)=FSI_hi(dirloc)
   idxR(dirloc)=FSI_hi(dirloc)
  enddo
  idxR(dir)=idxR(dir)+1
  xhi(dir)=half*(xdata3D(idxL(1),idxL(2),idxL(3),dir)+ &
                 xdata3D(idxR(1),idxR(2),idxR(3),dir))
  dxBB(dir)=(xhi(dir)-xlo(dir))/(FSI_hi(dir)-FSI_lo(dir)+1)
 enddo ! dir=1..3

! NINT = round to nearest whole number
! x=(i-lo)dx+xlo+dx/2=(i-lo+1/2)dx+xlo
! i=(x-xlo)/dx+lo-1/2

 do dir=1,3

  if (xnot(dir).lt.xlo(dir)-VOFTOL*dxBB(dir)) then
   print *,"node should be within grid interior"
   stop
  endif
  if (xnot(dir).gt.xhi(dir)+VOFTOL*dxBB(dir)) then
   print *,"node should be within grid interior"
   stop
  endif
  if (xnot(dir).le.xlo(dir)+VOFTOL*dxBB(dir)) then
   gridlo(dir)=FSI_lo(dir)
  else if (xnot(dir).ge.xhi(dir)-VOFTOL*dxBB(dir)) then
   gridlo(dir)=FSI_hi(dir)
  else if ((xnot(dir).ge.xlo(dir)).and. &
           (xnot(dir).le.xhi(dir))) then
   gridlo(dir)=NINT( (xnot(dir)-xlo(dir))/dxBB(dir)-half+FSI_lo(dir) )
   if ((gridlo(dir).lt.FSI_lo(dir)).or. &
       (gridlo(dir).gt.FSI_hi(dir))) then
    print *,"node should be within grid interior"
    stop
   endif
   do dirloc=1,3
    idx(dirloc)=FSI_lo(dirloc)
   enddo
   idx(dir)=gridlo(dir)
   dist=abs(xdata3D(idx(1),idx(2),idx(3),dir)-xnot(dir)) 
   change=1
   iter=0
   do while (change.eq.1)
    do dirloc=1,3
     idxL(dirloc)=idx(dirloc)
     idxR(dirloc)=idx(dirloc)
    enddo
    idxL(dir)=idx(dir)-1
    idxR(dir)=idx(dir)+1
    distL=abs(xdata3D(idxL(1),idxL(2),idxL(3),dir)-xnot(dir))
    distR=abs(xdata3D(idxR(1),idxR(2),idxR(3),dir)-xnot(dir))
    if ((distL.le.dist).and.(distR.le.dist)) then
     print *,"distL or distR invalid"
     stop
    endif
    if (distL.lt.dist) then
     change=1
     dist=distL
     idx(dir)=idx(dir)-1
    else if (distR.lt.dist) then
     change=1
     dist=distR
     idx(dir)=idx(dir)+1
    else if ((distL.ge.dist).and. &
             (distR.ge.dist)) then
     change=0
    else
     print *,"distL or distR invalid"
     stop
    endif
    iter=iter+1
    if (iter.gt.FSI_hi(dir)-FSI_lo(dir)+1) then
     print *,"iter.gt.FSI_hi(dir)-FSI_lo(dir)+1"
     stop
    endif
   enddo ! while (change==1)
   gridlo(dir)=idx(dir)
  else
   print *,"xnot(dir) invalid"
   stop
  endif

  gridhi(dir)=gridlo(dir)+interp_support 
  gridlo(dir)=gridlo(dir)-interp_support 

 enddo  ! dir=1..sdim

return
end subroutine find_grid_bounding_box_node


! called from FORT_HEADERMSG with FSI_operation=1
! isout==1 => verbose
subroutine CLSVOF_ReadNodes( &
  FSI_refine_factor, &
  FSI_bounding_box_ngrow, &
  CLSVOF_curtime,CLSVOF_dt, &
  h_small, &
  problo,probhi,current_step,plot_interval,ioproc,isout)
  use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

  INTEGER_T :: current_step,plot_interval
  INTEGER_T :: initflag,ioproc,isout,part_id
  REAL_T :: CLSVOF_curtime,CLSVOF_dt
  REAL_T :: h_small
  REAL_T problo(3),probhi(3)
  INTEGER_T nmat
  INTEGER_T node_factor 
  INTEGER_T ctml_part_id 
  INTEGER_T fsi_part_id 
  INTEGER_T :: inode_crit,inode
  INTEGER_T FSI_refine_factor(num_materials)
  INTEGER_T FSI_bounding_box_ngrow(num_materials)
  INTEGER_T im_sanity_check

  nmat=num_materials

  do im_sanity_check=1,nmat
   if ((FSI_refine_factor(im_sanity_check).lt.0).or. &
       (FSI_refine_factor(im_sanity_check).gt.100)) then
    print *,"FSI_refine_factor(im_sanity_check) invalid"
    stop
   endif
   if (FSI_bounding_box_ngrow(im_sanity_check).ne.3) then
    print *,"FSI_bounding_box_ngrow(im_sanity_check) invalid"
    stop
   endif
  enddo

  if (h_small.le.zero) then
   print *,"h_small invalid"
   stop
  endif
  if (current_step.lt.0) then
   print *,"current_step invalid"
   stop
  endif
  if (plot_interval.lt.-1) then
   print *,"plot_interval invalid"
   stop
  endif

  if ((TOTAL_NPARTS.ge.1).and.(TOTAL_NPARTS.le.MAX_PARTS)) then

   if (CTML_FSI_flagF(nmat).eq.1) then ! FSI_flag==4
#ifdef MVAHABFSI
    call CTML_SOLVE_SOLID( &
     CLSVOF_curtime, &
     CLSVOF_dt, &
     current_step, &
     isout, &
     plot_interval, &
     ioproc)
#else
    print *,"define MVAHABFSI"
    stop
#endif
   else if (CTML_FSI_flagF(nmat).eq.0) then
    ! do nothing
   else
    print *,"CTML_FSI_flagF(nmat) invalid"
    stop
   endif

   do part_id=1,TOTAL_NPARTS

    FSI(part_id)%PartID=part_id

    ctml_part_id=CTML_partid_map(part_id)
    fsi_part_id=FSI_partid_map(part_id)

    if (((ctml_part_id.ge.1).and. &
         (ctml_part_id.le.CTML_NPARTS)).or. &
        ((fsi_part_id.ge.1).and. &
         (fsi_part_id.le.FSI_NPARTS))) then

     if ((ctml_part_id.gt.0).and. &
         (ctml_part_id.le.CTML_NPARTS)) then

      if (AMREX_SPACEDIM.eq.3) then
       node_factor=1
       print *,"3D not supported yet"
      else if (AMREX_SPACEDIM.eq.2) then
       node_factor=2
      else
       print *,"dimension bust"
       stop
      endif
      if (FSI(part_id)%NumNodes.ne. &
          ctml_n_fib_active_nodes(ctml_part_id)*node_factor) then
       print *,"NumNodes is corrupt"
       stop
      endif
      if (FSI(part_id)%NumIntElems.ne. &
          (ctml_n_fib_active_nodes(ctml_part_id)-1)*node_factor) then
       print*,"NumIntElems is corrupt"
       stop
      endif
      if (FSI(part_id)%IntElemDim.ne.3) then
       print *,"FSI(part_id)%IntElemDim.ne.3"
       stop
      endif

      if ((ctml_part_id.ge.1).and. &
          (ctml_part_id.le.CTML_NPARTS)) then

       inode_crit=0 ! node index of first inactive node.
       inode=0
#ifdef MVAHABFSI
       call CTML_GET_POS_VEL_FORCE_WT( &
        ctml_fib_pst, &
        ctml_fib_vel, &
        ctml_fib_frc, &
        ctml_fib_mass, &
        ctml_n_fib_bodies, &
        ctml_max_n_fib_nodes, &
        ctml_part_id)

       do inode=1,ctml_n_fib_nodes(ctml_part_id)
        if (ctml_fib_mass(ctml_part_id,inode).gt.zero) then
         ! do nothing
        else if (ctml_fib_mass(ctml_part_id,inode).eq.zero) then 
         if (inode_crit.eq.0) then
          inode_crit=inode
         endif
        else
         print *,"ctml_fib_mass(ctml_part_id,inode) invalid"
         stop
        endif
       enddo ! inode=1,ctml_n_fib_nodes(ctml_part_id)

       if (inode_crit.eq.0) then
        ctml_n_fib_active_nodes(ctml_part_id)= &
         ctml_n_fib_nodes(ctml_part_id)
       else if (inode_crit.gt.1) then
        print *,"WARNING inode_crit>1  inode_crit=",inode_crit
         ctml_n_fib_active_nodes(ctml_part_id)=inode_crit-1
       else
        print *,"inode_crit invalid"
        stop
       endif
#else
       print *,"define MVAHABFSI"
       stop
#endif
      else
       print *,"ctml_part_id invalid"
       stop
      endif

      call CTML_init_sci_node(ioproc,part_id,isout)
     else if (ctml_part_id.eq.0) then
      ! do nothing
     else
      print *,"ctml_part_id invalid" 
      stop
     endif

     if ((fsi_part_id.ge.1).and. &
         (fsi_part_id.le.FSI_NPARTS)) then
      ! ifirst=0 (iread=0) internal to overall_solid_advance.
      ! since ifirst=0 (iread=0) when initinjector called (if probtype=538),
      ! code initializes FSI%solid_displ and FSI%solid_speed
      ! if probtype=701 then initflapping called, but nothing done.
      ! overall_solid_advance calls e.g. whale_geominit, or geominit, or
      !  initinjector, etc
      call overall_solid_advance(CLSVOF_curtime,CLSVOF_dt,part_id, &
       ioproc,isout)
     else if (fsi_part_id.eq.0) then
      ! do nothing
     else
      print *,"fsi_part_id invalid"
      stop
     endif

     initflag=0

     if (FSI(part_id)%deforming_part.eq.0) then
      ! do nothing
     else if (FSI(part_id)%deforming_part.eq.1) then
      call generate_new_triangles(initflag,problo,probhi, &
       part_id,ioproc,isout,h_small)
     else
      print *,"FSI(part_id)%deforming_part invalid"
      stop
     endif

    else if ((ctml_part_id.eq.0).and. &
             (fsi_part_id.eq.0)) then
     ! do nothing
    else
     print *,"ctml_part_id or fsi_part_id invalid"
     stop
    endif

   enddo ! part_id=1,TOTAL_NPARTS

  else
   print *,"TOTAL_NPARTS invalid: ",TOTAL_NPARTS
   stop
  endif
 
return
end subroutine CLSVOF_ReadNodes

!!!!!!!! Austen Duffy's routines here !!!!!!
! probtype=562 (whale)
! set axis_dir=4 (animated whalenormal.txt)

      subroutine getinfo(whalein,Nodes,whaleout)
      IMPLICIT NONE 
      
      Character (len=*) whalein, whaleout
      INTEGER_T Nodes
      
      whaleout=whalein
      
      print *, whaleout
      
      open(unit=4,file=whalein)
      
      read(4,*) Nodes
      
      close(4)
      
      END subroutine getinfo
   
      subroutine runonce(R,S,T,U,V,W,X,Y,Z,List,angle,timestep,spring, &
                   counter,counter_real,Nodes,Cells, &
                   X_init,Y_init,Z_init,whaleout)


!     creates files list.dat and connect.txt
!
!     ALL arguments are output variables used by get_new_geometry. 
!
!     R,S,T -> accel. , U,V,W -> vel. , X,Y,Z -> pos.
!
!     List -> data structure that stores node connectivity
!
!     angle and timestep -> arrays for tail angle data
!
!     spring -> array of spring constants 
!
!     DT_DUFFY * STEPS_DUFFY = 8


      IMPLICIT NONE

      CHARACTER*35 whaleout
      INTEGER Nodes, Cells, Shape
      INTEGER_T, DIMENSION(Nodes,20) :: List
      INTEGER i,  k, n1, n2, n3, garbage, counter
      REAL_T, Dimension(Nodes) ::  X, Y, Z, R, S, T, X_init, Y_init
      REAL_T, Dimension(Nodes) ::  spring, Z_init
      REAL_T, Dimension(Nodes) ::  U, V, W
      REAL_T, Dimension(22) :: angle, timestep
      REAL_T  a, b, c, value,temp
      REAL_T  counter_real, tailpos, L1,L2,L4,L5
      REAL_T  xtailpos

      tailpos=7.77
      xtailpos=0.5

      
      open(unit=2, file= whaleout)
      
      read(2,*) Nodes

      close(2)
      

      call nodelist(Nodes, 20, List,whaleout)

      open(unit=20, file='list.dat')

      DO i=1,Nodes
         write(20,*) (List(i,k), k=1,20)
      END DO

      Close(20)

      open(unit=2, file= whaleout)      
      
      read(2,*) Nodes
      read(2,*) Cells
      read(2,*) Shape

      

      DO i=1,Nodes
         read(2,*) garbage, a, b, c
         X(i)=a
         Y(i)=b
         Z(i)=c
         X_init(i)=a
         Y_init(i)=b
         Z_init(i)=c
      END DO   


      open(unit=16, file='connect.txt')

      DO i=1,3
         read(2,*) garbage
      END DO
  
      FSI(1)%IntElemDim=3 
      FSI(1)%NumIntElems=Cells
      allocate(FSI(1)%IntElem(FSI(1)%IntElemDim,FSI(1)%NumIntElems))

      DO i=1,Cells-1
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
         write(16,*) n1, n2, n3

         FSI(1)%IntElem(1,i)=n1
         FSI(1)%IntElem(2,i)=n2
         FSI(1)%IntElem(3,i)=n3

         read(2,*) garbage
         read(2,*) garbage
      END DO 

      read(2,*) n1
      read(2,*) n2
      read(2,*) n3
      write(16,*) n1, n2, n3

      FSI(1)%IntElem(1,Cells)=n1
      FSI(1)%IntElem(2,Cells)=n2
      FSI(1)%IntElem(3,Cells)=n3

      close(unit=16)         
      close(2)
 
      open(unit=18, file='angles3.txt')

      DO i=1,21
         read(18,*) a
         angle(i)=a*1.5
         timestep(i)=(i-1.0)*(STEPS_DUFFY/20.0)*4.0
      END DO

      close(18)


      call tailup(R,S,T,U,V,W,X,Y,Z,Nodes,tailpos)


      call springs(Z, spring, Nodes, tailpos)


       temp=36.28*0.0174539252/2.0

  
        DO i=1,Nodes

        If (Z_init(i).gt.Z_init(7355)) Then


        If (value.le.0.0) Then


          L1=sqrt((Z_init(i)-Z_init(7355))**2+(Y_init(7355)-Y_init(i))**2)


             L4=L1*cos(temp)

             
             L5=sqrt(L1**2-L4**2)


             L2=sqrt(2.0*L5**2)




             Y(i)=Y_init(i)-L2+(Y(7355)-Y_init(7355))
             Z(i)=Z_init(i)-L2+(Z(7355)-Z_init(7355))


        End If

        If (value.gt.0.0) Then

         L1=sqrt((Y_init(i)-Y_init(7355))**2+(Z_init(i)-Z_init(7355) )**2)

             L4=L1*cos(temp)
 
             L5=sqrt(L1**2-L4**2)

             L2=sqrt(2.0*L5**2)

             Y(i)=Y_init(i)+L2+(Y(7355)-Y_init(7355))
             Z(i)=Z_init(i)-L2+(Z(7355)-Z_init(7355))



        End If


        END IF

        END DO
        
        counter=1
        counter_real=1.0
   
      End Subroutine
      
      
   
      subroutine new_geometry(R,S,T,U,V,W,X,Y,Z,List,angle,timestep, &
                    spring,counter,counter_real,Nodes,Cells, &
                    X_init,Y_init,Z_init)

      IMPLICIT NONE

      INTEGER Nodes, Cells
      INTEGER_T, DIMENSION(Nodes,20) :: List
      INTEGER i, j, p, counter
      INTEGER_T count1, count2, count3,count4
      REAL_T, Dimension(Nodes) ::  X, Y, Z, R, S, T, X_init, Y_init
      REAL_T, Dimension(Nodes) ::  spring, Z_init
      REAL_T, Dimension(Nodes) ::  U, V, W
      REAL_T, Dimension(22) :: angle, timestep
      REAL_T  sum_x, sum_y, sum_z, mu, nu, value,temp
      REAL_T  counter_real, tailpos, mag, mag_init,L1,L2,L4,L5

      tailpos=7.77

          call plininterp(timestep, angle, counter_real, 21, value)
   
             temp=abs(value)*0.0174539252/2.0

  
        DO i=1,Nodes

        If (Z_init(i).gt.Z_init(7355)) Then


        If (value.le.0.0) Then



          L1=sqrt((Z_init(i)-Z_init(7355))**2+(Y_init(7355)-Y_init(i) )**2)


             L4=L1*cos(temp)

             
             L5=sqrt(L1**2-L4**2)


             L2=sqrt(2.0*L5**2)


             Y(i)=Y_init(i)-L2+(Y(7355)-Y_init(7355))
             Z(i)=Z_init(i)-L2+(Z(7355)-Z_init(7355))


          End If

        If (value.gt.0.0) Then

         L1=sqrt((Y_init(i)-Y_init(7355))**2+(Z_init(i)-Z_init(7355) )**2)

             L4=L1*cos(temp)
 
             L5=sqrt(L1**2-L4**2)

             L2=sqrt(2.0*L5**2)

             Y(i)=Y_init(i)+L2+(Y(7355)-Y_init(7355))
             Z(i)=Z_init(i)-L2+(Z(7355)-Z_init(7355))


          End If

        END IF

        END DO

       mu=0.05
       nu=0.0
       
        DO i=1,Nodes

          If (Z(i).lt.tailpos) Then
           sum_x=0.0
           sum_y=0.0
           sum_z=0.0
           DO j=1,20
              p=List(i,j)
              If (p.ne.0) Then
                mag=sqrt((X(p)-X(i))**2+(Y(p)-Y(i))**2+(Z(p)-Z(i))**2)
                mag_init=sqrt((X_init(p)-X_init(i))**2+ &
                 (Y_init(p)-Y_init(i))**2+(Z_init(p)-Z_init(i))**2)
                sum_x=sum_x+mu*(U(p)-U(i))+spring(i) &
                       *(X(p)-X(i))*(1.0-mag_init/mag)
                sum_y=sum_y+mu*(V(p)-V(i))+spring(i) &
                      *(Y(p)-Y(i))*(1.0-mag_init/mag)
                sum_z=sum_z+mu*(W(p)-W(i))+spring(i) &
                      *(Z(p)-Z(i))*(1.0-mag_init/mag)
              End If
           END DO
           R(i)=sum_x-nu*U(i)
           S(i)=sum_y-nu*V(i)
           T(i)=sum_z-nu*W(i)
          End If
        
        END DO
           
      DO i=1,Nodes
         U(i)=U(i)+DT_DUFFY*R(i)
         V(i)=V(i)+DT_DUFFY*S(i)
         W(i)=W(i)+DT_DUFFY*T(i)
      END DO


      DO i=1,Nodes
         X(i)=X(i)+DT_DUFFY*U(i)
         Y(i)=Y(i)+DT_DUFFY*V(i)
         Z(i)=Z(i)+DT_DUFFY*W(i)
      END DO


      count1=1*STEPS_DUFFY
      count2=2*STEPS_DUFFY
      count3=3*STEPS_DUFFY
      count4=4*STEPS_DUFFY

      counter=counter+1
      counter_real=counter_real+1.0

      if (counter.ge.count4) then
       counter=0
       counter_real=0.0
      endif

! W => horizontal direction
! V => vertical direction
! U => lateral direction

      If ((count1.eq.counter).or.(count3.eq.counter)) Then

      DO j=1,Nodes
         U(j)=0.0
         V(j)=-V(j)
         W(j)=0.0
         R(j)=0.0
         S(j)=0.0
         T(j)=0.0
      END DO

      End If

      If (count2.eq.counter) Then

         call taildown(R,S,T,U,V,W,X,Y,Z,Nodes,tailpos)

      END IF

      if (counter.eq.0) then
       call tailup(R,S,T,U,V,W,X,Y,Z,Nodes,tailpos)
      endif

      
      END subroutine new_geometry
    
      subroutine nodelist(Nodes, Ext, List,whaleout)

      IMPLICIT NONE

      CHARACTER*35 whaleout
      INTEGER Nodes, Ext
      INTEGER garbage, kk
      INTEGER N, Cells, Shape, i, j, k, n1, n2, n3
      REAL x, y, z

      INTEGER, DIMENSION(Nodes,Ext) :: List
      INTEGER, DIMENSION(Nodes) :: Counter

      open(unit=12, file=whaleout)

      read(12,*) garbage
      read(12,*) Cells
      read(12,*) Shape

      DO i=1,Nodes
         read(12,*) N, x, y, z
         List(i,1)=N
         Counter(i)=1
      END DO

      
      DO i=1,3 
         read(12,*) garbage
      END DO

      DO i=1,Cells-1
         read(12,*) n1
         read(12,*) n2
         read(12,*) n3
         
         j=Counter(n1)
         List(n1,j)=n2
         j=j+1
         List(n1,j)=n3
         Counter(n1)=j+1

         j=Counter(n2)
         List(n2,j)=n1
         j=j+1
         List(n2,j)=n3
         Counter(n2)=j+1

         j=Counter(n3)
         List(n3,j)=n1
         j=j+1
         List(n3,j)=n2
         Counter(n3)=j+1

         read(12,*) garbage
         read(12,*) garbage
      END DO


         read(12,*) n1
         read(12,*) n2
         read(12,*) n3

         j=Counter(n1)
         List(n1,j)=n2
         j=j+1
         List(n1,j)=n3
         Counter(n1)=j+1

         j=Counter(n2)
         List(n2,j)=n1
         j=j+1
         List(n2,j)=n3
         Counter(n2)=j+1

         j=Counter(n3)
         List(n3,j)=n1
         j=j+1
         List(n3,j)=n2
         Counter(n3)=j+1

      close(12)


      DO i=1,Nodes
         j=Counter(i)
         DO k=1,j-1
            n1=List(i,k)
            DO kk=k+1,j-1
               n2=List(i,kk)
               If (n1.eq.n2) Then
                  List(i,kk)=0
               End If
            END DO
         END DO
      END DO

      open(unit=10, file='testfile2.dat')

      DO i=1,Nodes
         j=Counter(i)
         write(10,*) i, (List(i,k), k=1,j-1)
      END DO

      close(10)

      END subroutine
         
      subroutine plininterp(xvect, yvect, x, n, value)

      IMPLICIT NONE

      INTEGER_T k, n
      REAL_T, Dimension(n) :: xvect, yvect
      REAL_T x, value

      
      DO k=1,n-1
      
         If (x.gt.xvect(k).and.x.lt.xvect(k+1)) Then
            value=yvect(k)+(yvect(k+1)-yvect(k))*(x-xvect(k)) &
                 /(xvect(k+1)-xvect(k))
         End If
         
         
         If (x.eq.xvect(k)) Then
            value=yvect(k)
         End if
       

      END DO


      END subroutine plininterp





      subroutine tailup(R,S,T,U,V,W,X,Y,Z, Nodes, tailpos)

      IMPLICIT NONE

      INTEGER_T Nodes, i
      REAL_T tailpos
      REAL_T, Dimension(Nodes) :: R,S,T,U,V,W,X,Y,Z

      
      DO i=1,Nodes
         U(i)=0.0
         V(i)=0.0
         V(i)=0.1+0.1*tanh(0.5*Z(i)-0.75)

              If (Z(i).lt.0.0) Then
                 V(i)=V(i)-.01*Z(i)
              End If

             If (Z(i).ge.0.0) Then
                V(i)=V(i)+.003*Z(i)**2
             End If
              
              If (Z(i).ge.tailpos) Then
                 V(i)=V(i)-.006*Z(i)+.006*tailpos
              End If   
      END DO

      DO i=1,Nodes
         W(i)=0.0
         R(i)=0.0
         S(i)=0.0
         T(i)=0.0
      END DO

      
      END subroutine tailup
      




      subroutine taildown(R,S,T,U,V,W,X,Y,Z,Nodes,tailpos)

      IMPLICIT NONE

      INTEGER_T Nodes, i
      REAL_T tailpos
      REAL_T, Dimension(Nodes) :: R,S,T,U,V,W,X,Y,Z

      
      DO i=1,Nodes
         U(i)=0.0
         V(i)=0.0
         V(i)=-(0.1+0.1*tanh(0.55*Z(i)-2.0))

               If (Z(i).lt.0.0) Then
                  V(i)=V(i)+0.015*Z(i)
               End If

               If (Z(i).gt.0.0) Then
                  V(i)=V(i)-.003*Z(i)**2
               End If

              If (Z(i).gt.tailpos) Then
                   V(i)=V(i)+.004*Z(i)-.004*tailpos
               End If

         W(i)=0.0
         R(i)=0.0
         S(i)=0.0
         T(i)=0.0

      END DO

      END subroutine taildown

      subroutine springs(Z, spring, Nodes,tailpos)
      IMPLICIT NONE

      INTEGER_T Nodes, i
      REAL_T tailpos
      REAL_T, Dimension(Nodes) :: spring, Z
      REAL_T springcons


      springcons=0.5

      DO i=1,Nodes
           spring(i)=0.1

          If (Z(i).gt.tailpos) Then
            spring(i)=spring(i)+springcons*(10.0-Z(i))**2 &
                     -springcons*(10.0-tailpos)**2
          End If

      END DO

      END subroutine springs

end module CLSVOFCouplerIO

#undef STANDALONE

