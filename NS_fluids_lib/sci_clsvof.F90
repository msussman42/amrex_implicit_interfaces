#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_FORT_INTEGER.H"
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "AMReX_ArrayLim_SUSSMAN.H"
#include "EXTRAP_COMP.H"

#define element_buffer_tol (0.0d0)
#define angle_tol (5.0d0)
#define max_plane_intersects 100
#define crossing_tol (0.001d0)

#define tecplot_post_process 1

! 10 seconds for tail to do a full period
#define WHALE_LENGTH 13.0
#define PERIOD_TAIL 10.0
#define DT_DUFFY 0.01d0
#define STEPS_DUFFY 800
#define injG 1
#define MAX_PARTS 100

#define flags_per_element 3
#define DOUBLYCOMP 3

#define CTMLoverflow (1.0D+20)
#define CTMLunderflow (1.0D-20)

#define BoundingBoxRadNode 2
#define BoundingBoxRadCell 3

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

INTEGER_T, PARAMETER :: sci_sdim=3

type lag_type
 INTEGER_T :: n_nodes,n_elems
 REAL_T, pointer :: nd(:,:)    ! nd(dir,node_id) dir=1..3
 REAL_T, pointer :: ndvel(:,:) ! ndvel(dir,node_id) dir=1..3
  ! ndmass=nddensity * ndvolume
 REAL_T, pointer :: ndmass(:) ! ndmass(node_id) 
 REAL_T, pointer :: nddensity(:) ! nddensity(node_id) 
 REAL_T, pointer :: ndforce(:,:) ! ndforce(dir,node_id)  dir=1..6
 REAL_T, pointer :: ndtemp(:)  ! ndtemp(node_id)
  ! number of nodes in element=elemdt(1,elemid)
  ! part number=elemdt(2,elemid)
  ! doubly wetted=elemdt(DOUBLYCOMP,elemid)
 INTEGER_T, pointer :: elemdt(:,:) 
  ! node_id=intelemdt(1..3,elemid)
  ! root (parent) element id = intelemdt(4,elemid)
 INTEGER_T, pointer :: intelemdt(:,:) 
end type lag_type

type mesh_type
 INTEGER_T :: part_id
 INTEGER_T :: flag_2D_to_3D
 INTEGER_T :: LS_FROM_SUBROUTINE
 INTEGER_T :: CTML_flag
 INTEGER_T :: refine_factor
 INTEGER_T :: bounding_box_ngrow
 REAL_T :: max_side_len
 REAL_T :: min_side_len
 REAL_T :: max_side_len_refined
 REAL_T :: min_side_len_refined
 INTEGER_T :: IntElemDim ! number of nodes (or edges) per element
 INTEGER_T :: IntElemDimPaddle
 INTEGER_T :: IntElemDimPool
 INTEGER_T :: NumNodes,NumNodesPaddle,NumNodesPool
 INTEGER_T :: NumIntElems,NumIntElemsPaddle,NumIntElemsPool
 INTEGER_T, pointer :: ElemData(:,:) !(nflags,nelements)
 INTEGER_T, pointer :: IntElem(:,:) !(nodes_per_elem,nelements)
 !average of adjoining element normals.
 REAL_T, pointer :: EdgeNormal(:,:) !(3*nodes_per_elem,nelements) 
 INTEGER_T, pointer :: EdgeElemId(:,:) !(nodes_per_elem,nelements)
 INTEGER_T, pointer :: EdgeElemIdNode(:,:) !(nodes_per_elem,nelements)
 INTEGER_T :: NumNodesBIG
 INTEGER_T :: NumIntElemsBIG
 INTEGER_T, pointer :: ElemNodeCountBIG(:)
 INTEGER_T, pointer :: ElemNodeCountEdgeBIG(:)
 REAL_T, pointer :: NodeBIG(:,:)  ! (3,node_id)
 REAL_T, pointer :: NodeVelBIG(:,:)
 REAL_T, pointer :: NodeForceBIG(:,:)  ! (NCOMP_FORCE_STRESS,node_id)
 REAL_T, pointer :: NodeDensityBIG(:)
 REAL_T, pointer :: NodeMassBIG(:)
 REAL_T, pointer :: NodeTempBIG(:)  
 REAL_T, pointer :: NodeNormalBIG(:,:)
 REAL_T, pointer :: NodeNormalEdgeBIG(:,:) ! sanity check purposes
 REAL_T, pointer :: ElemDataXnotBIG(:,:)
 INTEGER_T, pointer :: ElemDataBIG(:,:)
  ! root (parent) element id = IntElemBIG(4,elemid)
 INTEGER_T, pointer :: IntElemBIG(:,:) ! IntElemBIG(inode,ielem)
 !average of adjoining element normals.
 REAL_T, pointer :: EdgeNormalBIG(:,:) ! EdgeNormalBIG(3*(inode-1)+dir,ielem)
 INTEGER_T, pointer :: EdgeElemIdBIG(:,:)!EdgeElemIdBIG(inode,ielem)
 INTEGER_T, pointer :: EdgeElemIdNodeBIG(:,:)!EdgeElemIdNodeBIG(inode,ielem)
 REAL_T, pointer :: Node(:,:)  ! Node(dir,inode)
 REAL_T, pointer :: Node_old(:,:)
 REAL_T, pointer :: Node_new(:,:)
 REAL_T, pointer :: Node_current(:,:)
 REAL_T, pointer :: NodeVel(:,:)
 INTEGER_T, pointer :: ElemNodeCount(:)
 INTEGER_T, pointer :: ElemNodeCountEdge(:)
 REAL_T, pointer :: NodeNormal(:,:) !NodeNormal(dir,inode)
 REAL_T, pointer :: NodeNormalEdge(:,:) !sanity check purposes
 REAL_T, pointer :: NodeVel_old(:,:)
 REAL_T, pointer :: NodeVel_new(:,:)
 REAL_T, pointer :: NodeForce(:,:) ! NCOMP_FORCE_STRESS, NumNodes
 REAL_T, pointer :: NodeForce_old(:,:) ! NCOMP_FORCE_STRESS, NumNodes
 REAL_T, pointer :: NodeForce_new(:,:) ! NCOMP_FORCE_STRESS, NumNodes
 REAL_T, pointer :: NodeDensity(:)
 REAL_T, pointer :: NodeMass(:)
 REAL_T, pointer :: NodeTemp(:)
 REAL_T, pointer :: NodeTemp_old(:)
 REAL_T, pointer :: NodeTemp_new(:)
 REAL_T, pointer :: edge_endpoints(:,:)
 INTEGER_T, pointer :: edge_ielem(:)
 REAL_T soliddrop_displacement
 REAL_T soliddrop_speed
 REAL_T solid_displ(3)
 REAL_T solid_speed(3)
 REAL_T exterior_BB(3,2)
 REAL_T interior_BB(3,2)
 REAL_T center_BB(3)
 INTEGER_T deforming_part
 INTEGER_T normal_invert
 INTEGER_T exclusive_doubly_wetted
end type mesh_type


INTEGER_T :: use_temp
INTEGER_T :: istepB,sci_istop,sci_istep
REAL_T :: sci_curtime,sci_dt
REAL_T :: timeB,tstart,tfinish
REAL_T :: dtB

type(lag_type), dimension(:), allocatable :: multi_lag
type(mesh_type), dimension(MAX_PARTS) :: FSI

type(lag_type), dimension(:), allocatable :: aux_multi_lag
type(mesh_type), dimension(:), allocatable :: aux_FSI

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

REAL_T problo_ref(AMREX_SPACEDIM)
REAL_T probhi_ref(AMREX_SPACEDIM)
REAL_T problen_ref(AMREX_SPACEDIM)

REAL_T problo_act(3)
REAL_T probhi_act(3)
REAL_T problen_act(3)

INTEGER_T FSI_NPARTS
INTEGER_T CTML_NPARTS
INTEGER_T TOTAL_NPARTS
INTEGER_T im_solid_mapF(MAX_PARTS) ! type: 0..num_materials-1
INTEGER_T ctml_part_id_map(MAX_PARTS)
INTEGER_T fsi_part_id_map(MAX_PARTS)

INTEGER_T ctml_n_fib_bodies
INTEGER_T ctml_max_n_fib_nodes
INTEGER_T, dimension(:), allocatable :: ctml_n_fib_nodes
INTEGER_T, dimension(:), allocatable :: ctml_n_fib_active_nodes

REAL_T, dimension(:,:,:), allocatable :: ctml_fib_pst
REAL_T, dimension(:,:,:), allocatable :: ctml_fib_vel
REAL_T, dimension(:,:,:), allocatable :: ctml_fib_frc
REAL_T, dimension(:,:), allocatable :: ctml_fib_mass

REAL_T, dimension(:,:,:), allocatable :: ctml_fib_pst_prev
REAL_T, dimension(:,:,:), allocatable :: ctml_fib_vel_halftime_prev
REAL_T, dimension(:,:,:), allocatable :: ctml_fib_vel_prev
REAL_T, dimension(:,:,:), allocatable :: ctml_fib_frc_prev
REAL_T, dimension(:,:), allocatable :: ctml_fib_mass_prev

contains


! called from:
!   initinjector,initflapping,init_from_cas,init_gingerbread2D,
!   init_helix,initchannel,viorel_sphere_geominit,internal_inflow_geominit,
!   gearinit,initpaddle,initship
subroutine init2_FSI(part_id)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%part_id.ne.part_id) then
  print *,"FSI(part_id)%part_id.ne.part_id"
  stop
 endif

 call init2_FSI_mesh_type(FSI(part_id))

return
end subroutine init2_FSI


subroutine init3_FSI(part_id,ifirst,do_2nd_part,ioproc,isout)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: ifirst
INTEGER_T, INTENT(in) :: do_2nd_part
INTEGER_T, INTENT(in) :: ioproc
INTEGER_T, INTENT(in) :: isout
INTEGER_T inode,dir,it
REAL_T x,y,z,z0,z90,t,dt,t1,t2,inflowvel
REAL_T YK,ZK,lift0,lift90
REAL_T, dimension(3) :: displ1,displ2
REAL_T dilated_time

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id invalid"
   stop
  endif
  if (FSI(part_id)%part_id.ne.part_id) then
   print *,"FSI(part_id)%part_id.ne.part_id"
   stop
  endif

  call init3_FSI_mesh_type(FSI(part_id),ifirst)

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
! dist(x,y,z)= dist(x,y,z+0.01d0)
! IN FUTURE, DO NOT CALL "UNITE", INSTEAD
! USE BOTH FSI and FSI_NEEDLE when returning LS or VEL.
! ALSO IN FUTURE, DO NOT REPEATEDLY CALL GENERATE_NEW_TRIANGLES
! IN CLSVOF_ReadNodes.
!  
  else if ((probtype.eq.538).or.(probtype.eq.541)) then

   if (FSI(part_id)%part_id.eq.2) then

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
      displ1(3) = 0.5d0*(z0+z90)
     else if (t.gt.dilated_time) then
      t2 = t
      displ2(1) = x
      displ2(2) = y
      displ2(3) = 0.5d0*(z0+z90)
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

    if (levelrz.eq.COORDSYS_CARTESIAN) then
     ! do nothing
    else if (levelrz.eq.COORDSYS_RZ) then
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
   else if (FSI(part_id)%part_id.eq.1) then

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
! effect what happens in post_process_nodes_elements.

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

subroutine init3_FSI_mesh_type(FSI_mesh_type,ifirst)
type(mesh_type), INTENT(inout) :: FSI_mesh_type
INTEGER_T, INTENT(in) :: ifirst
INTEGER_T :: inode
INTEGER_T :: dir

 if (ifirst.eq.1) then
  allocate(FSI_mesh_type%Node(3,FSI_mesh_type%NumNodes))
  allocate(FSI_mesh_type%NodeVel(3,FSI_mesh_type%NumNodes))
  allocate(FSI_mesh_type%NodeForce(NCOMP_FORCE_STRESS,FSI_mesh_type%NumNodes))
  allocate(FSI_mesh_type%NodeNormal(3,FSI_mesh_type%NumNodes))
  allocate(FSI_mesh_type%NodeNormalEdge(3,FSI_mesh_type%NumNodes))
  allocate(FSI_mesh_type%ElemNodeCount(FSI_mesh_type%NumNodes))
  allocate(FSI_mesh_type%ElemNodeCountEdge(FSI_mesh_type%NumNodes))
  allocate(FSI_mesh_type%NodeTemp(FSI_mesh_type%NumNodes))
 else if (ifirst.eq.0) then
  ! do nothing
 else
  print *,"ifirst invalid"
  stop
 endif

 do inode=1,FSI_mesh_type%NumNodes
  do dir=1,3
   FSI_mesh_type%Node(dir,inode)=FSI_mesh_type%Node_current(dir,inode)
   FSI_mesh_type%NodeVel(dir,inode)=FSI_mesh_type%NodeVel_new(dir,inode)
  enddo
  do dir=1,NCOMP_FORCE_STRESS
   FSI_mesh_type%NodeForce(dir,inode)=FSI_mesh_type%NodeForce_new(dir,inode)
  enddo
  FSI_mesh_type%NodeTemp(inode)=FSI_mesh_type%NodeTemp_new(inode)
 enddo ! inode=1,FSI_mesh_type%NumNodes

 do dir=1,3
  FSI_mesh_type%solid_displ(dir)=zero
  FSI_mesh_type%solid_speed(dir)=zero
 enddo

 FSI_mesh_type%soliddrop_displacement=zero
 FSI_mesh_type%soliddrop_speed=zero

return
end subroutine init3_FSI_mesh_type


subroutine init2_FSI_mesh_type(FSI_mesh_type)
type(mesh_type), INTENT(inout) :: FSI_mesh_type
INTEGER_T :: inode
INTEGER_T :: dir

 do inode=1,FSI_mesh_type%NumNodes
  do dir=1,3
   FSI_mesh_type%NodeVel_old(dir,inode)=0.0
  enddo
  do dir=1,3
   FSI_mesh_type%NodeForce_old(dir,inode)=0.0
  enddo
  do dir=1,3
   FSI_mesh_type%Node_current(dir,inode)=FSI_mesh_type%Node_new(dir,inode)
   FSI_mesh_type%NodeVel_new(dir,inode)=FSI_mesh_type%NodeVel_old(dir,inode)
  enddo
  do dir=1,3
   FSI_mesh_type%NodeForce_new(dir,inode)=FSI_mesh_type%NodeForce_old(dir,inode)
  enddo
  FSI_mesh_type%NodeMass(inode)=one
  FSI_mesh_type%NodeDensity(inode)=one
 enddo  ! inode=1,NumNodes

return
end subroutine init2_FSI_mesh_type

subroutine init_FSI_mesh_type(FSI_mesh_type,allocate_intelem)
type(mesh_type), INTENT(inout) :: FSI_mesh_type
INTEGER_T, INTENT(in) :: allocate_intelem
INTEGER_T :: inode
INTEGER_T :: dir

 !(1,iface)=nodes per element (2,iface)=part num 
 !(DOUBLYCOMP,iface)=doubly wet flag
 allocate(FSI_mesh_type%ElemData(flags_per_element,FSI_mesh_type%NumIntElems))
 if (allocate_intelem.eq.1) then
  allocate(FSI_mesh_type%IntElem(FSI_mesh_type%IntElemDim, &
           FSI_mesh_type%NumIntElems))
  !average of adjoining element normals.
  allocate(FSI_mesh_type%EdgeNormal(3*FSI_mesh_type%IntElemDim, &
           FSI_mesh_type%NumIntElems))
  allocate(FSI_mesh_type%EdgeElemId(FSI_mesh_type%IntElemDim, &
           FSI_mesh_type%NumIntElems))
  allocate(FSI_mesh_type%EdgeElemIdNode(FSI_mesh_type%IntElemDim, &
           FSI_mesh_type%NumIntElems))
 else if (allocate_intelem.eq.0) then
  ! do nothing
 else
  print *,"allocate_intelem invalid"
  stop
 endif
 allocate(FSI_mesh_type%Node_old(3,FSI_mesh_type%NumNodes))
 allocate(FSI_mesh_type%Node_new(3,FSI_mesh_type%NumNodes))
 allocate(FSI_mesh_type%Node_current(3,FSI_mesh_type%NumNodes))
 allocate(FSI_mesh_type%NodeVel_old(3,FSI_mesh_type%NumNodes))
 allocate(FSI_mesh_type%NodeVel_new(3,FSI_mesh_type%NumNodes))
 allocate(FSI_mesh_type%NodeForce_old(3,FSI_mesh_type%NumNodes))
 allocate(FSI_mesh_type%NodeForce_new(3,FSI_mesh_type%NumNodes))
 allocate(FSI_mesh_type%NodeTemp_old(FSI_mesh_type%NumNodes))
 allocate(FSI_mesh_type%NodeTemp_new(FSI_mesh_type%NumNodes))

 allocate(FSI_mesh_type%NodeMass(FSI_mesh_type%NumNodes))
 allocate(FSI_mesh_type%NodeDensity(FSI_mesh_type%NumNodes))

 do inode=1,FSI_mesh_type%NumNodes

  FSI_mesh_type%NodeMass(inode)=one
  FSI_mesh_type%NodeDensity(inode)=one

  FSI_mesh_type%NodeTemp_old(inode)=0.0
  FSI_mesh_type%NodeTemp_new(inode)=0.0
  do dir=1,3
   FSI_mesh_type%Node_old(dir,inode)=0.0
   FSI_mesh_type%Node_new(dir,inode)=0.0
   FSI_mesh_type%Node_current(dir,inode)=0.0
   FSI_mesh_type%NodeVel_old(dir,inode)=0.0
   FSI_mesh_type%NodeVel_new(dir,inode)=0.0
  enddo
  do dir=1,3
   FSI_mesh_type%NodeForce_old(dir,inode)=0.0
   FSI_mesh_type%NodeForce_new(dir,inode)=0.0
  enddo
 enddo ! inode=1,FSI_mesh_type%NumNodes

return
end subroutine init_FSI_mesh_type


! called from:
!  CTML_init_sci (prior to CTML_init_sci_node),
!  initinjector,initflapping,init_from_cas,
!  init_gingerbread2D,init_helix,initchannel,geominit,
!  viorel_sphere_geominit,internal_inflow_geominit,
!  gearinit,whale_geominit,initpaddle,initship
subroutine init_FSI(part_id,allocate_intelem)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: allocate_intelem

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%part_id.ne.part_id) then
  print *,"FSI(part_id)%part_id.ne.part_id"
  stop
 endif

 call init_FSI_mesh_type(FSI(part_id),allocate_intelem)

return
end subroutine init_FSI

subroutine xdist_project(x1,x2, &
      FSI_mesh_type,part_id,max_part_id, &
      dist_project,dist_actual)
IMPLICIT NONE

type(mesh_type), INTENT(in) :: FSI_mesh_type
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
REAL_T, dimension(3),INTENT(in) :: x1,x2
REAL_T, INTENT(out) :: dist_project
REAL_T, INTENT(out) :: dist_actual
INTEGER_T sdim_local
INTEGER_T dir

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI_mesh_type%flag_2D_to_3D.eq.1) then
  sdim_local=2
 else if (FSI_mesh_type%flag_2D_to_3D.eq.0) then
  sdim_local=3
 else
  print *,"FSI_mesh_type%flag_2D_to_3D invalid"
  stop
 endif
 dist_project=zero
 do dir=1,sdim_local
  dist_project=dist_project+(x1(dir)-x2(dir))**2
 enddo
 if (dist_project.ge.zero) then
  dist_project=sqrt(dist_project)
 else
  print *,"dist_project bust"
  stop
 endif

 dist_actual=zero
 do dir=1,3
  dist_actual=dist_actual+(x1(dir)-x2(dir))**2
 enddo
 if (dist_actual.ge.zero) then
  dist_actual=sqrt(dist_actual)
 else
  print *,"dist_actual bust"
  stop
 endif


return
end subroutine xdist_project


subroutine xdistmin(x1,x2,dist)
IMPLICIT NONE

REAL_T, dimension(3),INTENT(in) :: x1,x2
REAL_T, INTENT(out) :: dist
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

subroutine get_new_half_vols(x1,x2,xsplit,volL,volR)
IMPLICIT NONE

REAL_T, INTENT(in), dimension(3) :: x1,x2,xsplit
REAL_T, INTENT(out) :: volL,volR
INTEGER_T :: dir

 volL=zero
 volR=zero
 do dir=1,3
  volL=volL+(xsplit(dir)-x1(dir))**2
  volR=volR+(xsplit(dir)-x2(dir))**2
 enddo
 if ((volL.gt.zero).and.(volR.gt.zero)) then
  volL=sqrt(volL)
  volR=sqrt(volR)
 else
  print *,"volL or volR invalid in get_new_half_vols"
  print *,"volL,volR= ",volL,volR
  stop
 endif

end subroutine get_new_half_vols


subroutine compare_core(nodej,nodejp1,coord_scale,compare_flag,ncore)
IMPLICIT NONE
INTEGER_T, INTENT(in) :: ncore
INTEGER_T, INTENT(out) :: compare_flag
REAL_T, INTENT(in) :: coord_scale
REAL_T, INTENT(in) :: nodej(ncore)
REAL_T, INTENT(in) :: nodejp1(ncore)
REAL_T :: mag
INTEGER_T :: dir

 if (coord_scale.gt.zero) then
  ! do nothing
 else
  print *,"coord_scale invalid"
  stop
 endif

 if ((ncore.eq.3).or.(ncore.eq.6)) then
  ! do nothing
 else
  print *,"ncore invalid"
  stop
 endif
 mag=zero
 do dir=1,ncore
  mag=mag+(nodej(dir)-nodejp1(dir))**2
 enddo
 mag=sqrt(mag)
 if (mag.le.ncore*VOFTOL*coord_scale) then
  compare_flag=0
 else if (mag.ge.ncore*VOFTOL*coord_scale) then
  compare_flag=0
  dir=1
  do while ((compare_flag.eq.0).and.(dir.le.ncore))

   if (nodej(dir).gt.nodejp1(dir)+VOFTOL*coord_scale) then
    compare_flag=1
   else if (nodej(dir).lt.nodejp1(dir)-VOFTOL*coord_scale) then
    compare_flag=-1
   else if (abs(nodej(dir)-nodejp1(dir)).le.VOFTOL*coord_scale) then
    ! do nothing
   else
    print *,"nodej or nodejp1 are NaN"
    stop
   endif
   dir=dir+1

  enddo !while ((compare_flag.eq.0).and.(dir.le.ncore))
 else
  print *,"mag invalid"
  stop
 endif

end subroutine compare_core


subroutine compare_edge(edgej,edgejp1,coord_scale,compare_flag,overlap_size)
use global_utility_module
IMPLICIT NONE
INTEGER_T, INTENT(out) :: compare_flag
REAL_T, INTENT(out) :: overlap_size
REAL_T, INTENT(in) :: coord_scale
REAL_T, INTENT(in) :: edgej(6)
REAL_T, INTENT(in) :: edgejp1(6)
REAL_T :: overlap_start,overlap_end
REAL_T :: mag
REAL_T :: map_mag
REAL_T :: mag_offline
REAL_T :: edgej_mag
REAL_T :: edgejp1_mag
REAL_T :: edgej_center(3)
REAL_T :: edgejp1_center(3)
INTEGER_T :: mincomp
INTEGER_T :: maxcomp
INTEGER_T :: dir
INTEGER_T :: imat,jmat
REAL_T :: mapx1(3)
REAL_T :: vec1(3)
REAL_T :: vec2(3)
REAL_T :: vec3(3)
REAL_T :: x1test(3)
REAL_T :: x2test(3)
REAL_T :: x1test_map(3)
REAL_T :: x2test_map(3)
REAL_T :: A(3,3)
REAL_T :: AINV(3,3)

 if (coord_scale.gt.zero) then
  ! do nothing
 else
  print *,"coord_scale invalid"
  stop
 endif
 compare_flag=1 ! default value

 edgej_mag=zero
 edgejp1_mag=zero
 do dir=1,3
  edgej_mag=edgej_mag+(edgej(dir)-edgej(dir+3))**2
  edgej_center(dir)=half*(edgej(dir)+edgej(dir+3))
  edgejp1_mag=edgejp1_mag+(edgejp1(dir)-edgejp1(dir+3))**2
  edgejp1_center(dir)=half*(edgejp1(dir)+edgejp1(dir+3))
 enddo
 edgej_mag=sqrt(edgej_mag)
 edgejp1_mag=sqrt(edgejp1_mag)
 if ((edgej_mag.gt.zero).and.(edgejp1_mag.gt.zero)) then

  call compare_core(edgej_center,edgejp1_center,coord_scale,compare_flag,3)
  if (compare_flag.eq.0) then
   overlap_size=one
  else if ((compare_flag.eq.1).or.(compare_flag.eq.-1)) then

   ! y=(x-x1)/||x2-x1|| gets mapped in such a way that 
   ! A y2 = (1 0 0)
   ! AINV (1 0 0)=y2
   ! first column of AINV is y2
   if (edgej_mag.ge.edgejp1_mag) then
    map_mag=edgej_mag
    do dir=1,3
     mapx1(dir)=edgej(dir)
     x1test(dir)=(edgejp1(dir)-mapx1(dir))/map_mag
     x2test(dir)=(edgejp1(dir+3)-mapx1(dir))/map_mag
     AINV(dir,1)=(edgej(dir+3)-mapx1(dir))/map_mag
    enddo
   else if (edgej_mag.le.edgejp1_mag) then
    map_mag=edgejp1_mag
    do dir=1,3
     mapx1(dir)=edgejp1(dir)
     x1test(dir)=(edgej(dir)-mapx1(dir))/map_mag
     x2test(dir)=(edgej(dir+3)-mapx1(dir))/map_mag
     AINV(dir,1)=(edgejp1(dir+3)-mapx1(dir))/map_mag
    enddo
   else
    print *,"edgej_mag or edgejp1_mag invalid"
    stop
   endif
   maxcomp=1
   mincomp=1
   do dir=2,3
    if (abs(AINV(dir,1)).gt.abs(AINV(maxcomp,1))) then
     maxcomp=dir
    endif
    if (abs(AINV(dir,1)).lt.abs(AINV(mincomp,1))) then
     mincomp=dir
    endif
   enddo
   if (abs(AINV(maxcomp,1)).gt.zero) then
    do dir=1,3
     AINV(dir,2)=zero
    enddo
    mag=sqrt(AINV(mincomp,1)**2+AINV(maxcomp,1)**2)
    if (mag.gt.zero) then
     AINV(mincomp,2)=AINV(maxcomp,1)/mag
     AINV(maxcomp,2)=-AINV(mincomp,1)/mag
     do dir=1,3
      vec1(dir)=AINV(dir,1)
      vec2(dir)=AINV(dir,2)
     enddo
     call crossprod(vec1,vec2,vec3)
     do dir=1,3
      AINV(dir,3)=vec3(dir)
     enddo
     do imat=1,3
     do jmat=1,3
      A(imat,jmat)=AINV(jmat,imat)
     enddo
     enddo
     mag_offline=zero
     do imat=1,3
      x1test_map(imat)=zero 
      x2test_map(imat)=zero 
      do jmat=1,3
       x1test_map(imat)=x1test_map(imat)+A(imat,jmat)*x1test(jmat)
       x2test_map(imat)=x2test_map(imat)+A(imat,jmat)*x2test(jmat)
      enddo
      if (imat.eq.1) then
       ! do nothing
      else if ((imat.eq.2).or.(imat.eq.3)) then
       mag_offline=mag_offline+x1test_map(imat)**2+x2test_map(imat)**2
      else
       print *,"imat invalid"
       stop
      endif
     enddo !imat=1..3
     mag_offline=sqrt(mag_offline)
     overlap_start=max(zero,min(x1test_map(1),x2test_map(1)))
     overlap_end=min(one,max(x1test_map(1),x2test_map(1)))
     overlap_size=max(zero,overlap_end-overlap_start)

     if ((mag_offline.le.VOFTOL).and. &
         (mag_offline.ge.zero).and. &
         (overlap_size.ge.VOFTOL).and. &
         (overlap_size.le.one+VOFTOL)) then
      if (overlap_size.ge.one-VOFTOL) then
       compare_flag=0
      else if ((overlap_size.le.one-VOFTOL).and. &
               (overlap_size.ge.VOFTOL)) then

       if (compare_flag.eq.1) then
        ! do nothing
       else if (compare_flag.eq.-1) then
        ! do nothing
       else
        print *,"compare_flag invalid"
        stop
       endif
  
       if (1.eq.0) then
        print *,"WARNING:overlap_size not equal to 1"
        print *,"check for hanging nodes"
        print *,"overlap_size=",overlap_size
        print *,"mag_offline=",mag_offline

        print *,"edgej_mag ",edgej_mag
        print *,"edgejp1_mag ",edgejp1_mag
        print *,"edgej ",edgej(1),edgej(2),edgej(3), &
         edgej(4),edgej(5),edgej(6)
        print *,"edgejp1 ",edgejp1(1),edgejp1(2),edgejp1(3), &
         edgejp1(4),edgejp1(5),edgejp1(6)
        print *,"x1test: ",x1test(1),x1test(2),x1test(3)
        print *,"x2test: ",x2test(1),x2test(2),x2test(3)
        print *,"x1test_map: ",x1test_map(1),x1test_map(2),x1test_map(3)
        print *,"x2test_map: ",x2test_map(1),x2test_map(2),x2test_map(3)
       endif

      else
       print *,"overlap_size invalid"
       stop
      endif
     else if ((mag_offline.ge.VOFTOL).or. &
              (overlap_size.le.VOFTOL)) then

      if (compare_flag.eq.1) then
       ! do nothing
      else if (compare_flag.eq.-1) then
       ! do nothing
      else
       print *,"compare_flag invalid"
       stop
      endif

      if (mag_offline.ge.VOFTOL) then
       ! do nothing
      else if ((mag_offline.le.VOFTOL).and.(mag_offline.ge.zero)) then 
       if ((overlap_size.ge.zero).and. &
           (overlap_size.le.VOFTOL)) then
        ! do nothing
       else
        print *,"overlap_size invalid1:",overlap_size
        print *,"overlap_size ",overlap_size
        print *,"mag_offline ",mag_offline
        print *,"edgej_mag ",edgej_mag
        print *,"edgejp1_mag ",edgejp1_mag
        print *,"edgej ",edgej(1),edgej(2),edgej(3), &
         edgej(4),edgej(5),edgej(6)
        print *,"edgejp1 ",edgejp1(1),edgejp1(2),edgejp1(3), &
         edgejp1(4),edgejp1(5),edgejp1(6)
        print *,"x1test: ",x1test(1),x1test(2),x1test(3)
        print *,"x2test: ",x2test(1),x2test(2),x2test(3)
        print *,"x1test_map: ",x1test_map(1),x1test_map(2),x1test_map(3)
        print *,"x2test_map: ",x2test_map(1),x2test_map(2),x2test_map(3)
        stop
       endif
      else
       print *,"mag_offline invalid"
       stop
      endif
     else
      print *,"mag_offline or overlap_size invalid:"
      print *,"overlap_size ",overlap_size
      print *,"mag_offline ",mag_offline
      stop
     endif
    else
     print *,"mag invalid"
     stop
    endif
   else
    print *,"abs(AINV(maxcomp,1)) invalid"
    stop
   endif
  else
   print *,"compare_flag invalid"
   stop
  endif

 else
  print *,"edgej_mag or edgejp1_mag invalid"
  stop
 endif

end subroutine compare_edge


subroutine compare_nodes(FSI_mesh_type, &
                jnode,jnodep1, &
                sorted_node_list, &
                coord_scale,compare_flag)
IMPLICIT NONE
type(mesh_type), INTENT(in) :: FSI_mesh_type
INTEGER_T, allocatable, INTENT(in) :: sorted_node_list(:)
REAL_T, INTENT(in) :: coord_scale
INTEGER_T, INTENT(in) :: jnode
INTEGER_T, INTENT(in) :: jnodep1
INTEGER_T, INTENT(out) :: compare_flag
REAL_T :: nodej(3)
REAL_T :: nodejp1(3)
INTEGER_T :: dir

 if (coord_scale.gt.zero) then
  ! do nothing
 else
  print *,"coord_scale invalid"
  stop
 endif

 do dir=1,3
  nodej(dir)=FSI_mesh_type%NodeBIG(dir,sorted_node_list(jnode))
  nodejp1(dir)=FSI_mesh_type%NodeBIG(dir,sorted_node_list(jnodep1))
 enddo
 call compare_core(nodej,nodejp1,coord_scale,compare_flag,3)

end subroutine compare_nodes

subroutine tecplot_normals(FSI_mesh_type,part_id,max_part_id,view_refined)
IMPLICIT NONE
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
INTEGER_T, INTENT(in) :: view_refined
type(mesh_type), INTENT(inout) :: FSI_mesh_type

character*2 partstr
character*20 auxfilename20
INTEGER_T :: stat
INTEGER_T :: nodes,cells
INTEGER_T :: dir
INTEGER_T :: i
INTEGER_T :: node_select
REAL_T :: xnode(3)
REAL_T :: nnode(3)
INTEGER_T :: inode(3)

 if (tecplot_post_process.eq.1) then
  write(partstr,'(I2)') part_id
  do i=1,2
   if (partstr(i:i).eq.' ') then
    partstr(i:i)='0'
   endif
  enddo

  do node_select=1,2

   if (view_refined.eq.1) then
    if (node_select.eq.1) then
     write(auxfilename20,'(A14,A2,A4)') 'auxfileBIGNODE',partstr,'.plt'
    else if (node_select.eq.2) then
     write(auxfilename20,'(A14,A2,A4)') 'auxfileBIGEDGE',partstr,'.plt'
    else
     print *,"node_select invalid"
     stop
    endif
   else if (view_refined.eq.0) then
    if (node_select.eq.1) then
     write(auxfilename20,'(A14,A2,A4)') 'auxfileREGNODE',partstr,'.plt'
    else if (node_select.eq.2) then
     write(auxfilename20,'(A14,A2,A4)') 'auxfileREGEDGE',partstr,'.plt'
    else
     print *,"node_select invalid"
     stop
    endif
   else
    print *,"view_refined invalid"
    stop
   endif

   print *,"auxfilename20 ",auxfilename20
   open(unit=4, file= auxfilename20,status='unknown',iostat=stat)
   if (stat.ne.0) then
    print *,auxfilename20," can not be opened"
    stop
   endif

   if (view_refined.eq.1) then
    nodes=FSI_mesh_type%NumNodesBIG
    cells=FSI_mesh_type%NumIntElemsBIG
   else if (view_refined.eq.0) then
    nodes=FSI_mesh_type%NumNodes
    cells=FSI_mesh_type%NumIntElems
   else
    print *,"view_refined invalid"
    stop
   endif

   write(4,*) 'TITLE = "3D surface" '
   write(4,*) 'VARIABLES = "X", "Y", "Z", "NX", "NY", "NZ" '
   write(4,*) 'ZONE T="TRIANGLES", N= ', Nodes, ', E= ', &
         Cells, ', DATAPACKING=POINT, '
   write(4,*) 'ZONETYPE=FETRIANGLE' 

   do i=1,Nodes
    if (view_refined.eq.1) then
     do dir=1,3
      xnode(dir)=FSI_mesh_type%NodeBIG(dir,i)
      if (node_select.eq.1) then
       nnode(dir)=FSI_mesh_type%NodeNormalBIG(dir,i)
      else if (node_select.eq.2) then
       nnode(dir)=FSI_mesh_type%NodeNormalEdgeBIG(dir,i)
      else
       print *,"node_select invalid"
       stop
      endif
     enddo
    else if (view_refined.eq.0) then
     do dir=1,3
      xnode(dir)=FSI_mesh_type%Node(dir,i)
      if (node_select.eq.1) then
       nnode(dir)=FSI_mesh_type%NodeNormal(dir,i)
      else if (node_select.eq.2) then
       nnode(dir)=FSI_mesh_type%NodeNormalEdge(dir,i)
      else
       print *,"node_select invalid"
       stop
      endif
     enddo ! dir=1..3
    else
     print *,"view_refined invalid"
     stop
    endif
    write(4,*) xnode(1),xnode(2),xnode(3),nnode(1),nnode(2),nnode(3)
   enddo ! i=1,Nodes

   do i=1,Cells
    if (view_refined.eq.1) then
     ! root (parent) element id = IntElemBIG(4,elemid)
     do dir=1,3
      inode(dir)=FSI_mesh_type%IntElemBIG(dir,i)
     enddo
    else if (view_refined.eq.0) then
     do dir=1,3
      inode(dir)=FSI_mesh_type%IntElem(dir,i)
     enddo
    else
     print *,"view_refined invalid"
     stop
    endif
    write(4,*) inode(1),inode(2),inode(3)
   enddo  

   close(4)
  enddo  !do node_select=1,2

 else if (tecplot_post_process.eq.0) then
  ! do nothing
 else
  print *,"tecplot_post_process invalid"
  stop
 endif

end subroutine tecplot_normals

subroutine TopDownMergeSort(FSI_mesh_type,coord_scale,A,B,n, &
      sort_nodes_flag)
IMPLICIT NONE
type(mesh_type), INTENT(in) :: FSI_mesh_type
REAL_T, INTENT(in) :: coord_scale
INTEGER_T, INTENT(in) :: sort_nodes_flag
INTEGER_T, INTENT(in) :: n
INTEGER_T, allocatable, INTENT(inout) :: A(:)
INTEGER_T, allocatable, INTENT(inout) :: B(:)

 print *,"in TopDownMergeSort: n,sort_nodes_flag=",n,sort_nodes_flag

 call CopyArray(FSI_mesh_type,coord_scale,A,0,n,B)
 call TopDownSplitMerge(FSI_mesh_type,coord_scale,B,0,n,A,sort_nodes_flag)

end subroutine TopDownMergeSort


recursive subroutine TopDownSplitMerge(FSI_mesh_type,coord_scale, &
 B,iBegin,iEnd,A,sort_nodes_flag)
IMPLICIT NONE
type(mesh_type), INTENT(in) :: FSI_mesh_type
REAL_T, INTENT(in) :: coord_scale
INTEGER_T, INTENT(in) :: sort_nodes_flag
INTEGER_T, INTENT(in) :: iBegin
INTEGER_T, INTENT(in) :: iEnd
INTEGER_T, allocatable, INTENT(inout) :: A(:)
INTEGER_T, allocatable, INTENT(inout) :: B(:)
INTEGER_T :: iMiddle

 if (iEnd-iBegin.le.1) then
  ! do nothing
 else
  iMiddle=(iEnd+iBegin)/2
  call TopDownSplitMerge(FSI_mesh_type,coord_scale,A,iBegin,iMiddle,B, &
    sort_nodes_flag)
  call TopDownSplitMerge(FSI_mesh_type,coord_scale,A,iMiddle,iEnd,B, &
    sort_nodes_flag)
  call TopDownMerge(FSI_mesh_type,coord_scale,B,iBegin,iMiddle,iEnd,A, &
    sort_nodes_flag)
 endif

end subroutine TopDownSplitMerge


subroutine TopDownMerge(FSI_mesh_type,coord_scale, &
 A,iBegin,iMiddle,iEnd,B,sort_nodes_flag)
IMPLICIT NONE
type(mesh_type), INTENT(in) :: FSI_mesh_type
REAL_T, INTENT(in) :: coord_scale
INTEGER_T, INTENT(in) :: sort_nodes_flag
INTEGER_T, INTENT(in) :: iBegin
INTEGER_T, INTENT(in) :: iMiddle
INTEGER_T, INTENT(in) :: iEnd
INTEGER_T, allocatable, INTENT(inout) :: A(:)
INTEGER_T, allocatable, INTENT(inout) :: B(:)
INTEGER_T :: i,j,k
INTEGER_T :: compare_flag
INTEGER_T :: dir
REAL_T :: edgej(6)
REAL_T :: edgejp1(6)
REAL_T :: overlap_size
INTEGER_T :: ielem,ielem_opp

 i=iBegin
 j=iMiddle
 k=iBegin

 do while (k.lt.iEnd)
   !Ai<Aj compare_flag=-1
   !Ai>Aj compare_flag=1

  compare_flag=0

  if ((i.lt.iMiddle).and.(j.lt.iEnd)) then

   if (sort_nodes_flag.eq.1) then
    call compare_nodes(FSI_mesh_type, &
     i+1,j+1, &
     A, &
     coord_scale,compare_flag)
   else if (sort_nodes_flag.eq.0) then
    do dir=1,6
     edgej(dir)=FSI_mesh_type%edge_endpoints(dir,A(i+1))
     edgejp1(dir)=FSI_mesh_type%edge_endpoints(dir,A(j+1))
    enddo
    call compare_edge(edgej,edgejp1,coord_scale,compare_flag,overlap_size)
    if ((compare_flag.eq.1).or. &
        (compare_flag.eq.-1)) then
     ! do nothing
    else if (compare_flag.eq.0) then
     ielem=FSI_mesh_type%edge_ielem(A(i+1))
     ielem_opp=FSI_mesh_type%edge_ielem(A(j+1))
     if (ielem.eq.ielem_opp) then
      if (i.ne.j) then
       print *,"cannot have two edges from the same element be equal"
       print *,"i,j = ",i,j
       stop
      endif 
     endif 
    else
     print *,"compare_flag invalid"
     stop
    endif
   else
    print *,"sort_nodes_flag invalid"
    stop
   endif
  else if ((i.ge.iMiddle).or.(j.ge.iEnd)) then
   ! do nothing
  else
   print *,"i,j bust"
   stop
  endif

  if ((i.lt.iMiddle).and. &
      ((j.ge.iEnd).or.(compare_flag.le.0))) then
   B(k+1)=A(i+1)
   i=i+1
  else
   B(k+1)=A(j+1)
   j=j+1
  endif
  k=k+1 
 enddo ! do while (k.lt.iEnd)

end subroutine TopDownMerge


subroutine CopyArray(FSI_mesh_type,coord_scale, &
 A,iBegin,iEnd,B)
IMPLICIT NONE
type(mesh_type), INTENT(in) :: FSI_mesh_type
REAL_T, INTENT(in) :: coord_scale
INTEGER_T, INTENT(in) :: iBegin
INTEGER_T, INTENT(in) :: iEnd
INTEGER_T, allocatable, INTENT(inout) :: A(:)
INTEGER_T, allocatable, INTENT(inout) :: B(:)
INTEGER_T :: k

 do k=iBegin,iEnd-1
  B(k+1)=A(k+1)
 enddo

end subroutine CopyArray

subroutine remove_duplicate_nodes(FSI_mesh_type,part_id,max_part_id)
IMPLICIT NONE
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
type(mesh_type), INTENT(inout) :: FSI_mesh_type
INTEGER_T, allocatable :: sorted_node_list(:)
INTEGER_T, allocatable :: B_list(:)
INTEGER_T, allocatable :: alternate_node_list(:)
INTEGER_T, allocatable :: new_node_list(:)
INTEGER_T, allocatable :: old_node_list(:)
REAL_T, allocatable :: NodeBIG_local(:,:)
REAL_T, allocatable :: NodeVelBIG_local(:,:)
REAL_T, allocatable :: NodeForceBIG_local(:,:)
REAL_T, allocatable :: NodeDensityBIG_local(:)
REAL_T, allocatable :: NodeMassBIG_local(:)
REAL_T, allocatable :: NodeTempBIG_local(:)
INTEGER_T :: inode,jnode,knode
INTEGER_T :: ilocal
INTEGER_T :: ielem
INTEGER_T :: nodes_per_elem
INTEGER_T :: compare_flag
INTEGER_T :: sort_nodes_flag
REAL_T :: min_coord
REAL_T :: max_coord
REAL_T :: coord_scale
REAL_T :: test_coord
INTEGER_T :: dir
INTEGER_T :: num_nodes_local 

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI_mesh_type%part_id.ne.part_id) then
  print *,"FSI_mesh_type%part_id.ne.part_id"
  stop
 endif

 print *,"removing duplicate nodes"

 allocate(sorted_node_list(FSI_mesh_type%NumNodesBIG))
 allocate(B_list(FSI_mesh_type%NumNodesBIG))
 allocate(alternate_node_list(FSI_mesh_type%NumNodesBIG))
 allocate(new_node_list(FSI_mesh_type%NumNodesBIG))
 allocate(old_node_list(FSI_mesh_type%NumNodesBIG))

 do inode=1,FSI_mesh_type%NumNodesBIG
  sorted_node_list(inode)=inode
  B_list(inode)=inode
  alternate_node_list(inode)=inode
  new_node_list(inode)=inode
  old_node_list(inode)=inode
 enddo

 min_coord=1.0D+20
 max_coord=-1.0D+20
 do inode=1,FSI_mesh_type%NumNodesBIG
  do dir=1,3
   test_coord=FSI_mesh_type%NodeBIG(dir,inode)
   if (test_coord.lt.min_coord) then
    min_coord=test_coord
   endif
   if (test_coord.gt.max_coord) then
    max_coord=test_coord
   endif
  enddo
 enddo
 coord_scale=max_coord-min_coord
 if (coord_scale.gt.zero) then
  ! do nothing
 else
  print *,"coord_scale invalid"
  stop
 endif

 print *,"min_coord: ",min_coord
 print *,"max_coord: ",max_coord
 print *,"coord_scale: ",coord_scale

 sort_nodes_flag=1
 call TopDownMergeSort(FSI_mesh_type,coord_scale, &
        sorted_node_list,B_list, &
        FSI_mesh_type%NumNodesBIG, &
        sort_nodes_flag)

  ! sanity check
 do inode=1,FSI_mesh_type%NumNodesBIG-1
  call compare_nodes(FSI_mesh_type, &
        inode,inode+1, &
        sorted_node_list, &
        coord_scale,compare_flag)

  if (compare_flag.eq.1) then ! (j) > (j+1)
   print *,"list not sorted properly"
   stop
  else if ((compare_flag.eq.0).or.(compare_flag.eq.-1)) then
   ! do nothing
  else
   print *,"compare_flag invalid"
   stop
  endif
 enddo ! do inode=1,FSI_mesh_type%NumNodesBIG-1 (sanity check)

 inode=1
 knode=1
 do while (inode.lt.FSI_mesh_type%NumNodesBIG)
  new_node_list(sorted_node_list(inode))=knode
  old_node_list(knode)=sorted_node_list(inode)
  jnode=inode
  call compare_nodes(FSI_mesh_type, &
    jnode,jnode+1, &
    sorted_node_list, &
    coord_scale,compare_flag)
  do while ((compare_flag.eq.0).and. &
            (jnode.lt.FSI_mesh_type%NumNodesBIG))
   jnode=jnode+1
   new_node_list(sorted_node_list(jnode))=knode
   alternate_node_list(sorted_node_list(jnode))=sorted_node_list(inode)
   if (jnode.lt.FSI_mesh_type%NumNodesBIG) then
    call compare_nodes(FSI_mesh_type, &
     jnode,jnode+1, &
     sorted_node_list, &
     coord_scale,compare_flag)
   else if (jnode.eq.FSI_mesh_type%NumNodesBIG) then
    compare_flag=0
   else 
    print *,"jnode invalid"
    stop
   endif
  enddo
  inode=jnode+1
  knode=knode+1
 enddo
 if (inode.eq.FSI_mesh_type%NumNodesBIG) then
  new_node_list(sorted_node_list(inode))=knode
  old_node_list(knode)=sorted_node_list(inode)
  inode=inode+1
  knode=knode+1
 else if (inode.eq.FSI_mesh_type%NumNodesBIG+1) then
  ! do nothing
 else
  print *,"inode invalid1: inode=",inode
  stop
 endif
 if (inode.eq.FSI_mesh_type%NumNodesBIG+1) then
  ! do nothing
 else
  print *,"inode invalid2: inode=",inode
  stop
 endif

 num_nodes_local=knode-1

 allocate(NodeBIG_local(3,num_nodes_local))
 allocate(NodeVelBIG_local(3,num_nodes_local))
 allocate(NodeForceBIG_local(3,num_nodes_local))
 allocate(NodeDensityBIG_local(num_nodes_local))
 allocate(NodeMassBIG_local(num_nodes_local))
 allocate(NodeTempBIG_local(num_nodes_local))

 do ielem=1,FSI_mesh_type%NumIntElemsBIG
  nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,ielem)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem.ne.3 not supported"
   stop
  endif
  ! root (parent) element id = IntElemBIG(4,elemid)
  do ilocal=1,nodes_per_elem
   inode=FSI_mesh_type%IntElemBIG(ilocal,ielem)
   FSI_mesh_type%IntElemBIG(ilocal,ielem)=new_node_list(inode)
  enddo
 enddo ! do ielem=1,FSI_mesh_type%NumIntElemsBIG

 print *,"old nodes ",FSI_mesh_type%NumNodesBIG
 print *,"new nodes ",num_nodes_local

 do knode=1,num_nodes_local
  inode=old_node_list(knode)

  do dir=1,3
   NodeBIG_local(dir,knode)=FSI_mesh_type%NodeBIG(dir,inode)
   NodeVelBIG_local(dir,knode)=FSI_mesh_type%NodeVelBIG(dir,inode)
   NodeForceBIG_local(dir,knode)=FSI_mesh_type%NodeForceBIG(dir,inode)
  enddo
  NodeDensityBIG_local(knode)=FSI_mesh_type%NodeDensityBIG(inode)
  NodeMassBIG_local(knode)=FSI_mesh_type%NodeMassBIG(inode)
  NodeTempBIG_local(knode)=FSI_mesh_type%NodeTempBIG(inode)
 enddo !knode=1,num_nodes_local

 deallocate(FSI_mesh_type%NodeBIG)
 deallocate(FSI_mesh_type%NodeVelBIG)
 deallocate(FSI_mesh_type%NodeForceBIG)
 deallocate(FSI_mesh_type%NodeDensityBIG)
 deallocate(FSI_mesh_type%NodeMassBIG)
 deallocate(FSI_mesh_type%NodeTempBIG)
 deallocate(FSI_mesh_type%NodeNormalBIG)
 deallocate(FSI_mesh_type%NodeNormalEdgeBIG)
 deallocate(FSI_mesh_type%ElemNodeCountBIG)
 deallocate(FSI_mesh_type%ElemNodeCountEdgeBIG)

 allocate(FSI_mesh_type%NodeBIG(3,num_nodes_local))
 allocate(FSI_mesh_type%NodeVelBIG(3,num_nodes_local))
 allocate(FSI_mesh_type%NodeForceBIG(3,num_nodes_local))
 allocate(FSI_mesh_type%NodeDensityBIG(num_nodes_local))
 allocate(FSI_mesh_type%NodeMassBIG(num_nodes_local))
 allocate(FSI_mesh_type%NodeTempBIG(num_nodes_local))
 allocate(FSI_mesh_type%NodeNormalBIG(3,num_nodes_local))
 allocate(FSI_mesh_type%NodeNormalEdgeBIG(3,num_nodes_local))
 allocate(FSI_mesh_type%ElemNodeCountBIG(num_nodes_local))
 allocate(FSI_mesh_type%ElemNodeCountEdgeBIG(num_nodes_local))

 do knode=1,num_nodes_local
  do dir=1,3
   FSI_mesh_type%NodeBIG(dir,knode)=NodeBIG_local(dir,knode)
   FSI_mesh_type%NodeVelBIG(dir,knode)=NodeVelBIG_local(dir,knode)
   FSI_mesh_type%NodeForceBIG(dir,knode)=NodeForceBIG_local(dir,knode)
   FSI_mesh_type%NodeNormalBIG(dir,knode)=zero
   FSI_mesh_type%NodeNormalEdgeBIG(dir,knode)=zero
  enddo
  FSI_mesh_type%NodeDensityBIG(knode)=NodeDensityBIG_local(knode)
  FSI_mesh_type%NodeMassBIG(knode)=NodeMassBIG_local(knode)
  FSI_mesh_type%NodeTempBIG(knode)=NodeTempBIG_local(knode)
  FSI_mesh_type%ElemNodeCountBIG(knode)=0
  FSI_mesh_type%ElemNodeCountEdgeBIG(knode)=0
 enddo !knode=1,num_nodes_local

 FSI_mesh_type%NumNodesBIG=num_nodes_local

 deallocate(NodeBIG_local)
 deallocate(NodeVelBIG_local)
 deallocate(NodeForceBIG_local)
 deallocate(NodeDensityBIG_local)
 deallocate(NodeMassBIG_local)
 deallocate(NodeTempBIG_local)
 
 deallocate(old_node_list)
 deallocate(new_node_list)
 deallocate(sorted_node_list)
 deallocate(B_list)
 deallocate(alternate_node_list)

end subroutine remove_duplicate_nodes

subroutine print_edge( &
 FSI_mesh_type,edit_refined_data,ielem,inode,edge_data)
IMPLICIT NONE
INTEGER_T, INTENT(in) :: edit_refined_data
type(mesh_type), INTENT(in) :: FSI_mesh_type
INTEGER_T, INTENT(in) :: ielem
INTEGER_T, INTENT(in) :: inode
REAL_T, INTENT(out) :: edge_data(6)
INTEGER_T :: inodep1
INTEGER_T :: local_nodes_per_elem
INTEGER_T :: dir
REAL_T :: x1(3)
REAL_T :: x2(3)

 inodep1=inode+1
 if (edit_refined_data.eq.0) then
  local_nodes_per_elem=FSI_mesh_type%ElemData(1,ielem)
  if (inode.eq.local_nodes_per_elem) then
   inodep1=1
  endif 
  do dir=1,3
   x1(dir)=FSI_mesh_type%Node(dir,FSI_mesh_type%IntElem(inode,ielem))
   x2(dir)=FSI_mesh_type%Node(dir,FSI_mesh_type%IntElem(inodep1,ielem))
  enddo
 else if (edit_refined_data.eq.1) then
  local_nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,ielem)
  if (inode.eq.local_nodes_per_elem) then
   inodep1=1
  endif 
  do dir=1,3
   x1(dir)=FSI_mesh_type%NodeBIG(dir,FSI_mesh_type%IntElemBIG(inode,ielem))
   x2(dir)=FSI_mesh_type%NodeBIG(dir,FSI_mesh_type%IntElemBIG(inodep1,ielem))
  enddo
 else
  print *,"edit_refined_data invalid"
  stop
 endif
 do dir=1,3
  edge_data(dir)=x1(dir)
  edge_data(dir+3)=x2(dir)
 enddo

 print *,"print_edge:"
 print *,"ielem,inode ",ielem,inode
 print *,"x1 ",x1(1),x1(2),x1(3)
 print *,"x2 ",x2(1),x2(2),x2(3)

end subroutine print_edge

! if an opposite element cannot be found for a given edge, then there
! is a "hanging node" which must be removed by either strategically
! splitting or merging selective elements.
subroutine init_EdgeNormal( &
   FSI_mesh_type,part_id,max_part_id,ioproc,isout,edit_refined_data, &
   generate_time)
IMPLICIT NONE
REAL_T, INTENT(in) :: generate_time
INTEGER_T, INTENT(in) :: edit_refined_data
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
type(mesh_type), INTENT(inout) :: FSI_mesh_type
INTEGER_T, INTENT(in) :: ioproc,isout
INTEGER_T, allocatable :: sorted_edge_list(:)
INTEGER_T, allocatable :: B_list(:)
INTEGER_T, allocatable :: edge_inode(:)
INTEGER_T :: sort_nodes_flag
INTEGER_T :: compare_flag
REAL_T :: min_coord
REAL_T :: max_coord
REAL_T :: coord_scale
REAL_T :: test_coord
REAL_T, allocatable :: xnode(:,:)
REAL_T :: normal(3)
REAL_T :: normal_opp(3)
REAL_T :: local_normal(3)
REAL_T :: mag
REAL_T :: overlap_size
INTEGER_T :: dir
INTEGER_T :: doubly_flag
INTEGER_T :: edge_id
INTEGER_T :: old_edge_id
INTEGER_T :: old_edge_id_opp
INTEGER_T :: old_edge_id_node
INTEGER_T :: old_edge_id_node_opp
INTEGER_T :: iedge
INTEGER_T :: jedge
INTEGER_T :: ielem
INTEGER_T :: ielem_opp
INTEGER_T :: inode
INTEGER_T :: inode_opp
INTEGER_T :: inodep1
INTEGER_T :: local_nodes_per_elem
INTEGER_T :: nodes_per_elem
INTEGER_T :: nelems
INTEGER_T :: num_equal
REAL_T :: edgej(6)
REAL_T :: edgejp1(6)
REAL_T :: old_edge_data(6)
REAL_T :: cur_edge_data(6)
REAL_T :: opp_edge_data(6)

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI_mesh_type%part_id.ne.part_id) then
  print *,"FSI_mesh_type%part_id.ne.part_id"
  stop
 endif
 if (edit_refined_data.eq.0) then
  nelems=FSI_mesh_type%NumIntElems
  nodes_per_elem=FSI_mesh_type%IntElemDim
  if (nelems.eq.UBOUND(FSI_mesh_type%IntElem,2)) then
   ! do nothing
  else
   print *,"nelems invalid"
   stop
  endif
  if (nodes_per_elem.eq.UBOUND(FSI_mesh_type%IntElem,1)) then
   ! do nothing
  else
   print *,"nodes_per_elem invalid"
   stop
  endif
 else if (edit_refined_data.eq.1) then
  nelems=FSI_mesh_type%NumIntElemsBIG
  nodes_per_elem=3
  if (nelems.eq.UBOUND(FSI_mesh_type%IntElemBIG,2)) then
   ! do nothing
  else
   print *,"nelems invalid"
   stop
  endif
  ! root (parent) element id = IntElemBIG(4,elemid)
  if (nodes_per_elem+1.eq.UBOUND(FSI_mesh_type%IntElemBIG,1)) then
   ! do nothing
  else
   print *,"nodes_per_elem invalid"
   stop
  endif
 else
  print *,"edit_refined_data invalid"
  stop
 endif

 allocate(FSI_mesh_type%edge_endpoints(6,nodes_per_elem*nelems))
 allocate(FSI_mesh_type%edge_ielem(nodes_per_elem*nelems))
 ! node id of first point on the edge.
 allocate(edge_inode(nodes_per_elem*nelems))
 allocate(sorted_edge_list(nodes_per_elem*nelems))
 allocate(B_list(nodes_per_elem*nelems))
 allocate(xnode(nodes_per_elem,3))

 min_coord=1.0D+20
 max_coord=-1.0D+20
 edge_id=0

 do ielem=1,nelems
  if (edit_refined_data.eq.0) then
   do dir=1,3*nodes_per_elem
     !average of adjoining element normals.
    FSI_mesh_type%EdgeNormal(dir,ielem)=zero
   enddo
   do dir=1,nodes_per_elem
    FSI_mesh_type%EdgeElemId(dir,ielem)=0
    FSI_mesh_type%EdgeElemIdNode(dir,ielem)=0
   enddo
   local_nodes_per_elem=FSI_mesh_type%ElemData(1,ielem)
   if (local_nodes_per_elem.le.nodes_per_elem) then
    ! do nothing
   else
    print *,"local_nodes_per_elem invalid"
    stop
   endif
   do inode=1,local_nodes_per_elem
    do dir=1,3
     xnode(inode,dir)= &
      FSI_mesh_type%Node(dir,FSI_mesh_type%IntElem(inode,ielem))
    enddo
   enddo
  else if (edit_refined_data.eq.1) then
   do dir=1,3*nodes_per_elem
    !average of adjoining element normals.
    FSI_mesh_type%EdgeNormalBIG(dir,ielem)=zero
   enddo
   do dir=1,nodes_per_elem
    FSI_mesh_type%EdgeElemIdBIG(dir,ielem)=0
    FSI_mesh_type%EdgeElemIdNodeBIG(dir,ielem)=0
   enddo
   local_nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,ielem)
   if (local_nodes_per_elem.eq.3) then
    do inode=1,local_nodes_per_elem
     do dir=1,3
      xnode(inode,dir)= &
       FSI_mesh_type%NodeBIG(dir,FSI_mesh_type%IntElemBIG(inode,ielem))
     enddo
    enddo
   else
    print *,"local_nodes_per_elem invalid"
    stop
   endif
  else
   print *,"edit_refined_data invalid"
   stop
  endif

  do inode=1,local_nodes_per_elem
   if ((inode.ge.1).and.(inode.lt.local_nodes_per_elem)) then
    inodep1=inode+1
   else if (inode.eq.local_nodes_per_elem) then
    inodep1=1
   else
    print *,"inode invalid3: inode=",inode
    stop
   endif
   edge_id=edge_id+1
   if ((edge_id.ge.1).and.(edge_id.le.nodes_per_elem*nelems)) then
    do dir=1,3
     FSI_mesh_type%edge_endpoints(dir,edge_id)=xnode(inode,dir)
     FSI_mesh_type%edge_endpoints(dir+3,edge_id)=xnode(inodep1,dir)
    enddo
    FSI_mesh_type%edge_ielem(edge_id)=ielem
    edge_inode(edge_id)=inode
    sorted_edge_list(edge_id)=edge_id
    B_list(edge_id)=edge_id
   else
    print *,"edge_id invalid"
    stop
   endif

   do dir=1,3
    test_coord=xnode(inode,dir)
    if (test_coord.lt.min_coord) then
     min_coord=test_coord
    endif
    if (test_coord.gt.max_coord) then
     max_coord=test_coord
    endif
   enddo
  enddo !inode=1,local_nodes_per_elem
 enddo !ielem=1,nelems

 coord_scale=max_coord-min_coord
 if (coord_scale.gt.zero) then
  ! do nothing
 else
  print *,"coord_scale invalid"
  stop
 endif
 print *,"min_coord(edgelist): ",min_coord
 print *,"max_coord(edgelist): ",max_coord
 print *,"coord_scale(edgelist): ",coord_scale

 sort_nodes_flag=0
 call TopDownMergeSort(FSI_mesh_type,coord_scale, &
        sorted_edge_list,B_list, &
        edge_id, &
        sort_nodes_flag)

  ! sanity check
 do iedge=1,edge_id-1
  do dir=1,6
   edgej(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(iedge))
   edgejp1(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(iedge+1))
  enddo
  call compare_edge(edgej,edgejp1,coord_scale,compare_flag,overlap_size)

  if (compare_flag.eq.1) then ! (j) > (j+1)
   print *,"list not sorted properly"
   stop
  else if ((compare_flag.eq.0).or.(compare_flag.eq.-1)) then
   ! do nothing
  else
   print *,"compare_flag invalid"
   stop
  endif
 enddo ! do iedge=1,edge_id-1 (sanity check)

 iedge=1
 do while (iedge.lt.edge_id)
  num_equal=0
  jedge=iedge
  compare_flag=0

  do while ((compare_flag.eq.0).and.(jedge.lt.edge_id)) 

   do dir=1,6
    edgej(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(jedge))
    edgejp1(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(jedge+1))
   enddo
   call compare_edge(edgej,edgejp1,coord_scale,compare_flag,overlap_size)
   if (compare_flag.eq.0) then
    jedge=jedge+1
    num_equal=num_equal+1
   else if (compare_flag.eq.-1)then
    ! do nothing
   else
    print *,"compare_flag invalid"
    stop
   endif

  enddo ! do while ((compare_flag.eq.0).and.(jedge.lt.edge_id)) 

  if (num_equal.eq.0) then
   ielem=FSI_mesh_type%edge_ielem(sorted_edge_list(iedge))
   ielem_opp=ielem
   inode=edge_inode(sorted_edge_list(iedge))
   inode_opp=inode

   if (edit_refined_data.eq.0) then
    old_edge_id=FSI_mesh_type%EdgeElemId(inode,ielem)
    old_edge_id_opp=old_edge_id
    old_edge_id_node=FSI_mesh_type%EdgeElemIdNode(inode,ielem)
    old_edge_id_node_opp=old_edge_id_node
    doubly_flag=FSI_mesh_type%ElemData(DOUBLYCOMP,ielem)
   else if (edit_refined_data.eq.1) then
    old_edge_id=FSI_mesh_type%EdgeElemIdBIG(inode,ielem)
    old_edge_id_opp=old_edge_id
    old_edge_id_node=FSI_mesh_type%EdgeElemIdNodeBIG(inode,ielem)
    old_edge_id_node_opp=old_edge_id_node
    doubly_flag=FSI_mesh_type%ElemDataBIG(DOUBLYCOMP,ielem)
   else
    print *,"edit_refined_data invalid"
    stop
   endif
   if (doubly_flag.eq.1) then
    ! do nothing
   else if (doubly_flag.eq.0) then
    if (1.eq.0) then
     print *,"Warning:cannot have a hanging edge for a singly wetted element"
    endif
   else
    print *,"doubly_flag invalid"
    stop
   endif
  else if (num_equal.eq.1) then
   ielem=FSI_mesh_type%edge_ielem(sorted_edge_list(iedge))
   ielem_opp=FSI_mesh_type%edge_ielem(sorted_edge_list(iedge+1))
   inode=edge_inode(sorted_edge_list(iedge))
   inode_opp=edge_inode(sorted_edge_list(iedge+1))
   if (edit_refined_data.eq.0) then
    old_edge_id=FSI_mesh_type%EdgeElemId(inode,ielem)
    old_edge_id_opp=FSI_mesh_type%EdgeElemId(inode_opp,ielem_opp)
    old_edge_id_node=FSI_mesh_type%EdgeElemIdNode(inode,ielem)
    old_edge_id_node_opp=FSI_mesh_type%EdgeElemIdNode(inode_opp,ielem_opp)
   else if (edit_refined_data.eq.1) then
    old_edge_id=FSI_mesh_type%EdgeElemIdBIG(inode,ielem)
    old_edge_id_opp=FSI_mesh_type%EdgeElemIdBIG(inode_opp,ielem_opp)
    old_edge_id_node=FSI_mesh_type%EdgeElemIdNodeBIG(inode,ielem)
    old_edge_id_node_opp=FSI_mesh_type%EdgeElemIdNodeBIG(inode_opp,ielem_opp)
   else
    print *,"edit_refined_data invalid"
    stop
   endif
  else
   print *,"num_equal invalid"
   stop
  endif

  if (edit_refined_data.eq.0) then
   if ((old_edge_id.eq.0).and. &
       (old_edge_id_opp.eq.0).and. &
       (old_edge_id_node.eq.0).and. &
       (old_edge_id_node_opp.eq.0)) then
    FSI_mesh_type%EdgeElemId(inode,ielem)=ielem_opp
    FSI_mesh_type%EdgeElemId(inode_opp,ielem_opp)=ielem
    FSI_mesh_type%EdgeElemIdNode(inode,ielem)=inode_opp
    FSI_mesh_type%EdgeElemIdNode(inode_opp,ielem_opp)=inode
    call scinormal(ielem,normal, &
      FSI_mesh_type,part_id,max_part_id, &
      generate_time)
    call scinormal(ielem_opp,normal_opp, &
      FSI_mesh_type,part_id,max_part_id, &
      generate_time)
    do dir=1,3
     local_normal(dir)=half*(normal(dir)+normal_opp(dir))
    enddo

    mag=zero
    do dir=1,3
     mag=mag+local_normal(dir)**2
    enddo
    mag=sqrt(mag)
    if (mag.gt.zero) then
     do dir=1,3
      local_normal(dir)=local_normal(dir)/mag
     enddo
    else if (mag.eq.zero) then
     do dir=1,3
      local_normal(dir)=zero
     enddo
    else
     print *,"mag bust"
     stop
    endif

    do dir=1,3
     FSI_mesh_type%EdgeNormal(3*(inode-1)+dir,ielem)= &
             local_normal(dir)
     FSI_mesh_type%EdgeNormal(3*(inode_opp-1)+dir,ielem_opp)= &
             local_normal(dir)
    enddo

   else
    print *,"old_edge_id,old_edge_id_opp,old_edge_id_node,node_opp bad0"
    print *,"old_edge_id ",old_edge_id
    print *,"old_edge_id_opp ",old_edge_id_opp
    print *,"old_edge_id_node ",old_edge_id_node
    print *,"old_edge_id_node_opp ",old_edge_id_node_opp
    print *,"ielem ",ielem
    print *,"ielem_opp ",ielem_opp
    print *,"inode ",inode
    print *,"inode_opp ",inode_opp
    call print_edge(FSI_mesh_type,edit_refined_data,old_edge_id, &
     old_edge_id_node,old_edge_data)
    call print_edge(FSI_mesh_type,edit_refined_data,ielem,inode,cur_edge_data)
    call print_edge(FSI_mesh_type,edit_refined_data,ielem_opp,inode_opp, &
     opp_edge_data)
    call compare_edge(old_edge_data,cur_edge_data,coord_scale, &
      compare_flag,overlap_size)
    print *,"compare_flag old,cur ",compare_flag
    call compare_edge(cur_edge_data,old_edge_data,coord_scale, &
      compare_flag,overlap_size)
    print *,"compare_flag cur,old ",compare_flag
    call compare_edge(opp_edge_data,cur_edge_data,coord_scale, &
      compare_flag,overlap_size)
    print *,"compare_flag opp,cur ",compare_flag
    print *,"iedge=",iedge
    do dir=1,6
     edgej(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(iedge))
     edgejp1(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(iedge+1))
     print *,"dir,edgej,edgejp1 ",dir,edgej(dir),edgejp1(dir)
    enddo
    stop
   endif

  else if (edit_refined_data.eq.1) then

   if ((old_edge_id.eq.0).and. &
       (old_edge_id_opp.eq.0).and. &
       (old_edge_id_node.eq.0).and. &
       (old_edge_id_node_opp.eq.0)) then
    FSI_mesh_type%EdgeElemIdBIG(inode,ielem)=ielem_opp
    FSI_mesh_type%EdgeElemIdBIG(inode_opp,ielem_opp)=ielem
    FSI_mesh_type%EdgeElemIdNodeBIG(inode,ielem)=inode_opp
    FSI_mesh_type%EdgeElemIdNodeBIG(inode_opp,ielem_opp)=inode
    call scinormalBIG(ielem,normal, &
      FSI_mesh_type,part_id,max_part_id, &
      generate_time)
    call scinormalBIG(ielem_opp,normal_opp, &
      FSI_mesh_type,part_id,max_part_id, &
      generate_time)
    do dir=1,3
     local_normal(dir)=half*(normal(dir)+normal_opp(dir))
    enddo

    mag=zero
    do dir=1,3
     mag=mag+local_normal(dir)**2
    enddo
    mag=sqrt(mag)
    if (mag.gt.zero) then
     do dir=1,3
      local_normal(dir)=local_normal(dir)/mag
     enddo
    else if (mag.eq.zero) then
     do dir=1,3
      local_normal(dir)=zero
     enddo
    else
     print *,"mag bust"
     stop
    endif

    do dir=1,3
     FSI_mesh_type%EdgeNormalBIG(3*(inode-1)+dir,ielem)= &
             local_normal(dir)
     FSI_mesh_type%EdgeNormalBIG(3*(inode_opp-1)+dir,ielem_opp)= &
             local_normal(dir)
    enddo
   else
    print *,"old_edge_id,old_edge_id_opp,old_edge_id_node,node_opp bad1"
    print *,"old_edge_id ",old_edge_id
    print *,"old_edge_id_opp ",old_edge_id_opp
    print *,"old_edge_id_node ",old_edge_id_node
    print *,"old_edge_id_node_opp ",old_edge_id_node_opp
    print *,"ielem ",ielem
    print *,"ielem_opp ",ielem_opp
    print *,"inode ",inode
    print *,"inode_opp ",inode_opp
    do dir=1,6
     edgej(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(iedge))
     edgejp1(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(iedge+1))
     print *,"dir,edgej,edgejp1 ",dir,edgej(dir),edgejp1(dir)
    enddo
    stop
   endif

  else
   print *,"edit_refined_data invalid"
   stop
  endif

  if ((1.eq.0).and.(num_equal.eq.1)) then
   print *,"START PAIRING INFO ---------------------------------"
   print *,"iedge: ",iedge
   print *,"ielem ",ielem
   print *,"ielem_opp ",ielem_opp
   print *,"inode ",inode
   print *,"inode_opp ",inode_opp
   call print_edge(FSI_mesh_type,edit_refined_data,ielem,inode,cur_edge_data)
   call print_edge(FSI_mesh_type,edit_refined_data,ielem_opp,inode_opp, &
    opp_edge_data)
   call compare_edge(opp_edge_data,cur_edge_data,coord_scale, &
    compare_flag,overlap_size)
   print *,"compare_flag opp,cur ",compare_flag
   do dir=1,6
    edgej(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(iedge))
    edgejp1(dir)=FSI_mesh_type%edge_endpoints(dir,sorted_edge_list(iedge+1))
    print *,"dir,edgej,edgejp1 ",dir,edgej(dir),edgejp1(dir)
   enddo
   print *,"END PAIRING INFO ---------------------------------"
  endif 

  iedge=iedge+num_equal+1

 enddo !while (iedge.lt.edge_id)

 deallocate(FSI_mesh_type%edge_endpoints)
 deallocate(FSI_mesh_type%edge_ielem)
 deallocate(edge_inode)
 deallocate(sorted_edge_list)
 deallocate(B_list)
 deallocate(xnode)

end subroutine init_EdgeNormal

! split triangles so that size is no bigger than "h_small"
subroutine post_process_nodes_elements(initflag, &
  problo,probhi, &
  FSI_mesh_type, &
  part_id, &
  max_part_id, &
  ioproc,isout,h_small)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
type(mesh_type), INTENT(inout) :: FSI_mesh_type

REAL_T, INTENT(in) :: problo(3),probhi(3)
INTEGER_T, INTENT(in) :: initflag
INTEGER_T, INTENT(in) :: ioproc,isout
REAL_T, INTENT(in) :: h_small
INTEGER_T :: edit_refined_data
INTEGER_T :: ielem,nodes_per_elem,dir
INTEGER_T :: ilev_lag
INTEGER_T :: inode_list
INTEGER_T :: inode_elem
INTEGER_T :: inode_elem_p1
INTEGER_T :: new_NumIntElems
REAL_T :: generate_time
REAL_T, dimension(3) :: x1,x2,x3
REAL_T, dimension(3) :: vel1,vel2,vel3
REAL_T, dimension(6) :: force1,force2,force3,forcesplit
REAL_T, dimension(3) :: xsplit
REAL_T, dimension(3) :: velsplit
REAL_T :: biggest_h
REAL_T :: smallest_h
INTEGER_T :: first_measure
REAL_T :: temp_h
REAL_T :: mass1,mass2,mass3,mass_split
REAL_T :: new_massL,new_massR
REAL_T :: volL,volR
REAL_T :: den1,den2,den3,den_split
REAL_T :: d12_2D,d23_2D,d13_2D
REAL_T :: d12_3D,d23_3D,d13_3D
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
REAL_T    :: mag
INTEGER_T :: view_refined
INTEGER_T, allocatable :: DoublyWettedNode(:)

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI_mesh_type%part_id.ne.part_id) then
  print *,"FSI_mesh_type%part_id.ne.part_id"
  stop
 endif
 if (h_small.gt.zero) then
  ! do nothing
 else
  print *,"h_small invalid"
  stop
 endif

 if ((ioproc.ne.1).and.(ioproc.ne.0)) then
  print *,"ioproc invalid"
  stop
 endif

 generate_time=zero

 local_refine_factor=FSI_mesh_type%refine_factor

 if ((initflag.eq.1).and.(isout.eq.1)) then
  if (ioproc.eq.1) then
   print *,"in post_process_nodes_elements..."
  endif
 endif

 edit_refined_data=0
 call init_EdgeNormal(FSI_mesh_type,part_id,max_part_id,ioproc,isout, &
    edit_refined_data,generate_time)

! UPDATE NODENORMAL ON ORIGINAL LAGRANGIAN GRID --------------

 allocate(DoublyWettedNode(FSI_mesh_type%NumNodes))

 do inode_list=1,FSI_mesh_type%NumNodes
  do dir=1,3
   FSI_mesh_type%NodeNormal(dir,inode_list)=0.0
   FSI_mesh_type%NodeNormalEdge(dir,inode_list)=0.0
  enddo
  FSI_mesh_type%ElemNodeCount(inode_list)=0
  FSI_mesh_type%ElemNodeCountEdge(inode_list)=0
  DoublyWettedNode(inode_list)=1
 enddo ! inode_list=1,FSI_mesh_type%NumNodes

 do ielem=1,FSI_mesh_type%NumIntElems
  nodes_per_elem=FSI_mesh_type%ElemData(1,ielem)
  if (nodes_per_elem.lt.3) then
   print *,"nodes_per_elem<3 not supported"
   stop
  endif
  call scinormal(ielem,normal, &
     FSI_mesh_type,part_id,max_part_id, &
     generate_time)
  do inode_elem=1,nodes_per_elem
   inode_list=FSI_mesh_type%IntElem(inode_elem,ielem)

   if ((inode_elem.ge.1).and.(inode_elem.lt.nodes_per_elem)) then
    inode_elem_p1=inode_elem+1
   else if (inode_elem.eq.nodes_per_elem) then
    inode_elem_p1=1
   else
    print *,"inode_elem invalid4: inode_elem=",inode_elem
    stop
   endif

   if (FSI_mesh_type%ElemData(DOUBLYCOMP,ielem).eq.0) then 
    DoublyWettedNode(inode_list)=0
   else if (FSI_mesh_type%ElemData(DOUBLYCOMP,ielem).eq.1) then 
    ! do nothing
   else
    print *,"ElemData(DOUBLYCOMP) invalid"
    stop
   endif

   do dir=1,3
    FSI_mesh_type%NodeNormal(dir,inode_list)= &
      FSI_mesh_type%NodeNormal(dir,inode_list)+normal(dir)
    FSI_mesh_type%NodeNormalEdge(dir,inode_list)= &
      FSI_mesh_type%NodeNormalEdge(dir,inode_list)+ &
       FSI_mesh_type%EdgeNormal(3*(inode_elem-1)+dir,ielem)
    FSI_mesh_type%NodeNormalEdge(dir,inode_list)= &
      FSI_mesh_type%NodeNormalEdge(dir,inode_list)+ &
       FSI_mesh_type%EdgeNormal(3*(inode_elem_p1-1)+dir,ielem)
   enddo !dir=1...3
   FSI_mesh_type%ElemNodeCount(inode_list)= &
     FSI_mesh_type%ElemNodeCount(inode_list)+1
   FSI_mesh_type%ElemNodeCountEdge(inode_list)= &
     FSI_mesh_type%ElemNodeCountEdge(inode_list)+2
  enddo ! inode_elem=1,nodes_per_elem
 enddo ! do ielem=1,FSI_mesh_type%NumIntElems

 do inode_list=1,FSI_mesh_type%NumNodes

  normal_cnt=FSI_mesh_type%ElemNodeCount(inode_list)
  if (normal_cnt.gt.0) then

   if (normal_cnt.ge.3) then
    ! do nothing
   else

    if (DoublyWettedNode(inode_list).eq.1) then
     ! do nothing
    else
     print *,"Warning(crse): expecting normal_cnt>=3; normal_cnt=",normal_cnt
    endif

   endif

   mag=zero
   do dir=1,3
    FSI_mesh_type%NodeNormal(dir,inode_list)= &
     FSI_mesh_type%NodeNormal(dir,inode_list)/normal_cnt
    mag=mag+FSI_mesh_type%NodeNormal(dir,inode_list)**2
   enddo ! dir=1,3
   mag=sqrt(mag)
   if (mag.gt.zero) then
    do dir=1,3
     FSI_mesh_type%NodeNormal(dir,inode_list)= &
       FSI_mesh_type%NodeNormal(dir,inode_list)/mag
    enddo
   else if (mag.eq.zero) then
    ! do nothing
   else
    print *,"mag invalid"
    stop
   endif 
  else if (normal_cnt.eq.0) then
   print *,"node has no elements connected"
   stop

   do dir=1,3
    if (FSI_mesh_type%NodeNormal(dir,inode_list).eq.zero) then
     ! do nothing
    else
     print *,"FSI_mesh_type%NodeNormal(dir,inode_list) invalid"
     stop
    endif
   enddo
  else
   print *,"normal_cnt invalid"
   stop
  endif

  normal_cnt=FSI_mesh_type%ElemNodeCountEdge(inode_list)
  if (normal_cnt.gt.0) then

   if (normal_cnt.ge.6) then
    ! do nothing
   else

    if (DoublyWettedNode(inode_list).eq.1) then
     ! do nothing
    else
     print *,"Warning(crse): expecting normal_cnt>=6; normal_cnt=",normal_cnt
    endif

   endif

   mag=zero
   do dir=1,3
    FSI_mesh_type%NodeNormalEdge(dir,inode_list)= &
     FSI_mesh_type%NodeNormalEdge(dir,inode_list)/normal_cnt
    mag=mag+FSI_mesh_type%NodeNormalEdge(dir,inode_list)**2
   enddo ! dir=1,3
   mag=sqrt(mag)
   if (mag.gt.zero) then
    do dir=1,3
     FSI_mesh_type%NodeNormalEdge(dir,inode_list)= &
       FSI_mesh_type%NodeNormalEdge(dir,inode_list)/mag
    enddo
   else if (mag.eq.zero) then
    ! do nothing
   else
    print *,"mag invalid"
    stop
   endif 
  else if (normal_cnt.eq.0) then
   print *,"node has no elements connected"
   stop

   do dir=1,3
    if (FSI_mesh_type%NodeNormalEdge(dir,inode_list).eq.zero) then
     ! do nothing
    else
     print *,"FSI_mesh_type%NodeNormalEdge(dir,inode_list) invalid"
     stop
    endif
   enddo
  else
   print *,"normal_cnt invalid"
   stop
  endif

 enddo ! do inode_list=1,FSI_mesh_type%NumNodes

 deallocate(DoublyWettedNode)

 view_refined=0
 call tecplot_normals(FSI_mesh_type,part_id,max_part_id,view_refined)

 do dir=1,3
  FSI_mesh_type%center_BB(dir)=zero
  FSI_mesh_type%exterior_BB(dir,1)=zero
  FSI_mesh_type%exterior_BB(dir,2)=zero
  FSI_mesh_type%interior_BB(dir,1)=zero
  FSI_mesh_type%interior_BB(dir,2)=zero
 enddo

 do inode_list=1,FSI_mesh_type%NumNodes
  do dir=1,3
   x1(dir)=FSI_mesh_type%Node(dir,inode_list)
   FSI_mesh_type%center_BB(dir)=FSI_mesh_type%center_BB(dir)+x1(dir)
   if (inode_list.eq.1) then
    FSI_mesh_type%exterior_BB(dir,1)=x1(dir)
    FSI_mesh_type%exterior_BB(dir,2)=x1(dir)
   endif
   if (x1(dir).lt.FSI_mesh_type%exterior_BB(dir,1)) then
    FSI_mesh_type%exterior_BB(dir,1)=x1(dir)
   endif
   if (x1(dir).gt.FSI_mesh_type%exterior_BB(dir,2)) then
    FSI_mesh_type%exterior_BB(dir,2)=x1(dir)
   endif
  enddo !dir=1..3
 enddo ! inode_list=1,FSI_mesh_type%NumNodes

 if (FSI_mesh_type%NumNodes.gt.0) then
  do dir=1,3
   FSI_mesh_type%center_BB(dir)= &
       FSI_mesh_type%center_BB(dir)/FSI_mesh_type%NumNodes
   if ((FSI_mesh_type%center_BB(dir).ge. &
        FSI_mesh_type%exterior_BB(dir,1)).and. &
       (FSI_mesh_type%center_BB(dir).le. &
        FSI_mesh_type%exterior_BB(dir,2))) then
    ! do nothing
   else
    print *,"center_BB invalid: dir,center_BB: ", &
      dir,FSI_mesh_type%center_BB(dir)
    stop
   endif
   FSI_mesh_type%interior_BB(dir,1)=FSI_mesh_type%exterior_BB(dir,1)
   FSI_mesh_type%interior_BB(dir,2)=FSI_mesh_type%exterior_BB(dir,2)
  enddo !dir=1..3

  do inode_list=1,FSI_mesh_type%NumNodes
   do dir=1,3
    x1(dir)=FSI_mesh_type%Node(dir,inode_list)
    if (x1(dir).ge.FSI_mesh_type%center_BB(dir)) then
     FSI_mesh_type%interior_BB(dir,2)= &
          min(FSI_mesh_type%interior_BB(dir,2),x1(dir))
    endif
    if (x1(dir).le.FSI_mesh_type%center_BB(dir)) then
     FSI_mesh_type%interior_BB(dir,1)= &
          max(FSI_mesh_type%interior_BB(dir,1),x1(dir))
    endif
   enddo !dir=1..3
  enddo ! inode_list=1,FSI_mesh_type%NumNodes
 else if (FSI_mesh_type%NumNodes.eq.0) then
  ! do nothing
 else
  print *,"NumNodes invalid"
  stop
 endif

 if (ioproc.eq.1) then
  do dir=1,3
   print *,"part_id,dir,center_BB ",part_id,dir,FSI_mesh_type%center_BB(dir)
   print *,"part_id,dir,exterior_BB(dir,1) ", &
           part_id,dir,FSI_mesh_type%exterior_BB(dir,1)
   print *,"part_id,dir,exterior_BB(dir,2) ", &
           part_id,dir,FSI_mesh_type%exterior_BB(dir,2)
   print *,"part_id,dir,interior_BB(dir,1) ", &
           part_id,dir,FSI_mesh_type%interior_BB(dir,1)
   print *,"part_id,dir,interior_BB(dir,2) ", &
           part_id,dir,FSI_mesh_type%interior_BB(dir,2)
  enddo !dir=1..3
 else if (ioproc.eq.0) then
  ! do nothing
 else
  print *,"ioproc invalid"
  stop
 endif 
   
! START OF LAGRANGIAN REFINEMENT SECTION ---------------------

 biggest_h=0.0
 smallest_h=0.0
 first_measure=0

 new_NumIntElems=0
 do ielem=1,FSI_mesh_type%NumIntElems
  nodes_per_elem=FSI_mesh_type%ElemData(1,ielem)
  do isub=0,nodes_per_elem-3
   new_NumIntElems=new_NumIntElems+1
   node1=FSI_mesh_type%IntElem(1,ielem)
   node2=FSI_mesh_type%IntElem(nodes_per_elem-isub-1,ielem)
   node3=FSI_mesh_type%IntElem(nodes_per_elem-isub,ielem)
   if ((node1.gt.FSI_mesh_type%NumNodes).or. &
       (node2.gt.FSI_mesh_type%NumNodes).or. &
       (node3.gt.FSI_mesh_type%NumNodes).or. &
       (node1.lt.1).or. &
       (node2.lt.1).or. &
       (node3.lt.1)) then
    print *,"node1,node2, or node3 invalid ",node1,node2,node3
    print *,"ielem=",ielem
    print *,"FSI_mesh_type%NumNodes ",FSI_mesh_type%NumNodes
    print *,"FSI_mesh_type%NumIntElems ",FSI_mesh_type%NumIntElems
    print *,"nodes_per_elem=",nodes_per_elem
    print *,"part_id=",part_id
    stop
   endif
    
   do dir=1,3
    x1(dir)=FSI_mesh_type%Node(dir,node1)
    x2(dir)=FSI_mesh_type%Node(dir,node2)
    x3(dir)=FSI_mesh_type%Node(dir,node3)
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

   call xdist_project(x1,x2, &
      FSI_mesh_type,part_id,max_part_id, &
      d12_2D,d12_3D)
   call xdist_project(x2,x3, &
      FSI_mesh_type,part_id,max_part_id, &
      d23_2D,d23_3D)
   call xdist_project(x3,x1, &
      FSI_mesh_type,part_id,max_part_id, &
      d13_2D,d13_3D)
   if ((d12_3D.ge.CTMLunderflow).and. &
       (d12_3D.le.CTMLoverflow).and. &
       (d23_3D.ge.CTMLunderflow).and. &
       (d23_3D.le.CTMLoverflow).and. &
       (d13_3D.ge.CTMLunderflow).and. &
       (d13_3D.le.CTMLoverflow)) then
    ! do nothing
   else
    print *,"ielem ",ielem
    print *,"isub= ",isub
    print *,"d12_3D,d23_3D,d13_3D out of range: ",d12_3D,d23_3D,d13_3D
    print *,"part_id,node1,node2,node3 ", &
      part_id,node1,node2,node3
    do dir=1,3
     print *,"dir,x1,x2,x3 ",dir,x1(dir),x2(dir),x3(dir)
    enddo
    stop
   endif

   if (first_measure.eq.0) then
    biggest_h=d12_2D
    smallest_h=d12_3D
   endif
   first_measure=1

   if (biggest_h.lt.d12_2D) then
    biggest_h=d12_2D
   endif
   if (biggest_h.lt.d23_2D) then
    biggest_h=d23_2D
   endif
   if (biggest_h.lt.d13_2D) then
    biggest_h=d13_2D
   endif

   if (smallest_h.gt.d12_3D) then
    smallest_h=d12_3D
   endif
   if (smallest_h.gt.d23_3D) then
    smallest_h=d23_3D
   endif
   if (smallest_h.gt.d13_3D) then
    smallest_h=d13_3D
   endif

  enddo ! isub=0..nodes_per_elem-3
 enddo ! ielem=1..FSI_mesh_type%NumIntElems

 FSI_mesh_type%max_side_len=biggest_h  
 FSI_mesh_type%min_side_len=smallest_h  
 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"part_id,max_side_len,min_side_len ",part_id, &
   FSI_mesh_type%max_side_len,FSI_mesh_type%min_side_len
 endif

 if (local_refine_factor.eq.0) then
  n_lag_levels=2
 else if (local_refine_factor.ge.1) then
  n_lag_levels=1
  temp_h=biggest_h
  do while (temp_h.gt.(local_refine_factor-0.01d0)*h_small)
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
 multi_lag(1)%n_nodes=FSI_mesh_type%NumNodes
 multi_lag(1)%n_elems=new_NumIntElems
 allocate(multi_lag(1)%nd(3,FSI_mesh_type%NumNodes))
 allocate(multi_lag(1)%ndvel(3,FSI_mesh_type%NumNodes))
 allocate(multi_lag(1)%ndforce(6,FSI_mesh_type%NumNodes))
 allocate(multi_lag(1)%ndmass(FSI_mesh_type%NumNodes))
 allocate(multi_lag(1)%nddensity(FSI_mesh_type%NumNodes))
 allocate(multi_lag(1)%ndtemp(FSI_mesh_type%NumNodes))
 allocate(multi_lag(1)%elemdt(flags_per_element,new_NumIntElems))
  !root (parent) element id = intelemdt(4,elemid)
 allocate(multi_lag(1)%intelemdt(4,new_NumIntElems))
 do inode_list=1,FSI_mesh_type%NumNodes
  do dir=1,3
   multi_lag(1)%nd(dir,inode_list)=FSI_mesh_type%Node(dir,inode_list)
   multi_lag(1)%ndvel(dir,inode_list)=FSI_mesh_type%NodeVel(dir,inode_list)
  enddo
  multi_lag(1)%ndmass(inode_list)=FSI_mesh_type%NodeMass(inode_list)
  multi_lag(1)%nddensity(inode_list)=FSI_mesh_type%NodeDensity(inode_list)
  do dir=1,3
   multi_lag(1)%ndforce(dir,inode_list)=FSI_mesh_type%NodeForce(dir,inode_list)
  enddo
  multi_lag(1)%ndtemp(inode_list)=FSI_mesh_type%NodeTemp(inode_list)
 enddo
  
 new_NumIntElems=0
 do ielem=1,FSI_mesh_type%NumIntElems
  nodes_per_elem=FSI_mesh_type%ElemData(1,ielem)
  do isub=0,nodes_per_elem-3
   new_NumIntElems=new_NumIntElems+1
   multi_lag(1)%intelemdt(1,new_NumIntElems)=FSI_mesh_type%IntElem(1,ielem)
   multi_lag(1)%intelemdt(2,new_NumIntElems)= &
     FSI_mesh_type%IntElem(nodes_per_elem-isub-1,ielem)
   multi_lag(1)%intelemdt(3,new_NumIntElems)= &
     FSI_mesh_type%IntElem(nodes_per_elem-isub,ielem)
   ! root (parent) element id = intelemdt(4,elemid)
   multi_lag(1)%intelemdt(4,new_NumIntElems)=ielem
    ! FSI_mesh_type%ElemData(DOUBLYCOMP,ielem) is the doubly wetted flag
   do dir=1,flags_per_element
    multi_lag(1)%elemdt(dir,new_NumIntElems)=FSI_mesh_type%ElemData(dir,ielem)
   enddo
    ! always 3 nodes per element for refined surface
   multi_lag(1)%elemdt(1,new_NumIntElems)=3
  enddo ! isub=0..nodes_per_elem-3
 enddo ! ielem=1..FSI_mesh_type%NumIntElems

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"numintelems,numnodes ", &
      FSI_mesh_type%NumIntElems,FSI_mesh_type%NumNodes
  print *,"numintelems(3node),numnodes ", &
      new_NumIntElems,FSI_mesh_type%NumNodes
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
    allocate(multi_lag(ilevel+1)%nddensity(save_n_nodes))
    allocate(multi_lag(ilevel+1)%ndtemp(save_n_nodes))
    allocate(multi_lag(ilevel+1)%elemdt(flags_per_element,save_n_elems))
     ! root (parent) element id = intelemdt(4,elemid)
    allocate(multi_lag(ilevel+1)%intelemdt(4,save_n_elems))
    do inode_list=1,multi_lag(ilevel)%n_nodes
     do dir=1,3
      multi_lag(ilevel+1)%nd(dir,inode_list)= &
        multi_lag(ilevel)%nd(dir,inode_list)
      multi_lag(ilevel+1)%ndvel(dir,inode_list)= &
        multi_lag(ilevel)%ndvel(dir,inode_list)
     enddo
     multi_lag(ilevel+1)%ndmass(inode_list)= &
        multi_lag(ilevel)%ndmass(inode_list)
     multi_lag(ilevel+1)%nddensity(inode_list)= &
        multi_lag(ilevel)%nddensity(inode_list)
     do dir=1,3
      multi_lag(ilevel+1)%ndforce(dir,inode_list)= &
        multi_lag(ilevel)%ndforce(dir,inode_list)
     enddo
     multi_lag(ilevel+1)%ndtemp(inode_list)= &
       multi_lag(ilevel)%ndtemp(inode_list)
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
    den1=multi_lag(ilevel_current)%nddensity(node1)
    den2=multi_lag(ilevel_current)%nddensity(node2)
    den3=multi_lag(ilevel_current)%nddensity(node3)
    do dir=1,3
     force1(dir)=multi_lag(ilevel)%ndforce(dir,node1)
     force2(dir)=multi_lag(ilevel)%ndforce(dir,node2)
     force3(dir)=multi_lag(ilevel)%ndforce(dir,node3)
    enddo
    temp1=multi_lag(ilevel)%ndtemp(node1)
    temp2=multi_lag(ilevel)%ndtemp(node2)
    temp3=multi_lag(ilevel)%ndtemp(node3)
    call xdist_project(x1,x2, &
      FSI_mesh_type,part_id,max_part_id, &
      d12_2D,d12_3D)
    call xdist_project(x2,x3, &
      FSI_mesh_type,part_id,max_part_id, &
      d23_2D,d23_3D)
    call xdist_project(x3,x1, &
      FSI_mesh_type,part_id,max_part_id, &
      d13_2D,d13_3D)

    if ((local_refine_factor.gt.0).and. &
        (d12_2D.ge.(local_refine_factor-0.01d0)*h_small).and. &
        (d12_2D.ge.d23_2D).and. &
        (d12_2D.ge.d13_2D)) then
     do dir=1,3
      xsplit(dir)=0.5d0*(x1(dir)+x2(dir))
      velsplit(dir)=0.5d0*(vel1(dir)+vel2(dir))
     enddo
     do dir=1,3
      forcesplit(dir)=0.5d0*(force1(dir)+force2(dir))
     enddo

     tempsplit=0.5d0*(temp1+temp2)
     den_split=0.5d0*(den1+den2)
       !  x1     x_split     x2
       ! den1 * vol1R = mass1R   
       ! den2 * vol2L = mass2L   
       ! 0.5d0 * den12 * (vol1R+vol2L) = mass12
     call get_new_half_vols(x1,x2,xsplit,volL,volR)

     if (FSI_mesh_type%CTML_flag.eq.1) then
      mass_split=0.5d0*den_split*(volL+volR)
      new_massL=mass1-0.5d0*den_split*volL
      new_massR=mass2-0.5d0*den_split*volR
     else if (FSI_mesh_type%CTML_flag.eq.0) then
      mass_split=0.5d0*(mass1+mass2)
      new_massL=mass1
      new_massR=mass2
     else
      print *,"FSI_mesh_type%CTML_flag invalid"
      stop
     endif

     if ((new_massL.gt.0.0d0).and.(new_massR.gt.0.0d0)) then
      ! do nothing
     else
      print *,"new_massL or new_massR invalid"
      stop
     endif

     multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel+1)%n_nodes+1
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+2
     nsplit=multi_lag(ilevel+1)%n_nodes
     esplit=multi_lag(ilevel+1)%n_elems

     if (iter.eq.2) then
      do dir=1,3
       multi_lag(ilevel+1)%nd(dir,nsplit)=xsplit(dir)
       multi_lag(ilevel+1)%ndvel(dir,nsplit)=velsplit(dir)
      enddo 
      do dir=1,3
       multi_lag(ilevel+1)%ndforce(dir,nsplit)=forcesplit(dir)
      enddo
      multi_lag(ilevel+1)%nddensity(nsplit)=den_split
      multi_lag(ilevel+1)%ndmass(nsplit)=mass_split
      multi_lag(ilevel+1)%ndmass(node1)=new_massL
      multi_lag(ilevel+1)%ndmass(node2)=new_massR

      multi_lag(ilevel+1)%ndtemp(nsplit)=tempsplit

      ! root (parent) element id = intelemdt(4,elemid)
      base_ielem=multi_lag(ilevel)%intelemdt(4,ielem)
      multi_lag(ilevel+1)%intelemdt(1,esplit-1)=node1
      multi_lag(ilevel+1)%intelemdt(2,esplit-1)=nsplit
      multi_lag(ilevel+1)%intelemdt(3,esplit-1)=node3
      ! root (parent) element id = intelemdt(4,elemid)
      multi_lag(ilevel+1)%intelemdt(4,esplit-1)=base_ielem
      ! elemdt(DOUBLYCOMP,ielem) is the doubly wetted flag
      do dir=1,flags_per_element
       multi_lag(ilevel+1)%elemdt(dir,esplit-1)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo

      multi_lag(ilevel+1)%intelemdt(1,esplit)=nsplit
      multi_lag(ilevel+1)%intelemdt(2,esplit)=node2
      multi_lag(ilevel+1)%intelemdt(3,esplit)=node3
      ! root (parent) element id = intelemdt(4,elemid)
      multi_lag(ilevel+1)%intelemdt(4,esplit)=base_ielem
      ! elemdt(DOUBLYCOMP,ielem) is the doubly wetted flag
      do dir=1,flags_per_element
       multi_lag(ilevel+1)%elemdt(dir,esplit)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo
     endif ! iter.eq.2
    else if ((local_refine_factor.gt.0).and. &
             (d23_2D.ge.(local_refine_factor-0.01d0)*h_small).and. &
             (d23_2D.ge.d12_2D).and. &
             (d23_2D.ge.d13_2D)) then
     do dir=1,3
      xsplit(dir)=0.5d0*(x2(dir)+x3(dir))
      velsplit(dir)=0.5d0*(vel2(dir)+vel3(dir))
     enddo
     do dir=1,3
      forcesplit(dir)=0.5d0*(force2(dir)+force3(dir))
     enddo
 
     tempsplit=0.5d0*(temp2+temp3)
     den_split=0.5d0*(den2+den3)

     call get_new_half_vols(x2,x3,xsplit,volL,volR)

     if (FSI_mesh_type%CTML_flag.eq.1) then
      mass_split=0.5d0*den_split*(volL+volR)
      new_massL=mass2-0.5d0*den_split*volL
      new_massR=mass3-0.5d0*den_split*volR
     else if (FSI_mesh_type%CTML_flag.eq.0) then
      mass_split=0.5d0*(mass2+mass3)
      new_massL=mass2
      new_massR=mass3
     else
      print *,"FSI_mesh_type%CTML_flag invalid"
      stop
     endif

     if ((new_massL.gt.0.0d0).and.(new_massR.gt.0.0d0)) then
      ! do nothing
     else
      print *,"new_massL or new_massR invalid"
      stop
     endif

     multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel+1)%n_nodes+1
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+2
     nsplit=multi_lag(ilevel+1)%n_nodes
     esplit=multi_lag(ilevel+1)%n_elems

     if (iter.eq.2) then
      do dir=1,3
       multi_lag(ilevel+1)%nd(dir,nsplit)=xsplit(dir)
       multi_lag(ilevel+1)%ndvel(dir,nsplit)=velsplit(dir)
      enddo 
      do dir=1,3
       multi_lag(ilevel+1)%ndforce(dir,nsplit)=forcesplit(dir)
      enddo
      multi_lag(ilevel+1)%nddensity(nsplit)=den_split
      multi_lag(ilevel+1)%ndmass(nsplit)=mass_split
      multi_lag(ilevel+1)%ndmass(node2)=new_massL
      multi_lag(ilevel+1)%ndmass(node3)=new_massR

      multi_lag(ilevel+1)%ndtemp(nsplit)=tempsplit

      ! root (parent) element id = intelemdt(4,elemid)
      base_ielem=multi_lag(ilevel)%intelemdt(4,ielem)
      multi_lag(ilevel+1)%intelemdt(1,esplit-1)=node1
      multi_lag(ilevel+1)%intelemdt(2,esplit-1)=node2
      multi_lag(ilevel+1)%intelemdt(3,esplit-1)=nsplit
      multi_lag(ilevel+1)%intelemdt(4,esplit-1)=base_ielem
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit-1)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo

      ! root (parent) element id = intelemdt(4,elemid)
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
             (d13_2D.ge.(local_refine_factor-0.01d0)*h_small).and. &
             (d13_2D.ge.d12_2D).and. &
             (d13_2D.ge.d23_2D)) then
     do dir=1,3
      xsplit(dir)=0.5d0*(x1(dir)+x3(dir))
      velsplit(dir)=0.5d0*(vel1(dir)+vel3(dir))
     enddo
     do dir=1,3
      forcesplit(dir)=0.5d0*(force1(dir)+force3(dir))
     enddo

     tempsplit=0.5d0*(temp1+temp3)
     den_split=0.5d0*(den1+den3)

     call get_new_half_vols(x1,x3,xsplit,volL,volR)

     if (FSI_mesh_type%CTML_flag.eq.1) then
      mass_split=0.5d0*den_split*(volL+volR)
      new_massL=mass1-0.5d0*den_split*volL
      new_massR=mass3-0.5d0*den_split*volR
     else if (FSI_mesh_type%CTML_flag.eq.0) then
      mass_split=0.5d0*(mass1+mass3)
      new_massL=mass1
      new_massR=mass3
     else
      print *,"FSI_mesh_type%CTML_flag invalid"
      stop
     endif

     if ((new_massL.gt.0.0d0).and.(new_massR.gt.0.0d0)) then
      ! do nothing
     else
      print *,"new_massL or new_massR invalid"
      stop
     endif

     multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel+1)%n_nodes+1
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+2
     nsplit=multi_lag(ilevel+1)%n_nodes
     esplit=multi_lag(ilevel+1)%n_elems

     if (iter.eq.2) then
      do dir=1,3
       multi_lag(ilevel+1)%nd(dir,nsplit)=xsplit(dir)
       multi_lag(ilevel+1)%ndvel(dir,nsplit)=velsplit(dir)
      enddo 
      do dir=1,3
       multi_lag(ilevel+1)%ndforce(dir,nsplit)=forcesplit(dir)
      enddo
      multi_lag(ilevel+1)%nddensity(nsplit)=den_split
      multi_lag(ilevel+1)%ndmass(nsplit)=mass_split
      multi_lag(ilevel+1)%ndmass(node1)=new_massL
      multi_lag(ilevel+1)%ndmass(node3)=new_massR

      multi_lag(ilevel+1)%ndtemp(nsplit)=tempsplit

      ! root (parent) element id = intelemdt(4,elemid)
      base_ielem=multi_lag(ilevel)%intelemdt(4,ielem)
      multi_lag(ilevel+1)%intelemdt(1,esplit-1)=node1
      multi_lag(ilevel+1)%intelemdt(2,esplit-1)=node2
      multi_lag(ilevel+1)%intelemdt(3,esplit-1)=nsplit
      multi_lag(ilevel+1)%intelemdt(4,esplit-1)=base_ielem
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit-1)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo

      ! root (parent) element id = intelemdt(4,elemid)
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
      ! root (parent) element id = intelemdt(4,elemid)
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

 FSI_mesh_type%NumNodesBIG=multi_lag(n_lag_levels)%n_nodes
 FSI_mesh_type%NumIntElemsBIG=multi_lag(n_lag_levels)%n_elems

  ! in: post_process_nodes_elements
 if (initflag.eq.0) then 
  deallocate(FSI_mesh_type%NodeBIG)
  deallocate(FSI_mesh_type%NodeVelBIG)
  deallocate(FSI_mesh_type%NodeForceBIG)
  deallocate(FSI_mesh_type%NodeMassBIG)
  deallocate(FSI_mesh_type%NodeDensityBIG)
  deallocate(FSI_mesh_type%NodeTempBIG)
  deallocate(FSI_mesh_type%NodeNormalBIG)
  deallocate(FSI_mesh_type%NodeNormalEdgeBIG)
  deallocate(FSI_mesh_type%ElemNodeCountBIG)
  deallocate(FSI_mesh_type%ElemNodeCountEdgeBIG)

  deallocate(FSI_mesh_type%ElemDataXnotBIG)
  deallocate(FSI_mesh_type%ElemDataBIG)
  deallocate(FSI_mesh_type%IntElemBIG)
  deallocate(FSI_mesh_type%EdgeNormalBIG)
  deallocate(FSI_mesh_type%EdgeElemIdBIG)
  deallocate(FSI_mesh_type%EdgeElemIdNodeBIG)

 endif

  ! in: post_process_nodes_elements
 allocate(FSI_mesh_type%NodeBIG(3,FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%NodeNormalBIG(3,FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%NodeNormalEdgeBIG(3,FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%NodeVelBIG(3,FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%NodeForceBIG(3,FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%NodeMassBIG(FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%NodeDensityBIG(FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%NodeTempBIG(FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%ElemNodeCountBIG(FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%ElemNodeCountEdgeBIG(FSI_mesh_type%NumNodesBIG))
 allocate(FSI_mesh_type%ElemDataXnotBIG(3,FSI_mesh_type%NumIntElemsBIG))
 allocate(FSI_mesh_type%ElemDataBIG(flags_per_element, &
          FSI_mesh_type%NumIntElemsBIG))
  ! root (parent) element id = IntElemBIG(4,elemid)
 allocate(FSI_mesh_type%IntElemBIG(4,FSI_mesh_type%NumIntElemsBIG))
 !average of adjoining element normals.
 allocate(FSI_mesh_type%EdgeNormalBIG(9,FSI_mesh_type%NumIntElemsBIG))
 allocate(FSI_mesh_type%EdgeElemIdBIG(3,FSI_mesh_type%NumIntElemsBIG))
 allocate(FSI_mesh_type%EdgeElemIdNodeBIG(3,FSI_mesh_type%NumIntElemsBIG))

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"NumNodes, NumIntElems ",FSI_mesh_type%NumNodes, &
   FSI_mesh_type%NumIntElems
  print *,"NumNodesBIG, NumIntElemsBIG ",FSI_mesh_type%NumNodesBIG, &
   FSI_mesh_type%NumIntElemsBIG
 endif

  ! in: post_process_nodes_elements
 do ielem=1,FSI_mesh_type%NumIntElemsBIG
  ! root (parent) element id = intelemdt(4,elemid)
  do dir=1,4
   FSI_mesh_type%IntElemBIG(dir,ielem)= &
     multi_lag(n_lag_levels)%intelemdt(dir,ielem)
  enddo
   ! ElemDataBIG(DOUBLYCOMP,ielem) is the doubly wetted flag
  do dir=1,flags_per_element
   FSI_mesh_type%ElemDataBIG(dir,ielem)= &
       multi_lag(n_lag_levels)%elemdt(dir,ielem)
  enddo
 enddo
  ! in: post_process_nodes_elements
 do inode_list=1,FSI_mesh_type%NumNodesBIG
  do dir=1,3
   FSI_mesh_type%NodeNormalBIG(dir,inode_list)=0.0
   FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)=0.0
   FSI_mesh_type%NodeBIG(dir,inode_list)= &
      multi_lag(n_lag_levels)%nd(dir,inode_list)
   FSI_mesh_type%NodeVelBIG(dir,inode_list)= &
      multi_lag(n_lag_levels)%ndvel(dir,inode_list)
  enddo
  do dir=1,3
   FSI_mesh_type%NodeForceBIG(dir,inode_list)= &
     multi_lag(n_lag_levels)%ndforce(dir,inode_list)
  enddo
  FSI_mesh_type%NodeMassBIG(inode_list)= &
     multi_lag(n_lag_levels)%ndmass(inode_list)
  FSI_mesh_type%NodeDensityBIG(inode_list)= &
     multi_lag(n_lag_levels)%nddensity(inode_list)
  FSI_mesh_type%NodeTempBIG(inode_list)= &
     multi_lag(n_lag_levels)%ndtemp(inode_list)
  FSI_mesh_type%ElemNodeCountBIG(inode_list)=0
  FSI_mesh_type%ElemNodeCountEdgeBIG(inode_list)=0
 enddo ! inode_list=1,FSI_mesh_type%NumNodesBIG

 do ilev_lag=1,n_lag_levels
  deallocate(multi_lag(ilev_lag)%intelemdt)
  deallocate(multi_lag(ilev_lag)%elemdt)
  deallocate(multi_lag(ilev_lag)%nd)
  deallocate(multi_lag(ilev_lag)%ndvel)
  deallocate(multi_lag(ilev_lag)%ndforce)
  deallocate(multi_lag(ilev_lag)%ndmass)
  deallocate(multi_lag(ilev_lag)%nddensity)
  deallocate(multi_lag(ilev_lag)%ndtemp)
 enddo
 deallocate(multi_lag)

 if ((ioproc.eq.1).and.(isout.eq.1)) then 
  print *,"prior to remove_duplicate_nodes"
 endif
 
 call remove_duplicate_nodes(FSI_mesh_type,part_id,max_part_id)

 edit_refined_data=1
 call init_EdgeNormal(FSI_mesh_type,part_id,max_part_id,ioproc,isout, &
          edit_refined_data,generate_time)

 if ((ioproc.eq.1).and.(isout.eq.1)) then 
  print *,"creating normals for the nodes and Xnot"
 endif

 ! UPDATE NODENORMALBIG ON REFINED LAGRANGIAN GRID --------------

 allocate(DoublyWettedNode(FSI_mesh_type%NumNodesBIG))

 do inode_list=1,FSI_mesh_type%NumNodesBIG
  do dir=1,3
   FSI_mesh_type%NodeNormalBIG(dir,inode_list)=0.0
   FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)=0.0
  enddo
  FSI_mesh_type%ElemNodeCountBIG(inode_list)=0
  FSI_mesh_type%ElemNodeCountEdgeBIG(inode_list)=0
  DoublyWettedNode(inode_list)=1
 enddo ! inode_list=1,FSI_mesh_type%NumNodesBIG

 do ielem=1,FSI_mesh_type%NumIntElemsBIG
  nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,ielem)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem invalid"
   stop
  endif
  do dir=1,3
   FSI_mesh_type%ElemDataXnotBIG(dir,ielem)=0.0
  enddo
  call scinormalBIG(ielem,normal, &
     FSI_mesh_type,part_id,max_part_id, &
     generate_time)
  do inode_elem=1,nodes_per_elem
   inode_list=FSI_mesh_type%IntElemBIG(inode_elem,ielem)

   if ((inode_elem.ge.1).and.(inode_elem.lt.nodes_per_elem)) then
    inode_elem_p1=inode_elem+1
   else if (inode_elem.eq.nodes_per_elem) then
    inode_elem_p1=1
   else
    print *,"inode_elem invalid5: inode_elem=",inode_elem
    stop
   endif

   if (FSI_mesh_type%ElemDataBIG(DOUBLYCOMP,ielem).eq.0) then 
    DoublyWettedNode(inode_list)=0
   else if (FSI_mesh_type%ElemDataBIG(DOUBLYCOMP,ielem).eq.1) then 
    ! do nothing
   else
    print *,"ElemDataBIG(DOUBLYCOMP) invalid"
    stop
   endif

   do dir=1,3
    FSI_mesh_type%NodeNormalBIG(dir,inode_list)= &
      FSI_mesh_type%NodeNormalBIG(dir,inode_list)+normal(dir)
    FSI_mesh_type%ElemDataXnotBIG(dir,ielem)= &
      FSI_mesh_type%ElemDataXnotBIG(dir,ielem)+ &
      FSI_mesh_type%NodeBIG(dir,inode_list)

    FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)= &
      FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)+ &
       FSI_mesh_type%EdgeNormalBIG(3*(inode_elem-1)+dir,ielem)
    FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)= &
      FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)+ &
       FSI_mesh_type%EdgeNormalBIG(3*(inode_elem_p1-1)+dir,ielem)
   enddo
   FSI_mesh_type%ElemNodeCountBIG(inode_list)= &
     FSI_mesh_type%ElemNodeCountBIG(inode_list)+1
   FSI_mesh_type%ElemNodeCountEdgeBIG(inode_list)= &
     FSI_mesh_type%ElemNodeCountEdgeBIG(inode_list)+2
  enddo ! inode_elem=1,nodes_per_elem
  do dir=1,3
   FSI_mesh_type%ElemDataXnotBIG(dir,ielem)= &
    FSI_mesh_type%ElemDataXnotBIG(dir,ielem)/3.0
  enddo
 enddo ! do ielem=1,FSI_mesh_type%NumIntElemsBIG

 do inode_list=1,FSI_mesh_type%NumNodesBIG
  normal_cnt=FSI_mesh_type%ElemNodeCountBIG(inode_list)
  if (normal_cnt.gt.0) then

   if (normal_cnt.ge.3) then
    ! do nothing
   else

    if (DoublyWettedNode(inode_list).eq.1) then
     ! do nothing
    else
     print *,"warning(BIG): expecting normal_cnt>=3; normal_cnt=",normal_cnt
    endif

   endif

   mag=zero
   do dir=1,3
    FSI_mesh_type%NodeNormalBIG(dir,inode_list)= &
     FSI_mesh_type%NodeNormalBIG(dir,inode_list)/normal_cnt
    mag=mag+FSI_mesh_type%NodeNormalBIG(dir,inode_list)**2
   enddo ! dir=1,3
   mag=sqrt(mag)
   if (mag.gt.zero) then
    do dir=1,3
     FSI_mesh_type%NodeNormalBIG(dir,inode_list)= &
       FSI_mesh_type%NodeNormalBIG(dir,inode_list)/mag
    enddo
   else if (mag.eq.zero) then
    ! do nothing
   else
    print *,"mag invalid"
    stop
   endif 
  else if (normal_cnt.eq.0) then
   do dir=1,3
    if (FSI_mesh_type%NodeNormalBIG(dir,inode_list).eq.zero) then
     ! do nothing
    else
     print *,"FSI_mesh_type%NodeNormalBIG(dir,inode_list) invalid"
     stop
    endif
   enddo
  else
   print *,"normal_cnt invalid"
   stop
  endif

  normal_cnt=FSI_mesh_type%ElemNodeCountEdgeBIG(inode_list)
  if (normal_cnt.gt.0) then

   if (normal_cnt.ge.6) then
    ! do nothing
   else

    if (DoublyWettedNode(inode_list).eq.1) then
     ! do nothing
    else
     print *,"warning(BIG): expecting normal_cnt>=6, normal_cnt=",normal_cnt
    endif

   endif

   mag=zero
   do dir=1,3
    FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)= &
     FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)/normal_cnt
    mag=mag+FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)**2
   enddo ! dir=1,3
   mag=sqrt(mag)
   if (mag.gt.zero) then
    do dir=1,3
     FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)= &
       FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list)/mag
    enddo
   else if (mag.eq.zero) then
    ! do nothing
   else
    print *,"mag invalid"
    stop
   endif 
  else if (normal_cnt.eq.0) then
   print *,"node has no elements connected"
   stop

   do dir=1,3
    if (FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list).eq.zero) then
     ! do nothing
    else
     print *,"FSI_mesh_type%NodeNormalEdgeBIG(dir,inode_list) invalid"
     stop
    endif
   enddo
  else
   print *,"normal_cnt invalid"
   stop
  endif

 enddo ! do inode_list=1,FSI_mesh_type%NumNodesBIG

 deallocate(DoublyWettedNode)

 biggest_h=0.0
 smallest_h=0.0
 first_measure=0

 do ielem=1,FSI_mesh_type%NumIntElemsBIG
  nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,ielem)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem invalid"
   stop
  endif
  node1=FSI_mesh_type%IntElemBIG(1,ielem)
  node2=FSI_mesh_type%IntElemBIG(2,ielem)
  node3=FSI_mesh_type%IntElemBIG(3,ielem)
  do dir=1,3
   x1(dir)=FSI_mesh_type%NodeBIG(dir,node1) 
   x2(dir)=FSI_mesh_type%NodeBIG(dir,node2)
   x3(dir)=FSI_mesh_type%NodeBIG(dir,node3)

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

  call xdist_project(x1,x2, &
    FSI_mesh_type,part_id,max_part_id, &
    d12_2D,d12_3D)
  call xdist_project(x2,x3, &
    FSI_mesh_type,part_id,max_part_id, &
    d23_2D,d23_3D)
  call xdist_project(x3,x1, &
    FSI_mesh_type,part_id,max_part_id, &
    d13_2D,d13_3D)

  if (first_measure.eq.0) then
   biggest_h=d12_2D
   smallest_h=d12_3D
  endif
  first_measure=1

  if (biggest_h.lt.d12_2D) then
   biggest_h=d12_2D
  endif
  if (biggest_h.lt.d23_2D) then
   biggest_h=d23_2D
  endif
  if (biggest_h.lt.d13_2D) then
   biggest_h=d13_2D
  endif

  if (smallest_h.gt.d12_3D) then
   smallest_h=d12_3D
  endif
  if (smallest_h.gt.d23_3D) then
   smallest_h=d23_3D
  endif
  if (smallest_h.gt.d13_3D) then
   smallest_h=d13_3D
  endif

 enddo ! ielem=1..FSI_mesh_type%NumIntElemsBIG

 FSI_mesh_type%max_side_len_refined=biggest_h  
 FSI_mesh_type%min_side_len_refined=smallest_h  

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"part_id,flag_2D_to_3D ",part_id,FSI_mesh_type%flag_2D_to_3D
  print *,"part_id,max_side_len_refined,min_side_len_refined ",part_id, &
   FSI_mesh_type%max_side_len_refined,FSI_mesh_type%min_side_len_refined
  print *,"local_refine_factor ",local_refine_factor
  print *,"h_small ",h_small
 endif

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print*,"in allocate: part ID ",FSI_mesh_type%part_id
 endif

 view_refined=1
 call tecplot_normals(FSI_mesh_type,part_id,max_part_id,view_refined)

! END OF LAGRANGIAN REFINEMENT SECTION -----------------------

return
end subroutine post_process_nodes_elements


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

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: isout
INTEGER_T :: inode
INTEGER_T, INTENT(in) :: ioproc
REAL_T, dimension(3) :: maxnode,minnode
REAL_T, dimension(3) :: xval,xval1,xval2,xvalm1,xvalp1
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
INTEGER_T :: dir

REAL_T, dimension(3) :: xxblob1,newxxblob1,xxblob2,newxxblob2
REAL_T, dimension(3) :: vel_local
REAL_T, dimension(NCOMP_FORCE_STRESS) :: force_local
REAL_T :: mass_local
REAL_T :: density_local
REAL_T :: volm1,volp1
REAL_T :: radradblob1,radradblob2
INTEGER_T :: stand_alone_flag
INTEGER_T :: orig_nodes
INTEGER_T :: ctml_part_id

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id out of range, part_id, TOTAL_NPARTS:",part_id,TOTAL_NPARTS
   stop
  endif
  if (FSI(part_id)%part_id.ne.part_id) then
   print *,"FSI(part_id)%part_id.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
   print *,"ioproc invalid"
   stop
  endif

  ctml_part_id=ctml_part_id_map(part_id)

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
    maxnode(dir)=-1.0e+10
    minnode(dir)=1.0e+10
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   if (FSI(part_id)%flag_2D_to_3D.eq.1) then
    orig_nodes=FSI(part_id)%NumNodes/2
    if (orig_nodes*2.eq.FSI(part_id)%NumNodes) then
     ! do nothing
    else
     print *,"orig_nodes invalid"
     stop
    endif
   else if (FSI(part_id)%flag_2D_to_3D.eq.0) then
    orig_nodes=FSI(part_id)%NumNodes
   else
    print *,"FSI(part_id)%flag_2D_to_3D invalid"
    stop
   endif


   do inode=1,orig_nodes

    do dir=1,3
     xval(dir)=zero
     xvalm1(dir)=zero
     xvalp1(dir)=zero
     vel_local(dir)=zero
    enddo

    do dir=1,AMREX_SPACEDIM
     xval(dir)=ctml_fib_pst(ctml_part_id,inode,dir)
     if (inode.gt.1) then
      xvalm1(dir)=ctml_fib_pst(ctml_part_id,inode-1,dir)
     endif
     if (inode.lt.orig_nodes) then
      xvalp1(dir)=ctml_fib_pst(ctml_part_id,inode+1,dir)
     endif
     vel_local(dir)=ctml_fib_vel(ctml_part_id,inode,dir)
    enddo
    volm1=zero
    volp1=zero
    if (inode.gt.1) then
     volm1=zero
     do dir=1,AMREX_SPACEDIM
      volm1=volm1+(xval(dir)-xvalm1(dir))**2
     enddo
     volm1=half*sqrt(volm1)
    endif 
    if (inode.lt.orig_nodes) then
     volp1=zero
     do dir=1,AMREX_SPACEDIM
      volp1=volp1+(xval(dir)-xvalp1(dir))**2
     enddo
     volp1=half*sqrt(volp1)
    endif

    do dir=1,NCOMP_FORCE_STRESS
     force_local(dir)=zero
    enddo
    do dir=1,AMREX_SPACEDIM
     force_local(dir)=ctml_fib_frc(ctml_part_id,inode,dir)
    enddo
    mass_local=ctml_fib_mass(ctml_part_id,inode)
    if (mass_local.gt.zero) then
     ! do nothing
    else
     print *,"mass_local invalid, mass_local=",mass_local
     stop
    endif
    if (volm1+volp1.gt.zero) then
     if (FSI(part_id)%flag_2D_to_3D.eq.1) then
      density_local=mass_local/(volm1+volp1)
     else
      print *,"do not know how to fine density_local if 3d"
      stop
     endif
    else
     print *,"volm1 or volp1 invalid",volm1,volp1
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
    enddo ! dir=1..3

    do dir=1,3
     FSI(part_id)%Node_old(dir,inode)=xval1(dir)
     FSI(part_id)%Node_new(dir,inode)=xval1(dir)
     FSI(part_id)%Node_current(dir,inode)=xval1(dir)
     FSI(part_id)%Node(dir,inode)=xval1(dir)
     FSI(part_id)%NodeVel(dir,inode)=vel_local(dir)
     FSI(part_id)%NodeVel_old(dir,inode)=vel_local(dir)
     FSI(part_id)%NodeVel_new(dir,inode)=vel_local(dir)
    enddo ! dir=1,3
    do dir=1,NCOMP_FORCE_STRESS
     FSI(part_id)%NodeForce(dir,inode)=force_local(dir)
     FSI(part_id)%NodeForce_old(dir,inode)=force_local(dir)
     FSI(part_id)%NodeForce_new(dir,inode)=force_local(dir)
    enddo

     ! in: CTML_init_sci_node
    FSI(part_id)%NodeMass(inode)=mass_local
    FSI(part_id)%NodeDensity(inode)=density_local

    FSI(part_id)%NodeTemp(inode)=zero
    FSI(part_id)%NodeTemp_new(inode)=zero

    if (FSI(part_id)%flag_2D_to_3D.eq.1) then
     stand_alone_flag=0
     call convert_2D_to_3D_nodes_FSI(part_id,inode,stand_alone_flag)
    else if (FSI(part_id)%flag_2D_to_3D.eq.0) then
     ! do nothing
    else
     print *,"FSI(part_id)%flag_2D_to_3D invalid"
     stop
    endif
       
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

! called from overall_solid_init
! overall_solid_init is called from CLSVOF_ReadHeader
! CLSVOF_ReadHeader is called from fort_headermsg when 
!  FSI_operation.eq.OP_FSI_INITIALIZE_NODES.
! fort_headermsg, when FSI_operation.eq.OP_FSI_INITIALIZE_NODES, is called from
! NavierStokes::ns_header_msg_level
! ns_header_msg_level is called from:
!   NavierStokes::FSI_make_distance 
!   (FSI_operation.eq.OP_FSI_MAKE_DISTANCE(SIGN))
!   NavierStokes::post_restart (FSI_operation.eq.OP_FSI_INITIALIZE_NODES)
!   NavierStokes::initData (FSI_operation.eq.OP_FSI_INITIALIZE_NODES)
!   NavierStokes::nonlinear_advection 
!    (FSI_operation.eq.OP_FSI_LAG_STRESS or OP_FSI_UPDATE_NODES)
! SUMMARY (typical time step):
!   0. for k=0,1,2,3,...
!   1. NavierStokes::nonlinear_advection
!   2. CLSMOF advection
!   3. copy Eulerian force to Lagrangian
!   4. Lagrangian (structure) evolution
!   5. copy Lagrangian velocity to Eulerian.
!   6. viscous solve
!   7. pressure solve
!   8. velocity and pressure at t^{n+1}(k) differ too much from
!      t^{n+1}(k-1)? if yes, go back to step 0.
subroutine CTML_init_sci(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  part_id,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: sdim
INTEGER_T, INTENT(in) :: ifirst
INTEGER_T, INTENT(in) :: isout
INTEGER_T :: iface
INTEGER_T, INTENT(in) :: ioproc
REAL_T, INTENT(in) :: curtime,dt
INTEGER_T :: dir
INTEGER_T, INTENT(in) :: istep
INTEGER_T, INTENT(in) :: istop
INTEGER_T :: ctml_part_id
INTEGER_T :: inode_crit,inode
INTEGER_T :: orig_nodes
INTEGER_T :: orig_elements
INTEGER_T :: local_nodes
INTEGER_T :: local_elements

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id out of range, part_id, TOTAL_NPARTS:",part_id,TOTAL_NPARTS
   stop
  endif
  if (FSI(part_id)%part_id.ne.part_id) then
   print *,"FSI(part_id)%part_id.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
  endif

  ctml_part_id=ctml_part_id_map(part_id)

  if ((ctml_part_id.ge.1).and. &
      (ctml_part_id.le.CTML_NPARTS)) then

   FSI(part_id)%CTML_flag=1

   inode_crit=0  ! node index of first inactive node.
   inode=0
#ifdef MVAHABFSI 
    ! declared in: CTMLFSI.F90
   call CTML_GET_POS_VEL_WT( &
    ctml_fib_pst, & ! =coord_fib
    ctml_fib_vel, & ! =vel_fib
    ctml_fib_mass,& ! =ds_fib
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
     ! all the node masses are positive
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
   print *,"in: CTML_init_sci; define MVAHABFSI"
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
    FSI(part_id)%flag_2D_to_3D=0
    print *,"3D not supported yet"
   else if (AMREX_SPACEDIM.eq.2) then
    FSI(part_id)%flag_2D_to_3D=1
   else
    print *,"dimension bust"
    stop
   endif

   orig_nodes=ctml_n_fib_active_nodes(ctml_part_id)
   orig_elements=ctml_n_fib_active_nodes(ctml_part_id)-1

   if (FSI(part_id)%flag_2D_to_3D.eq.1) then
    FSI(part_id)%NumNodes=orig_nodes*2
    FSI(part_id)%NumIntElems=orig_elements*2
   else if (FSI(part_id)%flag_2D_to_3D.eq.0) then
    FSI(part_id)%NumNodes=orig_nodes
    FSI(part_id)%NumIntElems=orig_elements
   else
    print *,"FSI(part_id)%flag_2D_to_3D invalid"
    stop
   endif
   local_nodes=FSI(part_id)%NumNodes
   local_elements=FSI(part_id)%NumIntElems
   FSI(part_id)%IntElemDim=3
 
    ! in: CTML_init_sci 
   if (ifirst.eq.1) then

    ! allocates and inits: ElemData,EdgeNormal(just allocates),
    !  EdgeElemId(just allocates),IntElem,Node_old,
    !  Node_new,Node_current,NodeVel_old,NodeVel_new,NodeForce_old,
    !  NodeForce_new,NodeTemp_old,NodeTemp_new,NodeMass,NodeDensity
    ! allocate_intelem=1
    ! allocates NodeMass,NodeDensity
    call init_FSI(part_id,1)  

    allocate(FSI(part_id)%Node(3,FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%NodeVel(3,FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%NodeNormal(3,FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%NodeNormalEdge(3,FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%ElemNodeCount(FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%ElemNodeCountEdge(FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%NodeForce(3,FSI(part_id)%NumNodes))
    allocate(FSI(part_id)%NodeTemp(FSI(part_id)%NumNodes))
   else
    print *,"something wrong, ifirst should be 1 here"
    stop
   endif  ! ifirst.eq.1

    ! initialize nodes, elements,... from CTML data.
   call CTML_init_sci_node(ioproc,part_id,isout)

   do iface=1,orig_elements

    if (FSI(part_id)%flag_2D_to_3D.eq.0) then
     print *,"3d not ready yet"
     stop
    else if (FSI(part_id)%flag_2D_to_3D.eq.1) then
     FSI(part_id)%IntElem(1,iface)=iface
     FSI(part_id)%IntElem(2,iface)=iface+1
     FSI(part_id)%IntElem(3,iface)=iface+1+orig_nodes

     FSI(part_id)%ElemData(1,iface)=3   ! number of nodes in element
     FSI(part_id)%ElemData(2,iface)=1   ! part number
     FSI(part_id)%ElemData(DOUBLYCOMP,iface)=1 !doubly wetted (0=singly wetted)

     call convert_2D_to_3D_elements_FSI(part_id,iface)

    else
     print *,"dimension bust"
     stop
    endif

   enddo ! iface=1,orig_elements

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

subroutine abort_sci_clsvof()
IMPLICIT NONE

print *,"normal points from solid to fluid (consistent with tecplot)"
print *,"solid nodes are ordered counter clockwise when viewed from fluid."
print *,"n=v1 x v2   v1=node2-node1  v2=node3-node2"
print *,"initially, hitsign=1 in the fluid but the sign will be"
print *,"switched so that LS>0 in the solid"
stop
return
end subroutine abort_sci_clsvof

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
  if (FSI(part_id)%part_id.ne.part_id) then
   print *,"FSI(part_id)%part_id.ne.part_id"
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
       ! IF probtype=538 and part_id=2, then add 0.01d0 in order to shift the
       ! needle all the way in the domain (newxxblob1(3)=0.01d0)
       ! this shift will be taken into account in the line that has:
       ! "xfoot(3)=xfoot(3)+0.01d0" (get_foot_from_target, get_target_from_foot)
     xxblob1(1)=0.0
     xxblob1(2)=0.0
     xxblob1(3)=0.0
     newxxblob1(1)=0.0
     newxxblob1(2)=0.0
     newxxblob1(3)=0.0

     if ((probtype.eq.538).or.(probtype.eq.541)) then
      if (FSI(part_id)%part_id.eq.2) then !want the whole needle in the domain.
       if ((AMREX_SPACEDIM.eq.3).or.(probtype.eq.538)) then
        newxxblob1(3)=0.01d0
       else if ((AMREX_SPACEDIM.eq.2).and.(probtype.eq.541)) then
        newxxblob1(3)=0.0
       endif
      else if (FSI(part_id)%part_id.eq.1) then
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
    if (FSI(part_id)%part_id.eq.1) then

      if (probtype.eq.5501) then
       dwave="rough.cas"
      else if ((probtype.eq.53).and.(axis_dir.eq.100)) then
       dwave="flat_fan_s.cas"
      else
       dwave="injectorgeom.dat"
      endif

    else if (FSI(part_id)%part_id.eq.2) then
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
     ! allocates NodeMass,NodeDensity
     ! allocate_intelem=1
     call init_FSI(part_id,1)  
    else
     print *,"something wrong, ifirst should be 1 here"
     stop
    endif  ! ifirst.eq.1

    do dir=1,3
     maxnode(dir)=-1.0e+10
     minnode(dir)=1.0e+10
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
     FSI(part_id)%ElemData(DOUBLYCOMP,iface)=0 !singly wetted (1=doubly wetted)
     call abort_sci_clsvof()
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
  if (FSI(part_id)%part_id.ne.part_id) then
   print *,"FSI(part_id)%part_id.ne.part_id"
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

   if (FSI(part_id)%part_id.eq.1 ) then

    dwave="foreWing.sci"
    newxxblob1(1)=0.0  ! 1st wing at a "default" location.
    print *,"init_flapping part_id=1"

   else if (FSI(part_id)%part_id.eq.2) then

    if ((axis_dir.eq.0).or.(axis_dir.eq.2)) then
     print *,"part id invalid"
     stop
    else if (axis_dir.eq.1) then
!     dwave="naca0012.sci"
     dwave="hindWing.sci"
     ! position 2nd wing in "default" location too.
     ! kinematics will shift wing in the x direction.
     newxxblob1(1)=0.0  
     print *,"init_flapping part_id=2"
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
    maxnode(dir)=-1.0e+10
    minnode(dir)=1.0e+10
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
     FSI(part_id)%ElemData(DOUBLYCOMP,iface+quad_counter)=0 ! singly wetted
     call abort_sci_clsvof()
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

!called from: overall_solid_advance and overall_solid_init
!"init_from_cas" is the default routine that is called when
!a probtype is not specifically supported.
subroutine init_from_cas(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  part_id,isout)
use global_utility_module

IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: sdim,ifirst,isout
INTEGER_T, INTENT(in) :: ioproc
REAL_T, INTENT(in) :: curtime,dt
INTEGER_T, INTENT(in) :: istep,istop
INTEGER_T :: inode,iface
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
INTEGER_T :: dir
INTEGER_T :: file_format

REAL_T, dimension(3) :: xxblob1,newxxblob1
REAL_T :: radradblob1
INTEGER_T localElem(3)

character(80) :: discard
character(80) :: points_line
INTEGER_T :: ivtk,dummy_num_nodes_per_elem

REAL_T :: local_nodes(3,3)  ! dir,node num
INTEGER_T :: raw_num_nodes
INTEGER_T :: raw_num_elements
REAL_T, allocatable :: raw_nodes(:,:)
INTEGER_T, allocatable :: raw_elements(:,:)

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id invalid"
   stop
  endif
  if (FSI(part_id)%part_id.ne.part_id) then
   print *,"FSI(part_id)%part_id.ne.part_id"
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
   xxblob1(1)=0.0d0
   xxblob1(2)=0.0d0
   xxblob1(3)=0.0d0

   newxxblob1(1)=0.0d0
   newxxblob1(2)=0.0d0
   newxxblob1(3)=0.0d0
   radradblob1=1.0d0

   denpaddle=one
   dampingpaddle=zero

    !file_format=0 cas format
    !file_format=1 vtk format
   call SUB_OPEN_CASFILE(part_id,14,file_format)

   if (ifirst.ne.1) then
    print *,"ifirst bust"
    stop
   endif
   FSI(part_id)%IntElemDim=3

   do dir=1,3
    maxnode(dir)=-1.0e+10
    minnode(dir)=1.0e+10
    maxnodebefore(dir)=-1.0d+10
    minnodebefore(dir)=1.0d+10
   enddo

   if (file_format.eq.0) then ! cas file

    READ(14,*) raw_num_nodes,raw_num_elements
    allocate(raw_nodes(raw_num_nodes,3))
    allocate(raw_elements(raw_num_elements,3))
    do inode=1,raw_num_nodes
     READ(14,*) raw_nodes(inode,1), &
          raw_nodes(inode,2),raw_nodes(inode,3)
    enddo
    do iface=1,raw_num_elements
     READ(14,*) raw_elements(iface,1),raw_elements(iface,2), &
         raw_elements(iface,3)
    enddo

   else if (file_format.eq.1) then ! vtk file

    do ivtk=1,4
     read(14,*) discard
    enddo
    read(14,'(a6)',advance='no') points_line
    read(14,*) raw_num_nodes

    allocate(raw_nodes(raw_num_nodes,3))
    do inode=1,raw_num_nodes
     READ(14,*) raw_nodes(inode,1), &
          raw_nodes(inode,2),raw_nodes(inode,3)
    enddo

    read(14,*) discard,raw_num_elements

    allocate(raw_elements(raw_num_elements,3))

    do iface=1,raw_num_elements
     READ(14,*) dummy_num_nodes_per_elem, &
         raw_elements(iface,1), &
         raw_elements(iface,2), &
         raw_elements(iface,3)
     do dir=1,3
      raw_elements(iface,dir)=raw_elements(iface,dir)+1
     enddo
    enddo

   else
    print *,"file_format invalid"
    stop
   endif

   close(14)

   FSI(part_id)%NumNodes=raw_num_nodes
   FSI(part_id)%NumIntElems=raw_num_elements
   print *,"NumNodes ",FSI(part_id)%NumNodes
   print *,"NumIntElems ",FSI(part_id)%NumIntElems

   call init_FSI(part_id,1)  ! allocate_intelem=1

   do inode=1,FSI(part_id)%NumNodes
    do dir=1,3
     xval(dir)=raw_nodes(inode,dir)
    enddo
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
    do dir=1,3
     localElem(dir)=raw_elements(iface,dir)
    enddo
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

    call SUB_ORDER_NODES(local_nodes,localElem)

    FSI(part_id)%ElemData(1,iface)=3 ! number of nodes in element
    FSI(part_id)%ElemData(2,iface)=1 ! part number
    FSI(part_id)%ElemData(DOUBLYCOMP,iface)=0 ! singly wetted
    call abort_sci_clsvof()
    do dir=1,3
     FSI(part_id)%IntElem(dir,iface)=localElem(dir)
    enddo
   enddo  ! iface, looping faces

   deallocate(raw_nodes)
   deallocate(raw_elements)

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

subroutine convert_2D_to_3D_nodes_FSI(part_id,inode,stand_alone_flag)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: stand_alone_flag
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: inode
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
  !  new nodes: z=1   x          x
  !  old nodes: z=-1  x          x
 local_nodes=FSI(part_id)%NumNodes
 orig_nodes=local_nodes/2
 if (orig_nodes*2.eq.local_nodes) then
  if ((inode.ge.1).and.(inode.le.orig_nodes)) then
   FSI(part_id)%Node_old(3,inode)=-one ! z coordinate
   do dir=1,3
    FSI(part_id)%Node_old(dir,inode+orig_nodes)= &
            FSI(part_id)%Node_old(dir,inode)
   enddo
   FSI(part_id)%Node_old(3,inode+orig_nodes)=one ! z coordinate
   do dir=1,3
    FSI(part_id)%Node_new(dir,inode)= &
      FSI(part_id)%Node_old(dir,inode)
    FSI(part_id)%Node_new(dir,inode+orig_nodes)= &
      FSI(part_id)%Node_old(dir,inode+orig_nodes)
   enddo

   if (stand_alone_flag.eq.1) then
    ! do nothing
   else if (stand_alone_flag.eq.0) then

    FSI(part_id)%NodeVel(3,inode)=zero

    do dir=1,3
     FSI(part_id)%Node_current(dir,inode)= &
             FSI(part_id)%Node_old(dir,inode)
     FSI(part_id)%Node_current(dir,inode+orig_nodes)= &
             FSI(part_id)%Node_old(dir,inode+orig_nodes)

     FSI(part_id)%Node(dir,inode)= &
             FSI(part_id)%Node_old(dir,inode)
     FSI(part_id)%Node(dir,inode+orig_nodes)= &
             FSI(part_id)%Node_old(dir,inode+orig_nodes)

     FSI(part_id)%NodeVel(dir,inode+orig_nodes)= &
        FSI(part_id)%NodeVel(dir,orig_nodes)

     FSI(part_id)%NodeVel_old(dir,inode+orig_nodes)= &
        FSI(part_id)%NodeVel(dir,orig_nodes)
     FSI(part_id)%NodeVel_new(dir,inode+orig_nodes)= &
        FSI(part_id)%NodeVel(dir,orig_nodes)
    enddo ! dir=1..3
    do dir=1,3
     FSI(part_id)%NodeForce(dir,inode+orig_nodes)= &
        FSI(part_id)%NodeForce(dir,inode)
     FSI(part_id)%NodeForce_old(dir,inode+orig_nodes)= &
        FSI(part_id)%NodeForce_old(dir,inode)
     FSI(part_id)%NodeForce_new(dir,inode+orig_nodes)= &
        FSI(part_id)%NodeForce_new(dir,inode)
    enddo

    FSI(part_id)%NodeMass(inode+orig_nodes)= &
       FSI(part_id)%NodeMass(inode)
    FSI(part_id)%NodeDensity(inode+orig_nodes)= &
       FSI(part_id)%NodeDensity(inode)

    FSI(part_id)%NodeTemp(inode+orig_nodes)= &
       FSI(part_id)%NodeTemp(inode)
    FSI(part_id)%NodeTemp_new(inode+orig_nodes)= &
       FSI(part_id)%NodeTemp_new(inode)

   else
    print *,"stand_alone_flag invalid"
    stop
   endif
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

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: iface
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

     !           3x         2y     1y
     !   1x      2x         3y
   FSI(part_id)%IntElem(1,iface+orig_elements)= &
     FSI(part_id)%IntElem(3,iface)
   FSI(part_id)%IntElem(2,iface+orig_elements)= &
     FSI(part_id)%IntElem(1,iface)+orig_nodes
   FSI(part_id)%IntElem(3,iface+orig_elements)= &
     FSI(part_id)%IntElem(1,iface)
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

INTEGER_T, INTENT(in) :: part_id
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
INTEGER_T :: stand_alone_flag

  if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
   print *,"part_id invalid"
   stop
  endif
  if (FSI(part_id)%part_id.ne.part_id) then
   print *,"FSI(part_id)%part_id.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
   print *,"ioproc invalid"
   stop
  endif

  FSI(part_id)%flag_2D_to_3D=1
  FSI(part_id)%CTML_flag=0

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

    ! gingerbread man and Xue: 1x1 box  0<x,y<1
   if (probtype.eq.400) then
    newxxblob1(1)=2.0/30.0
    newxxblob1(2)=2.0/30.0
    newxxblob1(3)=0.0
    radradblob1=30.0 
   else if (probtype.eq.404) then
    newxxblob1(1)=0.5
    newxxblob1(2)=0.5
    newxxblob1(3)=0.0
    radradblob1=4000.0 
   else if (probtype.eq.406) then
    newxxblob1(1)=0.0
    newxxblob1(2)=0.0
    newxxblob1(3)=0.0
    radradblob1=1.0 
   else
    print *,"probtype invalid"
    stop
   endif


   denpaddle=one
   dampingpaddle=zero

   if (probtype.eq.400) then ! gingerbread
    dwave="gingeroutline"
    print *,"opening ",dwave
    OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
   else if (probtype.eq.404) then ! Xue
    dwave="xueoutline"
    print *,"opening ",dwave
    OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
   else if (probtype.eq.406) then ! Fractal
    dwave="kochsnowoutline"
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

    ! allocates and inits: ElemData,IntElem,Node_old,
    !  Node_new,Node_current,NodeVel_old,NodeVel_new,NodeForce_old,
    !  NodeForce_new,NodeTemp_old,NodeTemp_new,NodeMass,NodeDensity
   call init_FSI(part_id,1)  ! allocate_intelem=1

   do dir=1,3
    maxnode(dir)=-1.0e+10
    minnode(dir)=1.0e+10
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
    enddo ! dir=1..3
    
    do dir=1,3
     FSI(part_id)%Node_old(dir,inode)=xval1(dir)
     FSI(part_id)%Node_new(dir,inode)=xval1(dir)
    enddo
   
    stand_alone_flag=1 
    call convert_2D_to_3D_nodes_FSI(part_id,inode,stand_alone_flag)

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

     !          3x
     !  1x      2x
    localElem(3)=localElem(2)+orig_nodes

    if ((localElem(1).eq.localElem(2)).or. &
        (localElem(1).eq.localElem(3)).or. &
        (localElem(2).eq.localElem(3))) then
     print *,"duplicate nodes for triangle"
     stop
    endif

    FSI(part_id)%ElemData(1,iface)=3 ! number of nodes in element
    FSI(part_id)%ElemData(2,iface)=1 ! part number
    FSI(part_id)%ElemData(DOUBLYCOMP,iface)=0 ! singly wetted
    call abort_sci_clsvof()
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

INTEGER_T, INTENT(in) :: part_id
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
  if (FSI(part_id)%part_id.ne.part_id) then
   print *,"FSI(part_id)%part_id.ne.part_id"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
   print *,"ioproc invalid"
   stop
  endif

  FSI(part_id)%flag_2D_to_3D=0
  FSI(part_id)%CTML_flag=0

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

   if ((probtype.eq.401).or. &
       (probtype.eq.415)) then ! helix or shock sphere interaction
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
    maxnode(dir)=-1.0e+10
    minnode(dir)=1.0e+10
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
    FSI(part_id)%ElemData(DOUBLYCOMP,iface)=0 ! singly wetted
    call abort_sci_clsvof()
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
    maxnode(dir)=-1.0e+10
    minnode(dir)=1.0e+10
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
     FSI(1)%ElemData(DOUBLYCOMP,iface+quad_counter)=0   ! singly wetted
     call abort_sci_clsvof()
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
     newxxblob(3)=0.5d0
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

     newxxblob(1)=0.5d0
     newxxblob(2)=0.5d0
     newxxblob(3)=0.5d0

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
      newxxblob(3)=0.5d0
      radradblob=7.5     ! 12.5/1.67=7.5

! was 10 for whale file dated December, 2007  (whalenormal)
! xval(1)=xtemp(3) xval(2)=xtemp(1) xval(3)=xtemp(2)
     else if (axis_dir.eq.2) then
      xxblob(1)=0.0
      xxblob(2)=0.0
      xxblob(3)=0.0
      newxxblob(1)=0.0
      newxxblob(2)=0.5d0
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
      newxxblob(2)=0.5d0
      radradblob=10.87   ! was 12.5, now 10.87
! whaletailup 
     else if (axis_dir.eq.0) then
      xxblob(1)=0.0 ! will become y
      xxblob(2)=-1.0 ! will become -x
      invertfactor(2)=-1.0
      xxblob(3)=0.0  ! will become z
      newxxblob(1)=0.0
      newxxblob(2)=1.0
      newxxblob(3)=0.5d0
      radradblob=10.87   ! was 12.5, now 10.87
! whaletaildown
     else if (axis_dir.eq.5) then
      xxblob(1)=0.0 ! will become y
      xxblob(2)=-1.0 ! will become -x
      invertfactor(2)=-1.0
      xxblob(3)=0.0  ! will become z
      newxxblob(1)=0.0
      newxxblob(2)=1.0
      newxxblob(3)=0.5d0
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
     xxblob(3)=0.5d0
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
      xxblob(3)=-0.5d0
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
      ! allocates NodeMass, NodeDensity + other variables.
      ! it is assumed that the number of nodes and elements does not
      ! change from frame to frame.
     call init_FSI(local_part_id,1)
    endif

    do dir=1,3
     maxnode(dir)=-1.0e+10
     minnode(dir)=1.0e+10
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
      !singly wetted, but do not call "fill" for these?
      FSI(1)%ElemData(DOUBLYCOMP,i)=2 
      call abort_sci_clsvof()

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
      print *,"k invalid geominit k=",k
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
      FSI(1)%ElemData(DOUBLYCOMP,i+quad_counter)=0   ! singly wetted
      if (probtype.eq.57) then  ! 57 heart    56 swimmer
       FSI(1)%ElemData(DOUBLYCOMP,i+quad_counter)=1   ! doubly wetted
      endif
      if ((probtype.eq.56).and.(1.eq.1)) then  ! BOXSWIMMER
       FSI(1)%ElemData(DOUBLYCOMP,i+quad_counter)=1   ! doubly wetted
      endif
      if (probtype.eq.562) then  ! 562 whale
       FSI(1)%ElemData(DOUBLYCOMP,i+quad_counter)=0 ! 0=singly 1=doubly wetted
       if (axis_dir.eq.6) then
        FSI(1)%ElemData(DOUBLYCOMP,i+quad_counter)=0 
       endif
       call abort_sci_clsvof()
      endif
      if (probtype.eq.5600) then ! dog
       FSI(1)%ElemData(DOUBLYCOMP,i+quad_counter)=0 ! 0=singly 1=doubly wetted
       call abort_sci_clsvof()
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
    maxnode(dir)=-1.0e+10
    minnode(dir)=1.0e+10
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
     print *,"k invalid viorel_sphere_geominit k=",k
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
     FSI(1)%ElemData(DOUBLYCOMP,i+quad_counter)=0 ! singly (=1 doubly wetted)
     call abort_sci_clsvof()
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
    maxnode(dir)=-1.0e+10
    minnode(dir)=1.0e+10
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
     print *,"k invalid internal_inflow_geominit k=",k
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
     FSI(1)%ElemData(DOUBLYCOMP,i+quad_counter)=0 ! singly wetted (=1 doubly)
     call abort_sci_clsvof()
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


! overall_solid_advance and post_process_nodes_elements 
! called every FSI_interval steps.
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
  radradblob=0.01d0 ! units are expected in cm

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
   maxnode(dir)=-1.0e+10
   minnode(dir)=1.0e+10
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
   FSI(1)%ElemData(DOUBLYCOMP,i)=0   ! singly wetted
   call abort_sci_clsvof()
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

  REAL_T, INTENT(in) :: cur_time
  REAL_T, INTENT(out) :: value
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
   newxxblob(3)=0.5d0
   radradblob=7.5     ! 12.5/1.67=7.5

! was 10 for whale file dated December, 2007  (whalenormal)
! xval(1)=xtemp(3) xval(2)=xtemp(1) xval(3)=xtemp(2)
  else if (whale_type.eq.2) then
   xxblob(1)=0.0
   xxblob(2)=0.0
   xxblob(3)=0.0
   newxxblob(1)=0.0
   newxxblob(2)=0.5d0
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
   newxxblob(2)=0.5d0
   radradblob=10.87   ! was 12.5, now 10.87
! whaletailup 
  else if (whale_type.eq.0) then
   xxblob(1)=0.0 ! will become y
   xxblob(2)=-1.0 ! will become -x
   invertfactor(2)=-1.0
   xxblob(3)=0.0  ! will become z
   newxxblob(1)=0.0
   newxxblob(2)=1.0
   newxxblob(3)=0.5d0
   radradblob=10.87   ! was 12.5, now 10.87
! whaletaildown
  else if (whale_type.eq.5) then
   xxblob(1)=0.0 ! will become y
   xxblob(2)=-1.0 ! will become -x
   invertfactor(2)=-1.0
   xxblob(3)=0.0  ! will become z
   newxxblob(1)=0.0
   newxxblob(2)=1.0
   newxxblob(3)=0.5d0
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
   maxnode(dir)=-1.0e+10
   minnode(dir)=1.0e+10
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
   FSI(1)%ElemData(DOUBLYCOMP,i)=0   ! singly wetted
   call abort_sci_clsvof()
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
REAL_T, INTENT(in) :: paddle_pos,paddle_vel
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
  dampingpaddle=0.01d0

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
   maxnode(dir)=-1.0e+10
   minnode(dir)=1.0e+10
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
   FSI(1)%ElemData(DOUBLYCOMP,iface)=0   ! singly wetted
   FSI(1)%ElemData(DOUBLYCOMP,iface)=1  ! doubly wetted
   if (iface.gt.FSI(1)%NumIntElemsPaddle) then
    FSI(1)%ElemData(DOUBLYCOMP,iface)=1  ! doubly wetted
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
REAL_T, INTENT(in) :: paddle_pos,paddle_vel
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
  dampingpaddle=0.01d0

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
   maxnode(dir)=-1.0e+10
   minnode(dir)=1.0e+10
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
    FSI(1)%ElemData(DOUBLYCOMP,iface)=0   ! singly wetted
    call abort_sci_clsvof()
   else if (axis_dir.eq.2) then
    FSI(1)%ElemData(DOUBLYCOMP,iface)=1   ! doubly wetted
   else
    print *,"bad axis_dir initship2 probtype,axis_dir ",probtype,axis_dir
    stop
   endif
!   FSI(1)%ElemData(DOUBLYCOMP,iface)=1   ! doubly wetted

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
193    FORMAT(E20.11,E20.11,E20.11)
194    FORMAT(I12,I12,I12,I12)

return
end subroutine initship

! overall_solid_advance is called from:
!  CLSVOF_ReadNodes
! if probtype==701,538,541 then post_process_nodes_elements is not needed.
subroutine overall_solid_advance(CLSVOF_curtime,CLSVOF_dt, &
  part_id,ioproc,isout)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: ioproc
INTEGER_T, INTENT(in) :: isout
REAL_T, INTENT(in) :: CLSVOF_curtime,CLSVOF_dt
INTEGER_T :: ifirst
REAL_T :: STEPSPERIOD,LL_CLSVOF,UU_CLSVOF,TT_CLSVOF,whale_dt

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%part_id.ne.part_id) then
  print *,"FSI(part_id)%part_id.ne.part_id"
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
 else if ((probtype.eq.400).or. &
          (probtype.eq.406).or. &
          (probtype.eq.404)) then
  call init_gingerbread2D(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout) 
 else if ((probtype.eq.401).or.(probtype.eq.415)) then
  call init_helix(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout) 
 else
   ! ifirst=0 (in overall_solid_advance)
  call init_from_cas(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,part_id,isout) 
 endif

return
end subroutine overall_solid_advance

! called from CLSVOF_ReadHeader
subroutine overall_solid_init(CLSVOFtime,ioproc,part_id,isout)
use global_utility_module

IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: ioproc
INTEGER_T, INTENT(in) :: isout
REAL_T, INTENT(in) :: CLSVOFtime
INTEGER_T :: ifirst
REAL_T :: paddle_pos,paddle_vel,CLSVOF_dt
REAL_T :: STEPSPERIOD,LL_CLSVOF,UU_CLSVOF,TT_CLSVOF,whale_dt
INTEGER_T :: ctml_part_id
INTEGER_T :: fsi_part_id

 if ((ioproc.ne.1).and.(ioproc.ne.0)) then
  print *,"ioproc invalid"
  stop
 endif
 ctml_part_id=ctml_part_id_map(part_id)
 fsi_part_id=fsi_part_id_map(part_id)

 if (((ctml_part_id.ge.1).and. &
      (ctml_part_id.le.CTML_NPARTS)).or. &
     ((fsi_part_id.ge.1).and. &
      (fsi_part_id.le.FSI_NPARTS))) then

  if ((ctml_part_id.ge.1).and. &
      (ctml_part_id.le.CTML_NPARTS)) then

   FSI(part_id)%exclusive_doubly_wetted=1
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
    FSI(part_id)%exclusive_doubly_wetted=1
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
     FSI(part_id)%normal_invert=0
     call abort_sci_clsvof()
    endif
    call geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)
   else if (probtype.eq.563) then
    call gearinit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)
   else if (probtype.eq.5600) then ! dog
    call geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)
   else if (probtype.eq.5601) then ! viorel sphere
    FSI(part_id)%normal_invert=0
    call abort_sci_clsvof()
    call viorel_sphere_geominit(sci_curtime,sci_dt,ifirst, &
      sci_sdim,sci_istop,sci_istep)
   else if (probtype.eq.5602) then ! internal inflow
    FSI(part_id)%normal_invert=1
    call abort_sci_clsvof()
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
   else if ((probtype.eq.400).or. &
            (probtype.eq.406).or. &
            (probtype.eq.404)) then
    call init_gingerbread2D(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop, &
     sci_istep,ioproc,part_id,isout) 
   else if ((probtype.eq.401).or.(probtype.eq.415)) then
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

 if ((probtype.eq.563).or. & ! gear
     (probtype.eq.50).or. &  ! paddle
     (probtype.eq.9)) then   ! ship
  ! do nothing
 else
  print *,"advance_solid only for gear, paddle or ship"
  stop
 endif

 do i=1,FSI(part_id)%NumNodes
  do dir=1,3
   FSI(part_id)%Node_old(dir,i)=FSI(part_id)%Node_new(dir,i)
   FSI(part_id)%NodeVel_old(dir,i)=FSI(part_id)%NodeVel_new(dir,i)
  enddo
  do dir=1,3
   FSI(part_id)%NodeForce_old(dir,i)=FSI(part_id)%NodeForce_new(dir,i)
  enddo
 enddo

 do i=1,FSI(part_id)%NumNodes
  do dir=1,3
   FSI(part_id)%Node(dir,i)=FSI(part_id)%Node_new(dir,i)
   FSI(part_id)%NodeVel(dir,i)=FSI(part_id)%NodeVel_new(dir,i)
  enddo
  do dir=1,3
   FSI(part_id)%NodeForce(dir,i)=FSI(part_id)%NodeForce_new(dir,i)
  enddo

   ! in: advance_solid
  FSI(part_id)%NodeMass(i)=one
  FSI(part_id)%NodeDensity(i)=one
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
   do dir=1,3
    FSI(part_id)%NodeForce_new(dir,i)=0.0
   enddo
  enddo

 endif

return
end subroutine advance_solid





! --------------------  SOLID ADVANCE STUFF ENDS HERE --------------


subroutine checkinpointBIG( &
  xclosest, &
  normal_closest, &
  inode,elemnum, &
  unsigned_mindist, &
  xc, & ! target point at which the signed distance is sought.
  inplane, &
  FSI_mesh_type, &
  part_id, &
  max_part_id, &
  time,dx)
use global_utility_module
IMPLICIT NONE

type(mesh_type), INTENT(in) :: FSI_mesh_type
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
INTEGER_T, INTENT(in) :: inode,elemnum
INTEGER_T, INTENT(inout) :: inplane
REAL_T, INTENT(in) :: time 
REAL_T, dimension(3), INTENT(inout) :: xclosest
REAL_T, dimension(3), INTENT(inout) :: normal_closest
REAL_T, dimension(3), INTENT(in) :: dx
REAL_T, dimension(3), INTENT(in) :: xc
REAL_T, INTENT(inout) :: unsigned_mindist
REAL_T :: curdist,mag
INTEGER_T :: dir
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3) :: xfoot
REAL_T, dimension(3) :: xfoot_pert
REAL_T, dimension(3) :: xtarget
REAL_T, dimension(3) :: xtarget_pert
REAL_T, dimension(3) :: ntarget
REAL_T, dimension(3) :: velparm

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI_mesh_type%part_id.ne.part_id) then
  print *,"FSI_mesh_type%part_id.ne.part_id"
  stop
 endif
 if (time.ge.zero) then
  ! do nothing
 else
  print *,"time invalid"
  stop
 endif

 nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,elemnum)

 if (nodes_per_elem.ne.3) then
  print *,"nodes_per_elem invalid"
  stop
 endif

 do dir=1,3
  if (dx(dir).gt.zero) then
   ! do nothing
  else
   print *,"dx invalid"
   stop
  endif
  xfoot(dir)= &
    FSI_mesh_type%NodeBIG(dir,FSI_mesh_type%IntElemBIG(inode,elemnum))
  xfoot_pert(dir)=xfoot(dir)+0.1*dx(1)* &
    FSI_mesh_type%NodeNormalBIG(dir,FSI_mesh_type%IntElemBIG(inode,elemnum))
  velparm(dir)=zero
 enddo
 call get_target_from_foot(xfoot,xtarget, &
   velparm,time, &
   FSI_mesh_type, &
   part_id, &
   max_part_id)
 call get_target_from_foot(xfoot_pert,xtarget_pert, &
   velparm,time, &
   FSI_mesh_type, &
   part_id, &
   max_part_id)
 mag=zero
 do dir=1,3
  ntarget(dir)=xtarget_pert(dir)-xtarget(dir)
  mag=mag+ntarget(dir)**2
 enddo
 mag=sqrt(mag)
 if (mag.gt.zero) then
  ! do nothing
 else
  print *,"mag invalid checkinpointBIG 0"
  stop
 endif
 do dir=1,3
  ntarget(dir)=ntarget(dir)/mag
 enddo

 call global_xdist(xtarget,xc,curdist)
 if (curdist.ge.zero) then
  ! do nothing
 else
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
 else if ((curdist.ge.unsigned_mindist).and. &
          (inplane.eq.1)) then
  ! do nothing
 else
  print *,"curdist or inplane invalid"
  stop
 endif

return
end subroutine checkinpointBIG


subroutine checkinlineBIG( &
  eul_over_lag_scale, &
  xclosest, &
  normal_closest, &
  inode, &
  elemnum, &
  unsigned_mindist, &
  xc, &  ! target point at which the signed distance is sought.
  inplane, &
  FSI_mesh_type, &
  part_id, &
  max_part_id, &
  time,dx)
use global_utility_module
IMPLICIT NONE

type(mesh_type), INTENT(in) :: FSI_mesh_type
REAL_T, INTENT(in) :: eul_over_lag_scale
REAL_T :: adjusted_tol
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
INTEGER_T, INTENT(in) :: inode
INTEGER_T, INTENT(in) :: elemnum
INTEGER_T, INTENT(inout) :: inplane
REAL_T, INTENT(in) :: time 
REAL_T, dimension(3), INTENT(inout) :: xclosest
REAL_T, dimension(3), INTENT(inout) :: normal_closest
REAL_T, dimension(3), INTENT(in) :: dx
REAL_T, dimension(3), INTENT(in) :: xc
INTEGER_T :: inodep1
REAL_T :: local_normal
REAL_T, dimension(2,3) :: xnode,nnode
REAL_T, INTENT(inout) :: unsigned_mindist
REAL_T :: mag
INTEGER_T :: dir
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3) :: xfoot
REAL_T, dimension(3) :: xfoot_pert
REAL_T, dimension(3) :: xtarget
REAL_T, dimension(3) :: xtarget_pert
REAL_T, dimension(3) :: ntarget
REAL_T, dimension(3) :: velparm

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI_mesh_type%part_id.ne.part_id) then
  print *,"FSI_mesh_type%part_id.ne.part_id"
  stop
 endif
 if (time.ge.zero) then
  ! do nothing
 else
  print *,"time invalid"
  stop
 endif

 if ((inode.ge.1).and.(inode.le.3)) then

  if ((inode.ge.1).and.(inode.le.2)) then
   inodep1=inode+1
  else if (inode.eq.3) then
   inodep1=1
  else
   print *,"inode invalid6: inode=",inode
   stop
  endif

  nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,elemnum)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem invalid"
   stop
  endif

  do dir=1,3
   if (dx(dir).gt.zero) then
    ! do nothing
   else
    print *,"dx invalid"
    stop
   endif

   xfoot(dir)= &
     FSI_mesh_type%NodeBIG(dir,FSI_mesh_type%IntElemBIG(inode,elemnum))

   if (1.eq.0) then
    local_normal= &
     FSI_mesh_type%NodeNormalBIG(dir,FSI_mesh_type%IntElemBIG(inode,elemnum))
   else
    local_normal= &
     FSI_mesh_type%EdgeNormalBIG(3*(inode-1)+dir,elemnum)
   endif

   xfoot_pert(dir)=xfoot(dir)+0.1d0*dx(1)*local_normal
   velparm(dir)=zero
  enddo ! dir=1..3

  if ((eul_over_lag_scale.gt.zero).and. &
      (eul_over_lag_scale.le.one)) then
   adjusted_tol=eul_over_lag_scale*element_buffer_tol
  else
   print *,"eul_over_lag_scale invalid"
   stop
  endif

  call get_target_from_foot(xfoot,xtarget, &
    velparm,time, &
    FSI_mesh_type, &
    part_id, &
    max_part_id)
  call get_target_from_foot(xfoot_pert,xtarget_pert, &
    velparm,time, &
    FSI_mesh_type, &
    part_id, &
    max_part_id)
  mag=zero
  do dir=1,3
   ntarget(dir)=xtarget_pert(dir)-xtarget(dir)
   mag=mag+ntarget(dir)**2
  enddo
  mag=sqrt(mag)
  if (mag.gt.zero) then
   ! do nothing
  else
   print *,"mag invalid checkinlineBIG 0"
   stop
  endif
  do dir=1,3
   nnode(1,dir)=ntarget(dir)/mag
   xnode(1,dir)=xtarget(dir)
  enddo

  do dir=1,3
   xfoot(dir)= &
    FSI_mesh_type%NodeBIG(dir,FSI_mesh_type%IntElemBIG(inodep1,elemnum))

   if (1.eq.0) then
    local_normal= &
     FSI_mesh_type%NodeNormalBIG(dir,FSI_mesh_type%IntElemBIG(inodep1,elemnum))
   else
    local_normal= &
     FSI_mesh_type%EdgeNormalBIG(3*(inode-1)+dir,elemnum)
   endif

   xfoot_pert(dir)=xfoot(dir)+0.1d0*dx(1)*local_normal
   velparm(dir)=zero
  enddo ! dir=1..3

  call get_target_from_foot(xfoot,xtarget, &
    velparm,time, &
    FSI_mesh_type, &
    part_id, &
    max_part_id)
  call get_target_from_foot(xfoot_pert,xtarget_pert, &
    velparm,time, &
    FSI_mesh_type, &
    part_id, &
    max_part_id)
  mag=zero
  do dir=1,3
   ntarget(dir)=xtarget_pert(dir)-xtarget(dir)
   mag=mag+ntarget(dir)**2
  enddo
  mag=sqrt(mag)
  if (mag.gt.zero) then
   ! do nothing
  else
   print *,"mag invalid checkinlineBIG 1"
   stop
  endif
  do dir=1,3
   nnode(2,dir)=ntarget(dir)/mag
   xnode(2,dir)=xtarget(dir)
  enddo

  call global_checkinline(nnode,xnode,adjusted_tol,xc, &
        inplane, &  ! INTENT(inout)
        unsigned_mindist, & ! INTENT(inout)
        xclosest, & ! INTENT(inout)
        normal_closest) ! INTENT(inout)

 else
  print *,"inode invalid in checkinlineBIG"
  stop
 endif

return
end subroutine checkinlineBIG

subroutine checkinplaneBIG( &
  eul_over_lag_scale, &
  xc, &  ! target point at which the signed distance is sought.
  xclosest, &
  xclosest_project, &
  normal, &
  normal_project, &
  elemnum, &
  inplane, &
  FSI_mesh_type, &
  part_id, &
  max_part_id, &
  time)
use global_utility_module
IMPLICIT NONE

type(mesh_type), INTENT(in) :: FSI_mesh_type
REAL_T, INTENT(in) :: eul_over_lag_scale
REAL_T :: adjusted_tol
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
INTEGER_T, INTENT(in) :: elemnum
REAL_T, INTENT(in) :: time
REAL_T, dimension(3), INTENT(in) :: xc
REAL_T, dimension(3), INTENT(in) :: xclosest
REAL_T, dimension(3), INTENT(out) :: xclosest_project
REAL_T, dimension(3), INTENT(in) :: normal
REAL_T, dimension(3), INTENT(out) :: normal_project
INTEGER_T, INTENT(out) :: inplane
INTEGER_T :: dir,i
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3) :: xfoot
REAL_T, dimension(3) :: xtarget
REAL_T, dimension(3) :: velparm
REAL_T, dimension(3,3) :: xnode ! (ipoint,dir)

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI_mesh_type%part_id.ne.part_id) then
  print *,"FSI_mesh_type%part_id.ne.part_id"
  stop
 endif
 if (time.ge.zero) then
  ! do nothing
 else
  print *,"time invalid"
  stop
 endif

 nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,elemnum)
 if (nodes_per_elem.ne.3) then
  print *,"nodes_per_elem invalid"
  stop
 endif

 do dir=1,3
  xclosest_project(dir)=xclosest(dir)
  normal_project(dir)=normal(dir)
 enddo

 do i=1,3
  do dir=1,3
   xfoot(dir)=FSI_mesh_type%NodeBIG(dir,FSI_mesh_type%IntElemBIG(i,elemnum))
   velparm(dir)=zero
  enddo
  call get_target_from_foot(xfoot,xtarget, &
     velparm,time, &
     FSI_mesh_type, &
     part_id, &
     max_part_id)
  do dir=1,3
   xnode(i,dir)=xtarget(dir)
  enddo
 enddo  ! i=1..3

 if ((eul_over_lag_scale.ge.zero).and. &
     (eul_over_lag_scale.le.one)) then
  adjusted_tol=eul_over_lag_scale*element_buffer_tol
 else
  print *,"eul_over_lag_scale invalid"
  stop
 endif

 call global_checkinplane(xnode,xclosest,adjusted_tol, &
         xclosest_project,inplane)

return
end subroutine checkinplaneBIG



! Note: meshlab is a software that can fix normal orientation.
! Note: Tecplot normals face out of object if oriented 
! counter-clockwise looking from outside:
! example: ZHI FACE:
! (0,0,1), (1,0,1), (1,1,1)  (counter clockwise looking from the top)
! v1=x2-x1=(1,0,0)  v2=x3-x2=(0,1,0)
! v1 x v2 = i   j    k = (0,0,1)  
!           1   0    0
!           0   1    0
! n=v1xv2=(0,0,1)
!
! example: YLO FACE:
! (0,0,0), (1,0,0), (1,0,1)
! v1=x2-x1=(1,0,0)  v2=x3-x2=(0,0,1)
! v1xv2= i   j  k
!        1   0  0
!        0   0  1  = (0, -1,0)
!
! for 2d problems, it is assumed that the 3rd node is equal to the 2nd
! node, except that the 3rd node extends OUT of the paper. (positive z)
! If the nodes are counter clockwise in the 2d plane, then the extruded object
! will have normals pointing out of the object:
! e.g. (1,1), (0,1) -> (1,1,0), (0,1,0), (0,1,1)
! v1=(-1,0,0)  v2=(0,0,1)
! v1xv2= i  j   k
!        -1 0   0
!        0  0   1  = (0,1,0)
!
subroutine scinormalBIG(elemnum,normal, &
     FSI_mesh_type,part_id,max_part_id, &
     time)
IMPLICIT NONE

type(mesh_type), INTENT(in) :: FSI_mesh_type
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
INTEGER_T, INTENT(in) :: elemnum
REAL_T, dimension(3), INTENT(out) :: normal
REAL_T, INTENT(in) :: time
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

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI_mesh_type%part_id.ne.part_id) then
  print *,"FSI_mesh_type%part_id.ne.part_id"
  stop
 endif
 if (time.ge.zero) then
  ! do nothing
 else
  print *,"time invalid"
  stop
 endif

 nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,elemnum)
 if (nodes_per_elem.gt.3) then
  print *,"nodes_per_elem>3 not supported"
  stop
 endif

 do dir=1,3
  nodeavg(dir)=0.0
 enddo

 do i=1,3
  do dir=1,3
   xfoot(dir)=FSI_mesh_type%NodeBIG(dir,FSI_mesh_type%IntElemBIG(i,elemnum))
   velparm(dir)=zero
  enddo
  call get_target_from_foot(xfoot,xtarget, &
      velparm,time, &
      FSI_mesh_type, &
      part_id, &
      max_part_id)

  do dir=1,3
   nodesave(i,dir)=xtarget(dir)
   nodeavg(dir)=nodeavg(dir)+nodesave(i,dir)
  enddo
 enddo ! i=1..3

 do dir=1,3
  nodeavg(dir)=nodeavg(dir)/3.0d0
 enddo

  ! Note: meshlab is a software that can fix normal orientation.
  ! Note: Tecplot normals face out of object if oriented 
  ! counter-clockwise looking from outside:
  ! example: ZHI FACE:
  ! (0,0,1), (1,0,1), (1,1,1)  (counter clockwise looking from the top)
  ! v1=x2-x1=(1,0,0)  v2=x3-x2=(0,1,0)
  ! v1 x v2 = i   j    k = (0,0,1)  
  !           1   0    0
  !           0   1    0
  ! n=v1xv2=(0,0,1)
  !
  ! example: YLO FACE:
  ! (0,0,0), (1,0,0), (1,0,1)
  ! v1=x2-x1=(1,0,0)  v2=x3-x2=(0,0,1)
  ! v1xv2= i   j  k
  !        1   0  0
  !        0   0  1  = (0, -1,0)
 do i=1,2
  do dir=1,3
   vec(i,dir)=nodesave(i+1,dir)-nodesave(i,dir)
  enddo
 enddo

 normal(1)=vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2)
 normal(2)=vec(1,3)*vec(2,1)-vec(1,1)*vec(2,3)
 normal(3)=vec(1,1)*vec(2,2)-vec(2,1)*vec(1,2)

 dist=sqrt(normal(1)**2+normal(2)**2+normal(3)**2)

 if (dist.gt.0.0d0) then
  do dir=1,3
   normal(dir)=normal(dir)/dist
   normal(dir)=min(normal(dir),1.0d0)
   normal(dir)=max(normal(dir),-1.0d0)
  enddo
 else if (dist.eq.0.0d0) then
  do dir=1,3
   normal(dir)=0.0d0
  enddo
 else
  print *,"dist is NaN"
  stop
 endif

 local_normal_invert=FSI_mesh_type%normal_invert

 if (local_normal_invert.eq.1) then
  do dir=1,3
   normal(dir)=-normal(dir)
  enddo
 else if (local_normal_invert.ne.0) then
  print *,"local_normal_invert must be 0 or 1"
  stop
 endif

return
end subroutine scinormalBIG

! Note: meshlab is a software that can fix normal orientation.
! Note: Tecplot normals face out of object if oriented 
! counter-clockwise looking from outside:
! example: ZHI FACE:
! (0,0,1), (1,0,1), (1,1,1)  (counter clockwise looking from the top)
! v1=x2-x1=(1,0,0)  v2=x3-x2=(0,1,0)
! v1 x v2 = i   j    k = (0,0,1)  
!           1   0    0
!           0   1    0
! n=v1xv2=(0,0,1)
!
! example: YLO FACE:
! (0,0,0), (1,0,0), (1,0,1)
! v1=x2-x1=(1,0,0)  v2=x3-x2=(0,0,1)
! v1xv2= i   j  k
!        1   0  0
!        0   0  1  = (0, -1,0)
!
! for 2d problems, it is assumed that the 3rd node is equal to the 2nd
! node, except that the 3rd node extends OUT of the paper. (positive z)
! If the nodes are counter clockwise in the 2d plane, then the extruded object
! will have normals pointing out of the object:
! e.g. (1,1), (0,1) -> (1,1,0), (0,1,0), (0,1,1)
! v1=(-1,0,0)  v2=(0,0,1)
! v1xv2= i  j   k
!        -1 0   0
!        0  0   1  = (0,1,0)
!
subroutine scinormal(elemnum,normal, &
    FSI_mesh_type,part_id,max_part_id, &
    time)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
type(mesh_type), INTENT(in) :: FSI_mesh_type
INTEGER_T, INTENT(in) :: elemnum
REAL_T, dimension(3), INTENT(out) :: normal
REAL_T, INTENT(in) :: time
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

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI_mesh_type%part_id.ne.part_id) then
  print *,"FSI_mesh_type%part_id.ne.part_id"
  stop
 endif
 if (time.ge.zero) then
  ! do nothing
 else
  print *,"time invalid"
  stop
 endif

 nodes_per_elem=FSI_mesh_type%ElemData(1,elemnum)
 if (nodes_per_elem.lt.3) then
  print *,"nodes_per_elem<3 not supported"
  stop
 endif

 do dir=1,3
  nodeavg(dir)=0.0
 enddo

 do i=1,3
  do dir=1,3
   xfoot(dir)=FSI_mesh_type%Node(dir,FSI_mesh_type%IntElem(i,elemnum))
   velparm(dir)=zero
  enddo

  call get_target_from_foot(xfoot,xtarget, &
      velparm,time, &
      FSI_mesh_type, &
      part_id, &
      max_part_id)

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
 if (dist.ge.1.0e-15) then
  do dir=1,3
   normal(dir)=normal(dir)/dist
  enddo
 else if ((dist.ge.0.0d0).and.(dist.le.1.0E-15)) then
  do dir=1,3
   normal(dir)=0.0d0
  enddo
 else
  print *,"dist is NaN"
  stop
 endif

 local_normal_invert=FSI_mesh_type%normal_invert

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
INTEGER_T, INTENT(in) :: elemnum
REAL_T, INTENT(out) :: area
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3,3) :: nodesave
INTEGER_T :: i,j
REAL_T :: aa,bb,cc

 if ((part_id.lt.1).or.(part_id.gt.TOTAL_NPARTS)) then
  print *,"part_id invalid"
  stop
 endif
 if (FSI(part_id)%part_id.ne.part_id) then
  print *,"FSI(part_id)%part_id.ne.part_id"
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
 area=0.5d0*sqrt(aa**2+bb**2+cc**2)

 if ((area.eq.0.0).and.(1.eq.0)) then
   print *,"area=0 x1,y1,x2,y2 ",nodesave(1,1),nodesave(1,2), &
    nodesave(2,1),nodesave(2,2)
   print *,"nodes_per_elem,elemnum ",nodes_per_elem,elemnum
   stop
 endif

return
end subroutine sciarea

! called from fort_headermsg when 
!  FSI_operation.eq.OP_FSI_LAG_STRESS and 
!  FSI_sub_operation.eq.SUB_OP_FSI_CLEAR_LAG_DATA
! isout==1 => verbose
subroutine CLSVOF_clear_lag_data(ioproc,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

INTEGER_T, INTENT(in) :: ioproc,isout
INTEGER_T :: part_id
INTEGER_T :: ctml_part_id
INTEGER_T :: fsi_part_id
INTEGER_T :: dir,inode,num_nodes

 if (TOTAL_NPARTS.ge.1) then

  if (CTML_FSI_flagF().eq.1) then 
#ifdef MVAHABFSI
   call CTML_RESET_ARRAYS() ! vel_fib=zero  force_fib=zero
#else
   print *,"define MEHDI_VAHAB_FSI in GNUmakefile"
   stop
#endif
  else if (CTML_FSI_flagF().eq.0) then
   ! do nothing
  else 
   print *,"CTML_FSI_flagF() invalid"
   stop
  endif

  do part_id=1,TOTAL_NPARTS

   ctml_part_id=ctml_part_id_map(part_id)
   fsi_part_id=fsi_part_id_map(part_id)

   if ((ctml_part_id.gt.0).or. &
       (fsi_part_id.gt.0)) then

    if ((ctml_part_id.gt.0).and. &
        (ctml_part_id.le.CTML_NPARTS)) then
     do inode=1,ctml_n_fib_nodes(ctml_part_id)
      do dir=1,AMREX_SPACEDIM
       ctml_fib_vel(ctml_part_id,inode,dir)=zero
      enddo
      do dir=1,AMREX_SPACEDIM
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
    if (num_nodes.gt.0) then
     !do nothing
    else
     print *,"num_nodes invalid"
     stop
    endif
    do inode=1,num_nodes
     do dir=1,3
      FSI(part_id)%NodeVel(dir,inode)=zero
      FSI(part_id)%NodeVel_old(dir,inode)=zero
      FSI(part_id)%NodeVel_new(dir,inode)=zero
     enddo ! dir=1,3
     FSI(part_id)%NodeMass(inode)=one
     FSI(part_id)%NodeDensity(inode)=one
     do dir=1,3
      FSI(part_id)%NodeForce(dir,inode)=zero
      FSI(part_id)%NodeForce_old(dir,inode)=zero
      FSI(part_id)%NodeForce_new(dir,inode)=zero
     enddo ! dir=1,3
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


! called from fort_headermsg when FSI_operation.eq.OP_FSI_LAG_STRESS and 
! FSI_sub_operation.eq.2
! isout==1 => verbose
! NOTE: headermsg when FSI_operation.eq.OP_FSI_LAG_STRESS and FSI_sub_operation.eq.1
! (which calls CLSVOF_Copy_To_LAG)
! is called only for those blocks associated with a given node.
! It is the job of this routine to insure that all nodes have all
! the Lagrangian information.
subroutine CLSVOF_sync_lag_data(ioproc,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif
#if (mpi_activate==1)
use mpi
#endif

IMPLICIT NONE

INTEGER_T, INTENT(in) :: ioproc,isout
INTEGER_T :: ierr1
double precision, dimension(:), allocatable :: sync_velocity
double precision, dimension(:), allocatable :: temp_velocity
double precision, dimension(:), allocatable :: sync_force
double precision, dimension(:), allocatable :: temp_force

INTEGER_T part_id
INTEGER_T ctml_part_id
INTEGER_T fsi_part_id
INTEGER_T num_nodes,sync_dim,inode,inode_fiber,dir


 if (TOTAL_NPARTS.ge.1) then

  do part_id=1,TOTAL_NPARTS

   ctml_part_id=ctml_part_id_map(part_id)
   fsi_part_id=fsi_part_id_map(part_id)

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

    allocate(sync_force(sync_dim))
    allocate(temp_force(sync_dim))

    do inode=1,num_nodes
    do dir=1,3
     sync_velocity(3*(inode-1)+dir)=FSI(part_id)%NodeVel(dir,inode)
     temp_velocity(3*(inode-1)+dir)=zero

     sync_force(3*(inode-1)+dir)=FSI(part_id)%NodeForce(dir,inode)
     temp_force(3*(inode-1)+dir)=zero
    enddo
    enddo
    ierr1=0
#if (mpi_activate==1)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr1)
     ! sync_velocity is input
     ! temp_velocity is output
    call MPI_ALLREDUCE(sync_velocity,temp_velocity,sync_dim, &
     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
     ! sync_force is input
     ! temp_force is output
    call MPI_ALLREDUCE(sync_force,temp_force,sync_dim, &
     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr1)
    do inode=1,num_nodes
    do dir=1,3
     sync_velocity(3*(inode-1)+dir)=temp_velocity(3*(inode-1)+dir)
     sync_force(3*(inode-1)+dir)=temp_force(3*(inode-1)+dir)
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

     FSI(part_id)%NodeForce(dir,inode)=sync_force(3*(inode-1)+dir)
     FSI(part_id)%NodeForce_old(dir,inode)=sync_force(3*(inode-1)+dir)
     FSI(part_id)%NodeForce_new(dir,inode)=sync_force(3*(inode-1)+dir)
    enddo
    enddo

    deallocate(sync_velocity)
    deallocate(temp_velocity)

    deallocate(sync_force)
    deallocate(temp_force)

    if ((ctml_part_id.ge.1).and. &
        (ctml_part_id.le.CTML_NPARTS)) then
     do inode_fiber=1,ctml_n_fib_active_nodes(ctml_part_id)
      inode=inode_fiber
      do dir=1,AMREX_SPACEDIM
       ctml_fib_vel(ctml_part_id,inode_fiber,dir)= &
        FSI(part_id)%NodeVel(dir,inode)
       ctml_fib_frc(ctml_part_id,inode_fiber,dir)= &
        FSI(part_id)%NodeForce(dir,inode)
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

 else
  print *,"TOTAL_NPARTS invalid"
  stop
 endif

return
end subroutine CLSVOF_sync_lag_data

subroutine CLSVOF_Init_aux_Box(FSI_operation,iter,auxcomp, &
                FSI_touch_flag,ioproc,aux_isout)
use global_utility_module

IMPLICIT NONE

INTEGER_T, INTENT(in) :: FSI_operation
INTEGER_T, INTENT(in) :: iter
INTEGER_T, INTENT(in) :: auxcomp
INTEGER_T, INTENT(in) :: ioproc
INTEGER_T, INTENT(in) :: aux_isout
INTEGER_T, INTENT(inout) :: FSI_touch_flag
INTEGER_T :: lev77_local
INTEGER_T :: tid_local
INTEGER_T :: tilenum_local
INTEGER_T :: ngrow_local
INTEGER_T :: nFSI_local
REAL_T :: time_local
REAL_T :: dt_local
INTEGER_T :: xmap3D_local(3)
REAL_T :: dx3D_local(3)
INTEGER_T :: sdim_AMR_local
INTEGER_T :: dir
INTEGER_T :: LSLO(3),LSHI(3)
INTEGER_T :: i,j,k
REAL_T :: aux_xpos(3)

REAL_T, dimension(:,:,:,:), pointer :: aux_xdata3D_ptr
REAL_T, dimension(:,:,:,:), pointer :: aux_FSIdata3D_ptr
REAL_T, dimension(:,:,:,:), pointer :: aux_masknbr3D_ptr

 aux_xdata3D_ptr=>aux_xdata3D
 aux_FSIdata3D_ptr=>aux_FSIdata3D
 aux_masknbr3D_ptr=>aux_masknbr3D

 sdim_AMR_local=AMREX_SPACEDIM
 lev77_local=-1
 tid_local=0
 tilenum_local=0
 ngrow_local=3
 nFSI_local=NCOMP_FSI
 time_local=0.0d0
 dt_local=1.0d0
 do dir=1,3
  xmap3D_local(dir)=dir
 enddo

 do dir=1,3
  LSLO(dir)=contain_aux(auxcomp)%lo3D(dir)- &
          aux_FSI(auxcomp)%bounding_box_ngrow
  LSHI(dir)=contain_aux(auxcomp)%hi3D(dir)+ &
          aux_FSI(auxcomp)%bounding_box_ngrow
  dx3D_local(dir)=contain_aux(auxcomp)%dx3D
 enddo

 if (aux_FSI(auxcomp)%LS_FROM_SUBROUTINE.eq.0) then

  call CLSVOF_InitBox( &
   iter, &
   sdim_AMR_local, &
   lev77_local, &
   tid_local, &
   tilenum_local, &
   auxcomp, &
   fort_num_local_aux_grids, &
   auxcomp, &
   ngrow_local, &
   nFSI_local, &
   FSI_operation, &
   FSI_touch_flag, &
   time_local, &
   dt_local, &
   contain_aux(auxcomp)%xlo3D, &
   contain_aux(auxcomp)%xhi3D, &
   xmap3D_local, &
   dx3D_local, &
   contain_aux(auxcomp)%xlo3D, &
   contain_aux(auxcomp)%xhi3D, &
   contain_aux(auxcomp)%lo3D, &
   contain_aux(auxcomp)%hi3D, &
   LSLO,LSHI, &
   LSLO,LSHI, &
   aux_xdata3D_ptr, &
   aux_FSIdata3D_ptr, &
   aux_masknbr3D_ptr, &
   ioproc, &
   aux_isout)

  if (FSI_operation.eq.OP_FSI_MAKE_DISTANCE) then
   FSI_touch_flag=1
  else if (FSI_operation.eq.OP_FSI_MAKE_SIGN) then
   ! do nothing
  else
   print *,"FSI_operation invalid"
   stop
  endif

 else if (aux_FSI(auxcomp)%LS_FROM_SUBROUTINE.eq.1) then

  FSI_touch_flag=0
  do i=LSLO(1),LSHI(1)
  do j=LSLO(2),LSHI(2)
  do k=LSLO(3),LSHI(3)
   do dir=1,3
    aux_xpos(dir)=aux_xdata3D(i,j,k,dir)
   enddo
   call SUB_AUX_DATA(auxcomp,aux_xpos,aux_FSIdata3D(i,j,k,FSI_LEVELSET+1))
  enddo
  enddo
  enddo
 else
  print *,"aux_FSI(auxcomp)%LS_FROM_SUBROUTINE invalid"
  stop
 endif

return
end subroutine CLSVOF_Init_aux_Box

subroutine CLSVOF_Read_aux_Header(auxcomp,ioproc,aux_isout)
use global_utility_module

IMPLICIT NONE

INTEGER_T, INTENT(in) :: auxcomp
INTEGER_T, INTENT(in) :: ioproc
INTEGER_T, INTENT(in) :: aux_isout

INTEGER_T :: dir
INTEGER_T :: inode
INTEGER_T :: iface
INTEGER_T :: ifirst
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
REAL_T, dimension(3) :: xxblob1,newxxblob1
REAL_T, dimension(3) :: side_len
REAL_T :: max_side_len
INTEGER_T :: dir_max_side
REAL_T :: radradblob1
INTEGER_T localElem(3)
INTEGER_T aux_ncells
REAL_T :: local_midpoint
REAL_T :: local_sidelen
INTEGER_T :: local_ncells
INTEGER_T :: LSLO(3)
INTEGER_T :: LSHI(3)
INTEGER_T :: initflag
INTEGER_T :: file_format
INTEGER_T :: i,j,k
INTEGER_T, dimension(3) :: idx
character(80) :: discard
character(80) :: points_line
INTEGER_T :: ivtk,dummy_num_nodes_per_elem
INTEGER_T :: raw_num_nodes
INTEGER_T :: raw_num_elements
REAL_T, allocatable :: raw_nodes(:,:)
INTEGER_T, allocatable :: raw_elements(:,:)
INTEGER_T, PARAMETER :: aux_unit_id=14

 initflag=1

 if (aux_data_allocated.eq.0) then

  if (auxcomp.eq.1) then
   allocate(aux_FSI(fort_num_local_aux_grids))
  else if ((auxcomp.gt.1).and. &
           (auxcomp.le.fort_num_local_aux_grids)) then
   ! do nothing
  else
   print *,"auxcomp invalid"
   stop
  endif

  xxblob1(1)=0.0d0
  xxblob1(2)=0.0d0
  xxblob1(3)=0.0d0
  newxxblob1(1)=0.0d0
  newxxblob1(2)=0.0d0
  newxxblob1(3)=0.0d0
  radradblob1=1.0d0

  do dir=1,3
   maxnode(dir)=-1.0e+10
   minnode(dir)=1.0e+10
   maxnodebefore(dir)=-1.0d+10
   minnodebefore(dir)=1.0d+10
  enddo

  call SUB_BOUNDING_BOX_AUX(auxcomp,minnode,maxnode, &
          aux_FSI(auxcomp)%LS_FROM_SUBROUTINE, &
          contain_aux(auxcomp)%aux_ncells_max_side)

  aux_ncells=contain_aux(auxcomp)%aux_ncells_max_side

  aux_FSI(auxcomp)%IntElemDim=3
  aux_FSI(auxcomp)%part_id=auxcomp
  aux_FSI(auxcomp)%flag_2D_to_3D=0

   ! The normals in Cody's vtk files point out of the object (as seen
   ! from Tecplot: 1. Analyze 2. Calculate Variables 3. Grid K unit normal
   ! vector.  4. Calculate 5. Close 6. Vector
   ! Therefore, since the tecplot normals point out of the object,
   ! the node ordering in elements is counter-clockwise as viewed from 
   ! the outside.
   ! one will have the following:
   ! normal_invert=0 => sign of distance function is positive in the object
   ! normal_invert=1 => sign of distance function is negative in the object
  aux_FSI(auxcomp)%normal_invert=0
  aux_FSI(auxcomp)%exclusive_doubly_wetted=0
  aux_FSI(auxcomp)%deforming_part=0
  aux_FSI(auxcomp)%CTML_flag=0
    !refine_factor=1 => refine the Lagrangian mesh as necessary.
    !refine_factor=0 => n_lag_levels=2
  aux_FSI(auxcomp)%refine_factor=1
  aux_FSI(auxcomp)%bounding_box_ngrow=3

  if (aux_FSI(auxcomp)%LS_FROM_SUBROUTINE.eq.0) then

   print *,"calling SUB_OPEN_AUXFILE; auxcomp,aux_unit_id = ", &
     auxcomp,aux_unit_id

   call SUB_OPEN_AUXFILE(auxcomp,aux_unit_id,file_format)

   if (file_format.eq.0) then ! cas file

    READ(aux_unit_id,*) raw_num_nodes,raw_num_elements
    allocate(raw_nodes(raw_num_nodes,3))
    allocate(raw_elements(raw_num_elements,3))
    do inode=1,raw_num_nodes
     READ(aux_unit_id,*) raw_nodes(inode,1), &
          raw_nodes(inode,2),raw_nodes(inode,3)
    enddo
    do iface=1,raw_num_elements
     READ(aux_unit_id,*) raw_elements(iface,1),raw_elements(iface,2), &
         raw_elements(iface,3)
    enddo

   else if (file_format.eq.1) then ! vtk file

    do ivtk=1,4
     read(aux_unit_id,*) discard
    enddo
    read(aux_unit_id,'(a6)',advance='no') points_line
    read(aux_unit_id,*) raw_num_nodes

    allocate(raw_nodes(raw_num_nodes,3))
    do inode=1,raw_num_nodes
     READ(aux_unit_id,*) raw_nodes(inode,1), &
          raw_nodes(inode,2),raw_nodes(inode,3)
    enddo

    read(aux_unit_id,*) discard,raw_num_elements

    allocate(raw_elements(raw_num_elements,3))

    do iface=1,raw_num_elements
     READ(aux_unit_id,*) dummy_num_nodes_per_elem, &
         raw_elements(iface,1), &
         raw_elements(iface,2), &
         raw_elements(iface,3)
     do dir=1,3
      raw_elements(iface,dir)=raw_elements(iface,dir)+1
     enddo
    enddo

   else
    print *,"file_format invalid"
    stop
   endif

   close(aux_unit_id)

   aux_FSI(auxcomp)%NumNodes=raw_num_nodes
   aux_FSI(auxcomp)%NumIntElems=raw_num_elements

   if (ioproc.eq.1) then
    print *,"auxcomp,NumNodes ",auxcomp,aux_FSI(auxcomp)%NumNodes
    print *,"auxcomp,NumIntElems ",auxcomp,aux_FSI(auxcomp)%NumIntElems
   else if (ioproc.eq.0) then
    ! do nothing
   else
    print *,"ioproc invalid"
    stop
   endif

   call init_FSI_mesh_type(aux_FSI(auxcomp),1)  ! allocate_intelem=1

   do inode=1,aux_FSI(auxcomp)%NumNodes
    do dir=1,3
     xval(dir)=raw_nodes(inode,dir)
    enddo
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
     aux_FSI(auxcomp)%Node_old(dir,inode)=xval1(dir)
     aux_FSI(auxcomp)%Node_new(dir,inode)=xval1(dir)
    enddo
     
   enddo  ! inode=1,NumNodes
  
   do iface=1,aux_FSI(auxcomp)%NumIntElems
    do dir=1,3
     localElem(dir)=raw_elements(iface,dir)
    enddo
    do inode=1,3
     if ((localElem(inode).lt.1).or. &
         (localElem(inode).gt.aux_FSI(auxcomp)%NumNodes)) then
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

    aux_FSI(auxcomp)%ElemData(1,iface)=3 ! number of nodes in element
    aux_FSI(auxcomp)%ElemData(2,iface)=1 ! part number
    aux_FSI(auxcomp)%ElemData(DOUBLYCOMP,iface)=0 ! singly wetted
    do dir=1,3
     aux_FSI(auxcomp)%IntElem(dir,iface)=localElem(dir)
    enddo
   enddo  ! iface, looping faces

   deallocate(raw_nodes)
   deallocate(raw_elements)

   if (ioproc.eq.1) then
    do dir=1,3 
     print *,"(before)auxcomp,dir,min,max ",auxcomp, &
             dir,minnodebefore(dir),maxnodebefore(dir)
    enddo
    do dir=1,3 
     print *,"(after)auxcomp,dir,min,max ",auxcomp, &
             dir,minnode(dir),maxnode(dir)
    enddo
   else if (ioproc.eq.0) then
    ! do nothing
   else
    print *,"ioproc invalid"
    stop
   endif

   call init2_FSI_mesh_type(aux_FSI(auxcomp)) 

   ifirst=1
   call init3_FSI_mesh_type(aux_FSI(auxcomp),ifirst) 

  else if (aux_FSI(auxcomp)%LS_FROM_SUBROUTINE.eq.1) then
   ! do nothing
  else
   print *,"aux_FSI(auxcomp)%LS_FROM_SUBROUTINE invalid"
   stop
  endif 

  max_side_len=0.0d0
  dir_max_side=0
  do dir=1,3
   side_len(dir)=maxnode(dir)-minnode(dir)
   if (side_len(dir).gt.0.0d0) then
    if (side_len(dir).gt.max_side_len) then
     dir_max_side=dir
     max_side_len=side_len(dir)
    endif
   else
    print *,"side_len invalid"
    stop
   endif
  enddo ! dir=1..3

  if (max_side_len.gt.0.0d0) then
   contain_aux(auxcomp)%lo3D(dir_max_side)=0
   contain_aux(auxcomp)%hi3D(dir_max_side)=aux_ncells-1
   contain_aux(auxcomp)%xlo3D(dir_max_side)= &
           minnode(dir_max_side)-0.2d0*max_side_len
   contain_aux(auxcomp)%xhi3D(dir_max_side)= &
           maxnode(dir_max_side)+0.2d0*max_side_len
   contain_aux(auxcomp)%dx3D=(contain_aux(auxcomp)%xhi3D(dir_max_side)- &
                     contain_aux(auxcomp)%xlo3D(dir_max_side))/aux_ncells
   do dir=1,3
    if (dir.ne.dir_max_side) then
     local_ncells=NINT(side_len(dir)*aux_ncells/max_side_len)
     if (local_ncells.le.aux_ncells) then
      if (local_ncells.lt.aux_ncells) then
       local_ncells=local_ncells+1
      endif
      if (local_ncells.ge.1) then
       if (local_ncells.lt.4) then
        local_ncells=4
       endif
       contain_aux(auxcomp)%lo3D(dir)=0
       contain_aux(auxcomp)%hi3D(dir)=local_ncells-1
       local_sidelen=local_ncells*contain_aux(auxcomp)%dx3D
       local_midpoint=0.5d0*(minnode(dir)+maxnode(dir))
       contain_aux(auxcomp)%xlo3D(dir)=local_midpoint-0.5d0*local_sidelen 
       contain_aux(auxcomp)%xhi3D(dir)=local_midpoint+0.5d0*local_sidelen 
      else
       print *,"local_ncells invalid"
       stop
      endif
     else
      print *,"local_ncells invalid"
      stop
     endif
    else if (dir.eq.dir_max_side) then
     ! do nothing
    else
     print *,"dir invalid"
     stop
    endif
   enddo ! dir=1..3
   do dir=1,3
    LSLO(dir)=contain_aux(auxcomp)%lo3D(dir)- &
            aux_FSI(auxcomp)%bounding_box_ngrow
    LSHI(dir)=contain_aux(auxcomp)%hi3D(dir)+ &
            aux_FSI(auxcomp)%bounding_box_ngrow
   enddo
   allocate(contain_aux(auxcomp)%LS3D(LSLO(1):LSHI(1), &
           LSLO(2):LSHI(2),LSLO(3):LSHI(3),1))
   allocate(aux_xdata3D(LSLO(1):LSHI(1), &
           LSLO(2):LSHI(2),LSLO(3):LSHI(3),3))
   allocate(aux_FSIdata3D(LSLO(1):LSHI(1), &
           LSLO(2):LSHI(2),LSLO(3):LSHI(3),NCOMP_FSI))
   allocate(aux_masknbr3D(LSLO(1):LSHI(1), &
           LSLO(2):LSHI(2),LSLO(3):LSHI(3),2))
  else
   print *,"max_side_len invalid"
   stop
  endif

  do i=LSLO(1),LSHI(1)
  do j=LSLO(2),LSHI(2)
  do k=LSLO(3),LSHI(3)
   idx(1)=i
   idx(2)=j
   idx(3)=k
   do dir=1,3
    aux_xdata3D(i,j,k,dir)=contain_aux(auxcomp)%xlo3D(dir)+ &
          (idx(dir)+0.5d0)*contain_aux(auxcomp)%dx3D
   enddo
   ! mask2==1 if (i,j,k) interior to tile
   ! mask1==0 if (i,j,k) is exterior to tile and is a
   !  (coarse/fine ghost cell or EXT_DIR ghost cell).
   aux_masknbr3D(i,j,k,1)=0
   aux_masknbr3D(i,j,k,2)=0

   do dir=1,3
    aux_FSIdata3D(i,j,k,FSI_VELOCITY+dir)=0.0d0
   enddo
   do dir=1,NCOMP_FORCE_STRESS
    aux_FSIdata3D(i,j,k,FSI_FORCE+dir)=0.0d0
   enddo
   aux_FSIdata3D(i,j,k,FSI_LEVELSET+1)=-99999.0
   ! bit 0=1 if +sign hits
   ! bit 1=1 if -sign hits
   aux_FSIdata3D(i,j,k,FSI_SIGN_CONFLICT+1)=0.0d0
   aux_FSIdata3D(i,j,k,FSI_TEMPERATURE+1)=0.0d0
   aux_FSIdata3D(i,j,k,FSI_EXTRAP_FLAG+1)=FSI_NOTHING_VALID
   aux_FSIdata3D(i,j,k,FSI_AREA_PER_VOL+1)=0.0d0
  enddo
  enddo
  enddo
  do i=contain_aux(auxcomp)%lo3D(1),contain_aux(auxcomp)%hi3D(1)
  do j=contain_aux(auxcomp)%lo3D(2),contain_aux(auxcomp)%hi3D(2)
  do k=contain_aux(auxcomp)%lo3D(3),contain_aux(auxcomp)%hi3D(3)
   aux_masknbr3D(i,j,k,1)=1
   aux_masknbr3D(i,j,k,2)=1
  enddo
  enddo
  enddo

  if (ioproc.eq.1) then
   print *,"auxcomp,dir_max_side ",auxcomp,dir_max_side
   print *,"auxcomp,dx3D ",auxcomp,contain_aux(auxcomp)%dx3D
   do dir=1,3
    print *,"auxcomp,dir,lo3D,hi3D ", &
     auxcomp,dir, &
     contain_aux(auxcomp)%lo3D(dir), &
     contain_aux(auxcomp)%hi3D(dir)
    print *,"auxcomp,dir,xlo3D,xhi3D ", &
     auxcomp,dir, &
     contain_aux(auxcomp)%xlo3D(dir), &
     contain_aux(auxcomp)%xhi3D(dir)
    print *,"auxcomp,dir,LSLO,LSHI ",auxcomp,dir,LSLO(dir),LSHI(dir)
   enddo
  else if (ioproc.eq.0) then
   ! do nothing
  else
   print *,"ioproc invalid"
   stop
  endif

  if (aux_FSI(auxcomp)%LS_FROM_SUBROUTINE.eq.0) then

    !EdgeNormal(dir,inode) and EdgeNormalBIG initialized here.
    !NodeNormal(dir,inode) and NodeNormalBIG initialized here.
   call post_process_nodes_elements(initflag, &
          contain_aux(auxcomp)%xlo3D, &
          contain_aux(auxcomp)%xhi3D, &
          aux_FSI(auxcomp),auxcomp,fort_num_local_aux_grids, &
          ioproc,aux_isout, &
          contain_aux(auxcomp)%dx3D)
  else if (aux_FSI(auxcomp)%LS_FROM_SUBROUTINE.eq.1) then
   ! do nothing
  else
   print *,"aux_FSI(auxcomp)%LS_FROM_SUBROUTINE invalid"
   stop
  endif 
 else 
  print *,"aux_data_allocated invalid"
  stop
 endif

return
end subroutine CLSVOF_Read_aux_Header

! called from fort_headermsg when FSI_operation.eq.OP_FSI_INITIALIZE_NODES
! fort_headermsg is called from NavierStokes::ns_header_msg_level
! ns_header_msg_level is called from:
!   NavierStokes::FSI_make_distance 
!     (  FSI_operation.eq.OP_FSI_MAKE_DISTANCE(SIGN)  )
!   NavierStokes::post_restart 
!     (  FSI_operation.eq.OP_FSI_INITIALIZE_NODES  )
!   NavierStokes::initData 
!     (  FSI_operation.eq.OP_FSI_INITIALIZE_NODES  )
!   NavierStokes::nonlinear_advection 
!     (  FSI_operation.eq.OP_FSI_LAG_STRESS  )
!     (  FSI_operation.eq.OP_FSI_UPDATE_NODES )
! isout==1 => verbose
subroutine CLSVOF_ReadHeader( &
  im_critical, &
  max_num_nodes_list, &
  max_num_elements_list, &
  num_nodes_list, &
  num_elements_list, &
  FSI_input_max_num_nodes, &
  FSI_input_max_num_elements, &
  FSI_input_num_nodes, &
  FSI_input_num_elements, &
  FSI_input_node_list, &
  FSI_input_element_list, &
  FSI_input_displacement_list, &
  FSI_input_velocity_halftime_list, &
  FSI_input_velocity_list, &
  FSI_input_force_list, &
  FSI_input_mass_list, &
  FSI_input_temperature_list, &
  FSI_output_max_num_nodes, &
  FSI_output_max_num_elements, &
  FSI_output_num_nodes, &
  FSI_output_num_elements, &
  FSI_output_node_list, &
  FSI_output_element_list, &
  FSI_output_displacement_list, &
  FSI_output_velocity_halftime_list, &
  FSI_output_velocity_list, &
  FSI_output_force_list, &
  FSI_output_mass_list, &
  FSI_output_temperature_list, &
  FSI_refine_factor, &
  FSI_bounding_box_ngrow, &
  nparts_in, &
  im_solid_map_in, &  ! im_part=im_solid_map_in(part_id)+1
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

INTEGER_T, INTENT(in) :: im_critical
INTEGER_T, INTENT(inout) :: max_num_nodes_list(num_materials)
INTEGER_T, INTENT(inout) :: max_num_elements_list(num_materials)
INTEGER_T, INTENT(inout) :: num_nodes_list(num_materials)
INTEGER_T, INTENT(inout) :: num_elements_list(num_materials)
INTEGER_T, INTENT(in) :: FSI_input_max_num_nodes
INTEGER_T, INTENT(in) :: FSI_input_max_num_elements
INTEGER_T, INTENT(in) :: FSI_input_num_nodes
INTEGER_T, INTENT(in) :: FSI_input_num_elements
REAL_T, INTENT(inout) :: FSI_input_node_list(3*FSI_input_max_num_nodes)
INTEGER_T, INTENT(inout) :: FSI_input_element_list(4*FSI_input_max_num_elements)
REAL_T, INTENT(inout) :: FSI_input_displacement_list(3*FSI_input_max_num_nodes)
REAL_T, INTENT(inout) ::  &
        FSI_input_velocity_halftime_list(3*FSI_input_max_num_nodes)
REAL_T, INTENT(inout) ::  &
        FSI_input_velocity_list(3*FSI_input_max_num_nodes)
REAL_T, INTENT(inout) ::  &
        FSI_input_force_list(NCOMP_FORCE_STRESS*FSI_input_max_num_nodes)
REAL_T, INTENT(inout) :: FSI_input_mass_list(FSI_input_max_num_nodes)
REAL_T, INTENT(inout) :: FSI_input_temperature_list(FSI_input_max_num_nodes)

INTEGER_T, INTENT(inout) :: FSI_output_max_num_nodes
INTEGER_T, INTENT(inout) :: FSI_output_max_num_elements
INTEGER_T, INTENT(inout) :: FSI_output_num_nodes
INTEGER_T, INTENT(inout) :: FSI_output_num_elements
REAL_T, INTENT(inout) :: FSI_output_node_list(3*FSI_output_max_num_nodes)
INTEGER_T, INTENT(inout) :: &
        FSI_output_element_list(4*FSI_output_max_num_elements)
REAL_T, INTENT(inout) :: &
        FSI_output_displacement_list(3*FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: &
        FSI_output_velocity_halftime_list(3*FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: &
        FSI_output_velocity_list(3*FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: &
        FSI_output_force_list(NCOMP_FORCE_STRESS*FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: FSI_output_mass_list(FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: FSI_output_temperature_list(FSI_output_max_num_nodes)

INTEGER_T, INTENT(in) :: nparts_in
INTEGER_T, INTENT(in) :: im_solid_map_in(nparts_in)
INTEGER_T :: initflag
INTEGER_T, INTENT(in) :: ioproc,isout
INTEGER_T, INTENT(in) :: CTML_FSI_INIT
REAL_T, INTENT(in) :: dx_max_level(AMREX_SPACEDIM)
REAL_T, INTENT(in) :: h_small
REAL_T, INTENT(in) :: CLSVOFtime
REAL_T, INTENT(in) :: problo(3),probhi(3)
INTEGER_T :: test_NPARTS
INTEGER_T :: part_id
INTEGER_T :: dir
INTEGER_T :: im_part
INTEGER_T :: ctml_part_id,fsi_part_id
INTEGER_T FSI_refine_factor(num_materials)
INTEGER_T FSI_bounding_box_ngrow(num_materials)
INTEGER_T im_sanity_check
INTEGER_T idir,ielem,inode

  do im_sanity_check=1,num_materials
   if ((FSI_refine_factor(im_sanity_check).lt.0).or. &
       (FSI_refine_factor(im_sanity_check).gt.100)) then
    print *,"FSI_refine_factor(im_sanity_check) invalid"
    stop
   endif
   if (FSI_bounding_box_ngrow(im_sanity_check).ne.3) then
    print *,"FSI_bounding_box_ngrow(im_sanity_check) invalid"
    stop
   endif
  enddo ! do im_sanity_check=1,num_materials

  if ((nparts_in.lt.1).or.(nparts_in.ge.num_materials)) then
   print *,"nparts_in invalid"
   stop
  endif

  if (im_critical.eq.0) then

   TOTAL_NPARTS=nparts_in
   ctml_part_id=0
   fsi_part_id=0

   do part_id=1,TOTAL_NPARTS

    im_solid_mapF(part_id)=im_solid_map_in(part_id)
    im_part=im_solid_mapF(part_id)+1

    if (CTML_FSI_mat(im_part).eq.1) then 
     ctml_part_id=ctml_part_id+1
     ctml_part_id_map(part_id)=ctml_part_id
    else if (CTML_FSI_mat(im_part).eq.0) then
     ctml_part_id_map(part_id)=0
    else
     print *,"CTML_FSI_mat(im_part) invalid"
     stop
    endif

    if ((FSI_flag(im_part).eq.FSI_PRESCRIBED_NODES).or. & 
        (FSI_flag(im_part).eq.FSI_ICE_NODES_INIT).or. & 
        (FSI_flag(im_part).eq.FSI_FLUID_NODES_INIT)) then 
     fsi_part_id=fsi_part_id+1
     fsi_part_id_map(part_id)=fsi_part_id
    else if ((FSI_flag(im_part).eq.FSI_PRESCRIBED_PROBF90).or. & 
             (FSI_flag(im_part).eq.FSI_SHOELE_CTML)) then 
     fsi_part_id_map(part_id)=0
    else
     print *,"FSI_flag(im_part) invalid in CLSVOF_ReadHeader"
     print *,"part_id,im_part,FSI_flag(im_part) ", &
             part_id,im_part,FSI_flag(im_part)
     stop
    endif

    FSI(part_id)%flag_2D_to_3D=0
    FSI(part_id)%normal_invert=0
    FSI(part_id)%exclusive_doubly_wetted=0

   enddo ! part_id=1,TOTAL_NPARTS

   CTML_NPARTS=ctml_part_id
   FSI_NPARTS=fsi_part_id
   if (FSI_NPARTS+CTML_NPARTS.gt.TOTAL_NPARTS) then
    print *,"FSI_NPARTS+CTML_NPARTS.gt.TOTAL_NPARTS"
    stop
   endif

   if (h_small.gt.zero) then
    ! do nothing
   else
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
    if (problen_ref(dir).gt.zero) then
     ! do nothing
    else
     print *,"problen_ref(dir).le.zero"
     stop
    endif
    if (dx_max_level(dir).gt.zero) then
     ! do nothing
    else
     print *,"dx_max_level(dir).le.zero"
     stop
    endif
   enddo ! dir=1..AMREX_SPACEDIM

   use_temp=0

   if (CTML_FSI_flagF().eq.1) then 
#ifdef MVAHABFSI
    if (CTML_FSI_INIT.eq.0) then
       ! CTML_INIT_SOLID is declared in CTMLFSI.F90
       ! vel_fib and force_fib are initialized to 0.0d0
     call CTML_INIT_SOLID( &
      dx_max_level, &
      problo_ref,probhi_ref, &
      ioproc, &
      ctml_n_fib_bodies, &
      ctml_max_n_fib_nodes)
     if (ctml_n_fib_bodies.eq.CTML_NPARTS) then
      allocate(ctml_n_fib_nodes(ctml_n_fib_bodies))
      allocate(ctml_n_fib_active_nodes(ctml_n_fib_bodies))
       ! CTML_GET_FIB_NODE_COUNT is declared in CTMLFSI.F90
      call CTML_GET_FIB_NODE_COUNT( &
       ctml_n_fib_bodies, & !intent(in)
       ctml_n_fib_nodes) ! intent(out)

      allocate(ctml_fib_pst(ctml_n_fib_bodies,ctml_max_n_fib_nodes, &
             AMREX_SPACEDIM))
      allocate(ctml_fib_vel(ctml_n_fib_bodies,ctml_max_n_fib_nodes, &
             AMREX_SPACEDIM))
      allocate(ctml_fib_frc(ctml_n_fib_bodies,ctml_max_n_fib_nodes, &
             AMREX_SPACEDIM))
      allocate(ctml_fib_mass(ctml_n_fib_bodies,ctml_max_n_fib_nodes))

      allocate(ctml_fib_pst_prev(ctml_n_fib_bodies,ctml_max_n_fib_nodes, &
             AMREX_SPACEDIM))
      allocate(ctml_fib_vel_halftime_prev( &
             ctml_n_fib_bodies,ctml_max_n_fib_nodes, &
             AMREX_SPACEDIM))
      allocate(ctml_fib_vel_prev( &
             ctml_n_fib_bodies,ctml_max_n_fib_nodes, &
             AMREX_SPACEDIM))
      allocate(ctml_fib_frc_prev(ctml_n_fib_bodies,ctml_max_n_fib_nodes, &
             AMREX_SPACEDIM))
      allocate(ctml_fib_mass_prev(ctml_n_fib_bodies,ctml_max_n_fib_nodes))

     else
      print *,"ctml_n_fib_bodies.eq.CTML_NPARTS failed: ", &
         ctml_n_fib_bodies,CTML_NPARTS
      stop
     endif
    else if (CTML_FSI_INIT.eq.1) then
     ! do nothing
    else
     print *,"CTML_FSI_INIT invalid"
     stop
    endif

    ctml_part_id=0

    do part_id=1,TOTAL_NPARTS

     im_part=im_solid_mapF(part_id)+1
     if (CTML_FSI_mat(im_part).eq.1) then 

      ctml_part_id=ctml_part_id+1
      if (ctml_part_id.eq.ctml_part_id_map(part_id)) then
       !do nothing
      else
       print *,"ctml_part_id invalid"
       stop
      endif

      FSI(part_id)%deforming_part=1
      FSI(part_id)%CTML_flag=1

      num_nodes_list(im_part)=ctml_n_fib_nodes(ctml_part_id)
      num_elements_list(im_part)=ctml_n_fib_nodes(ctml_part_id)-1

      if (max_num_nodes_list(im_part).eq.ctml_max_n_fib_nodes) then
       ! do nothing
      else
       print *,"max_num_nodes_list(im_part).eq.ctml_max_n_fib_nodes failed: "
       print *,"im_part=",im_part
       print *,"max_num_nodes_list(im_part)=",max_num_nodes_list(im_part)
       print *,"ctml_max_n_fib_nodes=",ctml_max_n_fib_nodes
       stop
      endif
      if (max_num_elements_list(im_part).eq.ctml_max_n_fib_nodes-1) then
       ! do nothing
      else
       print *,"max_num_elements_list(im_part).eq.ctml_max_n_fib_nodes-1 fail"
       print *,"im_part=",im_part
       print *,"max_num_elements_list(im_part)=",max_num_elements_list(im_part)
       print *,"ctml_max_n_fib_nodes=",ctml_max_n_fib_nodes
       stop
      endif
     else if (CTML_FSI_mat(im_part).eq.0) then
      FSI(part_id)%deforming_part=0
      FSI(part_id)%CTML_flag=0
     else
      print *,"CTML_FSI_mat(im_part) invalid"
      stop
     endif
    enddo ! part_id=1,TOTAL_NPARTS
   
#else
    print *,"define MVAHABFSI"
    stop
#endif

   else if (CTML_FSI_flagF().eq.0) then

    part_id=1

    if((probtype.eq.538).or.(probtype.eq.541)) then ! needle, housing
     test_NPARTS=2
     FSI(1)%deforming_part=0
     FSI(2)%deforming_part=0
     FSI(1)%CTML_flag=0
     FSI(2)%CTML_flag=0
    else if (probtype.eq.701) then  ! flapping wing - ReadHeader
     if ((axis_dir.eq.0).or.(axis_dir.eq.2)) then
      test_NPARTS=1
      FSI(1)%deforming_part=0
      FSI(1)%CTML_flag=0
     else if (axis_dir.eq.1) then
      test_NPARTS=2
      FSI(1)%deforming_part=0
      FSI(2)%deforming_part=0
      FSI(1)%CTML_flag=0
      FSI(2)%CTML_flag=0
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
     FSI(part_id)%CTML_flag=0
    else if (probtype.eq.9) then ! ship
     test_NPARTS=1
     FSI(part_id)%deforming_part=0
     FSI(part_id)%CTML_flag=0
    else
     test_NPARTS=TOTAL_NPARTS
     do part_id=1,TOTAL_NPARTS
      FSI(part_id)%deforming_part=0
      FSI(part_id)%CTML_flag=0
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

     if ((fsi_part_id_map(part_id).gt.0).or. &
         (ctml_part_id_map(part_id).gt.0)) then
  
      FSI(part_id)%part_id=part_id
      FSI(part_id)%LS_FROM_SUBROUTINE=0

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
      ! if CTML materials exists, then,
      !   overall_solid_init calls CTML_init_sci 
      !   CTML_init_sci calls CTML_GET_POS_VEL_FORCE_WT 
      !   (declared in CTMLFSI.F90)
      ! ctml_fib_pst =coord_fib
      ! ctml_fib_vel =vel_fib
      ! ctml_fib_frc =force_fib
      ! ctml_fib_mass=ds_fib
      call overall_solid_init(CLSVOFtime,ioproc,part_id,isout)  

      ! ReadHeader
      initflag=1
       !EdgeNormal(dir,inode) and EdgeNormalBIG initialized here.
       !NodeNormal(dir,inode) and NodeNormalBIG initialized here.
      call post_process_nodes_elements(initflag,problo,probhi, &
       FSI(part_id),part_id,TOTAL_NPARTS, &
       ioproc,isout,h_small)

     else if ((fsi_part_id_map(part_id).eq.0).and. &
              (ctml_part_id_map(part_id).eq.0)) then
      ! do nothing
     else
      print *,"fsi_part_id_map or ctml_part_id_map invalid"
      stop
     endif

    enddo ! part_id=1..TOTAL_NPARTS

   else
    print *,"TOTAL_NPARTS invalid: ",TOTAL_NPARTS
    stop
   endif

  else if ((im_critical.ge.1).and. &
           (im_critical.lt.num_materials)) then
   !do nothing
  else
   print *,"im_critical invalid"
   stop
  endif

  if ((im_critical.ge.0).and. &
      (im_critical.lt.num_materials)) then
   ctml_part_id=0
   do part_id=1,TOTAL_NPARTS
    im_part=im_solid_mapF(part_id)+1
    if (CTML_FSI_mat(im_part).eq.1) then 
     ctml_part_id=ctml_part_id+1
     if (ctml_part_id_map(part_id).eq.ctml_part_id) then
      if (im_part.eq.im_critical+1) then

       FSI_output_max_num_nodes=max_num_nodes_list(im_part)
       FSI_output_max_num_elements=max_num_elements_list(im_part)
       FSI_output_num_nodes=num_nodes_list(im_part)
       FSI_output_num_elements=num_elements_list(im_part)

       do inode=1,ctml_max_n_fib_nodes

        do idir=1,NCOMP_FORCE_STRESS
         FSI_output_force_list((inode-1)*NCOMP_FORCE_STRESS+idir)=zero
        enddo

        do idir=1,AMREX_SPACEDIM
         FSI_output_node_list((inode-1)*3+idir)= &
           ctml_fib_pst(ctml_part_id,inode,idir)
         FSI_output_displacement_list((inode-1)*3+idir)=zero
         FSI_output_velocity_halftime_list((inode-1)*3+idir)= &
           ctml_fib_vel(ctml_part_id,inode,idir)
         FSI_output_velocity_list((inode-1)*3+idir)= &
           ctml_fib_vel(ctml_part_id,inode,idir)
         FSI_output_force_list((inode-1)*NCOMP_FORCE_STRESS+idir)= &
           ctml_fib_frc(ctml_part_id,inode,idir)
        enddo !idir=1,sdim

        FSI_output_mass_list(inode)=ctml_fib_mass(ctml_part_id,inode)
        FSI_output_temperature_list(inode)=293.0d0
       enddo ! inode=1,ctml_max_n_fib_nodes 

       if (FSI_output_max_num_elements.eq. &
           ctml_max_n_fib_nodes-1) then
        do ielem=1,ctml_max_n_fib_nodes-1
         do inode=1,AMREX_SPACEDIM
          FSI_output_element_list(4*(ielem-1)+inode)=ielem+inode-1
         enddo
        enddo
       else
        print *,"mismatch:"
        print *,"im_part=",im_part
        print *,"FSI_output_max_num_elements=",  &
          FSI_output_max_num_elements
        print *,"ctml_max_n_fib_nodes ", &
          ctml_max_n_fib_nodes
        stop
       endif
      else if ((im_part.ge.1).and.(im_part.le.num_materials)) then
       ! do nothing
      else
       print *,"im_part invalid"
       stop
      endif
     else
      print *,"ctml_part_id_map(part_id) invalid"
      stop
     endif
    else if (CTML_FSI_mat(im_part).eq.0) then
     ! do nothing
    else
     print *,"CTML_FSI_mat(im_part) invalid"
     stop
    endif
   enddo ! part_id=1,TOTAL_NPARTS

  else
   print *,"im_critical invalid"
   stop
  endif
 
return
end subroutine CLSVOF_ReadHeader

INTEGER_T function sign_valid(mask)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: mask

if ((mask.eq.FSI_FINE_SIGN_VEL_VALID).or. &  
    (mask.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID).or. &  
    (mask.eq.FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID)) then
 sign_valid=1
else if ((mask.eq.FSI_NOTHING_VALID ).or. &
         (mask.eq.FSI_FINE_VEL_VALID).or. &
         (mask.eq.FSI_COARSE_LS_SIGN_VEL_VALID).or. &
         (mask.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID)) then
 sign_valid=0
else
 print *,"mask invalid"
 stop
endif

end function sign_valid


INTEGER_T function vel_valid(mask)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: mask

if ((mask.eq.FSI_FINE_VEL_VALID).or. &  
    (mask.eq.FSI_FINE_SIGN_VEL_VALID).or. & 
    (mask.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID).or. & 
    (mask.eq.FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID).or. &
    (mask.eq.FSI_COARSE_LS_SIGN_VEL_VALID).or. & 
    (mask.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID)) then 
 vel_valid=1
else if (mask.eq.FSI_NOTHING_VALID) then
 vel_valid=0
else
 print *,"mask invalid"
 stop
endif

end function vel_valid


subroutine check_overlap(part_id,ielem,time,minnode,maxnode, &
 tid,tilenum,dx3D,lev77,interior_flag,overlap,isweep)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: tid,tilenum
INTEGER_T, INTENT(in) :: ielem
INTEGER_T, INTENT(in) :: lev77
INTEGER_T, INTENT(out) :: interior_flag
INTEGER_T, INTENT(out) :: overlap
INTEGER_T, INTENT(in) :: isweep
REAL_T, INTENT(in) :: time
REAL_T, INTENT(in) :: dx3D(3)
REAL_T, INTENT(in) :: minnode(3)
REAL_T, INTENT(in) :: maxnode(3)
INTEGER_T local_nelems
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
 if ((local_iband.eq.BoundingBoxRadCell).and. &
     (local_max_side.gt.zero)) then
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
   if (xhi.gt.xlo) then
    ! do nothing
   else
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
   if (minnode(dir).le.maxnode(dir)) then
    ! do nothing
   else
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

subroutine check_overlap_nodeBIG(part_id,inode,time, &
 minnode, &
 tid,tilenum,dx3D,lev77,overlap)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: tid,tilenum
INTEGER_T, INTENT(in) :: inode
INTEGER_T, INTENT(in) :: lev77
INTEGER_T, INTENT(out) :: overlap
REAL_T, INTENT(in) :: time
REAL_T, INTENT(in) :: dx3D(3)
REAL_T, INTENT(in) :: minnode(3)
INTEGER_T local_nnodes
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
     (inode.le.FSI(part_id)%NumNodesBIG)) then

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
    if (xhi.gt.xlo) then
     ! do nothing
    else
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
  print *,"inode invalid7: inode=",inode
  stop
 endif

return
end subroutine check_overlap_nodeBIG


subroutine get_minmax_nodeBIG( &
     FSI_mesh_type, &
     part_id, &
     max_part_id, &
     ielem,time,minnode,maxnode)
IMPLICIT NONE

type(mesh_type), INTENT(in) :: FSI_mesh_type
INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: max_part_id
INTEGER_T, INTENT(in) :: ielem
REAL_T, INTENT(in) :: time
REAL_T, INTENT(out) :: minnode(3)
REAL_T, INTENT(out) :: maxnode(3)
INTEGER_T inode,dir,node_id
REAL_T nodetest
REAL_T xtarget(3)
REAL_T xfoot(3)
REAL_T velparm(3)

 if ((ielem.ge.1).and. &
     (ielem.le.FSI_mesh_type%NumIntElemsBIG)) then
  do inode=1,3
   do dir=1,3
    node_id=FSI_mesh_type%IntElemBIG(inode,ielem)
    if ((node_id.ge.1).and. &
        (node_id.le.FSI_mesh_type%NumNodesBIG)) then
     xfoot(dir)=FSI_mesh_type%NodeBIG(dir,node_id)
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
    velparm,time, &
    FSI_mesh_type, &
    part_id, &
    max_part_id)
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
end subroutine get_minmax_nodeBIG


subroutine get_contained_nodeBIG(part_id,inode,time,minnode)
IMPLICIT NONE

INTEGER_T, INTENT(in) :: part_id
INTEGER_T, INTENT(in) :: inode
REAL_T, INTENT(in) :: time
REAL_T, INTENT(out) :: minnode(3)
INTEGER_T dir
REAL_T xtarget(3)
REAL_T xfoot(3)
REAL_T velparm(3)


 if ((inode.ge.1).and. &
     (inode.le.FSI(part_id)%NumNodesBIG)) then

  do dir=1,3
   xfoot(dir)=FSI(part_id)%NodeBIG(dir,inode)
   velparm(dir)=zero
  enddo
  call get_target_from_foot(xfoot,xtarget, &
   velparm,time, &
   FSI(part_id), &
   part_id, &
   TOTAL_NPARTS)
  do dir=1,3
   minnode(dir)=xtarget(dir)
  enddo ! dir=1..3
 else
  print *,"inode invalid8: inode=",inode
  stop
 endif

return
end subroutine get_contained_nodeBIG

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
 cur_time, &
 dt)
use global_utility_module

IMPLICIT NONE

 INTEGER_T, INTENT(in) :: lev77
 INTEGER_T, INTENT(in) :: sci_max_level
 INTEGER_T, INTENT(in) :: nthread_parm
 REAL_T, INTENT(in) :: dx3D(3)
 INTEGER_T, INTENT(in) :: part_id
 INTEGER_T, INTENT(in) :: im_part
 REAL_T, INTENT(in) :: cur_time
 REAL_T, INTENT(in) :: dt

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
 if ((part_id.lt.1).or.(part_id.ge.num_materials)) then
  print *,"part_id invalid"
  stop
 endif
 if ((im_part.lt.1).or.(im_part.gt.num_materials)) then
  print *,"im_part invalid"
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
 num_nodes=FSI(part_id)%NumNodesBIG
 if (num_nodes.le.0) then
  print *,"num_nodes invalid"
  stop
 endif

 ctml_part_id=ctml_part_id_map(part_id)
 fsi_part_id=fsi_part_id_map(part_id)

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
        call get_minmax_nodeBIG( &
                FSI(part_id), &
                part_id, &
                TOTAL_NPARTS, &
                ielem,cur_time,minnode,maxnode) 
        do tid_loop=1,nthread_parm
        do tilenum_loop=1, &
               contain_elem(lev77)%num_tiles_on_thread3D_proc(tid_loop)
         if ((tid_loop.ne.tid).or.(tilenum_loop.ne.tilenum)) then
           ! isweep==2
          call check_overlap(part_id,ielem,cur_time,minnode,maxnode, &
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

         ! traverse FSI(part_id)%IntElemBIG(inode,ielem)
         ! look at FSI(part_id)%NodeBIG(dir,node_id)
      call get_minmax_nodeBIG( &
              FSI(part_id), &
              part_id, &
              TOTAL_NPARTS, &
              ielem,cur_time,minnode,maxnode)
     
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
      call check_overlap(part_id,ielem,cur_time,minnode,maxnode, &
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
         call check_overlap(part_id,ielem,cur_time,minnode,maxnode, &
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
       print *,"inode invalid9: inode=",inode
       stop
      endif

      call get_contained_nodeBIG(part_id,inode,cur_time,minnode)

       ! isweep==1
       ! check_overlap_nodeBIG increments
       !  contain_elem(lev77)%level_node_data(tid,part_id,tilenum)%numNodes
      call check_overlap_nodeBIG(part_id,inode,cur_time, &
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
          call check_overlap_nodeBIG(part_id,inode,cur_time, &
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
      print *,"FSI(part_id)%NumNodesBIG ",FSI(part_id)%NumNodesBIG
      print *,"FSI(part_id)%NumIntElemsBIG ",FSI(part_id)%NumIntElemsBIG
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

INTEGER_T function is_less_than_list(t1,s1,t2,s2)
IMPLICIT NONE

REAL_T, INTENT(in) :: t1,t2
INTEGER_T, INTENT(in) :: s1,s2

 is_less_than_list=0

 if ((t1.ge.zero).and.(t1.le.one).and. &
     (t2.ge.zero).and.(t2.le.one).and. &
     ((s1.eq.-1).or.(s1.eq.1)).and. &
     ((s2.eq.-1).or.(s2.eq.1))) then
  if (abs(t1-t2).le.crossing_tol) then
   if (s1.eq.s2) then
    ! do nothing
   else if (s1.lt.s2) then
    is_less_than_list=1
   else if (s1.gt.s2) then
    ! do nothing
   else
    print *,"s1 or s2 bust"
    stop
   endif
  else if (t1.lt.t2) then
   is_less_than_list=1
  else if (t1.gt.t2) then
   ! do nothing
  else
   print *,"t1 or t2 is NaN"
   stop
  endif
 else 
  print *,"t1,t2,s1, or s2 invalid"
  stop
 endif

end function is_less_than_list


INTEGER_T function is_equal_list(t1,s1,t2,s2)
IMPLICIT NONE

REAL_T, INTENT(in) :: t1,t2
INTEGER_T, INTENT(in) :: s1,s2

 is_equal_list=0

 if ((t1.ge.zero).and.(t1.le.one).and. &
     (t2.ge.zero).and.(t2.le.one).and. &
     ((s1.eq.-1).or.(s1.eq.1)).and. &
     ((s2.eq.-1).or.(s2.eq.1))) then
  if (abs(t1-t2).le.crossing_tol) then
   if (s1.eq.s2) then
    is_equal_list=1
   else if (s1.lt.s2) then
    ! do nothing
   else if (s1.gt.s2) then
    ! do nothing
   else
    print *,"s1 or s2 bust"
    stop
   endif
  else if (t1.lt.t2) then
   ! do nothing
  else if (t1.gt.t2) then
   ! do nothing
  else
   print *,"t1 or t2 is NaN"
   stop
  endif
 else 
  print *,"t1,t2,s1, or s2 invalid"
  stop
 endif

end function is_equal_list


! nFSI==nparts*NCOMP_FSI
! mask=FSI_NOTHING_VALID or
! mask=FSI_FINE_VEL_VALID or
! mask=FSI_FINE_SIGN_VEL_VALID or
! mask=FSI_DOUBLY_WETTED_SIGN_VEL_VALID or
! mask=FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID or 
! mask=FSI_COARSE_LS_SIGN_VEL_VALID or
! mask=FSI_COARSE_LS_SIGN_FINE_VEL_VALID

! vel=temp=force=0 if no interfaces in cell's neighborhood
! called from: SOLIDFLUID.F90
subroutine CLSVOF_InitBox(  &
  iter, &
  sdim_AMR, & ! =AMREX_SPACEDIM for aux (this variable is not used)
  lev77, &    ! =-1 for aux
  tid, &      ! =0 for aux
  tilenum, &  ! =0 for aux
  im_part, &  ! =1...fort_num_local_aux_grids for aux (auxcomp)
  nparts, &   ! =fort_num_local_aux_grids for aux
  part_id, &  ! =(auxcomp)
  ngrow_make_distance_in, & ! =3 for all cases
  nFSI, & ! =NCOMP_FSI (aux)
  FSI_operation, &
  touch_flag, &
  time, & ! =0.0 for aux
  dt, &   ! =1.0 for aux
  problo3D,probhi3D, & ! unused
  xmap3D, & ! unused
  dx3D, &
  xlo3D_tile, &
  xhi3D_tile, &
  FSI_lo3D,FSI_hi3D, & ! =contain_aux(part_id)%lo3D, hi3D for aux
  FSI_growlo3D,FSI_growhi3D, &
  growlo3D,growhi3D, & ! =FSI_growlo3D,FSI_growhi3D for aux
  xdata3D, &
  FSIdata3D, &
  masknbr3D, &
  ioproc,isout)
use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

  INTEGER_T, INTENT(in) :: iter
  INTEGER_T, INTENT(in) :: sdim_AMR ! not used
  INTEGER_T, INTENT(in) :: lev77
  INTEGER_T, INTENT(in) :: tid
  INTEGER_T, INTENT(in) :: tilenum
  INTEGER_T, INTENT(in) :: im_part
  INTEGER_T, INTENT(in) :: nparts
  INTEGER_T, INTENT(in) :: part_id
  INTEGER_T, INTENT(in) :: ngrow_make_distance_in
  INTEGER_T, INTENT(in) :: nFSI
  INTEGER_T, INTENT(in) :: FSI_operation
  INTEGER_T, INTENT(inout) :: touch_flag
  INTEGER_T :: numtouch
  REAL_T, INTENT(in) :: time,dt
  REAL_T, INTENT(in) :: problo3D(3)
  REAL_T, INTENT(in) :: probhi3D(3)
  INTEGER_T, INTENT(in) :: xmap3D(3)
  REAL_T, INTENT(in) :: dx3D(3)
  REAL_T, INTENT(in) :: xlo3D_tile(3)
  REAL_T, INTENT(in) :: xhi3D_tile(3)
  INTEGER_T, INTENT(in) :: FSI_lo3D(3),FSI_hi3D(3)
  INTEGER_T, INTENT(in) :: FSI_growlo3D(3),FSI_growhi3D(3)
  INTEGER_T, INTENT(in) :: growlo3D(3),growhi3D(3)
  REAL_T, INTENT(in), pointer :: xdata3D(:,:,:,:)
  REAL_T, INTENT(in), pointer :: FSIdata3D(:,:,:,:)
  REAL_T, INTENT(in), pointer :: masknbr3D(:,:,:,:)
  REAL_T, allocatable, target :: old_FSIdata(:,:,:,:)
  REAL_T, pointer :: old_FSIdata_ptr(:,:,:,:)

  INTEGER_T, INTENT(in) :: ioproc,isout

  REAL_T override_LS
  REAL_T override_VEL(3)
  REAL_T override_TEMP
  INTEGER_T override_MASK

  REAL_T dxBB(3) ! set in find_grid_bounding_box
  REAL_T dxBB_min
  REAL_T FSI_delta_cutoff

  INTEGER_T :: ielem
  INTEGER_T :: ielem_container
  INTEGER_T :: nodes_per_elem,inode,nodeptr
  REAL_T, dimension(3) :: xc
  REAL_T, dimension(3) :: xclosest
  REAL_T, dimension(3) :: element_xclosest
  REAL_T, dimension(3) :: xclosest_project
  REAL_T, dimension(3) :: element_normal
  REAL_T, dimension(3) :: normal
  REAL_T, dimension(3) :: normal_closest
  REAL_T, dimension(3) :: ncrit
  REAL_T, dimension(3) :: ncrit_closest
  REAL_T, dimension(3) :: xnot
  REAL_T, dimension(3) :: xfoot
  REAL_T, dimension(3) :: xelem
  REAL_T, dimension(3) :: xtarget
  INTEGER_T, dimension(3) :: gridloBB,gridhiBB,gridlenBB
  INTEGER_T :: i,j,k
  INTEGER_T :: element_inplane
  INTEGER_T :: element_node_edge_inplane
  REAL_T, dimension(3) :: velparm
  REAL_T :: massparm
  INTEGER_T :: dir
  REAL_T :: dotprod
  REAL_T :: unsigned_mindist  ! unsigned
  REAL_T :: element_unsigned_mindist  ! unsigned
  REAL_T :: weighttotal,distwt,weight
  REAL_T, dimension(NCOMP_FSI) :: weight_top 
  REAL_T :: weight_bot
  REAL_T, dimension(3) :: minnode,maxnode
  REAL_T, dimension(3) :: xx
  REAL_T, dimension(3) :: x_outside
  REAL_T, dimension(3) :: xcen
  REAL_T, dimension(3) :: xleft
  REAL_T, dimension(3) :: xright
  REAL_T, dimension(3) :: xside
  REAL_T, dimension(3) :: xcrit
  REAL_T, dimension(3) :: xcrit_project
  INTEGER_T :: ii,jj,kk
  INTEGER_T :: hitflag
  REAL_T :: phiside
  REAL_T :: phicen
  REAL_T :: testdist
  REAL_T :: hitsign
  REAL_T :: totaldist

  REAL_T mag_n,mag_n_test
  REAL_T n_dot_x
  REAL_T mag_ncrit
  REAL_T mag_x
  REAL_T sign_conflict
  REAL_T sign_conflict_local
  INTEGER_T ibase
  REAL_T ls_local
  INTEGER_T mask_local,mask_node
  INTEGER_T new_mask_local
  REAL_T, dimension(3) :: vel_local
  REAL_T, dimension(NCOMP_FORCE_STRESS) :: force_local
  REAL_T temp_local
  INTEGER_T nc
  INTEGER_T ii_current

  REAL_T inner_band_size
  REAL_T LS_sum,LS_SIDE,dx_SIDE
  REAL_T total_variation_sum_plus
  REAL_T total_variation_sum_minus
  INTEGER_T weight_total_variation

  INTEGER_T num_plane_intersects
  INTEGER_T num_plane_intersects_new
  INTEGER_T cur_ptr
  INTEGER_T cur_sign
  REAL_T t_top,t_bottom,t_crit,swap_data
  INTEGER_T swap_sign
  INTEGER_T num_sign_changes,local_corner_count,local_smooth_count
  REAL_T far_field_sign
  REAL_T near_field_sign
  INTEGER_T less_than_flag,equal_flag
  REAL_T plane_diff

  REAL_T plane_intersect_list(max_plane_intersects)
  ! sign on the far field side
  ! sign=-n dot (xfar-x0)
  INTEGER_T plane_intersect_list_sign(max_plane_intersects)

  INTEGER_T sign_status_changed
  INTEGER_T num_elements
  INTEGER_T num_elements_container

  REAL_T element_scale
  REAL_T test_scale
  REAL_T test_scale_max
  REAL_T eul_over_lag_scale

  INTEGER_T debug_all
  INTEGER_T mask1,mask2
  INTEGER_T null_intersection
  INTEGER_T ctml_part_id
  INTEGER_T fsi_part_id
  INTEGER_T local_iband
  INTEGER_T CTML_DEBUG_Mass
  type(mesh_type) :: FSI_mesh_type

  call FLUSH(6)  ! 6=screen

  debug_all=0

  CTML_DEBUG_Mass=0

  if ((part_id.lt.1).or.(part_id.gt.nparts)) then
   print *,"part_id invalid"
   stop
  endif

  if ((tid.lt.0).or.(tilenum.lt.0)) then
   print *,"tid or tilenum invalid"
   stop
  endif

  if (lev77.eq.-1) then ! aux Lagrangian data

   if (nparts.ne.fort_num_local_aux_grids) then
    print *,"nparts.ne.fort_num_local_aux_grids"
    stop
   endif

   FSI_mesh_type=aux_FSI(part_id)

  else if (lev77.ge.1) then

   if (nparts.ne.TOTAL_NPARTS) then
    print *,"nparts.ne.TOTAL_NPARTS"
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

   FSI_mesh_type=FSI(part_id)

  else
   print *,"lev77 invalid"
   stop
  endif

  local_iband=FSI_mesh_type%bounding_box_ngrow
  if (local_iband.eq.BoundingBoxRadCell) then
   ! do nothing
  else
   print *,"local_iband invalid"
   stop
  endif

  do dir=1,3

   if (lev77.eq.-1) then

    if (contain_aux(part_id)%lo3D(dir).ne.FSI_lo3D(dir)) then
     print *,"contain_aux(part_id)%lo3D(dir).ne.FSI_lo3D(dir)"
     stop
    endif
    if (contain_aux(part_id)%hi3D(dir).ne.FSI_hi3D(dir)) then
     print *,"contain_aux(part_id)%hi3D(dir).ne.FSI_hi3D(dir)"
     stop
    endif

   else if (lev77.ge.1) then

    if (contain_elem(lev77)%tilelo3D(tid+1,tilenum+1,dir).ne. &
        FSI_lo3D(dir)) then
     print *,"tilelo3D(tid+1,tilenum+1,dir).ne.FSI_lo3D(dir)"
     stop
    endif 
    if (contain_elem(lev77)%tilehi3D(tid+1,tilenum+1,dir).ne. &
        FSI_hi3D(dir)) then
     print *,"tilehi3D(tid+1,tilenum+1,dir).ne.FSI_hi3D(dir)"
     stop
    endif 

   else
    print *,"lev77 invalid"
    stop
   endif

   if (FSI_growlo3D(dir).lt.FSI_lo3D(dir)) then
    ! do nothing
   else
    print *,"FSI_growlo3D(dir) invalid"
    stop
   endif
   if (FSI_growhi3D(dir).gt.FSI_hi3D(dir)) then
    ! do nothing
   else
    print *,"FSI_growhi3D(dir) invalid"
    stop
   endif
   if (xlo3D_tile(dir).lt.xhi3D_tile(dir)) then
    ! do nothing
   else
    print *,"xlo3D_tile(dir).lt.xhi3D_tile(dir) violated"
    stop
   endif

  enddo ! dir=1..3

  if (lev77.eq.-1) then
   ctml_part_id=0
   fsi_part_id=part_id
  else if (lev77.ge.1) then
   ctml_part_id=ctml_part_id_map(part_id)
   fsi_part_id=fsi_part_id_map(part_id)
  else
   print *,"lev77 invalid"
   stop
  endif

  if (((ctml_part_id.ge.1).and. &
       (ctml_part_id.le.CTML_NPARTS)).or. &
      ((fsi_part_id.ge.1).and. &
       (fsi_part_id.le.nparts))) then

   num_elements=FSI_mesh_type%NumIntElemsBIG

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
   if ((part_id.lt.1).or.(part_id.gt.nparts)) then
    print *,"part_id invalid"
    stop
   endif
   if ((num_materials.lt.1).or.(num_materials.gt.50)) then
    print *,"num_materials out of range"
    stop
   endif

   if (lev77.eq.-1) then

    if (im_part.eq.part_id) then
     ! do nothing
    else
     print *,"im_part invalid"
     stop
    endif
    if (nFSI.ne.NCOMP_FSI) then
     print *,"nFSI invalid lev77=-1"
     stop
    endif
    ibase=0

   else if (lev77.ge.1) then

    if ((im_part.lt.1).or.(im_part.gt.num_materials)) then
     print *,"im_part invalid"
     stop
    endif
    if (nFSI.ne.nparts*NCOMP_FSI) then
     print *,"nFSI invalid lev77>=1"
     stop
    endif
    if ((nparts.lt.1).or.(nparts.ge.num_materials)) then
     print *,"nparts invalid"
     stop
    endif
    ibase=(part_id-1)*NCOMP_FSI

   else
    print *,"lev77 invalid"
    stop
   endif

   if (isout.eq.0) then
    ! do nothing
   else if (isout.eq.1) then
    print *,"in (START): CLSVOF_InitBox"
    print *,"num_elements(NumIntElemsBIG)=",num_elements
    print *,"NumIntElems=",FSI_mesh_type%NumIntElems
    print *,"sdim_AMR=",sdim_AMR 
    print *,"im_part=",im_part
    print *,"part_id=",part_id
    print *,"num_materials=",num_materials
    print *,"FSI_operation=",FSI_operation
    print *,"iter= ",iter
    print *,"touch_flag=",touch_flag
    print *,"time,dt ",time,dt
    print *,"ioproc=",ioproc
    print '(A9,3(f12.6))',"problo3D ",problo3D(1),problo3D(2),problo3D(3)
    print '(A9,3(f12.6))',"probhi3D ",probhi3D(1),probhi3D(2),probhi3D(3)
    print '(A7,3(I10))',"FSI_lo3D ",FSI_lo3D(1),FSI_lo3D(2),FSI_lo3D(3)
    print '(A7,3(I10))',"FSI_hi3D ",FSI_hi3D(1),FSI_hi3D(2),FSI_hi3D(3)
    print '(A11,3(I10))',"FSI_growlo3D ",FSI_growlo3D(1),FSI_growlo3D(2),FSI_growlo3D(3)
    print '(A11,3(I10))',"FSI_growhi3D ",FSI_growhi3D(1),FSI_growhi3D(2),FSI_growhi3D(3)
    print '(A11,3(I10))',"growlo3D ",growlo3D(1),growlo3D(2),growlo3D(3)
    print '(A11,3(I10))',"growhi3D ",growhi3D(1),growhi3D(2),growhi3D(3)
   else
    print *,"isout invalid1: ",isout
    stop
   endif

   if ((touch_flag.ne.0).and.(touch_flag.ne.1)) then
    print *,"touch_flag invalid"
    stop
   endif
   if ((FSI_operation.ne.2).and.(FSI_operation.ne.3)) then
    print *,"FSI_operation invalid"
    stop
   endif
   if (ngrow_make_distance_in.ne.3) then
    print *,"ngrow_make_distance_in invalid"
    stop
   endif

   call checkbound3D_array(FSI_lo3D,FSI_hi3D, &
    FSIdata3D, &
    ngrow_make_distance_in,-1)
   call checkbound3D_array(FSI_lo3D,FSI_hi3D, &
    xdata3D, &
    ngrow_make_distance_in,-1)
   call checkbound3D_array(FSI_lo3D,FSI_hi3D, &
    masknbr3D, &
    ngrow_make_distance_in,-1)

   if (lev77.eq.-1) then
    num_elements_container=FSI_mesh_type%NumIntElemsBIG
   else if (lev77.ge.1) then           
    num_elements_container=contain_elem(lev77)% &
         level_elem_data(tid+1,part_id,tilenum+1)% &
         numElems
   else
    print *,"lev77 invalid"
    stop
   endif

   if ((num_elements_container.lt.0).or. &
       (num_elements_container.gt.num_elements)) then
    print *,"num_elements_container invalid"
    stop
   endif

   local_corner_count=0
   local_smooth_count=0

   if (FSI_operation.eq.OP_FSI_MAKE_DISTANCE) then 

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

    do ielem_container=1,num_elements_container

     if (lev77.eq.-1) then
      ielem=ielem_container
     else if (lev77.ge.1) then           
      ielem=contain_elem(lev77)% &
           level_elem_data(tid+1,part_id,tilenum+1)% &
           ElemData(ielem_container)
     else
      print *,"lev77 invalid"
      stop
     endif

     if ((ielem.lt.1).or. &
         (ielem.gt.num_elements)) then
      print *,"ielem invalid"
      stop
     endif

     do dir=1,3
      xelem(dir)=FSI_mesh_type%ElemDataXnotBIG(dir,ielem)
      velparm(dir)=zero
     enddo 
     call get_target_from_foot(xelem,xnot, &
       velparm,time, &
       FSI_mesh_type, &
       part_id, &
       nparts)

     nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,ielem)
     if (nodes_per_elem.lt.3) then
      print *,"elem,nodes_per_elem ",ielem,nodes_per_elem   
      stop
     endif
     ! normal points from solid to fluid
     ! phi=n dot (x-xnot)
     ! phi>0 in the fluid (the sign will be switched later)
     ! this is the element normal (in contrast to the node normal)
     call scinormalBIG(ielem, &
       normal, &
       FSI_mesh_type, &
       part_id, &
       nparts, &
       time)

     test_scale=sqrt(normal(1)**2+normal(2)**2+normal(3)**2)

     if (abs(test_scale-one).le.VOFTOL) then

      do dir=1,3
       element_normal(dir)=normal(dir)
      enddo

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
       ! minnode and maxnode are needed in order to find the
       ! stencil of surrounding Eulerian cells to lagrangian element (triangle)
      call get_minmax_nodeBIG( &
        FSI_mesh_type, &
        part_id, &
        nparts, &
        ielem, &
        time, &
        minnode,maxnode)

      test_scale_max=zero
      do dir=1,3
       test_scale=maxnode(dir)-minnode(dir)
       if (test_scale.ge.zero) then
        test_scale_max=max(test_scale_max,test_scale)
       else
        print *,"test_scale.lt.zero"
        print *,"dir,minnode,maxnode,test_scale ",dir,minnode(dir), &
         maxnode(dir),test_scale
        print *,"part_id,ielem,time ",part_id,ielem,time
        stop
       endif
      enddo ! dir=1..3

      element_scale=FSI_mesh_type%min_side_len_refined

      if (element_scale.gt.zero) then
       ! do nothing
      else
       print *,"element_scale.le.zero: ",element_scale
       print *,"part_id= ",part_id
       print *,"ielem= ",ielem
       print *,"time= ",time
       do dir=1,3
        print *,"dir,minnode,maxnode ",dir,minnode(dir),maxnode(dir)
       enddo
       stop
      endif
     
       ! stencil of surrounding Eulerian cells to lagrangian element (triangle)
       ! BoundingBoxRadCell determines buffer zone.
      call find_grid_bounding_box( &
        FSI_mesh_type, &
        part_id, &
        nparts, &
        null_intersection, &
        minnode,maxnode, &
        FSI_lo3D,FSI_hi3D, &
        FSI_growlo3D,FSI_growhi3D, &
        growlo3D,growhi3D, &
        xdata3D, &
        gridloBB,gridhiBB, &
        dxBB) ! for spectral methods, the element size is bfact*dxBB

      dxBB_min=dxBB(1)
      do dir=1,3
       if (dxBB(dir).lt.dxBB_min) then
        dxBB_min=dxBB(dir)
       endif
      enddo
      FSI_delta_cutoff=3.0d0*dxBB_min

      if (test_scale_max.gt.zero) then
       eul_over_lag_scale=min(one,dxBB_min/test_scale_max)
      else
       print *,"test_scale_max must be positive"
       stop
      endif

      do dir=1,3
       if (abs(dxBB(dir)-dx3D(dir)).le.VOFTOL*dxBB(dir)) then
        ! do nothing
       else
        print *,"abs(dxBB(dir)-dx3D(dir)).gt.VOFTOL*dxBB(dir)"
        stop
       endif
      enddo ! dir=1..3

      if (1.eq.0) then
       print *,"ielem=",ielem
       do dir=1,3
        print *,"dir,gridloBB,gridhiBB ",dir,gridloBB(dir),gridhiBB(dir)
       enddo
      endif

      if (null_intersection.eq.0) then

        ! sanity check
       do dir=1,3
        gridlenBB(dir)=gridhiBB(dir)-gridloBB(dir)+1
        if ((gridlenBB(dir).ge.1).and. &
            (gridlenBB(dir).le.2048)) then
         ! do nothing
        else
         print *,"gridlenBB(dir) invalid"
         stop
        endif
       enddo ! dir=1..3

        ! LOOP through bounding box of the element.
        ! this code is thread safe
        ! gridloBB,gridhiBB restricted to growlo3D and growhi3D 
       do i=gridloBB(1),gridhiBB(1)
       do j=gridloBB(2),gridhiBB(2)
       do k=gridloBB(3),gridhiBB(3)

        if ((i.ge.FSI_lo3D(1)).and.(i.le.FSI_hi3D(1)).and. &
            (j.ge.FSI_lo3D(2)).and.(j.le.FSI_hi3D(2)).and. &
            (k.ge.FSI_lo3D(3)).and.(k.le.FSI_hi3D(3))) then

         mask1=NINT(masknbr3D(i,j,k,1))
         mask2=NINT(masknbr3D(i,j,k,2))

          ! masknbr3D derived from localMF[MASK_NBR_MF]
          ! (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
          ! (2) =1 interior  =0 otherwise
          ! (3) =1 interior+ngrow-1  =0 otherwise
          ! (4) =1 interior+ngrow    =0 otherwise

         if ((mask1.eq.1).and. &
             (mask2.eq.1)) then
          ! do nothing
         else
          print *,"expecting mask1=mask2=1"
          print *,"mask1: ",mask1
          print *,"mask2: ",mask2
          stop
         endif

         do dir=1,3
          xx(dir)=xdata3D(i,j,k,dir)
         enddo

         ! normal points from solid to fluid
         ! -normal points from fluid to solid
         dotprod=zero
         do dir=1,3
          dotprod=dotprod+normal(dir)*(xx(dir)-xnot(dir))
         enddo

         do dir=1,3
          xclosest(dir)=xx(dir)-dotprod*normal(dir)

          element_xclosest(dir)=xclosest(dir)

          normal_closest(dir)=normal(dir)
          xclosest_project(dir)=xclosest(dir)
         enddo ! dir=1..3

         unsigned_mindist=abs(dotprod)

         element_unsigned_mindist=unsigned_mindist

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

         call checkinplaneBIG( &
          eul_over_lag_scale, &
          xx, & ! target point at which the signed distance is sought.
          xclosest, &
          xclosest_project, & !INTENT(out)
          normal, & ! INTENT(in)
          normal_closest, & ! INTENT(out)
          ielem, &
          element_node_edge_inplane, & !INTENT(out)
          FSI_mesh_type, &
          part_id, &
          nparts, &
          time)

         element_inplane=element_node_edge_inplane

! investigate using NodeNormalBIG (normal defined at nodes)
         do inode=1,nodes_per_elem
          ! check distance to the edges of a triangular element.
          call checkinlineBIG( &
           eul_over_lag_scale, &
           xclosest_project, & !INTENT(inout)
           normal_closest, & ! INTENT(inout), initially normal of element.
           inode,ielem, &
           unsigned_mindist, & !INTENT(inout)
           xx, & ! target point at which the signed distance is sought.
           element_node_edge_inplane, &  !INTENT(inout)
           FSI_mesh_type, &
           part_id, &
           nparts, &
           time,dxBB)
          ! check distance to the nodes of a triangular element.
          ! normal_closest is the element normal.
          call checkinpointBIG( &
           xclosest_project, & !INTENT(inout)
           normal_closest, & ! INTENT(inout)
           inode,ielem, &
           unsigned_mindist, & ! INTENT(inout)
           xx, & ! target point at which the signed distance is sought.
           element_node_edge_inplane, &  !INTENT(inout)
           FSI_mesh_type, &
           part_id, &
           nparts, &
           time,dxBB)
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
          if (xright(dir).gt.xleft(dir)) then
           ! do nothing
          else
           print *,"xright(dir).le.xleft(dir)"
           stop
          endif
         enddo ! dir=1..3

         ! normal points from solid to fluid
         ! phi>0 in the fluid (sign will be switched later)
         ! xclosest_project=xx-phi n
         ! phi n = xx-xclosest_project
         ! phi = n dot (xx-xclosest_project) (sign will be switched later)
         if (element_node_edge_inplane.eq.1) then

          n_dot_x=zero
          mag_n=zero
          mag_n_test=zero
          mag_x=zero

          do dir=1,3
           mag_n=mag_n+normal_closest(dir)**2
           mag_n_test=mag_n_test+normal(dir)**2
           !xx=cell center where LS needed
           mag_x=mag_x+(xx(dir)-xclosest_project(dir))**2 
           
           n_dot_x=n_dot_x+ &
             normal_closest(dir)*(xx(dir)-xclosest_project(dir))
      
          enddo ! dir=1..3

          mag_n_test=sqrt(mag_n_test)
          mag_n=sqrt(mag_n)
          mag_x=sqrt(mag_x)

          if (abs(mag_x).lt.unsigned_mindist) then
           unsigned_mindist=abs(mag_x)
          endif

           ! at this stage, normal points from solid to fluid.
           ! mag_x=||xx-xclosest_project||
           ! mag_n=||normal_closest||
           ! n dot x=normal_closest dot (xx-xclosest_project)
           ! n dot x = mag_x * mag_n * cos(theta)

          if (mag_n.ge.zero) then
           if (mag_n_test.ge.zero) then
            if (mag_x.ge.zero) then

             if (element_inplane.eq.1) then

              n_dot_x=zero
              mag_n=zero
              mag_x=zero
              do dir=1,3
               mag_n=mag_n+element_normal(dir)**2
               mag_x=mag_x+(xx(dir)-element_xclosest(dir))**2 
               n_dot_x=n_dot_x+ &
                element_normal(dir)*(xx(dir)-element_xclosest(dir))

               normal_closest(dir)=element_normal(dir)
              enddo ! dir=1..3

              mag_n=sqrt(mag_n)
              mag_x=sqrt(mag_x)

              if (abs(mag_x).le.unsigned_mindist) then
               xclosest_project(dir)=element_xclosest(dir)
               unsigned_mindist=abs(mag_x)
              endif

             else if (element_inplane.eq.0) then

              ! do nothing

             else
              print *,"element_inplane invalid"
              stop
             endif

             hitflag=1
             if (n_dot_x.eq.zero) then
              hitsign=zero
             else if (n_dot_x.gt.zero) then ! fluid(sign switched later)
              hitsign=one
             else if (n_dot_x.lt.zero) then ! solid(sign switched later)
              hitsign=-one
             else
              print *,"n_dot_x bust"
              stop
             endif
    
            else
             print *,"mag_x invalid"
             stop
            endif

           else
            print *,"mag_n_test invalid"
            stop
           endif

          else
           print *,"mag_n invalid"
           stop
          endif

         else if (element_node_edge_inplane.eq.0) then
          print *,"checkinpointBIG should guarantee inplane=1, aborting"
          stop
         else
          print *,"element_node_edge_inplane invalid"
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
           ! phi>0 in the fluid (sign will be switched later)
           phicen=zero
           phiside=zero
           do dir=1,3
             ! levelset at the center cell
            phicen=phicen+normal(dir)*(xx(dir)-xnot(dir))
             ! levelset at the side cell
            phiside=phiside+normal(dir)*(xside(dir)-xnot(dir))
           enddo
           if ((phicen.eq.zero).and. &
               (phiside.eq.zero)) then
            ! do nothing: either normal=0 or interface parallel to the edge
           else if (phicen*phiside.gt.zero) then
            ! do nothing (no crossing)
           else if (phicen*phiside.le.zero) then
            mag_ncrit=zero
            do dir=1,3
             if (phicen.gt.phiside) then ! i.e. phicen>=0 (fluid), phiside<=0
              ncrit(dir)=(xx(dir)-xside(dir)) 
             else if (phicen.lt.phiside) then ! phicen<=0, phiside>=0 (fluid)
              ncrit(dir)=(xside(dir)-xx(dir))
             else
              print *,"phicen or phiside NaN"
              stop
             endif
             mag_ncrit=mag_ncrit+ncrit(dir)**2
            enddo
            mag_ncrit=sqrt(mag_ncrit)
            if (mag_ncrit.gt.zero) then
             do dir=1,3
              ncrit(dir)=ncrit(dir)/mag_ncrit
             enddo
            else
             print *,"mag_ncrit is NaN"
             stop
            endif
            do dir=1,3
             if (phiside.eq.zero) then
              xcrit(dir)=xside(dir)
             else if (phicen.eq.zero) then
              xcrit(dir)=xx(dir)
             else if ((phiside.ne.zero).and. &
                      (phicen.ne.zero)) then
              xcrit(dir)=(abs(phiside)*xx(dir)+ &
                           abs(phicen)*xside(dir))/  &
                          (abs(phicen)+abs(phiside))
             else
              print *,"phiside or phicen is NaN"
              stop
             endif
            enddo  ! dir=1..3

            call checkinplaneBIG( &
             eul_over_lag_scale, &
             xx, & ! target point at which the signed distance is sought.
             xcrit, & !INTENT(in)
             xcrit_project, & ! INTENT(out)
             ncrit, &         ! INTENT(in)
             ncrit_closest, & ! INTENT(out)
             ielem, &
             element_node_edge_inplane, &
             FSI_mesh_type, &
             part_id, &
             nparts, &
             time)
! totaldist is the distance between xx and xside
! testdist  is the distance between xx and xcrit 
! xx=center point  xside=stencil point xcrit=crossing point
            if (element_node_edge_inplane.eq.1) then
             testdist=zero
             totaldist=zero
             do dir=1,3
              testdist=testdist+(xcrit(dir)-xx(dir))**2
              totaldist=totaldist+(xside(dir)-xx(dir))**2
             enddo
             testdist=sqrt(testdist)
             totaldist=sqrt(totaldist)
             if (testdist.le.unsigned_mindist) then
              do dir=1,3
               normal_closest(dir)=ncrit(dir)
               xclosest_project(dir)=xcrit(dir)
              enddo
              unsigned_mindist=testdist
              hitflag=1
              if (phicen.eq.zero) then
               hitsign=zero
               ! n_dot_x>0 in the fluid (sign will be switched later)
              else if (phicen.gt.zero) then
               hitsign=one
               ! n_dot_x<0 in the solid (sign will be switched later)
              else if (phicen.lt.zero) then
               hitsign=-one
              else
               print *,"n_dot_x bust"
               stop
              endif
             else if (testdist.gt.unsigned_mindist) then
              ! do nothing
             else
              print *,"testdist or unsigned_mindist bust"
              stop
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

            else if (element_node_edge_inplane.eq.0) then
             ! do nothing
            else
             print *,"inplane invalid"
             stop
            endif  ! inplane
           else
            print *,"phicen or phiside is NaN"
            stop
           endif 
          endif ! abs(ii)+abs(jj)+abs(kk)>0
         enddo 
         enddo 
         enddo  ! ii,jj,kk

         ls_local=FSIdata3D(i,j,k,ibase+FSI_LEVELSET+1)

         sign_conflict_local=FSIdata3D(i,j,k,ibase+FSI_SIGN_CONFLICT+1)
         sign_conflict=sign_conflict_local

         mask_local=NINT(FSIdata3D(i,j,k,ibase+FSI_EXTRAP_FLAG+1))
         do dir=1,3
          vel_local(dir)=FSIdata3D(i,j,k,ibase+FSI_VELOCITY+dir)
         enddo
         temp_local=FSIdata3D(i,j,k,ibase+FSI_TEMPERATURE+1)

         if (hitflag.eq.1) then
          ! =0 no hits
          ! =1 positive
          ! =2 negative
          ! =3 inconclusive
          if (abs(unsigned_mindist).le.VOFTOL*dx3D(1)) then
           unsigned_mindist=zero
           sign_conflict_local=one
          else if (abs(unsigned_mindist).ge.VOFTOL*dx3D(1)) then
           if (abs(unsigned_mindist).le.0.01d0*dx3D(1))  then
            if (hitsign.ge.zero) then
             sign_conflict_local=one
            else if (hitsign.lt.zero) then
             sign_conflict_local=two
            else
             print *,"hitsign invalid"
             stop
            endif
           else if (abs(unsigned_mindist).ge.0.01d0*dx3D(1)) then
            mag_n=zero
            mag_n_test=zero
            n_dot_x=zero
            do dir=1,3
!mag_n=mag_n+normal_closest(dir)**2
             mag_n=mag_n+element_normal(dir)**2
             mag_n_test=mag_n_test+(xclosest_project(dir)-xx(dir))**2
!n_dot_x=n_dot_x+normal_closest(dir)*(xclosest_project(dir)-xx(dir))
             n_dot_x=n_dot_x+element_normal(dir)*(xclosest_project(dir)-xx(dir))
            enddo
            mag_n_test=sqrt(mag_n_test)
            mag_n=sqrt(mag_n)
            if (mag_n_test.gt.zero) then
             if (mag_n.eq.zero) then
              sign_conflict_local=three
             else if (mag_n.gt.zero) then
              n_dot_x=n_dot_x/(mag_n*mag_n_test)

              if (abs(n_dot_x).lt.cos(angle_tol*Pi/180.0d0)) then
               sign_conflict_local=three
              else if (hitsign.ge.zero) then
               sign_conflict_local=one
              else if (hitsign.lt.zero) then
               sign_conflict_local=two
              else
               print *,"hitsign invalid"
               stop
              endif

             else
              print *,"mag_n is corrupt"
              stop
             endif
            else
             print *,"mag_n_test invalid"
             stop
            endif
           else
            print *,"unsigned_mindist is corrupt"
            stop
           endif
          else
           print *,"unsigned_mindist is corrupt"
           stop
          endif

         else
          print *,"hitflag invalid"
          stop
         endif 

         if ((mask_local.eq.FSI_NOTHING_VALID).or. &  
             (mask_local.eq.FSI_COARSE_LS_SIGN_VEL_VALID).or. & 
             (unsigned_mindist.lt.abs(ls_local))) then

          if ((mask_local.eq.FSI_NOTHING_VALID).or. &
              (mask_local.eq.FSI_COARSE_LS_SIGN_VEL_VALID)) then
           ! do nothing
          else if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID) then 
           ! do nothing
          else if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID) then 
           print *,"mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID invalid here"
           stop
          else if (mask_local.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID) then 
           ! do nothing
          else if (mask_local.eq.FSI_FINE_VEL_VALID) then 
           ! do nothing
          else if (mask_local.eq.FSI_FINE_SIGN_VEL_VALID) then 
           print *,"mask_local.eq.FSI_FINE_SIGN_VEL_VALID invalid here"
           stop
          else
           print *,"mask_local invalid"
           stop
          endif

          if ((ctml_part_id.ge.1).and. &
              (ctml_part_id.le.CTML_NPARTS)) then
           if (FSI_mesh_type%ElemDataBIG(DOUBLYCOMP,ielem).eq.1) then 
            ! do nothing
           else
            print *,"expecting doubly wetted"
            print *,"ctml_part_id=",ctml_part_id
            print *,"ielem=",ielem
            print *,"FSI_mesh_type%ElemDataBIG(DOUBLYCOMP,ielem)=", &
              FSI_mesh_type%ElemDataBIG(DOUBLYCOMP,ielem)
            stop
           endif
          else if (ctml_part_id.eq.0) then
           ! do nothing
          else
           print *,"ctml_part_id invalid"
           stop
          endif

          if (FSI_mesh_type%ElemDataBIG(DOUBLYCOMP,ielem).eq.1) then 

           mask_local=FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID
           ls_local=-unsigned_mindist
           sign_conflict=sign_conflict_local

          else if (FSI_mesh_type%ElemDataBIG(DOUBLYCOMP,ielem).eq.0) then

           if (FSI_mesh_type%exclusive_doubly_wetted.eq.0) then
            ! do nothing
           else
            print *,"doubly wetted flags are inconsistent"
            stop
           endif

           if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID) then 
            print *, &
              "mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID invalid here"
            stop
           endif

           if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID) then 
            mask_local=FSI_NOTHING_VALID
           endif

           if (mask_local.eq.FSI_NOTHING_VALID) then
            mask_local=FSI_FINE_VEL_VALID
           else if (mask_local.eq.FSI_COARSE_LS_SIGN_VEL_VALID) then
            mask_local=FSI_COARSE_LS_SIGN_FINE_VEL_VALID
           else if (mask_local.eq.FSI_FINE_SIGN_VEL_VALID) then 
            print *,"mask_local.eq.FSI_FINE_SIGN_VEL_VALID invalid here"
            stop
           else if ((mask_local.eq.FSI_FINE_VEL_VALID).or. &
                    (mask_local.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID)) then
            ! do nothing
           else
            print *,"mask_local invalid"
            stop
           endif

           if ((mask_local.eq.FSI_NOTHING_VALID).or. &
               (mask_local.eq.FSI_FINE_VEL_VALID)) then
            ls_local=unsigned_mindist
            sign_conflict=sign_conflict_local
           else if (mask_local.eq.FSI_FINE_SIGN_VEL_VALID) then 
            print *,"mask_local.eq.FSI_FINE_SIGN_VEL_VALID invalid here"
            stop
           else if ((mask_local.eq.FSI_COARSE_LS_SIGN_VEL_VALID).or. &
                    (mask_local.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID)) then
            if (ls_local.le.zero) then
             sign_conflict=one
             ls_local=-unsigned_mindist 
            else if (ls_local.gt.zero) then
             sign_conflict=two
             ls_local=unsigned_mindist 
            else
             print *,"ls_local invalid"
             stop
            endif
           else
            print *,"mask_local invalid"
            stop
           endif 

           if (hitflag.eq.1) then

            if ((sign_conflict_local.eq.one).or. &
                (sign_conflict_local.eq.two).or. &
                (sign_conflict_local.eq.three)) then
             ! do nothing
            else
             print *,"sign_conflict_local invalid"
             stop
            endif

             ! The "-" below asserts that ls_local now points into the solid
            ls_local=-hitsign*abs(ls_local)
            sign_conflict=sign_conflict_local

            if ((mask_local.eq.FSI_NOTHING_VALID).or. &
                (mask_local.eq.FSI_COARSE_LS_SIGN_VEL_VALID).or. &
                (mask_local.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID).or. &
                (mask_local.eq.FSI_FINE_VEL_VALID)) then
             mask_local=FSI_COARSE_LS_SIGN_FINE_VEL_VALID
            else if (mask_local.eq.FSI_FINE_SIGN_VEL_VALID) then 
             print *,"mask_local.eq.FSI_FINE_SIGN_VEL_VALID invalid here"
             stop
            else if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID) then 
             print *,"not a doubly wetted element"
             stop
            else if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID) then 
             print *,"not a doubly wetted element"
             stop
            else 
             print *,"mask_local invalid"
             stop
            endif
           else if (hitflag.eq.0) then
            print *,"expecting hitflag=1"
            stop
           else
            print *,"hitflag invalid"
            stop
           endif 
          else
           print *,"wetted flag invalid"
           stop
          endif

          do dir=1,3
           vel_local(dir)=0.0d0
          enddo
          do dir=1,NCOMP_FORCE_STRESS
           force_local(dir)=0.0d0
          enddo
          temp_local=0.0d0

          weighttotal=0.0d0
          do inode=1,nodes_per_elem
           nodeptr=FSI_mesh_type%IntElemBIG(inode,ielem)
           do dir=1,3
            xfoot(dir)=FSI_mesh_type%NodeBIG(dir,nodeptr)
            velparm(dir)=FSI_mesh_type%NodeVelBIG(dir,nodeptr)
           enddo
           massparm=FSI_mesh_type%NodeMassBIG(nodeptr)
           if (massparm.gt.zero) then
            ! do nothing
           else
            print *,"massparm invalid"
            stop
           endif
           call get_target_from_foot(xfoot,xtarget, &
             velparm,time, &
             FSI_mesh_type, &
             part_id, &
             nparts)
   
             ! xtarget is Lagrangian coordinate
             ! xx is grid coordinate 
           if (CTML_DEBUG_Mass.eq.1) then
            print *,"inode,ielem,xtarget,xx,massparm ", &
              inode,ielem, &
              xtarget(1),xtarget(2),xtarget(3), &  
              xx(1),xx(2),xx(3),massparm
           endif
           distwt=0.0d0
           do dir=1,3
            distwt=distwt+(xtarget(dir)-xx(dir))**2
           enddo
           weight=1.0d0/( (distwt+(1.0E-10)**2)**4 )
           weight=weight*massparm
           do dir=1,3
            vel_local(dir)=vel_local(dir)+weight*velparm(dir)
           enddo
           do dir=1,NCOMP_FORCE_STRESS
            force_local(dir)=force_local(dir)+ &
             weight*FSI_mesh_type%NodeForceBIG(dir,nodeptr)
           enddo
           temp_local=temp_local+ & 
            weight*FSI_mesh_type%NodeTempBIG(nodeptr)
           weighttotal=weighttotal+weight
          enddo ! inode=1..nodes_per_elem

          do dir=1,3
           vel_local(dir)=vel_local(dir)/weighttotal
          enddo
          do dir=1,NCOMP_FORCE_STRESS
           force_local(dir)=force_local(dir)/weighttotal
          enddo

          if ((probtype.eq.9).and.(axis_dir.gt.1)) then
           if ( xx(1) .gt. -0.4 .and. xx(1) .lt. -0.2 ) then
            if ( xx(3) .lt. 0.045 .and. xx(3) .gt. 0.03 ) then
             vel_local(3) = -1.0
            endif
           endif
          endif

          temp_local=temp_local/weighttotal
         else if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID) then 

          print *,"mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID invalid here"
          stop

         else if ( ((mask_local.eq.FSI_FINE_VEL_VALID).and. &
                    (unsigned_mindist.ge.abs(ls_local))).or. &
                   ((mask_local.eq.FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID).and. &
                    (unsigned_mindist.ge.abs(ls_local))).or. &
                   ((mask_local.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID).and. &
                    (unsigned_mindist.ge.abs(ls_local))) ) then
          ! do nothing
         else
          print *,"mask_local invalid"
          stop
         endif 

         FSIdata3D(i,j,k,ibase+FSI_LEVELSET+1)=ls_local
         FSIdata3D(i,j,k,ibase+FSI_SIGN_CONFLICT+1)=sign_conflict
         FSIdata3D(i,j,k,ibase+FSI_EXTRAP_FLAG+1)=mask_local
         do dir=1,3
          FSIdata3D(i,j,k,ibase+FSI_VELOCITY+dir)=vel_local(dir)
         enddo 
          ! the Eulerian force will be corrected later so that:
          ! a) integral 1 dA = integral delta dV
          !  delta_cor=
          !   delta * (integral 1 dA)/
          !           (integral 1 delta dV)
          ! b) integral F_lag dA=integral delta_cor F_cor dV
          !    F_cor=F+c
         do dir=1,NCOMP_FORCE_STRESS
          FSIdata3D(i,j,k,ibase+FSI_FORCE+dir)=force_local(dir)*dt
         enddo

          ! FIX ME: delta function must be corrected so that
          ! perimeter measured from Eulerian and Lagrangian perspectives
          ! match.  (total forces should match too?)
         FSIdata3D(i,j,k,ibase+FSI_AREA_PER_VOL+1)= &
                 hsprime(ls_local,FSI_delta_cutoff)
         FSIdata3D(i,j,k,ibase+FSI_TEMPERATURE+1)=temp_local

        else if ((i.ge.FSI_growlo3D(1)).and.(i.le.FSI_growhi3D(1)).and. &
                 (j.ge.FSI_growlo3D(2)).and.(j.le.FSI_growhi3D(2)).and. &
                 (k.ge.FSI_growlo3D(3)).and.(k.le.FSI_growhi3D(3))) then
         ! do nothing
        else
         print *,"i,j,k outside of FSI_growlo3D,FSI_growhi3D range"
         stop
        endif 

       enddo
       enddo
       enddo ! i,j,k=gridloBB..gridhiBB

      else if (null_intersection.eq.1) then
       ! do nothing
      else
       print *,"null_intersection invalid"
       stop
      endif

     else if ((test_scale.ge.zero).and.(test_scale.le.one-VOFTOL)) then
      ! do nothing
     else
      print *,"test_scale invalid"
      stop
     endif

    enddo ! ielem_container=1,num_elements_container

    do dir=1,3
     if (dx3D(dir).gt.zero) then
      ! do nothing
     else
      print *,"dx3D(dir).le.zero"
      stop
     endif
    enddo ! dir=1..3

    do i=FSI_lo3D(1),FSI_hi3D(1)
    do j=FSI_lo3D(2),FSI_hi3D(2)
    do k=FSI_lo3D(3),FSI_hi3D(3)

     mask1=NINT(masknbr3D(i,j,k,1))
     mask2=NINT(masknbr3D(i,j,k,2))

     if ((mask1.eq.1).and. &
         (mask2.eq.1)) then
      ! do nothing
     else
      print *,"expecting mask1=mask2=1"
      print *,"mask1: ",mask1
      print *,"mask2: ",mask2
      stop
     endif

     ls_local=FSIdata3D(i,j,k,ibase+FSI_LEVELSET+1)
     mask_local=NINT(FSIdata3D(i,j,k,ibase+FSI_EXTRAP_FLAG+1))

     if (mask_local.eq.FSI_NOTHING_VALID) then
      ls_local=8.0*dx3D(1)
     else if (mask_local.eq.FSI_COARSE_LS_SIGN_VEL_VALID ) then
      ! do nothing, just use the value from fill coarse patch
     else if (mask_local.eq.FSI_FINE_SIGN_VEL_VALID) then 
      print *,"mask_local.eq.FSI_FINE_SIGN_VEL_VALID invalid here"
      stop
     else if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID) then 
      print *,"mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID invalid here"
      stop
     else if ((mask_local.eq.FSI_FINE_VEL_VALID ).or. &
              (mask_local.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID).or. &
              (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID)) then
      ! do nothing
     else
      print *,"mask_local invalid"
      stop
     endif

     if ((ctml_part_id.ge.1).and. &
         (ctml_part_id.le.CTML_NPARTS)) then
      if (FSI(part_id)%exclusive_doubly_wetted.eq.1) then
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

     if (FSI_mesh_type%exclusive_doubly_wetted.eq.0) then

      ! do nothing

     else if (FSI_mesh_type%exclusive_doubly_wetted.eq.1) then

      if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID) then
       print *,"mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID invalid here"
       stop
      else if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID) then
       ! do nothing
      else if (mask_local.eq.FSI_FINE_SIGN_VEL_VALID) then 
       print *,"mask_local.eq.FSI_FINE_SIGN_VEL_VALID invalid here"
       stop
      else if ((mask_local.eq.FSI_COARSE_LS_SIGN_VEL_VALID).or. &
               (mask_local.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID).or. &
               (mask_local.eq.FSI_NOTHING_VALID).or. &
               (mask_local.eq.FSI_FINE_VEL_VALID)) then
       mask_local=FSI_DOUBLY_WETTED_SIGN_VEL_VALID
      else
       print *,"mask_local invalid"
       stop
      endif

      ls_local=-abs(ls_local)

     else
      print *,"exclusive_doubly_wetted invalid"
      stop
     endif

     FSIdata3D(i,j,k,ibase+FSI_LEVELSET+1)=ls_local
     FSIdata3D(i,j,k,ibase+FSI_EXTRAP_FLAG+1)=mask_local

    enddo   
    enddo   
    enddo   

    do i=FSI_lo3D(1),FSI_hi3D(1)
    do j=FSI_lo3D(2),FSI_hi3D(2)
    do k=FSI_lo3D(3),FSI_hi3D(3)

     mask1=NINT(masknbr3D(i,j,k,1))
     mask2=NINT(masknbr3D(i,j,k,2))

     if ((mask1.eq.1).and. &
         (mask2.eq.1)) then
      ! do nothing
     else
      print *,"expecting mask1=mask2=1"
      print *,"mask1: ",mask1
      print *,"mask2: ",mask2
      stop
     endif

     ls_local=FSIdata3D(i,j,k,ibase+FSI_LEVELSET+1)
     mask_local=NINT(FSIdata3D(i,j,k,ibase+FSI_EXTRAP_FLAG+1))

     if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID) then
       !ls_local<0 in the fluid
      if (ls_local.le.zero) then
       ls_local=ls_local+dx3D(1)
      else
       print *,"ls_local invalid"
       stop
      endif
     else if (mask_local.eq.FSI_DOUBLY_WETTED_SIGN_VEL_VALID) then 
       !ls_local<0 in the fluid
      if (ls_local.lt.zero) then
       ! do nothing
      else
       print *,"ls_local invalid"
       stop
      endif
     else if (mask_local.eq.FSI_FINE_SIGN_VEL_VALID) then 
      print *,"mask_local.eq.FSI_FINE_SIGN_VEL_VALID invalid here"
      stop
     else if ((mask_local.eq.FSI_NOTHING_VALID).or. &
              (mask_local.eq.FSI_COARSE_LS_SIGN_VEL_VALID).or. &
              (mask_local.eq.FSI_FINE_VEL_VALID).or. &
              (mask_local.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID)) then
      ! do nothing
     else
      print *,"mask_local invalid"
      stop
     endif

     FSIdata3D(i,j,k,ibase+FSI_LEVELSET+1)=ls_local

    enddo
    enddo
    enddo

   else if (FSI_operation.eq.OP_FSI_MAKE_SIGN) then

     !FSI_growlo3D(dir)=FSI_lo3D(dir)-ngrow_make_distance_in
     !FSI_growhi3D(dir)=FSI_hi3D(dir)+ngrow_make_distance_in
    allocate(old_FSIdata( &
         FSI_growlo3D(1):FSI_growhi3D(1), &
         FSI_growlo3D(2):FSI_growhi3D(2), &
         FSI_growlo3D(3):FSI_growhi3D(3), &
         nFSI))

    old_FSIdata_ptr=>old_FSIdata

    call checkbound3D_array(FSI_lo3D,FSI_hi3D, &
     old_FSIdata_ptr, &
     ngrow_make_distance_in,-1)

    do i=FSI_growlo3D(1),FSI_growhi3D(1)
    do j=FSI_growlo3D(2),FSI_growhi3D(2)
    do k=FSI_growlo3D(3),FSI_growhi3D(3)
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

    if (iter.ge.0) then
     ! do nothing
    else
     print *,"iter invalid"
     stop
    endif

    do i=growlo3D(1),growhi3D(1)
    do j=growlo3D(2),growhi3D(2)
    do k=growlo3D(3),growhi3D(3)

     mask1=NINT(masknbr3D(i,j,k,1))
     mask2=NINT(masknbr3D(i,j,k,2))

     ! masknbr3D derived from localMF[MASK_NBR_MF]
     ! (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
     ! (2) =1 interior  =0 otherwise
     ! (3) =1 interior+ngrow-1  =0 otherwise
     ! (4) =1 interior+ngrow    =0 otherwise
     if ((mask1.eq.0).or. &
         (mask2.eq.1)) then

      mask_local=NINT(old_FSIdata(i,j,k,ibase+FSI_EXTRAP_FLAG+1))
      ls_local=old_FSIdata(i,j,k,ibase+FSI_LEVELSET+1)
      ! =0 no hits
      ! =1 positive
      ! =2 negative
      ! =3 inconclusive
      sign_conflict_local=old_FSIdata(i,j,k,ibase+FSI_SIGN_CONFLICT+1)

      new_mask_local=mask_local

      do dir=1,3
       xcen(dir)=xdata3D(i,j,k,dir)
      enddo

      if (sign_valid(mask_local).eq.0) then

       sign_status_changed=0

       call SUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP( &
         FSI_mesh_type%exterior_BB, &
         FSI_mesh_type%interior_BB, &
         xcen,time, &
         override_LS, &
         override_VEL, &
         override_TEMP, &
         override_MASK, &
         lev77, & !lev77=-1 for aux,>=1 otherwise
         im_part, &
         part_id)

       if (sign_valid(override_MASK).eq.1) then

        ! induces "new_mask_local=FSI_FINE_SIGN_VEL_VALID" below.
        sign_status_changed=1 
        if (sign_conflict_local.eq.zero) then
         ls_local=override_LS
        else if ((sign_conflict_local.eq.one).or. &
                 (sign_conflict_local.eq.two).or. &
                 (sign_conflict_local.eq.three)) then
         if (override_LS.lt.zero) then
          ls_local=-abs(ls_local)
         else if (override_LS.gt.zero) then
          ls_local=abs(ls_local)
         else if (override_LS.eq.zero) then
          ls_local=zero
         else
          print *,"override_LS invalid"
          stop
         endif
        else
         print *,"sign_conflict_local invalid"
         stop
        endif
        do dir=1,3
         FSIdata3D(i,j,k,ibase+FSI_VELOCITY+dir)=override_VEL(dir)
        enddo
        FSIdata3D(i,j,k,ibase+FSI_TEMPERATURE+1)=override_TEMP

       else if (sign_valid(override_MASK).eq.0) then

        inner_band_size=1.0d0

        if (((sign_conflict_local.eq.one).or. &
             (sign_conflict_local.eq.two)).and. &
            (abs(ls_local).le.inner_band_size*dx3D(1))) then
         ! induces "new_mask_local=FSI_FINE_SIGN_VEL_VALID" below.
         local_smooth_count=local_smooth_count+1
         sign_status_changed=1
         if (sign_conflict_local.eq.one) then
          ls_local=-abs(ls_local)
         else if (sign_conflict_local.eq.two) then
          ls_local=abs(ls_local)
         else
          print *,"sign_conflict_local invalid"
          stop
         endif
        else if ((sign_conflict_local.eq.zero).and. &
                 ((mask_local.eq.FSI_COARSE_LS_SIGN_VEL_VALID).or. &
                  (mask_local.eq.FSI_COARSE_LS_SIGN_FINE_VEL_VALID))) then
         ! induces "new_mask_local=FSI_FINE_SIGN_VEL_VALID" below.
         sign_status_changed=1
        else if ((sign_conflict_local.eq.zero).or. &
                 (sign_conflict_local.eq.three).or. &
                 (abs(ls_local).ge.inner_band_size*dx3D(1))) then

         if ((sign_conflict_local.eq.three).and. &
             (abs(ls_local).le.inner_band_size*dx3D(1))) then

          call SUB_GET_OUTSIDE_POINT( &
            FSI_mesh_type%exterior_BB, &
            xcen,time, &
            x_outside, &
            im_part, &
            part_id)

           ! x(t)=xcen + t xright
          do dir=1,3
           xright(dir)=x_outside(dir)-xcen(dir)
          enddo
          num_plane_intersects=0

          do ielem_container=1,num_elements_container

           if (lev77.eq.-1) then
            ielem=ielem_container
           else if (lev77.ge.1) then           
            ielem=contain_elem(lev77)% &
             level_elem_data(tid+1,part_id,tilenum+1)% &
             ElemData(ielem_container)
           else
            print *,"lev77 invalid"
            stop
           endif

           if ((ielem.lt.1).or. &
               (ielem.gt.num_elements)) then
            print *,"ielem invalid"
            stop
           endif

           do dir=1,3
            xelem(dir)=FSI_mesh_type%ElemDataXnotBIG(dir,ielem)
           enddo 
           call get_target_from_foot(xelem,xnot, &
            velparm,time, &
            FSI_mesh_type, &
            part_id, &
            nparts)

           nodes_per_elem=FSI_mesh_type%ElemDataBIG(1,ielem)
           if (nodes_per_elem.lt.3) then
            print *,"elem,nodes_per_elem ",ielem,nodes_per_elem   
            stop
           endif
           call scinormalBIG(ielem, &
            normal, &
            FSI_mesh_type, &
            part_id, &
            nparts, &
            time)

           test_scale=sqrt(normal(1)**2+normal(2)**2+normal(3)**2)

           if (abs(test_scale-one).le.VOFTOL) then

            eul_over_lag_scale=zero
            t_top=zero
            t_bottom=zero
            far_field_sign=zero
            near_field_sign=zero
            do dir=1,3
             far_field_sign=far_field_sign- &
                normal(dir)*(x_outside(dir)-xnot(dir))
             near_field_sign=near_field_sign- &
                normal(dir)*(xcen(dir)-xnot(dir))
             t_top=t_top-normal(dir)*(xcen(dir)-xnot(dir))
             t_bottom=t_bottom+normal(dir)*xright(dir)
            enddo
            if (abs(t_bottom).gt.zero) then
             t_crit=t_top/t_bottom
             if ((t_crit.ge.zero).and.(t_crit.le.one)) then

              do dir=1,3
               xclosest(dir)=xcen(dir)+t_crit*xright(dir)
              enddo 

              call checkinplaneBIG( &
               eul_over_lag_scale, &
               xcen, & ! not used
               xclosest, &
               xclosest_project, & !INTENT(out)
               normal, & ! INTENT(in)
               normal_closest, & ! INTENT(out)
               ielem, &
               element_node_edge_inplane, & !INTENT(out)
               FSI_mesh_type, &
               part_id, &
               nparts, &
               time)

              if (element_node_edge_inplane.eq.1) then
               num_plane_intersects=num_plane_intersects+1
               if ((num_plane_intersects.ge.1).and. &
                   (num_plane_intersects.le.max_plane_intersects)) then
                if (far_field_sign*near_field_sign.lt.zero) then
                 plane_intersect_list(num_plane_intersects)=t_crit
                 if (far_field_sign.lt.zero) then
                  plane_intersect_list_sign(num_plane_intersects)=-1
                 else if (far_field_sign.gt.zero) then
                  plane_intersect_list_sign(num_plane_intersects)=1
                 else
                  print *,"far_field_sign invalid"
                  stop
                 endif
                else
                 print *,"far_field_sign or near_field_sign invalid"
                 stop
                endif
               else
                print *,"num_plane_intersects invalid"
                stop
               endif
              else if (element_node_edge_inplane.eq.0) then
               ! do nothing
              else
               print *,"element_node_edge_inplane invalid"
               stop
              endif
             else if ((t_crit.lt.zero).or.(t_crit.gt.one)) then
              ! do nothing
             else
              print *,"t_crit is NaN"
              stop
             endif
            else if (abs(t_bottom).eq.zero) then
             ! do nothing
            else
             print *,"t_bottom is NaN"
             stop
            endif

           else if ((test_scale.ge.zero).and. &
                    (test_scale.le.one-VOFTOL)) then
            ! do nothing
           else
            print *,"test_scale invalid"
            stop
           endif

          enddo ! ielem_container=1,num_elements_container

           ! sort from highest to lowest
          do ii=1,num_plane_intersects-1
           do jj=1,num_plane_intersects-ii
             ! x(t)=xcen + t xright
            less_than_flag=is_less_than_list( &
             plane_intersect_list(jj), &
             plane_intersect_list_sign(jj), &
             plane_intersect_list(jj+1), &
             plane_intersect_list_sign(jj+1))

            if (less_than_flag.eq.1) then
             swap_data=plane_intersect_list(jj+1)
             swap_sign=plane_intersect_list_sign(jj+1)
             plane_intersect_list(jj+1)=plane_intersect_list(jj)
             plane_intersect_list_sign(jj+1)=plane_intersect_list_sign(jj)
             plane_intersect_list(jj)=swap_data
             plane_intersect_list_sign(jj)=swap_sign
            else if (less_than_flag.eq.0) then
             ! do nothing
            else
             print *,"less_than_flag invalid"
             stop
            endif

           enddo
          enddo

          num_plane_intersects_new=0
          cur_ptr=1
          do while (cur_ptr.le.num_plane_intersects)
           if (cur_ptr+1.gt.num_plane_intersects) then
            num_plane_intersects_new=num_plane_intersects_new+1
            plane_intersect_list(num_plane_intersects_new)= &
               plane_intersect_list(cur_ptr)
            plane_intersect_list_sign(num_plane_intersects_new)= &
               plane_intersect_list_sign(cur_ptr)
            cur_ptr=cur_ptr+1
           else if (cur_ptr+1.le.num_plane_intersects) then
            equal_flag=is_equal_list( &
             plane_intersect_list(cur_ptr), &
             plane_intersect_list_sign(cur_ptr), &
             plane_intersect_list(cur_ptr+1), &
             plane_intersect_list_sign(cur_ptr+1))
            if (equal_flag.eq.0) then
             num_plane_intersects_new=num_plane_intersects_new+1
             plane_intersect_list(num_plane_intersects_new)= &
               plane_intersect_list(cur_ptr)
             plane_intersect_list_sign(num_plane_intersects_new)= &
               plane_intersect_list_sign(cur_ptr)
             cur_ptr=cur_ptr+1
            else if (equal_flag.eq.1) then
             cur_ptr=cur_ptr+1
            else
             print *,"equal_flag invalid"
             stop
            endif
           else
            print *,"cur_ptr invalid"
            stop
           endif
          enddo !while (cur_ptr.le.num_plane_intersects)

          num_plane_intersects=num_plane_intersects_new

           ! list sorted from highest to lowest
          ls_local=-abs(ls_local)
          cur_sign=-1
          cur_ptr=1
          do while (cur_ptr.le.num_plane_intersects)
           if (cur_ptr+1.gt.num_plane_intersects) then
            if (plane_intersect_list_sign(cur_ptr).eq.cur_sign) then
             ! do nothing
            else
             print *,"inconsistent sign(1)"
             print *,"im_part, part_id ",im_part,part_id
             print *,"num_plane_intersects ",num_plane_intersects
             print *,"xcen ",xcen(1),xcen(2),xcen(3)
             print *,"x_outside ",x_outside(1),x_outside(2),x_outside(3)
             do ii=1,num_plane_intersects
              print *,"ii,t,(outside) sign ",ii, &
                plane_intersect_list(ii), &
                plane_intersect_list_sign(ii)
             enddo
             ls_local=-ls_local
             cur_sign=-cur_sign
             stop
            endif

            cur_ptr=cur_ptr+1
            ls_local=-ls_local
            cur_sign=-cur_sign
           else if (cur_ptr+1.le.num_plane_intersects) then
            plane_diff=abs(plane_intersect_list(cur_ptr)- &
                           plane_intersect_list(cur_ptr+1))
            if (plane_diff.le.crossing_tol) then
             if (plane_intersect_list_sign(cur_ptr)* &
                 plane_intersect_list_sign(cur_ptr+1).eq.-1) then
              ! do nothing
             else
              print *,"not all duplicates deleted"
              stop
             endif
             cur_ptr=cur_ptr+2
            else if (plane_diff.gt.crossing_tol) then
             if (plane_intersect_list_sign(cur_ptr).eq.cur_sign) then
              ! do nothing
             else
              print *,"inconsistent sign(2)"
              print *,"cur_ptr,cur_sign,plane_diff ",cur_ptr,cur_sign,plane_diff

              print *,"im_part, part_id ",im_part,part_id
              print *,"num_plane_intersects ",num_plane_intersects
              print *,"xcen ",xcen(1),xcen(2),xcen(3)
              print *,"x_outside ",x_outside(1),x_outside(2),x_outside(3)
              do ii=1,num_plane_intersects
               print *,"ii,t,(outside) sign ",ii, &
                 plane_intersect_list(ii), &
                 plane_intersect_list_sign(ii)
              enddo
              ls_local=-ls_local
              cur_sign=-cur_sign
              stop
             endif

             cur_ptr=cur_ptr+1
             ls_local=-ls_local
             cur_sign=-cur_sign
            else
             print *,"plane_diff is NaN"
             stop
            endif
           else
            print *,"cur_ptr invalid"
            stop
           endif
          enddo !while (cur_ptr.le.num_plane_intersects)

          local_corner_count=local_corner_count+1
          if (ioproc.eq.1) then
           print *,"local_corner_count,part_id,xcen,num_sign_changes ", &
             local_corner_count,part_id,xcen(1),xcen(2),xcen(3), &
             num_sign_changes
          endif

          sign_status_changed=1

         else if ((sign_conflict_local.eq.zero).or. &
                  (sign_conflict_local.eq.one).or. &
                  (sign_conflict_local.eq.two).or. &
                  (abs(ls_local).ge.inner_band_size*dx3D(1))) then

          LS_sum=zero
          total_variation_sum_plus=zero
          total_variation_sum_minus=zero
          do dir=1,NCOMP_FSI
           weight_top(dir)=zero
          enddo
          weight_bot=zero
          weight_total_variation=0
      
          do ii=-1,1
          do jj=-1,1
          do kk=-1,1
           if ((i+ii.ge.FSI_growlo3D(1)).and.(i+ii.le.FSI_growhi3D(1)).and. &
               (j+jj.ge.FSI_growlo3D(2)).and.(j+jj.le.FSI_growhi3D(2)).and. &
               (k+kk.ge.FSI_growlo3D(3)).and.(k+kk.le.FSI_growhi3D(3))) then

            mask_node=NINT(old_FSIdata(i+ii,j+jj,k+kk,ibase+FSI_EXTRAP_FLAG+1))

            weight=VOFTOL*dx3D(1)
            do dir=1,3
             xc(dir)=xdata3D(i+ii,j+jj,k+kk,dir)
             weight=weight+(xc(dir)-xcen(dir))**2
            enddo
            if (weight.gt.zero) then
             ! do nothing
            else
             print *,"weight invalid"
             stop
            endif
            weight=one/weight

            if (vel_valid(mask_node).eq.1) then
             do dir=1,3
              weight_top(FSI_VELOCITY+dir)=weight_top(FSI_VELOCITY+dir)+ &
               old_FSIdata(i+ii,j+jj,k+kk,ibase+FSI_VELOCITY+dir)*weight
             enddo
              ! temperature
             weight_top(FSI_TEMPERATURE+1)=weight_top(FSI_TEMPERATURE+1)+ &
               old_FSIdata(i+ii,j+jj,k+kk,ibase+FSI_TEMPERATURE+1)*weight
       
             weight_bot=weight_bot+weight
            else if (vel_valid(mask_node).eq.0) then
             ! do nothing
            else
             print *,"vel_valid(mask_node) invalid"
             stop
            endif

            ! sign_valid==1 if mask=
            !   FSI_FINE_SIGN_VEL_VALID, 
            !   FSI_DOUBLY_WETTED_SIGN_VEL_VALID, 
            !   FSI_DOUBLY_WETTED_SIGN_LS_VEL_VALID.
            ! sign_valid==0 if 
            !   mask=FSI_NOTHING_VALID,
            !   FSI_FINE_VEL_VALID,
            !   FSI_COARSE_LS_SIGN_VEL_VALID, 
            !   FSI_COARSE_LS_SIGN_FINE_VEL_VALID
            if (sign_valid(mask_node).eq.1) then

             weight_total_variation=weight_total_variation+1

              !ibase=(part_id-1)*NCOMP_FSI
             LS_SIDE=old_FSIdata(i+ii,j+jj,k+kk,ibase+FSI_LEVELSET+1)
             dx_SIDE=zero
             do dir=1,3
              dx_SIDE=dx_SIDE+(xc(dir)-xcen(dir))**2
             enddo
             dx_SIDE=sqrt(dx_SIDE)

             LS_sum=LS_sum+LS_SIDE

             if (dx_SIDE.gt.zero) then
              total_variation_sum_plus=total_variation_sum_plus+ &
                 (abs(LS_SIDE-abs(ls_local)))/dx_SIDE
              total_variation_sum_minus=total_variation_sum_minus+ &
                 (abs(LS_SIDE+abs(ls_local)))/dx_SIDE
             else
              print *,"dx_SIDE invalid"
              stop
             endif

            else if (sign_valid(mask_node).eq.0) then
             ! do nothing
            else
             print *,"sign_valid(mask_node) invalid"
             stop
            endif
           endif ! ii,jj,kk in grid
          enddo
          enddo
          enddo ! ii,jj,kk

          if ((weight_total_variation.ge.1).and. &
              (weight_total_variation.le.3*3*3-1)) then
           if (sign_conflict_local.eq.zero) then
            ls_local=LS_sum/weight_total_variation
           else if ((sign_conflict_local.eq.one).or. &
                    (sign_conflict_local.eq.two).or. & 
                    (sign_conflict_local.eq.three)) then

            if ((total_variation_sum_plus.ge.zero).and. &
                (total_variation_sum_minus.ge.zero)) then 
             if (total_variation_sum_plus.le. &
                 total_variation_sum_minus) then
              ls_local=abs(ls_local)
             else if (total_variation_sum_plus.ge. &
                      total_variation_sum_minus) then
              ls_local=-abs(ls_local)
             else
              print *,"total_variation_sum plus or minus invalid(1)"
              stop
             endif
            else
             print *,"total_variation_sum plus or minus invalid(2)"
             stop
            endif
           else
            print *,"sign_conflict_local invalid"
            stop
           endif
           ! induces "new_mask_local=FSI_FINE_SIGN_VEL_VALID" below.
           sign_status_changed=1
          else if (weight_total_variation.eq.0) then
           ! do nothing
          else
           print *,"weight_total_variation invalid: ",weight_total_variation
           stop
          endif

          if (vel_valid(mask_local).eq.0) then 

           if (weight_bot.gt.zero) then
            do dir=1,3
             FSIdata3D(i,j,k,ibase+FSI_VELOCITY+dir)= &
                    weight_top(FSI_VELOCITY+dir)/weight_bot
            enddo
            FSIdata3D(i,j,k,ibase+FSI_TEMPERATURE+1)= &
                   weight_top(FSI_TEMPERATURE+1)/weight_bot

           else if (weight_bot.eq.zero) then
            ! do nothing
           else
            print *,"weight_bot is NaN"
            stop
           endif 

          else if (vel_valid(mask_local).eq.1) then
           ! do nothing
          else
           print *,"vel_valid(mask_local) invalid"
           stop
          endif

         else
          print *,"sign_conflict_local or ls_local invalid (1) "
          stop
         endif

        else
         print *,"sign_conflict_local or ls_local invalid (2) "
         stop
        endif
       else
        print *,"sign_valid(override_MASK) invalid"
        stop
       endif

       if (sign_status_changed.eq.1) then 
        new_mask_local=FSI_FINE_SIGN_VEL_VALID
        touch_flag=1
        numtouch=numtouch+1
       else if (sign_status_changed.eq.0) then 
        ! do nothing
       else
        print *,"sign_status_changed invalid"
        stop
       endif

      else if (sign_valid(mask_local).eq.1) then
       ! do nothing
      else
       print *,"mask_local invalid"
       stop
      endif

      FSIdata3D(i,j,k,ibase+FSI_EXTRAP_FLAG+1)=new_mask_local 
      FSIdata3D(i,j,k,ibase+FSI_LEVELSET+1)=ls_local 

     else if ((mask1.eq.1).and.(mask2.eq.0)) then
      ! do nothing
     else
      print *,"mask1 or mask2 invalid"
      stop
     endif

    enddo 
    enddo 
    enddo  ! i,j,k

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
    print *,"END: CLSVOF_InitBox"
    print *,"FSI_operation=",FSI_operation
    if ((local_corner_count.gt.0).or. &
        (local_smooth_count.gt.0)) then
     print *,"local_corner_count=",local_corner_count
     print *,"local_smooth_count=",local_smooth_count
    endif
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

        ! called from fort_headermsg when FSI_operation.eq.OP_FSI_LAG_STRESS and
        ! FSI_sub_operation.eq.1.
        ! This routine is called only for those blocks associated with a 
        ! given node.
        ! CLSVOF_sync_lag_data must be called later in order to synchronize
        ! all of the Lagrangian data across all of the nodes.
      subroutine CLSVOF_Copy_To_LAG(  &
       sdim_AMR, &
       lev77, &
       tid, &
       tilenum, &
       im_part, & ! 1..num_materials
       nparts, &
       part_id, &
       ngrow_make_distance_in, &
       nFSI, &
       FSI_operation, &
       time, &
       problo3D,probhi3D, &
       xmap3D, &
       dx3D, &
       xlo3D_tile, &
       xhi3D_tile, &
       FSI_lo,FSI_hi, &
       FSI_growlo,FSI_growhi, &
       growlo3D,growhi3D, &
       xdata3D, &
       veldata3D, &
       stressdata3D, &
       stressflag3D, &
       masknbr3D, &
       maskfiner3D, &
       ioproc,isout)
       use global_utility_module
#ifdef MVAHABFSI
       use CTML_module
#endif

       IMPLICIT NONE

      INTEGER_T, INTENT(in) :: sdim_AMR
      INTEGER_T, INTENT(in) :: lev77
      INTEGER_T, INTENT(in) :: tid
      INTEGER_T, INTENT(in) :: tilenum
      INTEGER_T, INTENT(in) :: im_part ! 1..num_materials
      INTEGER_T, INTENT(in) :: nparts
      INTEGER_T, INTENT(in) :: part_id
      INTEGER_T, INTENT(in) :: ngrow_make_distance_in
      INTEGER_T, INTENT(in) :: nFSI
      INTEGER_T, INTENT(in) :: FSI_operation
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(in) :: problo3D(3)
      REAL_T, INTENT(in) :: probhi3D(3)
      INTEGER_T, INTENT(in) :: xmap3D(3)
      REAL_T, INTENT(in) :: dx3D(3)
      REAL_T, INTENT(in) :: xlo3D_tile(3)
      REAL_T, INTENT(in) :: xhi3D_tile(3)
      INTEGER_T, INTENT(in) :: FSI_lo(3),FSI_hi(3)
      INTEGER_T, INTENT(in) :: FSI_growlo(3),FSI_growhi(3)
      INTEGER_T, INTENT(in) :: growlo3D(3),growhi3D(3)
      REAL_T, INTENT(in), pointer :: xdata3D(:,:,:,:)
      REAL_T, INTENT(in), pointer :: veldata3D(:,:,:,:)
      REAL_T, INTENT(in), pointer :: stressdata3D(:,:,:,:)
      REAL_T, INTENT(in), pointer :: stressflag3D(:,:,:,:)
      REAL_T, INTENT(in), pointer :: masknbr3D(:,:,:,:)
      REAL_T, INTENT(in), pointer :: maskfiner3D(:,:,:,:)

      INTEGER_T, INTENT(in) :: ioproc,isout

      REAL_T dxBB(3) ! set in find_grid_bounding_box_node
      REAL_T dxBB_probe(3) ! set in find_grid_bounding_box_node

      INTEGER_T :: inode
      INTEGER_T :: inode_container
      REAL_T, dimension(3) :: xnot
      REAL_T, dimension(3) :: xprobe
      REAL_T, dimension(3) :: xnode
      INTEGER_T, dimension(3) :: gridloBB,gridhiBB
      INTEGER_T, dimension(3) :: gridloBB_probe,gridhiBB_probe
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
      INTEGER_T inside_interior_box
      REAL_T dx3D_min
      REAL_T probe_size
      REAL_T null_probe_size
      REAL_T sign_normal
      REAL_T data_out
      REAL_T total_stress(6)
      REAL_T NodeStress(6)
      REAL_T local_node_normal(3)
      INTEGER_T idoubly
      INTEGER_T istress,jstress
      REAL_T stress_3d(3,3)
      REAL_T local_force(3)

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

      dx3D_min=dx3D(1)
      do dir=1,3
       if (dx3D(dir).lt.dx3D_min) then
        dx3D_min=dx3D(dir)
       endif
      enddo
      probe_size=one
      null_probe_size=zero

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

      ctml_part_id=ctml_part_id_map(part_id)
      fsi_part_id=fsi_part_id_map(part_id)

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
       if ((num_materials.lt.1).or.(num_materials.gt.50)) then
        print *,"num_materials out of range"
        stop
       endif
       if ((im_part.lt.1).or.(im_part.gt.num_materials)) then
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
        print *,"num_materials=",num_materials
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

       if (nFSI.ne.nparts*NCOMP_FSI) then
        print *,"nFSI invalid"
        stop
       endif
       if ((nparts.lt.1).or.(nparts.ge.num_materials)) then
        print *,"nparts invalid"
        stop
       endif

       if (FSI_operation.ne.4) then
        print *,"FSI_operation invalid"
        stop
       endif
       if (ngrow_make_distance_in.ne.3) then
        print *,"ngrow_make_distance_in invalid"
        stop
       endif

       call checkbound3D_array(FSI_lo,FSI_hi, &
        xdata3D, &
        ngrow_make_distance_in,-1)
       call checkbound3D_array(FSI_lo,FSI_hi, &
        veldata3D, &
        ngrow_make_distance_in,-1)
       call checkbound3D_array(FSI_lo,FSI_hi, &
        stressdata3D, &
        ngrow_make_distance_in,-1)
       call checkbound3D_array(FSI_lo,FSI_hi, &
        stressflag3D, &
        ngrow_make_distance_in,-1)
       call checkbound3D_array(FSI_lo,FSI_hi, &
        masknbr3D, &
        ngrow_make_distance_in,-1)
       call checkbound3D_array(FSI_lo,FSI_hi, &
        maskfiner3D, &
        ngrow_make_distance_in,-1)

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
         print *,"inode invalid10: inode=",inode
         stop
        endif

        do dir=1,3
         xnode(dir)=FSI(part_id)%Node(dir,inode)
         velparm(dir)=zero
        enddo 
        call get_target_from_foot(xnode,xnot, &
          velparm,time, &
          FSI(part_id), &
          part_id, &
          TOTAL_NPARTS)

        if (debug_all.eq.1) then
         print *,"inode=",inode
         do dir=1,3
          print *,"dir,xnode ",dir,xnode(dir)
          print *,"dir,xnot  ",dir,xnot(dir)
         enddo
        endif

        call find_grid_bounding_box_node( &
         ngrow_make_distance_in, &
         null_probe_size, &
         xnot, &
         FSI_lo,FSI_hi, &
         FSI_growlo,FSI_growhi, &
         xdata3D, &
         gridloBB,gridhiBB,dxBB) 

        do dir=1,3
         if (abs(dxBB(dir)-dx3D(dir)).le.VOFTOL*dxBB(dir)) then
          ! do nothing
         else
          print *,"abs(dxBB(dir)-dx3D(dir)).gt.VOFTOL*dxBB(dir)"
          stop
         endif
        enddo ! dir=1..3

        if (1.eq.0) then
         print *,"inode=",inode
         do dir=1,3
          print *,"dir,gridloBB,gridhiBB ",dir,gridloBB(dir),gridhiBB(dir)
         enddo
        endif

        do dir=1,3
         idx(dir)=(gridhiBB(dir)+gridloBB(dir))/2
         if (gridhiBB(dir)-gridloBB(dir).ne.2*BoundingBoxRadNode) then
          print *,"gridhiBB(dir)-gridloBB(dir).ne.2*BoundingBoxRad"
          stop
         endif
        enddo

        local_mask=NINT(maskfiner3D(idx(1),idx(2),idx(3),1))

        if (local_mask.eq.1) then

         inside_interior_box=1
         do dir=1,3
          if (sdim_AMR.eq.3) then
           if ((xnot(dir).lt.xhi3D_tile(dir)).and. &
               (xnot(dir).ge.xlo3D_tile(dir))) then
            ! do nothing
           else if ((xnot(dir).ge.xhi3D_tile(dir)).or. &
                    (xnot(dir).lt.xlo3D_tile(dir))) then
            inside_interior_box=0
           else
            print *,"xnot,xlo3d, or xhibc NaN"
            stop
           endif
          else if (sdim_AMR.eq.2) then
           if (xmap3D(dir).eq.0) then
            ! do nothing
           else if ((xmap3D(dir).eq.1).or. &
                    (xmap3D(dir).eq.2)) then
            if ((xnot(dir).lt.xhi3D_tile(dir)).and. &
                (xnot(dir).ge.xlo3D_tile(dir))) then
             ! do nothing
            else if ((xnot(dir).ge.xhi3D_tile(dir)).or. &
                     (xnot(dir).lt.xlo3D_tile(dir))) then
             inside_interior_box=0
            else
             print *,"xnot,xlo3d, or xhibc NaN"
             stop
            endif
           else
            print *,"xmap3D invalid"
            stop
           endif
          else
           print *,"sdim_AMR invalid"
           stop
          endif
         enddo ! dir=1,3

         if (inside_interior_box.eq.1) then

          do dir=1,3
           local_node_normal(dir)=FSI(part_id)%NodeNormal(dir,inode)
          enddo

          do dir=1,3
           FSI(part_id)%NodeForce(dir,inode)=zero
          enddo

          do idoubly=1,2

           if (idoubly.eq.1) then
            sign_normal=-one
           else if (idoubly.eq.2) then 
            sign_normal=one
           else 
            print *,"idoubly invalid"
            stop
           endif

           xprobe(dir)=xnot(dir)+ &
             sign_normal*probe_size*dx3D_min*local_node_normal(dir)

           call find_grid_bounding_box_node( &
            ngrow_make_distance_in, &
            probe_size, &
            xprobe, &
            FSI_lo,FSI_hi, &
            FSI_growlo,FSI_growhi, &
            xdata3D, &
            gridloBB_probe,gridhiBB_probe,dxBB_probe) 

           total_weight=zero
           do dir=1,6
            total_stress(dir)=zero
           enddo
           do i=gridloBB_probe(1),gridhiBB_probe(1)
           do j=gridloBB_probe(2),gridhiBB_probe(2)
           do k=gridloBB_probe(3),gridhiBB_probe(3)
            wt=one
            do dir=1,3
             call safe_data3D(i,j,k,dir,xdata3D,data_out)
             dist_scale=abs(data_out-xprobe(dir))/dxBB_probe(dir)
             support_size=two
             df=hsprime(dist_scale,support_size)
             if (sdim_AMR.eq.3) then
              wt=wt*df
             else if (sdim_AMR.eq.2) then
              if (xmap3D(dir).eq.0) then
               ! do nothing
              else if ((xmap3D(dir).eq.1).or. &
                       (xmap3D(dir).eq.2)) then
               wt=wt*df
              else
               print *,"xmap3D invalid"
               stop
              endif
             else
              print *,"sdim_AMR invalid"
              stop
             endif

            enddo ! dir=1..3

            call safe_data3D(i,j,k,im_part,stressflag3D,data_out)
            if (data_out.eq.zero) then
             wt=zero
            else if (data_out.eq.one) then ! drag initialized
             ! do nothing
            else if (data_out.eq.two) then ! drag extended
             ! do nothing
            else
             print *,"data_out (stressflag3D) invalid"
             stop
            endif

            total_weight=total_weight+wt
            do dir=1,6
             call safe_data3D(i,j,k,6*(im_part-1)+dir,stressdata3D,data_out)
             total_stress(dir)=total_stress(dir)+wt*data_out
            enddo
           enddo
           enddo
           enddo
           if (total_weight.gt.zero) then
            do dir=1,6
             NodeStress(dir)=total_stress(dir)/total_weight
             call stress_index(dir,istress,jstress)
             stress_3d(istress,jstress)=NodeStress(dir)
            enddo
            stress_3d(2,1)=stress_3d(1,2)
            stress_3d(3,1)=stress_3d(1,3)
            stress_3d(3,2)=stress_3d(2,3)
            do dir=1,3
             local_force(dir)=zero
            enddo
            do istress=1,3
             do jstress=1,3
              local_force(istress)=local_force(istress)+ &
                sign_normal*local_node_normal(jstress)* &
                stress_3d(istress,jstress) 
             enddo
            enddo

             !NodeForce used in CLSVOF_sync_lag_data 
            do dir=1,3
             FSI(part_id)%NodeForce(dir,inode)= &
               FSI(part_id)%NodeForce(dir,inode)+ &
               local_force(dir)
            enddo

           else
            print *,"total_weight invalid"
            stop
           endif

          enddo ! idoubly=1,2

          total_weight=zero
          do dir=1,3
           total_vel(dir)=zero
          enddo
          do i=gridloBB(1),gridhiBB(1)
          do j=gridloBB(2),gridhiBB(2)
          do k=gridloBB(3),gridhiBB(3)
           wt=one
           do dir=1,3
            call safe_data3D(i,j,k,dir,xdata3D,data_out)
            dist_scale=abs(data_out-xnot(dir))/dxBB(dir)
            support_size=two
            df=hsprime(dist_scale,support_size)
            if (sdim_AMR.eq.3) then
             wt=wt*df
            else if (sdim_AMR.eq.2) then
             if (xmap3D(dir).eq.0) then
              ! do nothing
             else if ((xmap3D(dir).eq.1).or. &
                      (xmap3D(dir).eq.2)) then
              wt=wt*df
             else
              print *,"xmap3D invalid"
              stop
             endif
            else
             print *,"sdim_AMR invalid"
             stop
            endif

           enddo ! dir=1..3

           total_weight=total_weight+wt
           do dir=1,3
            call safe_data3D(i,j,k,dir,veldata3D,data_out)
            total_vel(dir)=total_vel(dir)+wt*data_out
           enddo
          enddo
          enddo
          enddo
          if (total_weight.gt.zero) then
           !NodeVel used in CLSVOF_sync_lag_data 
           do dir=1,3
            FSI(part_id)%NodeVel(dir,inode)=total_vel(dir)/total_weight
           enddo
          else
           print *,"total_weight invalid"
           stop
          endif

         else if (inside_interior_box.eq.0) then
          ! do nothing
         else
          print *,"inside_interior_box invalid"
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
      INTEGER_T, INTENT(in) :: numMotion
      REAL_T, INTENT(in) :: motionPara(11,numMotion)
      REAL_T, INTENT(out) :: r(3,4)
      REAL_T, INTENT(in) :: t
      REAL_T xPoint(3,numMotion),vTan(3,numMotion)
      REAL_T x0(3),v(3)
      REAL_T vNorm,theta,thetaMag,fTheta,phiTheta,theta0
      REAL_T hMag,fH,phiH,h0,ct,st
      REAL_T r1(3,4),r2(3,4)
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
          r(1,3) = v(1)*v(3)*(1.-ct)+st*v(2)
          r(2,1) = v(1)*v(2)*(1.-ct)+st*v(3)  
          r(2,2) = v(2)*v(2)*(1.-ct)+ct
          r(2,3) = v(2)*v(3)*(1.-ct)-st*v(1)
          r(3,1) = v(1)*v(3)*(1.-ct)-st*v(2)  
          r(3,2) = v(3)*v(2)*(1.-ct)+st*v(1)  
          r(3,3) = v(3)*v(3)*(1.-ct)+ct     

          r(1,4) = (1.-r(1,1))*x0(1)    -r(1,2) *x0(2)    -r(1,3) *x0(3)
          r(2,4) =    -r(2,1) *x0(1)+(1.-r(2,2))*x0(2)    -r(2,3) *x0(3)
          r(3,4) =    -r(3,1) *x0(1)    -r(3,2) *x0(2)+(1.-r(3,3))*x0(3)

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

          r(1,4) = x0(1) + v(1)*theta
          r(2,4) = x0(2) + v(2)*theta
          r(3,4) = x0(3) + v(3)*theta

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
      if (FSI(part_id)%part_id.ne.part_id) then
       print *,"FSI(part_id)%part_id.ne.part_id"
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
        RPM=abs(vinletgas)
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

        if (levelrz.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else if (levelrz.eq.COORDSYS_RZ) then ! place inlet at center of domain
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
        if (levelrz.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else if (levelrz.eq.COORDSYS_RZ) then ! place inlet at center of domain
         xfoot(1)=xfoot(1)+0.0015
         xfoot(2)=xfoot(2)-0.0045
        else
         print *,"levelrz invalid get_foot_from_target 2"
         stop
        endif 

        xfoot(3)=xfoot(3)+0.01d0
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
         motionPara(2,2)=0.5d0 ! hMag (was 0.0)
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

         dt_flapping=0.01d0
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
       velparm,time, &
       FSI_mesh_type, &
       part_id, &
       max_part_id)
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: part_id
      INTEGER_T, INTENT(in) :: max_part_id
      type(mesh_type), INTENT(in) :: FSI_mesh_type
      REAL_T, INTENT(out) :: velparm(3)
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(out) :: xtarget(3)
      REAL_T, INTENT(inout) :: xfoot(3)
      INTEGER_T dir
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


      if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
       print *,"part_id invalid"
       stop
      endif
      if (FSI_mesh_type%part_id.ne.part_id) then
       print *,"FSI_mesh_type%part_id.ne.part_id"
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
        RPM=abs(vinletgas)
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

        if (levelrz.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else if (levelrz.eq.COORDSYS_RZ) then ! place inlet at center of domain
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
        if (levelrz.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else if (levelrz.eq.COORDSYS_RZ) then ! place inlet at center of domain
         xtarget(1)=xtarget(1)-0.0015
         xtarget(2)=xtarget(2)+0.0045
        else
         print *,"levelrz invalid get_target_from_foot 2"
         stop
        endif 

        xtarget(3)=xtarget(3)-0.01d0
        do dir=1,3
         xtarget(dir)=xtarget(dir)+FSI_mesh_type%solid_displ(dir)
        enddo

        do dir=1,3
         velparm(dir)=FSI_mesh_type%solid_speed(dir)
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
         motionPara(2,2)=0.5d0 ! hMag (was 0.0)
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

         dt_flapping=0.01d0
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

      else if ((part_id.gt.2).and.(part_id.le.max_part_id)) then
       ! do nothing
      else
       print *,"part_id invalid"
       stop
      endif

      return
      end subroutine get_target_from_foot


subroutine find_grid_bounding_box( &
 FSI_mesh_type, &
 part_id, &
 max_part_id, &
 null_intersection, &
 minnode,maxnode, &
 FSI_lo,FSI_hi, &
 FSI_growlo,FSI_growhi,  &
 growlo3D,growhi3D,  &
 xdata3D, &
 gridloBB,gridhiBB,dxBB)

use global_utility_module
IMPLICIT NONE

 type(mesh_type), INTENT(in) :: FSI_mesh_type
 INTEGER_T, INTENT(in) :: part_id
 INTEGER_T, INTENT(in) :: max_part_id
 INTEGER_T, INTENT(out) :: null_intersection
 REAL_T, INTENT(in) :: minnode(3),maxnode(3)
 INTEGER_T, INTENT(in) :: FSI_lo(3),FSI_hi(3)
 INTEGER_T, INTENT(in) :: FSI_growlo(3),FSI_growhi(3)
 INTEGER_T, INTENT(in) :: growlo3D(3),growhi3D(3)
 REAL_T, INTENT(in), pointer :: xdata3D(:,:,:,:)
 INTEGER_T, INTENT(out) :: gridloBB(3),gridhiBB(3)
 REAL_T, INTENT(out) :: dxBB(3)
 INTEGER_T dir
 INTEGER_T ii,jj,kk
 INTEGER_T i,j,k,incr,iter
 REAL_T xcontrol,xcost
 REAL_T xlo(3),xhi(3)
 INTEGER_T ngrow,ngrowtest
 INTEGER_T local_iband

 if ((part_id.lt.1).or.(part_id.gt.max_part_id)) then
  print *,"part_id invalid"
  stop
 endif

 local_iband=FSI_mesh_type%bounding_box_ngrow
 if (local_iband.ne.BoundingBoxRadCell) then
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

 call checkbound3D_array(FSI_lo,FSI_hi, &
  xdata3D, &
  ngrow,-1)

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

   gridloBB(dir)=NINT( (minnode(dir)-xlo(dir))/dxBB(dir)-half+FSI_lo(dir) )
   if (gridloBB(dir).lt.growlo3D(dir)) then
    gridloBB(dir)=growlo3D(dir)
   endif
   if (gridloBB(dir).gt.growhi3D(dir)) then
    gridloBB(dir)=growhi3D(dir)
    null_intersection=1
   endif

   if (null_intersection.eq.0) then

    incr=gridloBB(dir)-FSI_lo(dir)
    xcontrol=xdata3D(i+ii*incr,j+jj*incr,k+kk*incr,dir)
    xcost=minnode(dir)-local_iband*dxBB(dir)
    iter=0
    do while ((gridloBB(dir).gt.growlo3D(dir)).and. &
              (xcontrol.gt.xcost))
     gridloBB(dir)=gridloBB(dir)-1
     incr=gridloBB(dir)-FSI_lo(dir)
     xcontrol=xdata3D(i+ii*incr,j+jj*incr,k+kk*incr,dir)
     iter=iter+1
     if (iter.gt.FSI_growhi(dir)-FSI_growlo(dir)+1) then
      print *,"failure to find xcontrol"
      stop
     endif
    enddo
  
    gridhiBB(dir)=NINT( (maxnode(dir)-xlo(dir))/dxBB(dir)-half+FSI_lo(dir) )
    if (gridhiBB(dir).gt.growhi3D(dir)) then
     gridhiBB(dir)=growhi3D(dir)
    endif
    if (gridhiBB(dir).lt.growlo3D(dir)) then
     gridhiBB(dir)=growlo3D(dir)
     null_intersection=1
    endif

    if (null_intersection.eq.0) then

     incr=gridhiBB(dir)-FSI_lo(dir)
     xcontrol=xdata3D(i+ii*incr,j+jj*incr,k+kk*incr,dir)
     xcost=maxnode(dir)+local_iband*dxBB(dir)
     iter=0
     do while ((gridhiBB(dir).lt.growhi3D(dir)).and. &
               (xcontrol.lt.xcost))
      gridhiBB(dir)=gridhiBB(dir)+1
      incr=gridhiBB(dir)-FSI_lo(dir)
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
 ngrow_make_distance_in, &
 probe_size, &
 xnot, &
 FSI_lo,FSI_hi, &
 FSI_growlo,FSI_growhi,  &
 xdata3D, &
 gridloBB,gridhiBB,dxBB)

use global_utility_module
IMPLICIT NONE

 INTEGER_T, INTENT(in) :: ngrow_make_distance_in
 REAL_T, INTENT(in) :: probe_size
 REAL_T, INTENT(in) :: xnot(3)
 INTEGER_T, INTENT(in) :: FSI_lo(3),FSI_hi(3)
 INTEGER_T, INTENT(in) :: FSI_growlo(3),FSI_growhi(3)
 REAL_T, INTENT(in), pointer :: xdata3D(:,:,:,:)
 INTEGER_T, INTENT(out) :: gridloBB(3),gridhiBB(3)
 REAL_T, INTENT(out) :: dxBB(3)
 INTEGER_T dir,dirloc
 INTEGER_T idx(3),idxL(3),idxR(3)
 INTEGER_T iter,change
 REAL_T dist,distL,distR
 REAL_T xlo(3),xhi(3)
 INTEGER_T ngrowtest
 INTEGER_T interp_support
 INTEGER_T i_probe_size

 if (ngrow_make_distance_in.ne.3) then
  print *,"ngrow_make_distance_in invalid"
  stop
 endif

 i_probe_size=NINT(probe_size)
 if ((i_probe_size.eq.0).or. &
     (i_probe_size.eq.1)) then
  ! do nothing
 else
  print *,"i_probe_size invalid"
  stop
 endif
 if ((probe_size.eq.zero).or. &
     (probe_size.eq.one)) then
  ! do nothing
 else
  print *,"probe_size invalid"
  stop
 endif

 call checkbound3D_array(FSI_lo,FSI_hi, &
   xdata3D, &
   ngrow_make_distance_in,-1)

 interp_support=BoundingBoxRadNode

 do dir=1,3
  ngrowtest=FSI_lo(dir)-FSI_growlo(dir)
  if (ngrowtest.ne.ngrow_make_distance_in) then
   print *,"ngrowtest invalid1 ",ngrowtest
   stop
  endif
  ngrowtest=FSI_growhi(dir)-FSI_hi(dir)
  if (ngrowtest.ne.ngrow_make_distance_in) then
   print *,"ngrowtest invalid2 ",ngrowtest
   stop
  endif
 enddo ! dir=1..3

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

  if (xnot(dir).ge.xlo(dir)-(VOFTOL+probe_size)*dxBB(dir)) then
   ! do nothing
  else
   print *,"node should be within grid interior"
   stop
  endif
  if (xnot(dir).le.xhi(dir)+(VOFTOL+probe_size)*dxBB(dir)) then
   ! do nothing
  else
   print *,"node should be within grid interior"
   stop
  endif
  if (xnot(dir).le.xlo(dir)+(VOFTOL-probe_size)*dxBB(dir)) then
   gridloBB(dir)=FSI_lo(dir)-i_probe_size
  else if (xnot(dir).ge.xhi(dir)-(VOFTOL-probe_size)*dxBB(dir)) then
   gridloBB(dir)=FSI_hi(dir)+i_probe_size
  else if ((xnot(dir).ge.xlo(dir)-probe_size*dxBB(dir)).and. &
           (xnot(dir).le.xhi(dir)+probe_size*dxBB(dir))) then
      ! for evenly spaced points:
      ! x=xlo+(i-ilo+1/2)*dx
      ! (x-xlo)/dx -1/2 = i-ilo
      ! i=ilo+(x-xlo)/dx -1/2
   gridloBB(dir)=NINT( (xnot(dir)-xlo(dir))/dxBB(dir)-half+FSI_lo(dir) )
   if ((gridloBB(dir).lt.FSI_lo(dir)-i_probe_size).or. &
       (gridloBB(dir).gt.FSI_hi(dir)+i_probe_size)) then
    print *,"node should be within (probe biased) grid interior"
    stop
   endif
   do dirloc=1,3
    idx(dirloc)=FSI_lo(dirloc)
   enddo
   idx(dir)=gridloBB(dir)
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
    if ((idxL(dir).ge.FSI_growlo(dir)).and. &
        (idxL(dir).le.FSI_growhi(dir)).and. &
        (idxR(dir).ge.FSI_growlo(dir)).and. &
        (idxR(dir).le.FSI_growhi(dir))) then
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
    else
     print *,"idxL or idxR out of range"
     stop
    endif
    iter=iter+1
    if (iter.gt.FSI_hi(dir)-FSI_lo(dir)+1+2*i_probe_size) then
     print *,"iter.gt.FSI_hi(dir)-FSI_lo(dir)+1+2*i_probe_size"
     stop
    endif
   enddo ! while (change==1)
   gridloBB(dir)=idx(dir)
  else
   print *,"xnot(dir) invalid"
   stop
  endif

  gridhiBB(dir)=gridloBB(dir)+interp_support 
  gridloBB(dir)=gridloBB(dir)-interp_support 

 enddo  ! dir=1..sdim

return
end subroutine find_grid_bounding_box_node

! ns_header_msg_level is called from NavierStokes::nonlinear_advection
! fort_headermsg is called from NavierStokes::ns_header_msg_level
! CLSVOF_ReadNodes is called from fort_headermsg 
! (FSI_operation.eq.OP_FSI_UPDATE_NODES, 
!  FSI_sub_operation.eq.SUB_OP_FSI_DEFAULT)
! isout==1 => verbose
subroutine CLSVOF_ReadNodes( &
  im_critical, & !0..2*num_materials-1
  max_num_nodes_list, &
  max_num_elements_list, &
  num_nodes_list, &
  num_elements_list, &
  FSI_input_max_num_nodes, &
  FSI_input_max_num_elements, &
  FSI_input_num_nodes, &
  FSI_input_num_elements, &
  FSI_input_node_list, &
  FSI_input_element_list, &
  FSI_input_displacement_list, &
  FSI_input_velocity_halftime_list, &
  FSI_input_velocity_list, &
  FSI_input_force_list, &
  FSI_input_mass_list, &
  FSI_input_temperature_list, &
  FSI_output_max_num_nodes, &
  FSI_output_max_num_elements, &
  FSI_output_num_nodes, &
  FSI_output_num_elements, &
  FSI_output_node_list, &
  FSI_output_element_list, &
  FSI_output_displacement_list, &
  FSI_output_velocity_halftime_list, &
  FSI_output_velocity_list, &
  FSI_output_force_list, &
  FSI_output_mass_list, &
  FSI_output_temperature_list, &
  FSI_refine_factor, &
  FSI_bounding_box_ngrow, &
  CLSVOF_curtime, & ! t^{n+1}
  CLSVOF_dt, &
  h_small, &
  problo,probhi, &
  current_step,plot_interval, &
  ioproc,isout)
  use global_utility_module
#ifdef MVAHABFSI
use CTML_module
#endif

IMPLICIT NONE

INTEGER_T, INTENT(in) :: im_critical
INTEGER_T, INTENT(inout) :: max_num_nodes_list(num_materials)
INTEGER_T, INTENT(inout) :: max_num_elements_list(num_materials)
INTEGER_T, INTENT(inout) :: num_nodes_list(num_materials)
INTEGER_T, INTENT(inout) :: num_elements_list(num_materials)
INTEGER_T, INTENT(in) :: FSI_input_max_num_nodes
INTEGER_T, INTENT(in) :: FSI_input_max_num_elements
INTEGER_T, INTENT(in) :: FSI_input_num_nodes
INTEGER_T, INTENT(in) :: FSI_input_num_elements
REAL_T, INTENT(inout) :: FSI_input_node_list(3*FSI_input_max_num_nodes)
INTEGER_T, INTENT(inout) :: &
        FSI_input_element_list(4*FSI_input_max_num_elements)
REAL_T, INTENT(inout) :: &
        FSI_input_displacement_list(3*FSI_input_max_num_nodes)
REAL_T, INTENT(inout) :: &
        FSI_input_velocity_halftime_list(3*FSI_input_max_num_nodes)
REAL_T, INTENT(inout) :: &
        FSI_input_velocity_list(3*FSI_input_max_num_nodes)
REAL_T, INTENT(inout) :: &
        FSI_input_force_list(NCOMP_FORCE_STRESS*FSI_input_max_num_nodes)
REAL_T, INTENT(inout) :: FSI_input_mass_list(FSI_input_max_num_nodes)
REAL_T, INTENT(inout) :: FSI_input_temperature_list(FSI_input_max_num_nodes)

INTEGER_T, INTENT(in) :: FSI_output_max_num_nodes
INTEGER_T, INTENT(in) :: FSI_output_max_num_elements
INTEGER_T, INTENT(in) :: FSI_output_num_nodes
INTEGER_T, INTENT(in) :: FSI_output_num_elements
REAL_T, INTENT(inout) :: FSI_output_node_list(3*FSI_output_max_num_nodes)
INTEGER_T, INTENT(inout) :: &
        FSI_output_element_list(4*FSI_output_max_num_elements)
REAL_T, INTENT(inout) :: &
        FSI_output_displacement_list(3*FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: &
        FSI_output_velocity_halftime_list(3*FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: &
        FSI_output_velocity_list(3*FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: &
        FSI_output_force_list(NCOMP_FORCE_STRESS*FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: FSI_output_mass_list(FSI_output_max_num_nodes)
REAL_T, INTENT(inout) :: FSI_output_temperature_list(FSI_output_max_num_nodes)


INTEGER_T, INTENT(in) :: current_step,plot_interval
INTEGER_T :: initflag
INTEGER_T, INTENT(in) :: ioproc,isout
INTEGER_T :: part_id
REAL_T, INTENT(in) :: CLSVOF_curtime,CLSVOF_dt
REAL_T, INTENT(in) :: h_small
REAL_T, INTENT(in) :: problo(3),probhi(3)
INTEGER_T node_factor 
INTEGER_T ctml_part_id 
INTEGER_T fsi_part_id 
INTEGER_T :: inode_crit,inode
INTEGER_T, INTENT(in) :: FSI_refine_factor(num_materials)
INTEGER_T, INTENT(in) :: FSI_bounding_box_ngrow(num_materials)
INTEGER_T im_sanity_check
INTEGER_T :: idir,ielem,im_part

  do im_sanity_check=1,num_materials
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

  if (h_small.gt.zero) then
   ! do nothing
  else
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

   if ((im_critical.ge.0).and. &
       (im_critical.lt.num_materials)) then

    ctml_part_id=0
    do part_id=1,TOTAL_NPARTS
     im_part=im_solid_mapF(part_id)+1
     if (CTML_FSI_mat(im_part).eq.1) then 
      ctml_part_id=ctml_part_id+1
      if (ctml_part_id_map(part_id).eq.ctml_part_id) then
       if (im_part.eq.im_critical+1) then
        do inode=1,ctml_max_n_fib_nodes
         do idir=1,AMREX_SPACEDIM
          ctml_fib_pst_prev(ctml_part_id,inode,idir)= &
            FSI_input_node_list((inode-1)*3+idir)
          ctml_fib_vel_halftime_prev(ctml_part_id,inode,idir)= &
            FSI_input_velocity_halftime_list((inode-1)*3+idir)
          ctml_fib_vel_prev(ctml_part_id,inode,idir)= &
            FSI_input_velocity_list((inode-1)*3+idir)
          ctml_fib_frc(ctml_part_id,inode,idir)= &
            FSI_input_force_list((inode-1)*NCOMP_FORCE_STRESS+idir)
         enddo !idir=1,sdim
         ctml_fib_mass_prev(ctml_part_id,inode)= &
            FSI_input_mass_list(inode)
        enddo ! inode=1,ctml_max_n_fib_nodes 
       else if ((im_part.ge.1).and.(im_part.le.num_materials)) then
        ! do nothing
       else
        print *,"im_part invalid"
        stop
       endif
      else
       print *,"ctml_part_id_map(part_id) invalid"
       stop
      endif
     else if (CTML_FSI_mat(im_part).eq.0) then
      ! do nothing
     else
      print *,"CTML_FSI_mat(im_part) invalid"
      stop
     endif
    enddo ! part_id=1,TOTAL_NPARTS

   else if ((im_critical.ge.num_materials).and. &
            (im_critical.lt.2*num_materials)) then
    ! do nothing
   else
    print *,"im_critical invalid"
    stop
   endif

   if (im_critical.eq.num_materials) then

    if (CTML_FSI_flagF().eq.1) then 
#ifdef MVAHABFSI

     do part_id=1,TOTAL_NPARTS

      ctml_part_id=ctml_part_id_map(part_id)

      if ((ctml_part_id.ge.1).and. &
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
        print *,"NumIntElems is corrupt"
        stop
       endif
       if (FSI(part_id)%IntElemDim.ne.3) then
        print *,"FSI(part_id)%IntElemDim.ne.3"
        stop
       endif

       inode_crit=0 ! node index of first inactive node.
       call CTML_PUT_FORCE( &
         ctml_fib_frc, &
         ctml_n_fib_bodies, &
         ctml_max_n_fib_nodes, &
         ctml_part_id)

       do inode=1,ctml_n_fib_nodes(ctml_part_id)
        if (ctml_fib_mass_prev(ctml_part_id,inode).gt.zero) then
         ! do nothing
        else if (ctml_fib_mass_prev(ctml_part_id,inode).eq.zero) then 
         if (inode_crit.eq.0) then
          inode_crit=inode
         endif
        else
         print *,"ctml_fib_mass_prev(ctml_part_id,inode) invalid"
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

      else if (ctml_part_id.eq.0) then
       ! do nothing
      else
       print *,"ctml_part_id invalid" 
       stop
      endif

     enddo ! part_id=1,TOTAL_NPARTS

     ! declared in: CTMLFSI.F90
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
    else if (CTML_FSI_flagF().eq.0) then
     ! do nothing
    else
     print *,"CTML_FSI_flagF() invalid"
     stop
    endif

    do part_id=1,TOTAL_NPARTS

     FSI(part_id)%part_id=part_id

     ctml_part_id=ctml_part_id_map(part_id)
     fsi_part_id=fsi_part_id_map(part_id)

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
        print *,"NumIntElems is corrupt"
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
        call CTML_GET_POS_VEL_WT( &
         ctml_fib_pst, &
         ctml_fib_vel, &
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
        !NodeNormal(dir,inode) and NodeNormalBIG initialized here.
       call post_process_nodes_elements(initflag,problo,probhi, &
        FSI(part_id),part_id,TOTAL_NPARTS, &
        ioproc,isout,h_small)
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

   endif

   if ((im_critical.ge.num_materials).and. &
       (im_critical.lt.2*num_materials)) then

    ctml_part_id=0
    do part_id=1,TOTAL_NPARTS
     im_part=im_solid_mapF(part_id)+1
     if (CTML_FSI_mat(im_part).eq.1) then 
      ctml_part_id=ctml_part_id+1
      if (ctml_part_id_map(part_id).eq.ctml_part_id) then
       if (im_part.eq.im_critical-num_materials+1) then
        do inode=1,ctml_max_n_fib_nodes

         do idir=1,NCOMP_FORCE_STRESS
          FSI_output_force_list((inode-1)*NCOMP_FORCE_STRESS+idir)=zero
         enddo

         do idir=1,AMREX_SPACEDIM
          FSI_output_node_list((inode-1)*3+idir)= &
           ctml_fib_pst(ctml_part_id,inode,idir)

          FSI_output_velocity_halftime_list((inode-1)*3+idir)= &
             (FSI_output_node_list((inode-1)*3+idir)- &
              FSI_input_node_list((inode-1)*3+idir))/CLSVOF_dt

          FSI_output_displacement_list((inode-1)*3+idir)= &
            FSI_input_displacement_list((inode-1)*3+idir)+ &
            CLSVOF_dt* &
            FSI_output_velocity_halftime_list((inode-1)*3+idir)

          FSI_output_velocity_list((inode-1)*3+idir)= &
           ctml_fib_vel(ctml_part_id,inode,idir)

          FSI_output_force_list((inode-1)*NCOMP_FORCE_STRESS+idir)= &
           ctml_fib_frc(ctml_part_id,inode,idir)

         enddo !idir=1,sdim

         FSI_output_mass_list(inode)=ctml_fib_mass(ctml_part_id,inode)
         FSI_output_temperature_list(inode)=293.0d0
        enddo ! inode=1,ctml_max_n_fib_nodes 
        do ielem=1,ctml_max_n_fib_nodes
         do inode=1,4
          FSI_output_element_list(4*(ielem-1)+inode)=ielem
         enddo
        enddo
       else if ((im_part.ge.1).and.(im_part.le.num_materials)) then
        ! do nothing
       else
        print *,"im_part invalid"
        stop
       endif
      else
       print *,"ctml_part_id_map(part_id) invalid"
       stop
      endif
     else if (CTML_FSI_mat(im_part).eq.0) then
      ! do nothing
     else
      print *,"CTML_FSI_mat(im_part) invalid"
      stop
     endif
    enddo ! part_id=1,TOTAL_NPARTS

   endif

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
      xtailpos=0.5d0

      
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
         V(i)=0.1+0.1*tanh(0.5d0*Z(i)-0.75)

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


      springcons=0.5d0

      DO i=1,Nodes
           spring(i)=0.1

          If (Z(i).gt.tailpos) Then
            spring(i)=spring(i)+springcons*(10.0-Z(i))**2 &
                     -springcons*(10.0-tailpos)**2
          End If

      END DO

      END subroutine springs

end module CLSVOFCouplerIO

