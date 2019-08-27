#undef BL_LANG_CC
#define BL_LANG_FORT
#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"

#if (BL_SPACEDIM==3)
#include "LEVEL_F.H"
#include "NAVIERSTOKES_F.H"
#endif

#define SDIM 3

! 10 seconds for tail to do a full period
#define WHALE_LENGTH 13.0
#define PERIOD_TAIL 10.0
#define DT_DUFFY 0.01
#define STEPS_DUFFY 800
#define SCIGRIDSIZE 256.0
#define DEBUG_ISOSURFACE 0
#define MARCO 0



! crossing_ls is closest distance using "height distance"
! if a Lagrangian face cuts inbetween a cell and one of its 26
! neighbors, then the sign distance of the cell is updated.

! lsface first made magnitude of height fraction.
! lsface aka dataface aka CellLSFACE
! lsface is not used

! in older versions, quadratic interpolation transferred data from
! auxiliary grid to AMR grid, now slope limited linear interpolation
! transfers data.

module CLSVOFCouplerIO

implicit none

type lag_type
 INTEGER_T :: n_nodes,n_elems
 REAL_T, pointer :: nd(:,:)
 REAL_T, pointer :: ndvel(:,:)
 REAL_T, pointer :: ndtemp(:)
 INTEGER_T, pointer :: elemdt(:,:)
 INTEGER_T, pointer :: intelemdt(:,:)
end type lag_type

type mesh_type
 INTEGER_T :: PartID
 INTEGER_T :: IntElemDim,IntElemDimPaddle,IntElemDimPool
 INTEGER_T :: NumNodes,NumNodesPaddle,NumNodesPool
 INTEGER_T :: NumIntElems,NumIntElemsPaddle,NumIntElemsPool
 INTEGER_T, pointer :: ElemData(:,:)
 INTEGER_T, pointer :: IntElem(:,:)
 INTEGER_T, pointer :: Eul2IntNode(:)
 INTEGER_T :: NumNodesBIG,NumIntElemsBIG
 INTEGER_T, pointer :: ElemNodeCount(:)
 REAL_T, pointer :: NodeBIG(:,:)
 REAL_T, pointer :: NodeVelBIG(:,:)
 REAL_T, pointer :: NodeTempBIG(:)
 REAL_T, pointer :: NodeNormalBIG(:,:)
 REAL_T, pointer :: ElemDataXnotBIG(:,:)
 INTEGER_T, pointer :: ElemDataBIG(:,:)
 INTEGER_T, pointer :: IntElemBIG(:,:)
 REAL_T, pointer :: Node(:,:)
 REAL_T, pointer :: Node_old(:,:)
 REAL_T, pointer :: Node_new(:,:)
 REAL_T, pointer :: Node_current(:,:)
 REAL_T, pointer :: NodeVel(:,:)
 REAL_T, pointer :: NodeVel_old(:,:)
 REAL_T, pointer :: NodeVel_new(:,:)
 REAL_T, pointer :: NodeTemp(:)
 REAL_T, pointer :: NodeTemp_old(:)
 REAL_T, pointer :: NodeTemp_new(:)
 logical(4), pointer :: ActiveIntElem(:)
 REAL_T soliddrop_displacement
 REAL_T soliddrop_speed
 REAL_T solid_displ(SDIM)
 REAL_T solid_speed(SDIM)
end type mesh_type

type auxgrid_type
 REAL_T, pointer :: CellMASK(:,:,:)
 REAL_T, pointer :: CellLS(:,:,:)
 REAL_T, pointer :: CellLSFACE(:,:,:,:)
 REAL_T, pointer :: CellVEL(:,:,:,:)
 REAL_T, pointer :: CellTEMP(:,:,:)
 REAL_T celldx(SDIM)
 INTEGER_T cell_lo(SDIM),cell_hi(SDIM)
 REAL_T Cellproblo(SDIM),Cellprobhi(SDIM)
end type auxgrid_type

! initialized when init_dimensions(SCI_in) is called
INTEGER_T DIMDEC3D(CellLS)

logical(4) :: DebugCouple
INTEGER_T :: FlagDim
INTEGER_T :: UnitCouple 
INTEGER_T :: UnitOut
INTEGER_T :: UnitDebug
INTEGER_T :: StdOut
INTEGER_T :: use_temp
INTEGER_T, dimension(4) :: NumFlags
INTEGER_T, dimension(:), pointer :: HdrFlag
INTEGER_T :: istepB,sci_sdim,sci_istop,sci_istep
REAL_T :: sci_curtime,sci_dt
REAL_T :: timeB,tstart,tfinish
REAL_T :: dtB
INTEGER_T, dimension(:), pointer :: NptFlag
INTEGER_T, dimension(:), pointer :: ForFlag

type(lag_type), dimension(:), allocatable :: multi_lag

type(mesh_type) :: FSI
type(mesh_type) :: FSI_needle
type(mesh_type) :: viorel_FSI

type(auxgrid_type) :: SCI
type(auxgrid_type) :: SCI_loc
type(auxgrid_type) :: SCI_needle
type(auxgrid_type) :: SCI_needle_loc

INTEGER_T :: lag_dim

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

contains



      subroutine checkbound3D(lo,hi, &
      DIMS3D(data), &
      ngrow,dir,id)
      IMPLICIT NONE

      INTEGER_T    lo(SDIM), hi(SDIM)
      INTEGER_T    DIMDEC3D(data)
      INTEGER_T    ngrow,dir,id

      INTEGER_T ii,jj,kk

      INTEGER_T    hidata(SDIM)

      hidata(1)=ARG3D_H1(data)
      if (ARG3D_L1(data).gt.hidata(1)) then
       print *,"swapped bounds in checkbound 3d id=",id
       stop
      endif
      hidata(2)=ARG3D_H2(data)
      if (ARG3D_L2(data).gt.hidata(2)) then
       print *,"swapped bounds in checkbound 3d id=",id
       stop
      endif
      hidata(3)=ARG3D_H3(data)
      if (ARG3D_L3(data).gt.hidata(3)) then
       print *,"swapped bounds in checkbound 3d id=",id
       stop
      endif

 
      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if (dir.eq.2) then
       kk=1
      else if (dir.ne.-1) then
       print *,"dir out of range in checkbound 3d id=",id
       stop
      endif

      if ((lo(1).lt.0).or.(lo(2).lt.0).or.(lo(3).lt.0)) then
       print *,"lo invalid in checkbound 3d id=",id
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

      if (ARG3D_L1(data).gt.lo(1)-ngrow) then
       print *,"lo mismatch id=",id
       stop
      endif
      if (ARG3D_L2(data).gt.lo(2)-ngrow) then
       print *,"lo mismatch id=",id
       stop
      endif
      if (ARG3D_L3(data).gt.lo(3)-ngrow) then
       print *,"lo mismatch id=",id
       stop
      endif

      if (ARG3D_H1(data).lt.hi(1)+ngrow+ii) then
       print *,"hi mismatch id=",id
       stop
      endif
      if (ARG3D_H2(data).lt.hi(2)+ngrow+jj) then
       print *,"hi mismatch id=",id
       stop
      endif
      if (ARG3D_H3(data).lt.hi(3)+ngrow+kk) then
       print *,"hi mismatch id=",id
       stop
      endif

      return
      end subroutine checkbound3D


subroutine init2_FSI(FSI_in)
type(mesh_type) :: FSI_in
INTEGER_T inode,dir,iface

 do inode=1,FSI_in%NumNodes
  do dir=1,3
   FSI_in%NodeVel_old(dir,inode)=0.0
  enddo
  do dir=1,3
   FSI_in%Node_current(dir,inode)=FSI_in%Node_new(dir,inode)
   FSI_in%NodeVel_new(dir,inode)=FSI_in%NodeVel_old(dir,inode)
  enddo
 enddo  ! inode=1,NumNodes
return
end subroutine init2_FSI


subroutine init3_FSI(FSI_in,ifirst,do_2nd_part,ioproc,isout)
type(mesh_type) :: FSI_in
INTEGER_T inode,dir,iface,ifirst,do_2nd_part,it,ioproc,isout
REAL_T x,y,z,z0,z90,t,dt,t1,t2,inflowvel
REAL_T YK,ZK,lift0,lift90
REAL_T, dimension(3) :: displ1,displ2
REAL_T dilated_time

#include "probdataf95.H"

  if (ifirst.eq.1) then
   allocate(FSI_in%Node(3,FSI_in%NumNodes))
   allocate(FSI_in%NodeVel(3,FSI_in%NumNodes))
   allocate(FSI_in%NodeTemp(FSI_in%NumNodes))
  else if (ifirst.eq.0) then
   ! do nothing
  else
   print *,"ifirst invalid"
   stop
  endif

  do inode=1,FSI_in%NumNodes
   do dir=1,3
    FSI_in%Node(dir,inode)=FSI_in%Node_current(dir,inode)
    FSI_in%NodeVel(dir,inode)=FSI_in%NodeVel_new(dir,inode)
   enddo
   FSI_in%NodeTemp(inode)=FSI_in%NodeTemp_new(inode)
  enddo

  do dir=1,3
   FSI_in%solid_displ(dir)=0.
   FSI_in%solid_speed(dir)=0.
  enddo

  FSI_in%soliddrop_displacement=0.0
  FSI_in%soliddrop_speed=0.0

  if (probtype.eq.531) then

   if (timeB.le.0.03) then  ! actual time x 10^3 ?
    FSI_in%soliddrop_displacement=-0.115*timeB
    FSI_in%soliddrop_speed=-0.115
   else
    FSI_in%soliddrop_displacement=-0.115*0.03 ! actual time x 10^3 ?
    FSI_in%soliddrop_speed=0.0
   endif

! MARK:
! if probtype=538 and partid=1 and RZ, then 
! dist(x,y,z)=dist(x+0.0015,y-0.0045,z)
! if probtype=538 and partid=2, then
! dist(x,y,z)= dist(x,y,z+0.01)
! IN FUTURE, DO NOT CALL "UNITE", INSTEAD
! USE BOTH FSI and FSI_NEEDLE when returning LS or VEL.
! ALSO IN FUTURE, DO NOT REPEATEDLY CALL GENERATE_NEW_TRIANGLES
! IN CLSVOF_ReadNodes.
!  
  else if (probtype.eq.538) then

   if (FSI_in%PartID.eq.2) then

    if (MARCO.eq.0) then
     dilated_time=timeB
     dilated_time = dilated_time*51.
    else
     dilated_time=timeB-1.460e-3 ! skip first 200 mus
    endif

    ! load trajectory of needle tip

    OPEN(unit=15,file="needle_motion.dat",access='sequential', &
     form="formatted",status='old')
    dt = 0.
    do it=1,251

     if (MARCO.eq.0) then
      READ(15,*) t,z0,z90,x,y
     else
      READ(15,*) t,lift0,lift90,YK,ZK
     endif

     if (MARCO.eq.0) then
      t = t*1e-6 ! mus to s
      x = (x+24.32)*1e-4 ! mum to cm (IJFM 2013)
      y = -(y+33.56)*1e-4
      z0 = (z0-468.0)*1e-4 ! at t = 910 mus
      z90 = (z90-470.75)*1e-4 
     else
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
     stop
    endif
    do dir=1,3
     FSI_in%solid_displ(dir)= &
       displ1(dir)+(dilated_time-t1)/dt*(displ2(dir)-displ1(dir))
     FSI_in%solid_speed(dir)= &
       (displ2(dir)-displ1(dir))/dt
    enddo
    if (levelrz.eq.0) then
     ! do nothing
    else if (levelrz.eq.1) then
     FSI_in%solid_speed(1)=zero
     FSI_in%solid_speed(2)=zero
    else
     print *,"levelrz invalid"
     stop
    endif 

    if ((ioproc.eq.1).and.(isout.eq.1)) then
     print*,"sci_clsvof: interpolation at main time ",timeB
     print*,"sci_clsvof: interpolation at dilated time ",dilated_time, &
       " between ",t1," and ",t2
     print*,"solid displacement: ",FSI_in%solid_displ
     print*,"solid velocity: ",FSI_in%solid_speed
    endif
   else if (FSI_in%PartID.eq.1) then

    do dir=1,3
     FSI_in%solid_displ(dir)=0.
     FSI_in%solid_speed(dir)=0.
    enddo

   else
    print *,"part id invalid"
    stop
   endif
    
  endif  ! probtype.eq.538

  if (do_2nd_part.eq.1) then

    ! solid body translation is default

   do inode=1,FSI_in%NumNodes
    x=FSI_in%Node(1,inode)
    y=FSI_in%Node(2,inode)
    z=FSI_in%Node(3,inode)

     ! viorel's sphere problem, cannot specify velocity at outer shell
     ! nodes because they are too far apart.
    if (probtype.eq.5601) then  ! viorel's sphere problem
     inflowvel=0.0
     FSI_in%NodeVel_old(1,inode)=inflowvel
     FSI_in%NodeVel_new(1,inode)=inflowvel
     FSI_in%NodeVel(1,inode)=inflowvel
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
     FSI_in%NodeVel_old(3,inode)=inflowvel
     FSI_in%NodeVel_new(3,inode)=inflowvel
     FSI_in%NodeVel(3,inode)=inflowvel
    else if (1.eq.0) then
     FSI_in%NodeVel_old(3,inode)=FSI_in%soliddrop_speed
     FSI_in%NodeVel_new(3,inode)=FSI_in%soliddrop_speed
     FSI_in%NodeVel(3,inode)=FSI_in%soliddrop_speed
     FSI_in%Node(3,inode)=FSI_in%Node_old(3,inode)+ &
       FSI_in%soliddrop_displacement
     FSI_in%Node_new(3,inode)=FSI_in%Node(3,inode)
    else
     do dir=1,3
      FSI_in%NodeVel_old(dir,inode)=FSI_in%solid_speed(dir)
      FSI_in%NodeVel_new(dir,inode)=FSI_in%solid_speed(dir)
      FSI_in%NodeVel(dir,inode)=FSI_in%solid_speed(dir)

! MARK: displacement taken into account elsewhere.
! do not want to change the node positions since it will
! effect what happens in generate_new_triangles.

      if (probtype.eq.538) then
       FSI_in%Node(dir,inode)=FSI_in%Node_old(dir,inode)
      else
       FSI_in%Node(dir,inode)=FSI_in%Node_old(dir,inode)+ &
         FSI_in%solid_displ(dir)
      endif
      FSI_in%Node_new(dir,inode)=FSI_in%Node(dir,inode)
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



subroutine init_FSI(FSI_in,allocate_intelem)
IMPLICIT NONE

type(mesh_type) :: FSI_in
INTEGER_T inode,dir,iface,allocate_intelem

 allocate(FSI_in%Eul2IntNode(FSI_in%NumNodes))
 allocate(FSI_in%ElemData(NumFlags(2),FSI_in%NumIntElems))
 if (allocate_intelem.eq.1) then
  allocate(FSI_in%IntElem(FSI_in%IntElemDim,FSI_in%NumIntElems))
 endif
 allocate(FSI_in%ActiveIntElem(FSI_in%NumIntElems))
 allocate(FSI_in%Node_old(3,FSI_in%NumNodes))
 allocate(FSI_in%Node_new(3,FSI_in%NumNodes))
 allocate(FSI_in%Node_current(3,FSI_in%NumNodes))
 allocate(FSI_in%NodeVel_old(3,FSI_in%NumNodes))
 allocate(FSI_in%NodeVel_new(3,FSI_in%NumNodes))
 allocate(FSI_in%NodeTemp_old(FSI_in%NumNodes))
 allocate(FSI_in%NodeTemp_new(FSI_in%NumNodes))

 do inode=1,FSI_in%NumNodes
  FSI_in%Eul2IntNode(inode)=inode
  FSI_in%NodeTemp_old(inode)=0.0
  FSI_in%NodeTemp_new(inode)=0.0
  do dir=1,3
   FSI_in%Node_old(dir,inode)=0.0
   FSI_in%Node_new(dir,inode)=0.0
   FSI_in%Node_current(dir,inode)=0.0
   FSI_in%NodeVel_old(dir,inode)=0.0
   FSI_in%NodeVel_new(dir,inode)=0.0
  enddo
 enddo
 do iface=1,FSI_in%NumIntElems
  FSI_in%ActiveIntElem(iface)=.true.
 enddo

return
end subroutine init_FSI

subroutine init_dimensions(SCI_in)
IMPLICIT NONE

type(auxgrid_type) :: SCI_in

 ARG3D_L1(CellLS)=SCI_in%cell_lo(1)-1
 ARG3D_L2(CellLS)=SCI_in%cell_lo(2)-1
 ARG3D_L3(CellLS)=SCI_in%cell_lo(3)-1
 ARG3D_H1(CellLS)=SCI_in%cell_hi(1)+1
 ARG3D_H2(CellLS)=SCI_in%cell_hi(2)+1
 ARG3D_H3(CellLS)=SCI_in%cell_hi(3)+1

return
end subroutine init_dimensions

subroutine xdist(x1,x2,dist)
IMPLICIT NONE
 
REAL_T, dimension(3),intent(in) :: x1,x2
REAL_T, intent(out) :: dist

 dist=sqrt( (x1(1)-x2(1))**2+ &
            (x1(2)-x2(2))**2+ &
            (x1(3)-x2(3))**2 )

return
end subroutine xdist

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

subroutine isect(boxlo,boxhi,ielem,iflag,FSI_in)
IMPLICIT NONE

type(mesh_type) :: FSI_in
INTEGER_T :: ielem
REAL_T :: dist,eps
REAL_T, dimension(3) :: boxlo,boxhi
INTEGER_T, intent(out) :: iflag
INTEGER_T :: nodes_per_elem,inode,dir,istat
REAL_T, dimension(3) :: xnode

 nodes_per_elem=FSI_in%ElemDataBIG(1,ielem)
 if (nodes_per_elem.ne.3) then
  print *,"nodes_per_elem invalid"
  stop
 endif
 call xdist(boxlo,boxhi,dist)
 eps=0.0
 iflag=0
 do inode=1,nodes_per_elem
  do dir=1,3
   xnode(dir)=FSI_in%NodeBIG(dir,FSI_in%IntElemBIG(inode,ielem))
  enddo
  call checkinbox(xnode,boxlo,boxhi,istat,eps)
  if (istat.eq.1) then
   iflag=1
  endif
 enddo 
 
return
end subroutine isect

subroutine simpleBC(data,DIMS3D(data),scomp,clo,chi)
IMPLICIT NONE

INTEGER_T clo(SDIM)
INTEGER_T chi(SDIM)
INTEGER_T DIMDEC3D(data)
INTEGER_T scomp
REAL_T data(DIMV3D(data),scomp)
INTEGER_T iprime,jprime,kprime,i,j,k

 if ((scomp.lt.1).or.(scomp.gt.SDIM)) then
  print *,"scomp invalid"
  stop
 endif

 call checkbound3D(clo,chi,DIMS3D(data),1,-1,333)
 do i=clo(1)-1,chi(1)+1
 do j=clo(2)-1,chi(2)+1
 do k=clo(3)-1,chi(3)+1
  iprime=i
  jprime=j
  kprime=k
  if (i.lt.clo(1)) then
   iprime=clo(1)
  endif
  if (i.gt.chi(1)) then
   iprime=chi(1)
  endif
  if (j.lt.clo(2)) then
   jprime=clo(2)
  endif
  if (j.gt.chi(2)) then
   jprime=chi(2)
  endif
  if (k.lt.clo(3)) then
   kprime=clo(3)
  endif
  if (k.gt.chi(3)) then
   kprime=chi(3)
  endif
  data(i,j,k,scomp)=data(iprime,jprime,kprime,scomp)
 enddo
 enddo
 enddo

return
end subroutine simpleBC

! mask=0 velocity not init and sign not init.
! mask=1 velocity init sign not init
! mask=2 velocity and sign init
subroutine extrap_iter(data,DIMS3D(data), &
  scomp,clo,chi,SCI_in)
IMPLICIT NONE

type(auxgrid_type) :: SCI_in
INTEGER_T clo(SDIM)
INTEGER_T chi(SDIM)
INTEGER_T DIMDEC3D(data)
INTEGER_T scomp
REAL_T data(DIMV3D(data),scomp)

REAL_T, dimension(:,:,:), allocatable :: mask_local

INTEGER_T i,j,k
INTEGER_T i1,j1,k1,dir
INTEGER_T ngrow
REAL_T weight_top,weight_bot,weight
REAL_T, dimension(3) :: maxnodedomain,minnodedomain

#include "probdataf95.H"

 if ((scomp.lt.1).or.(scomp.gt.SDIM)) then
  print *,"scomp invalid"
  stop
 endif

 call init_dimensions(SCI_in)
 call checkbound3D(clo,chi,DIMS3D(data),1,-1,333)

 allocate(mask_local(DIMV3D(CellLS)))
 do i=clo(1)-1,chi(1)+1
 do j=clo(2)-1,chi(2)+1
 do k=clo(3)-1,chi(3)+1
  mask_local(i,j,k)=SCI_in%CellMASK(i,j,k)
 enddo
 enddo
 enddo

 do i=clo(1),chi(1)
 do j=clo(2),chi(2)
 do k=clo(3),chi(3)

   if (mask_local(i,j,k).eq.zero) then
    weight_top=zero
    weight_bot=zero
    do i1=-3,3
    do j1=-3,3
    do k1=-3,3
     if ((i+i1.ge.clo(1)).and.(i+i1.le.chi(1)).and. &
         (j+j1.ge.clo(2)).and.(j+j1.le.chi(2)).and. &
         (k+k1.ge.clo(3)).and.(k+k1.le.chi(3))) then
      if (mask_local(i+i1,j+j1,k+k1).gt.zero) then
       weight=i1**2+j1**2+k1**2
       if (weight.le.zero) then
        print *,"weight invalid"
        stop
       endif
       weight=one/weight
       weight_top=weight_top+data(i+i1,j+j1,k+k1,scomp)*weight
       weight_bot=weight_bot+weight
      endif ! mask>0
     endif ! i1,j1,k1 in grid
    enddo
    enddo
    enddo
    if (weight_bot.gt.zero) then
     data(i,j,k,scomp)=weight_top/weight_bot
    endif 
   endif ! mask=0
 enddo
 enddo
 enddo ! i,j,k

    ! neumann bc
 call simpleBC(data,DIMS3D(data),scomp,clo,chi)

 deallocate(mask_local)

return
end subroutine extrap_iter


! mask=0 velocity not init and sign not init.
! mask=1 velocity init sign not init
! mask=2 velocity and sign init
! mask=3 velocity and sign init doubly wetted
subroutine sign_extrap_iter(data,dataface, &
  DIMS3D(data), &
  clo,chi,FSI_in,SCI_in,problo,dx,ioproc)
IMPLICIT NONE

INTEGER_T ioproc
type(mesh_type) :: FSI_in
type(auxgrid_type) :: SCI_in
REAL_T problo(SDIM)
REAL_T dx(SDIM)
INTEGER_T clo(SDIM)
INTEGER_T chi(SDIM)
INTEGER_T DIMDEC3D(data)
REAL_T data(DIMV3D(data))
REAL_T dataface(DIMV3D(data),SDIM)

REAL_T, dimension(:,:,:), allocatable :: mask_local
REAL_T, dimension(:,:,:), allocatable :: old_state
REAL_T, dimension(:,:,:), allocatable :: old_mask

INTEGER_T i,j,k,dir,side,max_iter,inode
INTEGER_T touch_flag
INTEGER_T presbc(SDIM,2)
REAL_T, dimension(3) :: maxnodedomain,minnodedomain,xx
REAL_T, dimension(3) :: COM
REAL_T :: maxnoderadius
REAL_T maxside
REAL_T mcen,mplus,mminus,dcen,dplus,dminus
INTEGER_T ii,jj,kk
INTEGER_T niter,numtouch
REAL_T dist

#include "probdataf95.H"

 if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
 endif
 call GetDomainScale(maxnodedomain,minnodedomain,maxside, &
   maxnoderadius,COM,FSI_in)

 if (invert_solid_levelset.eq.1) then
  do i=clo(1)-1,chi(1)+1
  do j=clo(2)-1,chi(2)+1
  do k=clo(3)-1,chi(3)+1
   data(i,j,k)=-data(i,j,k)
  enddo
  enddo
  enddo
 else if (invert_solid_levelset.ne.0) then
  print *,"invert_solid_levelset invalid"
  stop
 endif

 call init_dimensions(SCI_in)
 call checkbound3D(clo,chi,DIMS3D(data),1,-1,333)

 allocate(mask_local(DIMV3D(CellLS)))
 do i=clo(1)-1,chi(1)+1
 do j=clo(2)-1,chi(2)+1
 do k=clo(3)-1,chi(3)+1
  mask_local(i,j,k)=SCI_in%CellMASK(i,j,k)
 enddo
 enddo
 enddo

 touch_flag=1
 niter=0
 do while (touch_flag.eq.1)
   touch_flag=0
   allocate(old_state(DIMV3D(CellLS)))
   allocate(old_mask(DIMV3D(CellLS)))
   do i=clo(1)-1,chi(1)+1
   do j=clo(2)-1,chi(2)+1
   do k=clo(3)-1,chi(3)+1
    old_state(i,j,k)=data(i,j,k)
    old_mask(i,j,k)=mask_local(i,j,k) 
   enddo
   enddo
   enddo

   if (niter.eq.0) then

     ! speed up whale or ship conversion
    if (MARCO.eq.0) then
     if ((probtype.eq.562).or.(probtype.eq.9)) then

      dir=2
      ii=0
      jj=1
      kk=0
      do i=clo(1),chi(1)
      do k=clo(3),chi(3)
       j=chi(dir)
       do while ((j.gt.clo(dir)).and.(abs(dataface(i,j,k,dir)).gt.one))
        old_mask(i,j,k)=two
        mask_local(i,j,k)=two
        old_state(i,j,k)=abs(old_state(i,j,k))
        data(i,j,k)=abs(data(i,j,k))
        j=j-1
       enddo  ! j
      enddo
      enddo ! i,k

     endif  ! probtype=562 or 9
    endif ! speed up the sign for the whale or ship problem 

     ! flapping wing problem (since sides are open)
     ! also accelerate finding sign
     ! auxiliary distance is negative in the solid
    if (probtype.eq.701) then  
      do i=clo(1),chi(1)
      do j=clo(2),chi(2)
      do k=clo(3),chi(3)
       xx(1)=SCI_in%Cellproblo(1)+(i-clo(1)+half)*SCI_in%celldx(1)
       xx(2)=SCI_in%Cellproblo(2)+(j-clo(2)+half)*SCI_in%celldx(2)
       xx(3)=SCI_in%Cellproblo(3)+(k-clo(3)+half)*SCI_in%celldx(3)

       do dir=1,SDIM
        if (dir.ne.2) then
         if ((xx(dir).lt.minnodedomain(dir)).or. &
             (xx(dir).gt.maxnodedomain(dir))) then
          old_mask(i,j,k)=two
          mask_local(i,j,k)=two
          old_state(i,j,k)=abs(old_state(i,j,k))
          data(i,j,k)=abs(data(i,j,k))
         endif
        endif
       enddo
      enddo
      enddo
      enddo
    endif  ! probtype=701

      ! speed up injector conversion
    if (MARCO.eq.1) then
     if (probtype.eq.537) then  
      do i=clo(1),chi(1)
      do j=clo(2),chi(2)
      do k=clo(3),chi(3)
       xx(1)=SCI_in%Cellproblo(1)+(i-clo(1)+half)*SCI_in%celldx(1)
       xx(2)=SCI_in%Cellproblo(2)+(j-clo(2)+half)*SCI_in%celldx(2)
       xx(3)=SCI_in%Cellproblo(3)+(k-clo(3)+half)*SCI_in%celldx(3)
       do dir=1,SDIM
        if ((xx(dir).lt.minnodedomain(dir)).or. &
            (xx(dir).gt.maxnodedomain(dir))) then
         old_mask(i,j,k)=two
         mask_local(i,j,k)=two
         old_state(i,j,k)=abs(old_state(i,j,k))
         data(i,j,k)=abs(data(i,j,k))
        endif
       enddo
      enddo
      enddo
      enddo
     endif  ! probtype=537
    endif

     ! speed up gear conversion (will not speed it up much, since most
     ! of time is spent initializing the sign within the gear. (-) )
    if (MARCO.eq.0) then
     if ((probtype.eq.563).and.(axis_dir.eq.2)) then

      do i=clo(1),chi(1)
      do j=clo(2),chi(2)
      do k=clo(3),chi(3)
       xx(1)=SCI_in%Cellproblo(1)+(i-clo(1)+half)*SCI_in%celldx(1)
       xx(2)=SCI_in%Cellproblo(2)+(j-clo(2)+half)*SCI_in%celldx(2)
       xx(3)=SCI_in%Cellproblo(3)+(k-clo(3)+half)*SCI_in%celldx(3)
       dist=zero
       do dir=1,SDIM
        dist=dist+(xx(dir)-COM(dir))**2
       enddo
       dist=sqrt(dist)
       if (dist.gt.maxnoderadius) then
        old_mask(i,j,k)=two
        mask_local(i,j,k)=two
        old_state(i,j,k)=abs(old_state(i,j,k))
        data(i,j,k)=abs(data(i,j,k))
       endif
      enddo
      enddo
      enddo

      dir=1
      ii=1
      jj=0
      kk=0
      do j=clo(2),chi(2)
      do k=clo(3),chi(3)
       i=chi(dir)
       do while ((i.gt.clo(dir)).and.(abs(dataface(i,j,k,dir)).gt.one))
        old_mask(i,j,k)=two
        mask_local(i,j,k)=two
        old_state(i,j,k)=abs(old_state(i,j,k))
        data(i,j,k)=abs(data(i,j,k))
        i=i-1
       enddo  ! i
       i=clo(dir)
       do while ((i.lt.chi(dir)).and. &
        (abs(dataface(i+ii,j+jj,k+kk,dir)).gt.one))
        old_mask(i,j,k)=two
        mask_local(i,j,k)=two
        old_state(i,j,k)=abs(old_state(i,j,k))
        data(i,j,k)=abs(data(i,j,k))
        i=i+1
       enddo  ! i
      enddo
      enddo ! j,k

      dir=3
      ii=0
      jj=0
      kk=1
      do i=clo(1),chi(1)
      do j=clo(2),chi(2)
       k=chi(dir)
       do while ((k.gt.clo(dir)).and.(abs(dataface(i,j,k,dir)).gt.one))
        old_mask(i,j,k)=two
        mask_local(i,j,k)=two
        old_state(i,j,k)=abs(old_state(i,j,k))
        data(i,j,k)=abs(data(i,j,k))
        k=k-1
       enddo  ! k
       k=clo(dir)
       do while ((k.lt.chi(dir)).and. &
          (abs(dataface(i+ii,j+jj,k+kk,dir)).gt.one))
        old_mask(i,j,k)=two
        mask_local(i,j,k)=two
        old_state(i,j,k)=abs(old_state(i,j,k))
        data(i,j,k)=abs(data(i,j,k))
        k=k+1
       enddo  ! k
      enddo
      enddo ! i,j

     endif  ! probtype=563 and axis_dir=2
    endif


    if (1.eq.0) then  

     dir=2
     ii=0
     jj=1
     kk=0
     do i=clo(1),chi(1)
     do k=clo(3),chi(3)
      j=chi(dir)
      do while ((j.gt.clo(dir)).and.(abs(dataface(i,j,k,dir)).gt.one))
       old_mask(i,j,k)=two
       mask_local(i,j,k)=two
       old_state(i,j,k)=abs(old_state(i,j,k))
       data(i,j,k)=abs(data(i,j,k))
       j=j-1
      enddo  ! j
      j=clo(dir)
      do while ((j.lt.chi(dir)).and. &
        (abs(dataface(i+ii,j+jj,k+kk,dir)).gt.one))
       old_mask(i,j,k)=two
       mask_local(i,j,k)=two
       old_state(i,j,k)=abs(old_state(i,j,k))
       data(i,j,k)=abs(data(i,j,k))
       j=j+1
      enddo  ! j
     enddo
     enddo ! i,k

     dir=1
     ii=1
     jj=0
     kk=0
     do j=clo(2),chi(2)
     do k=clo(3),chi(3)
      i=chi(dir)
      do while ((i.gt.clo(dir)).and.(abs(dataface(i,j,k,dir)).gt.one))
       old_mask(i,j,k)=two
       mask_local(i,j,k)=two
       old_state(i,j,k)=abs(old_state(i,j,k))
       data(i,j,k)=abs(data(i,j,k))
       i=i-1
      enddo  ! i
      i=clo(dir)
      do while ((i.lt.chi(dir)).and. &
       (abs(dataface(i+ii,j+jj,k+kk,dir)).gt.one))
       old_mask(i,j,k)=two
       mask_local(i,j,k)=two
       old_state(i,j,k)=abs(old_state(i,j,k))
       data(i,j,k)=abs(data(i,j,k))
       i=i+1
      enddo  ! i
     enddo
     enddo ! j,k

     dir=3
     ii=0
     jj=0
     kk=1
     do i=clo(1),chi(1)
     do j=clo(2),chi(2)
      k=chi(dir)
      do while ((k.gt.clo(dir)).and.(abs(dataface(i,j,k,dir)).gt.one))
       old_mask(i,j,k)=two
       mask_local(i,j,k)=two
       old_state(i,j,k)=abs(old_state(i,j,k))
       data(i,j,k)=abs(data(i,j,k))
       k=k-1
      enddo  ! k
      k=clo(dir)
      do while ((k.lt.chi(dir)).and. &
         (abs(dataface(i+ii,j+jj,k+kk,dir)).gt.one))
       old_mask(i,j,k)=two
       mask_local(i,j,k)=two
       old_state(i,j,k)=abs(old_state(i,j,k))
       data(i,j,k)=abs(data(i,j,k))
       k=k+1
      enddo  ! k
     enddo
     enddo ! i,j

    endif  ! initializing sign to be positive in exterior regions

   endif ! niter=0

   numtouch=0

   do i=clo(1),chi(1)
   do j=clo(2),chi(2)
   do k=clo(3),chi(3)
    mcen=old_mask(i,j,k)
    if (mcen.lt.two) then
     do dir=1,SDIM
      ii=0
      jj=0
      kk=0
      if (dir.eq.1) then
       ii=1
      else if (dir.eq.2) then
       jj=1
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir invalid"
       stop
      endif
      dcen=old_state(i,j,k)
      dplus=old_state(i+ii,j+jj,k+kk)
      dminus=old_state(i-ii,j-jj,k-kk)
      mplus=old_mask(i+ii,j+jj,k+kk)
      mminus=old_mask(i-ii,j-jj,k-kk)
      if (mplus.ge.two) then
       numtouch=numtouch+1
       if (dplus.gt.zero) then
        data(i,j,k)=abs(data(i,j,k))
        mask_local(i,j,k)=two
        touch_flag=1
       else if (dplus.lt.zero) then
        data(i,j,k)=-abs(data(i,j,k))
        mask_local(i,j,k)=two
        touch_flag=1
       endif
      endif
      if (mminus.ge.two) then
       numtouch=numtouch+1
       if (dminus.gt.zero) then
        data(i,j,k)=abs(data(i,j,k))
        mask_local(i,j,k)=two
        touch_flag=1
       else if (dminus.lt.zero) then
        data(i,j,k)=-abs(data(i,j,k))
        mask_local(i,j,k)=two
        touch_flag=1
       endif
      endif
     enddo !dir
    endif !mcen<2
   enddo 
   enddo 
   enddo  ! i,j,k

   call simpleBC(mask_local,DIMS3D(CellLS),1,clo,chi)
   call simpleBC(data,DIMS3D(data),1,clo,chi)

   deallocate(old_mask)
   deallocate(old_state)
  
   niter=niter+1 
   if (ioproc.eq.1) then
    print *,"niter=",niter
    print *,"numtouch ",numtouch
   endif
 enddo  ! touch_flag.eq.1

 deallocate(mask_local)

return
end subroutine sign_extrap_iter



! split triangles so that size is no bigger than "h_small"
subroutine generate_new_triangles(initflag,problo,probhi, &
  FSI_in,SCI_in,islocal,ioproc,isout)
IMPLICIT NONE

type(auxgrid_type) :: SCI_in
type(mesh_type) :: FSI_in
REAL_T problo(SDIM),probhi(SDIM)
INTEGER_T :: initflag,ioproc,isout,islocal
INTEGER_T :: ielem,nodes_per_elem,inode,jnode,i,j,k,dir
INTEGER_T :: scount,ncount,total_new_count,nodes_new_count,iupdate
INTEGER_T :: new_NumIntElems
REAL_T :: mindist,tempdist,h_small,maxdist
REAL_T, dimension(3) :: x1,x2,x3,maxnodedomain,minnodedomain
REAL_T, dimension(3) :: COM
REAL_T :: maxnoderadius
REAL_T, dimension(3) :: vel1,vel2,vel3,xsplit,velsplit
REAL_T :: diag,biggest_h,temp_h,d1,d2,d3,temp1,temp2,temp3,tempsplit
INTEGER_T :: ibox,ii,jj,kk,ihit,iflag,iter
INTEGER_T :: temp_leaves,temp_elements,temp_ptr
INTEGER_T :: current_number_leaves,current_number_elements
REAL_T, dimension(3) :: cboxlo,cboxhi,boxlo,boxhi,normal
INTEGER_T :: ilevel,nsplit,esplit,isub,normal_cnt,base_ielem
INTEGER_T :: node1,node2,node3,lag_dim,save_n_elems,save_n_nodes
REAL_T :: h_buffer
INTEGER_T n_cell(SDIM)
INTEGER_T scomp,redistflag
REAL_T local_sci_grid_size
INTEGER_T adapt_factor(SDIM)
REAL_T buffer_factor(SDIM)
INTEGER_T clo(SDIM)
INTEGER_T chi(SDIM)
INTEGER_T levelcomp,lvaptype,lgastype,lsoltype,level,gridno
INTEGER_T grids_per_level,finest_level,nsteps,arrdim
REAL_T dxs ! uniform grid spacing of the local auxiliary 
           ! grid: should be use input


REAL_T, pointer :: masklocal(:,:,:)

#include "probdataf95.H"

 if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
 endif
 local_sci_grid_size=SCIGRIDSIZE
 do dir=1,SDIM
  adapt_factor(dir)=1
  buffer_factor(dir)=20.0  ! buffer is 1/20 the size of the interior
 enddo
 if (probtype.eq.5700) then  ! microfluidics simulations
  if ((axis_dir.eq.0).or.(axis_dir.eq.1).or.(axis_dir.eq.3).or. &
      (axis_dir.eq.4)) then
   local_sci_grid_size=128.0
   adapt_factor(3)=4
  else if (axis_dir.eq.2) then
   local_sci_grid_size=128.0
   adapt_factor(3)=3
   adapt_factor(2)=3
   buffer_factor(3)=4.0
  else
   if (isout.eq.1) then
    print *,"axis_dir invalid"
   endif
   stop
  endif
 else if (probtype.eq.563) then
  local_sci_grid_size=384.0  ! gear
 else if (probtype.eq.5602) then
  local_sci_grid_size=64.0  ! internal inflow
 else if (probtype.eq.5601) then
  local_sci_grid_size=64.0  ! viorel sphere
 else if (probtype.eq.537) then 
  local_sci_grid_size=128.0  ! injector 128 for production
 else if ((probtype.eq.53).and.(axis_dir.eq.100)) then 
  local_sci_grid_size=96.0  ! "fan" injector
 else if (probtype.eq.5501) then
  local_sci_grid_size=96.0  ! rough surface (0.08^3 domain)
!  local_sci_grid_size=256.0  ! rough surface (0.64^3 domain)
 else if (probtype.eq.538) then
  if (MARCO.eq.0) then
!  local_sci_grid_size=144.0    ! injector with moving parts
   local_sci_grid_size=12.0     ! injector with moving parts
  else
!  local_sci_grid_size=196.0    ! injector with moving parts
   local_sci_grid_size=34.0    ! injector with moving parts
  endif

 else if (probtype.eq.539) then
! local_sci_grid_size=128    ! sup injector
! local_sci_grid_size=724    ! sup injector
  local_sci_grid_size=762    ! sup injector
 else if (probtype.eq.701) then ! flapping wing
  local_sci_grid_size=256
 endif

 h_buffer=0.8

 if ((initflag.eq.1).and.(isout.eq.1)) then
  if (ioproc.eq.1) then
   print *,"in generate_new_triangles..."
  endif
 endif

 call GetDomainScale(maxnodedomain,minnodedomain,mindist, &
   maxnoderadius,COM,FSI_in);
 maxdist=0.0
 do dir=1,SDIM
  SCI_in%Cellproblo(dir)=minnodedomain(dir)
  SCI_in%Cellprobhi(dir)=maxnodedomain(dir)

  ! ysl 09/02/14 should we comment off the following commands?
  !if(FSI_in%PartID.eq.2) then                      
  ! ! use hack to keep part 2 nested inside part 1
  ! if (islocal.eq.0) then
  !  SCI_in%Cellproblo(dir)=SCI%Cellproblo(dir)
  !  SCI_in%Cellprobhi(dir)=SCI%Cellprobhi(dir)
  ! else
  !   ! part 1 local auxgrid already generated.
  !  SCI_in%Cellproblo(dir)=SCI_loc%Cellproblo(dir)
  !  SCI_in%Cellprobhi(dir)=SCI_loc%Cellprobhi(dir)
  ! endif
  !endif

   ! problo and probhi might be a subset of the overall domain (islocal=1)
  if (SCI_in%Cellproblo(dir).lt.problo(dir)) then
   SCI_in%Cellproblo(dir)=problo(dir)
  endif
  if (SCI_in%Cellprobhi(dir).gt.probhi(dir)) then
   SCI_in%Cellprobhi(dir)=probhi(dir)
  endif
  if (SCI_in%Cellprobhi(dir).le.SCI_in%Cellproblo(dir)) then
   ! possible if the solid body is outside of the local domain
   if(islocal.eq.1) then
    print *,"solid body not included in the sub-domain, so no local LS"
    return
   else
    print *,"solid body not included in the whole domain: exiting ..."
    stop
   endif
  endif
 enddo

 maxdist=0.0
 do dir=1,SDIM
  diag=(SCI_in%Cellprobhi(dir)-SCI_in%Cellproblo(dir))/buffer_factor(dir)
  SCI_in%Cellproblo(dir)=SCI_in%Cellproblo(dir)-diag
  SCI_in%Cellprobhi(dir)=SCI_in%Cellprobhi(dir)+diag
  diag=SCI_in%Cellprobhi(dir)-SCI_in%Cellproblo(dir)
  if (diag.gt.maxdist) then
   maxdist=diag
  endif
 enddo

 h_small=1.0E+10
 if (MARCO.eq.0) then
  dxs = 6.0e-3  ! 6.0e-3 looks ok
 else
  dxs = 10e-4
 endif
 do dir=1,SDIM
  diag=SCI_in%Cellprobhi(dir)-SCI_in%Cellproblo(dir)
     ! NINT->round to nearest whole number
  if (islocal.eq.0) then !based on sci_grid_size for global auxiliary grid
   n_cell(dir)=NINT(adapt_factor(dir)*local_sci_grid_size*diag/maxdist)
  else if (islocal.eq.1) then ! based on specified dxs for local grid
   n_cell(dir)=NINT(diag/dxs)+1
  else
   print*,"generate_new_triangles: islocal bust ",islocal
   stop
  endif
  SCI_in%cell_lo(dir)=0
  SCI_in%cell_hi(dir)=n_cell(dir)-1
  SCI_in%celldx(dir)=diag/n_cell(dir)
  if (SCI_in%celldx(dir).lt.h_small) then
   h_small=SCI_in%celldx(dir)
  endif
 enddo ! dir

  !output when global grid and output processor or when local grid
 if (((ioproc.eq.1).and.(isout.eq.1)).or.(islocal.eq.1)) then

  print *,"COM: ",COM(1),COM(2),COM(3)
  print *,"maxnoderadius ",maxnoderadius
  print *,"mindist (scale of geom), h_small ",mindist,h_small
  do dir=1,SDIM
   print *,"islocal,dir,n_cell,dx,problo,probhi ",islocal,dir,n_cell(dir), &
    SCI_in%celldx(dir),SCI_in%Cellproblo(dir), &
    SCI_in%Cellprobhi(dir)
  enddo

 endif

 biggest_h=0.0
 new_NumIntElems=0
 do ielem=1,FSI_in%NumIntElems
  nodes_per_elem=FSI_in%ElemData(1,ielem)
  do isub=0,nodes_per_elem-3
   new_NumIntElems=new_NumIntElems+1
   node1=FSI_in%IntElem(1,ielem)
   node2=FSI_in%IntElem(nodes_per_elem-isub-1,ielem)
   node3=FSI_in%IntElem(nodes_per_elem-isub,ielem)
   do dir=1,3
    x1(dir)=FSI_in%Node(dir,node1)
    x2(dir)=FSI_in%Node(dir,node2)
    x3(dir)=FSI_in%Node(dir,node3)
   enddo
   call xdist(x1,x2,d1)
   call xdist(x2,x3,d2)
   call xdist(x3,x1,d3)
   if (biggest_h.lt.d1) then
    biggest_h=d1
   endif
   if (biggest_h.lt.d2) then
    biggest_h=d2
   endif
   if (biggest_h.lt.d3) then
    biggest_h=d3
   endif
  enddo ! isub
 enddo ! ielem
 
 lag_dim=1
 temp_h=biggest_h
 do while (temp_h.gt.h_buffer*h_small)
  temp_h=0.5*temp_h
  lag_dim=lag_dim+2  ! two sides might exceed limit
 enddo
 lag_dim=2

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"original biggest_h ",biggest_h
  print *,"number of refinements needed for nodes/elements ",lag_dim
 endif

 allocate(multi_lag(lag_dim))
 multi_lag(1)%n_nodes=FSI_in%NumNodes
 multi_lag(1)%n_elems=new_NumIntElems
 allocate(multi_lag(1)%nd(3,FSI_in%NumNodes))
 allocate(multi_lag(1)%ndvel(3,FSI_in%NumNodes))
 allocate(multi_lag(1)%ndtemp(FSI_in%NumNodes))
 allocate(multi_lag(1)%elemdt(3,new_NumIntElems))
 allocate(multi_lag(1)%intelemdt(4,new_NumIntElems))
 do i=1,FSI_in%NumNodes
  do dir=1,3
   multi_lag(1)%nd(dir,i)=FSI_in%Node(dir,i)
   multi_lag(1)%ndvel(dir,i)=FSI_in%NodeVel(dir,i)
  enddo
  multi_lag(1)%ndtemp(i)=FSI_in%NodeTemp(i)
 enddo
  
 new_NumIntElems=0
 do ielem=1,FSI_in%NumIntElems
  nodes_per_elem=FSI_in%ElemData(1,ielem)
  do isub=0,nodes_per_elem-3
   new_NumIntElems=new_NumIntElems+1
   multi_lag(1)%intelemdt(1,new_NumIntElems)=FSI_in%IntElem(1,ielem)
   multi_lag(1)%intelemdt(2,new_NumIntElems)= &
     FSI_in%IntElem(nodes_per_elem-isub-1,ielem)
   multi_lag(1)%intelemdt(3,new_NumIntElems)= &
     FSI_in%IntElem(nodes_per_elem-isub,ielem)
   multi_lag(1)%intelemdt(4,new_NumIntElems)=ielem
   do dir=1,3
    multi_lag(1)%elemdt(dir,new_NumIntElems)=FSI_in%ElemData(dir,ielem)
   enddo
   multi_lag(1)%elemdt(1,new_NumIntElems)=3
  enddo
 enddo

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"numintelems,numnodes ",FSI_in%NumIntElems,FSI_in%NumNodes
  print *,"numintelems(3node),numnodes ",new_NumIntElems,FSI_in%NumNodes
 endif

 do ilevel=1,lag_dim-1
  do iter=1,2
   if (iter.eq.2) then
    save_n_elems=multi_lag(ilevel+1)%n_elems
    save_n_nodes=multi_lag(ilevel+1)%n_nodes 
    allocate(multi_lag(ilevel+1)%nd(3,save_n_nodes))
    allocate(multi_lag(ilevel+1)%ndvel(3,save_n_nodes))
    allocate(multi_lag(ilevel+1)%ndtemp(save_n_nodes))
    allocate(multi_lag(ilevel+1)%elemdt(3,save_n_elems))
    allocate(multi_lag(ilevel+1)%intelemdt(4,save_n_elems))
    do i=1,multi_lag(ilevel)%n_nodes
     do dir=1,3
      multi_lag(ilevel+1)%nd(dir,i)=multi_lag(ilevel)%nd(dir,i)
      multi_lag(ilevel+1)%ndvel(dir,i)=multi_lag(ilevel)%ndvel(dir,i)
     enddo
     multi_lag(ilevel+1)%ndtemp(i)=multi_lag(ilevel)%ndtemp(i)
    enddo
   endif  ! iter.eq.2

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
    temp1=multi_lag(ilevel)%ndtemp(node1)
    temp2=multi_lag(ilevel)%ndtemp(node2)
    temp3=multi_lag(ilevel)%ndtemp(node3)
    call xdist(x1,x2,d1)
    call xdist(x2,x3,d2)
    call xdist(x3,x1,d3)

    if ((d1.ge.h_buffer*h_small).and.(d1.ge.d2).and.(d1.ge.d3)) then
     do dir=1,3
      xsplit(dir)=0.5*(x1(dir)+x2(dir))
      velsplit(dir)=0.5*(vel1(dir)+vel2(dir))
     enddo
     tempsplit=0.5*(temp1+temp2)

     multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel+1)%n_nodes+1
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+2
     nsplit=multi_lag(ilevel+1)%n_nodes
     esplit=multi_lag(ilevel+1)%n_elems

     if (iter.eq.2) then
      do dir=1,3
       multi_lag(ilevel+1)%nd(dir,nsplit)=xsplit(dir)
       multi_lag(ilevel+1)%ndvel(dir,nsplit)=velsplit(dir)
      enddo 
      multi_lag(ilevel+1)%ndtemp(nsplit)=tempsplit

      base_ielem=multi_lag(ilevel)%intelemdt(4,ielem)
      multi_lag(ilevel+1)%intelemdt(1,esplit-1)=node1
      multi_lag(ilevel+1)%intelemdt(2,esplit-1)=nsplit
      multi_lag(ilevel+1)%intelemdt(3,esplit-1)=node3
      multi_lag(ilevel+1)%intelemdt(4,esplit-1)=base_ielem
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit-1)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo

      multi_lag(ilevel+1)%intelemdt(1,esplit)=nsplit
      multi_lag(ilevel+1)%intelemdt(2,esplit)=node2
      multi_lag(ilevel+1)%intelemdt(3,esplit)=node3
      multi_lag(ilevel+1)%intelemdt(4,esplit)=base_ielem
      do dir=1,3
       multi_lag(ilevel+1)%elemdt(dir,esplit)= &
         multi_lag(ilevel)%elemdt(dir,ielem)
      enddo
     endif ! iter.eq.2
    else if ((d2.ge.h_buffer*h_small).and.(d2.ge.d1).and.(d2.ge.d3)) then
     do dir=1,3
      xsplit(dir)=0.5*(x2(dir)+x3(dir))
      velsplit(dir)=0.5*(vel2(dir)+vel3(dir))
     enddo
     tempsplit=0.5*(temp2+temp3)

     multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel+1)%n_nodes+1
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+2
     nsplit=multi_lag(ilevel+1)%n_nodes
     esplit=multi_lag(ilevel+1)%n_elems

     if (iter.eq.2) then
      do dir=1,3
       multi_lag(ilevel+1)%nd(dir,nsplit)=xsplit(dir)
       multi_lag(ilevel+1)%ndvel(dir,nsplit)=velsplit(dir)
      enddo 
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
    else if ((d3.ge.h_buffer*h_small).and.(d3.ge.d1).and.(d3.ge.d2)) then
     do dir=1,3
      xsplit(dir)=0.5*(x1(dir)+x3(dir))
      velsplit(dir)=0.5*(vel1(dir)+vel3(dir))
     enddo
     tempsplit=0.5*(temp1+temp3)

     multi_lag(ilevel+1)%n_nodes=multi_lag(ilevel+1)%n_nodes+1
     multi_lag(ilevel+1)%n_elems=multi_lag(ilevel+1)%n_elems+2
     nsplit=multi_lag(ilevel+1)%n_nodes
     esplit=multi_lag(ilevel+1)%n_elems

     if (iter.eq.2) then
      do dir=1,3
       multi_lag(ilevel+1)%nd(dir,nsplit)=xsplit(dir)
       multi_lag(ilevel+1)%ndvel(dir,nsplit)=velsplit(dir)
      enddo 
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

 FSI_in%NumNodesBIG=multi_lag(lag_dim)%n_nodes
 FSI_in%NumIntElemsBIG=multi_lag(lag_dim)%n_elems

 if (initflag.eq.0) then 
  deallocate(FSI_in%NodeBIG)
  deallocate(FSI_in%NodeVelBIG)
  deallocate(FSI_in%NodeTempBIG)
  deallocate(FSI_in%NodeNormalBIG)
  deallocate(FSI_in%ElemNodeCount)

  deallocate(FSI_in%ElemDataXnotBIG)
  deallocate(FSI_in%ElemDataBIG)
  deallocate(FSI_in%IntElemBIG)

  deallocate(SCI_in%CellLS)
  deallocate(SCI_in%CellVEL)
  deallocate(SCI_in%CellTEMP)
 endif

 allocate(FSI_in%NodeBIG(3,FSI_in%NumNodesBIG))
 allocate(FSI_in%NodeNormalBIG(3,FSI_in%NumNodesBIG))
 allocate(FSI_in%NodeVelBIG(3,FSI_in%NumNodesBIG))
 allocate(FSI_in%NodeTempBIG(FSI_in%NumNodesBIG))
 allocate(FSI_in%ElemNodeCount(FSI_in%NumNodesBIG))
 allocate(FSI_in%ElemDataXnotBIG(3,FSI_in%NumIntElemsBIG))
 allocate(FSI_in%ElemDataBIG(3,FSI_in%NumIntElemsBIG))
 allocate(FSI_in%IntElemBIG(4,FSI_in%NumIntElemsBIG))

 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"NumNodes, NumIntElems ",FSI_in%NumNodes,FSI_in%NumIntElems
  print *,"NumNodesBIG, NumIntElemsBIG ",FSI_in%NumNodesBIG, &
   FSI_in%NumIntElemsBIG
 endif

 do ielem=1,FSI_in%NumIntElemsBIG
  do dir=1,4
   FSI_in%IntElemBIG(dir,ielem)=multi_lag(lag_dim)%intelemdt(dir,ielem)
  enddo
  do dir=1,3
   FSI_in%ElemDataBIG(dir,ielem)=multi_lag(lag_dim)%elemdt(dir,ielem)
  enddo
 enddo
 do inode=1,FSI_in%NumNodesBIG
  do dir=1,3
   FSI_in%NodeNormalBIG(dir,inode)=0.0
   FSI_in%NodeBIG(dir,inode)=multi_lag(lag_dim)%nd(dir,inode)
   FSI_in%NodeVelBIG(dir,inode)=multi_lag(lag_dim)%ndvel(dir,inode)
  enddo
  FSI_in%NodeTempBIG(inode)=multi_lag(lag_dim)%ndtemp(inode)
  FSI_in%ElemNodeCount(inode)=0
 enddo

 do i=1,lag_dim
  deallocate(multi_lag(i)%intelemdt)
  deallocate(multi_lag(i)%elemdt)
  deallocate(multi_lag(i)%nd)
  deallocate(multi_lag(i)%ndvel)
  deallocate(multi_lag(i)%ndtemp)
 enddo
 deallocate(multi_lag)
 
 if ((ioproc.eq.1).and.(isout.eq.1)) then 
  print *,"creating normals for the nodes and Xnot"
 endif

 do ielem=1,FSI_in%NumIntElemsBIG
  nodes_per_elem=FSI_in%ElemDataBIG(1,ielem)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem invalid"
   stop
  endif
  do dir=1,3
   FSI_in%ElemDataXnotBIG(dir,ielem)=0.0
  enddo
  call scinormal(ielem,normal,FSI_in)
  do i=1,3
   inode=FSI_in%IntElemBIG(i,ielem)
   do dir=1,3
    FSI_in%NodeNormalBIG(dir,inode)=FSI_in%NodeNormalBIG(dir,inode)+normal(dir)
    FSI_in%ElemDataXnotBIG(dir,ielem)=FSI_in%ElemDataXnotBIG(dir,ielem)+ &
      FSI_in%NodeBIG(dir,inode)
   enddo
   FSI_in%ElemNodeCount(inode)=FSI_in%ElemNodeCount(inode)+1
  enddo
  do dir=1,3
   FSI_in%ElemDataXnotBIG(dir,ielem)=FSI_in%ElemDataXnotBIG(dir,ielem)/3.0
  enddo
 enddo

 do inode=1,FSI_in%NumNodesBIG
  normal_cnt=FSI_in%ElemNodeCount(inode)
  if (normal_cnt.gt.0) then
   do dir=1,3
    FSI_in%NodeNormalBIG(dir,inode)=FSI_in%NodeNormalBIG(dir,inode)/normal_cnt
   enddo
  endif
 enddo

 biggest_h=0.0
 do ielem=1,FSI_in%NumIntElemsBIG
  nodes_per_elem=FSI_in%ElemDataBIG(1,ielem)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem invalid"
   stop
  endif
  do i=1,2
   node1=FSI_in%IntElemBIG(i,ielem)
   node2=FSI_in%IntElemBIG(i+1,ielem)
   do dir=1,3
    x1(dir)=FSI_in%NodeBIG(dir,node1) 
    x2(dir)=FSI_in%NodeBIG(dir,node2)
   enddo
   call xdist(x1,x2,temp_h)
   if (temp_h.gt.biggest_h) then
    biggest_h=temp_h
   endif
  enddo
 enddo
 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"biggest_h=",biggest_h
  print *,"h_buffer ",h_buffer
  print *,"h_small ",h_small
 endif

 call init_dimensions(SCI_in)
 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print*,"in allocate: part ID ",FSI_in%PartID
 endif

 allocate(SCI_in%CellVEL(DIMV3D(CellLS),SDIM))
 allocate(SCI_in%CellLS(DIMV3D(CellLS)))
 allocate(SCI_in%CellLSFACE(DIMV3D(CellLS),SDIM))
 allocate(SCI_in%CellTEMP(DIMV3D(CellLS)))

 allocate(SCI_in%CellMASK(DIMV3D(CellLS)))

 do i=1,SDIM
  clo(i)=SCI_in%cell_lo(i)
  chi(i)=SCI_in%cell_hi(i)
 enddo

 do i=clo(1)-1,chi(1)+1
 do j=clo(2)-1,chi(2)+1
 do k=clo(3)-1,chi(3)+1
  SCI_in%CellMASK(i,j,k)=0.0
  SCI_in%CellLS(i,j,k)=0.0
  SCI_in%CellTEMP(i,j,k)=0.0
  do dir=1,SDIM
   SCI_in%CellVEL(i,j,k,dir)=0.0
   SCI_in%CellLSFACE(i,j,k,dir)=1.0E+10  ! height fraction and slope
  enddo
 enddo
 enddo
 enddo

  ! mask=1 velocity is init, but sign is not
  ! mask=2 both sign and velocity are init
  ! mask=3 both sign and velocity are init doubly wetted
  ! mask=0 neither velocity or sign are valid.  LS has a large value.
 call CLSVOF_InitBox( &
   SCI_in%CellMASK,DIMS3D(CellLS), &
   SCI_in%Cellproblo,SCI_in%celldx, &
   SCI_in%CellLS,DIMS3D(CellLS), &
   SCI_in%CellLSFACE,DIMS3D(CellLS), &
   SCI_in%CellVEL,DIMS3D(CellLS), &
   SCI_in%CellTEMP,DIMS3D(CellLS), &
   SCI_in%cell_lo,SCI_in%cell_hi, &
   FSI_in,ioproc)


  ! neumann bc
 call simpleBC(SCI_in%CellMASK,DIMS3D(CellLS),1,clo,chi)
 call simpleBC(SCI_in%CellLS,DIMS3D(CellLS),1,clo,chi)
 call simpleBC(SCI_in%CellTEMP,DIMS3D(CellLS),1,clo,chi)
 call simpleBC(SCI_in%CellVEL,DIMS3D(CellLS),SDIM,clo,chi)
 
 call sign_extrap_iter(SCI_in%CellLS,SCI_in%CellLSFACE, &
   DIMS3D(CellLS),clo,chi,FSI_in,SCI_in, &
   SCI_in%Cellproblo,SCI_in%celldx,ioproc)

 scomp=1
 call extrap_iter(SCI_in%CellTEMP,DIMS3D(CellLS),scomp,clo,chi,SCI_in)
 do dir=1,SDIM
  call extrap_iter(SCI_in%CellVEL,DIMS3D(CellLS),dir,clo,chi,SCI_in)
 enddo 

  ! debugging isosurface
#if (BL_SPACEDIM==3)

 if (DEBUG_ISOSURFACE.eq.1) then
  print *,"writing solid LS isosurface from separate buffer"
  level=0
  gridno=0
  allocate(masklocal(DIMV3D(CellLS)))

   ! clo(i)=SCI_in%cell_lo(i)
   ! chi(i)=SCI_in%cell_hi(i)
  do i=clo(1)-1,chi(1)+1
  do j=clo(2)-1,chi(2)+1
  do k=clo(3)-1,chi(3)+1
   masklocal(i,j,k)=one
  enddo
  enddo
  enddo
  call FORT_ISOGRIDSINGLE( &
   SCI_in%CellLS,DIMS3D(CellLS), & 
   SCI_in%Cellproblo,SCI_in%celldx, &
   masklocal,DIMS3D(CellLS), &
   clo,chi, &
   level,gridno)

  grids_per_level=1
  finest_level=0
  nsteps=0
  arrdim=1

  call FORT_COMBINETRIANGLESSINGLE( &
   grids_per_level,finest_level,nsteps,arrdim)

  deallocate(masklocal)
  print *,"end writing solid LS isosurface from separate buffer"
 endif

#elif (BL_SPACEDIM==2)

 if (DEBUG_ISOSURFACE.eq.1) then
  print *,"cannot debug isosurface in 2d"
  stop
 endif

#else

 print *,"dimension bust"
 stop

#endif
 
 deallocate(SCI_in%CellMASK)
 deallocate(SCI_in%CellLSFACE)

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

subroutine initinjector(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  FSI_read,isout)

type(mesh_type) :: FSI_read
INTEGER_T :: j,k,i1,j1,k1,itimecount,sdim,ifirst,isout
INTEGER_T :: inode,iface,ioproc
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1,xval2,xtemp
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
REAL_T :: tper,tcrit,theta
INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
INTEGER_T :: filler
character(40) :: dwave,dwave2
REAL_T :: ofs,xx,zz

REAL_T, dimension(3) :: xxblob1,newxxblob1,xxblob2,newxxblob2
REAL_T :: radradblob1,radradblob2

#include "probdataf95.H"

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
  endif

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=0.0
  TorquePos=0.0
  TorqueVel=0.0

  iread=0
  if (ifirst.eq.1) then
   allocate(HdrFlag(NumFlags(1)))
   allocate(NptFlag(1))
   allocate(ForFlag(NumFlags(4))) 
   iread=1
  else if (ifirst.eq.0) then
   iread=0
  else
   print *,"ifirst invalid"
   stop
  endif

  if (iread.eq.1) then
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
   else if ((probtype.eq.536).or.(probtype.eq.537).or. &
            (probtype.eq.538).or.(probtype.eq.539)) then

      ! MARK:
      ! IF probtype=538 and partID=2, then add 0.01 in order to shift the
      ! needle all the way in the domain (newxxblob1(3)=0.01)
      ! this shift will be taken into account later.
    xxblob1(1)=0.0
    xxblob1(2)=0.0
    xxblob1(3)=0.0
    newxxblob1(1)=0.0
    newxxblob1(2)=0.0
    newxxblob1(3)=0.0

    if (probtype.eq.538) then
     if (FSI_read%partID.eq.2) then  ! want the whole needle in the domain.
      newxxblob1(3)=0.01
     else if (FSI_read%partID.eq.1) then
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
   if (FSI_read%partID.eq.1 ) then

     if (probtype.eq.5501) then
      dwave="rough.cas"
     else if ((probtype.eq.53).and.(axis_dir.eq.100)) then
      dwave="flat_fan_s.cas"
     else
      dwave="injectorgeom.dat"
     endif

   else if (FSI_read%partID.eq.2) then
     dwave="injectorgeom_needle.dat"
   else
     print *,"part id invalid"
     stop
   endif

   if ((ioproc.eq.1).and.(isout.eq.1)) then
    print *,"opening ",dwave
   endif
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

   READ(14,*) FSI_read%NumNodes,FSI_read%NumIntElems
   FSI_read%NumNodes=FSI_read%NumNodes*2
   FSI_read%NumIntElems=FSI_read%NumIntElems*2
   FSI_read%IntElemDim=3
 
   if (ifirst.eq.1) then
    call init_FSI(FSI_read,1)  ! allocate_intelem=1
   else
    print *,"something wrong, ifirst should be 1 here"
    stop
   endif  ! ifirst.eq.1

   HdrFlag(1)=sdim
   HdrFlag(2)=0     ! stationary body

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   do inode=1,FSI_read%NumNodes/2
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
     FSI_read%Node_old(dir,inode)=xval1(dir)
     FSI_read%Node_new(dir,inode)=xval1(dir)
     FSI_read%Node_old(dir,inode+FSI_read%NumNodes/2)=xval2(dir)
     FSI_read%Node_new(dir,inode+FSI_read%NumNodes/2)=xval2(dir)
    enddo
      
   enddo  ! inode=1,NumNodes
    
   do iface=1,FSI_read%NumIntElems/2
    READ(14,*) FSI_read%IntElem(1,iface),FSI_read%IntElem(2,iface), &
               FSI_read%IntElem(3,iface)
    do dir=3,1,-1
     FSI_read%IntElem(dir,iface+FSI_read%NumIntElems/2)= &
       FSI_read%IntElem(dir,iface)+FSI_read%NumNodes/2 
    enddo
   enddo

   do iface=1,FSI_read%NumIntElems
    FSI_read%ElemData(1,iface)=3   ! number of nodes in element
    FSI_read%ElemData(2,iface)=1   ! part number
    FSI_read%ElemData(3,iface)=0   ! singly wetted (1=doubly wetted)
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

   call init2_FSI(FSI_read) 
  else if (iread.eq.0) then
   ! do nothing
  else
   print *,"iread invalid"
   stop
  endif 

  use_temp=0

    ! this routine initializes FSI_read%solid_displ and
    ! FSI_read%solid_speed
  call init3_FSI(FSI_read,ifirst,1,ioproc,isout) ! do_2nd_part=1 (initinjector)
  
191    FORMAT(I12,I12)
91     FORMAT(I7)
92     FORMAT(f11.3)
193    FORMAT(E15.11,E15.11,E15.11)
194    FORMAT(I12,I12,I12,I12)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine initinjector



subroutine initflapping(curtime,dt,ifirst,sdim,istop,istep,ioproc, &
  FSI_read,isout)

type(mesh_type) :: FSI_read
INTEGER_T :: j,k,i1,j1,k1,itimecount,sdim,ifirst,isout
INTEGER_T :: inode,iface,ioproc
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1,xval2,xtemp
REAL_T, dimension(3) :: xvalbefore
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
REAL_T :: tper,tcrit,theta
INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
INTEGER_T :: filler
character(40) :: dwave,dwave2
REAL_T :: ofs,xx,zz

REAL_T, dimension(3) :: xxblob1,newxxblob1,xxblob2,newxxblob2
REAL_T :: radradblob1,radradblob2
INTEGER_T tempelem1
INTEGER_T tempelem2
INTEGER_T tempelem3
INTEGER_T tempelem4
INTEGER_T numquads,quad_counter
REAL_T localElem(4)

#include "probdataf95.H"

  if (probtype.ne.701) then
   print *,"probtype invalid"
   stop
  endif

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
  endif

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=0.0
  TorquePos=0.0
  TorqueVel=0.0

  iread=0
  if (ifirst.eq.1) then
   allocate(HdrFlag(NumFlags(1)))
   allocate(NptFlag(1))
   allocate(ForFlag(NumFlags(4))) 
   iread=1
  else if (ifirst.eq.0) then
   iread=0
  else
   print *,"ifirst invalid"
   stop
  endif

  if (iread.eq.1) then
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

   if (FSI_read%partID.eq.1 ) then

    dwave="foreWing.sci"
    newxxblob1(1)=0.0  ! 1st wing at a "default" location.
    print *,"init_flapping partID=1"

   else if (FSI_read%partID.eq.2) then

    if (axis_dir.eq.0) then
     print *,"part id invalid"
     stop
    else if (axis_dir.eq.1) then
     dwave="hindWing.sci"
     newxxblob1(1)=0.0  
     print *,"init_flapping partID=2"
    else
     print *,"axis_dir invalid"
     stop
    endif
   else
     print*,"invalid partID"
     stop
   endif

   numquads=0

   print *,"opening to check to see how many quads to convert ",dwave 
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
   READ(14,*) FSI_read%NumNodes
   print *,"NumNodes(check) ",FSI_read%NumNodes
   READ(14,*) FSI_read%NumIntElems
   print *,"NumIntElems(check) ",FSI_read%NumIntElems
   do inode=1,FSI_read%NumNodes
    READ(14,*) xvalbefore(1),xvalbefore(2),xvalbefore(3)
   enddo
   do iface=1,FSI_read%NumIntElems
    numquads=numquads+1
    READ(14,*) tempelem1,tempelem2,tempelem3,tempelem4
   enddo
   close(14)
   print *,"numquads=",numquads

   print *,"opening ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

   READ(14,*) FSI_read%NumNodes
   print *,"NumNodes ",FSI_read%NumNodes
   READ(14,*) FSI_read%NumIntElems
   print *,"NumIntElems ",FSI_read%NumIntElems
   FSI_read%IntElemDim=3
   FSI_read%NumIntElems=FSI_read%NumIntElems+numquads

   if (ifirst.ne.1) then
    print *,"ifirst bust"
    stop
   endif

   call init_FSI(FSI_read,1)  ! allocate_intelem=1

   HdrFlag(1)=sdim
   HdrFlag(2)=0     ! stationary body

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   do inode=1,FSI_read%NumNodes
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
     FSI_read%Node_old(dir,inode)=xval1(dir)
     FSI_read%Node_new(dir,inode)=xval1(dir)
    enddo
      
   enddo  ! inode=1,NumNodes
   
   quad_counter=0 
   do iface=1,FSI_read%NumIntElems-numquads
    READ(14,*) localElem(4),localElem(3),localElem(2),localElem(1)
    do j1=0,1
     if (j1.eq.1) then
      quad_counter=quad_counter+1
     endif
     FSI_read%ElemData(1,iface+quad_counter)=3   ! number of nodes in element
     FSI_read%ElemData(2,iface+quad_counter)=1   ! part number
     FSI_read%ElemData(3,iface+quad_counter)=0   ! singly wetted
     if (j1.eq.0) then
      do dir=1,3
       FSI_read%IntElem(dir,iface+quad_counter)=localElem(dir)
      enddo
     else if (j1.eq.1) then
       FSI_read%IntElem(1,iface+quad_counter)=localElem(3)
       FSI_read%IntElem(2,iface+quad_counter)=localElem(4)
       FSI_read%IntElem(3,iface+quad_counter)=localElem(1)
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

   call init2_FSI(FSI_read)
  else if (iread.eq.0) then
   ! do nothing
  else
   print *,"iread invalid"
   stop
  endif 

  use_temp=0
   ! do_2nd_part=0  
   ! even if 2 wings, we want do_2nd_part=0 because
   ! the velocity is prescribed in the kinematic routines later on.
  call init3_FSI(FSI_read,ifirst,0,ioproc,isout) 

191    FORMAT(I12,I12)
91     FORMAT(I7)
92     FORMAT(f11.3)
193    FORMAT(E15.11,E15.11,E15.11)
194    FORMAT(I12,I12,I12,I12)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine initflapping




subroutine initchannel(curtime,dt,ifirst,sdim,istop,istep)
INTEGER_T :: j,k,i1,j1,k1,itimecount,sdim,ifirst
INTEGER_T :: inode,iface
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xval1,xval2,xtemp
REAL_T, dimension(3) :: maxnodebefore,minnodebefore,xvalbefore
REAL_T :: tper,tcrit,theta
INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
INTEGER_T :: numquads,quad_counter
INTEGER_T :: tempelem1,tempelem2,tempelem3,tempelem4
INTEGER_T, dimension(4) :: localElem
character(40) :: dwave
REAL_T :: ofs,xx,zz

REAL_T, dimension(3) :: xxblob1,newxxblob1
REAL_T :: radradblob1

#include "probdataf95.H"

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=0.0
  TorquePos=0.0
  TorqueVel=0.0

  iread=0
  if (ifirst.eq.1) then
   allocate(HdrFlag(NumFlags(1)))
   allocate(NptFlag(1))
   allocate(ForFlag(NumFlags(4))) 
   iread=1
  endif

  if (probtype.ne.5700) then
   print *,"probtype invalid"
   stop
  endif
  if (iread.eq.1) then
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
   READ(14,*) FSI%NumNodes
   print *,"NumNodes(check) ",FSI%NumNodes
   READ(14,*) FSI%NumIntElems
   print *,"NumIntElems(check) ",FSI%NumIntElems
   do inode=1,FSI%NumNodes
    READ(14,*) xvalbefore(1),xvalbefore(2),xvalbefore(3)
   enddo
   do iface=1,FSI%NumIntElems
    numquads=numquads+1
    READ(14,*) tempelem1,tempelem2,tempelem3,tempelem4
   enddo
   close(14)
   print *,"numquads=",numquads

   print *,"opening ",dwave
   OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

   READ(14,*) FSI%NumNodes
   print *,"NumNodes ",FSI%NumNodes
   READ(14,*) FSI%NumIntElems
   print *,"NumIntElems ",FSI%NumIntElems
   FSI%IntElemDim=3
   FSI%NumIntElems=FSI%NumIntElems+numquads

   if (ifirst.eq.1) then
    call init_FSI(FSI,1)
   endif  ! ifirst.eq.1

   HdrFlag(1)=sdim
   HdrFlag(2)=0     ! stationary body

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   do inode=1,FSI%NumNodes
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
     FSI%Node_old(dir,inode)=xval1(dir)
     FSI%Node_new(dir,inode)=xval1(dir)
    enddo
      
   enddo  ! inode=1,NumNodes
   
   quad_counter=0 
   do iface=1,FSI%NumIntElems-numquads
    READ(14,*) localElem(1),localElem(2),localElem(3),localElem(4)
    do j1=0,1
     if (j1.eq.1) then
      quad_counter=quad_counter+1
     endif
     FSI%ElemData(1,iface+quad_counter)=3   ! number of nodes in element
     FSI%ElemData(2,iface+quad_counter)=1   ! part number
     FSI%ElemData(3,iface+quad_counter)=0   ! singly wetted
     if (j1.eq.0) then
      do dir=1,3
       FSI%IntElem(dir,iface+quad_counter)=localElem(dir)
      enddo
     else if (j1.eq.1) then
       FSI%IntElem(1,iface+quad_counter)=localElem(3)
       FSI%IntElem(2,iface+quad_counter)=localElem(4)
       FSI%IntElem(3,iface+quad_counter)=localElem(1)
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

   call init2_FSI(FSI)
 
  endif ! iread=1

  use_temp=0
  call init3_FSI(FSI,ifirst,1,1,1)  ! do_2nd_part=1  isout=1 (initchannel)

191    FORMAT(I12,I12)
91     FORMAT(I7)
92     FORMAT(f11.3)
193    FORMAT(E15.11,E15.11,E15.11)
194    FORMAT(I12,I12,I12,I12)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine initchannel




subroutine geominit(curtime,dt,ifirst,sdim,istop,istep)
INTEGER_T :: i,j,k,i1,j1,k1,m1,tempelem,itimecount,ifirst,sdim
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
INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
INTEGER_T :: gwave
INTEGER_T :: shift_from_zero_node
INTEGER_T :: shift_from_zero_face
INTEGER_T :: override_IntElemDim
character(100) :: dwave,poolname
REAL_T :: ofs
REAL_T :: plungerfreq,plungeramp

#include "probdataf95.H"

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=dt

  if (ifirst.eq.1) then
   allocate(HdrFlag(NumFlags(1)))
   allocate(NptFlag(1))
   allocate(ForFlag(NumFlags(4))) 
  endif

  if (probtype.eq.58) then
   iread=1
  else if (probtype.eq.55) then
   if (ifirst.eq.1) then
    iread=1
   else
    iread=0
   endif
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

   FSI%NumNodesPool=0
   FSI%NumIntElemsPool=0
   FSI%IntElemDimPool=0

   do i1=1,2 
    if ((probtype.eq.55).or.(probtype.eq.58)) then
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
     FSI%NumNodesPool=0
     FSI%NumIntElemsPool=0
     FSI%IntElemDimPool=0

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
      READ(141,*) FSI%NumNodesPool
      print *,"NumNodesPool ",FSI%NumNodesPool
      READ(141,*) FSI%NumIntElemsPool
      print *,"NumIntElemsPool ",FSI%NumIntElemsPool
      READ(141,*) FSI%IntElemDimPool
      print *,"IntElemDimPool ",FSI%IntElemDimPool
     else
      FSI%NumNodesPool=0
      FSI%NumIntElemsPool=0
      FSI%IntElemDimPool=0
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
    if (FSI%IntElemDimPool.eq.0) then
     print *,"opening to check to see how many quads to convert ",dwave
     OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
     READ(14,91) FSI%NumNodes
     print *,"NumNodes(check) ",FSI%NumNodes
     READ(14,91) FSI%NumIntElems
     print *,"NumIntElems(check) ",FSI%NumIntElems
     READ(14,91) FSI%IntElemDim
     print *,"IntElemDim(check) ",FSI%IntElemDim
     do i=1,FSI%NumNodes
      READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
     enddo
     READ(14,*) j
     do i=1,FSI%NumIntElems
      READ(14,*) j
      READ(14,*) k
      if (k.eq.FSI%NumIntElems) then
       k=3
      endif
      if (k.eq.4) then
       numquads=numquads+1
      else if (k.ne.3) then
       print *,"elements must be triangles or quadrilaterals k=",k
       print *,"i,NumIntElems ",i,FSI%NumIntElems
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

    READ(14,91) FSI%NumNodes
    print *,"NumNodes ",FSI%NumNodes
    READ(14,91) FSI%NumIntElems
    print *,"NumIntElems ",FSI%NumIntElems
    READ(14,91) FSI%IntElemDim
    print *,"IntElemDim ",FSI%IntElemDim
    print *,"overriding IntElemDim to be 3"
    FSI%IntElemDim=3
    FSI%NumNodes=FSI%NumNodes+FSI%NumNodesPool
    FSI%NumIntElems=FSI%NumIntElems+FSI%NumIntElemsPool

    FSI%NumIntElems=FSI%NumIntElems+numquads

    if (FSI%IntElemDimPool.ne.0) then
     if (FSI%IntElemDim.ne.FSI%IntElemDimPool) then
      print *,"IntElemDim inconsistent"
      stop
     endif
    endif

    if ((ifirst.eq.1).and.(i1.eq.1)) then
     call init_FSI(FSI,1)
    endif

    HdrFlag(1)=sdim
    HdrFlag(2)=0     ! stationary body

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
     do i=1,FSI%NumNodesPool
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
        FSI%Node_old(dir,i)=xval(dir)
       enddo
      else if (i1.eq.2) then
       do dir=1,3
        FSI%Node_new(dir,i)=xval(dir)
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
     if (j.ne.FSI%NumIntElemsPool) then
      print *,"face mismatch reading pool file"
      stop
     endif

     do i=1,FSI%NumIntElemsPool
      READ(141,*) j
      if (i.ne.j) then
       print *,"face mismatch reading pool file"
       stop
      endif
      READ(141,*) k
      if (k.gt.FSI%IntElemDim) then
       print *,"too many vertices k,IntElemDim ",k,FSI%IntElemDim
       stop
      endif
      FSI%ElemData(1,i)=k   ! number of nodes in element
      FSI%ElemData(2,i)=1   ! part number
      FSI%ElemData(3,i)=2   ! singly wetted, but do not call "fill" for these...
      do j1=1,k
       READ(141,*) FSI%IntElem(j1,i) 
      enddo
     enddo  ! i, looping faces
     close(141)
    endif ! probtype.eq.561 


    shift_from_zero_node=0
    shift_from_zero_face=0
    override_IntElemDim=0
    do i=FSI%NumNodesPool+1,FSI%NumNodes
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
       FSI%Node_old(dir,i)=xval(dir)
      enddo
     else if (i1.eq.2) then
      do dir=1,3
       FSI%Node_new(dir,i)=xval(dir)
      enddo
     else
      print *,"i1 invalid"
      stop
     endif
     if (i-FSI%NumNodesPool.ne.j) then
      print *,"vertex mismatch"
      print *,"NumNodesPool=",FSI%NumNodesPool
      print *,"expected node index ",i-FSI%NumNodesPool
      print *,"node index read from file ",j
      stop
     endif
      
    enddo
    READ(14,*) j
    if (j.ne.FSI%NumIntElems-FSI%NumIntElemsPool-numquads) then
     print *,"face mismatch"
     stop
    endif

    quad_counter=0
    do i=1+FSI%NumIntElemsPool,FSI%NumIntElems-numquads
     READ(14,*) j
     if ((j.eq.0).and.(shift_from_zero_face.eq.0)) then
      shift_from_zero_face=1
      print *,"faces in file start at 0; shifting to start at 1"
     endif
     if (shift_from_zero_face.eq.1) then
      j=j+1
     endif

     if (i-FSI%NumIntElemsPool.ne.j) then
      print *,"face mismatch"
      print *,"NumIntElemsPool=",FSI%NumIntElemsPool
      print *,"expected face index ",i-FSI%NumIntElemsPool
      print *,"face index read from file ",j
      stop
     endif
     READ(14,*) k
     if (k.eq.FSI%NumIntElems) then
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
      FSI%ElemData(1,i+quad_counter)=3   ! number of nodes in element
      FSI%ElemData(2,i+quad_counter)=1   ! part number
      FSI%ElemData(3,i+quad_counter)=0   ! singly wetted
      if (probtype.eq.57) then  ! 57 heart    56 swimmer
       FSI%ElemData(3,i+quad_counter)=1   ! doubly wetted
      endif
      if ((probtype.eq.56).and.(1.eq.1)) then  ! BOXSWIMMER
       FSI%ElemData(3,i+quad_counter)=1   ! doubly wetted
      endif
      if (probtype.eq.562) then  ! 562 whale
       FSI%ElemData(3,i+quad_counter)=0   ! 0=singly wetted 1=doubly wetted
       if (axis_dir.eq.6) then
        FSI%ElemData(3,i+quad_counter)=0 
       endif
      endif
      if (probtype.eq.5600) then ! dog
       FSI%ElemData(3,i+quad_counter)=0   ! 0=singly wetted 1=doubly wetted
      endif
      if (j1.eq.0) then
       do dir=1,3
        FSI%IntElem(dir,i+quad_counter)=localElem(dir)
       enddo
      else if (j1.eq.1) then
       FSI%IntElem(1,i+quad_counter)=localElem(3)
       FSI%IntElem(2,i+quad_counter)=localElem(4)
       FSI%IntElem(3,i+quad_counter)=localElem(1)
      else
       print *,"j1 invalid"
       stop
      endif
      do m1=1,3
       FSI%IntElem(m1,i+quad_counter)=FSI%IntElem(m1,i+quad_counter)+ &
         FSI%NumNodesPool
       if (shift_from_zero_node.eq.1) then
        FSI%IntElem(m1,i+quad_counter)=FSI%IntElem(m1,i+quad_counter)+1
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

   if ((probtype.eq.55).or.(probtype.eq.58)) then
    theta=1.0
   else if (tcrit.le.tstart+1.0e-8) then
    theta=0.0
   else if (tcrit.ge.tfinish-1.0e-8) then
    theta=1.0
   else
    theta=(tcrit-tstart)/(tfinish-tstart)
   endif

   do i=1,FSI%NumNodes
    if (probtype.eq.58) then 
     do dir=1,3
      FSI%NodeVel_old(dir,i)=0.0
     enddo
     FSI%NodeVel_old(1,i)=2.0*3.14159*plungeramp*plungerfreq*  &
       cos(2.0*3.14159*plungerfreq*curtime)
    else if (tfinish-tstart.gt.1.0e-8) then
     do dir=1,3
      FSI%NodeVel_old(dir,i)= &
       (FSI%Node_new(dir,i)-FSI%Node_old(dir,i))/(tfinish-tstart)
     enddo
    else
     do dir=1,3
      FSI%NodeVel_old(dir,i)=0.0
     enddo
    endif
    do dir=1,3
     FSI%Node_current(dir,i)=theta*FSI%Node_new(dir,i)+ &
      (1.0-theta)*FSI%Node_old(dir,i)
     FSI%NodeVel_new(dir,i)=FSI%NodeVel_old(dir,i)
    enddo
   enddo  ! i=1,NumNodes
  endif ! iread=1



  use_temp=0
  call init3_FSI(FSI,ifirst,0,1,1) ! do_2nd_part=0 isout=1 (geominit)

  print *,"after geominit  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt,Theta_Dot,Theta ",timeB,dtB,TorqueVel,TorquePos
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
      xxblob(1),xxblob(3),newxxblob(1),newxxblob(3),radradblob
  print *," xnew=cos(theta)*(x-x_0)+sin(theta)*(z-z_0)+x_0"
  print *," znew=-sin(theta)*(x-x_0)+cos(theta)*(z-z_0)+z_0"

91     FORMAT(I7)
92     FORMAT(f11.3)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine geominit



! sphere centered almost at the origin with a radius of 5
! inflow at xlo, outflow at xhi
subroutine viorel_sphere_geominit(curtime,dt,ifirst,sdim,istop,istep)
INTEGER_T :: i,j,k,i1,j1,k1,m1,tempelem,itimecount,ifirst,sdim
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

INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
INTEGER_T :: gwave

INTEGER_T :: shift_from_zero_node
INTEGER_T :: shift_from_zero_face
INTEGER_T :: override_IntElemDim
character(100) :: dwave,poolname
REAL_T :: ofs
REAL_T :: plungerfreq,plungeramp

#include "probdataf95.H"

  if (probtype.ne.5601) then
   print *,"probtype should be 5601"
   stop
  endif

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=dt

  iread=0
  if (ifirst.eq.1) then
   allocate(HdrFlag(NumFlags(1)))
   allocate(NptFlag(1))
   allocate(ForFlag(NumFlags(4))) 
   iread=1
  endif

  if (iread.eq.1) then
   dwave="viorel_sphere.txt"

   FSI%NumNodesPool=0
   FSI%NumIntElemsPool=0
   FSI%IntElemDimPool=0

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
   READ(14,91) FSI%NumNodes
   print *,"NumNodes(check) ",FSI%NumNodes
   READ(14,91) FSI%NumIntElems
   print *,"NumIntElems(check) ",FSI%NumIntElems
   READ(14,91) FSI%IntElemDim
   print *,"IntElemDim(check) ",FSI%IntElemDim
   do i=1,FSI%NumNodes
    READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
   enddo
   READ(14,*) j
   do i=1,FSI%NumIntElems
    READ(14,*) j
    READ(14,*) k
    if (k.eq.FSI%NumIntElems) then ! use default 3 
     k=3
    endif
    if (k.eq.4) then
     numquads=numquads+1
    else if (k.ne.3) then
     print *,"elements must be triangles or quadrilaterals k=",k
     print *,"i,NumIntElems ",i,FSI%NumIntElems
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

   READ(14,91) FSI%NumNodes
   print *,"NumNodes ",FSI%NumNodes
   READ(14,91) FSI%NumIntElems
   print *,"NumIntElems ",FSI%NumIntElems
   READ(14,91) FSI%IntElemDim
   print *,"IntElemDim ",FSI%IntElemDim
   print *,"overriding IntElemDim to be 3"
   FSI%IntElemDim=3
 
   FSI%NumIntElems=FSI%NumIntElems+numquads

   if (ifirst.eq.1) then
    call init_FSI(FSI,1)
   endif

   HdrFlag(1)=sdim
   HdrFlag(2)=0     ! stationary body

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   shift_from_zero_node=0
   shift_from_zero_face=0
   override_IntElemDim=0
   do i=1,FSI%NumNodes
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
     FSI%Node_old(dir,i)=xval(dir)
    enddo
    do dir=1,3
     FSI%Node_new(dir,i)=xval(dir)
    enddo
    if (i.ne.j) then
     print *,"vertex mismatch"
     print *,"NumNodesPool=",FSI%NumNodesPool
     print *,"expected node index ",i-FSI%NumNodesPool
     print *,"node index read from file ",j
     stop
    endif
   enddo !  reading nodes
      
   READ(14,*) j
   if (j.ne.FSI%NumIntElems-numquads) then
    print *,"face mismatch"
    stop
   endif

   quad_counter=0
   do i=1,FSI%NumIntElems-numquads
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
     print *,"NumIntElemsPool=",FSI%NumIntElemsPool
     print *,"expected face index ",i-FSI%NumIntElemsPool
     print *,"face index read from file ",j
     stop
    endif
    READ(14,*) k
    if (k.eq.FSI%NumIntElems) then  ! use default 3
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
     FSI%ElemData(1,i+quad_counter)=3   ! number of nodes in element
     FSI%ElemData(2,i+quad_counter)=1   ! part number
     FSI%ElemData(3,i+quad_counter)=0   ! singly wetted (set =1 doubly)
     if (j1.eq.0) then
      do dir=1,SDIM
       FSI%IntElem(dir,i+quad_counter)=localElem(dir)
      enddo
     else if (j1.eq.1) then
      FSI%IntElem(1,i+quad_counter)=localElem(3)
      FSI%IntElem(2,i+quad_counter)=localElem(4)
      FSI%IntElem(3,i+quad_counter)=localElem(1)
     else
      print *,"j1 invalid"
      stop
     endif
     do m1=1,3
      if (shift_from_zero_node.eq.1) then
       FSI%IntElem(m1,i+quad_counter)=FSI%IntElem(m1,i+quad_counter)+1
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

   call init2_FSI(FSI)
  endif ! iread=1

  use_temp=0 
  call init3_FSI(FSI,ifirst,1,1,1) !do 2nd part=1 isout=1,viorel sphere geominit

  print *,"after viorel_sphere geominit  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt,Theta_Dot,Theta ",timeB,dtB,TorqueVel,TorquePos
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
      xxblob(1),xxblob(3),newxxblob(1),newxxblob(3),radradblob
  print *," xnew=cos(theta)*(x-x_0)+sin(theta)*(z-z_0)+x_0"
  print *," znew=-sin(theta)*(x-x_0)+cos(theta)*(z-z_0)+z_0"

91     FORMAT(I7)
92     FORMAT(f11.3)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine viorel_sphere_geominit




subroutine internal_inflow_geominit(curtime,dt,ifirst,sdim,istop,istep)
INTEGER_T :: i,j,k,i1,j1,k1,m1,tempelem,itimecount,ifirst,sdim
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

INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
INTEGER_T :: gwave

INTEGER_T :: shift_from_zero_node
INTEGER_T :: shift_from_zero_face
INTEGER_T :: override_IntElemDim
character(100) :: dwave,poolname
REAL_T :: ofs
REAL_T :: plungerfreq,plungeramp

#include "probdataf95.H"

  if (probtype.ne.5602) then
   print *,"probtype should be 5602"
   stop
  endif

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=dt

  iread=0
  if (ifirst.eq.1) then
   allocate(HdrFlag(NumFlags(1)))
   allocate(NptFlag(1))
   allocate(ForFlag(NumFlags(4))) 
   iread=1
  endif

  if (iread.eq.1) then
   dwave="internal_inflow.txt"

   FSI%NumNodesPool=0
   FSI%NumIntElemsPool=0
   FSI%IntElemDimPool=0

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
   READ(14,91) FSI%NumNodes
   print *,"NumNodes(check) ",FSI%NumNodes
   READ(14,91) FSI%NumIntElems
   print *,"NumIntElems(check) ",FSI%NumIntElems
   READ(14,91) FSI%IntElemDim
   print *,"IntElemDim(check) ",FSI%IntElemDim
   do i=1,FSI%NumNodes
    READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
   enddo
   READ(14,*) j
   do i=1,FSI%NumIntElems
    READ(14,*) j
    READ(14,*) k
    if (k.eq.FSI%NumIntElems) then ! use default 3 
     k=3
    endif
    if (k.eq.4) then
     numquads=numquads+1
    else if (k.ne.3) then
     print *,"elements must be triangles or quadrilaterals k=",k
     print *,"i,NumIntElems ",i,FSI%NumIntElems
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

   READ(14,91) FSI%NumNodes
   print *,"NumNodes ",FSI%NumNodes
   READ(14,91) FSI%NumIntElems
   print *,"NumIntElems ",FSI%NumIntElems
   READ(14,91) FSI%IntElemDim
   print *,"IntElemDim ",FSI%IntElemDim
   print *,"overriding IntElemDim to be 3"
   FSI%IntElemDim=3
 
   FSI%NumIntElems=FSI%NumIntElems+numquads

   if (ifirst.eq.1) then
    call init_FSI(FSI,1)
   endif

   HdrFlag(1)=sdim
   HdrFlag(2)=0     ! stationary body

   do dir=1,3
    maxnode(dir)=0.0
    minnode(dir)=0.0
    maxnodebefore(dir)=-1.0e+10
    minnodebefore(dir)=1.0e+10
   enddo

   shift_from_zero_node=0
   shift_from_zero_face=0
   override_IntElemDim=0
   do i=1,FSI%NumNodes
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
     FSI%Node_old(dir,i)=xval(dir)
    enddo
    do dir=1,3
     FSI%Node_new(dir,i)=xval(dir)
    enddo
    if (i.ne.j) then
     print *,"vertex mismatch"
     print *,"NumNodesPool=",FSI%NumNodesPool
     print *,"expected node index ",i-FSI%NumNodesPool
     print *,"node index read from file ",j
     stop
    endif
   enddo !  reading nodes
      
   READ(14,*) j
   if (j.ne.FSI%NumIntElems-numquads) then
    print *,"face mismatch"
    stop
   endif

   quad_counter=0
   do i=1,FSI%NumIntElems-numquads
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
     print *,"NumIntElemsPool=",FSI%NumIntElemsPool
     print *,"expected face index ",i-FSI%NumIntElemsPool
     print *,"face index read from file ",j
     stop
    endif
    READ(14,*) k
    if (k.eq.FSI%NumIntElems) then  ! use default 3
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
     FSI%ElemData(1,i+quad_counter)=3   ! number of nodes in element
     FSI%ElemData(2,i+quad_counter)=1   ! part number
     FSI%ElemData(3,i+quad_counter)=0   ! singly wetted (set =1 doubly)
     if (j1.eq.0) then
      do dir=1,SDIM
       FSI%IntElem(dir,i+quad_counter)=localElem(dir)
      enddo
     else if (j1.eq.1) then
      FSI%IntElem(1,i+quad_counter)=localElem(3)
      FSI%IntElem(2,i+quad_counter)=localElem(4)
      FSI%IntElem(3,i+quad_counter)=localElem(1)
     else
      print *,"j1 invalid"
      stop
     endif
     do m1=1,3
      if (shift_from_zero_node.eq.1) then
       FSI%IntElem(m1,i+quad_counter)=FSI%IntElem(m1,i+quad_counter)+1
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

   call init2_FSI(FSI)
  endif ! iread=1

  use_temp=0 
  call init3_FSI(FSI,ifirst,1,1,1)  ! do_2nd_part=1 isout=1

  print *,"after viorel_sphere geominit  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt,Theta_Dot,Theta ",timeB,dtB,TorqueVel,TorquePos
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
      xxblob(1),xxblob(3),newxxblob(1),newxxblob(3),radradblob
  print *," xnew=cos(theta)*(x-x_0)+sin(theta)*(z-z_0)+x_0"
  print *," znew=-sin(theta)*(x-x_0)+cos(theta)*(z-z_0)+z_0"

91     FORMAT(I7)
92     FORMAT(f11.3)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine internal_inflow_geominit


! overall_solid_advance and generate_new_triangles called every geom_interval
! time.
! This routine will be called just once at the very beginning or upon
! restart.
subroutine gearinit(curtime,dt,ifirst,sdim,istop,istep)
INTEGER_T :: i,j,k,i1,j1,k1,m1,tempelem,itimecount,ifirst,sdim
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
INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
INTEGER_T :: gwave
INTEGER_T :: shift_from_zero_node
INTEGER_T :: shift_from_zero_face
INTEGER_T :: override_IntElemDim
character(100) :: dwave,poolname
REAL_T :: ofs

#include "probdataf95.H"

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=dt

  allocate(HdrFlag(NumFlags(1)))
  allocate(NptFlag(1))
  allocate(ForFlag(NumFlags(4))) 

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

  READ(14,*) FSI%NumNodes,FSI%NumIntElems
  FSI%IntElemDim=3

  call init_FSI(FSI,1)

  HdrFlag(1)=sdim
  HdrFlag(2)=0     ! stationary body

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
   maxnodebefore(dir)=-1.0e+10
   minnodebefore(dir)=1.0e+10
  enddo

  do i=1,FSI%NumNodes
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
    FSI%Node_old(dir,i)=xval(dir)
   enddo
   do dir=1,3
    FSI%Node_new(dir,i)=xval(dir)
   enddo
      
  enddo

  do i=1,FSI%NumIntElems
   READ(14,*) FSI%IntElem(3,i),FSI%IntElem(2,i),FSI%IntElem(1,i)

   FSI%ElemData(1,i)=3   ! number of nodes in element
   FSI%ElemData(2,i)=1   ! part number
   FSI%ElemData(3,i)=0   ! singly wetted
  enddo  ! i, looping faces
  close(14)

  do dir=1,3 
   print *,"(before)dir,min,max ",dir,minnodebefore(dir),maxnodebefore(dir)
  enddo
  do dir=1,3 
   print *,"(after)dir,min,max ",dir,minnode(dir),maxnode(dir)
  enddo

  call init2_FSI(FSI)

  use_temp=0 
  call init3_FSI(FSI,1,0,1,1)  ! do_2nd_part=0  isout=1 (gearinit)
  do i=1,FSI%NumNodes
   do dir=1,3
    FSI%Node(dir,i)=FSI%Node_new(dir,i)
   enddo
  enddo

91     FORMAT(I7)
92     FORMAT(f11.3)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine gearinit




! this subroutine generates a speed vs. time profile for 
! use with the first type of inflow 
subroutine timefluct(cur_time,value)
IMPLICIT NONE

  REAL_T, intent(in) :: cur_time
  REAL_T, intent(out) :: value
  REAL_T inittime, medtime, endtime, cur_timeW

#include "probdataf95.H"
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
INTEGER_T :: i,j,k,i1,j1,k1,itimecount,ifirst,sdim
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xtemp,vtemp
REAL_T, dimension(3) :: maxnodebefore,minnodebefore
REAL_T, dimension(3) :: xvalbefore
REAL_T, dimension(3) :: invertfactor
REAL_T :: tper,tcrit,theta,radradblobpool
INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
INTEGER_T :: shift_from_zero_node
INTEGER_T :: shift_from_zero_face
INTEGER_T :: override_IntElemDim,whale_type
REAL_T :: ofs,CLSVOF_curtime,CLSVOF_dt
REAL_T :: STEPSPERIOD,LL_CLSVOF,UU_CLSVOF,TT_CLSVOF
REAL_T :: LL_DUFFY,TT_DUFFY,UU_DUFFY

#include "probdataf95.H"

  if (probtype.ne.562) then
   print *,"probtype must be 562 for animated whale"
   stop
  endif

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=dt

  if (ifirst.eq.1) then
   allocate(HdrFlag(NumFlags(1)))
   allocate(NptFlag(1))
   allocate(ForFlag(NumFlags(4))) 
  endif

  FSI%NumNodesPool=0
  FSI%NumIntElemsPool=0
  FSI%IntElemDimPool=0

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

  FSI%NumNodes=whale_nodes 
  FSI%NumIntElems=whale_cells
  FSI%IntElemDim=3

  if (ifirst.eq.1) then
   call init_FSI(FSI,0)
  endif ! ifirst.eq.1

  HdrFlag(1)=sdim
  HdrFlag(2)=0     ! stationary body

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
   maxnodebefore(dir)=-1.0e+10
   minnodebefore(dir)=1.0e+10
  enddo

  shift_from_zero_node=0
  shift_from_zero_face=0
  override_IntElemDim=0

  do i=1,FSI%NumNodes
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
    FSI%Node_old(dir,i)=xval(dir)
    FSI%Node_new(dir,i)=xval(dir)
   enddo
  enddo ! i=1,NumNodes

  do i=1,FSI%NumIntElems
   FSI%ElemData(1,i)=3   ! number of nodes in element
   FSI%ElemData(2,i)=1   ! part number
   FSI%ElemData(3,i)=0   ! singly wetted
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

  do i=1,FSI%NumNodes
   vtemp(1)=UU(i)
   vtemp(2)=VV(i)
   vtemp(3)=WW(i)
   do dir=1,3
    vtemp(dir)=vtemp(dir)*invertfactor(dir)*UU_DUFFY/UU_CLSVOF
   enddo
   FSI%NodeVel_old(1,i)=vtemp(3)
   FSI%NodeVel_old(2,i)=vtemp(1)
   FSI%NodeVel_old(3,i)=vtemp(2)
   do dir=1,3
    FSI%Node_current(dir,i)=FSI%Node_new(dir,i)
    FSI%NodeVel_new(dir,i)=FSI%NodeVel_old(dir,i)
   enddo
  enddo  ! i=1,NumNodes

  use_temp=0 
  call init3_FSI(FSI,ifirst,0,1,1)  ! do_2nd_part=0  isout=1 (whalegeominit)

  print *,"CLSVOF_curtime,CLSVOF_dt ",CLSVOF_curtime,CLSVOF_dt
  print *,"after whale_geominit  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt ",timeB,dtB
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
      xxblob(1),xxblob(3),newxxblob(1),newxxblob(3),radradblob

91     FORMAT(I7)
92     FORMAT(f11.3)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine whale_geominit



subroutine initpaddle(curtime,dt,sdim,istop,istep, &
  paddle_pos,paddle_vel)
INTEGER_T :: i,j,k,i1,j1,k1,itimecount,sdim
INTEGER_T :: inode,iface
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xtemp
REAL_T :: tper,tcrit,theta
INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
character(20) :: dwave,dwave2
REAL_T :: ofs,xx,zz
REAL_T, intent(in) :: paddle_pos,paddle_vel

#include "probdataf95.H"

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=0.0
  TorquePos=paddle_pos
  TorqueVel=paddle_vel

  allocate(HdrFlag(NumFlags(1)))
  allocate(NptFlag(1))
  allocate(ForFlag(NumFlags(4))) 

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

  READ(14,91) FSI%NumNodesPaddle
  print *,"NumNodesPaddle ",FSI%NumNodesPaddle
  READ(14,91) FSI%NumIntElemsPaddle
  print *,"NumIntElemsPaddle ",FSI%NumIntElemsPaddle
  READ(14,91) FSI%IntElemDimPaddle
  print *,"IntElemDimPaddle ",FSI%IntElemDimPaddle

  READ(15,91) FSI%NumNodes
  print *,"NumNodes ",FSI%NumNodes
  READ(15,91) FSI%NumIntElems
  print *,"NumIntElems ",FSI%NumIntElems
  READ(15,91) FSI%IntElemDim
  print *,"IntElemDim ",FSI%IntElemDim

  FSI%NumNodes=FSI%NumNodes+FSI%NumNodesPaddle
  FSI%NumIntElems=FSI%NumIntElems+FSI%NumIntElemsPaddle
  if (FSI%IntElemDim.lt.FSI%IntElemDimPaddle) then
   FSI%IntElemDim=FSI%IntElemDimPaddle
  endif

  call init_FSI(FSI,1)

  HdrFlag(1)=sdim
  HdrFlag(2)=0     ! stationary body

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
  enddo

  do inode=1,FSI%NumNodes
   if (inode.le.FSI%NumNodesPaddle) then
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
    if (inode.ne.j+FSI%NumNodesPaddle) then
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
    FSI%Node_old(dir,inode)=xval(dir)
    FSI%Node_new(dir,inode)=xval(dir)
   enddo
      
  enddo  ! inode=1,NumNodes
 
  READ(14,91) j
  if (j.ne.FSI%NumIntElemsPaddle) then
   print *,"face mismatch"
   stop
  endif

  READ(15,91) j
  if (j.ne.FSI%NumIntElems-FSI%NumIntElemsPaddle) then
   print *,"face mismatch"
   stop
  endif

  do iface=1,FSI%NumIntElems
   if (iface.le.FSI%NumIntElemsPaddle) then
    READ(14,91) j
    if (iface.ne.j) then
     print *,"face mismatch"
     stop
    endif
    READ(14,91) k
    if (k.gt.FSI%IntElemDim) then
     print *,"14 too many vertices k,IntElemDim ",k,FSI%IntElemDim
     stop
    endif
    do j1=1,k
     READ(14,91) FSI%IntElem(j1,iface) 
    enddo
   else
    READ(15,91) j
    if (iface.ne.j+FSI%NumIntElemsPaddle) then
     print *,"face mismatch"
     stop
    endif
    READ(15,91) k
    if (k.gt.FSI%IntElemDim) then
     print *,"15 too many vertices k,IntElemDim ",k,FSI%IntElemDim
     stop
    endif
    do j1=1,k
     READ(15,91) FSI%IntElem(j1,iface) 
     FSI%IntElem(j1,iface)=FSI%IntElem(j1,iface)+FSI%NumNodesPaddle
    enddo
   endif

   FSI%ElemData(1,iface)=k   ! number of nodes in element
   FSI%ElemData(2,iface)=1   ! part number
   FSI%ElemData(3,iface)=0   ! singly wetted
   FSI%ElemData(3,iface)=1  ! doubly wetted
   if (iface.gt.FSI%NumIntElemsPaddle) then
    FSI%ElemData(3,iface)=1  ! doubly wetted
   endif
  enddo  ! iface, looping faces

  close(14)
  close(15)
 
  print *,"minnode,maxnode ",minnode(1),minnode(2),minnode(3), &
      maxnode(1),maxnode(2),maxnode(3)

  call init2_FSI(FSI)

  do inode=1,FSI%NumNodesPaddle
   xx=FSI%Node_current(1,inode)-newxxblob(1)
   zz=FSI%Node_current(3,inode)-newxxblob(3)
   FSI%Node_new(1,inode)=xx*cos(TorquePos)+zz*sin(TorquePos)+newxxblob(1)
   FSI%Node_new(3,inode)=-xx*sin(TorquePos)+zz*cos(TorquePos)+newxxblob(3)
   xx=FSI%Node_new(1,inode)-newxxblob(1)
   zz=FSI%Node_new(3,inode)-newxxblob(3)
   FSI%NodeVel_new(1,inode)=TorqueVel*zz
   FSI%NodeVel_new(3,inode)=-TorqueVel*xx
  enddo

  use_temp=0

  call init3_FSI(FSI,1,0,1,1)  ! do_2nd_part=0 isout=1 (initpaddle)
  do i=1,FSI%NumNodes
   do dir=1,3
    FSI%Node(dir,i)=FSI%Node_new(dir,i)
   enddo
  enddo

91     FORMAT(I7)
92     FORMAT(f11.3)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine initpaddle


subroutine initship(curtime,dt,sdim,istop,istep, &
  paddle_pos,paddle_vel)
INTEGER_T :: i,j,k,i1,j1,k1,itimecount,sdim
INTEGER_T :: inode,iface
REAL_T :: curtime,dt
REAL_T, dimension(3) :: maxnode,minnode,xval,xtemp
REAL_T :: tper,tcrit,theta
INTEGER_T :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
INTEGER_T :: filler
character(40) :: dwave,dwave2
REAL_T :: ofs,xx,zz
REAL_T, intent(in) :: paddle_pos,paddle_vel

#include "probdataf95.H"

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=0.0
  TorquePos=paddle_pos
  TorqueVel=paddle_vel

  allocate(HdrFlag(NumFlags(1)))
  allocate(NptFlag(1))
  allocate(ForFlag(NumFlags(4))) 

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

  READ(14,191) FSI%NumNodes,FSI%NumIntElems
  FSI%IntElemDim=3

  call init_FSI(FSI,1)

  HdrFlag(1)=sdim
  HdrFlag(2)=0     ! stationary body

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
  enddo

  do inode=1,FSI%NumNodes
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
    FSI%Node_old(dir,inode)=xval(dir)
    FSI%Node_new(dir,inode)=xval(dir)
   enddo
      
  enddo  ! inode=1,NumNodes
 
  do iface=1,FSI%NumIntElems
   READ(14,194) FSI%IntElem(3,iface),FSI%IntElem(2,iface), &
     FSI%IntElem(1,iface),filler

   FSI%ElemData(1,iface)=3   ! number of nodes in element
   FSI%ElemData(2,iface)=1   ! part number
   if ((axis_dir.eq.1).or.(axis_dir.eq.3)) then
    FSI%ElemData(3,iface)=0   ! singly wetted
   else if (axis_dir.eq.2) then
    FSI%ElemData(3,iface)=1   ! doubly wetted
   else
    print *,"bad axis_dir initship2 probtype,axis_dir ",probtype,axis_dir
    stop
   endif
!   FSI%ElemData(3,iface)=1   ! doubly wetted

  enddo  ! iface, looping faces

  close(14)
 
  print *,"minnode,maxnode ",minnode(1),minnode(2),minnode(3), &
      maxnode(1),maxnode(2),maxnode(3)

  call init2_FSI(FSI)

  use_temp=0
  call init3_FSI(FSI,1,0,1,1)  ! do_2nd_part=0 isout=1 (initship)
  do i=1,FSI%NumNodes
   do dir=1,3
    FSI%Node(dir,i)=FSI%Node_new(dir,i)
   enddo
  enddo

191    FORMAT(I12,I12)
91     FORMAT(I7)
92     FORMAT(f11.3)
193    FORMAT(E15.11,E15.11,E15.11)
194    FORMAT(I12,I12,I12,I12)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)

return
end subroutine initship

subroutine overall_solid_advance(CLSVOF_curtime,CLSVOF_dt, &
  FSI_advance,ioproc,isout)
    
type(mesh_type) :: FSI_advance
INTEGER_T :: ifirst,ioproc,isout
REAL_T :: CLSVOF_curtime,CLSVOF_dt
REAL_T :: STEPSPERIOD,LL_CLSVOF,UU_CLSVOF,TT_CLSVOF,whale_dt

#include "probdataf95.H"

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
  call whale_geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,CLSVOF_whale_time,CLSVOF_dt)
  ! animated heart
 else if (probtype.eq.57) then
  call geominit(CLSVOF_curtime,CLSVOF_dt,ifirst,sci_sdim,sci_istop,sci_istep)
 else if ((probtype.eq.52).or. &
     ((probtype.eq.56).and.  &
      ((axis_dir.eq.0).or.(axis_dir.eq.2).or.(axis_dir.eq.3))).or.  &
     (probtype.eq.561).or. &
     ((probtype.eq.55).and.(raddust.eq.0.0)).or. &
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
 else if ((probtype.eq.53).and.(MARCO.eq.0)) then
  call initinjector(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,FSI_advance,isout)
 else if ((probtype.eq.53).and.(axis_dir.gt.0).and.(MARCO.eq.1)) then
  call initinjector(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,FSI_advance,isout)
 else if ((probtype.eq.531).or. &
          (probtype.eq.5501)) then
  call initinjector(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,FSI_advance,isout)

   ! ifirst=0 (in overall_solid_advance)
 else if ((probtype.eq.536).or.(probtype.eq.537).or. &
          (probtype.eq.538).or.(probtype.eq.539)) then
  call initinjector(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,FSI_advance,isout)    
   ! ifirst=0 (in overall_solid_advance)
 else if (probtype.eq.701) then
  call initflapping(CLSVOF_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,ioproc,FSI_advance,isout)    
 else
  call advance_solid(sci_sdim,sci_curtime,sci_dt,sci_istop,sci_istep)
 endif

return
end subroutine overall_solid_advance

subroutine overall_solid_init(CLSVOFtime,ioproc,FSI_init,isout)
type(mesh_type) :: FSI_init
INTEGER_T :: ifirst,ioproc,isout
REAL_T :: paddle_pos,paddle_vel,CLSVOFtime,CLSVOF_dt
REAL_T :: STEPSPERIOD,LL_CLSVOF,UU_CLSVOF,TT_CLSVOF,whale_dt

#include "probdataf95.H"

 if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
 endif

 normal_invert=0
 sci_sdim=3    
 exclusive_doubly_wetted=0

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
 ifirst=1  ! which means iread=1 in initinjector

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
 else if (probtype.eq.5700) then
  call initchannel(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep)    
 else if ((probtype.eq.53).or.(probtype.eq.531).or. &
          (probtype.eq.5501)) then
  call initinjector(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep, &
    ioproc,FSI_init,isout) 

   ! ifirst=1 (overall_solid_init)
 else if ((probtype.eq.536).or.(probtype.eq.537).or. &
          (probtype.eq.538).or.(probtype.eq.539)) then
  call initinjector(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep, &
    ioproc,FSI_init,isout) 
   ! ifirst=1 (overall_solid_init)
 else if (probtype.eq.701) then
  call initflapping(CLSVOFtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep, &
    ioproc,FSI_init,isout) 
 else if (probtype.eq.50) then
  call initpaddle(sci_curtime,sci_dt,sci_sdim,sci_istop,sci_istep, &
    paddle_pos,paddle_vel)
 else if (probtype.eq.9) then
  call initship(sci_curtime,sci_dt,sci_sdim,sci_istop,sci_istep, &
    paddle_pos,paddle_vel)
 else
  print *,"do not know how to initialize the solid"
  stop
 endif
 if ((ioproc.eq.1).and.(isout.eq.1)) then
  print *,"after initialize solid dt=",sci_dt
 endif

return
end subroutine overall_solid_init


subroutine advance_solid(sdim,curtime,dt,istop,istep)
INTEGER_T :: sdim
REAL_T :: curtime,dt,totaltorque,xx,zz
INTEGER_T :: istop
INTEGER_T :: i,j,k,dir,istep,istepcfd

#include "probdataf95.H"

 DebugCouple=.true.
 
 do i=1,FSI%NumNodes
  do dir=1,3
   FSI%Node_old(dir,i)=FSI%Node_new(dir,i)
   FSI%NodeVel_old(dir,i)=FSI%NodeVel_new(dir,i)
  enddo
 enddo

 do i=1,FSI%NumNodes
  do dir=1,3
   FSI%Node(dir,i)=FSI%Node_new(dir,i)
   FSI%NodeVel(dir,i)=FSI%NodeVel_new(dir,i)
  enddo
  FSI%NodeTemp(i)=FSI%NodeTemp_new(i)
 enddo

 dtB=0.0
 if (curtime.gt.timeB) then
  dtB=curtime-timeB
  timeB=curtime
 endif

  ! this first case will never happen since initinjector is called 
  ! instead of this routine. 
 if (probtype.eq.538) then  

  do i=1,FSI%NumNodes
   do dir=1,3
    FSI%Node_new(dir,i)=FSI%Node_new(dir,i)+dtB*FSI%NodeVel_new(dir,i)
   enddo
  enddo

  ! this second case will never happen since initflapping is called 
  ! instead of this routine. 
 else if (probtype.eq.701) then

  do i=1,FSI%NumNodes
   do dir=1,3
    FSI%Node_new(dir,i)=FSI%Node_new(dir,i)+dtB*FSI%NodeVel_new(dir,i)
   enddo
  enddo

 else

   ! default for ship, gear, ...
  do i=1,FSI%NumNodes
   do dir=1,3
    FSI%NodeVel_new(dir,i)=0.0
   enddo
  enddo

 endif

return
end subroutine advance_solid





! --------------------  SOLID ADVANCE STUFF ENDS HERE --------------

subroutine checkinbox(xc,boxlo,boxhi,istat,eps)
IMPLICIT NONE

REAL_T, dimension(3), intent(in) :: xc(3),boxlo(3),boxhi(3)
REAL_T :: eps
INTEGER_T, intent(out) :: istat
INTEGER_T :: dir

 istat=1
 do dir=1,3
  if ((xc(dir).lt.boxlo(dir)-eps).or.(xc(dir).gt.boxhi(dir)+eps)) then
   istat=0
  endif
 enddo

return
end subroutine checkinbox



subroutine checkinpoint(xclosest,inode,elemnum, &
  mindist,xc,inplane,FSI_in)
IMPLICIT NONE

type(mesh_type) :: FSI_in
INTEGER_T, intent(in) :: inode,elemnum
INTEGER_T, intent(inout) :: inplane
REAL_T, dimension(3), intent(inout) :: xclosest
REAL_T, dimension(3), intent(in) :: xc
REAL_T, intent(inout) :: mindist
REAL_T, dimension(3) :: xnot,normal
REAL_T :: curdist,dotprod
INTEGER_T :: dir,i,j,k
INTEGER_T :: nodes_per_elem


 nodes_per_elem=FSI_in%ElemDataBIG(1,elemnum)

 if (nodes_per_elem.ne.3) then
  print *,"nodes_per_elem invalid"
  stop
 endif

 do dir=1,3
  xnot(dir)=FSI_in%NodeBIG(dir,FSI_in%IntElemBIG(inode,elemnum))
  normal(dir)=FSI_in%NodeNormalBIG(dir,FSI_in%IntElemBIG(inode,elemnum))
 enddo

 dotprod=0.0
 do dir=1,3
  dotprod=dotprod+normal(dir)*(xc(dir)-xnot(dir))
 enddo

 call xdist(xnot,xc,curdist)
 if (dotprod.lt.0.0) then
  curdist=-curdist
 endif

 if ((abs(curdist).lt.abs(mindist)).or.(inplane.eq.0)) then
  inplane=1
  mindist=curdist
  do dir=1,3
   xclosest(dir)=xnot(dir)
  enddo
 endif

return
end subroutine checkinpoint



subroutine checkinline(xclosest,inode,elemnum, &
  tol,mindist,xc,inplane,FSI_in)
IMPLICIT NONE

type(mesh_type) :: FSI_in
INTEGER_T, intent(in) :: inode,elemnum
INTEGER_T, intent(inout) :: inplane
REAL_T, intent(in) :: tol
REAL_T, dimension(3), intent(inout) :: xclosest
REAL_T, dimension(3), intent(in) :: xc
REAL_T, dimension(2,3) :: xnode,nnode
REAL_T, intent(inout) :: mindist
REAL_T, dimension(3) :: xnot,normal
REAL_T :: dottop,dotbot,t,curdist,dotprod
INTEGER_T :: dir,i,j,k
INTEGER_T :: nodes_per_elem

 if (inode.le.2) then

  nodes_per_elem=FSI_in%ElemDataBIG(1,elemnum)
  if (nodes_per_elem.ne.3) then
   print *,"nodes_per_elem invalid"
   stop
  endif

  do dir=1,3
   xnode(1,dir)=FSI_in%NodeBIG(dir,FSI_in%IntElemBIG(inode,elemnum))
   xnode(2,dir)=FSI_in%NodeBIG(dir,FSI_in%IntElemBIG(inode+1,elemnum))
   nnode(1,dir)=FSI_in%NodeNormalBIG(dir,FSI_in%IntElemBIG(inode,elemnum))
   nnode(2,dir)=FSI_in%NodeNormalBIG(dir,FSI_in%IntElemBIG(inode+1,elemnum))
  enddo

  dottop=0.0
  dotbot=0.0
  do dir=1,3
   dottop=dottop+(xnode(2,dir)-xnode(1,dir))*(xnode(2,dir)-xc(dir))
   dotbot=dotbot+(xnode(2,dir)-xnode(1,dir))**2
  enddo
  dotbot=dotbot+1.0E-14
  t=dottop/dotbot
  if ((t.ge.-tol).and.(t.le.1.0+tol)) then
   do dir=1,3
    xnot(dir)=t*xnode(1,dir)+(1.0-t)*xnode(2,dir)
    normal(dir)=t*nnode(1,dir)+(1.0-t)*nnode(2,dir)
   enddo

   dotprod=0.0
   do dir=1,3
    dotprod=dotprod+normal(dir)*(xc(dir)-xnot(dir))
   enddo

   call xdist(xnot,xc,curdist)
   if (dotprod.lt.0.0) then
    curdist=-curdist
   endif

   if ((abs(curdist).lt.abs(mindist)).or.(inplane.eq.0)) then
    inplane=1
    mindist=curdist
    do dir=1,3
     xclosest(dir)=xnot(dir)
    enddo
   endif
  endif ! 0<=t<=1

 endif

return
end subroutine checkinline


subroutine checkinplane(xclosest,elemnum,inplane, &
  tol,minnode,maxnode,FSI_in)
IMPLICIT NONE

type(mesh_type) :: FSI_in
INTEGER_T, intent(in) :: elemnum
REAL_T, intent(in) :: tol
REAL_T, dimension(3), intent(in) :: minnode,maxnode,xclosest
INTEGER_T, intent(out) :: inplane
REAL_T, dimension(3,3) :: xnode,AA,AI
INTEGER_T :: dir,i,j,k
INTEGER_T :: nodes_per_elem
REAL_T :: det
REAL_T, dimension(3) :: tx
REAL_T, dimension(3) :: v1,v2,v1xv2

 nodes_per_elem=FSI_in%ElemDataBIG(1,elemnum)
 if (nodes_per_elem.ne.3) then
  print *,"nodes_per_elem invalid"
  stop
 endif

 inplane=1
 do dir=1,3
  if ((xclosest(dir).lt.minnode(dir)-tol).or. &
      (xclosest(dir).gt.maxnode(dir)+tol)) then
   inplane=0
  endif
 enddo

 if (inplane.eq.1) then

  do i=1,3
   do dir=1,3
    xnode(i,dir)=FSI_in%NodeBIG(dir,FSI_in%IntElemBIG(i,elemnum))
   enddo
  enddo  ! i

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
  
   if ((tx(1).lt.-tol).or.(tx(1).gt.1.0+tol).or. &
       (tx(2).lt.-tol).or.(tx(2).gt.1.0+tol).or. &
       (tx(1)+tx(2).gt.1.0+tol)) then
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
subroutine scinormal(elemnum,normal,FSI_in)
IMPLICIT NONE

type(mesh_type) :: FSI_in
INTEGER_T, intent(in) :: elemnum
REAL_T, dimension(3), intent(out) :: normal
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3,3) :: nodesave
REAL_T, dimension(3) :: nodeavg
REAL_T, dimension(2,3) :: vec
REAL_T :: dist
INTEGER_T :: i,j
INTEGER_T :: local_normal_invert

#include "probdataf95.H"

  nodes_per_elem=FSI_in%ElemDataBIG(1,elemnum)
  if (nodes_per_elem.gt.3) then
   print *,"nodes_per_elem>3 not supported"
   stop
  endif

  do j=1,3
   nodeavg(j)=0.0
   do i=1,3
    nodesave(i,j)=FSI_in%NodeBIG(j,FSI_in%IntElemBIG(i,elemnum))
    nodeavg(j)=nodeavg(j)+nodesave(i,j)
   enddo
   nodeavg(j)=nodeavg(j)/3.0
  enddo

  do i=1,2
   do j=1,3
    vec(i,j)=nodesave(i+1,j)-nodesave(i,j)
   enddo
  enddo

  normal(1)=vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2)
  normal(2)=vec(1,3)*vec(2,1)-vec(1,1)*vec(2,3)
  normal(3)=vec(1,1)*vec(2,2)-vec(2,1)*vec(1,2)

  dist=sqrt(normal(1)**2+normal(2)**2+normal(3)**2)
  if (dist.gt.1.0e-15) then
   do i=1,3
    normal(i)=normal(i)/dist
   enddo
  endif

  local_normal_invert=normal_invert

  if (local_normal_invert.eq.1) then
   do i=1,3
    normal(i)=-normal(i)
   enddo
  else if (local_normal_invert.ne.0) then
   print *,"local_normal_invert must be 0 or 1"
   stop
  endif

return
end subroutine scinormal


subroutine sciarea(elemnum,area,FSI_in)
IMPLICIT NONE

type(mesh_type) :: FSI_in
INTEGER_T, intent(in) :: elemnum
REAL_T, intent(out) :: area
INTEGER_T :: nodes_per_elem
REAL_T, dimension(3,3) :: nodesave
REAL_T, dimension(2,3) :: vec
INTEGER_T :: i,j
REAL_T :: aa,bb,cc

#include "probdataf95.H"

  nodes_per_elem=FSI_in%ElemDataBIG(1,elemnum)
  if (nodes_per_elem.ne.3) then
   print *,"nodes in element?  elem,nodes_per_elem ",elemnum,nodes_per_elem   
   stop
  endif

  do i=1,3
   do j=1,3
    nodesave(i,j)=FSI_in%NodeBIG(j,FSI_in%IntElemBIG(i,elemnum))
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




subroutine CLSVOF_ReadHeader(h_small,CLSVOFtime,problo,probhi, &
  xproblo,xprobhi,ioproc,isout)
IMPLICIT NONE

INTEGER_T :: i,j,k,initflag,ioproc,isout
REAL_T :: h_small,CLSVOFtime
REAL_T problo(SDIM),probhi(SDIM)
REAL_T xproblo(SDIM),xprobhi(SDIM) ! local processor domain

#include "probdataf95.H"

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=10
  UnitOut=20
  UnitDebug=30
  StdOut=6
  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4
  use_temp=0

  FSI%PartID = 1

   ! ifirst=1 (=>iread=1 in initinjector or initflapping)
  call overall_solid_init(CLSVOFtime,ioproc,FSI,isout)  

    ! ReadHeader
  initflag=1
  call generate_new_triangles(initflag,problo,probhi,FSI, &
   SCI,0,ioproc,isout)

  if(probtype.eq.538) then ! moving parts
    FSI_needle%PartID = 2
    call overall_solid_init(CLSVOFtime,ioproc,FSI_needle, &
      isout)  ! ifirst=1

    call generate_new_triangles(initflag,problo,probhi,FSI_needle, &
      SCI_needle,0,ioproc,isout)

    ! repeat for a subset of the domain
    call generate_new_triangles(initflag,xproblo,xprobhi,FSI,&
     SCI_loc,1,ioproc,isout)

    call generate_new_triangles(initflag,xproblo,xprobhi,FSI_needle,&
     SCI_needle_loc,1,ioproc,isout)

  else if (probtype.eq.701) then  ! flapping wing - ReadHeader

   if (axis_dir.eq.0) then
    ! do nothing - just 1 wing
   else if (axis_dir.eq.1) then
    FSI_needle%PartID = 2
    call overall_solid_init(CLSVOFtime,ioproc,FSI_needle, &
      isout)  ! ifirst=1
    call generate_new_triangles(initflag,problo,probhi,FSI_needle, &
      SCI_needle,0,ioproc,isout)
   else
    print *,"axis_dir invalid"
    stop
   endif 

  endif

 
return
end subroutine CLSVOF_ReadHeader

subroutine GetDomainScale(maxnodedomain,minnodedomain,maxside, &
  maxnoderadius,COM,FSI_in)
IMPLICIT NONE

type(mesh_type) :: FSI_in
REAL_T, dimension(3), intent(out) :: maxnodedomain,minnodedomain
REAL_T, dimension(3), intent(out) :: COM
REAL_T, intent(out) :: maxnoderadius
REAL_T, intent(out) :: maxside
REAL_T :: tempside
INTEGER_T :: dir,ifirst,ielem,nodes_per_elem,i,inode
INTEGER_T :: nodes_counted
REAL_T, dimension(3) :: xtarget(3)
REAL_T :: noderadius

 do dir=1,3
  maxnodedomain(dir)=0.0
  minnodedomain(dir)=0.0
  COM(dir)=0.0
 enddo
 nodes_counted=0
 maxnoderadius=0.0

 ifirst=1
 do ielem=1,FSI_in%NumIntElems
  if (FSI_in%ActiveIntElem(ielem)) then
   nodes_per_elem=FSI_in%ElemData(1,ielem)

   do i=1,nodes_per_elem
    inode=FSI_in%IntElem(i,ielem) 

    do dir=1,3
     if ((minnodedomain(dir).gt.FSI_in%Node(dir,inode)).or.(ifirst.eq.1)) then
      minnodedomain(dir)=FSI_in%Node(dir,inode)
     endif
     if ((maxnodedomain(dir).lt.FSI_in%Node(dir,inode)).or.(ifirst.eq.1)) then
      maxnodedomain(dir)=FSI_in%Node(dir,inode)
     endif
     COM(dir)=COM(dir)+FSI_in%Node(dir,inode)
    enddo
    nodes_counted=nodes_counted+1

    ifirst=0
   enddo  ! i
  endif ! ActiveIntElem
 enddo ! ielem

 maxside=0.0
 do dir=1,3
  tempside=maxnodedomain(dir)-minnodedomain(dir)
  if (maxside.lt.tempside) then
   maxside=tempside
  endif
  if (nodes_counted.gt.0) then
   COM(dir)=COM(dir)/nodes_counted
  endif
 enddo

 do ielem=1,FSI_in%NumIntElems
  if (FSI_in%ActiveIntElem(ielem)) then
   nodes_per_elem=FSI_in%ElemData(1,ielem)

   do i=1,nodes_per_elem
    inode=FSI_in%IntElem(i,ielem)

    noderadius=zero
    do dir=1,3
     xtarget(dir)=FSI_in%Node(dir,inode)
     noderadius=noderadius+(xtarget(dir)-COM(dir))**2
    enddo
    noderadius=sqrt(noderadius)
    if (noderadius.gt.maxnoderadius) then
     maxnoderadius=noderadius
    endif
   enddo  ! i
  endif ! ActiveIntElem
 enddo ! ielem


return
end subroutine GetDomainScale

! mask=0.0 prior to entry
! mask=1 velocity is init, but sign is not
! mask=2 both sign and velocity are init
! mask=3 both sign and velocity are init (doubly wetted)
! mask=0 neither velocity or sign are valid.  LS has a large value.

! vel=temp=0 if no interfaces in cell's neighborhood
subroutine CLSVOF_InitBox(  &
  mask,DIMS3D(mask), &
  cellproblo,celldx, &
  ls,DIMS3D(ls), &
  lsface,DIMS3D(lsface), &
  vel,DIMS3D(vel), &
  temp,DIMS3D(temp), &
  lo,hi,FSI_in,ioproc)

IMPLICIT NONE

  type(mesh_type) :: FSI_in
  INTEGER_T xpos_ngrow,ioproc
  INTEGER_T lo(SDIM),hi(SDIM)
  REAL_T cellproblo(SDIM)
  REAL_T celldx(SDIM)
  INTEGER_T DIMDEC3D(mask)
  INTEGER_T DIMDEC3D(ls)
  INTEGER_T DIMDEC3D(lsface)
  INTEGER_T DIMDEC3D(vel)
  INTEGER_T DIMDEC3D(temp)

  REAL_T mask(DIMV3D(mask))
  REAL_T ls(DIMV3D(ls))
  REAL_T lsface(DIMV3D(lsface),SDIM)
  REAL_T vel(DIMV3D(vel),SDIM)
  REAL_T temp(DIMV3D(temp))

  INTEGER_T :: ielem,nodes_per_elem,inode,nodeptr
  INTEGER_T :: jelem,elem1,elem2,node1,node2,nd,ngrowc,nn
  REAL_T, dimension(3) :: xc,xproj,velnode,xclosest
  REAL_T, dimension(3) :: boxlo,boxhi,normal,xnot
  REAL_T :: tempnode
  REAL_T, dimension(3) :: minnodedomain,maxnodedomain
  REAL_T, dimension(3) :: COM
  REAL_T :: maxnoderadius
  INTEGER_T, dimension(3) :: gridlo,gridhi
  INTEGER_T :: i,j,k,ifirst
  INTEGER_T :: inbox
  INTEGER_T :: ic,jc,kc
  INTEGER_T :: inplane,bandwidth
  REAL_T :: cutoff,tol,wallthick,maxside,tempside
  REAL_T, dimension(3) :: velcell
  REAL_T :: distcell,maskcell,tempcell
  REAL_T :: del,delmax
  INTEGER_T :: istat
  INTEGER_T :: depth,dir,dir1,dir2,side
  INTEGER_T :: ecount
  INTEGER_T :: local_element_count,iband
  REAL_T :: diag_depth,nodetest,dotprod,mindist,crossing_mindist
  REAL_T :: weighttotal,distwt,weight,eps,out_dist
  REAL_T, dimension(3) :: minnode,maxnode,xx,xside,xcrit
  INTEGER_T :: ii,jj,kk,hitflag
  REAL_T :: phiside,phicenter,testdist,hitsign,totaldist

  REAL_T, dimension(:,:,:), allocatable :: crossing_ls
  INTEGER_T modify_vel
  
#include "probdataf95.H"

  if ((ioproc.ne.1).and.(ioproc.ne.0)) then
    print *,"ioproc invalid"
    stop
  endif
  call checkbound3D(lo,hi,DIMS3D(mask),1,-1,521)
  call checkbound3D(lo,hi,DIMS3D(ls),1,-1,522)
  call checkbound3D(lo,hi,DIMS3D(lsface),1,-1,522)
  call checkbound3D(lo,hi,DIMS3D(vel),1,-1,523)
  call checkbound3D(lo,hi,DIMS3D(temp),1,-1,524)

  iband=3

  call GetDomainScale(maxnodedomain,minnodedomain,maxside, &
    maxnoderadius,COM,FSI_in)
  tol=maxside/10000.0
  if (ioproc.eq.1) then
   print *,"InitBox tol,maxside ",tol,maxside
  endif

  allocate(crossing_ls(DIMV3D(ls)))

  do i=lo(1)-1,hi(1)+1
  do j=lo(2)-1,hi(2)+1
  do k=lo(3)-1,hi(3)+1
   mask(i,j,k)=0.0
   temp(i,j,k)=0.0
   do dir=1,3
    vel(i,j,k,dir)=0.0
    lsface(i,j,k,dir)=1.0E+10
   enddo
   ls(i,j,k)=1.0E+10
   crossing_ls(i,j,k)=1.0E+10
  enddo   
  enddo   
  enddo   

  do ielem=1,FSI_in%NumIntElemsBIG
   do dir=1,3
    xnot(dir)=FSI_in%ElemDataXnotBIG(dir,ielem)
   enddo 

   nodes_per_elem=FSI_in%ElemDataBIG(1,ielem)
   if (nodes_per_elem.lt.3) then
    print *,"elem,nodes_per_elem ",ielem,nodes_per_elem   
    stop
   endif
    ! points from solid to fluid
    ! phi=n dot (x-xnot)
   call scinormal(ielem,normal,FSI_in)

   do inode=1,3
    do dir=1,3
     nodetest=FSI_in%NodeBIG(dir,FSI_in%IntElemBIG(inode,ielem))
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
    enddo
   enddo  ! inode=1,3
   
   call find_grid_bounding_box(minnode,maxnode,cellproblo,celldx, &
      lo,hi,gridlo,gridhi,iband) 

   do i=gridlo(1),gridhi(1)
   do j=gridlo(2),gridhi(2)
   do k=gridlo(3),gridhi(3)

    if ((i.ge.lo(1)-1).and.(i.le.hi(1)+1).and.(j.ge.lo(2)-1).and. &
        (j.le.hi(2)+1).and.(k.ge.lo(3)-1).and.(k.le.hi(3)+1)) then
     xx(1)=cellproblo(1)+(i-lo(1)+half)*celldx(1)
     xx(2)=cellproblo(2)+(j-lo(2)+half)*celldx(2)
     xx(3)=cellproblo(3)+(k-lo(3)+half)*celldx(3)

     dotprod=0.0 
     do dir=1,3
      dotprod=dotprod+normal(dir)*(xx(dir)-xnot(dir))
     enddo
     do dir=1,3
      xclosest(dir)=xx(dir)-dotprod*normal(dir)
     enddo
     mindist=dotprod
     phicenter=mindist

     call checkinplane(xclosest,ielem,inplane,tol,minnode,maxnode,FSI_in)

! investigate using NodeNormalBIG
     do inode=1,nodes_per_elem
      call checkinline(xclosest,inode,ielem,tol,mindist,xx,inplane,FSI_in)
      call checkinpoint(xclosest,inode,ielem,mindist,xx,inplane,FSI_in)
     enddo

     hitflag=0
     hitsign=zero
     crossing_mindist=1.0E+10

! check crossing between cells
! this is where the local sign is initialized.

     do ii=-1,1
     do jj=-1,1
     do kk=-1,1
       if (abs(ii)+abs(jj)+abs(kk).gt.0) then

        do dir2=1,3
         xside(dir2)=xx(dir2)
        enddo
        xside(1)=xside(1)+ii*celldx(1)
        xside(2)=xside(2)+jj*celldx(2)
        xside(3)=xside(3)+kk*celldx(3)

        phiside=zero
        do dir2=1,3
         phiside=phiside+normal(dir2)*(xside(dir2)-xnot(dir2))
        enddo
        if (phiside*phicenter.le.zero) then
         do dir2=1,3
          if (phiside.eq.zero) then
           xcrit(dir2)=xside(dir2)
          else if (phicenter.eq.zero) then
           xcrit(dir2)=xx(dir2)
          else
           xcrit(dir2)=(abs(phiside)*xx(dir2)+ &
                        abs(phicenter)*xside(dir2))/  &
                       (abs(phicenter)+abs(phiside))
          endif
         enddo  ! dir2

         call checkinplane(xcrit,ielem,inplane,tol,minnode,maxnode,FSI_in)
! totaldist is the distance between xx and xside
! testdist  is the distance between xx and xcrit
         if (inplane.eq.1) then
          testdist=zero
          totaldist=zero
          do dir2=1,3
           testdist=testdist+(xcrit(dir2)-xx(dir2))**2
           totaldist=totaldist+(xside(dir2)-xx(dir2))**2
          enddo
          testdist=sqrt(testdist)
          totaldist=sqrt(totaldist)
          if (abs(testdist).le.abs(mindist)) then
           mindist=testdist
          endif
          if (abs(testdist).le.abs(crossing_mindist)) then
           crossing_mindist=testdist
           hitflag=1
           if (phicenter.eq.zero) then
            hitsign=zero
           else if (phicenter.gt.zero) then
            hitsign=one
           else
            hitsign=-one
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

          if (abs(ii)+abs(jj)+abs(kk).eq.1) then
           if (ii.eq.-1) then
            dir=1
            side=1
           else if (ii.eq.1) then
            dir=1
            side=2
           else if (jj.eq.-1) then
            dir=2
            side=1
           else if (jj.eq.1) then
            dir=2
            side=2
           else if (kk.eq.-1) then
            dir=3
            side=1
           else if (kk.eq.1) then
            dir=3
            side=2
           else
            print *,"bust"
            stop
           endif
           if (side.eq.1) then
            lsface(i,j,k,dir)=testdist  ! positive height fraction
           else if (side.eq.2) then
            if ((i+ii.le.hi(1)+1).and. &
                (j+jj.le.hi(2)+1).and. &
                (k+kk.le.hi(3)+1)) then
             lsface(i+ii,j+jj,k+kk,dir)=testdist  ! positive height fraction
            endif
           else
            print *,"side invalid"
            stop
           endif
          endif ! abs(ii)+abs(jj)+abs(kk)=1

         else if (inplane.ne.0) then
          print *,"inplane invalid"
          stop
         endif  ! inplane
        endif ! phiside x phicenter <= 0
       endif ! abs(ii)+abs(jj)+abs(kk)>0
     enddo 
     enddo 
     enddo  ! ii,jj,kk

     modify_vel=0
     
     if ( (abs(mindist).lt.abs(ls(i,j,k))).or. &
          (abs(crossing_mindist).lt.abs(crossing_ls(i,j,k))) ) then

      if (abs(mindist).lt.abs(ls(i,j,k))) then
       modify_vel=1
      endif

      if (FSI_in%ElemDataBIG(3,ielem).eq.1) then ! doubly wetted
       mask(i,j,k)=3.0
       if (abs(mindist).lt.abs(ls(i,j,k))) then
        ls(i,j,k)=abs(mindist)
       endif
       if (abs(crossing_mindist).lt.abs(crossing_ls(i,j,k))) then
        if (hitflag.eq.1) then
         crossing_ls(i,j,k)=abs(crossing_mindist)
        endif
       endif
      else
       if (mask(i,j,k).lt.one) then
        mask(i,j,k)=one
       endif
       if (abs(mindist).lt.abs(ls(i,j,k))) then
        if (mask(i,j,k).ne.two) then
         ls(i,j,k)=abs(mindist)
        else
         if (ls(i,j,k).lt.zero) then
          ls(i,j,k)=-abs(mindist)
         else
          ls(i,j,k)=abs(mindist)
         endif
        endif 
       endif
       if (abs(crossing_mindist).lt.abs(crossing_ls(i,j,k))) then
        if (hitflag.eq.1) then
         mask(i,j,k)=2.0
         ls(i,j,k)=hitsign*abs(ls(i,j,k))
         crossing_ls(i,j,k)=hitsign*abs(crossing_mindist)
        endif
       endif 
      endif

      if (modify_vel.eq.1) then

       do dir2=1,3
        vel(i,j,k,dir2)=0.0d0
       enddo
       temp(i,j,k)=0.0d0

       weighttotal=0.0d0
       do inode=1,nodes_per_elem
        nodeptr=FSI_in%IntElemBIG(inode,ielem)
        distwt=0.0d0
        do dir1=1,3
         distwt=distwt+(FSI_in%NodeBIG(dir1,nodeptr)-xx(dir1))**2
        enddo
        weight=1.0d0/( (distwt+(1.0E-10)**2)**4 )
        do dir2=1,3
         vel(i,j,k,dir2)=vel(i,j,k,dir2)+ &
          weight*FSI_in%NodeVelBIG(dir2,nodeptr)
        enddo
        temp(i,j,k)=temp(i,j,k)+weight*FSI_in%NodeTempBIG(nodeptr)
        weighttotal=weighttotal+weight
       enddo ! inode

       do dir2=1,3
        vel(i,j,k,dir2)=vel(i,j,k,dir2)/weighttotal
       enddo

       if ((probtype.eq.9).and.(axis_dir.gt.1)) then
        if ( xx(1) .gt. -0.4 .and. xx(1) .lt. -0.2 ) then
         if ( xx(3) .lt. 0.045 .and. xx(3) .gt. 0.03 ) then
          vel(i,j,k,3) = -1.0
         endif
        endif
       endif

       temp(i,j,k)=temp(i,j,k)/weighttotal
      endif  ! modify_vel=1
     endif ! mindist.lt.abs(ls) or crossing_mindist.lt.abs(crossing_ls)
    endif ! i,j,k in domain
   enddo
   enddo
   enddo ! i,j,k
  enddo ! ielem

  do i=lo(1)-1,hi(1)+1
  do j=lo(2)-1,hi(2)+1
  do k=lo(3)-1,hi(3)+1
   if (mask(i,j,k).eq.zero) then
    ls(i,j,k)=8.0*celldx(1)
   endif

   if (exclusive_doubly_wetted.eq.0) then

   else if (exclusive_doubly_wetted.eq.1) then
    mask(i,j,k)=3.0
   else
    print *,"exclusive_doubly_wetted invalid"
    stop
   endif
  enddo   
  enddo   
  enddo   

! fix for doubly wetted 

  wallthick=1.0
  
  do i=lo(1)-1,hi(1)+1
  do j=lo(2)-1,hi(2)+1
  do k=lo(3)-1,hi(3)+1
   if (mask(i,j,k).eq.3.0) then
    ls(i,j,k)=ls(i,j,k)-wallthick*celldx(1)
   endif
  enddo
  enddo
  enddo

  deallocate(crossing_ls)

return
end subroutine CLSVOF_InitBox


      subroutine MinModGridInterp3D(data,xpos,x,H)
      IMPLICIT NONE

      REAL_T data(3,3,3)
      REAL_T xpos(3,3,3,SDIM)
      REAL_T x(SDIM)
      REAL_T x0(SDIM)
      REAL_T H
      REAL_T intercept
      REAL_T slopes(SDIM)
      INTEGER_T dir,i1,j1,k1
      REAL_T mindata,maxdata

      do dir=1,SDIM
       x0(dir)=xpos(2,2,2,dir)
      enddo
      call minmod3D(data,xpos,slopes,intercept,SDIM)
      call distfunc3D(intercept,slopes,x0,x,H,SDIM)
      mindata=data(2,2,2)
      maxdata=data(2,2,2)
      do i1=1,3
      do j1=1,3
      do k1=1,3
       if (data(i1,j1,k1).lt.mindata) then
        mindata=data(i1,j1,k1)
       endif
       if (data(i1,j1,k1).gt.maxdata) then
        maxdata=data(i1,j1,k1)
       endif
      enddo 
      enddo 
      enddo 
      if (mindata.gt.H) then
       H=mindata
      endif
      if (maxdata.lt.H) then
       H=maxdata
      endif
 
      return 
      end subroutine MinModGridInterp3D



! clo,chi are dimensions for the secondary grid.
subroutine interp_from_grid(zz,xtarget,scomp,lsflag, &
  data,DIMS3D(data),clo,chi,SCI_in)
IMPLICIT NONE

 type(auxgrid_type) :: SCI_in
 INTEGER_T clo(SDIM)
 INTEGER_T chi(SDIM)
 REAL_T zz
 REAL_T xtarget(SDIM)
 INTEGER_T scomp,lsflag
 INTEGER_T DIMDEC3D(data)
 REAL_T data(DIMV3D(data),scomp)
 INTEGER_T inbox,ic,jc,kc,ii,jj,kk,dir
 REAL_T dataquad(3,3,3)
 REAL_T xposquad(3,3,3,SDIM)
 REAL_T diag,xcell,ycell,zcell

 if ((scomp.lt.1).or.(scomp.gt.SDIM)) then
  print *,"scomp invalid"
  stop
 endif
 if ((lsflag.ne.0).and.(lsflag.ne.1)) then
  print *,"lsflag invalid"
  stop
 endif
 if ((lsflag.eq.1).and.(scomp.ne.1)) then
  print *,"scomp invalid"
  stop
 endif

 call init_dimensions(SCI_in)
 call checkbound3D(clo,chi,DIMS3D(data),1,-1,51)
 inbox=1

  ! NINT->round to nearest whole number
  ! x=(i-lo+half)dx+problo    i=(x-problo)/dx+lo-half
 dir=1
 if (xtarget(dir).lt.SCI_in%Cellproblo(dir)) then
  inbox=0
  ic=clo(dir)-1
 else if (xtarget(dir).gt.SCI_in%Cellprobhi(dir)) then
  inbox=0
  ic=chi(dir)+1
 else
  ic=NINT((xtarget(dir)-SCI_in%Cellproblo(dir))/SCI_in%celldx(dir)- &
    half+clo(dir))
  if (ic.lt.clo(dir)) then
   ic=clo(dir)
  endif
  if (ic.gt.chi(dir)) then
   ic=chi(dir)
  endif
 endif
 dir=2
 if (xtarget(dir).lt.SCI_in%Cellproblo(dir)) then
  inbox=0
  jc=clo(dir)-1
 else if (xtarget(dir).gt.SCI_in%Cellprobhi(dir)) then
  inbox=0
  jc=chi(dir)+1
 else
  jc=NINT((xtarget(dir)-SCI_in%Cellproblo(dir))/SCI_in%celldx(dir)- &
    half+clo(dir)) 
  if (jc.lt.clo(dir)) then
   jc=clo(dir)
  endif
  if (jc.gt.chi(dir)) then
   jc=chi(dir)
  endif
 endif
 dir=3
 if (xtarget(dir).lt.SCI_in%Cellproblo(dir)) then
  inbox=0
  kc=clo(dir)-1
 else if (xtarget(dir).gt.SCI_in%Cellprobhi(dir)) then
  inbox=0
  kc=chi(dir)+1
 else
  kc=NINT((xtarget(dir)-SCI_in%Cellproblo(dir))/SCI_in%celldx(dir)- &
    half+clo(dir))
  if (kc.lt.clo(dir)) then
   kc=clo(dir)
  endif
  if (kc.gt.chi(dir)) then
   kc=chi(dir)
  endif
 endif

 if (inbox.eq.0) then
  zz=data(ic,jc,kc,scomp)
  if (lsflag.eq.1) then
   xcell=SCI_in%Cellproblo(1)+(ic-clo(1)+half)*SCI_in%celldx(1)
   ycell=SCI_in%Cellproblo(2)+(jc-clo(2)+half)*SCI_in%celldx(2)
   zcell=SCI_in%Cellproblo(3)+(kc-clo(3)+half)*SCI_in%celldx(3)

   diag=(xtarget(1)-xcell)**2 + &
        (xtarget(2)-ycell)**2 + &
        (xtarget(3)-zcell)**2
   diag=sqrt(diag)
   if (zz.lt.zero) then
    zz=zz-diag
   else 
    zz=zz+diag
   endif
  endif  
 else
  do ii=1,3
  do jj=1,3
  do kk=1,3
   dataquad(ii,jj,kk)=data(ic+ii-2,jc+jj-2,kc+kk-2,scomp)
   xcell=SCI_in%Cellproblo(1)+(ic+ii-2-clo(1)+half)*SCI_in%celldx(1)
   ycell=SCI_in%Cellproblo(2)+(jc+jj-2-clo(2)+half)*SCI_in%celldx(2)
   zcell=SCI_in%Cellproblo(3)+(kc+kk-2-clo(3)+half)*SCI_in%celldx(3)
   xposquad(ii,jj,kk,1)=xcell
   xposquad(ii,jj,kk,2)=ycell
   xposquad(ii,jj,kk,3)=zcell
  enddo
  enddo
  enddo
  call MinModGridInterp3D(dataquad,xposquad,xtarget,zz)
 endif

return
end subroutine interp_from_grid

subroutine SCI_soliddist(x,y,z,time,dist)
IMPLICIT NONE

#include "probdataf95.H"

INTEGER_T islocal,dir
REAL_T x,y,z,time,dist
REAL_T dummy_vel(SDIM)
REAL_T xtarget(SDIM)

 if (SDIM.ne.3) then
  print *,"dimension bust"
  stop
 endif

 xtarget(1)=x
 xtarget(2)=y
 xtarget(3)=z

 if (probtype.eq.538) then  ! diesel injector with moving parts

  islocal = 1  ! assume local auxiliary grid is available
  do dir=1,SDIM
   if(SCI_loc%Cellproblo(dir).ge.SCI_loc%Cellprobhi(dir)) then
    ! outside the solid domain: revert to global
    islocal = 0
   else if((xtarget(dir).le.SCI_loc%Cellproblo(dir)) .or. &
           (xtarget(dir).ge.SCI_loc%Cellprobhi(dir))) then 
    ! outside the processor domain: revert to global
    islocal = 0
   endif
  enddo  ! dir

 else
  islocal = 0
 endif


 if(islocal.eq.0) then
  call SCI_soliddist_velocity(x,y,z,time,FSI_needle,SCI,SCI_needle, &
   dist,dummy_vel)
 else if(islocal.eq.1) then
  call SCI_soliddist_velocity(x,y,z,time,FSI_needle,SCI_loc,SCI_needle_loc, &
   dist,dummy_vel)
 else
  print*,"SCI_soliddist: islocal bust ",islocal
  stop
 endif

return
end subroutine SCI_soliddist


      subroutine flappingKinematics(numMotion,motionPara,r,t)
      IMPLICIT NONE
      INTEGER_T numMotion
      REAL_T xPoint(3,numMotion),vTan(3,numMotion)
      REAL_T x0(3),v(3),motionPara(11,numMotion)
      REAL_T vNorm,theta,thetaMag,fTheta,phiTheta,theta0
      REAL_T t,h,hMag,fH,phiH,h0,ct,st
      REAL_T r(3,4),r1(3,4),r2(3,4),r0(3,4)
      INTEGER_T motionType
      INTEGER_T rotateAlongALine,translateAlongALine
      INTEGER_T i,j,k,iMotion

!     pitching motion theta(t)=theta_mag*cos(2*pi*f_theta*t+phi_theta)+theta_0
!     plunging motion h(t)    =h_mag    *cos(2*pi*f_h    *t+phi_h    )+h_0
      rotateAlongALine   =0
      translateAlongALine=1 

      do iMotion=1,numMotion 
        do i=1,3
           x0(i)=xPoint(i,iMotion)
           v(i) =vTan(i,iMotion)
        enddo
        vNorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
        v(1)=v(1)/vNorm
        v(2)=v(2)/vNorm
        v(3)=v(3)/vNorm

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


! FSI_displ is the needle geometry if probtype=538
! if probtype=701,
!  SCI_1=>SCI
!  SCI_2=>SCI_needle
!  FSI_displ is not used for 701.
subroutine SCI_soliddist_velocity(x,y,z,time, &
  FSI_displ,SCI_1,SCI_2, &
  dist,velparm)

IMPLICIT NONE

REAL_T velparm(SDIM)
REAL_T velparm2(SDIM)
REAL_T x,y,z,time,dist,dist2
INTEGER_T scomp,lsflag,dir
REAL_T xtarget(SDIM)
REAL_T xtarget2(SDIM)
REAL_T xfoot(SDIM)
REAL_T xfoot2(SDIM)
REAL_T RPM,alpha,RR,radgear,theta
type(mesh_type) :: FSI_displ
type(auxgrid_type) :: SCI_1,SCI_2

INTEGER_T numMotion
REAL_T r(3,4)
REAL_T rinv(3,4)
REAL_T rplus(3,4)
REAL_T det
REAL_T motionPara(11,2)
REAL_T flapping_time
REAL_T dt_flapping
REAL_T flapping_time_plus
INTEGER_T itime

#include "probdataf95.H"

 if (SDIM.ne.3) then
  print *,"dimension bust"
  stop
 endif

 do dir=1,SDIM
  velparm(dir)=0.0
 enddo

 dist=1.0E+10
 call init_dimensions(SCI_1)  ! ysl 09/02/14 init dims(CellLS) should we initialize it? 

 xtarget(1)=x
 xtarget(2)=y
 xtarget(3)=z
 do dir=1,SDIM
  xfoot(dir)=xtarget(dir)
  xfoot2(dir)=xfoot(dir)
 enddo

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
  ! below 538=diesel injector
 else if (probtype.eq.538) then

  if (levelrz.eq.0) then
   ! do nothing
  else if (levelrz.eq.1) then ! place inlet at center of domain
   xfoot(1)=xfoot(1)+0.0015
   xfoot(2)=xfoot(2)-0.0045
  else
   print *,"levelrz invalid"
   stop
  endif 

 else if (probtype.eq.701) then ! SCI_soliddist_velocity: first flapping wing

  numMotion=2
      
  motionPara(1,1)=0.0             ! 0: rotation; 1: translation
  motionPara(2,1)=30./180.*Pi     ! amplitude
  motionPara(3,1)=1.              ! f=frequency 2*pi*f
  motionPara(4,1)=0.              ! offset alpha0
  motionPara(5,1)=0./180.*Pi      ! phase change
  motionPara(6,1)=0.25            ! pivoting point
  motionPara(7,1)=0.        
  motionPara(8,1)=0.
  motionPara(9,1)=0.              ! rotation direction
  motionPara(10,1)=0.
  motionPara(11,1)=1.

  motionPara(1,2)=1             ! motionType 
  motionPara(2,2)=0.0           ! hMag
  motionPara(3,2)=1.            ! frequency 
  motionPara(4,2)=0.            ! hPhi
  motionPara(5,2)=0.            ! h0 
  motionPara(6,2)=0.0           ! x(1)
  motionPara(7,2)=0.            ! x(2)
  motionPara(8,2)=0.            ! x(3)
  motionPara(9,2)=0.            ! v(1)
  motionPara(10,2)=1.           ! v(2) 
  motionPara(11,2)=0.           ! v(3)

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
  rinv(1,1)=r(2,2)/det
  rinv(1,2)=-r(1,2)/det
  rinv(2,1)=-r(2,1)/det
  rinv(2,2)=r(1,1)/det
  xfoot2(1)=xtarget(1)-r(1,4)
  xfoot2(3)=xtarget(3)-r(2,4)
  xfoot(1)=rinv(1,1)*xfoot2(1)+rinv(1,2)*xfoot2(3)
  xfoot(3)=rinv(2,1)*xfoot2(1)+rinv(2,2)*xfoot2(3)

  flapping_time_plus=flapping_time+dt_flapping
  call flappingKinematics(numMotion,motionPara,rplus,flapping_time_plus)
  velparm(1)=(rplus(1,1)-r(1,1))*xfoot(1)+ &
             (rplus(1,2)-r(1,2))*xfoot(3)+ &
             rplus(1,4)-r(1,4)
  velparm(2)=0.0
  velparm(3)=(rplus(2,1)-r(2,1))*xfoot(1)+ &
             (rplus(2,2)-r(2,2))*xfoot(3)+ &
             rplus(2,4)-r(2,4)
  do dir=1,SDIM
   velparm(dir)=velparm(dir)/(2.0*dt_flapping)
  enddo
  
 endif   ! probtype.eq.701

 scomp=1
 lsflag=1
 call interp_from_grid(dist,xfoot,scomp,lsflag, &
    SCI_1%CellLS,DIMS3D(CellLS),SCI_1%cell_lo,SCI_1%cell_hi,SCI_1)

  ! needle for diesel injector (part 2)
 if (probtype.eq.538) then

   ! if injector is shifted, then needle should be shifted too.
  if (levelrz.eq.0) then
   ! do nothing
  else if (levelrz.eq.1) then ! place inlet at center of domain
   xfoot2(1)=xfoot2(1)+0.0015
   xfoot2(2)=xfoot2(2)-0.0045
  else
   print *,"levelrz invalid"
   stop
  endif 

  xfoot2(3)=xfoot2(3)+0.01
  do dir=1,SDIM
   xfoot2(dir)=xfoot2(dir)-FSI_displ%solid_displ(dir)
  enddo

  call init_dimensions(SCI_2)  ! ysl 09/02/14; should we initialize it?
  call interp_from_grid(dist2,xfoot2,scomp,lsflag, &
     SCI_2%CellLS,DIMS3D(CellLS), &
     SCI_2%cell_lo,SCI_2%cell_hi,SCI_2)

  if (dist2.lt.dist) then
   dist=dist2
   do dir=1,SDIM
    velparm(dir)=FSI_displ%solid_speed(dir)
   enddo
  endif

 else if (probtype.eq.701) then  ! 2nd flapping wing

  if (axis_dir.eq.0) then
   ! do nothing - no part 2 for flapping wing
  else if (axis_dir.eq.1) then

   do dir=1,SDIM
    xtarget2(dir)=xtarget(dir)
   enddo
   !xtarget2(1)=xtarget2(1)-one  ! translate 2nd wing to the right 1 unit.
   do dir=1,SDIM
    xfoot(dir)=xtarget2(dir)
    xfoot2(dir)=xfoot(dir)
   enddo

   numMotion=2
       
   motionPara(1,1)=0.0
   motionPara(2,1)=0./180.*Pi
   motionPara(3,1)=1.
   motionPara(4,1)=0.
   motionPara(5,1)=0./180.*Pi
   motionPara(6,1)=.0
   motionPara(7,1)=0.
   motionPara(8,1)=0.
   motionPara(9,1)=0.
   motionPara(10,1)=0.
   motionPara(11,1)=1.

   motionPara(1,2)=1             ! motionType 
   motionPara(2,2)=0.5          ! hMag
   motionPara(3,2)=1.            ! hF
   motionPara(4,2)=0.            ! hPhi
   motionPara(5,2)=0.            ! h0 
   motionPara(6,2)=0.0           ! x(1)  the default translation point should be (0,0,0) 
   motionPara(7,2)=0.            ! x(2)  otherwise it shifts to the new location
   motionPara(8,2)=0.            ! x(3)  (x1,x2,x3)
   motionPara(9,2)=0.            ! v(1)
   motionPara(10,2)=1.           ! v(2) 
   motionPara(11,2)=0.           ! v(3)

     ! NINT = round to nearest whole number 
     ! for 2nd wing motion, reverse direction of flapping motion.
   !itime=NINT(time-0.49999999999999)
   flapping_time=time

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
   xfoot2(1)=xtarget2(1)-r(1,4)
   xfoot2(3)=xtarget2(3)-r(2,4)
   xfoot(1)=rinv(1,1)*xfoot2(1)+rinv(1,2)*xfoot2(3)
   xfoot(3)=rinv(2,1)*xfoot2(1)+rinv(2,2)*xfoot2(3)

   dt_flapping=0.01
   flapping_time_plus=flapping_time+dt_flapping
   call flappingKinematics(numMotion,motionPara,rplus,flapping_time_plus)
   velparm2(1)=(rplus(1,1)-r(1,1))*xfoot(1)+ &
               (rplus(1,2)-r(1,2))*xfoot(3)+ &
               rplus(1,4)-r(1,4)
   velparm2(2)=0.0
   velparm2(3)=(rplus(2,1)-r(2,1))*xfoot(1)+ &
               (rplus(2,2)-r(2,2))*xfoot(3)+ &
               rplus(2,4)-r(2,4)
   do dir=1,SDIM
    velparm2(dir)=velparm2(dir)/dt_flapping
   enddo

   scomp=1
   lsflag=1
   call init_dimensions(SCI_2)  ! ysl 09/02/14 should we initialize it? 
   call interp_from_grid(dist2,xfoot,scomp,lsflag, &
    SCI_2%CellLS,DIMS3D(CellLS), &
    SCI_2%cell_lo,SCI_2%cell_hi,SCI_2)

!   if (xfoot(1).gt.1.45.and.xfoot(1).lt.1.5.and.abs(xfoot(3)).le.0.05) then
!       print*, xfoot(1), xfoot(3), dist2,dist
!       read(*,*) dt_flapping 
!       endif     

   if (dist2.lt.dist) then
    dist=dist2
    do dir=1,SDIM
     velparm(dir)=velparm2(dir)
    enddo
   endif

  else
   print *,"axis_dir invalid"
   stop
  endif
 endif  

return
end subroutine SCI_soliddist_velocity

subroutine SCI_velsolid(x,y,z,time,vel)
IMPLICIT NONE

REAL_T x,y,z,time,inflowvel
INTEGER_T lsflag,dir,islocal
REAL_T xtarget(SDIM)
REAL_T vel(SDIM)
REAL_T RPM,alpha,RR,radgear,theta
REAL_T dummy_dist

#include "probdataf95.H"


 do dir=1,SDIM
  vel(dir)=0.0
 enddo

 call init_dimensions(SCI)  ! the dimensions for global auxgrid
 call init_dimensions(SCI_needle) ! ysl; do we need to initialize SCI_needle? 

 xtarget(1)=x
 xtarget(2)=y
 xtarget(3)=z

 islocal = 0
 if (probtype.eq.701) then
  call SCI_soliddist_velocity(x,y,z,time,FSI_needle,&
    SCI,SCI_needle,dummy_dist,vel)
 else if (probtype.eq.538) then  ! diesel injector with moving parts

  islocal = 1  ! assume local auxiliary grid is available
  do dir=1,SDIM
   if(SCI_loc%Cellproblo(dir).ge.SCI_loc%Cellprobhi(dir)) then
    ! outside the solid domain: revert to global
    islocal = 0
   else if((xtarget(dir).le.SCI_loc%Cellproblo(dir)) .or. &
           (xtarget(dir).ge.SCI_loc%Cellprobhi(dir))) then 
    ! outside the processor domain: revert to global
    islocal = 0
   endif
  enddo  ! dir

  if(islocal.eq.0) then
   call SCI_soliddist_velocity(x,y,z,time,FSI_needle,&
    SCI,SCI_needle,dummy_dist,vel)
  else if(islocal.eq.1) then
   call SCI_soliddist_velocity(x,y,z,time,FSI_needle,&
    SCI_loc,SCI_needle_loc,dummy_dist,vel)
  else
   print*,"SCI_velsolid: islocal bust ",islocal
   stop
  endif

 else if (probtype.ne.538) then
   ! SCI is the global auxgrid
  call init_dimensions(SCI)
  
  lsflag=0
  do dir=1,SDIM
   call interp_from_grid(vel(dir),xtarget,dir,lsflag, &
     SCI%CellVEL,DIMS3D(CellLS),SCI%cell_lo,SCI%cell_hi,SCI)
  enddo

   ! viorel sphere problem
  if (probtype.eq.5601) then
   inflowvel=0.0
   if ((sqrt(xtarget(2)**2+xtarget(3)**2).le.2.0).and. &
       (abs(xtarget(1)+10.0).le.1.0)) then
    inflowvel=1.0
   endif
   if ((sqrt(xtarget(2)**2+(xtarget(3)-3.0)**2).le.2.0).and. &
       (abs(xtarget(1)-10.0).le.1.0)) then
    inflowvel=1.0
   endif
   vel(1)=inflowvel
  endif  

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

    vel(1)=alpha*xtarget(3)
    vel(2)=zero
    vel(3)=-alpha*xtarget(1)
  endif  ! probtype=563

 else
  print *,"probtype bust"
  stop
 endif   

return
end subroutine SCI_velsolid


subroutine SCI_tempsolid(x,y,z,time,temperature)
IMPLICIT NONE

REAL_T x,y,z,time
INTEGER_T scomp,lsflag
REAL_T xtarget(SDIM)
REAL_T temperature

#include "probdataf95.H"

 if (SDIM.ne.3) then
  print *,"dimension bust"
  stop
 endif

 temperature=0.0

   ! SCI is the global auxgrid
 call init_dimensions(SCI)

 xtarget(1)=x
 xtarget(2)=y
 xtarget(3)=z

 lsflag=0
 scomp=1
 call interp_from_grid(temperature,xtarget,scomp,lsflag, &
    SCI%CellTEMP,DIMS3D(CellLS),SCI%cell_lo,SCI%cell_hi,SCI)

return
end subroutine SCI_tempsolid



! This routine only works for dimension by dimension stretching.
! The cell locations are monotonically increasing.


subroutine box_bounds(cellproblo,celldx, &
  lo,hi,xlo,xhi)
IMPLICIT NONE

  INTEGER_T dir
  INTEGER_T lo(SDIM),hi(SDIM)
  REAL_T xlo(SDIM),xhi(SDIM)
  REAL_T cellproblo(SDIM)
  REAL_T celldx(SDIM)

  do dir=1,3
   xlo(dir)=cellproblo(dir)
   xhi(dir)=cellproblo(dir)+(hi(dir)-lo(dir)+one)*celldx(dir)
  enddo
  

return
end subroutine box_bounds


subroutine find_grid_bounding_box(minnode,maxnode, &
  cellproblo,celldx, &
  lo,hi,gridlo,gridhi,iband)
IMPLICIT NONE

 REAL_T minnode(SDIM),maxnode(SDIM)
 INTEGER_T lo(SDIM),hi(SDIM)
 INTEGER_T gridlo(SDIM),gridhi(SDIM)
 INTEGER_T iband,dir
 REAL_T cellproblo(SDIM)
 REAL_T celldx(SDIM)
 REAL_T xlo(SDIM),xhi(SDIM)

 if (iband.lt.0) then
  print *,"iband bust"
  stop
 endif

 do dir=1,SDIM
  xlo(dir)=cellproblo(dir)
  xhi(dir)=cellproblo(dir)+(hi(dir)-lo(dir)+one)*celldx(dir)
 enddo

! NINT = round to nearest whole number
! x=(i-lo)dx+xlo+dx/2=(i-lo+1/2)dx+xlo
! i=(x-xlo)/dx+lo-1/2

 do dir=1,SDIM
  gridlo(dir)=NINT( (minnode(dir)-xlo(dir))/celldx(dir)-half+lo(dir) )-iband
  gridhi(dir)=NINT( (maxnode(dir)-xlo(dir))/celldx(dir)-half+lo(dir) )+iband
 enddo 

return
end subroutine find_grid_bounding_box


subroutine CLSVOF_ReadNodes(CLSVOF_curtime,CLSVOF_dt,h_small, &
  problo,probhi,ioproc,isout)
IMPLICIT NONE

  INTEGER_T :: initflag,ioproc,isout
  REAL_T :: CLSVOF_curtime,CLSVOF_dt,h_small
  REAL_T problo(SDIM),probhi(SDIM)

#include "probdataf95.H"

  DebugCouple=.true.
  FSI%PartID = 1

   ! ifirst=0 (iread=0) when initinjector called (if probtype=538)
   ! initializes FSI%solid_displ and FSI%solid_speed
   ! if probtype=701 then initflapping called, but nothing done.
  call overall_solid_advance(CLSVOF_curtime,CLSVOF_dt,FSI, &
     ioproc,isout)

  initflag=0

   ! for neither probtype=538 or probtype=701, is "generate_new_triangles"
   ! called from this routine; it would be inefficient since solid
   ! motion is known to be rigid body motion in these two cases.
  if (probtype.eq.701) then

   if (axis_dir.eq.0) then
    ! do nothing (flapping wing - just 1 part)
   else if (axis_dir.eq.1) then ! 2nd wing - ReadNodes
    FSI_needle%PartID = 2
    call overall_solid_advance(CLSVOF_curtime,CLSVOF_dt,FSI_needle, &
      ioproc,isout)
   else
    print *,"axis_dir invalid"
    stop
   endif

  else if (probtype.ne.538) then

   call generate_new_triangles(initflag,problo,probhi,FSI,SCI, &
    0,ioproc,isout)

  else if (probtype.eq.538) then ! moving parts (needle)

    FSI_needle%PartID = 2

      ! ifirst=0 (iread=0)
      ! initializes FSI_needle%solid_displ and 
      ! FSI_needle%solid_speed
      ! routines that use local auxgrid must use the global
      ! FSI_needle%solid_displ and
      ! FSI_needle%solid_speed
    call overall_solid_advance(CLSVOF_curtime,CLSVOF_dt,FSI_needle, &
      ioproc,isout)

  else
   print *,"probtype bust"
   stop
  endif
 
return
end subroutine CLSVOF_ReadNodes

!!!!!!!! Austen Duffy's routines here !!!!!!
! probtype=562 (whale)
! set axis_dir=4 (animated whalenormal.txt)

      SUBROUTINE getinfo(whalein,Nodes,whaleout)
      
      
      Character (len=*) whalein, whaleout
      INTEGER_T Nodes
      
      whaleout=whalein
      
      print *, whaleout
      
      open(unit=4,file=whalein)
      
      read(4,*) Nodes
      
      close(4)
      
      END SUBROUTINE
   
      SUBROUTINE runonce(R,S,T,U,V,W,X,Y,Z,List,angle,timestep,spring, &
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

      CHARACTER*35 filename, whaleout
      INTEGER Nodes, Cells, Shape, cycles
      INTEGER_T, DIMENSION(Nodes,20) :: List
      INTEGER i, j, k, n1, n2, n3, garbage, p, counter, count, n
      REAL_T, Dimension(Nodes) ::  X, Y, Z, R, S, T, X_init, Y_init
      REAL_T, Dimension(Nodes) ::  spring, Z_init
      REAL_T, Dimension(Nodes) ::  U, V, W, V_init, temp_S, temp_T
      REAL_T, Dimension(Nodes) ::  temp_U,temp_V,temp_W,temp_R
      REAL_T, Dimension(22) :: angle, timestep
      REAL_T  a, b, c, sum_x, sum_y, sum_z, mu, nu, value,temp
      REAL_T  counter_real, tailpos, mag, mag_init,L1,L2,L4,L5
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
  
      FSI%IntElemDim=3 
      FSI%NumIntElems=Cells
      allocate(FSI%IntElem(FSI%IntElemDim,FSI%NumIntElems))

      DO i=1,Cells-1
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
         write(16,*) n1, n2, n3

         FSI%IntElem(1,i)=n1
         FSI%IntElem(2,i)=n2
         FSI%IntElem(3,i)=n3

         read(2,*) garbage
         read(2,*) garbage
      END DO 

      read(2,*) n1
      read(2,*) n2
      read(2,*) n3
      write(16,*) n1, n2, n3

      FSI%IntElem(1,Cells)=n1
      FSI%IntElem(2,Cells)=n2
      FSI%IntElem(3,Cells)=n3

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
      
      
   
      SUBROUTINE new_geometry(R,S,T,U,V,W,X,Y,Z,List,angle,timestep, &
                    spring,counter,counter_real,Nodes,Cells, &
                    X_init,Y_init,Z_init)

      IMPLICIT NONE

      CHARACTER*35 filename
      INTEGER Nodes, Cells, Shape, cycles
      INTEGER_T, DIMENSION(Nodes,20) :: List
      INTEGER i, j, k, n1, n2, n3, garbage, p, counter, count, n
      INTEGER_T count1, count2, count3,count4
      REAL_T, Dimension(Nodes) ::  X, Y, Z, R, S, T, X_init, Y_init
      REAL_T, Dimension(Nodes) ::  spring, Z_init
      REAL_T, Dimension(Nodes) ::  U, V, W, V_init, temp_S, temp_T
      REAL_T, Dimension(Nodes) ::  temp_U,temp_V,temp_W,temp_R
      REAL_T, Dimension(22) :: angle, timestep
      REAL_T  a, b, c, sum_x, sum_y, sum_z, mu, nu, value,temp
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

      
      END SUBROUTINE
    
      SUBROUTINE nodelist(Nodes, Ext, List,whaleout)

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

      END SUBROUTINE
         
      SUBROUTINE plininterp(xvect, yvect, x, n, value)

      IMPLICIT NONE

      INTEGER_T i, j, k, n, m
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


      END SUBROUTINE





      SUBROUTINE tailup(R,S,T,U,V,W,X,Y,Z, Nodes, tailpos)

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

      
      END SUBROUTINE
      




      SUBROUTINE taildown(R,S,T,U,V,W,X,Y,Z,Nodes,tailpos)

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

      END SUBROUTINE





      SUBROUTINE springs(Z, spring, Nodes,tailpos)


      IMPLICIT NONE

      INTEGER_T Nodes, i
      REAL_T tailpos
      REAL_T, Dimension(Nodes) :: spring, Z
      REAL_T const1, const2,const3,springcons


      springcons=0.5

      DO i=1,Nodes
           spring(i)=0.1

          If (Z(i).gt.tailpos) Then
            spring(i)=spring(i)+springcons*(10.0-Z(i))**2 &
                     -springcons*(10.0-tailpos)**2
          End If

      END DO
  
            


      END SUBROUTINE

end module CLSVOFCouplerIO


