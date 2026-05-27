#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"

#include "EXTRAP_COMP.H"

#define NOTUS_REAL (8)

#define MOFDEB (0)

#define MAXTET (5)
#define MAXAREA (5)
#define USERAND (1)

! Author: Mark Sussman sussman@math.fsu.edu
! Department of Mathematics
! Florida State University
! Tallahassee, FL 32306
!

module geometry_intersect_module
use amrex_fort_module, only : amrex_real

implicit none

integer, PARAMETER :: INTERCEPT_MAXITER=100
integer, PARAMETER :: INTERCEPT_MAXITER_NEWTON=25

integer, PARAMETER :: maxfacelist=20
integer, PARAMETER :: maxnodelist=20,maxtetlist=15,maxcapfacelist=5
integer, PARAMETER :: maxmappednodes=8
integer, PARAMETER :: n_vol_listmax=400,n_area_listmax=200

integer :: geom_nthreads
integer :: geom_nmax
 ! 4,3,geom_nmax,geom_nthreads
real(amrex_real), ALLOCATABLE :: geom_xtetlist_local(:,:,:,:)
real(amrex_real), ALLOCATABLE :: geom_xtetlist(:,:,:,:)
real(amrex_real), ALLOCATABLE :: geom_xtetlist_old(:,:,:,:)
real(amrex_real), ALLOCATABLE :: geom_xtetlist_uncapt(:,:,:,:)
integer, ALLOCATABLE :: mof_calls(:,:)
integer, ALLOCATABLE :: mof_iterations(:,:)
real(amrex_real), ALLOCATABLE :: mof_errors(:,:)

real(amrex_real), ALLOCATABLE :: intercept_error_history(:,:)

! the 3rd component of a facelist is 0 in 2D
type intersect_type
  integer :: n_nodes,n_tet,n_faces,n_capfaces
  integer :: n_pos_nodes
  integer :: nodelistmap(maxnodelist,3)  
  integer :: nodelist(maxnodelist)
  integer :: tetlist(maxtetlist,4)
  integer :: facelist(maxfacelist,3)
  integer :: aligned(maxfacelist)
  integer :: capfacelist(maxcapfacelist,3)
end type intersect_type

type mapping_type
  integer :: mapped_nodes(maxmappednodes)
  type(intersect_type) :: intersect_geometry
end type mapping_type

! check sum and headnode determine geometry map to use
type(mapping_type), dimension(256,8)  :: hexahedron_maps
type(mapping_type), dimension(16,4)    :: tetrahedron_maps
type(mapping_type), dimension(16,4)    :: rectangle_maps
type(mapping_type), dimension(8,3)    :: triangle_maps

type(intersect_type), dimension(9) :: template_hex_plane
type(intersect_type), dimension(4) :: template_tet_plane
type(intersect_type), dimension(4) :: template_rec_plane
type(intersect_type), dimension(3) :: template_tri_plane


 !e.g. phi=z-eta(x,y)
type levelset_parm_type
integer :: ROOT_SDIM ! 1,2,3,.... ROOT_SDIM>=LOCAL_SDIM
integer :: LOCAL_SDIM ! 1,2,3,....
integer :: height_dir !=0,1,2,3,....  height_dir<ROOT_SDIM
integer, pointer :: frozen_parms_flag(:) ! 1...ROOT_SDIM
real(amrex_real), pointer :: frozen_parms(:) ! 1...ROOT_SDIM
integer :: signLS
integer :: activeLS
integer :: order
!e.g. sum a_ij x^i y^j i,j=0..order
real(amrex_real), pointer, dimension(:) :: LSCOEFF_FLATTEN 
end type levelset_parm_type

type levelset_array_parm_type
integer :: nLS
integer :: xkL
integer :: xkU
integer :: k_reduce
integer :: S_quad_type
type(levelset_parm_type), pointer, dimension(:) :: LS_array !1..nLS
end type levelset_array_parm_type

type function_parm_type
integer :: ROOT_SDIM ! 1,2,3,.... ROOT_SDIM>=LOCAL_SDIM
integer :: LOCAL_SDIM ! 1,2,3,....
integer :: order
!e.g. sum a_{i,j,k} x^i y^j z^k i,j,k=0..order
real(amrex_real), pointer, dimension(:) :: FNCOEFF_FLATTEN
! 1..ROOT_SDIM
type(levelset_array_parm_type), pointer, dimension(:) :: LS_FN_array 
end type function_parm_type

contains

subroutine copy_intersect_type(source,dest)
IMPLICIT NONE

type(intersect_type), INTENT(in) :: source
type(intersect_type), INTENT(out) :: dest

integer i,dir

 dest%n_nodes=source%n_nodes
 dest%n_tet=source%n_tet
 dest%n_faces=source%n_faces
 dest%n_capfaces=source%n_capfaces
 dest%n_pos_nodes=source%n_pos_nodes
 do i=1,maxnodelist
  do dir=1,3
   dest%nodelistmap(i,dir)=0
  enddo
  dest%nodelist(i)=0
 enddo
 do i=1,maxtetlist
  do dir=1,4
   dest%tetlist(i,dir)=0
  enddo
 enddo
 do i=1,maxfacelist
  do dir=1,3
   dest%facelist(i,dir)=0
  enddo
  dest%aligned(i)=0
 enddo
 do i=1,maxcapfacelist
  do dir=1,3
   dest%capfacelist(i,dir)=0
  enddo
 enddo

 do i=1,source%n_nodes
  do dir=1,3
   dest%nodelistmap(i,dir)=source%nodelistmap(i,dir)
  enddo
  dest%nodelist(i)=source%nodelist(i)
 enddo
 do i=1,source%n_tet
  do dir=1,4
   dest%tetlist(i,dir)=source%tetlist(i,dir)
  enddo
 enddo
 do i=1,source%n_faces
  do dir=1,3
   dest%facelist(i,dir)=source%facelist(i,dir)
  enddo
  dest%aligned(i)=source%aligned(i)
 enddo
 do i=1,source%n_capfaces
  do dir=1,3
   dest%capfacelist(i,dir)=source%capfacelist(i,dir)
  enddo
 enddo

return
end subroutine copy_intersect_type


subroutine fast_copy_intersect_type(source,dest,sdim)
IMPLICIT NONE

type(intersect_type), INTENT(in) :: source
type(intersect_type), INTENT(out) :: dest
integer, INTENT(in) :: sdim
integer :: i,dir

 if ((sdim.ne.2).and.(sdim.ne.3)) then
  print *,"sdim invalid"
  stop
 endif

 dest%n_nodes=source%n_nodes
 dest%n_tet=source%n_tet
 dest%n_capfaces=source%n_capfaces

 do i=1,source%n_nodes
  do dir=2,3
   dest%nodelistmap(i,dir)=source%nodelistmap(i,dir)
  enddo
 enddo
 do i=1,source%n_tet
  do dir=1,sdim+1
   dest%tetlist(i,dir)=source%tetlist(i,dir)
  enddo
 enddo
 do i=1,source%n_capfaces
  do dir=1,sdim
   dest%capfacelist(i,dir)=source%capfacelist(i,dir)
  enddo
 enddo

return
end subroutine fast_copy_intersect_type



subroutine add_to_hex(gridmap,template_geom)
IMPLICIT NONE

type(intersect_type), INTENT(in) :: template_geom
integer, INTENT(in) :: gridmap(2,2,2)
integer :: mapped_nodes(maxmappednodes)
integer :: power2(maxmappednodes)
integer :: nn,n_nodes,inode,checksum
integer :: igrid_node
integer :: jgrid_node
integer :: kgrid_node
integer :: i_power
integer :: i_map

 power2(1)=1
 do i_power=2,maxmappednodes
  power2(i_power)=2*power2(i_power-1)
 enddo

 do i_map=1,maxmappednodes
  mapped_nodes(i_map)=0
 enddo

 n_nodes=template_geom%n_nodes
 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n_nodes invalid"
  stop
 endif

 inode=1
 checksum=0
 do kgrid_node=1,2
 do jgrid_node=1,2
 do igrid_node=1,2
  mapped_nodes(inode)=gridmap(igrid_node,jgrid_node,kgrid_node)
  if ((mapped_nodes(inode).lt.1).or.(mapped_nodes(inode).gt.8)) then
   print *,"mapped nodes bust"
   stop
  endif

  do nn=1,n_nodes
   if ((template_geom%nodelistmap(nn,2).eq.inode).and. &
       (template_geom%nodelistmap(nn,3).eq.inode)) then
    checksum=checksum+power2(mapped_nodes(inode))
   endif
  enddo

  inode=inode+1
 enddo
 enddo
 enddo
 if ((checksum.lt.1).or.(checksum.gt.255)) then
  print *,"checksum invalid"
  stop
 endif
 call copy_intersect_type(template_geom, &
   hexahedron_maps(checksum,mapped_nodes(1))%intersect_geometry)
 do i_map=1,maxmappednodes
  hexahedron_maps(checksum,mapped_nodes(1))%mapped_nodes(i_map)= &
          mapped_nodes(i_map)
 enddo

return
end subroutine add_to_hex


subroutine add_to_tet(linemap,template_geom)
IMPLICIT NONE

type(intersect_type), INTENT(in) :: template_geom
integer, INTENT(in) :: linemap(4)
integer :: mapped_nodes(maxmappednodes)
integer :: power2(maxmappednodes)
integer :: nn,n_nodes,inode,checksum
integer :: i_power
integer :: i_map

 power2(1)=1
 do i_power=2,maxmappednodes
  power2(i_power)=2*power2(i_power-1)
 enddo

 do i_map=1,maxmappednodes
  mapped_nodes(i_map)=0
 enddo

 n_nodes=template_geom%n_nodes
 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n_nodes invalid"
  stop
 endif

 inode=1
 checksum=0
 do inode=1,4
  mapped_nodes(inode)=linemap(inode)
  if ((mapped_nodes(inode).lt.1).or.(mapped_nodes(inode).gt.4)) then
   print *,"mapped nodes bust"
   stop
  endif

  do nn=1,n_nodes
   if ((template_geom%nodelistmap(nn,2).eq.inode).and. &
       (template_geom%nodelistmap(nn,3).eq.inode)) then
    checksum=checksum+power2(mapped_nodes(inode))
   endif
  enddo

 enddo

 if ((checksum.lt.1).or.(checksum.gt.15)) then
  print *,"checksum invalid"
  stop
 endif
 call copy_intersect_type(template_geom, &
   tetrahedron_maps(checksum,mapped_nodes(1))%intersect_geometry)
 do i_map=1,maxmappednodes
  tetrahedron_maps(checksum,mapped_nodes(1))%mapped_nodes(i_map)= &
          mapped_nodes(i_map)
 enddo

return
end subroutine add_to_tet


subroutine add_to_tri(linemap,template_geom)
IMPLICIT NONE

type(intersect_type), INTENT(in) :: template_geom
integer, INTENT(in) :: linemap(3)
integer :: mapped_nodes(maxmappednodes)
integer :: power2(maxmappednodes)
integer :: nn,n_nodes,inode,checksum
integer :: i_power
integer :: i_map

 power2(1)=1
 do i_power=2,maxmappednodes
  power2(i_power)=2*power2(i_power-1)
 enddo

 do i_map=1,maxmappednodes
  mapped_nodes(i_map)=0
 enddo

 n_nodes=template_geom%n_nodes
 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n_nodes invalid"
  stop
 endif

 inode=1
 checksum=0
 do inode=1,3
  mapped_nodes(inode)=linemap(inode)
  if ((mapped_nodes(inode).lt.1).or.(mapped_nodes(inode).gt.3)) then
   print *,"mapped nodes bust"
   stop
  endif

  do nn=1,n_nodes
   if ((template_geom%nodelistmap(nn,2).eq.inode).and. &
       (template_geom%nodelistmap(nn,3).eq.inode)) then
    checksum=checksum+power2(mapped_nodes(inode))
   endif
  enddo

 enddo

 if ((checksum.lt.1).or.(checksum.gt.7)) then
  print *,"checksum invalid"
  stop
 endif
 call copy_intersect_type(template_geom, &
   triangle_maps(checksum,mapped_nodes(1))%intersect_geometry)
 do i_map=1,maxmappednodes
  triangle_maps(checksum,mapped_nodes(1))%mapped_nodes(i_map)= &
          mapped_nodes(i_map)
 enddo

return
end subroutine add_to_tri



subroutine add_to_rec(gridmap,template_geom)
IMPLICIT NONE

type(intersect_type), INTENT(in) :: template_geom
integer, INTENT(in) :: gridmap(2,2)
integer :: mapped_nodes(maxmappednodes)
integer :: power2(maxmappednodes)
integer :: nn,n_nodes,inode,checksum
integer :: igrid_node
integer :: jgrid_node
integer :: i_power
integer :: i_map

 if (maxmappednodes.lt.4) then
  print *,"maxmappednodes invalid"
  stop
 endif

 power2(1)=1
 do i_power=2,maxmappednodes
  power2(i_power)=2*power2(i_power-1)
 enddo

 do i_map=1,maxmappednodes
  mapped_nodes(i_map)=0
 enddo

 n_nodes=template_geom%n_nodes
 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n_nodes invalid"
  stop
 endif

 inode=1
 checksum=0
 do jgrid_node=1,2
 do igrid_node=1,2
  mapped_nodes(inode)=gridmap(igrid_node,jgrid_node)
  if ((mapped_nodes(inode).lt.1).or.(mapped_nodes(inode).gt.4)) then
   print *,"mapped nodes bust"
   stop
  endif

  do nn=1,n_nodes
   if ((template_geom%nodelistmap(nn,2).eq.inode).and. &
       (template_geom%nodelistmap(nn,3).eq.inode)) then
    checksum=checksum+power2(mapped_nodes(inode))
   endif
  enddo

  inode=inode+1
 enddo
 enddo

 if ((checksum.lt.1).or.(checksum.gt.15)) then
  print *,"checksum invalid"
  stop
 endif
 if (template_geom%n_nodes.le.0) then
  print *,"n_nodes bust"
  stop
 endif
 if (1.eq.0) then
  print *,"init rectangle checksum,map(1) ",checksum,mapped_nodes(1)
  print *,"n_pos_nodes ",template_geom%n_pos_nodes
 endif

 call copy_intersect_type(template_geom, &
   rectangle_maps(checksum,mapped_nodes(1))%intersect_geometry)
 do i_map=1,maxmappednodes
  rectangle_maps(checksum,mapped_nodes(1))%mapped_nodes(i_map)= &
          mapped_nodes(i_map)
 enddo

return
end subroutine add_to_rec

subroutine init_intersect_type(template_geom,n_nodes,n_faces, &
  n_capfaces,n_pos_nodes,node_array,face_array, &
  aligned_array,capface_array,sdim)
IMPLICIT NONE

integer, INTENT(in) :: n_nodes,n_faces,n_capfaces,n_pos_nodes,sdim
integer, INTENT(in) :: node_array(n_nodes)
integer, INTENT(in) :: face_array(3*n_faces)
integer, INTENT(in) :: aligned_array(n_faces)
integer, INTENT(in) :: capface_array(3*n_capfaces)
integer :: i,j,icomp,dir,firstnode,secondnode
integer :: rawnode,found,itet

type(intersect_type), INTENT(out) :: template_geom

 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n nodes invalid"
  stop
 endif
 if ((n_faces.lt.1).or.(n_faces.gt.maxfacelist)) then
  print *,"n faces invalid"
  stop
 endif
 if ((n_capfaces.lt.0).or.(n_capfaces.gt.maxcapfacelist)) then
  print *,"n faces invalid"
  stop
 endif
 if ((n_pos_nodes.lt.1).or.(n_pos_nodes.gt.n_nodes)) then
  print *,"n nodes invalid"
  stop
 endif


 template_geom%n_nodes=n_nodes
 template_geom%n_faces=n_faces
 template_geom%n_capfaces=n_capfaces
 template_geom%n_pos_nodes=n_pos_nodes
 template_geom%n_tet=0

 if ((sdim.ne.2).and.(sdim.ne.3)) then
  print *,"sdim invalid"
  stop
 endif

 if (node_array(1).ne.11) then
  print *,"node_array(1) should be 11"
  stop
 endif

 do i=1,n_nodes
  template_geom%nodelist(i)=node_array(i)
  template_geom%nodelistmap(i,1)=node_array(i)
  firstnode=node_array(i)/10
  secondnode=node_array(i)-10*firstnode
  template_geom%nodelistmap(i,2)=firstnode
  template_geom%nodelistmap(i,3)=secondnode
 enddo

 icomp=1
 do i=1,n_faces
  do dir=1,3

   template_geom%facelist(i,dir)=0

   if (dir.le.sdim) then

    rawnode=face_array(icomp)
    found=0
    do j=1,n_nodes
     if (rawnode.eq.node_array(j)) then
      found=1
      template_geom%facelist(i,dir)=j
     endif
    enddo
    if (found.ne.1) then
     print *,"could not find node in list"
     stop
    endif

   endif ! dir<=sdim

   icomp=icomp+1
  enddo  ! dir
 enddo 

 icomp=1
 do i=1,n_capfaces
  do dir=1,3

   template_geom%capfacelist(i,dir)=0

   if (dir.le.sdim) then

    rawnode=capface_array(icomp)
    found=0
    do j=1,n_nodes
     if (rawnode.eq.node_array(j)) then
      found=1
      template_geom%capfacelist(i,dir)=j
     endif
    enddo
    if (found.ne.1) then
     print *,"could not find node in list"
     stop
    endif

   endif ! dir<=sdim

   icomp=icomp+1
  enddo  ! dir
 enddo ! i

 if (n_nodes.gt.0) then
  template_geom%n_tet=0
  itet=0

   ! first add all the tets with a capface base
  do i=1,n_capfaces
   itet=itet+1
   if ((itet.lt.1).or.(itet.gt.maxtetlist)) then
    print *,"itet invalid"
    stop
   endif
   template_geom%tetlist(itet,1)=1
   do dir=1,3
    template_geom%tetlist(itet,dir+1)=0
   enddo
   do dir=1,sdim
    template_geom%tetlist(itet,dir+1)=template_geom%capfacelist(i,dir)
   enddo
  enddo

  do i=1,n_faces
   template_geom%aligned(i)=aligned_array(i)
   if (aligned_array(i).eq.0) then
    itet=itet+1
    if ((itet.lt.1).or.(itet.gt.maxtetlist)) then
     print *,"itet invalid"
     stop
    endif
    template_geom%tetlist(itet,1)=1
    do dir=1,3
     template_geom%tetlist(itet,dir+1)=0
    enddo
    do dir=1,sdim
     template_geom%tetlist(itet,dir+1)=template_geom%facelist(i,dir)
    enddo
   else if (aligned_array(i).ne.1) then
    print *,"aligned array invalid"
    stop
   endif
  enddo  ! i=1,n_faces

  template_geom%n_tet=itet
 else 
  print *,"n_nodes invalid"
  stop
 endif

return 
end subroutine init_intersect_type
 

subroutine init_geometry_tables()
IMPLICIT NONE

integer i,j,k,n,sdim,n_nodes,n_faces,n_capfaces,n_pos_nodes
integer ii,jj,kk
integer node_array3(3)
integer node_array4(4)
integer node_array5(5)
integer node_array6(6)
integer node_array8(8)
integer node_array10(10)
integer node_array11(11)
integer face_array2(6)
integer face_array3(9)
integer face_array4(12)
integer face_array6(18)
integer face_array7(21)
integer face_array10(30)
integer face_array12(36)
integer face_array13(39)
integer face_array14(42)
integer face_array15(45)
integer capface_array1(3)
integer capface_array2(6)
integer capface_array3(9)
integer capface_array4(12)
integer aligned_array2(2)
integer aligned_array3(3)
integer aligned_array4(4)
integer aligned_array6(6)
integer aligned_array7(7)
integer aligned_array10(10)
integer aligned_array12(12)
integer aligned_array13(13)
integer aligned_array14(14)
integer aligned_array15(15)

integer inode
integer basegridmap(2,2,2)
integer basegridmap2d(2,2)
integer baselinemap(4)
integer rotmap(2,2,2,0:3)
integer rotmap2d(2,2,0:3)
integer flipx,flipy,flipz
integer rotatex,rotatey,rotatez
integer i1,i2,i3,i4,npos,nrot

 do i=1,256
 do j=1,8
  hexahedron_maps(i,j)%intersect_geometry%n_nodes=0
  do n=1,maxmappednodes
   hexahedron_maps(i,j)%mapped_nodes(n)=0
  enddo
 enddo
 enddo
 do i=1,16
 do j=1,4
  tetrahedron_maps(i,j)%intersect_geometry%n_nodes=0
  rectangle_maps(i,j)%intersect_geometry%n_nodes=0
  do n=1,maxmappednodes
   tetrahedron_maps(i,j)%mapped_nodes(n)=0
   rectangle_maps(i,j)%mapped_nodes(n)=0
  enddo
 enddo
 enddo
 do i=1,8
 do j=1,3
  triangle_maps(i,j)%intersect_geometry%n_nodes=0
  do n=1,maxmappednodes
   triangle_maps(i,j)%mapped_nodes(n)=0
  enddo
 enddo
 enddo
 
 sdim=3
   
 n_nodes=4
 n_faces=3
 n_capfaces=1
 n_pos_nodes=1
 node_array4=(/ 11,12,13,15 /)
 face_array3=(/ 11,13,15, 11,12,15, 11,12,13 /)
 aligned_array3=(/ 1,1,1 /)  
 capface_array1=(/ 12,13,15 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array3, &
   aligned_array3,capface_array1,sdim)

 n_nodes=6
 n_faces=6
 n_capfaces=2
 n_pos_nodes=2
 node_array6=(/ 11,22,13,15,24,26 /)
 face_array6=(/ 22,26,24, 11,15,13, 11,22,24, 11,13,24, &
  11,26,15, 11,26,22 /)
 aligned_array6=(/ 0,1,1,1,1,1 /)
 capface_array2=(/ 26,15,24, 15,24,13 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array6,face_array6, &
   aligned_array6,capface_array2,sdim)

 n_nodes=8
 n_faces=6
 n_capfaces=3
 n_pos_nodes=3
 node_array8=(/ 11,22,33,15,24,34,26,37 /)
 face_array6=(/ 11,26,22, 11,26,15, 11,37,33, 11,37,15, &
   33,34,37, 24,26,22 /)
 aligned_array6=(/ 1,1,1,1,0,0 /)
 capface_array3=(/ 24,26,15, 24,34,15, 34,37,15 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array8,face_array6, &
   aligned_array6,capface_array3,sdim)

 n_nodes=8
 n_faces=10
 n_capfaces=2
 n_pos_nodes=4
 node_array8=(/ 11,22,33,44,26,15,37,48 /)
 face_array10=(/ 11,26,22, 11,26,15, 11,37,33, 11,37,15, &
                 11,44,22, 11,44,33, 33,48,44, 33,48,37, &
                 22,48,44, 22,48,26 /)
 aligned_array10=(/ 1,1,1,1,1,1,0,0,0,0 /)
 capface_array2=(/ 15,48,26, 15,48,37 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array8,face_array10, &
   aligned_array10,capface_array2,sdim)

 n_nodes=10
 n_faces=12
 n_capfaces=4
 n_pos_nodes=4
 node_array10=(/ 11,22,33,55,34,37,57,24,26,56 /)
 face_array12=(/ 11,34,33, 11,34,24, 11,24,22, &
     11,37,33, 11,37,57, 11,57,55, &
     11,56,55, 11,56,26, 11,26,22, &
     55,56,57, 33,34,37, 24,22,26 /)
 aligned_array12=(/ 1,1,1,1,1,1,1,1,1,0,0,0 /)
 capface_array4=(/ 34,37,57, 34,57,24, 24,26,57, 57,26,56 /)

   ! 9th element will have this extra shape
 call init_intersect_type(template_hex_plane(9),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array10,face_array12, &
   aligned_array12,capface_array4,sdim)

 n_nodes=11
 n_faces=13
 n_capfaces=3
 n_pos_nodes=5
 node_array11=(/ 11,22,33,44,55,37,48,26,56,57,37 /)
 face_array13=(/ 11,44,22, 11,44,33, 11,57,55, 11,57,37, 11,37,33, &
                 11,56,55, 11,56,26, 11,26,22, 44,37,48, 44,37,33, &
                 44,26,48, 44,26,22, 55,57,56 /)

 aligned_array13=(/ 1,1,1,1,1,1,1,1,0,0,0,0,0 /)
 capface_array3=(/ 48,57,37, 48,57,56, 48,56,26 /)
 
 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array11,face_array13, &
   aligned_array13,capface_array3,sdim)

 n_nodes=10
 n_faces=14
 n_capfaces=2
 n_pos_nodes=6
 node_array10=(/ 11,22,33,44,55,66,68,48,57,37 /)
 face_array14=(/ 11,44,33, 11,44,22, 11,66,55, 11,66,22, &
                 44,37,33, 44,37,48, 11,37,33, 11,37,57, 11,57,55, &
                 22,48,44, 22,48,68, 22,68,66, 55,68,66, 55,68,57 /)

 aligned_array14=(/ 1,1,1,1,0,0,1,1,1,0,0,0,0,0 /)
 capface_array2=(/ 48,57,37, 48,57,68 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array10,face_array14, &
   aligned_array14,capface_array2,sdim)

 n_nodes=10
 n_faces=15
 n_capfaces=1
 n_pos_nodes=7
 node_array10=(/ 11,22,33,44,55,66,77,68,48,78 /)
 face_array15=(/ 11,44,22, 11,44,33, 55,68,66, 55,68,78, 55,78,77, &
                 11,66,55, 11,66,22, 11,77,55, 11,77,33, &
                 22,68,66, 22,68,48, 22,48,44, &
                 33,78,77, 33,78,48, 33,48,44 /)

 aligned_array15=(/ 1,1,0,0,0,1,1,1,1,0,0,0,0,0,0 /)
 capface_array1=(/ 68,48,78 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array10,face_array15, &
   aligned_array15,capface_array1,sdim)

 n_nodes=8
 n_faces=12
 n_capfaces=0
 n_pos_nodes=8
 node_array8=(/ 11,22,33,44,55,66,77,88 /)
 face_array12=(/ 11,44,22, 11,44,33, 55,88,77, 55,88,66, &
   11,66,55, 11,66,22, 33,88,77, 33,88,44,  &
   11,77,33, 11,77,55, 22,88,44, 22,88,66 /)
 aligned_array12=(/ 1,1,0,0,1,1,0,0,1,1,0,0 /)
 capface_array1=(/ 0,0,0 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array8,face_array12, &
   aligned_array12,capface_array1,sdim)

  ! consider all possible orderings of the 8 nodes relative to each other
  ! that does not change the "shape" of a regular hexahedron.
 do flipx=0,1
 do flipy=0,1
 do flipz=0,1

 do rotatez=0,3
 do rotatey=0,3
 do rotatex=0,3

  inode=1
  do k=1,2
  do j=1,2
  do i=1,2
   basegridmap(i,j,k)=inode
   inode=inode+1
  enddo
  enddo
  enddo

  do k=1,2
  do j=1,2
  do i=1,2
   if (flipx.eq.0) then
    ii=i
   else 
    ii=3-i
   endif
   if (flipy.eq.0) then
    jj=j
   else 
    jj=3-j
   endif
   if (flipz.eq.0) then
    kk=k
   else 
    kk=3-k
   endif
 
 
   rotmap(i,j,k,0)=basegridmap(ii,jj,kk)
  enddo
  enddo
  enddo

  ! rotate about z axis
  do k=1,2
   do nrot=1,rotatez
    rotmap(1,1,k,nrot)=rotmap(2,1,k,nrot-1)
    rotmap(2,1,k,nrot)=rotmap(2,2,k,nrot-1)
    rotmap(2,2,k,nrot)=rotmap(1,2,k,nrot-1)
    rotmap(1,2,k,nrot)=rotmap(1,1,k,nrot-1)
   enddo
  enddo

  do k=1,2
  do j=1,2
  do i=1,2
   rotmap(i,j,k,0)=rotmap(i,j,k,rotatez)
  enddo
  enddo
  enddo

  ! rotate about y axis
  do j=1,2
   do nrot=1,rotatey
    rotmap(1,j,1,nrot)=rotmap(2,j,1,nrot-1)
    rotmap(2,j,1,nrot)=rotmap(2,j,2,nrot-1)
    rotmap(2,j,2,nrot)=rotmap(1,j,2,nrot-1)
    rotmap(1,j,2,nrot)=rotmap(1,j,1,nrot-1)
   enddo
  enddo


  do k=1,2
  do j=1,2
  do i=1,2
   rotmap(i,j,k,0)=rotmap(i,j,k,rotatey)
  enddo
  enddo
  enddo

  ! rotate about x axis
  do i=1,2
   do nrot=1,rotatex
    rotmap(i,1,1,nrot)=rotmap(i,2,1,nrot-1)
    rotmap(i,2,1,nrot)=rotmap(i,2,2,nrot-1)
    rotmap(i,2,2,nrot)=rotmap(i,1,2,nrot-1)
    rotmap(i,1,2,nrot)=rotmap(i,1,1,nrot-1)
   enddo
  enddo

  do k=1,2
  do j=1,2
  do i=1,2
   basegridmap(i,j,k)=rotmap(i,j,k,rotatex)
  enddo
  enddo
  enddo

  do npos=1,9
   call add_to_hex(basegridmap,template_hex_plane(npos))
  enddo

 enddo
 enddo
 enddo  ! rotatex,y,z

 enddo
 enddo
 enddo  ! flipx,y,z

  ! now initialize the tetrahedron intersection table

 sdim=3

 n_nodes=4
 n_faces=3
 n_capfaces=1
 n_pos_nodes=1
 node_array4=(/ 11,12,13,14 /)
 face_array3=(/ 11,13,14, 11,12,14, 11,12,13 /)
 aligned_array3=(/ 1,1,1 /)  
 capface_array1=(/ 12,13,14 /)
 
 call init_intersect_type(template_tet_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array3, &
   aligned_array3,capface_array1,sdim)

 n_nodes=6
 n_faces=6
 n_capfaces=2
 n_pos_nodes=2
 node_array6=(/ 11,22,14,13,24,23 /)
 face_array6=(/ 11,24,22, 11,24,14, 11,23,22, 11,23,13, &
    11,13,14, 22,23,24 /)
 aligned_array6=(/ 1,1,1,1,1,0 /)  
 capface_array2=(/ 14,23,13, 14,23,24 /)

 call init_intersect_type(template_tet_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array6,face_array6, &
   aligned_array6,capface_array2,sdim)

 n_nodes=6
 n_faces=7
 n_capfaces=1
 n_pos_nodes=3
 node_array6=(/ 11,22,33,14,34,24 /)
 face_array7=(/ 11,22,33, 22,34,33, 22,34,24, &
   11,34,33, 11,34,14, 11,24,22, 11,24,14 /)
 aligned_array7=(/ 1,0,0,1,1,1,1 /)
 capface_array1=(/ 14,34,24 /)

 call init_intersect_type(template_tet_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array6,face_array7, &
   aligned_array7,capface_array1,sdim)

 n_nodes=4
 n_faces=4
 n_capfaces=0
 n_pos_nodes=4
 node_array4=(/ 11,22,33,44 /)
 face_array4=(/ 11,22,33, 11,22,44, 11,33,44, 22,33,44 /)
 aligned_array4=(/ 1,1,1,0 /)
 capface_array1=(/ 0,0,0 /)

 call init_intersect_type(template_tet_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array4, &
   aligned_array4,capface_array1,sdim)

  ! all possible orderings of the nodes with respect to each other
  ! does not change the "shape" of a tetrahedra.
 do i1=1,4
 do i2=1,4
 do i3=1,4
 do i4=1,4
  baselinemap(1)=i1
  if (i2.ne.i1) then
   baselinemap(2)=i2
   if ((i3.ne.i2).and.(i3.ne.i1)) then
    baselinemap(3)=i3
    if ((i4.ne.i3).and.(i4.ne.i2).and.(i4.ne.i1)) then
     baselinemap(4)=i4

     do npos=1,4
      call add_to_tet(baselinemap,template_tet_plane(npos))
     enddo
    endif
   endif
  endif
 enddo
 enddo
 enddo
 enddo  ! i1,i2,i3,i4
    
   ! now initialize the rectangle intersection table

 sdim=2

 n_nodes=3
 n_faces=2
 n_capfaces=1
 n_pos_nodes=1
 node_array3=(/ 11,13,12 /)
 face_array2=(/ 11,13,0, 11,12,0 /)
 aligned_array2=(/ 1,1 /)
 capface_array1=(/ 12,13,0 /)

 call init_intersect_type(template_rec_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array3,face_array2, &
   aligned_array2,capface_array1,sdim)

 n_nodes=4
 n_faces=3
 n_capfaces=1
 n_pos_nodes=2
 node_array4=(/ 11,22,13,24 /)
 face_array3=(/ 11,22,0, 11,13,0, 22,24,0 /)
 aligned_array3=(/ 1,1,0 /)
 capface_array1=(/ 13,24,0 /) 

 call init_intersect_type(template_rec_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array3, &
   aligned_array3,capface_array1,sdim)

 n_nodes=5
 n_faces=4
 n_capfaces=1
 n_pos_nodes=3
 node_array5=(/ 11,22,33,34,24 /)
 face_array4=(/ 11,22,0, 11,33,0, 33,34,0, 22,24,0 /)
 aligned_array4=(/ 1,1,0,0 /)
 capface_array1=(/ 34,24,0 /)  

 call init_intersect_type(template_rec_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array5,face_array4, &
   aligned_array4,capface_array1,sdim)

 n_nodes=4
 n_faces=4
 n_capfaces=0
 n_pos_nodes=4
 node_array4=(/ 11,22,33,44 /)
 face_array4=(/ 11,22,0, 22,44,0, 33,44,0, 11,33,0 /)
 aligned_array4=(/ 1,0,0,1 /)
 capface_array1=(/ 0,0,0 /)

 call init_intersect_type(template_rec_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array4, &
   aligned_array4,capface_array1,sdim)


  ! consider all possible orderings of the 4 nodes relative to each other
  ! that does not change the "shape" of a rectangle.
  !   3 4
  !   1 2
 do flipx=0,1
 do flipy=0,1

 do rotatez=0,3

  inode=1
  do j=1,2
  do i=1,2
   basegridmap2d(i,j)=inode
   inode=inode+1
  enddo
  enddo

  do j=1,2
  do i=1,2
   if (flipx.eq.0) then
    ii=i
   else 
    ii=3-i
   endif
   if (flipy.eq.0) then
    jj=j
   else 
    jj=3-j
   endif
 
   rotmap2d(i,j,0)=basegridmap2d(ii,jj)
  enddo
  enddo

  ! rotate about z axis
  do nrot=1,rotatez
    rotmap2d(1,1,nrot)=rotmap2d(2,1,nrot-1)
    rotmap2d(2,1,nrot)=rotmap2d(2,2,nrot-1)
    rotmap2d(2,2,nrot)=rotmap2d(1,2,nrot-1)
    rotmap2d(1,2,nrot)=rotmap2d(1,1,nrot-1)
  enddo

  do j=1,2
  do i=1,2
   basegridmap2d(i,j)=rotmap2d(i,j,rotatez)
  enddo
  enddo

  do npos=1,4
   call add_to_rec(basegridmap2d,template_rec_plane(npos))
  enddo

 enddo  ! rotatez

 enddo
 enddo  ! flipx,y


  ! now initialize the triangle intersection table

 sdim=2

 n_nodes=3
 n_faces=2
 n_capfaces=1
 n_pos_nodes=1
 node_array3=(/ 11,12,13 /)
 face_array2=(/ 11,13,0, 11,12,0 /)
 aligned_array2=(/ 1,1 /)  
 capface_array1=(/ 12,13,0 /)
 
 call init_intersect_type(template_tri_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array3,face_array2, &
   aligned_array2,capface_array1,sdim)

 n_nodes=4
 n_faces=3
 n_capfaces=1
 n_pos_nodes=2
 node_array4=(/ 11,22,13,23 /)
 face_array3=(/ 11,22,0, 11,13,0, 22,23,0 /)
 aligned_array3=(/ 1,1,0 /)  
 capface_array1=(/ 13,23,0 /)

 call init_intersect_type(template_tri_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array3, &
   aligned_array3,capface_array1,sdim)

 n_nodes=3
 n_faces=3
 n_capfaces=0
 n_pos_nodes=3
 node_array3=(/ 11,22,33 /)
 face_array3=(/ 11,22,0, 11,33,0, 22,33,0 /)
 aligned_array3=(/ 1,1,0 /)      
 capface_array1=(/ 0,0,0 /)

 call init_intersect_type(template_tri_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array3,face_array3, &
   aligned_array3,capface_array1,sdim)

  ! all possible orderings of the nodes with respect to each other
  ! does not change the "shape" of a triangle.
 do i1=1,3
 do i2=1,3
 do i3=1,3
  baselinemap(1)=i1
  if (i2.ne.i1) then
   baselinemap(2)=i2
   if ((i3.ne.i2).and.(i3.ne.i1)) then
    baselinemap(3)=i3

    do npos=1,3
     call add_to_tri(baselinemap,template_tri_plane(npos))
    enddo
   endif
  endif
 enddo
 enddo
 enddo  ! i1,i2,i3

return
end subroutine init_geometry_tables

subroutine create_xnodelist( &
  n_vol,n_area, &
  cum_volume,cum_area,cum_centroid, &
  xnode,phinode,checksum,maxnode, &
  shapeflag, & !0=regular hexahedron  1=tetrahedron
  nodedomain,sdim)

IMPLICIT NONE

integer, INTENT(in) :: sdim

integer, INTENT(inout) :: n_vol
integer, INTENT(inout) :: n_area
real(amrex_real), INTENT(inout) :: cum_volume
real(amrex_real), INTENT(inout) :: cum_area
real(amrex_real), INTENT(inout) :: cum_centroid(sdim)

integer, INTENT(in) :: shapeflag
integer, INTENT(in) :: nodedomain

integer, INTENT(in) :: checksum,maxnode
real(amrex_real), INTENT(in) :: xnode(nodedomain,sdim)
real(amrex_real), INTENT(in) :: phinode(nodedomain)
integer mapped_nodes(nodedomain)
integer n_nodes,dir,index1,index2,n_tet,n_capfaces
integer i
integer j_tet_node
real(amrex_real) x1(sdim),x2(sdim)
real(amrex_real) phi1,phi2
real(amrex_real) xtet(sdim+1,sdim)
real(amrex_real) xtri(sdim,sdim)
real(amrex_real) local_volume,local_area
real(amrex_real) local_centroid(sdim)
real(amrex_real) xnodelist_array(maxnodelist,sdim)
integer local_mapped_node
type(intersect_type) :: template_geom

 if (sdim.eq.3) then

  if (shapeflag.eq.0) then
   call fast_copy_intersect_type( &
    hexahedron_maps(checksum,maxnode)%intersect_geometry, &
    template_geom, &
    sdim)
   do i=1,nodedomain
    mapped_nodes(i)=hexahedron_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else if (shapeflag.eq.1) then
   call fast_copy_intersect_type( &
    tetrahedron_maps(checksum,maxnode)%intersect_geometry, &
    template_geom, &
    sdim)
   do i=1,nodedomain
    mapped_nodes(i)=tetrahedron_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else
   print *,"shapeflag invalid"
   stop
  endif

 else if (sdim.eq.2) then

  if (shapeflag.eq.0) then
   call fast_copy_intersect_type( &
    rectangle_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=rectangle_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else if (shapeflag.eq.1) then
   call fast_copy_intersect_type( &
    triangle_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=triangle_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else
   print *,"shapeflag invalid"
   stop
  endif

 else
  print *,"sdim invalid"
  stop
 endif

 n_nodes=template_geom%n_nodes
 if (n_nodes.le.0) then
  print *,"no support for this intersection shape"
  stop
 else if (n_nodes.gt.0) then
  do i=1,n_nodes
   index1=template_geom%nodelistmap(i,2) 
   index2=template_geom%nodelistmap(i,3) 
   if (index1.eq.index2) then
    if ((index1.lt.1).or.(index1.gt.nodedomain)) then
     print *,"index1 invalid"
     stop
    endif
    do dir=1,sdim
     xnodelist_array(i,dir)=xnode(mapped_nodes(index1),dir)
    enddo
   else if (index1.lt.index2) then
    if ((index1.lt.1).or.(index2.gt.nodedomain)) then
     print *,"index1 or index2 invalid"
     stop
    endif
    do dir=1,sdim
     x1(dir)=xnode(mapped_nodes(index1),dir)
     x2(dir)=xnode(mapped_nodes(index2),dir)
    enddo
    phi1=phinode(mapped_nodes(index1))
    phi2=phinode(mapped_nodes(index2))

    if (((phi1.ge.zero).and.(phi2.le.zero)).or. &
        ((phi1.le.zero).and.(phi2.ge.zero))) then
     do dir=1,sdim
      if (phi1.eq.zero) then
       xnodelist_array(i,dir)=x1(dir)
      else if (phi2.eq.zero) then
       xnodelist_array(i,dir)=x2(dir)
      else if (((phi1.gt.zero).and.(phi2.lt.zero)).or. &
               ((phi1.lt.zero).and.(phi2.gt.zero))) then
       xnodelist_array(i,dir)= &
        (abs(phi1)*x2(dir)+abs(phi2)*x1(dir))/ &
        (abs(phi1)+abs(phi2))
      else
       print *,"phi1 or phi2 bust"
       print *,"phi1,phi2 = ",phi1,phi2
       stop
      endif
     enddo ! dir=1..sdim
    else
     print *,"phi does not change sign"
     print *,"phi1,phi2 = ",phi1,phi2
     print *,"index1,index2 ",index1,index2
     print *,"mapped: ",mapped_nodes(index1),mapped_nodes(index2)
     print *,"x1 ",x1(1),x1(2),x1(sdim)
     print *,"x2 ",x2(1),x2(2),x2(sdim)
     print *,"i,n_nodes ",i,n_nodes
     print *,"(breakpoint) break point and gdb: "
     print *,"(1) compile with the -g option"
     print *,"(2) break MOF.F90:1339"
     stop
    endif

   else
    print *,"index1 should be < index2"
    stop
   endif  

   if (1.eq.0) then
    print *,"in create xnodelist checksum,maxnode,shapeflag,nodedomain ", &
      checksum,maxnode,shapeflag,nodedomain
    print *,"sdim,n_nodes ",sdim,n_nodes
    print *,"i,index1,index2 ",i,index1,index2 
    print *,"i,mapped(index1),mapped(index2) ",i, &
     mapped_nodes(index1),mapped_nodes(index2) 
    do dir=1,sdim
     print *,"i,dir,xnodelist_array ",i,dir,xnodelist_array(i,dir)
    enddo
   endif
  enddo ! initializing the nodes

  n_tet=template_geom%n_tet
  n_capfaces=template_geom%n_capfaces

  do i=1,n_tet
   n_vol=n_vol+1

   do j_tet_node=1,sdim+1
    do dir=1,sdim
     xtet(j_tet_node,dir)= &
        xnodelist_array(template_geom%tetlist(i,j_tet_node),dir)
    enddo
   enddo

   call tetrahedron_volume(xtet,local_volume,local_centroid,sdim)

   cum_volume=cum_volume+local_volume
   do dir=1,sdim
    cum_centroid(dir)=cum_centroid(dir)+local_volume*local_centroid(dir)
   enddo

   if (1.eq.0) then
    print *,"n_tet,n_vol,cum_volume ",n_tet,n_vol,cum_volume
    do j_tet_node=1,sdim+1
     print *,"i,j_tet_node,tetnode ",i,j_tet_node, &
             template_geom%tetlist(i,j_tet_node)
     do dir=1,sdim
      print *,"i,j_tet_node,dir,xtetnode ",i,j_tet_node,dir, &
       xnodelist_array(template_geom%tetlist(i,j_tet_node),dir)
     enddo
    enddo
   endif

  enddo ! looping through all tets 
   
  if ((n_capfaces.ge.0).and.(n_capfaces.le.maxcapfacelist)) then
  
   do i=1,n_capfaces
    n_area=n_area+1

    do j_tet_node=1,sdim

     local_mapped_node=template_geom%capfacelist(i,j_tet_node)

     if ((local_mapped_node.ge.1).and. &
         (local_mapped_node.le.maxnodelist)) then

      do dir=1,sdim
       xtri(j_tet_node,dir)=xnodelist_array(local_mapped_node,dir)
      enddo

     else
      print *,"local_mapped_node invalid: ",local_mapped_node
      print *,"maxnodelist=",maxnodelist
      stop
     endif

    enddo !j_tet_node=1..sdim

    call surface_area(xtri,local_area,sdim)

    cum_area=cum_area+local_area
   enddo ! i=1..n_capfaces

  else
   print *,"n_capfaces invalid: ",n_capfaces
   print *,"maxcapfacelist=",maxcapfacelist
   stop
  endif

 else
  print *,"n_nodes invalid: ",n_nodes
  stop
 endif 

return
end subroutine create_xnodelist


subroutine create_xnodelist_simple( &
  n_vol, &
  cum_volume,cum_centroid, &
  xnode,phinode,checksum,maxnode, &
  shapeflag, & !0=regular hexahedron 1=tetrahedron
  nodedomain,sdim)

IMPLICIT NONE

integer, INTENT(in) :: sdim

integer, INTENT(inout) :: n_vol
real(amrex_real), INTENT(inout) :: cum_volume
real(amrex_real), INTENT(inout) :: cum_centroid(sdim)

integer, INTENT(in) :: checksum,maxnode,shapeflag,nodedomain
real(amrex_real), INTENT(in) :: xnode(nodedomain,sdim)
real(amrex_real), INTENT(in) :: phinode(nodedomain)
integer mapped_nodes(nodedomain)
integer n_nodes,dir,index1,index2,n_tet
integer i
integer j_tet_node
real(amrex_real) x1(sdim),x2(sdim)
real(amrex_real) phi1,phi2
real(amrex_real) xtet(sdim+1,sdim)
real(amrex_real) local_volume
real(amrex_real) local_centroid(sdim)
real(amrex_real) xnodelist_array(maxnodelist,sdim)
type(intersect_type) :: template_geom

 if (sdim.eq.3) then

  if (shapeflag.eq.0) then
   call fast_copy_intersect_type( &
    hexahedron_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=hexahedron_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else if (shapeflag.eq.1) then
   call fast_copy_intersect_type( &
    tetrahedron_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=tetrahedron_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else
   print *,"shapeflag invalid"
   stop
  endif

 else if (sdim.eq.2) then

  if (shapeflag.eq.0) then
   call fast_copy_intersect_type( &
    rectangle_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=rectangle_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else if (shapeflag.eq.1) then
   call fast_copy_intersect_type( &
    triangle_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=triangle_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else
   print *,"shapeflag invalid"
   stop
  endif

 else
  print *,"sdim invalid"
  stop
 endif

 n_nodes=template_geom%n_nodes
 if (n_nodes.le.0) then
  print *,"no support for this intersection shape"
  stop
 else if (n_nodes.gt.0) then
  do i=1,n_nodes
   index1=template_geom%nodelistmap(i,2) 
   index2=template_geom%nodelistmap(i,3) 
   if (index1.eq.index2) then
    if ((index1.lt.1).or.(index1.gt.nodedomain)) then
     print *,"index1 invalid"
     stop
    endif
    do dir=1,sdim
     xnodelist_array(i,dir)=xnode(mapped_nodes(index1),dir)
    enddo
   else if (index1.lt.index2) then
    if ((index1.lt.1).or.(index2.gt.nodedomain)) then
     print *,"index1 or index2 invalid"
     stop
    endif
    do dir=1,sdim
     x1(dir)=xnode(mapped_nodes(index1),dir)
     x2(dir)=xnode(mapped_nodes(index2),dir)
    enddo
    phi1=phinode(mapped_nodes(index1))
    phi2=phinode(mapped_nodes(index2))
    if (((phi1.lt.zero).and.(phi2.lt.zero)).or. &
        ((phi1.ge.zero).and.(phi2.ge.zero))) then
     print *,"phi does not change sign"
     print *,"phi1,phi2 = ",phi1,phi2
     stop
    endif
    do dir=1,sdim
     xnodelist_array(i,dir)= &
      (abs(phi1)*x2(dir)+abs(phi2)*x1(dir))/ &
      (abs(phi1)+abs(phi2))
    enddo 
   else
    print *,"index1 should be < index2"
    stop
   endif  

   if (1.eq.0) then
    print *,"in create xnodelist checksum,maxnode,shapeflag,nodedomain ", &
      checksum,maxnode,shapeflag,nodedomain
    print *,"sdim,n_nodes ",sdim,n_nodes
    print *,"i,index1,index2 ",i,index1,index2 
    print *,"i,mapped(index1),mapped(index2) ",i, &
     mapped_nodes(index1),mapped_nodes(index2) 
    do dir=1,sdim
     print *,"i,dir,xnodelist_array ",i,dir,xnodelist_array(i,dir)
    enddo
   endif
  enddo ! initializing the nodes

  n_tet=template_geom%n_tet

  do i=1,n_tet
   n_vol=n_vol+1

   do j_tet_node=1,sdim+1
    do dir=1,sdim
     xtet(j_tet_node,dir)= &
        xnodelist_array(template_geom%tetlist(i,j_tet_node),dir)
    enddo
   enddo

   call tetrahedron_volume(xtet,local_volume,local_centroid,sdim)

   cum_volume=cum_volume+local_volume
   do dir=1,sdim
    cum_centroid(dir)=cum_centroid(dir)+local_volume*local_centroid(dir)
   enddo

   if (1.eq.0) then
    print *,"n_tet,n_vol,cum_volume ",n_tet,n_vol,cum_volume
    do j_tet_node=1,sdim+1
     print *,"i,j_tet_node,tetnode ",i,j_tet_node, &
             template_geom%tetlist(i,j_tet_node)
     do dir=1,sdim
      print *,"i,j_tet_node,dir,xtetnode ",i,j_tet_node,dir, &
       xnodelist_array(template_geom%tetlist(i,j_tet_node),dir)
     enddo
    enddo
   endif

  enddo ! looping through all tets 
 
 else
  print *,"n_nodes invalid"
  stop  
 endif 

return
end subroutine create_xnodelist_simple


subroutine create_xnodelist_and_map( &
  normdir,coeff, &
  n_vol, &
  cum_volume,cum_centroid, &
  cum_volume_map,cum_centroid_map, &
  xnode,phinode,checksum,maxnode, &
  shapeflag, & !0=regular hexahedron  1=tetrahedron
  nodedomain,sdim)

IMPLICIT NONE

integer, INTENT(in) :: sdim
integer, INTENT(in) :: normdir
real(amrex_real), INTENT(in) :: coeff(2)
integer, INTENT(inout) :: n_vol
real(amrex_real), INTENT(inout) :: cum_volume
real(amrex_real), INTENT(inout) :: cum_centroid(sdim)
real(amrex_real), INTENT(inout) :: cum_volume_map
real(amrex_real), INTENT(inout) :: cum_centroid_map(sdim)

integer, INTENT(in) :: checksum,maxnode,shapeflag,nodedomain
real(amrex_real), INTENT(in) :: xnode(nodedomain,sdim)
real(amrex_real), INTENT(in) :: phinode(nodedomain)
integer mapped_nodes(nodedomain)
integer n_nodes,dir,index1,index2,n_tet
integer i
integer j_tet_node
real(amrex_real) x1(sdim),x2(sdim)
real(amrex_real) phi1,phi2
real(amrex_real) xtet(sdim+1,sdim)
real(amrex_real) local_volume
real(amrex_real) local_centroid(sdim)
real(amrex_real) local_volume_map
real(amrex_real) local_centroid_map(sdim)
real(amrex_real) xnodelist_array(maxnodelist,sdim)
type(intersect_type) :: template_geom

 if ((normdir.lt.0).or.(normdir.ge.sdim)) then
  print *,"normdir invalid"
  stop
 endif
 if (coeff(1).le.zero) then
  print *,"coeff(1) invalid"
  stop
 endif

 if (sdim.eq.3) then

  if (shapeflag.eq.0) then
   call fast_copy_intersect_type( &
    hexahedron_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=hexahedron_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else if (shapeflag.eq.1) then
   call fast_copy_intersect_type( &
    tetrahedron_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=tetrahedron_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else
   print *,"shapeflag invalid"
   stop
  endif

 else if (sdim.eq.2) then

  if (shapeflag.eq.0) then
   call fast_copy_intersect_type( &
    rectangle_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=rectangle_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else if (shapeflag.eq.1) then
   call fast_copy_intersect_type( &
    triangle_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=triangle_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else
   print *,"shapeflag invalid"
   stop
  endif

 else
  print *,"sdim invalid"
  stop
 endif

 n_nodes=template_geom%n_nodes
 if (n_nodes.le.0) then
  print *,"no support for this intersection shape"
  stop
 else if (n_nodes.gt.0) then
  do i=1,n_nodes
   index1=template_geom%nodelistmap(i,2) 
   index2=template_geom%nodelistmap(i,3) 
   if (index1.eq.index2) then
    if ((index1.lt.1).or.(index1.gt.nodedomain)) then
     print *,"index1 invalid"
     stop
    endif
    do dir=1,sdim
     xnodelist_array(i,dir)=xnode(mapped_nodes(index1),dir)
    enddo
   else if (index1.lt.index2) then
    if ((index1.lt.1).or.(index2.gt.nodedomain)) then
     print *,"index1 or index2 invalid"
     stop
    endif
    do dir=1,sdim
     x1(dir)=xnode(mapped_nodes(index1),dir)
     x2(dir)=xnode(mapped_nodes(index2),dir)
    enddo
    phi1=phinode(mapped_nodes(index1))
    phi2=phinode(mapped_nodes(index2))
    if (((phi1.lt.zero).and.(phi2.lt.zero)).or. &
        ((phi1.ge.zero).and.(phi2.ge.zero))) then
     print *,"phi does not change sign"
     print *,"phi1,phi2 = ",phi1,phi2
     stop
    endif
    do dir=1,sdim
     xnodelist_array(i,dir)= &
      (abs(phi1)*x2(dir)+abs(phi2)*x1(dir))/ &
      (abs(phi1)+abs(phi2))
    enddo 
   else
    print *,"index1 should be < index2"
    stop
   endif  

   if (1.eq.0) then
    print *,"in create xnodelist checksum,maxnode,shapeflag,nodedomain ", &
      checksum,maxnode,shapeflag,nodedomain
    print *,"sdim,n_nodes ",sdim,n_nodes
    print *,"i,index1,index2 ",i,index1,index2 
    print *,"i,mapped(index1),mapped(index2) ",i, &
     mapped_nodes(index1),mapped_nodes(index2) 
    do dir=1,sdim
     print *,"i,dir,xnodelist_array ",i,dir,xnodelist_array(i,dir)
    enddo
   endif
  enddo ! initializing the nodes

  n_tet=template_geom%n_tet

  do i=1,n_tet
   n_vol=n_vol+1

   do j_tet_node=1,sdim+1
    do dir=1,sdim
     xtet(j_tet_node,dir)= &
       xnodelist_array(template_geom%tetlist(i,j_tet_node),dir)
    enddo
   enddo

   call tetrahedron_volume_and_map(normdir,coeff, &
     xtet, &
     local_volume,local_centroid, &
     local_volume_map,local_centroid_map, &
     sdim)

   cum_volume=cum_volume+local_volume
   cum_volume_map=cum_volume_map+local_volume_map
   do dir=1,sdim
    cum_centroid(dir)=cum_centroid(dir)+local_volume*local_centroid(dir)
    cum_centroid_map(dir)=cum_centroid_map(dir)+ &
      local_volume_map*local_centroid_map(dir)
   enddo

   if (1.eq.0) then
    print *,"n_tet,n_vol,cum_volume ",n_tet,n_vol,cum_volume
    do j_tet_node=1,sdim+1
     print *,"i,j_tet_node,tetnode ",i,j_tet_node, &
             template_geom%tetlist(i,j_tet_node)
     do dir=1,sdim
      print *,"i,j_tet_node,dir,xtetnode ",i,j_tet_node,dir, &
       xnodelist_array(template_geom%tetlist(i,j_tet_node),dir)
     enddo
    enddo
   endif

  enddo ! looping through all tets 
 
 else
  print *,"n_nodes invalid"
  stop  
 endif

return
end subroutine create_xnodelist_and_map

subroutine increment_volume( &
 n_vol,n_area, &
 cum_volume,cum_area,cum_centroid, &
 phinode,xnode,nodedomain,sdim)

IMPLICIT NONE

integer, INTENT(in) :: sdim

integer, INTENT(inout) :: n_vol,n_area
real(amrex_real), INTENT(inout) :: cum_volume,cum_area
real(amrex_real), INTENT(inout) :: cum_centroid(sdim)

integer, INTENT(in) :: nodedomain
real(amrex_real), INTENT(in) :: phinode(nodedomain)
real(amrex_real), INTENT(in) :: xnode(nodedomain,sdim)
integer maxchecksum,checksum,power2,n
real(amrex_real) phimax
integer maxnode,n_nodes
integer, PARAMETER :: shapeflag=1  !tetrahedron

 if (nodedomain.ne.sdim+1) then
  print *,"nodedomain invalid"
  stop
 endif

 if (sdim.eq.2) then
  maxchecksum=7
 else if (sdim.eq.3) then
  maxchecksum=15
 else
  print *,"sdim invalid"
  stop
 endif

 checksum=0 
 phimax=zero
 maxnode=0
 power2=1
 do n=1,nodedomain
  if (phinode(n).ge.zero) then
   if (phinode(n).gt.phimax) then
    maxnode=n
    phimax=phinode(n)
   endif
   checksum=checksum+power2
  endif
  power2=2*power2
 enddo

 if (phimax.le.zero) then   

  ! do nothing

 else 

   if (sdim.eq.3) then
    n_nodes=tetrahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
    if (n_nodes.le.0) then
     print *,"tetra maps incomplete"
     stop
    endif
   else if (sdim.eq.2) then
    n_nodes=triangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
    if (n_nodes.le.0) then
     print *,"triangle maps incomplete"
     stop
    endif
   else
    print *,"sdim invalid"
    stop
   endif

   call create_xnodelist( &
    n_vol,n_area, &
    cum_volume,cum_area,cum_centroid, &
    xnode,phinode,checksum,maxnode,shapeflag, &
    nodedomain,sdim)

 endif ! phimax>0

return
end subroutine increment_volume


subroutine increment_volume_simple( &
 n_vol, &
 cum_volume,cum_centroid, &
 phinode,xnode,nodedomain,sdim)

IMPLICIT NONE

integer, INTENT(in) :: sdim

integer, INTENT(inout) ::  n_vol
real(amrex_real), INTENT(inout) :: cum_volume
real(amrex_real), INTENT(inout) :: cum_centroid(sdim)

integer, INTENT(in) :: nodedomain
real(amrex_real), INTENT(in) :: phinode(nodedomain)
real(amrex_real), INTENT(in) :: xnode(nodedomain,sdim)
integer maxchecksum,checksum,power2,n
real(amrex_real) phimax
integer maxnode,n_nodes
integer, PARAMETER :: shapeflag=1 !tetrahedron

 if (nodedomain.ne.sdim+1) then
  print *,"nodedomain invalid"
  stop
 endif

 if (sdim.eq.2) then
  maxchecksum=7
 else if (sdim.eq.3) then
  maxchecksum=15
 else
  print *,"sdim invalid"
  stop
 endif

 checksum=0 
 phimax=zero
 maxnode=0
 power2=1
 do n=1,nodedomain
  if (phinode(n).ge.zero) then
   if (phinode(n).gt.phimax) then
    maxnode=n
    phimax=phinode(n)
   endif
   checksum=checksum+power2
  endif
  power2=2*power2
 enddo

 if (phimax.le.zero) then   

  ! do nothing

 else 

   if (sdim.eq.3) then
    n_nodes=tetrahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
    if (n_nodes.le.0) then
     print *,"tetra maps incomplete"
     stop
    endif
   else if (sdim.eq.2) then
    n_nodes=triangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
    if (n_nodes.le.0) then
     print *,"triangle maps incomplete"
     stop
    endif
   else
    print *,"sdim invalid"
    stop
   endif

   call create_xnodelist_simple( &
    n_vol, &
    cum_volume,cum_centroid, &
    xnode,phinode,checksum,maxnode,shapeflag, &
    nodedomain,sdim)

 endif ! phimax>0

return
end subroutine increment_volume_simple



subroutine increment_volume_and_map( &
 normdir,coeff, &
 n_vol, &
 cum_volume,cum_centroid, &
 cum_volume_map,cum_centroid_map, &
 phinode,xnode,nodedomain,sdim)

IMPLICIT NONE

integer, INTENT(in) :: sdim
integer, INTENT(in) :: normdir
real(amrex_real), INTENT(in) :: coeff(2)

integer, INTENT(inout) :: n_vol
real(amrex_real), INTENT(inout) :: cum_volume
real(amrex_real), INTENT(inout) :: cum_centroid(sdim)
real(amrex_real), INTENT(inout) :: cum_volume_map
real(amrex_real), INTENT(inout) :: cum_centroid_map(sdim)

integer, INTENT(in) :: nodedomain
real(amrex_real), INTENT(in) :: phinode(nodedomain)
real(amrex_real), INTENT(in) :: xnode(nodedomain,sdim)
integer maxchecksum,checksum,power2,n
real(amrex_real) phimax
integer maxnode,n_nodes
integer, PARAMETER :: shapeflag=1 !tetrahedron

 if ((normdir.lt.0).or.(normdir.ge.sdim)) then
  print *,"normdir invalid"
  stop
 endif
 if (coeff(1).le.zero) then
  print *,"coeff(1) invalid"
  stop
 endif
 if (nodedomain.ne.sdim+1) then
  print *,"nodedomain invalid"
  stop
 endif

 if (sdim.eq.2) then
  maxchecksum=7
 else if (sdim.eq.3) then
  maxchecksum=15
 else
  print *,"sdim invalid"
  stop
 endif

 checksum=0 
 phimax=zero
 maxnode=0
 power2=1
 do n=1,nodedomain
  if (phinode(n).ge.zero) then
   if (phinode(n).gt.phimax) then
    maxnode=n
    phimax=phinode(n)
   endif
   checksum=checksum+power2
  endif
  power2=2*power2
 enddo

 if (phimax.le.zero) then   

  ! do nothing

 else 

   if (sdim.eq.3) then
    n_nodes=tetrahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
    if (n_nodes.le.0) then
     print *,"tetra maps incomplete"
     stop
    endif
   else if (sdim.eq.2) then
    n_nodes=triangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
    if (n_nodes.le.0) then
     print *,"triangle maps incomplete"
     stop
    endif
   else
    print *,"sdim invalid"
    stop
   endif

   call create_xnodelist_and_map( &
    normdir,coeff, &
    n_vol, &
    cum_volume,cum_centroid, &
    cum_volume_map,cum_centroid_map, &
    xnode,phinode,checksum,maxnode,shapeflag, &
    nodedomain,sdim)

 endif ! phimax>0

return
end subroutine increment_volume_and_map



! linearcut=1 => input always intersection with plane
! linearcut=0 => input might be intersection with plane
! linearcut=-1 => always break up initial domain into tets/tris
subroutine intersection_volume( &
  cum_volume,cum_area,cum_centroid, &
  phinode,xnode,nodedomain, &
  sdim,fullelementfast,linearcut)

IMPLICIT NONE

integer, INTENT(in) :: sdim

integer :: n_vol,n_area
real(amrex_real), INTENT(out) :: cum_volume,cum_area
real(amrex_real), INTENT(out) :: cum_centroid(sdim)

integer, INTENT(in) :: nodedomain
integer, INTENT(in) :: fullelementfast
integer, INTENT(in) :: linearcut
integer, PARAMETER :: bfact=1
real(amrex_real), INTENT(in) :: phinode(nodedomain)
real(amrex_real), INTENT(in) :: xnode(nodedomain,sdim)
integer maxchecksum,power2
real(amrex_real) xsten_grid(-3:3,sdim)
real(amrex_real) dxgrid(sdim)
integer shapeflag
integer, PARAMETER :: symmetry_flag=0
integer ntetbox,dir
real(amrex_real) phimax
integer checksum,maxnode,n,n_nodes,id,sub_nodedomain
real(amrex_real) xx(sdim+1,sdim)
real(amrex_real) ls(sdim+1)
integer, PARAMETER :: nhalf=3

 call get_ntetbox(ntetbox,symmetry_flag,sdim)

 if (sdim.eq.2) then
  if (nodedomain.eq.4) then ! domain is a box
   maxchecksum=15
   shapeflag=0
  else if (nodedomain.eq.3) then  ! domain is a triangle
   maxchecksum=7
   shapeflag=1
  else
   print *,"nodedomain invalid"
   stop
  endif
 else if (sdim.eq.3) then
  if (nodedomain.eq.8) then  ! domain is a regular hexahedron
   maxchecksum=255
   shapeflag=0
  else if (nodedomain.eq.4) then  ! domain is a tetrahedron
   maxchecksum=15
   shapeflag=1
  else
   print *,"nodedomain invalid"
   stop
  endif
 else
  print *,"sdim invalid"
  stop
 endif

 checksum=0 
 phimax=zero
 maxnode=0
 power2=1
 do n=1,nodedomain
  if (phinode(n).ge.zero) then
   if (phinode(n).gt.phimax) then
    maxnode=n
    phimax=phinode(n)
   endif
   checksum=checksum+power2
  endif
  power2=2*power2
 enddo

 cum_volume=zero
 cum_area=zero
 do dir=1,sdim
   cum_centroid(dir)=zero
 enddo
 n_vol=0
 n_area=0
 
 if (phimax.le.zero) then   

  ! do nothing

 else 

  if ((checksum.eq.maxchecksum).and. &
      (fullelementfast.eq.1)) then
    if (shapeflag.eq.0) then
     do dir=1,sdim
      xsten_grid(-1,dir)=xnode(1,dir)
      xsten_grid(1,dir)=xnode(nodedomain,dir)
      dxgrid(dir)=xsten_grid(1,dir)-xsten_grid(-1,dir)
      xsten_grid(0,dir)=half*(xsten_grid(-1,dir)+xsten_grid(1,dir))
     enddo  ! dir
     call Box_volumeFAST(bfact,dxgrid,xsten_grid,nhalf, &
       cum_volume,cum_centroid,sdim)
    else if (shapeflag.eq.1) then
     call tetrahedron_volume(xnode,cum_volume,cum_centroid,sdim) 
    else
     print *,"shapeflag invalid"
     stop
    endif
  else
   if (sdim.eq.3) then
    if (shapeflag.eq.0) then
     n_nodes=hexahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
    else if (shapeflag.eq.1) then
     n_nodes=tetrahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
     if (n_nodes.le.0) then
      print *,"tetra maps incomplete"
      stop
     endif
    else
     print *,"shapeflag invalid"
     stop
    endif
   else if (sdim.eq.2) then
    if (shapeflag.eq.0) then
     n_nodes=rectangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
    else if (shapeflag.eq.1) then
     n_nodes=triangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
     if (n_nodes.le.0) then
      print *,"triangle maps incomplete"
      stop
     endif
    else
     print *,"shapeflag invalid" 
     stop
    endif
   else
    print *,"sdim invalid"
    stop
   endif

! linearcut=1 => input always intersection with plane
! linearcut=0 => input might be intersection with plane
! linearcut=-1 => always break up initial domain into tets/tris

   if ((n_nodes.eq.0).or.(linearcut.eq.-1)) then

    if ((linearcut.eq.1).and.(1.eq.0)) then
     print *,"WARNING: table entry missing"
     print *,"nodedomain=",nodedomain
     print *,"sdim=",sdim
     print *,"fullelementfast=",fullelementfast
     print *,"checksum,maxnode,phimax ",checksum,maxnode,phimax
     print *,"ntetbox,maxchecksum,shapeflag ",ntetbox,maxchecksum, &
       shapeflag
     do n=1,nodedomain 
      print *,"n,phinode ",n,phinode(n)
      do dir=1,sdim
       print *,"n,dir,xnode ",n,dir,xnode(n,dir)
      enddo
     enddo
    endif

    sub_nodedomain=sdim+1
    do id=1,ntetbox
     call extract_tet(xnode,phinode,xx,ls,id,symmetry_flag,sdim)
     call increment_volume( &
      n_vol,n_area, &
      cum_volume,cum_area,cum_centroid, &
      ls,xx,sub_nodedomain,sdim)
    enddo 

   else if (n_nodes.gt.0) then

    call create_xnodelist( &
     n_vol,n_area, &
     cum_volume,cum_area,cum_centroid, &
     xnode,phinode,checksum,maxnode,shapeflag, &
     nodedomain,sdim)

   else
    print *,"n_nodes invalid"
    stop
   endif 

   if (cum_volume.gt.zero) then
    do dir=1,sdim
     cum_centroid(dir)=cum_centroid(dir)/cum_volume
    enddo
   else
    do dir=1,sdim
     cum_centroid(dir)=zero
    enddo
   endif
  endif ! not a full element or a full element that is "tetrahedralized"

 endif ! phimax>0

return
end subroutine intersection_volume

! nodedomain is 4*(sdim-1)
subroutine intersection_volume_simple( &
  cum_volume,cum_centroid, &
  phinode,xnode,nodedomain, &
  sdim)

IMPLICIT NONE

integer, INTENT(in) :: sdim

integer :: n_vol
real(amrex_real), INTENT(inout) :: cum_volume
real(amrex_real), INTENT(inout) :: cum_centroid(sdim)

integer, INTENT(in) :: nodedomain
integer :: bfact
real(amrex_real), INTENT(in) :: phinode(nodedomain)
real(amrex_real), INTENT(in) :: xnode(nodedomain,sdim)
integer maxchecksum,power2
real(amrex_real) xsten_grid(-3:3,sdim)
real(amrex_real) dxgrid(sdim)
integer symmetry_flag,ntetbox,dir
real(amrex_real) phimax
integer checksum,maxnode,n,n_nodes,id,sub_nodedomain
real(amrex_real) xx(sdim+1,sdim)
real(amrex_real) ls(sdim+1)
integer nhalf
integer, PARAMETER :: shapeflag=0 !regular hexahedron

 nhalf=3

 bfact=1
 symmetry_flag=0
 call get_ntetbox(ntetbox,symmetry_flag,sdim)

 if (sdim.eq.2) then
  if (nodedomain.eq.4) then ! domain is a box
   maxchecksum=15
  else
   print *,"nodedomain invalid"
   stop
  endif
 else if (sdim.eq.3) then
  if (nodedomain.eq.8) then  ! domain is a regular hexahedron
   maxchecksum=255
  else
   print *,"nodedomain invalid"
   stop
  endif
 else
  print *,"sdim invalid"
  stop
 endif

 checksum=0 
 phimax=zero
 maxnode=0
 power2=1
 do n=1,nodedomain
  if (phinode(n).ge.zero) then
   if (phinode(n).gt.phimax) then
    maxnode=n
    phimax=phinode(n)
   endif
   checksum=checksum+power2
  endif
  power2=2*power2
 enddo

 cum_volume=zero
 do dir=1,sdim
   cum_centroid(dir)=zero
 enddo
 n_vol=0
 
 if (phimax.le.zero) then   

  ! do nothing

 else 

  if (checksum.eq.maxchecksum) then
    do dir=1,sdim
     xsten_grid(-1,dir)=xnode(1,dir)
     xsten_grid(1,dir)=xnode(nodedomain,dir)
     dxgrid(dir)=xsten_grid(1,dir)-xsten_grid(-1,dir)
     xsten_grid(0,dir)=half*(xsten_grid(-1,dir)+xsten_grid(1,dir))
    enddo  ! dir
    call Box_volumeFAST(bfact,dxgrid,xsten_grid,nhalf, &
      cum_volume,cum_centroid,sdim)
  else
   if (sdim.eq.3) then
    n_nodes=hexahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
   else if (sdim.eq.2) then
    n_nodes=rectangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
   else
    print *,"sdim invalid"
    stop
   endif

   if (n_nodes.eq.0) then

    if (1.eq.0) then
     print *,"WARNING: table entry missing"
     print *,"nodedomain=",nodedomain
     print *,"sdim=",sdim
     print *,"checksum,maxnode,phimax ",checksum,maxnode,phimax
     print *,"ntetbox,maxchecksum,shapeflag ",ntetbox,maxchecksum, &
       shapeflag
     do n=1,nodedomain 
      print *,"n,phinode ",n,phinode(n)
      do dir=1,sdim
       print *,"n,dir,xnode ",n,dir,xnode(n,dir)
      enddo
     enddo
    endif

    sub_nodedomain=sdim+1
    do id=1,ntetbox
     call extract_tet(xnode,phinode,xx,ls,id,symmetry_flag,sdim)
     call increment_volume_simple( &
      n_vol, &
      cum_volume,cum_centroid, &
      ls,xx,sub_nodedomain,sdim)
    enddo 

   else if (n_nodes.gt.0) then

    call create_xnodelist_simple( &
     n_vol, &
     cum_volume,cum_centroid, &
     xnode,phinode,checksum,maxnode,shapeflag, &
     nodedomain,sdim)

   else
    print *,"n_nodes invalid"
    stop
   endif 

   if (cum_volume.gt.zero) then
    do dir=1,sdim
     cum_centroid(dir)=cum_centroid(dir)/cum_volume
    enddo
   else
    do dir=1,sdim
     cum_centroid(dir)=zero
    enddo
   endif

  endif ! not a full element or a full element that is "tetrahedralized"

 endif ! phimax>0

return
end subroutine intersection_volume_simple


! nodedomain is 4*(sdim-1)
subroutine intersection_volume_and_map( &
  normdir, &
  coeff, &
  cum_volume,cum_centroid, &
  cum_volume_map,cum_centroid_map, &
  phinode,xnode,nodedomain, &
  sdim)

IMPLICIT NONE

integer, INTENT(in) :: sdim
integer, INTENT(in) :: normdir
real(amrex_real), INTENT(in) :: coeff(2)
integer :: n_vol
real(amrex_real), INTENT(inout) :: cum_volume
real(amrex_real), INTENT(inout) :: cum_centroid(sdim)
real(amrex_real), INTENT(inout) :: cum_volume_map
real(amrex_real), INTENT(inout) :: cum_centroid_map(sdim)

integer, INTENT(in) :: nodedomain
integer :: bfact
real(amrex_real), INTENT(in) :: phinode(nodedomain)
real(amrex_real), INTENT(in) :: xnode(nodedomain,sdim)
integer maxchecksum,power2
real(amrex_real) xsten_grid(-3:3,sdim)
real(amrex_real) dxgrid(sdim)
integer symmetry_flag,ntetbox,dir
real(amrex_real) phimax
integer checksum,maxnode,n,n_nodes,id,sub_nodedomain
real(amrex_real) xx(sdim+1,sdim)
real(amrex_real) ls(sdim+1)
integer nhalf
integer, PARAMETER :: shapeflag=0 !regular hexahedron

 if ((normdir.lt.0).or.(normdir.ge.sdim)) then
  print *,"normdir invalid"
  stop
 endif
 if (coeff(1).le.zero) then
  print *,"coeff(1) invalid"
  stop
 endif

 nhalf=3

 bfact=1
 symmetry_flag=0
 call get_ntetbox(ntetbox,symmetry_flag,sdim)

 if (sdim.eq.2) then
  if (nodedomain.eq.4) then ! domain is a box
   maxchecksum=15
  else
   print *,"nodedomain invalid"
   stop
  endif
 else if (sdim.eq.3) then
  if (nodedomain.eq.8) then  ! domain is a regular hexahedron
   maxchecksum=255
  else
   print *,"nodedomain invalid"
   stop
  endif
 else
  print *,"sdim invalid"
  stop
 endif

 checksum=0 
 phimax=zero
 maxnode=0
 power2=1
 do n=1,nodedomain
  if (phinode(n).ge.zero) then
   if (phinode(n).gt.phimax) then
    maxnode=n
    phimax=phinode(n)
   endif
   checksum=checksum+power2
  endif
  power2=2*power2
 enddo

 cum_volume=zero
 cum_volume_map=zero
 do dir=1,sdim
   cum_centroid(dir)=zero
   cum_centroid_map(dir)=zero
 enddo
 n_vol=0
 
 if (phimax.le.zero) then   

  ! do nothing

 else 

  if (checksum.eq.maxchecksum) then
    do dir=1,sdim
     xsten_grid(-1,dir)=xnode(1,dir)
     xsten_grid(1,dir)=xnode(nodedomain,dir)
     dxgrid(dir)=xsten_grid(1,dir)-xsten_grid(-1,dir)
     xsten_grid(0,dir)=half*(xsten_grid(-1,dir)+xsten_grid(1,dir))
    enddo  ! dir
    call Box_volumeFAST_and_map( &
      normdir, &
      coeff, &
      bfact, &
      dxgrid, &
      xsten_grid, &
      nhalf, &
      cum_volume, &
      cum_centroid, &
      cum_volume_map, &
      cum_centroid_map, &
      sdim)
  else
   if (sdim.eq.3) then
    n_nodes=hexahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
   else if (sdim.eq.2) then
    n_nodes=rectangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
   else
    print *,"sdim invalid"
    stop
   endif

   if (n_nodes.eq.0) then

    if (1.eq.0) then
     print *,"WARNING: table entry missing"
     print *,"nodedomain=",nodedomain
     print *,"sdim=",sdim
     print *,"checksum,maxnode,phimax ",checksum,maxnode,phimax
     print *,"ntetbox,maxchecksum ",ntetbox,maxchecksum
     do n=1,nodedomain 
      print *,"n,phinode ",n,phinode(n)
      do dir=1,sdim
       print *,"n,dir,xnode ",n,dir,xnode(n,dir)
      enddo
     enddo
    endif

    sub_nodedomain=sdim+1
    do id=1,ntetbox
     call extract_tet(xnode,phinode,xx,ls,id,symmetry_flag,sdim)
     call increment_volume_and_map( &
      normdir,coeff, &
      n_vol, &
      cum_volume,cum_centroid, &
      cum_volume_map,cum_centroid_map, &
      ls,xx,sub_nodedomain,sdim)
    enddo 

   else if (n_nodes.gt.0) then

    call create_xnodelist_and_map( &
     normdir,coeff, &
     n_vol, &
     cum_volume,cum_centroid, &
     cum_volume_map,cum_centroid_map, &
     xnode,phinode,checksum,maxnode,shapeflag, &
     nodedomain,sdim)

   else
    print *,"n_nodes invalid"
    stop
   endif 

   if ((cum_volume.gt.zero).and.(cum_volume_map.gt.zero)) then
    do dir=1,sdim
     cum_centroid(dir)=cum_centroid(dir)/cum_volume
     cum_centroid_map(dir)=cum_centroid_map(dir)/cum_volume_map
    enddo
   else if ((cum_volume.eq.zero).or.(cum_volume_map.eq.zero)) then
    do dir=1,sdim
     cum_centroid(dir)=zero
     cum_centroid_map(dir)=zero
    enddo
   else
    print *,"cum_volume or cum_volume_map invalid"
    stop
   endif

  endif ! not a full element or a full element that is "tetrahedralized"

 endif ! phimax>0

return
end subroutine intersection_volume_and_map



      subroutine angle_to_slope(angle,nslope,sdim)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(out) :: nslope(sdim)
      real(amrex_real), INTENT(in) :: angle(sdim-1)
      integer :: dir

      if (MOF_PI.eq.zero) then
       MOF_PI=four*atan(one)
      endif

      do dir=1,sdim-1
       if ((angle(dir).ge.-four*MOF_PI).and. &
           (angle(dir).le.four*MOF_PI)) then
        !do nothing
       else
        print *,"angle(dir) invalid: ",dir,angle(dir)
        stop
       endif
      enddo !dir=1..sdim-1

      if (sdim.eq.3) then
       nslope(sdim)=cos(angle(sdim-1))
       nslope(1)=sin(angle(sdim-1))*cos(angle(1))
       nslope(2)=sin(angle(sdim-1))*sin(angle(1))
       if (abs(one- &
               sqrt(nslope(1)**2+nslope(2)**2+nslope(sdim)**2)).le.EPS2) then
        !do nothing
       else
        print *,"nslope not unit mag"
        do dir=1,sdim
         print *,"dir,nslope ",dir,nslope(dir)
        enddo
        stop
       endif

      else if (sdim.eq.2) then
       nslope(1)=cos(angle(1))
       nslope(2)=sin(angle(1))
       if (abs(one- &
               sqrt(nslope(1)**2+nslope(2)**2)).le.EPS2) then
        !do nothing
       else
        print *,"nslope not unit mag"
        do dir=1,sdim
         print *,"dir,nslope ",dir,nslope(dir)
        enddo
        stop
       endif
      else
       print *,"sdim invalid in angle_to_slope"
       stop
      endif

      return
      end subroutine angle_to_slope


        ! returns centroid in absolute coordinate system

      subroutine cell_intersection_grid( &
        bfact,dxgrid,xgrid,nhalf, &
        lnode, &
        volumedark,centroiddark, &
        area,volall,cenall,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: bfact

      integer, INTENT(in) :: sdim,nhalf
      real(amrex_real), INTENT(in) :: lnode(4*(sdim-1))
      real(amrex_real) xnode(4*(sdim-1),sdim)
      real(amrex_real), INTENT(in) :: xgrid(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dxgrid(sdim)
      real(amrex_real), INTENT(out) :: volumedark
      real(amrex_real), INTENT(out) :: centroiddark(sdim)
      real(amrex_real), INTENT(out) :: volall,area
      real(amrex_real), INTENT(out) :: cenall(sdim)
      real(amrex_real) centroididdark(sdim)
      real(amrex_real) xx(sdim+1,sdim)
      real(amrex_real) lsdark(sdim+1)
      real(amrex_real) volumeiddark
      real(amrex_real) areaiddark
      integer j_dir
      integer i_grid_node
      integer j_grid_node
      integer k_grid_node
      integer id,symmetry_flag,ntetbox
      integer power2,inode,nplus,nminus,ndark,sumdark

      if (nhalf.lt.1) then
       print *,"nhalf invalid cell_intersection_grid"
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust cell intersection grid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid120"
       stop
      endif
 
      do j_dir=1,sdim
       centroiddark(j_dir)=zero
      enddo
      volumedark=zero
      area=zero

      inode=1
      power2=1
      nplus=0
      nminus=0
      sumdark=0
      ndark=0

      call Box_volumeFAST(bfact,dxgrid,xgrid,nhalf,volall,cenall,sdim)

      if (volall.gt.zero) then
    
       symmetry_flag=0
       call get_ntetbox(ntetbox,symmetry_flag,sdim)
 
       if (sdim.eq.3) then
        do k_grid_node=-1,1,2
        do j_grid_node=-1,1,2
        do i_grid_node=-1,1,2 
         xnode(inode,1)=xgrid(i_grid_node,1)
         xnode(inode,2)=xgrid(j_grid_node,2)
         xnode(inode,sdim)=xgrid(k_grid_node,sdim)

         if (lnode(inode).ge.zero) then
          nplus=nplus+1
          sumdark=sumdark+power2
          ndark=ndark+1
         endif
         if (lnode(inode).le.zero) then
          nminus=nminus+1
         endif

         inode=inode+1
         power2=power2*2
        enddo
        enddo
        enddo
       else if (sdim.eq.2) then
        do j_grid_node=-1,1,2
        do i_grid_node=-1,1,2 
         xnode(inode,1)=xgrid(i_grid_node,1)
         xnode(inode,2)=xgrid(j_grid_node,2)

         if (lnode(inode).ge.zero) then
          nplus=nplus+1
          sumdark=sumdark+power2
          ndark=ndark+1
         endif
         if (lnode(inode).le.zero) then
          nminus=nminus+1
         endif

         inode=inode+1
         power2=power2*2
        enddo
        enddo
       else 
        print *,"sdim invalid"
        stop
       endif

       do id=1,ntetbox

        call extract_tet(xnode,lnode,xx,lsdark,id,symmetry_flag,sdim)

        if (sdim.eq.3) then
         call intersection_volumeXYZ(lsdark,xx,volumeiddark, &
           centroididdark,areaiddark,sdim)
        else if (sdim.eq.2) then
         call int_volumeXYorRZ(lsdark,xx,volumeiddark, &
           centroididdark,areaiddark,sdim)
        else
         print *,"sdim invalid"
         stop
        endif

        volumedark=volumedark+volumeiddark
        area=area+areaiddark
        do j_dir=1,sdim
         centroiddark(j_dir)=centroiddark(j_dir)+ &
                 centroididdark(j_dir)*volumeiddark
        enddo

       enddo ! id
      else if (volall.eq.zero) then
       print *,"WARNING:volall cannot be zero in cell_intersection_grid"
      else if (volall.lt.zero) then
       print *,"volall invalid"
       stop
      endif

      if (volumedark.gt.zero) then
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)/volumedark
       enddo
      else
       volumedark=zero
       do j_dir=1,sdim
        centroiddark(j_dir)=zero
       enddo
      endif

      return
      end subroutine cell_intersection_grid


! finds the volume of intersection between the region phi>=0 and
! a column with triangular cross section.
! Also, obtain the centroid and surface area of intersection.

      subroutine int_volumeXYorRZ(phi,x,volumedark, &
       centroiddark,area,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim

      real(amrex_real) xtrilist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(sdim,sdim,MAXAREA)
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real) xint(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: volumedark,area
      real(amrex_real) volumelistdark,arealist
      real(amrex_real), INTENT(out) :: centroiddark(sdim)
      real(amrex_real) centroidlistdark(sdim)

      integer i_tet_node
      integer n,nlist,narea
      integer j_dir

      if (sdim.ne.2) then
       print *,"sdim invalid"
       stop
      endif

      call list_tris(phi,x,xtrilist,MAXTET,nlist,xarealist,narea,sdim)

      volumedark=zero
      area=zero
      do j_dir=1,sdim
       centroiddark(j_dir)=zero
      enddo
      do n=1,nlist
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=xtrilist(i_tet_node,j_dir,n)
       enddo
       enddo
       call tetrahedron_volume(xint,volumelistdark,centroidlistdark,sdim)
       volumedark=volumedark+volumelistdark
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)+ &
                centroidlistdark(j_dir)*volumelistdark
       enddo
      enddo
      if (volumedark.gt.zero) then
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)/volumedark
       enddo
      else
       do j_dir=1,sdim
        centroiddark(j_dir)=zero
       enddo
      endif

      do n=1,narea
       do i_tet_node=1,sdim
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=xarealist(i_tet_node,j_dir,n)
       enddo
       enddo
       do j_dir=1,sdim
        xint(sdim+1,j_dir)=zero
       enddo
       call areaXYorRZ(xint,1,2,arealist)
       area=area+arealist
      enddo

      return
      end subroutine int_volumeXYorRZ


! finds the volume of intersection between the region phi>=0 and
! a column with triangular cross section.
! Also, obtain the centroid of intersection.

      subroutine int_volumeXYorRZ_simple(phi,x,volumedark, &
       centroiddark,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim

      real(amrex_real) xtrilist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(sdim,sdim,MAXAREA)
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real) xint(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: volumedark
      real(amrex_real) volumelistdark
      real(amrex_real), INTENT(out) :: centroiddark(sdim)
      real(amrex_real) centroidlistdark(sdim)

      integer i_tet_node
      integer j_dir
      integer n,nlist,narea

      if (sdim.ne.2) then
       print *,"sdim invalid"
       stop
      endif

      call list_tris(phi,x,xtrilist,MAXTET,nlist,xarealist,narea,sdim)

      volumedark=zero
      do j_dir=1,sdim
       centroiddark(j_dir)=zero
      enddo
      do n=1,nlist
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=xtrilist(i_tet_node,j_dir,n)
       enddo
       enddo
       call tetrahedron_volume(xint,volumelistdark,centroidlistdark,sdim)
       volumedark=volumedark+volumelistdark
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)+ &
                centroidlistdark(j_dir)*volumelistdark
       enddo
      enddo
      if (volumedark.gt.zero) then
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)/volumedark
       enddo
      else
       do j_dir=1,sdim
        centroiddark(j_dir)=zero
       enddo
      endif

      return
      end subroutine int_volumeXYorRZ_simple


! finds the volume of intersection between the region phi>=0 and
! a column with triangular cross section.
! Also, obtain the centroid of intersection.

      subroutine int_volumeXYorRZ_and_map( &
       normdir,coeff, &
       phi,x, &
       volumedark, &
       volumedark_map, &
       centroiddark, &
       centroiddark_map, &
       sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: normdir
      real(amrex_real), INTENT(in) :: coeff(2)

      real(amrex_real) xtrilist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(sdim,sdim,MAXAREA)
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real) xint(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: volumedark
      real(amrex_real) volumelistdark
      real(amrex_real), INTENT(out) :: volumedark_map
      real(amrex_real) volumelistdark_map
      real(amrex_real), INTENT(out) :: centroiddark(sdim)
      real(amrex_real) centroidlistdark(sdim)
      real(amrex_real), INTENT(out) :: centroiddark_map(sdim)
      real(amrex_real) centroidlistdark_map(sdim)

      integer i_tet_node
      integer j_dir
      integer n,nlist,narea

      if (sdim.ne.2) then
       print *,"sdim invalid"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.sdim)) then
       print *,"normdir invalid"
       stop
      endif
      if (coeff(1).le.zero) then
       print *,"coeff(1) invalid"
       stop
      endif

      call list_tris(phi,x,xtrilist,MAXTET,nlist,xarealist,narea,sdim)

      volumedark=zero
      volumedark_map=zero
      do j_dir=1,sdim
       centroiddark(j_dir)=zero
       centroiddark_map(j_dir)=zero
      enddo
      do n=1,nlist
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=xtrilist(i_tet_node,j_dir,n)
       enddo
       enddo
       call tetrahedron_volume_and_map(normdir,coeff, &
        xint,volumelistdark,centroidlistdark, &
        volumelistdark_map,centroidlistdark_map, &
        sdim)
       volumedark=volumedark+volumelistdark
       volumedark_map=volumedark_map+volumelistdark_map
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)+ &
                centroidlistdark(j_dir)*volumelistdark
        centroiddark_map(j_dir)=centroiddark_map(j_dir)+ &
         centroidlistdark_map(j_dir)*volumelistdark_map
       enddo
      enddo
      if ((volumedark.gt.zero).and.(volumedark_map.gt.zero)) then
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)/volumedark
        centroiddark_map(j_dir)=centroiddark_map(j_dir)/volumedark_map
       enddo
      else if ((volumedark.eq.zero).or.(volumedark_map.eq.zero)) then
       do j_dir=1,sdim
        centroiddark(j_dir)=zero
        centroiddark_map(j_dir)=zero
       enddo
      else
       print *,"volumedark or volumedark_map invalid"
       stop
      endif

      return
      end subroutine int_volumeXYorRZ_and_map




! find triangles that make up the intersection of a plane (phi>0) with
! a triangle
      subroutine list_tris(phi,x,xtrilist, &
                    nlist_alloc,nlist,xarealist,narea,sdim)
      IMPLICIT NONE 
  
      integer, INTENT(in) :: nlist_alloc 
      integer, INTENT(in) :: sdim 
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: xtrilist(sdim+1,sdim,nlist_alloc)
      real(amrex_real), INTENT(out) :: xarealist(2,2,MAXAREA)
      integer, INTENT(out) :: nlist
      integer, INTENT(out) :: narea
      integer i_tet_node
      integer j_dir

      if (sdim.ne.2) then
       print *,"sdim invalid list_tris"
       stop
      endif
      if (nlist_alloc.ge.1) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      nlist=0
      narea=0
      if ((phi(1).le.zero).and.(phi(2).le.zero).and. &
          (phi(3).le.zero)) then
       nlist=0
      else if ((phi(1).ge.zero).and.(phi(2).ge.zero).and. &
               (phi(3).ge.zero)) then
       nlist=1
       do i_tet_node=1,3
       do j_dir=1,2
        xtrilist(i_tet_node,j_dir,1)=x(i_tet_node,j_dir)
       enddo 
       enddo 
      else if ((phi(1).lt.zero).and.(phi(2).ge.zero).and. &
               (phi(3).ge.zero)) then
       call shrink_list_tri(phi,x,1,2,3,1,xtrilist, &
               nlist_alloc,nlist,xarealist,narea)
      else if ((phi(2).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(3).ge.zero)) then
       call shrink_list_tri(phi,x,2,1,3,1,xtrilist, &
               nlist_alloc,nlist,xarealist,narea)
      else if ((phi(3).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(2).ge.zero)) then
       call shrink_list_tri(phi,x,3,1,2,1,xtrilist, &
               nlist_alloc,nlist,xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(2).lt.zero).and. &
               (phi(3).lt.zero)) then
       call shrink_list_tri(phi,x,1,2,3,0,xtrilist, &
               nlist_alloc,nlist,xarealist,narea)
      else if ((phi(2).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(3).lt.zero)) then
       call shrink_list_tri(phi,x,2,1,3,0,xtrilist, &
               nlist_alloc,nlist,xarealist,narea)
      else if ((phi(3).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(2).lt.zero)) then
       call shrink_list_tri(phi,x,3,1,2,0,xtrilist, &
               nlist_alloc,nlist,xarealist,narea)
      else
       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break MOF.F90:3767"
       print *,"bust list_tris"
       print *,"phi : ",phi(1),phi(2),phi(3)
       print *,"sdim= ",sdim
       print *,"nlist_alloc= ",nlist_alloc
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        print *,"i_tet_node,j_dir,x ",i_tet_node,j_dir,x(i_tet_node,j_dir)
       enddo 
       enddo 
       stop
      endif

      return
      end subroutine list_tris



      subroutine get_xbounds(x,xmin,xmax,dxmax,npoints,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: npoints,sdim
      real(amrex_real), INTENT(in) :: x(npoints,sdim)
      real(amrex_real), INTENT(out) :: xmin(sdim)
      real(amrex_real), INTENT(out) :: xmax(sdim)
      real(amrex_real), INTENT(out) :: dxmax
      real(amrex_real) :: dxtrial
      integer i,dir,dir2

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif

      dxmax=zero

      do dir=1,sdim
       xmin(dir)=x(1,dir)
       xmax(dir)=x(1,dir)
      enddo
      do dir=1,sdim
       do i=1,npoints
        if (x(i,dir).lt.xmin(dir)) then
         xmin(dir)=x(i,dir)
        endif
        if (x(i,dir).gt.xmax(dir)) then
         xmax(dir)=x(i,dir)
        endif
       enddo !i=1,npoints
       dxtrial=xmax(dir)-xmin(dir)
       if (dxtrial.ge.zero) then
        ! do nothing
       else
        print *,"(breakpoint) break point and gdb: "
        print *,"(1) compile with the -g option"
        print *,"(2) break MOF.F90:3239"
        print *,"dxtrial is corrupt: ",dxtrial
        print *,"sdim ",sdim
        print *,"npoints ",npoints
        do i=1,npoints
        do dir2=1,sdim
         print *,"i,dir2,x ",i,dir2,x(i,dir)
        enddo
        enddo
        print *,"xmin ",xmin
        print *,"xmax ",xmax
        stop
       endif
       if (dxtrial.gt.dxmax) then
        dxmax=dxtrial
       else if (dxtrial.le.dxmax) then
        ! do nothing
       else
        print *,"dxtrial invalid"
        stop
       endif
      enddo !dir=1,sdim

      return
      end subroutine get_xbounds

! find the length of a side of a triangle      
      subroutine areaXYorRZ(x,i1,i2,area)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x(3,2)
      real(amrex_real) xx(2,2)
      real(amrex_real), INTENT(out) :: area
      integer, INTENT(in) :: i1,i2
      integer dir,sdim

      sdim=2

      if ((i1.lt.1).or.(i1.gt.3).or.(i2.lt.1).or.(i2.gt.3).or. &
          (i1.eq.i2)) then
       print *,"index invalid area xy or rz"
       print *,"i1,i2 ",i1,i2
       stop
      endif

      do dir=1,sdim
       xx(1,dir)=x(i1,dir) 
       xx(2,dir)=x(i2,dir) 
      enddo
      call surface_area(xx,area,sdim)

      return
      end subroutine areaXYorRZ

       ! in=1  (1,1,1) node
       ! in=2  (2,1,1) node
       ! in=3  (1,2,1) node
       ! in=4  (2,2,1) node
       ! in=5  (1,1,2) node
       ! in=6  (2,1,2) node
       ! in=7  (1,2,2) node
       ! in=8  (2,2,2) node
       ! in=1
       ! do k=1,2
       ! do j=1,2
       ! do i=1,2      
       !  in=in+1
       !  ....
      subroutine extract_tet(xnode,datanode,xtet,datatet, &
       id,symmetry_flag,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(out) :: xtet(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: datatet(sdim+1)
      integer, INTENT(in) :: id
      integer inode,jnode
      integer j_tet_node
      integer k_dir
      integer, INTENT(in) :: symmetry_flag

      real(amrex_real), INTENT(in) :: datanode(4*(sdim-1))
      real(amrex_real), INTENT(in) :: xnode(4*(sdim-1),sdim)
      integer nodelist(sdim+1)
      real(amrex_real) xavg,davg

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust extract tet"
       stop
      endif

      if (symmetry_flag.eq.0) then

       if (sdim.eq.3) then

        if (id.eq.1) then
         nodelist(1)=3 
         nodelist(2)=1 
         nodelist(3)=2 
         nodelist(4)=5
        else if (id.eq.2) then
         nodelist(1)=5 
         nodelist(2)=6 
         nodelist(3)=2 
         nodelist(4)=8
        else if (id.eq.3) then
         nodelist(1)=3 
         nodelist(2)=5 
         nodelist(3)=7 
         nodelist(4)=8
        else if (id.eq.4) then
         nodelist(1)=3 
         nodelist(2)=4 
         nodelist(3)=2 
         nodelist(4)=8
        else if (id.eq.5) then
         nodelist(1)=5 
         nodelist(2)=8 
         nodelist(3)=3 
         nodelist(4)=2
        else
         print *,"id invalid extract_tet"
         stop
        endif

       else if (sdim.eq.2) then

        if (id.eq.1) then
         nodelist(1)=1 
         nodelist(2)=2 
         nodelist(3)=3
        else if (id.eq.2) then 
         nodelist(1)=2 
         nodelist(2)=3 
         nodelist(3)=4
        else
         print *,"id invalid"
         stop
        endif

       else
        print *,"sdim invalid"
        stop
       endif

       do j_tet_node=1,sdim+1
        inode=nodelist(j_tet_node)
        do k_dir=1,sdim
         xtet(j_tet_node,k_dir)=xnode(inode,k_dir)
        enddo
        datatet(j_tet_node)=datanode(inode)
       enddo 

      else if (symmetry_flag.eq.1) then

        !  3 4
        !  1 2
       if (sdim.eq.2) then

        if (id.eq.1) then  ! dir=2 side=1
         nodelist(1)=1
         nodelist(2)=2
         nodelist(3)=0
        else if (id.eq.2) then ! dir=1 side=1
         nodelist(1)=1
         nodelist(2)=3
         nodelist(3)=0
        else if (id.eq.3) then ! dir=2 side=2
         nodelist(1)=3
         nodelist(2)=4
         nodelist(3)=0
        else if (id.eq.4) then ! dir=1 side=2
         nodelist(1)=4
         nodelist(2)=2
         nodelist(3)=0
        else
         print *,"id invalid"
         stop
        endif

        do j_tet_node=1,sdim+1
         inode=nodelist(j_tet_node)
         do k_dir=1,sdim
          if (inode.eq.0) then
           xavg=zero
           do jnode=1,4
            xavg=xavg+xnode(jnode,k_dir)
           enddo
           xtet(j_tet_node,k_dir)=xavg/four
          else
           xtet(j_tet_node,k_dir)=xnode(inode,k_dir)
          endif
         enddo
         if (inode.eq.0) then
          davg=zero
          do jnode=1,4
           davg=davg+datanode(jnode)
          enddo
          datatet(j_tet_node)=davg/four
         else
          datatet(j_tet_node)=datanode(inode)
         endif
        enddo  ! j_tet_node

       else if (sdim.eq.3) then
        print *,"symmetric tetrahedrazation for 3d not complete"
        stop
       else
        print *,"sdim invalid"
        stop
       endif

      else
       print *,"symmetry_flag invalid"
       stop
      endif

      return
      end subroutine extract_tet 



! find area of a face of a tetrahedron
! find the centroid of a face of a tetrahedron
      subroutine areaXYZ(x,i1,i2,i3,area)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x(4,3)
      real(amrex_real) xx(3,3)
      real(amrex_real), INTENT(out) :: area
      integer, INTENT(in) :: i1,i2,i3 
      integer dir
      integer sdim

      sdim=3

      if ((i1.lt.1).or.(i1.gt.4).or.(i2.lt.1).or.(i2.gt.4).or. &
          (i3.lt.1).or.(i3.gt.4).or.(i1.eq.i2).or.(i1.eq.i3).or. &
          (i2.eq.i3)) then
       print *,"index invalid areaXYZ"
       print *,"i1,i2,i3 ",i1,i2,i3
       stop
      endif

      do dir=1,sdim
       xx(1,dir)=x(i1,dir) 
       xx(2,dir)=x(i2,dir) 
       xx(3,dir)=x(i3,dir) 
      enddo
      call surface_area(xx,area,sdim)

      return
      end subroutine areaXYZ

! find the area of a triangular planar element in 2d or 3d
! Acknowledgement: JOHN BURKHARDT
      subroutine surface_area(x,area,sdim)
      use probcommon_module
      use LegendreNodes
      use triangle_fekete_module, only : fekete_degree,fekete_order_num, &
                                         fekete_rule
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: x(sdim,sdim)
      real(amrex_real), INTENT(out) :: area
      real(amrex_real) xx(sdim-1,sdim)
      integer i,dir,dircrit,itan,jtan
      real(amrex_real), dimension(:,:), allocatable :: xy
      real(amrex_real), dimension(:), allocatable :: w
      real(amrex_real) y1_cross_y2(sdim)
      real(amrex_real) coeff(sdim)
      real(amrex_real) mag,DA
      integer rule,degree,order_num  ! rule=1 degree=3
      real(amrex_real) rho,jac
      real(amrex_real) dxpos(sdim)

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust surface area"
       stop
      endif

      do i=1,sdim-1
      do dir=1,sdim
       xx(i,dir)=x(i,dir)-x(sdim,dir)
      enddo
      enddo

      if (sdim.eq.3) then
       y1_cross_y2(1)=xx(1,2)*xx(sdim-1,sdim)-xx(sdim-1,2)*xx(1,sdim)
       y1_cross_y2(2)=-(xx(1,1)*xx(sdim-1,sdim)-xx(sdim-1,1)*xx(1,sdim))
       y1_cross_y2(sdim)=xx(1,1)*xx(sdim-1,2)-xx(sdim-1,1)*xx(1,2)
      else if (sdim.eq.2) then
       y1_cross_y2(1)=xx(1,2)
       y1_cross_y2(2)=-xx(1,1)
      else
       print *,"dimension bust"
       stop
      endif
      
      mag=zero
      do dir=1,sdim
       mag=mag+y1_cross_y2(dir)**2
      enddo
      mag=sqrt(mag)

      if (sdim.eq.3) then
       area=half*mag  ! mag is area of parallelogram
      else if (sdim.eq.2) then
       area=mag
      else
       print *,"dimension bust"
       stop
      endif

      if (mag.gt.zero) then

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then

        dircrit=1
        if (abs(y1_cross_y2(2)).gt.abs(y1_cross_y2(dircrit))) then
         dircrit=2
        endif
        if (sdim.eq.3) then
         if (abs(y1_cross_y2(sdim)).gt.abs(y1_cross_y2(dircrit))) then
          dircrit=sdim
         endif
        endif
        if (abs(y1_cross_y2(dircrit)).gt.zero) then
         ! do nothing
        else
         print *,"y1_cross_y2 bust"
         stop
        endif

        do dir=1,sdim
         if (dir.eq.dircrit) then
          if (abs(y1_cross_y2(dir)).gt.zero) then
           coeff(dir)=one
          else
           print *,"y1_cross_y2(dircrit) bad, value=",y1_cross_y2(dir)
           stop
          endif
         else if (dir.ne.dircrit) then
          coeff(dir)=y1_cross_y2(dir)/y1_cross_y2(dircrit)
          if (abs(coeff(dir)).le.one) then
           ! do nothing
          else
           print *,"coeff(dir) out of range, dir,value=",dir,coeff(dir)
           stop
          endif
         else
          print *,"dir or dircrit bust"
          stop
         endif
        enddo !dir=1..sdim

        if (dircrit.eq.1) then
         itan=2
         jtan=3
        else if (dircrit.eq.2) then
         itan=1
         jtan=3
        else if ((dircrit.eq.3).and.(sdim.eq.3)) then
         itan=1
         jtan=2
        else
         print *,"dircrit invalid"
         stop
        endif

       else
        print *,"levelrz invalid"
        stop
       endif

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then

        area=zero
        jac=abs(y1_cross_y2(dircrit))

        if (sdim.ne.2) then
         print *,"dimension bust"
         stop
        endif
        degree=3
        order_num=2
        allocate(xy(sdim-1,order_num))
        allocate(w(order_num))

        if (GQTYPE.eq.0) then
         ! do nothing (Legendre Gauss)
        else if (GQTYPE.eq.1) then
         print *,"Clenshaw Curtis has problems"
         stop
        else
         print *,"GQTYPE invalid"
         stop
        endif

        do i=1,order_num
         xy(1,i)=(cache_gauss(order_num,i-1,GQTYPE)+one)/two
         w(i)=cache_gauss_w(order_num,i-1,GQTYPE)*half
        
         dxpos(itan)=xx(1,itan)*xy(1,i)
         dxpos(dircrit)=-(coeff(itan)*dxpos(itan))
         rho=x(sdim,1)+dxpos(1)
        
         ! F(rho,theta,z)=0
         ! area=integral_R |grad F|/|grad F dot p| dR
         ! p is normal to R.   
         ! dR=2 pi rho dz or
         ! dR=2 pi rho drho
         DA=two*Pi*abs(rho)*sqrt( one+coeff(itan)**2 )
 
         area=area+DA*jac*w(i)

        enddo ! i

        if (area.gt.zero) then
         ! do nothing
        else
         print *,"area became 0 even though mag>0, area=",area
         stop
        endif

        deallocate(xy)
        deallocate(w)

       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then ! in: surface_area

        area=zero
        jac=abs(y1_cross_y2(dircrit))

        if (sdim.eq.2) then
         degree=3
         order_num=2
         allocate(xy(sdim-1,order_num))
         allocate(w(order_num))

         do i=1,order_num
          xy(1,i)=(cache_gauss(order_num,i-1,GQTYPE)+one)/two
          w(i)=cache_gauss_w(order_num,i-1,GQTYPE)*half
        
          dxpos(itan)=xx(1,itan)*xy(1,i)
          dxpos(dircrit)=-(coeff(itan)*dxpos(itan))
          rho=x(sdim,1)+dxpos(1)
        
          ! F(rho,theta,z)=0
          ! area=integral_R |grad F|/|grad F dot p| dR
          ! p is normal to R.   
          ! dR=rho dtheta or
          ! dR=drho
          if (dircrit.eq.1) then
           DA=sqrt( rho**2+coeff(itan)**2 )
          else if (dircrit.eq.2) then
           DA=sqrt( one+(coeff(itan)*rho)**2 )
          else
           print *,"dircrit invalid"
           stop
          endif

          area=area+DA*jac*w(i)

         enddo ! i

         if (area.gt.zero) then
          ! do nothing
         else
          print *,"area became 0 even though mag>0, area=",area
          stop
         endif

         deallocate(xy)
         deallocate(w)

        else if (sdim.eq.3) then

         rule=1  ! degree=3
         call fekete_degree(rule,degree)
         if (degree.ne.3) then
          print *,"degree invalid"
          stop
         endif
         call fekete_order_num(rule,order_num)
         allocate(xy(sdim-1,order_num))
         allocate(w(order_num))
         call fekete_rule(rule,order_num,xy,w)

         do i=1,order_num
          dxpos(itan)=xx(1,itan)*xy(1,i)+xx(2,itan)*xy(2,i)
          dxpos(jtan)=xx(1,jtan)*xy(1,i)+xx(2,jtan)*xy(2,i)
          dxpos(dircrit)=-(coeff(itan)*dxpos(itan)+coeff(jtan)*dxpos(jtan))
          rho=x(sdim,1)+dxpos(1)
        
          ! F(rho,theta,z)=0
          ! area=integral_R |grad F|/|grad F dot p| dR
          ! p is normal to R.   
          ! dR=rho dtheta dz or
          ! dR=drho dz or
          ! dR=rho drho dtheta
          if (dircrit.eq.1) then 
           DA=sqrt( (coeff(jtan)*rho)**2+rho**2+coeff(itan)**2 )
          else if (dircrit.eq.2) then
           DA=sqrt( (coeff(jtan)*rho)**2+one+(coeff(itan)*rho)**2 )
          else if (dircrit.eq.3) then
           DA=sqrt( (coeff(jtan))**2+rho**2+(coeff(itan)*rho)**2 )
          else
           print *,"dircrit invalid"
           stop
          endif
 
          area=area+DA*jac*w(i)

         enddo ! i

         if (area.gt.zero) then
          ! do nothing
         else
          print *,"area became 0 even though mag>0, area=",area
          stop
         endif

         area=area*half

         deallocate(xy)
         deallocate(w)
        else
         print *,"dimension bust"
         stop
        endif
       else
        print *,"levelrz invalid ",levelrz
        stop
       endif
      else if (mag.eq.zero) then
       ! do nothing
      else
       print *,"mag invalid in surface_area; mag=",mag
       print *,"sdim=",sdim
       do i=1,sdim
       do dir=1,sdim
        print *,"triangle index i=",i
        print *,"dir=",dir
        print *,"x(i,dir)=",x(i,dir)
       enddo
       enddo
       stop
      endif

      return
      end subroutine surface_area

! find the volume of a tetrahedron (in 3d) or area of triangle (2d)
      subroutine tetrahedron_volume(x,volume,centroid,sdim)
      use probcommon_module
      use tetrahedron_keast_module, only : keast_degree,keast_order_num, &
                                           keast_rule
      use triangle_fekete_module, only : fekete_degree,fekete_order_num, &
                                         fekete_rule

      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: volume
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real) centroid_def(sdim)
      real(amrex_real) xx(sdim,sdim)
      integer i_tet_node
      integer i_order
      integer j_dir,j_vec,dir
      real(amrex_real), dimension(:,:), allocatable :: xyz
      real(amrex_real), dimension(:), allocatable :: w
      integer rule,degree,order_num
      real(amrex_real) rho,volzero
      real(amrex_real) xmin(sdim)
      real(amrex_real) xmax(sdim)
      real(amrex_real) :: dxmax
      real(amrex_real) dxpos(sdim)
      real(amrex_real) total_weight
      real(amrex_real) centroid_eps

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid"
       stop
      endif

      call get_xbounds(x,xmin,xmax,dxmax,sdim+1,sdim)
  
      do i_tet_node=1,sdim
      do dir=1,sdim
       xx(i_tet_node,dir)=x(i_tet_node,dir)-x(sdim+1,dir)
      enddo
      enddo

      if (sdim.eq.3) then
       volzero= &
        xx(2,1)*(xx(1,2)*xx(sdim,sdim)-xx(sdim,2)*xx(1,sdim))- &
        xx(1,1)*(xx(2,2)*xx(sdim,sdim)-xx(sdim,2)*xx(2,sdim))+ &
        xx(sdim,1)*(xx(2,2)*xx(1,sdim)-xx(1,2)*xx(2,sdim))
      else if (sdim.eq.2) then
       volzero=xx(1,1)*xx(2,2)-xx(2,1)*xx(1,2)
      else
       print *,"dimension bust"
       stop
      endif
      volzero=abs(volzero)

      do dir=1,sdim
       centroid(dir)=zero
       do i_tet_node=1,sdim+1
        centroid(dir)=centroid(dir)+x(i_tet_node,dir)
       enddo
       centroid(dir)=centroid(dir)/(sdim+one)
       centroid_def(dir)=centroid(dir)
      enddo ! dir

      if (levelrz.eq.COORDSYS_CARTESIAN) then

       if (sdim.eq.3) then
        volume=volzero/six
       else if (sdim.eq.2) then
        volume=volzero/two
       else
        print *,"sdim invalid"
        stop
       endif
  
      else if (levelrz.eq.COORDSYS_RZ) then

       if (sdim.ne.2) then
        print *,"dimension bust"
        stop
       endif

       rule=1  ! degree=3
       call fekete_degree(rule,degree)
       if (degree.ne.3) then
        print *,"degree invalid"
        stop
       endif
       call fekete_order_num(rule,order_num)
       allocate(xyz(sdim,order_num))
       allocate(w(order_num))
       call fekete_rule(rule,order_num,xyz,w)

       if (1.eq.0) then
        print *,"---------------------------"
        total_weight=zero
        do i_order=1,order_num
         print *,"i_order,x,y,w ", &
           i_order,xyz(1,i_order),xyz(2,i_order),w(i_order)
         total_weight=total_weight+w(i_order)
        enddo
        print *,"total_weight: ",total_weight
        print *,"---------------------------"
       endif

       volume=zero
       do dir=1,sdim
        centroid(dir)=zero
       enddo
       do i_order=1,order_num
        do dir=1,sdim
         dxpos(dir)=zero
         do j_vec=1,sdim
          dxpos(dir)=dxpos(dir)+xx(j_vec,dir)*xyz(j_vec,i_order)
         enddo
        enddo
        rho=x(sdim+1,1)+dxpos(1)
        volume=volume+two*Pi*abs(rho)*w(i_order)
        do dir=1,sdim
         centroid(dir)=centroid(dir)+(x(sdim+1,dir)+dxpos(dir))* &
           two*Pi*abs(rho)*w(i_order)
        enddo
       enddo ! i_order=1,order_num
       if (volume.gt.zero) then
        do dir=1,sdim
         centroid(dir)=centroid(dir)/volume
        enddo
       else
        do dir=1,sdim
         centroid(dir)=centroid_def(dir)
        enddo
       endif

       volume=abs(volume)*volzero*half

       deallocate(xyz)
       deallocate(w)

      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then

       if (sdim.eq.2) then

        rule=1  ! degree=3
        call fekete_degree(rule,degree)
        if (degree.ne.3) then
         print *,"degree invalid"
         stop
        endif
        call fekete_order_num(rule,order_num)
        allocate(xyz(sdim,order_num))
        allocate(w(order_num))
        call fekete_rule(rule,order_num,xyz,w)
        volume=zero
        do dir=1,sdim
         centroid(dir)=zero
        enddo
        do i_order=1,order_num
         do dir=1,sdim
          dxpos(dir)=zero
          do j_vec=1,sdim
           dxpos(dir)=dxpos(dir)+xx(j_vec,dir)*xyz(j_vec,i_order)
          enddo
         enddo
         rho=x(sdim+1,1)+dxpos(1)
         volume=volume+abs(rho)*w(i_order)
         do dir=1,sdim
          centroid(dir)=centroid(dir)+ &
             (x(sdim+1,dir)+dxpos(dir))*abs(rho)*w(i_order)
         enddo
        enddo ! i=1,order_num
        if (volume.gt.zero) then
         do dir=1,sdim
          centroid(dir)=centroid(dir)/volume
         enddo
        else
         do dir=1,sdim
          centroid(dir)=centroid_def(dir)
         enddo
        endif

        volume=abs(volume)*volzero*half

        deallocate(xyz)
        deallocate(w)

       else if (sdim.eq.3) then

        rule=2  ! linear integrands exactly
        call keast_degree(rule,degree)
        if (degree.ne.1) then
         print *,"degree invalid"
         stop
        endif
        call keast_order_num(rule,order_num)
        allocate(xyz(sdim,order_num))
        allocate(w(order_num))
        call keast_rule(rule,order_num,xyz,w)
        volume=zero
        do dir=1,sdim
         centroid(dir)=zero
        enddo
        do i_order=1,order_num
         do dir=1,sdim
          dxpos(dir)=zero
          do j_vec=1,sdim
           dxpos(dir)=dxpos(dir)+xx(j_vec,dir)*xyz(j_vec,i_order)
          enddo
         enddo
         rho=x(sdim+1,1)+dxpos(1)
         volume=volume+abs(rho)*w(i_order)
         do dir=1,sdim
          centroid(dir)=centroid(dir)+(x(sdim+1,dir)+dxpos(dir))* &
            abs(rho)*w(i_order)
         enddo ! dir
        enddo ! i_order=1,...,order_num
        if (volume.gt.zero) then
         do dir=1,sdim
          centroid(dir)=centroid(dir)/volume
         enddo
        else
         do dir=1,sdim
          centroid(dir)=centroid_def(dir)
         enddo
        endif

        volume=abs(volume)*volzero*sixth

        deallocate(xyz)
        deallocate(w)

       else
        print *,"dimension bust"
        stop
       endif

      else
       print *,"levelrz invalid ",levelrz
       stop
      endif

      do j_dir=1,sdim

       if (xmin(j_dir).gt.xmax(j_dir)) then
        print *,"xmin(j_dir).gt.xmax(j_dir)"
        stop
       else if (xmin(j_dir).le.xmax(j_dir)) then
        ! do nothing
       else
        print *,"xmin or xmax is NaN"
        stop
       endif

       if (dxmax.ge.zero) then
        !do nothing
       else
        print *,"dxmax invalid: ",dxmax
        stop
       endif

       centroid_eps=EPS_13_5*max(dxmax,one)

       if (centroid(j_dir)+centroid_eps.lt.xmin(j_dir)) then
        print *,"WARN centroid(j_dir)+centroid_eps.lt.xmin(j_dir) XYZ"
        print *,"j_dir,centroid,EPS_13_5,dxmax,xmin ",j_dir, &
         centroid(j_dir),EPS_13_5,dxmax,xmin(j_dir)
        centroid(j_dir)=xmin(j_dir)
       else if (centroid(j_dir)-centroid_eps.gt.xmax(j_dir)) then
        print *,"WARN centroid(j_dir)-centroid_eps.gt.xmax(j_dir) XYZ"
        print *,"j_dir,centroid,EPS_13_5,dxmax,xmax ",j_dir, &
         centroid(j_dir),EPS_13_5,dxmax,xmax(j_dir)
        centroid(j_dir)=xmax(j_dir)
       endif

       if ((centroid(j_dir)+centroid_eps.lt.xmin(j_dir)).or. &
           (centroid(j_dir)-centroid_eps.gt.xmax(j_dir))) then
        print *,"centroid still invalid XYZ"
        stop
       endif

      enddo ! j=1,sdim

      return
      end subroutine tetrahedron_volume


! find the volume of a tetrahedron (in 3d) or area of triangle (2d)
      subroutine tetrahedron_volume_and_map(normdir,coeff, &
       x,volume,centroid,volume_map,centroid_map,sdim)
      use probcommon_module
      use tetrahedron_keast_module, only : keast_degree,keast_order_num, &
                                           keast_rule
      use triangle_fekete_module, only : fekete_degree,fekete_order_num, &
                                         fekete_rule

      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: normdir
      real(amrex_real), INTENT(in) :: coeff(2)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: volume
      real(amrex_real), INTENT(out) :: volume_map
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real), INTENT(out) :: centroid_map(sdim)
      real(amrex_real) centroid_def(sdim)
      real(amrex_real) centroid_def_map(sdim)
      real(amrex_real) xx(sdim,sdim)
      integer i_tet_node
      integer i_order
      integer j_dir,j_vec,dir
      real(amrex_real), dimension(:,:), allocatable :: xyz
      real(amrex_real), dimension(:), allocatable :: w
      integer rule,degree,order_num
      real(amrex_real) rho
      real(amrex_real) volzero
      real(amrex_real) volzero_map
      real(amrex_real) xmin(sdim)
      real(amrex_real) xmax(sdim)
      real(amrex_real) :: dxmax
      real(amrex_real) xmin_map(sdim)
      real(amrex_real) xmax_map(sdim)
      real(amrex_real) dxpos(sdim)
      real(amrex_real) total_weight
      real(amrex_real) xfactor
      real(amrex_real) centroid_eps

      if ((normdir.lt.0).or.(normdir.ge.sdim)) then
       print *,"normdir invalid"
       stop
      endif
      if (coeff(1).le.zero) then
       print *,"coeff(1) invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid"
       stop
      endif

      call get_xbounds(x,xmin,xmax,dxmax,sdim+1,sdim)
  
      do i_tet_node=1,sdim
      do dir=1,sdim
       xx(i_tet_node,dir)=x(i_tet_node,dir)-x(sdim+1,dir)
      enddo
      enddo

      if (sdim.eq.3) then
       volzero= &
        xx(2,1)*(xx(1,2)*xx(sdim,sdim)-xx(sdim,2)*xx(1,sdim))- &
        xx(1,1)*(xx(2,2)*xx(sdim,sdim)-xx(sdim,2)*xx(2,sdim))+ &
        xx(sdim,1)*(xx(2,2)*xx(1,sdim)-xx(1,2)*xx(2,sdim))
       volzero_map=volzero*coeff(1)
      else if (sdim.eq.2) then
       volzero=xx(1,1)*xx(2,2)-xx(2,1)*xx(1,2)
       volzero_map=volzero*coeff(1)
      else
       print *,"dimension bust"
       stop
      endif
      volzero=abs(volzero)
      volzero_map=abs(volzero_map)

      do dir=1,sdim
       centroid(dir)=zero
       do i_tet_node=1,sdim+1
        centroid(dir)=centroid(dir)+x(i_tet_node,dir)
       enddo
       centroid(dir)=centroid(dir)/(sdim+one)
       centroid_def(dir)=centroid(dir)
       if (dir.eq.normdir+1) then
        centroid_map(dir)=coeff(1)*centroid(dir)+coeff(2)
       else
        centroid_map(dir)=centroid(dir)
       endif
       centroid_def_map(dir)=centroid_map(dir)
      enddo ! dir

      if (levelrz.eq.COORDSYS_CARTESIAN) then

       if (sdim.eq.3) then
        volume=volzero/six
        volume_map=volzero_map/six
       else if (sdim.eq.2) then
        volume=volzero/two
        volume_map=volzero_map/two
       else
        print *,"sdim invalid"
        stop
       endif
  
      else if (levelrz.eq.COORDSYS_RZ) then

       if (sdim.ne.2) then
        print *,"dimension bust"
        stop
       endif

       rule=1  ! degree=3
       call fekete_degree(rule,degree)
       if (degree.ne.3) then
        print *,"degree invalid"
        stop
       endif
       call fekete_order_num(rule,order_num)
       allocate(xyz(sdim,order_num))
       allocate(w(order_num))
       call fekete_rule(rule,order_num,xyz,w)

       if (1.eq.0) then
        print *,"---------------------------"
        total_weight=zero
        do i_order=1,order_num
         print *,"i_order,x,y,w ", &
           i_order,xyz(1,i_order),xyz(2,i_order),w(i_order)
         total_weight=total_weight+w(i_order)
        enddo
        print *,"total_weight: ",total_weight
        print *,"---------------------------"
       endif

       volume=zero
       volume_map=zero
       do dir=1,sdim
        centroid(dir)=zero
        centroid_map(dir)=zero
       enddo
       do i_order=1,order_num
        do dir=1,sdim
         dxpos(dir)=zero
         do j_vec=1,sdim
          dxpos(dir)=dxpos(dir)+xx(j_vec,dir)*xyz(j_vec,i_order)
         enddo
        enddo

        rho=x(sdim+1,1)+dxpos(1)
        volume=volume+two*Pi*abs(rho)*w(i_order)
        do dir=1,sdim
         xfactor=x(sdim+1,dir)+dxpos(dir)
         centroid(dir)=centroid(dir)+xfactor*two*Pi*abs(rho)*w(i_order)
        enddo

        if (normdir.eq.0) then
         rho=coeff(1)*rho+coeff(2)
        endif
        volume_map=volume_map+two*Pi*abs(rho)*w(i_order)
        do dir=1,sdim
         xfactor=x(sdim+1,dir)+dxpos(dir)
         if (dir.eq.normdir+1) then
          xfactor=coeff(1)*xfactor+coeff(2)
         endif
         centroid_map(dir)=centroid_map(dir)+xfactor*two*Pi*abs(rho)*w(i_order)
        enddo

       enddo ! i_order=1,order_num
       if ((volume.gt.zero).and.(volume_map.gt.zero)) then
        do dir=1,sdim
         centroid(dir)=centroid(dir)/volume
         centroid_map(dir)=centroid_map(dir)/volume_map
        enddo
       else if ((volume.eq.zero).or.(volume_map.eq.zero)) then
        do dir=1,sdim
         centroid(dir)=centroid_def(dir)
         centroid_map(dir)=centroid_def_map(dir)
        enddo
       else
        print *,"volume or volume_map invalid"
        stop
       endif

       volume=abs(volume)*volzero*half
       volume_map=abs(volume_map)*volzero_map*half

       deallocate(xyz)
       deallocate(w)

      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then

       if (sdim.eq.2) then

        rule=1  ! degree=3
        call fekete_degree(rule,degree)
        if (degree.ne.3) then
         print *,"degree invalid"
         stop
        endif
        call fekete_order_num(rule,order_num)
        allocate(xyz(sdim,order_num))
        allocate(w(order_num))
        call fekete_rule(rule,order_num,xyz,w)
        volume=zero
        volume_map=zero
        do dir=1,sdim
         centroid(dir)=zero
         centroid_map(dir)=zero
        enddo
        do i_order=1,order_num
         do dir=1,sdim
          dxpos(dir)=zero
          do j_vec=1,sdim
           dxpos(dir)=dxpos(dir)+xx(j_vec,dir)*xyz(j_vec,i_order)
          enddo
         enddo

         rho=x(sdim+1,1)+dxpos(1)
         volume=volume+abs(rho)*w(i_order)
         do dir=1,sdim
          xfactor=x(sdim+1,dir)+dxpos(dir)
          centroid(dir)=centroid(dir)+xfactor*abs(rho)*w(i_order)
         enddo

         if (normdir.eq.0) then
          rho=coeff(1)*rho+coeff(2)
         endif
         volume_map=volume_map+abs(rho)*w(i_order)
         do dir=1,sdim
          xfactor=x(sdim+1,dir)+dxpos(dir)
          if (dir.eq.normdir+1) then
           xfactor=coeff(1)*xfactor+coeff(2)
          endif
          centroid_map(dir)=centroid_map(dir)+xfactor*abs(rho)*w(i_order)
         enddo

        enddo ! i_order=1,order_num
        if ((volume.gt.zero).and.(volume_map.gt.zero)) then
         do dir=1,sdim
          centroid(dir)=centroid(dir)/volume
          centroid_map(dir)=centroid_map(dir)/volume_map
         enddo
        else if ((volume.eq.zero).or.(volume_map.eq.zero)) then
         do dir=1,sdim
          centroid(dir)=centroid_def(dir)
          centroid_map(dir)=centroid_def_map(dir)
         enddo
        else
         print *,"volume or volume_map invalid"
         stop
        endif

        volume=abs(volume)*volzero*half
        volume_map=abs(volume_map)*volzero_map*half

        deallocate(xyz)
        deallocate(w)

       else if (sdim.eq.3) then

        rule=2  ! linear integrands exactly
        call keast_degree(rule,degree)
        if (degree.ne.1) then
         print *,"degree invalid"
         stop
        endif
        call keast_order_num(rule,order_num)
        allocate(xyz(sdim,order_num))
        allocate(w(order_num))
        call keast_rule(rule,order_num,xyz,w)
        volume=zero
        volume_map=zero
        do dir=1,sdim
         centroid(dir)=zero
         centroid_map(dir)=zero
        enddo
        do i_order=1,order_num
         do dir=1,sdim
          dxpos(dir)=zero
          do j_vec=1,sdim
           dxpos(dir)=dxpos(dir)+xx(j_vec,dir)*xyz(j_vec,i_order)
          enddo
         enddo

         rho=x(sdim+1,1)+dxpos(1)
         volume=volume+abs(rho)*w(i_order)
         do dir=1,sdim
          xfactor=x(sdim+1,dir)+dxpos(dir)
          centroid(dir)=centroid(dir)+xfactor*abs(rho)*w(i_order)
         enddo ! dir

         if (normdir.eq.0) then
          rho=coeff(1)*rho+coeff(2)
         endif
         volume_map=volume_map+abs(rho)*w(i_order)
         do dir=1,sdim
          xfactor=x(sdim+1,dir)+dxpos(dir)
          if (dir.eq.normdir+1) then
           xfactor=coeff(1)*xfactor+coeff(2)
          endif
          centroid_map(dir)=centroid_map(dir)+xfactor*abs(rho)*w(i_order)
         enddo

        enddo ! i_order=1,...,order_num
        if ((volume.gt.zero).and.(volume_map.gt.zero)) then
         do dir=1,sdim
          centroid(dir)=centroid(dir)/volume
          centroid_map(dir)=centroid_map(dir)/volume_map
         enddo
        else if ((volume.eq.zero).or.(volume_map.eq.zero)) then
         do dir=1,sdim
          centroid(dir)=centroid_def(dir)
          centroid_map(dir)=centroid_def_map(dir)
         enddo
        else
         print *,"volume or volume_map invalid"
         stop
        endif

        volume=abs(volume)*volzero*sixth
        volume_map=abs(volume_map)*volzero_map*sixth

        deallocate(xyz)
        deallocate(w)

       else
        print *,"dimension bust"
        stop
       endif

      else
       print *,"levelrz invalid ",levelrz
       stop
      endif

      do j_dir=1,sdim

       if (xmin(j_dir).gt.xmax(j_dir)) then
        print *,"xmin(j_dir).gt.xmax(j_dir)"
        stop
       else if (xmin(j_dir).le.xmax(j_dir)) then
        ! do nothing
       else
        print *,"xmin or xmax is NaN"
        stop
       endif

       if (dxmax.ge.zero) then
        !do nothing
       else
        print *,"dxmax invalid"
        stop
       endif

       centroid_eps=EPS_13_5*max(dxmax,one)

       if (centroid(j_dir)+centroid_eps.lt.xmin(j_dir)) then
        print *,"WARN centroid(j_dir)+centroid_eps.lt.xmin(j_dir) XYZ"
        print *,"j_dir,centroid,EPS_13_5,dxmax,xmin ",j_dir, &
         centroid(j_dir),EPS_13_5,dxmax,xmin(j_dir)
        centroid(j_dir)=xmin(j_dir)
       else if (centroid(j_dir)-centroid_eps.gt.xmax(j_dir)) then
        print *,"WARN centroid(j_dir)-centroid_eps.gt.xmax(j_dir) XYZ"
        print *,"j_dir,centroid,EPS_13_5,dxmax,xmax ",j_dir, &
         centroid(j_dir),EPS_13_5,dxmax,xmax(j_dir)
        centroid(j_dir)=xmax(j_dir)
       endif

       if ((centroid(j_dir)+centroid_eps.lt.xmin(j_dir)).or. &
           (centroid(j_dir)-centroid_eps.gt.xmax(j_dir))) then
        print *,"centroid still invalid XYZ"
        stop
       endif

       xmin_map(j_dir)=xmin(j_dir)
       xmax_map(j_dir)=xmax(j_dir)
       if (j_dir.eq.normdir+1) then
        xmin_map(j_dir)=coeff(1)*xmin_map(j_dir)+coeff(2)
        xmax_map(j_dir)=coeff(1)*xmax_map(j_dir)+coeff(2)
       endif

       if (xmin_map(j_dir).gt.xmax_map(j_dir)) then
        print *,"xmin_map(j_dir).gt.xmax_map(j_dir)"
        stop
       endif

       if (centroid_map(j_dir)+EPS_13_5*dxmax.lt.xmin_map(j_dir)) then
        print *,"WARN: centroid_map(j_dir)+EPS_13_5*dxmax.lt.xmin_map(j_dir) XYZ"
        centroid_map(j_dir)=xmin_map(j_dir)
       else if (centroid_map(j_dir)-EPS_13_5*dxmax.gt.xmax_map(j_dir)) then
        print *,"WARN: centroid_map(j_dir)-EPS_13_5*dxmax.gt.xmax_map(j_dir) XYZ"
        centroid_map(j_dir)=xmax_map(j_dir)
       endif

       if ((centroid_map(j_dir)+EPS_13_5*dxmax.lt.xmin_map(j_dir)).or. &
           (centroid_map(j_dir)-EPS_13_5*dxmax.gt.xmax_map(j_dir))) then
        print *,"centroid_map still invalid XYZ"
        stop
       endif

      enddo ! j=1,sdim

      return
      end subroutine tetrahedron_volume_and_map


! internal routine
      subroutine shrink2D(xint,x,phi,isrc,itarg,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(out) :: xint(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      integer, INTENT(in) :: isrc,itarg
      integer j_tet_node
      integer j_dir

      if (sdim.ne.2) then
       print *,"sdim bust shring 2d"
       stop
      endif

      if ((itarg.lt.1).or.(itarg.gt.3).or.(isrc.lt.1).or. &
          (isrc.gt.3).or.(itarg.eq.isrc)) then
       print *,"index invalid shrink2d "
       do j_tet_node=1,3
        print *,"j_tet_node,xint ", &
           j_tet_node,xint(j_tet_node,1),xint(j_tet_node,2)
        print *,"j_tet_node,x ",j_tet_node,x(j_tet_node,1),x(j_tet_node,2)
        print *,"j_tet_node,phi ",j_tet_node,phi(j_tet_node)
       enddo
       stop
      endif
      if (((phi(itarg).ge.zero).and.(phi(isrc).ge.zero)).or. &
          ((phi(itarg).lt.zero).and.(phi(isrc).lt.zero))) then
       print *,"levelset does not change sign"
       stop
      endif

      do j_dir=1,2
       xint(itarg,j_dir)=(abs(phi(itarg))*x(isrc,j_dir)+ &
                          abs(phi(isrc))*x(itarg,j_dir))/ &
                         (abs(phi(itarg))+abs(phi(isrc)))
      enddo

      return
      end subroutine shrink2D

! internal routine
      subroutine shrink3D(xint,x,phi,isrc,itarg,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(out) :: xint(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      integer, INTENT(in) :: isrc,itarg
      integer j_tet_node
      integer j_dir

      if (sdim.ne.3) then
       print *,"sdim bust shrink 3d"
       stop
      endif

      if ((itarg.lt.1).or.(itarg.gt.4).or.(isrc.lt.1).or. &
          (isrc.gt.4).or.(itarg.eq.isrc)) then
       print *,"index invalid shrink3d "
       do j_tet_node=1,4
        print *,"j_tet_node,xint ",j_tet_node, &
            xint(j_tet_node,1),xint(j_tet_node,2),xint(j_tet_node,3)
        print *,"j_tet_node,x ",j_tet_node, &
            x(j_tet_node,1),x(j_tet_node,2),x(j_tet_node,3)
        print *,"j_tet_node,phi ",j_tet_node,phi(j_tet_node)
       enddo
       stop
      endif
      if (((phi(itarg).ge.zero).and.(phi(isrc).ge.zero)).or. &
          ((phi(itarg).lt.zero).and.(phi(isrc).lt.zero))) then
       print *,"levelset does not change sign"
       stop
      endif

      do j_dir=1,3
       xint(itarg,j_dir)=(abs(phi(itarg))*x(isrc,j_dir)+ &
                          abs(phi(isrc))*x(itarg,j_dir))/ &
                         (abs(phi(itarg))+abs(phi(isrc)))
      enddo

      return
      end subroutine shrink3D

 
! internal routine, do not call
      subroutine shrink_list_tri(phi,x,i1,i2,i3,complement, &
        xtrilist, &
        nlist_alloc,nlist,xarealist,narea)
      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      real(amrex_real), INTENT(out) :: xtrilist(3,2,nlist_alloc)
      real(amrex_real), INTENT(out) :: xarealist(2,2,MAXAREA)
      integer, INTENT(inout) :: nlist,narea
      real(amrex_real), INTENT(in) :: phi(3)
      real(amrex_real), INTENT(in) :: x(3,2)
      real(amrex_real) xint(3,2)
      integer i_tet_node
      integer j_dir
      integer, INTENT(in) :: i1,i2,i3,complement
      integer sdim

      sdim=2

      if (nlist_alloc.ge.nlist+1) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif 

      if (complement.eq.0) then
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
       enddo
       call shrink2D(xint,x,phi,i1,i2,sdim)
       call shrink2D(xint,x,phi,i1,i3,sdim)
       nlist=nlist+1 
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif 
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtrilist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
       enddo

       narea=narea+1
       do j_dir=1,sdim
        xarealist(1,j_dir,narea)=xint(i2,j_dir)
        xarealist(2,j_dir,narea)=xint(i3,j_dir)
       enddo
      else if (complement.eq.1) then
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
       enddo
       call shrink2D(xint,x,phi,i2,i1,sdim)
       nlist=nlist+1  
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif 
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtrilist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
       enddo

       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
       enddo
       call shrink2D(xint,x,phi,i1,i2,sdim)
       call shrink2D(xint,x,phi,i3,i1,sdim)
       nlist=nlist+1  
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif 
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtrilist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
       enddo

       narea=narea+1
       if (narea.gt.MAXAREA) then
        print *,"nlist invalid"
        stop
       endif 
       do j_dir=1,sdim
        xarealist(1,j_dir,narea)=xint(i2,j_dir)
        xarealist(2,j_dir,narea)=xint(i1,j_dir)
       enddo
      else
       print *,"complement invalid"
       stop
      endif

      return
      end subroutine shrink_list_tri



! internal routine, do not call 
! was: nlist_alloc=MAXTET
      subroutine shrink_list_tet(phi,x,i1,i2,i3,i4,complement, &
        xtetlist,nlist_alloc,nlist,xarealist,narea)
      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real), INTENT(out) :: xarealist(3,3,MAXAREA)
      integer, INTENT(inout) :: nlist,narea
      real(amrex_real), INTENT(in) :: phi(4)
      real(amrex_real), INTENT(in) :: x(4,3)
      real(amrex_real) xint(4,3)
      integer i_tet_node
      integer j_dir
      integer, INTENT(in) :: i1,i2,i3,i4,complement
      integer sdim
 
      sdim=3
    
      if (nlist_alloc.ge.nlist+1) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif 

      if (complement.eq.0) then
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
       enddo
       call shrink3D(xint,x,phi,i1,i2,sdim)  ! xint,x,phi,isrc,itarg
       call shrink3D(xint,x,phi,i1,i3,sdim)  
       call shrink3D(xint,x,phi,i1,i4,sdim)
       nlist=nlist+1  
       if (nlist.gt.nlist_alloc) then
        print *,"nlist invalid"
        stop
       endif
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtetlist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
       enddo
       narea=narea+1
       do j_dir=1,sdim
        xarealist(1,j_dir,narea)=xint(i2,j_dir)
        xarealist(2,j_dir,narea)=xint(i3,j_dir)
        xarealist(3,j_dir,narea)=xint(i4,j_dir)
       enddo
      else if (complement.eq.1) then
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
       enddo
       call shrink3D(xint,x,phi,i3,i1,sdim)
       nlist=nlist+1
       if (nlist.gt.nlist_alloc) then
        print *,"nlist invalid"
        stop
       endif
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtetlist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
       enddo

       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
       enddo
       call shrink3D(xint,x,phi,i4,i1,sdim)
       call shrink3D(xint,x,phi,i1,i3,sdim)
       nlist=nlist+1
       if (nlist.gt.nlist_alloc) then
        print *,"nlist invalid"
        stop
       endif
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtetlist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
       enddo

       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
       enddo
       call shrink3D(xint,x,phi,i1,i3,sdim)
       call shrink3D(xint,x,phi,i1,i4,sdim)
       call shrink3D(xint,x,phi,i2,i1,sdim)
       nlist=nlist+1
       if (nlist.gt.nlist_alloc) then
        print *,"nlist invalid"
        stop
       endif
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtetlist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
       enddo
       narea=narea+1
       if (narea.gt.MAXAREA) then
        print *,"narea invalid"
        stop
       endif
       do j_dir=1,sdim
        xarealist(1,j_dir,narea)=xint(i3,j_dir)
        xarealist(2,j_dir,narea)=xint(i4,j_dir)
        xarealist(3,j_dir,narea)=xint(i1,j_dir)
       enddo
      else
       print *,"complement invalid"
       stop
      endif

      return
      end subroutine shrink_list_tet

! internal routine, do not call
! i1,i2 nodes are positive
! i3,i4 nodes are negative
      subroutine shrink_gableroof_list(phi,x,i1,i2,i3,i4, &
        xtetlist,nlist_alloc,nlist,xarealist,narea)
      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real), INTENT(out) :: xarealist(3,3,MAXAREA)
      integer, INTENT(inout) :: nlist,narea
      real(amrex_real), INTENT(in) :: phi(4)
      real(amrex_real), INTENT(in) :: x(4,3)
      real(amrex_real) xint(4,3)
      integer, INTENT(in) :: i1,i2,i3,i4
      integer i_tet_node
      integer j_dir

      integer sdim

      sdim=3

      if (nlist_alloc.ge.nlist+1) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif 

      do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
      enddo
      call shrink3D(xint,x,phi,i1,i3,sdim)  ! xint,x,phi,isrc,itarg
      call shrink3D(xint,x,phi,i1,i4,sdim)  
      call shrink3D(xint,x,phi,i4,i2,sdim)  
      nlist=nlist+1
      if (nlist.gt.nlist_alloc) then
        print *,"nlist invalid"
        stop
      endif
      do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtetlist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
      enddo
      narea=narea+1
      if (narea.gt.MAXAREA) then
        print *,"narea invalid"
        stop
      endif
      do j_dir=1,sdim
       xarealist(1,j_dir,narea)=xint(i2,j_dir)
       xarealist(2,j_dir,narea)=xint(i3,j_dir)
       xarealist(3,j_dir,narea)=xint(i4,j_dir)
      enddo

      do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
      enddo
      call shrink3D(xint,x,phi,i1,i3,sdim)
      call shrink3D(xint,x,phi,i2,i4,sdim)
      call shrink3D(xint,x,phi,i3,i2,sdim)
      nlist=nlist+1
      if (nlist.gt.nlist_alloc) then
        print *,"nlist invalid"
        stop
      endif
      do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtetlist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
      enddo
      narea=narea+1
      if (narea.gt.MAXAREA) then
        print *,"narea invalid"
        stop
      endif
      do j_dir=1,sdim
       xarealist(1,j_dir,narea)=xint(i2,j_dir)
       xarealist(2,j_dir,narea)=xint(i3,j_dir)
       xarealist(3,j_dir,narea)=xint(i4,j_dir)
      enddo

      do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=x(i_tet_node,j_dir)
       enddo
      enddo
      call shrink3D(xint,x,phi,i2,i4,sdim)
      call shrink3D(xint,x,phi,i2,i3,sdim)
      nlist=nlist+1
      if (nlist.gt.nlist_alloc) then
        print *,"nlist invalid"
        stop
      endif
      do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtetlist(i_tet_node,j_dir,nlist)=xint(i_tet_node,j_dir)
       enddo
      enddo

      return
      end subroutine shrink_gableroof_list



   
subroutine volume_sanity_check()

use probcommon_module
use global_utility_module

IMPLICIT NONE

integer bfact
real(amrex_real) xsanity(3),dxgrid(3)
real(amrex_real) xsten0(-3:3,3)
real(amrex_real) hangle,hintercept,intercept
real(amrex_real) angle(2)
real(amrex_real) nslope(3)
integer Nangle,sdim,Nangle2,nodedomain
integer i_int,a1,a2,inode,dir
integer isten_loop
integer i_grid_node
integer j_grid_node
integer k_grid_node
real(amrex_real) volslow,areaslow,volall
real(amrex_real) cenall(3)
real(amrex_real) censlow(3)
real(amrex_real) xnode3d(8,3)
real(amrex_real) xnode2d(4,2)
real(amrex_real) phinode(8)
real(amrex_real) xtarget(3)
integer fullelementfast,linearcut,nhalf
real(amrex_real) t1,t2

real(amrex_real) cum_volume,cum_area
real(amrex_real) cum_centroid(3)
real(amrex_real) local_scale

 nhalf=3
 bfact=2
 xsanity(1)=0.125
 xsanity(2)=0.25
 xsanity(3)=0.32
 dxgrid(1)=0.013
 dxgrid(2)=0.031
 dxgrid(3)=0.072
 do isten_loop=-nhalf,nhalf
  do dir=1,3
   xsten0(isten_loop,dir)=xsanity(dir)+isten_loop*dxgrid(dir)*half
  enddo
 enddo

 Nangle=512
 hangle=2*Pi/Nangle
 hintercept=two*dxgrid(3)/Nangle

 call cpu_time(t1)

 print *,"sanity check N=",Nangle
 do sdim=2,3
  print *,"sanity check sdim= ",sdim
  if (sdim.eq.2) then
   Nangle2=0
  endif

  nodedomain=4*(sdim-1)
  print *,"sanity check levelrz=",levelrz
  do i_int=0,Nangle
   intercept=-dxgrid(3)+i_int*hintercept
   do a1=0,Nangle
    angle(1)=a1*hangle

    do a2=0,Nangle2
     angle(2)=a2*hangle
     call angle_to_slope(angle,nslope,sdim)

     inode=1
     if (sdim.eq.3) then
      do k_grid_node=-1,1,2
      do j_grid_node=-1,1,2
      do i_grid_node=-1,1,2
       do dir=1,sdim
        if (dir.eq.1) then
         xnode3d(inode,dir)=xsten0(i_grid_node,dir)
        else if (dir.eq.2) then
         xnode3d(inode,dir)=xsten0(j_grid_node,dir)
        else if (dir.eq.sdim) then
         xnode3d(inode,dir)=xsten0(k_grid_node,dir)
        else
         print *,"dir invalid volume sanity check"
         stop
        endif
        xtarget(dir)=xnode3d(inode,dir)
       enddo  ! dir
       call distfunc(bfact,dxgrid,xsten0,nhalf, &
          intercept,nslope,xtarget, &
          phinode(inode),sdim)
       
       inode=inode+1
      enddo
      enddo
      enddo  ! i,j,k
     else if (sdim.eq.2) then
      do j_grid_node=-1,1,2
      do i_grid_node=-1,1,2
       do dir=1,sdim
        if (dir.eq.1) then
         xnode2d(inode,dir)=xsten0(i_grid_node,dir)
        else if (dir.eq.2) then
         xnode2d(inode,dir)=xsten0(j_grid_node,dir)
        else
         print *,"dir invalid volume sanity check 2"
         stop
        endif
        xtarget(dir)=xnode2d(inode,dir)
       enddo  ! dir
       call distfunc(bfact,dxgrid,xsten0,nhalf, &
          intercept,nslope,xtarget, &
          phinode(inode),sdim)
       
       inode=inode+1
      enddo
      enddo ! i,j
     else
      print *,"sdim invalid"
      stop
     endif

     if (inode.ne.nodedomain+1) then
      print *,"inode invalid1"
      print *,"inode=",inode
      print *,"nodedomain=",nodedomain
      stop
     endif

     call cell_intersection_grid( &
       bfact,dxgrid,xsten0,nhalf, &
       phinode, &
       volslow,censlow,areaslow, &
       volall,cenall,sdim)

     if (1.eq.1) then

      fullelementfast=1
        ! 1 if input is a plane, 0 if input might be plane, 
        ! -1 if force break up of domain.
      linearcut=1
      if (sdim.eq.2) then
       call intersection_volume( &
        cum_volume,cum_area,cum_centroid, &
        phinode,xnode2d,nodedomain, &
        sdim,fullelementfast,linearcut)
      else if (sdim.eq.3) then
       call intersection_volume( &
        cum_volume,cum_area,cum_centroid, &
        phinode,xnode3d,nodedomain, &
        sdim,fullelementfast,linearcut)
      else
       print *,"sdim invalid"
       stop
      endif
    
     else
      call cell_intersection_grid( &
         bfact,dxgrid,xsten0,nhalf, &
         phinode, &
         cum_volume,cum_centroid,cum_area, &
         volall,cenall,sdim)

     endif

     local_scale=max(cum_volume,volslow)

     if (abs(cum_volume-volslow).gt.EPS_8_4*local_scale) then
      print *,"volume incorrect"
      print *,"cum_volume,volslow ",cum_volume,volslow
      print *,"nodedomain ",nodedomain
      do inode=1,nodedomain
       print *,"inode,phi ",inode,phinode(inode)
       do dir=1,sdim
        if (sdim.eq.2) then
         print *,"inode,dir,xnode2d ",inode,dir,xnode2d(inode,dir)
        else
         print *,"inode,dir,xnode3d ",inode,dir,xnode3d(inode,dir)
        endif
       enddo 
      enddo
      stop
     else if (abs(cum_volume-volslow).le.EPS_8_4*local_scale) then
      ! do nothing
     else
      print *,"cum_volume corrupt: ",cum_volume,volslow
      stop
     endif

     local_scale=max(cum_area,areaslow)

     if (abs(cum_area-areaslow).gt.EPS_8_4*local_scale) then
      print *,"area incorrect"
      print *,"cum_area,areaslow ",cum_area,areaslow
      print *,"nodedomain ",nodedomain
      do inode=1,nodedomain
       print *,"inode,phi ",inode,phinode(inode)
       do dir=1,sdim
        if (sdim.eq.2) then
         print *,"inode,dir,xnode2d ",inode,dir,xnode2d(inode,dir)
        else
         print *,"inode,dir,xnode3d ",inode,dir,xnode3d(inode,dir)
        endif
       enddo 
      enddo
      stop
     else if (abs(cum_area-areaslow).le.EPS_8_4*local_scale) then
      ! do nothing
     else
      print *,"cum_area corrupt: ",cum_area,areaslow
      stop
     endif

     do dir=1,sdim
      local_scale=max(abs(cum_centroid(dir)),abs(censlow(dir)))
      if (abs(cum_centroid(dir)-censlow(dir)).gt.EPS_8_4*local_scale) then
       print *,"centroid incorrect"
       stop
      else if (abs(cum_centroid(dir)-censlow(dir)).le. &
               EPS_8_4*local_scale) then
       ! do nothing
      else
       print *,"cum_centroid corrupt"
       stop
      endif
     enddo !dir=1,sdim

    enddo ! a2=0,Nangle
   enddo ! a1=0,Nangle
  enddo ! i_int=0,Nangle
   
 enddo ! sdim

 call cpu_time(t2)
 print *,"elapsed time for sanity check: ",t2-t1

return
end subroutine volume_sanity_check


! find line segment that make up the intersection of the region where
! phi>0 and the given segment
      subroutine list_segment(phi,x,xsegmentlist,nseg,sdim)
      IMPLICIT NONE

      integer, INTENT(out) :: nseg
      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: phi(2)
      real(amrex_real), INTENT(in) :: x(2)
      real(amrex_real), INTENT(out) :: xsegmentlist(2)

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif

      if ((phi(1).ge.zero).and.(phi(2).ge.zero)) then
       nseg=1
       xsegmentlist(1)=x(1) 
       xsegmentlist(2)=x(2) 
      else if ((phi(1).lt.zero).and.(phi(2).lt.zero)) then
       nseg=0
      else if ((phi(1).ge.zero).and.(phi(2).lt.zero)) then
       nseg=1
       xsegmentlist(1)=x(1) 
       xsegmentlist(2)=x(1)+phi(1)*(x(2)-x(1))/(phi(1)-phi(2))
      else if ((phi(1).lt.zero).and.(phi(2).ge.zero)) then
       nseg=1
       xsegmentlist(2)=x(2) 
       xsegmentlist(1)=x(2)-phi(2)*(x(2)-x(1))/(phi(2)-phi(1))
      else
       print *,"phi invalid in list segment"
       stop
      endif

      return
      end subroutine list_segment

      subroutine list_tets(phi,x,xtetlist,nlist_alloc, &
                      nlist,xarealist,narea,sdim)
      IMPLICIT NONE 
   
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: sdim 
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: xtetlist(sdim+1,sdim,nlist_alloc)
      real(amrex_real), INTENT(out) :: xarealist(3,3,MAXAREA)
      integer, INTENT(out) :: nlist
      integer i_tet_node
      integer j_dir
      integer, INTENT(out) :: narea

      if (sdim.ne.3) then
       print *,"sdim invalid list_tets"
       stop
      endif

      if (nlist_alloc.ge.1) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      nlist=0
      narea=0
      if ((phi(1).le.zero).and.(phi(2).le.zero).and. &
          (phi(3).le.zero).and.(phi(4).le.zero)) then
       nlist=0
      else if ((phi(1).ge.zero).and.(phi(2).ge.zero).and. &
               (phi(3).ge.zero).and.(phi(4).ge.zero)) then
       nlist=1
       do i_tet_node=1,4
       do j_dir=1,3
        xtetlist(i_tet_node,j_dir,1)=x(i_tet_node,j_dir)
       enddo 
       enddo 
      else if ((phi(1).lt.zero).and.(phi(2).ge.zero).and. &
               (phi(3).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_list_tet(phi,x,1,2,3,4,1,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(2).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(3).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_list_tet(phi,x,2,1,3,4,1,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(3).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(2).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_list_tet(phi,x,3,1,2,4,1,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(4).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(2).ge.zero).and.(phi(3).ge.zero)) then
       call shrink_list_tet(phi,x,4,1,2,3,1,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(2).lt.zero).and. &
               (phi(3).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_list_tet(phi,x,1,2,3,4,0,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(2).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(3).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_list_tet(phi,x,2,1,3,4,0,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(3).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(2).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_list_tet(phi,x,3,1,2,4,0,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(4).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(2).lt.zero).and.(phi(3).lt.zero)) then
       call shrink_list_tet(phi,x,4,1,2,3,0,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(1).lt.zero).and.(phi(2).lt.zero).and. &
               (phi(3).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_gableroof_list(phi,x,3,4,1,2,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(1).lt.zero).and.(phi(3).lt.zero).and. &
               (phi(2).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_gableroof_list(phi,x,2,4,1,3,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(2).ge.zero).and. &
               (phi(3).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_gableroof_list(phi,x,1,2,3,4,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(3).ge.zero).and. &
               (phi(2).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_gableroof_list(phi,x,1,3,2,4,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(3).ge.zero).and.(phi(2).ge.zero).and. &
               (phi(1).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_gableroof_list(phi,x,3,2,1,4,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(4).ge.zero).and. &
               (phi(2).lt.zero).and.(phi(3).lt.zero)) then
       call shrink_gableroof_list(phi,x,1,4,2,3,xtetlist, &
        nlist_alloc,nlist, &
        xarealist,narea)
      else
       print *,"bust list_tets"
       print *,"phi : ",phi(1),phi(2),phi(3),phi(4)
       stop
      endif

      return
      end subroutine list_tets





! find the triangles that make up the intersection of two triangles.
! xtetlist_old is a scratch variable
      subroutine intersect_tri(x1,x2,xtetlist_old,xtetlist, &
                      nlist_alloc,nlist,nmax,sdim)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: sdim
      integer, INTENT(out) :: nlist
      integer, INTENT(in) :: nmax
      integer nlist_old
      integer i_tet_node
      integer j_dir
      integer iplane,n,nsub,narea,n2
      real(amrex_real), INTENT(in) :: x1(sdim+1,sdim)
      real(amrex_real) x1old(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: x2(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: xtetlist_old(4,3,nlist_alloc)
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real) xsublist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(2,2,MAXAREA)
      real(amrex_real) phi1(sdim+1)
      real(amrex_real) coeff(sdim)
      real(amrex_real) voltest
      real(amrex_real) centroidtest(sdim)
      real(amrex_real) xcandidate(sdim+1,sdim)

      integer itan(sdim)
      real(amrex_real) sign,mag,maxside,testside1,testside2
      
      if (sdim.ne.2) then
       print *,"sdim bust intersect tri"
       stop
      endif

      if (nlist_alloc.ge.1) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      maxside=zero
      do i_tet_node=2,sdim+1
       testside1=zero
       testside2=zero
       do j_dir=1,sdim
        testside1=testside1+(x1(1,j_dir)-x1(i_tet_node,j_dir))**2
        testside2=testside2+(x2(1,j_dir)-x2(i_tet_node,j_dir))**2
       enddo
       testside1=sqrt(testside1)
       testside2=sqrt(testside2)
       if (testside1.gt.maxside) then
        maxside=testside1
       endif
       if (testside2.gt.maxside) then
        maxside=testside2
       endif
      enddo ! i

      if (maxside.gt.one) then
       maxside=one
      endif
      if (maxside.le.zero) then
       print *,"maxside invalid"
       stop
      endif


      nlist_old=1
      do i_tet_node=1,sdim+1
      do j_dir=1,sdim
       xtetlist_old(i_tet_node,j_dir,1)=x1(i_tet_node,j_dir)
      enddo
      enddo

! intersect the members of xtetlist_old with the planes of x2
! phi = s(n dot (x-x0)) where x0 is a point on the plane.
! s=n dot (xp-x0) where xp is not on the plane.

      do iplane=1,sdim+1

       nlist=0
       if (iplane.eq.1) then
        itan(1)=2
        itan(2)=3
       else if (iplane.eq.2) then
        itan(1)=1
        itan(2)=3
       else if (iplane.eq.3) then
        itan(1)=1
        itan(2)=2
       else
        print *,"iplane invalid"
        stop
       endif
 
       coeff(1)=-(x2(itan(2),2)-x2(itan(1),2))
       coeff(2)=(x2(itan(2),1)-x2(itan(1),1))

       mag=sqrt(coeff(1)**2+coeff(2)**2)

       if (mag.gt.EPS_14_7*maxside) then
        do j_dir=1,sdim
         coeff(j_dir)=coeff(j_dir)/mag
        enddo
        sign=zero
        do j_dir=1,sdim
         sign=sign+coeff(j_dir)*(x2(iplane,j_dir)-x2(itan(1),j_dir))
        enddo
       
        if (abs(sign).gt.EPS_14_7*maxside) then

         do n=1,nlist_old

          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           x1old(i_tet_node,j_dir)=xtetlist_old(i_tet_node,j_dir,n)
          enddo
          enddo

          do i_tet_node=1,sdim+1
           phi1(i_tet_node)=zero
           do j_dir=1,sdim
            phi1(i_tet_node)=phi1(i_tet_node)+ &
              coeff(j_dir)*(x1old(i_tet_node,j_dir)-x2(itan(1),j_dir))
           enddo
           phi1(i_tet_node)=phi1(i_tet_node)*sign
          enddo
          call list_tris(phi1,x1old,xsublist,MAXTET,nsub,xarealist,narea,sdim)
          do n2=1,nsub

           do i_tet_node=1,sdim+1
           do j_dir=1,sdim
            xcandidate(i_tet_node,j_dir)=xsublist(i_tet_node,j_dir,n2)
           enddo
           enddo
           call tetrahedron_volume(xcandidate,voltest, &
            centroidtest,sdim)

           if (voltest.gt.zero) then
            nlist=nlist+1
            if (nlist.gt.nmax) then
             print *,"nlist overflow in intersect tri"
             print *,"nlist,nmax ",nlist,nmax
             stop
            endif
            do i_tet_node=1,sdim+1
            do j_dir=1,sdim
             xtetlist(i_tet_node,j_dir,nlist)=xcandidate(i_tet_node,j_dir)
            enddo
            enddo
           endif

          enddo  ! n2
         enddo  ! n
         nlist_old=nlist
         do n=1,nlist
          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           xtetlist_old(i_tet_node,j_dir,n)=xtetlist(i_tet_node,j_dir,n)
          enddo
          enddo
         enddo
        else
         print *,"sign=0 in intersect_tri"
         stop
        endif 
       else
        print *,"mag=0 in intersect_tri: ",mag
        stop
       endif
      enddo  ! iplane

      return
      end subroutine intersect_tri

! find the tetrahedra that make up the intersection of two 
! tetrahedrons.  xtetlist_old is a scratch variable.
      subroutine intersect_tet(x1,x2,xtetlist_old,xtetlist, &
                      nlist_alloc,nlist,nmax,sdim)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: sdim
      integer, INTENT(out) :: nlist
      integer, INTENT(in) :: nmax
      integer nlist_old,iplane,n,nsub,narea,n2
      integer i_tet_node
      integer i_tan_idx
      integer j_dir
      real(amrex_real), INTENT(in) :: x1(sdim+1,sdim)
      real(amrex_real) x1old(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: x2(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: xtetlist_old(4,3,nlist_alloc)
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real) xsublist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(3,3,MAXAREA)
      real(amrex_real) phi1(sdim+1)
      real(amrex_real) coeff(sdim)
      real(amrex_real) voltest
      real(amrex_real) centroidtest(sdim)
      real(amrex_real) xcandidate(sdim+1,sdim)

      integer itan(sdim)
      real(amrex_real) sign,mag
      real(amrex_real) vec(2,sdim)
      real(amrex_real) maxside,testside1,testside2

      if (sdim.ne.3) then
       print *,"sdim bust intersect tet"
       stop
      endif

      if (nlist_alloc.ge.1) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      maxside=zero
      do i_tet_node=2,sdim+1
       testside1=zero
       testside2=zero
       do j_dir=1,sdim
        testside1=testside1+(x1(1,j_dir)-x1(i_tet_node,j_dir))**2
        testside2=testside2+(x2(1,j_dir)-x2(i_tet_node,j_dir))**2
       enddo
       testside1=sqrt(testside1)
       testside2=sqrt(testside2)
       if (testside1.gt.maxside) then
        maxside=testside1
       endif
       if (testside2.gt.maxside) then
        maxside=testside2
       endif
      enddo ! i

      if (maxside.gt.one) then
       maxside=one
      endif
      if (maxside.le.zero) then
       print *,"maxside invalid"
       stop
      endif


      nlist_old=1
      do i_tet_node=1,sdim+1
      do j_dir=1,sdim
       xtetlist_old(i_tet_node,j_dir,1)=x1(i_tet_node,j_dir)
      enddo
      enddo

! intersect the members of xtetlist_old with the planes of x2
! phi = s(n dot (x-x0)) where x0 is a point on the plane.
! s=n dot (xp-x0) where xp is not on the plane.

      do iplane=1,sdim+1

       nlist=0
       if (iplane.eq.1) then
        itan(1)=2
        itan(2)=3
        itan(3)=4
       else if (iplane.eq.2) then
        itan(1)=1
        itan(2)=3
        itan(3)=4
       else if (iplane.eq.3) then
        itan(1)=1
        itan(2)=2
        itan(3)=4
       else if (iplane.eq.4) then
        itan(1)=1
        itan(2)=2
        itan(3)=3
       else
        print *,"iplane invalid"
        stop
       endif

       do i_tan_idx=1,2
        if ((iplane.eq.itan(i_tan_idx)).or.(iplane.eq.itan(3))) then
         print *,"bust intersect_tet"
         print *,"iplane=",iplane
         print *,"i_tan_idx=",i_tan_idx
         print *,"itan(i_tan_idx)=",itan(i_tan_idx)
         print *,"itan(3)=",itan(3)
         print *,"sdim=",sdim
         stop
        endif
        do j_dir=1,sdim
         vec(i_tan_idx,j_dir)=x2(itan(i_tan_idx),j_dir)-x2(itan(3),j_dir)
        enddo
       enddo
       coeff(1)=vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2)
       coeff(2)=vec(1,3)*vec(2,1)-vec(1,1)*vec(2,3)
       coeff(3)=vec(1,1)*vec(2,2)-vec(2,1)*vec(1,2)
       mag=sqrt(coeff(1)**2+coeff(2)**2+coeff(3)**2)

       if (mag.gt.EPS_14_7*maxside) then
        do j_dir=1,sdim
         coeff(j_dir)=coeff(j_dir)/mag
        enddo
        sign=zero
        do j_dir=1,sdim
         sign=sign+coeff(j_dir)*(x2(iplane,j_dir)-x2(itan(1),j_dir))
        enddo
       
        if (abs(sign).gt.EPS_14_7*maxside) then

         do n=1,nlist_old

          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           x1old(i_tet_node,j_dir)=xtetlist_old(i_tet_node,j_dir,n)
          enddo
          enddo

          do i_tet_node=1,sdim+1
           phi1(i_tet_node)=zero
           do j_dir=1,sdim
            phi1(i_tet_node)=phi1(i_tet_node)+ &
             coeff(j_dir)*(x1old(i_tet_node,j_dir)-x2(itan(1),j_dir))
           enddo
           phi1(i_tet_node)=phi1(i_tet_node)*sign
          enddo
          call list_tets(phi1,x1old,xsublist,MAXTET,nsub,xarealist,narea,sdim)
          do n2=1,nsub

           do i_tet_node=1,sdim+1
           do j_dir=1,sdim
            xcandidate(i_tet_node,j_dir)=xsublist(i_tet_node,j_dir,n2)
           enddo
           enddo
           call tetrahedron_volume(xcandidate,voltest, &
            centroidtest,sdim)

           if (voltest.gt.zero) then
            nlist=nlist+1
            if (nlist.gt.nmax) then
             print *,"nlist overflow in intersect_tet"
             print *,"nlist,nmax ",nlist,nmax
             stop
            endif
            do i_tet_node=1,sdim+1
            do j_dir=1,sdim
             xtetlist(i_tet_node,j_dir,nlist)=xcandidate(i_tet_node,j_dir)
            enddo
            enddo
           endif

          enddo  ! n2
         enddo  ! n
         nlist_old=nlist
         do n=1,nlist
          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           xtetlist_old(i_tet_node,j_dir,n)=xtetlist(i_tet_node,j_dir,n)
          enddo
          enddo
         enddo
        else
         print *,"sign=0 in intersect_tet"
         stop
        endif 
       else
        print *,"mag=0 in intersect_tet: ",mag
        stop
       endif
      enddo  ! iplane

      return
      end subroutine intersect_tet

! find the tets of intersection of a cube with a tet.
! xtetlist_old is a scratch variable.
      subroutine intersect_cube( &
        x1, &
        xsten,nhalf, &
        xtetlist_old, &
        xtetlist, &
        nlist_alloc, &
        nlist,nmax,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: sdim,nhalf
      integer, INTENT(out) :: nlist
      integer, INTENT(in) :: nmax
      integer nlist_old,iplane,n,nsub,narea,n2
      integer i_tet_node
      integer j_dir
      real(amrex_real), INTENT(in) :: x1(sdim+1,sdim)
      real(amrex_real) x1old(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real) x0(sdim)
      real(amrex_real), INTENT(out) :: xtetlist_old(4,3,nlist_alloc)
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real) xsublist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(sdim,sdim,MAXAREA)
      real(amrex_real) phi1(sdim+1)
      real(amrex_real) coeff(sdim)
      real(amrex_real) voltest
      real(amrex_real) centroidtest(sdim)
      real(amrex_real) xcandidate(sdim+1,sdim)

      if (nhalf.lt.1) then
       print *,"nhalf invalid intersect cube"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim bust intersect cube"
       stop
      endif
      if (nlist_alloc.ge.1) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      nlist_old=1
      do i_tet_node=1,sdim+1
      do j_dir=1,sdim
       xtetlist_old(i_tet_node,j_dir,1)=x1(i_tet_node,j_dir)
      enddo
      enddo

! intersect the members of xtetlist_old with the planes of x2
! phi = (n dot (x-x0)) where x0 is a point on the plane.

      do iplane=1,2*sdim

       nlist=0

       do j_dir=1,sdim
        coeff(j_dir)=zero
        x0(j_dir)=xsten(-1,j_dir)
       enddo
 
       if (iplane.eq.1) then  ! left side
        coeff(1)=one
        x0(1)=xsten(-1,1)
       else if (iplane.eq.2) then  ! right side
        coeff(1)=-one
        x0(1)=xsten(1,1)
       else if (iplane.eq.3) then  ! front
        coeff(2)=one
        x0(2)=xsten(-1,2)
       else if (iplane.eq.4) then  ! back
        coeff(2)=-one
        x0(2)=xsten(1,2)
       else if (iplane.eq.5) then  ! bottom
        if (sdim.ne.3) then
         print *,"dimension bust"
         stop
        endif
        coeff(sdim)=one
        x0(sdim)=xsten(-1,sdim)
       else if (iplane.eq.6) then  ! top
        if (sdim.ne.3) then
         print *,"dimension bust"
         stop
        endif
        coeff(sdim)=-one
        x0(sdim)=xsten(1,sdim)
       else
        print *,"iplane invalid"
        stop
       endif

       do n=1,nlist_old

        do i_tet_node=1,sdim+1
        do j_dir=1,sdim
         if ((n.ge.1).and.(n.le.nlist_alloc)) then
          x1old(i_tet_node,j_dir)=xtetlist_old(i_tet_node,j_dir,n)
         else
          print *,"n out of range1"
          stop
         endif
        enddo
        enddo

        do i_tet_node=1,sdim+1
         phi1(i_tet_node)=zero
         do j_dir=1,sdim
          phi1(i_tet_node)=phi1(i_tet_node)+ &
            coeff(j_dir)*(x1old(i_tet_node,j_dir)-x0(j_dir))
         enddo
        enddo
     
         ! in: intersect_cube 
        if (sdim.eq.3) then
         call list_tets(phi1,x1old,xsublist,MAXTET,nsub,xarealist,narea,sdim)
        else if (sdim.eq.2) then
         call list_tris(phi1,x1old,xsublist,MAXTET,nsub,xarealist,narea,sdim)
        else
         print *,"sdim invalid"
         stop
        endif


        do n2=1,nsub

         do i_tet_node=1,sdim+1
         do j_dir=1,sdim
          xcandidate(i_tet_node,j_dir)=xsublist(i_tet_node,j_dir,n2)
         enddo
         enddo
         call tetrahedron_volume(xcandidate,voltest, &
          centroidtest,sdim)

         if (voltest.gt.zero) then
          nlist=nlist+1
          if (nlist.gt.nmax) then
           print *,"nlist overflow in intersect_cube"
           print *,"nlist,nmax ",nlist,nmax
           stop
          endif
          if ((nlist.ge.1).and.(nlist.le.nlist_alloc)) then
           do i_tet_node=1,sdim+1
           do j_dir=1,sdim
            xtetlist(i_tet_node,j_dir,nlist)=xcandidate(i_tet_node,j_dir)
           enddo
           enddo
          else
           print *,"nlist out of range2"
           stop
          endif
         endif

        enddo  ! n2
       enddo  ! n
       nlist_old=nlist
       do n=1,nlist
        if ((n.ge.1).and.(n.le.nlist_alloc)) then
         do i_tet_node=1,sdim+1
         do j_dir=1,sdim
          xtetlist_old(i_tet_node,j_dir,n)=xtetlist(i_tet_node,j_dir,n)
         enddo
         enddo
        else
         print *,"n out of range3"
         stop
        endif
       enddo
      enddo  ! iplane

      return
      end subroutine intersect_cube



! get volume/centroid of a triangulated region
      subroutine get_cut_geom3D( &
        xtetlist, & ! intent(in)
        nlist_alloc, &
        nlist, &
        nmax, &
        volcut, &
        cencut,sdim)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: nlist,nmax,sdim
      integer i_tet_node
      integer j_dir
      integer n
      real(amrex_real), INTENT(out) :: volcut
      real(amrex_real), INTENT(out) :: cencut(sdim)
      real(amrex_real), INTENT(in) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real) volumelist
      real(amrex_real) centroidlist(sdim)
      real(amrex_real) xint(sdim+1,sdim)

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid get_cut_geom3D"
       stop
      endif 
      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      volcut=zero
      do j_dir=1,sdim
       cencut(j_dir)=zero
      enddo
      do n=1,nlist
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=xtetlist(i_tet_node,j_dir,n)
       enddo
       enddo
     
       call tetrahedron_volume(xint,volumelist,centroidlist,sdim)

       volcut=volcut+volumelist
       do j_dir=1,sdim
        cencut(j_dir)=cencut(j_dir)+centroidlist(j_dir)*volumelist
       enddo
      enddo
      if (volcut.gt.zero) then
       do j_dir=1,sdim
        cencut(j_dir)=cencut(j_dir)/volcut
       enddo
      else
       do j_dir=1,sdim
        cencut(j_dir)=zero
       enddo
      endif

      return
      end subroutine get_cut_geom3D

      integer function get_material_rank(im,order_algorithm_in)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, intent(in) :: im
      integer, intent(in) :: order_algorithm_in
      integer :: prefix

      if ((im.ge.1).and.(im.le.num_materials)) then
       !do nothing
      else
       print *,"im invalid ",im
       stop
      endif
      if ((order_algorithm_in.ge.1).and. &
          (order_algorithm_in.le.num_materials+1)) then
       !do nothing
      else
       print *,"order_algorithm_in invalid ",order_algorithm_in
       stop
      endif
      if (is_rigid(im).eq.1) then
       prefix=0
      else if (is_elastic(im).eq.1) then
       prefix=num_materials+1
      else if ((is_rigid(im).eq.0).and.(is_elastic(im).eq.0)) then
       prefix=2*(num_materials+1)
      else
       print *,"im invalid ",im
       print *,"or is_rigid invalid ",is_rigid(im)
       print *,"or is_elastic invalid ",is_elastic(im)
       stop
      endif
      get_material_rank=prefix+order_algorithm_in

      end function get_material_rank

      integer function is_proper_layer(im,layer_iter)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, intent(in) :: im
      integer, intent(in) :: layer_iter

      if (layer_iter.eq.RIGID_LAYER_INDEX) then
       is_proper_layer=is_rigid(im)
      else if (layer_iter.eq.ELASTIC_LAYER_INDEX) then
       is_proper_layer=is_elastic(im)
      else if (layer_iter.eq.FLUID_LAYER_INDEX) then
       if ((is_rigid(im).eq.0).and. &
           (is_elastic(im).eq.0)) then
        is_proper_layer=1
       else 
        is_proper_layer=0
       endif
      else
       print *,"layer_iter invalid ",layer_iter
       stop
      endif

      end function is_proper_layer

! normal points from phi<0 to phi>0   phi=n dot (x-x0plane) + intercept
! vof, ref centroid, order,slope,intercept  x num_materials
! tetrahedra representing unfilled region
! (where LS < 0)
! routine starts with tetrahedralization of a cell,
! then intersects this region with the plane (LS<0 side)
! of already initialized materials.
!
! tets_box_planes is called from 
!  individual_MOF,
!  multi_get_volume_grid,
!  multi_get_volume_grid_simple,
!  multi_get_volume_grid_and_map,
!  multi_get_area_pairs
      subroutine tets_box_planes( &
       layer_flag, &
       bfact,dx, &
       xsten0,nhalf0, &!tet domain,tet LS x0=centroid(tet domain),hex LS x0
       xsten_box,nhalf_box, &
       mofdata, &
       xtetlist, &
       nlist_alloc, &
       nlist, &
       nmax, &
       sdim)

      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: layer_flag
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: sdim,bfact
      integer, INTENT(in) :: nhalf0
      integer, INTENT(in) :: nhalf_box
      integer symmetry_flag,ntetbox
      integer, INTENT(out) :: nlist
      integer, INTENT(in) :: nmax
      integer nlist_old,iplane,n,nsub,narea,n2
      integer j_dir
      integer i_tet_node
      integer i_grid_node
      integer j_grid_node
      integer k_grid_node
      integer id,icrit,im,vofcomp,iorder
      integer dir
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: xsten_box(-nhalf_box:nhalf_box,sdim)
      real(amrex_real) :: xsten0_LS(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: mofdata(num_materials*(2*sdim+3))
      real(amrex_real) xtetlist_old(4,3,nlist_alloc)
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real) xsublist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(sdim,sdim,MAXAREA)
      real(amrex_real) phi1(sdim+1),dummyphi(sdim+1)
      real(amrex_real) voltest
      real(amrex_real) centroidtest(sdim)
      real(amrex_real) xcandidate(sdim+1,sdim)
      real(amrex_real) xx(sdim+1,sdim)
      real(amrex_real) x1old(sdim+1,sdim)
      real(amrex_real) xnode(4*(sdim-1),sdim)
      real(amrex_real) phinode(4*(sdim-1))
      integer inode
      integer layer_iter
      real(amrex_real) nn(sdim)
      real(amrex_real) intercept
      integer is_masked(num_materials)

      call init_mask_flag(layer_flag,is_masked)

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid tets box planes"
       stop
      endif
      if (nhalf_box.lt.1) then
       print *,"nhalf_box invalid tets box planes"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid tets_box_planes"
       stop
      endif
      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid tets box planes"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid121"
       stop
      endif

      do i_tet_node=-nhalf0,nhalf0
      do j_dir=1,sdim
       xsten0_LS(i_tet_node,j_dir)=xsten0(i_tet_node,j_dir)
      enddo
      enddo

      symmetry_flag=0
      call get_ntetbox(ntetbox,symmetry_flag,sdim)

      inode=1

      if (sdim.eq.3) then
       do k_grid_node=-1,1,2
       do j_grid_node=-1,1,2
       do i_grid_node=-1,1,2 
        xnode(inode,1)=xsten_box(i_grid_node,1)
        xnode(inode,2)=xsten_box(j_grid_node,2)
        xnode(inode,sdim)=xsten_box(k_grid_node,sdim)
        phinode(inode)=one
        inode=inode+1
       enddo
       enddo
       enddo
      else if (sdim.eq.2) then
       do j_grid_node=-1,1,2
       do i_grid_node=-1,1,2 
        xnode(inode,1)=xsten_box(i_grid_node,1)
        xnode(inode,2)=xsten_box(j_grid_node,2)
        phinode(inode)=one
        inode=inode+1
       enddo
       enddo
      else
       print *,"sdim invalid"
       stop
      endif

      do id=1,ntetbox
       call extract_tet(xnode,phinode,xx,dummyphi,id,symmetry_flag,sdim)
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xtetlist_old(i_tet_node,j_dir,id)=xx(i_tet_node,j_dir)
        xtetlist(i_tet_node,j_dir,id)=xx(i_tet_node,j_dir)
       enddo 
       enddo 
      enddo !id=1,ntetbox

      nlist_old=ntetbox
      nlist=ntetbox
    
      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       do iplane=1,num_materials

        icrit=0

        do im=1,num_materials
         vofcomp=(im-1)*(2*sdim+3)+1
         if ((is_masked(im).eq.0).and. &
             (is_proper_layer(im,layer_iter).eq.1)) then
          iorder=NINT(mofdata(vofcomp+sdim+1))
          if ((iorder.ge.0).and.(iorder.le.num_materials)) then
           if (iorder.eq.iplane) then
            icrit=im
           endif
          else
           print *,"iorder invalid ",iorder
           stop
          endif
         else if ((is_masked(im).eq.1).or. &
                  (is_proper_layer(im,layer_iter).eq.0)) then
          ! do nothing
         else
          print *,"is_masked invalid: ",is_masked
          print *,"or is_proper_layer invalid: ", &
               im,layer_iter,is_proper_layer(im,layer_iter)
          stop
         endif
        enddo ! im=1..num_materials

        if ((icrit.gt.0).and.(icrit.le.num_materials)) then
         nlist=0
         vofcomp=(icrit-1)*(2*sdim+3)+1
          ! want intersection where phi_original<0 but routine
          ! checks if phi_input>0.  So negate phi_original.
         do dir=1,sdim
          nn(dir)=-mofdata(vofcomp+sdim+1+dir)  
         enddo
         intercept=-mofdata(vofcomp+2*sdim+2)
 
         do n=1,nlist_old

          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           x1old(i_tet_node,j_dir)=xtetlist_old(i_tet_node,j_dir,n)
          enddo
          enddo

            ! LEVELSET FOR CUTTING PLANE
          do i_tet_node=1,sdim+1
           phi1(i_tet_node)=intercept
           do j_dir=1,sdim
            phi1(i_tet_node)=phi1(i_tet_node)+ &
              nn(j_dir)*(x1old(i_tet_node,j_dir)-xsten0_LS(0,j_dir))
           enddo
          enddo
           ! tetrahedras representing intersection of region where phi1>0
           ! (phi_original<0) and the tetrahedra x1old
          if (sdim.eq.3) then
           call list_tets(phi1,x1old,xsublist,MAXTET,nsub,xarealist,narea,sdim)
          else if (sdim.eq.2) then
           call list_tris(phi1,x1old,xsublist,MAXTET,nsub,xarealist,narea,sdim)
          else
           print *,"sdim invalid"
           stop
          endif

          do n2=1,nsub

           do i_tet_node=1,sdim+1
           do j_dir=1,sdim
            xcandidate(i_tet_node,j_dir)=xsublist(i_tet_node,j_dir,n2)
           enddo
           enddo
           call tetrahedron_volume(xcandidate,voltest, &
            centroidtest,sdim)

           if (voltest.gt.zero) then
            nlist=nlist+1
            if (nlist.gt.nmax) then
             print *,"nlist overflow in tets_box_planes"
             print *,"nlist,nmax ",nlist,nmax
             stop
            endif
            do i_tet_node=1,sdim+1
            do j_dir=1,sdim
             xtetlist(i_tet_node,j_dir,nlist)=xcandidate(i_tet_node,j_dir)
            enddo
            enddo
           endif

          enddo  ! n2
         enddo  ! n
         nlist_old=nlist
         do n=1,nlist
          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           xtetlist_old(i_tet_node,j_dir,n)=xtetlist(i_tet_node,j_dir,n)
          enddo
          enddo
         enddo
        else if (icrit.eq.0) then
         nlist=nlist_old
         do n=1,nlist
          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           xtetlist(i_tet_node,j_dir,n)=xtetlist_old(i_tet_node,j_dir,n)
          enddo
          enddo
         enddo
        else
         print *,"icrit invalid"
         stop
        endif 
       enddo ! iplane=1,num_materials

      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      return
      end subroutine tets_box_planes

      subroutine init_mask_flag( &
        layer_flag, &
        is_masked)
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      integer, intent(in) :: layer_flag
      integer, intent(out) :: is_masked(num_materials)
      integer im

      do im=1,num_materials

       if (layer_flag.eq.FLUIDS_ELASTIC_RIGID_LAYER) then
        is_masked(im)=0
       else if (layer_flag.eq.NULL_LAYER) then
        is_masked(im)=1
       else if (layer_flag.eq.FLUIDS_RIGID_LAYER) then
        is_masked(im)=is_elastic(im)
       else if (layer_flag.eq.FLUIDS_ELASTIC_LAYER) then
        is_masked(im)=is_rigid(im)
       else if (layer_flag.eq.ELASTIC_RIGID_LAYER) then
        if ((is_rigid(im).eq.1).or. &
            (is_elastic(im).eq.1)) then
         is_masked(im)=0
        else if ((is_rigid(im).eq.0).and. &
                 (is_elastic(im).eq.0)) then
         is_masked(im)=1
        else
         print *,"is_rigid or is_elastic invalid: ", &
           im,is_rigid(im),is_elastic(im)
         stop
        endif
       else if (layer_flag.eq.RIGID_LAYER) then
        is_masked(im)=1-is_rigid(im)
       else if (layer_flag.eq.ELASTIC_LAYER) then
        is_masked(im)=1-is_elastic(im)
       else if (layer_flag.eq.FLUIDS_LAYER) then
        if ((is_rigid(im).eq.1).or. &
            (is_elastic(im).eq.1)) then
         is_masked(im)=1
        else if ((is_rigid(im).eq.0).and. &
                 (is_elastic(im).eq.0)) then
         is_masked(im)=0
        else
         print *,"is_rigid or is_elastic invalid ",im, &
           is_rigid(im),is_elastic(im)
         stop
        endif
       else
        print *,"layer_flag invalid multi_get_volume_grid ",layer_flag
        stop
       endif

      enddo !im=1,num_materials

      return
      end subroutine init_mask_flag

      subroutine tets_tet_planes( &
        layer_flag, &
        bfact,dx, &
        xsten0,nhalf0, &
        xtet,mofdata, &
        xtetlist, &
        nlist_alloc, &
        nlist, &
        nmax, &
        sdim)

      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: layer_flag
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: sdim,bfact,nhalf0
      integer, INTENT(out) :: nlist
      integer, INTENT(in) :: nmax
      integer nlist_old,iplane,n,nsub,narea,n2
      integer i_tet_node
      integer j_dir
      integer icrit,im,vofcomp,iorder
      integer dir
      real(amrex_real), INTENT(in) :: xtet(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: mofdata(num_materials*(2*sdim+3))
      real(amrex_real) xtetlist_old(4,3,nlist_alloc)
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real) xsublist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(sdim,sdim,MAXAREA)
      real(amrex_real) phi1(sdim+1)
      real(amrex_real) voltest
      real(amrex_real) centroidtest(sdim)
      real(amrex_real) xcandidate(sdim+1,sdim)
      real(amrex_real) x1old(sdim+1,sdim)
      integer layer_iter
      real(amrex_real) nn(sdim)
      real(amrex_real) intercept

      integer is_masked(num_materials)

      call init_mask_flag(layer_flag,is_masked)

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid124"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid tets_tet_planes"
       stop
      endif
      if ((num_materials.lt.1).or.(num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid tets tet planes"
       stop
      endif

      nlist_old=1
      nlist=1
      do i_tet_node=1,sdim+1
      do j_dir=1,sdim
        xtetlist_old(i_tet_node,j_dir,1)=xtet(i_tet_node,j_dir)
        xtetlist(i_tet_node,j_dir,1)=xtet(i_tet_node,j_dir)
      enddo 
      enddo 

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       do iplane=1,num_materials

        icrit=0

        do im=1,num_materials
         vofcomp=(im-1)*(2*sdim+3)+1
         if ((is_masked(im).eq.0).and. &
             (is_proper_layer(im,layer_iter).eq.1)) then
          iorder=NINT(mofdata(vofcomp+sdim+1))
          if ((iorder.ge.0).and.(iorder.le.num_materials)) then
           if (iorder.eq.iplane) then
            icrit=im
           endif
          else
           print *,"iorder invalid ",iorder
           stop
          endif
         else if ((is_masked(im).eq.1).or. &
                  (is_proper_layer(im,layer_iter).eq.0)) then
          ! do nothing
         else
          print *,"is_masked invalid: ",is_masked
          print *,"or is_proper_layer invalid: ", &
               im,layer_iter,is_proper_layer(im,layer_iter)
          stop
         endif
        enddo ! im=1..num_materials

        if ((icrit.gt.0).and.(icrit.le.num_materials)) then

         nlist=0
         vofcomp=(icrit-1)*(2*sdim+3)+1
         ! want intersection where phi_original<0 but routine
         ! checks if phi_input>0.  So negate phi_original.
         do dir=1,sdim
          nn(dir)=-mofdata(vofcomp+sdim+1+dir)  
         enddo
         intercept=-mofdata(vofcomp+2*sdim+2)
 
         do n=1,nlist_old

          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           x1old(i_tet_node,j_dir)=xtetlist_old(i_tet_node,j_dir,n)
          enddo
          enddo

          do i_tet_node=1,sdim+1
           phi1(i_tet_node)=intercept
           do j_dir=1,sdim
            phi1(i_tet_node)=phi1(i_tet_node)+ &
              nn(j_dir)*(x1old(i_tet_node,j_dir)-xsten0(0,j_dir))
           enddo
          enddo
           ! triangles representing intersection of region where phi1>0
           ! (phi_original<0) and the triangle x1old
          if (sdim.eq.3) then
           call list_tets(phi1,x1old,xsublist,MAXTET,nsub,xarealist,narea,sdim)
          else if (sdim.eq.2) then
           call list_tris(phi1,x1old,xsublist,MAXTET,nsub,xarealist,narea,sdim)
          else
           print *,"sdim invalid"
           stop
          endif


          do n2=1,nsub

           do i_tet_node=1,sdim+1
           do j_dir=1,sdim
            xcandidate(i_tet_node,j_dir)=xsublist(i_tet_node,j_dir,n2)
           enddo
           enddo

           call tetrahedron_volume(xcandidate,voltest, &
            centroidtest,sdim)

           if (voltest.gt.zero) then
            nlist=nlist+1
            if (nlist.gt.nmax) then
             print *,"nlist overflow in tets_tet_planes"
             print *,"nlist,nmax ",nlist,nmax
             stop
            endif
            do i_tet_node=1,sdim+1
            do j_dir=1,sdim
             xtetlist(i_tet_node,j_dir,nlist)=xcandidate(i_tet_node,j_dir)
            enddo
            enddo
           endif

          enddo  ! n2
         enddo  ! n
         nlist_old=nlist
         do n=1,nlist
          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           xtetlist_old(i_tet_node,j_dir,n)=xtetlist(i_tet_node,j_dir,n)
          enddo
          enddo
         enddo
        else if (icrit.eq.0) then
         nlist=nlist_old
         do n=1,nlist
          do i_tet_node=1,sdim+1
          do j_dir=1,sdim
           xtetlist(i_tet_node,j_dir,n)=xtetlist_old(i_tet_node,j_dir,n)
          enddo
          enddo
         enddo
        else
         print *,"icrit invalid"
         stop
        endif 
       enddo ! iplane=1,num_materials

      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      return
      end subroutine tets_tet_planes
 

! find the volume/area/etc. of intersection of a plane with
! a tetrahedron
      subroutine intersection_volumeXYZ(phi,x,volumedark, &
       centroiddark,area,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real) xtetlist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(sdim,sdim,MAXAREA)
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real) xint(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: centroiddark(sdim)
      real(amrex_real) centroidlistdark(sdim)

      real(amrex_real), INTENT(out) :: volumedark,area
      real(amrex_real) volumelistdark,arealist

      integer n,nlist,narea
      integer i_tet_node
      integer j_dir

      if (sdim.ne.3) then
       print *,"sdim invalid"
       stop
      endif

      call list_tets(phi,x,xtetlist,MAXTET,nlist,xarealist,narea,sdim)

      volumedark=zero
      area=zero
      do j_dir=1,sdim
       centroiddark(j_dir)=zero
      enddo

      do n=1,nlist
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=xtetlist(i_tet_node,j_dir,n)
       enddo
       enddo
       call tetrahedron_volume(xint,volumelistdark,centroidlistdark,sdim)
       volumedark=volumedark+volumelistdark
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)+ &
                centroidlistdark(j_dir)*volumelistdark
       enddo
      enddo
      if (volumedark.gt.zero) then
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)/volumedark
       enddo
      else
       do j_dir=1,sdim
        centroiddark(j_dir)=zero
       enddo
      endif

      do n=1,narea
       do i_tet_node=1,sdim
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=xarealist(i_tet_node,j_dir,n)
       enddo
       enddo
       do j_dir=1,sdim
        xint(sdim+1,j_dir)=zero
       enddo
       call areaXYZ(xint,1,2,3,arealist)
       area=area+arealist
      enddo

      return
      end subroutine intersection_volumeXYZ


! find the volume/centroid of intersection of a plane with
! a tetrahedron
      subroutine intersection_volumeXYZ_simple(phi,x,volumedark, &
       centroiddark,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real) xtetlist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(sdim,sdim,MAXAREA)
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real) xint(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: centroiddark(sdim)
      real(amrex_real) centroidlistdark(sdim)

      real(amrex_real), INTENT(out) :: volumedark
      real(amrex_real) volumelistdark

      integer i_tet_node
      integer j_dir
      integer n,nlist,narea

      if (sdim.ne.3) then
       print *,"sdim invalid"
       stop
      endif

      call list_tets(phi,x,xtetlist,MAXTET,nlist,xarealist,narea,sdim)

      volumedark=zero
      do j_dir=1,sdim
       centroiddark(j_dir)=zero
      enddo

      do n=1,nlist
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=xtetlist(i_tet_node,j_dir,n)
       enddo
       enddo
       call tetrahedron_volume(xint,volumelistdark,centroidlistdark,sdim)
       volumedark=volumedark+volumelistdark
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)+ &
                centroidlistdark(j_dir)*volumelistdark
       enddo
      enddo
      if (volumedark.gt.zero) then
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)/volumedark
       enddo
      else
       do j_dir=1,sdim
        centroiddark(j_dir)=zero
       enddo
      endif

      return
      end subroutine intersection_volumeXYZ_simple



! find the volume/centroid of intersection of a plane with
! a tetrahedron
      subroutine intersection_volumeXYZ_and_map( &
       normdir,coeff, &
       phi,x, &
       volumedark, &
       volumedark_map, &
       centroiddark, &
       centroiddark_map, &
       sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: normdir
      real(amrex_real), INTENT(in) :: coeff(2)

      real(amrex_real) xtetlist(sdim+1,sdim,MAXTET)
      real(amrex_real) xarealist(sdim,sdim,MAXAREA)
      real(amrex_real), INTENT(in) :: phi(sdim+1)
      real(amrex_real), INTENT(in) :: x(sdim+1,sdim)
      real(amrex_real) xint(sdim+1,sdim)
      real(amrex_real), INTENT(out) :: centroiddark(sdim)
      real(amrex_real) centroidlistdark(sdim)
      real(amrex_real), INTENT(out) :: centroiddark_map(sdim)
      real(amrex_real) centroidlistdark_map(sdim)

      real(amrex_real), INTENT(out) :: volumedark
      real(amrex_real) volumelistdark
      real(amrex_real), INTENT(out) :: volumedark_map
      real(amrex_real) volumelistdark_map

      integer n,nlist,narea
      integer i_tet_node
      integer j_dir

      if (sdim.ne.3) then
       print *,"sdim invalid"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.sdim)) then
       print *,"normdir invalid"
       stop
      endif
      if (coeff(1).le.zero) then
       print *,"coeff(1) invalid"
       stop
      endif
      call list_tets(phi,x,xtetlist,MAXTET,nlist,xarealist,narea,sdim)

      volumedark=zero
      volumedark_map=zero
      do j_dir=1,sdim
       centroiddark(j_dir)=zero
       centroiddark_map(j_dir)=zero
      enddo

      do n=1,nlist
       do i_tet_node=1,sdim+1
       do j_dir=1,sdim
        xint(i_tet_node,j_dir)=xtetlist(i_tet_node,j_dir,n)
       enddo
       enddo
       call tetrahedron_volume_and_map(normdir,coeff, &
        xint,volumelistdark,centroidlistdark, &
        volumelistdark_map,centroidlistdark_map, &
        sdim)
       volumedark=volumedark+volumelistdark
       volumedark_map=volumedark_map+volumelistdark_map
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)+ &
                centroidlistdark(j_dir)*volumelistdark
        centroiddark_map(j_dir)=centroiddark_map(j_dir)+ &
         centroidlistdark_map(j_dir)*volumelistdark_map
       enddo
      enddo
      if ((volumedark.gt.zero).and.(volumedark_map.gt.zero)) then
       do j_dir=1,sdim
        centroiddark(j_dir)=centroiddark(j_dir)/volumedark
        centroiddark_map(j_dir)=centroiddark_map(j_dir)/volumedark_map
       enddo
      else if ((volumedark.eq.zero).or.(volumedark_map.eq.zero)) then
       do j_dir=1,sdim
        centroiddark(j_dir)=zero
        centroiddark_map(j_dir)=zero
       enddo
      else
       print *,"volumedark or volumedark_map invalid"
       stop
      endif

      return
      end subroutine intersection_volumeXYZ_and_map


        ! distance from line(im) to point(im_opp)
      subroutine dist_centroid_line( &
       xsten,nhalf, &
       im,im_opp, &
       mofdata, &
       cencell, &
       dist_tol, &
       sdim, &
       dist_to_line)
      
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: nhalf,im,im_opp
      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: mofdata(num_materials*ngeom_recon)
      real(amrex_real), INTENT(in) :: cencell(sdim)
      real(amrex_real) slope(sdim)
      real(amrex_real) xtarget(sdim)
      real(amrex_real), INTENT(in) :: dist_tol
      real(amrex_real) mag
      real(amrex_real), INTENT(out) :: dist_to_line
      integer ibase,ibase_opp,dir
      real(amrex_real) intercept
     

      if (nhalf.lt.1) then
       print *,"nhalf invalid"
       stop
      endif
      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon invalid"
       stop
      endif
      if ((num_materials.lt.1).or.(num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials out of range"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      if ((im.lt.1).or.(im.gt.num_materials).or. &
          (im_opp.lt.1).or.(im_opp.gt.num_materials)) then
       print *,"im or im_opp invalid"
       stop
      endif
      if (dist_tol.gt.zero) then
       !do nothing
      else
       print *,"dist_tol invalid: ",dist_tol
       stop
      endif

       ! vfrac,centroid,order,slope,intercept x num_materials
      ibase=(im-1)*ngeom_recon+1+sdim+1
      ibase_opp=(im_opp-1)*ngeom_recon+1

      mag=zero
      do dir=1,sdim
       slope(dir)=mofdata(ibase+dir)
       mag=mag+slope(dir)**2
       xtarget(dir)=mofdata(ibase_opp+dir)+cencell(dir)
      enddo
      mag=sqrt(mag)
      intercept=mofdata(ibase+sdim+1)
      if (mag.gt.zero) then
       do dir=1,sdim
        slope(dir)=slope(dir)/mag
       enddo
       intercept=intercept/mag
      else
       print *,"slope has 0 magnitude"
       stop
      endif
      dist_to_line=intercept
      do dir=1,sdim
       dist_to_line=dist_to_line+ &
         slope(dir)*(xtarget(dir)-xsten(0,dir))
      enddo

      if ((levelrz.eq.COORDSYS_CARTESIAN).or.(levelrz.eq.COORDSYS_RZ)) then
       ! do nothing
      else
       print *,"levelrz invalid ",levelrz
       stop
      endif

      dist_to_line=abs(dist_to_line)
      if (dist_to_line.lt.dist_tol) then
       dist_to_line=dist_tol
      endif

      end subroutine dist_centroid_line

 
      subroutine CISBOX( &
       xsten,nhalf, &
       xlo,dx,i,j,k, &
       bfact,level, &
       vol,cen,sdim)
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: bfact,level
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(out) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xlo(sdim)
      integer, INTENT(in) :: i,j,k
      real(amrex_real), INTENT(out) :: vol
      real(amrex_real), INTENT(out) :: cen(sdim)

      if (nhalf.lt.1) then
       print *,"nhalf invalid cisbox"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid125"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      call gridsten_level(xsten,i,j,k,level,nhalf)
      call Box_volumeFAST(bfact,dx,xsten,nhalf,vol,cen,sdim)

      return
      end subroutine CISBOX 


      subroutine CISBOXFINE( &
       xsten,nhalf, &
       xlo,dx, &
       i,j,k, &
       ifine,jfine,kfine, &
       bfact,level, &
       vol,cen,sdim)
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: bfact,level
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(out) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xlo(sdim)
      integer, INTENT(in) :: i,j,k
      integer, INTENT(in) :: ifine,jfine,kfine
      real(amrex_real), INTENT(out) :: vol
      real(amrex_real), INTENT(out) :: cen(sdim)
      integer :: dir2

      if (nhalf.lt.1) then
       print *,"nhalf invalid cisboxfine: ",nhalf
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid cisboxfine"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid cisboxfine"
       stop
      endif
      call gridsten_level(xsten,i,j,k,level,nhalf)

      dir2=1
      if (ifine.eq.0) then
       xsten(1,dir2)=xsten(0,dir2)
      else if (ifine.eq.1) then
       xsten(-1,dir2)=xsten(0,dir2)
      else
       print *,"ifine invalid"
       stop
      endif
      dir2=2
      if (jfine.eq.0) then
       xsten(1,dir2)=xsten(0,dir2)
      else if (jfine.eq.1) then
       xsten(-1,dir2)=xsten(0,dir2)
      else
       print *,"jfine invalid"
       stop
      endif
      if (SDIM.eq.3) then
       dir2=SDIM
       if (kfine.eq.0) then
        xsten(1,dir2)=xsten(0,dir2)
       else if (kfine.eq.1) then
        xsten(-1,dir2)=xsten(0,dir2)
       else
        print *,"kfine invalid: ",kfine
        stop
       endif
      endif
      do dir2=1,SDIM
       xsten(0,dir2)=half*(xsten(-1,dir2)+xsten(1,dir2))
      enddo

        ! bfact is not used.
      call Box_volumeFAST(bfact,dx,xsten,nhalf,vol,cen,sdim)

      return
      end subroutine CISBOXFINE 





      subroutine CISBOXHALF( &
       xsten,nhalf, &
       xlo,dx,i,j,k,iside,veldir, &
       bfact,level, &
       vol,cen,sdim)
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: bfact,iside,veldir,level
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(out) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xlo(sdim)
      integer, INTENT(in) :: i,j,k
      real(amrex_real), INTENT(out) :: vol
      real(amrex_real), INTENT(out) :: cen(sdim)

      if ((iside.ne.-1).and.(iside.ne.1)) then
       print *,"iside invalid"
       stop
      endif
      if ((veldir.lt.1).or.(veldir.gt.sdim)) then
       print *,"veldir invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid126"
       stop
      endif
      if (nhalf.ne.1) then
       print *,"nhalf invalid cisboxhalf"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      call CISBOX(xsten,nhalf, &
       xlo,dx,i,j,k,bfact,level,vol,cen,sdim)
     
       ! if iside=-1
       !  xstennew(-1)=xsten(-1)
       !  xstennew(0)=(xsten(0)+xsten(-1))/2
       !  xstennew(1)=xsten(0)
       ! if iside=1
       !  xstennew(-1)=xsten(0)
       !  xstennew(0)=(xsten(0)+xsten(1))/2
       !  xstennew(1)=xsten(1)
      xsten(-iside,veldir)=xsten(0,veldir)
      xsten(0,veldir)=half*(xsten(1,veldir)+xsten(-1,veldir))

      call Box_volumeFAST(bfact,dx,xsten,nhalf,vol,cen,sdim)

      return
      end subroutine CISBOXHALF 



       ! centroid is in absolute coordinate system 
      subroutine Box_volumeFAST( &
        bfact,dx, &
        xsten,nhalf, &
        volume,centroid,sdim)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: sdim,bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(out) :: volume
      real(amrex_real), INTENT(out) :: centroid(sdim)
      integer dir
      real(amrex_real) rval,dr
    
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid Box_volumeFAST: ",sdim
       stop
      endif 
      if (nhalf.lt.1) then
       print *,"nhalf invalid Box_volumeFAST: ",nhalf
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid127 Box_volumeFAST: ",bfact
       stop
      endif

      if (sdim.eq.2) then

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        volume=one
        do dir=1,sdim
         volume=volume*(xsten(1,dir)-xsten(-1,dir))
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
        enddo
       else if (levelrz.eq.COORDSYS_RZ) then
        volume=one
        do dir=1,sdim
         volume=volume*(xsten(1,dir)-xsten(-1,dir))
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
        enddo
        rval=centroid(1)
        dr=xsten(1,1)-xsten(-1,1)
        if (rval.ne.zero) then
         volume=volume*two*Pi*abs(rval)
         centroid(1)=rval+dr*dr/(12.0*rval)
        else if (rval.eq.zero) then
         volume=zero
         centroid(1)=zero
        else
         print *,"rval is NaN (Box_volumeFAST): ",rval
         stop
        endif
       else
        print *,"levelrz invalid Box_volumeFAST(1): ",levelrz
        stop
       endif

      else if (sdim.eq.3) then
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        volume=one
        do dir=1,sdim
         volume=volume*(xsten(1,dir)-xsten(-1,dir))
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
        enddo
       else
        print *,"levelrz invalid Box_volumeFAST(2): ",levelrz
        stop
       endif
      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine Box_volumeFAST


       ! centroid is in absolute coordinate system 
      subroutine Box_volumeFAST_and_map( &
        normdir, &
        coeff, &
        bfact, &
        dx, &
        xsten, &
        nhalf, &
        volume, &
        centroid, &
        volume_map, &
        centroid_map, &
        sdim)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: sdim,bfact,nhalf
      integer, INTENT(in) :: normdir
      real(amrex_real), INTENT(in) :: coeff(2)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(out) :: volume
      real(amrex_real), INTENT(out) :: volume_map
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real), INTENT(out) :: centroid_map(sdim)
      integer dir
      real(amrex_real) rval,dr
    
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid Box_volumeFAST"
       stop
      endif 
      if (nhalf.lt.1) then
       print *,"nhalf invalid boxvolumefast"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.sdim)) then
       print *,"normdir invalid"
       stop
      endif
      if (coeff(1).le.zero) then
       print *,"coeff(1) invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid128"
       stop
      endif

      if (sdim.eq.2) then

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        volume=one
        do dir=1,sdim
         dr=xsten(1,dir)-xsten(-1,dir)
         volume=volume*dr
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
         if (dir.eq.normdir+1) then
          centroid_map(dir)=coeff(1)*centroid(dir)+coeff(2)
         else
          centroid_map(dir)=centroid(dir)
         endif
        enddo ! dir=1..sdim
        volume_map=volume*coeff(1)
       else if (levelrz.eq.COORDSYS_RZ) then
        volume=one
        do dir=1,sdim
         dr=xsten(1,dir)-xsten(-1,dir)
         volume=volume*dr
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
         if (dir.eq.normdir+1) then
          centroid_map(dir)=coeff(1)*centroid(dir)+coeff(2)
         else
          centroid_map(dir)=centroid(dir)
         endif
        enddo ! dir=1..sdim
        volume_map=volume*coeff(1)

        rval=centroid(1)
        dr=xsten(1,1)-xsten(-1,1)
        if (rval.ne.zero) then
         volume=volume*two*Pi*abs(rval)
         centroid(1)=rval+dr*dr/(12.0*rval)
        else
         volume=zero
         centroid(1)=zero
        endif

        rval=centroid_map(1)
        if (normdir.eq.0) then
         dr=dr*coeff(1)
        endif
        if (rval.ne.zero) then
         volume_map=volume_map*two*Pi*abs(rval)
         centroid_map(1)=rval+dr*dr/(12.0*rval)
        else
         volume_map=zero
         centroid_map(1)=zero
        endif
       else
        print *,"levelrz invalid ",levelrz
        stop
       endif

      else if (sdim.eq.3) then
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        volume=one
        do dir=1,sdim
         dr=xsten(1,dir)-xsten(-1,dir)
         volume=volume*dr
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
         if (dir.eq.normdir+1) then
          centroid_map(dir)=coeff(1)*centroid(dir)+coeff(2)
         else
          centroid_map(dir)=centroid(dir)
         endif
        enddo ! dir=1..sdim
        volume_map=volume*coeff(1)

       else
        print *,"levelrz invalid ",levelrz
        stop
       endif
      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine Box_volumeFAST_and_map


! f(x,y,z)=a+b(x-x0)+c(y-y0)+d(z-z0)
!   3  4  
!   1  2
!
!   7 8
!   5 6
!
!  the 5 tetrahedra:
!   3,1,2,5
!   5,6,2,8
!   3,5,7,8
!   3,4,2,8
!   5,8,3,2
!
! in 2d
! 1,2,3
! 2,3,4

      subroutine get_ntetbox(ntetbox,symmetry_flag,sdim)
      IMPLICIT NONE

      integer, INTENT(out) :: ntetbox
      integer, INTENT(in) :: symmetry_flag,sdim

      if (sdim.eq.2) then
       if (symmetry_flag.eq.0) then
        ntetbox=2
       else if (symmetry_flag.eq.1) then
        ntetbox=4
       else
        print *,"symmetry_flag invalid"
        stop
       endif
      else if (sdim.eq.3) then
       if (symmetry_flag.eq.0) then
        ntetbox=5
       else if (symmetry_flag.eq.1) then
        ntetbox=24
       else
        print *,"symmetry_flag invalid"
        stop
       endif
      else
       print *,"sdim invalid"
       stop
      endif

      return
      end subroutine get_ntetbox


       ! centroid in absolute coordinate system
       ! returns a volume fraction
      subroutine getvolume( &
        bfact,dxgrid,xsten,nhalf, &
        ldata,volume,facearea, &
        centroid,EBVOFTOL,sdim)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: sdim,bfact,nhalf
      real(amrex_real), INTENT(in) :: EBVOFTOL
      real(amrex_real), INTENT(in) :: ldata(D_DECL(3,3,3))
      real(amrex_real) lnode(4*(sdim-1))
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dxgrid(sdim)
      real(amrex_real), INTENT(out) :: volume,facearea
      real(amrex_real) volcell
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real) cenall(sdim)
      integer dir

      if (bfact.lt.1) then
       print *,"bfact invalid130"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid getvolume"
       stop
      endif
      if (EBVOFTOL.gt.zero) then
       ! do nothing
      else
       print *,"getvolume: EBVOFTOL too small EBVOFTOL=",EBVOFTOL
       stop
      endif
      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (sdim.ne.2) then
        print *,"levelrz invalid get volume"
        stop
       endif
      else
       print *,"levelrz invalid get volume 2 ",levelrz
       stop
      endif
      if (AMREX_SPACEDIM.ne.sdim) then
       print *,"dimension mismatch"
       stop
      endif
      call data_to_node(ldata,lnode,1,xsten,nhalf,sdim)

      call fast_cell_intersection_grid( &
       bfact,dxgrid,xsten,nhalf, &       
       lnode, &
       volume,centroid,facearea, &
       volcell,cenall,sdim)

      if (volcell.le.zero) then
       volume=zero
      else
       volume=volume/volcell
      endif

      if (volume.le.EBVOFTOL) then
       volume=zero
       do dir=1,sdim
        centroid(dir)=cenall(dir)
       enddo
      endif
      if (volume.ge.one-EBVOFTOL) then
       volume=one
       do dir=1,sdim
        centroid(dir)=cenall(dir)
       enddo
      endif

      return
      end subroutine getvolume



         ! "centroid" in absolute coordinate system      
      subroutine getvolumebatch(bfact,dxgrid,xsten,nhalf, &
       ldata,volume,facearea, &
       centroid,EBVOFTOL,sdim)
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim,bfact,nhalf
      real(amrex_real), INTENT(in) :: EBVOFTOL
      integer im
      real(amrex_real), INTENT(in) :: ldata(D_DECL(3,3,3),num_materials)
      real(amrex_real) lnode(4*(sdim-1),num_materials)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dxgrid(sdim)
      real(amrex_real), INTENT(out) :: volume(num_materials)
      real(amrex_real), INTENT(out) :: facearea(num_materials)
      real(amrex_real) volcell
      real(amrex_real), INTENT(out) :: centroid(num_materials,sdim)
      real(amrex_real) cenall(sdim)

      if (bfact.lt.1) then
       print *,"bfact invalid131"
       stop
      endif

      if (EBVOFTOL.gt.zero) then
       !do nothing
      else
       print *,"getvolumebatch: EBVOFTOL too small EBVOFTOL=",EBVOFTOL
       stop
      endif
      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (sdim.ne.2) then
        print *,"sdim invalid"
        stop
       endif
      else
       print *,"levelrz invalid get volume batch ",levelrz
       stop
      endif 
      if (AMREX_SPACEDIM.ne.sdim) then
       print *,"dimension mismatch"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid getvolumebatch"
       stop
      endif

      call data_to_node(ldata,lnode,num_materials,xsten,nhalf,sdim)

      call fast_cell_intersection_grid_batch( &
        bfact,dxgrid,xsten,nhalf, &
        lnode, &
        volume, &
        centroid,facearea, &
        volcell,cenall,sdim)

      do im=1,num_materials
       if (volcell.le.zero) then
        volume(im)=zero
       else
        volume(im)=volume(im)/volcell
       endif

       if (volume(im).le.EBVOFTOL) then
        volume(im)=zero
       endif
       if (volume(im).ge.one-EBVOFTOL) then
        volume(im)=one
       endif
      enddo ! im

      return
      end subroutine getvolumebatch

      subroutine data_to_node(datasten,datanode,ncomp,xsten,nhalf,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: ncomp,nhalf
      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: datasten(D_DECL(3,3,3),ncomp)
      real(amrex_real), INTENT(out) :: datanode(4*(sdim-1),ncomp)
      integer i,j,k,i1,j1,k1,inode,im,klo,khi,dir
      real(amrex_real) xlonode(sdim)
      real(amrex_real) xhinode(sdim)
      real(amrex_real) xlocell(sdim)
      real(amrex_real) xhicell(sdim)
      real(amrex_real) xlo,xhi,wtprod,wtsum

      if (nhalf.lt.3) then
       print *,"nhalf invalid data to node"
       stop
      endif
      if (sdim.eq.3) then
       klo=-1
       khi=1
      else if (sdim.eq.2) then
       klo=0
       khi=0
      else
       print *,"sdim bust data_to_node"
       stop
      endif

      if (ncomp.lt.1) then
       print *,"ncomp invalid data_to_node"
       stop
      endif

       ! wt=|omega_cell int omega_node|/|omega_node|

      inode=1
      do k=klo,khi,2
      do j=-1,1,2
      do i=-1,1,2
       dir=1
       xlonode(dir)=xsten(i-1,dir)
       xhinode(dir)=xsten(i+1,dir)
       dir=2
       xlonode(dir)=xsten(j-1,dir)
       xhinode(dir)=xsten(j+1,dir)
       if (sdim.eq.3) then
        dir=sdim
        xlonode(dir)=xsten(k-1,dir)
        xhinode(dir)=xsten(k+1,dir)
       endif
       do im=1,ncomp
        datanode(inode,im)=zero
       enddo
       wtsum=zero
       do k1=klo,khi
       do j1=-1,1
       do i1=-1,1
        dir=1
        xlocell(dir)=xsten(2*i1-1,dir)
        xhicell(dir)=xsten(2*i1+1,dir)
        dir=2
        xlocell(dir)=xsten(2*j1-1,dir)
        xhicell(dir)=xsten(2*j1+1,dir)
        if (sdim.eq.3) then
         dir=sdim
         xlocell(dir)=xsten(2*k1-1,dir)
         xhicell(dir)=xsten(2*k1+1,dir)
        endif
        wtprod=one
        do dir=1,sdim
         xlo=max(xlocell(dir),xlonode(dir))
         xhi=min(xhicell(dir),xhinode(dir))
         if (xhi.gt.xlo) then
          wtprod=wtprod*(xhi-xlo)
         else
          wtprod=zero
         endif
        enddo
        do im=1,ncomp
         datanode(inode,im)=datanode(inode,im)+ &
           wtprod*datasten(D_DECL(i1+2,j1+2,k1+2),im)
        enddo
        wtsum=wtsum+wtprod
       enddo
       enddo
       enddo ! i1,j1,k1
       if (wtsum.le.zero) then
        print *,"wtsum invalid in data_to_node"
        stop
       endif
       do im=1,ncomp
        datanode(inode,im)=datanode(inode,im)/wtsum
       enddo
       inode=inode+1
      enddo  
      enddo  
      enddo  ! i,j,k

      return
      end subroutine data_to_node


        ! reconstruction relative to xsten0
        ! volume domain: xsten_grid
        ! shapeflag=0 find volumes within xsten_grid
        ! shapeflag=1 find volumes within xtet
      subroutine fast_cut_cell_intersection( &
        bfact,dx,xsten0,nhalf0, &
        slope,intercept, &
        volume, &
        centroid, &
        area, &
        xsten_grid,nhalf_grid, &
        xtet,shapeflag, &
        sdim)

      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim

      real(amrex_real) cum_volume,cum_area
      real(amrex_real) cum_centroid(sdim)

      integer, INTENT(in) :: bfact,nhalf0,nhalf_grid

      integer, INTENT(in) :: shapeflag
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      real(amrex_real), INTENT(in) :: xtet(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(in) :: intercept

      real(amrex_real), INTENT(out) :: volume,area
      real(amrex_real), INTENT(out) :: centroid(sdim)
      integer dir,inode
      integer j_dir
      integer i_grid_node
      integer j_grid_node
      integer k_grid_node
      integer i_tet_node
      real(amrex_real) xtarget(sdim)
      real(amrex_real) ls(sdim+1)

!   3  4  
!   1  2
!
!   7 8
!   5 6

      real(amrex_real) xnode(4*(sdim-1),sdim)
      real(amrex_real) phinode(4*(sdim-1))

      integer linearcut,fullelementfast,nodedomain

      linearcut=1
      fullelementfast=1

      if ((nhalf0.lt.1).or.(nhalf_grid.lt.1)) then
       print *,"nhalf invalid fast cut cell intersection"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid fast_cut_cell_intersection"
       stop
      endif

      if (shapeflag.eq.1) then

       do i_tet_node=1,sdim+1
        do j_dir=1,sdim
         xtarget(j_dir)=xtet(i_tet_node,j_dir)
        enddo
        call distfunc(bfact,dx,xsten0,nhalf0, &
         intercept,slope,xtarget,ls(i_tet_node),sdim)
       enddo  ! i_tet_node=1,sdim+1

       nodedomain=sdim+1
       call intersection_volume( &
        cum_volume,cum_area,cum_centroid, &
        ls,xtet,nodedomain, &
        sdim,fullelementfast,linearcut)

      else if (shapeflag.eq.0) then

       inode=1

       if (sdim.eq.3) then
        do k_grid_node=-1,1,2
        do j_grid_node=-1,1,2
        do i_grid_node=-1,1,2 
         do dir=1,sdim
          if (dir.eq.1) then
           xnode(inode,dir)=xsten_grid(i_grid_node,dir)
          else if (dir.eq.2) then
           xnode(inode,dir)=xsten_grid(j_grid_node,dir)
          else if (dir.eq.sdim) then
           xnode(inode,dir)=xsten_grid(k_grid_node,dir)
          else
           print *,"dir invalid face cut cell intersection"
           stop
          endif
          xtarget(dir)=xnode(inode,dir)
         enddo  ! dir
         call distfunc(bfact,dx,xsten0,nhalf0, &
           intercept,slope,xtarget, &
           phinode(inode),sdim)
         
         inode=inode+1
        enddo
        enddo
        enddo  ! i,j,k

       else if (sdim.eq.2) then

        do j_grid_node=-1,1,2
        do i_grid_node=-1,1,2 
         do dir=1,sdim
          if (dir.eq.1) then
           xnode(inode,dir)=xsten_grid(i_grid_node,dir)
          else if (dir.eq.2) then
           xnode(inode,dir)=xsten_grid(j_grid_node,dir)
          else
           print *,"dir invalid face cut cell intersection 2"
           stop
          endif
          xtarget(dir)=xnode(inode,dir)
         enddo  ! dir
         call distfunc(bfact,dx,xsten0,nhalf0, &
           intercept,slope,xtarget, &
           phinode(inode),sdim)
         
         inode=inode+1
        enddo
        enddo ! i,j
       else
        print *,"sdim invalid"
        stop
       endif

       nodedomain=4*(sdim-1)

       if (inode.ne.nodedomain+1) then
        print *,"inode invalid2"
        print *,"inode=",inode
        print *,"nodedomain=",nodedomain
        stop
       endif
       call intersection_volume( &
        cum_volume,cum_area,cum_centroid, &
        phinode,xnode,nodedomain, &
        sdim,fullelementfast,linearcut)
      else
       print *,"shapeflag invalid"
       stop
      endif

      volume=cum_volume
      area=cum_area
      do dir=1,sdim
       centroid(dir)=cum_centroid(dir)
      enddo

      return
      end subroutine fast_cut_cell_intersection


        ! reconstruction relative to xsten0
        ! volume domain: xsten_grid
      subroutine fast_cut_cell_intersection_simple( &
        bfact,dx,xsten0,nhalf0, &
        slope,intercept, &
        volume,centroid, &
        xsten_grid,nhalf_grid, &
        sdim)

      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim

      real(amrex_real) cum_volume
      real(amrex_real) cum_centroid(sdim)

      integer, INTENT(in) :: bfact,nhalf0,nhalf_grid

      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(in) :: intercept

      real(amrex_real), INTENT(out) :: volume
      real(amrex_real), INTENT(out) :: centroid(sdim)
      integer dir,inode
      real(amrex_real) xtarget(sdim)
      integer i_grid_node
      integer j_grid_node
      integer k_grid_node

!   3  4  
!   1  2
!
!   7 8
!   5 6

      real(amrex_real) xnode(4*(sdim-1),sdim)
      real(amrex_real) phinode(4*(sdim-1))

      integer nodedomain

      if ((nhalf0.lt.1).or.(nhalf_grid.lt.1)) then
       print *,"nhalf invalid fast cut cell intersection"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid fast_cut_cell_intersection"
       stop
      endif

      inode=1

      if (sdim.eq.3) then
       do k_grid_node=-1,1,2
       do j_grid_node=-1,1,2
       do i_grid_node=-1,1,2 
        do dir=1,sdim
         if (dir.eq.1) then
          xnode(inode,dir)=xsten_grid(i_grid_node,dir)
         else if (dir.eq.2) then
          xnode(inode,dir)=xsten_grid(j_grid_node,dir)
         else if (dir.eq.sdim) then
          xnode(inode,dir)=xsten_grid(k_grid_node,dir)
         else
          print *,"dir invalid face cut cell intersection"
          stop
         endif
         xtarget(dir)=xnode(inode,dir)
        enddo  ! dir
        call distfunc(bfact,dx,xsten0,nhalf0, &
          intercept,slope,xtarget, &
          phinode(inode),sdim)
         
        inode=inode+1
       enddo
       enddo
       enddo  ! i,j,k

      else if (sdim.eq.2) then

       do j_grid_node=-1,1,2
       do i_grid_node=-1,1,2 
        do dir=1,sdim
         if (dir.eq.1) then
          xnode(inode,dir)=xsten_grid(i_grid_node,dir)
         else if (dir.eq.2) then
          xnode(inode,dir)=xsten_grid(j_grid_node,dir)
         else
          print *,"dir invalid face cut cell intersection 2"
          stop
         endif
         xtarget(dir)=xnode(inode,dir)
        enddo  ! dir
        call distfunc(bfact,dx,xsten0,nhalf0, &
          intercept,slope,xtarget, &
          phinode(inode),sdim)
         
        inode=inode+1
       enddo
       enddo ! i,j
      else
       print *,"sdim invalid"
       stop
      endif

      nodedomain=4*(sdim-1)

      if (inode.ne.nodedomain+1) then
       print *,"inode invalid3"
       print *,"inode=",inode
       print *,"nodedomain=",nodedomain
       stop
      endif
      call intersection_volume_simple( &
        cum_volume,cum_centroid, &
        phinode,xnode,nodedomain, &
        sdim)

      volume=cum_volume
      do dir=1,sdim
       centroid(dir)=cum_centroid(dir)
      enddo

      return
      end subroutine fast_cut_cell_intersection_simple


        ! reconstruction relative to xsten0
        ! volume domain: xsten_grid
      subroutine fast_cut_cell_intersection_and_map( &
        normdir, &
        coeff, &
        bfact,dx,xsten0,nhalf0, &
        slope,intercept, &
        volume,centroid, &
        volume_map,centroid_map, &
        xsten_grid,nhalf_grid, &
        sdim)

      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim

      integer, INTENT(in) :: normdir
      real(amrex_real), INTENT(in) :: coeff(2)
      real(amrex_real) cum_volume
      real(amrex_real) cum_volume_map
      real(amrex_real) cum_centroid(sdim)
      real(amrex_real) cum_centroid_map(sdim)

      integer, INTENT(in) :: bfact,nhalf0,nhalf_grid

      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(in) :: intercept

      real(amrex_real), INTENT(out) :: volume
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real), INTENT(out) :: volume_map
      real(amrex_real), INTENT(out) :: centroid_map(sdim)
      integer dir,inode
      real(amrex_real) xtarget(sdim)
      integer i_grid_node
      integer j_grid_node
      integer k_grid_node

!   3  4  
!   1  2
!
!   7 8
!   5 6

      real(amrex_real) xnode(4*(sdim-1),sdim)
      real(amrex_real) phinode(4*(sdim-1))

      integer nodedomain

      if ((nhalf0.lt.1).or.(nhalf_grid.lt.1)) then
       print *,"nhalf invalid fast cut cell intersection_and_map"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid fast_cut_cell_intersection_and_map"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.sdim)) then
       print *,"normdir invalid"
       stop
      endif
      if (coeff(1).le.zero) then
       print *,"coeff(1) invalid"
       stop
      endif

      inode=1

      if (sdim.eq.3) then
       do k_grid_node=-1,1,2
       do j_grid_node=-1,1,2
       do i_grid_node=-1,1,2 
        do dir=1,sdim
         if (dir.eq.1) then
          xnode(inode,dir)=xsten_grid(i_grid_node,dir)
         else if (dir.eq.2) then
          xnode(inode,dir)=xsten_grid(j_grid_node,dir)
         else if (dir.eq.sdim) then
          xnode(inode,dir)=xsten_grid(k_grid_node,dir)
         else
          print *,"dir invalid face cut cell intersection and map"
          stop
         endif
         xtarget(dir)=xnode(inode,dir)
        enddo  ! dir
        call distfunc(bfact,dx,xsten0,nhalf0, &
          intercept,slope,xtarget, &
          phinode(inode),sdim)
         
        inode=inode+1
       enddo
       enddo
       enddo  ! i,j,k

      else if (sdim.eq.2) then

       do j_grid_node=-1,1,2
       do i_grid_node=-1,1,2 
        do dir=1,sdim
         if (dir.eq.1) then
          xnode(inode,dir)=xsten_grid(i_grid_node,dir)
         else if (dir.eq.2) then
          xnode(inode,dir)=xsten_grid(j_grid_node,dir)
         else
          print *,"dir invalid face cut cell intersection and map 2"
          stop
         endif
         xtarget(dir)=xnode(inode,dir)
        enddo  ! dir
        call distfunc(bfact,dx,xsten0,nhalf0, &
          intercept,slope,xtarget, &
          phinode(inode),sdim)
        
        inode=inode+1
       enddo
       enddo ! i,j
      else
       print *,"sdim invalid"
       stop
      endif

      nodedomain=4*(sdim-1)

      if (inode.ne.nodedomain+1) then
       print *,"inode invalid4"
       print *,"inode=",inode
       print *,"nodedomain=",nodedomain
       stop
      endif
      call intersection_volume_and_map( &
       normdir, &
       coeff, &
       cum_volume,cum_centroid, &
       cum_volume_map,cum_centroid_map, &
       phinode,xnode,nodedomain, &
       sdim)

      volume=cum_volume
      volume_map=cum_volume_map
      do dir=1,sdim
       centroid(dir)=cum_centroid(dir)
       centroid_map(dir)=cum_centroid_map(dir)
      enddo

      return
      end subroutine fast_cut_cell_intersection_and_map


        ! returns centroid in absolute coordinate system

      subroutine fast_cell_intersection_grid_batch( &
        bfact,dxgrid,xsten0,nhalf0, &
        lnodebatch, &
        volumedark,centroiddark, &
        area,volall,cenall,sdim)

      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: nhalf0

      real(amrex_real) cum_volume,cum_area
      real(amrex_real) cum_centroid(sdim)

      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: lnodebatch(4*(sdim-1),num_materials)
      real(amrex_real) lnode(4*(sdim-1))
      real(amrex_real) xnode(4*(sdim-1),sdim)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dxgrid(sdim)
      real(amrex_real), INTENT(out) :: volumedark(num_materials)
      real(amrex_real), INTENT(out) :: centroiddark(num_materials,sdim)
      real(amrex_real), INTENT(out) :: volall
      real(amrex_real), INTENT(out) :: area(num_materials)
      real(amrex_real), INTENT(out) :: cenall(sdim)
      integer j_dir
      integer i_grid_node
      integer j_grid_node
      integer k_grid_node
      integer inode,im
      integer linearcut,fullelementfast,nodedomain

      if (bfact.lt.1) then
       print *,"bfact invalid135"
       stop
      endif
      if ((num_materials.lt.1).or.(num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid fast cell intersection grid batch"
       stop
      endif

      linearcut=0
      fullelementfast=1
      nodedomain=4*(sdim-1)

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust fast cell intersection grid batch"
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif

      inode=1

      call Box_volumeFAST(bfact,dxgrid,xsten0,nhalf0,volall,cenall,sdim)

      if (volall.gt.zero) then
     
       if (sdim.eq.3) then
        do k_grid_node=-1,1,2
        do j_grid_node=-1,1,2
        do i_grid_node=-1,1,2 
         xnode(inode,1)=xsten0(i_grid_node,1)
         xnode(inode,2)=xsten0(j_grid_node,2)
         xnode(inode,sdim)=xsten0(k_grid_node,sdim)
         inode=inode+1
        enddo
        enddo
        enddo
       else if (sdim.eq.2) then
        do j_grid_node=-1,1,2
        do i_grid_node=-1,1,2 
         xnode(inode,1)=xsten0(i_grid_node,1)
         xnode(inode,2)=xsten0(j_grid_node,2)
         inode=inode+1
        enddo
        enddo
       else 
        print *,"sdim invalid"
        stop
       endif

       if (inode.ne.nodedomain+1) then
        print *,"inode invalid5"
        print *,"inode=",inode
        print *,"nodedomain=",nodedomain
        stop
       endif

       do im=1,num_materials
        do inode=1,nodedomain
         lnode(inode)=lnodebatch(inode,im)
        enddo        
        call intersection_volume( &
         cum_volume,cum_area,cum_centroid, &
         lnode,xnode,nodedomain, &
         sdim,fullelementfast,linearcut)

        volumedark(im)=cum_volume
        area(im)=cum_area
        do j_dir=1,sdim
         centroiddark(im,j_dir)=cum_centroid(j_dir)
        enddo
       enddo  ! im

      else if (volall.eq.zero) then

       print *,"WARNING:volall cannot be zero: fast_cell_intersection_grid"
       do im=1,num_materials
        volumedark(im)=zero
        area(im)=zero
        do j_dir=1,sdim
         centroiddark(im,j_dir)=zero
        enddo
       enddo  ! im

      else if (volall.lt.zero) then
       print *,"volall invalid"
       stop
      endif

      return
      end subroutine fast_cell_intersection_grid_batch

        ! returns centroid in absolute coordinate system

      subroutine fast_cell_intersection_grid( &
        bfact,dxgrid,xsten0,nhalf0, &
        lnode, &
        volumedark,centroiddark, &
        area,volall,cenall,sdim)

      IMPLICIT NONE

      integer, INTENT(in) :: sdim

      real(amrex_real) cum_volume,cum_area
      real(amrex_real) cum_centroid(sdim)

      integer, INTENT(in) :: bfact,nhalf0
      real(amrex_real), INTENT(in) :: lnode(4*(sdim-1))
      real(amrex_real) xnode(4*(sdim-1),sdim)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dxgrid(sdim)
      real(amrex_real), INTENT(out) :: volumedark
      real(amrex_real), INTENT(out) :: centroiddark(sdim)
      real(amrex_real), INTENT(out) :: volall,area
      real(amrex_real), INTENT(out) :: cenall(sdim)
      integer i_grid_node
      integer j_grid_node
      integer k_grid_node
      integer j_dir
      integer inode
      integer linearcut,fullelementfast,nodedomain

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      linearcut=0
      fullelementfast=1
      nodedomain=4*(sdim-1)

      if (bfact.lt.1) then
       print *,"bfact invalid135"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust fast cell intersection grid"
       stop
      endif

      inode=1

      call Box_volumeFAST(bfact,dxgrid,xsten0,nhalf0,volall,cenall,sdim)

      if (volall.gt.zero) then
     
       if (sdim.eq.3) then
        do k_grid_node=-1,1,2
        do j_grid_node=-1,1,2
        do i_grid_node=-1,1,2 
         xnode(inode,1)=xsten0(i_grid_node,1)
         xnode(inode,2)=xsten0(j_grid_node,2)
         xnode(inode,sdim)=xsten0(k_grid_node,sdim)
         inode=inode+1
        enddo
        enddo
        enddo
       else if (sdim.eq.2) then
        do j_grid_node=-1,1,2
        do i_grid_node=-1,1,2 
         xnode(inode,1)=xsten0(i_grid_node,1)
         xnode(inode,2)=xsten0(j_grid_node,2)
         inode=inode+1
        enddo
        enddo
       else 
        print *,"sdim invalid"
        stop
       endif

       if (inode.ne.nodedomain+1) then
        print *,"inode invalid6"
        print *,"inode=",inode
        print *,"nodedomain=",nodedomain
        stop
       endif
       call intersection_volume( &
        cum_volume,cum_area,cum_centroid, &
        lnode,xnode,nodedomain, &
        sdim,fullelementfast,linearcut)

       volumedark=cum_volume
       area=cum_area
       do j_dir=1,sdim
        centroiddark(j_dir)=cum_centroid(j_dir)
       enddo

      else if (volall.eq.zero) then
       print *,"WARNING:volall cannot be zero: fast_cell_intersection_grid"
       volumedark=zero
       area=zero
       do j_dir=1,sdim
        centroiddark(j_dir)=zero
       enddo
      else if (volall.lt.zero) then
       print *,"volall invalid"
       stop
      endif

      return
      end subroutine fast_cell_intersection_grid

end module geometry_intersect_module


module MOF_routines_module
use amrex_fort_module, only : amrex_real

      integer :: MOF_DEBUG_RECON_COUNT
       ! 1=>output errors when fastflag=0 or 1
       ! 2=>output errors when fastflag=0 
      integer :: MOF_DEBUG_RECON=0

contains

      subroutine set_MOFITERMAX(MOFITERMAX_in)
      IMPLICIT NONE

      integer, INTENT(in) :: MOFITERMAX_in

#include "mofdata.H"

      MOFITERMAX=MOFITERMAX_in
      if ((MOFITERMAX.lt.0).or. &
          (MOFITERMAX.gt.MOFITERMAX_LIMIT)) then
       print *,"MOFITERMAX invalid in set mofitermax"
       print *,"MOFITERMAX ",MOFITERMAX
       print *,"MOFITERMAX_LIMIT ",MOFITERMAX_LIMIT
       stop
      endif

      return
      end subroutine set_MOFITERMAX


       ! order_algorithm=0 => try different combinations and
       ! choose combination with smallest MOF error
      subroutine set_order_algorithm( &
         order_algorithm_in, &
         renormalize_order_algorithm_in)
      use probcommon_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: order_algorithm_in(num_materials)
      integer, INTENT(in) :: renormalize_order_algorithm_in(num_materials)
      integer im
      integer im_debug
      integer irank

#include "mofdata.H"

      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid set order algorithm"
       print *,"num_materials= ",num_materials
       stop
      endif

      do im=1,num_materials
       order_algorithm(im)=order_algorithm_in(im)
       if (order_algorithm(im).lt.0) then
        print *,"order_alg bust im,order_algorithm ",im,order_algorithm(im)
        stop
       endif
       renormalize_order_algorithm(im)=renormalize_order_algorithm_in(im)
       if (renormalize_order_algorithm(im).lt.0) then
        print *,"renormalize_order_alg bust im,renormalize_order_algorithm ", &
           im,renormalize_order_algorithm(im)
        stop
       endif
       irank=renormalize_order_algorithm(im)
       if (irank.ne.0) then
        if (rank_algorithm(irank).eq.im) then
         !do nothing
        else
         print *,"im=",im
         print *,"irank=",irank
         print *,"rank_algorithm invalid: "
         do im_debug=1,num_materials
          print *,"im_debug,rank_algorithm ", &
             im_debug,rank_algorithm(im_debug)
          print *,"im_debug,renormalize_order_algorithm ", &
             im_debug,renormalize_order_algorithm(im_debug)
         enddo
         stop
        endif
       endif
      enddo !im=1..num_materials

      return
      end subroutine set_order_algorithm


! intercept=-maxphi  => refvfrac=0
! intercept=-minphi  => refvfrac=1
      subroutine multi_phi_bounds( &
        bfact,dx, &
        xsten,nhalf, &
        slope, &
        xtetlist, &
        nlist_alloc, &
        nlist, &
        nmax, &
        minphi,maxphi, &
        sdim)

      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: nlist,nmax,sdim,bfact,nhalf
      real(amrex_real), INTENT(in) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real) :: xsten_local(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real) xtarget(sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      integer dir
      integer i_tet_node
      integer i_stencil_node
      integer n
      real(amrex_real), INTENT(out) :: minphi,maxphi
      real(amrex_real) intercept,dist

      if (nhalf.lt.1) then
       print *,"nhalf invalid multi phi bounds 3d"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid multi_phi_bounds"
       stop
      endif
      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif


       ! xsten_local used for LS dist.
      do i_stencil_node=-nhalf,nhalf
      do dir=1,sdim
       xsten_local(i_stencil_node,dir)=xsten(i_stencil_node,dir)
      enddo !dir
      enddo !i_stencil_node

      intercept=zero
      minphi=1.0D+10
      maxphi=-1.0D+10
      do n=1,nlist
       do i_tet_node=1,sdim+1
        do dir=1,sdim
         xtarget(dir)=xtetlist(i_tet_node,dir,n)
        enddo

        call distfunc(bfact,dx,xsten_local,nhalf, &
         intercept,slope,xtarget,dist,sdim)

        if (dist.lt.minphi) then
         minphi=dist
        endif
        if (dist.gt.maxphi) then
         maxphi=dist
        endif
       enddo  ! sweeping nodes of triangle
      enddo ! sweeping triangles

      if (((minphi.eq.zero).and. &
           (maxphi.eq.zero)).or. &
          (minphi.ge.maxphi)) then
       print *,"cannot have zero slope"
       print *,"minphi=",minphi
       print *,"maxphi=",maxphi
       print *,"slopexyz=",slope(1),slope(2),slope(sdim)
       print *,"xcell xyz=",xsten(0,1),xsten(0,2),xsten(0,sdim)
       stop
      else if (((minphi.ne.zero).or. &
                (maxphi.ne.zero)).and. &
               (minphi.lt.maxphi)) then
       ! do nothing
      else
       print *,"minphi or maxphi is NaN: multi_phi_bounds: ",minphi,maxphi
       stop
      endif

      return
      end subroutine multi_phi_bounds 
 
! phi>0 in the material.  n points into the phi>0 region.
      subroutine closest(xint,x,nn,phi,sdim)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(out) :: xint(sdim)
      real(amrex_real), INTENT(in) :: x(sdim)
      real(amrex_real), INTENT(in) :: nn(sdim)
      real(amrex_real), INTENT(in) :: phi
      integer dir

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (sdim.ne.2) then
        print *,"dimension bust"
       endif
      else
       print *,"levelrz invalid ",levelrz
       stop
      endif

      do dir=1,sdim
       xint(dir)=x(dir)-nn(dir)*phi
      enddo

      return
      end subroutine closest

 
      subroutine dist2or3D(x1,x2,dist,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer dir
      real(amrex_real), INTENT(in) :: x1(sdim),x2(sdim)
      real(amrex_real), INTENT(out) :: dist

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"dimension bust"
       stop
      endif
      dist=zero
      do dir=1,sdim
       dist=dist+(x1(dir)-x2(dir))**2
      enddo
      dist=sqrt(dist)

      return
      end subroutine dist2or3D


! find volume, centroid and area of intersection of a line with a 
! collection of triangles.
      subroutine multi_cell_intersection( &
        bfact,dx, &
        xsten,nhalf, &
        slope,intercept, &
        volume,centroid,area, &
        xtetlist, & !intent(in)
        nlist_alloc, &
        nlist, &
        nmax, &
        sdim)

      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: bfact,nhalf

      integer, INTENT(in) :: nlist,nmax,sdim
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(in) :: intercept

      real(amrex_real), INTENT(out) :: volume,area
      real(amrex_real) volumelist,arealist
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real) centroidlist(sdim)
      integer i_tet_node
      integer j_dir
      integer n
      real(amrex_real) xx(sdim+1,sdim)
      real(amrex_real) xtarget(sdim)
      real(amrex_real) ls(sdim+1)

      if (nhalf.lt.1) then
       print *,"nhalf invalid multi cell intersection"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_cell_intersection"
       stop
      endif
      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      do j_dir=1,sdim
       centroid(j_dir)=zero
      enddo
      area=zero
      volume=zero

      do n=1,nlist
       do i_tet_node=1,sdim+1
        do j_dir=1,sdim
         xx(i_tet_node,j_dir)=xtetlist(i_tet_node,j_dir,n)
         xtarget(j_dir)=xx(i_tet_node,j_dir)
        enddo
        call distfunc(bfact,dx,xsten,nhalf, &
         intercept,slope,xtarget,ls(i_tet_node),sdim)
       enddo

       if (sdim.eq.3) then
        call intersection_volumeXYZ(ls,xx,volumelist, &
         centroidlist,arealist,sdim)
       else if (sdim.eq.2) then
        call int_volumeXYorRZ(ls,xx,volumelist, &
         centroidlist,arealist,sdim)
       else
        print *,"sdim invalid"
        stop
       endif

       volume=volume+volumelist
       area=area+arealist
       do j_dir=1,sdim
        centroid(j_dir)=centroid(j_dir)+centroidlist(j_dir)*volumelist
       enddo
      enddo ! n

      if (volume.gt.zero) then
       do j_dir=1,sdim
        centroid(j_dir)=centroid(j_dir)/volume
       enddo
      else
       volume=zero
       do j_dir=1,sdim
        centroid(j_dir)=zero
       enddo
      endif

      return
      end subroutine multi_cell_intersection

! find volume, centroid of intersection of a line with a 
! collection of triangles.
      subroutine multi_cell_intersection_simple( &
        bfact,dx,xsten,nhalf, &
        slope,intercept, &
        volume,centroid, &
        xtetlist, &
        nlist_alloc, &
        nlist, &
        nmax, &
        sdim)

      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: bfact,nhalf

      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: nlist,nmax,sdim
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(in) :: intercept

      real(amrex_real), INTENT(out) :: volume
      real(amrex_real) volumelist
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real) centroidlist(sdim)
      integer i_tet_node
      integer j_dir
      integer n
      real(amrex_real) xx(sdim+1,sdim)
      real(amrex_real) xtarget(sdim)
      real(amrex_real) ls(sdim+1)

      if (nhalf.lt.1) then
       print *,"nhalf invalid multi cell intersection"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_cell_intersection"
       stop
      endif
      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      do j_dir=1,sdim
       centroid(j_dir)=zero
      enddo
      volume=zero

      do n=1,nlist
       do i_tet_node=1,sdim+1
        do j_dir=1,sdim
         xx(i_tet_node,j_dir)=xtetlist(i_tet_node,j_dir,n)
         xtarget(j_dir)=xx(i_tet_node,j_dir)
        enddo
        call distfunc(bfact,dx,xsten,nhalf, &
         intercept,slope,xtarget,ls(i_tet_node),sdim)
       enddo

       if (sdim.eq.3) then
        call intersection_volumeXYZ_simple(ls,xx,volumelist, &
         centroidlist,sdim)
       else if (sdim.eq.2) then
        call int_volumeXYorRZ_simple(ls,xx,volumelist, &
         centroidlist,sdim)
       else
        print *,"sdim invalid"
        stop
       endif

       volume=volume+volumelist
       do j_dir=1,sdim
        centroid(j_dir)=centroid(j_dir)+centroidlist(j_dir)*volumelist
       enddo
      enddo ! n

      if (volume.gt.zero) then
       do j_dir=1,sdim
        centroid(j_dir)=centroid(j_dir)/volume
       enddo
      else
       volume=zero
       do j_dir=1,sdim
        centroid(j_dir)=zero
       enddo
      endif

      return
      end subroutine multi_cell_intersection_simple



! find volume, centroid and area of intersection of a line with a 
! collection of triangles.
      subroutine multi_cell_intersection_and_map( &
        normdir, &
        coeff, &
        bfact,dx,xsten,nhalf, &
        slope,intercept, &
        volume,centroid, &
        volume_map,centroid_map, &
        xtetlist, &
        nlist_alloc, &
        nlist, &
        nmax, &
        sdim)

      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: normdir
      real(amrex_real), INTENT(in) :: coeff(2)

      integer, INTENT(in) :: bfact,nhalf

      integer, INTENT(in) :: nlist,nmax,sdim
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(in) :: intercept

      real(amrex_real), INTENT(out) :: volume,volume_map
      real(amrex_real) volumelist
      real(amrex_real) volumelist_map
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real), INTENT(out) :: centroid_map(sdim)
      real(amrex_real) centroidlist(sdim)
      real(amrex_real) centroidlist_map(sdim)
      integer i_tet_node
      integer j_dir
      integer n
      real(amrex_real) xx(sdim+1,sdim)
      real(amrex_real) xtarget(sdim)
      real(amrex_real) ls(sdim+1)

      if (nhalf.lt.1) then
       print *,"nhalf invalid multi cell intersection"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_cell_intersection"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.sdim)) then
       print *,"normdir invalid"
       stop
      endif
      if (coeff(1).le.zero) then
       print *,"coeff(1) invalid"
       stop
      endif
      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      do j_dir=1,sdim
       centroid(j_dir)=zero
       centroid_map(j_dir)=zero
      enddo
      volume=zero
      volume_map=zero

      do n=1,nlist
       do i_tet_node=1,sdim+1
        do j_dir=1,sdim
         xx(i_tet_node,j_dir)=xtetlist(i_tet_node,j_dir,n)
         xtarget(j_dir)=xx(i_tet_node,j_dir)
        enddo
        call distfunc(bfact,dx,xsten,nhalf, &
         intercept,slope,xtarget,ls(i_tet_node),sdim)
       enddo

       if (sdim.eq.3) then
        call intersection_volumeXYZ_and_map( &
         normdir,coeff, &
         ls,xx, &
         volumelist, &
         volumelist_map, &
         centroidlist, &
         centroidlist_map, &
         sdim)
       else if (sdim.eq.2) then
        call int_volumeXYorRZ_and_map( &
         normdir,coeff, &
         ls,xx, &
         volumelist, &
         volumelist_map, &
         centroidlist, &
         centroidlist_map, &
         sdim)
       else
        print *,"sdim invalid"
        stop
       endif

       volume=volume+volumelist
       volume_map=volume_map+volumelist_map
       do j_dir=1,sdim
        centroid(j_dir)=centroid(j_dir)+centroidlist(j_dir)*volumelist
        centroid_map(j_dir)=centroid_map(j_dir)+ &
                centroidlist_map(j_dir)*volumelist_map
       enddo
      enddo ! n

      if ((volume.gt.zero).and.(volume_map.gt.zero)) then
       do j_dir=1,sdim
        centroid(j_dir)=centroid(j_dir)/volume
        centroid_map(j_dir)=centroid_map(j_dir)/volume_map
       enddo
      else if ((volume.eq.zero).or.(volume_map.eq.zero)) then
       volume=zero
       volume_map=zero
       do j_dir=1,sdim
        centroid(j_dir)=zero
        centroid_map(j_dir)=zero
       enddo
      else
       print *,"volume or volume_map invalid"
       stop
      endif

      return
      end subroutine multi_cell_intersection_and_map


       ! centroid relative to absolute coordinate system

      subroutine multi_ff(bfact,dx, &
        xsten0,nhalf0, &
        ff,slope, &
        intercept, & !intent(in)
        arean, &
        vtarget, &
        xtetlist, & !intent(in)
        nlist_alloc, &
        centroid, & !intent(out)
        nlist, &
        nmax, &
        fastflag,sdim)

      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: bfact,nhalf0
      integer, INTENT(in) :: nlist,nmax,fastflag
      real(amrex_real), INTENT(in) :: xtetlist(4,3,nlist_alloc)
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(in) :: intercept
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)

      real(amrex_real), INTENT(out) :: ff
      real(amrex_real) voln
      real(amrex_real), INTENT(out) :: arean
      real(amrex_real), INTENT(in) :: vtarget
      real(amrex_real) volcell
      real(amrex_real) cencell(sdim)
      real(amrex_real), INTENT(out) :: centroid(sdim)
      integer, PARAMETER :: shapeflag=0 !regular hexahedron
      real(amrex_real) xtet(sdim+1,sdim)
      real(amrex_real) xsten_local(-nhalf0:nhalf0,sdim)
      integer i,dir

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_ff"
       stop
      endif
      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif

      call Box_volumeFAST( &
        bfact,dx,xsten0,nhalf0, &
        volcell, &
        cencell, &
        sdim)

      if (volcell.gt.zero) then
       ! do nothing
      else
       print *,"volcell bust: ",volcell
       stop
      endif

      if (fastflag.eq.0) then

         ! xsten_local used for LS dist.
       do i=-nhalf0,nhalf0
       do dir=1,sdim
         xsten_local(i,dir)=xsten0(i,dir)
       enddo
       enddo

         !calling from "multi_ff"
       call multi_cell_intersection( &
         bfact,dx, &
         xsten_local,nhalf0, &
         slope, &
         intercept, &
         voln, &
         centroid, &
         arean, &
         xtetlist, & !intent(in)
         nlist_alloc, &
         nlist, &
         nmax, &
         sdim)

      else if (fastflag.eq.1) then

        call fast_cut_cell_intersection( &
         bfact,dx,xsten0,nhalf0, &
         slope,intercept, &
         voln,centroid, &
         arean, &
         xsten0,nhalf0, &
         xtet,shapeflag, &
         sdim) 

      else
       print *,"fastflag invalid multi_ff"
       stop
      endif

      ff=(voln-vtarget)/volcell

      return
      end subroutine multi_ff


      subroutine single_ff(bfact,dx, &
        xsten0,nhalf0, &
        ff,slope,intercept, &
        arean, &
        vtarget, &
        centroid,sdim)

      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: bfact,nhalf0
      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(in) :: intercept
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)

      real(amrex_real), INTENT(out) :: ff
      real(amrex_real) voln
      real(amrex_real), INTENT(out) :: arean
      real(amrex_real), INTENT(in) :: vtarget
      real(amrex_real) volcell
      real(amrex_real) cencell(sdim)
      real(amrex_real), INTENT(out) :: centroid(sdim)
      integer, PARAMETER :: shapeflag=0 !regular hexahedron
      real(amrex_real) xtet(sdim+1,sdim)

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid single_ff"
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif

      call Box_volumeFAST(bfact,dx,xsten0,nhalf0,volcell, &
       cencell,sdim)

      if (volcell.gt.zero) then
       !do nothing
      else
       print *,"volcell bust: ",volcell
       stop
      endif

      call fast_cut_cell_intersection( &
       bfact,dx,xsten0,nhalf0, &
       slope,intercept, &
       voln,centroid, &
       arean, &
       xsten0,nhalf0, &
       xtet,shapeflag, &
       sdim) 

      ff=(voln-vtarget)/volcell

      return
      end subroutine single_ff


       ! centroid in absolute coordinate system
      subroutine multi_find_intercept( &
       tid_in, &
       nDOF, & ! sdim-1
       nEQN, & ! sdim
       bfact,dx, &
       xsten0,nhalf0, &
       slope, & !intent(in)
       intercept, &
       xtetlist, & !intent(in)
       nlist_alloc, & !intent(in)
       nlist, & !intent(in)
       nmax, & !intent(in)
       vfrac, &
       use_initial_guess, &
       centroid, & !intent(out)
       fastflag,sdim)

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: nDOF ! sdim-1
      integer, INTENT(in) :: nEQN ! sdim 
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: bfact,nhalf0
      integer, INTENT(in):: nlist
      integer, INTENT(in):: nlist_alloc
      integer, INTENT(in) :: nmax,fastflag
      real(amrex_real),INTENT(in):: xtetlist(4,3,nlist_alloc)
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(inout) :: intercept
      real(amrex_real), INTENT(in) :: vfrac
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)

      integer niter,maxiter
      integer i,j,k,dir
      real(amrex_real) minphi,maxphi
      real(amrex_real) volcell
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real) cencell(sdim)
      real(amrex_real) arean
      real(amrex_real) err
      real(amrex_real) moftol
      real(amrex_real) min_err
      real(amrex_real) vtarget,fc
      integer debug_root
      integer, INTENT(in) :: use_initial_guess
      real(amrex_real) err_default,fc_default,arean_default,intercept_default
      real(amrex_real) volcut
      real(amrex_real) cencut(sdim)
      real(amrex_real) xtarget(sdim)
      real(amrex_real) vfrac_normalize
      real(amrex_real) null_intercept,dist
      integer klo_stencil,khi_stencil
      integer nn
      real(amrex_real) intercept_upper,intercept_lower
      real(amrex_real) intercept_test,aa,bb,fa,fb

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid (multi_find_intercept) ",tid_in
       stop
      endif 

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_find_intercept"
       stop
      endif
      if ((nlist.ge.0).and.(nlist_alloc.ge.1)) then
       ! do nothing
      else
       print *,"nlist or nlist_alloc invalid"
       stop
      endif
      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      if (sdim.eq.2) then
       klo_stencil=0
       khi_stencil=0
      else if (sdim.eq.3) then
       klo_stencil=-1
       khi_stencil=1
      else
       print *,"sdim invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid 9514 ",bfact
       stop
      endif

      debug_root=0

      if (INTERCEPT_MAXITER_NEWTON.lt.INTERCEPT_MAXITER) then
       ! do nothing
      else
       print *,"INTERCEPT_MAXITER_NEWTON invalid"
       stop
      endif

      maxiter=INTERCEPT_MAXITER
       ! INTERCEPT_TOL is declared in PROBCOMMON.F90
      moftol=INTERCEPT_TOL
      min_err=-one

! phi=n dot (x-x0)+int
! find max,min n dot (x-x0)
! if fastflag=0, 
!  search the vertices of all triangles that make up the "cut" domain.

      call Box_volumeFAST( &
        bfact,dx, &
        xsten0,nhalf0, &
        volcell, &
        cencell, &
        sdim)

       ! at each node x_i one solves:
       !  n dot (x_i-x0) + b_i = 0
       !  b_i=-n dot (x_i-x0)
       !  b_Lower_bound=min b_i = - maxphi
       !  b_Upper_bound=max b_i = - minphi
       !  In otherwords:
       !  if intercept>max b_i, then F=1
       !  if intercept<min b_i, then F=0
       !  min b_i <= intercept <= max b_i
       !  i.e.
       !  -maxphi <=intercept <= -minphi
       !  i.e.
       !  intercept_lower=-maxphi
       !  intercept_upper=-minphi
       !  Also define,
       !  maxphi=max n dot (x_i -x0)=-min -(n dot (x_i-x0))=-min b_i
       !    =>min b_i = -maxphi
       !  minphi=min n dot (x_i -x0)=-max -(n dot (x_i-x0))=-max b_i
       !    =>max b_i = -minphi
      if (fastflag.eq.0) then

       call get_cut_geom3D( &
         xtetlist, & !intent(in)
         nlist_alloc,nlist,nmax,volcut, &
         cencut,sdim)

       call multi_phi_bounds( &
         bfact,dx, &
         xsten0,nhalf0, &
         slope, &
         xtetlist, & !intent(in) 
         nlist_alloc, &
         nlist, &
         nmax, &
         minphi,maxphi,sdim)

      else if (fastflag.eq.1) then

       volcut=volcell
       do dir=1,sdim
        cencut(dir)=cencell(dir)
       enddo

       minphi=1.0D+10
       maxphi=-1.0D+10
       null_intercept=zero

       do k=klo_stencil,khi_stencil,2
       do j=-1,1,2
       do i=-1,1,2

        dir=1
        xtarget(dir)=xsten0(i,dir)
        dir=2
        xtarget(dir)=xsten0(j,dir)

        if (sdim.eq.3) then
         dir=sdim
         xtarget(dir)=xsten0(k,dir)
        endif

        call distfunc(bfact,dx, &
         xsten0,nhalf0, &
         null_intercept,slope, &
         xtarget,dist,sdim)

        if (dist.lt.minphi) then
         minphi=dist
        endif
        if (dist.gt.maxphi) then
         maxphi=dist
        endif
        
       enddo
       enddo
       enddo  ! i,j,k

       if (((minphi.eq.zero).and. &
            (maxphi.eq.zero)).or. &
           (minphi.ge.maxphi)) then
        print *,"cannot have zero slope"
        print *,"fastflag=",fastflag
        print *,"minphi=",minphi
        print *,"maxphi=",maxphi
        print *,"slopexyz=",slope(1),slope(2),slope(sdim)
        print *,"xsten(0) xyz=",xsten0(0,1),xsten0(0,2),xsten0(0,sdim)
        stop
       else if (((minphi.ne.zero).or. &
                 (maxphi.ne.zero)).and. &
                (minphi.lt.maxphi)) then
        ! do nothing
       else
        print *,"minphi or maxphi is NaN"
        stop
       endif

      else
       print *,"fastflag invalid multi_find_intercept 2"
       stop
      endif

      do dir=1,sdim
       centroid(dir)=cencut(dir)
      enddo

      intercept_lower=-maxphi
      intercept_upper=-minphi

       ! volcut is the uncaptured volume of the cell.
      if ((volcut.le.zero).and.(vfrac.gt.EPS_14_7)) then
       print *,"ERROR: volcut<=0 and vfrac>EPS_14_7"
       stop
      else if ((volcut.gt.zero).or.(vfrac.le.EPS_14_7)) then
       ! do nothing
      else
       print *,"volcut or vfrac is NaN: ",volcut,vfrac
       stop
      endif

      if (vfrac.le.EPS_14_7) then
       intercept=intercept_lower
      else if (vfrac.ge.volcut/volcell-EPS_14_7) then
       intercept=intercept_upper
      else if ((vfrac.ge.EPS_14_7).and. &
               (vfrac.le.volcut/volcell-EPS_14_7)) then

       vtarget=volcell*vfrac

! solve f(xx)=0 where f(xx)=(V(n dot (x-x0)+intercept-xx)-Vtarget)/volcell
! maxphi -> (maxphi-minphi)vfrac 
! minphi -> (minphi-maxphi)(1-vfrac)

       vfrac_normalize=vfrac*volcell/volcut
       if ((vfrac_normalize.le.zero).or. &
           (vfrac_normalize.ge.one)) then
        print *,"ERROR: vfrac_normalize out of range"
        stop
       else if ((vfrac_normalize.gt.zero).and. &
                (vfrac_normalize.lt.one)) then
        !do nothing
       else
        print *,"vfrac_normalize invalid 9714: ",vfrac_normalize
        stop
       endif

       intercept_default=intercept_lower*(one-vfrac_normalize)+ &
           intercept_upper*vfrac_normalize

         ! fc_default=(voln-vtarget)/volcell
       call multi_ff( &
        bfact,dx, &
        xsten0,nhalf0, &
        fc_default,slope, &
        intercept_default, & !intent(in)
        arean_default,vtarget, &
        xtetlist, & !intent(in)
        nlist_alloc, &
        centroid, & !intent(out) 
        nlist, &
        nmax, &
        fastflag,sdim)
       err_default=abs(fc_default)

       if (use_initial_guess.eq.0) then

        intercept=intercept_default
        arean=arean_default
        err=err_default
        fc=fc_default

       else if (use_initial_guess.eq.1) then

        if ((intercept.lt.intercept_lower).or. &
            (intercept.gt.intercept_upper)) then

         intercept=intercept_default
         arean=arean_default
         err=err_default
         fc=fc_default

        else if ((intercept.ge.intercept_lower).and. &
                 (intercept.le.intercept_upper)) then

          ! fc=(voln-vtarget)/volcell
         call multi_ff( &
          bfact,dx, &
          xsten0,nhalf0, &
          fc,slope, &
          intercept, & !intent(in)
          arean, &
          vtarget, &
          xtetlist, & !intent(in)
          nlist_alloc, &
          centroid, & !intent(out)
          nlist, &
          nmax, &
          fastflag,sdim)
         err=abs(fc)
         if ((err.ge.err_default).or.(arean.eq.zero)) then
          intercept=intercept_default
          arean=arean_default
          err=err_default
          fc=fc_default
         endif

        else
         print *,"(breakpoint) break point and gdb: "
         print *,"(1) compile with the -g option"
         print *,"(2) break MOF.F90:9740"
         print *,"intercept is NaN (multi_find_intercept): ",intercept
         stop 
        endif

       else
        print *,"use_initial_guess invalid multi_find_intercept"
        print *,"use_initial_guess = ",use_initial_guess
        stop
       endif

        ! err=abs(fc)=(voln-vtarget)/volcell
       if (err.gt.moftol) then

        niter=0
        do while ((niter.lt.maxiter).and.(err.gt.moftol)) 

         if (min_err.eq.-one) then
          min_err=err
         else if (min_err.gt.err) then
          min_err=err
         else if ((min_err.le.err).and. &
                  (min_err.ge.zero)) then
          ! do nothing
         else
          print *,"min_err invalid"
          stop
         endif
         intercept_error_history(niter+1,tid_in+1)=err

         if (arean.ge.zero) then

          if (arean.gt.zero) then
           intercept_test=intercept-fc*volcell/arean
          else if (arean.eq.zero) then
           intercept_test=intercept_upper+1.0D+20
          else
           print *,"arean should not be negative: ",arean
           stop
          endif

          if ((intercept_test.le.intercept_lower).or. &
              (intercept_test.ge.intercept_upper).or. &
              (niter.ge.INTERCEPT_MAXITER_NEWTON)) then
           aa=intercept_lower
           bb=intercept_upper
           niter=0
           do while ((niter.lt.maxiter).and.(err.gt.moftol))

            if (niter.eq.0) then
             call multi_ff( &
              bfact,dx, &
              xsten0,nhalf0, &
              fa,slope, &
              aa, & !intent(in)
              arean,vtarget, &
              xtetlist, & !intent(in)
              nlist_alloc, &
              centroid, & !intent(out)
              nlist, &
              nmax, &
              fastflag,sdim)
             call multi_ff( &
              bfact,dx,xsten0,nhalf0, &
              fb,slope, &
              bb, & !intent(in)
              arean,vtarget, &
              xtetlist, & !intent(in)
              nlist_alloc, &
              centroid, & !intent(out)
              nlist, &
              nmax, &
              fastflag,sdim)
            else if (niter.gt.0) then
             ! do nothing
            else
             print *,"niter invalid"
             stop
            endif

             !fa=(voln-vtarget)/volcell
            if (abs(fa).le.moftol) then
             intercept=aa
             err=zero
            else if (abs(fb).le.moftol) then
             intercept=bb
             err=zero
            else if (fa*fb.lt.zero) then
             intercept=half*(aa+bb)
             call multi_ff( &
              bfact,dx, &
              xsten0,nhalf0, &
              fc,slope, &
              intercept, & !intent(in)
              arean,vtarget, &
              xtetlist, & !intent(in)
              nlist_alloc, & 
              centroid, & !intent(out)
              nlist, &
              nmax, &
              fastflag,sdim)
             err=abs(fc)
             if (fa*fc.lt.zero) then
              bb=intercept
              fb=fc
             else if (fb*fc.lt.zero) then
              aa=intercept
              fa=fc
             else if (fc.eq.zero) then
              aa=intercept
              fa=fc
             else
              print *,"fa,fb, or fc bust"
              print *,"fa=",fa
              print *,"fb=",fb
              print *,"fc=",fc
              stop
             endif
            else
             print *,"(breakpoint) break point and gdb: "
             print *,"(1) compile with the -g option"
             print *,"(2) break MOF.F90:9859"
             print *,"signs of fa and fb are inconsistent"
             stop
            endif

            niter=niter+1
            if (debug_root.eq.1) then
             print *,"bisection: niter,intercept,fc ", &
                niter,intercept,fc
            endif  
           enddo ! bisection while

           if ((niter.ge.1).and.(niter.lt.maxiter)) then
            if (err.le.moftol) then
             ! do nothing
            else
             print *,"err invalid err,moftol=",err,moftol
             stop
            endif
           else if (niter.eq.maxiter) then
            ! 10^15 < 2^x   x=15 log 10/log 2=50
            if (niter.gt.50) then
             err=zero
            endif
           else
            print *,"niter invalid"
            stop
           endif

          else if ((intercept_test.ge.intercept_lower).and. &
                   (intercept_test.le.intercept_upper).and. &
                   (niter.lt.INTERCEPT_MAXITER_NEWTON)) then
            !intercept_test=intercept-fc*volcell/arean
           intercept=intercept_test
           call multi_ff( &
            bfact,dx, &
            xsten0,nhalf0, &
            fc,slope, &
            intercept, & !intent(in)
            arean,vtarget, &
            xtetlist, & !intent(in)
            nlist_alloc, &
            centroid, & !intent(out)
            nlist, &
            nmax, &
            fastflag,sdim)
           err=abs(fc)
          else
           print *,"intercept_test or intercept_lower(upper) invalid"
           print *,"intercept_test=",intercept_test
           print *,"intercept_lower=",intercept_lower
           print *,"intercept_upper=",intercept_upper
           print *,"niter ",niter
           print *,"INTERCEPT_MAXITER_NEWTON ",INTERCEPT_MAXITER_NEWTON
           stop
          endif

          niter=niter+1
          if (debug_root.eq.1) then
           print *,"newton: niter,intercept,fc ",niter,intercept,fc
          endif  
         else
          print *,"multi_find_intercept: "
          print *,"arean should not be negative"
          print *,"fastflag=",fastflag
          print *,"use_initial_guess=",use_initial_guess
          print *,"niter,vfrac,arean,volcut,fc ",niter,vfrac,arean, &
           volcut,fc
          print *,"volcell ",volcell
          print *,"vfrac_normalize ",vfrac_normalize
          print *,"minphi ",minphi
          print *,"maxphi ",maxphi
          print *,"intercept_default ",intercept_default
          stop
         endif

         if (niter.gt.maxiter-2) then
          if (abs(min_err-err).le.moftol) then
           if ((abs(err-intercept_error_history(niter,tid_in+1)).le. &
                moftol).and. &
               (abs(err-intercept_error_history(niter-1,tid_in+1)).le.  &
                moftol)) then
            err=zero
           endif
          else if (abs(min_err-err).gt.moftol) then
           ! do nothing
          else
           print *,"min_err or err invalid"
           stop
          endif
         else if ((niter.ge.1).and.(niter.le.maxiter-2)) then
          ! do nothing
         else
          print *,"niter invalid in the newtons method"
          stop
         endif

         ! outer while statement: Newton's method:
         ! do while((niter.lt.maxiter).and.(err.gt.moftol))

        enddo 

        if ((niter.ge.maxiter).and.(err.gt.moftol)) then
         print *,"vof recon failed in multi_find_intercept"
         print *,"niter,maxiter ",niter,maxiter
         print *,"bfact ",bfact
         print *,"vfrac ",vfrac
         do dir=1,sdim
          print *,"dir, dx, xsten0(0) ",dir,dx(dir),xsten0(0,dir)
         enddo
         do dir=1,sdim
          print *,"dir,slope ",dir,slope(dir)
         enddo
         print *,"moftol ",moftol
         print *,"error history: "
         do nn=1,maxiter
          print *,"nn,tid_in,error ",nn,tid_in, &
             intercept_error_history(nn,tid_in+1)
         enddo
         stop
        else if ((niter.lt.maxiter).or.(err.le.moftol)) then
         ! do nothing
        else
         print *,"niter or err invalid"
         stop
        endif
       else if ((err.le.moftol).and.(err.ge.zero)) then
        ! do nothing
       else
        print *,"err invalid"
        stop       
       endif 
      else
       print *,"vfrac invalid vfrac= ",vfrac
       stop
      endif 

      return
      end subroutine multi_find_intercept


       ! centroid in absolute coordinate system
      subroutine single_find_intercept( &
       tid_in, &
       nDOF, & ! sdim-1
       nEQN, & ! sdim
       bfact,dx, &
       xsten0,nhalf0, &
       slope, &
       intercept, &
       vfrac, &
       centroid,sdim)

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: nDOF ! sdim-1
      integer, INTENT(in) :: nEQN ! sdim 
      integer, INTENT(in) :: bfact,nhalf0
      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: slope(sdim)
      real(amrex_real), INTENT(inout) :: intercept
      real(amrex_real), INTENT(in) :: vfrac
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)

      integer niter,maxiter
      integer i,j,k,dir
      real(amrex_real) minphi,maxphi
      real(amrex_real) volcell
      real(amrex_real), INTENT(out) :: centroid(sdim)
      real(amrex_real) cencell(sdim)
      real(amrex_real) arean
      real(amrex_real) err
      real(amrex_real) moftol
      real(amrex_real) min_err
      real(amrex_real) vtarget,fc
      integer debug_root
      real(amrex_real) err_default,fc_default,arean_default,intercept_default
      real(amrex_real) volcut
      real(amrex_real) cencut(sdim)
      real(amrex_real) xtarget(sdim)
      real(amrex_real) vfrac_normalize
      real(amrex_real) null_intercept,dist
      integer klo_stencil,khi_stencil
      integer nn
      real(amrex_real) intercept_upper,intercept_lower
      real(amrex_real) intercept_test,aa,bb,fa,fb

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid (single_find_intercept) ",tid_in
       stop
      endif 

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid single_find_intercept"
       stop
      endif

      if (sdim.eq.2) then
       klo_stencil=0
       khi_stencil=0
      else if (sdim.eq.3) then
       klo_stencil=-1
       khi_stencil=1
      else
       print *,"sdim invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid10116 ",bfact
       stop
      endif

      debug_root=0

      if (INTERCEPT_MAXITER_NEWTON.lt.INTERCEPT_MAXITER) then
       ! do nothing
      else
       print *,"INTERCEPT_MAXITER_NEWTON invalid"
       stop
      endif

      maxiter=INTERCEPT_MAXITER
       ! INTERCEPT_TOL is declared in PROBCOMMON.F90
      moftol=INTERCEPT_TOL
      min_err=-one

! phi=n dot (x-x0)+int
! find max,min n dot (x-x0)
! if fastflag=0, 
!  search the vertices of all triangles that make up the "cut" domain.

      call Box_volumeFAST( &
       bfact,dx, &
       xsten0,nhalf0, &
       volcell, &
       cencell, &
       sdim)

       ! at each node x_i one solves:
       !  n dot (x_i-x0) + b_i = 0
       !  b_i=-n dot (x_i-x0)
       !  b_Lower_bound=min b_i = - maxphi
       !  b_Upper_bound=max b_i = - minphi
      volcut=volcell
      do dir=1,sdim
       cencut(dir)=cencell(dir)
      enddo

      minphi=1.0D+10
      maxphi=-1.0D+10
      null_intercept=zero

      do k=klo_stencil,khi_stencil,2
      do j=-1,1,2
      do i=-1,1,2
       dir=1
       xtarget(dir)=xsten0(i,dir)
       dir=2
       xtarget(dir)=xsten0(j,dir)

       if (sdim.eq.3) then
        dir=sdim
        xtarget(dir)=xsten0(k,dir)
       endif

       call distfunc(bfact,dx, &
        xsten0,nhalf0, &
        null_intercept,slope, &
        xtarget,dist,sdim)

       if (dist.lt.minphi) then
        minphi=dist
       endif
       if (dist.gt.maxphi) then
        maxphi=dist
       endif
         
      enddo
      enddo
      enddo  ! i,j,k

      if (((minphi.eq.zero).and. &
           (maxphi.eq.zero)).or. &
          (minphi.ge.maxphi)) then
       print *,"cannot have zero slope"
       print *,"minphi=",minphi
       print *,"maxphi=",maxphi
       print *,"slopexyz=",slope(1),slope(2),slope(sdim)
       print *,"xsten(0) xyz=",xsten0(0,1),xsten0(0,2),xsten0(0,sdim)
       stop
      else if (((minphi.ne.zero).or. &
                (maxphi.ne.zero)).and. &
               (minphi.lt.maxphi)) then
       ! do nothing
      else
       print *,"minphi or maxphi is NaN"
       stop
      endif

      do dir=1,sdim
       centroid(dir)=cencut(dir)
      enddo

      intercept_lower=-maxphi
      intercept_upper=-minphi

       ! volcut is the uncaptured volume of the cell.
      if ((volcut.le.zero).and.(vfrac.gt.EPS_14_7)) then
       print *,"ERROR: volcut<=0 and vfrac>EPS_14_7"
       stop
      else if ((volcut.gt.zero).or.(vfrac.le.EPS_14_7)) then
       !do nothing
      else
       print *,"volcut or vfrac corrupt: ",volcut,vfrac
       stop
      endif

      if (vfrac.le.EPS_14_7) then
       intercept=intercept_lower
      else if (vfrac.ge.volcut/volcell-EPS_14_7) then
       intercept=intercept_upper
      else if ((vfrac.ge.EPS_14_7).and. &
               (vfrac.le.volcut/volcell-EPS_14_7)) then

       vtarget=volcell*vfrac

! solve f(xx)=0 where f(xx)=(V(n dot (x-x0)+intercept-xx)-Vtarget)/volcell
! maxphi -> (maxphi-minphi)vfrac 
! minphi -> (minphi-maxphi)(1-vfrac)

       vfrac_normalize=vfrac*volcell/volcut
       if ((vfrac_normalize.le.zero).or. &
           (vfrac_normalize.ge.one)) then
        print *,"ERROR: vfrac_normalize out of range"
        stop
       else if ((vfrac_normalize.gt.zero).and. &
                (vfrac_normalize.lt.one)) then
        ! do nothing
       else
        print *,"vfrac_normalize is invalid 10247: ",vfrac_normalize
        stop
       endif

       intercept_default=intercept_lower*(one-vfrac_normalize)+ &
           intercept_upper*vfrac_normalize

         ! fc_default=(voln-vtarget)/volcell
       call single_ff( &
        bfact,dx, &
        xsten0,nhalf0, &
        fc_default,slope, &
        intercept_default, &
        arean_default,vtarget, &
        centroid, &
        sdim)
       err_default=abs(fc_default)

       intercept=intercept_default
       arean=arean_default
       err=err_default
       fc=fc_default

        ! err=abs(fc)=(voln-vtarget)/volcell
       if (err.gt.moftol) then

        niter=0
        do while ((niter.lt.maxiter).and.(err.gt.moftol)) 

         if (min_err.eq.-one) then
          min_err=err
         else if (min_err.gt.err) then
          min_err=err
         else if ((min_err.le.err).and. &
                  (min_err.ge.zero)) then
          ! do nothing
         else
          print *,"min_err invalid"
          stop
         endif
         intercept_error_history(niter+1,tid_in+1)=err

         if (arean.ge.zero) then

          if (arean.gt.zero) then
           intercept_test=intercept-fc*volcell/arean
          else if (arean.eq.zero) then
           intercept_test=intercept_upper+1.0D+20
          else
           print *,"arean should not be negative: ",arean
           stop
          endif

          if ((intercept_test.le.intercept_lower).or. &
              (intercept_test.ge.intercept_upper).or. &
              (niter.ge.INTERCEPT_MAXITER_NEWTON)) then
           aa=intercept_lower
           bb=intercept_upper
           niter=0
           do while ((niter.lt.maxiter).and.(err.gt.moftol))

            if (niter.eq.0) then
             call single_ff( &
              bfact,dx, &
              xsten0,nhalf0, &
              fa,slope, &
              aa, &
              arean,vtarget, &
              centroid, &
              sdim)
             call single_ff( &
              bfact,dx,xsten0,nhalf0, &
              fb,slope, &
              bb, &
              arean,vtarget, &
              centroid, &
              sdim)
            else if (niter.gt.0) then
             ! do nothing
            else
             print *,"niter invalid"
             stop
            endif

             !fa=(voln-vtarget)/volcell
            if (abs(fa).le.moftol) then
             intercept=aa
             err=zero
            else if (abs(fb).le.moftol) then
             intercept=bb
             err=zero
            else if (fa*fb.lt.zero) then
             intercept=half*(aa+bb)
             call single_ff( &
              bfact,dx, &
              xsten0,nhalf0, &
              fc,slope, &
              intercept, &
              arean,vtarget, &
              centroid, &
              sdim)
             err=abs(fc)
             if (fa*fc.lt.zero) then
              bb=intercept
              fb=fc
             else if (fb*fc.lt.zero) then
              aa=intercept
              fa=fc
             else
              print *,"fa,fb, or fc bust"
              print *,"fa=",fa
              print *,"fb=",fb
              print *,"fc=",fc
              stop
             endif
            else
             print *,"(breakpoint) break point and gdb: "
             print *,"(1) compile with the -g option"
             print *,"(2) break MOF.F90:10252"
             print *,"signs of fa and fb are inconsistent"
             stop
            endif

            niter=niter+1
            if (debug_root.eq.1) then
             print *,"bisection(single): niter,intercept,fc ", &
                niter,intercept,fc
            endif  
           enddo ! bisection while

           if ((niter.ge.1).and.(niter.lt.maxiter)) then
            if (err.le.moftol) then
             ! do nothing
            else
             print *,"err invalid err,moftol=",err,moftol
             stop
            endif
           else if (niter.eq.maxiter) then
            ! 10^15 < 2^x   x=15 log 10/log 2=50
            if (niter.gt.50) then
             err=zero
            endif
           else
            print *,"niter invalid"
            stop
           endif

          else if ((intercept_test.ge.intercept_lower).and. &
                   (intercept_test.le.intercept_upper).and. &
                   (niter.lt.INTERCEPT_MAXITER_NEWTON)) then
            !intercept_test=intercept-fc*volcell/arean
           intercept=intercept_test
           call single_ff( &
            bfact,dx, &
            xsten0,nhalf0, &
            fc,slope, &
            intercept, &
            arean,vtarget, &
            centroid, &
            sdim)
           err=abs(fc)
          else
           print *,"intercept_test or intercept_lower(upper) invalid"
           print *,"intercept_test=",intercept_test
           print *,"intercept_lower=",intercept_lower
           print *,"intercept_upper=",intercept_upper
           print *,"niter ",niter
           print *,"INTERCEPT_MAXITER_NEWTON ",INTERCEPT_MAXITER_NEWTON
           stop
          endif

          niter=niter+1
          if (debug_root.eq.1) then
           print *,"newton(single): niter,intercept,fc ",niter,intercept,fc
          endif  
         else
          print *,"single_find_intercept: "
          print *,"arean should not be negative"
          print *,"niter,vfrac,arean,volcut,fc ",niter,vfrac,arean, &
           volcut,fc
          print *,"volcell ",volcell
          print *,"vfrac_normalize ",vfrac_normalize
          print *,"minphi ",minphi
          print *,"maxphi ",maxphi
          print *,"intercept_default ",intercept_default
          stop
         endif

         if (niter.gt.maxiter-2) then
          if (abs(min_err-err).le.moftol) then
           if ((abs(err-intercept_error_history(niter,tid_in+1)).le.moftol).and. &
               (abs(err-intercept_error_history(niter-1,tid_in+1)).le.moftol)) then
            err=zero
           endif
          else if (abs(min_err-err).gt.moftol) then
           ! do nothing
          else
           print *,"min_err or err invalid"
           stop
          endif
         else if ((niter.ge.1).and.(niter.le.maxiter-2)) then
          ! do nothing
         else
          print *,"niter invalid in the newtons method"
          stop
         endif

         ! outer while statement: Newton's method:
         ! do while((niter.lt.maxiter).and.(err.gt.moftol))

        enddo ! outer while statement: Newton's method

        if ((niter.ge.maxiter).and.(err.gt.moftol)) then
         print *,"vof recon failed in single_find_intercept"
         print *,"niter,maxiter ",niter,maxiter
         print *,"bfact ",bfact
         print *,"vfrac ",vfrac
         do dir=1,sdim
          print *,"dir, dx, xsten0(0) ",dir,dx(dir),xsten0(0,dir)
         enddo
         do dir=1,sdim
          print *,"dir,slope ",dir,slope(dir)
         enddo
         print *,"moftol ",moftol
         print *,"error history: "
         do nn=1,maxiter
          print *,"nn,tid_in,error ",nn,tid_in, &
             intercept_error_history(nn,tid_in+1)
         enddo
         stop
        else if ((niter.lt.maxiter).or.(err.le.moftol)) then
         ! do nothing
        else
         print *,"niter or err invalid"
         stop
        endif
       else if ((err.le.moftol).and.(err.ge.zero)) then
        ! do nothing
       else
        print *,"err invalid"
        stop       
       endif 
      else
       print *,"vfrac invalid vfrac= ",vfrac
       stop
      endif  

      return
      end subroutine single_find_intercept

      real(NOTUS_REAL) function atan_verify_NOTUS(x)
      IMPLICIT NONE
      real(NOTUS_REAL), INTENT(in) :: x
      real(NOTUS_REAL) NOTUS_PI

      NOTUS_PI=4.0d0*atan(1.0d0)

      atan_verify_NOTUS=atan(x)
      if ((atan_verify_NOTUS.ge.-0.5d0*NOTUS_PI-0.01d0).and. &
          (atan_verify_NOTUS.le.0.5d0*NOTUS_PI+0.01d0)) then
       !do nothing
      else
       print *,"atan out of range NOTUS"
       stop
      endif
 
      return
      end function atan_verify_NOTUS

      subroutine arctan2_NOTUS(y,x,angle)
      IMPLICIT NONE

      real(NOTUS_REAL), INTENT(in) :: y,x
      real(NOTUS_REAL), INTENT(out) :: angle
      real(NOTUS_REAL) NOTUS_PI

      NOTUS_PI=4.0d0*atan(1.0d0)

      if ((y.gt.0.0d0).and.(x.eq.0.0d0)) then
       angle=0.5d0*NOTUS_PI
      else if ((y.lt.0.0d0).and.(x.eq.0.0d0)) then
       angle=-0.5d0*NOTUS_PI
        ! pi/4 <angle<pi/2 (first quadrant)
      else if ((y.gt.0.0d0).and.(x.gt.0.0d0).and.(y.ge.x)) then
       angle=atan_verify_NOTUS(y/x)
        ! pi/2 < angle < 3pi/4 (second quadrant)
      else if ((y.gt.0.0d0).and.(x.lt.0.0d0).and.(y.ge.abs(x))) then
       angle=atan_verify_NOTUS(y/x)+NOTUS_PI
        ! -pi/2<angle<-pi/4 (fourth quadrant)
      else if ((y.lt.0.0d0).and.(x.gt.0.0d0).and.(abs(y).ge.x)) then
       angle=atan_verify_NOTUS(y/x)
        ! -3pi/4 < angle < -pi/2  (3rd quadrant)
      else if ((y.lt.0.0d0).and.(x.lt.0.0d0).and.(abs(y).ge.abs(x))) then
       angle=atan_verify_NOTUS(y/x)-NOTUS_PI
      else if ((y.eq.0.0d0).and.(x.gt.0.0d0)) then
       angle=0.0d0
      else if ((y.eq.0.0d0).and.(x.lt.0.0d0)) then
       angle=NOTUS_PI
       ! 0<angle<pi/4 (1st quadrant)
      else if ((x.gt.0.0d0).and.(y.gt.0.0d0).and.(y.le.x)) then
       angle=atan_verify_NOTUS(y/x)
       ! -pi/4<angle<0 (4th quadrant)
      else if ((x.gt.0.0d0).and.(y.lt.0.0d0).and.(abs(y).le.x)) then
       angle=atan_verify_NOTUS(y/x)
       ! 3pi/4<angle<pi (second quadrant)
      else if ((x.lt.0.0d0).and.(y.gt.0.0d0).and.(y.le.abs(y))) then
       angle=atan_verify_NOTUS(y/x)+NOTUS_PI
       ! -pi<angle<-3pi/4 (3rd quadrant)
      else if ((x.lt.0.0d0).and.(y.lt.0.0d0).and.(abs(y).le.abs(x))) then
       angle=atan_verify_NOTUS(y/x)-NOTUS_PI
      else
       angle=0.0d0
      endif

      return
      end subroutine arctan2_NOTUS


        ! "Find the actual centroid given the angle" 
        ! xcell is center of cell, not the cell centroid
        ! refcentroid is passed into this routine.
        ! refcentroid is relative to cell centroid of the cell.
      subroutine multi_rotatefunc( &
        tid_in, &
        uncaptured_volume_vof, & !intent(in)
        mofdata, &
        nDOF, & ! sdim-1 
        nEQN, & ! sdim 
        use_MilcentLemoine, &
        bfact,dx, &
        xsten0,nhalf0, &
        xtetlist_vof, & !intent(in)
        nlist_vof, & !intent(in)
        nlist_alloc, & !intent(in)
        nmax, &
        refcentroid, & !relative to cell centroid
        refvfrac, &
        angle, & !intent(in)
        ff, & !intent(out)
        intercept, & !intent(inout)
        testcen, &  !intent(out)
        use_initial_guess, &
        fastflag,sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module
      use mod_mof3d_tetra_analytic_centroid
      use mod_mof3d_analytic_centroid
      use mod_mof2d_analytic_centroid

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: nlist_alloc

      real(amrex_real), INTENT(in) :: uncaptured_volume_vof

      real(amrex_real), INTENT(inout) :: mofdata(num_materials*ngeom_recon)

      integer, INTENT(in) :: nDOF ! sdim-1
      integer, INTENT(in) :: nEQN ! sdim 
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: use_MilcentLemoine
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: nhalf0
      integer,INTENT(in):: nlist_vof

      integer, INTENT(in) :: nmax,fastflag
      real(amrex_real), INTENT(in) :: xtetlist_vof(4,3,nlist_alloc)

      real(amrex_real), INTENT(in) :: refcentroid(nEQN)
      real(amrex_real) refcentroidT(nEQN)
      real(amrex_real), INTENT(in) :: refvfrac
      integer, INTENT(in) :: use_initial_guess 

      real(amrex_real), INTENT(out) :: testcen(nEQN)
      real(amrex_real) testcenT(nEQN)
      real(amrex_real), INTENT(in) :: angle(nDOF)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real) xsten_local(-nhalf0:nhalf0,sdim)
      integer isten
      integer, PARAMETER :: nhalf2=1
      real(amrex_real), INTENT(out) :: ff(nEQN) 
      real(amrex_real) volume_cut,facearea
      integer dir
      integer dir_local
      real(amrex_real) :: nslope(nEQN)
      real(amrex_real), INTENT(inout) :: intercept

      real(amrex_real) volcell_vof
      real(amrex_real) cencell_vof(sdim)
      real(amrex_real) xtet(sdim+1,sdim)
      integer, PARAMETER :: shapeflag=0 !regular hexahedron
      integer ksten_low,ksten_high
      real(NOTUS_REAL) local_angles_NOTUS(2)
      real(NOTUS_REAL) local_volume_NOTUS
      real(NOTUS_REAL) local_cell_size_NOTUS(3)
      real(NOTUS_REAL) local_nslope_NOTUS(3)
      real(NOTUS_REAL) local_ref_centroid_NOTUS(3)
      real(NOTUS_REAL) local_centroid_NOTUS(3)
      real(NOTUS_REAL) local_gradient_NOTUS(2)
      real(NOTUS_REAL) local_refvfrac_NOTUS
      real(NOTUS_REAL) NOTUS_PI

      integer :: NOTUS_worked

      NOTUS_PI=4.0d0*atan(1.0d0)

      if ((tid_in.lt.0).or.(tid_in.ge.geom_nthreads)) then
       print *,"tid_in invalid (multi_rotatefunc) ",tid_in
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_rotatefunc"
       stop
      endif
      if (nhalf0.lt.3) then
       print *,"nhalf0 invalid"
       stop
      endif

      if ((nlist_vof.ge.0).and.(nlist_vof.le.nlist_alloc)) then
       ! do nothing
      else
       print *,"nlist_vof invalid"
       print *,"nlist_vof: ",nlist_vof
       print *,"nlist_alloc: ",nlist_alloc
       print *,"nmax: ",nmax
       stop
      endif
      if (fastflag.eq.1) then
       if (nlist_vof.eq.0) then
        !do nothing
       else
        print *,"(breakpoint) break point and gdb: "
        print *,"(1) compile with the -g option"
        print *,"(2) break MOF.F90:10852"
        print *,"expecting nlist_vof=0"
        print *,"nlist_vof: ",nlist_vof
        print *,"nlist_alloc: ",nlist_alloc
        print *,"nmax: ",nmax
        print *,"fastflag: ",fastflag
        stop
       endif
      else if (fastflag.eq.0) then
       if ((nlist_vof.ge.1).and.(nlist_vof.le.nlist_alloc)) then
        !do nothing
       else
        print *,"expecting nlist_vof >= 1"
        stop
       endif
      else
       print *,"fastflag invalid"
       stop
      endif

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      if ((nDOF.eq.sdim-1).and. &
          (nEQN.eq.sdim)) then
       ! do nothing
      else
       print *,"invalid nDOF or nEQN"
       stop
      endif

      if ((use_initial_guess.ne.0).and. &
          (use_initial_guess.ne.1)) then
       print *,"use_initial_guess invalid multirotatefunc"
       stop
      endif

      if (uncaptured_volume_vof.gt.zero) then
       ! do nothing
      else
       print *,"uncaptured_volume_vof invalid"
       stop
      endif

      if (refvfrac**2.ge.zero) then
       !do nothing
      else
       print *,"corrupt: refvfrac ",refvfrac
       stop
      endif

      do dir=1,sdim
       if (refcentroid(dir)**2.ge.zero) then
        !do nothing
       else
        print *,"corrupt: dir,refcentroid(dir) ",dir,refcentroid(dir)
        stop
       endif
      enddo

      do dir=1,nEQN
       ff(dir)=zero
      enddo

       ! if RZ, cencell can be negative

      call Box_volumeFAST( &
       bfact,dx,xsten0,nhalf0, &
       volcell_vof, &
       cencell_vof, &
       sdim)

      call angle_to_slope( &
        angle, &
        nslope, &
        sdim)

      if (sdim.eq.3) then
       ksten_low=-1
       ksten_high=1
      else if (sdim.eq.2) then
       ksten_low=0
       ksten_high=0
      else
       print *,"sdim invalid"
       stop
      endif

      if (use_MilcentLemoine.eq.0) then

        ! inside of multi_rotatefunc
        ! testcen in absolute coordinate system
        ! (testcen is the centroid of the intersection of
        !  the material region with the center cell 
       call multi_find_intercept( &
        tid_in, &
        nDOF, & ! sdim-1
        nEQN, & ! sdim 
        bfact,dx, &
        xsten0,nhalf0, &
        nslope, &  ! intent(in)
        intercept, & ! intent(inout)
        xtetlist_vof, & !intent(in)
        nlist_alloc, &
        nlist_vof, &
        nmax, &
        refvfrac, &
        use_initial_guess, &
        testcen, & !intent(out)
        fastflag,sdim)

         ! testcen in absolute coordinate system
       if (fastflag.eq.0) then
          ! xsten_local used for LS dist.
        do isten=-nhalf0,nhalf0
        do dir=1,sdim
         xsten_local(isten,dir)=xsten0(isten,dir)
        enddo
        enddo

          ! (testcen is the centroid of the intersection of
          !  the material region with the center cell)
          ! calling from "mult_rotatefunc"
        call multi_cell_intersection( &
          bfact,dx, &
          xsten_local,nhalf0, &
          nslope, &
          intercept, &
          volume_cut, & !intent(out)
          testcen, & !intent(out)
          facearea, & !intent(out)
          xtetlist_vof, & !intent(in)
          nlist_alloc, &
          nlist_vof, &
          nmax, &
          sdim)

       else if (fastflag.eq.1) then

          ! (testcen is the centroid of the intersection of
          !  the material region with the center cell)
        call fast_cut_cell_intersection( &
          bfact,dx,xsten0,nhalf0, &
          nslope, &
          intercept, &
          volume_cut,testcen, &
          facearea, &
          xsten0,nhalf0, &
          xtet,shapeflag, &
          sdim) 

       else
        print *,"fastflag invalid multi rotatefunc: ",fastflag
        stop
       endif

       do dir=1,sdim
        testcen(dir)=testcen(dir)-cencell_vof(dir)
       enddo

        !refcentroid and testcen are both "local"
        !refcentroidT and testcenT are both "global"
       call RT_transform_offset(refcentroid,cencell_vof,refcentroidT)
       call RT_transform_offset(testcen,cencell_vof,testcenT)

      else if (use_MilcentLemoine.eq.1) then

       if (fastflag.eq.1) then
        !do nothing
       else
        print *,"expecting fastflag=1 if use_MilcentLemoine==1"
        stop
       endif

       intercept=zero

        ! MilcentLemoine slope points away from the material,
        ! so we have to reverse the normal and angles
        ! nslope comes from "angle_to_slope"
       do dir=1,sdim
        local_nslope_NOTUS(dir)=-nslope(dir)
       enddo
       call slope_to_angle_NOTUS(local_nslope_NOTUS,local_angles_NOTUS,sdim)

       local_refvfrac_NOTUS=refvfrac
       do dir=1,sdim
        local_ref_centroid_NOTUS(dir)=refcentroid(dir)
       enddo

       if ((refvfrac.ge.0.5d0).and. &
           (refvfrac.lt.1.0d0)) then
        local_refvfrac_NOTUS=1.0d0-local_refvfrac_NOTUS

        if (local_refvfrac_NOTUS.gt.0.0d0) then
         do dir=1,sdim
          local_ref_centroid_NOTUS(dir)= &
            -refvfrac*refcentroid(dir)/local_refvfrac_NOTUS
          local_nslope_NOTUS(dir)=-local_nslope_NOTUS(dir)
         enddo
        else if (local_refvfrac_NOTUS.eq.0.0d0) then
         do dir=1,sdim
          local_ref_centroid_NOTUS(dir)=0.0d0
          local_nslope_NOTUS(dir)=-local_nslope_NOTUS(dir)
         enddo
        else
         print *,"local_refvfrac_NOTUS corrupt: ",local_refvfrac_NOTUS
         stop
        endif

        call slope_to_angle_NOTUS(local_nslope_NOTUS, &
                local_angles_NOTUS,sdim)

       else if ((refvfrac.gt.0.0d0).and. &
                (refvfrac.le.0.5d0)) then
        ! do nothing
       else
        print *,"local_refvfrac_NOTUS invalid: ",local_refvfrac_NOTUS
        print *,"refvfrac Invalid: ",refvfrac
        stop
       endif

       local_volume_NOTUS=local_refvfrac_NOTUS*volcell_vof

       do dir=1,sdim
        local_cell_size_NOTUS(dir)=xsten0(1,dir)-xsten0(-1,dir)
        if (local_cell_size_NOTUS(dir).gt.0.0d0) then
         if ((local_refvfrac_NOTUS.gt.0.0d0).and. &
             (local_refvfrac_NOTUS.le.0.5d0)) then
          local_ref_centroid_NOTUS(dir)= &
               local_ref_centroid_NOTUS(dir)+ &
               0.5d0*local_cell_size_NOTUS(dir)
         else
          print *,"local_refvfrac_NOTUS invalid: ",local_refvfrac_NOTUS
          stop
         endif
        else
         print *,"local_cell_size_NOTUS(dir) invalid: ", &
                 dir,local_cell_size_NOTUS(dir)
         stop
        endif 
       enddo ! dir=1..sdim

       if (sdim.eq.3) then
        ! do nothing
       else if (sdim.eq.2) then
        local_angles_NOTUS(2)=0.5d0*NOTUS_PI
        local_nslope_NOTUS(3)=0.0d0
        local_cell_size_NOTUS(3)=0.0d0
        local_ref_centroid_NOTUS(3)=0.5d0
       else
        print *,"sdim invalid"
        stop
       endif

       if (sdim.eq.3) then
        call mof3d_compute_analytic_gradient( &
         local_angles_NOTUS, & !intent(in)
         local_ref_centroid_NOTUS, & !intent(in) relative to LL corner.
         local_volume_NOTUS, & !intent(in)
         local_cell_size_NOTUS, & !intent(in)
         local_centroid_NOTUS, & !intent(out)
         local_gradient_NOTUS) !intent(out) relative to LL corner
       else if (sdim.eq.2) then
        call mof2d_compute_analytic_gradient( &
         local_angles_NOTUS, & !intent(in)
         local_volume_NOTUS, &
         local_cell_size_NOTUS, &
         local_centroid_NOTUS)  !intent(out) relative to LL corner
       else
        print *,"sdim invalid"
        stop
       endif

       NOTUS_worked=1

       do dir=1,sdim
        local_centroid_NOTUS(dir)= &
                local_centroid_NOTUS(dir)- &
                0.5d0*local_cell_size_NOTUS(dir)
        testcen(dir)=local_centroid_NOTUS(dir) !relative of cell centroid !
        if (testcen(dir)**2.ge.zero) then
         !do nothing
        else
         print *,"WARNING"
         print *,"testcen(dir) corrupt (multi_rotatefunc): ",dir,testcen(dir)
         print *,"local_centroid_NOTUS(dir): ",dir, &
                 local_centroid_NOTUS(dir)
         print *,"volcell_vof: ",volcell_vof
         print *,"local_volume_NOTUS: ",local_volume_NOTUS
         print *,"refvfrac=",refvfrac
         print *,"local_refvfrac_NOTUS=",local_refvfrac_NOTUS
         print *,"refcentroid(dir)=",dir,refcentroid(dir)
         print *,"local_ref_centroid_NOTUS(dir)=", &
                 dir,local_ref_centroid_NOTUS(dir)
         print *,"uncaptured_volume_vof=",uncaptured_volume_vof
         print *,"dx(dir)=",dir,dx(dir)
         print *,"bfact=",bfact
         print *,"local_angles_NOTUS: ", &
           local_angles_NOTUS(1),local_angles_NOTUS(2)
         print *,"local_volume_NOTUS: ",local_volume_NOTUS
         do dir_local=1,sdim
          print *,"dir_local,local_ref_centroid_NOTUS ", &
            dir_local,local_ref_centroid_NOTUS(dir_local)
          print *,"dir_local,local_cell_size_NOTUS ", &
            dir_local,local_cell_size_NOTUS(dir_local)
         enddo !dir_local=1,sdim
         print *,"local_gradient_NOTUS: ",local_gradient_NOTUS
!        stop
         NOTUS_worked=0
        endif
       enddo !dir=1,sdim

       if (NOTUS_worked.eq.1) then

        if ((refvfrac.ge.0.5d0).and. &
            (refvfrac.lt.1.0d0)) then
         do dir=1,sdim
          testcen(dir)=-local_refvfrac_NOTUS*testcen(dir)/refvfrac
         enddo
        else if ((refvfrac.gt.0.0d0).and. &
                 (refvfrac.le.0.5d0)) then
         ! do nothing
        else
         print *,"refvfrac invalid: ",refvfrac
         stop
        endif

       else if (NOTUS_worked.eq.0) then

        ! inside of multi_rotatefunc
        ! testcen in absolute coordinate system
        ! (testcen is the centroid of the intersection of
        !  the material region with the center cell 
        call multi_find_intercept( &
         tid_in, &
         nDOF, & ! sdim-1
         nEQN, & ! sdim 
         bfact,dx, &
         xsten0,nhalf0, &
         nslope, &  ! intent(in)
         intercept, & ! intent(inout)
         xtetlist_vof, & !intent(in)
         nlist_alloc, &
         nlist_vof, &
         nmax, &
         refvfrac, &
         use_initial_guess, &
         testcen, & !intent(out)
         fastflag,sdim)

        ! testcen in absolute coordinate system
        if (fastflag.eq.0) then
         print *,"expecting fastflag=1"
         stop
        else if (fastflag.eq.1) then

         ! (testcen is the centroid of the intersection of
         !  the material region with the center cell)
         call fast_cut_cell_intersection( &
           bfact,dx,xsten0,nhalf0, &
           nslope, &
           intercept, &
           volume_cut,testcen, &
           facearea, &
           xsten0,nhalf0, &
           xtet,shapeflag, &
           sdim) 

        else
         print *,"fastflag invalid multi rotatefunc: ",fastflag
         stop
        endif

         ! testcen relative to cell centroid
        do dir=1,sdim
         testcen(dir)=testcen(dir)-cencell_vof(dir)
        enddo

        if (levelrz.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else 
         print *,"expecting levelrz.eq.COORDSYS_CARTESIAN: ",levelrz
         stop
        endif

       else
        print *,"NOTUS_worked invalid"
        stop
       endif

        ! refcentroid is relative to the cell centroid
        ! testcen is relative to the cell centroid
       do dir=1,sdim
        refcentroidT(dir)=refcentroid(dir)
        testcenT(dir)=testcen(dir)
       enddo

      else
       print *,"use_MilcentLemoine invalid: ",use_MilcentLemoine
       stop
      endif

      do dir=1,nEQN
        !either both refcentroidT and testcenT are global or both are
        !local.
       ff(dir)=(refcentroidT(dir)-testcenT(dir))
       if (ff(dir)**2.ge.zero) then
        ! do nothing
       else
        print *,"ff(dir) bust"
        print *,"use_MilcentLemoine=",use_MilcentLemoine
        print *,"dir,ff(dir) ",dir,ff(dir)
        do dir_local=1,nEQN
         print *,"dir_local,refcentroidT,testcenT ", &
            refcentroidT(dir_local),testcenT(dir_local)
        enddo
        stop
       endif
      enddo ! do dir=1,nEQN

      return
      end subroutine multi_rotatefunc


      subroutine advance_angle(angle,delangle)
      use global_utility_module
      IMPLICIT NONE
      real(amrex_real), INTENT(inout) :: angle
      real(amrex_real), INTENT(in) :: delangle

      if (MOF_PI.eq.zero) then
       MOF_PI=four*atan(one)
      endif

      angle=angle+delangle
      if (angle.lt.-MOF_PI) then
       angle=angle+two*MOF_PI
      endif
      if (angle.gt.MOF_PI) then
       angle=angle-two*MOF_PI
      endif

      return
      end subroutine advance_angle


      subroutine slope_to_angle_NOTUS(nn,angle,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer dir
      real(NOTUS_REAL), INTENT(in) :: nn(sdim)
      real(NOTUS_REAL), INTENT(out) :: angle(sdim-1)
      real(NOTUS_REAL) local_nn(sdim)
      real(NOTUS_REAL) NOTUS_PI
      real(NOTUS_REAL) mag
      real(NOTUS_REAL) x,y,z

      NOTUS_PI=4.0d0*atan(1.0d0)

      mag=0.0d0
      do dir=1,sdim
       mag=mag+nn(dir)**2
      enddo
      mag=sqrt(mag)
      if (mag.gt.0.0d0) then
       !do nothing
      else
       print *,"slope_to_angle_NOTUS: invalid slope mag=",mag
       do dir=1,sdim
        print *,"dir,nn(dir) ",dir,nn(dir)
       enddo
       stop
      endif
      do dir=1,sdim
       local_nn(dir)=nn(dir)/mag
      enddo
      x=local_nn(1)
      y=local_nn(2)
      z=local_nn(sdim)

      if (sdim.eq.3) then

        ! tan(phi)=y/x   sin(phi)/cos(phi)=y/x
       call arctan2_NOTUS(y,x,angle(1))
       if ((y.eq.0.0d0).and.(x.eq.0.0d0)) then
        if (z.eq.0.0d0) then
         print *,"z cannot be 0.0d0"
         stop
        else if (z.gt.0.0d0) then
         angle(sdim-1)=0.0d0
        else if (z.lt.0.0d0) then
         angle(sdim-1)=NOTUS_PI
        else
         print *,"bust slope_to_angle_NOTUS"
         print *,"x,y,z ",x,y,z
         print *,"sdim=",sdim
         stop
        endif
       else if (abs(x).ge.abs(y)) then
        ! tan(theta)=x/(z cos(phi))
        ! sin(theta)/cos(theta)=x/(z cos(phi))
        ! z=cos(theta)  x=sin(theta)cos(phi) y=x tan(phi)=sin(theta)sin(phi)
        call arctan2_NOTUS(x/cos(angle(1)),z,angle(sdim-1))
       else if (abs(x).le.abs(y)) then
        ! tan(theta)=y/(z sin(phi)
        ! sin(theta)/cos(theta)=y/(z sin(phi))
        ! z=cos(theta) y=sin(theta)sin(phi)  x=y/tan(phi)=cos(phi)sin(theta)
        call arctan2_NOTUS(y/sin(angle(1)),z,angle(sdim-1))
       else
        print *,"x or y bust: ",x,y
        stop
       endif 

      else if (sdim.eq.2) then
        ! tan(phi)=y/x   sin(phi)/cos(phi)=y/x
       call arctan2_NOTUS(y,x,angle(1))
 
      else
       print *,"slope_to_angle_NOTUS: sdim invalid"
       stop
      endif

      return
      end subroutine slope_to_angle_NOTUS



! angle=(phi,theta)
! 3D: x=cos(phi)sin(theta)  y=sin(phi)sin(theta)  z=cos(theta)
! 2D: (cos(phi),sin(phi)) 
! arctan2: -pi < angle < pi (declared in GLOBALUTIL.F90)
      subroutine slope_to_angle(nn,angle,sdim)

      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer dir
      real(amrex_real), INTENT(in) :: nn(sdim)
      real(amrex_real), INTENT(out) :: angle(sdim-1)
      real(amrex_real) local_nn(sdim)
      real(amrex_real) mag
      real(amrex_real) x,y,z

      if (MOF_PI.eq.zero) then
       MOF_PI=four*atan(one)
      endif

      mag=zero
      do dir=1,sdim
       mag=mag+nn(dir)**2
      enddo
      mag=sqrt(mag)
      if (mag.gt.zero) then
       !do nothing
      else
       print *,"slope_to_angle: invalid slope mag=",mag
       do dir=1,sdim
        print *,"dir,nn(dir) ",dir,nn(dir)
       enddo
       stop
      endif
      do dir=1,sdim
       local_nn(dir)=nn(dir)/mag
      enddo
      x=local_nn(1)
      y=local_nn(2)
      z=local_nn(sdim)

      if (sdim.eq.3) then

        ! tan(phi)=y/x   sin(phi)/cos(phi)=y/x
       call arctan2(y,x,angle(1))
       if ((y.eq.zero).and.(x.eq.zero)) then
        if (z.eq.zero) then
         print *,"z cannot be zero"
         stop
        else if (z.gt.zero) then
         angle(sdim-1)=zero
        else if (z.lt.zero) then
         angle(sdim-1)=MOF_PI
        else
         print *,"bust slope_to_angle"
         print *,"x,y,z ",x,y,z
         print *,"sdim=",sdim
         stop
        endif
       else if (abs(x).ge.abs(y)) then
        ! tan(theta)=x/(z cos(phi))
        ! sin(theta)/cos(theta)=x/(z cos(phi))
        ! z=cos(theta)  x=sin(theta)cos(phi) y=x tan(phi)=sin(theta)sin(phi)
        call arctan2(x/cos(angle(1)),z,angle(sdim-1))
       else if (abs(x).le.abs(y)) then
        ! tan(theta)=y/(z sin(phi)
        ! sin(theta)/cos(theta)=y/(z sin(phi))
        ! z=cos(theta) y=sin(theta)sin(phi)  x=y/tan(phi)=cos(phi)sin(theta)
        call arctan2(y/sin(angle(1)),z,angle(sdim-1))
       else
        print *,"x or y bust: ",x,y
        stop
       endif 

      else if (sdim.eq.2) then
        ! tan(phi)=y/x   sin(phi)/cos(phi)=y/x
       call arctan2(y,x,angle(1))
 
      else
       print *,"slope_to_angle: sdim invalid"
       stop
      endif

      return
      end subroutine slope_to_angle


        ! refcentroid and centroidA relative to cell centroid of the
        ! cell.
        ! xsten0(0,dir) is center of cell, not the cell centroid
        ! output: intercept,centroidA,nslope
        ! called from: individual_MOF and multimaterial_MOF
      subroutine find_cut_geom_slope( &
        tid_in, &
        uncaptured_volume_vof, &
        mofdata, &
        ls_mof, &
        lsnormal, &
        lsnormal_valid, &
        override_normal, &
        override_normal_valid, &
        bfact,dx, &
        xsten0,nhalf0, &
        refcentroid, & ! relative to cell centroid of the cell.
        refvfrac, &
        npredict, & !intent(in)
        nslope, &
        intercept, &
        xtetlist_vof, & !intent(in)
        nlist_vof, &  !intent(in)
        nlist_alloc, & !intent(in)
        centroidA, &
        nmax, &
        critical_material, & !INTENT(in)
        fastflag, &
        sdim, &
        nDOF, & ! sdim-1
        nEQN)   ! sdim 

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

#include "mofdata.H"

      integer, INTENT(in) :: tid_in

      integer, INTENT(in) :: nDOF ! sdim-1
      integer, INTENT(in) :: nEQN ! sdim 

      real(amrex_real), INTENT(in) :: uncaptured_volume_vof
      real(amrex_real), INTENT(inout) :: mofdata(num_materials*ngeom_recon)

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: bfact,nhalf0
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: nlist_vof
      integer, INTENT(in) :: nmax
      integer, INTENT(in) :: critical_material
      integer, INTENT(in) :: fastflag
      real(amrex_real), INTENT(in) :: xtetlist_vof(4,3,nlist_alloc)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: refcentroid(nEQN)
      real(amrex_real), INTENT(in) :: refvfrac
      real(amrex_real), INTENT(in) :: npredict(nEQN)
      real(amrex_real), INTENT(out) :: intercept
      real(amrex_real), INTENT(out) :: nslope(nEQN) 

      real(amrex_real) new_angle(nDOF)
      real(amrex_real) angle_opt(nDOF)

      real(amrex_real) angle_base(nDOF)
      real(amrex_real) angle_plus(nDOF)
      real(amrex_real) angle_minus(nDOF)
      real(amrex_real) angle_plus_archive(nDOF,nDOF)
      real(amrex_real) angle_minus_archive(nDOF,nDOF)
      real(amrex_real) angle_gauss_newton(nDOF)
      real(amrex_real) angle_steepest(nDOF)

      real(amrex_real) angle_init(nDOF)

      real(amrex_real) intercept_init
      real(amrex_real) cen_derive_init(nEQN)

      real(amrex_real) intercept_array(MOFITERMAX+1)
      real(amrex_real) cen_array(nEQN,MOFITERMAX+1)
      real(amrex_real) angle_array(nDOF,MOFITERMAX+1)
      real(amrex_real) f_array(nEQN,MOFITERMAX+1)  
      real(amrex_real) err_array(MOFITERMAX+1)
      real(amrex_real) delangle_array(MOFITERMAX+1)
      real(amrex_real) delta_theta_local

      integer dir,iter
      real(amrex_real) best_error
      real(amrex_real) finit(nEQN)
      real(amrex_real) fp(nEQN)
      real(amrex_real) fm(nEQN)
      real(amrex_real) fopt(nEQN)
      real(amrex_real) f_gauss_newton(nEQN)
      real(amrex_real) f_steepest(nEQN)
      real(amrex_real) fbase(nEQN)
      real(amrex_real) f_plus(nEQN,nDOF)
      real(amrex_real) f_minus(nEQN,nDOF)

      real(amrex_real) local_int

      real(amrex_real) intp(nDOF)
      real(amrex_real) intm(nDOF)
      real(amrex_real) int_gauss_newton
      real(amrex_real) int_steepest
      real(amrex_real) cen_gauss_newton(nEQN)
      real(amrex_real) cen_steepest(nEQN)
      real(amrex_real) intopt
      real(amrex_real) cenp(nEQN)
      real(amrex_real) cenm(nEQN)
      real(amrex_real) cenopt(nEQN)
      real(amrex_real) cen_plus(nEQN,nDOF)
      real(amrex_real) cen_minus(nEQN,nDOF)

      real(amrex_real) delta_theta
      real(amrex_real) delta_theta_max
      real(amrex_real) err
      real(amrex_real) err_gauss_newton
      real(amrex_real) err_steepest
      real(amrex_real) fgrad(nEQN,nDOF)  
      integer ii,iicrit
      integer i_angle,j_angle
      real(amrex_real) delangle(nDOF)
      real(amrex_real) RHS(nDOF)
      real(amrex_real) JTJ(nDOF,nDOF)
      real(amrex_real) JTJINV(nDOF,nDOF)
      real(amrex_real) tol,local_tol,DET
      real(amrex_real) err_local_min
      real(amrex_real) errinit

      real(amrex_real) err_plus(nDOF)
      real(amrex_real) err_minus(nDOF)

      real(amrex_real), INTENT(out) :: centroidA(nEQN)
      integer use_initial_guess
      real(amrex_real) dx_normalize
      integer nguess
      real(amrex_real) nLS(sdim)
      real(amrex_real) magLS

      real(amrex_real), INTENT(in) ::  &
         ls_mof(D_DECL(-1:1,-1:1,-1:1),num_materials)

      real(amrex_real), INTENT(in) :: lsnormal(num_materials,sdim)
      integer, INTENT(in) :: lsnormal_valid(num_materials)

      real(amrex_real), INTENT(in) :: override_normal(num_materials,sdim)
      integer, INTENT(in) :: override_normal_valid(num_materials)

      integer use_only_override_data

      integer local_MOFITERMAX

      integer :: use_MilcentLemoine

      if ((tid_in.lt.0).or.(tid_in.ge.geom_nthreads)) then
       print *,"tid_in invalid find_cut_geom_slope: ",tid_in
       stop
      endif

      if (MOF_PI.eq.zero) then
       MOF_PI=four*atan(one)
      endif

      if (nhalf0.lt.1) then
       print *,"expecting nhalf0>=1: ",nhalf0
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid135"
       stop
      endif
      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif
      if (nlist_vof.le.nlist_alloc) then
       !do nothing
      else
       print *,"nlist_vof corrupt"
       print *,"nlist_vof=",nlist_vof
       print *,"nlist_alloc=",nlist_alloc
       stop
      endif

      if ((MOFITERMAX.lt.num_materials+3).or. &
          (MOFITERMAX.gt.MOFITERMAX_LIMIT)) then
       print *,"MOFITERMAX out of range find cut geom slope"
       print *,"MOFITERMAX: ",MOFITERMAX
       print *,"MOFITERMAX_LIMIT: ",MOFITERMAX_LIMIT
       print *,"num_materials: ",num_materials

       stop
      endif

      local_MOFITERMAX=MOFITERMAX

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid find_cut_geom_slope"
       stop
      endif 

      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid find cut geom slope: ",num_materials
       stop
      endif
      if ((critical_material.lt.1).or. &
          (critical_material.gt.num_materials)) then
       print *,"critical_material invalid 11653: ",critical_material
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"expecting nhalf0>=1: ",nhalf0
       stop
      endif

      if (uncaptured_volume_vof.gt.zero) then
       ! do nothing
      else
       print *,"uncaptured_volume_vof invalid: ",uncaptured_volume_vof
       stop
      endif

      if ((nDOF.eq.sdim-1).and. &
          (nEQN.eq.sdim)) then
       ! do nothing
      else
       print *,"invalid nDOF or nEQN:",nDOF,nEQN
       stop
      endif

      if ((fastflag.eq.1).and. &
          (levelrz.eq.COORDSYS_CARTESIAN)) then
       use_MilcentLemoine=1 ! "1" option available here.
      else if ((fastflag.eq.0).or. &
               (levelrz.eq.COORDSYS_RZ)) then
       use_MilcentLemoine=0 ! "0" only allowed here.
      else
       print *,"parameters invalid"
       print *,"fastflag ",fastflag
       print *,"levelrz ",levelrz
       stop
      endif

      dx_normalize=dx(1)
      if (dx_normalize.gt.one) then
       dx_normalize=one
      else if ((dx_normalize.gt.zero).and. &
               (dx_normalize.le.one)) then
       ! do nothing
      else
       print *,"dx_normalize invalid: ",dx_normalize
       stop
      endif
  
      tol=dx(1)*EPS_10_5
      local_tol=dx(1)*tol*EPS2

      nguess=0
      use_only_override_data=0

      if (override_normal_valid(critical_material).eq.1) then

       use_only_override_data=1

       magLS=zero
       do dir=1,sdim
        nLS(dir)=override_normal(critical_material,dir)
        magLS=magLS+nLS(dir)**2
       enddo
       magLS=sqrt(magLS)
        ! -pi < angle < pi
       if ((magLS.ge.one-EPS2).and. &
           (magLS.le.one+EPS2)) then
        call slope_to_angle(nLS,angle_init,sdim)
       else
        print *,"magLS invalid (override): ",magLS
        do dir=1,sdim
         print *,"dir,nLS (override) ",dir,nLS(dir)
        enddo
        stop
       endif

       nguess=nguess+1 

       do dir=1,sdim-1
        angle_array(dir,nguess)=angle_init(dir)
       enddo

      else if (override_normal_valid(critical_material).eq.0) then

       use_only_override_data=0

       if (lsnormal_valid(critical_material).eq.1) then
        magLS=zero
        do dir=1,sdim
         nLS(dir)=lsnormal(critical_material,dir)
         magLS=magLS+nLS(dir)**2
        enddo
        magLS=sqrt(magLS)
         ! -pi < angle < pi
        if (magLS.gt.zero) then
         call slope_to_angle(nLS,angle_init,sdim)
        else
         print *,"magLS invalid: ",magLS
         do dir=1,sdim
          print *,"dir,nLS ",dir,nLS(dir)
         enddo
         stop
        endif
        nguess=nguess+1 
        do dir=1,sdim-1
         angle_array(dir,nguess)=angle_init(dir)
        enddo
       else if (lsnormal_valid(critical_material).eq.0) then
        ! do nothing
       else
        print *,"LSNORMAL_valid invalid1"
        print *,"critical_material,flag ",critical_material, &
         lsnormal_valid(critical_material) 
        stop
       endif

        ! -pi < angle < pi
       call slope_to_angle(npredict,angle_init,sdim)
       nguess=nguess+1
       do dir=1,nDOF
        angle_array(dir,nguess)=angle_init(dir)
       enddo

       if (num_materials+4.le.MOFITERMAX) then
        ! do nothing
       else
        print *,"no room for initial guess"
        print *,"MOFITERMAX ",MOFITERMAX
        print *,"num_materials ",num_materials
        stop
       endif

       if (sdim.eq.2) then
        if (nDOF.eq.1) then
         ! do nothing
        else
         print *,"nDOF invalid"
         stop
        endif
       else if (sdim.eq.3) then
        if (nDOF.eq.2) then
         ! do nothing
        else
         print *,"nDOF invalid"
         stop
        endif
       else
        print *,"sdim invalid" 
        stop
       endif

      else
       print *,"override_normal_valid incorrect: ",override_normal_valid
       print *,"critical_material ",critical_material
       stop
      endif

      iicrit=0
      do iter=1,nguess 

       do dir=1,nDOF
        angle_init(dir)=angle_array(dir,iter)
       enddo

       ! find finit=xref-xact for cut domain cut by a line.
       use_initial_guess=0

       intercept_init=zero

        ! calling from:
        !   find_cut_geom_slope
       call multi_rotatefunc( &
        tid_in, &
        uncaptured_volume_vof, &
        mofdata, &
        nDOF, & ! sdim-1
        nEQN, & ! sdim
        use_MilcentLemoine, &
        bfact,dx, &
        xsten0,nhalf0, &
        xtetlist_vof, & !intent(in)
        nlist_vof, & !intent(in)
        nlist_alloc, & !intent(in)
        nmax, &
        refcentroid, &
        refvfrac, &
        angle_init, &
        finit, &
        intercept_init, &
        cen_derive_init, &
        use_initial_guess, &
        fastflag,sdim)

       errinit=zero
       do dir=1,nEQN
        errinit=errinit+finit(dir)**2
       enddo
       errinit=sqrt(errinit)
       err=errinit
       best_error=errinit*refvfrac

       do dir=1,nEQN
        f_array(dir,iter)=finit(dir)
       enddo
       err_array(iter)=err
       delangle_array(iter)=zero

       if (MOFDEB.eq.1) then
        print *,"initial guess iter,angle,del,err ",iter, &
          angle_array(1,iter),angle_array(2,iter), &
          delangle_array(iter), &
          err_array(iter)
       endif

       intercept_array(iter)=intercept_init

       do dir=1,nEQN
        cen_array(dir,iter)=cen_derive_init(dir)
       enddo 
       if (iicrit.eq.0) then
        iicrit=iter
       else if (err.lt.err_array(iicrit)) then
        iicrit=iter
       endif
      enddo ! iter=1..nguess

      if ((iicrit.lt.1).or.(iicrit.gt.nguess)) then
       print *,"iicrit invalid:",iicrit
       stop
      endif

      do dir=1,nDOF
       angle_array(dir,1)=angle_array(dir,iicrit)
      enddo
      do dir=1,nEQN
       f_array(dir,1)=f_array(dir,iicrit)
      enddo
      err_array(1)=err_array(iicrit)
      delangle_array(1)=zero
      intercept_array(1)=intercept_array(iicrit)
      do dir=1,nEQN
       cen_array(dir,1)=cen_array(dir,iicrit)
      enddo 
  
      err_local_min=err_array(1)
      err=err_array(1)
      best_error=err*refvfrac

      iter=0

      if (use_only_override_data.eq.1) then
       local_MOFITERMAX=0
      else if (use_only_override_data.eq.0) then
       !do nothing
      else
       print *,"use_only_override_data invalid: ",use_only_override_data
       stop
      endif

      delta_theta=Pi/180.0d0  ! 1 degree=pi/180
      delta_theta_max=10.0d0*Pi/180.0d0  ! 10 degrees
      delta_theta_local=delta_theta

      do while ((iter.lt.local_MOFITERMAX).and. &
                (err.gt.tol).and. &
                (err_local_min.gt.local_tol))

       do dir=1,nEQN
        fbase(dir)=f_array(dir,iter+1)
       enddo

       do i_angle=1,nDOF
        angle_base(i_angle)=angle_array(i_angle,iter+1)
       enddo

       do i_angle=1,nDOF

        do j_angle=1,nDOF
         angle_plus(j_angle)=angle_base(j_angle)
         angle_minus(j_angle)=angle_base(j_angle)
        enddo

        call advance_angle(angle_plus(i_angle),delta_theta_local)
        call advance_angle(angle_minus(i_angle),-delta_theta_local)

        do j_angle=1,nDOF
         angle_plus_archive(i_angle,j_angle)=angle_plus(j_angle)
         angle_minus_archive(i_angle,j_angle)=angle_minus(j_angle)
        enddo

        intp(i_angle)=intercept_array(iter+1)
        intm(i_angle)=intercept_array(iter+1)

        use_initial_guess=1

        local_int=intp(i_angle)

         ! fp=xref-cenp
         ! calling from:
         !   find_cut_geom_slope
        call multi_rotatefunc( &
         tid_in, &
         uncaptured_volume_vof, &
         mofdata, &
         nDOF, & ! sdim-1
         nEQN, & ! sdim 
         use_MilcentLemoine, &
         bfact,dx,xsten0,nhalf0, &
         xtetlist_vof, & !intent(in)
         nlist_vof, & !intent(in)
         nlist_alloc, & !intent(in)
         nmax, &
         refcentroid, &
         refvfrac, &
         angle_plus, & !intent(in)
         fp, & !intent(out)
         local_int, & !intent(inout)
         cenp, & !intent(out)
         use_initial_guess, &
         fastflag,sdim)

        intp(i_angle)=local_int

        err_plus(i_angle)=zero
        do dir=1,nEQN
         err_plus(i_angle)=err_plus(i_angle)+fp(dir)**2
         f_plus(dir,i_angle)=fp(dir)
         cen_plus(dir,i_angle)=cenp(dir)
        enddo
        err_plus(i_angle)=sqrt(err_plus(i_angle))

        local_int=intm(i_angle)

         ! fm=xref-cenm
         ! calling from:
         !   find_cut_geom_slope
        call multi_rotatefunc( &
         tid_in, &
         uncaptured_volume_vof, &
         mofdata, &
         nDOF, & ! sdim-
         nEQN, & ! sdim or 3 * sdim
         use_MilcentLemoine, &
         bfact,dx,xsten0,nhalf0, &
         xtetlist_vof, & !intent(in)
         nlist_vof, & !intent(in)
         nlist_alloc, & !intent(in)
         nmax, &
         refcentroid,refvfrac, &
         angle_minus, & !intent(in)
         fm, & !intent(out)
         local_int, & !intent(inout)
         cenm, & !intent(out)
         use_initial_guess, &
         fastflag,sdim)

        intm(i_angle)=local_int

        err_minus(i_angle)=zero
        do dir=1,nEQN
         err_minus(i_angle)=err_minus(i_angle)+fm(dir)**2
         f_minus(dir,i_angle)=fm(dir)
         cen_minus(dir,i_angle)=cenm(dir)
        enddo
        err_minus(i_angle)=sqrt(err_minus(i_angle))
                
! jacobian matrix has:
!   f1_1  f1_2
!   f2_1  f2_2
!   f3_1  f3_2 

        do dir=1,nEQN
         fgrad(dir,i_angle)=(fp(dir)-fm(dir))/(two*delta_theta_local)
        enddo !do dir=1,nEQN

       enddo  ! i_angle=1..sdim-1

       do i_angle=1,nDOF
        do j_angle=1,nDOF
         JTJ(i_angle,j_angle)=zero
         do dir=1,nEQN
          JTJ(i_angle,j_angle)=JTJ(i_angle,j_angle)+ &
           fgrad(dir,i_angle)*fgrad(dir,j_angle)
         enddo
        enddo !j_angle
       enddo !i_angle

! jacobian matrix has:
!   f1_1  f1_2
!   f2_1  f2_2
!   f3_1  f3_2 
! cost=f dot f
! grad cost = 2 sum_{ij} fi,j * fi   i=1,sdim j=1,sdim-1
! -JT * f = -grad cost/2
       do i_angle=1,nDOF  
        RHS(i_angle)=zero
        do dir=1,nEQN
         RHS(i_angle)=RHS(i_angle)-fgrad(dir,i_angle)*fbase(dir)
        enddo
       enddo !i_angle=1,nDOF

       err_local_min=zero
       do i_angle=1,nDOF  
        err_local_min=err_local_min+RHS(i_angle)**2
       enddo
       err_local_min=sqrt(err_local_min)

       if (sdim.eq.3) then
        DET=JTJ(1,1)*JTJ(2,2)-JTJ(1,2)*JTJ(2,1)
       else if (sdim.eq.2) then
        DET=JTJ(1,1)
       else
        print *,"sdim invalid"
        stop
       endif 

        ! DET has dimensions of length squared
       if (abs(DET).ge.EPS_13_5*(dx_normalize**2)) then 
        !do nothing
       else if (abs(DET).le.EPS_13_5*(dx_normalize**2)) then
        err_local_min=zero
       else
        print *,"DET invalid: ",DET
        stop
       endif

        !err_local_min=||grad COST(angle(s))||
       if (err_local_min.le.local_tol) then

        do dir=1,nEQN
         fopt(dir)=fbase(dir)
         cenopt(dir)=cen_array(dir,iter+1)
        enddo
        do i_angle=1,nDOF
         angle_opt(i_angle)=angle_array(i_angle,iter+1)
        enddo
        intopt=intercept_array(iter+1)

        err=zero 
        do dir=1,nEQN
         err=err+fopt(dir)**2
        enddo
        err=sqrt(err)

        !err_local_min=||grad COST(angles)||
       else if (err_local_min.gt.local_tol) then

        if (sdim.eq.3) then
         JTJINV(1,1)=JTJ(2,2)
         JTJINV(2,2)=JTJ(1,1)
         JTJINV(1,2)=-JTJ(1,2)
         JTJINV(2,1)=-JTJ(2,1)
        else if (sdim.eq.2) then
         JTJINV(1,1)=one
        else
         print *,"sdim invalid"
         stop
        endif

        do i_angle=1,nDOF
         do j_angle=1,nDOF
          JTJINV(i_angle,j_angle)=JTJINV(i_angle,j_angle)/DET
         enddo
        enddo

        do i_angle=1,nDOF  ! compute JTJ^-1 (RHS)
         delangle(i_angle)=zero
         do j_angle=1,nDOF
          delangle(i_angle)=delangle(i_angle)+ &
           JTJINV(i_angle,j_angle)*RHS(j_angle)
         enddo
         if (delangle(i_angle).gt.delta_theta_max) then
          delangle(i_angle)=delta_theta_max
         else if (delangle(i_angle).lt.-delta_theta_max) then
          delangle(i_angle)=-delta_theta_max
         endif
        enddo !i_angle=1,nDOF

         ! -pi<angle<pi 
        do i_angle=1,nDOF
         angle_gauss_newton(i_angle)=angle_base(i_angle)
         call advance_angle(angle_gauss_newton(i_angle),delangle(i_angle))
        enddo

        int_gauss_newton=intercept_array(iter+1)

        use_initial_guess=1

         ! calling from:
         !   find_cut_geom_slope
        call multi_rotatefunc( &
         tid_in, &
         uncaptured_volume_vof, &
         mofdata, &
         nDOF, & ! sdim-1
         nEQN, & ! sdim
         use_MilcentLemoine, &
         bfact,dx,xsten0,nhalf0, &
         xtetlist_vof, & !intent(in)
         nlist_vof, & !intent(in)
         nlist_alloc, & !intent(in)
         nmax, &
         refcentroid, &
         refvfrac, &
         angle_gauss_newton, & !intent(in)
         f_gauss_newton, & !intent(out)
         int_gauss_newton, & !intent(inout)
         cen_gauss_newton, & !real(amrex_real), INTENT(out) :: testcen(nEQN)
         use_initial_guess, &
         fastflag,sdim)

        err_gauss_newton=zero
        do dir=1,nEQN
         err_gauss_newton=err_gauss_newton+f_gauss_newton(dir)**2
        enddo
        err_gauss_newton=sqrt(err_gauss_newton)

        do i_angle=1,nDOF  !  delta_theta_local*(RHS)/||RHS||
         delangle(i_angle)=delta_theta_local*RHS(i_angle)/err_local_min
        enddo !i_angle=1,nDOF

         ! -pi<angle<pi 
        do i_angle=1,nDOF
         angle_steepest(i_angle)=angle_base(i_angle)
         call advance_angle(angle_steepest(i_angle),delangle(i_angle))
        enddo

        int_steepest=intercept_array(iter+1)

        use_initial_guess=1

         ! calling from:
         !   find_cut_geom_slope
        call multi_rotatefunc( &
         tid_in, &
         uncaptured_volume_vof, &
         mofdata, &
         nDOF, & ! sdim-1
         nEQN, & ! sdim
         use_MilcentLemoine, &
         bfact,dx,xsten0,nhalf0, &
         xtetlist_vof, & !intent(in)
         nlist_vof, & !intent(in)
         nlist_alloc, & !intent(in)
         nmax, &
         refcentroid, &
         refvfrac, &
         angle_steepest, & !intent(in)
         f_steepest, & !intent(out)
         int_steepest, &  !intent(inout)
         cen_steepest, & !real(amrex_real), INTENT(out) :: testcen(nEQN)
         use_initial_guess, &
         fastflag,sdim)

        err_steepest=zero
        do dir=1,nEQN
         err_steepest=err_steepest+f_steepest(dir)**2
        enddo
        err_steepest=sqrt(err_steepest)

        do dir=1,nEQN
         fopt(dir)=fbase(dir)
         cenopt(dir)=cen_array(dir,iter+1)
        enddo
        do i_angle=1,nDOF
         angle_opt(i_angle)=angle_array(i_angle,iter+1)
        enddo
        intopt=intercept_array(iter+1)
     
        err=zero 
        do dir=1,nEQN
         err=err+fopt(dir)**2
        enddo
        err=sqrt(err)

        do i_angle=1,nDOF

         if (err.le.err_plus(i_angle)) then
          ! do nothing
         else if (err.ge.err_plus(i_angle)) then

          err=err_plus(i_angle)
          do dir=1,nEQN
           fopt(dir)=f_plus(dir,i_angle)
           cenopt(dir)=cen_plus(dir,i_angle)
          enddo 
          intopt=intp(i_angle)
          do j_angle=1,nDOF
           angle_opt(j_angle)=angle_plus_archive(i_angle,j_angle)
          enddo
         else
          print *,"err or err_plus bad:",i_angle,err,err_plus(i_angle)
          stop
         endif

         if (err.le.err_minus(i_angle)) then
          ! do nothing
         else if (err.ge.err_minus(i_angle)) then
          err=err_minus(i_angle)
          do dir=1,nEQN
           fopt(dir)=f_minus(dir,i_angle)
           cenopt(dir)=cen_minus(dir,i_angle)
          enddo 
          intopt=intm(i_angle)
          do j_angle=1,nDOF
           angle_opt(j_angle)=angle_minus_archive(i_angle,j_angle)
          enddo
         else
          print *,"err or err_minus bad:",i_angle,err,err_minus(i_angle)
          stop
         endif

        enddo ! i_angle=1..nDOF
   
        if (err.le.err_gauss_newton) then
         ! do nothing
        else if (err.ge.err_gauss_newton) then

         err=err_gauss_newton
         do dir=1,nEQN
          fopt(dir)=f_gauss_newton(dir)
          cenopt(dir)=cen_gauss_newton(dir)
         enddo 
         intopt=int_gauss_newton
         do j_angle=1,nDOF
          angle_opt(j_angle)=angle_gauss_newton(j_angle)
         enddo

        else
         print *,"err or err_gauss_newton bad:",err,err_gauss_newton
         stop
        endif

        if (err.le.err_steepest) then
         ! do nothing
        else if (err.ge.err_steepest) then

         err=err_steepest
         do dir=1,nEQN
          fopt(dir)=f_steepest(dir)
          cenopt(dir)=cen_steepest(dir)
         enddo 
         intopt=int_steepest
         do j_angle=1,nDOF
          angle_opt(j_angle)=angle_steepest(j_angle)
         enddo

        else
         print *,"err or err_steepest bad:",err,err_steepest
         stop
        endif

        if (MOFDEB.eq.1) then
         print *,"iter,err_plus,err_minus,err_gauss_newton,err_steepest ", &
          iter,err_plus,err_minus,err_gauss_newton,err_steepest
        endif

       else
        print *,"err_local_min invalid ",err_local_min
        stop
       endif

       do dir=1,nEQN
        f_array(dir,iter+2)=fopt(dir)
        cen_array(dir,iter+2)=cenopt(dir)
       enddo

       do i_angle=1,nDOF
        angle_array(i_angle,iter+2)=angle_opt(i_angle) 
       enddo
       delangle_array(iter+2)=delta_theta_local

       err_array(iter+2)=err

       intercept_array(iter+2)=intopt

       if (iter.ge.1) then
        if (err_array(iter+2).lt.err_array(iter+1)) then
         delta_theta_local=min(two*delta_theta_local,delta_theta)
        else
         delta_theta_local=half*delta_theta_local
        endif
       endif
         
       iter=iter+1
      enddo ! while error>tol and iter<local_MOFITERMAX

      iicrit=iter
      do ii=0,iter

       if ((MOF_DEBUG_RECON.eq.1).or. &
           (MOFDEB.eq.1).or. &
           ((MOF_DEBUG_RECON.eq.2).and.(fastflag.eq.0))) then
        print *,"iimof,eemof,delta ",ii,err_array(ii+1), &
           delangle_array(ii+1)
       endif

       if (err_array(ii+1).le.err_array(iicrit+1)) then
        iicrit=ii
       endif
       if (MOFDEB.eq.1) then
        print *,"regular iter,angle,err,del ",ii+1, &
          angle_array(1,ii+1),angle_array(2,ii+1), &
          err_array(ii+1),delangle_array(ii+1)
       endif

      enddo ! ii=0..iter

      if ((MOF_DEBUG_RECON.eq.1).or. &
          ((MOF_DEBUG_RECON.eq.2).and.(fastflag.eq.0))) then
       MOF_DEBUG_RECON_COUNT=MOF_DEBUG_RECON_COUNT+1
       print *,"MOF_DEBUG_RECON_COUNT=",MOF_DEBUG_RECON_COUNT
       print *,"AFTER------------------------------- "
      endif

      best_error=err_array(iicrit+1)*refvfrac

      do dir=1,nDOF
       new_angle(dir)=angle_array(dir,iicrit+1)
      enddo

      intercept=intercept_array(iicrit+1)

      call angle_to_slope( &
        new_angle,nslope,sdim)

      do dir=1,nEQN
       centroidA(dir)=cen_array(dir,iicrit+1)
      enddo

      if (use_MilcentLemoine.eq.0) then
       ! do nothing
      else if (use_MilcentLemoine.eq.1) then

       if ((nDOF.eq.sdim-1).and. &
           (nEQN.eq.sdim)) then
        ! do nothing
       else
        print *,"invalid nDOF or nEQN"
        stop
       endif

       use_initial_guess=0

       call multi_find_intercept( &
        tid_in, &
        nDOF, & ! sdim-1
        nEQN, & ! sdim
        bfact,dx, &
        xsten0,nhalf0, &
        nslope, &
        intercept, &
        xtetlist_vof, & !intent(in)
        nlist_alloc, &
        nlist_vof, &
        nmax, &
        refvfrac, &
        use_initial_guess, &
        cen_derive_init, & !intent(out)
        fastflag,sdim)
      else
       print *,"use_MilcentLemoine invalid"
       stop
      endif

      mof_iterations(tid_in+1,critical_material)= &
       mof_iterations(tid_in+1,critical_material)+iter
      mof_errors(tid_in+1,critical_material)= &
       mof_errors(tid_in+1,critical_material)+best_error
      mof_calls(tid_in+1,critical_material)= &
       mof_calls(tid_in+1,critical_material)+1

      return
      end subroutine find_cut_geom_slope

      subroutine find_predict_slope( &
       slope, & !intent(out)
       mag_vec, & !intent(out)
       vof_scale, &
       vof_free, &
       cen_free, &
       vof_ref, &
       cen_ref, &
       bfact,dx, &
       xsten,nhalf,sdim)
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: nhalf
      integer :: dir
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(out) :: slope(sdim)
      real(amrex_real) local_slope(sdim)
      real(amrex_real), INTENT(in) :: vof_scale
      real(amrex_real), INTENT(in) :: vof_free
      real(amrex_real), INTENT(in) :: cen_free(sdim)
      real(amrex_real), INTENT(in) :: vof_ref
      real(amrex_real), INTENT(in) :: cen_ref(sdim)
      real(amrex_real) :: local_ref(sdim)
      real(amrex_real) :: local_free(sdim)
      real(amrex_real), INTENT(out) :: mag_vec

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid135"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid"
       stop
      endif
      if ((levelrz.eq.COORDSYS_RZ).and.(sdim.ne.2)) then
       print *,"dimension bust"
       stop
      endif
      if ((vof_ref.ge.zero).and. &
          (vof_scale.ge.zero).and. &
          (vof_free.ge.zero)) then
       !do nothing
      else
       print *,"vof_ref, vof_free, or vof_scale less than zero"
       print *,"vof_ref: ",vof_ref
       print *,"vof_scale: ",vof_scale
       print *,"vof_free: ",vof_free
       stop
      endif

      do dir=1,sdim
       local_ref(dir)=cen_ref(dir)
       local_free(dir)=cen_free(dir)
      enddo !dir=1,sdim
           
      if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
          (levelrz.eq.COORDSYS_RZ)) then
       do dir=1,sdim
        local_slope(dir)=local_ref(dir)-local_free(dir)
       enddo
       ! mag_vec=|slope|
       call prepare_normal(local_slope,mag_vec,sdim)

       if (mag_vec.gt.zero) then
        !do nothing
       else if (mag_vec.eq.zero) then
        do dir=1,sdim
         local_slope(dir)=zero
        enddo
        local_slope(1)=one
       else
        print *,"mag_vec invalid find_predict_slope: ",mag_vec
        print *,"levelrz= ",levelrz
        stop
       endif

       do dir=1,sdim
        slope(dir)=local_slope(dir)
       enddo

      else
       print *,"find_predict_slope: levelrz invalid:",levelrz
       stop
      endif

      return
      end subroutine find_predict_slope

! this routine is called at most num_materials times; after uncaptured_volume=0,
! routine is not called anymore.
! First time routine is called, the flag for all materials is 0
! reference and actual centroid relative to cell centroid.
! xcell is cell center (not cell centroid)

      subroutine individual_MOF( &
        layer_flag, &
        tid_in, &
        ls_mof, &
        lsnormal, &
        lsnormal_valid, &
        override_normal, &
        override_normal_valid, &
        bfact,dx, &
        xsten0,nhalf0, &
        order_algorithm_in, &
        xtetlist_vof, & !intent(out)
        nlist_alloc, & !intent(in)
        nmax, &
        mofdata, &
        loop_counter, &
        uncaptured_volume_fraction, &
        vfrac_sum_local, &
        num_processed_total, &
        num_processed, &
        num_override_history, &
        override_target, &
        recon_history, &
        material_used, &
        uncaptured_volume_vof, &
        multi_centroidA, &
        sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(IN) :: layer_flag
      integer, INTENT(IN) :: nlist_alloc 
      integer, INTENT(IN) :: tid_in
      integer, INTENT(in) :: sdim 
      integer, INTENT (IN) :: nhalf0
      integer, INTENT (IN) :: bfact
      real(amrex_real), INTENT (IN), DIMENSION(sdim) :: dx
      real(amrex_real), INTENT (IN), DIMENSION(-nhalf0:nhalf0,sdim) :: xsten0
      real(amrex_real), DIMENSION(-nhalf0:nhalf0,sdim) :: xsten_local
      integer dir
      integer i_stencil_node
      integer, INTENT(in) :: nmax
      integer, INTENT(in) :: order_algorithm_in(num_materials)
      real(amrex_real), INTENT(inout) :: uncaptured_volume_fraction
      real(amrex_real), INTENT(in) :: vfrac_sum_local
      real(amrex_real), INTENT(inout) :: uncaptured_volume_vof 
      real(amrex_real), INTENT(inout) :: mofdata(num_materials*ngeom_recon)
      real(amrex_real) centroidA(num_materials*sdim)
      real(amrex_real), INTENT(out) :: multi_centroidA(num_materials,sdim)
      integer im
      integer vofcomp
      integer vofcomp_single
      integer, INTENT(inout) :: num_processed_total
      integer, INTENT(inout) :: num_processed
      integer, INTENT(inout) :: material_used(num_materials)
      integer, INTENT(in) :: loop_counter
      integer, INTENT(in) :: num_override_history
      integer, INTENT(in) :: override_target(num_materials)
      integer, INTENT(inout) :: recon_history(num_materials)
      real(amrex_real) distmax
      integer ordermax
      integer rank_min
      integer rank_current
      integer critical_material
      real(amrex_real) volcut_vof,volcell_vof
      real(amrex_real) mag_vec
      real(amrex_real) cencut_vof(sdim)
      real(amrex_real) cencell_vof(sdim)
      real(amrex_real), INTENT(out) :: xtetlist_vof(4,3,nlist_alloc)
      integer nlist_vof
      integer single_material_takes_all
      integer single_material_im
      integer test_order
      real(amrex_real) test_vfrac,max_vfrac
      integer use_initial_guess
      real(amrex_real) npredict(sdim)
      real(amrex_real) refcentroid(num_materials*sdim)
      real(amrex_real) centroid_ref(sdim)
      real(amrex_real) centroid_free(sdim)
      real(amrex_real) refvfrac
      real(amrex_real) single_volume
      real(amrex_real) nslope(num_materials*sdim)
      real(amrex_real) intercept
      integer mat_before,vofcomp_before
      integer fastflag
      real(amrex_real) vofmain(num_materials)

      real(amrex_real), INTENT(in) :: &
        ls_mof(D_DECL(-1:1,-1:1,-1:1),num_materials)

      real(amrex_real), INTENT(in) :: lsnormal(num_materials,sdim)
      integer, INTENT(in) :: lsnormal_valid(num_materials)

      real(amrex_real), INTENT(in) :: override_normal(num_materials,sdim)
      integer, INTENT(in) :: override_normal_valid(num_materials)

      integer :: nDOF_standard
      integer :: nEQN_standard
      integer is_masked(num_materials)

#include "mofdata.H"

      nDOF_standard=sdim-1
      nEQN_standard=sdim

      fastflag=1
      nlist_vof=0

      call init_mask_flag(layer_flag,is_masked)

      if ((tid_in.lt.0).or.(tid_in.ge.geom_nthreads)) then
       print *,"tid_in invalid individual_MOF: ",tid_in
       stop
      endif

      if (nhalf0.lt.3) then
       print *,"expecting nhalf0>=3: ",nhalf0
       stop
      endif
      if ((uncaptured_volume_fraction.ge.one-vfrac_sum_local).and. &
          (uncaptured_volume_fraction.le.one+VOFTOL_MATERIAL)) then
       ! do nothing
      else
       print *,"uncaptured_volume_fraction invalid: ", &
         uncaptured_volume_fraction
       print *,"vfrac_sum_local= ", &
         vfrac_sum_local
       print *,"VOFTOL_MATERIAL= ",VOFTOL_MATERIAL
       stop
      endif
      if ((vfrac_sum_local.ge.zero).and. &
          (vfrac_sum_local.le.one+VOFTOL_MATERIAL)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ", &
         vfrac_sum_local
       print *,"VOFTOL_MATERIAL= ",VOFTOL_MATERIAL
       stop
      endif

      if (uncaptured_volume_vof.gt.zero) then
       ! do nothing
      else
       print *,"uncaptured_volume_vof invalid: ",uncaptured_volume_vof
       stop
      endif
      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon invalid(individual_MOF) ",ngeom_recon
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid individual_MOF ",sdim
       stop
      endif
      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid individual_MOF ",num_materials
       stop
      endif
      if ((loop_counter.lt.0).or. &
          (loop_counter.ge.num_materials)) then
       print *,"loop_counter invalid ",loop_counter
       stop
      endif
      if (nmax.lt.10) then
       print *,"nmax too small ",nmax
       stop
      endif

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid: ",nlist_alloc
       stop
      endif

      do i_stencil_node=-nhalf0,nhalf0
      do dir=1,sdim
       xsten_local(i_stencil_node,dir)=xsten0(i_stencil_node,dir)
      enddo !dir
      enddo !i_stencil_node


       ! cencell_vof is in absolute coordinate system

      call Box_volumeFAST( &
        bfact,dx,xsten0,nhalf0, &
        volcell_vof, &
        cencell_vof, &
        sdim)

      if ((num_processed_total.ge.1).and. &
          (num_processed_total.lt.num_materials)) then
       !do nothing
      else if (num_processed_total.eq.0) then
       !do nothing
      else 
       print *,"num_processed_total invalid: ",num_processed_total
       stop
      endif

      if ((num_processed.ge.1).and. &
          (num_processed.lt.num_materials)) then
       fastflag=0
      else if (num_processed.eq.0) then
       fastflag=1
      else 
       print *,"num_processed invalid: ",num_processed
       stop
      endif

       ! cencut is in absolute coordinate system

      if (fastflag.eq.0) then

         ! get triangulation of uncaptured space in the cell
         ! in: individual_MOF
       call tets_box_planes( &
         layer_flag, &
         bfact,dx, &
         xsten0,nhalf0, &
         xsten0,nhalf0, &
         mofdata, &
         xtetlist_vof, &
         nlist_alloc, &
         nlist_vof, &
         nmax, &
         sdim)

       call get_cut_geom3D(xtetlist_vof, &
         nlist_alloc,nlist_vof,nmax, &
         volcut_vof,cencut_vof,sdim)

      else if (fastflag.eq.1) then
       volcut_vof=volcell_vof
       do dir=1,sdim
        cencut_vof(dir)=cencell_vof(dir)  
       enddo
      else
       print *,"fastflag invalid individual MOF ",fastflag
       stop
      endif      

      if (abs(volcut_vof-uncaptured_volume_vof).le. &
          EPS_8_2*volcell_vof) then
        !do nothing
      else
        print *,"(breakpoint) break point and gdb: "
        print *,"(1) compile with the -g option"
        print *,"(2) break MOF.F90:13110"
        print *,"volcut_vof invalid individual mof"
        print *,"num_processed ",num_processed
        print *,"order_algorithm_in ",order_algorithm_in
        print *,"bfact ",bfact
        print *,"dx ",dx

        do dir=1,sdim
        do im=-nhalf0,nhalf0
         print *,"i,dir,xsten0(i,dir) ",im,dir,xsten0(im,dir)
        enddo
        enddo

        print *,"xsten0 ",xsten0
        print *,"nhalf0 ",nhalf0
        print *,"fastflag ",fastflag
        print *,"nlist_alloc ",nlist_alloc
        print *,"nlist_vof ",nlist_vof
        print *,"volcell_vof ",volcell_vof
        print *,"cencell_vof ",cencell_vof
        print *,"volcut_vof ",volcut_vof
        print *,"cencut_vof ",cencut_vof
        print *,"uncaptured_volume_fraction ",uncaptured_volume_fraction
        print *,"vfrac_sum_local ",vfrac_sum_local
        print *,"uncaptured_volume_vof ",uncaptured_volume_vof
        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         print *,"im,vof,centroid ",im,mofdata(vofcomp),mofdata(vofcomp+1), &
           mofdata(vofcomp+2),mofdata(vofcomp+sdim)
         print *,"im,flag ",mofdata(vofcomp+sdim+1)
         print *,"im,slope ",im,mofdata(vofcomp+sdim+2), &
            mofdata(vofcomp+sdim+3), &
            mofdata(vofcomp+sdim+sdim+1)
         print *,"im,intercept ",im,mofdata(vofcomp+2*sdim+2)
        enddo
        stop
      endif

         ! figure out the next material to fill the unoccupied region.
         ! the order of the next material will be "ordermax+1"
      distmax=-one
      rank_min=9999
      ordermax=0

      do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        if (is_masked(im).eq.1) then
         ! do nothing
        else if (is_masked(im).eq.0) then
         if (NINT(mofdata(vofcomp+sdim+1)).gt.ordermax) then
          ordermax=NINT(mofdata(vofcomp+sdim+1))
         endif
        else
         print *,"is_masked invalid: ",im,is_masked(im)
         stop
        endif
        vofmain(im)=mofdata(vofcomp)
      enddo ! im=1..num_materials

      if (ordermax.ge.num_materials) then
        print *,"all the materials already initialized"
        print *,"ordermax= ",ordermax
        print *,"num_materials= ",num_materials
        stop
      endif
      if (ordermax.lt.0) then
        print *,"ordermax invalid: ",ordermax
        print *,"num_materials= ",num_materials
        stop
      endif

      critical_material=0
      single_material_takes_all=0
      single_material_im=0
      do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        if (is_masked(im).eq.1) then
         ! do nothing
        else if (is_masked(im).eq.0) then
         if (NINT(mofdata(vofcomp+sdim+1)).eq.0) then
          single_volume=mofdata(vofcomp)*volcell_vof
          if (single_volume.ge.(one-EPS_UNCAPTURED)*uncaptured_volume_vof) then

           if (single_material_im.eq.0) then
            single_material_im=im
            single_material_takes_all=1
            critical_material=im
            distmax=-one  ! tells code below not to search for a slope
           else if ((single_material_im.ge.1).and. &
                    (single_material_im.le.num_materials)) then
            vofcomp_single=(single_material_im-1)*ngeom_recon+1
            if (mofdata(vofcomp_single).lt. &
                mofdata(vofcomp)) then
             single_material_im=im
             critical_material=im
            else if (mofdata(vofcomp_single).ge. &
                     mofdata(vofcomp)) then
             !do nothing
            else
             print *,"mofdata invalid ",mofdata
             stop
            endif
           else
            print *,"single_material_im invalid ",single_material_im
            stop
           endif

          endif
         endif  ! material not already processed.
        else
         print *,"is_masked invalid MOF.F90: ",im,is_masked(im)
         stop
        endif
      enddo ! do im=1,num_materials (init single_material,crit_mat,distmax)

      if (single_material_takes_all.eq.1) then
        uncaptured_volume_vof=zero
        uncaptured_volume_fraction=zero
      else if (single_material_takes_all.eq.0) then
        ! do nothing
      else 
        print *,"bust individual_MOF (1)"
        print *,"sdim,num_materials ",sdim,num_materials
        print *,"single_material_takes_all ",single_material_takes_all
        stop
      endif
        
          ! if uncaptured_volume_vof=0, then there is no need to find the
          ! slope since "single_material_takes_all=1". 
      if ((uncaptured_volume_vof.gt.zero).and. &
          (uncaptured_volume_fraction.gt.zero)) then

         ! find unprocessed material whose moment is furthest from cencut. 

        do im=1,num_materials

         if (is_masked(im).eq.1) then
          ! do nothing
         else if (is_masked(im).eq.0) then

          vofcomp=(im-1)*ngeom_recon+1

          single_volume=mofdata(vofcomp)*volcell_vof

          ! only find the slope if refvfrac<>0 and refvfrac<available vol.
          if ((NINT(mofdata(vofcomp+sdim+1)).eq.0).and. & ! order==0?
              (single_volume.ge.EPS_UNCAPTURED*volcell_vof).and. & ! 0<V<Vavail?
              (single_volume.le. &
               (one-EPS_UNCAPTURED)*uncaptured_volume_vof)) then

           do dir=1,sdim
            centroid_free(dir)=cencut_vof(dir)
            centroid_ref(dir)=mofdata(vofcomp+dir)+cencell_vof(dir)
           enddo

           ! centroid_ref-centroid_free
           ! normal points from light to dark
           call find_predict_slope( &
             npredict, & !intent(out)
             mag_vec, & !intent(out)
             volcell_vof, &
             volcut_vof, &
             centroid_free, & !centroid of uncaptured region
             single_volume, &
             centroid_ref, &  !absolute coordinate system
             bfact,dx, &
             xsten_local,nhalf0,sdim)
          
           if (mag_vec.gt.zero) then

            rank_current=get_material_rank(im,order_algorithm_in(im))

             ! rank_min initialized to be 9999.
            if ((order_algorithm_in(im).lt.1).or. &
                (order_algorithm_in(im).gt.num_materials+1).or. &
                (rank_current.lt.1)) then
             print *,"order_algorithm_in invalid: ", &
               im,order_algorithm_in(im)
             stop
            else if ((order_algorithm_in(im).ge.1).and. &
                     (order_algorithm_in(im).le.num_materials+1).and. &
                     (rank_current.ge.1).and. &
                     (num_override_history.gt.0).and. &
                     (num_processed_total.lt.num_override_history).and. &
                     (override_target(num_processed_total+1).eq.im)) then
             distmax=1.0e+10
             rank_min=0
             critical_material=im
            else if (rank_current.lt.rank_min) then
             distmax=mag_vec
             rank_min=rank_current
             critical_material=im
            else if (rank_current.eq.rank_min) then
             if (mag_vec.gt.distmax) then
              distmax=mag_vec
              critical_material=im
             endif
            else if (rank_current.gt.rank_min) then
             !do nothing
            else
             print *,"rank_current invalid ",rank_current
             print *,"or rank_min invalid ",rank_min
             stop
            endif

           else if (mag_vec.eq.zero) then
            ! do nothing
           else
            print *,"mag_vec invalid: ",mag_vec
            stop
           endif

          else if ((NINT(mofdata(vofcomp+sdim+1)).ge.1).or. & ! order>=1?
                   (abs(single_volume).le.EPS_UNCAPTURED*volcell_vof).or. &
                   (single_volume.ge. &
                    (one-EPS_UNCAPTURED)*uncaptured_volume_vof)) then
           ! do nothing
          else
           print *,"order or single_volume invalid"
           print *,"single_volume=",single_volume
           print *,"order=",NINT(mofdata(vofcomp+sdim+1))
           stop
          endif 
         else
          print *,"is_masked invalid MOF.F90: ",is_masked
          print *,"im=",im
          stop
         endif

        enddo ! im=1..num_materials

      else if ((uncaptured_volume_vof.eq.zero).or. &
               (uncaptured_volume_fraction.eq.zero)) then

        if (distmax.eq.-one) then
         ! do nothing
        else
         print *,"distmax should be -one here: ",distmax
         stop
        endif

      else
        print *,"uncaptured_volume_vof invalid: ",uncaptured_volume_vof
        print *,"or uncaptured_volume_fraction invalid: ", &
           uncaptured_volume_fraction
        stop
      endif 

        ! find MOF slope 
      if (distmax.gt.EPS_8_4*dx(1)) then

         if ((critical_material.lt.1).or. &
             (critical_material.gt.num_materials)) then
          print *,"bust individual_MOF (2)"
          print *,"critical_material=",critical_material
          print *,"num_materials=",num_materials
          stop
         endif

         if (is_masked(critical_material).ne.0) then
          print *,"is_masked invalid MOF.F90: ",is_masked
          print *,"critical_material=",critical_material
          stop
         endif

         vofcomp=(critical_material-1)*ngeom_recon+1

         do dir=1,sdim
          centroid_free(dir)=cencut_vof(dir)
          centroid_ref(dir)=mofdata(vofcomp+dir)+cencell_vof(dir)
         enddo

         single_volume=mofdata(vofcomp)*volcell_vof

          ! centroid_ref-centroid_free
          ! normal points from light to dark
         call find_predict_slope( &
           npredict, & !intent(out)
           mag_vec, & !intent(out)
           volcell_vof, &
           volcut_vof, &
           centroid_free, & !centroid of uncaptured region
           single_volume, &
           centroid_ref, &  !absolute coordinate system
           bfact,dx, &
           xsten_local,nhalf0,sdim)

         if (mag_vec.gt.zero) then
          !do nothing
         else
          print *,"mag_vec underflow: ",mag_vec
          stop
         endif

         do dir=1,sdim
          refcentroid(dir)=mofdata(vofcomp+dir)
         enddo
         refvfrac=mofdata(vofcomp)

           ! centroidA and refcentroid relative to cell centroid of the cell.
           ! find_cut_geom_slope called from: individual_MOF
         call find_cut_geom_slope( &
           tid_in, &
           uncaptured_volume_vof, &
           mofdata, &
           ls_mof, &
           lsnormal, &
           lsnormal_valid, &
           override_normal, &
           override_normal_valid, &
           bfact,dx, &
           xsten0,nhalf0, &
           refcentroid, &
           refvfrac, &
           npredict, & !intent(in)
           nslope, &
           intercept, &
           xtetlist_vof, & !intent(in)
           nlist_vof, & !intent(in)
           nlist_alloc, & !intent(in)
           centroidA, &
           nmax, &
           critical_material, & !INTENT(in)
           fastflag, &
           sdim, &
           nDOF_standard, & !sdim-1
           nEQN_standard)   !sdim

         material_used(critical_material)=ordermax+1
         num_processed=num_processed+1
         num_processed_total=num_processed_total+1
         recon_history(num_processed_total)=critical_material

         mofdata(vofcomp+sdim+1)=ordermax+1
         mofdata(vofcomp+2*sdim+2)=intercept
         do dir=1,sdim
          mofdata(vofcomp+sdim+1+dir)=nslope(dir)
          multi_centroidA(critical_material,dir)=centroidA(dir)
         enddo 
         uncaptured_volume_vof=uncaptured_volume_vof-refvfrac*volcell_vof
         if (uncaptured_volume_vof.le.volcell_vof*EPS_UNCAPTURED) then
          uncaptured_volume_vof=zero
         endif
         uncaptured_volume_fraction=uncaptured_volume_fraction-refvfrac
         if (uncaptured_volume_fraction.le. &
             one-vfrac_sum_local+EPS_UNCAPTURED) then
          uncaptured_volume_fraction=one-vfrac_sum_local
         endif

         ! above MOF reconstruct, below default slopes.
      else if ((distmax.eq.-one).or. &
               ((distmax.ge.zero).and. &
                (distmax.le.EPS_8_4*dx(1)))) then

          ! if single_material_takes_all=1, then
          !  distmax<0
          !  critical_material>0

         if (single_material_takes_all.eq.1) then
          critical_material=single_material_im
         else if (single_material_takes_all.eq.0) then

          max_vfrac=zero
          do im=1,num_materials 
           vofcomp=(im-1)*ngeom_recon+1
           test_vfrac=mofdata(vofcomp)
           test_order=NINT(mofdata(vofcomp+sdim+1))
           if (is_masked(im).eq.1) then
            ! do nothing
           else if (is_masked(im).eq.0) then
            if ((test_order.eq.0).and. &
                (test_vfrac.ge.max_vfrac)) then
             max_vfrac=test_vfrac
             critical_material=im
            endif
           else
            print *,"is_masked invalid: ",im,is_masked
            stop
           endif
          enddo ! im=1..num_materials
        
         else 
          print *,"single_material_takes_all invalid: ", &
            single_material_takes_all
          stop
         endif

         if ((critical_material.lt.1).or. &
             (critical_material.gt.num_materials)) then
          print *,"bust individual_MOF (3)"
          print *,"sdim,num_materials,critical_material ", &
                  sdim,num_materials,critical_material
          print *,"ngeom_recon= ",ngeom_recon
          print *,"vfrac_sum_local ",vfrac_sum_local
          print *,"uncaptured_volume_vof ",uncaptured_volume_vof
          do im=1,num_materials
           vofcomp=(im-1)*ngeom_recon+1
           print *,"im,vf ",im,mofdata(vofcomp)
           print *,"im,status ",im,NINT(mofdata(vofcomp+sdim+1))
          enddo
          stop
         endif

         if (is_masked(critical_material).ne.0) then
          print *,"is_masked invalid MOF.F90: ",is_masked
          print *,"critical_material=",critical_material
          stop
         endif
          
         vofcomp=(critical_material-1)*ngeom_recon+1

         test_order=NINT(mofdata(vofcomp+sdim+1))
         if (test_order.ne.0) then
          print *,"test_order invalid: ",critical_material,test_order
          stop
         endif

         mat_before=0
         if (ordermax.gt.0) then
          do im=1,num_materials
           vofcomp_before=(im-1)*ngeom_recon+1
           if (is_masked(im).eq.1) then
            ! do nothing
           else if (is_masked(im).eq.0) then
            if (NINT(mofdata(vofcomp_before+sdim+1)).eq.ordermax) then
             mat_before=im
            endif
           else
            print *,"is_masked invalid MOF.F90: ",is_masked
            stop
           endif
          enddo ! im=1..num_materials
         else if (ordermax.eq.0) then
          ! do nothing
         else
          print *,"ordermax invalid: ",ordermax
          stop
         endif

         if ((mat_before.ge.1).and. &
             (mat_before.le.num_materials)) then

          if (is_masked(mat_before).ne.0) then
           print *,"is_masked invalid MOF.F90: ",is_masked
           print *,"mat_before=",mat_before
           stop
          endif

          vofcomp_before=(mat_before-1)*ngeom_recon+1
          do dir=1,sdim
           npredict(dir)=-mofdata(vofcomp_before+sdim+1+dir)
          enddo
          intercept=-mofdata(vofcomp_before+2*sdim+2)

         else if (mat_before.eq.0) then

          do dir=1,sdim
           centroid_free(dir)=cencell_vof(dir)
           centroid_ref(dir)=cencut_vof(dir)
          enddo

          ! centroid_ref-centroid_free
          call find_predict_slope( &
            npredict, & !intent(out)
            mag_vec, & !intent(out)
            volcell_vof, &
            volcell_vof, &
            centroid_free, &
            volcut_vof, &
            centroid_ref, &
            bfact,dx, &
            xsten_local,nhalf0,sdim)

          if (mag_vec.gt.zero) then
           ! do nothing
          else if (mag_vec.eq.zero) then
           do dir=1,sdim
            npredict(dir)=zero
           enddo
           npredict(1)=one
          else
           print *,"mag_vec invalid: ",mag_vec
           stop
          endif
          intercept=mofdata(vofcomp)-half
         else
          print *,"mat_before invalid: ",mat_before
          stop
         endif
 
         if ((single_material_takes_all.eq.1).and. &
             (mat_before.ge.1).and. &
             (mat_before.le.num_materials)) then

          material_used(critical_material)=ordermax+1
          num_processed=num_processed+1
          num_processed_total=num_processed_total+1
          recon_history(num_processed_total)=critical_material

          mofdata(vofcomp+sdim+1)=ordermax+1

          mofdata(vofcomp+2*sdim+2)=intercept
          do dir=1,sdim
           mofdata(vofcomp+sdim+1+dir)=npredict(dir)
            ! cencut_vof is the uncaptured centroid in absolute frame 
           centroidA(dir)=cencut_vof(dir)
           multi_centroidA(critical_material,dir)= &
            centroidA(dir)-cencell_vof(dir)
          enddo 

         else if ((single_material_takes_all.eq.0).or. &
                  (mat_before.eq.0)) then

          refvfrac=mofdata(vofcomp)
          use_initial_guess=0

           ! inside of individual MOF
           ! centroidA in absolute coordinate system.

          call multi_find_intercept( &
           tid_in, &
           nDOF_standard, & !sdim-1
           nEQN_standard, & !sdim
           bfact,dx, &
           xsten0,nhalf0, &
           npredict, & !intent(in)
           intercept, &
           xtetlist_vof, & !intent(in)
           nlist_alloc, &  !intent(in)
           nlist_vof, &    !intent(in)
           nmax, &         !intent(in)
           refvfrac, &
           use_initial_guess,centroidA,fastflag,sdim)

          uncaptured_volume_vof=uncaptured_volume_vof-refvfrac*volcell_vof
          if (uncaptured_volume_vof.le.volcell_vof*EPS_UNCAPTURED) then
           uncaptured_volume_vof=zero
          endif
          uncaptured_volume_fraction=uncaptured_volume_fraction-refvfrac
          if (uncaptured_volume_fraction.le. &
              one-vfrac_sum_local+EPS_UNCAPTURED) then
           uncaptured_volume_fraction=one-vfrac_sum_local
          endif
  
          if (single_material_takes_all.eq.1) then
            ! centroidA is the uncaptured centroid (absolute frame)
           do dir=1,sdim
            centroidA(dir)=cencut_vof(dir)
           enddo
          else if (single_material_takes_all.eq.0) then
           ! do nothing
          else
           print *,"single material takes all invalid: ", &
            single_material_takes_all
           stop
          endif

          material_used(critical_material)=ordermax+1
          num_processed=num_processed+1
          num_processed_total=num_processed_total+1
          recon_history(num_processed_total)=critical_material
  
          mofdata(vofcomp+sdim+1)=ordermax+1
          mofdata(vofcomp+2*sdim+2)=intercept
           ! cencell is the cell centroid
          do dir=1,sdim
           mofdata(vofcomp+sdim+1+dir)=npredict(dir)
           multi_centroidA(critical_material,dir)= &
            centroidA(dir)-cencell_vof(dir)
          enddo 
         else
          print *,"single_material_takes_all or mat_before invalid"
          print *,"single_material_takes_all=",single_material_takes_all
          print *,"mat_before=",mat_before
          stop
         endif

         if (is_masked(critical_material).ne.0) then
          print *,"is_masked invalid MOF.F90: ",is_masked
          stop
         endif

      else
        print *,"distmax invalid: ",distmax
        stop
      endif  

      return
      end subroutine individual_MOF

        ! input: n1d
        ! n1d=1 => im material on top
        ! n1d=-1 => im material on bottom
      subroutine get_col_ht_LS( &
       vof_height_function, &
       crossing_status, &
       bfact, &
       dx, &
       xsten0, &
       dx_col, &
       x_col, &
       x_col_avg, &
       lsdata, &
       vofdata, &
       ht_from_LS, &
       ht_from_VOF, &
       dircrit, & ! dircrit=1..sdim
       n1d, &
       sdim)
      use probcommon_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, parameter :: nhalf_height=9
      integer, parameter :: ngrow_tower=4

      integer, INTENT(in) :: vof_height_function
      integer, INTENT(out) :: crossing_status
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten0( &
        -nhalf_height:nhalf_height, &
        sdim)
      real(amrex_real), INTENT(in) :: dx_col(sdim)
      real(amrex_real), INTENT(in) :: x_col(sdim)
      real(amrex_real), INTENT(in) :: x_col_avg(sdim)
      real(amrex_real), INTENT(in) :: lsdata( &
        -ngrow_tower:ngrow_tower)
      real(amrex_real), INTENT(in) :: vofdata( &
        -ngrow_tower:ngrow_tower)
      real(amrex_real), INTENT(out) :: ht_from_LS
      real(amrex_real), INTENT(out) :: ht_from_VOF
      real(amrex_real), INTENT(in) :: n1d
      integer, INTENT(in) :: dircrit

      integer l
      integer l_vof
      integer lmin,lmax
      integer lvof_min,lvof_max
      integer lcrit
      real(amrex_real) xbottom_stencil
      real(amrex_real) xtop_stencil
      real(amrex_real) xbottom_adjusted
      integer vof_ratio_ht_power
      real(amrex_real) X_AT_ABS_LSMIN,ABS_LSMIN,LSTEST

      real(amrex_real) ls1,ls2,x1,x2,slope
      real(amrex_real) charfn(-ngrow_tower:ngrow_tower)
      real(amrex_real) LS
      real(amrex_real) vof_top_sum,vof_bot_sum
      real(amrex_real) dr,dz
      real(amrex_real) volcell
      real(amrex_real) vof_crit
      real(amrex_real) dx_norm,dx_tan1,dx_tan2

      if (ngrow_tower.ne.4) then
       print *,"expecting ngrow_tower=4 ",ngrow_tower
       stop
      endif
      if (nhalf_height.ne.9) then
       print *,"expecting nhalf_height=9 ",nhalf_height
       stop
      endif

      if (ngrow_distance.lt.4) then
       print *,"ngrow_distance invalid: ",ngrow_distance
       stop
      endif
      if (ngrow_make_distance.ne.ngrow_distance-1) then
       print *,"ngrow_make_distance invalid: ",ngrow_make_distance
       stop
      endif
      if ((num_materials.lt.1).or.(num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid get col ht LS"
       print *,"num_materials= ",num_materials
       stop
      endif

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       if (abs(x_col(1)).ge.zero) then
        !do nothing
       else
        print *,"x_col(1) undef ",x_col(1)
        stop
       endif
       if (abs(x_col_avg(1)).ge.zero) then
        !do nothing
       else
        print *,"x_col_avg(1) undef ",x_col_avg(1)
        stop
       endif
      else if (levelrz.eq.COORDSYS_RZ) then
       if (sdim.ne.2) then
        print *,"sdim invalid"
        stop
       endif
       if (x_col(1).gt.zero) then
        ! do nothing
       else
        print *,"cannot have r<0"
        stop
       endif
       if (x_col_avg(1).gt.zero) then
        ! do nothing
       else
        print *,"cannot have r<0"
        stop
       endif
      else 
       print *,"levelrz invalid get col ht ls"
       stop
      endif

      if ((dircrit.lt.1).or.(dircrit.gt.sdim)) then
       print *,"dircrit invalid get_col_ht_LS dircrit=",dircrit
       stop
      endif
      if ((n1d.ne.one).and.(n1d.ne.-one)) then
       print *,"n1d invalid"
       stop
      endif

      lmin=-ngrow_tower
      lmax=ngrow_tower
      if (levelrz.eq.COORDSYS_RZ) then
       if (dircrit.eq.1) then ! horizontal column
        do while (xsten0(2*lmin,dircrit).lt.zero)
         lmin=lmin+1
         if ((2*lmin.gt.2*ngrow_tower+1).or. &
             (lmin.gt.lmax)) then
          print *,"lmin too big"
          stop
         endif
        enddo ! while (xsten0(2*lmin,dircrit).lt.zero)
        if (xsten0(2*lmin,dircrit).gt.zero) then
         ! do nothing
        else
         print *,"xsten0(2*lmin,dircrit).gt.zero not true"
         stop
        endif
       endif
      else if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else
       print *,"levelrz invalid get col ht ls 2"
       stop
      endif 

      ABS_LSMIN=abs(lsdata(lmin))
      X_AT_ABS_LSMIN=xsten0(2*lmin,dircrit)
      do l=lmin+1,lmax
       LSTEST=abs(lsdata(l))
       if (LSTEST.lt.ABS_LSMIN) then
        ABS_LSMIN=LSTEST
        X_AT_ABS_LSMIN=xsten0(2*l,dircrit)
       else if (LSTEST.ge.ABS_LSMIN) then
        ! do nothing
       else
        print *,"LSTEST or ABS_LSMIN is NaN"
        stop
       endif
      enddo

      xbottom_stencil=xsten0(2*lmin-1,dircrit)
      xtop_stencil=xsten0(2*lmax+1,dircrit)

      lcrit=0
      crossing_status=0

      do l=lmin,lmax

       LS=lsdata(l)
       if (LS.ge.zero) then
        charfn(l)=one
       else if (LS.le.zero) then
        charfn(l)=-one
       else
        print *,"LS is NaN: ",l,LS
        stop
       endif
        
      enddo  ! l=lmin,lmax

       ! n1d>0 if LS>0 material on top and LS<0 material on bottom.
       ! n1d<0 otherwise
      do l=0,ngrow_tower-1

        ! first check upper half of stencil for a crossing
       if ((crossing_status.eq.0).and.(l+1.le.lmax)) then 
        ls1=charfn(l)
        ls2=charfn(l+1)
        ! n1d=1 => im material on top
        ! n1d=-1 => im material on bottom
        if ((ls1*ls2.lt.zero).and. &
            ((ls2-ls1)*n1d.gt.zero)) then
         lcrit=l
         crossing_status=1
         ls1=lsdata(l)
         ls2=lsdata(l+1)
         x1=xsten0(2*l,dircrit)
         x2=xsten0(2*l+2,dircrit)
         if (x1.lt.x2) then
          ! LS=LS1+slope(x-x1)  xzero=-LS1/slope+x1
          if (ls1.eq.zero) then
           ht_from_LS=x1
          else if (ls2.eq.zero) then 
           ht_from_LS=x2
          else if ((ls1.ne.zero).and. &
                   (ls2.ne.zero)) then
           slope=(ls1-ls2)/(x1-x2)
           if (slope.ne.zero) then
            ht_from_LS=x1-ls1/slope
           else
            print *,"slope invalid"
            stop
           endif
          else
           print *,"ls1 or ls2 invalid"
           stop
          endif
         else
          print *,"x1 or x2 invalid"
          stop
         endif
        else if ((ls1*ls2.ge.zero).or. &
                 ((ls2-ls1)*n1d.le.zero)) then
         ! do nothing
        else
         print *,"ls1 or ls2 invalid"
         stop
        endif
       else if ((crossing_status.ne.0).or.(l+1.gt.lmax)) then
        ! do nothing
       else
        print *,"crossing_status or l+1 problem"
        stop
       endif


        ! second: check lower half of stencil for a crossing
       if ((crossing_status.eq.0).and. &
           (-(l+1).ge.lmin)) then 
        ls1=charfn(-l)
        ls2=charfn(-(l+1))
        if ((ls1*ls2.lt.zero).and.((ls1-ls2)*n1d.gt.zero)) then
         lcrit=-(l+1)
         crossing_status=1
         ls1=lsdata(-l)
         ls2=lsdata(-(l+1))
         x1=xsten0(-2*l,dircrit)
         x2=xsten0(-(2*l+2),dircrit)
         if (x1.gt.x2) then
          ! LS=LS1+slope(x-x1)  xzero=-LS1/slope+x1
          if (ls1.eq.zero) then
           ht_from_LS=x1
          else if (ls2.eq.zero) then
           ht_from_LS=x2
          else if ((ls1.ne.zero).and. &
                   (ls2.ne.zero)) then
           slope=(ls1-ls2)/(x1-x2)
           if (slope.ne.zero) then
            ht_from_LS=x1-ls1/slope
           else
            print *,"slope invalid"
            stop
           endif
          else
           print *,"ls1 or ls2 invalid"
           stop
          endif
         else
          print *,"x1 or x2 invalid"
          stop
         endif
        else if ((ls1*ls2.ge.zero).or. &
                 ((ls1-ls2)*n1d.le.zero)) then
         ! do nothing
        else
         print *,"ls1 or ls2 invalid"
         stop
        endif  

       else if ((crossing_status.ne.0).or. &
                (-(l+1).lt.lmin)) then
        ! do nothing
       else
        print *,"crossing_status or -(l+1) problem"
        stop
       endif

      enddo ! l=0,ngrow_tower-1

         ! given volume for column, find the interface.
         ! for RZ, if dircrit=1, then column extends to r=0

      if (crossing_status.eq.1) then ! crossing found
       if ((lcrit.lt.lmin).or.(lcrit+1.gt.lmax)) then
        print *,"lcrit invalid"
        stop
       endif
      else if (crossing_status.eq.0) then ! crossing not found
       ! do nothing
      else 
       print *,"crossing_status invalid"
       stop
      endif

      ht_from_VOF=ht_from_LS

      if (crossing_status.eq.1) then
       if (vof_height_function.eq.1) then
        lvof_min=max(lcrit-1,lmin)
        lvof_max=min(lcrit+2,lmax)

        if ((lvof_min.ge.lmin).and. &
            (lvof_max.le.lmax).and. &
            (lcrit.ge.lmin).and. &
            (lcrit+1.le.lmax)) then

         vof_top_sum=zero
         vof_bot_sum=zero

         do l_vof=lmin,lmax

          if (n1d.eq.1) then ! n1d=1 => im material on top
           vof_crit=one-vofdata(l_vof)
          else if (n1d.eq.-1) then ! n1d=-1 => im material on bottom
           vof_crit=vofdata(l_vof)
          else
           print *,"n1d invalid"
           stop
          endif

          if (l_vof.lt.lvof_min) then
           vof_crit=one
          else if ((l_vof.ge.lvof_min).and. &
                   (l_vof.le.lvof_max)) then
           ! do nothing
          else if (l_vof.gt.lvof_max) then
           vof_crit=zero
          else
           print *,"l_vof invalid"
           stop
          endif

          if ((vof_crit.ge.-EPS1).and. &
              (vof_crit.le.VOFTOL_MATERIAL)) then
           vof_crit=zero
          else if ((vof_crit.ge.one-VOFTOL_MATERIAL).and. &
                   (vof_crit.le.one+EPS1)) then
           vof_crit=one
          else if ((vof_crit.gt.zero).and. &
                   (vof_crit.lt.one)) then
           ! do nothing
          else
           print *,"vof_crit invalid: ",vof_crit
           stop
          endif

          if (levelrz.eq.COORDSYS_RZ) then
           if (dircrit.eq.1) then ! horizontal column

            dr=xsten0(2*l_vof+1,dircrit)-xsten0(2*l_vof-1,dircrit)
            dz=dx_col(2)
            if ((dz.gt.zero).and.(dr.gt.zero)) then
              !2 pi r * dr * dz
             volcell=Pi*(xsten0(2*l_vof-1,dircrit)+ &
                         xsten0(2*l_vof+1,dircrit))*dz*dr
            else
             print *,"dz or dr invalid: ",dz,dr
             stop
            endif
             
           else if (dircrit.eq.2) then ! vertical column

            dz=xsten0(2*l_vof+1,dircrit)-xsten0(2*l_vof-1,dircrit)
            dr=dx_col(1)
            if ((dz.gt.zero).and.(dr.gt.zero)) then
             volcell=Pi*two*x_col_avg(1)*dz*dr
            else
             print *,"dz or dr invalid: ",dz,dr
             stop
            endif

           else
            print *,"dircrit invalid: ",dircrit
            stop
           endif

          else if (levelrz.eq.COORDSYS_CARTESIAN) then

           if (SDIM.eq.2) then

            if (dircrit.eq.1) then ! horizontal column
             dr=xsten0(2*l_vof+1,dircrit)-xsten0(2*l_vof-1,dircrit)
             dz=dx_col(2)
             if ((dz.gt.zero).and.(dr.gt.zero)) then
              volcell=dz*dr
             else
              print *,"dz or dr invalid: ",dz,dr
              stop
             endif
            else if (dircrit.eq.2) then ! vertical column

             dz=xsten0(2*l_vof+1,dircrit)-xsten0(2*l_vof-1,dircrit)
             dr=dx_col(1)
             if ((dz.gt.zero).and.(dr.gt.zero)) then
              volcell=dz*dr
             else
              print *,"dz or dr invalid: ",dz,dr
              stop
             endif
            else
             print *,"dircrit invalid: ",dircrit
             stop
            endif

           else if (SDIM.eq.3) then

            dx_norm=xsten0(2*l_vof+1,dircrit)-xsten0(2*l_vof-1,dircrit)
            if (dx_norm.gt.zero) then
             if (dircrit.eq.1) then ! horizontal column
              dx_tan1=dx_col(2)
              dx_tan2=dx_col(SDIM)
             else if (dircrit.eq.2) then ! vertical column
              dx_tan1=dx_col(1)
              dx_tan2=dx_col(SDIM)
             else if ((dircrit.eq.3).and.(SDIM.eq.3)) then
              dx_tan1=dx_col(1)
              dx_tan2=dx_col(2)
             else
              print *,"dircrit invalid: ",dircrit
              stop
             endif
             if ((dx_tan1.gt.zero).and.(dx_tan2.gt.zero)) then
              volcell=dx_norm*dx_tan1*dx_tan2
             else
              print *,"dx_tan1 or dx_tan2 invalid: ",dx_tan1,dx_tan2
              stop
             endif
            else
             print *,"dx_norm invalid: ",dx_norm
             stop
            endif

           else
            print *,"dimension bust"
            stop
           endif
          else
           print *,"levelrz invalid: ",levelrz
           stop
          endif

          vof_top_sum=vof_top_sum+vof_crit*volcell
          vof_bot_sum=vof_bot_sum+volcell
         enddo !l_vof=lmin,lmax

         xbottom_adjusted=xbottom_stencil
         vof_ratio_ht_power=1

         if (levelrz.eq.COORDSYS_RZ) then

          if (dircrit.eq.1) then ! horizontal column

           if (problox.ge.zero) then
            vof_ratio_ht_power=2
            xbottom_adjusted=zero
            dr=xsten0(2*lmin-1,dircrit)-xbottom_adjusted
            dz=dx_col(2)
            if ((dz.gt.zero).and.(dr.ge.zero)) then
             volcell=Pi*(xsten0(2*lmin-1,dircrit)+ &
                         xbottom_adjusted)*dr*dz
            else
             print *,"dz or dr invalid: ",dz,dr
             stop
            endif
           else
            print *,"expecting  problox>=0: ",problox
            stop
           endif  

          else if (dircrit.eq.2) then ! vertical column

           dz=xsten0(2*lmin-1,dircrit)-xbottom_adjusted
           dr=dx_col(1)
           if ((dz.ge.zero).and.(dr.gt.zero)) then
            volcell=two*Pi*x_col_avg(1)*dz*dr
           else
            print *,"dz or dr invalid: ",dz,dr
            stop
           endif

          else
           print *,"dircrit invalid: ",dircrit
           stop
          endif

         else if (levelrz.eq.COORDSYS_CARTESIAN) then

          if (SDIM.eq.2) then

           if (dircrit.eq.1) then ! horizontal column
            dr=xsten0(2*lmin-1,dircrit)-xbottom_adjusted
            dz=dx_col(2)
            if ((dz.gt.zero).and.(dr.ge.zero)) then
             volcell=dr*dz
            else
             print *,"dz or dr invalid: ",dz,dr
             stop
            endif
           else if (dircrit.eq.2) then ! vertical column
            dz=xsten0(2*lmin-1,dircrit)-xbottom_adjusted
            dr=dx_col(1)
            if ((dz.ge.zero).and.(dr.gt.zero)) then
             volcell=dz*dr
            else
             print *,"dz or dr invalid: ",dz,dr
             stop
            endif
           else
            print *,"dircrit invalid: ",dircrit
            stop
           endif

          else if (SDIM.eq.3) then

           dx_norm=xsten0(2*lmin-1,dircrit)-xbottom_adjusted
           if (dx_norm.ge.zero) then
            if (dircrit.eq.1) then ! horizontal column
             dx_tan1=dx_col(2)
             dx_tan2=dx_col(SDIM)
            else if (dircrit.eq.2) then ! vertical column
             dx_tan1=dx_col(1)
             dx_tan2=dx_col(SDIM)
            else if ((dircrit.eq.3).and.(SDIM.eq.3)) then
             dx_tan1=dx_col(1)
             dx_tan2=dx_col(2)
            else
             print *,"dircrit invalid: ",dircrit
             stop
            endif
            if ((dx_tan1.gt.zero).and.(dx_tan2.gt.zero)) then
             volcell=dx_norm*dx_tan1*dx_tan2
            else
             print *,"dx_tan1 or dx_tan2 invalid: ",dx_tan1,dx_tan2
             stop
            endif
           else
            print *,"dx_norm invalid: ",dx_norm
            stop
           endif

          else
           print *,"dimension bust"
           stop
          endif
         else
          print *,"levelrz invalid: ",levelrz
          stop
         endif

         vof_top_sum=vof_top_sum+volcell
         vof_bot_sum=vof_bot_sum+volcell

         if (vof_bot_sum.gt.zero) then
          if (vof_ratio_ht_power.eq.1) then
           if (abs(xbottom_adjusted-xbottom_stencil).le. &
               EPS_8_4*dx(dircrit)) then
            if (vof_top_sum.le.vof_bot_sum) then
             if (xtop_stencil.gt.xbottom_stencil) then
              ht_from_VOF= &
               vof_top_sum*(xtop_stencil-xbottom_adjusted)/vof_bot_sum+ &
                 xbottom_adjusted
              if (abs(ht_from_VOF-xbottom_stencil).le. &
                  EPS_8_4*dx(dircrit)) then
               ht_from_VOF=xbottom_stencil
              else if (abs(ht_from_VOF-xtop_stencil).le. &
                       EPS_8_4*dx(dircrit)) then
               ht_from_VOF=xtop_stencil
              else if ((ht_from_VOF.ge.xbottom_stencil).and. &
                       (ht_from_VOF.le.xtop_stencil)) then
               ! do nothing
              else
               print*,"ht_from_VOF invalid(1): ",ht_from_VOF
               stop
              endif
             else
              print *,"xtop_stencil or xbottom_stencil bad"
              print *,"xtop_stencil ",xtop_stencil
              print *,"xbottom_stencil ",xbottom_stencil
              stop
             endif
            else
             print *,"vof_top_sum or vof_bot_sum bad"
             print *,"vof_top_sum : ",vof_top_sum
             print *,"vof_bot_sum : ",vof_bot_sum
             stop
            endif
           else
            print *,"expecting xbottom_adjusted=xbottom_stencil"
            print *,"xbottom_adjusted ",xbottom_adjusted
            print *,"xbottom_stencil ",xbottom_stencil
            stop
           endif
          else if (vof_ratio_ht_power.eq.2) then
           dr=dx_col(1)
           ht_from_VOF=(xtop_stencil**2)*vof_top_sum/vof_bot_sum
           if (ht_from_VOF.ge.zero) then
            ht_from_VOF=sqrt(ht_from_VOF)
            if (abs(ht_from_VOF-xbottom_stencil).le.EPS_8_4*dr) then
             ht_from_VOF=xbottom_stencil
            else if (abs(ht_from_VOF-xtop_stencil).le.EPS_8_4*dr) then
             ht_from_VOF=xtop_stencil
            else if ((ht_from_VOF.ge.xbottom_stencil).and. &
                     (ht_from_VOF.le.xtop_stencil)) then
             ! do nothing
            else
             print*,"ht_from_VOF invalid: ",ht_from_VOF
             stop
            endif
                    
           else
            print *,"ht_from_VOF invalid: ",ht_from_VOF
            stop
           endif

          else 
           print *,"vof_ratio_ht_power invalid: ",vof_ratio_ht_power
           stop
          endif
         else
          print *,"vof_bot_sum invalid: ",vof_bot_sum
          stop
         endif

        else
         print *,"lcrit,lvof_min, or lvof_max invalid: ", &
          lcrit,lvof_min,lvof_max
         stop
        endif
                ! do nothing
       else if (vof_height_function.eq.0) then
        ! do nothing
       else
        print *,"vof_height_function invalid"
        stop
       endif
      else if (crossing_status.eq.0) then
       ht_from_LS=X_AT_ABS_LSMIN
       ht_from_VOF=X_AT_ABS_LSMIN
      else
       print *,"crossing_status invalid"
       stop
      endif

      return
      end subroutine get_col_ht_LS


! normal points from (-) to (+)   phi=n dot (x-x0) + intercept
! xsten0(0,dir) is cell center (not cell centroid)
! finds grad phi/|grad phi| where grad=(d/dx,d/dy,d/dz) or
! grad=(d/dr,d/dz) or
! grad=(d/dr,d/dtheta,d/dz)
      subroutine find_cut_geom_slope_CLSVOF( &
        ls_mof, &
        lsnormal, &
        lsnormal_valid, &
        ls_intercept, &
        bfact,dx, &
        xsten0,nhalf0, &
        im, &
        dxmax_volume_constraint, & 
        sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: nhalf0
      integer, INTENT(in) :: im
      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: dxmax_volume_constraint
      real(amrex_real) xpoint(sdim)
      real(amrex_real) cutoff,m1,m2
      real(amrex_real) nsimple(sdim)
      real(amrex_real) nn(sdim)
      real(amrex_real) distsimple,dist
      real(amrex_real) LSWT
      real(amrex_real) LSWT_grid
      real(amrex_real) w(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) aa(sdim+1,sdim+1)
      real(amrex_real) xx(sdim+1)
      real(amrex_real) bb(sdim+1)
      integer matstatus
      integer ii,jj,kk
      integer i,j,k
      integer i1,j1,k1
      integer klo,khi
      integer dir

      real(amrex_real), INTENT(in) :: &
            ls_mof(D_DECL(-1:1,-1:1,-1:1),num_materials)

      real(amrex_real), INTENT(out)    :: lsnormal(num_materials,sdim)
      integer, INTENT(out) :: lsnormal_valid(num_materials)
      real(amrex_real), INTENT(out)    :: ls_intercept(num_materials)
      real(amrex_real) dxplus,dxminus
      real(amrex_real) LS_cen,LS_plus,LS_minus
      real(amrex_real) slope_plus,slope_minus
      real(amrex_real) theta_plus,theta_minus
      integer, parameter :: always_use_default=0

#include "mofdata.H"

      if (nhalf0.lt.3) then
       print *,"expecting nhalf0>=3: ",nhalf0
       stop
      endif

      if (dxmax_volume_constraint.gt.zero) then
       ! do nothing
      else
       print *,"dxmax_volume_constraint invalid: ", &
               dxmax_volume_constraint
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid find_cut_geom_slope_CLSVOF ",bfact
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid find_cut_geom_slope_CLSVOF: ",sdim
       stop
      endif

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid find_cut_geom_slope_CLSVOF ",im
       stop
      endif
      if (num_materials.lt.1) then
       print *,"num_materials invalid"
       stop
      endif

      if (sdim.eq.3) then
       klo=-1
       khi=1
       cutoff=four*dxmax_volume_constraint
      else if (sdim.eq.2) then
       klo=0
       khi=0
       cutoff=four*dxmax_volume_constraint
      else
       print *,"dimension bust"
       stop
      endif

      lsnormal_valid(im)=1

      LS_cen=ls_mof(D_DECL(0,0,0),im)

       ! if 0<F<1 then |LS_{center}|<dx*sqrt(3)/2 (3D)
      if (abs(LS_cen).gt.two*dxmax_volume_constraint) then
       lsnormal_valid(im)=0
      else if (abs(LS_cen).le.two*dxmax_volume_constraint) then
       !do nothing
      else
       print *,"LS_cen invalid: ",LS_cen
       stop
      endif

      if (lsnormal_valid(im).eq.1) then

       do dir=1,sdim
        ii=0
        jj=0
        kk=0
        if (dir.eq.1) then
         ii=1
        else if (dir.eq.2) then
         jj=1
        else if ((dir.eq.3).and.(sdim.eq.3)) then
         kk=1
        else
         print *,"dir invalid CLSVOF slope: ",dir
         stop
        endif
        dxplus=xsten0(2,dir)-xsten0(0,dir)
        dxminus=xsten0(0,dir)-xsten0(-2,dir)

        LS_plus=ls_mof(D_DECL(ii,jj,kk),im)
        LS_minus=ls_mof(D_DECL(-ii,-jj,-kk),im)

        ! if 0<F<1 then |LS_{center}|<dx*sqrt(3)/2 (3D)
        if ((abs(LS_plus).le.three*dxmax_volume_constraint).and. &
            (abs(LS_minus).le.three*dxmax_volume_constraint)) then
         !do nothing
        else if ((abs(LS_plus).ge.three*dxmax_volume_constraint).or. &
                 (abs(LS_minus).ge.three*dxmax_volume_constraint)) then
         lsnormal_valid(im)=0
        else
         print *,"LS_plus or LS_minus invalid: ",LS_plus,LS_minus
         stop
        endif
  
        if ((dxplus.gt.zero).and.(dxminus.gt.zero)) then

         slope_plus=(LS_plus-LS_cen)/dxplus
         slope_minus=(LS_cen-LS_minus)/dxminus

         if ((LS_plus*LS_cen.le.zero).and. &
             (LS_minus*LS_cen.gt.zero)) then
          nsimple(dir)=slope_plus
         else if ((LS_plus*LS_cen.gt.zero).and. &
                  (LS_minus*LS_cen.le.zero)) then
          nsimple(dir)=slope_minus
         else if ((LS_plus*LS_cen.le.zero).and. &
                  (LS_minus*LS_cen.le.zero)) then
          if ((LS_plus.eq.zero).and.(LS_cen.eq.zero)) then
           nsimple(dir)=zero
          else if ((LS_minus.eq.zero).and.(LS_cen.eq.zero)) then
           nsimple(dir)=zero
          else if ((abs(LS_plus).gt.zero).or. &
                   (abs(LS_minus).gt.zero).or. &
                   (abs(LS_cen).gt.zero)) then
           theta_plus=abs(LS_plus/(LS_plus-LS_cen)-half)
           theta_minus=abs(LS_minus/(LS_minus-LS_cen)-half)
           if ((theta_plus.le.half).and.(theta_minus.le.half)) then
            if (theta_plus.le.theta_minus) then
             nsimple(dir)=slope_plus
            else if (theta_minus.le.theta_plus) then
             nsimple(dir)=slope_minus
            else
             print *,"theta_plus or theta_minus invalid1"
             print *,"LS_plus,LS_minus,LS_cen ", &
              LS_plus,LS_minus,LS_cen
             print *,"theta_plus,theta_minus ", &
              theta_plus,theta_minus
             stop
            endif
           else
            print *,"theta_plus or theta_minus invalid2"
            print *,"LS_plus,LS_minus,LS_cen ", &
             LS_plus,LS_minus,LS_cen
            print *,"theta_plus,theta_minus ", &
             theta_plus,theta_minus
            stop
           endif
          else
           print *,"LS_plus, LS_minus, or LS_cen invalid: ", &
              LS_plus,LS_minus,LS_cen
           stop
          endif
         else if ((LS_plus*LS_cen.gt.zero).and. &
                  (LS_minus*LS_cen.gt.zero)) then
          if (abs(LS_plus).le.abs(LS_minus)) then
           nsimple(dir)=slope_plus
          else if (abs(LS_plus).ge.abs(LS_minus)) then
           nsimple(dir)=slope_minus
          else
           print *,"LS_plus or LS_minus invalid: ",LS_plus,LS_minus
           stop
          endif
         else
          print *,"LS_plus, LS_minus, or LS_cen invalid: ", &
             LS_plus,LS_minus,LS_cen
          stop
         endif

        else
         print *,"dxplus or dxminus invalid: ",dxminus,dxplus
         stop
        endif

       enddo ! dir=1..sdim

       if (lsnormal_valid(im).eq.1) then

        distsimple=zero
        do dir=1,sdim
         distsimple=distsimple+nsimple(dir)**2
        enddo
        distsimple=sqrt(distsimple)
        if (distsimple.gt.zero) then
         do dir=1,sdim
          nsimple(dir)=nsimple(dir)/distsimple
         enddo
         ls_intercept(im)=LS_cen/distsimple
        else if (distsimple.eq.zero) then
         lsnormal_valid(im)=0
        else
         print *,"distsimple invalid"
         stop
        endif

        if (lsnormal_valid(im).eq.1) then

         if (always_use_default.eq.1) then

          do dir=1,sdim
           lsnormal(im,dir)=nsimple(dir)
          enddo

         else if (always_use_default.eq.0) then

          do k=klo,khi
          do j=-1,1
          do i=-1,1

           LSWT=abs(ls_mof(D_DECL(i,j,k),im))
           if (LSWT.le.three*dxmax_volume_constraint) then
            !do nothing
           else if (LSWT.gt.three*dxmax_volume_constraint) then
            lsnormal_valid(im)=0
           else
            print *,"LSWT invalid: ",LSWT
            stop
           endif

           LSWT_grid=(xsten0(2*i,1)-xsten0(0,1))**2
           LSWT_grid=LSWT_grid+(xsten0(2*j,2)-xsten0(0,2))**2
           if (sdim.eq.3) then
            LSWT_grid=LSWT_grid+(xsten0(2*k,sdim)-xsten0(0,sdim))**2
           endif
           LSWT_grid=sqrt(LSWT_grid)

           w(D_DECL(i,j,k))=hsprime(LSWT_grid,cutoff)*hsprime(LSWT,cutoff)
          enddo
          enddo
          enddo ! i,j,k

          if (lsnormal_valid(im).eq.1) then

           do i=1,sdim+1
            do j=1,sdim+1
             aa(i,j)=zero
            enddo
            bb(i)=zero
           enddo

           do k1=klo,khi
           do j1=-1,1
           do i1=-1,1
            xpoint(1)=xsten0(2*i1,1)-xsten0(0,1)
            xpoint(2)=xsten0(2*j1,2)-xsten0(0,2)
            if (sdim.eq.3) then
             xpoint(sdim)=xsten0(2*k1,sdim)-xsten0(0,sdim)
            endif
     
            do i=1,sdim+1
             if (i.eq.sdim+1) then
              m1=one
             else
              m1=xpoint(i)
             endif

             do j=1,sdim+1
              if (j.eq.sdim+1) then
               m2=one
              else
               m2=xpoint(j)
              endif
              aa(i,j)=aa(i,j)+w(D_DECL(i1,j1,k1))*m1*m2
             enddo ! j=1,sdim+1 
             bb(i)=bb(i)+ &
              w(D_DECL(i1,j1,k1))*m1*ls_mof(D_DECL(i1,j1,k1),im)
            enddo ! i=1,sdim+1
           enddo
           enddo
           enddo ! i1,j1,k1

           call matrix_solve(aa,xx,bb,matstatus,sdim+1)
           if (matstatus.eq.1) then
            dist=zero
            do dir=1,sdim
             nn(dir)=xx(dir)
             dist=dist+nn(dir)**2
            enddo
            dist=sqrt(dist)
            if (dist.gt.zero) then

             do dir=1,sdim
              nn(dir)=nn(dir)/dist
             enddo
             ls_intercept(im)=xx(sdim+1)/dist 

            else if (dist.eq.zero) then

             do dir=1,sdim
              nn(dir)=nsimple(dir)
             enddo

            else
             print *,"dist invalid: ",dist
             stop
            endif

           else if (matstatus.eq.0) then
            do dir=1,sdim
             nn(dir)=nsimple(dir)
            enddo
           else
            print *,"matstatus invalid ",matstatus
            stop
           endif
           do dir=1,sdim
            lsnormal(im,dir)=nn(dir)
           enddo

          else if (lsnormal_valid(im).eq.0) then
           ! do nothing
          else
           print *,"lsnormal_valid(im) invalid"
           stop
          endif

         else 
          print *,"always_use_default invalid ",always_use_default
          stop
         endif

        else if (lsnormal_valid(im).eq.0) then
         ! do nothing
        else
         print *,"lsnormal_valid(im) invalid 15498: ",im,lsnormal_valid(im)
         stop
        endif

       else if (lsnormal_valid(im).eq.0) then
        ! do nothing
       else
        print *,"lsnormal_valid(im) invalid 15505: ",im,lsnormal_valid(im)
        stop
       endif
      else if (lsnormal_valid(im).eq.0) then
       ! do nothing
      else
       print *,"lsnormal_valid(im) invalid 15511: ",im,lsnormal_valid(im)
       stop
      endif

      return
      end subroutine find_cut_geom_slope_CLSVOF

      subroutine nfact(n,nf)
      IMPLICIT NONE

      integer, INTENT(in) :: n
      integer, INTENT(out) :: nf
      integer i

      if (n.le.0) then
       print *,"n invalid"
       stop
      else if (n.ge.1) then
       nf=1
       do i=2,n
        nf=nf*i
       enddo
      else
       print *,"n invalid ",n
       stop
      endif
 
      return
      end subroutine nfact

      subroutine permute_sanity( &
        permutation_array,i,j,ncomp_i,ncomp_j)
      IMPLICIT NONE

      integer, intent(in) :: i,j,ncomp_i,ncomp_j
      integer, intent(in) :: permutation_array(ncomp_i,ncomp_j)
      integer :: same_flag,jj,kk

      do jj=1,ncomp_j
       if ((permutation_array(i,jj).ge.1).and. &
           (permutation_array(i,jj).le.ncomp_j).and. &
           (permutation_array(j,jj).ge.1).and. &
           (permutation_array(j,jj).le.ncomp_j)) then
        !do nothing
       else
        print *,"permutation_array failure"
        stop
       endif
      enddo
   
      do jj=1,ncomp_j-1
      do kk=jj+1,ncomp_j
       if ((permutation_array(i,jj).ne.permutation_array(i,kk)).and. &
           (permutation_array(j,jj).ne.permutation_array(j,kk))) then
        !do nothing
       else
        print *,"permutation_array failure"
        stop
       endif
      enddo
      enddo
      same_flag=1
      do jj=1,ncomp_j
       if (permutation_array(i,jj).ne.permutation_array(j,jj)) then
        same_flag=0
       endif
      enddo
      if (same_flag.eq.0) then
       !do nothing
      else
       print *,"permutation_array failure"
       stop
      endif
     
      return
      end subroutine permute_sanity

       
      subroutine multimaterial_MOF( &
        tessellate, &
        tid_in, &
        bfact,dx, &
        xsten0,nhalf0, &
        mof_verbose, &
        use_ls_data, &
        LS_stencil, &
        xtetlist_vof, &
        nlist_alloc, &
        nmax, &
        mofdata, & !intent(inout)
        vof_super, &
        multi_centroidA, &
        sdim)

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT (IN) :: tessellate
      integer, INTENT (IN) :: tid_in
      integer, INTENT (IN) :: bfact,nhalf0
      integer, INTENT (IN) :: use_ls_data
      integer, INTENT (IN) :: mof_verbose
      integer, INTENT (IN) :: sdim
      integer, INTENT (IN) :: nlist_alloc
      integer, INTENT (IN) :: nmax

       ! D_DECL(i,j,k) = i,j  in 2D
       !               = i,j,k in 3D
      real(amrex_real), INTENT (IN) :: &
            LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)

      real(amrex_real), INTENT (IN), DIMENSION(sdim) :: dx
      real(amrex_real), DIMENSION(sdim) :: dx_extend

      real(amrex_real), INTENT (INOUT), DIMENSION(4,3,nlist_alloc) :: &
              xtetlist_vof
      real(amrex_real), INTENT (IN), DIMENSION(-nhalf0:nhalf0,sdim) :: xsten0
      real(amrex_real), INTENT (INOUT),DIMENSION(num_materials*ngeom_recon) :: &
              mofdata
      real(amrex_real), DIMENSION(num_materials*ngeom_recon) :: &
              mofdata_init
      real(amrex_real), DIMENSION(num_materials*ngeom_recon) :: &
              mofdata_current
      real(amrex_real), DIMENSION(num_materials*ngeom_recon) :: &
              mofdata_save
      real(amrex_real), INTENT (IN), DIMENSION(num_materials) :: vof_super
      real(amrex_real), INTENT (OUT), DIMENSION(num_materials,sdim) :: &
              multi_centroidA
      real(amrex_real), DIMENSION(num_materials,sdim) :: &
              multi_centroidA_current
      real(amrex_real), DIMENSION(num_materials,sdim) :: &
              multi_centroidA_save

      integer imaterial_prev
      integer imaterial
      integer im_opp
      integer vofcomp
      integer vofcomp_prev
      integer dir
      integer loop_counter

      real(amrex_real) vfrac_sum_local( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) vfrac_sum( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer layer_flag
      integer layer_iter

      integer local_num_materials( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer override_target(num_materials)
      integer num_override_history
      integer recon_history(num_materials)
      integer recon_history_save(num_materials)
      integer outer_sweeps
      integer num_outer_sweeps

      integer material_used(num_materials)
      integer num_processed_total
      integer sanity_num_processed_total
      integer layer_sanity

      integer num_processed( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) &
          uncaptured_volume_fraction(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) &
          uncaptured_volume_vof(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) &
          uncaptured_centroid_vof(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX,sdim)

      real(amrex_real) uncaptured_volume_local
      real(amrex_real) uncaptured_centroid_local(sdim)

      integer single_material
      real(amrex_real) remaining_vfrac

      real(amrex_real) nrecon(sdim)
      integer order_algorithm_local(num_materials)
      real(amrex_real) voftest(num_materials)
      integer i1,j1,k1,k1lo,k1hi
      real(amrex_real) dxmax
      real(amrex_real) dxmax_volume_constraint
      real(amrex_real) null_intercept
      real(amrex_real), dimension(D_DECL(:,:,:),:), allocatable :: ls_mof

      real(amrex_real), dimension(:,:), allocatable :: override_normal
      integer, dimension(:), allocatable :: override_normal_valid
      integer, dimension(:), allocatable :: override_order

      real(amrex_real), dimension(:,:), allocatable :: lsnormal
      integer, dimension(:), allocatable :: lsnormal_valid

      real(amrex_real), dimension(:), allocatable :: ls_intercept

      integer, dimension(:,:), allocatable :: permutation_array
      integer permute_c(num_materials)
      integer repeat_count
      integer ihistory
      integer leading_rank
      integer compare_rank
      integer hold_swap
      integer ifirst
      integer ipermute
      integer ipermute_prev
      real(amrex_real) leading_vfrac_sum
      real(amrex_real) mof_err,mof_err_save
      real(amrex_real) xref_mat(sdim)
      real(amrex_real) xact_mat(sdim)
      real(amrex_real) xref_matT(sdim)
      real(amrex_real) xact_matT(sdim)

      real(amrex_real) mag_vec

      integer, PARAMETER :: use_initial_guess=0
      integer :: nDOF_standard
      integer :: nEQN_standard
      integer is_masked(num_materials)
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer at_least_one(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      integer order_sanity(num_materials)
      integer test_order

!      logical :: get_dyn=.false.
      integer :: get_dyn

#ifdef _OPENMP
      integer omp_get_dynamic
#endif

#include "mofdata.H"

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid (multimaterial_MOF): ",tid_in
       stop
      endif

      get_dyn=0
#ifdef _OPENMP
      get_dyn=omp_get_dynamic()
#endif
      if (get_dyn.eq.0) then
       !do nothing
      else
       print *,"get_dyn invalid: ",get_dyn
       stop
      endif

      nDOF_standard=sdim-1
      nEQN_standard=sdim

      do dir=1,sdim
       dx_extend(dir)=xsten0(1,dir)-xsten0(-1,dir)
       if (dx_extend(dir).gt.zero) then
        !do nothing
       else
        print *,"dx_extend invalid dir,dx_extend(dir) ", &
          dir,dx_extend(dir)
        stop
       endif
      enddo !dir=1,sdim

      do imaterial=1,num_materials

       if (is_rigid(imaterial).eq.1) then
        ! order_algorithm is declared in mofdata.H
        if (order_algorithm(imaterial).eq.1) then
         !do nothing
        else
         print *,"expecting order_algorithm(imaterial)=1 if is_rigid"
         stop
        endif
       else if (is_rigid(imaterial).eq.0) then
        !do nothing
       else
        print *,"is_rigid invalid: ", &
           imaterial,is_rigid(imaterial)
        stop
       endif

       if (is_elastic(imaterial).eq.1) then

        ! order_algorithm is declared in mofdata.H
        if ((order_algorithm(imaterial).ge.1).and. &
            (order_algorithm(imaterial).le.num_materials-1)) then
         !do nothing
        else
         print *,"expecting order_algorithm(imaterial)>=1 if is_elastic"
         print *,"and expecting order_algorithm(imaterial)<=nmat-1"
         stop
        endif

        do im_opp=1,num_materials
         if (im_opp.ne.imaterial) then
          if ((is_elastic(im_opp).eq.0).and. &
              (is_rigid(im_opp).eq.0)) then
           if (order_algorithm(imaterial).ge.order_algorithm(im_opp)) then
            print *,"expecting elastic ordering to lead the fluids"
            print *,"order_algorithm=",order_algorithm
            print *,"is_elastic(imaterial) ",imaterial,is_elastic(imaterial)
            print *,"im_opp=",im_opp
            stop
           endif
          endif
         endif
        enddo !im_opp=1,num_materials

       else if (is_elastic(imaterial).eq.0) then
        !do nothing
       else
        print *,"is_elastic invalid: ", &
          imaterial,is_elastic(imaterial)
        stop
       endif

       if (tessellate.eq.TESSELLATE_IGNORE_ISRIGID) then
        if (use_ls_data.eq.0) then
         !do nothing
        else
         print *,"expecting use_ls_data==0"
         stop
        endif
       else if (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC) then
        if (use_ls_data.eq.0) then
         !do nothing
        else
         print *,"expecting use_ls_data==0"
         stop
        endif
       else if (tessellate.eq.TESSELLATE_FLUIDS) then
        ! do nothing; fluids tessellate, rigid|elastic embedded.
       else
        print *,"tessellate invalid8 make_vfrac_sum_ok_base: ",tessellate
        stop
       endif

      enddo ! imaterial=1..num_materials

      call init_local_material_vars( & !subroutine multimaterial_MOF
       is_rigid_local, &
       is_elastic_local, &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
       tessellate, &
       is_proper_layer_local)

      if (bfact.lt.1) then
       print *,"bfact invalid135 multimaterial_mof ",bfact
       stop
      endif

      if (nhalf0.lt.3) then
       print *,"expecting nhalf0>=3: ",nhalf0
       stop
      endif

      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon invalid multimaterial_MOF ",ngeom_recon
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multimaterial_MOF"
       stop
      endif

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid: ",nlist_alloc
       stop
      endif

      if (nmax.eq.geom_nmax) then
       ! do nothing
      else
       print *,"multimaterial_MOF: nmax<>geom_nmax"
       print *,"nmax= ",nmax
       print *,"geom_nmax= ",geom_nmax
       stop
      endif

      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid multimaterial mof"
       print *,"num_materials= ",num_materials
       stop
      endif
      if (nmax.lt.10) then
       print *,"nmax too small in multimaterial_MOF"
       stop
      endif
      if ((use_ls_data.ne.0).and. &
          (use_ls_data.ne.1)) then
       print *,"use_ls_data invalid: ",use_ls_data
       stop
      endif

      call get_dxmax(dx_extend,bfact,dxmax)

      dxmax_volume_constraint=dxmax
      null_intercept=two*bfact*dxmax_volume_constraint

      if (sdim.eq.2) then
       k1lo=0
       k1hi=0
      else if (sdim.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      do imaterial=1,num_materials

       vofcomp=(imaterial-1)*ngeom_recon+1

       voftest(imaterial)=mofdata(vofcomp)

       if ((vof_super(imaterial).ge.-EPS1).and. &
           (vof_super(imaterial).le.one+EPS1)) then
        ! do nothing
       else
        print *,"vof_super out of range: ",imaterial,vof_super(imaterial)
        stop
       endif
       if ((voftest(imaterial).ge.-EPS1).and. &
           (voftest(imaterial).le.one+EPS1)) then
        ! do nothing
       else
        print *,"voftest out of range: ",imaterial,voftest(imaterial)
        stop
       endif

       if ((is_rigid(imaterial).eq.1).or. &
           (is_elastic(imaterial).eq.1)) then

        if (abs(voftest(imaterial)-vof_super(imaterial)).le.EPS12) then
         !do nothing
        else
         print *,"vof_super mismatch with voftest"
         print *,"imaterial=",imaterial
         print *,"is_rigid(imaterial)=",is_rigid(imaterial)
         print *,"is_elastic(imaterial)=",is_elastic(imaterial)
         print *,"voftest=",voftest(imaterial)
         print *,"vof_super=",vof_super(imaterial)
         stop
        endif

       else if ((is_rigid(imaterial).eq.0).and. &
                (is_elastic(imaterial).eq.0)) then

        if (abs(voftest(imaterial)-vof_super(imaterial)).le.EPS12) then
         !do nothing
        else
         print *,"vof_super mismatch with voftest(2)"
         print *,"imaterial=",imaterial
         print *,"is_rigid(imaterial)=",is_rigid(imaterial)
         print *,"is_elastic(imaterial)=",is_elastic(imaterial)
         print *,"voftest=",voftest(imaterial)
         print *,"vof_super=",vof_super(imaterial)
         stop
        endif

       else
        print *,"is_rigid(imaterial) invalid ",is_rigid(imaterial)
        print *,"or is_elastic(imaterial) invalid ", &
           is_elastic(imaterial)
        stop
       endif

      enddo ! imaterial=1,num_materials

      allocate(ls_mof(D_DECL(-1:1,-1:1,-1:1),num_materials))

      allocate(override_normal(num_materials,sdim))
      allocate(override_normal_valid(num_materials))
      allocate(override_order(num_materials))

      allocate(lsnormal(num_materials,sdim))
      allocate(lsnormal_valid(num_materials))

      allocate(ls_intercept(num_materials))

       ! if F<eps or F>1-eps, then moments and vfracs are truncated.
       ! sum of F_fluid=1
       ! sum of F_rigid<=1
      call make_vfrac_sum_ok_base( &
        xsten0, &
        nhalf0, &
        bfact,dx, &
        tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
        mofdata, &
        sdim)

       ! clear flag for all num_materials materials.
       ! vfrac,centroid,order,slope,intercept x num_materials

      do imaterial=1,num_materials

       vofcomp=(imaterial-1)*ngeom_recon+1

       mof_iterations(tid_in+1,imaterial)=0
       mof_errors(tid_in+1,imaterial)=0.0d0
       mof_calls(tid_in+1,imaterial)=0

       if ((mofdata(vofcomp).ge.-EPS1).and. &
           (mofdata(vofcomp).le.zero)) then
        mofdata(vofcomp)=zero
        do dir=1,sdim
         mofdata(vofcomp+dir)=zero  ! centroid=0
        enddo
       else if (mofdata(vofcomp).le.one+EPS1) then
        !do nothing
       else
        print *,"mofdata(vofcomp) invalid ",mofdata
        stop
       endif

       mofdata(vofcomp+sdim+1)=zero  ! order=0
       do dir=1,sdim
        mofdata(vofcomp+sdim+1+dir)=zero  ! slope=0 
       enddo
       mofdata(vofcomp+2*sdim+2)=zero  ! intercept=0

       if (vofcomp+2*sdim+2.eq.imaterial*ngeom_recon) then
        ! do nothing
       else
        print *,"ngeom_recon,vofcomp, or imaterial invalid: ", &
           ngeom_recon,vofcomp,imaterial
        stop
       endif

      enddo !imaterial=1,num_materials

      call Box_volumeFAST( &
       bfact,dx, &
       xsten0,nhalf0, &
       uncaptured_volume_local, &
       uncaptured_centroid_local, &
       sdim)

      if ((VOFTOL_LAYER.gt.VOFTOL_MATERIAL).and. &
          (VOFTOL_MATERIAL.gt.EPS_UNCAPTURED).and. &
          (EPS_UNCAPTURED.gt.INTERCEPT_TOL)) then
       ! do nothing
      else
       print *,"tolerances incorrect"
       print *,"VOFTOL_LAYER ",VOFTOL_LAYER
       print *,"VOFTOL_MATERIAL ",VOFTOL_MATERIAL
       print *,"EPS_UNCAPTURED ",EPS_UNCAPTURED
       print *,"INTERCEPT_TOL ",INTERCEPT_TOL
       stop
      endif

      do imaterial=1,num_materials
       lsnormal_valid(imaterial)=0
       override_normal_valid(imaterial)=0
       override_order(imaterial)=0
      enddo

       !update override_normal_valid

      do imaterial=1,num_materials

       if ((voftest(imaterial).le.VOFTOL_MATERIAL).or. &
           (voftest(imaterial).ge.one-VOFTOL_MATERIAL)) then

        !do nothing

       else if ((voftest(imaterial).ge.VOFTOL_MATERIAL).and. &
                (voftest(imaterial).le.one-VOFTOL_MATERIAL)) then

        mag_vec=zero
        !refvfrac,refcen,order,slope,intercept
        do dir=1,sdim
         override_normal(imaterial,dir)= &
              mofdata((imaterial-1)*ngeom_recon+sdim+2+dir)
         mag_vec=mag_vec+override_normal(imaterial,dir)**2
        enddo
        mag_vec=sqrt(mag_vec)
        if (mag_vec.eq.zero) then
         override_normal_valid(imaterial)=0
         override_order(imaterial)=0
        else if ((mag_vec.ge.one-EPS2).and. &
                 (mag_vec.le.one+EPS2)) then
         do dir=1,sdim
          override_normal(imaterial,dir)=override_normal(imaterial,dir)/mag_vec
         enddo
         override_normal_valid(imaterial)=1
         override_order(imaterial)= &
            NINT(mofdata((imaterial-1)*ngeom_recon+sdim+2))
        else
         print *,"mag_vec invalid: ",mag_vec
         stop
        endif

       else
        print *,"voftest(imaterial) invalid: ",imaterial,voftest(imaterial)
        stop
       endif

      enddo ! imaterial=1,num_materials

      if (use_ls_data.eq.1) then

       do imaterial=1,num_materials

        if ((voftest(imaterial).le.VOFTOL_MATERIAL).or. &
            (voftest(imaterial).ge.one-VOFTOL_MATERIAL)) then

         !do nothing

        else if ((voftest(imaterial).ge.VOFTOL_MATERIAL).and. &
                 (voftest(imaterial).le.one-VOFTOL_MATERIAL)) then

         do k1=k1lo,k1hi
         do j1=-1,1
         do i1=-1,1
          ls_mof(D_DECL(i1,j1,k1),imaterial)= &
           LS_stencil(D_DECL(i1,j1,k1),imaterial) 
         enddo
         enddo
         enddo

         ! in multimaterial_MOF
         ! find n=grad phi/|grad phi| corresponding to "imaterial"
         call find_cut_geom_slope_CLSVOF( &
           ls_mof, & !intent(in)
           lsnormal, & !intent(out)
           lsnormal_valid, & !intent(out)
           ls_intercept, & !intent(out)
           bfact,dx, & !intent(in)
           xsten0,nhalf0, & !intent(in)
           imaterial, & !intent(in)
           dxmax_volume_constraint, & !intent(in)
           sdim) !intent(in)

        else
         print *,"voftest(imaterial) invalid: ",imaterial,voftest(imaterial)
         stop
        endif

       enddo !imaterial=1 ... num_materials

      else if (use_ls_data.eq.0) then

       !do nothing

      else
       print *,"use_ls_data invalid: ",use_ls_data
       stop
      endif


      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       local_num_materials(layer_iter)=0
       vfrac_sum_local(layer_iter)=zero
       vfrac_sum(layer_iter)=zero

       do imaterial=1,num_materials
        vofcomp=(imaterial-1)*ngeom_recon+1
        if (is_proper_layer_local(imaterial,layer_iter).eq.1) then
         local_num_materials(layer_iter)=local_num_materials(layer_iter)+1
         vfrac_sum_local(layer_iter)= &
            vfrac_sum_local(layer_iter)+mofdata(vofcomp)
        endif
        if (is_proper_layer(imaterial,layer_iter).eq.1) then
         vfrac_sum(layer_iter)= &
            vfrac_sum(layer_iter)+mofdata(vofcomp)
        endif
       enddo !imaterial=1,num_materials

       if ((vfrac_sum_local(layer_iter).ge.zero).and. &
           (vfrac_sum_local(layer_iter).le.one+VOFTOL_MATERIAL)) then
        !do nothing
       else
        print *,"layer_iter=",layer_iter
        print *,"vfrac_sum_local invalid: ",vfrac_sum_local
        stop
       endif
       if ((vfrac_sum(layer_iter).ge.zero).and. &
           (vfrac_sum(layer_iter).le.one+VOFTOL_MATERIAL)) then
        !do nothing
       else
        print *,"layer_iter=",layer_iter
        print *,"vfrac_sum invalid: ",vfrac_sum
        stop
       endif

      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      if (local_num_materials(FLUID_LAYER_INDEX)+ &
          local_num_materials(ELASTIC_LAYER_INDEX)+ &
          local_num_materials(RIGID_LAYER_INDEX).ne.num_materials) then
       print *,"local_num_materials invalid ",local_num_materials
       stop
      endif
      if (abs(one-vfrac_sum_local(FLUID_LAYER_INDEX)).le.EPS1) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif

      do imaterial=1,num_materials

       ! order_algorithm is declared in mofdata.H
       order_algorithm_local(imaterial)=order_algorithm(imaterial)

       if (is_elastic(imaterial).eq.1) then
        if ((order_algorithm(imaterial).ge.1).and. &
            (order_algorithm(imaterial).le.num_materials)) then
         ! do nothing
        else
         print *,"order_algorithm(imaterial) invalid (elastic): ",imaterial, &
          order_algorithm(imaterial)
         stop
        endif
       else if (is_rigid(imaterial).eq.1) then
        if ((order_algorithm(imaterial).ge.1).and. &
            (order_algorithm(imaterial).le.1)) then
         ! do nothing
        else
         print *,"order_algorithm(imaterial) invalid (rigid): ",imaterial, &
          order_algorithm(imaterial)
         stop
        endif
       else if ((is_elastic(imaterial).eq.0).and. &
                (is_rigid(imaterial).eq.0)) then
        if ((order_algorithm(imaterial).ge.1).and. &
            (order_algorithm(imaterial).le.num_materials+1)) then
         ! do nothing
        else
         print *,"order_algorithm(imaterial) invalid (fluid): ",imaterial, &
          order_algorithm(imaterial)
         stop
        endif
       else
        print *,"is_rigid(imaterial) invalid ", &
               imaterial,is_rigid(imaterial)
        print *,"or is_elastic(imaterial) invalid ", &
               imaterial,is_elastic(imaterial)
        stop
       endif

      enddo !imaterial=1,num_materials

      do imaterial=1,num_materials

       vofcomp=(imaterial-1)*ngeom_recon+1

       if (override_normal_valid(imaterial).eq.1) then
        if (override_order(imaterial).eq.0) then
         !do nothing
        else if ((override_order(imaterial).ge.1).and. &
                 (override_order(imaterial).le.num_materials+1)) then
         order_algorithm_local(imaterial)=override_order(imaterial)
        else
         print *,"override_order invalid: ",override_order
         print *,"imaterial=",imaterial
         stop
        endif
       else if (override_normal_valid(imaterial).eq.0) then
        !do nothing
       else
        print *,"override_normal_valid invalid ",override_normal_valid
        print *,"imaterial=",imaterial
        stop
       endif 
        
       if (mofdata(vofcomp).ge.one-VOFTOL_MATERIAL) then
        order_algorithm_local(imaterial)=num_materials+1
       else if (mofdata(vofcomp).le.VOFTOL_MATERIAL) then
        order_algorithm_local(imaterial)=num_materials+1  
       else if ((mofdata(vofcomp).gt.VOFTOL_MATERIAL).and. &
                (mofdata(vofcomp).lt.one-VOFTOL_MATERIAL)) then
        ! do nothing
       else
        print *,"mofdata(vofcomp) invalid: ",vofcomp,mofdata(vofcomp)
        stop
       endif

      enddo ! imaterial=1,num_materials

      do imaterial=1,num_materials*ngeom_recon
       mofdata_init(imaterial)=mofdata(imaterial)
      enddo

      num_outer_sweeps=1
      outer_sweeps=0
      num_override_history=0
      do imaterial=1,num_materials
       override_target(imaterial)=0
      enddo

      do while (outer_sweeps.lt.num_outer_sweeps)

       do imaterial=1,num_materials*ngeom_recon
        mofdata(imaterial)=mofdata_init(imaterial)
       enddo

       do imaterial=1,num_materials
        material_used(imaterial)=0
        recon_history(imaterial)=0
        do dir=1,sdim
         multi_centroidA(imaterial,dir)=zero
        enddo
       enddo
       num_processed_total=0
       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
        num_processed(layer_iter)=0
        uncaptured_volume_fraction(layer_iter)=one
        uncaptured_volume_vof(layer_iter)=uncaptured_volume_local
        do dir=1,sdim
         uncaptured_centroid_vof(layer_iter,dir)=uncaptured_centroid_local(dir)
        enddo !dir=1,sdim
       enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX


       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

        at_least_one(layer_iter)=0
        do imaterial=1,num_materials
         if (is_proper_layer_local(imaterial,layer_iter).eq.1) then
          at_least_one(layer_iter)=at_least_one(layer_iter)+1
         endif
        enddo

        if ((at_least_one(layer_iter).ge.1).and. &
            (at_least_one(layer_iter).le.num_materials).and. &
            (vfrac_sum_local(layer_iter).gt.zero)) then

         call init_layer_flag( &  !multimaterial_MOF
          at_least_one(layer_iter), &
          vfrac_sum, &
          vfrac_sum_local, &
          layer_flag, & !intent(out)
          layer_iter, & !intent(in)
          tessellate, &
          tessellate)

         call init_mask_flag( &
          layer_flag, &
          is_masked)

         loop_counter=0
         do while (continue_material_processing( &
              loop_counter, &
              local_num_materials(layer_iter), &
              num_processed(layer_iter), &
              uncaptured_volume_fraction(layer_iter), &
              vfrac_sum_local(layer_iter), &
              uncaptured_volume_vof(layer_iter)).eq.1) 

          call check_for_single_material( & !multimaterial_MOF
           single_material, & !intent(out)
           remaining_vfrac, & !intent(out)
           num_processed_total, & !intent(out)
           material_used, & !intent(in)
           is_masked, & !intent(in)
           mofdata, &
           mofdata, &  !intent(inout)
           EPS_UNCAPTURED, &
           uncaptured_volume_fraction, & !intent(in)
           layer_iter, &
           tessellate, &
           tessellate, &
           sdim)

          sanity_num_processed_total=0
          do layer_sanity=RIGID_LAYER_INDEX,layer_iter
           sanity_num_processed_total=sanity_num_processed_total+ &
             num_processed(layer_sanity)
          enddo
          if (sanity_num_processed_total.eq.num_processed_total) then
           !do nothing
          else
           print *,"num_processed_total invalid ",num_processed_total
           print *,"sanity_num_processed_total ",sanity_num_processed_total
           stop
          endif

          if ((single_material.ge.1).and. &
              (single_material.le.num_materials).and. &
              (uncaptured_volume_fraction(layer_iter).eq.one).and. &
              (remaining_vfrac.le.EPS_UNCAPTURED)) then

           if (num_processed(layer_iter).eq.0) then
            !do nothing
           else
            print *,"expecting num_processed(layer_iter)=0 ", &
              layer_iter,num_processed(layer_iter)
            stop
           endif
 
           uncaptured_volume_vof(layer_iter)=zero
           uncaptured_volume_fraction(layer_iter)=zero

           num_processed(layer_iter)=num_processed(layer_iter)+1
           num_processed_total=num_processed_total+1
           material_used(single_material)=1
           recon_history(num_processed_total)=single_material

           vofcomp=(single_material-1)*ngeom_recon+1
           mofdata(vofcomp+sdim+1)=one  ! order=1
           do dir=1,sdim
            nrecon(dir)=zero  
           enddo
           nrecon(1)=one ! null slope=(1 0 0) 

           !null_intercept=two*bfact*dxmax_volume_constraint
           mofdata(vofcomp+2*sdim+2)=null_intercept
           do dir=1,sdim
            mofdata(vofcomp+sdim+1+dir)=nrecon(dir)
            multi_centroidA(single_material,dir)=zero  ! cell is full
           enddo 

          else if ((single_material.eq.0).or. &
                   (remaining_vfrac.ge.EPS_UNCAPTURED).or. &
                   (uncaptured_volume_fraction(layer_iter).lt.one)) then

           call individual_MOF( &
            layer_flag, &
            tid_in, &
            ls_mof, &
            lsnormal, &
            lsnormal_valid, &
            override_normal, &
            override_normal_valid, &
            bfact,dx, &
            xsten0,nhalf0, &
            order_algorithm_local, &
            xtetlist_vof, & !intent(out)
            nlist_alloc, & !intent(in)
            nmax, &
            mofdata, &
            loop_counter, & 
            uncaptured_volume_fraction(layer_iter), &
            vfrac_sum_local(layer_iter), &
            num_processed_total, &
            num_processed(layer_iter), &
            num_override_history, &
            override_target, &
            recon_history, &
            material_used, &
            uncaptured_volume_vof(layer_iter), &
            multi_centroidA, &
            sdim)

          else
           print *,"single_material invalid ",single_material
           print *,"or remaining_vfrac invalid: ",remaining_vfrac
           print *,"or uncaptured_volume_fraction invalid: ", &
            uncaptured_volume_fraction
           stop
          endif

          loop_counter=loop_counter+1
         enddo !while continue_material_processing==1

        else if ((at_least_one(layer_iter).eq.0).or. &
                 (vfrac_sum_local(layer_iter).eq.zero)) then
         !do nothing
        else
         print *,"at_least_one or vfrac_sum_local invalid ", &
           at_least_one,vfrac_sum_local
         stop
        endif

       enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

        do imaterial=1,num_materials
         order_sanity(imaterial)=0
        enddo
        do imaterial=1,num_materials
         if (is_proper_layer_local(imaterial,layer_iter).eq.1) then
          vofcomp=(imaterial-1)*ngeom_recon+1
          if ((mofdata(vofcomp).ge.EPS2).and. &
              (mofdata(vofcomp).le.one-EPS2)) then
           test_order=NINT(mofdata(vofcomp+sdim+1))
           if ((test_order.ge.1).and. &
               (test_order.le.at_least_one(layer_iter))) then
            if (material_used(imaterial).eq.test_order) then
             !do nothing
            else
             print *,"material_used invalid ",material_used
             print *,"imaterial= ",imaterial
             print *,"test_order= ",test_order
             stop
            endif
            if (order_sanity(test_order).eq.0) then
             order_sanity(test_order)=imaterial
            else
             print *,"expecting order_sanity(test_order).eq.0"
             stop
            endif
           else
            print *,"expecting test_order in proper range: ",test_order
            stop
           endif
          else if ((mofdata(vofcomp).ge.-EPS1).and. &
                   (mofdata(vofcomp).le.EPS2)) then
           !do nothing
          else if ((mofdata(vofcomp).ge.one-EPS2).and. &
                   (mofdata(vofcomp).le.one+EPS1)) then
           !do nothing
          else
           print *,"mofdata(vofcomp) invalid ",mofdata(vofcomp)
           stop
          endif
         else if (is_proper_layer_local(imaterial,layer_iter).eq.0) then
          !do nothing
         else
          print *,"is_proper_layer_local(imaterial,layer_iter) invalid"
          stop
         endif
        enddo !imaterial=1,num_materials

       enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       do imaterial=1,ngeom_recon*num_materials
        mofdata_current(imaterial)=mofdata(imaterial)
       enddo
       do imaterial=1,num_materials
       do dir=1,sdim
        multi_centroidA_current(imaterial,dir)=multi_centroidA(imaterial,dir)
       enddo
       enddo

       repeat_count=0
       imaterial=recon_history(1)

       if ((imaterial.ge.1).and.(imaterial.le.num_materials)) then
        !do nothing
       else
        print *,"imaterial invalid ",imaterial
        print *,"recon_history ",recon_history
        print *,"mofdata_init ",mofdata_init
        print *,"mofdata ",mofdata
        print *,"mofdata_current ",mofdata_current
        print *,"outer_sweeps ",outer_sweeps
        print *,"num_outer_sweeps ",num_outer_sweeps
        do imaterial=1,num_materials
         vofcomp=(imaterial-1)*ngeom_recon+1
         print *,"vofcomp ",vofcomp
         print *,"mofdata_init(vofcomp) ",mofdata_init(vofcomp)
         print *,"mofdata(vofcomp) ",mofdata(vofcomp)
         print *,"mofdata_current(vofcomp) ",mofdata_current(vofcomp)
        enddo
        stop
       endif

       vofcomp=(imaterial-1)*ngeom_recon+1
       leading_vfrac_sum=mofdata_current(vofcomp)

       do ihistory=2,num_processed_total
        imaterial_prev=recon_history(ihistory-1)
        vofcomp_prev=(imaterial_prev-1)*ngeom_recon+1

        leading_rank=get_material_rank(imaterial_prev, &
           order_algorithm_local(imaterial_prev))

        imaterial=recon_history(ihistory)
        vofcomp=(imaterial-1)*ngeom_recon+1

        compare_rank=get_material_rank(imaterial, &
           order_algorithm_local(imaterial))

        if (compare_rank.eq.leading_rank) then
         if (override_normal_valid(imaterial).eq.0) then
          !do nothing
         else
          print *,"expecting override_normal_valid(imaterial).eq.0"
          print *,"override_normal_valid=",override_normal_valid
          print *,"order_algorithm_local= ",order_algorithm_local
          stop
         endif
         if (override_normal_valid(imaterial_prev).eq.0) then
          !do nothing
         else
          print *,"expecting override_normal_valid(imaterial_prev).eq.0"
          print *,"order_algorithm_local= ",order_algorithm_local
          stop
         endif

         if ((mofdata_current(vofcomp).ge.VOFTOL_ORDERING).and. &
             (mofdata_current(vofcomp_prev).ge.VOFTOL_ORDERING)) then

          if (repeat_count.eq.0) then
           leading_vfrac_sum=leading_vfrac_sum+mofdata_current(vofcomp)
          endif

         else if ((mofdata_current(vofcomp).lt.VOFTOL_ORDERING).or. &
                  (mofdata_current(vofcomp_prev).lt.VOFTOL_ORDERING)) then

          if (repeat_count.eq.0) then
           repeat_count=ihistory-1
          endif

         else
          print *,"mofdata_current invalid ",mofdata_current
          stop
         endif

        else if (compare_rank.gt.leading_rank) then
         if (repeat_count.eq.0) then
          repeat_count=ihistory-1
         endif
        else
         print *,"leading_rank or compare_rank invalid ", &
          leading_rank,compare_rank
         stop
        endif
       enddo ! ihistory=2,num_processed_total

       if (repeat_count.eq.0) then
        repeat_count=num_processed_total
       else if ((repeat_count.ge.1).and. &
                (repeat_count.lt.num_processed_total)) then
        !do nothing
       else
        print *,"repeat_count invalid ",repeat_count
        stop
       endif

       mof_err=zero
       do ihistory=1,repeat_count
        imaterial=recon_history(ihistory)
        if ((imaterial.ge.1).and. &
            (imaterial.le.num_materials)) then
         vofcomp=(imaterial-1)*ngeom_recon+1
         do dir=1,sdim
          xref_mat(dir)=mofdata_current(vofcomp+dir)
          xact_mat(dir)=multi_centroidA_current(imaterial,dir)
         enddo
         call RT_transform_offset(xref_mat,uncaptured_centroid_local,xref_matT)
         call RT_transform_offset(xact_mat,uncaptured_centroid_local,xact_matT)
         do dir=1,sdim
          mof_err = mof_err + &
            mofdata_current(vofcomp)*((xref_matT(dir)-xact_matT(dir))**2)
         enddo ! dir
        else
         print *,"imaterial invalid ",imaterial
         print *,"ihistory ",ihistory
         stop
        endif
       enddo !ihistory=1,repeat_count

       if (outer_sweeps.eq.0) then

        do imaterial=1,num_materials
         recon_history_save(imaterial)=recon_history(imaterial)
        enddo
        mof_err_save=mof_err
        do imaterial=1,ngeom_recon*num_materials
         mofdata_save(imaterial)=mofdata_current(imaterial)
        enddo
        do imaterial=1,num_materials
        do dir=1,sdim
         multi_centroidA_save(imaterial,dir)= &
            multi_centroidA_current(imaterial,dir)
        enddo
        enddo

        num_override_history=repeat_count
        if ((leading_vfrac_sum.ge.zero).and. &
            (leading_vfrac_sum.le.one-VOFTOL_LAYER)) then
         num_override_history=0
        else if ((leading_vfrac_sum.ge.one-VOFTOL_LAYER).and. &
                 (leading_vfrac_sum.le.one+EPS1)) then
         !do nothing
        else
         print *,"leading_vfrac_sum invalid ",leading_vfrac_sum
         stop
        endif

        if ((num_override_history.ge.0).and. &
            (num_override_history.le.2)) then
         !do nothing
        else if ((num_override_history.ge.3).and. &
                 (num_override_history.le.num_materials)) then
         call nfact(num_override_history,num_outer_sweeps)
         if (num_outer_sweeps.ge.6) then
          !do nothing
         else
          print *,"nfact invalid"
          print *,"num_outer_sweeps ",num_outer_sweeps
          print *,"num_override_history ",num_override_history
          stop
         endif
         allocate(permutation_array(num_outer_sweeps,num_override_history))
         do ihistory=1,num_override_history
          do ipermute=1,2
           permutation_array(ipermute,ihistory)=ihistory
          enddo
         enddo 
          !Heap's non recursive algorithm (B.R. Heap 1963)
         do imaterial=1,num_materials
          permute_c(imaterial)=0
         enddo
         ihistory=1
         ipermute=2
         do while (ihistory.lt.num_override_history)
          if (permute_c(ihistory+1).lt.ihistory) then

           if ((ipermute.ge.2).and.(ipermute.le.num_outer_sweeps)) then
            !do nothing
           else
            print *,"ipermute out of range ",ipermute
            print *,"num_outer_sweeps ",num_outer_sweeps
            stop
           endif
           
           if (2*(ihistory/2).eq.ihistory) then
            ifirst=1
            hold_swap=permutation_array(ipermute,ifirst)
            permutation_array(ipermute,ifirst)= &
               permutation_array(ipermute,ihistory+1)
            permutation_array(ipermute,ihistory+1)= &
               hold_swap
           else
            ifirst=permute_c(ihistory+1)+1
            hold_swap=permutation_array(ipermute,ifirst)
            permutation_array(ipermute,ifirst)= &
               permutation_array(ipermute,ihistory+1)
            permutation_array(ipermute,ihistory+1)= &
               hold_swap
           endif
           do ipermute_prev=1,ipermute-1
            call permute_sanity( &
              permutation_array, &
              ipermute_prev,ipermute, &
              num_outer_sweeps,num_override_history)
           enddo
           ipermute=ipermute+1
           if (ipermute.le.num_outer_sweeps) then
            do imaterial=1,num_override_history
             permutation_array(ipermute,imaterial)= &
               permutation_array(ipermute-1,imaterial)
            enddo
           else if (ipermute.eq.num_outer_sweeps+1) then
            !do nothing
           else
            print *,"ipermute out of range ",ipermute
            print *,"num_outer_sweeps ",num_outer_sweeps
            stop
           endif
           permute_c(ihistory+1)=permute_c(ihistory+1)+1
           ihistory=1
          else
           permute_c(ihistory+1)=0
           ihistory=ihistory+1
          endif
         enddo ! while (ihistory.lt.num_override_history)
        else
         print *,"num_override_history invalid ",num_override_history
         stop
        endif

       else if (outer_sweeps.gt.0) then

        do ihistory=1,num_override_history
         if (override_target(ihistory).eq.recon_history(ihistory)) then
          !do nothing
         else
          print *,"override_target invalid ",override_target
          print *,"or recon_history invalid ",recon_history
          print *,"repeat_count ",repeat_count
          print *,"num_override_history ",num_override_history
          stop
         endif
        enddo !ihistory=1,num_override_history

        if (repeat_count.eq.num_override_history) then
         !do nothing
        else
         print *,"repeat_count changed ",repeat_count
         print *,"num_override_history ",num_override_history
         stop
        endif

        if (mof_err.lt.mof_err_save) then
         mof_err_save=mof_err
         do imaterial=1,ngeom_recon*num_materials
          mofdata_save(imaterial)=mofdata_current(imaterial)
         enddo
         do imaterial=1,num_materials
         do dir=1,sdim
          multi_centroidA_save(imaterial,dir)= &
            multi_centroidA_current(imaterial,dir)
         enddo
         enddo
        else if (mof_err.ge.mof_err_save) then
         !do nothing
        else
         print *,"mof_err invalid ",mof_err
         print *,"mof_err_save ",mof_err_save
         stop
        endif
 
       else 
        print *,"outer_sweeps invalid ",outer_sweeps
        stop
       endif

       outer_sweeps=outer_sweeps+1

       if (num_outer_sweeps.ge.6) then
        do imaterial=1,num_materials
         override_target(imaterial)=0
        enddo 
        do ihistory=1,num_override_history
         override_target(ihistory)= &
          recon_history_save(permutation_array(outer_sweeps,ihistory))
        enddo
       else if (num_outer_sweeps.eq.1) then
        !do nothing
       else
        print *,"num_outer_sweeps invalid ",num_outer_sweeps
        stop
       endif

      enddo !while(outer_sweeps.lt.num_outer_sweeps)

      if (num_outer_sweeps.eq.1) then
       !do nothing
      else if (num_outer_sweeps.ge.6) then
       deallocate(permutation_array)
      else
       print *,"num_outer_sweeps invalid ",num_outer_sweeps
       stop
      endif

      do imaterial=1,ngeom_recon*num_materials
       mofdata(imaterial)=mofdata_save(imaterial)
      enddo
      do imaterial=1,num_materials
      do dir=1,sdim
       multi_centroidA(imaterial,dir)= &
            multi_centroidA_save(imaterial,dir)
      enddo
      enddo

      if (mof_verbose.eq.1) then
       print *,"AFTER AFTER"
       print *,"nmax = ",nmax
       print *,"levelrz = ",levelrz
       print *,"num_materials = ",num_materials
       print *,"sdim = ",sdim
       print *,"ngeom_recon = ",ngeom_recon
       do imaterial=1,num_materials*ngeom_recon
        print *,"i,mofdata ",imaterial,mofdata(imaterial)
       enddo
       do imaterial=1,num_materials
        print *,"imaterial,order_algorithm ",imaterial, &
         order_algorithm(imaterial)
        print *,"imaterial,order_algorithm_local ",imaterial, &
         order_algorithm_local(imaterial)
       enddo
       do dir=1,sdim
        print *,"dir,xsten0(0,dir) ",dir,xsten0(0,dir)
       enddo
       print *,"MOFITERMAX ",MOFITERMAX
      else if (mof_verbose.eq.0) then
       ! do nothing
      else
       print *,"mof_verbose invalid in multimaterial_MOF 2"
       print *,"mof_verbose= ",mof_verbose
       stop
      endif

      deallocate(ls_mof)

      deallocate(override_normal)
      deallocate(override_normal_valid)
      deallocate(override_order)

      deallocate(lsnormal)
      deallocate(lsnormal_valid)
      deallocate(ls_intercept)

      return
      end subroutine multimaterial_MOF


      subroutine diagnostic_MOF(sdim,nmax)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module


      IMPLICIT NONE

      integer, parameter :: tid_diagnostic=0

      integer, INTENT(in) :: nmax

      integer, INTENT(in) :: sdim
      integer :: nDOF_standard

      real(amrex_real) xtetlist(4,3,nmax)
      real(amrex_real) dx(sdim) 

      real(amrex_real) mofdata(num_materials*(2*sdim+3))
      real(amrex_real) vof_super(num_materials)

      real(amrex_real) multi_centroidA(num_materials,sdim)
      real(amrex_real) angle(sdim-1)
      real(amrex_real) xpoint(sdim)
      integer nrecon
      integer, parameter :: nsamples=10000 !10000 
      integer im,ntry,iangle,vofcomp
      integer inode
      integer dir,dir2
      integer seed
      real(amrex_real) max_mof_error,moferror
      real(amrex_real) nslope(sdim)
      real(amrex_real) nslope2(sdim)
      real(amrex_real) intercept,intercept2
      real(amrex_real) volcut,volcut2,areacut,areacut2
      real(amrex_real) cencut(sdim)
      real(amrex_real) cencut2(sdim)
      real(amrex_real) xtet_domain(sdim+1,sdim)
      real(amrex_real) volcell_hex
      real(amrex_real) cencell_hex(sdim)
      real(amrex_real) xref_mat(sdim)
      real(amrex_real) xact_mat(sdim)
      real(amrex_real) xref_matT(sdim)
      real(amrex_real) xact_matT(sdim)
      real(amrex_real) avgiter
      real(amrex_real) avgerror
      real(amrex_real) total_calls(num_materials)
      real(amrex_real) total_iterations(num_materials)
      real(amrex_real) total_errors(num_materials)
      integer max_iterations(num_materials)
      integer, parameter :: mof_verbose=0
      integer, parameter :: use_ls_data=0
      integer isten
      integer, parameter :: bfact=1
      real(amrex_real), dimension(D_DECL(:,:,:),:), allocatable :: LS_stencil
      integer i1,j1,k1,k1lo,k1hi
      integer, PARAMETER :: nhalf0=3
      real(amrex_real), PARAMETER :: shrink_factor=0.333333333333d0
      real(amrex_real) xsten0_LS(-nhalf0:nhalf0,sdim) 
      real(amrex_real) xsten0_hex(-nhalf0:nhalf0,sdim) 
      real(amrex_real) xsten0_recon(-nhalf0:nhalf0,sdim) 
      real(amrex_real) volcell_recon
      real(amrex_real) cencell_recon(sdim)
      real(amrex_real) LS_local
      real(amrex_real) x_stencil(sdim)

      integer, parameter :: tessellate=TESSELLATE_FLUIDS
      integer, parameter :: shapeflag=0

#include "mofdata.H"

      nDOF_standard=sdim-1

      print *,"in diagnostic_MOF"
      print *,"DO NOT RUN THIS TEST WITH MULTIPLE THREADS"

      if (geom_nthreads.eq.1) then
       !do nothing
      else
       print *,"expecting geom_nthreads=1"
       stop
      endif

      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon invalid"
       print *,"sdim (parameter) = ",sdim
       print *,"ngeom_recon (global) = ",ngeom_recon
       stop
      endif
      nrecon=num_materials*ngeom_recon
      allocate(LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials))
      if (sdim.eq.2) then
       k1lo=0
       k1hi=0
      else if (sdim.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif
      do k1=k1lo,k1hi
      do j1=-1,1
      do i1=-1,1
      do im=1,num_materials
       LS_stencil(D_DECL(i1,j1,k1),im)=zero
      enddo
      enddo
      enddo
      enddo

      if (nmax.ne.POLYGON_LIST_MAX) then
       print *,"nmax invalid diagnostic MOF nmax=",nmax
       print *,"POLYGON_LIST_MAX=",POLYGON_LIST_MAX
       stop
      endif
      if ((num_materials.ge.2).and. &
          (num_materials.le.MAX_NUM_MATERIALS)) then
       ! do nothing
      else
       print *,"num_materials not initialized properly"
       print *,"num_materials=",num_materials
       stop
      endif

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif

      if (levelrz.ne.COORDSYS_CARTESIAN) then
       print *,"levelrz should be 0 for this sanity check"
       stop
      endif

      do dir=1,sdim
       dx(dir)=one
      enddo

      ! tetrahedra "origin node" mapped to (-1/2,-1/2,-1/2)
      do inode=1,sdim+1
      do dir=1,sdim
       xtet_domain(inode,dir)=zero
      enddo
      enddo
      do dir=1,sdim
       xtet_domain(dir+1,dir)=one
      enddo
      do inode=1,sdim+1
      do dir=1,sdim
       xtet_domain(inode,dir)=xtet_domain(inode,dir)-half
      enddo
      enddo

      do dir=1,sdim
      do isten=-nhalf0,nhalf0
       xsten0_hex(isten,dir)=isten*half*dx(dir)
      enddo
      enddo ! dir

      call Box_volumeFAST( &
       bfact,dx, &
       xsten0_hex,nhalf0, &
       volcell_hex, &
       cencell_hex, &
       sdim)

      if (volcell_hex.ge.zero) then
       !do nothing
      else
       print *,"volcell_hex invalid: ",volcell_hex
       stop
      endif

      volcell_recon=volcell_hex

      ! do test for [-1/2,1/2]^sdim cube
      do dir=1,sdim
       cencell_recon(dir)=cencell_hex(dir)
       do isten=-nhalf0,nhalf0
        xsten0_LS(isten,dir)=xsten0_hex(isten,dir)
        xsten0_recon(isten,dir)=xsten0_LS(isten,dir)
       enddo
      enddo ! dir

      do im=1,num_materials
       total_calls(im)=zero
       total_iterations(im)=zero
       total_errors(im)=zero
       max_iterations(im)=0
      enddo

      seed=86456
#if (USERAND==1)
      call srand(seed)
#else
      print *,"set userand (caps) = 1"
      stop
#endif

      max_mof_error=zero

      do ntry=1,nsamples

       ! 0<=rand()<=1
       do iangle=1,nDOF_standard
#if (USERAND==1)
        angle(iangle)=rand()
#else
        print *,"set userand (caps) = 1"
        stop
#endif
        angle(iangle)=angle(iangle)*two*Pi 
       enddo ! iangle

       call angle_to_slope(angle,nslope,sdim)

        ! phi=n dot (x-xpoint)= n dot (x-xcell) +intercept
        ! intercept=n dot (xcell-xpoint)

       do dir=1,sdim
        nslope2(dir)=-nslope(dir)

        ! 0<=rand()<=1
#if (USERAND==1)
        xpoint(dir)=rand()
#else
        print *,"set userand (caps) = 1"
        stop
#endif
       enddo !dir=1..sdim

       ! 0<=rand()<=1

       do dir=1,sdim
        ! make: -shrink_factor/2 <= xpoint <= shrink_factor/2
        xpoint(dir)=(xpoint(dir)-half)*shrink_factor
       enddo

       ! phi=n dot (x-xpoint)= n dot (x-xcell) +intercept
       ! intercept=n dot (xcell-xpoint)

       intercept=zero
       do dir=1,sdim
        intercept=intercept+nslope(dir)*(xsten0_LS(0,dir)-xpoint(dir))
       enddo  ! dir
       intercept2=-intercept

       call fast_cut_cell_intersection( &
         bfact,dx, &
         xsten0_LS,nhalf0, &
         nslope,intercept, &
         volcut,cencut, &
         areacut, &
         xsten0_LS,nhalf0, &
         xtet_domain,shapeflag, &
         sdim)
       call fast_cut_cell_intersection( &
        bfact,dx, &
        xsten0_LS,nhalf0, &
        nslope2,intercept2, &
        volcut2,cencut2, &
        areacut2, &
        xsten0_LS,nhalf0, &
        xtet_domain,shapeflag, &
        sdim)

       do dir2=1,nrecon
        mofdata(dir2)=zero 
       enddo
       do im=1,num_materials
        vof_super(im)=zero
       enddo

       im=1
       vofcomp=(im-1)*(2*sdim+3)+1
       mofdata(vofcomp)=volcut/volcell_recon
       vof_super(im)=mofdata(vofcomp)

       do dir=1,sdim
        mofdata(vofcomp+dir)=cencut(dir)-cencell_recon(dir)
       enddo

       im=2
       vofcomp=(im-1)*(2*sdim+3)+1
       mofdata(vofcomp)=volcut2/volcell_recon
       vof_super(im)=mofdata(vofcomp)

       do dir=1,sdim
        mofdata(vofcomp+dir)=cencut2(dir)-cencell_recon(dir)
       enddo

       do k1=k1lo,k1hi
       do j1=-1,1
       do i1=-1,1
        dir=1
        x_stencil(dir)=xsten0_recon(i1,dir)
        dir=2
        x_stencil(dir)=xsten0_recon(j1,dir)
        if (sdim.eq.2) then
         ! do nothing
        else if (sdim.eq.3) then
         dir=sdim
         x_stencil(dir)=xsten0_recon(k1,dir)
        else
         print *,"sdim invalid"
         stop
        endif

        LS_local=intercept
        do dir=1,sdim
         LS_local=LS_local+nslope(dir)*(x_stencil(dir)-xsten0_LS(0,dir))
        enddo

        im=1
        LS_stencil(D_DECL(i1,j1,k1),im)=LS_local
        im=2
        LS_stencil(D_DECL(i1,j1,k1),im)=-LS_local
       enddo
       enddo
       enddo

       call multimaterial_MOF( &
        tessellate , &  !TESSELLATE_FLUIDS
        tid_diagnostic, &
        bfact,dx, &
        xsten0_recon,nhalf0, &
        mof_verbose, &
        use_ls_data, &
        LS_stencil, &
        xtetlist, &
        nmax, &
        nmax, &
        mofdata, & !intent(inout)
        vof_super, &
        multi_centroidA, &
        sdim)

       moferror=zero
       do im=1,num_materials
        vofcomp=(im-1)*(2*sdim+3)+1

        do dir=1,sdim
         xref_mat(dir)=mofdata(vofcomp+dir)
         xact_mat(dir)=multi_centroidA(im,dir)
        enddo
        call RT_transform_offset(xref_mat,cencell_recon,xref_matT)
        call RT_transform_offset(xact_mat,cencell_recon,xact_matT)

        do dir=1,sdim
         moferror=moferror+ &
          mofdata(vofcomp)*((xref_matT(dir)-xact_matT(dir))**2)
        enddo
       enddo  ! im
       moferror=sqrt(moferror)/dx(1)

       if (moferror.gt.max_mof_error) then
        max_mof_error=moferror
       endif

        ! diagnostic_MOF
       do im=1,num_materials
        total_calls(im)=total_calls(im)+mof_calls(1,im)
        total_iterations(im)= &
            total_iterations(im)+mof_iterations(1,im)
        total_errors(im)= &
            total_errors(im)+mof_errors(1,im)
        if (mof_iterations(1,im).gt.max_iterations(im)) then
         max_iterations(im)=mof_iterations(1,im)
        endif
       enddo

      enddo ! ntry

      print *,"sdim= ",sdim
      print *,"MOFITERMAX= ",MOFITERMAX
      print *,"nsamples= ",nsamples
      print *,"shrink_factor=",shrink_factor
      do im=1,num_materials
       print *,"im= ",im
       print *,"total calls= ",total_calls(im)
       if (total_calls(im).gt.zero) then
        avgiter=total_iterations(im)/total_calls(im)
        avgerror=total_errors(im)/total_calls(im)
        print *,"avgiter= ",avgiter
        print *,"avgerror= ",avgerror
       endif
       print *,"max iterations= ",max_iterations(im)
      enddo
      print *,"max_mof_error=",max_mof_error


      deallocate(LS_stencil)

      return
      end subroutine diagnostic_MOF


! vof,ref centroid,order,slope,intercept  x num_materials
      subroutine make_vfrac_sum_ok_base( &
        xsten, &
        nhalf, &
        bfact,dx, &
        tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
        mofdata, &
        sdim)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      integer, INTENT(in) :: tessellate
      real(amrex_real), INTENT(inout) :: mofdata(num_materials*ngeom_recon)

      integer irank
      integer last
      integer last_elastic
      real(amrex_real) uncapt
      real(amrex_real) uncapt_elastic
      real(amrex_real) test_vof

      integer im
      integer dir
      integer vofcomp
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) volcell
      real(amrex_real) cencell(sdim)
      real(amrex_real) xtarget(sdim)
      real(amrex_real) boxlo,boxhi
      integer layer_iter
      real(amrex_real) vfrac_sum(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) vfrac_sum_truncate(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) vfrac_sum_local(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      integer num_materials_local(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

#include "mofdata.H"

      if (nhalf.ge.1) then
       ! do nothing
      else
       print *,"nhalf invalid in make_vfrac_sum_ok_base"
       print *,"nhalf ",nhalf
       stop
      endif

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid: ",bfact
       stop
      endif

      if ((tessellate.eq.TESSELLATE_FLUIDS).or. &
          (tessellate.eq.TESSELLATE_IGNORE_ISRIGID).or. &
          (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC)) then
       !do nothing
      else
       print *,"in make_vfrac_sum_ok_base"
       print *,"expecting tessellate=TESSELLATE_FLUIDS or"
       print *,"tessellate=TESSELLATE_IGNORE_ISRIGID or"
       print *,"tessellate=TESSELLATE_IGNORE_ISELASTIC"
       print *,"tessellate=",tessellate
       stop
      endif

      call init_local_material_vars( & !subroutine make_vfrac_sum_ok_base
       is_rigid_local, &
       is_elastic_local, &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
       tessellate, &
       is_proper_layer_local)

      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials bust: ",num_materials
       stop
      endif
      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon invalid: ",ngeom_recon
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid: ",sdim
       stop
      endif

      call Box_volumeFAST( &
        bfact,dx, &
        xsten,nhalf, &
        volcell, &
        cencell, &
        sdim)

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       vfrac_sum(layer_iter)=zero

       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        if (is_proper_layer(im,layer_iter).eq.1) then

         if ((mofdata(vofcomp).ge.-EPS1).and. &
             (mofdata(vofcomp).le.zero)) then
          mofdata(vofcomp)=zero
          do dir=1,sdim
           mofdata(vofcomp+dir)=zero
          enddo
         else if ((mofdata(vofcomp).le.one+EPS1).and. &
                  (mofdata(vofcomp).ge.one)) then
          mofdata(vofcomp)=one
          do dir=1,sdim
           mofdata(vofcomp+dir)=zero
          enddo
         else if ((mofdata(vofcomp).gt.zero).and. &
                (mofdata(vofcomp).lt.one)) then
          ! do nothing
         else
          print *,"make_vfrac_sum_ok_base"
          print *,"dx=",dx
          print *,"xsten= ",xsten
          print *,"xpoint=",xsten(0,1),xsten(0,2),xsten(0,SDIM)
          print *,"tessellate=",tessellate
          print *,"mofdata(vofcomp) invalid 2"
          print *,"mofdata(vofcomp)=",mofdata(vofcomp)
          print *,"im,num_materials,ngeom_recon,sdim ", &
                  im,num_materials,ngeom_recon,sdim
          print *,"mofdata= ",mofdata
          print *,"put breakpoint here to see the caller"
          stop
         endif

         vfrac_sum(layer_iter)= &
            vfrac_sum(layer_iter)+mofdata(vofcomp)
        endif
       enddo !im=1,num_materials

       if (layer_iter.eq.RIGID_LAYER_INDEX) then

        if ((vfrac_sum(layer_iter).ge.-EPS1).and. &
            (vfrac_sum(layer_iter).le.VOFTOL_LAYER)) then
         vfrac_sum(layer_iter)=zero
         do im=1,num_materials
          vofcomp=(im-1)*ngeom_recon+1
          if (is_proper_layer(im,layer_iter).eq.1) then
           mofdata(vofcomp)=zero
           do dir=1,sdim
            mofdata(vofcomp+dir)=zero
           enddo
          endif
         enddo
        else if ((vfrac_sum(layer_iter).le.one+EPS1).and. &
                 (vfrac_sum(layer_iter).ge.one-VOFTOL_LAYER)) then
         vfrac_sum(layer_iter)=one 
        else if ((vfrac_sum(layer_iter).ge.VOFTOL_LAYER).and. &
                 (vfrac_sum(layer_iter).le.one-VOFTOL_LAYER)) then
         !do nothing
        else
         print *,"vfrac_sum(layer_iter) invalid ",vfrac_sum(layer_iter)
         stop
        endif

       else if (layer_iter.eq.ELASTIC_LAYER_INDEX) then

        if ((vfrac_sum(layer_iter).ge.-EPS1).and. &
            (vfrac_sum(layer_iter).le.VOFTOL_LAYER)) then
         vfrac_sum(layer_iter)=zero
         do im=1,num_materials
          vofcomp=(im-1)*ngeom_recon+1
          if (is_proper_layer(im,layer_iter).eq.1) then
           mofdata(vofcomp)=zero
           do dir=1,sdim
            mofdata(vofcomp+dir)=zero
           enddo
          endif
         enddo
        else if ((vfrac_sum(layer_iter).le.one+EPS1).and. &
                 (vfrac_sum(layer_iter).ge.one-VOFTOL_LAYER)) then
         if (vfrac_sum(RIGID_LAYER_INDEX).eq.zero) then
          vfrac_sum(layer_iter)=one 
         else if (tessellate.eq.TESSELLATE_FLUIDS) then
          vfrac_sum(layer_iter)=one 
         else if ((vfrac_sum(RIGID_LAYER_INDEX).gt.zero).or. &
                  (tessellate.eq.TESSELLATE_IGNORE_ISRIGID).or. &
                  (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC)) then
          !do nothing
         else
          print *,"tessellate invalid ",tessellate
          print *,"or vfrac_sum invalid ",vfrac_sum
          stop
         endif
        else if ((vfrac_sum(layer_iter).ge.VOFTOL_LAYER).and. &
                 (vfrac_sum(layer_iter).le.one-VOFTOL_LAYER)) then
         !do nothing
        else
         print *,"vfrac_sum(layer_iter) invalid ",vfrac_sum(layer_iter)
         stop
        endif

       else if (layer_iter.eq.FLUID_LAYER_INDEX) then

        if ((vfrac_sum(layer_iter).ge.-EPS1).and. &
            (vfrac_sum(layer_iter).le.VOFTOL_LAYER)) then
         vfrac_sum(layer_iter)=zero
         do im=1,num_materials
          vofcomp=(im-1)*ngeom_recon+1
          if (is_proper_layer(im,layer_iter).eq.1) then
           mofdata(vofcomp)=zero
           do dir=1,sdim
            mofdata(vofcomp+dir)=zero
           enddo
          endif
         enddo
        else if ((vfrac_sum(layer_iter).le.one+EPS1).and. &
                 (vfrac_sum(layer_iter).ge.VOFTOL_LAYER)) then

         if ((vfrac_sum(RIGID_LAYER_INDEX).eq.zero).and. &
             (vfrac_sum(ELASTIC_LAYER_INDEX).eq.zero)) then
          vfrac_sum(layer_iter)=one 
         else if (tessellate.eq.TESSELLATE_FLUIDS) then
          vfrac_sum(layer_iter)=one 
         else if ((tessellate.eq.TESSELLATE_IGNORE_ISELASTIC).and. &
                  (vfrac_sum(ELASTIC_LAYER_INDEX).eq.zero)) then
          vfrac_sum(layer_iter)=one 
         else if ((vfrac_sum(RIGID_LAYER_INDEX).gt.zero).or. &
                  (vfrac_sum(ELASTIC_LAYER_INDEX).gt.zero).or. &
                  (tessellate.eq.TESSELLATE_IGNORE_ISRIGID).or. &
                  (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC)) then
          !do nothing
         else
          print *,"tessellate invalid ",tessellate
          print *,"or vfrac_sum invalid ",vfrac_sum
          stop
         endif

        else
         print *,"vfrac_sum(layer_iter) invalid ",vfrac_sum(layer_iter)
         stop
        endif

       else
        print *,"layer_iter invalid"
        stop
       endif

      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1

       if ((mofdata(vofcomp).ge.zero).and. &
           (mofdata(vofcomp).le.VOFTOL_MATERIAL)) then
        mofdata(vofcomp)=zero
        do dir=1,sdim
         mofdata(vofcomp+dir)=zero
        enddo
       else if ((mofdata(vofcomp).le.one).and. &
                (mofdata(vofcomp).ge.one-VOFTOL_MATERIAL)) then
        mofdata(vofcomp)=one
        do dir=1,sdim
         mofdata(vofcomp+dir)=zero
        enddo
       else if ((mofdata(vofcomp).gt.zero).and. &
                (mofdata(vofcomp).lt.one)) then
        ! do nothing
       else
        print *,"make_vfrac_sum_ok_base"
        print *,"dx=",dx
        print *,"xsten= ",xsten
        print *,"xpoint=",xsten(0,1),xsten(0,2),xsten(0,SDIM)
        print *,"tessellate=",tessellate
        print *,"mofdata(vofcomp) invalid 2"
        print *,"mofdata(vofcomp)=",mofdata(vofcomp)
        print *,"im,num_materials,ngeom_recon,sdim ", &
                im,num_materials,ngeom_recon,sdim
        print *,"mofdata= ",mofdata
        print *,"put breakpoint here to see the caller"
        stop
       endif
      enddo ! im=1,num_materials

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       vfrac_sum_truncate(layer_iter)=zero
       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        if (is_proper_layer(im,layer_iter).eq.1) then
         vfrac_sum_truncate(layer_iter)= &
            vfrac_sum_truncate(layer_iter)+mofdata(vofcomp)
        endif
       enddo !im=1,num_materials
       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        if (vfrac_sum_truncate(layer_iter).gt.zero) then
         mofdata(vofcomp)=mofdata(vofcomp)*vfrac_sum(layer_iter)/ &
               vfrac_sum_truncate(layer_iter)
        else if (vfrac_sum_truncate(layer_iter).eq.zero) then
         !do nothing
        else
         print *,"vfrac_sum_truncate invalid ",vfrac_sum_truncate
         stop
        endif
       enddo !im=1,num_materials
      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       vfrac_sum_local(layer_iter)=zero
       num_materials_local(layer_iter)=0
       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        if (is_proper_layer_local(im,layer_iter).eq.1) then
         num_materials_local(layer_iter)=num_materials_local(layer_iter)+1
         vfrac_sum_local(layer_iter)= &
            vfrac_sum_local(layer_iter)+mofdata(vofcomp)
        endif
       enddo !im=1,num_materials
      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       ! for all of the cases,
       !  TESSELLATE_IGNORE_ISRIGID
       !  TESSELLATE_IGNORE_ISELASTIC
       !  TESSELLATE_FLUIDS
       ! the "fluids" should tessellate the domain.
      if (vfrac_sum_local(FLUID_LAYER_INDEX).gt.zero) then
       ! do nothing
      else if (vfrac_sum_local(FLUID_LAYER_INDEX).le.zero) then
       print *,"vacuum bust in make_vfrac_sum_ok_base"
       print *,"put breakpoint here to see the caller"
       print *,"break MOF.F90:16201 (compile -g option)"
       print *,"num_materials= ",num_materials
       print *,"sdim= ",sdim
       print *,"vfrac_sum_local= ",vfrac_sum_local
       print *,"xsten= ",xsten
       print *,"xpoint=",xsten(0,1),xsten(0,2),xsten(0,SDIM)
       print *,"nhalf= ",nhalf
       print *,"bfact=",bfact
       print *,"dx=",dx
       print *,"tessellate=",tessellate
       print *,"mofdata=",mofdata
       print *,"is_rigid_local=",is_rigid_local
       print *,"is_elastic_local=",is_elastic_local
       stop
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1
       if (is_rigid_local(im).eq.1) then
        if (vfrac_sum_local(RIGID_LAYER_INDEX).gt.one) then
         mofdata(vofcomp)=mofdata(vofcomp)/vfrac_sum_local(RIGID_LAYER_INDEX)
        else if ((vfrac_sum_local(RIGID_LAYER_INDEX).ge.zero).and. &
                 (vfrac_sum_local(RIGID_LAYER_INDEX).le.one)) then
         ! do nothing
        else
         print *,"vfrac_sum_local invalid: ",vfrac_sum_local
         stop
        endif
       else if (is_rigid_local(im).eq.0) then
        !do nothing
       else
        print *,"is_rigid_local invalid MOF.F90: ",im,is_rigid_local(im)
        stop
       endif
      enddo  ! im=1..num_materials

      if (tessellate.eq.TESSELLATE_IGNORE_ISRIGID) then

       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        if ((is_rigid_local(im).eq.1).or. &
            (is_elastic_local(im).eq.1)) then
         print *,"expecting is_rigid_local(im)==0"
         print *,"and expecting is_elastic_local(im)==0"
         print *,"tessellate=",tessellate
         stop
        else if ((is_rigid_local(im).eq.0).and. &
                 (is_elastic_local(im).eq.0)) then
         mofdata(vofcomp)=mofdata(vofcomp)/vfrac_sum_local(FLUID_LAYER_INDEX)
        else
         print *,"is_rigid_local invalid MOF.F90: ",is_rigid_local
         print *,"or is_elastic_local invalid MOF.F90: ",is_elastic_local
         stop
        endif
       enddo  ! im=1..num_materials

       ! fluids tessellate, is_elastic and is_rigid embedded.
      else if ((tessellate.eq.TESSELLATE_FLUIDS).or. &
               (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC)) then

       if ((nonzero_ranks.eq. &
            num_materials_local(FLUID_LAYER_INDEX)+ &
            num_materials_local(ELASTIC_LAYER_INDEX)).and. &
           (num_materials_local(FLUID_LAYER_INDEX)+ &
            num_materials_local(ELASTIC_LAYER_INDEX).ge.1)) then
        !do nothing
       else
        print *,"nonzero_ranks invalid ",nonzero_ranks
        print *,"num_materials_local ",num_materials_local
        stop
       endif

       last=0
       last_elastic=0
       uncapt=one
       uncapt_elastic=one

       do irank=1,nonzero_ranks

        im=rank_algorithm(irank)

        if (renormalize_order_algorithm(im).eq.irank) then
         !do nothing
        else
         print *,"irank inconsistent"
         print *,"irank=",irank
         print *,"im=",im
         print *,"renormalize_order_algorithm(im)=", &
             renormalize_order_algorithm(im)
         stop
        endif

        vofcomp=(im-1)*ngeom_recon+1

        if ((im.ge.1).and.(im.le.num_materials)) then

         if (is_rigid_local(im).eq.0) then
          !do nothing
         else
          print *,"expecting is_rigid_local(im).eq.0 if ranked"
          print *,"tessellate=",tessellate
          print *,"im=",im
          print *,"irank=",irank
          stop
         endif 

        else
         print *,"im invalid ",im
         stop
        endif

        if (is_elastic_local(im).eq.0) then

          !uncapt=1 initially
         if (uncapt.gt.zero) then
          if ((im.ge.1).and.(im.le.num_materials)) then
           if (is_rigid_local(im).eq.0) then
            !do nothing
           else
            print *,"expecting is_rigid_local(im).eq.0"
            print *,"im=",im
            print *,"irank=",irank
            stop
           endif 
           test_vof=mofdata(vofcomp)
           if (test_vof.gt.zero) then
            last=im
           endif
           uncapt=uncapt-test_vof
          else
           print *,"im invalid ",im
           stop
          endif
         else if (uncapt.le.zero) then
          mofdata(vofcomp)=zero
          do dir=1,sdim
           mofdata(vofcomp+dir)=zero
          enddo
         else
          print *,"uncapt=NaN: ",uncapt
          stop
         endif

        else if (is_elastic_local(im).eq.1) then
       
          !uncapt_elastic=1 initially
         if (uncapt_elastic.gt.zero) then
          if ((im.ge.1).and.(im.le.num_materials)) then
           if (is_rigid_local(im).eq.0) then
            !do nothing
           else
            print *,"expecting is_rigid_local(im).eq.0"
            print *,"im=",im
            print *,"irank=",irank
            stop
           endif 
           test_vof=mofdata(vofcomp)
           if (test_vof.gt.zero) then
            last_elastic=im
           endif
           uncapt_elastic=uncapt_elastic-test_vof
          else
           print *,"im invalid ",im
           stop
          endif
         else if (uncapt_elastic.le.zero) then
          mofdata(vofcomp)=zero
          do dir=1,sdim
           mofdata(vofcomp+dir)=zero
          enddo
         else
          print *,"uncapt_elastic=NaN: ",uncapt_elastic
          stop
         endif

        else
         print *,"is_elastic_local(im) invalid: ",im,is_elastic_local(im)
         stop
        endif

       enddo !irank=1,nonzero_ranks

       if (last.eq.0) then
        print *,"all volume vanished make_vfrac_sum_ok_base"
        stop
       else if ((last.ge.1).and.(last.le.num_materials)) then
        im=last
        vofcomp=(im-1)*ngeom_recon+1
        mofdata(vofcomp)=mofdata(vofcomp)+uncapt

        if (uncapt_elastic.lt.zero) then
         if ((last_elastic.ge.1).and.(last_elastic.le.num_materials)) then
          im=last_elastic
          vofcomp=(im-1)*ngeom_recon+1
          mofdata(vofcomp)=mofdata(vofcomp)+uncapt_elastic
         else
          print *,"expecting 1<=last_elastic<=num_materials: ",last_elastic
          stop
         endif
        else if ((uncapt_elastic.ge.zero).and.(uncapt_elastic.le.one)) then
         !do nothing
        else
         print *,"uncapt_elastic invalid ",uncapt_elastic
         stop
        endif

       else 
        print *,"last invalid: ",last
        stop
       endif
       
      else
       print *,"tessellate invalid make_vfrac_sum_ok_base: ",tessellate
       stop
      endif

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1
       if ((mofdata(vofcomp).eq.zero).or. &
           (mofdata(vofcomp).eq.one)) then
        do dir=1,sdim
         mofdata(vofcomp+dir)=zero
        enddo
       else if ((mofdata(vofcomp).gt.zero).and. &
                (mofdata(vofcomp).lt.one)) then

        do dir=1,sdim
         xtarget(dir)=mofdata(vofcomp+dir)+cencell(dir)
        enddo

        do dir=1,sdim
         boxlo=xsten(-1,dir)
         boxhi=xsten(1,dir)
         if (boxhi.gt.boxlo) then
          ! do nothing
         else
          print *,"boxhi>boxlo condition failed: ",boxlo,boxhi
          stop
         endif

         if (xtarget(dir).le.boxlo) then
          mofdata(vofcomp+dir)=boxlo-cencell(dir)+EPS_13_5*dx(dir)
         else if (xtarget(dir).ge.boxhi) then
          mofdata(vofcomp+dir)=boxhi-cencell(dir)-EPS_13_5*dx(dir)
         else if ((xtarget(dir).gt.boxlo).and. &
                  (xtarget(dir).lt.boxhi)) then
          ! do nothing
         else
          print *,"xtarget(dir) invalid: ",xtarget(dir)
          stop
         endif
        enddo !dir=1..sdim

       else
        print *,"mofdata(vofcomp) invalid"
        print *,"vofcomp=",vofcomp
        print *,"mofdata(vofcomp)=",mofdata(vofcomp)
        stop
       endif
      enddo ! im=1..num_materials

      return
      end subroutine make_vfrac_sum_ok_base


! vof,ref centroid,order,slope,intercept  x num_materials
      subroutine make_vfrac_sum_ok_copy( &
        xsten,nhalf, &
        bfact,dx, &
        tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
        mofdata,mofdatavalid, &
        sdim)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)

      integer, INTENT(in) :: tessellate
      real(amrex_real), INTENT(in) :: mofdata(num_materials*ngeom_recon)
      real(amrex_real), INTENT(out) :: mofdatavalid(num_materials*ngeom_recon)
      integer im
      integer dir
      integer vofcomp
      real(amrex_real) vfrac_sum_local(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) vof_test
      integer layer_iter
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      if (nhalf.ge.1) then
       ! do nothing
      else
       print *,"nhalf invalid in make_vfrac_sum_ok_copy"
       print *,"nhalf ",nhalf
       stop
      endif

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid"
       stop
      endif

      call init_local_material_vars( & !subroutine make_vfrac_sum_ok_copy
       is_rigid_local, &
       is_elastic_local, &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
       tessellate, &
       is_proper_layer_local)

      if ((num_materials.ge.1).and.(num_materials.le.MAX_NUM_MATERIALS)) then
       ! do nothing
      else
       print *,"num_materials bust: ",num_materials
       stop
      endif
      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon invalid"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       vfrac_sum_local(layer_iter)=zero

       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        vof_test=mofdata(vofcomp)
        if ((vof_test.ge.-EPS1).and. &
            (vof_test.le.VOFTOL_MATERIAL)) then
         vof_test=zero
        else if ((vof_test.le.one+EPS1).and. &
                 (vof_test.ge.one-VOFTOL_MATERIAL)) then
         vof_test=one
        else if ((vof_test.gt.zero).and. &
                 (vof_test.lt.one)) then
         ! do nothing
        else
         print *,"vof_test invalid: ",vof_test
         stop
        endif
        if (is_proper_layer_local(im,layer_iter).eq.1) then
         vfrac_sum_local(layer_iter)= &
            vfrac_sum_local(layer_iter)+mofdata(vofcomp)
        endif
        do dir=1,ngeom_recon
         mofdatavalid(dir+vofcomp-1)=mofdata(dir+vofcomp-1)
        enddo
       enddo ! im=1..num_materials
      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      if (vfrac_sum_local(FLUID_LAYER_INDEX).gt.zero) then
       !do nothing
      else
       print *,"vacuum bust in make_vfrac_sum_ok_copy"
       print *,"put breakpoint here to see the caller"
       print *,"num_materials= ",num_materials
       print *,"sdim= ",sdim
       print *,"vfrac_sum_local ",vfrac_sum_local
       print *,"tessellate=",tessellate
       print *,"dx=",dx
       print *,"xsten=",xsten
       print *,"xpoint=",xsten(0,1),xsten(0,2),xsten(0,SDIM)
       print *,"mofdata=",mofdata
       print *,"is_rigid_local=",is_rigid_local
       print *,"is_elastic_local=",is_elastic_local
       stop
      endif

      call make_vfrac_sum_ok_base( &
       xsten, &
       nhalf, &
       bfact,dx, &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
       mofdatavalid, &
       sdim)

      return
      end subroutine make_vfrac_sum_ok_copy


      subroutine project_slopes_to_face( &
       bfact,dx,xsten0,nhalf0, &
       mofdata,mofdataproject, &
       sdim,dir_side,side)

      use probcommon_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in)  :: sdim,dir_side,side
      integer, INTENT(in)  :: bfact,nhalf0
      real(amrex_real), INTENT(in)  :: mofdata(num_materials*ngeom_recon)
      real(amrex_real), INTENT(out) :: mofdataproject(num_materials*ngeom_recon)
      real(amrex_real), INTENT(in)  :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in)  :: dx(sdim)
      integer vofcomp,im
      integer dir_local
      real(amrex_real) x0_face(sdim)
      real(amrex_real) slope(sdim)
      real(amrex_real) slope_project(sdim)
      real(amrex_real) intercept_project
      real(amrex_real) mag

      if (ngeom_recon.eq.2*sdim+3) then
       ! do nothing
      else
       print *,"ngeom_recon.ne.2*sdim+3"
       stop
      endif

      if ((dir_side.ge.1).and.(dir_side.le.sdim)) then
       ! do nothing
      else
       print *,"dir_side invalid project slopes to face"
       stop
      endif
      if ((side.eq.1).or.(side.eq.2)) then
       ! do nothing
      else
       print *,"side invalid"
       stop
      endif

      if (nhalf0.ge.1) then
       ! do nothing
      else
       print *,"nhalf0 invalid project_slopes_to_face"
       stop
      endif
      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid135"
       stop
      endif
      if ((sdim.eq.3).or.(sdim.eq.2)) then
       ! do nothing
      else
       print *,"sdim invalid project_slopes_to_face"
       stop
      endif
      if ((num_materials.ge.1).and.(num_materials.le.MAX_NUM_MATERIALS)) then
       ! do nothing
      else
       print *,"num_materials invalid project_slopes_to_face"
       stop
      endif

       ! F,X,order,slope,intercept
      do im=1,num_materials

       vofcomp=(im-1)*ngeom_recon+1

       do dir_local=1,ngeom_recon
        mofdataproject(vofcomp+dir_local-1)=mofdata(vofcomp+dir_local-1)
       enddo

       do dir_local=1,sdim
        slope(dir_local)=mofdata(vofcomp+sdim+1+dir_local)
        slope_project(dir_local)=slope(dir_local)
       enddo
       intercept_project=mofdata(vofcomp+2*sdim+2)
        ! Pr(n) = (I-nf nf^T)n  nf=normal to face
       slope_project(dir_side)=zero  
       mag=zero
       do dir_local=1,sdim
        mag=mag+slope_project(dir_local)*slope_project(dir_local)
       enddo
       mag=sqrt(mag)
       if (mag.gt.zero) then
        do dir_local=1,sdim
         slope_project(dir_local)=slope_project(dir_local)/mag
         x0_face(dir_local)=xsten0(0,dir_local)
        enddo 
        if (side.eq.1) then
         x0_face(dir_side)=xsten0(-1,dir_side)
        else if (side.eq.2) then
         x0_face(dir_side)=xsten0(1,dir_side)
        else
         print *,"side invalid"
         stop
        endif
        intercept_project=(intercept_project+ &
         slope(dir_side)*(x0_face(dir_side)-xsten0(0,dir_side)))/mag

         ! Pr(n)=(I-nf nf^T)n
         ! ntilde=Pr(n)/||Pr(n)||
         ! intercept_project=(n dot (xf-x0) + intercept_old)/||Pr(n)||
         ! projected representation:
         ! ntilde dot (x-x0) + intercept_project =
         ! (Pr(n) dot (x-x0) + n dot (xf-x0) + intercept_old)/||Pr(n)||
         ! Since Pr(n) has no nf component, Pr(n) dot (x-x0)=
         ! Pr(n) dot (x-xf).
         ! If x is on the face, then Pr(n) dot (x-xf)=n dot (x-xf)
         ! so that
         ! ntilde dot (x-x0)+ intercept_project=
         ! (n dot (x-xf) + n dot (xf-x0) + intercept_old)/||Pr(n)||=
         ! (n dot (x-x0) + intercept_old)/||Pr(n)||.

        do dir_local=1,sdim
         mofdataproject(vofcomp+sdim+1+dir_local)=slope_project(dir_local)
        enddo
        mofdataproject(vofcomp+2*sdim+2)=intercept_project
       else if (mag.eq.zero) then
        ! do not modify the slope or intercept
       else
        print *,"mag invalid in project_slopes_to_face"
        print *,"mag=",mag
        print *,"dir_side=",dir_side
        print *,"intercept_project= ",intercept_project
        print *,"slope_project=",slope_project
        stop
       endif

      enddo ! im=1..num_materials

      return
      end subroutine project_slopes_to_face

      subroutine init_local_material_vars( &
         is_rigid_local, &
         is_elastic_local, &
         tessellate_source, &
         tessellate_dest, &
         is_proper_layer_local)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: tessellate_source
      integer, INTENT(in) :: tessellate_dest
      integer, intent(out) :: is_rigid_local(num_materials)
      integer, intent(out) :: is_elastic_local(num_materials)
      integer, intent(out) :: is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      integer :: im

      do im=1,num_materials
       is_rigid_local(im)=is_rigid(im)
       is_elastic_local(im)=is_elastic(im)
       if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then
        is_rigid_local(im)=0
        is_elastic_local(im)=0
        if (tessellate_source.eq.tessellate_dest) then
         !do nothing
        else
         print *,"expecting tessellate_source==tessellate_dest"
         print *,"tessellate_source=",tessellate_source
         print *,"tessellate_dest=",tessellate_dest
         stop
        endif
       else if (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC) then
        is_elastic_local(im)=0
        if ((tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC).or. &
            (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
         !do nothing
        else
         print *,"tessellate_source or tessellate_dest invalid"
         print *,"tessellate_source=",tessellate_source
         print *,"tessellate_dest=",tessellate_dest
         stop
        endif
       else if (tessellate_source.eq.TESSELLATE_FLUIDS) then
        if ((tessellate_dest.eq.TESSELLATE_ALL).or. &
            (tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
            (tessellate_dest.eq.TESSELLATE_ALL_RASTER).or. &
            (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC)) then
         ! do nothing
        else
         print *,"tessellate_dest invalid: ",tessellate_dest
         stop
        endif
       else
        print *,"tessellate_source invalid10: ",tessellate_source
        stop
       endif
      enddo ! im=1..num_materials

      do im=1,num_materials
       is_proper_layer_local(im,RIGID_LAYER_INDEX)= &
               is_rigid_local(im)
       is_proper_layer_local(im,ELASTIC_LAYER_INDEX)= &
               is_elastic_local(im)
       if ((is_rigid_local(im).eq.0).and. &
           (is_elastic_local(im).eq.0)) then
        is_proper_layer_local(im,FLUID_LAYER_INDEX)=1
       else if ((is_rigid_local(im).eq.1).or. &
                (is_elastic_local(im).eq.1)) then
        is_proper_layer_local(im,FLUID_LAYER_INDEX)=0
       else
        print *,"is_rigid_local or is_elastic_local invalid: ", &
          is_rigid_local,is_elastic_local
        stop
       endif
      enddo !im=1,num_materials

      end subroutine init_local_material_vars

      subroutine tiny_fraction_case( &
        mofdatasave, &
        uncaptured_volume_START, &
        uncaptured_centroid, &
        multi_volume, &
        multi_cen, &
        is_rigid_local, &
        is_elastic_local, &
        tessellate_source, &
        tessellate_dest, &
        vfrac_sum_local, &
        sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      real(amrex_real), intent(in) :: mofdatasave(num_materials*(2*sdim+3))
      real(amrex_real), intent(in) :: uncaptured_volume_START
      real(amrex_real), intent(in) :: &
          uncaptured_centroid(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX,sdim)
      real(amrex_real), INTENT(inout) :: multi_volume(num_materials)
      real(amrex_real), INTENT(inout) :: multi_cen(sdim,num_materials)
      integer, intent(in) :: is_rigid_local(num_materials)
      integer, intent(in) :: is_elastic_local(num_materials)
      integer, INTENT(in) :: tessellate_source
      integer, INTENT(in) :: tessellate_dest
      real(amrex_real), intent(in) :: vfrac_sum_local( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      integer, INTENT(in) :: sdim
      integer im,vofcomp,dir

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1

       multi_volume(im)=uncaptured_volume_START*mofdatasave(vofcomp)

       if (is_rigid_local(im).eq.1) then
        !do nothing
       else if ((is_rigid_local(im).eq.0).and. &
                (is_elastic_local(im).eq.1)) then

        if (tessellate_source.eq.TESSELLATE_FLUIDS) then

         if (tessellate_dest.eq.TESSELLATE_FLUIDS) then
          !do nothing
         else if (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC) then
          !do nothing
         else if (tessellate_dest.eq.TESSELLATE_ALL) then
          multi_volume(im)=multi_volume(im)* &
            abs(one-vfrac_sum_local(RIGID_LAYER_INDEX))
         else if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then
          !do nothing
         else
          print *,"tessellate_dest invalid: ",tessellate_dest
          stop
         endif

        else if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then

         if (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID) then
          !do nothing
         else
          print *,"tessellate_dest invalid: ",tessellate_dest
          stop
         endif

        else
         print *,"tessellate_source invalid ",tessellate_source
         stop
        endif

       else if ((is_rigid_local(im).eq.0).and. &
                (is_elastic_local(im).eq.0)) then

        if ((tessellate_source.eq.TESSELLATE_FLUIDS).and. &
            (tessellate_dest.eq.TESSELLATE_FLUIDS)) then
         !do nothing
        else if ((tessellate_source.eq.TESSELLATE_FLUIDS).and. &
                 (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC)) then
         multi_volume(im)=multi_volume(im)* &
           abs(one-vfrac_sum_local(ELASTIC_LAYER_INDEX))
        else if ((tessellate_source.eq.TESSELLATE_FLUIDS).and. &
                 (tessellate_dest.eq.TESSELLATE_ALL)) then
         multi_volume(im)=multi_volume(im)* &
            abs(one-vfrac_sum_local(RIGID_LAYER_INDEX))* &
            abs(one-vfrac_sum_local(ELASTIC_LAYER_INDEX))
        else if ((tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC).and. &
                 (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
         !do nothing
        else if ((tessellate_source.eq.TESSELLATE_FLUIDS).and. &
                 (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
         multi_volume(im)=multi_volume(im)* &
            abs(one-vfrac_sum_local(ELASTIC_LAYER_INDEX))
        else if ((tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID).and. &
                 (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID)) then 
         !do nothing
        else if ((tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC).and. &
                 (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC)) then 
         !do nothing
        else
         print *,"tessellate_source invalid: ",tessellate_source
         print *,"or tessellate_dest invalid: ",tessellate_dest
         stop
        endif
                
       else
        print *,"is_rigid_local invalid: ",is_rigid_local
        print *,"or is_elastic_local invalid: ",is_elastic_local
        stop
       endif

       do dir=1,sdim
        multi_cen(dir,im)=uncaptured_centroid(FLUID_LAYER_INDEX,dir)
       enddo
      enddo ! im=1..num_materials

      return
      end subroutine tiny_fraction_case
       

      subroutine init_layer_flag( &
         at_least_one, &
         vfrac_sum, &
         vfrac_sum_local, &
         layer_flag, & !intent(out)
         layer_iter, & !intent(in)
         tessellate_source, &
         tessellate_dest)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, intent(in) :: at_least_one

      real(amrex_real), intent(in) :: vfrac_sum( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real), intent(in) :: vfrac_sum_local( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      integer, intent(out) :: layer_flag
      integer, intent(in) :: layer_iter
      integer, intent(in) :: tessellate_source
      integer, intent(in) :: tessellate_dest

      if (at_least_one.eq.0) then
       !do nothing
      else if ((at_least_one.ge.1).and. &
               (at_least_one.le.num_materials)) then
       !do nothing
      else
       print *,"at_least_one invalid ",at_least_one
       stop
      endif

      if (layer_iter.eq.RIGID_LAYER_INDEX) then

       if (tessellate_source.eq.TESSELLATE_FLUIDS) then
        if ((tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
            (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC).or. &
            (tessellate_dest.eq.TESSELLATE_ALL)) then
         layer_flag=RIGID_LAYER
        else if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then
         layer_flag=RIGID_LAYER
         if (at_least_one.eq.0) then
          !do nothing
         else if ((vfrac_sum(layer_iter).eq.zero).and. &
                  (vfrac_sum_local(layer_iter).eq.zero)) then
          !do nothing
         else
          print *,"break MOF.F90:16992"
          print *,"expecting at_least_one=0 ",at_least_one
          print *,"or expecting vfrac_sum=vfrac_sum_local=0 ", &
           vfrac_sum(layer_iter),vfrac_sum_local(layer_iter)
          print *,"tessellate_source=",tessellate_source
          print *,"tessellate_dest=",tessellate_dest
          print *,"layer_iter ",layer_iter
          stop
         endif
        else
         print *,"tessellate_dest invalid in init_layer_flag ",tessellate_dest
         print *,"tessellate_source ",tessellate_source
         stop
        endif
       else if (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC) then
        if (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC) then
         layer_flag=RIGID_LAYER
        else if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then
         layer_flag=RIGID_LAYER
         if (at_least_one.eq.0) then
          !do nothing
         else if ((vfrac_sum(layer_iter).eq.zero).and. &
                  (vfrac_sum_local(layer_iter).eq.zero)) then
          !do nothing
         else
          print *,"expecting at_least_one=0 ",at_least_one
          print *,"or expecting vfrac_sum=vfrac_sum_local=0 ", &
           vfrac_sum(layer_iter),vfrac_sum_local(layer_iter)
          print *,"tessellate_source=",tessellate_source
          print *,"tessellate_dest=",tessellate_dest
          print *,"layer_iter ",layer_iter
          stop
         endif
        else
         print *,"tessellate_dest invalid ",tessellate_dest
         print *,"tessellate_source ",tessellate_source
         print *,"layer_iter ",layer_iter
         stop
        endif
       else if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then
        layer_flag=RIGID_LAYER
        if (at_least_one.eq.0) then
         !do nothing
        else
         print *,"expecting at_least_one=0 ",at_least_one
         print *,"tessellate_source=",tessellate_source
         print *,"tessellate_dest=",tessellate_dest
         print *,"layer_iter ",layer_iter
         stop
        endif
       else
        print *,"tessellate_source invalid in init_layer_flag ", &
          tessellate_source
        print *,"tessellate_dest ",tessellate_dest
        stop
       endif

      else if (layer_iter.eq.ELASTIC_LAYER_INDEX) then
         
       if (tessellate_source.eq.TESSELLATE_FLUIDS) then
        if ((tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
            (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC).or. &
            (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
         layer_flag=ELASTIC_LAYER
        else if (tessellate_dest.eq.TESSELLATE_ALL) then
         layer_flag=ELASTIC_RIGID_LAYER
        else
         print *,"tessellate_dest invalid in init_layer_flag ",tessellate_dest
         print *,"tessellate_source ",tessellate_source
         stop
        endif
       else if ((tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID).or. &
                (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC)) then
        layer_flag=ELASTIC_LAYER
        if (at_least_one.eq.0) then
         !do nothing
        else
         print *,"expecting at_least_one=0 ",at_least_one
         print *,"tessellate_source=",tessellate_source
         print *,"tessellate_dest=",tessellate_dest
         print *,"layer_iter ",layer_iter
         stop
        endif
       else
        print *,"tessellate_source invalid ",tessellate_source
        print *,"tessellate_dest ",tessellate_dest
        print *,"layer_iter ",layer_iter
        stop
       endif

      else if (layer_iter.eq.FLUID_LAYER_INDEX) then

       if (tessellate_source.eq.TESSELLATE_FLUIDS) then

        if (tessellate_dest.eq.TESSELLATE_FLUIDS) then
         layer_flag=FLUIDS_LAYER
        else if ((tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC).or. &
                 (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
         if (vfrac_sum_local(ELASTIC_LAYER_INDEX).eq. &
             vfrac_sum(ELASTIC_LAYER_INDEX)) then
          !do nothing
         else 
          print *,"expecting vfrac_sum=vfrac_sum_local: ", &
            vfrac_sum,vfrac_sum_local
          stop
         endif
         if (vfrac_sum_local(ELASTIC_LAYER_INDEX).gt.zero) then
          layer_flag=FLUIDS_ELASTIC_LAYER
         else if (vfrac_sum_local(ELASTIC_LAYER_INDEX).eq.zero) then
          layer_flag=FLUIDS_LAYER
         else
          print *,"vfrac_sum_local(ELASTIC_LAYER_INDEX) invalid"
          stop
         endif
        else if (tessellate_dest.eq.TESSELLATE_ALL) then

         if (vfrac_sum_local(RIGID_LAYER_INDEX).eq. &
             vfrac_sum(RIGID_LAYER_INDEX)) then
          !do nothing
         else 
          print *,"expecting vfrac_sum=vfrac_sum_local: ", &
            vfrac_sum,vfrac_sum_local
          stop
         endif
         if (vfrac_sum_local(ELASTIC_LAYER_INDEX).eq. &
             vfrac_sum(ELASTIC_LAYER_INDEX)) then
          !do nothing
         else 
          print *,"expecting vfrac_sum=vfrac_sum_local: ", &
            vfrac_sum,vfrac_sum_local
          stop
         endif
         if ((vfrac_sum_local(RIGID_LAYER_INDEX).gt.zero).and. &
             (vfrac_sum_local(ELASTIC_LAYER_INDEX).gt.zero)) then
          layer_flag=FLUIDS_ELASTIC_RIGID_LAYER
         else if ((vfrac_sum_local(RIGID_LAYER_INDEX).eq.zero).and. &
                  (vfrac_sum_local(ELASTIC_LAYER_INDEX).gt.zero)) then
          layer_flag=FLUIDS_ELASTIC_LAYER
         else if ((vfrac_sum_local(RIGID_LAYER_INDEX).gt.zero).and. &
                  (vfrac_sum_local(ELASTIC_LAYER_INDEX).eq.zero)) then
          layer_flag=FLUIDS_RIGID_LAYER
         else if ((vfrac_sum_local(RIGID_LAYER_INDEX).eq.zero).and. &
                  (vfrac_sum_local(ELASTIC_LAYER_INDEX).eq.zero)) then
          layer_flag=FLUIDS_LAYER
         else
          print *,"vfrac_sum_local invalid ",vfrac_sum_local
          stop
         endif
        else
         print *,"tessellate_dest invalid ",tessellate_dest
         stop
        endif

       else if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then

        if (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID) then
         if (vfrac_sum_local(RIGID_LAYER_INDEX).eq.zero) then
          !do nothing
         else
          print *,"expecting vfrac_sum_local=0 ",vfrac_sum_local
          stop
         endif
         if (vfrac_sum_local(ELASTIC_LAYER_INDEX).eq.zero) then
          !do nothing
         else
          print *,"expecting vfrac_sum_local=0 ",vfrac_sum_local
          stop
         endif
         if ((vfrac_sum(RIGID_LAYER_INDEX).gt.zero).and. &
             (vfrac_sum(ELASTIC_LAYER_INDEX).gt.zero)) then
          layer_flag=FLUIDS_ELASTIC_RIGID_LAYER
         else if ((vfrac_sum(RIGID_LAYER_INDEX).eq.zero).and. &
                  (vfrac_sum(ELASTIC_LAYER_INDEX).gt.zero)) then
          layer_flag=FLUIDS_ELASTIC_LAYER
         else if ((vfrac_sum(RIGID_LAYER_INDEX).gt.zero).and. &
                  (vfrac_sum(ELASTIC_LAYER_INDEX).eq.zero)) then
          layer_flag=FLUIDS_RIGID_LAYER
         else if ((vfrac_sum(RIGID_LAYER_INDEX).eq.zero).and. &
                  (vfrac_sum(ELASTIC_LAYER_INDEX).eq.zero)) then
          layer_flag=FLUIDS_LAYER
         else
          print *,"vfrac_sum invalid ",vfrac_sum
          stop
         endif
        else
         print *,"tessellate_dest invalid ",tessellate_dest
         stop
        endif

       else if (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC) then

        if (vfrac_sum_local(ELASTIC_LAYER_INDEX).eq.zero) then
         !do nothing
        else
         print *,"expecting vfrac_sum_local=0 ",vfrac_sum_local
         stop
        endif

        if ((tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC).or. &
            (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
         if (vfrac_sum(ELASTIC_LAYER_INDEX).gt.zero) then
          layer_flag=FLUIDS_ELASTIC_LAYER
         else if (vfrac_sum(ELASTIC_LAYER_INDEX).eq.zero) then
          layer_flag=FLUIDS_LAYER
         else
          print *,"vfrac_sum invalid ",vfrac_sum
          stop
         endif
        else
         print *,"tessellate_dest invalid ",tessellate_dest
         stop
        endif

       else
        print *,"tessellate_source invalid ",tessellate_source
        stop
       endif

      else
       print *,"layer_iter invalid ",layer_iter
       stop
      endif

      return
      end subroutine init_layer_flag

      subroutine derive_critical_material( &
        mofdatalocal, &
        mofdatasave, &
        material_used, &
        is_masked, &
        critical_material, &
        sdim)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module
      IMPLICIT NONE
 
      integer, intent(in) :: sdim
      integer, intent(out) :: critical_material
      real(amrex_real), intent(in) :: mofdatasave(num_materials*(2*sdim+3))
      real(amrex_real), intent(in) :: mofdatalocal(num_materials*(2*sdim+3))
      integer, intent(in) :: material_used(num_materials)
      integer, intent(in) :: is_masked(num_materials)
      integer im,vofcomp,testflag,testflag_save
      integer vofcomp_compare,testflag_compare
      integer old_material_rank,current_material_rank

      critical_material=0
      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1

       if (is_masked(im).eq.0) then
        testflag=NINT(mofdatalocal(vofcomp+sdim+1)) !"progress"
        testflag_save=NINT(mofdatasave(vofcomp+sdim+1)) !original
        if ((testflag_save.ge.1).and. & !old_flag=processed
            (testflag_save.le.num_materials)) then
   
         if ((testflag.eq.0).and. &
             (material_used(im).eq.0)) then

          if (critical_material.eq.0) then
           critical_material=im
          else if ((critical_material.ge.1).and. &
                   (critical_material.le.num_materials)) then

           vofcomp_compare=(critical_material-1)*ngeom_recon+1
           testflag_compare=NINT(mofdatasave(vofcomp_compare+sdim+1))

           old_material_rank=get_material_rank(critical_material, &
             testflag_compare)
           current_material_rank=get_material_rank(im, &
             testflag_save)
           if (current_material_rank.lt.old_material_rank) then
            critical_material=im
           else if (current_material_rank.gt.old_material_rank) then
            !do nothing
           else
            print *,"expecting current_material_rank<>old_material_rank"
            stop
           endif
          else
           print *,"critical_material invalid ",critical_material
           stop
          endif

         else if ((testflag.ge.1).and. &
                  (testflag.le.num_materials).and. &
                  (material_used(im).ge.1).and. &
                  (material_used(im).le.num_materials)) then
          !do nothing
         else
          print *,"testflag invalid ",testflag
          print *,"or material_used invalid: ",material_used
          stop
         endif
 
        else if (testflag_save.eq.0) then ! old flag=unprocessed
         ! do nothing
        else
         print *,"testflag_save invalid? ",testflag_save
         stop         
        endif 
       else if (is_masked(im).eq.1) then
        ! do nothing
       else
        print *,"is_masked invalid (derive_critical_material) ", &
          im,is_masked(im)
        stop
       endif
      enddo ! im=1..num_materials

      return
      end subroutine derive_critical_material

      subroutine advance_uncaptured_vars( &
         at_least_one, &
         tessellate_source, &
         tessellate_dest, &
         layer_iter, &
         uncaptured_volume, &
         uncaptured_centroid, &
         uncaptured_volume_fraction, &
         sdim)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, intent(in) :: at_least_one
      integer, intent(in) :: sdim
      integer, intent(in) :: tessellate_source
      integer, intent(in) :: tessellate_dest
      integer, intent(in) :: layer_iter
      real(amrex_real), intent(inout) :: &
          uncaptured_volume(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real), intent(inout) :: &
          uncaptured_centroid(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX,sdim)
      real(amrex_real), intent(in) :: &
          uncaptured_volume_fraction( &
             RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      integer dir

      if (layer_iter.eq.RIGID_LAYER_INDEX) then
       ! do nothing
      else if (layer_iter.eq.ELASTIC_LAYER_INDEX) then

       if ((tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
           (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC)) then
        ! do nothing; uncaptured_volume_elastic remains to represent
        ! the original uncaptured space.
        if (tessellate_source.eq.TESSELLATE_FLUIDS) then
         !do nothing
        else
         print *,"expecting tessellate_source=TESSELLATE_FLUIDS ", &
            tessellate_source
         stop
        endif
       else if (tessellate_dest.eq.TESSELLATE_ALL) then
        ! modify the uncaptured regions to recognize the presence
        ! of the is_rigid==1 materials.
        uncaptured_volume(layer_iter)= &
            uncaptured_volume(RIGID_LAYER_INDEX)
        uncaptured_volume(FLUID_LAYER_INDEX)= &
            uncaptured_volume(RIGID_LAYER_INDEX)
        do dir=1,sdim
         uncaptured_centroid(layer_iter,dir)= &
             uncaptured_centroid(RIGID_LAYER_INDEX,dir)
         uncaptured_centroid(FLUID_LAYER_INDEX,dir)= &
             uncaptured_centroid(RIGID_LAYER_INDEX,dir)
        enddo
        if (tessellate_source.eq.TESSELLATE_FLUIDS) then
         !do nothing
        else
         print *,"expecting tessellate_source=TESSELLATE_FLUIDS ", &
            tessellate_source
         stop
        endif
       else if (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID) then
        ! do nothing (is_rigid_local(im)=0 for all materials)
        if (uncaptured_volume(layer_iter).eq. &
            uncaptured_volume(RIGID_LAYER_INDEX)) then
         !do nothing
        else
         print *,"expect uncaptured_volume(layer_iter).eq.(RIGID_LAYER_INDEX)"
         stop
        endif
        if (uncaptured_volume_fraction(RIGID_LAYER_INDEX).eq.one) then
         !do nothing
        else
         print *,"expect uncaptured_volume_fraction(RIGID_LAYER_INDEX).eq.one"
         stop
        endif
        if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then
         !do nothing
        else
         print *,"expect tessellate_source=TESSELLATE_IGNORE_ISRIGID ", &
            tessellate_source
         stop
        endif
       else if (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC) then 
        ! do nothing; uncaptured_volume_elastic remains to represent
        ! the original uncaptured space.
        if (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC) then
         !do nothing
        else
         print *,"expect tessellate_source=TESSELLATE_IGNORE_ISELASTIC ", &
            tessellate_source
         stop
        endif
       else if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then
        ! do nothing (is_rigid_local(im)=1 materials all zapped out)
        if (uncaptured_volume(layer_iter).eq. &
            uncaptured_volume(RIGID_LAYER_INDEX)) then
         !do nothing
        else
         print *,"expect uncaptured_volume(layer_iter).eq.(RIGID_LAYER_INDEX)"
         stop
        endif
        if (uncaptured_volume_fraction(RIGID_LAYER_INDEX).eq.one) then
         !do nothing
        else
         print *,"expecting uncaptured_volume_fraction(RIGID_LAYER_INDEX).eq.one"
         stop
        endif
        if ((tessellate_source.eq.TESSELLATE_FLUIDS).or. &
            (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC)) then
         !do nothing
        else
         print *,"expecting tessellate_source=TESSELLATE_IGNORE_ISELASTIC ", &
            tessellate_source
         print *,"or expecting tessellate_source=TESSELLATE_FLUIDS ", &
            tessellate_source
         stop
        endif

       else
        print *,"tessellate_dest bad(multi_get_volume_grid) ", &
         tessellate_dest
        stop
       endif

      else if (layer_iter.eq.FLUID_LAYER_INDEX) then

       if (tessellate_dest.eq.TESSELLATE_FLUIDS) then
        ! do nothing; uncaptured_volume_fluid remains to represent
        ! the original uncaptured space.
        if (tessellate_source.eq.TESSELLATE_FLUIDS) then
         !do nothing
        else
         print *,"expecting tessellate_source=TESSELLATE_FLUIDS ", &
            tessellate_source
         stop
        endif
       else if ((tessellate_dest.eq.TESSELLATE_ALL).or. &
                (tessellate_dest.eq.TESSELLATE_ALL_RASTER).or. &
                (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC)) then

        uncaptured_volume(layer_iter)= &
              uncaptured_volume(ELASTIC_LAYER_INDEX)
        do dir=1,sdim
         uncaptured_centroid(layer_iter,dir)= &
            uncaptured_centroid(ELASTIC_LAYER_INDEX,dir)
        enddo
        if ((tessellate_source.eq.TESSELLATE_FLUIDS).or. &
            (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC)) then
         !do nothing
        else
         print *,"expecting tessellate_source=TESSELLATE_IGNORE_ISELASTIC ", &
            tessellate_source
         print *,"or expecting tessellate_source=TESSELLATE_FLUIDS ", &
            tessellate_source
         stop
        endif

       else if (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID) then
        ! do nothing (is_rigid_local(im)=0 for all materials)
        if (uncaptured_volume(layer_iter).eq. &
            uncaptured_volume(ELASTIC_LAYER_INDEX)) then
         !do nothing
        else
         print *,"expect uncaptured_volume(layer_iter).eq.(ELASTIC_LAYER_INDEX)"
         stop
        endif
        if (uncaptured_volume_fraction(ELASTIC_LAYER_INDEX).eq.one) then
         !do nothing
        else
         print *,"expect uncaptured_volume_fraction(ELASTIC_LAYER_INDEX).eq.one"
         stop
        endif
        if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then
         !do nothing
        else
         print *,"expecting tessellate_source=TESSELLATE_IGNORE_ISRIGID ", &
            tessellate_source
         stop
        endif

       else if (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC) then 
        ! do nothing; uncaptured_volume_fluid remains to represent
        ! the original uncaptured space.
        if (uncaptured_volume(layer_iter).eq. &
            uncaptured_volume(ELASTIC_LAYER_INDEX)) then
         !do nothing
        else
         print *,"expect uncaptured_volume(layer_iter).eq.(ELASTIC_LAYER_INDEX)"
         stop
        endif
        if (uncaptured_volume_fraction(ELASTIC_LAYER_INDEX).eq.one) then
         !do nothing
        else
         print *,"expect uncaptured_volume_fraction(ELASTIC_LAYER_INDEX).eq.one"
         stop
        endif
        if (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC) then
         !do nothing
        else
         print *,"expecting tessellate_source=TESSELLATE_IGNORE_ISELASTIC ", &
            tessellate_source
         stop
        endif
       else
        print *,"tessellate_dest bad(multi_get_volume_grid) ", &
         tessellate_dest
        stop
       endif

      else
       print *,"layer_iter invalid(multi_get_volume_grid): ",layer_iter
       stop
      endif

      return 
      end subroutine advance_uncaptured_vars

      integer function continue_material_processing( &
           loop_counter, &
           local_num_materials, &
           num_processed, &
           uncaptured_volume_fraction, &
           vfrac_sum_local, &
           uncaptured_volume)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, intent(in) :: loop_counter
      integer, intent(in) :: local_num_materials
      integer, intent(in) :: num_processed
      real(amrex_real), intent(in) :: uncaptured_volume_fraction
      real(amrex_real), intent(in) :: vfrac_sum_local
      real(amrex_real), intent(in) :: uncaptured_volume

      if ((loop_counter.ge.0).and. &
          (loop_counter.lt. &
           local_num_materials).and. &
          (num_processed.ge.0).and. &
          (num_processed.lt. &
           local_num_materials).and. &
          (uncaptured_volume_fraction.le.one).and. &
          ((one-VOFTOL_MATERIAL)* &
           uncaptured_volume_fraction.gt. &
           one-vfrac_sum_local).and. &
          (uncaptured_volume.gt.zero)) then
       continue_material_processing=1
      else if ((loop_counter.eq.local_num_materials).or. &
               (num_processed.eq.local_num_materials).or. &
               ((one-VOFTOL_MATERIAL)* &
                uncaptured_volume_fraction.le. &
                one-vfrac_sum_local).or. &
               (uncaptured_volume.eq.zero)) then
       continue_material_processing=0
      else
       print *,"continue_material_processing parameters invalid"
       print *,"loop_counter ",loop_counter
       print *,"local_num_materials ",local_num_materials
       print *,"num_processed ",num_processed
       print *,"uncaptured_volume_fraction ",uncaptured_volume_fraction
       print *,"vfrac_sum_local ",vfrac_sum_local
       print *,"uncaptured_volume ",uncaptured_volume
       stop
      endif

      end function continue_material_processing

      subroutine check_for_single_material( &
         single_material, & !intent(out)
         remaining_vfrac, & !intent(out)
         num_processed_total, & !intent(out)
         material_used, &
         is_masked, &
         mofdatasave, &
         mofdatalocal, & !intent(inout)
         EPS_LOCAL, &
         uncaptured_volume_fraction, &
         layer_iter, &
         tessellate_source, &
         tessellate_dest, &
         sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, intent(out) :: single_material
      real(amrex_real), intent(out) :: remaining_vfrac
      integer, intent(out) :: num_processed_total
      real(amrex_real), intent(in) :: EPS_LOCAL
      integer, intent(in) :: material_used(num_materials)
      integer, intent(in) :: is_masked(num_materials)
      integer, intent(in) :: sdim
      integer, intent(in) :: tessellate_source
      integer, intent(in) :: tessellate_dest
      integer, intent(in) :: layer_iter
      real(amrex_real), intent(in) :: mofdatasave(num_materials*(2*sdim+3))
      real(amrex_real), intent(inout) :: mofdatalocal(num_materials*(2*sdim+3))
      real(amrex_real), intent(in) :: &
          uncaptured_volume_fraction(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      integer im_test,vofcomp,vofcomp_single


      remaining_vfrac=zero
      single_material=0
      num_processed_total=0

      do im_test=1,num_materials
       vofcomp=(im_test-1)*ngeom_recon+1

       if ((material_used(im_test).eq.0).and. &
           (is_masked(im_test).eq.0)) then
        !cut this material from avail if used

        if (mofdatasave(vofcomp).gt. &
            (one-EPS_LOCAL)*uncaptured_volume_fraction(layer_iter)) then

         if (single_material.eq.0) then
          single_material=im_test
         else if ((single_material.ge.1).and. &
                  (single_material.le.num_materials)) then
          vofcomp_single=(single_material-1)*ngeom_recon+1
          if (mofdatasave(vofcomp_single).lt. &
              mofdatasave(vofcomp)) then
           single_material=im_test
          else if (mofdatasave(vofcomp_single).ge. &
                   mofdatasave(vofcomp)) then
           !do nothing
          else
           print *,"mofdatasave bad (multi_get_volume_grid_simple) ", &
             mofdatasave
           stop
          endif
         else
          print *,"single_material bad(multi_get_volume_grid_simple)", &
           single_material
          stop
         endif

        else
         remaining_vfrac=remaining_vfrac+mofdatasave(vofcomp)
        endif

       else if (((material_used(im_test).ge.1).and. &
                 (material_used(im_test).le.num_materials)).or. &
                (is_masked(im_test).eq.1)) then
        ! do nothing
       else
        print *,"material used bust: ",material_used
        print *,"im_test: ",im_test
        print *,"is_masked ",is_masked
        stop
       endif

       mofdatalocal(vofcomp+sdim+1)=zero ! order=0
             
       if ((is_masked(im_test).eq.0).or. &!cut this material from avail if used
           (is_masked(im_test).eq.1)) then !do not cut this
        if ((material_used(im_test).ge.1).and. &
            (material_used(im_test).le.num_materials)) then
         mofdatalocal(vofcomp+sdim+1)=material_used(im_test)
         num_processed_total=num_processed_total+1
        else if (material_used(im_test).eq.0) then
         ! do nothing
        else
         print *,"material_used invalid: ",material_used
         stop
        endif
       else
        print *,"is_masked invalid ",is_masked
        stop
       endif

      enddo  ! im_test=1..num_materials

      return 
      end subroutine check_for_single_material


      subroutine sanity_check_multi_get_volume( &
       tessellate_source, &
       tessellate_dest, &
       multi_volume, &
       uncaptured_volume_START, &
       uncaptured_volume, &
       volcell, &
       EPS_SANITY, &
       EPS_SINGLE, &
       caller_id, &
       mofdata, &
       uncaptured_volume_fraction, &
       vfrac_sum_local, &
       xsten0, &
       nhalf0, &
       xsten_grid, &
       nhalf_grid, &
       is_rigid_local, &
       is_elastic_local, &
       sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, intent(in) :: tessellate_source
      integer, intent(in) :: tessellate_dest
      real(amrex_real), intent(in) :: multi_volume(num_materials)
      real(amrex_real), intent(in) :: uncaptured_volume_START
      real(amrex_real), intent(in) :: &
          uncaptured_volume(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real), intent(in) :: &
         uncaptured_volume_fraction( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real), intent(in) :: &
         vfrac_sum_local( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real), intent(in) :: volcell
      real(amrex_real) :: volcell_compare
      real(amrex_real), intent(in) :: EPS_SANITY,EPS_SINGLE
      integer, intent(in) :: caller_id
      real(amrex_real), INTENT(in) :: mofdata(num_materials*(2*sdim+3))
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      integer, INTENT(in) :: nhalf0,nhalf_grid,sdim
      integer, intent(in) :: is_rigid_local(num_materials)
      integer, intent(in) :: is_elastic_local(num_materials)
      integer :: layer_iter,im,vofcomp
 
      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
  
       if ((vfrac_sum_local(layer_iter).ge.one-VOFTOL_MATERIAL).or. &
           (layer_iter.eq.FLUID_LAYER_INDEX)) then
 
        if (uncaptured_volume(layer_iter).le.four*EPS_SANITY*volcell) then
         !do nothing
        else
         print *,"not all volume accounted for multi_get_volume_grid"
         print *,"caller_id=",caller_id
         print *,"sdim=",sdim
         print *,"layer_iter=",layer_iter
         do im=1,num_materials
          vofcomp=(im-1)*ngeom_recon+1
          print *,"im,vofcomp,mofdata ",im,vofcomp,mofdata(vofcomp)
         enddo
         print *,"mofdata=",mofdata
         print *,"EPS_SINGLE=",EPS_SINGLE
         print *,"uncaptured_volume ",uncaptured_volume
         print *,"uncaptured_volume_fraction ", &
           uncaptured_volume_fraction
         print *,"vfrac_sum_local ",vfrac_sum_local
         print *,"volcell ",volcell
         print *,"fraction of uncapt volume ", &
           uncaptured_volume(layer_iter)/volcell
         print *,"tolerance: ",four*EPS_SANITY
         print *,"xsten0= ",xsten0
         print *,"nhalf0= ",nhalf0
         print *,"xsten_grid= ",xsten_grid
         print *,"nhalf_grid= ",nhalf_grid
         print *,"is_rigid_local=",is_rigid_local
         print *,"is_elastic_local=",is_elastic_local
         stop
        endif
       else if ((vfrac_sum_local(layer_iter).le.one-VOFTOL_MATERIAL).and. &
                (vfrac_sum_local(layer_iter).ge.zero).and. &
                (layer_iter.ne.FLUID_LAYER_INDEX)) then
        !do nothing
       else
        print *,"parameter bust: ",vfrac_sum_local
        print *,"layer_iter ",layer_iter
        stop
       endif

      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      if ((tessellate_source.eq.TESSELLATE_FLUIDS).or. &
          (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID).or. &
          (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC)) then
       volcell_compare=zero
       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        if ((is_rigid(im).eq.0).and.(is_elastic(im).eq.0)) then
         volcell_compare=volcell_compare+multi_volume(im)
        else if (is_rigid(im).eq.1) then
         if (tessellate_dest.eq.TESSELLATE_FLUIDS) then
          !do nothing
         else if (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC) then
          !do nothing
         else if (tessellate_dest.eq.TESSELLATE_ALL) then
          volcell_compare=volcell_compare+multi_volume(im)
         else if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then
          volcell_compare=volcell_compare+multi_volume(im)
         else if (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID) then
          volcell_compare=volcell_compare+multi_volume(im)
         else if (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC) then
          !do nothing
         else
          print *,"tessellate_dest invalid ",tessellate_dest
          stop
         endif
        else if (is_elastic(im).eq.1) then
         if (tessellate_dest.eq.TESSELLATE_FLUIDS) then
          !do nothing
         else if (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC) then
          volcell_compare=volcell_compare+multi_volume(im)
         else if (tessellate_dest.eq.TESSELLATE_ALL) then
          volcell_compare=volcell_compare+multi_volume(im)
         else if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then
          volcell_compare=volcell_compare+multi_volume(im)
         else if (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID) then
          volcell_compare=volcell_compare+multi_volume(im)
         else if (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC) then
          volcell_compare=volcell_compare+multi_volume(im)
         else
          print *,"tessellate_dest invalid ",tessellate_dest
          stop
         endif
        else
         print *,"is_rigid(im) invalid ",im,is_rigid(im)
         print *,"or is_elastic(im) invalid ",im,is_elastic(im)
         stop
        endif
       enddo !im=1,num_materials

       if (abs(uncaptured_volume_START-volcell_compare).le. &
           EPS2*uncaptured_volume_START) then
        !do nothing
       else
        print *,"tessellate_source ",tessellate_source
        print *,"tessellate_dest ",tessellate_dest
        print *,"uncaptured_volume_fraction ", &
         uncaptured_volume_fraction
        print *,"multi_volume ",multi_volume
        print *,"uncaptured_volume ",uncaptured_volume
        print *,"is_rigid_local ",is_rigid_local
        print *,"is_elastic_local ",is_elastic_local
        print *,"volcell= ",volcell
        print *,"uncaptured_volume_START or volcell_compare invalid "
        print *,"uncaptured_volume_START= ",uncaptured_volume_START
        print *,"volcell_compare= ",volcell_compare
        print *,"relative error: ", &
          (uncaptured_volume_START-volcell_compare)/uncaptured_volume_START
        stop
       endif
      endif

      return

      end subroutine sanity_check_multi_get_volume

      subroutine multi_get_volume_grid_simple( &
       tid_in, &
       EPS_SINGLE, &
       tessellate_source, & 
       tessellate_dest, & 
       bfact,dx, &
       xsten0,nhalf0, & ! phi = n dot (x-x0) + intercept
       mofdata, &
       xsten_grid,nhalf_grid, & ! find volumes within xsten_grid or,
       multi_volume, &
       multi_cen, &
       xtetlist, &
       nlist_alloc, &
       nmax, &
       sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, PARAMETER :: caller_id_local=999
      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: nmax
      integer, INTENT(in) :: tessellate_source
      integer, INTENT(in) :: tessellate_dest
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: nhalf0,nhalf_grid
      real(amrex_real), INTENT(in) :: EPS_SINGLE
      real(amrex_real), INTENT(in) :: mofdata(num_materials*(2*sdim+3))
      real(amrex_real) mofdatalocal(num_materials*(2*sdim+3))
      real(amrex_real) mofdatasave(num_materials*(2*sdim+3))
      real(amrex_real) mofdatavalid(num_materials*(2*sdim+3))
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      real(amrex_real), INTENT(out) :: multi_volume(num_materials)
      real(amrex_real), INTENT(out) :: multi_cen(sdim,num_materials)
      integer dir
      integer vofcomp
      integer im
      integer layer_iter

      real(amrex_real) :: uncaptured_volume_START

      real(amrex_real) &
          uncaptured_volume(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) &
          uncaptured_centroid(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX,sdim)
      real(amrex_real) uncaptured_centroid_local(sdim)
      real(amrex_real) uncaptured_volume_fraction( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer at_least_one( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) volcell
      real(amrex_real) cencell(sdim)
      real(amrex_real) volcut,cencut(sdim)
      integer nlist
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      integer critical_material
      real(amrex_real) nrecon(sdim)
      real(amrex_real) intercept
      real(amrex_real) voltemp,centemp(sdim)
      integer single_material
      real(amrex_real) remaining_vfrac
      real(amrex_real) uncaptured_volume_save

      integer material_used(num_materials)

      integer im_raster_solid
      integer return_raster_info
      integer fastflag
      integer layer_flag

      integer local_num_materials( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) vfrac_sum_local( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) vfrac_sum( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer num_processed( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer num_processed_total

      integer sanity_num_processed_total
      integer layer_sanity

      integer loop_counter
      integer is_masked(num_materials)
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) vfrac_raster_solid

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid (multi_get_volume_grid_simple) ",tid_in
       stop
      endif 
      if ((EPS_SINGLE.gt.zero).and. &
          (EPS_SINGLE.lt.half)) then
       !do nothing
      else
       print *,"EPS_SINGLE invalid: ",EPS_SINGLE
       stop
      endif

      if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then
       if (tessellate_source.eq.tessellate_dest) then
        !do nothing
       else
        print *,"expecting tessellate_source==tessellate_dest"
        print *,"tessellate_source=",tessellate_source
        print *,"tessellate_dest=",tessellate_dest
        stop
       endif
      else if (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC) then
       if ((tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC).or. &
           (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
        !do nothing
       else
        print *,"tessellate_source or tessellate_dest invalid"
        print *,"tessellate_source=",tessellate_source
        print *,"tessellate_dest=",tessellate_dest
        stop
       endif
      else if (tessellate_source.eq.TESSELLATE_FLUIDS) then
       if ((tessellate_dest.eq.TESSELLATE_ALL).or. &
           (tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
           (tessellate_dest.eq.TESSELLATE_ALL_RASTER).or. &
           (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC)) then
        ! do nothing
       else
        print *,"tessellate_dest invalid: ",tessellate_dest
        stop
       endif
      else
       print *,"tessellate_source invalid10: ",tessellate_source
       stop
      endif

      call init_local_material_vars( & !subroutine multi_get_volume_grid_simple
         is_rigid_local, &
         is_elastic_local, &
         tessellate_source, &
         tessellate_dest, &
         is_proper_layer_local)

      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon.ne.2*sdim+3: ",ngeom_recon
       stop
      endif

      if (nmax.lt.4) then
       print *,"nmax invalid multi_get_volume_grid_simple nmax=",nmax
       stop
      endif
      if ((nhalf0.lt.1).or.(nhalf_grid.lt.1)) then
       print *,"nhalf invalid multi get volume grid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid135"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_volume_grid_simple"
       stop
      endif
      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid multi get volume grid"
       stop
      endif

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif
 
      do im=1,num_materials
       multi_volume(im)=zero
       do dir=1,sdim
        multi_cen(dir,im)=zero
       enddo
      enddo

       ! sum Frigid <=1
       ! sum Ffluid = 1
       ! multi_get_volume_grid_simple
      call make_vfrac_sum_ok_copy( &
        xsten0,nhalf0, &
        bfact,dx, &
        tessellate_source, &!TESSELLATE_FLUIDS|IGNORE_ISELASTIC|IGNORE_ISRIGID
        mofdata, &
        mofdatavalid,sdim)

      do dir=1,num_materials*ngeom_recon
       mofdatalocal(dir)=mofdatavalid(dir)
       mofdatasave(dir)=mofdatavalid(dir)
      enddo

      call Box_volumeFAST( &
       bfact, &
       dx, &
       xsten0, &
       nhalf0, &
       volcell, &
       cencell, &
       sdim)

      call Box_volumeFAST( &
         bfact, &
         dx, &
         xsten_grid, &
         nhalf_grid, &
         uncaptured_volume_START, &
         uncaptured_centroid_local,sdim)

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       uncaptured_volume(layer_iter)=uncaptured_volume_START
       do dir=1,sdim
        uncaptured_centroid(layer_iter,dir)= &
            uncaptured_centroid_local(dir)
       enddo !dir=1,sdim
       if (uncaptured_volume(layer_iter).ge.zero) then
        !do nothing
       else
        print *,"uncaptured_volume(layer_iter) invalid ", &
         layer_iter,uncaptured_volume(layer_iter)
        stop
       endif
      enddo ! layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      if (volcell.gt.zero) then
       !do nothing
      else
       print *,"volcell invalid multigetvolume grid: ",volcell
       stop
      endif

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       vfrac_sum(layer_iter)=zero
       vfrac_sum_local(layer_iter)=zero
       local_num_materials(layer_iter)=0
      enddo ! layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      im_raster_solid=0
      vfrac_raster_solid=zero

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1

       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

        if (is_proper_layer(im,layer_iter).eq.1) then
         vfrac_sum(layer_iter)= &
            vfrac_sum(layer_iter)+mofdatasave(vofcomp)
        endif

        if (is_proper_layer_local(im,layer_iter).eq.1) then

         vfrac_sum_local(layer_iter)= &
            vfrac_sum_local(layer_iter)+mofdatasave(vofcomp)

         local_num_materials(layer_iter)=local_num_materials(layer_iter)+1

         if (layer_iter.eq.RIGID_LAYER_INDEX) then

          if (im_raster_solid.eq.0) then
           im_raster_solid=im
           vfrac_raster_solid=mofdatasave(vofcomp)
          else if ((im_raster_solid.ge.1).and. &
                   (im_raster_solid.le.num_materials).and. &
                   (is_rigid_local(im_raster_solid).eq.1)) then
           if (vfrac_raster_solid.lt.mofdatasave(vofcomp)) then
            im_raster_solid=im
            vfrac_raster_solid=mofdatasave(vofcomp)
           endif
          else
           print *,"im_raster_solid invalid: ",im_raster_solid
           stop
          endif
      
         endif !if (layer_iter.eq.RIGID_LAYER_INDEX) 
        else if (is_proper_layer_local(im,layer_iter).eq.0) then
         !do nothing
        else
         print *,"is_proper_layer_local(im,layer_iter) invalid"
         stop
        endif
       enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
      enddo ! im=1..num_materials

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       if ((vfrac_sum(layer_iter).ge.zero).and. &
           (vfrac_sum(layer_iter).le.one+VOFTOL_MATERIAL)) then
        !do nothing
       else
        print *,"vfrac_sum invalid: ",vfrac_sum
        print *,"layer_iter= ",layer_iter
        stop
       endif
       !tessellate_source=TESSELLATE_FLUIDS|IGNORE_ISELASTIC|IGNORE_ISRIGID
       if ((vfrac_sum_local(layer_iter).ge.zero).and. &
           (vfrac_sum_local(layer_iter).le.one+VOFTOL_MATERIAL)) then
        !do nothing
       else
        print *,"vfrac_sum_local invalid: ",vfrac_sum
        print *,"layer_iter= ",layer_iter
        stop
       endif
      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      if (local_num_materials(FLUID_LAYER_INDEX)+ &
          local_num_materials(ELASTIC_LAYER_INDEX)+ &
          local_num_materials(RIGID_LAYER_INDEX).ne.num_materials) then
       print *,"local_num_materials invalid ",local_num_materials
       stop
      endif

      if (abs(one-vfrac_sum_local(FLUID_LAYER_INDEX)).le.EPS1) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif

      if ((vfrac_sum_local(FLUID_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
          (vfrac_sum_local(FLUID_LAYER_INDEX).ge.zero)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif
      if ((vfrac_sum_local(ELASTIC_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
          (vfrac_sum_local(ELASTIC_LAYER_INDEX).ge.zero)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif
      if ((vfrac_sum_local(RIGID_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
          (vfrac_sum_local(RIGID_LAYER_INDEX).ge.zero)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif

      return_raster_info=0

       ! if majority=solid => all solid.
      if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then 

       if ((tessellate_source.eq.TESSELLATE_FLUIDS).or. &
           (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC)) then
        !do nothing
       else
        print *,"tessellate_source invalid ",tessellate_source
        stop
       endif

       if (vfrac_sum(RIGID_LAYER_INDEX).eq. &
           vfrac_sum_local(RIGID_LAYER_INDEX)) then
        !do nothing
       else
        print *,"vfrac_sum<>vfrac_sum_local"
        print *,"vfrac_sum ",vfrac_sum
        print *,"vfrac_sum_local ",vfrac_sum_local
        stop
       endif

       if (vfrac_sum_local(RIGID_LAYER_INDEX).ge.half) then

        return_raster_info=1

        if ((im_raster_solid.ge.1).and. &
            (im_raster_solid.le.num_materials)) then
         vofcomp=(im_raster_solid-1)*ngeom_recon+1
         multi_volume(im_raster_solid)=uncaptured_volume(RIGID_LAYER_INDEX)
         do dir=1,sdim
          multi_cen(dir,im_raster_solid)= &
             uncaptured_centroid(RIGID_LAYER_INDEX,dir)
         enddo
        else
         print *,"im_raster_solid invalid: ",im_raster_solid
         stop
        endif

       else if (vfrac_sum_local(RIGID_LAYER_INDEX).le.half) then

        vfrac_sum_local(RIGID_LAYER_INDEX)=zero
        vfrac_sum(RIGID_LAYER_INDEX)=zero
        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         if (is_rigid_local(im).eq.0) then
          ! do nothing
         else if (is_rigid_local(im).eq.1) then
          do dir=0,sdim
           mofdatasave(vofcomp+dir)=zero
           mofdatavalid(vofcomp+dir)=zero
           mofdatalocal(vofcomp+dir)=zero
          enddo
         else
          print *,"is_rigid_local invalid: ",im,is_rigid_local(im)
          stop
         endif
        enddo ! im=1..num_materials

       else
        print *,"vfrac_sum_local invalid: ",vfrac_sum_local
        stop
       endif

      else if ((tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
               (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC).or. &
               (tessellate_dest.eq.TESSELLATE_ALL).or. &
               (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC).or. &
               (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID)) then
       ! do nothing
      else
       print *,"tessellate_dest bad (multi_get_volume_grid_simple): ", &
               tessellate_dest
       stop
      endif
      
      if (return_raster_info.eq.1) then

       ! do nothing (in this case F_solid>1/2 and TESSELLATE_ALL_RASTER)

      else if (return_raster_info.eq.0) then
 
       num_processed_total=0

        !multi_get_volume_grid_simple
       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
        uncaptured_volume_fraction(layer_iter)=one
        num_processed(layer_iter)=0
       enddo

       do im=1,num_materials
        material_used(im)=0
       enddo ! im=1..num_materials

       fastflag=1

        ! intersection of departure region with grid cell is 
        ! very small.
       if (uncaptured_volume_START.le.EPS_12_6*volcell) then

        call tiny_fraction_case( &
          mofdatasave, &
          uncaptured_volume_START, &
          uncaptured_centroid, &
          multi_volume, &
          multi_cen, &
          is_rigid_local, &
          is_elastic_local, &
          tessellate_source, &
          tessellate_dest, &
          vfrac_sum_local, &
          sdim)

       else if (uncaptured_volume_START.ge.EPS_12_6*volcell) then

        do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

         at_least_one(layer_iter)=0
         do im=1,num_materials
          if (is_proper_layer_local(im,layer_iter).eq.1) then
           at_least_one(layer_iter)=at_least_one(layer_iter)+1
          endif
         enddo

         call init_layer_flag( & !multi_get_volume_grid_simple
          at_least_one(layer_iter), &
          vfrac_sum, &
          vfrac_sum_local, &
          layer_flag, & !intent(out)
          layer_iter, & !intent(in)
          tessellate_source, &
          tessellate_dest)
            
         call init_mask_flag( &
            layer_flag, &
            is_masked)

         call advance_uncaptured_vars( & !multi_get_volume_grid_simple
          at_least_one(layer_iter), &
          tessellate_source, &
          tessellate_dest, &
          layer_iter, &
          uncaptured_volume, &
          uncaptured_centroid, &
          uncaptured_volume_fraction, &
          sdim)


         if ((at_least_one(layer_iter).ge.1).and. &
             (at_least_one(layer_iter).le.num_materials).and. &
             (vfrac_sum_local(layer_iter).gt.zero)) then

           !uncaptured_volume(layer_iter) initialized to the box volume
           !for all "layer_iter"
           ! 
           !uncaptured_volume_fraction(layer_iter) initialized to 1.0
           !for all "layer_iter"
          loop_counter=0
          do while (continue_material_processing( &
                loop_counter, &
                local_num_materials(layer_iter), &
                num_processed(layer_iter), &
                uncaptured_volume_fraction(layer_iter), &
                vfrac_sum_local(layer_iter), &
                uncaptured_volume(layer_iter)).eq.1) 

           call check_for_single_material( & !multi_get_volume_grid_simple
            single_material, &
            remaining_vfrac, &
            num_processed_total, &
            material_used, &
            is_masked, &
            mofdatasave, &
            mofdatalocal, &
            EPS_SINGLE, &
            uncaptured_volume_fraction, &
            layer_iter, &
            tessellate_source, &
            tessellate_dest, &
            sdim)

           sanity_num_processed_total=0
           do layer_sanity=RIGID_LAYER_INDEX,layer_iter
            sanity_num_processed_total=sanity_num_processed_total+ &
             num_processed(layer_sanity)
           enddo
           if (sanity_num_processed_total.eq.num_processed_total) then
            !do nothing
           else
            print *,"num_processed_total invalid ",num_processed_total
            print *,"sanity_num_processed_total ",sanity_num_processed_total
            stop
           endif

           if ((single_material.ge.1).and. &
               (single_material.le.num_materials).and. &
               (remaining_vfrac.le.EPS_SINGLE)) then

            vofcomp=(single_material-1)*ngeom_recon+1
            multi_volume(single_material)=uncaptured_volume(layer_iter)
            do dir=1,sdim
             multi_cen(dir,single_material)= &
              uncaptured_centroid(layer_iter,dir)
            enddo

            uncaptured_volume(layer_iter)=zero
            uncaptured_volume_fraction(layer_iter)=zero

            num_processed(layer_iter)=num_processed(layer_iter)+1
            num_processed_total=num_processed_total+1
            material_used(single_material)=num_processed_total

           else if ((single_material.eq.0).or. &
                    (remaining_vfrac.ge.EPS_SINGLE)) then

            if ((num_processed_total.gt.0).and. &
                (num_processed_total.lt.num_materials)) then
             fastflag=0
            else if (num_processed_total.eq.0) then
             fastflag=1
            else          
             print *,"num_processed_total invalid: ",num_processed_total
             stop
            endif

            if (fastflag.eq.0) then

             ! only xsten0(0,dir) dir=1..sdim used
             ! in: multi_volume_grid
             call tets_box_planes( &
               layer_flag, &
               bfact,dx, &
               xsten0,nhalf0, &
               xsten_grid,nhalf_grid, &
               mofdatalocal, & !only slope,intercept,order used.
               xtetlist, &
               nlist_alloc, &
               nlist, &
               nmax, &
               sdim)

             call get_cut_geom3D(xtetlist, &
               nlist_alloc,nlist,nmax, &
               volcut,cencut,sdim)

             if (abs(volcut-uncaptured_volume(layer_iter)).gt. &
                 EPS1*volcell) then
              print *,"volcut invalid multi_get_volume_grid_simple 1"
              print *,"volcut= ",volcut
              print *,"uncaptured_volume=",uncaptured_volume
              print *,"volcell= ",volcell
              print *,"VOFTOL_MATERIAL= ",VOFTOL_MATERIAL
              print *,"EPS_11_4= ",EPS_11_4
              print *,"EPS_12_6= ",EPS_12_6
              print *,"EPS_8_4= ",EPS_8_4
              print *,"EPS_SINGLE= ",EPS_SINGLE
              print *,"xsten0 ",xsten0(0,1),xsten0(0,2),xsten0(0,sdim)
              print *,"xsten_grid ",xsten_grid(0,1),xsten_grid(0,2), &
                xsten_grid(0,sdim)
              do im=1,num_materials
               vofcomp=(im-1)*ngeom_recon+1
               print *,"im,mofdatavalid(vofcomp) ",im,mofdatavalid(vofcomp)
              enddo 
              print *,"is_rigid_local=",is_rigid_local
              print *,"is_elastic_local=",is_elastic_local
              print *,"is_masked=",is_masked
              stop
             endif

            else if (fastflag.eq.1) then

             ! do nothing; unnecessary to intersect the original box with
             ! the compliment of materials already processed.

            else 
             print *,"fastflag bad multi_get_volume_grid_simple: ",fastflag
             stop
            endif

            call derive_critical_material( & !multi_get_volume_grid_simple
             mofdatalocal, &
             mofdatasave, &
             material_used, &
             is_masked, &
             critical_material, &
             sdim)

            if ((critical_material.ge.1).and. &
                (critical_material.le.num_materials)) then        
             vofcomp=(critical_material-1)*ngeom_recon+1
             do dir=1,sdim
              nrecon(dir)=mofdatalocal(vofcomp+sdim+1+dir)
             enddo
             intercept=mofdatalocal(vofcomp+2*sdim+2)

             if (fastflag.eq.0) then
              !only xsten0(0,dir) dir=1..sdim used
              !calling from "multi_get_volume_grid_simple"
              call multi_cell_intersection_simple( &
               bfact,dx, &
               xsten0,nhalf0, &
               nrecon,intercept, &
               voltemp,centemp, &
               xtetlist, &
               nlist_alloc, &
               nlist, &
               nmax, &
               sdim) 
             else if (fastflag.eq.1) then
              ! only xsten0(0,dir) dir=1..sdim used
              call fast_cut_cell_intersection_simple( &
               bfact,dx, &
               xsten0,nhalf0, &
               nrecon,intercept, &
               voltemp,centemp, &
               xsten_grid,nhalf_grid, &
               sdim) 
             else 
              print *,"fastflag invalid multi get volume grid 2"
              stop
             endif

             multi_volume(critical_material)=voltemp
             do dir=1,sdim
              if (voltemp.gt.zero) then
               multi_cen(dir,critical_material)=centemp(dir)
              else
               multi_cen(dir,critical_material)=zero
              endif
             enddo ! dir=1..sdim

             uncaptured_volume_save=uncaptured_volume(layer_iter)
             uncaptured_volume(layer_iter)=uncaptured_volume(layer_iter)- &
                voltemp
             if (uncaptured_volume(layer_iter).lt. &
                 EPS_SINGLE*uncaptured_volume_START) then
              uncaptured_volume(layer_iter)=zero
             endif

             ! V^{uncapt,k}=V+V^{uncapt,k+1}
             ! V^{uncapt,k}x^{uncapt,k}=V x+V^{uncapt,k+1}x^{uncapt,k+1}

             do dir=1,sdim
              if (uncaptured_volume(layer_iter).le.zero) then
               uncaptured_centroid(layer_iter,dir)=zero
              else
               uncaptured_centroid(layer_iter,dir)= &
                (uncaptured_volume_save*uncaptured_centroid(layer_iter,dir)- &
                 voltemp*centemp(dir))/uncaptured_volume(layer_iter)
              endif
             enddo ! dir=1..sdim
  
             !uncaptured_volume_fraction_solid=1.0 initially. 
             uncaptured_volume_fraction(layer_iter)= &
              uncaptured_volume_fraction(layer_iter)- &
              mofdatalocal(vofcomp)
             if (uncaptured_volume_fraction(layer_iter).lt. &
                 one-vfrac_sum_local(layer_iter)+EPS_SINGLE) then
              uncaptured_volume_fraction(layer_iter)= &
                one-vfrac_sum_local(layer_iter)
             endif

             num_processed(layer_iter)=num_processed(layer_iter)+1
             num_processed_total=num_processed_total+1
             material_used(critical_material)=num_processed_total

            else if (critical_material.eq.0) then
             ! do nothing
            else
             print *,"critical_material invalid 18455: ",critical_material
             stop
            endif
  
           else
            print *,"single_material or remaining_vfrac invalid ", &
             single_material,remaining_vfrac
            stop
           endif

           loop_counter=loop_counter+1
          enddo !while continue_material_processing==1

         else if ((at_least_one(layer_iter).eq.0).or. &
                  (vfrac_sum_local(layer_iter).eq.zero)) then
          !do nothing
         else
          print *,"at_least_one or vfrac_sum_local invalid ", &
           at_least_one,vfrac_sum_local
          stop
         endif
         
        enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

        call sanity_check_multi_get_volume( & !multi_get_volume_grid_simple
         tessellate_source, &
         tessellate_dest, &
         multi_volume, &
         uncaptured_volume_START, &
         uncaptured_volume, &
         volcell, &
         EPS2, &
         EPS_SINGLE, &
         caller_id_local, &
         mofdata, &
         uncaptured_volume_fraction, &
         vfrac_sum_local, &
         xsten0, &
         nhalf0, &
         xsten_grid, &
         nhalf_grid, &
         is_rigid_local, &
         is_elastic_local, &
         sdim)

       else
        print *,"uncaptured_volume_START invalid ",uncaptured_volume_START
        stop
       endif

      else
       print *,"return_raster_info invalid ",return_raster_info
       stop
      endif

      return
      end subroutine multi_get_volume_grid_simple



      subroutine multi_get_volume_grid( &
       caller_id, &
       tid_in, &
       EPS_SINGLE, &
       tessellate_source, & 
       tessellate_dest, & 
       bfact,dx, &
       xsten0,nhalf0, & ! phi = n dot (x-x0) + intercept
       mofdata, &
       xsten_grid,nhalf_grid, & ! find volumes within xsten_grid or,
       xtet, &                  ! within xtet
       multi_volume, &
       multi_cen, &
       multi_area, &
       xtetlist, &
       nlist_alloc, &
       nmax, &
       sdim, &
       shapeflag)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: caller_id
      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: nmax
      integer, INTENT(in) :: tessellate_source
      integer, INTENT(in) :: tessellate_dest
      integer, INTENT(in) :: shapeflag
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: nhalf0,nhalf_grid
      real(amrex_real), INTENT(in) :: EPS_SINGLE
      real(amrex_real), INTENT(in) :: xtet(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: mofdata(num_materials*(2*sdim+3))
      real(amrex_real) mofdatalocal(num_materials*(2*sdim+3))
      real(amrex_real) mofdatasave(num_materials*(2*sdim+3))
      real(amrex_real) mofdatavalid(num_materials*(2*sdim+3))
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      real(amrex_real), INTENT(out) :: multi_volume(num_materials)
      real(amrex_real), INTENT(out) :: multi_cen(sdim,num_materials)
      real(amrex_real), INTENT(out) :: multi_area(num_materials)
      integer dir
      integer vofcomp
      integer im
      integer layer_iter

      real(amrex_real) :: uncaptured_volume_START

      real(amrex_real) &
          uncaptured_volume(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) &
          uncaptured_centroid(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX,sdim)
      real(amrex_real) uncaptured_centroid_local(sdim)
      real(amrex_real) uncaptured_volume_fraction( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer at_least_one( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) volcell
      real(amrex_real) cencell(sdim)
      real(amrex_real) volcut,cencut(sdim)
      integer nlist
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      integer critical_material
      real(amrex_real) nrecon(sdim)
      real(amrex_real) intercept
      real(amrex_real) voltemp,centemp(sdim)
      real(amrex_real) areatemp
      integer single_material
      real(amrex_real) remaining_vfrac
      real(amrex_real) uncaptured_volume_save

      integer material_used(num_materials)

      integer im_raster_solid
      integer return_raster_info
      integer fastflag
      integer layer_flag

      integer local_num_materials( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) vfrac_sum_local( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) vfrac_sum( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer num_processed( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer num_processed_total

      integer sanity_num_processed_total
      integer layer_sanity

      integer loop_counter
      integer is_masked(num_materials)
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) vfrac_raster_solid

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid (multi_get_volume_grid) ",tid_in
       stop
      endif 
      if ((EPS_SINGLE.gt.zero).and. &
          (EPS_SINGLE.lt.half)) then
       !do nothing
      else
       print *,"EPS_SINGLE invalid: ",EPS_SINGLE
       stop
      endif

      if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then
       if (tessellate_source.eq.tessellate_dest) then
        !do nothing
       else
        print *,"expecting tessellate_source==tessellate_dest"
        print *,"tessellate_source=",tessellate_source
        print *,"tessellate_dest=",tessellate_dest
        stop
       endif
      else if (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC) then
       if ((tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC).or. &
           (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
        !do nothing
       else
        print *,"tessellate_source or tessellate_dest invalid"
        print *,"tessellate_source=",tessellate_source
        print *,"tessellate_dest=",tessellate_dest
        stop
       endif
      else if (tessellate_source.eq.TESSELLATE_FLUIDS) then
       if ((tessellate_dest.eq.TESSELLATE_ALL).or. &
           (tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
           (tessellate_dest.eq.TESSELLATE_ALL_RASTER).or. &
           (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC)) then
        ! do nothing
       else
        print *,"tessellate_dest invalid: ",tessellate_dest
        stop
       endif
      else
       print *,"tessellate_source invalid10: ",tessellate_source
       stop
      endif

      call init_local_material_vars( & !subroutine multi_get_volume_grid
         is_rigid_local, &
         is_elastic_local, &
         tessellate_source, &
         tessellate_dest, &
         is_proper_layer_local)

      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon.ne.2*sdim+3: ",ngeom_recon
       stop
      endif

      if (nmax.lt.4) then
       print *,"nmax invalid multi_get_volume_grid nmax=",nmax
       stop
      endif
      if ((nhalf0.lt.1).or.(nhalf_grid.lt.1)) then
       print *,"nhalf invalid multi get volume grid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid multi_get_volume_grid ",bfact
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_volume_grid"
       stop
      endif
      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid multi get volume grid"
       stop
      endif

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif
 
      do im=1,num_materials
       multi_volume(im)=zero
       do dir=1,sdim
        multi_cen(dir,im)=zero
       enddo
       multi_area(im)=zero
      enddo

       ! sum Frigid <=1
       ! sum Ffluid = 1
       ! multi_get_volume_grid
      call make_vfrac_sum_ok_copy( &
        xsten0,nhalf0, &
        bfact,dx, &
        tessellate_source, &!TESSELLATE_FLUIDS|IGNORE_ISELASTIC|IGNORE_ISRIGID
        mofdata, &
        mofdatavalid,sdim)

      do dir=1,num_materials*ngeom_recon
       mofdatalocal(dir)=mofdatavalid(dir)
       mofdatasave(dir)=mofdatavalid(dir)
      enddo

      call Box_volumeFAST( &
       bfact, &
       dx, &
       xsten0, &
       nhalf0, &
       volcell, &
       cencell, &
       sdim)

      if (shapeflag.eq.0) then
 
       call Box_volumeFAST(bfact,dx, &
         xsten_grid,nhalf_grid, &
         uncaptured_volume_START, &
         uncaptured_centroid_local,sdim)

      else if (shapeflag.eq.1) then

       call tetrahedron_volume(xtet, &
         uncaptured_volume_START, &
         uncaptured_centroid_local,sdim)

      else
       print *,"shapeflag invalid: ",shapeflag
       stop
      endif

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       uncaptured_volume(layer_iter)=uncaptured_volume_START
       do dir=1,sdim
        uncaptured_centroid(layer_iter,dir)= &
            uncaptured_centroid_local(dir)
       enddo !dir=1,sdim
       if (uncaptured_volume(layer_iter).ge.zero) then
        !do nothing
       else
        print *,"uncaptured_volume(layer_iter) invalid ", &
         layer_iter,uncaptured_volume(layer_iter)
        stop
       endif
      enddo ! layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      if (volcell.gt.zero) then
       !do nothing
      else
       print *,"volcell invalid multigetvolume grid: ",volcell
       stop
      endif

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       vfrac_sum(layer_iter)=zero
       vfrac_sum_local(layer_iter)=zero
       local_num_materials(layer_iter)=0
      enddo ! layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      im_raster_solid=0
      vfrac_raster_solid=zero

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1

       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

        if (is_proper_layer(im,layer_iter).eq.1) then
         vfrac_sum(layer_iter)= &
            vfrac_sum(layer_iter)+mofdatasave(vofcomp)
        endif

        if (is_proper_layer_local(im,layer_iter).eq.1) then

         vfrac_sum_local(layer_iter)= &
            vfrac_sum_local(layer_iter)+mofdatasave(vofcomp)

         local_num_materials(layer_iter)=local_num_materials(layer_iter)+1

         if (layer_iter.eq.RIGID_LAYER_INDEX) then

          if (im_raster_solid.eq.0) then
           im_raster_solid=im
           vfrac_raster_solid=mofdatasave(vofcomp)
          else if ((im_raster_solid.ge.1).and. &
                   (im_raster_solid.le.num_materials).and. &
                   (is_rigid_local(im_raster_solid).eq.1)) then
           if (vfrac_raster_solid.lt.mofdatasave(vofcomp)) then
            im_raster_solid=im
            vfrac_raster_solid=mofdatasave(vofcomp)
           endif
          else
           print *,"im_raster_solid invalid: ",im_raster_solid
           stop
          endif
      
         endif !if (layer_iter.eq.RIGID_LAYER_INDEX) 
        else if (is_proper_layer_local(im,layer_iter).eq.0) then
         !do nothing
        else
         print *,"is_proper_layer_local(im,layer_iter) invalid"
         stop
        endif
       enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
      enddo ! im=1..num_materials

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       if ((vfrac_sum(layer_iter).ge.zero).and. &
           (vfrac_sum(layer_iter).le.one+VOFTOL_MATERIAL)) then
        !do nothing
       else
        print *,"vfrac_sum invalid: ",vfrac_sum
        print *,"layer_iter= ",layer_iter
        stop
       endif
       !tessellate_source=TESSELLATE_FLUIDS|IGNORE_ISELASTIC|IGNORE_ISRIGID
       if ((vfrac_sum_local(layer_iter).ge.zero).and. &
           (vfrac_sum_local(layer_iter).le.one+VOFTOL_MATERIAL)) then
        !do nothing
       else
        print *,"vfrac_sum_local invalid: ",vfrac_sum
        print *,"layer_iter= ",layer_iter
        stop
       endif
      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      if (local_num_materials(FLUID_LAYER_INDEX)+ &
          local_num_materials(ELASTIC_LAYER_INDEX)+ &
          local_num_materials(RIGID_LAYER_INDEX).ne.num_materials) then
       print *,"local_num_materials invalid ",local_num_materials
       stop
      endif

      if (abs(one-vfrac_sum_local(FLUID_LAYER_INDEX)).le.EPS1) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif

      if ((vfrac_sum_local(FLUID_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
          (vfrac_sum_local(FLUID_LAYER_INDEX).ge.zero)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif
      if ((vfrac_sum_local(ELASTIC_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
          (vfrac_sum_local(ELASTIC_LAYER_INDEX).ge.zero)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif
      if ((vfrac_sum_local(RIGID_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
          (vfrac_sum_local(RIGID_LAYER_INDEX).ge.zero)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif

      return_raster_info=0

       ! if majority=solid => all solid.
      if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then 

       if ((tessellate_source.eq.TESSELLATE_FLUIDS).or. &
           (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC)) then
        !do nothing
       else
        print *,"tessellate_source invalid ",tessellate_source
        stop
       endif

       if (vfrac_sum(RIGID_LAYER_INDEX).eq. &
           vfrac_sum_local(RIGID_LAYER_INDEX)) then
        !do nothing
       else
        print *,"vfrac_sum<>vfrac_sum_local"
        print *,"vfrac_sum ",vfrac_sum
        print *,"vfrac_sum_local ",vfrac_sum_local
        stop
       endif

       if (vfrac_sum_local(RIGID_LAYER_INDEX).ge.half) then

        return_raster_info=1

        if ((im_raster_solid.ge.1).and. &
            (im_raster_solid.le.num_materials)) then
         vofcomp=(im_raster_solid-1)*ngeom_recon+1
         multi_volume(im_raster_solid)=uncaptured_volume(RIGID_LAYER_INDEX)
         do dir=1,sdim
          multi_cen(dir,im_raster_solid)= &
             uncaptured_centroid(RIGID_LAYER_INDEX,dir)
         enddo
        else
         print *,"im_raster_solid invalid: ",im_raster_solid
         stop
        endif

       else if (vfrac_sum_local(RIGID_LAYER_INDEX).le.half) then

        vfrac_sum_local(RIGID_LAYER_INDEX)=zero
        vfrac_sum(RIGID_LAYER_INDEX)=zero
        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         if (is_rigid_local(im).eq.0) then
          ! do nothing
         else if (is_rigid_local(im).eq.1) then
          do dir=0,sdim
           mofdatasave(vofcomp+dir)=zero
           mofdatavalid(vofcomp+dir)=zero
           mofdatalocal(vofcomp+dir)=zero
          enddo
         else
          print *,"is_rigid_local invalid: ",im,is_rigid_local(im)
          stop
         endif
        enddo ! im=1..num_materials

       else
        print *,"vfrac_sum_local invalid: ",vfrac_sum_local
        stop
       endif

      else if ((tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
               (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC).or. &
               (tessellate_dest.eq.TESSELLATE_ALL).or. &
               (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC).or. &
               (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID)) then
       ! do nothing
      else
       print *,"tessellate_dest bad (multi_get_volume_grid): ",tessellate_dest
       stop
      endif
      
      if (return_raster_info.eq.1) then

       ! do nothing (in this case F_solid>1/2 and TESSELLATE_ALL_RASTER)

      else if (return_raster_info.eq.0) then
 
       num_processed_total=0

        ! multi_get_volume_grid
       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
        uncaptured_volume_fraction(layer_iter)=one
        num_processed(layer_iter)=0
       enddo

       do im=1,num_materials
        material_used(im)=0
       enddo ! im=1..num_materials

       fastflag=1

        ! intersection of departure region with grid cell is 
        ! very small.
       if (uncaptured_volume_START.le.EPS_12_6*volcell) then

        call tiny_fraction_case( &
          mofdatasave, &
          uncaptured_volume_START, &
          uncaptured_centroid, &
          multi_volume, &
          multi_cen, &
          is_rigid_local, &
          is_elastic_local, &
          tessellate_source, &
          tessellate_dest, &
          vfrac_sum_local, &
          sdim)

       else if (uncaptured_volume_START.ge.EPS_12_6*volcell) then

        do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

         at_least_one(layer_iter)=0
         do im=1,num_materials
          if (is_proper_layer_local(im,layer_iter).eq.1) then
           at_least_one(layer_iter)=at_least_one(layer_iter)+1
          endif
         enddo

         call init_layer_flag( &
          at_least_one(layer_iter), &
          vfrac_sum, &
          vfrac_sum_local, &
          layer_flag, & !intent(out)
          layer_iter, & !intent(in)
          tessellate_source, &
          tessellate_dest)

         call init_mask_flag( &
          layer_flag, &
          is_masked)

         call advance_uncaptured_vars( & !multi_get_volume_grid
          at_least_one(layer_iter), &
          tessellate_source, &
          tessellate_dest, &
          layer_iter, &
          uncaptured_volume, &
          uncaptured_centroid, &
          uncaptured_volume_fraction, &
          sdim)

         if ((at_least_one(layer_iter).ge.1).and. &
             (at_least_one(layer_iter).le.num_materials).and. &
             (vfrac_sum_local(layer_iter).gt.zero)) then


           !uncaptured_volume(layer_iter) initialized to the box volume
           !for all "layer_iter"
           ! 
           !uncaptured_volume_fraction(layer_iter) initialized to 1.0
           !for all "layer_iter"
          loop_counter=0
          do while (continue_material_processing( &
                loop_counter, &
                local_num_materials(layer_iter), &
                num_processed(layer_iter), &
                uncaptured_volume_fraction(layer_iter), &
                vfrac_sum_local(layer_iter), &
                uncaptured_volume(layer_iter)).eq.1) 

           call check_for_single_material( & !multi_get_volume_grid
            single_material, &
            remaining_vfrac, &
            num_processed_total, &
            material_used, &
            is_masked, &
            mofdatasave, &
            mofdatalocal, &
            EPS_SINGLE, &
            uncaptured_volume_fraction, &
            layer_iter, &
            tessellate_source, &
            tessellate_dest, &
            sdim)

           sanity_num_processed_total=0
           do layer_sanity=RIGID_LAYER_INDEX,layer_iter
            sanity_num_processed_total=sanity_num_processed_total+ &
              num_processed(layer_sanity)
           enddo
           if (sanity_num_processed_total.eq.num_processed_total) then
            !do nothing
           else
            print *,"num_processed_total invalid ",num_processed_total
            print *,"sanity_num_processed_total ",sanity_num_processed_total
            stop
           endif

           if ((single_material.ge.1).and. &
               (single_material.le.num_materials).and. &
               (remaining_vfrac.le.EPS_SINGLE)) then

            vofcomp=(single_material-1)*ngeom_recon+1
            multi_volume(single_material)=uncaptured_volume(layer_iter)
            do dir=1,sdim
             multi_cen(dir,single_material)= &
              uncaptured_centroid(layer_iter,dir)
            enddo

            uncaptured_volume(layer_iter)=zero
            uncaptured_volume_fraction(layer_iter)=zero

            num_processed(layer_iter)=num_processed(layer_iter)+1
            num_processed_total=num_processed_total+1
            material_used(single_material)=num_processed_total

           else if ((single_material.eq.0).or. &
                    (remaining_vfrac.ge.EPS_SINGLE)) then

            if ((num_processed_total.gt.0).and. &
                (num_processed_total.lt.num_materials)) then
             fastflag=0
            else if (num_processed_total.eq.0) then
             fastflag=1
            else          
             print *,"num_processed_total invalid: ",num_processed_total
             stop
            endif

            if (fastflag.eq.0) then

             if (shapeflag.eq.0) then ! volumes in a box
              ! only xsten0(0,dir) dir=1..sdim used
              ! in: multi_volume_grid
              call tets_box_planes( &
               layer_flag, &
               bfact,dx, &
               xsten0,nhalf0, &
               xsten_grid,nhalf_grid, &
               mofdatalocal, & !only slope,intercept,order used.
               xtetlist, &
               nlist_alloc, &
               nlist, &
               nmax, &
               sdim)

             else if (shapeflag.eq.1) then ! volumes in a tet.

              ! only xsten0(0,dir) dir=1..sdim used
              ! xtetlist=xtet - highest order material
              call tets_tet_planes( &
                layer_flag, &
                bfact,dx, &
                xsten0,nhalf0, &
                xtet, &
                mofdatalocal, &
                xtetlist, &
                nlist_alloc, &
                nlist, &
                nmax, &
                sdim)

             else
              print *,"shapeflag invalid: ",shapeflag
              stop
             endif

             call get_cut_geom3D(xtetlist, &
               nlist_alloc,nlist,nmax, &
               volcut,cencut,sdim)

             if (abs(volcut-uncaptured_volume(layer_iter)).gt. &
                 EPS1*volcell) then
              print *,"volcut invalid multi_get_volume_grid 1"
              print *,"volcut= ",volcut
              print *,"uncaptured_volume=",uncaptured_volume
              print *,"volcell= ",volcell
              print *,"VOFTOL_MATERIAL= ",VOFTOL_MATERIAL
              print *,"EPS_11_4= ",EPS_11_4
              print *,"EPS_12_6= ",EPS_12_6
              print *,"EPS_8_4= ",EPS_8_4
              print *,"EPS_SINGLE= ",EPS_SINGLE
              print *,"xsten0 ",xsten0(0,1),xsten0(0,2),xsten0(0,sdim)
              print *,"xsten_grid ",xsten_grid(0,1),xsten_grid(0,2), &
                xsten_grid(0,sdim)
              do im=1,num_materials
               vofcomp=(im-1)*ngeom_recon+1
               print *,"im,mofdatavalid(vofcomp) ",im,mofdatavalid(vofcomp)
              enddo 
              print *,"is_rigid_local=",is_rigid_local
              print *,"is_elastic_local=",is_elastic_local
              print *,"is_masked=",is_masked
              stop
             endif

            else if (fastflag.eq.1) then

             ! do nothing; unnecessary to intersect the original box with
             ! the compliment of materials already processed.

            else 
             print *,"fastflag invalid multi get volume grid: ",fastflag
             stop
            endif

            call derive_critical_material( & !multi_get_volume_grid
             mofdatalocal, &
             mofdatasave, &
             material_used, &
             is_masked, &
             critical_material, &
             sdim)
          
            if ((critical_material.ge.1).and. &
                (critical_material.le.num_materials)) then        
             vofcomp=(critical_material-1)*ngeom_recon+1
             do dir=1,sdim
              nrecon(dir)=mofdatalocal(vofcomp+sdim+1+dir)
             enddo
             intercept=mofdatalocal(vofcomp+2*sdim+2)

             if (fastflag.eq.0) then
              !only xsten0(0,dir) dir=1..sdim used
              !calling from "multi_get_volume_grid"
              call multi_cell_intersection( &
               bfact,dx, &
               xsten0,nhalf0, &
               nrecon,intercept, &
               voltemp,centemp, &
               areatemp, &
               xtetlist, &
               nlist_alloc, &
               nlist, &
               nmax, &
               sdim) 
             else if (fastflag.eq.1) then
              ! only xsten0(0,dir) dir=1..sdim used
              call fast_cut_cell_intersection( &
               bfact,dx, &
               xsten0,nhalf0, &
               nrecon,intercept, &
               voltemp,centemp, &
               areatemp, &
               xsten_grid,nhalf_grid, &
               xtet,shapeflag, &
               sdim) 
             else 
              print *,"fastflag invalid multi get volume grid 2"
              stop
             endif

             multi_volume(critical_material)=voltemp
             do dir=1,sdim
              if (voltemp.gt.zero) then
               multi_cen(dir,critical_material)=centemp(dir)
              else
               multi_cen(dir,critical_material)=zero
              endif
             enddo ! dir=1..sdim
             multi_area(critical_material)=areatemp

             uncaptured_volume_save=uncaptured_volume(layer_iter)
             uncaptured_volume(layer_iter)=uncaptured_volume(layer_iter)- &
                voltemp
             if (uncaptured_volume(layer_iter).lt. &
                 EPS_SINGLE*uncaptured_volume_START) then
              uncaptured_volume(layer_iter)=zero
             endif

             ! V^{uncapt,k}=V+V^{uncapt,k+1}
             ! V^{uncapt,k}x^{uncapt,k}=V x+V^{uncapt,k+1}x^{uncapt,k+1}

             do dir=1,sdim
              if (uncaptured_volume(layer_iter).le.zero) then
               uncaptured_centroid(layer_iter,dir)=zero
              else
               uncaptured_centroid(layer_iter,dir)= &
                (uncaptured_volume_save*uncaptured_centroid(layer_iter,dir)- &
                 voltemp*centemp(dir))/uncaptured_volume(layer_iter)
              endif
             enddo ! dir=1..sdim
  
             !uncaptured_volume_fraction(layer_iter)=1.0 initially. 
             uncaptured_volume_fraction(layer_iter)= &
              uncaptured_volume_fraction(layer_iter)- &
              mofdatalocal(vofcomp)
             if (uncaptured_volume_fraction(layer_iter).lt. &
                 one-vfrac_sum_local(layer_iter)+EPS_SINGLE) then
              uncaptured_volume_fraction(layer_iter)= &
                one-vfrac_sum_local(layer_iter)
             endif

             num_processed(layer_iter)=num_processed(layer_iter)+1
             num_processed_total=num_processed_total+1
             material_used(critical_material)=num_processed_total

            else if (critical_material.eq.0) then
             ! do nothing
            else
             print *,"critical_material invalid 18455: ",critical_material
             stop
            endif
  
           else
            print *,"single_material or remaining_vfrac invalid ", &
             single_material,remaining_vfrac
            stop
           endif

           loop_counter=loop_counter+1
          enddo !while continue_material_processing==1

         else if ((at_least_one(layer_iter).eq.0).or. &
                  (vfrac_sum_local(layer_iter).eq.zero)) then
          !do nothing
         else
          print *,"at_least_one or vfrac_sum_local invalid ", &
           at_least_one,vfrac_sum_local
          stop
         endif
         
        enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

        call sanity_check_multi_get_volume( & !multi_get_volume_grid
         tessellate_source, &
         tessellate_dest, &
         multi_volume, &
         uncaptured_volume_START, &
         uncaptured_volume, &
         volcell, &
         EPS2, &
         EPS_SINGLE, &
         caller_id, &
         mofdata, &
         uncaptured_volume_fraction, &
         vfrac_sum_local, &
         xsten0, &
         nhalf0, &
         xsten_grid, &
         nhalf_grid, &
         is_rigid_local, &
         is_elastic_local, &
         sdim)

       else
        print *,"uncaptured_volume_START invalid ",uncaptured_volume_START
        stop
       endif

      else
       print *,"return_raster_info invalid ",return_raster_info
       stop
      endif

      return
      end subroutine multi_get_volume_grid

        !tessellate_dest=TESSELLATE_FLUIDS for advecting fluids
        !tessellate_dest=TESSELLATE_FLUIDS_ELASTIC for advecting elastics
      subroutine multi_get_volume_grid_and_map( &
       tessellate_source, & !TESSELLATE_FLUIDS
       tessellate_dest, & !TESSELLATE_FLUIDS|FLUIDS_ELASTIC
       tid_in, &
       normdir, &
       coeff, &
       bfact,dx, &
       xsten0,nhalf0, & ! phi = n dot (x-x0) + intercept
       mofdata, &
       xsten_grid,nhalf_grid, & ! find volumes within xsten_grid or,
       multi_volume, &
       multi_cen, &
       multi_volume_map, &
       multi_cen_map, &
       xtetlist, &
       nlist_alloc, &
       nmax, &
       sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, PARAMETER :: caller_id_local=999
      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: nmax
      integer, INTENT(in) :: tessellate_source
      integer, INTENT(in) :: tessellate_dest
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: nhalf0,nhalf_grid
      integer, INTENT(in) :: normdir
      real(amrex_real), INTENT(in) :: coeff(2)
      real(amrex_real), INTENT(in) :: mofdata(num_materials*(2*sdim+3))
      real(amrex_real) mofdatalocal(num_materials*(2*sdim+3))
      real(amrex_real) mofdatasave(num_materials*(2*sdim+3))
      real(amrex_real) mofdatavalid(num_materials*(2*sdim+3))
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      real(amrex_real), INTENT(out) :: multi_volume(num_materials)
      real(amrex_real), INTENT(out) :: multi_cen(sdim,num_materials)
      real(amrex_real), INTENT(out) :: multi_volume_map(num_materials)
      real(amrex_real), INTENT(out) :: multi_cen_map(sdim,num_materials)
      integer dir
      integer vofcomp
      integer im
      integer layer_iter

      real(amrex_real) uncaptured_volume_START
      real(amrex_real) uncaptured_volume_map_START

      real(amrex_real) &
          uncaptured_volume(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) &
          uncaptured_volume_map(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) &
          uncaptured_centroid(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX,sdim)
      real(amrex_real) &
          uncaptured_centroid_map(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX,sdim)
      real(amrex_real) uncaptured_centroid_local(sdim)
      real(amrex_real) uncaptured_centroid_local_map(sdim)
      real(amrex_real) uncaptured_volume_fraction( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer at_least_one( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) volcell
      real(amrex_real) cencell(sdim)
      real(amrex_real) volcut,cencut(sdim)
      integer nlist
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      integer critical_material
      real(amrex_real) nrecon(sdim)
      real(amrex_real) intercept
      real(amrex_real) voltemp,centemp(sdim)
      real(amrex_real) voltemp_map,centemp_map(sdim)
      integer single_material
      real(amrex_real) remaining_vfrac
      real(amrex_real) uncaptured_volume_save
      real(amrex_real) uncaptured_volume_save_map

      integer material_used(num_materials)

      integer fastflag
      integer layer_flag

      integer local_num_materials( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) vfrac_sum_local( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      real(amrex_real) vfrac_sum( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer num_processed( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer num_processed_total

      integer sanity_num_processed_total
      integer layer_sanity

      integer loop_counter
      integer is_masked(num_materials)
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid (multi_get_volume_grid_and_map) ",tid_in
       stop
      endif 

      if (tessellate_source.eq.TESSELLATE_FLUIDS) then
       if ((tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
           (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC)) then
        ! do nothing
       else
        print *,"tessellate_dest invalid: ",tessellate_dest
        stop
       endif
      else
       print *,"tessellate_source invalid10: ",tessellate_source
       stop
      endif

      call init_local_material_vars( &!subroutine multi_get_volume_grid_and_map
         is_rigid_local, &
         is_elastic_local, &
         tessellate_source, &
         tessellate_dest, &
         is_proper_layer_local)

      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon.ne.2*sdim+3: ",ngeom_recon
       stop
      endif

      if (nmax.lt.4) then
       print *,"nmax invalid multi_get_volume_grid_and_map nmax=",nmax
       stop
      endif
      if ((nhalf0.lt.1).or.(nhalf_grid.lt.1)) then
       print *,"nhalf invalid multi get volume grid and map"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid135 multi_get_volume_grid_and_map ",bfact
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_volume_grid_and_map"
       stop
      endif
      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid multi get volume grid and map"
       stop
      endif

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid: ",nlist_alloc
       stop
      endif

      if ((normdir.lt.0).or.(normdir.ge.sdim)) then
       print *,"normdir invalid multi_get_volume_grid_and_map ",normdir
       stop
      endif
      if (coeff(1).gt.zero) then
       !do nothing
      else
       print *,"coeff(1) invalid multi_get_volume_grid_and_map ",coeff(1)
       stop
      endif
 
      do im=1,num_materials
       multi_volume(im)=zero
       multi_volume_map(im)=zero
       do dir=1,sdim
        multi_cen(dir,im)=zero
        multi_cen_map(dir,im)=zero
       enddo
      enddo

       ! sum Frigid <=1
       ! sum Ffluid = 1
       ! multi_get_volume_grid
      call make_vfrac_sum_ok_copy( &
        xsten0,nhalf0, &
        bfact,dx, &
        tessellate_source, & !TESSELLATE_FLUIDS
        mofdata, &
        mofdatavalid,sdim)

      do dir=1,num_materials*ngeom_recon
       mofdatalocal(dir)=mofdatavalid(dir)
       mofdatasave(dir)=mofdatavalid(dir)
      enddo

      call Box_volumeFAST(bfact,dx,xsten0,nhalf0,volcell,cencell,sdim)

      call Box_volumeFAST_and_map( &
        normdir, &
        coeff, &
        bfact, &
        dx, &
        xsten_grid, &
        nhalf_grid, &
        uncaptured_volume_START, &
        uncaptured_centroid_local, &
        uncaptured_volume_map_START, &
        uncaptured_centroid_local_map, &
        sdim)

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       uncaptured_volume(layer_iter)=uncaptured_volume_START
       do dir=1,sdim
        uncaptured_centroid(layer_iter,dir)= &
            uncaptured_centroid_local(dir)
       enddo !dir=1,sdim
       if (uncaptured_volume(layer_iter).ge.zero) then
        !do nothing
       else
        print *,"uncaptured_volume(layer_iter) invalid ", &
         layer_iter,uncaptured_volume(layer_iter)
        stop
       endif

       uncaptured_volume_map(layer_iter)=uncaptured_volume_map_START
       do dir=1,sdim
        uncaptured_centroid_map(layer_iter,dir)= &
            uncaptured_centroid_local_map(dir)
       enddo !dir=1,sdim
       if (uncaptured_volume_map(layer_iter).ge.zero) then
        !do nothing
       else
        print *,"uncaptured_volume_map(layer_iter) invalid ", &
         layer_iter,uncaptured_volume_map(layer_iter)
        stop
       endif

      enddo ! layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      if (volcell.gt.zero) then
       !do nothing
      else
       print *,"volcell invalid multi_getvolume_grid_and_map : ",volcell
       stop
      endif

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       vfrac_sum(layer_iter)=zero
       vfrac_sum_local(layer_iter)=zero
       local_num_materials(layer_iter)=0
      enddo ! layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1

       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

        if (is_proper_layer(im,layer_iter).eq.1) then
         vfrac_sum(layer_iter)= &
            vfrac_sum(layer_iter)+mofdatasave(vofcomp)
        endif

        if (is_proper_layer_local(im,layer_iter).eq.1) then

         vfrac_sum_local(layer_iter)= &
            vfrac_sum_local(layer_iter)+mofdatasave(vofcomp)

         local_num_materials(layer_iter)=local_num_materials(layer_iter)+1

        else if (is_proper_layer_local(im,layer_iter).eq.0) then
         !do nothing
        else
         print *,"is_proper_layer_local(im,layer_iter) invalid"
         stop
        endif
       enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
      enddo ! im=1..num_materials

      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       if ((vfrac_sum(layer_iter).ge.zero).and. &
           (vfrac_sum(layer_iter).le.one+VOFTOL_MATERIAL)) then
        !do nothing
       else
        print *,"vfrac_sum invalid: ",vfrac_sum
        print *,"layer_iter= ",layer_iter
        stop
       endif
       !tessellate_source=TESSELLATE_FLUIDS|IGNORE_ISELASTIC|IGNORE_ISRIGID
       if ((vfrac_sum_local(layer_iter).ge.zero).and. &
           (vfrac_sum_local(layer_iter).le.one+VOFTOL_MATERIAL)) then
        !do nothing
       else
        print *,"vfrac_sum_local invalid: ",vfrac_sum
        print *,"layer_iter= ",layer_iter
        stop
       endif
      enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      if (local_num_materials(FLUID_LAYER_INDEX)+ &
          local_num_materials(ELASTIC_LAYER_INDEX)+ &
          local_num_materials(RIGID_LAYER_INDEX).ne.num_materials) then
       print *,"local_num_materials invalid ",local_num_materials
       stop
      endif

      if (abs(one-vfrac_sum_local(FLUID_LAYER_INDEX)).le.EPS1) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif

      if ((vfrac_sum_local(FLUID_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
          (vfrac_sum_local(FLUID_LAYER_INDEX).ge.zero)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif
      if ((vfrac_sum_local(ELASTIC_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
          (vfrac_sum_local(ELASTIC_LAYER_INDEX).ge.zero)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif
      if ((vfrac_sum_local(RIGID_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
          (vfrac_sum_local(RIGID_LAYER_INDEX).ge.zero)) then
       ! do nothing
      else
       print *,"vfrac_sum_local invalid: ",vfrac_sum_local
       stop
      endif

      num_processed_total=0

       !multi_get_volume_grid_and_map
      do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       uncaptured_volume_fraction(layer_iter)=one
       num_processed(layer_iter)=0
      enddo

      do im=1,num_materials
       material_used(im)=0
      enddo ! im=1..num_materials

      fastflag=1

       ! intersection of departure region with grid cell is 
       ! very small.
      if (uncaptured_volume_START.le.EPS_12_6*volcell) then

       call tiny_fraction_case( &
          mofdatasave, &
          uncaptured_volume_START, &
          uncaptured_centroid, &
          multi_volume, &
          multi_cen, &
          is_rigid_local, &
          is_elastic_local, &
          tessellate_source, &
          tessellate_dest, &
          vfrac_sum_local, &
          sdim)

       call tiny_fraction_case( &
          mofdatasave, &
          uncaptured_volume_map_START, &
          uncaptured_centroid_map, &
          multi_volume_map, &
          multi_cen_map, &
          is_rigid_local, &
          is_elastic_local, &
          tessellate_source, &
          tessellate_dest, &
          vfrac_sum_local, &
          sdim)

      else if (uncaptured_volume_START.ge.EPS_12_6*volcell) then

       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

        at_least_one(layer_iter)=0
        do im=1,num_materials
         if (is_proper_layer_local(im,layer_iter).eq.1) then
          at_least_one(layer_iter)=at_least_one(layer_iter)+1
         endif
        enddo

        call init_layer_flag( &
          at_least_one(layer_iter), &
          vfrac_sum, &
          vfrac_sum_local, &
          layer_flag, & !intent(out)
          layer_iter, & !intent(in)
          tessellate_source, &
          tessellate_dest)

        call init_mask_flag( &
           layer_flag, &
           is_masked)

        call advance_uncaptured_vars( & !multi_get_volume_grid_and_map
          at_least_one(layer_iter), &
          tessellate_source, &
          tessellate_dest, &
          layer_iter, &
          uncaptured_volume, &
          uncaptured_centroid, &
          uncaptured_volume_fraction, &
          sdim)
        call advance_uncaptured_vars( & !multi_get_volume_grid_and_map
          at_least_one(layer_iter), &
          tessellate_source, &
          tessellate_dest, &
          layer_iter, &
          uncaptured_volume_map, &
          uncaptured_centroid_map, &
          uncaptured_volume_fraction, &
          sdim)


        if ((at_least_one(layer_iter).ge.1).and. &
            (at_least_one(layer_iter).le.num_materials).and. &
            (vfrac_sum_local(layer_iter).gt.zero)) then

          !uncaptured_volume(layer_iter) initialized to the box volume
          !for all "layer_iter"
          ! 
          !uncaptured_volume_fraction(layer_iter) initialized to 1.0
          !for all "layer_iter"
         loop_counter=0
         do while (continue_material_processing( &
               loop_counter, &
               local_num_materials(layer_iter), &
               num_processed(layer_iter), &
               uncaptured_volume_fraction(layer_iter), &
               vfrac_sum_local(layer_iter), &
               uncaptured_volume(layer_iter)).eq.1) 

          call check_for_single_material( & !multi_get_volume_grid_and_map
            single_material, &
            remaining_vfrac, &
            num_processed_total, &
            material_used, &
            is_masked, &
            mofdatasave, &
            mofdatalocal, &
            EPS_UNCAPTURED, &
            uncaptured_volume_fraction, &
            layer_iter, &
            tessellate_source, &
            tessellate_dest, &
            sdim)

          sanity_num_processed_total=0
          do layer_sanity=RIGID_LAYER_INDEX,layer_iter
           sanity_num_processed_total=sanity_num_processed_total+ &
             num_processed(layer_sanity)
          enddo
          if (sanity_num_processed_total.eq.num_processed_total) then
           !do nothing
          else
           print *,"num_processed_total invalid ",num_processed_total
           print *,"sanity_num_processed_total ",sanity_num_processed_total
           stop
          endif

          if ((single_material.ge.1).and. &
              (single_material.le.num_materials).and. &
              (remaining_vfrac.le.EPS_UNCAPTURED)) then

           vofcomp=(single_material-1)*ngeom_recon+1
           multi_volume(single_material)=uncaptured_volume(layer_iter)
           multi_volume_map(single_material)=uncaptured_volume_map(layer_iter)
           do dir=1,sdim
            multi_cen(dir,single_material)= &
             uncaptured_centroid(layer_iter,dir)
            multi_cen_map(dir,single_material)= &
             uncaptured_centroid_map(layer_iter,dir)
           enddo

           uncaptured_volume(layer_iter)=zero
           uncaptured_volume_map(layer_iter)=zero
           uncaptured_volume_fraction(layer_iter)=zero

           num_processed(layer_iter)=num_processed(layer_iter)+1
           num_processed_total=num_processed_total+1
           material_used(single_material)=num_processed_total

          else if ((single_material.eq.0).or. &
                   (remaining_vfrac.ge.EPS_UNCAPTURED)) then

           if ((num_processed_total.gt.0).and. &
               (num_processed_total.lt.num_materials)) then
            fastflag=0
           else if (num_processed_total.eq.0) then
            fastflag=1
           else          
            print *,"num_processed_total invalid: ",num_processed_total
            stop
           endif

           if (fastflag.eq.0) then

            ! only xsten0(0,dir) dir=1..sdim used
            ! in: multi_volume_grid
            call tets_box_planes( &
              layer_flag, &
              bfact,dx, &
              xsten0,nhalf0, &
              xsten_grid,nhalf_grid, &
              mofdatalocal, & !only slope,intercept,order used.
              xtetlist, &
              nlist_alloc, &
              nlist, &
              nmax, &
              sdim)

            call get_cut_geom3D(xtetlist, &
              nlist_alloc,nlist,nmax, &
              volcut,cencut,sdim)

            if (abs(volcut-uncaptured_volume(layer_iter)).gt. &
                EPS1*volcell) then
             print *,"volcut invalid multi_get_volume_grid 1"
             print *,"volcut= ",volcut
             print *,"uncaptured_volume=",uncaptured_volume
             print *,"volcell= ",volcell
             print *,"VOFTOL_MATERIAL= ",VOFTOL_MATERIAL
             print *,"EPS_11_4= ",EPS_11_4
             print *,"EPS_12_6= ",EPS_12_6
             print *,"EPS_8_4= ",EPS_8_4
             print *,"xsten0 ",xsten0(0,1),xsten0(0,2),xsten0(0,sdim)
             print *,"xsten_grid ",xsten_grid(0,1),xsten_grid(0,2), &
               xsten_grid(0,sdim)
             do im=1,num_materials
              vofcomp=(im-1)*ngeom_recon+1
              print *,"im,mofdatavalid(vofcomp) ",im,mofdatavalid(vofcomp)
             enddo 
             print *,"is_rigid_local=",is_rigid_local
             print *,"is_elastic_local=",is_elastic_local
             print *,"is_masked=",is_masked
             stop
            endif

           else if (fastflag.eq.1) then

            ! do nothing; unnecessary to intersect the original box with
            ! the compliment of materials already processed.

           else 
            print *,"fastflag invalid multi get volume grid map: ",fastflag
            stop
           endif

           call derive_critical_material( &
             mofdatalocal, &
             mofdatasave, &
             material_used, &
             is_masked, &
             critical_material, &
             sdim)

           if ((critical_material.ge.1).and. &
               (critical_material.le.num_materials)) then        
            vofcomp=(critical_material-1)*ngeom_recon+1
            do dir=1,sdim
             nrecon(dir)=mofdatalocal(vofcomp+sdim+1+dir)
            enddo
            intercept=mofdatalocal(vofcomp+2*sdim+2)

            if (fastflag.eq.0) then
             !only xsten0(0,dir) dir=1..sdim used
             !calling from "multi_get_volume_grid_and_map"
             call multi_cell_intersection_and_map( &
              normdir, &
              coeff, &
              bfact,dx, &
              xsten0,nhalf0, &
              nrecon,intercept, &
              voltemp,centemp, &
              voltemp_map,centemp_map, &
              xtetlist, &
              nlist_alloc, &
              nlist, &
              nmax, &
              sdim) 
            else if (fastflag.eq.1) then
             ! only xsten0(0,dir) dir=1..sdim used
             call fast_cut_cell_intersection_and_map( &
              normdir, &
              coeff, &
              bfact,dx, &
              xsten0,nhalf0, &
              nrecon,intercept, &
              voltemp,centemp, &
              voltemp_map,centemp_map, &
              xsten_grid,nhalf_grid, &
              sdim) 
            else 
             print *,"fastflag invalid multi get volume grid and map 2"
             stop
            endif

            multi_volume(critical_material)=voltemp
            multi_volume_map(critical_material)=voltemp_map
            do dir=1,sdim
             if (voltemp.gt.zero) then
              multi_cen(dir,critical_material)=centemp(dir)
              multi_cen_map(dir,critical_material)=centemp_map(dir)
             else
              multi_cen(dir,critical_material)=zero
              multi_cen_map(dir,critical_material)=zero
             endif
            enddo ! dir=1..sdim

            uncaptured_volume_save=uncaptured_volume(layer_iter)
            uncaptured_volume_save_map=uncaptured_volume_map(layer_iter)
            uncaptured_volume(layer_iter)= &
               uncaptured_volume(layer_iter)- &
               voltemp
            uncaptured_volume_map(layer_iter)= &
               uncaptured_volume_map(layer_iter)- &
               voltemp_map
            if ((uncaptured_volume(layer_iter).lt. &
                 EPS_UNCAPTURED*uncaptured_volume_START).or. &
                (uncaptured_volume_map(layer_iter).lt. &
                 EPS_UNCAPTURED*uncaptured_volume_map_START)) then
             uncaptured_volume(layer_iter)=zero
             uncaptured_volume_map(layer_iter)=zero
            endif

            ! V^{uncapt,k}=V+V^{uncapt,k+1}
            ! V^{uncapt,k}x^{uncapt,k}=V x+V^{uncapt,k+1}x^{uncapt,k+1}

            do dir=1,sdim
             if ((uncaptured_volume(layer_iter).le.zero).or. &
                 (uncaptured_volume_map(layer_iter).le.zero)) then
              uncaptured_centroid(layer_iter,dir)=zero
              uncaptured_centroid_map(layer_iter,dir)=zero
             else
              uncaptured_centroid(layer_iter,dir)= &
               (uncaptured_volume_save*uncaptured_centroid(layer_iter,dir)- &
                voltemp*centemp(dir))/uncaptured_volume(layer_iter)
              uncaptured_centroid_map(layer_iter,dir)= &
               (uncaptured_volume_save_map* &
                uncaptured_centroid_map(layer_iter,dir)- &
                voltemp_map*centemp_map(dir))/ &
               uncaptured_volume_map(layer_iter)
             endif
            enddo ! dir=1..sdim
 
            !uncaptured_volume_fraction(layer_iter)=1.0 initially. 
            uncaptured_volume_fraction(layer_iter)= &
             uncaptured_volume_fraction(layer_iter)- &
             mofdatalocal(vofcomp)
            if (uncaptured_volume_fraction(layer_iter).lt. &
                one-vfrac_sum_local(layer_iter)+EPS_UNCAPTURED) then
             uncaptured_volume_fraction(layer_iter)= &
               one-vfrac_sum_local(layer_iter)
            endif

            num_processed(layer_iter)=num_processed(layer_iter)+1
            num_processed_total=num_processed_total+1
            material_used(critical_material)=num_processed_total

           else if (critical_material.eq.0) then
            ! do nothing
           else
            print *,"critical_material invalid 18455: ",critical_material
            stop
           endif
 
          else
           print *,"single_material or remaining_vfrac invalid ", &
            single_material,remaining_vfrac
           stop
          endif

          loop_counter=loop_counter+1
         enddo !while continue_material_processing==1

        else if ((at_least_one(layer_iter).eq.0).or. &
                 (vfrac_sum_local(layer_iter).eq.zero)) then
         !do nothing
        else
         print *,"at_least_one or vfrac_sum_local invalid ", &
          at_least_one,vfrac_sum_local
         stop
        endif
        
       enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       call sanity_check_multi_get_volume( & !multi_get_volume_grid_and_map
         tessellate_source, &
         tessellate_dest, &
         multi_volume, &
         uncaptured_volume_START, &
         uncaptured_volume, &
         volcell, &
         EPS2, &
         EPS2, &
         caller_id_local, &
         mofdata, &
         uncaptured_volume_fraction, &
         vfrac_sum_local, &
         xsten0, &
         nhalf0, &
         xsten_grid, &
         nhalf_grid, &
         is_rigid_local, &
         is_elastic_local, &
         sdim)

      else
       print *,"uncaptured_volume_START invalid ",uncaptured_volume_START
       stop
      endif

      return
      end subroutine multi_get_volume_grid_and_map


      subroutine multi_get_volume_tetlist( &
       tid_in, &
       EPS_SINGLE, &
       bfact,dx, &
       xsten0,nhalf0, & ! phi = n dot (x-x0) + intercept
       mofdata, &
       xtetlist_in, & 
       nlist_alloc_in, &
       nlist_in, &
       multi_volume, &
       multi_cen, &
       multi_area, &
       xtetlist, &
       nlist_alloc, &
       nmax, &
       sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: nlist_alloc_in
      integer, INTENT(in) :: nlist_in
      integer, INTENT(in) :: nmax
      integer, parameter :: tessellate=TESSELLATE_IGNORE_ISRIGID
      integer, INTENT(in) :: sdim,bfact
      integer, INTENT(in) :: nhalf0
      real(amrex_real), INTENT(in) :: EPS_SINGLE
      real(amrex_real), INTENT(in) :: xtetlist_in(4,3,nlist_alloc_in)
      real(amrex_real) :: xtet(sdim+1,sdim)
      real(amrex_real), INTENT(in) :: mofdata(num_materials*(2*sdim+3))
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(out) :: multi_volume(num_materials)
      real(amrex_real), INTENT(out) :: multi_cen(sdim,num_materials)
      real(amrex_real), INTENT(out) :: multi_area(num_materials)
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)
      integer dir_local
      integer im
      integer, PARAMETER :: shapeflag=1 !tetrahedron
      real(amrex_real) :: multi_volume_sub(num_materials)
      real(amrex_real) :: multi_cen_sub(sdim,num_materials)
      real(amrex_real) :: multi_area_sub(num_materials)
      integer i_list
      integer itet_node
      integer, parameter :: caller_id=1

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid (multi_get_volume_tetlist) ",tid_in
       stop
      endif 
      if ((EPS_SINGLE.gt.zero).and. &
          (EPS_SINGLE.lt.half)) then
       !do nothing
      else
       print *,"EPS_SINGLE invalid: ",EPS_SINGLE
       stop
      endif
      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon.ne.2*sdim+3"
       stop
      endif
      if (nmax.lt.4) then
       print *,"nmax invalid multi_get_volume_tetlist nmax=",nmax
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid multi get volume tetlist"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid135"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_volume_tetlist"
       stop
      endif
      if ((num_materials.lt.1).or.(num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid multi get volume tetlist"
       stop
      endif

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif
      if ((nlist_alloc_in.ge.1).and.(nlist_alloc_in.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif
 
      do im=1,num_materials
       multi_volume(im)=zero
       do dir_local=1,sdim
        multi_cen(dir_local,im)=zero
       enddo
       multi_area(im)=zero
      enddo

      if (nlist_in.gt.0) then
       do i_list=1,nlist_in
        do itet_node=1,sdim+1
         do dir_local=1,sdim
          xtet(itet_node,dir_local)=xtetlist_in(itet_node,dir_local,i_list)
         enddo
        enddo
         ! multi_cen_sub is "absolute" (not relative to cell centroid)
        call multi_get_volume_grid( &
          caller_id, &
          tid_in, &
          EPS_SINGLE, &
          tessellate, &  !TESSELLATE_IGNORE_ISRIGID
          tessellate, &  !TESSELLATE_IGNORE_ISRIGID
          bfact,dx, &
          xsten0,nhalf0, &
          mofdata, &
          xsten0,nhalf0, &
          xtet, &
          multi_volume_sub, &
          multi_cen_sub, &
          multi_area_sub, &
          xtetlist, &
          nlist_alloc, &
          nmax, &
          sdim, &
          shapeflag)
        do im=1,num_materials
         multi_volume(im)=multi_volume(im)+multi_volume_sub(im)
         multi_area(im)=multi_area(im)+multi_area_sub(im)
         do dir_local=1,sdim
          multi_cen(dir_local,im)=multi_cen(dir_local,im)+ &
              multi_cen_sub(dir_local,im)*multi_volume_sub(im)
         enddo
        enddo ! im=1..num_materials
       enddo ! i_list=1..nlist_in

       do im=1,num_materials
        if (multi_volume(im).eq.0) then
         ! do nothing
        else if (multi_volume(im).gt.zero) then
         do dir_local=1,sdim
          multi_cen(dir_local,im)=multi_cen(dir_local,im)/multi_volume(im)
         enddo
        else
         print *,"multi_volume invalid: ",multi_volume
         stop
        endif
       enddo ! im=1..num_materials
         
      else if (nlist_in.eq.0) then
       ! do nothing
      else
       print *,"nlist_in invalid"
       stop
      endif

      return
      end subroutine multi_get_volume_tetlist

      subroutine project_centroid_box(cen_in,xsten,nhalf,sdim)
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,sdim)
      real(amrex_real), INTENT(inout) :: cen_in(sdim)

      integer dir_local

      do dir_local=1,sdim
       if (cen_in(dir_local).lt.xsten(-1,dir_local)) then
        cen_in(dir_local)=xsten(-1,dir_local)
       else if (cen_in(dir_local).gt.xsten(1,dir_local)) then
        cen_in(dir_local)=xsten(1,dir_local)
       else if ((cen_in(dir_local).ge.xsten(-1,dir_local)).and. &
                (cen_in(dir_local).le.xsten(1,dir_local))) then
        ! do nothing
       else
        print *,"cen_in invalid"
        stop
       endif
      enddo

      return
      end subroutine project_centroid_box

      subroutine multi_get_area_pairs( &
       tid_in, &
       tessellate_source, & !TESSELLATE_FLUIDS
       tessellate_dest, & !TESSELLATE_ALL,TESSELLATE_ALL_RASTER
       bfact,dx, &
       xsten0_plus, &
       xsten0_minus, & !phi = n dot (x-x0) + intercept (phi>0 in omega_m)
       nhalf0, & 
       mofdata_plus, & ! vfrac,centroid(unused),order,slope,intercept
       mofdata_minus, & ! vfrac,centroid(unused),order,slope,intercept
       dir_plus, & ! 1..sdim
       multi_area_pair, & ! (num_materials,num_materials) (left,right)
       sdim, &
       xtetlist_plus, &
       nlist_alloc_plus, &
       xtetlist_minus, &
       nlist_alloc_minus, &
       nmax)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: tessellate_source !TESSELLATE_FLUIDS
      integer, INTENT(in) :: tessellate_dest !TESSELLATE_ALL|ALL_RASTER
      integer, INTENT(in) :: nlist_alloc_plus
      integer, INTENT(in) :: nlist_alloc_minus
      integer, INTENT(in) :: nmax
      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: nhalf0
      integer, INTENT(in) :: dir_plus
      real(amrex_real), INTENT(in) :: mofdata_plus(num_materials*ngeom_recon)
      real(amrex_real), INTENT(in) :: mofdata_minus(num_materials*ngeom_recon)
      real(amrex_real), INTENT(in) :: xsten0_plus(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: xsten0_minus(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(out) :: xtetlist_plus(4,3,nlist_alloc_plus)
      real(amrex_real), INTENT(out) :: xtetlist_minus(4,3,nlist_alloc_minus)
      real(amrex_real), INTENT(out) :: &
        multi_area_pair(num_materials,num_materials)
      real(amrex_real) :: mofdatavalid_plus(num_materials*ngeom_recon)
      real(amrex_real) :: mofdatavalid_minus(num_materials*ngeom_recon)
      real(amrex_real) :: mofdataproject_plus(num_materials*ngeom_recon)
      real(amrex_real) :: mofdataproject_minus(num_materials*ngeom_recon)

      integer im,im_opp,im_test
      integer side
      integer dir_local
      integer, PARAMETER :: nhalf_thin=1
      integer isten 
      real(amrex_real) dxthin
      real(amrex_real) xsten_thin(-nhalf_thin:nhalf_thin,sdim)
      real(amrex_real) xtet(sdim+1,sdim)
      integer, PARAMETER :: shapeflag=0 !regular hexahedron

      real(amrex_real) multi_volume_plus_thin(num_materials)
      real(amrex_real) multi_area_plus_thin(num_materials)
      real(amrex_real) multi_cen_plus_thin(sdim,num_materials)

      real(amrex_real) multi_volume_plus_thin_shrink(num_materials)
      real(amrex_real) multi_area_plus_thin_shrink(num_materials)
      real(amrex_real) multi_cen_plus_thin_shrink(sdim,num_materials)

      real(amrex_real) multi_volume_pair(num_materials,num_materials)
      real(amrex_real) multi_volume_cen_pair(num_materials,num_materials,sdim)

      real(amrex_real) &
          uncaptured_volume(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) &
          uncaptured_centroid(RIGID_LAYER_INDEX:FLUID_LAYER_INDEX,sdim)
      real(amrex_real) uncaptured_volume_fraction( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) uncaptured_volume_save
      real(amrex_real) uncaptured_volume_START
      real(amrex_real) uncaptured_centroid_START(sdim)
      real(amrex_real) uncaptured_area

      real(amrex_real) volume_plus
      real(amrex_real) centroid_plus(sdim)

      real(amrex_real) voltemp
      real(amrex_real) areatemp
      real(amrex_real) centemp(sdim)

      integer num_processed_total
      integer num_processed( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      integer material_used(num_materials)

      integer, PARAMETER :: IGNORE_ISRIGID_tessellate=TESSELLATE_IGNORE_ISRIGID

      integer nlist
      integer loop_counter 
      integer single_material 
      integer critical_material 
      integer fastflag
      integer vofcomp
      real(amrex_real) remaining_vfrac
      integer at_least_one( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) vfrac_sum_local( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) vfrac_sum( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      real(amrex_real) volcut
      real(amrex_real) cencut(sdim)
      real(amrex_real) cen_diff(sdim)
      real(amrex_real) vol_old,vol_new,vol_diff

      real(amrex_real) nrecon(sdim)
      real(amrex_real) intercept

      integer layer_iter,layer_flag
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_masked(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)
      integer local_num_materials( &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      integer, parameter :: caller_id=2

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid (multi_get_area_pairs) ",tid_in
       stop
      endif 
      if (ngeom_recon.eq.2*sdim+3) then
       ! do nothing
      else
       print *,"ngeom_recon.ne.2*sdim+3: ",ngeom_recon
       stop
      endif
      if (nmax.lt.4) then
       print *,"nmax invalid multi_get_area pairs nmax=",nmax
       stop
      endif

      if (nhalf0.ge.1) then
       ! do nothing
      else
       print *,"nhalf0 invalid multi get area pairs: ",nhalf0
       stop
      endif

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid multi_get_area_pairs: ",bfact
       stop
      endif
      if ((sdim.eq.3).or.(sdim.eq.2)) then
       ! do nothing
      else
       print *,"sdim invalid multi_get_area_pairs: ",sdim
       stop
      endif
      if ((num_materials.ge.1).and. &
          (num_materials.le.MAX_NUM_MATERIALS)) then
       ! do nothing
      else
       print *,"num_materials invalid multi get area pairs: ",num_materials
       stop
      endif

      if ((dir_plus.ge.1).and.(dir_plus.le.sdim)) then
       ! do nothing
      else
       print *,"dir_plus invalid: ",dir_plus
       stop
      endif 

      if (tessellate_source.eq.TESSELLATE_FLUIDS) then
       !do nothing
      else
       print *,"tessellate_source invalid: ",tessellate_source
       stop
      endif
      
      if ((tessellate_dest.eq.TESSELLATE_ALL).or. &
          (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
       ! do nothing
      else
       print *,"tessellate_dest invalid: ",tessellate_dest
       stop
      endif

       ! left index: mofdata_minus
       ! right index: mofdata_plus
      do im=1,num_materials
       do im_opp=1,num_materials
        multi_area_pair(im,im_opp)=zero
        multi_volume_pair(im,im_opp)=zero
        do dir_local=1,sdim
         multi_volume_cen_pair(im,im_opp,dir_local)=zero
        enddo
       enddo
      enddo  ! im=1..num_materials

      call Box_volumeFAST(bfact,dx, &
         xsten0_plus,nhalf0, &
         volume_plus, &
         centroid_plus,sdim)

      dxthin=EPS_3_2*half* &
          (xsten0_plus(1,dir_plus)-xsten0_minus(-1,dir_plus))
      if (dxthin.gt.zero) then
       ! do nothing
      else
       print *,"dxthin invalid: ",dxthin
       stop
      endif

      do isten=-nhalf_thin,nhalf_thin
       do dir_local=1,sdim
        xsten_thin(isten,dir_local)=xsten0_plus(isten,dir_local)
       enddo
      enddo ! isten
      xsten_thin(1,dir_plus)=xsten0_plus(-1,dir_plus)+dxthin
      xsten_thin(-1,dir_plus)=xsten0_minus(1,dir_plus)-dxthin
      xsten_thin(0,dir_plus)=half*(xsten0_plus(-1,dir_plus)+ &
                                   xsten0_minus(1,dir_plus))

       !multi_get_area_pairs
      call make_vfrac_sum_ok_copy( &
        xsten0_plus,nhalf0, &
        bfact,dx, &
        tessellate_source, &  ! =TESSELLATE_FLUIDS
        mofdata_plus, &
        mofdatavalid_plus, &
        sdim)
       !multi_get_area_pairs
      call make_vfrac_sum_ok_copy( &
        xsten0_minus,nhalf0, &
        bfact,dx, &
        tessellate_source, & ! =TESSELLATE_FLUIDS
        mofdata_minus, &
        mofdatavalid_minus, &
        sdim)

      if ((tessellate_dest.eq.TESSELLATE_ALL).or. &
          (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then

       call multi_get_volume_tessellate( &
        tid_in, &
        tessellate_source, & !TESSELLATE_FLUIDS 
        tessellate_dest, & !=TESSELLATE_ALL or TESSELLATE_ALL_RASTER 
        bfact,dx, &
        xsten0_plus,nhalf0, &
        mofdatavalid_plus, &
        xtetlist_plus, &
        nlist_alloc_plus, &
        nmax, &
        sdim)

       call multi_get_volume_tessellate( &
        tid_in, &
        tessellate_source, & !TESSELLATE_FLUIDS 
        tessellate_dest, & ! =TESSELLATE_ALL or TESSELLATE_ALL_RASTER
        bfact,dx, &
        xsten0_minus,nhalf0, &
        mofdatavalid_minus, &
        xtetlist_minus, &
        nlist_alloc_minus, &
        nmax, &
        sdim)

      else
       print *,"tessellate_dest invalid: ",tessellate_dest
       stop
      endif

       ! only slopes and intercepts modified.
      side=1
      call project_slopes_to_face( &
        bfact,dx,xsten0_plus,nhalf0, &
        mofdatavalid_plus, &
        mofdataproject_plus, &
        sdim,dir_plus,side)
      side=2
      call project_slopes_to_face( &
        bfact,dx,xsten0_minus,nhalf0, &
        mofdatavalid_minus, &
        mofdataproject_minus, &
        sdim,dir_plus,side)

      call multi_get_volume_grid( &
       caller_id, &
       tid_in, &
       EPS_FULL_WEAK, & !tolerance for "single material" criterion
       IGNORE_ISRIGID_tessellate, & !TESSELLATE_IGNORE_ISRIGID
       IGNORE_ISRIGID_tessellate, & !TESSELLATE_IGNORE_ISRIGID
       bfact,dx, &
       xsten0_plus,nhalf0, &
       mofdataproject_plus, &
       xsten_thin,nhalf_thin, &
       xtet, &
       multi_volume_plus_thin, &
       multi_cen_plus_thin, &
       multi_area_plus_thin, &
       xtetlist_plus, &
       nlist_alloc_plus, &
       nmax, &
       sdim, &
       shapeflag)

       !  xmid-dxthin < x < xmid+dxthin
      call Box_volumeFAST(bfact,dx, &
         xsten_thin,nhalf_thin, &
         uncaptured_volume_START, &
         uncaptured_centroid_START,sdim)

       ! can use "subroutine gridarea" instead (declared in GLOBALUTIL.F90).
      uncaptured_area=uncaptured_volume_START/(2.0*dxthin)

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (dir_plus.eq.1) then
        if ((xsten0_minus(0,dir_plus).lt.zero).and. &
            (xsten0_plus(0,dir_plus).gt.zero)) then
         uncaptured_area=zero 
        else if ((xsten0_minus(0,dir_plus).gt.zero).and. &
                 (xsten0_plus(0,dir_plus).gt.zero)) then
         ! do nothing
        else
         print *,"xsten0_minus or xsten0_plus invalid"
         print *,"xsten0_minus ",xsten0_minus
         print *,"xsten0_plus ",xsten0_plus
         stop
        endif
       else if ((dir_plus.eq.2).or.(dir_plus.eq.sdim)) then
        ! do nothing
       else
        print *,"dir_plus invalid: ",dir_plus
        stop
       endif
      else
       print *,"levelrz invalid: ",levelrz
       stop
      endif

      if (uncaptured_volume_START.gt.zero) then
       ! do nothing
      else
       print *,"uncaptured_volume_START invalid: ",uncaptured_volume_START
       stop
      endif
      if (uncaptured_area.ge.zero) then
       ! do nothing
      else
       print *,"uncaptured_area invalid: ",uncaptured_area
       stop
      endif

      if (uncaptured_area.gt.zero) then

       call init_local_material_vars( & !subroutine multi_get_area_pairs
         is_rigid_local, &
         is_elastic_local, &
         IGNORE_ISRIGID_tessellate, &
         IGNORE_ISRIGID_tessellate, &
         is_proper_layer_local)

       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
        uncaptured_volume(layer_iter)=uncaptured_volume_START
        do dir_local=1,sdim
         uncaptured_centroid(layer_iter,dir_local)= &
            uncaptured_centroid_START(dir_local)
        enddo !dir_local=1,sdim
        if (uncaptured_volume(layer_iter).ge.zero) then
         !do nothing
        else
         print *,"uncaptured_volume(layer_iter) invalid ", &
          layer_iter,uncaptured_volume(layer_iter)
         stop
        endif
       enddo ! layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
        vfrac_sum(layer_iter)=zero
        vfrac_sum_local(layer_iter)=zero
        local_num_materials(layer_iter)=0
       enddo ! layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1

        do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

         if (is_proper_layer(im,layer_iter).eq.1) then
          vfrac_sum(layer_iter)= &
            vfrac_sum(layer_iter)+mofdataproject_minus(vofcomp)
         endif
 
         if (is_proper_layer_local(im,layer_iter).eq.1) then

          vfrac_sum_local(layer_iter)= &
            vfrac_sum_local(layer_iter)+mofdataproject_minus(vofcomp)

          local_num_materials(layer_iter)=local_num_materials(layer_iter)+1

         else if (is_proper_layer_local(im,layer_iter).eq.0) then
          !do nothing
         else
          print *,"is_proper_layer_local(im,layer_iter) invalid"
          stop
         endif
        enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
       enddo ! im=1..num_materials

       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
        if ((vfrac_sum(layer_iter).ge.zero).and. &
            (vfrac_sum(layer_iter).le.one+VOFTOL_MATERIAL)) then
         !do nothing
        else
         print *,"vfrac_sum invalid: ",vfrac_sum
         print *,"layer_iter= ",layer_iter
         stop
        endif
        !IGNORE_ISRIGID_tessellate
        if ((vfrac_sum_local(layer_iter).ge.zero).and. &
            (vfrac_sum_local(layer_iter).le.one+VOFTOL_MATERIAL)) then
         !do nothing
        else
         print *,"vfrac_sum_local invalid: ",vfrac_sum
         print *,"layer_iter= ",layer_iter
         stop
        endif
       enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

       if (local_num_materials(FLUID_LAYER_INDEX)+ &
           local_num_materials(ELASTIC_LAYER_INDEX)+ &
           local_num_materials(RIGID_LAYER_INDEX).ne.num_materials) then
        print *,"local_num_materials invalid ",local_num_materials
        stop
       endif

       if (abs(one-vfrac_sum_local(FLUID_LAYER_INDEX)).le.EPS1) then
        ! do nothing
       else
        print *,"vfrac_sum_local invalid: ",vfrac_sum_local
        stop
       endif
       if ((vfrac_sum_local(ELASTIC_LAYER_INDEX).eq.zero).and. &
           (vfrac_sum_local(RIGID_LAYER_INDEX).eq.zero)) then
        !do nothing
       else
        print *,"in: multi_get_area_pairs"
        print *,"expecting vfrac_sum_local(ELASTIC_LAYER_INDEX)=0 ", &
         vfrac_sum_local
        print *,"or expecting vfrac_sum_local(RIGID_LAYER_INDEX)=0 ", &
         vfrac_sum_local
        stop
       endif

       if ((vfrac_sum_local(FLUID_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
           (vfrac_sum_local(FLUID_LAYER_INDEX).ge.zero)) then
        ! do nothing
       else
        print *,"vfrac_sum_local invalid: ",vfrac_sum_local
        stop
       endif
       if ((vfrac_sum_local(ELASTIC_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
           (vfrac_sum_local(ELASTIC_LAYER_INDEX).ge.zero)) then
        ! do nothing
       else
        print *,"vfrac_sum_local invalid: ",vfrac_sum_local
        stop
       endif
       if ((vfrac_sum_local(RIGID_LAYER_INDEX).le.one+VOFTOL_MATERIAL).and. &
           (vfrac_sum_local(RIGID_LAYER_INDEX).ge.zero)) then
        ! do nothing
       else
        print *,"vfrac_sum_local invalid: ",vfrac_sum_local
        stop
       endif

        ! multi_get_area_pairs
       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX
        uncaptured_volume_fraction(layer_iter)=one
        num_processed(layer_iter)=0
       enddo

       do im=1,num_materials
        material_used(im)=0 !i.e. "material_subtracted"
       enddo ! im=1..num_materials

       critical_material=0

       fastflag=1

       do layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

        at_least_one(layer_iter)=0
        do im=1,num_materials
         if (is_proper_layer_local(im,layer_iter).eq.1) then
          at_least_one(layer_iter)=at_least_one(layer_iter)+1
         endif
        enddo


        call init_layer_flag( &
         at_least_one(layer_iter), &
         vfrac_sum, &
         vfrac_sum_local, &
         layer_flag, & !intent(out)
         layer_iter, & !intent(in)
         IGNORE_ISRIGID_tessellate, &
         IGNORE_ISRIGID_tessellate)

        call init_mask_flag( &
           layer_flag, &
           is_masked)

        call advance_uncaptured_vars( & !multi_get_area_pairs
         at_least_one(layer_iter), &
         IGNORE_ISRIGID_tessellate, &
         IGNORE_ISRIGID_tessellate, &
         layer_iter, &
         uncaptured_volume, &
         uncaptured_centroid, &
         uncaptured_volume_fraction, &
         sdim)


        if ((at_least_one(layer_iter).ge.1).and. &
            (at_least_one(layer_iter).le.num_materials).and. &
            (vfrac_sum_local(layer_iter).gt.zero)) then


         loop_counter=0

         do while (continue_material_processing( &
               loop_counter, &
               local_num_materials(layer_iter), &
               num_processed(layer_iter), &
               uncaptured_volume_fraction(layer_iter), &
               vfrac_sum_local(layer_iter), &
               uncaptured_volume(layer_iter)).eq.1) 

          if (layer_iter.eq.FLUID_LAYER_INDEX) then
           !do nothing
          else
           print *,"expecting layer_iter=FLUID_LAYER_INDEX"
           stop
          endif

          call check_for_single_material( & !multi_get_area_pairs
            single_material, & !intent(out)
            remaining_vfrac, & !intent(out)
            num_processed_total, & !intent(out)
            material_used, &
            is_masked, &
            mofdataproject_minus, &
            mofdataproject_minus, &
            EPS_FULL_WEAK, &
            uncaptured_volume_fraction, &
            layer_iter, &
            IGNORE_ISRIGID_tessellate, &
            IGNORE_ISRIGID_tessellate, &
            sdim)


          ! num_processed_total=number of materials subtracted from the
          ! "minus" side.
          if ((num_processed_total.gt.0).and. &
              (num_processed_total.lt.num_materials)) then
           fastflag=0
          else if (num_processed_total.eq.0) then
           fastflag=1
          else          
           print *,"num_processed_total invalid: ",num_processed_total
           stop
          endif

          if (fastflag.eq.0) then ! num_processed_total>=1

           ! Only xsten0(0,dir) dir=1..sdim used.
           ! tets_box_planes intersects xsten_thin with the compliments 
           ! of already initialized (subtracted) materials. 
           ! (i.e. materials with 1<=material_used(im)<=num_materials)

           call tets_box_planes( &
            layer_flag, &
            bfact,dx, &
            xsten0_minus,nhalf0, &
            xsten_thin,nhalf_thin, &
            mofdataproject_minus, &
            xtetlist_minus, &
            nlist_alloc_minus, &
            nlist, &
            nmax, &
            sdim)

           call get_cut_geom3D(xtetlist_minus, &
            nlist_alloc_minus,nlist,nmax, &
            volcut,cencut,sdim)

           if (abs(volcut-uncaptured_volume(layer_iter)).gt. &
               EPS_FULL_WEAK*volume_plus) then
            print *,"volcut invalid multi get area pairs 2 "
            print *,"volcut= ",volcut
            print *,"volume_plus=",volume_plus
            print *,"uncaptured_volume=",uncaptured_volume
            print *,"uncaptured_volume_START= ",uncaptured_volume_START
            if (uncaptured_volume_START.gt.zero) then
             print *,"abs(volcut-uncapt_vol)/uncaptured_volume_START=", &
              abs(volcut-uncaptured_volume(layer_iter))/ &
              uncaptured_volume_START
            endif
            print *,"VOFTOL_MATERIAL= ",VOFTOL_MATERIAL
            print *,"EPS_12_6= ",EPS_12_6
            print *,"EPS_8_4= ",EPS_8_4
            print *,"xsten0_minus ", &
              xsten0_minus(0,1),xsten0_minus(0,2),xsten0_minus(0,sdim)
            stop
           endif

           !"critical_material" is the most recent material to be 
           !subtracted from the minus side.
           if ((critical_material.ge.1).and. &
               (critical_material.le.num_materials)) then

            if ((material_used(critical_material).ge.1).and. &
                (material_used(critical_material).le.num_materials)) then

              !intersection of plus-side data with the 
              !minus side uncaptured space.
             call multi_get_volume_tetlist( &
              tid_in, &
              EPS_FULL_WEAK, & !tolerance for "single material" criterion
              bfact,dx, &
              xsten0_plus,nhalf0, &
              mofdataproject_plus, &
              xtetlist_minus, &
              nlist_alloc_minus, &
              nlist, &
              multi_volume_plus_thin_shrink, &
              multi_cen_plus_thin_shrink, &
              multi_area_plus_thin_shrink, &
              xtetlist_plus, &
              nlist_alloc_plus, &
              nmax, &
              sdim)

             do im_opp=1,num_materials

              vol_old=multi_volume_plus_thin(im_opp)
              vol_new=multi_volume_plus_thin_shrink(im_opp)
              if (vol_old-vol_new.ge.-EPS_FULL_WEAK*volume_plus) then
               if (vol_old-vol_new.le.zero) then
                vol_diff=zero
               else if (vol_old-vol_new.ge.zero) then
                vol_diff=vol_old-vol_new
               else
                print *,"vol_diff bust: ",vol_diff
                print *,"vol_old, vol_new ",vol_old,vol_new
                stop
               endif
               multi_volume_pair(critical_material,im_opp)=vol_diff
               if (vol_diff.gt.zero) then
                do dir_local=1,sdim
                 cen_diff(dir_local)= &
                  (multi_cen_plus_thin(dir_local,im_opp)*vol_old- &
                   multi_cen_plus_thin_shrink(dir_local,im_opp)*vol_new)/ &
                  vol_diff
                enddo

                call project_centroid_box(cen_diff,xsten_thin,nhalf_thin,sdim)

                multi_volume_plus_thin(im_opp)=vol_new

                do dir_local=1,sdim
                 multi_volume_cen_pair(critical_material,im_opp,dir_local)= &
                         cen_diff(dir_local)
                 multi_cen_plus_thin(dir_local,im_opp)= &
                     multi_cen_plus_thin_shrink(dir_local,im_opp)
                enddo
               else if (vol_diff.eq.zero) then
                ! do nothing
               else
                print *,"vol_diff invalid:",vol_diff
                stop
               endif
              else
               print *,"vol_old or vol_new invalid"
               print *,"dir_plus ",dir_plus
               print *,"dxthin=",dxthin
               print *,"uncaptured_volume_START=",uncaptured_volume_START
               print *,"uncaptured_area=",uncaptured_area
               print *,"volume_plus=",volume_plus
               print *,"vol_old=",vol_old
               print *,"vol_new=",vol_new
               print *,"critical_material=",critical_material
               print *,"im_opp=",im_opp
               print *,"loop_counter=",loop_counter
               print *,"num_processed=",num_processed
               print *,"num_materials=",num_materials
               print *,"uncaptured_volume_fraction=",uncaptured_volume_fraction
               print *,"uncaptured_volume=",uncaptured_volume
               do im_test=1,num_materials*ngeom_recon
                print *,"im_test,mofdata_plus,mofdata_minus ", &
                  im_test,mofdata_plus(im_test),mofdata_minus(im_test)
               enddo
               do im_test=-nhalf0,nhalf0
               do dir_local=1,sdim
                print *,"ii,dir,xsten0_plus ",im_test,dir_local, &
                        xsten0_plus(im_test,dir_local)
                print *,"ii,dir,xsten0_minus ",im_test,dir_local, &
                        xsten0_minus(im_test,dir_local)
               enddo
               enddo
               stop
              endif
             enddo ! im_opp=1..num_materials

            else 
             print *,"expecting 1<=material_used(critical_material)<=nmat ", &
               material_used(critical_material)
             stop
            endif
           else
            print *,"warning critical_material invalid 19694:",critical_material
            do im_test=1,num_materials
             print *,"im_test,material_used ",im_test,material_used(im_test)
            enddo
            print *,"tessellate_source=",tessellate_source
            print *,"tessellate_dest=",tessellate_dest
            do im_test=1,num_materials*ngeom_recon
             print *,"i,mofdatavalid_plus ",im_test,mofdatavalid_plus(im_test)
            enddo
            do im_test=1,num_materials*ngeom_recon
             print *,"i,mofdatavalid_minus ",im_test,mofdatavalid_minus(im_test)
            enddo
            do im_test=1,num_materials*ngeom_recon
             print *,"i,mofdataproject_plus ",im_test, &
               mofdataproject_plus(im_test)
            enddo
            do im_test=1,num_materials*ngeom_recon
             print *,"i,mofdataproject_minus ", &
              im_test,mofdataproject_minus(im_test)
            enddo
            print *,"uncaptured_volume: ",uncaptured_volume
            print *,"uncaptured_volume_fraction: ", &
                  uncaptured_volume_fraction
            print *,"loop_counter: ",loop_counter
            print *,"num_processed: ",num_processed
            print *,"dir_plus=",dir_plus
!           stop
           endif

          else if (fastflag.eq.1) then

           if (critical_material.eq.0) then
            ! do nothing
           else
            print *,"expecting critical_material=0 22618: ",critical_material
            stop
           endif

          else 
           print *,"fastflag invalid multi get area pairs:",fastflag
           stop
          endif

          !determine the next ("minus") material to be subtracted.
          critical_material=0

          if ((single_material.ge.1).and. &
              (single_material.le.num_materials).and. &
              (remaining_vfrac.le.EPS_FULL_WEAK)) then

           !we leave critical_material==0 since the "pair" variables
           !are updated right here and only "single_material" is
           !subtracted.
           vofcomp=(single_material-1)*ngeom_recon+1
           do im_opp=1,num_materials
             multi_volume_pair(single_material,im_opp)= &
                  multi_volume_plus_thin(im_opp)
             do dir_local=1,sdim
              multi_volume_cen_pair(single_material,im_opp,dir_local)= &
               multi_cen_plus_thin(dir_local,im_opp)
             enddo
           enddo ! im_opp=1..num_materials

           uncaptured_volume(layer_iter)=zero
           uncaptured_volume_fraction(layer_iter)=zero

           num_processed(layer_iter)=num_processed(layer_iter)+1
           num_processed_total=num_processed_total+1

           if (num_processed(layer_iter).eq.num_processed_total) then
            !do nothing
           else
            print *,"expect num_processed(layer_iter).eq.num_processed_total"
            stop
           endif

           material_used(single_material)=num_processed_total

          else if ((single_material.eq.0).or. &
                   (remaining_vfrac.ge.EPS_FULL_WEAK)) then

           call derive_critical_material( & !multi_get_area_pairs
             mofdataproject_minus, &
             mofdatavalid_minus, &
             material_used, &
             is_masked, &
             critical_material, &
             sdim)

           !"critical_material" is the next material to be subtracted
           !from the "minus" side.
           if ((critical_material.ge.1).and. &
               (critical_material.le.num_materials)) then        

            vofcomp=(critical_material-1)*ngeom_recon+1
            do dir_local=1,sdim
             nrecon(dir_local)=mofdataproject_minus(vofcomp+sdim+1+dir_local)
            enddo
            intercept=mofdataproject_minus(vofcomp+2*sdim+2)

            if (fastflag.eq.0) then
             !calling from "multi_get_area_pairs"
             call multi_cell_intersection( &
              bfact,dx, &
              xsten0_minus,nhalf0, &
              nrecon,intercept, &
              voltemp,centemp, &
              areatemp, &
              xtetlist_minus, &
              nlist_alloc_minus, &
              nlist, &
              nmax, &
              sdim) 
            else if (fastflag.eq.1) then
             call fast_cut_cell_intersection( &
              bfact,dx, &
              xsten0_minus,nhalf0, &
              nrecon,intercept, &
              voltemp,centemp, &
              areatemp, &
              xsten_thin,nhalf_thin, &
              xtet,shapeflag, &
              sdim) 
            else 
             print *,"fastflag invalid multi get area pairs 2:",fastflag
             stop
            endif

            uncaptured_volume_save=uncaptured_volume(layer_iter)
            uncaptured_volume(layer_iter)=uncaptured_volume(layer_iter)- &
                voltemp
            if (uncaptured_volume(layer_iter).lt. &
                EPS_FULL_WEAK*uncaptured_volume_START) then
             uncaptured_volume(layer_iter)=zero
            endif

            uncaptured_volume_fraction(layer_iter)= &
              uncaptured_volume_fraction(layer_iter)- &
              mofdataproject_minus(vofcomp)
            if (uncaptured_volume_fraction(layer_iter).lt.EPS_FULL_WEAK) then
             uncaptured_volume_fraction(layer_iter)=zero
            endif

            num_processed(layer_iter)=num_processed(layer_iter)+1
            num_processed_total=num_processed_total+1
            material_used(critical_material)=num_processed_total

            if (num_processed(layer_iter).eq.num_processed_total) then
             !do nothing
            else
             print *,"expect num_processed(layer_iter).eq.num_processed_total"
             stop
            endif

           else if (critical_material.eq.0) then
            ! do nothing
           else
            print *,"critical_material invalid 19840: ",critical_material
            stop
           endif 

          else
           print *,"single_material or remaining_vfrac invalid"
           print *,"single_material: ",single_material
           print *,"remaining_vfrac: ",remaining_vfrac
           stop
          endif

          loop_counter=loop_counter+1
         enddo !while continue_material_processing==1

         if ((critical_material.ge.1).and. &
             (critical_material.le.num_materials)) then        

          vofcomp=(critical_material-1)*ngeom_recon+1
          do im_opp=1,num_materials

           multi_volume_pair(critical_material,im_opp)= &
                multi_volume_plus_thin(im_opp)
           do dir_local=1,sdim
            multi_volume_cen_pair(critical_material,im_opp,dir_local)= &
             multi_cen_plus_thin(dir_local,im_opp)
           enddo

          enddo ! im_opp=1..num_materials

         else if (critical_material.eq.0) then
          ! do nothing (pair variables already updated)
         else
          print *,"critical_material invalid 19890: ",critical_material
          stop
         endif

         xsten_thin(-1,dir_plus)=xsten_thin(0,dir_plus)
         xsten_thin(1,dir_plus)=xsten_thin(0,dir_plus)
         call project_centroid_box(uncaptured_centroid_START, &
               xsten_thin,nhalf_thin,sdim)

         voltemp=zero

         do im=1,num_materials

          do im_opp=1,num_materials
 
           voltemp=voltemp+multi_volume_pair(im,im_opp)
           do dir_local=1,sdim
            cen_diff(dir_local)=multi_volume_cen_pair(im,im_opp,dir_local)
           enddo
           call project_centroid_box(cen_diff,xsten_thin,nhalf_thin,sdim)

          enddo !im_opp=1..num_materials

         enddo ! im=1..num_materials

         if (voltemp.gt.zero) then
          do im=1,num_materials
           do im_opp=1,num_materials
            if (multi_volume_pair(im,im_opp).le.voltemp) then
             multi_area_pair(im,im_opp)= &
              uncaptured_area*multi_volume_pair(im,im_opp)/voltemp
            else
             print *,"im,im_opp,multi_volume_pair: ", &
               im,im_opp,multi_volume_pair(im,im_opp)
             print *,"warning: voltemp underflow: ",voltemp
!             stop
            endif
           enddo ! im_opp=1..num_materials
          enddo ! im=1..num_materials
         else
          print *,"warning: voltemp invalid: ",voltemp
          print *,"tessellate_source=",tessellate_source
          print *,"tessellate_dest=",tessellate_dest
          do im_test=1,num_materials*ngeom_recon
           print *,"i,mofdatavalid_plus ",im_test,mofdatavalid_plus(im_test)
          enddo
          do im_test=1,num_materials*ngeom_recon
           print *,"i,mofdatavalid_minus ",im_test,mofdatavalid_minus(im_test)
          enddo
          do im_test=1,num_materials*ngeom_recon
           print *,"i,mofdataproject_plus ",im_test,mofdataproject_plus(im_test)
          enddo
          do im_test=1,num_materials*ngeom_recon
           print *,"i,mofdataproject_minus ", &
            im_test,mofdataproject_minus(im_test)
          enddo
!         stop
         endif

         dxthin=half*(xsten0_plus(1,dir_plus)-xsten0_minus(-1,dir_plus))

         if (abs(voltemp-uncaptured_volume_START).le. &
             two*EPS_FULL_WEAK*volume_plus) then
          ! do nothing
         else
          print *,"warning"
          print *,"voltemp or uncaptured_volume_START invalid"
          print *,"voltemp: ",voltemp
          print *,"uncaptured_volume_START: ",uncaptured_volume_START
!         stop
         endif

         if (uncaptured_volume(layer_iter).le. &
             four*EPS_FULL_WEAK*volume_plus) then
          !do nothing
         else
          print *,"warning"
          print *,"not all volume accounted for multi get area pairs"
          print *,"uncaptured_volume ",uncaptured_volume
          print *,"uncaptured_volume_START ",uncaptured_volume_START
          print *,"fraction of uncapt volume ", &
            uncaptured_volume(layer_iter)/uncaptured_volume_START
          print *,"tolerance: ",four*EPS_FULL_WEAK
!         stop
         endif

        else if ((at_least_one(layer_iter).eq.0).or. &
                 (vfrac_sum_local(layer_iter).eq.zero)) then
         !do nothing
        else
         print *,"at_least_one or vfrac_sum_local invalid ", &
          at_least_one,vfrac_sum_local
         stop
        endif

       enddo !layer_iter=RIGID_LAYER_INDEX,FLUID_LAYER_INDEX

      else if (uncaptured_area.eq.zero) then
       ! do nothing
      else
       print *,"uncaptured_area invalid: ",uncaptured_area
       stop
      endif

      return
      end subroutine multi_get_area_pairs



      subroutine FIX_LS_tessellate(LS,LS_new)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: LS(num_materials)
      real(amrex_real), INTENT(out) :: LS_new(num_materials)
      integer im,im_opp,im_tessellate
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)

      do im=1,num_materials
       is_rigid_local(im)=is_rigid(im)
       is_elastic_local(im)=is_elastic(im)
      enddo ! im=1..num_materials

      do im=1,num_materials
       LS_new(im)=LS(im)
      enddo
      do im=1,num_materials
       if ((is_rigid_local(im).eq.0).and. &
           (is_elastic_local(im).eq.0)) then
         ! im_tessellate=argmax_{im_opp<>im} LS(im_opp)
        im_tessellate=0
        do im_opp=1,num_materials
         if ((is_rigid_local(im_opp).eq.0).and. &
             (is_elastic_local(im_opp).eq.0)) then
          if (im_opp.ne.im) then
           if (im_tessellate.eq.0) then
            im_tessellate=im_opp
           else if ((im_tessellate.ge.1).and. &
                    (im_tessellate.le.num_materials)) then
            if (LS(im_tessellate).lt.LS(im_opp)) then
             im_tessellate=im_opp
            endif
           else
            print *,"im_tessellate invalid: ",im_tessellate
            stop
           endif
          else if (im.eq.im_opp) then
           !do nothing
          else
           print *,"im or im_opp invalid: ",im,im_opp
           stop
          endif
         else if ((is_rigid_local(im_opp).eq.1).or. &
                  (is_elastic_local(im_opp).eq.1)) then
          ! do nothing
         else
          print *,"is_rigid_local(im_opp) invalid: ", &
             im_opp,is_rigid_local(im_opp)
          print *,"or is_elastic_local(im_opp) invalid: ", &
             im_opp,is_elastic_local(im_opp)
          stop
         endif
        enddo !im_opp=1..num_materials
        if (im_tessellate.eq.0) then
         if (LS(im).gt.zero) then
          !do nothing
         else
          print *,"im_tessellate invalid: ",im_tessellate
          print *,"im,LS(im) ",im,LS(im)
          stop
         endif
        else if ((im_tessellate.ge.1).and. &
                 (im_tessellate.le.num_materials)) then
         LS_new(im)=(LS(im)-LS(im_tessellate))/two
        else
         print *,"im_tessellate invalid: ",im_tessellate
         stop
        endif
       else if ((is_rigid_local(im).eq.1).or. &
                (is_elastic_local(im).eq.1)) then
        ! do nothing
       else
        print *,"is_rigid_local(im) invalid: ",im,is_rigid_local(im)
        print *,"or is_elastic_local(im) invalid: ",im,is_elastic_local(im)
        stop
       endif
      enddo ! im=1..num_materials

      end subroutine FIX_LS_tessellate

      subroutine LS_tessellate(dx,LS,LS_new,sdim,tessellate_source)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: tessellate_source
      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: LS(num_materials)
      real(amrex_real), INTENT(out) :: LS_new(num_materials)
      real(amrex_real) :: LS_new_hold(num_materials)
      real(amrex_real) LS_primary
      real(amrex_real) LS_max_rest
      integer im,im_primary,im_rest
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      if ((tessellate_source.eq.TESSELLATE_FLUIDS).or. &
          (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC)) then
       !do nothing
      else
       print *,"expecting tessellate_source.eq.TESSELLATE_FLUIDS ", &
         tessellate_source
       print *,"or expect tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC ", &
         tessellate_source
       stop
      endif
      call init_local_material_vars( & !subroutine LS_tessellate
       is_rigid_local, &
       is_elastic_local, &
       tessellate_source, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
       tessellate_source, &
       is_proper_layer_local)

      call get_primary_material(dx,LS,im_primary)

      if ((is_rigid_local(im_primary).eq.0).and. &
          (is_elastic_local(im_primary).eq.0)) then
       do im=1,num_materials
        LS_new(im)=LS(im)
       enddo
      else if (is_rigid_local(im_primary).eq.1) then
       LS_primary=max(LS(im_primary),zero)
       if (LS_primary.ge.0.0d0) then 
        !do nothing
       else
        print *,"LS_primary invalid: ",LS_primary
        stop
       endif
       do im=1,num_materials
        if (im.eq.im_primary) then
         LS_new(im)=LS_primary
        else if (is_rigid_local(im).eq.1) then
         LS_new(im)=LS(im)
        else if (is_rigid_local(im).eq.0) then
         if (LS(im).gt.-LS_primary) then
          LS_new(im)=-LS_primary
         else
          LS_new(im)=LS(im)
         endif
        else
         print *,"is_rigid_local invalid MOF.F90: ",is_rigid_local
         stop
        endif
       enddo ! im=1..num_materials
      else if ((is_rigid_local(im_primary).eq.0).and. &
               (is_elastic_local(im_primary).eq.1)) then
       LS_primary=max(LS(im_primary),zero)
       if (LS_primary.ge.0.0d0) then 
        !do nothing
       else
        print *,"LS_primary invalid: ",LS_primary
        stop
       endif
       do im=1,num_materials
        if (im.eq.im_primary) then
         LS_new(im)=LS_primary
        else if (is_elastic_local(im).eq.1) then
         LS_new(im)=LS(im)
        else if (is_rigid_local(im).eq.1) then
         LS_new(im)=LS(im)
        else if ((is_rigid_local(im).eq.0).and. &
                 (is_elastic_local(im).eq.0)) then
         if (LS(im).gt.-LS_primary) then
          LS_new(im)=-LS_primary
         else
          LS_new(im)=LS(im)
         endif
        else
         print *,"is_rigid_local invalid MOF.F90: ",is_rigid_local
         print *,"or is_elastic_local invalid MOF.F90: ",is_elastic_local
         stop
        endif
       enddo ! im=1..num_materials
      else
       print *,"is_rigid_local invalid MOF.F90: ",is_rigid_local
       print *,"or is_elastic_local invalid MOF.F90: ",is_elastic_local
       stop
      endif

      do im=1,num_materials
       LS_new_hold(im)=LS_new(im)
      enddo

      do im=1,num_materials
       LS_max_rest=-99999.0d0
       do im_rest=1,num_materials
        if (im_rest.ne.im) then
         if (LS_new_hold(im_rest).gt.LS_max_rest) then
          LS_max_rest=LS_new_hold(im_rest)
         else if (LS_new_hold(im_rest).le.LS_max_rest) then
          ! do nothing
         else
          print *,"LS_max_rest bust"
          stop
         endif
        endif
       enddo ! im_rest=1..num_materials
       LS_new(im)=half*(LS_new(im)-LS_max_rest)
      enddo ! im=1..num_materials

      end subroutine LS_tessellate

      subroutine multi_get_volume_tessellate( &
       tid_in, &
       tessellate_source, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
       tessellate_dest, & !TESSELLATE_ALL|ALL_RASTER|FLUIDS_ELASTIC
       bfact,dx, &
       xsten0,nhalf0, &
       mofdata, &
       xtetlist, &
       nlist_alloc, &
       nmax, &
       sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: nlist_alloc
      integer, INTENT(in) :: nmax
      integer, INTENT(in) :: sdim
      integer, PARAMETER :: shapeflag=0 !regular hexahedron
      integer, INTENT(in) :: bfact,nhalf0
      integer, INTENT(in) :: tessellate_source
      integer, INTENT(in) :: tessellate_dest
      real(amrex_real) xtet(sdim+1,sdim)
      real(amrex_real), INTENT(inout) :: mofdata(num_materials*(2*sdim+3))
      real(amrex_real) :: mofdata_convert(num_materials*(2*sdim+3))
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)

      integer :: local_tessellate

      integer, parameter :: mof_verbose=0
      integer, parameter :: use_ls_data=0
       ! D_DECL(i,j,k) = i,j  in 2D
       !               = i,j,k in 3D
      real(amrex_real) :: &
            LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)

      real(amrex_real) multi_volume(num_materials)
      real(amrex_real) multi_cen(sdim,num_materials)
      real(amrex_real) multi_area(num_materials)
      real(amrex_real), INTENT(out) :: xtetlist(4,3,nlist_alloc)

      integer im
      integer vofcomp
      real(amrex_real) volcell
      real(amrex_real) volcell_compare
      real(amrex_real) cencell(sdim)
      integer dir
      real(amrex_real) vof_super(num_materials)
      real(amrex_real) multi_centroidA(num_materials,sdim)

      integer, parameter :: caller_id=3

      if ((tid_in.ge.geom_nthreads).or.(tid_in.lt.0)) then
       print *,"tid_in invalid (multi_get_volume_tessellate) ",tid_in
       stop
      endif 

      if ((tessellate_source.eq.TESSELLATE_FLUIDS).or. &
          (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC)) then
       !do nothing
      else
       print *,"tessellate_source invalid ",tessellate_source
       stop
      endif

      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon.ne.2*sdim+3"
       stop
      endif

      if (nmax.lt.4) then
       print *,"nmax invalid multi_get_volume_tessellate nmax=",nmax
       stop
      endif

      if ((nlist_alloc.ge.1).and.(nlist_alloc.le.nmax)) then
       ! do nothing
      else
       print *,"nlist_alloc invalid"
       stop
      endif

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid multi get volume tessellate"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid135"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_volume_tessellate"
       stop
      endif
      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid multi get volume tessellate"
       stop
      endif
 
       ! sum Frigid <=1
       ! sum Ffluid = 1
      call make_vfrac_sum_ok_base( &
       xsten0, &
       nhalf0, &
       bfact,dx, &
       tessellate_source, & !=TESSELLATE_FLUIDS|IGNORE_ISELASTIC
       mofdata, &
       sdim)

      call multi_get_volume_grid( &
       caller_id, &
       tid_in, &
       EPS_FULL_WEAK, & ! tolerance for "single material" criterion
       tessellate_source, & !=TESSELLATE_FLUIDS|IGNORE_ISELASTIC
       tessellate_dest, & ! =TESSELLATE_ALL|ALL_RASTER|FLUIDS_ELASTIC
       bfact,dx, &
       xsten0,nhalf0, &
       mofdata, &
       xsten0,nhalf0, &
       xtet, &
       multi_volume, &
       multi_cen, & !(sdim,num_materials) absolute frame of ref.
       multi_area, &
       xtetlist, &
       nlist_alloc, &
       nmax, &
       sdim, &
       shapeflag) 

      call Box_volumeFAST(bfact,dx, &
        xsten0,nhalf0, &
        volcell,cencell,sdim)

      if (volcell.gt.zero) then
       !do nothing
      else
       print *,"expecting volcell>0.0 ",volcell
       stop
      endif

      volcell_compare=zero
      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1
       if ((tessellate_dest.eq.TESSELLATE_ALL).or. &
           (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
        volcell_compare=volcell_compare+multi_volume(im)
       else if (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC) then
        if (is_rigid(im).eq.1) then
         !do nothing
        else if (is_rigid(im).eq.0) then
         volcell_compare=volcell_compare+multi_volume(im)
        else
         print *,"is_rigid(im) invalid ",im,is_rigid(im)
         stop
        endif
       else
        print *,"tessellate_dest invalid ",tessellate_dest
        stop
       endif
      enddo !im=1,num_materials

      if (abs(volcell-volcell_compare).le.EPS2*volcell) then
       !do nothing
      else
       print *,"tessellate_source ",tessellate_source
       print *,"tessellate_dest ",tessellate_dest
       print *,"multi_volume ",multi_volume
       print *,"volcell or volcell_compare invalid: ", &
         volcell,volcell_compare, &
         (volcell-volcell_compare)/volcell 
       stop
      endif

       ! zap out override slope and order.
      do im=1,num_materials*ngeom_recon
       mofdata_convert(im)=zero
      enddo

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1
       mofdata_convert(vofcomp)=multi_volume(im)/volcell
       vof_super(im)=mofdata_convert(vofcomp)
       do dir=1,sdim
        mofdata_convert(vofcomp+dir)=multi_cen(dir,im)-cencell(dir)
       enddo
      enddo !im=1,num_materials

      if (tessellate_dest.eq.TESSELLATE_ALL) then
       local_tessellate=TESSELLATE_IGNORE_ISRIGID
      else if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then
       local_tessellate=TESSELLATE_IGNORE_ISRIGID
      else if (tessellate_dest.eq.TESSELLATE_FLUIDS_ELASTIC) then
       local_tessellate=TESSELLATE_IGNORE_ISELASTIC
      else
       print *,"tessellate_dest invalid: ",tessellate_dest
       stop
      endif

      call multimaterial_MOF( &
        local_tessellate, & !TESSELLATE_IGNORE_ISRIGID|ISELASTIC
        tid_in, &
        bfact,dx, &
        xsten0,nhalf0, &
        mof_verbose, &
        use_ls_data, &
        LS_stencil, &
        xtetlist, &
        nlist_alloc, &
        nmax, &
        mofdata_convert, & !intent(inout)
        vof_super, &
        multi_centroidA, &
        sdim)

      do im=1,num_materials*ngeom_recon
       mofdata(im)=mofdata_convert(im)
      enddo

      return
      end subroutine multi_get_volume_tessellate

      subroutine update_touchLS(newLS,minLS,maxLS,touch_hold,im,sdim)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: im,sdim
      real(amrex_real), INTENT(in) :: newLS(num_materials*(1+sdim)) 
      integer, INTENT(inout) :: touch_hold(num_materials)
      real(amrex_real), INTENT(inout) :: minLS(num_materials)
      real(amrex_real), INTENT(inout) :: maxLS(num_materials)
      integer ctouch

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif 
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid37"
       stop
      endif

      ctouch=touch_hold(im)
      if ((ctouch.eq.0).or.(ctouch.eq.1).or.(ctouch.eq.2)) then
       if (newLS(im).lt.minLS(im)) then
        minLS(im)=newLS(im)
       endif
       if (newLS(im).gt.maxLS(im)) then
        maxLS(im)=newLS(im)
       endif
       touch_hold(im)=2
      else
       print *,"ctouch invalid"
       stop
      endif

      return
      end subroutine update_touchLS

      subroutine compare_distance( &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
       bfact,dx, &
       xsten0,nhalf0, &
       xaccept, &
       xdonate, & !intent(in)
       newLS, &
       touch_hold, &
       minLS, &
       maxLS, &
       im_test,n_im, &
       slope,imslope, &
       imcell,sdim, &
       center_stencil, &
       donateflag) !intent(in)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: tessellate
      integer, INTENT(in) :: n_im
      integer, INTENT(in) :: bfact,nhalf0
      integer, INTENT(in) :: im_test(n_im)
      integer, INTENT(in) :: imslope,imcell
      integer, INTENT(in) :: sdim,center_stencil
      integer, INTENT(in) :: donateflag(num_materials+1)
      real(amrex_real), INTENT(in) :: xaccept(sdim) 
      real(amrex_real), INTENT(in) :: xdonate(sdim) 
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim) 
      real(amrex_real), INTENT(in) :: dx(sdim) 
      real(amrex_real), INTENT(in) :: slope(sdim) 
      real(amrex_real), INTENT(inout) :: newLS(num_materials*(1+sdim)) 
      integer, INTENT(inout) :: touch_hold(num_materials)
      real(amrex_real), INTENT(inout) :: minLS(num_materials)
      real(amrex_real), INTENT(inout) :: maxLS(num_materials)
      integer distzero,im_here,im_opp_here
      integer im,im3,dir,nc
      real(amrex_real) disttest,dist_compare,LSSIGN
      real(amrex_real) slopetest(sdim)
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      if ((tessellate.eq.TESSELLATE_FLUIDS).or. &
          (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC)) then
       !do nothing
      else
       print *,"expecting tessellate.eq.TESSELLATE_FLUIDS ", &
         tessellate
       print *,"or expect tessellate.eq.TESSELLATE_IGNORE_ISELASTIC ", &
         tessellate
       stop
      endif
      call init_local_material_vars( & !subroutine compare_distance
       is_rigid_local, &
       is_elastic_local, &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
       tessellate, &
       is_proper_layer_local)

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid135"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid compare_distance"
       stop
      endif
      if ((num_materials.lt.1).or.(num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid multi_get_distance"
       stop
      endif
      if ((n_im.lt.1).or.(n_im.gt.6)) then
       print *,"n_im invalid"
       stop
      endif

      do nc=1,n_im
       if ((im_test(nc).lt.1).or.(im_test(nc).gt.num_materials)) then
        print *,"im_test invalid: ",nc,im_test(nc)
        stop
       endif
       if ((is_rigid_local(im_test(nc)).eq.0).and. &
           (is_elastic_local(im_test(nc)).eq.0)) then
        !do nothing
       else
        print *,"is_rigid_local invalid MOF.F90: ", &
          nc,im_test(nc),is_rigid_local(im_test(nc))
        print *,"or is_elastic_local invalid MOF.F90: ", &
          nc,im_test(nc),is_elastic_local(im_test(nc))
        print *,"compare_distance"
        stop
       endif
      enddo ! nc=1..n_im

      if ((imcell.lt.1).or.(imcell.gt.num_materials)) then
       print *,"imcell invalid imcell=",imcell
       print *,"put breakpoint here to see caller"
       stop
      endif
      if ((imslope.lt.0).or.(imslope.gt.num_materials)) then
       print *,"imslope invalid: ",imslope
       stop
      endif
      if ((center_stencil.ne.0).and.(center_stencil.ne.1)) then
       print *,"center_stencil invalid: ",center_stencil
       stop
      endif

      disttest=zero
      do dir=1,sdim
       disttest=disttest+(xaccept(dir)-xdonate(dir))**2
       slopetest(dir)=xaccept(dir)-xdonate(dir)
      enddo
      disttest=sqrt(disttest)
      distzero=0
      if (disttest.lt.EPS_8_4*dx(1)) then
       distzero=1
       if (center_stencil.eq.1) then
        do dir=1,sdim
         slopetest(dir)=slope(dir)
        enddo
       else if (center_stencil.eq.0) then
        print *,"center_stencil invalid when disttest=0: ",center_stencil
        stop
       else
        print *,"center_stencil invalid: ",center_stencil
        stop
       endif
      else if (disttest.gt.zero) then
       ! xclosest=x -  d n
       ! n will point from xclosest to x.
       do dir=1,sdim
        slopetest(dir)=slopetest(dir)/disttest
       enddo
      else
       print *,"disttest invalid: ",disttest
       stop
      endif

      do im=1,num_materials
       if ((donateflag(im).ne.0).and. &
           (donateflag(im).ne.1)) then
        print *,"donateflag invalid: im,donateflag = ",im,donateflag(im)
        print *,"num_materials=",num_materials
        stop
       endif

       if ((is_rigid_local(im).eq.0).and. &
           (is_elastic_local(im).eq.0)) then

        im_here=0
        im_opp_here=0
        do nc=1,n_im
         im3=im_test(nc)
          
         if (donateflag(im3).eq.1) then
          if (im3.eq.im) then
           im_here=1
          endif
          if (im3.ne.im) then
           im_opp_here=1
          endif
         endif
        enddo  ! nc

        dist_compare=abs(newLS(im))
 
         ! disttest=|xaccept-xdonate|
        if ((dist_compare.gt.disttest).or. &
            (touch_hold(im).eq.0).or. &
            (touch_hold(im).eq.1)) then
         if (center_stencil.eq.1) then
          LSSIGN=zero
          if ((imcell.eq.im).and.(im_opp_here.eq.1)) then
           LSSIGN=one
          endif
          if ((imcell.ne.im).and.(im_here.eq.1)) then
           LSSIGN=-one
          endif
          if (LSSIGN.ne.zero) then
           newLS(im)=LSSIGN*disttest  
           if (distzero.eq.1) then
            if ((imslope.ge.1).and.(imslope.le.num_materials)) then
             if (imslope.eq.im) then
              LSSIGN=one
             else
              LSSIGN=-one
             endif
            else
             print *,"imslope invalid imslope= ",imslope
             print *,"num_materials ",num_materials
             print *,"xaccept ",xaccept(1),xaccept(2),xaccept(sdim) 
             print *,"xdonate ",xdonate(1),xdonate(2),xdonate(sdim) 
             print *,"xsten0 ",xsten0(0,1),xsten0(0,2),xsten0(0,sdim) 
             print *,"dx ",dx(1),dx(2),dx(sdim) 
             do nc=1,n_im
              print *,"nc,im_test ",nc,im_test(nc)
             enddo
             print *,"imcell,center_stencil ",imcell,center_stencil
             stop
            endif
           else if (distzero.eq.0) then
            ! do nothing
           else
            print *,"distzero invalid: ",distzero
            stop
           endif
           do dir=1,sdim
            newLS(num_materials+sdim*(im-1)+dir)=LSSIGN*slopetest(dir)
           enddo
           call update_touchLS(newLS,minLS,maxLS,touch_hold,im,sdim)
          else if (LSSIGN.eq.zero) then
           ! do nothing
          else
           print *,"LSSIGN BUST"
           stop
          endif
         else if (center_stencil.eq.0) then
          if (distzero.ne.0) then
           print *,"distzero invalid: ",distzero
           stop
          endif
          LSSIGN=zero
          if (newLS(im).ge.zero) then
           im_here=1
          else
           im_opp_here=1
          endif
          if ((newLS(im).ge.zero).and.(im_opp_here.eq.1)) then 
           LSSIGN=one
          endif
          if ((newLS(im).lt.zero).and.(im_here.eq.1)) then
           LSSIGN=-one
          endif
          if (LSSIGN.ne.zero) then
           newLS(im)=LSSIGN*disttest  
           do dir=1,sdim
            newLS(num_materials+sdim*(im-1)+dir)=LSSIGN*slopetest(dir)
           enddo
           call update_touchLS(newLS,minLS,maxLS,touch_hold,im,sdim)
          else if (LSSIGN.eq.zero) then
           ! do nothing
          else
           print *,"LSSIGN BUST: ",LSSIGN
           stop
          endif
         else
          print *,"center_stencil invalid"
          stop
         endif
        else if ((dist_compare.le.disttest).and. &
                 (touch_hold(im).eq.2)) then
         ! do nothing
        else
         print *,"dist_compare, disttest, or touch_hold invalid: ", &
                 dist_compare,disttest,im,touch_hold(im)
         stop
        endif   

       else if ((is_rigid_local(im).eq.1).or. &
                (is_elastic_local(im).eq.1)) then
        ! do nothing
       else
        print *,"is_rigid_local invalid MOF.F90 ",is_rigid_local
        print *,"or is_elastic_local invalid MOF.F90 ",is_elastic_local
        stop
       endif 
      enddo ! im=1..num_materials
  
      return
      end subroutine compare_distance 

! find the closest point (x_cp) on the intersection of two planes (in
! 3D) or two lines (2D) to a point (x_a)
! NOTE: Assumed the vectors slope_list are "unit" normals
      subroutine closestINT( &
         bfact,dx,xsten0,nhalf0, &
         x_cp,xplus,xminus, &
         xplus2,xminus2, &
         x_a,maxdx,ilist,jlist,nlist, &
         slope_list,intercept_list,im_list,sdim, &
         inboxflag)
       use probcommon_module
       use global_utility_module

       IMPLICIT NONE

       integer, INTENT(in) :: sdim,nlist,bfact,nhalf0
       real(amrex_real), INTENT(in) :: maxdx
       real(amrex_real), INTENT(out) :: x_cp(sdim)
       real(amrex_real), INTENT(out) :: xplus(sdim)
       real(amrex_real), INTENT(out) :: xminus(sdim)
       real(amrex_real), INTENT(out) :: xplus2(sdim)
       real(amrex_real), INTENT(out) :: xminus2(sdim)
       real(amrex_real), INTENT(in) :: x_a(sdim)
       integer, INTENT(in) :: ilist, jlist
       real(amrex_real), INTENT(in) :: slope_list(num_materials,sdim)
       real(amrex_real), INTENT(in) :: intercept_list(num_materials)
       integer, INTENT(in) :: im_list(num_materials)
       real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
       real(amrex_real), INTENT(in) :: dx(sdim)
       integer, INTENT(out) :: inboxflag

       integer dir,islope
       real(amrex_real) c_i,c_j,s_c
       real(amrex_real) n_i_d_n_j
       real(amrex_real) n_i_c_n_j(sdim)
       real(amrex_real) dtrmn,mag
       real(amrex_real) h(2)
       integer indexlist(2)

       if (nhalf0.lt.1) then
        print *,"nhalf0 invalid"
        stop
       endif
       if (bfact.lt.1) then
        print *,"bfact invalid135"
        stop
       endif
       if (maxdx.gt.zero) then
        !do nothing
       else
        print *,"maxdx invalid (closestINT): ",maxdx
        stop
       endif

       if ((nlist.lt.2).or.(nlist.gt.num_materials-1).or. &
           (ilist.lt.1).or.(ilist.gt.nlist).or. &
           (jlist.lt.1).or.(jlist.gt.nlist).or. &
           (ilist.eq.jlist)) then
        print *,"ilist,jlist, or nlist invalid"
        stop
       endif

       inboxflag = 0
       indexlist(1)=ilist
       indexlist(2)=jlist

        ! P_i and P_j are in the form n.(x-x_0)+d=0
        ! transform them to the form n.x=h
        ! h=n.x_0-d
       do islope=1,2
        h(islope)=zero
        do dir=1,sdim
         h(islope)=h(islope)+slope_list(indexlist(islope),dir)*xsten0(0,dir)
        enddo
        h(islope)=h(islope)-intercept_list(indexlist(islope))
       enddo ! islope

       if (sdim.eq.2) then
         ! solve the 2-by-2 linear system Ax=b
         ! A = [n_i1 n_i2]
         !     [n_j1 n_j2]
         ! x = [x_cp_1]
         !     [x_cp_2]
         ! h = [n_i1*x_01+n_i2*x_02-d_i]
         !     [n_j1*x_01+n_j2*x_02-d_j]

         dtrmn = slope_list(ilist,1)*slope_list(jlist,2)- &
                 slope_list(ilist,2)*slope_list(jlist,1)

         if (dtrmn.eq.zero) then
          inboxflag=0
         else
          x_cp(1)= &
           (slope_list(jlist,2)*h(1)- &
            slope_list(ilist,2)*h(2))/dtrmn
          x_cp(2)= &
           (-slope_list(jlist,1)*h(1)+ &
             slope_list(ilist,1)*h(2))/dtrmn
           ! in: closestINT
          call check_inbox(x_cp,xsten0,nhalf0,inboxflag)
         endif

       else if (sdim.eq.3) then

        ! The intersction line I is 
        ! I:{x|x=c_i n_i + c_j n_j + c_k (n_i cross n_j)} 
        ! c_i = (h_i-h_j(n_i.n_j))/(1-(n_i.n_j)^2)
        ! c_j = (h_j-h_i(n_i.n_j))/(1-(n_i.n_j)^2)

        n_i_d_n_j=zero
        do dir=1,sdim
         n_i_d_n_j=n_i_d_n_j+slope_list(ilist,dir)*slope_list(jlist,dir)
        enddo

        dtrmn=one-n_i_d_n_j*n_i_d_n_j
        if (dtrmn.le.zero) then
         inboxflag=0
        else
         c_i =(h(1)-h(2)*n_i_d_n_j)/dtrmn
         c_j =(h(2)-h(1)*n_i_d_n_j)/dtrmn

         n_i_c_n_j(1)=slope_list(ilist,2)*slope_list(jlist,3) &
                     -slope_list(ilist,3)*slope_list(jlist,2)
         n_i_c_n_j(2)=slope_list(ilist,3)*slope_list(jlist,1) &
                     -slope_list(ilist,1)*slope_list(jlist,3)
         n_i_c_n_j(3)=slope_list(ilist,1)*slope_list(jlist,2) &
                     -slope_list(ilist,2)*slope_list(jlist,1)


         ! The x_a in the {n_i,n_j,n_i cross n_j} basis form is
         ! x_a = s_i n_i + s_j n_j + s_c (n_i cross n_j)
         ! s_c = x_a.(n_i cross n_j)
         ! The x_cp would be
         ! x_cp = c_i n_i + c_j n_j + s_c (n_i cross n_j)

         mag=zero
         do dir=1,sdim
          mag=mag+n_i_c_n_j(dir)**2
         enddo
         mag=sqrt(mag)

         if (mag.gt.zero) then
          do dir=1,sdim
           n_i_c_n_j(dir)=n_i_c_n_j(dir)/mag
          enddo
 
          s_c=zero
          do dir=1,sdim
           s_c=s_c+x_a(dir)*n_i_c_n_j(dir)
          enddo

          do dir=1,sdim
           x_cp(dir)=c_i*slope_list(ilist,dir)+&
                     c_j*slope_list(jlist,dir)+&
                     s_c*n_i_c_n_j(dir)
          enddo
           ! in: closestINT
          call check_inbox(x_cp,xsten0,nhalf0,inboxflag)
         else
          inboxflag=0
         endif
        endif

       else
        print *,"MOF.F90::closestINT - sdim invalid!"
        stop
       endif ! sdim

       if (inboxflag.eq.0) then
        ! do nothing
       else if (inboxflag.eq.1) then
        do dir=1,sdim
         xplus(dir)=x_cp(dir)-maxdx*EPS4*half* &
           (slope_list(ilist,dir)+slope_list(jlist,dir)) 
         xminus(dir)=x_cp(dir)+maxdx*EPS4*half* &
           (slope_list(ilist,dir)+slope_list(jlist,dir)) 
         xplus2(dir)=x_cp(dir)-maxdx*EPS4*half* &
           (slope_list(ilist,dir)-slope_list(jlist,dir)) 
         xminus2(dir)=x_cp(dir)+maxdx*EPS4*half* &
           (slope_list(ilist,dir)-slope_list(jlist,dir)) 
        enddo
       else
        print *,"inboxflag invalid"
        stop
       endif

      end subroutine closestINT


      subroutine closestINT3( &
        bfact,dx,xsten0,nhalf0, &
        x_cp,xplus,xminus, &
        xplus2,xminus2, &
        x_a,maxdx,ilist,jlist,nlist, &
        slope_list,intercept_list,im_list, &
        gphi,&
        x_0side,sdim,inboxflag)
       use probcommon_module
       use global_utility_module

       IMPLICIT NONE

       integer, INTENT(in) :: sdim,nlist,bfact,nhalf0
       real(amrex_real), INTENT(in) :: maxdx
       real(amrex_real), INTENT(out) :: x_cp(sdim)
       real(amrex_real), INTENT(out) :: xplus(sdim)
       real(amrex_real), INTENT(out) :: xminus(sdim)
       real(amrex_real), INTENT(out) :: xplus2(sdim)
       real(amrex_real), INTENT(out) :: xminus2(sdim)
       real(amrex_real), INTENT(in) :: x_a(sdim)
       real(amrex_real), INTENT(in) :: gphi(sdim)
       real(amrex_real), INTENT(in) :: x_0side(sdim)
       integer, INTENT(in) :: ilist, jlist
       real(amrex_real), INTENT(in) :: slope_list(num_materials,sdim)
       real(amrex_real), INTENT(in) :: intercept_list(num_materials)
       integer, INTENT(in) :: im_list(num_materials)
       real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
       real(amrex_real), INTENT(in) :: dx(sdim)
       integer, INTENT(out) :: inboxflag

       integer dir,islope,dir_2
       real(amrex_real) c_i,c_j,c_k
       real(amrex_real) n_i_d_n_j
       real(amrex_real) n_i_c_n_j(sdim)
       real(amrex_real) dtrmn,mag
       real(amrex_real) h(2)
       integer cp_init
       integer indexlist(2)

       if (bfact.lt.1) then
        print *,"bfact invalid135"
        stop
       endif

       if (maxdx.gt.zero) then
        !do nothing
       else
        print *,"maxdx invalid (closestINT3): ",maxdx
        stop
       endif

       if ((nlist.lt.2).or.(nlist.gt.num_materials-1).or. &
           (ilist.lt.1).or.(ilist.gt.nlist).or. &
           (jlist.lt.1).or.(jlist.gt.nlist).or. &
           (ilist.eq.jlist)) then
        print *,"ilist,jlist, or nlist invalid"
        stop
       endif

       if (sdim.eq.3) then
        ! do nothing
       else
        print *,"MOF.F90::closestINT3 - sdim invalid!"
        stop
       endif

       inboxflag = 0
       indexlist(1)=ilist
       indexlist(2)=jlist

        ! P_i and P_j are in the form n.(x-x_0)+d=0
        ! transform them to the form n.x=h
        ! h=n.x_0-d
       do islope=1,2
        h(islope)=zero
        do dir=1,sdim
         h(islope)=h(islope)+slope_list(indexlist(islope),dir)*xsten0(0,dir)
        enddo
        h(islope)=h(islope)-intercept_list(indexlist(islope))
       enddo ! islope

       ! The intersction line I is 
       ! I:{x|x=c_i n_i + c_j n_j + c_k (n_i cross n_j)} 
       ! c_i = (h_i-h_j(n_i.n_j))/(1-(n_i.n_j)^2)
       ! c_j = (h_j-h_i(n_i.n_j))/(1-(n_i.n_j)^2)

       n_i_d_n_j=zero
       do dir=1,sdim
        n_i_d_n_j=n_i_d_n_j+slope_list(ilist,dir)*slope_list(jlist,dir)
       enddo

       dtrmn=one-n_i_d_n_j*n_i_d_n_j
       if (dtrmn.le.zero) then
        inboxflag=0
       else
        c_i =(h(1)-h(2)*n_i_d_n_j)/dtrmn
        c_j =(h(2)-h(1)*n_i_d_n_j)/dtrmn

        n_i_c_n_j(1)=slope_list(ilist,2)*slope_list(jlist,3) &
                    -slope_list(ilist,3)*slope_list(jlist,2)
        n_i_c_n_j(2)=slope_list(ilist,3)*slope_list(jlist,1) &
                    -slope_list(ilist,1)*slope_list(jlist,3)
        n_i_c_n_j(3)=slope_list(ilist,1)*slope_list(jlist,2) &
                    -slope_list(ilist,2)*slope_list(jlist,1)

        mag=zero
        do dir=1,sdim
         mag=mag+n_i_c_n_j(dir)**2
        enddo
        mag=sqrt(mag)

        if (mag.gt.zero) then
         do dir=1,sdim
          n_i_c_n_j(dir)=n_i_c_n_j(dir)/mag
         enddo

         ! Find the intersection of line I with a cell face 
         ! Solve for c_k
         ! I:{x|x=c_i n_i + c_j n_j + c_k (n_i cross n_j)} 
         !                 and
         ! x=x_0side(1) or y=x_0side(2) or z=x_0side(3) 
         !        {based on the direction of gphi}

         cp_init=0
         do dir=1,sdim
          if (gphi(dir).eq.zero) then
           ! do nothing
          else if (gphi(dir).eq.one) then
           if (cp_init.eq.1) then
            print *,"cp already init"
            stop
           else if (cp_init.eq.0) then
            cp_init=1
           else
            print *,"cp_init bust"
            stop
           endif
            ! Face is in this direction
           if(n_i_c_n_j(dir).ne.zero) then
            ! Line I is NOT parallel to the face
            c_k = (x_0side(dir)-c_i*slope_list(ilist,dir)&
                  -c_j*slope_list(jlist,dir))/n_i_c_n_j(dir)

            do dir_2=1,sdim
             x_cp(dir_2)=c_i*slope_list(ilist,dir_2)+&
                         c_j*slope_list(jlist,dir_2)+&
                         c_k*n_i_c_n_j(dir_2)
            enddo
             ! in: closestINT3
            call check_inbox(x_cp,xsten0,nhalf0,inboxflag)
           else
            inboxflag=0
           endif
          else
           print *,"gphi invalid"
           stop
          endif
         enddo !dir
         if (cp_init.ne.1) then
          print *,"face plane incorrectly specified"
          stop
         endif
        else
         inboxflag=0
        endif
       endif ! dtrmn<>0 ?

       if (inboxflag.eq.0) then
        ! do nothing
       else if (inboxflag.eq.1) then
        do dir=1,sdim
         xplus(dir)=x_cp(dir)-maxdx*EPS4*half* &
           (slope_list(ilist,dir)+slope_list(jlist,dir)) 
         xminus(dir)=x_cp(dir)+maxdx*EPS4*half* &
           (slope_list(ilist,dir)+slope_list(jlist,dir)) 
         xplus2(dir)=x_cp(dir)-maxdx*EPS4*half* &
           (slope_list(ilist,dir)-slope_list(jlist,dir)) 
         xminus2(dir)=x_cp(dir)+maxdx*EPS4*half* &
           (slope_list(ilist,dir)-slope_list(jlist,dir)) 
        enddo
       else
        print *,"inboxflag invalid"
        stop
       endif

       end subroutine closestINT3

         ! n dot (x-x0) + intercept = 0
         ! n dot x = d    d=n dot x0 - intercept
         ! P: ( x= a1 t1 + a2 t2 + a3 n )  a3=d/(n dot n)
         ! xaccept=b1 t1 + b2 t2 + b3 n  b3=xaccept dot n/(n dot n)
         ! xcp=b1 t1 + b2 t2 + a3 n=xaccept+(a3-b3)n
         ! if n dot n=1,
         ! xcp=xaccept+(d-xaccept dot n)n=
         !     xaccept+(n dot x0-int-xaccept dot n)n=
         !     xaccept-(n dot (xaccept-x0)+int)n
         ! xsten0 can be different from xstenbox since the LS
         ! might be projected to a lower dimension.
       subroutine closestPLANE( &
         bfact,dx,xsten0,nhalf0, &
         xcp,xcp_plus,xcp_minus,xaccept, &
         maxdx,slope_in,intercept_in, &
         xstenbox,nhalfbox,sdim,inboxflag)
       use probcommon_module
       use global_utility_module
       IMPLICIT NONE

       integer, INTENT(in) :: sdim
       integer, INTENT(out) :: inboxflag
       integer, INTENT(in) :: bfact,nhalf0,nhalfbox
       integer dir
       real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
       real(amrex_real), INTENT(in) :: xstenbox(-nhalfbox:nhalfbox,sdim)
       real(amrex_real), INTENT(in) :: dx(sdim)
       real(amrex_real), INTENT(out) :: xcp(sdim)
       real(amrex_real), INTENT(out) :: xcp_plus(sdim)
       real(amrex_real), INTENT(out) :: xcp_minus(sdim)
       real(amrex_real), INTENT(in) :: xaccept(sdim)
       real(amrex_real), INTENT(in) :: maxdx
       real(amrex_real), INTENT(in) :: slope_in(sdim)
       real(amrex_real) slope(sdim)
       real(amrex_real), INTENT(in) :: intercept_in
       real(amrex_real) intercept,dist

       if (nhalf0.lt.1) then
        print *,"nhalf0 invalid"
        stop
       endif
       if (nhalfbox.lt.1) then
        print *,"nhalfbox invalid"
        stop
       endif
       if ((sdim.ne.2).and.(sdim.ne.3)) then
        print *,"sdim invalid"
        stop
       endif

       if (maxdx.gt.zero) then
        !do nothing
       else
        print *,"maxdx invalid (closestPLANE): ",maxdx
        stop
       endif
       inboxflag=1
       dist=zero
       do dir=1,sdim
        dist=dist+slope_in(dir)**2
       enddo
       dist=sqrt(dist)
       if (dist.gt.zero) then
        do dir=1,sdim
         slope(dir)=slope_in(dir)/dist
        enddo
        intercept=intercept_in/dist
        dist=intercept
        do dir=1,sdim
         dist=dist+slope(dir)*(xaccept(dir)-xsten0(0,dir))
        enddo
        do dir=1,sdim
         xcp(dir)=xaccept(dir)-dist*slope(dir)
         xcp_plus(dir)=xcp(dir)-maxdx*EPS4*slope(dir)
         xcp_minus(dir)=xcp(dir)+maxdx*EPS4*slope(dir)
        enddo
         ! in: closestPLANE
        call check_inbox(xcp,xstenbox,nhalfbox,inboxflag)
       else
        inboxflag=0
       endif

       return
       end subroutine closestPLANE

        ! get the slope of the fluid material whose interface is closest to
        ! the center of the cell.
       subroutine get_primary_slope( &
         tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
         bfact,dx,xsten0,nhalf0, &
         mofdata, &
         slope,imslope,sdim)
       use probcommon_module
       use geometry_intersect_module
       use global_utility_module
       IMPLICIT NONE

       integer, INTENT(in) :: tessellate
       integer, INTENT(in) :: sdim
       integer, INTENT(out) :: imslope
       integer, INTENT(in) :: bfact,nhalf0
       real(amrex_real), INTENT(in) :: mofdata(num_materials*(2*sdim+3))
       real(amrex_real) mofdatavalid(num_materials*(2*sdim+3))
       real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
       real(amrex_real), INTENT(in) :: dx(sdim)
       real(amrex_real), INTENT(out) :: slope(sdim)
       real(amrex_real) slope_test(sdim)
       real(amrex_real) intercept,intercept_min
       real(amrex_real) vfrac_data(num_materials)
       integer sorted_list(num_materials)
       real(amrex_real) uncaptured_volume
       real(amrex_real) mag
       integer im,vofcomp,irank,testflag,dir
       integer is_rigid_local(num_materials)
       integer is_elastic_local(num_materials)
       integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      if ((tessellate.eq.TESSELLATE_FLUIDS).or. &
          (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC)) then
       !do nothing
      else
       print *,"expecting tessellate.eq.TESSELLATE_FLUIDS ", &
         tessellate
       print *,"or expect tessellate.eq.TESSELLATE_IGNORE_ISELASTIC ", &
         tessellate
       stop
      endif
       call init_local_material_vars( & !subroutine get_primary_slope
         is_rigid_local, &
         is_elastic_local, &
         tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
         tessellate, &
         is_proper_layer_local)

       if (bfact.lt.1) then
        print *,"bfact invalid135"
        stop
       endif
       if (nhalf0.lt.1) then
        print *,"nhalf0 invalid"
        stop
       endif

       if (ngeom_recon.ne.2*sdim+3) then
        print *,"ngeom_recon.ne.2*sdim+3"
        stop
       endif

       if ((sdim.ne.3).and.(sdim.ne.2)) then
        print *,"sdim invalid get_primary_slope"
        stop
       endif
     
       imslope=0
       do dir=1,sdim
        slope(dir)=zero
       enddo
       intercept_min=1.0e+20

        ! sum F_fluid=1  
        ! sum F_solid<=1 
       call make_vfrac_sum_ok_copy( &
         xsten0,nhalf0, &
         bfact,dx, &
         tessellate, & ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
         mofdata,mofdatavalid,sdim)

       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        vfrac_data(im)=mofdatavalid(vofcomp)
       enddo
       call sort_volume_fraction(vfrac_data,tessellate,sorted_list)
       im=sorted_list(1)
       if ((is_rigid_local(im).eq.0).and. &
           (is_elastic_local(im).eq.0)) then
        ! do nothing
       else
        print *,"is_rigid_local invalid MOF.F90 ",is_rigid_local
        print *,"or is_elastic_local invalid MOF.F90 ",is_elastic_local
        stop
       endif

       if (vfrac_data(im).ge.one-EPS_UNCAPTURED) then
        ! do nothing, there are no reconstructed interfaces in the cell.
       else if ((vfrac_data(im).ge.EPS_UNCAPTURED).and. &
                (vfrac_data(im).lt.one-EPS_UNCAPTURED)) then

        uncaptured_volume=one
        irank=1

        do while ((irank.le.num_materials).and. &
                  (uncaptured_volume.gt.zero))
         do im=1,num_materials
          if ((is_rigid_local(im).eq.0).and. &
              (is_elastic_local(im).eq.0)) then
           vofcomp=(im-1)*ngeom_recon+1
           testflag=NINT(mofdatavalid(vofcomp+sdim+1))
           if (testflag.eq.irank) then
            do dir=1,sdim
             slope_test(dir)=mofdatavalid(vofcomp+sdim+1+dir)
            enddo
            intercept=mofdatavalid(vofcomp+2*sdim+2)

            uncaptured_volume=uncaptured_volume-mofdatavalid(vofcomp)
            if (uncaptured_volume.lt.EPS_UNCAPTURED) then
             uncaptured_volume=zero
            endif
            if (uncaptured_volume.gt.zero) then ! we have a valid interface.
             if (abs(intercept).lt.abs(intercept_min)) then
              imslope=im
              intercept_min=intercept
              do dir=1,sdim
               slope(dir)=slope_test(dir)
              enddo
             endif
            endif
           endif ! testflag=irank?
          else if ((is_rigid_local(im).eq.1).or. &
                   (is_elastic_local(im).eq.1)) then
           ! do nothing
          else
           print *,"is_rigid_local invalid MOF.F90 ",is_rigid_local
           print *,"or is_elastic_local invalid MOF.F90 ",is_elastic_local
           print *,"get_primary_slope  MOF.F90"
           stop
          endif
         enddo ! im=1..num_materials
         irank=irank+1
        enddo ! while irank<=num_materials and uncaptured _vol>0
       else
        print *,"vfrac_data max out of range: ",im,vfrac_data(im)
        stop
       endif ! vfrac(im)>1-eps ?  (im==sorted_list(1))

       if (imslope.eq.0) then
        !do nothing
       else if ((imslope.ge.1).and.(imslope.le.num_materials)) then
        mag=zero
        do dir=1,sdim
         mag=mag+slope(dir)**2
        enddo 
        mag=sqrt(mag)
        if (mag.eq.zero) then
         imslope=0
        else if (mag.gt.zero) then
         do dir=1,sdim
          slope(dir)=slope(dir)/mag
         enddo
        else
         print *,"mag invalid: ",mag
         stop
        endif
       else
        print *,"imslope invalid: ",imslope
        stop
       endif

       return
       end subroutine get_primary_slope

       
      subroutine multi_get_distance( &
        tessellate, &
        bfact,dx, &
        xsten_recon,nhalf_recon, &
        xgrid, &
        mofdata, &
        multi_distance, &
        touch_hold, &
        minLS, &
        maxLS, &
        sdim, &
        center_stencil, &
        donateflag)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: tessellate
      integer, INTENT(in) :: nhalf_recon
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: sdim
      integer dir,side,l
      integer, INTENT(in) :: center_stencil
      integer, INTENT(in) :: donateflag(num_materials+1)
      real(amrex_real), INTENT(in) :: mofdata(num_materials*(2*sdim+3))
      real(amrex_real) mofdatavalid(num_materials*(2*sdim+3))
      real(amrex_real), INTENT(in) :: xsten_recon(-nhalf_recon:nhalf_recon,sdim)
      real(amrex_real) xstenface_recon(-nhalf_recon:nhalf_recon,sdim)
      real(amrex_real) x0side(sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(inout) :: multi_distance(num_materials*(1+sdim))
      integer, INTENT(inout) :: touch_hold(num_materials)
      real(amrex_real), INTENT(inout) :: minLS(num_materials)
      real(amrex_real), INTENT(inout) :: maxLS(num_materials)
      integer irank,vofcomp,im,im_plus,im_minus,im0
      integer im_plus2,im_minus2
      real(amrex_real) vfrac_data(num_materials)
      integer sorted_list(num_materials)
      real(amrex_real) uncaptured_volume
      integer testflag
      real(amrex_real) slopes(sdim)
      real(amrex_real) intercept,intercept_face
      real(amrex_real) maxdx,normgphi
      integer inboxflag
      real(amrex_real) xx(sdim)
      real(amrex_real) gphi(sdim)
      real(amrex_real) slope_list(num_materials,sdim)
      real(amrex_real) intercept_list(num_materials)
      integer im_list(num_materials)
      integer nlist
      real(amrex_real) x0(sdim)
      integer dir1,dir2,side1,ilist,jlist
      integer n_im
      integer im_test(6)
      real(amrex_real) x0face(sdim)
      real(amrex_real), INTENT(in) :: xgrid(sdim)
      real(amrex_real) xgrid_cen(sdim)
      real(amrex_real) xgrid_plus(sdim)
      real(amrex_real) xgrid_minus(sdim)
      real(amrex_real) xgrid_plus2(sdim)
      real(amrex_real) xgrid_minus2(sdim)
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      if ((tessellate.eq.TESSELLATE_FLUIDS).or. &
          (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC)) then
       !do nothing
      else
       print *,"expecting tessellate.eq.TESSELLATE_FLUIDS ", &
         tessellate
       print *,"or expect tessellate.eq.TESSELLATE_IGNORE_ISELASTIC ", &
         tessellate
       stop
      endif
      call init_local_material_vars( & !subroutine multi_get_distance
       is_rigid_local, &
       is_elastic_local, &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
       tessellate, &
       is_proper_layer_local)

      if ((tessellate.eq.TESSELLATE_FLUIDS).or. &
          (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC)) then
       !do nothing
      else
       print *,"expecting tessellate=TESSELLATE_FLUIDS|IGNORE_ISELASTIC"
       stop
      endif

      if (nhalf_recon.lt.1) then
       print *,"nhalf_recon invalid"
       stop
      endif
      if ((center_stencil.ne.0).and.(center_stencil.ne.1)) then
       print *,"center_stencil invalid: ",center_stencil
       stop
      endif
      maxdx=dx(1)
      if (dx(2).gt.maxdx) then
       maxdx=dx(2)
      endif
      if (dx(sdim).gt.maxdx) then
       maxdx=dx(sdim)
      endif

      if (maxdx.gt.zero) then
       !do nothing
      else
       print *,"maxdx invalid (multi_get_distance): ",maxdx
       stop
      endif

      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon.ne.2*sdim+3"
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_distance"
       stop
      endif
      if ((num_materials.lt.1).or.(num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid multi_get_distance"
       stop
      endif

       ! sum F_fluid = 1
       ! sum F_solid <= 1
      call make_vfrac_sum_ok_copy( &
        xsten_recon,nhalf_recon, &
        bfact,dx, &
        tessellate, &  ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
        mofdata,mofdatavalid,sdim)

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1
       vfrac_data(im)=mofdatavalid(vofcomp)
      enddo
      call sort_volume_fraction(vfrac_data, &
         tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
         sorted_list)
      im=sorted_list(1)
      if ((is_rigid_local(im).eq.0).and. &
          (is_elastic_local(im).eq.0)) then
       ! do nothing
      else
       print *,"is_rigid_local invalid MOF.F90 ",is_rigid_local
       print *,"or is_elastic_local invalid MOF.F90 ",is_elastic_local
       stop
      endif
       ! a full cell, so distance is either +bigdist or -bigdist,
       ! and a default normal is used.
      if (vfrac_data(im).ge.one-EPS_UNCAPTURED) then
       ! do nothing, there are no reconstructed interfaces in the cell.
      else if ((vfrac_data(im).ge.EPS_UNCAPTURED).and. &
               (vfrac_data(im).lt.one-EPS_UNCAPTURED)) then

       do dir=1,sdim
        x0(dir)=xsten_recon(0,dir)
       enddo
       call multi_get_volumePOINT( &
        tessellate, &  ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
        bfact,dx,xsten_recon,nhalf_recon, &
        mofdata,x0,im0,sdim)

       uncaptured_volume=one
       irank=1
       nlist=0
       do while ((irank.le.num_materials).and. &
                 (uncaptured_volume.gt.zero))
        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         testflag=NINT(mofdatavalid(vofcomp+sdim+1))

         if ((is_rigid_local(im).eq.0).and. &
             (is_elastic_local(im).eq.0)) then

          if (testflag.eq.irank) then
           do dir=1,sdim
            slopes(dir)=mofdatavalid(vofcomp+sdim+1+dir)
           enddo
           intercept=mofdatavalid(vofcomp+2*sdim+2)

           uncaptured_volume=uncaptured_volume-mofdatavalid(vofcomp)
           if (uncaptured_volume.lt.EPS_UNCAPTURED) then
            uncaptured_volume=zero
           endif
           if (uncaptured_volume.gt.zero) then ! we have a valid interface.
            nlist=nlist+1
            if ((nlist.ge.1).and.(nlist.lt.num_materials)) then
             do dir=1,sdim
              slope_list(nlist,dir)=slopes(dir)
             enddo
             intercept_list(nlist)=intercept
             im_list(nlist)=im
            else 
             print *,"nlist invalid"
             stop
            endif
            call closestPLANE( &
             bfact,dx,xsten_recon,nhalf_recon, &
             xgrid_cen,xgrid_plus,xgrid_minus, &
             xgrid,maxdx,slopes,intercept, &
             xsten_recon,nhalf_recon,sdim,inboxflag) 

            if (inboxflag.eq.1) then
             call multi_get_volumePOINT( &
              tessellate, & ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
              bfact,dx,xsten_recon,nhalf_recon, &
              mofdata,xgrid_plus, &
              im_plus,sdim)
             call multi_get_volumePOINT( &
              tessellate, & ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
              bfact,dx,xsten_recon,nhalf_recon, &
              mofdata,xgrid_minus, &
              im_minus,sdim)
             n_im=2
             im_test(1)=im_plus
             im_test(2)=im_minus
             if (center_stencil.eq.1) then
              n_im=n_im+1
              im_test(n_im)=im0
             endif
             if (irank.eq.1) then
              n_im=n_im+1
              im_test(n_im)=im
             endif
             ! xgrid_cen is closest point on the plane.
             call compare_distance( &
              tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
              bfact,dx, &
              xsten_recon,nhalf_recon, &
              xgrid,xgrid_cen, &
              multi_distance, &
              touch_hold, &
              minLS, &
              maxLS, &
              im_test,n_im, &
              slopes, &
              im,im0,sdim, &
              center_stencil, &
              donateflag)
            else if (inboxflag.eq.0) then
             ! do nothing
            else
             print *,"inboxflag invalid"
             stop
            endif

            do dir=1,sdim
            do side=-1,1,2
             do l=1,sdim
              xx(l)=xgrid(l)
              gphi(l)=slopes(l)
              x0face(l)=x0(l)
             enddo
             x0face(dir)=xsten_recon(side,dir)
             xx(dir)=x0face(dir)
             gphi(dir)=zero
             normgphi=zero
             do l=1,sdim
              normgphi=normgphi+gphi(l)**2
             enddo 
             normgphi=sqrt(normgphi)
  
             if (normgphi.ge.EPS12) then 
              !  n dot (x-x0)+intercept=0
              !  n dot (x-x0+x0face-x0face)+intercept=0
              !  n dot (x-x0face) + n dot (x0face-x0)+intercept=0
              intercept_face=intercept

              do dir2=1,sdim
               intercept_face=intercept_face+ &
                 slopes(dir2)*(x0face(dir2)-x0(dir2))
               xstenface_recon(0,dir2)=x0face(dir2)
              enddo
           
              call closestPLANE( &
               bfact,dx,xstenface_recon,nhalf_recon, &
               xgrid_cen,xgrid_plus,xgrid_minus, &
               xx,maxdx,gphi,intercept_face, &
               xsten_recon,nhalf_recon,sdim,inboxflag) 

              if (inboxflag.eq.1) then
               call multi_get_volumePOINT( &
                 tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
                 bfact,dx,xsten_recon,nhalf_recon, &
                 mofdata,xgrid_plus, &
                 im_plus,sdim)
               call multi_get_volumePOINT( &
                 tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
                 bfact,dx,xsten_recon,nhalf_recon, &
                 mofdata,xgrid_minus, &
                 im_minus,sdim)
               n_im=2
               im_test(1)=im_plus
               im_test(2)=im_minus
               if (center_stencil.eq.1) then
                n_im=n_im+1
                im_test(n_im)=im0
               endif
               if (irank.eq.1) then
                n_im=n_im+1
                im_test(n_im)=im
               endif

               call compare_distance( &
                tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
                bfact,dx, &
                xsten_recon,nhalf_recon, &
                xgrid,xgrid_cen, &
                multi_distance, &
                touch_hold, &
                minLS, &
                maxLS, &
                im_test,n_im, &
                slopes, &
                im,im0,sdim, &
                center_stencil, &
                donateflag)
              else if (inboxflag.eq.0) then
               ! do nothing
              else
               print *,"inboxflag invalid"
               stop
              endif

              if (sdim.eq.2) then
               ! do nothing
              else if (sdim.eq.3) then
  
               do dir1=1,sdim
               do side1=-1,1,2
                if (dir1.ne.dir) then
                 do l=1,sdim
                  xx(l)=xgrid(l)
                  gphi(l)=slopes(l)
                  x0face(l)=x0(l)
                 enddo
                 x0face(dir)=xsten_recon(side,dir)
                 x0face(dir1)=xsten_recon(side1,dir1)
                 xx(dir)=x0face(dir)
                 xx(dir1)=x0face(dir1)
                 gphi(dir)=zero
                 gphi(dir1)=zero
                 normgphi=zero
                 do l=1,sdim
                  normgphi=normgphi+gphi(l)**2
                 enddo 
                 normgphi=sqrt(normgphi)
   
                 if (normgphi.ge.EPS12) then
                  intercept_face=intercept
                  do dir2=1,sdim
                   intercept_face=intercept_face+ &
                    slopes(dir2)*(x0face(dir2)-x0(dir2))
                   xstenface_recon(0,dir2)=x0face(dir2)
                  enddo
                  call closestPLANE( &
                   bfact,dx,xstenface_recon,nhalf_recon, &
                   xgrid_cen,xgrid_plus,xgrid_minus, &
                   xx,maxdx,gphi,intercept_face, &
                   xsten_recon,nhalf_recon,sdim,inboxflag) 

                  if (inboxflag.eq.1) then
                   call multi_get_volumePOINT( &
                    tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
                    bfact,dx,xsten_recon,nhalf_recon, &
                    mofdata,xgrid_plus, &
                    im_plus,sdim)
                   call multi_get_volumePOINT( &
                    tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
                    bfact,dx,xsten_recon,nhalf_recon, &
                    mofdata,xgrid_minus, &
                    im_minus,sdim)
                   n_im=2
                   im_test(1)=im_plus
                   im_test(2)=im_minus
                   if (center_stencil.eq.1) then
                    n_im=n_im+1
                    im_test(n_im)=im0
                   endif
                   if (irank.eq.1) then
                    n_im=n_im+1
                    im_test(n_im)=im
                   endif

                   call compare_distance( &
                    tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
                    bfact,dx, &
                    xsten_recon,nhalf_recon, &
                    xgrid,xgrid_cen, &
                    multi_distance, &
                    touch_hold, &
                    minLS, &
                    maxLS, &
                    im_test,n_im, &
                    slopes, &
                    im,im0,sdim, &
                    center_stencil, &
                    donateflag)
                  else if (inboxflag.eq.0) then
                   ! do nothing
                  else
                   print *,"inboxflag invalid"
                   stop
                  endif
                 endif ! normgphi>0 ?
                endif ! dir1<>dir
               enddo ! side1
               enddo ! dir1

              else
               print *,"sdim invalid"
               stop
              endif
             endif ! normgphi>0 ?
            enddo ! side
            enddo ! dir

           endif ! uncaptured_volume>0
          endif  ! testflag=irank
         else if ((is_rigid_local(im).eq.1).or. &
                  (is_elastic_local(im).eq.1))  then
          ! do nothing
         else
          print *,"is_rigid_local invalid MOF.F90 ",is_rigid_local
          print *,"or is_elastic_local invalid MOF.F90 ",is_elastic_local
          print *,"multi_get_distance"
          stop
         endif
        enddo ! im=1,num_materials
        irank=irank+1
       enddo  ! while irank<=num_materials and uncaptured_volume>0 

       if (nlist.ge.2) then
        do ilist=1,nlist-1
        do jlist=ilist+1,nlist
      

          ! xgrid_cen=xgrid-d n
          ! n=(xgrid-xgrid_cen)/d 
          ! if no intersection, or not in the box, then
          ! inboxflag=0.
         call closestINT( &
          bfact,dx,xsten_recon,nhalf_recon, &
          xgrid_cen,xgrid_plus,xgrid_minus, &
          xgrid_plus2,xgrid_minus2, &
          xgrid,maxdx,ilist,jlist,nlist, &
          slope_list,intercept_list,im_list,sdim, &
          inboxflag)
         if (inboxflag.eq.1) then
          call multi_get_volumePOINT( &
            tessellate, &  ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
            bfact,dx,xsten_recon,nhalf_recon, &
            mofdata,xgrid_plus, &
            im_plus,sdim)
          call multi_get_volumePOINT( &
            tessellate, &  ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
            bfact,dx,xsten_recon,nhalf_recon, &
            mofdata,xgrid_minus, &
            im_minus,sdim)
          call multi_get_volumePOINT( &
            tessellate, &  ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
            bfact,dx,xsten_recon,nhalf_recon, &
            mofdata,xgrid_plus2, &
            im_plus2,sdim)
          call multi_get_volumePOINT( &
            tessellate, &  ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
            bfact,dx,xsten_recon,nhalf_recon, &
            mofdata,xgrid_minus2, &
            im_minus2,sdim)
          n_im=4
          im_test(1)=im_plus
          im_test(2)=im_minus
          im_test(3)=im_plus2
          im_test(4)=im_minus2
          if (center_stencil.eq.1) then
           n_im=n_im+1
           im_test(n_im)=im0
          endif
          do dir2=1,sdim
           slopes(dir2)=slope_list(ilist,dir2)
          enddo
          im=im_list(ilist)
          call compare_distance( &
            tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
            bfact,dx, &
            xsten_recon,nhalf_recon, &
            xgrid,xgrid_cen, &
            multi_distance, &
            touch_hold, &
            minLS, &
            maxLS, &
            im_test,n_im, &
            slopes, &
            im,im0,sdim, &
            center_stencil, &
            donateflag)
         else if (inboxflag.eq.0) then
          ! do nothing
         else
          print *,"inboxflag invalid"
          stop
         endif 

          ! this is case where 2 planes, and a cell wall plane have
          ! a single common intersection point.
         if (sdim.eq.3) then

          do dir=1,sdim
          do side=-1,1,2
           do l=1,sdim
            xx(l)=xgrid(l)
            gphi(l)=0.0
            x0side(l)=x0(l)
           enddo
           gphi(dir)=one
           xx(dir)=xsten_recon(side,dir)
           x0side(dir)=xx(dir)

           call closestINT3( &
            bfact,dx,xsten_recon,nhalf_recon, &
            xgrid_cen,xgrid_plus,xgrid_minus, &
            xgrid_plus2,xgrid_minus2, &
            xx,maxdx,ilist,jlist,nlist, &
            slope_list,intercept_list,im_list, &
            gphi, &
            x0side,sdim, &
            inboxflag)

           if (inboxflag.eq.1) then
            call multi_get_volumePOINT( &
             tessellate, & ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
             bfact,dx,xsten_recon,nhalf_recon, &
             mofdata,xgrid_plus, &
             im_plus,sdim)
            call multi_get_volumePOINT( &
             tessellate, & ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
             bfact,dx,xsten_recon,nhalf_recon, &
             mofdata,xgrid_minus, &
             im_minus,sdim)
            call multi_get_volumePOINT( &
             tessellate, & ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
             bfact,dx,xsten_recon,nhalf_recon, &
             mofdata,xgrid_plus2, &
             im_plus2,sdim)
            call multi_get_volumePOINT( &
             tessellate, & ! =TESSELLATE_FLUIDS|IGNORE_ISELASTIC
             bfact,dx,xsten_recon,nhalf_recon, &
             mofdata,xgrid_minus2, &
             im_minus2,sdim)

            n_im=4
            im_test(1)=im_plus
            im_test(2)=im_minus
            im_test(3)=im_plus2
            im_test(4)=im_minus2
            if (center_stencil.eq.1) then
             n_im=n_im+1
             im_test(n_im)=im0
            endif
            do dir2=1,sdim
             slopes(dir2)=slope_list(ilist,dir2)
            enddo
            im=im_list(ilist)

            call compare_distance( &
              tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
              bfact,dx, &
              xsten_recon,nhalf_recon, &
              xgrid,xgrid_cen, &
              multi_distance, &
              touch_hold, &
              minLS, &
              maxLS, &
              im_test,n_im, &
              slopes, &
              im,im0,sdim, &
              center_stencil, &
              donateflag)
           else if (inboxflag.eq.0) then
            ! do nothing
           else
            print *,"inboxflag invalid"
            stop
           endif 

          enddo ! side
          enddo ! dir
 
         else if (sdim.eq.2) then
           ! do nothing
         else
          print *,"sdim invalid"
          stop
         endif
        enddo 
        enddo  ! ilist,jlist
       endif ! nlist >=2
 
      else
       print *,"vfrac_data(im) out of range: ",im,vfrac_data(im)
       stop
      endif ! vfrac(im)>1-eps ?  im==sorted_list(1)

      return
      end subroutine multi_get_distance


! vof, ref centroid, order,slope,intercept  x num_materials
! phi=n dot (x-x0) + intercept
! x0 is center of cell (not centroid)
! in: MOF_routines_module
      subroutine multi_get_volumePOINT( &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
       bfact,dx, &
       xsten0,nhalf0, & ! absolute coordinate system.
       mofdata, &
       xgrid, &  ! absolute coordinate system.
       im_crit, &
       sdim)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      integer, INTENT(in) :: tessellate
      integer, INTENT(in) :: bfact,nhalf0
      real(amrex_real), INTENT(in) :: mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdatavalid(num_materials*ngeom_recon)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim)
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: xgrid(sdim)
      integer, INTENT(out) :: im_crit
      integer irank,vofcomp,im
      real(amrex_real) uncaptured_volume_fraction
      integer testflag,dir
      real(amrex_real) maxvof
      real(amrex_real) slopes(sdim)
      real(amrex_real) intercept,ls
      real(amrex_real) vfrac_data(num_materials)
      integer vfrac_checked(num_materials)
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)


#include "mofdata.H"

      call init_local_material_vars( & !subroutine multi_get_volumePOINT
       is_rigid_local, &
       is_elastic_local, &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
       tessellate, &
       is_proper_layer_local)

      if (ngeom_recon.ne.(2*sdim+3)) then
       print *,"ngeom_recon invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid135"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_volumePOINT"
       stop
      endif
      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid multi get volume point"
       stop
      endif

       ! sum voffluid=1 , sum vofsolid <= 1
      call make_vfrac_sum_ok_copy( &
        xsten0,nhalf0, &
        bfact,dx, &
        tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC|IGNORE_ISRIGID
        mofdata,mofdatavalid,sdim)

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1
       vfrac_data(im)=mofdatavalid(vofcomp)
      enddo ! im=1..num_materials

       !in: multi_get_volume_POINT
      call check_full_cell_vfrac( &
        vfrac_data, &
        tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC|IGNORE_ISRIGID
        tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC|IGNORE_ISRIGID
        im_crit, &
        EPS_8_4)

      if ((im_crit.ge.1).and. &
          (im_crit.le.num_materials)) then

       ! do nothing, im_crit is set.

      else if (im_crit.eq.0) then
      
       do im=1,num_materials
        vfrac_checked(im)=0
       enddo

       uncaptured_volume_fraction=one
       irank=1
       do while ((irank.le.num_materials).and. &
                 (uncaptured_volume_fraction.gt.zero))

        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         testflag=NINT(mofdatavalid(vofcomp+sdim+1))

         if ((is_rigid_local(im).eq.0).and. &
             (is_elastic_local(im).eq.0)) then

          if (testflag.eq.irank) then
           vfrac_checked(im)=1
           do dir=1,sdim
            slopes(dir)=mofdatavalid(vofcomp+sdim+1+dir)
           enddo
           intercept=mofdatavalid(vofcomp+2*sdim+2)
           call distfunc(bfact,dx,xsten0,nhalf0, &
            intercept,slopes,xgrid,ls,sdim)

           if ((ls.ge.zero).or. &
               (mofdatavalid(vofcomp).ge. &
                (one-EPS_FULL_WEAK)*uncaptured_volume_fraction)) then
            im_crit=im
            uncaptured_volume_fraction=zero
           else 
            uncaptured_volume_fraction=uncaptured_volume_fraction- &
              mofdatavalid(vofcomp)
           endif
          endif  ! testflag=irank
         else if ((is_rigid_local(im).eq.1).or. &
                  (is_elastic_local(im).eq.1)) then
          ! do nothing
         else
          print *,"is_rigid_local invalid MOF.F90: ",is_rigid_local
          print *,"or is_elastic_local invalid MOF.F90: ",is_elastic_local
          print *,"multi_get_volumePOINT"
          stop
         endif
           
        enddo ! im=1..num_materials
        irank=irank+1
       enddo  ! while irank<=num_materials and uncaptured_volume_fraction>0 

       if ((im_crit.ge.1).and.(im_crit.le.num_materials)) then

        if (uncaptured_volume_fraction.eq.zero) then
         !do nothing
        else
         print *,"uncaptured_volume_fraction invalid: ", &
           uncaptured_volume_fraction
         stop
        endif

       else if (im_crit.eq.0) then

        maxvof=zero
        do im=1,num_materials
         if ((is_rigid_local(im).eq.0).and. &
             (is_elastic_local(im).eq.0)) then
          if (vfrac_checked(im).eq.1) then
           ! do nothing
          else if (vfrac_checked(im).eq.0) then
           if (vfrac_data(im).gt.maxvof) then
            maxvof=vfrac_data(im)
            im_crit=im
           endif
          else
           print *,"vfrac_checked invalid"
           stop
          endif
         else if ((is_rigid_local(im).eq.1).or. &
                  (is_elastic_local(im).eq.1)) then
          ! do nothing
         else
          print *,"is_rigid_local invalid MOF.F90: ",is_rigid_local
          print *,"or is_elastic_local invalid MOF.F90: ",is_elastic_local
          print *,"multi_get_volumePOINT"
          stop
         endif
        enddo ! im=1..num_materials

        if ((im_crit.ge.1).and.(im_crit.le.num_materials)) then
         ! do nothing
        else if (im_crit.eq.0) then

         maxvof=zero
         do im=1,num_materials
          if ((is_rigid_local(im).eq.0).and. &
              (is_elastic_local(im).eq.0)) then
           if (vfrac_data(im).gt.maxvof) then
            maxvof=vfrac_data(im)
            im_crit=im
           endif
          else if ((is_rigid_local(im).eq.1).or. &
                   (is_elastic_local(im).eq.1)) then
           ! do nothing
          else
           print *,"is_rigid_local invalid MOF.F90: ",is_rigid_local
           print *,"or is_elastic_local invalid MOF.F90: ",is_elastic_local
           print *,"multi_get_volumePOINT"
           stop
          endif
         enddo ! im=1..num_materials

         if (maxvof.gt.zero) then
          ! do nothing
         else if (maxvof.le.zero) then
          print *,"failed to find material that covers point"
          do im=1,num_materials
           vofcomp=(im-1)*ngeom_recon+1
           print *,"im,vof,flag,int ",im,mofdatavalid(vofcomp), &
            NINT(mofdatavalid(vofcomp+sdim+1)), &
            mofdatavalid(vofcomp+2*sdim+2)
          enddo ! im=1..num_materials
          print *,"xgrid,xsten0 ",xgrid(1),xgrid(2),xsten0(0,1),xsten0(0,2)
          stop
         else
          print *,"maxvof invalid: ",maxvof
          stop
         endif 

        else 
         print *,"im_crit invalid 24338: ",im_crit
         stop
        endif

       else
        print *,"im_crit invalid 24343"
        stop
       endif

      else
       print *,"im_crit invalid 24348: ",im_crit
       stop
      endif

      return
      end subroutine multi_get_volumePOINT


      subroutine check_full_cell_vfrac( &
         vfrac, &
         tessellate_source, &
         tessellate_dest, &
         im_full,tol)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module
 
      IMPLICIT NONE

      real(amrex_real), intent(in) :: tol
      integer, INTENT(in) :: tessellate_source
      integer, INTENT(in) :: tessellate_dest
      real(amrex_real), INTENT(in) :: vfrac(num_materials)
      integer, INTENT(out) :: im_full
      integer im
      integer im_fluid_max,im_solid_max,im_elastic_max
      real(amrex_real) sum_solid_vfrac,sum_fluid_vfrac,sum_elastic_vfrac
      real(amrex_real) max_solid_vfrac,max_fluid_vfrac,max_elastic_vfrac
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      if ((tol.gt.zero).and.(tol.lt.one)) then
       !do nothing
      else
       print *,"tol invalid: ",tol
       stop
      endif

      if (tessellate_source.eq.TESSELLATE_FLUIDS) then
       if ((tessellate_dest.eq.TESSELLATE_FLUIDS).or. &
           (tessellate_dest.eq.TESSELLATE_ALL).or. &
           (tessellate_dest.eq.TESSELLATE_ALL_RASTER)) then
        !do nothing
       else
        print *,"tessellate_dest invalid: ",tessellate_dest
        stop
       endif
      else if (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC) then
       if (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC) then
        !do nothing
       else
        print *,"tessellate_dest invalid: ",tessellate_dest
        stop
       endif
      else if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then
       if (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID) then
        !do nothing
       else
        print *,"tessellate_dest invalid: ",tessellate_dest
        stop
       endif
      else
       print *,"tessellate_source invalid check_full_cell_vfrac: ", &
         tessellate_source
       stop
      endif

      call init_local_material_vars( & !subroutine check_full_cell_vfrac
       is_rigid_local, &
       is_elastic_local, &
       tessellate_source, & !TESSELLATE_FLUIDS|IGNORE_ISRIGID|IGNORE_ISELASTIC
       tessellate_dest, &
       is_proper_layer_local)

      if ((num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid check_full_cell_vfrac"
       print *,"num_materials= ",num_materials
       stop
      endif

      do im=1,num_materials
       if ((vfrac(im).ge.-tol).and. &
           (vfrac(im).le.one+tol)) then
        ! do nothing
       else
        print *,"check_full_cell_vfrac: "
        print *,"vfrac out of range im,vfrac(im): ",im,vfrac(im)
        stop
       endif
      enddo ! im=1..num_materials

      im_full=0

      sum_solid_vfrac=zero
      sum_fluid_vfrac=zero
      sum_elastic_vfrac=zero

      max_solid_vfrac=zero
      max_fluid_vfrac=zero
      max_elastic_vfrac=zero

      im_fluid_max=0
      im_solid_max=0
      im_elastic_max=0

      do im=1,num_materials

       if (is_rigid_local(im).eq.1) then

        sum_solid_vfrac=sum_solid_vfrac+vfrac(im)
        if (im_solid_max.eq.0) then
         im_solid_max=im
         max_solid_vfrac=vfrac(im)
        else if ((im_solid_max.ge.1).and. &
                 (im_solid_max.le.num_materials)) then
         if (vfrac(im).gt.vfrac(im_solid_max)) then
          im_solid_max=im
          max_solid_vfrac=vfrac(im)
         endif
        else
         print *,"im_solid_max invalid"
         stop
        endif

       else if (is_elastic_local(im).eq.1) then

        sum_elastic_vfrac=sum_elastic_vfrac+vfrac(im)
        if (im_elastic_max.eq.0) then
         im_elastic_max=im
         max_elastic_vfrac=vfrac(im)
        else if ((im_elastic_max.ge.1).and. &
                 (im_elastic_max.le.num_materials)) then
         if (vfrac(im).gt.vfrac(im_elastic_max)) then
          im_elastic_max=im
          max_elastic_vfrac=vfrac(im)
         endif
        else
         print *,"im_elastic_max invalid"
         stop
        endif

       else if ((is_rigid_local(im).eq.0).and. &
                (is_elastic_local(im).eq.0)) then

        sum_fluid_vfrac=sum_fluid_vfrac+vfrac(im)
        if (im_fluid_max.eq.0) then
         im_fluid_max=im
         max_fluid_vfrac=vfrac(im)
        else if ((im_fluid_max.ge.1).and.(im_fluid_max.le.num_materials)) then
         if (vfrac(im).gt.vfrac(im_fluid_max)) then
          im_fluid_max=im
          max_fluid_vfrac=vfrac(im)
         endif
        else
         print *,"im_fluid_max invalid (check_full_cell_vfrac) ",im_fluid_max
         stop
        endif

       else
        print *,"is_rigid_local or is_elastic_local invalid MOF.F90"
        stop
       endif

      enddo !im=1..num_materials

      if (tessellate_source.eq.TESSELLATE_FLUIDS) then

       if (tessellate_dest.eq.TESSELLATE_FLUIDS) then

        if (max_fluid_vfrac.ge.one-tol) then
         im_full=im_fluid_max
        endif

       else if (tessellate_dest.eq.TESSELLATE_ALL) then

        if (max_solid_vfrac.ge.one-tol) then
         im_full=im_solid_max
        endif

        if (sum_solid_vfrac.le.tol) then
         if (max_elastic_vfrac.ge.one-tol) then
          im_full=im_elastic_max
         endif
         if (sum_elastic_vfrac.le.tol) then
          if (max_fluid_vfrac.ge.one-tol) then
           im_full=im_fluid_max
          endif
         endif
        endif

       else if (tessellate_dest.eq.TESSELLATE_ALL_RASTER) then

        if (max_solid_vfrac.ge.half) then
         im_full=im_solid_max
        else if (max_solid_vfrac.le.half) then

         if (max_elastic_vfrac.ge.one-tol) then
          im_full=im_elastic_max
         endif
         if (sum_elastic_vfrac.le.tol) then
          if (max_fluid_vfrac.ge.one-tol) then
           im_full=im_fluid_max
          endif
         endif
        else
         print *,"max_solid_vfrac invalid: ",max_solid_vfrac
         stop
        endif

       else
        print *,"tessellate_dest invalid check_full_cell_vfrac: ", &
         tessellate_dest
        stop
       endif

      else if (tessellate_source.eq.TESSELLATE_IGNORE_ISELASTIC) then

       if (tessellate_dest.eq.TESSELLATE_IGNORE_ISELASTIC) then

        if (max_fluid_vfrac.ge.one-tol) then
         im_full=im_fluid_max
        endif

       else
        print *,"tessellate_dest invalid check_full_cell_vfrac: ", &
         tessellate_dest
        stop
       endif

      else if (tessellate_source.eq.TESSELLATE_IGNORE_ISRIGID) then

       if (tessellate_dest.eq.TESSELLATE_IGNORE_ISRIGID) then

        if (max_fluid_vfrac.ge.one-tol) then
         im_full=im_fluid_max
        endif

       else
        print *,"tessellate_dest invalid check_full_cell_vfrac: ", &
         tessellate_dest
        stop
       endif

      else
       print *,"tessellate_source invalid check_full_cell_vfrac: ", &
         tessellate_source
       stop
      endif
     
      end subroutine check_full_cell_vfrac

      subroutine get_secondary_material(dx,LS,im_primary,im_secondary,sdim)
      use probcommon_module
      use geometry_intersect_module
      use global_utility_module
 
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: dx(sdim)
      real(amrex_real), INTENT(in) :: LS(num_materials)
      integer, INTENT(in) :: im_primary
      integer, INTENT(out) :: im_secondary
      integer im
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      real(amrex_real) :: dist_cutoff
      real(amrex_real) :: second_dist_cutoff
      real(amrex_real) :: dx_sanity
      integer dir

      if ((sdim.eq.2).or.(sdim.eq.3)) then
       !do nothing
      else
       print *,"sdim invalid: ",sdim
       stop
      endif

      if ((im_primary.ge.1).and.(im_primary.le.num_materials)) then
       ! do nothing
      else
       print *,"im_primary invalid get_secondary_material: ",im_primary
       stop
      endif

      do dir=1,sdim
       dx_sanity=problen_array(dir)*half
       if ((dx(dir).gt.zero).and.(dx(dir).le.dx_sanity)) then
        !do nothing
       else
        print *,"dx invalid in get_secondary_material: ",dx
        print *,"problen_array=",problen_array
        stop
       endif
      enddo !dir=1,SDIM

      dist_cutoff=dx(1)*EPS3

      second_dist_cutoff=-abs(LS(im_primary))-dist_cutoff

      do im=1,num_materials
       is_rigid_local(im)=is_rigid(im)
       is_elastic_local(im)=is_elastic(im)
      enddo ! im=1..num_materials

      im_secondary=0

       ! check if rigid material is second.
      do im=1,num_materials
       if (im.ne.im_primary) then

        if (is_rigid_local(im).eq.1) then
         if (LS(im).ge.second_dist_cutoff) then
          im_secondary=im
         endif
        endif

       else if (im.eq.im_primary) then
        ! do nothing
       else
        print *,"im bust: ",im,im_primary
        stop
       endif
      enddo !im=1..num_materials

      if (im_secondary.eq.0) then

        ! check if elastic material is second.
       do im=1,num_materials
        if (im.ne.im_primary) then

         if (is_elastic_local(im).eq.1) then
          if (LS(im).ge.second_dist_cutoff) then
           im_secondary=im
          endif
         endif

        else if (im.eq.im_primary) then
         ! do nothing
        else
         print *,"im bust: ",im,im_primary
         stop
        endif

       enddo !im=1..num_materials

       if (im_secondary.eq.0) then
        do im=1,num_materials
         if (im.ne.im_primary) then

          if (im_secondary.eq.0) then
           im_secondary=im
          else if (LS(im).gt.LS(im_secondary)) then
           im_secondary=im
          else if (LS(im).le.LS(im_secondary)) then
           ! do nothing
          else
           print *,"LS bust in get_secondary_material: ",LS(im)
           stop
          endif
         else if (im.eq.im_primary) then
          ! do nothing
         else
          print *,"im bust: ",im,im_primary
          stop
         endif

        enddo !im=1..num_materials

       else if ((im_secondary.ge.1).and. &
                (im_secondary.le.num_materials).and. &
                (im_secondary.ne.im_primary)) then
        !do nothing
       else
        print *,"im_secondary invalid: ",im_secondary
        stop
       endif

      else if ((im_secondary.ge.1).and. &
               (im_secondary.le.num_materials).and. &
               (im_secondary.ne.im_primary)) then
       !do nothing
      else
       print *,"im_secondary invalid: ",im_secondary
       stop
      endif

      end subroutine get_secondary_material

      subroutine sort_volume_fraction( &
       vfrac_data, &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
       sorted_list)

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: tessellate
      real(amrex_real), INTENT(in) :: vfrac_data(num_materials)
      integer, INTENT(out) :: sorted_list(num_materials)
      integer im,changed,nsweeps,swap,do_swap
      integer is_rigid_local(num_materials)
      integer is_elastic_local(num_materials)
      integer is_proper_layer_local(num_materials, &
           RIGID_LAYER_INDEX:FLUID_LAYER_INDEX)

      if ((tessellate.eq.TESSELLATE_FLUIDS).or. &
          (tessellate.eq.TESSELLATE_IGNORE_ISELASTIC)) then
       !do nothing
      else
       print *,"expecting tessellate.eq.TESSELLATE_FLUIDS ", &
         tessellate
       print *,"or expect tessellate.eq.TESSELLATE_IGNORE_ISELASTIC ", &
         tessellate
       stop
      endif
      call init_local_material_vars( & !subroutine sort_volume_fraction
       is_rigid_local, &
       is_elastic_local, &
       tessellate, & !TESSELLATE_FLUIDS|IGNORE_ISELASTIC
       tessellate, &
       is_proper_layer_local)

      if ((num_materials.lt.1).or.(num_materials.gt.MAX_NUM_MATERIALS)) then
       print *,"num_materials invalid sort_volume_fraction"
       print *,"num_materials= ",num_materials
       stop
      endif

      do im=1,num_materials
       sorted_list(im)=im
      enddo
      changed=1
      nsweeps=0
      do while ((changed.eq.1).and.(nsweeps.lt.num_materials-1))
       changed=0
       do im=1,num_materials-nsweeps-1
        do_swap=0

        if (is_rigid_local(sorted_list(im)).eq.1) then
         do_swap=1
        else if (is_elastic_local(sorted_list(im)).eq.1) then
         do_swap=1
        else if (is_rigid_local(sorted_list(im+1)).eq.1) then
         do_swap=0
        else if (is_elastic_local(sorted_list(im+1)).eq.1) then
         do_swap=0
        else if ((is_rigid_local(sorted_list(im)).eq.0).and. &
                 (is_rigid_local(sorted_list(im+1)).eq.0).and. &
                 (is_elastic_local(sorted_list(im)).eq.0).and. &
                 (is_elastic_local(sorted_list(im+1)).eq.0)) then
         if (vfrac_data(sorted_list(im)).lt. &
             vfrac_data(sorted_list(im+1))) then
          do_swap=1
         else if (vfrac_data(sorted_list(im)).ge. &
                  vfrac_data(sorted_list(im+1))) then
          do_swap=0
         else
          print *,"vfrac_data invalid: ",vfrac_data
          print *,"sorted_list= ",sorted_list
          stop
         endif
        else
         print *,"is_rigid_local invalid MOF.F90 ",is_rigid_local
         print *,"or is_elastic_local invalid MOF.F90 ",is_elastic_local
         stop 
        endif
       
        if (do_swap.eq.1) then
         swap=sorted_list(im)
         sorted_list(im)=sorted_list(im+1)
         sorted_list(im+1)=swap
         changed=1
        endif

       enddo ! im=1,num_materials-nsweeps-1

       nsweeps=nsweeps+1
      enddo ! while ((changed.eq.1).and.(nsweeps.lt.num_materials-1))

      return
      end subroutine sort_volume_fraction
       

end module MOF_routines_module


module MOF_cpp_module

contains

        !fort_initmof is called from:
        !void fortran_parameters() (in NavierStokes.cpp)
      subroutine fort_initmof( &
       order_algorithm_in, &
       renormalize_order_algorithm_in, &
       denconst_local, &
       is_rigid_local, &
       MOFITERMAX_in, &
       MOF_DEBUG_RECON_in, &
       nthreads, &
       nmax_in) &
      bind(c,name='fort_initmof')

      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: nmax_in,nthreads
      integer, INTENT(in) :: order_algorithm_in(num_materials)
      integer, INTENT(in) :: renormalize_order_algorithm_in(num_materials)
      real(amrex_real), INTENT(in) :: denconst_local(num_materials)
      integer, INTENT(in) :: is_rigid_local(num_materials)
      integer, INTENT(in) :: MOFITERMAX_in
      integer, INTENT(in) :: MOF_DEBUG_RECON_in
      integer sdim,nmax_test
      integer im,sub_im
      integer default_flag
      integer irank

#include "mofdata.H"

      if ((num_materials.ge.2).and. &
          (num_materials.le.MAX_NUM_MATERIALS)) then
       ! do nothing
      else
       print *,"num_materials not initialized properly"
       print *,"num_materials=",num_materials
       stop
      endif

      MOF_DEBUG_RECON_COUNT=0
      MOF_DEBUG_RECON=MOF_DEBUG_RECON_in

      call set_MOFITERMAX(MOFITERMAX_in)

      if (nthreads.lt.1) then
       print *,"nthreads invalid"
       stop
      endif
      if (nmax_in.ne.POLYGON_LIST_MAX) then
       print *,"expecting nmax_in = POLYGON_LIST_MAX"
       print *,"nmax_in=",nmax_in
       print *,"POLYGON_LIST_MAX=",POLYGON_LIST_MAX
       stop
      endif

      geom_nthreads=nthreads
      geom_nmax=nmax_in
      allocate(geom_xtetlist_local(4,3,geom_nmax,geom_nthreads))
      allocate(geom_xtetlist(4,3,geom_nmax,geom_nthreads))
      allocate(geom_xtetlist_old(4,3,geom_nmax,geom_nthreads))
      allocate(geom_xtetlist_uncapt(4,3,geom_nmax,geom_nthreads))

      allocate(intercept_error_history(INTERCEPT_MAXITER,geom_nthreads))

      allocate(mof_calls(geom_nthreads,num_materials))
      allocate(mof_errors(geom_nthreads,num_materials))
      allocate(mof_iterations(geom_nthreads,num_materials))

      default_flag=1
      do im=1,num_materials
       if (renormalize_order_algorithm_in(im).eq.0) then
        !do nothing
       else if (renormalize_order_algorithm_in(im).gt.0) then
        default_flag=0
       else
        print *,"renormalize_order_algorithm_in(im) invalid"
        stop
       endif
      enddo
      if (default_flag.ne.0) then
       print *,"expecting default_flag==0"
       stop
      endif

      nonzero_ranks=0

      do im=1,num_materials
       if (denconst_local(im).gt.zero) then
        !do nothing
       else
        print *,"denconst_local(im) invalid"
        stop
       endif
       if (is_rigid_local(im).eq.1) then
        !do nothing
       else if (is_rigid_local(im).eq.0) then
        irank=1
        do sub_im=1,num_materials
         if (im.ne.sub_im) then
          if (is_rigid_local(sub_im).eq.1) then
           !do nothing
          else if (is_rigid_local(sub_im).eq.0) then
           if (renormalize_order_algorithm_in(sub_im).lt. &
               renormalize_order_algorithm_in(im)) then
            irank=irank+1
           else if (renormalize_order_algorithm_in(sub_im).eq. &
                    renormalize_order_algorithm_in(im)) then
            print *,"renormalize_order_algorithm_in invalid", &
                    renormalize_order_algorithm_in
            print *,"sub_im,im ",sub_im,im
            stop
           else if (renormalize_order_algorithm_in(sub_im).gt. &
                    renormalize_order_algorithm_in(im)) then
            !do nothing
           endif
          else
           print *,"is_rigid_local(sub_im) invalid"
           stop
          endif
         endif !im.ne.sub_im
        enddo !sub_im=1,num_materials

        do while (rank_algorithm(irank).ne.0)
         irank=irank+1
        enddo

        if (irank.gt.num_materials) then
         print *,"irank invalid: ",irank
         stop
        endif

        rank_algorithm(irank)=im
        nonzero_ranks=nonzero_ranks+1
        if (renormalize_order_algorithm_in(im).eq.irank) then
         !do nothing
        else
         print *,"renormalize_order_algorithm_in(im)<>irank"
         print *,"im=",im
         print *,"renormalize_order_algorithm_in(im) ", &
           renormalize_order_algorithm_in(im)
         print *,"renormalize_order_algorithm_in ", &
           renormalize_order_algorithm_in
         stop
        endif
       else
        print *,"is_rigid_local(im) invalid ",im,is_rigid_local(im)
        stop
       endif
      enddo ! im=1,num_materials

      call set_order_algorithm( &
        order_algorithm_in, &
        renormalize_order_algorithm_in)

      print *,"initializing geometry tables"

      call init_geometry_tables()
      if (1.eq.0) then
       call volume_sanity_check()
      endif
      if (1.eq.0) then
       if (geom_nthreads.eq.1) then
        sdim=3
        nmax_test=POLYGON_LIST_MAX
        call diagnostic_MOF(sdim,nmax_test)
        stop
       else
        print *,"expecting geom_nthreads.eq.1"
        stop
       endif
      endif

      return
      end subroutine fort_initmof

      subroutine delete_mof()

      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

#include "mofdata.H"

      deallocate(geom_xtetlist_local)
      deallocate(geom_xtetlist)
      deallocate(geom_xtetlist_old)
      deallocate(geom_xtetlist_uncapt)

      deallocate(intercept_error_history)

      deallocate(mof_calls)
      deallocate(mof_iterations)
      deallocate(mof_errors)

      return
      end subroutine delete_mof

end module MOF_cpp_module

