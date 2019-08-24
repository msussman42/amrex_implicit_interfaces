MODULE GeneralClass
USE probcommon_module
use MOF_routines_module
use geometry_intersect_module
IMPLICIT NONE


real(kind=8),parameter :: eps = 1e-10
real(kind=8),parameter :: tol = 1e-6
real(kind=8),parameter :: pi  = 4.0d0 * atan(1.0d0)
integer, parameter :: dp = selected_real_kind(p=8)

TYPE::POINTS
    REAL(KIND=dp)       :: VAL(2)
END TYPE

TYPE:: CUBE
   REAL(KIND=dp)            :: vt(8,3)
END TYPE

TYPE:: tetrahedron
    REAL(KIND=dp),allocatable :: vt(:,:)   ! 4,3 in 3D  3,2 in 2D
    INTEGER                   :: id
    INTEGER,allocatable       :: id_bc(:,:)    ! (3,2) in 3D, (2,2) in 2D, 0 for not touching the wall
    INTEGER,allocatable       :: id_face(:)    ! 4 faces in a tet. (3D), 3 faces in a triangle (2D) 
    INTEGER,allocatable       :: id_node(:,:)  ! index 1: 1..8 in 3D, 1..4 in 2D  index 2: 1..nmat
    INTEGER,allocatable       :: is_corner(:)  ! 4 vertices in a tet. (3D), 3 vert for triangle (2D)
    INTEGER,allocatable       :: is_face(:)    ! 4 faces in a tet. (3D), 3 faces in a triangle (2D)
END TYPE

TYPE:: tetra_data_list
   integer                         :: num
   TYPE(tetrahedron),allocatable   :: tetra(:) 
END TYPE

TYPE:: tetra_linked_list
   integer                         :: freehead,acthead,freerear,actrear         
   integer                         :: freenum,activenum
   integer,allocatable             :: next(:),prev(:)                     ! allocate list   
END TYPE


contains

!------------construct tetrahedron linked list----------------
subroutine init_freelist(list,maxlen)
implicit none

TYPE(tetra_linked_list)      :: list
integer                      :: maxlen
integer                      :: i

 if (maxlen.lt.10) then
  print *,"maxlen too small"
  stop
 endif
 if (maxlen.gt.10000) then
  print *,"maxlen too big"
  stop
 endif
 list%freenum = 0
 list%activenum=0
 list%freehead = 1
 list%acthead = 0
 list%actrear=0
 allocate(list%next(maxlen),list%prev(maxlen)) 
 do i = 1,maxlen 
   list%next(i) = i+1
   list%prev(i) = i-1
 enddo
   list%next(maxlen) = 0
   list%freenum=maxlen
   list%freerear=maxlen
end subroutine init_freelist

! TODO
!  - RETURN ALL ACTIVE MEMBERS BACK TO FREE 
!  - INIT FVM STENCIL:
!      -MOFDATA STENCIL AND XSTENCIL as input
!        MOFDATA_STENCIL(-1:1,SDIM,nmat*ngeom_recon)
!        XSTENCIL(-3:3,SDIM)   half cell increments.
!      -Dimension is input
!      -INITIALIZE WITH 5 tet. or 2 triangles
!      -Intersect the "objects" with all the planes.
!   |    .   |   .   |   .   |
!  -3   -2  -1   0   1   2   3


!---------------- reset active list----------------------
subroutine list_purge(list,maxlen) 
implicit none

TYPE(tetra_linked_list)         :: list
integer                         :: maxlen
integer                         :: i

if(list%actnum .eq. 0)then
  print *, "no need to purge, active list is empty."
  stop
endif

do while (list%actnum.ne.0)

 i=list%actnum
 call list_del(list,i)

enddo

end subroutine list_purge

!------------delete tetrahedron linked list---------
subroutine list_void(list)
implicit none

TYPE(tetra_linked_list)      :: list

  list%freenum = 0
  list%activenum = 0
  list%freehead = 0
  list%acthead  = 0
  list%actrear=0
  list%freerear=0
  deallocate(list%next,list%prev)

end subroutine list_void
!------------------------------------------------
subroutine list_add(list,id)
implicit none

TYPE(tetra_linked_list)         :: list
integer                         :: id  ! this is OUTPUT
integer                         :: i
integer                         :: save_act

if(list%acthead .eq. 0) then
   ! SANITY CHECK
  if(list%actrear .ne. 0 .or. list%actnum .ne. 0) then  
    print *,"error, active list should be empty" 
    stop
  endif
  if(list%freehead .eq. 0 .or. list%freenum .eq. 0)then
    print *,"error, freelist should not be empty"
    stop
  endif

  list%acthead  = list%freehead
  list%freehead = list%next(list%freehead)
  list%prev(list%freehead) = 0
  list%next(list%acthead) = 0
  list%prev(list%acthead) = 0
  list%actrear=list%acthead
  id=list%acthead
else if (list%acthead.gt.0) then
  ! SANITY CHECKS
  if(list%actnum .eq. 0) then  
    print *,"error, active list should not be empty" 
    stop
  endif
  if(list%freehead .eq. 0 .or. list%freenum .eq. 0)then
    print *,"error, freelist should not be empty"
    stop
  endif

  save_act = list%acthead

  list%acthead  = list%freehead
  list%freehead = list%next(list%freehead)
  list%prev(list%freehead) = 0 

  list%next(list%acthead) = save_act
  list%prev(save_act) = list%acthead
  list%prev(list%acthead) = 0
  id=list%acthead
else
  print *,"acthead invalid"
  stop
endif

  list%activenum = list%activenum + 1
  list%freenum = list%freenum - 1

end subroutine 
!----------------------------------------------
subroutine list_del(list,id)
implicit none

TYPE(tetra_linked_list)         :: list
integer, INTENT(IN)             :: id ! INPUT
integer                         :: i
integer                         :: save_free
integer                         :: old_prev,old_next

  ! SANITY CHECKS
if(list%actnum .eq. 0)then
  print *, "the active list should not be empty"
endif


old_prev = list%prev(id)
old_next = list%next(id)
save_free  = list%freehead
list%freehead = id
list%next(list%freehead) = save_free
list%prev(save_free) = list%freehead
list%prev(list%freehead) = 0

list%next(old_prev) = old_next
list%prev(old_next) = old_prev

list%activenum = list%activenum - 1
list%freenum = list%freenum + 1
if (id.eq.list%actrear) then
 list%actrear=old_prev
endif

end subroutine list_del

!--------------------------------------------------
!  - INIT FVM STENCIL:
!      -MOFDATA STENCIL AND XSTENCIL as input
!        MOFDATA_STENCIL(-1:1,SDIM,nmat*ngeom_recon)
!        XSTENCIL(-3:3,SDIM)   half cell increments.
!      -Dimension is input
!      -INITIALIZE WITH 5 tet. or 2 triangles
!      -Intersect the "objects" with all the planes.
!   |    .   |   .   |   .   |
!  -3   -2  -1   0   1   2   3
!
subroutine init_fvm_stencil(sdim,nmat,dx,mofdata,xstencil)
implicit none

integer,intent(in)          :: sdim,nmat
real(kind=dp),intent(in)    :: mofdata(-1:1,sdim,nmat*ngeom_recon)
                             ! ngeom_recon = 3+sdim*2
real(kind=dp),intent(in)    :: xstencil(-3:3,sdim)
real(kind=dp),intent(in)    :: dx(sdim)

real(kind=dp)               :: xnode(4*(sdim-1),sdim)
type(tetrahedron)           :: cc  !(center cell) 
type(tetra_data_list)       :: tetras 
integer                     :: i,j,k,i1,j1


IF(SDIM .EQ. 2) THEN
  CALL init_tetras(sdim,dx,xstencil,tetras)
  
  call intersect_tet(x1,x2,xtetlist_old,xtetlist,nlist,nmax,sdim) 


ELSEIF(SDIM .EQ. 3) THEN  
  CALL init_tetras(sdim,dx,xstencil,tetras)
  
  call intersect_tet(x1,x2,xtetlist_old,xtetlist,nlist,nmax,sdim) 

ELSE
  print*, "error, dimension must either be 2 or 3"
ENDIF




deallocate(tetras%tetra)

end subroutine init_fvm_stencil

!-------------------------------------------------------
subroutine alloc_tetra(tetra,sdim)
implicit none

type(tetrahedron)         :: tetra
integer                   :: sdim

!TYPE:: tetrahedron
!    REAL(KIND=dp),allocatable :: vt(:,:)   ! 4,3 in 3D  3,2 in 2D
!    INTEGER                   :: id
!    INTEGER,allocatable       :: id_bc(:,:)    ! (3,2) in 3D, (2,2) in 2D, 0 for not touching the wall
!    INTEGER,allocatable       :: id_face(:)    ! 4 faces in a tet. (3D), 3 faces in a triangle (2D) 
!    INTEGER,allocatable       :: id_node(:,:)  ! index 1: 1..8 in 3D, 1..4 in 2D  index 2: 1..nmat
!    INTEGER,allocatable       :: is_corner(:)  ! 4 vertices in a tet. (3D), 3 vert for triangle (2D)
!    INTEGER,allocatable       :: is_face(:)    ! 4 faces in a tet. (3D), 3 faces in a triangle (2D)
!END TYPE

if(sdim .eq. 2)then

  allocate(vt(3,2))
  allocate(id_bc(2,2))
  allocate(id_face(3))
  allocate(id_node(3))    !?????
  allocate(is_corner(3))  
  allocate(is_face(3))

elseif(sdim .eq. 3)then

  allocate(vt(4,3))
  allocate(id_bc(3,2))
  allocate(id_face(4))
  allocate(id_node(4))    !?????
  allocate(is_corner(4))  
  allocate(is_face(4))


else
  print *,"sdim invalid  in alloca_tetra"
endif



end subroutine alloc_tetra

!------------------------------------------------------
subroutine init_tetras(sdim,dx,xstencil,tetras)
implicit none

integer       ,intent(in)    :: sdim
real(kind=dp) ,intent(in)    :: xstencil(-3:3,sdim)
real(kind=dp) ,intent(in)    :: dx(sdim)
type(tetra_data_list)        :: xtets

!real(kind=dp),allocatable    :: xnode(4*(sdim-1),sdim) 
real(kind=dp)               :: xnode(4*(sdim-1),sdim)
integer                     :: nodelist(sdim+1)
integer                     :: id,i

real(kind=dp)               :: datanode(4*(sdim-1))
real(kind=dp)               :: datatet(sdim+1)

!TYPE:: tetra_data_list
!   integer                         :: num
!   TYPE(tetrahedron),allocatable   :: tetra(:) 
!END TYPE
!TYPE:: tetrahedron
!    REAL(KIND=dp),allocatable :: vt(:,:)   ! 4,3 in 3D  3,2 in 2D
!    INTEGER                   :: id
!    INTEGER,allocatable       :: id_bc(:,:)    ! (3,2) in 3D, (2,2) in 2D, 0 for not touching the wall
!    INTEGER,allocatable       :: id_face(:)    ! 4 faces in a tet. (3D), 3 faces in a triangle (2D) 
!    INTEGER,allocatable       :: id_node(:,:)  ! index 1: 1..8 in 3D, 1..4 in 2D  index 2: 1..nmat
!    INTEGER,allocatable       :: is_corner(:)  ! 4 vertices in a tet. (3D), 3 vert for triangle (2D)
!    INTEGER,allocatable       :: is_face(:)    ! 4 faces in a tet. (3D), 3 faces in a triangle (2D)
!END TYPE


IF(sdim .eq. 2) THEN
  xnode(1,1) = xstencil(0,1) - 0.5d0*dx(1)
  xnode(1,2) = xstencil(0,2) - 0.5d0*dx(2)
  xnode(2,1) = xstencil(0,1) + 0.5d0*dx(1)
  xnode(2,2) = xstencil(0,2) - 0.5d0*dx(2)
  xnode(3,1) = xstencil(0,1) + 0.5d0*dx(1)
  xnode(3,2) = xstencil(0,2) + 0.5d0*dx(2)
  xnode(4,1) = xstencil(0,1) - 0.5d0*dx(1)
  xnode(4,2) = xstencil(0,2) + 0.5d0*dx(2)

  xtets%num = 2
  allocate(xtets%tetra(2))
  do  i =1 ,2
     CALL alloc_tetra(xtets%tetra(i),sdim)    
  enddo
  ! square side order
  ! 1: 12   2: 24   3:34   4: 31 
  !
  !  3   4
  !    \
  !  1   2
  ! 123 234
  ! 
  do i1 = 1,2
    call extract_tet(xnode,datanode,xtets%tetra(i1)%vt,datatet,i1,0,sdim)
        if(i1 .eq. 1) then
              xtets%tetra(i1)%is_corner(1) = 1
              xtets%tetra(i1)%is_corner(2) = 2
              xtets%tetra(i1)%is_corner(3) = 3
              xtets%tetra(i1)%is_face(1) = 1
              xtets%tetra(i1)%is_face(2) = 0
              xtets%tetra(i1)%is_face(3) = 4
        elseif(i1 .eq. 2) then
              xtets%tetra(i1)%is_corner(1) = 2
              xtets%tetra(i1)%is_corner(2) = 3
              xtets%tetra(i1)%is_corner(3) = 4
              xtets%tetra(i1)%is_face(1) = 0
              xtets%tetra(i1)%is_face(2) = 3
              xtets%tetra(i1)%is_face(3) = 2
          
        else
          print *,"error, wrong number of tetra in init_tetra"
        endif
  enddo
   
 
ELSEIF(sdim .eq. 3) THEN
  xnode(1,1) = xstencil(0,1) - 0.5d0*dx(1)
  xnode(1,2) = xstencil(0,2) - 0.5d0*dx(2)
  xnode(1,3) = xstencil(0,3) - 0.5d0*dx(3)
  xnode(2,1) = xstencil(0,1) + 0.5d0*dx(1)
  xnode(2,2) = xstencil(0,2) - 0.5d0*dx(2)
  xnode(2,3) = xstencil(0,3) - 0.5d0*dx(3)
  xnode(3,1) = xstencil(0,1) - 0.5d0*dx(1)
  xnode(3,2) = xstencil(0,2) - 0.5d0*dx(2)
  xnode(3,3) = xstencil(0,3) + 0.5d0*dx(3)
  xnode(4,1) = xstencil(0,1) + 0.5d0*dx(1)
  xnode(4,2) = xstencil(0,2) - 0.5d0*dx(2)
  xnode(4,3) = xstencil(0,3) + 0.5d0*dx(3)
  xnode(5,1) = xstencil(0,1) - 0.5d0*dx(1)
  xnode(5,2) = xstencil(0,2) + 0.5d0*dx(2)
  xnode(5,3) = xstencil(0,3) - 0.5d0*dx(3)
  xnode(6,1) = xstencil(0,1) + 0.5d0*dx(1)
  xnode(6,2) = xstencil(0,2) + 0.5d0*dx(2)
  xnode(6,3) = xstencil(0,3) - 0.5d0*dx(3)
  xnode(7,1) = xstencil(0,1) - 0.5d0*dx(1)
  xnode(7,2) = xstencil(0,2) + 0.5d0*dx(2)
  xnode(7,3) = xstencil(0,3) + 0.5d0*dx(3)
  xnode(8,1) = xstencil(0,1) + 0.5d0*dx(1)
  xnode(8,2) = xstencil(0,2) + 0.5d0*dx(2)
  xnode(8,3) = xstencil(0,3) + 0.5d0*dx(3)

  xtets%num = 5
  allocate(xtets%tetra(5))
  do  i =1,5
     CALL alloc_tetra(xtets%tetra(i),sdim)    
  enddo

  do i1 = 1,5
    call extract_tet(xnode,datanode,xtets%tetra(i1)%vt,datatet,i1,0,sdim)

      ! cube face order
      ! 1:left  2:right  3: down 4: up 5: back 6: front
      ! tetra face order
      ! 123 124 134 234
      !            
           if(i1 .eq. 1) then
              xtets%tetra(i1)%is_corner(1) = 3
              xtets%tetra(i1)%is_corner(2) = 1
              xtets%tetra(i1)%is_corner(3) = 5
              xtets%tetra(i1)%is_corner(4) = 2
              xtets%tetra(i1)%is_face(1) = 1
              xtets%tetra(i1)%is_face(2) = 5
              xtets%tetra(i1)%is_face(3) = 0
              xtets%tetra(i1)%is_face(4) = 3
           elseif(i1 .eq. 2) then
              xtets%tetra(i1)%is_corner(1) = 5
              xtets%tetra(i1)%is_corner(2) = 6
              xtets%tetra(i1)%is_corner(3) = 2
              xtets%tetra(i1)%is_corner(4) = 8 
              xtets%tetra(i1)%is_face(1) = 3
              xtets%tetra(i1)%is_face(2) = 6
              xtets%tetra(i1)%is_face(3) = 0
              xtets%tetra(i1)%is_face(4) = 2    
           elseif(i1 .eq. 3) then
              xtets%tetra(i1)%is_corner(1) = 3
              xtets%tetra(i1)%is_corner(2) = 5
              xtets%tetra(i1)%is_corner(3) = 7
              xtets%tetra(i1)%is_corner(4) = 8
              xtets%tetra(i1)%is_face(1) = 1
              xtets%tetra(i1)%is_face(2) = 0
              xtets%tetra(i1)%is_face(3) = 4
              xtets%tetra(i1)%is_face(4) = 6
           elseif(i1 .eq. 4) then
              xtets%tetra(i1)%is_corner(1) = 3
              xtets%tetra(i1)%is_corner(2) = 4
              xtets%tetra(i1)%is_corner(3) = 2
              xtets%tetra(i1)%is_corner(4) = 8
              xtets%tetra(i1)%is_face(1) = 5
              xtets%tetra(i1)%is_face(2) = 4
              xtets%tetra(i1)%is_face(3) = 0
              xtets%tetra(i1)%is_face(4) = 2
           elseif(i1 .eq. 5) then
              xtets%tetra(i1)%is_corner(1) = 5
              xtets%tetra(i1)%is_corner(2) = 8
              xtets%tetra(i1)%is_corner(3) = 3
              xtets%tetra(i1)%is_corner(4) = 2
              xtets%tetra(i1)%is_face = 0
           else
              print *,"error, wrong number of tetra in init_tetra"
           
           endif
    
  enddo
  

ELSE
  PRINT *," WRONG DIMENSION"

ENDIF


end subroutine init_tetras




















!-----------------------------------------------
!-------------------------------------------------
subroutine tetra_copy(tetra_a,tetra_b)
implicit none

TYPE(tetrahedron),intent(in)   :: tetra_a
TYPE(tetrahedron)              :: tetra_b
integer                        :: i

do i = 1,4
  tetra_b%nodes(i)%val = tetra_a%nodes(i)%val
enddo

tetra_b%id = tetra_a%id

tetra_b%id_bc = tetra_a%id_bc

tetra_b%id_face = tetra_a%id_face

tetra_b%id_nodenum = tetra_a%id_nodenum

allocate(tetra_b%id_node(tetra_b%id_nodenum))

do i = 1,tetra_b%id_nodenum
   tetra_b%id_node(i) = tetra_a%id_node(i) 
enddo

tetra_b%is_corner = tetra_a%is_corner

tetra_b%is_face = tetra_a%is_face


end subroutine tetra_copy

!--------------------------------------------------
subroutine tetra_del(tetra_b)
implicit none

TYPE(tetrahedron)              :: tetra_b
integer                        :: i


do i = 1,4
  tetra_b%nodes(i)%val = 0.0d0
enddo


tetra_b%id = 0

tetra_b%is_boundary_tag = 0

tetra_b%id_bc = 0

tetra_b%id_face = 0

tetra_b%id_nodenum = 0

deallocate(tetra_b%id_node)

tetra_b%is_corner = 0

tetra_b%is_face = 0

end subroutine tetra_del
!-------------------------------------------------------



subroutine init_partition(sdim,nmat,lo,hi,numofpart,dx,datanode, &
                          nodx,nody,nodz,corx,cory,corz, &
                          dx,tetras)
implicit none

real(kind=dp),intent(in)            :: lo(sdim), hi(sdim)
integer,intent(in)                  :: sdim,nmat
integer,intent(in)                  :: numofpart(sdim)

integer                             :: N,M,L
integer                             :: i,j,k,i1,j1
real(kind=dp)                       :: dx(sdim)
real(kind=dp),allocatable           :: nodx(:),nody(:),nodz(:)
real(kind=dp),allocatable           :: corx(:),cory(:),corz(:)
TYPE(CUBE),   allocatable           :: CELL_3d(:,:,:)
type(tetra_list),allocatable        :: tetras(:,:,:)


REAL(kind=dp) :: xtet(sdim+1,sdim)
REAL(kind=dp) :: datatet(sdim+1)
!INTEGER       :: symmetry_flag   
REAL(kind=dp) :: datanode(4*(sdim-1))
REAL(kind=dp) :: xnode(4*(sdim-1),sdim)
INTEGER       :: nodelist(sdim+1)


 if(sdim .eq. 2) then
   N = numofpart(1)
   M = numofpart(2)
   dx(1) = (hi(1)-lo(1))/N
   dx(2) = (hi(2)-lo(2))/M   
   allocate(nodx(N))
   allocate(nody(M))
   allocate(corx(N+1))
   allocate(cory(M+1))
   do i = 1,N
      nodx(i) = 0.5_dp*dx(1) + (i-1)* dx(1)
   enddo
   do i = 1,M
      nody(i) = 0.5_dp*dx(2) + (i-1)* dx(2)
   enddo
   do i = 1,N+1
      corx(i) = (i-1)* dx(1)
   enddo
   do i = 1,M+1
      cory(i) = (i-1)* dx(2)
   enddo

 elseif(sdim .eq. 3)then
   N = numofpart(1)
   M = numofpart(2)
   L = numofpart(3)
   dx(1) = (hi(1)-lo(1))/N
   dx(2) = (hi(2)-lo(2))/M   
   dx(3) = (hi(3)-lo(3))/L 
   allocate(nodx(N))                             ! remember to deallocate
   allocate(nody(M))
   allocate(nodz(L))
   allocate(corx(N+1))
   allocate(cory(M+1))
   allocate(corz(L+1))
   do i = 1,N
      nodx(i) = 0.5_dp*dx(1) + (i-1)* dx(1)
   enddo
   do i = 1,M
      nody(i) = 0.5_dp*dx(2) + (i-1)* dx(2)
   enddo
   do i = 1,L
      nodz(i) = 0.5_dp*dx(3) + (i-1)* dx(3)
   enddo

   do i = 1,N+1
      corx(i) = (i-1)* dx(1)
   enddo
   do i = 1,M+1
      cory(i) = (i-1)* dx(2)
   enddo
   do i = 1,L+1
      corz(i) = (i-1)* dx(3)
   enddo
 else
   write(*,*) "sdim is invalid"
 endif


if (sdim .eq. 3) then
  allocate(cell_3d(N,M,L))                 ! deallocate
  do i = 1,N
    do j = 1,M
      do k = 1,L
   
           cell_3d(i,j,k)%vt(1,1) = corx(i)              
           cell_3d(i,j,k)%vt(1,2) = cory(j)
           cell_3d(i,j,k)%vt(1,3) = corz(k)

           cell_3d(i,j,k)%vt(2,1) = corx(i+1)              
           cell_3d(i,j,k)%vt(2,2) = cory(j)
           cell_3d(i,j,k)%vt(2,3) = corz(k)         

           cell_3d(i,j,k)%vt(3,1) = corx(i)              
           cell_3d(i,j,k)%vt(3,2) = cory(j+1)
           cell_3d(i,j,k)%vt(3,3) = corz(k) 

           cell_3d(i,j,k)%vt(4,1) = corx(i+1)              
           cell_3d(i,j,k)%vt(4,2) = cory(j+1)
           cell_3d(i,j,k)%vt(4,3) = corz(k) 

           cell_3d(i,j,k)%vt(5,1) = corx(i)              
           cell_3d(i,j,k)%vt(5,2) = cory(j)
           cell_3d(i,j,k)%vt(5,3) = corz(k+1) 

           cell_3d(i,j,k)%vt(6,1) = corx(i+1)              
           cell_3d(i,j,k)%vt(6,2) = cory(j)
           cell_3d(i,j,k)%vt(6,3) = corz(k+1)

           cell_3d(i,j,k)%vt(7,1) = corx(i)              
           cell_3d(i,j,k)%vt(7,2) = cory(j+1)
           cell_3d(i,j,k)%vt(7,3) = corz(k+1) 

           cell_3d(i,j,k)%vt(8,1) = corx(i+1)              
           cell_3d(i,j,k)%vt(8,2) = cory(j+1)
           cell_3d(i,j,k)%vt(8,3) = corz(k+1)  
      enddo
    enddo
  enddo

elseif(sdim .eq. 2)then





else
  write(*,*) "sdim invalid"
endif


if(sdim .eq. 3) then
  allocate (tetras(N,M,L))
  do i = 1,N
   do j = 1,M
     do k = 1,L
        tetras(i,j,k)%num = 5
        allocate(tetras(i,j,k)%tetra(5))
        do i1 = 1,5
! cube face order
! 1:left  2:right  3: down 4: up 5: back 6: front
! 
! 
! tetra face order
! 123 124 134 234
!            
           call extract_tet(cell_3d(i,j,k)%vt,datanode,xtet,datatet,i1,0,sdim)
           tetras(i,j,k)%tetra(i1)%vt  = xtet                                  ! init vertex
           tetras(i,j,k)%tetra(i1)%id  = 0                                     ! init id           
           tetras(i,j,k)%tetra(i1)%id_bc = 0             
           tetras(i,j,k)%tetra(i1)%id_face = 0
           allocate(tetras(i,j,k)%tetra(i1)%id_node(4*nmat))               ! remember to deallocate
           tetras(i,j,k)%tetra(i1)%id_node = 0
       
           if(i1 .eq. 1) then
              tetras(i,j,k)%tetra(i1)%is_corner(1) = 3
              tetras(i,j,k)%tetra(i1)%is_corner(2) = 1
              tetras(i,j,k)%tetra(i1)%is_corner(3) = 5
              tetras(i,j,k)%tetra(i1)%is_corner(4) = 2
              tetras(i,j,k)%tetra(i1)%is_face(1) = 1
              tetras(i,j,k)%tetra(i1)%is_face(2) = 5
              tetras(i,j,k)%tetra(i1)%is_face(3) = 0
              tetras(i,j,k)%tetra(i1)%is_face(4) = 3
           elseif(i1 .eq. 2) then
              tetras(i,j,k)%tetra(i1)%is_corner(1) = 5
              tetras(i,j,k)%tetra(i1)%is_corner(2) = 6
              tetras(i,j,k)%tetra(i1)%is_corner(3) = 2
              tetras(i,j,k)%tetra(i1)%is_corner(4) = 8 
              tetras(i,j,k)%tetra(i1)%is_face(1) = 3
              tetras(i,j,k)%tetra(i1)%is_face(2) = 6
              tetras(i,j,k)%tetra(i1)%is_face(3) = 0
              tetras(i,j,k)%tetra(i1)%is_face(4) = 2    
           elseif(i1 .eq. 3) then
              tetras(i,j,k)%tetra(i1)%is_corner(1) = 3
              tetras(i,j,k)%tetra(i1)%is_corner(2) = 5
              tetras(i,j,k)%tetra(i1)%is_corner(3) = 7
              tetras(i,j,k)%tetra(i1)%is_corner(4) = 8
              tetras(i,j,k)%tetra(i1)%is_face(1) = 1
              tetras(i,j,k)%tetra(i1)%is_face(2) = 0
              tetras(i,j,k)%tetra(i1)%is_face(3) = 4
              tetras(i,j,k)%tetra(i1)%is_face(4) = 6
           elseif(i1 .eq. 4) then
              tetras(i,j,k)%tetra(i1)%is_corner(1) = 3
              tetras(i,j,k)%tetra(i1)%is_corner(2) = 4
              tetras(i,j,k)%tetra(i1)%is_corner(3) = 2
              tetras(i,j,k)%tetra(i1)%is_corner(4) = 8
              tetras(i,j,k)%tetra(i1)%is_face(1) = 5
              tetras(i,j,k)%tetra(i1)%is_face(2) = 4
              tetras(i,j,k)%tetra(i1)%is_face(3) = 0
              tetras(i,j,k)%tetra(i1)%is_face(4) = 2
           elseif(i1 .eq. 5) then
              tetras(i,j,k)%tetra(i1)%is_corner(1) = 5
              tetras(i,j,k)%tetra(i1)%is_corner(2) = 8
              tetras(i,j,k)%tetra(i1)%is_corner(3) = 3
              tetras(i,j,k)%tetra(i1)%is_corner(4) = 2
              tetras(i,j,k)%tetra(i1)%is_face = 0
           endif

        enddo

     enddo
   enddo
  enddo

endif







end subroutine init_partition










end module




program superthing_3d
use generalclass
USE probcommon_module
use MOF_routines_module
use geometry_intersect_module
implicit none

integer,parameter                   :: N=100,M =100, L=100
integer                             :: sdim,nmat
integer                             :: numofpart(3)
real(kind=dp)                       :: lo(3), hi(3)
integer                             :: i,j,k,i1,j1
real(kind=dp)                       :: dx(3)
real(kind=dp)                       :: nodx(N),nody(M),nodz(L)
real(kind=dp)                       :: corx(N+1),cory(M+1),corz(L+1)
type(tetra_list)                    :: init_tetras(N,M,L)

sdim = 3
nmat = 2
lo = 0.0d0
hi = 1.0d0

numofpart(1) = N
numofpart(2) = M
numofpart(3) = L


call  init_partition(sdim,nmat,lo,hi,numofpart,dx,datanode, &
                          nodx,nody,nodz,corx,cory,corz, &
                          dx,init_tetras)







end program superthing_3d































