! --------------------------------------------------------------------
!  Flow Physics and Computation Lab
!   
!
!  FSI module for
!  VICAR3D, a viscous, Cartesian, 3D flow solver.
!
!  This is a contineously developing project.
!
!  Starting Developers:
!  Kourosh Shoele
!
!
!  Final Filename: UTIL_BOUNDARY_FORCE_FSI.F90
!  Latest Modification: June, 01 2016 
!  by Kourosh Shoele
! --------------------------------------------------------------------
#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "AMReX_SPACE.H"

#ifdef BL_USE_MPI
#define mpi_activate 1
#else
#ifdef BL_USE_MPI3
#define mpi_activate 1
#else
#define mpi_activate 0
#endif
#endif

# define L2GI(i)       myIs+i-1
# define L2GJ(j)       myJs+j-1
# define G2LI(i)       i-(myIs-1)
# define G2LJ(j)       j-(myJs-1)


MODULE dummy_module
implicit none 
include './distFSI/grid_def'
 integer:: Nout, nsecIBM, nsecIBMmax,NrIBM
 integer:: nIBM_fib, nrIBM_fib
 integer:: nqIBM_fsh, nIBM_fsh, nrIBM_fsh
 integer:: nIBM_esh, nrIBM_esh
 integer:: nIBM_fbc, nrIBM_fbc


 integer,SAVE :: dtypeDelta(3) 
 real*8,SAVE :: fx, fy, fz 
 real*8,SAVE :: dsecIBMy,dsecIBMz
 integer, SAVE, allocatable, dimension(:) :: nIBM_r,nIBM_q 
 integer, SAVE, allocatable, dimension(:) :: nIBM_r_fib 
 integer, SAVE, allocatable, dimension(:) :: nIBM_r_fsh,nIBM_q_fsh 
 integer, SAVE, allocatable, dimension(:) :: nIBM_r_esh 
 integer, SAVE, allocatable, dimension(:) :: nIBM_r_fbc 

 real*8, SAVE, allocatable, dimension(:,:,:) :: coord_fib
 real*8, SAVE, allocatable, dimension(:,:,:,:) :: coord_fsh
 real*8, SAVE, allocatable, dimension(:,:,:)   :: coord_esh
 real*8, SAVE, allocatable, dimension(:,:,:)   :: coord_fbc

 real*8, SAVE, allocatable, dimension(:,:,:) :: coord_fib_prev
 real*8, SAVE, allocatable, dimension(:,:,:,:) :: coord_fsh_prev
 real*8, SAVE, allocatable, dimension(:,:,:)   :: coord_esh_prev
 real*8, SAVE, allocatable, dimension(:,:,:)   :: coord_fbc_prev

 real*8, SAVE, allocatable, dimension(:,:,:,:) :: vel_fib, force_fib
 real*8, SAVE, allocatable, dimension(:,:,:,:) :: vel_fsh, force_fsh
 real*8, SAVE, allocatable, dimension(:,:,:)   :: vel_esh, force_esh
 real*8, SAVE, allocatable, dimension(:,:,:)   :: vel_fbc, force_fbc

 real*8, SAVE, allocatable, dimension(:,:,:,:) :: vel_fib_prev, force_fib_prev
 real*8, SAVE, allocatable, dimension(:,:,:,:) :: vel_fsh_prev, force_fsh_prev
 real*8, SAVE, allocatable, dimension(:,:,:)   :: vel_esh_prev, force_esh_prev
 real*8, SAVE, allocatable, dimension(:,:,:)   :: vel_fbc_prev, force_fbc_prev


 real*8, SAVE, allocatable, dimension(:,:)   :: ds_fib
 real*8, SAVE, allocatable, dimension(:,:,:) :: ds_fsh
 real*8, SAVE, allocatable, dimension(:,:)   :: ds_esh
 real*8, SAVE, allocatable, dimension(:,:)   :: ds_fbc

 real*8, SAVE, allocatable, dimension(:,:)   :: ds_fib_prev
 real*8, SAVE, allocatable, dimension(:,:,:) :: ds_fsh_prev
 real*8, SAVE, allocatable, dimension(:,:)   :: ds_esh_prev
 real*8, SAVE, allocatable, dimension(:,:)   :: ds_fbc_prev


 parameter(nqIBM_fsh=nq_IBM_fsh)

 parameter(nIBM_fib=ns_IBM_fib) ! Max number of points in each fiber
 parameter(nIBM_fsh=ns_IBM_fsh)
 parameter(nIBM_esh=ns_IBM_esh)
 parameter(nIBM_fbc=ns_IBM_fbc)

 parameter(nrIBM=nr_IBM)         ! Number of solid geometries
 parameter(nrIBM_fib=nr_IBM_fib) ! Number of fiber geometries
 parameter(nrIBM_fsh=nr_IBM_fsh)
 parameter(nrIBM_esh=nr_IBM_esh)
 parameter(nrIBM_fbc=nr_IBM_fbc)

 parameter(nsecIBMmax=Nsec_IBMmax)

END MODULE dummy_module


! interpolate velocities on the body marker
! inputs:
!   bodyIBM(nIBM,3) : 
!   nIBM : number of body marker points
!   dtype : delta function type
! output:
!   velIBM(nIBM,3) : velocities on the body marker points

! integer ndim1 : dimension, 2 or 3
! REAL*8  bodyIBM(nrIBM,nIBM,3) : body marker coordinates
! integer dtype(3) : delta function type for each direction
! REAL*8  dsecIBMy : section info (?)
! REAL*8  dsecIBMz : section info (?)
! integer nrIBM : number of fiber bodies(=nr_IBM_fib)
! integer nsecIBMmax : number of sections (?) (=Nsec_IBMmax)
! integer nsecIBM : section info (?)
! integer nIBM : max number of points in each body (=ns_ibm_fib)
! integer nIBM_r(nrIBM) : number of points in each fiber body



SUBROUTINE cal_velIBM_fib(&
 ndim1,&  
 bodyIBM,&      
 velIBM,&
 dtype,&
 dsecIBMy,&
 dsecIBMz,&
 nrIBM,&
 nsecIBMmax,&
 nsecIBM,&
 nIBM,&
 nIBM_r)


 USE global_parameters
 USE flow_parameters
 USE unstructured_surface_arrays 
 USE boundary_arrays
 USE flow_arrays
 USE grid_arrays
#if (mpi_activate==1)
 USE mpi
#endif
 
 implicit NONE
 integer ndim1
 REAL*8  bodyIBM(nrIBM,nIBM,3)
 REAL*8  velIBM(nrIBM,nsecIBMmax,nIBM,3)
 integer dtype(3)
 REAL*8  dsecIBMy
 REAL*8  dsecIBMz
 integer nrIBM
 integer nsecIBMmax
 integer nsecIBM
 integer nIBM
 integer nIBM_r(nrIBM)

 integer i, j, iCell, jCell, kCell, ii, jj, kk, ierr1,jq
 integer INTP_CORONA, ip, jp, kp
 REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: tempvel
 REAL*8 xp, yp, zp, gx, gy, gz, dfx, dfy, dfz
 integer imin, imax, jmin, jmax, kmin, kmax
 integer, allocatable, dimension (:,:) :: iCellMarker, jCellMarker, kCellMarker
 integer isecIBM,isec
 REAL*8  zsecIBM,ysecIBM

 

 
! print*, 'inside cal_velIBM'
!print*, dsecIBMy,dsecIBMz,nrIBM,nsecIBMmax,nsecIBM
 ALLOCATE(iCellMarker(nsecIBM,nIBM))
 ALLOCATE(jCellMarker(nsecIBM,nIBM))
 ALLOCATE(kCellMarker(nsecIBM,nIBM))
 
 ALLOCATE(tempvel(nrIBM,nsecIBMmax,nIBM))
  
 INTP_CORONA = 1
 
 velIBM = 0.0
 tempvel = 0.0



 ! find grid index
Do i = 1,nrIBM

  do isec=1,nsecIBM
   ysecIBM=real(isec*2-nsecIBM-1)/2.0*dsecIBMy
   zsecIBM=real(isec*2-nsecIBM-1)/2.0*dsecIBMz
! print*, 'Start with cell markers', nIBM_r(i),i
   Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)+ysecIBM
   zp = bodyIBM(i,j,3)+zsecIBM

   iCellMarker(isec,j) = -255
   jCellMarker(isec,j) = -255
   kCellMarker(isec,j) = -255
   
   select case ( ndim1 )
   
   case ( 2 )
   
   ! check the body marker is inside the current sub-domain
   IF( xp .ge. x(myIS) .and. xp .lt. x(myIe+1) .and. yp .ge. y(myJs) .and. yp .lt. y(myJe+1) ) THEN
   
   !-----+---------+---------
   !     |    x    |             
   !-----+---------+---------
   !    ip-1       ip
   !         iCell = ip
   
   ! ip, jp, kp are local index
   
   do ip = 1, nxc+1
    if( xp .ge. xc(L2GI(ip-1)) .and. xp .le. xc(L2GI(ip)))  iCellMarker(isec,j) = ip
   enddo
   
   do jp = 1, nyc+1
    if( yp .ge. yc(L2GJ(jp-1)) .and. yp .le. yc(L2GJ(jp)) ) jCellMarker(isec,j) = jp
   enddo   
   
   kCellMarker(isec,j) = 1
   
   ENDIF ! xp
   
   case ( 3 )
    
   IF( xp .ge. x(myIS) .and. xp .lt. x(myIe+1) .and. yp .ge. y(myJs) .and. yp .lt. y(myJe+1) .and. zp .ge. z(1) .and. zp .le. z(nz) ) THEN
   
   do ip = 1, nxc+1
    if( xp .ge. xc(L2GI(ip-1)) .and. xp .le. xc(L2GI(ip))) iCellMarker(isec,j) = ip
   enddo
   
   do jp = 1, nyc+1
    if( yp .ge. yc(L2GJ(jp-1)) .and. yp .le. yc(L2GJ(jp)) ) jCellMarker(isec,j) = jp
   enddo
   
   do kp = 1, nzc+1
    if( zp .ge. zc(kp-1) .and. zp .le. zc(kp) ) kCellMarker(isec,j) = kp
   enddo
   
   ENDIF ! xp
    
   end select
  ENDDO  !j
  ENDDO  !isection

  do isec=1,nsecibm
   ysecIBM=real(isec*2-nsecIBM-1)/2.0*dsecIBMy
   zsecIBM=real(isec*2-nsecIBM-1)/2.0*dsecIBMz
  Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)+ysecIBM
   zp = bodyIBM(i,j,3)+zsecIBM
   
   iCell = iCellMarker(isec,j)
   jCell = jCellMarker(isec,j)
   kCell = kCellMarker(isec,j)
    
   IF( iCell .ne. -255 .and. jCell .ne. -255 .and. kCell .ne. -255 ) THEN
   
   imin = max(iCell-1-INTP_CORONA,1-ngl)
   imax = min(iCell  +INTP_CORONA,nxc+ngl)
   jmin = max(jCell-1-INTP_CORONA,1-ngl)
   jmax = min(jCell  +INTP_CORONA,nyc+ngl)
   kmin = max(kCell-1-INTP_CORONA,0)
   kmax = min(kCell  +INTP_CORONA,nzc+1)
   
   select case ( ndim1 )
   
   case ( 2 )
  
   do ii=imin,imax
   do jj=jmin,jmax
   
   kk = kCell
   
   gx = (xp-xc(L2GI(ii)))/dxc(L2GI(iCell))
   gy = (yp-yc(L2GJ(jj)))/dyc(L2GJ(jCell))
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
 
   velIBM(i,isec,j,1) = velIBM(i,isec,j,1) + dfx*dfy*u(ii,jj,kk)
   velIBM(i,isec,j,2) = velIBM(i,isec,j,2) + dfx*dfy*v(ii,jj,kk)
    
   enddo
   enddo
   
   
   case ( 3 )
   
    
   do ii=imin, imax
   do jj=jmin, jmax
   do kk=kmin, kmax

   gx = (xp-xc(L2GI(ii)))/dxc(L2GI(iCell))
   gy = (yp-yc(L2GJ(jj)))/dyc(L2GJ(jCell))
   gz = (zp-zc(kk)      )/dzc(kcell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
   call deltao_fun(dtype(3),gz,dfz)
 
   velIBM(i,isec,j,1) = velIBM(i,isec,j,1) + dfx*dfy*dfz*u(ii,jj,kk)
   velIBM(i,isec,j,2) = velIBM(i,isec,j,2) + dfx*dfy*dfz*v(ii,jj,kk)   
   velIBM(i,isec,j,3) = velIBM(i,isec,j,3) + dfx*dfy*dfz*w(ii,jj,kk)
  
   enddo
   enddo
   enddo
   
   end select 
    
  ENDIF ! iCell
  
  ENDDO ! nIBM
  ENDDO ! nsecIBM
  ENDDO ! nrIBM
 
 tempvel=velIBM(:,:,:,1)
#if (mpi_activate==1)
 tempvel=0.0 
 CALL MPI_ALLREDUCE(velIBM(:,:,:,1),tempvel,nIBM*nsecIBMmax*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
 velIBM(:,:,:,1) = tempvel

 tempvel=velIBM(:,:,:,2)
#if (mpi_activate==1)
 tempvel=0.0
 CALL MPI_ALLREDUCE(velIBM(:,:,:,2),tempvel,nIBM*nsecIBMmax*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
 velIBM(:,:,:,2) = tempvel

 if (ndim1 .gt. 2) then
  tempvel=velIBM(:,:,:,3)
#if (mpi_activate==1)
  tempvel=0.0
  CALL MPI_ALLREDUCE(velIBM(:,:,:,3),tempvel,nIBM*nsecIBMmax*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
  velIBM(:,:,:,3) = tempvel
  tempvel=0.0
 ENDIF !ndim1   
 DEALLOCATE(tempvel)
 DEALLOCATE(iCellMarker,jCellMarker,kCellMarker)
 
END SUBROUTINE
!====================================================================================
SUBROUTINE cal_velIBM_fsh(ndim1, bodyIBM, velIBM, dtype,nrIBM,nIBM,nqIBM,nIBM_r,nIBM_q)
 USE global_parameters
 USE flow_parameters
 USE unstructured_surface_arrays 
 USE boundary_arrays
 USE flow_arrays
 USE grid_arrays
#if (mpi_activate==1)
 USE mpi
#endif
 
 implicit NONE
 
  integer i, j, iCell, jCell, kCell, ii, jj, kk, ierr1,jq
  integer ndim1, INTP_CORONA, ip, jp, kp, dtype(3)
  integer nqIBM,nIBM, nrIBM, nIBM_r(nrIBM),nIBM_q(nrIBM)
  REAL*8 bodyIBM(nrIBM,nqIBM,nIBM,3), velIBM(nrIBM,nqIBM,nIBM,3)
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: tempvel
  REAL*8 xp, yp, zp, gx, gy, gz, dfx, dfy, dfz
  integer imin, imax, jmin, jmax, kmin, kmax
  integer, allocatable, dimension (:,:) :: iCellMarker, jCellMarker, kCellMarker

 ALLOCATE(iCellMarker(nqIBM,nIBM))
 ALLOCATE(jCellMarker(nqIBM,nIBM))
 ALLOCATE(kCellMarker(nqIBM,nIBM))
 
 ALLOCATE(tempvel(nrIBM,nqIBM,nIBM))
  
 INTP_CORONA = 1
 
 velIBM = 0.0
 tempvel = 0.0



 ! find grid index
Do i = 1,nrIBM

   Do jq= 1, nIBM_q(i)
   Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,jq,j,1)
   yp = bodyIBM(i,jq,j,2)
   zp = bodyIBM(i,jq,j,3)

   iCellMarker(jq,j) = -255
   jCellMarker(jq,j) = -255
   kCellMarker(jq,j) = -255
   
   select case ( ndim1 )
   
   case ( 2 )
   
   ! check the body marker is inside the current sub-domain
   IF( xp .ge. x(myIS) .and. xp .lt. x(myIe+1) .and. yp .ge. y(myJs) .and. yp .lt. y(myJe+1) ) THEN
   
   !-----+---------+---------
   !     |    x    |             
   !-----+---------+---------
   !    ip-1       ip
   !         iCell = ip
   
   ! ip, jp, kp are local index
   
   do ip = 1, nxc+1
    if( xp .ge. xc(L2GI(ip-1)) .and. xp .le. xc(L2GI(ip)))  iCellMarker(jq,j) = ip
   enddo
   
   do jp = 1, nyc+1
    if( yp .ge. yc(L2GJ(jp-1)) .and. yp .le. yc(L2GJ(jp)) ) jCellMarker(jq,j) = jp
   enddo   
   
   kCellMarker(jq,j) = 1
   
   ENDIF ! xp
   
   case ( 3 )
    
   IF( xp .ge. x(myIS) .and. xp .lt. x(myIe+1) .and. yp .ge. y(myJs) .and. yp .lt. y(myJe+1) .and. zp .ge. z(1) .and. zp .le. z(nz) ) THEN
   
   do ip = 1, nxc+1
    if( xp .ge. xc(L2GI(ip-1)) .and. xp .le. xc(L2GI(ip))) iCellMarker(jq,j) = ip
   enddo
   
   do jp = 1, nyc+1
    if( yp .ge. yc(L2GJ(jp-1)) .and. yp .le. yc(L2GJ(jp)) ) jCellMarker(jq,j) = jp
   enddo
   
   do kp = 1, nzc+1
    if( zp .ge. zc(kp-1) .and. zp .le. zc(kp) ) kCellMarker(jq,j) = kp
   enddo
   
   ENDIF ! xp
    
   end select
  ENDDO 
  ENDDO

  Do jq= 1, nIBM_q(i)
  Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,jq,j,1)
   yp = bodyIBM(i,jq,j,2)
   zp = bodyIBM(i,jq,j,3)
   
   iCell = iCellMarker(jq,j)
   jCell = jCellMarker(jq,j)
   kCell = kCellMarker(jq,j)
    
   IF( iCell .ne. -255 .and. jCell .ne. -255 .and. kCell .ne. -255 ) THEN
   
   imin = max(iCell-1-INTP_CORONA,1-ngl)
   imax = min(iCell  +INTP_CORONA,nxc+ngl)
   jmin = max(jCell-1-INTP_CORONA,1-ngl)
   jmax = min(jCell  +INTP_CORONA,nyc+ngl)
   kmin = max(kCell-1-INTP_CORONA,0)
   kmax = min(kCell  +INTP_CORONA,nzc+1)
   
   select case ( ndim1 )
   
   case ( 2 )
  
   do ii=imin,imax
   do jj=jmin,jmax
   
   kk = kCell
   
   gx = (xp-xc(L2GI(ii)))/dxc(L2GI(iCell))
   gy = (yp-yc(L2GJ(jj)))/dyc(L2GJ(jCell))
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
 
   velIBM(i,jq,j,1) = velIBM(i,jq,j,1) + dfx*dfy*u(ii,jj,kk)
   velIBM(i,jq,j,2) = velIBM(i,jq,j,2) + dfx*dfy*v(ii,jj,kk)
    
   enddo
   enddo
   
   
   case ( 3 )
   
   do ii=imin, imax
   do jj=jmin, jmax
   do kk=kmin, kmax

   gx = (xp-xc(L2GI(ii)))/dxc(L2GI(iCell))
   gy = (yp-yc(L2GJ(jj)))/dyc(L2GJ(jCell))
   gz = (zp-zc(kk)      )/dzc(kcell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
   call deltao_fun(dtype(3),gz,dfz)
 
   velIBM(i,jq,j,1) = velIBM(i,jq,j,1) + dfx*dfy*dfz*u(ii,jj,kk)
   velIBM(i,jq,j,2) = velIBM(i,jq,j,2) + dfx*dfy*dfz*v(ii,jj,kk)   
   velIBM(i,jq,j,3) = velIBM(i,jq,j,3) + dfx*dfy*dfz*w(ii,jj,kk)
  
   enddo
   enddo
   enddo
   
   end select 
    
  ENDIF ! iCell
  
  ENDDO ! nIBM
  ENDDO ! nqIBM
  ENDDO ! nrIBM
  
 tempvel=velIBM(:,:,:,1)
#if (mpi_activate==1)
 tempvel=0.0 
 CALL MPI_ALLREDUCE(velIBM(:,:,:,1),tempvel,nIBM*nqIBM*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
 velIBM(:,:,:,1) = tempvel

 tempvel=velIBM(:,:,:,2)
#if (mpi_activate==1)
 tempvel=0.0
 CALL MPI_ALLREDUCE(velIBM(:,:,:,2),tempvel,nIBM*nqIBM*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
 velIBM(:,:,:,2) = tempvel

 if (ndim1 .gt. 2) then
  tempvel=velIBM(:,:,:,3)
#if (mpi_activate==1)
  tempvel=0.0
  CALL MPI_ALLREDUCE(velIBM(:,:,:,3),tempvel,nIBM*nqIBM*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
  velIBM(:,:,:,3) = tempvel
  tempvel=0.0
 ENDIF !ndim1   
 DEALLOCATE(tempvel)
 DEALLOCATE(iCellMarker,jCellMarker,kCellMarker)
 
END SUBROUTINE
!====================================================================================
SUBROUTINE cal_velIBM_esh(ndim1, bodyIBM, velIBM, dtype,nrIBM,nIBM,nIBM_r)
 USE global_parameters
 USE flow_parameters
 USE unstructured_surface_arrays 
 USE boundary_arrays
 USE flow_arrays
 USE grid_arrays
#if (mpi_activate==1)
 USE mpi
#endif
 
 implicit NONE
 
  integer i, j, iCell, jCell, kCell, ii, jj, kk, ierr1,jq
  integer ndim1, INTP_CORONA, ip, jp, kp, dtype(3)
  integer nIBM, nrIBM, nIBM_r(nrIBM)
  REAL*8 bodyIBM(nrIBM,nIBM,3), velIBM(nrIBM,nIBM,3)
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tempvel
  REAL*8 xp, yp, zp, gx, gy, gz, dfx, dfy, dfz
  integer imin, imax, jmin, jmax, kmin, kmax
  integer, allocatable, dimension (:) :: iCellMarker, jCellMarker, kCellMarker

 ALLOCATE(iCellMarker(nIBM))
 ALLOCATE(jCellMarker(nIBM))
 ALLOCATE(kCellMarker(nIBM))
 
 ALLOCATE(tempvel(nrIBM,nIBM))
  
 INTP_CORONA = 1
 
 velIBM = 0.0
 tempvel = 0.0



 ! find grid index
Do i = 1,nrIBM

   Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)
   zp = bodyIBM(i,j,3)

   iCellMarker(j) = -255
   jCellMarker(j) = -255
   kCellMarker(j) = -255
   
   select case ( ndim1 )
   
   case ( 2 )
   
   ! check the body marker is inside the current sub-domain
   IF( xp .ge. x(myIS) .and. xp .lt. x(myIe+1) .and. yp .ge. y(myJs) .and. yp .lt. y(myJe+1) ) THEN
   
   !-----+---------+---------
   !     |    x    |             
   !-----+---------+---------
   !    ip-1       ip
   !         iCell = ip
   
   ! ip, jp, kp are local index
   
   do ip = 1, nxc+1
    if( xp .ge. xc(L2GI(ip-1)) .and. xp .le. xc(L2GI(ip)))  iCellMarker(j) = ip
   enddo
   
   do jp = 1, nyc+1
    if( yp .ge. yc(L2GJ(jp-1)) .and. yp .le. yc(L2GJ(jp)) ) jCellMarker(j) = jp
   enddo   
   
   kCellMarker(j) = 1
   
   ENDIF ! xp
   
   case ( 3 )
    
   IF( xp .ge. x(myIS) .and. xp .lt. x(myIe+1) .and. yp .ge. y(myJs) .and. yp .lt. y(myJe+1) .and. zp .ge. z(1) .and. zp .le. z(nz) ) THEN
   
   do ip = 1, nxc+1
    if( xp .ge. xc(L2GI(ip-1)) .and. xp .le. xc(L2GI(ip))) iCellMarker(j) = ip
   enddo
   
   do jp = 1, nyc+1
    if( yp .ge. yc(L2GJ(jp-1)) .and. yp .le. yc(L2GJ(jp)) ) jCellMarker(j) = jp
   enddo
   
   do kp = 1, nzc+1
    if( zp .ge. zc(kp-1) .and. zp .le. zc(kp) ) kCellMarker(j) = kp
   enddo
   
   ENDIF ! xp
    
   end select
  ENDDO

  Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)
   zp = bodyIBM(i,j,3)
   
 
   iCell = iCellMarker(j)
   jCell = jCellMarker(j)
   kCell = kCellMarker(j)
   
   IF( iCell .ne. -255 .and. jCell .ne. -255 .and. kCell .ne. -255 ) THEN
   
   imin = max(iCell-1-INTP_CORONA,1-ngl)
   imax = min(iCell  +INTP_CORONA,nxc+ngl)
   jmin = max(jCell-1-INTP_CORONA,1-ngl)
   jmax = min(jCell  +INTP_CORONA,nyc+ngl)
   kmin = max(kCell-1-INTP_CORONA,0)
   kmax = min(kCell  +INTP_CORONA,nzc+1)
   
   select case ( ndim1 )
   
   case ( 2 )
  
   do ii=imin,imax
   do jj=jmin,jmax
   
   kk = kCell
   
   gx = (xp-xc(L2GI(ii)))/dxc(L2GI(iCell))
   gy = (yp-yc(L2GJ(jj)))/dyc(L2GJ(jCell))
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
 
   velIBM(i,j,1) = velIBM(i,j,1) + dfx*dfy*u(ii,jj,kk)
   velIBM(i,j,2) = velIBM(i,j,2) + dfx*dfy*v(ii,jj,kk)
    
   enddo
   enddo
   
   
   case ( 3 )
   
! print*, 'before delta'
    
   do ii=imin, imax
   do jj=jmin, jmax
   do kk=kmin, kmax

   gx = (xp-xc(L2GI(ii)))/dxc(L2GI(iCell))
   gy = (yp-yc(L2GJ(jj)))/dyc(L2GJ(jCell))
   gz = (zp-zc(kk)      )/dzc(kcell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
   call deltao_fun(dtype(3),gz,dfz)
 
   velIBM(i,j,1) = velIBM(i,j,1) + dfx*dfy*dfz*u(ii,jj,kk)
   velIBM(i,j,2) = velIBM(i,j,2) + dfx*dfy*dfz*v(ii,jj,kk)   
   velIBM(i,j,3) = velIBM(i,j,3) + dfx*dfy*dfz*w(ii,jj,kk)
  
   enddo
   enddo
   enddo
   
   end select 
    
  ENDIF ! iCell
  
  ENDDO ! nIBM
  ENDDO ! nrIBM
  
 tempvel=velIBM(:,:,1)
#if (mpi_activate==1)
 tempvel=0.0 
 CALL MPI_ALLREDUCE(velIBM(:,:,1),tempvel,nIBM*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
 velIBM(:,:,1) = tempvel

 tempvel=velIBM(:,:,2)
#if (mpi_activate==1)
 tempvel=0.0
 CALL MPI_ALLREDUCE(velIBM(:,:,2),tempvel,nIBM*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
 velIBM(:,:,2) = tempvel

 if (ndim1 .gt. 2) then
  tempvel=velIBM(:,:,3)
#if (mpi_activate==1)
  tempvel=0.0
  CALL MPI_ALLREDUCE(velIBM(:,:,3),tempvel,nIBM*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
  velIBM(:,:,3) = tempvel
  tempvel=0.0
 ENDIF !ndim1   
 DEALLOCATE(tempvel)
 DEALLOCATE(iCellMarker,jCellMarker,kCellMarker)
 
END SUBROUTINE
!====================================================================================
SUBROUTINE cal_velIBM_fbc(ndim1, bodyIBM, velIBM, dtype,nrIBM,nIBM,nIBM_r)
 USE global_parameters
 USE flow_parameters
 USE unstructured_surface_arrays 
 USE boundary_arrays
 USE flow_arrays
 USE grid_arrays
#if (mpi_activate==1)
 USE mpi
#endif
 
 implicit NONE
 
  integer i, j, iCell, jCell, kCell, ii, jj, kk, ierr1,jq
  integer ndim1, INTP_CORONA, ip, jp, kp, dtype(3)
  integer nIBM, nrIBM, nIBM_r(nrIBM)
  REAL*8 bodyIBM(nrIBM,nIBM,3), velIBM(nrIBM,nIBM,3)
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tempvel
  REAL*8 xp, yp, zp, gx, gy, gz, dfx, dfy, dfz
  integer imin, imax, jmin, jmax, kmin, kmax
  integer, allocatable, dimension (:) :: iCellMarker, jCellMarker, kCellMarker

 ALLOCATE(iCellMarker(nIBM))
 ALLOCATE(jCellMarker(nIBM))
 ALLOCATE(kCellMarker(nIBM))
 
 ALLOCATE(tempvel(nrIBM,nIBM))
  
 INTP_CORONA = 1
 
 velIBM = 0.0
 tempvel = 0.0



 ! find grid index
Do i = 1,nrIBM

   Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)
   zp = bodyIBM(i,j,3)

   iCellMarker(j) = -255
   jCellMarker(j) = -255
   kCellMarker(j) = -255
   
   select case ( ndim1 )
   
   case ( 2 )
   
   ! check the body marker is inside the current sub-domain
   IF( xp .ge. x(myIS) .and. xp .lt. x(myIe+1) .and. yp .ge. y(myJs) .and. yp .lt. y(myJe+1) ) THEN
   
   !-----+---------+---------
   !     |    x    |             
   !-----+---------+---------
   !    ip-1       ip
   !         iCell = ip
   
   ! ip, jp, kp are local index
   
   do ip = 1, nxc+1
    if( xp .ge. xc(L2GI(ip-1)) .and. xp .le. xc(L2GI(ip)))  iCellMarker(j) = ip
   enddo
   
   do jp = 1, nyc+1
    if( yp .ge. yc(L2GJ(jp-1)) .and. yp .le. yc(L2GJ(jp)) ) jCellMarker(j) = jp
   enddo   
   
   kCellMarker(j) = 1
   
   ENDIF ! xp
   
   case ( 3 )
    
   IF( xp .ge. x(myIS) .and. xp .lt. x(myIe+1) .and. yp .ge. y(myJs) .and. yp .lt. y(myJe+1) .and. zp .ge. z(1) .and. zp .le. z(nz) ) THEN
   
   do ip = 1, nxc+1
    if( xp .ge. xc(L2GI(ip-1)) .and. xp .le. xc(L2GI(ip))) iCellMarker(j) = ip
   enddo
   
   do jp = 1, nyc+1
    if( yp .ge. yc(L2GJ(jp-1)) .and. yp .le. yc(L2GJ(jp)) ) jCellMarker(j) = jp
   enddo
   
   do kp = 1, nzc+1
    if( zp .ge. zc(kp-1) .and. zp .le. zc(kp) ) kCellMarker(j) = kp
   enddo
   
   ENDIF ! xp
    
   end select
  ENDDO

  Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)
   zp = bodyIBM(i,j,3)
   
 
   iCell = iCellMarker(j)
   jCell = jCellMarker(j)
   kCell = kCellMarker(j)
   
   IF( iCell .ne. -255 .and. jCell .ne. -255 .and. kCell .ne. -255 ) THEN
   
   imin = max(iCell-1-INTP_CORONA,1-ngl)
   imax = min(iCell  +INTP_CORONA,nxc+ngl)
   jmin = max(jCell-1-INTP_CORONA,1-ngl)
   jmax = min(jCell  +INTP_CORONA,nyc+ngl)
   kmin = max(kCell-1-INTP_CORONA,0)
   kmax = min(kCell  +INTP_CORONA,nzc+1)
   
   select case ( ndim1 )
   
   case ( 2 )
  
   do ii=imin,imax
   do jj=jmin,jmax
   
   kk = kCell
   
   gx = (xp-xc(L2GI(ii)))/dxc(L2GI(iCell))
   gy = (yp-yc(L2GJ(jj)))/dyc(L2GJ(jCell))
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
 
   velIBM(i,j,1) = velIBM(i,j,1) + dfx*dfy*u(ii,jj,kk)
   velIBM(i,j,2) = velIBM(i,j,2) + dfx*dfy*v(ii,jj,kk)
    
   enddo
   enddo
   
   
   case ( 3 )
   
! print*, 'before delta'
    
   do ii=imin, imax
   do jj=jmin, jmax
   do kk=kmin, kmax

   gx = (xp-xc(L2GI(ii)))/dxc(L2GI(iCell))
   gy = (yp-yc(L2GJ(jj)))/dyc(L2GJ(jCell))
   gz = (zp-zc(kk)      )/dzc(kcell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
   call deltao_fun(dtype(3),gz,dfz)
 
   velIBM(i,j,1) = velIBM(i,j,1) + dfx*dfy*dfz*u(ii,jj,kk)
   velIBM(i,j,2) = velIBM(i,j,2) + dfx*dfy*dfz*v(ii,jj,kk)   
   velIBM(i,j,3) = velIBM(i,j,3) + dfx*dfy*dfz*w(ii,jj,kk)
  
   enddo
   enddo
   enddo
   
   end select 
    
  ENDIF ! iCell
  
  ENDDO ! nIBM
  ENDDO ! nrIBM
  
 tempvel=velIBM(:,:,1)
#if (mpi_activate==1)
 tempvel=0.0 
 CALL MPI_ALLREDUCE(velIBM(:,:,1),tempvel,nIBM*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
 velIBM(:,:,1) = tempvel

 tempvel=velIBM(:,:,2)
#if (mpi_activate==1)
 tempvel=0.0
 CALL MPI_ALLREDUCE(velIBM(:,:,2),tempvel,nIBM*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
 velIBM(:,:,2) = tempvel

 if (ndim1 .gt. 2) then
  tempvel=velIBM(:,:,3)
#if (mpi_activate==1)
  tempvel=0.0
  CALL MPI_ALLREDUCE(velIBM(:,:,3),tempvel,nIBM*nrIBM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
#endif
  velIBM(:,:,3) = tempvel
  tempvel=0.0
 ENDIF !ndim1   
 DEALLOCATE(tempvel)
 DEALLOCATE(iCellMarker,jCellMarker,kCellMarker)
 
END SUBROUTINE
!====================================================================================
!====================================================================================
!====================================================================================

! distribute the force on the body marker point to surrounding cartesian grids
! Fxyz = fIBM*d(x)*d(y)*d(z)*ds/dV
!
! inputs:
! 	ndim1 : dimension, 2 or 3
!   bodyIBM(nIBM,3) : body marker coordinates
!   fIBM(nIBM,3) : force vector on the body marker points
!   dsIBM(nIBM) : element area
!   nIBM : number of body marker points
!   dtype : delta function type
!

SUBROUTINE cal_fxyz_fib(ndim1, bodyIBM, fIBM, dsIBM, dtype,dsecIBMy,dsecIBMz,nrIBM,nsecIBMmax,nsecIBM,nIBM,nIBM_r)
 USE global_parameters
 USE flow_parameters
 USE unstructured_surface_arrays 
 USE boundary_arrays
 USE grid_arrays
 USE multiuse_arrays
#if (mpi_activate==1)
 USE mpi
#endif
 
 implicit NONE
 
 integer i, j, iCell, jCell, kCell, ii, jj, kk, ierr1,jq
 integer ndim1, INTP_CORONA, ip, jp, kp, dtype(3)
 integer nsecIBM,nsecIBMmax
 integer nIBM, nrIBM, nIBM_r(nrIBM)
 REAL*8 bodyIBM(nrIBM,nIBM,3), fIBM(nrIBM,nsecIBMmax,nIBM,3), dsIBM(nrIBM,nIBM)
 REAL*8 xp, yp, zp, gx, gy, gz, dfx, dfy, dfz
 integer imin, imax, jmin, jmax, kmin, kmax
 integer, allocatable, dimension (:,:) :: iCellMarker, jCellMarker, kCellMarker
 integer isecIBM,isec
 REAL*8  dsecIBMy,dsecIBMz,zsecIBM,ysecIBM,CoefsecIBM

 ALLOCATE(iCellMarker(nsecIBM,nIBM))
 ALLOCATE(jCellMarker(nsecIBM,nIBM))
 ALLOCATE(kCellMarker(nsecIBM,nIBM))
  
 INTP_CORONA = 1
 
 ! find grid index
  Do i = 1, nrIBM 
  Do isec = 1, nsecIBM
   zsecIBM=real(isec*2-nsecIBM-1)/2.0*dsecIBMz
   ysecIBM=real(isec*2-nsecIBM-1)/2.0*dsecIBMy

  Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)+ysecIBM
   zp = bodyIBM(i,j,3)+zsecIBM
   
   iCellMarker(isec,j) = -255
   jCellMarker(isec,j) = -255
   kCellMarker(isec,j) = -255
   
   select case ( ndim1 )
   
   case ( 2 )
   
   !-----+---------+---------
   !     |    x    |             
   !-----+---------+---------
   !    ip-1       ip
   !         iCell = ip

   ! ip,jp,kp are global index now 
   
   do ip = 1, nxc_GLBL+1
    if( xp .ge. xc(ip-1) .and. xp .le. xc(ip))  iCellMarker(isec,j) = ip
   enddo
   
   do jp = 1, nyc_GLBL+1
    if( yp .ge. yc(jp-1) .and. yp .le. yc(jp))  jCellMarker(isec,j) = jp
   enddo   
   
   kCellMarker(isec,j) = 1

   
   case ( 3 )
    
   do ip = 1, nxc_GLBL+1
    if( xp .ge. xc(ip-1) .and. xp .le. xc(ip)) iCellMarker(isec,j) = ip
   enddo
   
   do jp = 1, nyc_GLBL+1
    if( yp .ge. yc(jp-1) .and. yp .le. yc(jp)) jCellMarker(isec,j) = jp
   enddo
   
   do kp = 1, nzc+1
    if( zp .ge. zc(kp-1) .and. zp .le. zc(kp) ) kCellMarker(isec,j) = kp
   enddo
    
   end select
  
  ENDDO !nIBM
  ENDDO !nsecIBM
 
 
 ! distribute force
  Do isec = 1, nsecIBM
   ysecIBM=real(isec*2-nsecIBM-1)/2.0*dsecIBMy
   zsecIBM=real(isec*2-nsecIBM-1)/2.0*dsecIBMz

  if(nsecIBM .eq. 1) then
     CoefsecIBM=1.0/real(nsecIBM)
  elseif((isec .eq. 1) .or. (isec .eq. nsecIBM)) then
     CoefsecIBM=0.5*(dsecIBMy+dsecIBMz)
  else
     CoefsecIBM=1.0*(dsecIBMy+dsecIBMz)
  endif

  Do j = 1, nIBM_r(i)   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)+ysecIBM
   zp = bodyIBM(i,j,3)+zsecIBM
   

   iCell = iCellMarker(isec,j)
   jCell = jCellMarker(isec,j)
   kCell = kCellMarker(isec,j)
    
   IF( iCell .ne. -255 .and. jCell .ne. -255 .and. kCell .ne. -255 ) THEN
   
   imin = iCell-1-INTP_CORONA
   imax = iCell  +INTP_CORONA
   jmin = jCell-1-INTP_CORONA
   jmax = jCell  +INTP_CORONA
   kmin = kCell-1-INTP_CORONA
   kmax = kCell  +INTP_CORONA
   
   select case ( ndim1 )
   
   case ( 2 )
  
   do ii=imin,imax
   do jj=jmin,jmax
   
   kk = kCell
   
   ! if it is in my sub-domain
   IF( G2LI(ii) .ge. 1 .and. G2LI(ii) .le. nxc .and. G2LJ(jj) .ge. 1 .and. G2LJ(jj) .le. nyc ) THEN
   
   gx = (xp-xc(ii))/dxc(iCell)
   gy = (yp-yc(jj))/dyc(jCell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
 
   nlu(G2LI(ii),G2LJ(jj),kk) = nlu(G2LI(ii),G2LJ(jj),kk) + CoefsecIBM*dt*dfx*dfy*dsIBM(i,j)*fIBM(i,isec,j,1)/(dx(ii)*dy(jj))
   nlv(G2LI(ii),G2LJ(jj),kk) = nlv(G2LI(ii),G2LJ(jj),kk) + CoefsecIBM*dt*dfx*dfy*dsIBM(i,j)*fIBM(i,isec,j,2)/(dx(ii)*dy(jj))
   
   ENDIF
    
   enddo
   enddo
   
   
   case ( 3 )
   
   do ii=imin, imax
   do jj=jmin, jmax
   do kk=kmin, kmax

   ! if it is in my sub-domain
   IF( G2LI(ii) .ge. 1 .and. G2LI(ii) .le. nxc .and. G2LJ(jj) .ge. 1 .and. G2LJ(jj) .le. nyc .and. kk .ge. 1 .and. kk .le. nzc ) THEN
   
   gx = (xp-xc(ii))/dxc(iCell)
   gy = (yp-yc(jj))/dyc(jCell)
   gz = (zp-zc(kk))/dzc(kCell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
   call deltao_fun(dtype(3),gz,dfz)

   nlu(G2LI(ii),G2LJ(jj),kk) = nlu(G2LI(ii),G2LJ(jj),kk) + CoefsecIBM*dt*dfx*dfy*dfz*dsIBM(i,j)*fIBM(i,isec,j,1)/(dx(ii)*dy(jj)*dz(kk))
   nlv(G2LI(ii),G2LJ(jj),kk) = nlv(G2LI(ii),G2LJ(jj),kk) + CoefsecIBM*dt*dfx*dfy*dfz*dsIBM(i,j)*fIBM(i,isec,j,2)/(dx(ii)*dy(jj)*dz(kk))
   nlw(G2LI(ii),G2LJ(jj),kk) = nlw(G2LI(ii),G2LJ(jj),kk) + CoefsecIBM*dt*dfx*dfy*dfz*dsIBM(i,j)*fIBM(i,isec,j,3)/(dx(ii)*dy(jj)*dz(kk))
   ENDIF
  
   enddo
   enddo
   enddo
   
   end select 
    
  ENDIF ! iCell

  ENDDO ! nIBM
  ENDDO ! nsecIBM

  ENDDO !nrIBM
  
 DEALLOCATE(iCellMarker,jCellMarker,kCellMarker)
 
END SUBROUTINE
!====================================================================================
SUBROUTINE cal_fxyz_fsh(ndim1, bodyIBM, fIBM, dsIBM, dtype,nrIBM,nIBM,nqIBM,nIBM_r,nIBM_q)
 USE global_parameters
 USE flow_parameters
 USE unstructured_surface_arrays 
 USE boundary_arrays
 USE grid_arrays
 USE multiuse_arrays
#if (mpi_activate==1)
 USE mpi
#endif
 
 implicit NONE
 
 integer i, j, iCell, jCell, kCell, ii, jj, kk, ierr1,jq
 integer ndim1, INTP_CORONA, ip, jp, kp, dtype(3)
 integer nqIBM,nIBM, nrIBM, nIBM_r(nrIBM),nIBM_q(nrIBM)
 REAL*8 bodyIBM(nrIBM,nqIBM,nIBM,3), fIBM(nrIBM,nqIBM,nIBM,3), dsIBM(nrIBM,nqIBM,nIBM)
 REAL*8 xp, yp, zp, gx, gy, gz, dfx, dfy, dfz
 integer imin, imax, jmin, jmax, kmin, kmax
 integer, allocatable, dimension (:,:) :: iCellMarker, jCellMarker, kCellMarker

 ALLOCATE(iCellMarker(nqIBM,nIBM))
 ALLOCATE(jCellMarker(nqIBM,nIBM))
 ALLOCATE(kCellMarker(nqIBM,nIBM))
  
 INTP_CORONA = 1
 
 ! find grid index
  Do i = 1, nrIBM 


  Do jq= 1, nIBM_q(i)
  Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,jq,j,1)
   yp = bodyIBM(i,jq,j,2)
   zp = bodyIBM(i,jq,j,3)
   
   iCellMarker(jq,j) = -255
   jCellMarker(jq,j) = -255
   kCellMarker(jq,j) = -255
   
   select case ( ndim1 )
   
   case ( 2 )
   
   !-----+---------+---------
   !     |    x    |             
   !-----+---------+---------
   !    ip-1       ip
   !         iCell = ip

   ! ip,jp,kp are global index now 
   
   do ip = 1, nxc_GLBL+1
    if( xp .ge. xc(ip-1) .and. xp .le. xc(ip))  iCellMarker(jq,j) = ip
   enddo
   
   do jp = 1, nyc_GLBL+1
    if( yp .ge. yc(jp-1) .and. yp .le. yc(jp))  jCellMarker(jq,j) = jp
   enddo   
   
   kCellMarker(jq,j) = 1

   
   case ( 3 )
    
   do ip = 1, nxc_GLBL+1
    if( xp .ge. xc(ip-1) .and. xp .le. xc(ip)) iCellMarker(jq,j) = ip
   enddo
   
   do jp = 1, nyc_GLBL+1
    if( yp .ge. yc(jp-1) .and. yp .le. yc(jp)) jCellMarker(jq,j) = jp
   enddo
   
   do kp = 1, nzc+1
    if( zp .ge. zc(kp-1) .and. zp .le. zc(kp) ) kCellMarker(jq,j) = kp
   enddo
    
   end select
  
  ENDDO !nIBM
  ENDDO !nqIBM
 
 
 ! distribute force

  Do jq= 1, nIBM_q(i)
  Do j = 1, nIBM_r(i)   
   xp = bodyIBM(i,jq,j,1)
   yp = bodyIBM(i,jq,j,2)
   zp = bodyIBM(i,jq,j,3)
   

   iCell = iCellMarker(jq,j)
   jCell = jCellMarker(jq,j)
   kCell = kCellMarker(jq,j)
    
   IF( iCell .ne. -255 .and. jCell .ne. -255 .and. kCell .ne. -255 ) THEN
   
   imin = iCell-1-INTP_CORONA
   imax = iCell  +INTP_CORONA
   jmin = jCell-1-INTP_CORONA
   jmax = jCell  +INTP_CORONA
   kmin = kCell-1-INTP_CORONA
   kmax = kCell  +INTP_CORONA
   
   select case ( ndim1 )
   
   case ( 2 )
  
   do ii=imin,imax
   do jj=jmin,jmax
   
   kk = kCell
   
   ! if it is in my sub-domain
   IF( G2LI(ii) .ge. 1 .and. G2LI(ii) .le. nxc .and. G2LJ(jj) .ge. 1 .and. G2LJ(jj) .le. nyc ) THEN
   
   gx = (xp-xc(ii))/dxc(iCell)
   gy = (yp-yc(jj))/dyc(jCell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
 
   nlu(G2LI(ii),G2LJ(jj),kk) = nlu(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dsIBM(i,jq,j)*fIBM(i,jq,j,1)/(dx(ii)*dy(jj))
   nlv(G2LI(ii),G2LJ(jj),kk) = nlv(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dsIBM(i,jq,j)*fIBM(i,jq,j,2)/(dx(ii)*dy(jj))
   
   ENDIF
    
   enddo
   enddo
   
   
   case ( 3 )
   
   do ii=imin, imax
   do jj=jmin, jmax
   do kk=kmin, kmax

   ! if it is in my sub-domain
   IF( G2LI(ii) .ge. 1 .and. G2LI(ii) .le. nxc .and. G2LJ(jj) .ge. 1 .and. G2LJ(jj) .le. nyc .and. kk .ge. 1 .and. kk .le. nzc ) THEN
   
   gx = (xp-xc(ii))/dxc(iCell)
   gy = (yp-yc(jj))/dyc(jCell)
   gz = (zp-zc(kk))/dzc(kCell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
   call deltao_fun(dtype(3),gz,dfz)

   nlu(G2LI(ii),G2LJ(jj),kk) = nlu(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dfz*dsIBM(i,jq,j)*fIBM(i,jq,j,1)/(dx(ii)*dy(jj)*dz(kk))
   nlv(G2LI(ii),G2LJ(jj),kk) = nlv(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dfz*dsIBM(i,jq,j)*fIBM(i,jq,j,2)/(dx(ii)*dy(jj)*dz(kk))
   nlw(G2LI(ii),G2LJ(jj),kk) = nlw(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dfz*dsIBM(i,jq,j)*fIBM(i,jq,j,3)/(dx(ii)*dy(jj)*dz(kk))
   ENDIF
  
   enddo
   enddo
   enddo
   
   end select 
    
  ENDIF ! iCell

  ENDDO ! nqIBM  
  ENDDO ! nIBM

  ENDDO !nrIBM
  
 DEALLOCATE(iCellMarker,jCellMarker,kCellMarker)
 
END SUBROUTINE
!====================================================================================
SUBROUTINE cal_fxyz_esh(ndim1, bodyIBM, fIBM, dsIBM, dtype,nrIBM,nIBM,nIBM_r)
 USE global_parameters
 USE flow_parameters
 USE unstructured_surface_arrays 
 USE boundary_arrays
 USE grid_arrays
 USE multiuse_arrays
#if (mpi_activate==1)
 USE mpi
#endif
 
 implicit NONE
 
 integer i, j, iCell, jCell, kCell, ii, jj, kk, ierr1,jq
 integer ndim1, INTP_CORONA, ip, jp, kp, dtype(3)
 integer nIBM, nrIBM, nIBM_r(nrIBM)
 REAL*8 bodyIBM(nrIBM,nIBM,3), fIBM(nrIBM,nIBM,3), dsIBM(nrIBM,nIBM)
 REAL*8 xp, yp, zp, gx, gy, gz, dfx, dfy, dfz
 integer imin, imax, jmin, jmax, kmin, kmax
 integer, allocatable, dimension (:) :: iCellMarker, jCellMarker, kCellMarker

 ALLOCATE(iCellMarker(nIBM))
 ALLOCATE(jCellMarker(nIBM))
 ALLOCATE(kCellMarker(nIBM))
  
 INTP_CORONA = 1
 
 ! find grid index
  Do i = 1, nrIBM 


   Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)
   zp = bodyIBM(i,j,3)
   
   iCellMarker(j) = -255
   jCellMarker(j) = -255
   kCellMarker(j) = -255
   
   select case ( ndim1 )
   
   case ( 2 )
   
   !-----+---------+---------
   !     |    x    |             
   !-----+---------+---------
   !    ip-1       ip
   !         iCell = ip

   ! ip,jp,kp are global index now 
   
   do ip = 1, nxc_GLBL+1
    if( xp .ge. xc(ip-1) .and. xp .le. xc(ip))  iCellMarker(j) = ip
   enddo
   
   do jp = 1, nyc_GLBL+1
    if( yp .ge. yc(jp-1) .and. yp .le. yc(jp))  jCellMarker(j) = jp
   enddo   
   
   kCellMarker(j) = 1

   
   case ( 3 )
    
   do ip = 1, nxc_GLBL+1
    if( xp .ge. xc(ip-1) .and. xp .le. xc(ip)) iCellMarker(j) = ip
   enddo
   
   do jp = 1, nyc_GLBL+1
    if( yp .ge. yc(jp-1) .and. yp .le. yc(jp)) jCellMarker(j) = jp
   enddo
   
   do kp = 1, nzc+1
    if( zp .ge. zc(kp-1) .and. zp .le. zc(kp) ) kCellMarker(j) = kp
   enddo
    
   end select
  
  ENDDO !nIBM

 
 
 ! distribute force

  Do j = 1, nIBM_r(i)   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)
   zp = bodyIBM(i,j,3)
   

   iCell = iCellMarker(j)
   jCell = jCellMarker(j)
   kCell = kCellMarker(j)
    
   IF( iCell .ne. -255 .and. jCell .ne. -255 .and. kCell .ne. -255 ) THEN
   
   imin = iCell-1-INTP_CORONA
   imax = iCell  +INTP_CORONA
   jmin = jCell-1-INTP_CORONA
   jmax = jCell  +INTP_CORONA
   kmin = kCell-1-INTP_CORONA
   kmax = kCell  +INTP_CORONA
   
   select case ( ndim1 )
   
   case ( 2 )
  
   do ii=imin,imax
   do jj=jmin,jmax
   
   kk = kCell
   
   ! if it is in my sub-domain
   IF( G2LI(ii) .ge. 1 .and. G2LI(ii) .le. nxc .and. G2LJ(jj) .ge. 1 .and. G2LJ(jj) .le. nyc ) THEN
   
   gx = (xp-xc(ii))/dxc(iCell)
   gy = (yp-yc(jj))/dyc(jCell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
 
   nlu(G2LI(ii),G2LJ(jj),kk) = nlu(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dsIBM(i,j)*fIBM(i,j,1)/(dx(ii)*dy(jj))
   nlv(G2LI(ii),G2LJ(jj),kk) = nlv(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dsIBM(i,j)*fIBM(i,j,2)/(dx(ii)*dy(jj))
   
   ENDIF
    
   enddo
   enddo
   
   
   case ( 3 )
   
   do ii=imin, imax
   do jj=jmin, jmax
   do kk=kmin, kmax

   ! if it is in my sub-domain
   IF( G2LI(ii) .ge. 1 .and. G2LI(ii) .le. nxc .and. G2LJ(jj) .ge. 1 .and. G2LJ(jj) .le. nyc .and. kk .ge. 1 .and. kk .le. nzc ) THEN
   
   gx = (xp-xc(ii))/dxc(iCell)
   gy = (yp-yc(jj))/dyc(jCell)
   gz = (zp-zc(kk))/dzc(kCell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
   call deltao_fun(dtype(3),gz,dfz)

   nlu(G2LI(ii),G2LJ(jj),kk) = nlu(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dfz*dsIBM(i,j)*fIBM(i,j,1)/(dx(ii)*dy(jj)*dz(kk))
   nlv(G2LI(ii),G2LJ(jj),kk) = nlv(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dfz*dsIBM(i,j)*fIBM(i,j,2)/(dx(ii)*dy(jj)*dz(kk))
   nlw(G2LI(ii),G2LJ(jj),kk) = nlw(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dfz*dsIBM(i,j)*fIBM(i,j,3)/(dx(ii)*dy(jj)*dz(kk))
   ENDIF
  
   enddo
   enddo
   enddo
   
   end select 
    
  ENDIF ! iCell

  ENDDO ! nIBM

  ENDDO !nrIBM
  
 DEALLOCATE(iCellMarker,jCellMarker,kCellMarker)
 
END SUBROUTINE
!====================================================================================
SUBROUTINE cal_fxyz_fbc(ndim1, bodyIBM, fIBM, dsIBM, dtype,nrIBM,nIBM,nIBM_r)
 USE global_parameters
 USE flow_parameters
 USE unstructured_surface_arrays 
 USE boundary_arrays
 USE grid_arrays
 USE multiuse_arrays
#if (mpi_activate==1)
 USE mpi
#endif
 
 implicit NONE
 
 integer i, j, iCell, jCell, kCell, ii, jj, kk, ierr1,jq
 integer ndim1, INTP_CORONA, ip, jp, kp, dtype(3)
 integer nIBM, nrIBM, nIBM_r(nrIBM)
 REAL*8 bodyIBM(nrIBM,nIBM,3), fIBM(nrIBM,nIBM,3), dsIBM(nrIBM,nIBM)
 REAL*8 xp, yp, zp, gx, gy, gz, dfx, dfy, dfz
 integer imin, imax, jmin, jmax, kmin, kmax
 integer, allocatable, dimension (:) :: iCellMarker, jCellMarker, kCellMarker

 ALLOCATE(iCellMarker(nIBM))
 ALLOCATE(jCellMarker(nIBM))
 ALLOCATE(kCellMarker(nIBM))
  
 INTP_CORONA = 1
 
 ! find grid index
  Do i = 1, nrIBM 


   Do j = 1, nIBM_r(i)
   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)
   zp = bodyIBM(i,j,3)
   
   iCellMarker(j) = -255
   jCellMarker(j) = -255
   kCellMarker(j) = -255
   
   select case ( ndim1 )
   
   case ( 2 )
   
   !-----+---------+---------
   !     |    x    |             
   !-----+---------+---------
   !    ip-1       ip
   !         iCell = ip

   ! ip,jp,kp are global index now 
   
   do ip = 1, nxc_GLBL+1
    if( xp .ge. xc(ip-1) .and. xp .le. xc(ip))  iCellMarker(j) = ip
   enddo
   
   do jp = 1, nyc_GLBL+1
    if( yp .ge. yc(jp-1) .and. yp .le. yc(jp))  jCellMarker(j) = jp
   enddo   
   
   kCellMarker(j) = 1

   
   case ( 3 )
    
   do ip = 1, nxc_GLBL+1
    if( xp .ge. xc(ip-1) .and. xp .le. xc(ip)) iCellMarker(j) = ip
   enddo
   
   do jp = 1, nyc_GLBL+1
    if( yp .ge. yc(jp-1) .and. yp .le. yc(jp)) jCellMarker(j) = jp
   enddo
   
   do kp = 1, nzc+1
    if( zp .ge. zc(kp-1) .and. zp .le. zc(kp) ) kCellMarker(j) = kp
   enddo
    
   end select
  
  ENDDO !nIBM

 
 
 ! distribute force

  Do j = 1, nIBM_r(i)   
   xp = bodyIBM(i,j,1)
   yp = bodyIBM(i,j,2)
   zp = bodyIBM(i,j,3)
   

   iCell = iCellMarker(j)
   jCell = jCellMarker(j)
   kCell = kCellMarker(j)
    
   IF( iCell .ne. -255 .and. jCell .ne. -255 .and. kCell .ne. -255 ) THEN
   
   imin = iCell-1-INTP_CORONA
   imax = iCell  +INTP_CORONA
   jmin = jCell-1-INTP_CORONA
   jmax = jCell  +INTP_CORONA
   kmin = kCell-1-INTP_CORONA
   kmax = kCell  +INTP_CORONA
   
   select case ( ndim1 )
   
   case ( 2 )
  
   do ii=imin,imax
   do jj=jmin,jmax
   
   kk = kCell
   
   ! if it is in my sub-domain
   IF( G2LI(ii) .ge. 1 .and. G2LI(ii) .le. nxc .and. G2LJ(jj) .ge. 1 .and. G2LJ(jj) .le. nyc ) THEN
   
   gx = (xp-xc(ii))/dxc(iCell)
   gy = (yp-yc(jj))/dyc(jCell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
 
   nlu(G2LI(ii),G2LJ(jj),kk) = nlu(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dsIBM(i,j)*fIBM(i,j,1)/(dx(ii)*dy(jj))
   nlv(G2LI(ii),G2LJ(jj),kk) = nlv(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dsIBM(i,j)*fIBM(i,j,2)/(dx(ii)*dy(jj))
   
   ENDIF
    
   enddo
   enddo
   
   
   case ( 3 )
   
   do ii=imin, imax
   do jj=jmin, jmax
   do kk=kmin, kmax

   ! if it is in my sub-domain
   IF( G2LI(ii) .ge. 1 .and. G2LI(ii) .le. nxc .and. G2LJ(jj) .ge. 1 .and. G2LJ(jj) .le. nyc .and. kk .ge. 1 .and. kk .le. nzc ) THEN
   
   gx = (xp-xc(ii))/dxc(iCell)
   gy = (yp-yc(jj))/dyc(jCell)
   gz = (zp-zc(kk))/dzc(kCell)
   
   call deltao_fun(dtype(1),gx,dfx)
   call deltao_fun(dtype(2),gy,dfy)
   call deltao_fun(dtype(3),gz,dfz)

   nlu(G2LI(ii),G2LJ(jj),kk) = nlu(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dfz*dsIBM(i,j)*fIBM(i,j,1)/(dx(ii)*dy(jj)*dz(kk))
   nlv(G2LI(ii),G2LJ(jj),kk) = nlv(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dfz*dsIBM(i,j)*fIBM(i,j,2)/(dx(ii)*dy(jj)*dz(kk))
   nlw(G2LI(ii),G2LJ(jj),kk) = nlw(G2LI(ii),G2LJ(jj),kk) + dt*dfx*dfy*dfz*dsIBM(i,j)*fIBM(i,j,3)/(dx(ii)*dy(jj)*dz(kk))
   ENDIF
  
   enddo
   enddo
   enddo
   
   end select 
    
  ENDIF ! iCell

  ENDDO ! nIBM

  ENDDO !nrIBM
  
 DEALLOCATE(iCellMarker,jCellMarker,kCellMarker)
 
END SUBROUTINE
!====================================================================================
!====================================================================================
!====================================================================================

subroutine init_membrane_solver(&
 idim1,&
 nReadin,&
 theboss)
 
 USE dummy_module
 USE global_parameters
 USE flow_parameters
 USE grid_arrays
#if (mpi_activate==1)
 USE mpi
#endif

 Implicit none
 integer idim1
 integer nReadin
 logical theboss
 integer i,j,isec,jq,iErrin,ierr1
 real*8 min_grid_x,min_grid_y,min_grid_z
 WRITE(*,*) 'Nr_IBM', Nr_IBM
      ALLOCATE(nIBM_r(Nr_IBM),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
        'Allocate_memory: Memory Allocation Error for nIBM_r'
        STOP
      ENDIF ! iErrin
      ALLOCATE(nIBM_q(Nr_IBM),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for nIBM_q'
        STOP
      ENDIF ! iErrin

!%===============================================================
      ALLOCATE(nIBM_r_fib(Nr_IBM_fib),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
        'Allocate_memory: Memory Allocation Error for nIBM_r_fib'
        STOP
      ENDIF ! iErrin
!%===============================================================
      ALLOCATE(nIBM_r_fsh(Nr_IBM_fsh),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
        'Allocate_memory: Memory Allocation Error for nIBM_r_fsh'
        STOP
      ENDIF ! iErrin
      ALLOCATE(nIBM_q_fsh(Nr_IBM_fsh),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for nIBM_q_fsh'
        STOP
      ENDIF ! iErrin
!%===============================================================
      ALLOCATE(nIBM_r_esh(Nr_IBM_esh),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
        'Allocate_memory: Memory Allocation Error for nIBM_r_esh'
        STOP
      ENDIF ! iErrin!
!%===============================================================
      ALLOCATE(nIBM_r_fbc(Nr_IBM_esh),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
        'Allocate_memory: Memory Allocation Error for nIBM_r_fbc'
        STOP
      ENDIF ! iErrin

!%===============================================================
 

      ALLOCATE(coord_fib(nr_IBM_fib,ns_IBM_fib,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for coord_fib'
        STOP
      ENDIF ! iErrin

      ALLOCATE(vel_fib(nr_IBM_fib,Nsec_IBMmax,ns_IBM_fib,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for vel_fib'
        STOP
      ENDIF ! iErrin

      ALLOCATE(force_fib(nr_IBM_fib,Nsec_IBMmax,ns_IBM_fib,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for force_fib'
        STOP
      ENDIF ! iErrin

      ALLOCATE(ds_fib(nr_IBM_fib,ns_IBM_fib),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for ds_fib'
        STOP
      ENDIF ! iErrin


      ALLOCATE(coord_fib_prev(nr_IBM_fib,ns_IBM_fib,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for coord_fib_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(vel_fib_prev(nr_IBM_fib,Nsec_IBMmax,ns_IBM_fib,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for vel_fib_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(force_fib_prev(nr_IBM_fib,Nsec_IBMmax,ns_IBM_fib,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for force_fib_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(ds_fib_prev(nr_IBM_fib,ns_IBM_fib),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for ds_fib_prev'
        STOP
      ENDIF ! iErrin


!%===============================================================
      ALLOCATE(coord_fsh(nr_IBM_fsh,nq_IBM_fsh,ns_IBM_fsh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for coord_fsh'
        STOP
      ENDIF ! iErrin

      ALLOCATE(vel_fsh(nr_IBM_fsh,nq_IBM_fsh,ns_IBM_fsh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for vel_fsh'
        STOP
      ENDIF ! iErrin

      ALLOCATE(force_fsh(nr_IBM_fsh,nq_IBM_fsh,ns_IBM_fsh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for force_fsh'
        STOP
      ENDIF ! iErrin

      ALLOCATE(ds_fsh(nr_IBM_fsh,nq_IBM_fsh,ns_IBM_fsh),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for ds_fsh'
        STOP
      ENDIF ! iErrin


      ALLOCATE(coord_fsh_prev(nr_IBM_fsh,nq_IBM_fsh,ns_IBM_fsh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for coord_fsh_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(vel_fsh_prev(nr_IBM_fsh,nq_IBM_fsh,ns_IBM_fsh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for vel_fsh_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(force_fsh_prev(nr_IBM_fsh,nq_IBM_fsh,ns_IBM_fsh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for force_fsh_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(ds_fsh_prev(nr_IBM_fsh,nq_IBM_fsh,ns_IBM_fsh),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for ds_fsh_prev'
        STOP
      ENDIF ! iErrin


!%===============================================================

      ALLOCATE(coord_esh(nr_IBM_esh,ns_IBM_esh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for coord_esh'
        STOP
      ENDIF ! iErrin

      ALLOCATE(vel_esh(nr_IBM_esh,ns_IBM_esh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for vel_esh'
        STOP
      ENDIF ! iErrin

      ALLOCATE(force_esh(nr_IBM_esh,ns_IBM_esh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for force_esh'
        STOP
      ENDIF ! iErrin

      ALLOCATE(ds_esh(nr_IBM_esh,ns_IBM_esh),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for ds_esh'
        STOP
      ENDIF ! iErrin


      ALLOCATE(coord_esh_prev(nr_IBM_esh,ns_IBM_esh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for coord_esh_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(vel_esh_prev(nr_IBM_esh,ns_IBM_esh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for vel_esh_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(force_esh_prev(nr_IBM_esh,ns_IBM_esh,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for force_esh_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(ds_esh_prev(nr_IBM_esh,ns_IBM_esh),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for ds_esh_prev'
        STOP
      ENDIF ! iErrin

!%===============================================================

      ALLOCATE(coord_fbc(nr_IBM_fbc,ns_IBM_fbc,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for coord_fbc'
        STOP
      ENDIF ! iErrin

      ALLOCATE(vel_fbc(nr_IBM_fbc,ns_IBM_fbc,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for vel_fbc'
        STOP
      ENDIF ! iErrin

      ALLOCATE(force_fbc(nr_IBM_fbc,ns_IBM_fbc,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for force_fbc'
        STOP
      ENDIF ! iErrin

      ALLOCATE(ds_fbc(nr_IBM_fbc,ns_IBM_fbc),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for ds_fbc'
        STOP
      ENDIF ! iErrin


      ALLOCATE(coord_fbc_prev(nr_IBM_fbc,ns_IBM_fbc,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for coord_fbc_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(vel_fbc_prev(nr_IBM_fbc,ns_IBM_fbc,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for vel_fbc_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(force_fbc_prev(nr_IBM_fbc,ns_IBM_fbc,3),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for force_fbc_prev'
        STOP
      ENDIF ! iErrin

      ALLOCATE(ds_fbc_prev(nr_IBM_fbc,ns_IBM_fbc),STAT=iErrin)
      IF ( iErrin/= 0 ) THEN
        WRITE(*,*) &
       'Allocate_memory: Memory Allocation Error for ds_fbc_prev'
        STOP
      ENDIF ! iErrin


!%===============================================================




!   Loop on cells in X direction
!   ----------------------------
    min_grid_x=1000.0
    DO i=0,nxc_GLBL+1
      min_grid_x=min(min_grid_x,dx(i) )  
    ENDDO
   
!   Loop on cells in Y direction
!   ----------------------------
    min_grid_y=1000.0
    DO i=0,nyc_GLBL+1
      min_grid_y=min(min_grid_y,dy(i) )  
    ENDDO
    if(idim1 .gt. 2) then
!   Loop on cells in Z direction
!   ----------------------------
    min_grid_z=1000.0
    DO i=0,nzc+1
      min_grid_z=min(min_grid_z,dz(i) )  
    ENDDO
    else
        min_grid_z=1.0
    endif

 if(nrIBM_fib>=1) ds_fib=0.0
 if(nrIBM_fsh>=1) ds_fsh=0.0
 if(nrIBM_esh>=1) ds_esh=0.0
 if(nrIBM_fbc>=1) ds_fbc=0.0
!   write(*,*) nrIBM_fsh,nqIBM_fsh,nIBM_fsh
  
 call initialize_ibm(nIBM_q,nIBM_q_fsh & 
                    ,nIBM_r,nIBM_r_fib,nIBM_r_fsh &
                    ,nIBM_r_esh, nIBM_r_fbc, dtypeDelta  & 
                    ,ds_fib(1:nrIBM_fib,1:nIBM_fib) &
                    ,ds_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh) &
                    ,ds_esh(1:nrIBM_esh,1:nIBM_esh) &
                    ,ds_fbc(1:nrIBM_fbc,1:nIBM_fbc) &
                    ,coord_fib(1:nrIBM_fib,1:nIBM_fib,1) &
                    ,coord_fib(1:nrIBM_fib,1:nIBM_fib,2) &
                    ,coord_fib(1:nrIBM_fib,1:nIBM_fib,3) &
                    ,coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1) &
                    ,coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2) &
                    ,coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3) &
                    ,coord_esh(1:nrIBM_esh,1:nIBM_esh,1) &
                    ,coord_esh(1:nrIBM_esh,1:nIBM_esh,2) &
                    ,coord_esh(1:nrIBM_esh,1:nIBM_esh,3) &
                    ,coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,1) &
                    ,coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,2) &
                    ,coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,3) &
                    ,min_grid_x,min_grid_y,min_grid_z &
                    ,nx_GLBL,ny_GLBL,nz &
                    ,x(1:nx_GLBL),y(1:ny_GLBL),z(1:nz)  &
                    ,nsecIBM,dsecIBMy,dsecIBMz   &
                    ,idim1,nReadin &
                    ,theboss)
!                    ,nrIBM,nsecIBMmax &
!                    ,nrIBM_fib,nrIBM_fsh,nrIBM_esh &
!                    ,nqIBM_fsh &
!                    ,nIBM_fib,nIBM_fsh,nIBM_esh &

!   write(*,*) nrIBM_fsh,nqIBM_fsh,nIBM_fsh



 do i=1,nrIBM_fib
    do isec=1,nsecIBM
    do j=1,nIBM_fib
      force_fib(i,isec,j,1:3)=0.d0
      force_fib_prev(i,isec,j,1:3)=0.d0
    enddo
    enddo      
 enddo
 do i=1,nrIBM_fsh
    do jq=1,nqIBM_fsh
    do j=1,nIBM_fsh
      force_fsh(i,jq,j,1:3)=0.d0
      force_fsh_prev(i,jq,j,1:3)=0.d0
    enddo
    enddo      
 enddo
 do i=1,nrIBM_esh
    do j=1,nIBM_esh
      force_esh(i,j,1:3)=0.d0
      force_esh_prev(i,j,1:3)=0.d0
    enddo
 enddo
 do i=1,nrIBM_fbc
    do j=1,nIBM_fbc
      force_fbc(i,j,1:3)=0.d0
      force_fbc_prev(i,j,1:3)=0.d0
    enddo
 enddo
#if (mpi_activate==1)
 CALL MPI_BARRIER(MPI_COMM_WORLD, ierr1)
#endif
return
end subroutine

subroutine solve_membranetest(iter,nmonitor,theboss,dt,rtime,ntsave,irestart,idim1,nrestart,nReadtmp)
 USE dummy_module
#if (mpi_activate==1)
 USE mpi
#endif
 implicit none
 integer irestart,nrestart,ierr1,nReadtmp
 integer iter, i, kk, ntsave,idim1,j
 character*80 fname
 real*8 dt, uavg, vavg, xc, yc, rtime
 logical monitorON,theboss,temprestart
 integer nmonitor
 real*8 xp,yp,zp,time2tmp 


  real*8 testtime,temp_ibm2
  integer itesttime
   testtime=0.d0
      if (nReadtmp .eq. 1) then
!       temprestart=.true.
!       call write_restart_DISRIBM(time2tmp,temprestart,1)
!      testtime=time2tmp
       call adjustsub(testtime)
      endif
DO itesttime = 1,100000000        ! move solution from n --> n+1
    testtime = testtime + dt
    rtime=testtime
    iter=itesttime
  if(mod(iter,nmonitor)==0 .and. theboss) then
   monitorON=.true. 
  else
   monitorON=.false.
  endif
    IF(monitorON) write(*,*) '============================================='
     IF(monitorON) write(*,*) 'time is =', testtime,'at',itesttime
    IF(monitorON) write(*,*) ' INTERPOLATE VELOCITY'
!!!-----------------------------------------------
!    Membrane solver here
 IF(monitorON) write(*,*) ' SOLVING EQ'
 IF(monitorON) write(*,*) ' SOLVING EQ'
 call tick(rtime,dt,iter,monitorON,ntsave, &
      vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), &
      vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), &
      vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3),  &
      vel_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1), &
      vel_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2), &
      vel_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3),  &
      vel_esh(1:nrIBM_esh,1:nIBM_esh,1), &
      vel_esh(1:nrIBM_esh,1:nIBM_esh,2), &
      vel_esh(1:nrIBM_esh,1:nIBM_esh,3),  &
      vel_fbc(1:nrIBM_fbc,1:nIBM_fbc,1), &
      vel_fbc(1:nrIBM_fbc,1:nIBM_fbc,2), &
      vel_fbc(1:nrIBM_fbc,1:nIBM_fbc,3),  &
      force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), &
      force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), &
      force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3), &
      force_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1), &
      force_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2), &
      force_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3), &
      force_esh(1:nrIBM_esh,1:nIBM_esh,1), &
      force_esh(1:nrIBM_esh,1:nIBM_esh,2), &
      force_esh(1:nrIBM_esh,1:nIBM_esh,3), &
      force_fbc(1:nrIBM_fbc,1:nIBM_fbc,1), &
      force_fbc(1:nrIBM_fbc,1:nIBM_fbc,2), &
      force_fbc(1:nrIBM_fbc,1:nIBM_fbc,3), &
      coord_fib(1:nrIBM_fib,1:nIBM_fib,1),   &
      coord_fib(1:nrIBM_fib,1:nIBM_fib,2),   &
      coord_fib(1:nrIBM_fib,1:nIBM_fib,3),   &
      coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1),   &
      coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2),   &
      coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3),   &
      coord_esh(1:nrIBM_esh,1:nIBM_esh,1),   &
      coord_esh(1:nrIBM_esh,1:nIBM_esh,2),   &
      coord_esh(1:nrIBM_esh,1:nIBM_esh,3),   &
      coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,1),   &
      coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,2),   &
      coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,3),   &
      theboss)
!!!-----------------------------------------------

    irestart=0
    IF ( MOD(itesttime,nrestart) == 0 ) THEN
       irestart=1
    ENDIF   
 IF ( irestart == 1 .and. theboss) THEN
     write(*,*) 'Writing restart file by the boss '
     temp_ibm2=0.0
    call write_restart_DISRIBM(temp_ibm2,1,0)
 ENDIF
#if (mpi_activate==1)
 CALL mpi_barrier(MPI_COMM_WORLD, ierr1)
#endif
enddo
pause

end subroutine

subroutine solve_membrane(iter,monitorON,theboss,dt,rtime,ntsave,irestart,idim1)
 USE dummy_module
#if (mpi_activate==1)
 USE mpi
#endif
 implicit none
 integer irestart, ierr1
 logical monitorON ,theboss
 integer iter, i, kk, ntsave,idim1,j
 character*80 fname
 real*8 dt, uavg, vavg, xc, yc, rtime,temp_ibm2

 real*8 xp,yp,zp 
 IF(monitorON) write(*,*) ' INTERPOLATE VELOCITY'
 if(nrIBM_fib>=1) then
  call cal_velIBM_fib( &
   idim1,&
   coord_fib,&
   vel_fib,&
   dtypeDelta(1:3),&
   dsecIBMy,&
   dsecIBMz,&
   nrIBM_fib,&
   nsecIBMmax,&
   nsecIBM,&
   nIBM_fib,&
   nIBM_r_fib)
  end if
 if(nrIBM_fsh>=1) call cal_velIBM_fsh(idim1,coord_fsh,vel_fsh,dtypeDelta(1:3),nrIBM_fsh,nIBM_fsh,nqIBM_fsh,nIBM_r_fsh,nIBM_q_fsh)
 if(nrIBM_esh>=1) call cal_velIBM_esh(idim1,coord_esh,vel_esh,dtypeDelta(1:3),nrIBM_esh,nIBM_esh,nIBM_r_esh)
 if(nrIBM_fbc>=1) call cal_velIBM_fbc(idim1,coord_fbc,vel_fbc,dtypeDelta(1:3),nrIBM_fbc,nIBM_fbc,nIBM_r_fbc)
!!!-----------------------------------------------
!    Membrane solver here
 IF(monitorON) write(*,*) ' SOLVING EQ'
   !write(*,*) nrIBM_fsh,nqIBM_fsh,nIBM_fsh
 call tick(rtime,dt,iter,monitorON,ntsave, &
      vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), &
      vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), &
      vel_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3),  &
      vel_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1), &
      vel_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2), &
      vel_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3),  &
      vel_esh(1:nrIBM_esh,1:nIBM_esh,1), &
      vel_esh(1:nrIBM_esh,1:nIBM_esh,2), &
      vel_esh(1:nrIBM_esh,1:nIBM_esh,3),  &
      vel_fbc(1:nrIBM_fbc,1:nIBM_fbc,1), &
      vel_fbc(1:nrIBM_fbc,1:nIBM_fbc,2), &
      vel_fbc(1:nrIBM_fbc,1:nIBM_fbc,3),  &
      force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,1), &
      force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,2), &
      force_fib(1:nrIBM_fib,1:nsecIBMmax,1:nIBM_fib,3), &
      force_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1), &
      force_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2), &
      force_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3), &
      force_esh(1:nrIBM_esh,1:nIBM_esh,1), &
      force_esh(1:nrIBM_esh,1:nIBM_esh,2), &
      force_esh(1:nrIBM_esh,1:nIBM_esh,3), &
      force_fbc(1:nrIBM_fbc,1:nIBM_fbc,1), &
      force_fbc(1:nrIBM_fbc,1:nIBM_fbc,2), &
      force_fbc(1:nrIBM_fbc,1:nIBM_fbc,3), &
      coord_fib(1:nrIBM_fib,1:nIBM_fib,1),   &
      coord_fib(1:nrIBM_fib,1:nIBM_fib,2),   &
      coord_fib(1:nrIBM_fib,1:nIBM_fib,3),   &
      coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,1),   &
      coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,2),   &
      coord_fsh(1:nrIBM_fsh,1:nqIBM_fsh,1:nIBM_fsh,3),   &
      coord_esh(1:nrIBM_esh,1:nIBM_esh,1),   &
      coord_esh(1:nrIBM_esh,1:nIBM_esh,2),   &
      coord_esh(1:nrIBM_esh,1:nIBM_esh,3),   &
      coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,1),   &
      coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,2),   &
      coord_fbc(1:nrIBM_fbc,1:nIBM_fbc,3),   &
      theboss)
!!!-----------------------------------------------

 IF ( irestart == 1 .and. theboss) THEN
     temp_ibm2=0.0
    call write_restart_DISRIBM(temp_ibm2,1,0)
 ENDIF
#if (mpi_activate==1)
 CALL mpi_barrier(MPI_COMM_WORLD, ierr1)
#endif


 IF(monitorON) write(*,*) ' DISTRIBUTE FORCE '
 if(nrIBM_fib>=1) then
  call cal_fxyz_fib(&
   idim1,&
   coord_fib,&
   force_fib,&
   ds_fib,&
   dtypeDelta(1:3),&
   dsecIBMy,&
   dsecIBMz,&
   nrIBM_fib,&
   nsecIBMmax,&
   nsecIBM,&
   nIBM_fib,&
   nIBM_r_fib)
  end if
 if(nrIBM_fsh>=1) call cal_fxyz_fsh(idim1,coord_fsh,force_fsh,ds_fsh,dtypeDelta(1:3),nrIBM_fsh,nIBM_fsh,nqIBM_fsh,nIBM_r_fsh,nIBM_q_fsh)
 if(nrIBM_esh>=1) call cal_fxyz_esh(idim1,coord_esh,force_esh,ds_esh,dtypeDelta(1:3),nrIBM_esh,nIBM_esh,nIBM_r_esh)
 if(nrIBM_fbc>=1) call cal_fxyz_fbc(idim1,coord_fbc,force_fbc,ds_fbc,dtypeDelta(1:3),nrIBM_fbc,nIBM_fbc,nIBM_r_fbc)
end subroutine


subroutine deltao_fun(delta_type,r,f)

 implicit none
 integer delta_type 
 real*8 r,f,pi1,fi_ib_6,r2,r3,rp,K
 
 pi1 = 4. * atan(1.0)
 f=0.0

 select case (delta_type)

  case(1)
   r=abs(r)
   if (r .lt. 2.d0) f=0.25d0*(1.0d0+dcos(pi1*r/2.d0))

  case(2)
   r=abs(r)
   if (r .lt. 1.d0) then
    f=0.125*(3.-2.*r+sqrt(1.+4.*r-4.*r**2.))
   else if (r .le. 2.d0) then
    f=0.125*(5.-2.*r-sqrt(-7.+12.*r-4.*r**2.))
   end if

   case(3)
    if (r .le. 1.d0) then
             r2=r
             fi_ib_6=61./112.-11./42.*r2-11./56.*r2**2.+1./12.*r2**3 &
               +sqrt(3.)/336.*sqrt(243.+1584.*r2-748.*r2**2-1560.*r2**3 &
               +500.*r2**4+336.*r2**5-112.*r2**6)
              f=fi_ib_6
    else if (r .le. 2.d0) then
             r2=r-1.0
             fi_ib_6=61./112.-11./42.*r2-11./56.*r2**2.+1./12.*r2**3 &
               +sqrt(3.)/336.*sqrt(243.+1584.*r2-748.*r2**2-1560.*r2**3 &
               +500.*r2**4+336.*r2**5-112.*r2**6)
              f=21./16.+7./12.*r-7./8.*r**2+1./6.*r**3-1.5*fi_ib_6
    else if (r .le. 3.d0) then
             r2=r-2.0
             fi_ib_6=61./112.-11./42.*r2-11./56.*r2**2.+1./12.*r2**3 &
               +sqrt(3.)/336.*sqrt(243.+1584.*r2-748.*r2**2-1560.*r2**3 &
               +500.*r2**4+336.*r2**5-112.*r2**6)
             f=9./8.-23./12.*r+3./4.*r**2-1/12.*r**3+0.5*fi_ib_6  
    end if

  case(4)
   r=abs(r)
   if (r .lt. 1.d0) then
    f=1.-0.5*r-r**2+0.5*r**3
   else if (r .le. 2.d0) then
    f=1.-11./6.*r+r**2-1./6.*r**3
  end if

  case(6)
   K=0.714075092976608   !59.0d0/60.0d0-sqrt(29.0d0)/20.0d0
   if((-3.0 .ge. r) .and. (r .le. 3)) then
      Rp = r - ceiling(r) + 1.0  !!Rp between [0,1] 
      R2 = Rp * Rp; R3 = R2*Rp
      call deltaCK(f,rp,r2,r3,K,r)
    end if  
  case default
    write(*,*) 'wrong delta type'
    stop
 end select
return
end subroutine


subroutine deltaCK(phi,r,R2,R3,K,rc)
implicit none
real::  r,K,r2,r3,rc
real:: alpha,beta,gamma,discr,sgnK,phi

      alpha=28.0

      beta = 9./4. - 1.5 * (K + R2) + (22./3.-7.*K)*R - 7./3.*R3

      gamma = 0.25 * ( 0.5*(161./36. - 59./6.*K + 5.*K*K)*R2 + 1./3.*(-109./24. + 5.*K)*R2*R2 + 5./18.*R3*R3  )

      discr=alpha*gamma*(-4.0D0)+beta**2

      sgnK=-1.0
      if((1.5-K) .gt. 0) sgnK=1.0
    if ( (-3.0 < rc) .and. (rc <= -2.0) ) then
            phi = 1./(2*alpha) * ( -beta + sgnK * sqrt(discr) )
    else if ( (-2.0 < rc) .and. (rc <= -1.0) ) then
            phi = -3./(2*alpha) * ( -beta + sgnK * sqrt(discr) ) - 1./16 + 1./8*( K+(rc+2)*(rc+2) ) + 1./12*(3*K-1)*(rc+2) + 1./12*(rc+2)*(rc+2)*(rc+2)
    else if ( (-1.0 < rc) .and. (rc <= 0.0) ) then
            phi = 2./(2*alpha) * ( -beta + sgnK * sqrt(discr) ) + 1./4 + 1./6*(4-3*K)*(rc+1) - 1./6*(rc+1)*(rc+1)*(rc+1)
    else if ( (0.0 < rc) .and. (rc <= 1.0) ) then
            phi = 2./(2*alpha) * ( -beta + sgnK * sqrt(discr) ) +5./8 - 1./4 * ( K+rc*rc )
    else if ( (1.0 < rc) .and. (rc <= 2.0) ) then
            phi = -3./(2*alpha) * ( -beta + sgnK * sqrt(discr) ) +  1./4 - 1./6*(4-3*K)*(rc-1) + 1./6*(rc-1)*(rc-1)*(rc-1)
    else if ( (2.0 < rc) .and. (rc <= 3.0) ) then
            phi = 1./(2*alpha) * ( -beta + sgnK * sqrt(discr) ) -   1./16 + 1./8*(K+(rc-2)*(rc-2)) - 1./12*(3*K-1)*(rc-2) - 1./12*(rc-2)*(rc-2)*(rc-2)
    else
            phi=0.0
    endif     
return
end subroutine
