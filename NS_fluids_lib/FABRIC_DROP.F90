#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

#include "AMReX_ArrayLim.H"

! #define DEBUG_FABRIC_DROP
#define FABRIC_DROP_PROB_TYPE 7001


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

! probtype==FABRIC_DROP_PROB_TYPE (see run2d/inputs.FABRIC_DROP)
module FABRIC_DROP_MODULE
implicit none 
REAL_T, allocatable, dimension(:,:,:) :: thread_nodes
REAL_T, allocatable, dimension(:) :: thread_radius
INTEGER_T, allocatable, dimension(:) :: num_nodes
INTEGER_T num_threads

REAL_T, allocatable, dimension(:,:,:) :: internal_thread_ls
REAL_T internal_dx(3)
contains

! do any initial preparation needed
subroutine INIT_FABRIC_DROP_MODULE()
use probcommon_module
IMPLICIT NONE
REAL_T R
REAL_T internal_x(3)
INTEGER_T ithread
INTEGER_T i,j,k
INTEGER_T N, N_max
character*40 tr_file

#ifdef DEBUG_FABRIC_DROP
INTEGER_T inode
#endif

! Reading and storing thread files
100  FORMAT('../threads/thread_',I3.3,'.txt')

if ((probtype.ne.FABRIC_DROP_PROB_TYPE).or.(SDIM.ne.3)) then
 print *,"probtype or SDIM invalid!"
 stop
endif

num_threads = IDNINT(xblob2)
N_max=0

do ithread=1,num_threads
 WRITE(tr_file,100) ithread
 OPEN(unit=101,&
  file=tr_file,&
  form="formatted",&
  status='old',&
  action='read')
 READ(101,*) N, R

 if(N.gt.N_max) then
  N_max=N
 endif
 CLOSE(101) 
enddo ! ithread



allocate(thread_nodes(num_threads,N_max,3))
allocate(thread_radius(num_threads))
allocate(num_nodes(num_threads))

do ithread=1,num_threads
 WRITE(tr_file,100) ithread
 OPEN(unit=101,&
  file=tr_file,&
  form="formatted",&
  status='old',&
  action='read')
 READ(101,*) num_nodes(ithread), thread_radius(ithread)
 READ(101,*) ((thread_nodes(ithread,i,j), j=1,3), i=1,num_nodes(ithread))
 CLOSE(101) 
enddo ! ithread



! Evaluate level set function on internal fine mesh
allocate(internal_thread_ls(IDNINT(xblob5),IDNINT(yblob5),IDNINT(zblob5)))

internal_dx(1)=(xblob4-xblob3)/xblob5
internal_dx(2)=(yblob4-yblob3)/yblob5
internal_dx(3)=(zblob4-zblob3)/zblob5

print *,"Internal dx for thread LS calculation:", internal_dx
print *,"Calculating thread level set on internal mesh..."

do k=1,IDNINT(zblob5)
#ifdef DEBUG_FABRIC_DROP
 print *,"k=", k
#endif
 do j=1,IDNINT(yblob5)
  do i=1,IDNINT(xblob5)
   internal_x(1)=xblob3+(i-half)*internal_dx(1)
   internal_x(2)=yblob3+(j-half)*internal_dx(2)
   internal_x(3)=zblob3+(k-half)*internal_dx(3)

   internal_thread_ls(i,j,k) = DIST_THREADS(internal_x)
  enddo
 enddo
enddo

print *,"... done!"


return
end subroutine INIT_FABRIC_DROP_MODULE


! fluids tessellate the domain, solids are immersed. 
subroutine FABRIC_DROP_LS(x,t,LS,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(out) :: LS(nmat)

REAL_T x_d(3),x_0(3),x_p(3)
INTEGER_T ind(3)
REAL_T c_000, c_001, c_010, c_011, c_100, c_101, c_110, c_111
REAL_T c_00, c_01, c_10, c_11
REAL_T c_0, c_1, c

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif


LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)
LS(2)=-LS(1)

! Direct thread level set calculation 
! LS(3)= DIST_THREADS(x)

! Approximate level set calcualtion from internal fine mesh
if ((x(3).lt.(zblob3+internal_dx(3))).or.(x(3).gt.(zblob4-internal_dx(3)))) then
 ! use large level set value
 LS(3)=-abs(radblob3)
else
 ! assuming [xblob3,yblob3]:[xblob4,yblob4] is the repeating pattern
 x_p(1)=modulo(x(1),xblob4)
 x_p(2)=modulo(x(2),yblob4)
 x_p(3)=x(3)

 ! ! For 'x' as object point 
 ! ! interpolate from internal fine mesh
 ! ! low side cell index
 ! ind(1)=IDINT((x(1)-xblob3)/internal_dx(1))
 ! ind(2)=IDINT((x(2)-yblob3)/internal_dx(2)) 
 ! ind(3)=IDINT((x(3)-zblob3)/internal_dx(3))

 ! ! x is in the cube {ind, ind+1}
 ! ! trilinear interpolation
 ! ! (https://en.wikipedia.org/wiki/Trilinear_interpolation)
 ! x_0(1)=xblob3+(ind(1)-half)*internal_dx(1)
 ! x_0(2)=yblob3+(ind(2)-half)*internal_dx(2)
 ! x_0(3)=zblob3+(ind(3)-half)*internal_dx(3)
 
 ! x_d=(x-x_0)/internal_dx ! component-wise division

 ! c_000=internal_thread_ls(ind(1)+0, ind(2)+0, ind(3)+0)
 ! c_100=internal_thread_ls(ind(1)+1, ind(2)+0, ind(3)+0)
 ! c_010=internal_thread_ls(ind(1)+0, ind(2)+1, ind(3)+0)
 ! c_110=internal_thread_ls(ind(1)+1, ind(2)+1, ind(3)+0)
 ! c_001=internal_thread_ls(ind(1)+0, ind(2)+0, ind(3)+1)
 ! c_101=internal_thread_ls(ind(1)+1, ind(2)+0, ind(3)+1)
 ! c_011=internal_thread_ls(ind(1)+0, ind(2)+1, ind(3)+1)
 ! c_111=internal_thread_ls(ind(1)+1, ind(2)+1, ind(3)+1)

 ! c_00=c_000*(one-x_d(1)) + c_100*x_d(1)
 ! c_01=c_001*(one-x_d(1)) + c_101*x_d(1)
 ! c_10=c_010*(one-x_d(1)) + c_110*x_d(1)
 ! c_11=c_011*(one-x_d(1)) + c_111*x_d(1)

 ! c_0=c_00*(one-x_d(2)) + c_10*x_d(2)
 ! c_1=c_01*(one-x_d(2)) + c_11*x_d(2)

 ! c=c_0*(one-x_d(3)) + c_1*x_d(3)

 ! LS(3)=c
 
 ! For 'x_p' as object point 
 ! interpolate from internal fine mesh
 ! low side cell index
 ind(1)=IDINT((x_p(1)-xblob3)/internal_dx(1))
 ind(2)=IDINT((x_p(2)-yblob3)/internal_dx(2)) 
 ind(3)=IDINT((x_p(3)-zblob3)/internal_dx(3))

 ! x is in the cube {ind, ind+1}
 ! trilinear interpolation
 ! (https://en.wikipedia.org/wiki/Trilinear_interpolation)
 x_0(1)=xblob3+(ind(1)-half)*internal_dx(1)
 x_0(2)=yblob3+(ind(2)-half)*internal_dx(2)
 x_0(3)=zblob3+(ind(3)-half)*internal_dx(3)
 
 x_d=(x_p-x_0)/internal_dx ! component-wise division

 c_000=internal_thread_ls(ind(1)+0, ind(2)+0, ind(3)+0)
 c_100=internal_thread_ls(ind(1)+1, ind(2)+0, ind(3)+0)
 c_010=internal_thread_ls(ind(1)+0, ind(2)+1, ind(3)+0)
 c_110=internal_thread_ls(ind(1)+1, ind(2)+1, ind(3)+0)
 c_001=internal_thread_ls(ind(1)+0, ind(2)+0, ind(3)+1)
 c_101=internal_thread_ls(ind(1)+1, ind(2)+0, ind(3)+1)
 c_011=internal_thread_ls(ind(1)+0, ind(2)+1, ind(3)+1)
 c_111=internal_thread_ls(ind(1)+1, ind(2)+1, ind(3)+1)

 c_00=c_000*(one-x_d(1)) + c_100*x_d(1)
 c_01=c_001*(one-x_d(1)) + c_101*x_d(1)
 c_10=c_010*(one-x_d(1)) + c_110*x_d(1)
 c_11=c_011*(one-x_d(1)) + c_111*x_d(1)

 c_0=c_00*(one-x_d(2)) + c_10*x_d(2)
 c_1=c_01*(one-x_d(2)) + c_11*x_d(2)

 c=c_0*(one-x_d(3)) + c_1*x_d(3)

 LS(3)=c

endif

! print *,"X, LS", x, LS

return
end subroutine FABRIC_DROP_LS

subroutine FABRIC_DROP_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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

if (adv_dir.eq.SDIM) then
  ! material 3 is the fabric
  ! velsolid_flag==1 if initializing the solid velocity.
 if ((LS(3).ge.zero).or.(velsolid_flag.eq.1)) then
  ! in solid
  do dir=1,SDIM
   VEL(dir)=zero
  enddo
 else if ((LS(3).le.zero).and. &
          (velsolid_flag.eq.0)) then
  ! material 1 is the drop
  if ((LS(1).ge.zero).or. &
      (LS(1).ge.-two*dx(1))) then
   ! in drop
   do dir=1,SDIM
    VEL(dir)=zero
   enddo
   VEL(SDIM)=-abs(advbot)
   ! material 2 is the gas
  else if (LS(2).ge.zero) then 
   do dir=1,SDIM
    VEL(dir)=zero
   enddo
  else
   print *,"LS bust"
   stop
  endif
 else
  print *,"LS(3) or velsolid bust"
  stop
 endif

else
 print *,"expecting adv_dir = SDIM in FABRIC_DROP"
 stop
endif

return 
end subroutine FABRIC_DROP_VEL

function DIST_THREADS(P)
 ! Returns the signed distance function to the thrads' surfaces
 ! Inside the threads > 0
 ! Outside the threads < 0
 IMPLICIT NONE
 REAL_T DIST_THREADS
 REAL_T P(SDIM)
 REAL_T dist,seg_dist,thr_dist
 INTEGER_T ithread,inode

 ! Iterate over threads
 ! Iterrate over thread nodes-1
 ! Find distance to line segment/thread surface
 dist=-1.0d15
 do ithread = 1,num_threads
  thr_dist=1.0d15
  do inode = 1,num_nodes(ithread)-1
   seg_dist=DIST_SEGMENT(P,&
    thread_nodes(ithread,inode,:),&
    thread_nodes(ithread,inode+1,:))
   thr_dist=min(thr_dist,seg_dist)
  enddo ! num_nodes

  thr_dist=-(thr_dist-thread_radius(ithread))

  if(abs(thr_dist).lt.abs(dist)) then
   dist=thr_dist
  endif
 enddo ! num_threads

 DIST_THREADS=dist

 return
end function DIST_THREADS

function DIST_SEGMENT(P,M,N)
 ! Returns the distance o point P from line segment MN
 IMPLICIT NONE
 REAL_T DISt_SEGMENT
 REAL_T P(SDIM)
 REAL_T M(SDIM)
 REAL_T N(SDIM)
 REAL_T NM_unit(SDIM)
 REAL_T LMN,LMP,L,dist

 LMN=sqrt(dot_product(N-M,N-M))
 NM_unit=(N-M)/LMN
 L=dot_product(P-M,NM_unit)

 if ((L.ge.zero).and.(L.le.LMN)) then
  LMP=sqrt(dot_product(P-M,P-M))
  dist=sqrt(LMP**2-L**2)
 elseif (L.lt.zero) then
  dist=sqrt(dot_product(P-M,P-M))
 elseif (L.gt.LMN) then
  dist=sqrt(dot_product(P-N,P-N))
 endif

 DIST_SEGMENT=dist
 
 return
end function DIST_SEGMENT
  
  
 


! ! These next routines only used for compressible materials.
! !***********************************************
! ! compressible material functions for (ns.material_type = 24)
! subroutine EOS_FABRIC_DROP(rho,internal_energy,pressure, &
!   imattype,im)
!  IMPLICIT NONE
!  INTEGER_T, intent(in) :: imattype,im
!  REAL_T, intent(in) :: rho
!  REAL_T, intent(in) :: internal_energy
!  REAL_T, intent(out) :: pressure

!  if (imattype.eq.24) then
!   pressure=zero
!   print *,"Not expecting compressible material! "
!  else
!   print *,"imattype invalid EOS_FABRIC_DROP"
!   stop
!  endif

!  return
! end subroutine EOS_FABRIC_DROP

! subroutine SOUNDSQR_FABRIC_DROP(rho,internal_energy,soundsqr, &
!   imattype,im)
!  IMPLICIT NONE
!  INTEGER_T, intent(in) :: imattype,im
!  REAL_T, intent(in) :: rho
!  REAL_T, intent(in) :: internal_energy
!  REAL_T, intent(out) :: soundsqr
!  REAL_T pressure

!  if (imattype.eq.24) then
!   call EOS_FABRIC_DROP(rho,internal_energy,pressure,imattype,im)
!   soundsqr=zero
!   print *,"Not expecting compressible material! "
!   stop

!  else
!   print *,"imattype invalid SOUNDSQR_FABRIC_DROP"
!   stop
!  endif

!  return
! end subroutine SOUNDSQR_FABRIC_DROP

! subroutine INTERNAL_FABRIC_DROP(rho,temperature, &
!   local_internal_energy, &
!   imattype,im)
!  use global_utility_module
!  IMPLICIT NONE
!  INTEGER_T, intent(in) :: imattype,im
!  REAL_T, intent(in) :: rho
!  REAL_T, intent(in) :: temperature 
!  REAL_T, intent(out) :: local_internal_energy

!  call INTERNAL_default(rho,temperature,local_internal_energy, &
!         imattype,im)

!  return
! end subroutine INTERNAL_FABRIC_DROP

! subroutine TEMPERATURE_FABRIC_DROP(rho,temperature,internal_energy, &
!   imattype,im)
!  use global_utility_module
!  IMPLICIT NONE
!  INTEGER_T, intent(in) :: imattype,im
!  REAL_T, intent(in) :: rho
!  REAL_T, intent(out) :: temperature 
!  REAL_T, intent(in) :: internal_energy

!  call TEMPERATURE_default(rho,temperature,internal_energy, &
!         imattype,im)

!  return
! end subroutine TEMPERATURE_FABRIC_DROP

! This routine will not effect the simulation since
! all of the domain BC will be no-slip.


subroutine FABRIC_DROP_PRES(x,t,LS,PRES,nmat)
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
PRES=outflow_pressure

return 
end subroutine FABRIC_DROP_PRES


subroutine FABRIC_DROP_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
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

if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.FABRIC_DROP_PROB_TYPE)) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  STATE(ibase+1)=fort_denconst(im)
  if (t.eq.zero) then
   STATE(ibase+2)=fort_initial_temperature(im)
  else if (t.gt.zero) then
   STATE(ibase+2)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif
  do n=1,num_species_var
   STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine FABRIC_DROP_STATE

 ! dir=1..sdim  side=1..2
subroutine FABRIC_DROP_LS_BC(xwall,xghost,t,LS, &
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
REAL_T, intent(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call FABRIC_DROP_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine FABRIC_DROP_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine FABRIC_DROP_VEL_BC(xwall,xghost,t,LS, &
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

 call FABRIC_DROP_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine FABRIC_DROP_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine FABRIC_DROP_PRES_BC(xwall,xghost,t,LS, &
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

 call FABRIC_DROP_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine FABRIC_DROP_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine FABRIC_DROP_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T local_STATE(nmat*num_state_material)
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
 call FABRIC_DROP_STATE(xghost,t,LS,local_STATE, &
         local_bcflag,nmat,num_state_material)
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
end subroutine FABRIC_DROP_STATE_BC

subroutine FABRIC_DROP_HEATSOURCE( &
     im,VFRAC, &
     time, &
     x, &
     xsten, & ! xsten(-nhalf:nhalf,SDIM)
     nhalf, &
     temp, &
     heat_source,den,CV,dt, &
     nmat)
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

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((num_materials.eq.3).and.(probtype.eq.FABRIC_DROP_PROB_TYPE)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine FABRIC_DROP_HEATSOURCE

end module FABRIC_DROP_MODULE
