#undef BL_LANG_CC
#define BL_LANG_FORT

#include "AMReX_BC_TYPES.H"

module integrate_module 
USE GeneralClass ! vof_cisl.F90
USE probmain_module 
USE probcommon_module 
USE LegendreNodes
USE global_utility_module 
USE MOF_pair_module ! vfrac_pair.F90
USE mmat_FVM  ! multimat_FVM.F90
USE bicgstab_module ! BICGSTAB_Yang_MULTI.F90
USE supercooled_exact_sol
USE variable_temperature_drop
USE tsat_module

IMPLICIT NONE

contains

subroutine integrate_steps( &
 nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in, &
 y_fluxtest1, &
 y_fluxtest2, &
 TSTOP, &
 TDIFF_in, &
 subcycling_step, &
 rstefan, &
 M_START,max_front_vel, &
 iter,iter_average, &
 isink, &
 flxavg1,flxavg2, &
 Ts, &
 deltat_in,deltat_polar, &
 sdim_in,N_CURRENT, &
 plot_int, &
 mofdata_FAB_in,T_STATE,local_state_ncomp, &
 local_operator_internal,local_operator_external,local_linear_exact, &
 probtype_in,ngeom_recon_in, &
 h_in, &
 nmat_in, &
 T_new, &
 M_CURRENT,M_MAX_TIME_STEP,fixed_dt_main,dx_in, &
 local_nten,stefan_flag,xCC,yCC)
IMPLICIT NONE

INTEGER, intent(in) :: M_MAX_TIME_STEP
real(kind=8), intent(inout), dimension(1:M_MAX_TIME_STEP+1) :: Ts
INTEGER :: tm
INTEGER, intent(inout) ::  nx_in,ny_in
INTEGER, intent(in) ::  plot_int
INTEGER, intent(in) ::  N_CURRENT
INTEGER, intent(in) ::  sdim_in
REAL(kind=8), intent(in) :: TDIFF_in
REAL(kind=8), intent(inout) :: deltat_in
REAL(kind=8), intent(in) :: TSTOP
INTEGER, intent(in) :: nmat_in
INTEGER, intent(in) :: ngeom_recon_in
integer, intent(in) :: local_state_ncomp
! -1:N,-1:N,nmat*ngeom_recon
real(kind=8),dimension(-1:N_CURRENT,-1:N_CURRENT,ngeom_recon_in*nmat_in), &
  intent(inout) :: mofdata_FAB_in
! -1:N,-1:N,nmat
real(kind=8),dimension(-1:N_CURRENT,-1:N_CURRENT,local_state_ncomp), &
  intent(inout) :: T_STATE
integer, intent(in) :: local_operator_internal
integer, intent(in) :: local_operator_external
integer, intent(in) :: local_linear_exact
integer, intent(in) :: probtype_in
REAL(KIND=8), intent(in) :: h_in
REAL(KIND=8), intent(in) :: fixed_dt_main
real(kind=8), dimension(:,:,:), allocatable :: VFRAC_MOF_in
INTEGER, intent(in) :: M_CURRENT
real(kind=8), intent(in) :: dx_in(sdim_in)
real(kind=8), intent(inout) :: rstefan
integer, intent(in) :: local_nten
INTEGER, intent(in) :: subcycling_step
INTEGER, intent(inout) :: stefan_flag ! VARIABLE TSAT
!real(kind=8),dimension(-1:N) :: xCC,yCC       
real(kind=8),dimension(-1:N_CURRENT), intent(in) :: xCC,yCC       ! cell centers
real(kind=8),dimension(-1:N_CURRENT,-1:N_CURRENT,local_state_ncomp), &
  intent(inout) :: T_new

real(kind=8), dimension(:,:,:), allocatable :: UNEW_in
real(kind=8), dimension(:,:,:), allocatable :: UOLD_in
real(kind=8), dimension(:,:,:), allocatable :: beta_in
REAL(kind=8), intent(in) :: deltat_polar
integer :: finished_flag
real(kind=8), intent(inout) :: flxavg1,flxavg2
integer, intent(in) :: isink
integer, intent(in) :: M_START
integer, intent(inout) :: iter
real(kind=8), intent(inout) :: iter_average
real(kind=8), intent(inout) :: max_front_vel
real(kind=8), intent(inout) :: y_fluxtest1
real(kind=8), intent(inout) :: y_fluxtest2
integer, intent(inout) :: lox_in,loy_in,hix_in,hiy_in
REAL(kind=8) :: bicgstab_tol_in
REAL(kind=8) :: current_time_in
real(kind=8) :: eff_radius
real(kind=8) :: T_FIELD
real(kind=8) :: expect_radius
real(kind=8) :: flxtot1,flxtot2
real(kind=8) :: test_vel
REAL(KIND=8) :: time_n,time_np1
real(kind=8) :: TSAT
real(kind=8) :: xgrid,ygrid
real(kind=8) :: xcen,ycen
real(kind=8) :: xlo_fluxtest,xhi_fluxtest
integer hflag
integer i,j
integer imof
integer im,im1,im_opp
integer iten
integer ireverse
integer icen,jcen
integer j_fluxtest,ilo_fluxtest,ihi_fluxtest,isum
integer vofcomp
integer vofcomp2
real(kind=8) :: LL
real(kind=8) :: voltotal
real(kind=8) :: local_Pi
integer nsteps
integer total_nsteps_parm
INTEGER :: precond_type_in
integer scomp
real(kind=8) :: stefan_time
real(kind=8) :: sum_alpha
real(kind=8) :: sumvf,sumvf2
REAL(kind=8) :: alpha_in(100)
real(kind=8) :: local_dist

integer local_state_ncomp_test
integer allocate_flag

! temperature, velocity, interface reconstruction, level set
local_state_ncomp_test=nmat_in+local_nten*sdim_in+ &
    ngeom_recon_in*nmat_in+nmat_in*(sdim_in+1)
if (local_state_ncomp_test.eq.local_state_ncomp) then
 ! do nothing
else
 print *,"local_state_ncomp_test invalid"
 stop
endif

if (M_MAX_TIME_STEP.ge.M_CURRENT) then
 ! do nothing
else
 print *,"M_MAX_TIME_STEP or M_CURRENT invalid"
 stop
endif

local_Pi=4.0d0*atan(1.0d0)

tm=1
finished_flag=0
do while (finished_flag.eq.0)

 current_time_in=Ts(tm) ! t^{n} (Ts(i)=(i-1) * deltat)
 nsteps=tm-1 ! NSTEPS

 print *,"STEP (>=1), TIME, DT ",tm,current_time_in,deltat_in

 nx_in=N_CURRENT
 ny_in=N_CURRENT
 bicgstab_tol_in=1.0D-10
 precond_type_in=1 ! 0 M=I  1=Jacobi precond.
 hflag=0
 do im=1,nmat_in
  alpha_in(im)=fort_heatviscconst(im)
 enddo
 lox_in=0
 loy_in=0
 hix_in=nx_in-1
 hiy_in=ny_in-1
 allocate(UNEW_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,local_state_ncomp)) 
 allocate(UOLD_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,local_state_ncomp)) 
 allocate(beta_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)) 
 allocate(VFRAC_MOF_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1, &
    nmat_in)) 

 if (ngeom_recon_in.eq.2*sdim_in+3) then
  ! do nothing
 else
  print *,"ngeom_recon_in invalid"
  stop
 endif

 do i=lox_in-1,hix_in+1
 do j=loy_in-1,hiy_in+1
  do im=1,nmat_in
   xgrid=(i+0.5)*h_in
   ygrid=(j+0.5)*h_in

   sumvf=0.0
   sum_alpha=0.0
   do im1=1,nmat_in
    vofcomp=ngeom_recon_in*(im1-1)+1
    sumvf=sumvf+mofdata_FAB_in(i,j,vofcomp)
    sum_alpha=sum_alpha+mofdata_FAB_in(i,j,vofcomp)/ &
            (alpha_in(im1)+1.0E-10)
   enddo
   sum_alpha=sum_alpha/sumvf
   sum_alpha=1.0/sum_alpha 
   beta_in(i,j,im)=sum_alpha

   vofcomp=ngeom_recon_in*(im-1)+1
   VFRAC_MOF_in(i,j,im)=mofdata_FAB_in(i,j,vofcomp)
  enddo ! im=1..nmat_in
  do im=1,local_state_ncomp
   UNEW_in(i,j,im)=T_STATE(i,j,im)
   UOLD_in(i,j,im)=T_STATE(i,j,im)
  enddo
 enddo
 enddo

 if (tm.eq.1) then
  call init_tsatfab(N_CURRENT) ! VARIABLE TSAT, allocate swept too
 endif

  ! UOLD=UOLD_in  UNEW=UNEW_in
  ! in: BICGSTAB_Yang_MULTI.F90
 allocate_flag=1
 call INIT_GLOBALS( &
  allocate_flag, &
  nsteps, & ! NSTEPS
  local_state_ncomp, &
  local_operator_internal, &
  local_operator_external, &
  local_linear_exact, &
  probtype_in, &
  sdim_in,ngeom_recon_in, &
  nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in, &
  UNEW_in,UOLD_in, &
  beta_in,h_in,precond_type_in,bicgstab_tol_in, &
  VFRAC_MOF_in,nmat_in,alpha_in,deltat_in, &
  mofdata_FAB_in,current_time_in)

 time_n=current_time_in
 time_np1=current_time_in+deltat_in

 if (tm.eq.1) then

     ! output_solution declared in: BICGSTAB_Yang_MULTI.F90
  if (fixed_dt_main.eq.0.0d0) then
   total_nsteps_parm=M_MAX_TIME_STEP
  else
   total_nsteps_parm=M_CURRENT
  endif
  call output_solution(UNEW_in,time_n,nsteps,plot_int, &
          total_nsteps_parm, &
          fixed_dt_main)

 endif

  ! in: MOVE_INTERFACE_3D.F90
  ! interface updated here: 
  ! input: UOLD
  ! output: UNEW  (has reconstructed interface)
  ! update swept and interface temperaure tsatfab
 call update_interface(UOLD,UNEW,N_CURRENT,local_state_ncomp, &
   dx_in,time_n,deltat_in,nsteps,local_nten,stefan_flag)

  ! hflag==0
 call set_boundary(UNEW,0,local_state_ncomp)

 do i=lox_in-1,hix_in+1
 do j=loy_in-1,hiy_in+1
  do im=1,local_state_ncomp
   UNEW_in(i,j,im)=UNEW(i,j,im)
   UOLD_in(i,j,im)=UNEW(i,j,im)
   UOLD(i,j,im)=UOLD_in(i,j,im)
  enddo
 enddo
 enddo

 do i= -1,N_CURRENT
  do j= -1,N_CURRENT
   do imof=1,ngeom_recon_in*nmat_in
    scomp=nmat_in+local_nten*sdim_in+imof
    mofdata_FAB_in(i,j,imof)=UNEW_in(i,j,scomp)
    mofdata_FAB(i,j,imof)=UNEW_in(i,j,scomp)
   enddo
   do im=1,nmat_in
    vofcomp=ngeom_recon_in*(im-1)+1
    VFRAC_MOF_in(i,j,im)=mofdata_FAB_in(i,j,vofcomp)
    VFRAC_MOF(i,j,im)=mofdata_FAB_in(i,j,vofcomp)
   enddo
  enddo
 enddo

  ! Dirichlet BC use t^n+1 data
  ! polar solver called before bicgstab is called.
 if (probtype_in.eq.19) then
  print *,"deltat_polar",deltat_polar                            
  print *,"subcycling_step", subcycling_step    
  print *,"Np,Mp= ",Np,Mp
  do i=1,subcycling_step                               
   call polar_2d_heat(sdim_in,Np,Mp,fort_heatviscconst(2),deltat_polar, &
    r_polar,z_polar,dr_polar,dz_polar,upolar)
   if (i.eq.(i/1000)*1000) then
    print *,"subcycling_step number: i=",i
   endif
  enddo
 endif

  ! UOLD=UOLD_in  UNEW=UNEW_in
  ! in: BICGSTAB_Yang_MULTI.F90
  ! uses tsatfab and swept
 allocate_flag=0
 call INIT_GLOBALS( &
  allocate_flag, &
  nsteps, & ! NSTEPS
  local_state_ncomp, &
  local_operator_internal, &
  local_operator_external, &
  local_linear_exact, &
  probtype_in, &
  sdim_in,ngeom_recon_in, &
  nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in, &
  UNEW_in,UOLD_in, &
  beta_in,h_in,precond_type_in,bicgstab_tol_in, &
  VFRAC_MOF_in,nmat_in,alpha_in,deltat_in, &
  mofdata_FAB_in,current_time_in)

  ! uses tsatfab and swept
 call bicgstab(UNEW_in,hflag,iter)

 iter_average=iter_average+iter

  ! hflag==0
 call set_boundary(UNEW_in,0,local_state_ncomp)
 do i=lox_in-1,hix_in+1
 do j=loy_in-1,hiy_in+1
  do im=1,local_state_ncomp
   UNEW(i,j,im)=UNEW_in(i,j,im)
   UOLD_in(i,j,im)=UNEW_in(i,j,im)
   UOLD(i,j,im)=UOLD_in(i,j,im) !UOLD declared in bicgstab_module:INIT_GLOBALS
  enddo
 enddo
 enddo

 if (probtype_in.eq.0) then
  ! do nothing
 else if (probtype_in.eq.1) then
  ! do nothing
 else if (probtype_in.eq.2) then
  ! do nothing
 else if (probtype_in.eq.3) then
  ! do nothing
 else if (probtype_in.eq.4) then
  call axisymmetric_disk_advance(deltat_in)
 else if (probtype_in.eq.400) then
  ! do nothing
 else if (probtype_in.eq.401) then
  ! do nothing
 else if (probtype_in.eq.402) then
  ! do nothing
 else if (probtype_in.eq.403) then
  ! do nothing
 else if (probtype_in.eq.404) then
  ! do nothing
 else if (probtype_in.eq.5) then
  ! do nothing
 else if (probtype_in.eq.19) then   ! annulus cvg test
  ! do nothing (polar solver called before bicgstab is called)
 elseif (probtype_in.eq.13)then
  ! do nothing
 elseif(probtype_in.eq.14)then
  ! do nothing
 elseif(probtype_in.eq.15)then
  ! do nothing
 else if (probtype_in.eq.16) then
  call axisymmetric_disk_advance(deltat_in)
 elseif(probtype_in.eq.17)then
  ! do nothing
 elseif(probtype_in.eq.20)then
  ! do nothing
 else
  print *,"probtype_in invalid"
  stop
 endif 

 do i= -1,N_CURRENT
 do j= -1,N_CURRENT
  do im=1,nmat_in
   xcen=xCC(i)
   ycen=yCC(j)
   if (probtype_in.eq.0) then
    ! do nothing
   else if (probtype_in.eq.2) then
    ! do nothing
   else if (probtype_in.eq.1) then
    ! do nothing
   else if (probtype_in.eq.3) then
    rstefan=sqrt((xcen-xblob)**2+(ycen-yblob)**2)
    if (rstefan.ge.0.5d0-radblob2) then
     stefan_time=fort_time_radblob(2)+Ts(tm+1)
     call liquid_temperature_driver( &
      rstefan, &
      stefan_time, &
      T_FIELD)
     UNEW(i,j,im)=T_FIELD
    endif
   else if (probtype_in.eq.4) then
    ! do nothing - this is a shrinking disk, outside temperature
    ! is uniform, grad T dot n=0 on the outer walls.
   else if (probtype_in.eq.400) then
    ! do nothing
   else if (probtype_in.eq.401) then
    ! do nothing
   else if (probtype_in.eq.402) then
    ! do nothing
   else if (probtype_in.eq.403) then ! dendrite

    if (transition_region.gt.0.0d0) then
     local_dist=sqrt((xcen-xblob)**2+(ycen-yblob)**2)
     if (local_dist.ge.xblob-1.0d0/16.0d0) then
      UNEW(i,j,im)=fort_tempconst(1)
     endif
    else if (transition_region.eq.0.0d0) then
     ! do nothing
    else
     print *,"transition_region invalid"
     stop
    endif

   else if (probtype_in.eq.404) then
    ! do nothing
   else if (probtype_in.eq.5) then
    if (xcen.ge.1.0d0-1.0d0/16.0d0) then !buffer size:1 cell coarse grid
     T_FIELD=272.0d0+exp(-(xcen-0.1d0-Ts(tm+1)))
     UNEW(i,j,im)=T_FIELD
    endif 
   else if (probtype_in.eq.19) then   ! annulus cvg test
    ! do nothing
   elseif (probtype_in.eq.13)then
    ! do nothing
   elseif(probtype_in.eq.14)then
    ! do nothing
   elseif(probtype_in.eq.15)then
    ! do nothing
   elseif(probtype_in.eq.16)then
    ! do nothing
   elseif(probtype_in.eq.17)then
    ! do nothing
   elseif(probtype_in.eq.20)then
    ! do nothing
   else
    print *,"probtype_in invalid5 ",probtype_in
    stop
   endif
  enddo !im=1..nmat_in
 enddo ! j
 enddo ! i

  ! hflag==0
 call set_boundary(UNEW,0,local_state_ncomp)

 do i=lox_in-1,hix_in+1
 do j=loy_in-1,hiy_in+1
  do im=1,local_state_ncomp
   UNEW_in(i,j,im)=UNEW(i,j,im)
   UOLD_in(i,j,im)=UNEW(i,j,im)
   UOLD(i,j,im)=UOLD_in(i,j,im)
  enddo
 enddo
 enddo

 do i= -1,N_CURRENT
  do j= -1,N_CURRENT
   do imof=1,ngeom_recon_in*nmat_in
    scomp=nmat_in+local_nten*sdim_in+imof
    mofdata_FAB_in(i,j,imof)=UNEW_in(i,j,scomp)
    mofdata_FAB(i,j,imof)=UNEW_in(i,j,scomp)
   enddo
   do im=1,nmat_in
    vofcomp=ngeom_recon_in*(im-1)+1
    VFRAC_MOF_in(i,j,im)=mofdata_FAB_in(i,j,vofcomp)
    VFRAC_MOF(i,j,im)=mofdata_FAB_in(i,j,vofcomp)
   enddo
  enddo
 enddo

 nsteps=tm
     ! output_solution declared in: BICGSTAB_Yang_MULTI.F90
 if (fixed_dt_main.eq.0.0d0) then
  if (Ts(tm+1).ge.TSTOP-1.0D-14) then
   total_nsteps_parm=nsteps
  else
   total_nsteps_parm=M_MAX_TIME_STEP
  endif
 else
  total_nsteps_parm=M_CURRENT
 endif
 call output_solution(UNEW_in,time_np1,nsteps,plot_int, &
         total_nsteps_parm, &
         fixed_dt_main)

  ! deallocates UNEW, UOLD, VFRAC_MOF, mofdata_FAB
 call DEALLOCATE_GLOBALS()

 deallocate(UOLD_in) 
 deallocate(beta_in) 
 deallocate(VFRAC_MOF_in) 

 do i=lox_in-1,hix_in+1
 do j=loy_in-1,hiy_in+1
  do im=1,local_state_ncomp
   T_new(i,j,im)=UNEW_in(i,j,im)
  enddo
 enddo
 enddo
 icen=(hix_in+lox_in)/2
 jcen=(hiy_in+loy_in)/2
 if (isink.eq.1) then
  do im=1,nmat
   T_new(icen,jcen,im)=TDIFF_in+saturation_temp(1)
  enddo
 else if (isink.eq.0) then
  ! do nothing
 else
  print *,"isink invalid"
  stop
 endif

 deallocate(UNEW_in) 

 do i= -1,N_CURRENT
 do j= -1,N_CURRENT
 do im=1,local_state_ncomp
  T_STATE(i,j,im) = T_new(i,j,im)
 enddo
 enddo
 enddo

 do im=1,nmat_in
  voltotal=0.0d0
  do i= 0,N_CURRENT-1
   do j= 0,N_CURRENT-1
    vofcomp=nmat_in+local_nten*sdim_in+(im-1)*ngeom_recon_in+1
    sumvf=T_STATE(i,j,vofcomp)
    voltotal=voltotal+sumvf*h_in*h_in
   enddo
  enddo
  print *,"TIME= ",Ts(tm+1)," MAT= ",im," VOLUME= ",voltotal
  eff_radius=sqrt(voltotal/local_Pi)   ! pi r^2 = V  r=(V/pi)^(1/2)
  if (probtype_in.eq.5) then
   if (im.eq.1) then
    eff_radius=voltotal
   else if (im.eq.2) then
    eff_radius=1.0d0-voltotal
   else
    print *,"im invalid 113"
    stop
   endif
  endif
  print *,"TIME= ",Ts(tm+1)," MAT= ",im," EFF RADIUS= ",eff_radius
 enddo ! im=1..nmat_in

 if (probtype_in.eq.0) then
  ! do nothing
 else if (probtype_in.eq.2) then
  ! do nothing
 else if (probtype_in.eq.1) then
  ! do nothing
 else if (probtype_in.eq.3) then
  im=1
  stefan_time=fort_time_radblob(2)+Ts(tm+1)
  call solidification_front_radius_driver(stefan_time,expect_radius) 
  print *,"TIME= ",Ts(tm+1)," MAT= ",im," EXACT RADIUS= ",expect_radius
 else if (probtype_in.eq.4) then
  im=1
  expect_radius=axisymmetric_disk_radblob(2)
  print *,"TIME= ",Ts(tm+1)," MAT= ",im," EXACT RADIUS= ",expect_radius
 else if (probtype_in.eq.400) then
  ! do nothing
 else if (probtype_in.eq.401) then
  ! do nothing
 else if (probtype_in.eq.402) then
  ! do nothing
 else if (probtype_in.eq.403) then
  ! do nothing
 else if (probtype_in.eq.404) then
  ! do nothing
 else if (probtype_in.eq.5) then
  expect_radius=0.1d0+Ts(tm+1)
  print *,"TIME= ",Ts(tm+1)," MAT= ",im," EXACT RADIUS= ",expect_radius
 else if (probtype_in.eq.19) then   ! annulus cvg test
  ! do nothing
 elseif (probtype_in.eq.13)then
  ! do nothing
 elseif(probtype_in.eq.14)then
  ! do nothing
 elseif(probtype_in.eq.15)then
  ! do nothing
 elseif(probtype_in.eq.16)then

  flxtot1=0.0d0 
  flxtot2=0.0d0 

  if ((N_CURRENT.eq.16).or. &
      (N_CURRENT.eq.32).or. &
      (N_CURRENT.eq.64).or. &
      (N_CURRENT.eq.128).or. &
      (N_CURRENT.eq.256).or. &
      (N_CURRENT.eq.512).or. &
      (N_CURRENT.eq.1024)) then

    ! filament_test_type declared in vof_cisl.F90
   if (filament_test_type.eq.0) then ! irregular material 2
    y_fluxtest1=0.40625d0
   else if (filament_test_type.eq.1) then ! circular material 2
    y_fluxtest1=0.59375d0
   else
    print *,"filament_test_type invalid"
    stop
   endif

   j_fluxtest=NINT(y_fluxtest1/dx_in(2))
   xlo_fluxtest=0.375d0
   xhi_fluxtest=0.625d0
   ilo_fluxtest=NINT(xlo_fluxtest/dx_in(1))
   ihi_fluxtest=NINT(xhi_fluxtest/dx_in(1))

   isum=0
   do i=ilo_fluxtest,ihi_fluxtest-1
    flxtot1=flxtot1+(T_STATE(i,j_fluxtest,2)-T_STATE(i,j_fluxtest-1,2))/dx_in(2)
    isum=isum+1
   enddo
   flxtot1=flxtot1/real(isum,8)
   flxavg1=flxavg1+flxtot1
   print *,"TIME,flxtot1 ",Ts(tm+1),flxtot1

   y_fluxtest2=0.84375d0
   j_fluxtest=NINT(y_fluxtest2/dx_in(2))
   xlo_fluxtest=0.0d0
   xhi_fluxtest=1.0d0
   ilo_fluxtest=NINT(xlo_fluxtest/dx_in(1))
   ihi_fluxtest=NINT(xhi_fluxtest/dx_in(1))

   isum=0
   do i=ilo_fluxtest,ihi_fluxtest-1
    flxtot2=flxtot2+(T_STATE(i,j_fluxtest,3)-T_STATE(i,j_fluxtest-1,3))/dx_in(2)
    isum=isum+1
   enddo
   flxtot2=flxtot2/real(isum,8)
   flxavg2=flxavg2+flxtot2
   print *,"TIME,flxtot2 ",Ts(tm+1),flxtot2

  else
   print *,"N_CURRENT out of range N_CURRENT=",N_CURRENT
   stop
  endif

 elseif(probtype_in.eq.17)then
  ! do nothing
 elseif(probtype_in.eq.20)then
  ! do nothing
 else
  print *,"probtype_in invalid20 ",probtype_in
  stop
 endif

 tm=tm+1

 finished_flag=0
 if (Ts(tm).ge.TSTOP-1.0D-14) then
  finished_flag=1
 endif
 if (tm-1.ge.M_MAX_TIME_STEP) then
  finished_flag=1
 endif
 if (tm-1.ge.M_CURRENT) then
  if ((fixed_dt_main.eq.-1.0d0).or. &
      (fixed_dt_main.gt.0.0d0)) then
   if ((finished_flag.ne.1).or. &
       (abs(Ts(tm)-TSTOP).gt.1.0D-14)) then
    print *,"expecting Ts(M_CURRENT+1) == TSTOP"
    stop
   endif
  else if (fixed_dt_main.eq.0.0d0) then
   ! check nothing
  else
   print *,"fixed_dt_main invalid"
   stop
  endif
 else if (tm.ge.2) then
  !check nothing
 else
  print *,"tm invalid"
  stop
 endif

  ! gingerbread man, xue, ice melt, NASA bubble, or dendrite
 if ((probtype_in.eq.400).or. &
     (probtype_in.eq.401).or. &
     (probtype_in.eq.402).or. &
     (probtype_in.eq.403).or. &
     (probtype_in.eq.404)) then 

  max_front_vel=0.0
  do i= 0,N_CURRENT-1
  do j= 0,N_CURRENT-1
   do im=1,nmat_in
    vofcomp=nmat_in+local_nten*sdim_in+(im-1)*ngeom_recon_in+1
    sumvf=T_STATE(i,j,vofcomp)
    if ((sumvf.ge.VOFTOL_REDIST).and. &
        (sumvf.le.1.0d0-VOFTOL_REDIST)) then
     do im_opp=im+1,nmat_in
      vofcomp2=nmat_in+local_nten*sdim_in+(im_opp-1)*ngeom_recon_in+1
      sumvf2=T_STATE(i,j,vofcomp2)
      if ((sumvf2.ge.VOFTOL_REDIST).and. &
          (sumvf2.le.1.0d0-VOFTOL_REDIST)) then
       call get_iten(im,im_opp,iten,nmat_in)
       do ireverse=0,1
        LL=abs(latent_heat(iten+ireverse*local_nten))
        TSAT=saturation_temp(iten+ireverse*local_nten)
         
        if (LL.gt.0.0d0) then
         if (probtype.eq.403) then !dendrite
          TSAT=tsatfab(i,j,local_nten+1)
         endif
         test_vel=fort_heatviscconst(im)*abs(T_STATE(i,j,im)-TSAT)/LL+ &
                  fort_heatviscconst(im_opp)*abs(T_STATE(i,j,im_opp)-TSAT)/LL
         test_vel=test_vel*4.0d0/dx_in(1)
         if (test_vel.gt.max_front_vel) then
          max_front_vel=test_vel
         endif
        else if (LL.eq.0.0d0) then
         ! do nothing
        else
         print *,"LL invalid"
         stop
        endif
       enddo ! ireverse=0..1
      else if ((sumvf2.ge.-VOFTOL_REDIST).and. &
               (sumvf2.le.VOFTOL_REDIST)) then
       ! do nothing
      else if ((sumvf2.ge.1.0d0-VOFTOL_REDIST).and. &
               (sumvf2.le.1.0d0+VOFTOL_REDIST)) then
        ! do nothing
      else
       print *,"sumvf2 invalid"
       stop
      endif
     enddo ! im_opp=1..nmat_in
    else if ((sumvf.ge.-VOFTOL_REDIST).and. &
             (sumvf.le.VOFTOL_REDIST)) then
     ! do nothing
    else if ((sumvf.ge.1.0d0-VOFTOL_REDIST).and. &
             (sumvf.le.1.0d0+VOFTOL_REDIST)) then
     ! do nothing
    else
     print *,"sumvf invalid"
     stop
    endif
   enddo ! im=1..nmat_in
  enddo
  enddo

  if (max_front_vel.gt.0.0d0) then
   deltat_in=h_in*0.25d0/max_front_vel
    ! VERIFICATION
   deltat_in=TSTOP/M_START
   if (finished_flag.eq.0) then

    if (Ts(tm)+deltat_in.ge.TSTOP-1.0D-14) then
     deltat_in=TSTOP-Ts(tm)
    endif
    if (deltat_in.gt.0.0d0) then
     ! do nothing
    else
     print *,"deltat_in invalid"
     stop
    endif
    Ts(tm+1)=Ts(tm)+deltat_in
   else if (finished_flag.eq.1) then
    ! do nothing
   else
    print *,"finished_flag invalid"
    stop
   endif
  else
   print *,"max_front_vel invalid"
   stop
  endif

 else if ((probtype_in.ne.400).and. &
          (probtype_in.ne.404).and. &
          (probtype_in.ne.403)) then
  ! do not alter dt
 else
  print *,"probtype_in invalid"
  stop
 endif

enddo ! do while (finished_flag.eq.0)

end subroutine integrate_steps

end module integrate_module
