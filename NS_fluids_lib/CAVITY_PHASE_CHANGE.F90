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

! probtype==710
module CAVITY_PHASE_CHANGE_module
implicit none                   

REAL_T :: DEF_VAPOR_GAMMA
INTEGER_T, parameter :: max_sitesnum=1000
!REAL_T,parameter  :: tmax=116.0d0    ! end of the curve, currently should be 20 higher than tref 
!REAL_T,parameter  :: tinit=109.0d0
!REAL_T,parameter  :: tref=100.0d0  


INTEGER_T :: sitesnum,sitesnum2,sitesnum3
REAL_T,allocatable :: sites(:,:),sites2(:,:),sites3(:,:)  ! nucleate sites list
INTEGER_T,allocatable  :: active_flag(:),active_flag2(:),active_flag3(:)
INTEGER_T,allocatable  :: flagrecord(:)


contains

subroutine findmax_dist(flagnum,flaglist,distin,iout)
implicit none

REAL_T,intent(in) :: distin(8)
INTEGER_T,intent(in) :: flagnum
INTEGER_T,intent(in) :: flaglist(1000)
INTEGER_T            :: i,j
INTEGER_T,intent(out):: iout
REAL_T       :: tempdist

tempdist=0.0d0

do i=1,8
 if(flagnum.ne.0)then
  do j=1,flagnum
  if(distin(i).gt.tempdist.and.i.ne.flaglist(flagnum))then
   iout=i  
  endif
  enddo
 else
  if(distin(i).gt.tempdist)then
   iout=i
  endif
 endif
enddo

end subroutine findmax_dist


subroutine distcal(samp,pp,dist,vflag)
implicit none

REAL_T,intent(in) :: samp(2,8)
REAL_T,intent(in) :: pp(2)
REAL_T            :: dist(8)
INTEGER_T                 :: i
INTEGER_T,intent(out)     :: vflag
REAL_T,parameter  :: eps=0.01d0

vflag=0

do i=1,8
 dist(i)=sqrt((samp(1,i)-pp(1))**2.0+(samp(2,i)-pp(2))**2.0)
  if(dist(i).lt.eps)then
    vflag=i
  endif
enddo




end subroutine distcal

  ! do any initial preparation needed
subroutine INIT_CAVITY_PHASE_CHANGE_MODULE()
 use probcommon_module
IMPLICIT NONE
INTEGER_T :: i,j
REAL_T    :: t
REAL_T    :: r(2)
INTEGER_T :: isite
REAL_T    :: samp(2,8)
INTEGER_T :: maxi,n
REAL_T    :: a,b,fo
REAL_T    :: totdist(8),tempdist(8)
INTEGER_T :: vflag
INTEGER_T :: flagnum
REAL_T    :: tmax    ! end of the curve, currently should be 20 higher than tref 
REAL_T    :: tinit
REAL_T    :: tref  


DEF_VAPOR_GAMMA =  1.666666667D0
number_of_source_regions=1

tref=xblob4
tinit=yblob4
tmax=tref+zblob4

print *,"number_of_regions",number_of_source_regions

allocate(sites(3,max_sitesnum)) !x,y,temperature
allocate(active_flag(max_sitesnum))
allocate(flagrecord(max_sitesnum))
allocate(sites2(3,max_sitesnum)) !x,y,temperature
allocate(active_flag2(max_sitesnum))
allocate(sites3(3,max_sitesnum)) !x,y,temperature
allocate(active_flag3(max_sitesnum))

fo=5.0d0       ! function order
a=201.0/(zblob4**fo-(tinit-tref)**fo)   ! 20 superheat condition with 202 sites
b=1.0-a*((tinit-tref)**fo)            ! tinit-tref superheat condition for the first site
!print *,"a=",a,"b=",b

n=50
call random_seed(size=n)
isite=0
sites=0.0d0
t=0.0d0

do while(t.le.tmax)
 isite=isite+1
 if(isite.eq.1)then
  call random_number(r)
  if (SDIM.eq.3) then
   ! do nothing
  else if (SDIM.eq.2) then
    if(axis_dir.ne.12)then
     r(2)=zero
    else
     r(1)=zero
    endif
  else
   print *,"dimension invalid"
   stop
  endif
!  print *,r
  do i=1,2
   sites(i,isite)=r(i)
  enddo
 elseif(isite.gt.1)then
  do i=1,8

   call random_number(r)
   if (SDIM.eq.3) then
    ! do nothing
   else if (SDIM.eq.2) then
    if(axis_dir.ne.12)then
     r(2)=zero
    else
     r(1)=zero
    endif
   else
    print *,"dimension invalid"
    stop
   endif
   samp(1,i)=r(1)
   samp(2,i)=r(2)
!   print *,samp(:,i)
  enddo ! do i=1,8
!  print *,"-----------------"
  totdist=0.0
  flagrecord=0
  flagnum=0
  do i=1,isite-1
   call distcal(samp,sites(1:2,i),tempdist,vflag)
   do j=1,8
    totdist(j)=totdist(j)+tempdist(j)
   enddo
   if(vflag.ne.0)then
    flagnum=flagnum+1
    flagrecord(flagnum)=vflag
   endif
  enddo ! do i=1,isite-1
!  print *,"flagnum=",flagnum
  call findmax_dist(flagnum,flagrecord,totdist,maxi)

  do i=1,2
   sites(i,isite)=samp(i,maxi)
  enddo  
 else
  print *,"isite invalid"
  stop
 endif

 if(isite.eq.1)then
  sites(3,isite)=tinit
  t=tinit
 else
  sites(3,isite)=((real(isite,8)-b)/a)**(1.0d0/fo)+tref
  t=sites(3,isite)
 endif
! print *,"ith=",isite, "t=",t
 
enddo ! do while(t.le.tmax)
!do i=1,isite
!  sites(3,i)=sites(3,i)+273.0d0
!enddo

sitesnum=isite

!print *,"probx",problox,probhix

do j=1,sitesnum
do i=1,2
! sites(i,j)=(0.6d0-0.0d0)*sites(i,j)
! sites(i,j)=(6.0d0-0.0d0)*sites(i,j)
 if(axis_dir.eq.9)then
  sites(i,j)=(4.8d0-1.2d0)*sites(i,j)+1.2d0 
 elseif(axis_dir.eq.12)then
  if(i.eq.1)then
   sites(i,j)=(probhix-problox)*sites(i,j)
  elseif(i.eq.2)then
   sites(i,j)=0.105d0+(0.295d0-0.105d0)*sites(i,j)
  endif
 else
  sites(i,j)=(probhix-problox)*sites(i,j)
 endif
enddo
enddo

print *,"sitesnum",sitesnum,"sites list:"
do i=1,sitesnum
 print *, sites(:,i)
enddo

active_flag=0

! second group
n=50
call random_seed(size=n)
isite=0
sites2=0.0d0
t=0.0d0

do while(t.le.tmax)
 isite=isite+1
 if(isite.eq.1)then
  call random_number(r)
  if (SDIM.eq.3) then
   ! do nothing
  else if (SDIM.eq.2) then
    if(axis_dir.ne.12)then
     r(2)=zero
    else
     r(1)=zero
    endif
  else
   print *,"dimension invalid"
   stop
  endif
!  print *,r
  do i=1,2
   sites2(i,isite)=r(i)
  enddo
 elseif(isite.gt.1)then
  do i=1,8

   call random_number(r)
   if (SDIM.eq.3) then
    ! do nothing
   else if (SDIM.eq.2) then
    if(axis_dir.ne.12)then
     r(2)=zero
    else
     r(1)=zero
    endif
   else
    print *,"dimension invalid"
    stop
   endif
   samp(1,i)=r(1)
   samp(2,i)=r(2)
!   print *,samp(:,i)
  enddo ! do i=1,8
!  print *,"-----------------"
  totdist=0.0
  flagrecord=0
  flagnum=0
  do i=1,isite-1
   call distcal(samp,sites2(1:2,i),tempdist,vflag)
   do j=1,8
    totdist(j)=totdist(j)+tempdist(j)
   enddo
   if(vflag.ne.0)then
    flagnum=flagnum+1
    flagrecord(flagnum)=vflag
   endif
  enddo ! do i=1,isite-1
!  print *,"flagnum=",flagnum
  call findmax_dist(flagnum,flagrecord,totdist,maxi)

  do i=1,2
   sites2(i,isite)=samp(i,maxi)
  enddo  
 else
  print *,"isite invalid"
  stop
 endif

 if(isite.eq.1)then
  sites2(3,isite)=tinit
  t=tinit
 else
  sites2(3,isite)=((real(isite,8)-b)/a)**(1.0d0/fo)+tref
  t=sites2(3,isite)
 endif
! print *,"ith=",isite, "t=",t
 
enddo ! do while(t.le.tmax)

sitesnum2=isite


do j=1,sitesnum2
do i=1,2
 if(axis_dir.eq.9)then
  sites2(i,j)=(4.8d0-1.2d0)*sites2(i,j)+1.2d0 
 elseif(axis_dir.eq.12)then
  if(i.eq.1)then
   sites2(i,j)=(probhix-problox)*sites2(i,j)
  elseif(i.eq.2)then
   sites2(i,j)=0.1d0+(0.2d0-0.1d0)*sites2(i,j)
  endif
 else
  sites2(i,j)=(probhix-problox)*sites2(i,j)
 endif
enddo
enddo

print *,"sitesnum2",sitesnum2,"sites 2 list:"
do i=1,sitesnum2
 print *, sites2(:,i)
enddo

active_flag2=0

! thrid group
n=50
call random_seed(size=n)
isite=0
sites3=0.0d0
t=0.0d0

do while(t.le.tmax)
 isite=isite+1
 if(isite.eq.1)then
  call random_number(r)
  if (SDIM.eq.3) then
   ! do nothing
  else if (SDIM.eq.2) then
     if(axis_dir.ne.12)then
     r(2)=zero
    else
     r(1)=zero
    endif
  else
   print *,"dimension invalid"
   stop
  endif
!  print *,r
  do i=1,2
   sites3(i,isite)=r(i)
  enddo
 elseif(isite.gt.1)then
  do i=1,8

   call random_number(r)
   if (SDIM.eq.3) then
    ! do nothing
   else if (SDIM.eq.2) then
     if(axis_dir.ne.12)then
     r(2)=zero
    else
     r(1)=zero
    endif
   else
    print *,"dimension invalid"
    stop
   endif
   samp(1,i)=r(1)
   samp(2,i)=r(2)
!   print *,samp(:,i)
  enddo ! do i=1,8
!  print *,"-----------------"
  totdist=0.0
  flagrecord=0
  flagnum=0
  do i=1,isite-1
   call distcal(samp,sites3(1:2,i),tempdist,vflag)
   do j=1,8
    totdist(j)=totdist(j)+tempdist(j)
   enddo
   if(vflag.ne.0)then
    flagnum=flagnum+1
    flagrecord(flagnum)=vflag
   endif
  enddo ! do i=1,isite-1
!  print *,"flagnum=",flagnum
  call findmax_dist(flagnum,flagrecord,totdist,maxi)

  do i=1,2
   sites3(i,isite)=samp(i,maxi)
  enddo  
 else
  print *,"isite invalid"
  stop
 endif

 if(isite.eq.1)then
  sites3(3,isite)=tinit
  t=tinit
 else
  sites3(3,isite)=((real(isite,8)-b)/a)**(1.0d0/fo)+tref
  t=sites3(3,isite)
 endif
! print *,"ith=",isite, "t=",t
 
enddo ! do while(t.le.tmax)

sitesnum3=isite


do j=1,sitesnum3
do i=1,2
 if(axis_dir.eq.9)then
  sites3(i,j)=(4.8d0-1.2d0)*sites3(i,j)+1.2d0 
 elseif(axis_dir.eq.12)then
  if(i.eq.1)then
   sites3(i,j)=(probhix-problox)*sites3(i,j)
  elseif(i.eq.2)then
   sites3(i,j)=0.1d0+(0.2d0-0.1d0)*sites3(i,j)   ! z
  endif
 else
  sites3(i,j)=(probhix-problox)*sites3(i,j)
 endif
enddo
enddo

print *,"sitesnum3",sitesnum3,"sites 3 list:"
do i=1,sitesnum3
 print *, sites3(:,i)
enddo

active_flag3=0


return
end subroutine INIT_CAVITY_PHASE_CHANGE_MODULE

! this routine called from PROB.F90
subroutine Satomodel_nucleation(nucleate_in,xsten,nhalf,make_seed)
use probcommon_module_types
use probcommon_module
IMPLICIT NONE
INTEGER_T, INTENT(in) :: nhalf
REAL_T, dimension(-nhalf:nhalf,SDIM), INTENT(in) :: xsten
INTEGER_T, INTENT(inout) :: make_seed
type(nucleation_parm_type_input), INTENT(in) :: nucleate_in
INTEGER_T    :: i,j,k,im_l,im_s
INTEGER_T    :: temperature_component 
INTEGER_T    :: ii
REAL_T       :: tempt,vf_sol,ls_sol,ls_liq
REAL_T       :: tempvec1(SDIM),tempvec2(SDIM)
REAL_T       :: tempdist

!  print *,"in Satomodel_nucleation"
 if(axis_dir.eq.8.or.axis_dir.eq.9)then
  i=nucleate_in%i
  j=nucleate_in%j
  k=nucleate_in%k
! 1<=im<=num_materials
  if (num_materials.eq.3) then
   ! do nothing
  else
   print *,"num_materials invalid"
   stop
  endif
  im_l=1
  im_s=3
  temperature_component=(im_l-1)*num_state_material+ENUM_TEMPERATUREVAR+1
  tempt=nucleate_in%EOS(D_DECL(i,j,k),temperature_component)
  vf_sol=nucleate_in%Snew(D_DECL(i,j,k), &
      STATECOMP_MOF+(im_s-1)*ngeom_raw+1)
  ls_sol=nucleate_in%LSnew(D_DECL(i,j,k),im_s)
!  print *,"i=",i,"j=",j

!  print *,"ls_sol=",ls_sol,"temperature",tempt,"dx",nucleate_in%dx(SDIM)
!  print *,"xsten", xsten(-1,1),xsten(1,1)

!  ls_liq=nucleate_in%LSnew(D_DECL(i,j,k),1)
  tempvec1(1)=xsten(0,1)
  tempvec1(2)=xsten(0,2)
  tempvec1(SDIM)=xsten(0,SDIM)

    do ii=1,sitesnum
!     print *,"sitesnum=", ii
     if(tempt.ge.sites(3,ii).and.active_flag(ii).eq.0)then
!      print *,"tempt satisfied"
      if(abs(ls_sol).le.nucleate_in%dx(SDIM)+radblob3)then
          tempvec2(1)=sites(1,ii)
          tempvec2(2)=sites(2,ii)
          tempvec2(SDIM)=zblob
          call l2norm(tempvec1,tempvec2, tempdist)                               
        if(tempdist.le.radblob3)then
         print *,"make seed","temp","ls_sol","site"
         print *,"make seed",tempt,ls_sol,sites(1:2,ii),tempvec1(SDIM)
           make_seed=1
         active_flag(ii)=1
        elseif(tempdist.gt.radblob3)then
                ! do nothing
        else
          print *,"invalid tempdist",tempdist
          stop
        endif
      endif    
     endif
    enddo ! do ii=1,sitesnum

 elseif(axis_dir.eq.12)then
  i=nucleate_in%i
  j=nucleate_in%j
  k=nucleate_in%k
! 1<=im<=num_materials
  if (num_materials.eq.3) then
   ! do nothing
  else
   print *,"num_materials invalid"
   stop
  endif
  im_l=1
  im_s=3
  temperature_component=(im_l-1)*num_state_material+ENUM_TEMPERATUREVAR+1
  tempt=nucleate_in%EOS(D_DECL(i,j,k),temperature_component)
  vf_sol=nucleate_in%Snew(D_DECL(i,j,k), &
      STATECOMP_MOF+(im_s-1)*ngeom_raw+1)
  ls_sol=nucleate_in%LSnew(D_DECL(i,j,k),im_s)
 if(1.eq.0)then
  if(abs(ls_sol).le.nucleate_in%dx(SDIM)+radblob3)then
!   print *,"i=",i,"j=",j
   print *,"ls_sol=",ls_sol,"temperature",tempt,"dx",nucleate_in%dx(SDIM)
   print *,"xsten", xsten(0,1),xsten(0,2)
  endif
 endif
  ls_liq=nucleate_in%LSnew(D_DECL(i,j,k),1)
  tempvec1(1)=xsten(0,1)
  tempvec1(2)=xsten(0,2)
  tempvec1(SDIM)=xsten(0,SDIM)
    do ii=1,sitesnum
!     print *,"sitesnum=", ii
     if(tempt.ge.sites(3,ii).and.active_flag(ii).eq.0)then
!       print *,"tempt satisfied"
!       print *,"tempt=",tempt,"ls_sol=",ls_sol,"threshold",nucleate_in%dx(SDIM)+radblob3
      if(abs(ls_sol).le.nucleate_in%dx(SDIM)+radblob3)then
!       print *,"ls_sol satisfied"
!       print *,"tempt< ",sites(3,ii),"ls_sol< ",nucleate_in%dx(SDIM)," + ",radblob3
          tempvec2(1)=sites(1,ii)
          tempvec2(2)=sites(2,ii)
          tempvec2(SDIM)=zblob
          if(SDIM.eq.2)then
          tempvec2(1)=sites(2,ii)
          tempvec2(SDIM)=zblob 
          endif
          call l2norm(tempvec1,tempvec2, tempdist)
!          print *,"pointin",tempvec1
!          print *,"pointlist",tempvec2          
        if(tempdist.le.radblob3)then
         print *,"make seed ","temp ","ls_sol ","site "
         print *,"make seed",tempt,ls_sol,sites(:,ii)
           make_seed=1
         active_flag(ii)=1
        elseif(tempdist.gt.radblob3)then
                ! do nothing
        else
          print *,"invalid tempdist",tempdist
          stop
        endif
      endif    
     endif
    enddo ! do ii=1,sitesnum
    if(1.eq.1)then
    do ii=1,sitesnum2
!     print *,"sitesnum=", ii
     if(tempt.ge.sites2(3,ii).and.active_flag2(ii).eq.0)then
!      print *,"tempt satisfied"
      if(abs(ls_sol).le.nucleate_in%dx(SDIM)+radblob3)then
          tempvec2(1)=sites2(1,ii)
          tempvec2(2)=0.105d0
          tempvec2(SDIM)=sites2(2,ii)
          if(SDIM.eq.2)then
          tempvec2(1)=0.105d0
          tempvec2(SDIM)=sites2(2,ii) 
          endif
          call l2norm(tempvec1,tempvec2, tempdist)                               
        if(tempdist.le.radblob3)then
         print *,"make seed","temp","ls_sol","site"
         print *,"make seed",tempt,ls_sol,sites2(:,ii)
           make_seed=1
         active_flag2(ii)=1
        elseif(tempdist.gt.radblob3)then
                ! do nothing
        else
          print *,"invalid tempdist",tempdist
          stop
        endif
      endif    
     endif
    enddo ! do ii=1,sitesnum2
    do ii=1,sitesnum3
     if(tempt.ge.sites3(3,ii).and.active_flag3(ii).eq.0)then
      if(abs(ls_sol).le.nucleate_in%dx(SDIM)+radblob3)then
          tempvec2(1)=sites3(1,ii)
          tempvec2(2)=0.295d0
          tempvec2(SDIM)=sites3(2,ii)
          if(SDIM.eq.2)then
          tempvec2(1)=0.295d0
          tempvec2(SDIM)=sites3(2,ii) 
          endif
          call l2norm(tempvec1,tempvec2, tempdist)                               
        if(tempdist.le.radblob3)then
         print *,"make seed","temp","ls_sol","site"
         print *,"make seed",tempt,ls_sol,sites3(:,ii)
           make_seed=1
         active_flag3(ii)=1
        elseif(tempdist.gt.radblob3)then
                ! do nothing
        else
          print *,"invalid tempdist",tempdist
          stop
        endif
      endif    
     endif
    enddo ! do ii=1,sitesnum3
   endif   ! 1=0

 else
   make_seed=0 
 endif

return
end subroutine Satomodel_nucleation

subroutine CAVITY_BOILING_INIT_REGIONS_LIST(constant_density_all_time, &
      num_materials_in,num_threads_in)
use probcommon_module

IMPLICIT NONE

INTEGER_T, INTENT(in) :: num_materials_in
INTEGER_T, INTENT(in) :: num_threads_in
INTEGER_T, INTENT(in) :: constant_density_all_time(num_materials_in)
INTEGER_T :: im,iregion,dir

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
 print *,"num_threads_in",num_threads_in

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

 if (axis_dir.eq.8.or. &
     axis_dir.eq.9.or. &
     axis_dir.eq.10.or.&
     axis_dir.eq.12) then
  regions_list(1,0)%region_material_id=3 !heater
  regions_list(1,0)%region_energy_flux=xblob3 ! Watts=J/s
!  print *,"m_id",regions_list(1,0)%region_material_id
!  print *,"flux",regions_list(1,0)%region_energy_flux
 elseif(axis_dir.eq.11)then
  regions_list(1,0)%region_material_id=2 !heater
  regions_list(1,0)%region_volume_flux=xblob6
  regions_list(1,0)%region_mass_flux=xblob6*fort_denconst(2)
  regions_list(1,0)%region_temperature_prescribe=yblob6
 else
  print *,"axis_dir must equal to 8 or 9or10"
  stop
 endif
end subroutine CAVITY_BOILING_INIT_REGIONS_LIST

subroutine CAVITY_BOILING_CHARFN_REGION(region_id,x,cur_time,charfn_out)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: region_id
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: cur_time
REAL_T, INTENT(out) :: charfn_out

REAL_T           :: rtemp

  if(axis_dir.eq.8)then 
   if(x(SDIM).lt. xblob2.and.x(SDIM).gt.zblob2)then
    charfn_out=one
   else
    charfn_out=0.0d0
   endif
  elseif(axis_dir.eq.9)then
   if(x(SDIM).lt. xblob2.and.x(SDIM).gt.zblob2 .and.&
      x(1).le.5.0d0 .and. x(1).ge.1.0d0 .and. &     
      x(2).le.5.0d0 .and. x(2).ge.1.0d0)then
    charfn_out=one
   else
    charfn_out=0.0d0
   endif
  elseif(axis_dir.eq.10)then
   if(x(SDIM).lt. xblob2+1.0e-3.and.x(SDIM).gt.xblob2 .and.&
      x(1).le.5.0e-3 .and. x(1).ge.1.0e-3 .and. &     
      x(2).le.5.0e-3 .and. x(2).ge.1.0e-3)then
    charfn_out=one
   else
    charfn_out=0.0d0
   endif
  elseif(axis_dir.eq.11)then
   if(x(SDIM).le.probloz+0.000375d0 .and.&
      x(1).le.(probhix-problox)*0.5d0+problox+0.000375d0 .and. &
      x(1).ge.(probhix-problox)*0.5d0+problox)then      
    charfn_out=one
   else
    charfn_out=0.0d0
   endif
  elseif(axis_dir.eq.12)then
      ! circumfrence = 0.2  
    if(SDIM.eq.3)then
      rtemp=0.2d0/(2.0*4.0*atan(1.0d0))-&
              sqrt((x(2)-0.2d0)**2.0d0+(x(SDIM)-0.15d0)**2.0d0)
    elseif(SDIM.eq.2)then
      rtemp=0.2d0/(2.0*4.0*atan(1.0d0))-&
              sqrt((x(1)-0.2d0)**2.0d0+(x(SDIM)-0.15d0)**2.0d0)
    else
      print *,"SDIM INVALID"
      stop
    endif

   if(rtemp.ge.0.0d0)then
    charfn_out=one
   else
    charfn_out=0.0d0
   endif
  else
   print *,"axis_dir must equal to 8-11"
   stop
  endif

end subroutine CAVITY_BOILING_CHARFN_REGION

subroutine CAVITY_BOILING_THERMAL_K(x,dx,cur_time, &
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

INTEGER_T, INTENT(in) :: im
INTEGER_T, INTENT(in) :: im_solid
INTEGER_T, INTENT(in) :: near_interface
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: dx(SDIM)
REAL_T, INTENT(in) :: cur_time
REAL_T, INTENT(in) :: density
REAL_T, INTENT(in) :: temperature
REAL_T, INTENT(in) :: temperature_wall
REAL_T, INTENT(in) :: temperature_wall_max
REAL_T, INTENT(in) :: temperature_probe
REAL_T, INTENT(in) :: nrm(SDIM) ! nrm points from solid to fluid
REAL_T, INTENT(inout) :: thermal_k
INTEGER_T :: local1
REAL_T   :: rtemp

if(axis_dir.eq.8)then
  if(x(SDIM).lt.zblob2)then
   thermal_k=0.0d0
  endif
elseif(axis_dir.eq.9)then
  local1=0   
   if(x(SDIM).le.zblob.and.x(SDIM).ge.zblob2+0.2d0 .and. &
     x(1).le.5.0d0-0.2d0 .and. x(1).ge.1.0d0+0.2d0 .and. &     
     x(2).le.5.0d0-0.2d0 .and. x(2).ge.1.0d0+0.2d0)then
     local1=1
   else
     ! do nothing
   endif  

   if(x(SDIM).le.zblob.and.x(SDIM).ge.zblob2 .and.&
     x(1).le.5.0d0 .and. x(1).ge.1.0d0 .and. &     
     x(2).le.5.0d0 .and. x(2).ge.1.0d0 .and. &
     local1.eq.0) then
     thermal_k=0.0d0
   else
     ! do nothing
   endif  
elseif(axis_dir.eq.10)then
  local1=0   
   if(x(SDIM).le.zblob-0.2e-3 .and. x(SDIM).ge.zblob-1.0e-3 .and. &
     x(1).le.5.0e-3-0.2e-3 .and. x(1).ge.1.0e-3+0.2e-3 .and. &     
     x(2).le.5.0e-3-0.2e-3 .and. x(2).ge.1.0e-3+0.2e-3)then
     local1=1
   else
     ! do nothing
   endif  

   if(x(SDIM).le.zblob.and.x(SDIM).ge.zblob-1.0e-3 .and.&
     x(1).le.5.0e-3 .and. x(1).ge.1.0e-3 .and. &     
     x(2).le.5.0e-3 .and. x(2).ge.1.0e-3 .and. &
     local1.eq.0) then
     thermal_k=0.0d0
   else
     ! do nothing
   endif 
elseif(axis_dir.eq.12)then
    if(SDIM.eq.3)then
      rtemp=0.2d0/(2.0*4.0*atan(1.0d0))-& 
           sqrt((x(2)-0.2d0)**2.0d0+(x(SDIM)-0.15d0)**2.0d0)
    elseif(SDIM.eq.2)then
      rtemp=0.2d0/(2.0*4.0*atan(1.0d0))-&
           sqrt((x(1)-0.2d0)**2.0d0+(x(SDIM)-0.15d0)**2.0d0)
    else
      print *,"SDIM INVALID"
      stop
    endif
   if(rtemp.ge.0.0d0)then
    thermal_k=1.0e+8
   else
    ! do nothing
   endif
else
  ! do nothing
endif

end subroutine CAVITY_BOILING_THERMAL_K


! ----------
! YANG LIU
! ----------
!---------------------------------------------
subroutine find_p_plane_dist_wb(flag,nv,vts,p,dist)
! check a point if its projection of a plane fall inside the 
! close domain surrounded by the boudary 
! if it is, then dist = find_p_phase_dist
! if not, find the closest point of on the boundary 
implicit none

INTEGER_T,INTENT(in)  :: nv ! number of vertices 
REAL_T,INTENT(in) :: vts(nv,SDIM)                          
     !vertices corrd in clockwise or counterclockwise
REAL_T,INTENT(in) :: p(SDIM)
REAL_T        :: pp(SDIM)
INTEGER_T             :: i,j
INTEGER_T,INTENT(in)  :: flag
! if flag = +1, then counter clockwise
! if flag = -1, then clockwise
! the normal of the plane always point to positive material.
REAL_T         :: sign1,sign2,sign3
REAL_T         :: s1(SDIM),s2(SDIM)
INTEGER_T              :: isd_flag
REAL_T         :: dist
REAL_T         :: dist1(nv)
REAL_T         :: normal(SDIM)
REAL_T         :: d
REAL_T         :: x1temp(SDIM)
REAL_T         :: x2temp(SDIM)
REAL_T         :: x3temp(SDIM)
INTEGER_T      :: dir

if(SDIM .ne. 3)then
 print *,"invalid dimension"
 stop
endif
if ((nv.ne.3).and.(nv.ne.4)) then
 print *,"nv invalid"
 stop
endif

isd_flag = 1 ! inside the domain

do dir=1,SDIM
 x1temp(dir)=vts(1,dir)
 x2temp(dir)=vts(2,dir)
 x3temp(dir)=vts(3,dir)
enddo

! use the first three point in vts to make a plane
call make_3points_plane(x1temp,x2temp,x3temp,normal,d)


! find the projection point pp
call find_pp_plane(p,normal,d,pp)

! determine if the pp is inside the vts(nv)
 do i = 1, nv-1
  do j = 1,SDIM
   if(flag .eq. +1) then  ! vts in counter clockwise direction
    s1(j) = vts(i+1,j) - vts(i,j)
    s2(j) = pp(j) -vts(i,j)   
   elseif(flag .eq. -1)then
    s2(j) = vts(i+1,j) - vts(i,j)
    s1(j) = pp(j) -vts(i,j)
   else
    print *,"invalid flag"
    stop
   endif
  enddo
  sign1 = s1(2)*s2(SDIM) - s1(SDIM)*s2(2)
  sign2 = -s1(1)*s2(SDIM) + s1(SDIM)*s2(1)
  sign3 = s1(1)*s2(2) - s1(2)*s2(1)
  if(sign1*normal(1) .lt. 0.0d0 .or. &
     sign2*normal(2) .lt. 0.0d0 .or. &
     sign3*normal(SDIM) .lt. 0.0d0)then
   isd_flag = 0
   exit
  endif
 enddo
 if(isd_flag .eq. 1)then
  do j = 1,SDIM
   if(flag .eq. 1)then
    s1(j) = vts(1,j) - vts(nv,j)
    s2(j) = pp(j) -vts(nv,j)
   elseif(flag .eq. -1)then
    s2(j) = vts(1,j) - vts(nv,j)
    s1(j) = pp(j) -vts(nv,j)
   else
    print *,"invalid flag"
    stop
   endif
  enddo
  sign1 = s1(2)*s2(SDIM) - s1(SDIM)*s2(2)
  sign2 = -s1(1)*s2(SDIM) + s1(SDIM)*s2(1)
  sign3 = s1(1)*s2(2) - s1(2)*s2(1)
  if(sign1*normal(1) .lt. 0.0d0 .or. &
     sign2*normal(2) .lt. 0.0d0 .or. &
     sign3*normal(SDIM) .lt. 0.0d0)then
   isd_flag = 0
  endif
 endif


if(isd_flag .eq. 1)then
 call find_p_plane_dist(p,normal,d,dist) 
elseif(isd_flag .eq. 0)then
 do i = 1,nv-1
  do dir=1,SDIM
   x1temp(dir)=vts(i,dir)
   x2temp(dir)=vts(i+1,dir)
  enddo
  call dist_point_to_line(x1temp,x2temp,p,dist1(i))
 enddo
 do dir=1,SDIM
  x1temp(dir)=vts(nv,dir)
  x2temp(dir)=vts(1,dir)
 enddo
 call dist_point_to_line(x1temp,x2temp,p,dist1(nv))
 dist = minval(abs(dist1))
else
 print *,"invalid isd_flag"
 stop
endif


end subroutine find_p_plane_dist_wb
!----------------------------------------------
! find the shortest distance from a point to plane(normal,d)
subroutine find_p_plane_dist(p,normal,d,dist)
implicit none

REAL_T, INTENT(in) :: p(SDIM)
REAL_T, INTENT(in) :: normal(SDIM)
REAL_T, INTENT(in) :: d
REAL_T,INTENT(out)  :: dist
REAL_T             :: nn

if(SDIM .ne. 3)then
 print *,"invalid dimension"
 stop
endif

dist = 0.0d0
nn = sqrt(normal(1)**2.0d0 + normal(2)**2.0d0 &
          + normal(SDIM)**2.0d0)
if(nn .le. 0.0d0)then
 print *,"warning, the denomintor is zero"
 stop
endif

dist = abs(p(1)*normal(1) + p(2)*normal(2) + &
          p(SDIM)*normal(SDIM) + d)/nn

end subroutine find_p_plane_dist

!----------------------------------------------
subroutine plane_distf(p,normal,d,val)
implicit none

! plane level set function 
REAL_T, INTENT(in) :: p(SDIM)
REAL_T, INTENT(in) :: normal(SDIM)
REAL_T, INTENT(in) :: d

REAL_T             :: val

if(SDIM .ne. 3)then
 print *,"invalid dimension"
 stop
endif

VAL=normal(1)*p(1)+normal(2)*p(2)+normal(SDIM)*p(SDIM)+d

end subroutine plane_distf
!---------------------------------------------
! find the projection pp of a point p on a given plane (normal,d)
subroutine find_pp_plane(p,normal,d,pp)
implicit none

REAL_T,INTENT(in)   :: p(SDIM),normal(SDIM)
REAL_T,INTENT(in)   :: d
REAL_T              :: pp(SDIM)

INTEGER_T                   :: i
REAL_T              :: s

if(SDIM .ne. 3)then
 print *,"invalid dimension"
 stop
endif
if(normal(1) .eq. 0.0d0 .and. normal(2) .eq. 0.0d0 & 
   .and. normal(SDIM) .eq. 0.0d0 )then
 print *,"warning, 0 normal"
 stop
endif

s = (-d - normal(1)*pp(1) - normal(2)*pp(2) - &
    normal(SDIM)*pp(SDIM))/ &
    (normal(1)**2.0d0 + normal(2)**2.0d0 + &
     normal(SDIM)**2.0d0)

do i = 1,SDIM
 pp(i) = s*normal(i) + pp(i)
enddo

end subroutine find_pp_plane
!-----------------------------------------------
subroutine adjust_plane_sign(flag,p,normal,d)
implicit none
! for test point p
! if \phi(p) >= 0,  flag = 1
! if \phi(p) < 0,   flag = -1
INTEGER_T    ,INTENT(in)   :: flag
REAL_T, INTENT(in) :: p(SDIM)
REAL_T             :: normal(SDIM)
REAL_T             :: d

INTEGER_T                  :: i
REAL_T             :: val

if(SDIM .ne. 3)then
 print *,"invalid dimension"
 stop
endif

if(flag .ne. 1 .and. flag .ne. -1)then
 print *,"invalid test flag in adjust_plane_sign"
 stop
endif

VAL=normal(1)*p(1)+normal(2)*p(2)+normal(SDIM)*p(SDIM)+d

if(flag .eq. 1)then
 if(val .ge. 0.0d0)then
  ! do nothing
 else
  do i = 1,SDIM
   normal(i) = -1.0d0*normal(i)
   d = -1.0d0*d
  enddo
 endif
elseif(flag .eq. -1)then
 if(val .lt. 0.0d0) then
  ! do nothing
 else
  do i = 1,SDIM
   normal(i) = -1.0d0*normal(i)
   d = -1.0d0*d
  enddo
 endif
endif

end subroutine
!---------------------------------------------------------
subroutine make_3points_plane(p1,p2,p3,normal,d)
implicit none
! find a plane (normal,d) with 3 points p1,p2,p3

REAL_T, INTENT(in)     :: p1(SDIM),P2(SDIM),p3(SDIM)

REAL_T                 :: normal(SDIM)
REAL_T                 :: d

if(SDIM .ne. 3)then
 print *,"invalid call make_3points_plane"
 stop
endif

call crossproduct(p1,p2,p3,normal)
d = -normal(1)*p1(1) - normal(2)*p1(2) - normal(SDIM)*p1(SDIM)


end subroutine make_3points_plane
!------------------------------------------------------------
!       p1-------------------------p2
!               *
!               *  
!               x
!----------------------------------------------------------
subroutine dist_point_to_line(p1,p2,x,dist)
implicit none
! represent the line in parametric form,(v = f(s))
! v^x = x1 + (x2-x1)s
! v^y = y1 + (y2-y1)s
! v^z = z1 + (z2-z1)s
! 
! if the closest point on the line to point x  is outside p1 -- p2, 
! --------> return the distance from x either to p1 or p2 which is shorter.
! otherwise
! --------> return the distance from x to the cloest point

REAL_T,INTENT(in)      ::  p1(SDIM),p2(SDIM),x(SDIM)
REAL_T                 ::  dist,dist1,dist2

REAL_T                 :: diff10,diff21
REAL_T                 :: x10(SDIM), x21(SDIM)
INTEGER_T              :: i
REAL_T                 :: s


dist = 0.0d0

do i = 1,SDIM
 x10(i) = p1(i) - x(i)
 x21(i) = p2(i) - p1(i)

 if(abs(x10(i)) .lt. 10e-10)then
  x10(i) = 0.0d0
 endif
 if(abs(x21(i)) .lt. 10e-10)then
  x21(i) = 0.0d0
 endif

enddo

if (maxval(abs(x21)) .lt. 10d-8)then
 print *,"p1 and p2 are coincide with each other",p1,p2
 stop
endif
if (maxval(abs(x10)) .lt. 10d-8)then
 dist=zero
else

 call l2norm(p1, x, diff10)
 call l2norm(p2, p1, diff21)

 s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)

 if(s .le. 1.0d0 .and. s .ge. 0.0)then

  do i = 1,SDIM
   dist = dist + (x10(i) + x21(i)*s)**2.0d0
  enddo
  dist = sqrt(dist)
 else
  call l2norm(p1,x,dist1)
  call l2norm(p2,x,dist2)
  dist = min(dist1,dist2)
 endif

endif

end subroutine dist_point_to_line
!---------------------------------------------------------
subroutine crossproduct(p1,p2,p3,oter)
implicit none

REAL_T,INTENT(in)       :: p1(SDIM),p2(SDIM),p3(SDIM)
REAL_T,INTENT(out)      :: oter(SDIM)

INTEGER_T                      :: i
REAL_T                 :: p1p2(SDIM),p1p3(SDIM)

if(SDIM .ne. 3)then
 print *,"invalid call outerproduct"
 stop
endif

oter = 0.0d0

do i = 1,SDIM
 p1p2(i)=p2(i)-p1(i)
 p1p3(i)=p3(i)-p1(i)
enddo 

oter(1) = p1p2(2)*p1p3(SDIM) - p1p2(SDIM)*p1p3(2) 
oter(2) = -p1p2(1)*p1p3(SDIM)+ p1p2(SDIM)*p1p3(1)
oter(SDIM) = p1p2(1)*p1p3(2) - p1p2(2)*p1p3(1)

end subroutine crossproduct
!------------------------------------------------------------
subroutine l2norm(x1,x2, x1x2norm)
implicit none

REAL_T,INTENT(in)  :: x1(SDIM),x2(SDIM)
REAL_T             :: x1x2norm

INTEGER_T                  :: i
REAL_T,allocatable :: diff(:)

x1x2norm = 0.0d0
allocate(diff(SDIM))
do i = 1,SDIM
 diff(i) = x1(i)-x2(i)
enddo

do i = 1,SDIM
 x1x2norm = x1x2norm + diff(i)**2.0d0
enddo

x1x2norm = sqrt(x1x2norm)

deallocate(diff)

end subroutine l2norm

! negative outside the object, positive inside.
subroutine melting_ls(x1,x2,y1,y2, x_in,dist)
! 2d distant function
!             :------                    y1
!             :      |    
!             :      |
!             :      |
!             :      |----------|        y2
!             :                 | 
!             :                 |
!             :-----------------|                      
!  
!           (0,-y1)     x1       x2



! x_in(2):  the input x
! dist:  output signed distance

implicit none

REAL_T,INTENT(in)  :: x1,x2,y1,y2
! x1 x2 y1 y2 shape parameter
REAL_T,INTENT(in)  ::  x_in(SDIM)
REAL_T,INTENT(out) :: dist

REAL_T             :: x(SDIM)
REAL_T             :: dist1,dist2,dist3
REAL_T             :: dist4,dist5
REAL_T             :: p1(SDIM),p2(SDIM),p3(SDIM)
REAL_T             :: p4(SDIM),p5(SDIM),p6(SDIM)
REAL_T             :: ss    ! sign
INTEGER_T                  :: i

p1(1)=0.0d0
p1(2)=y1
p2(1)=x1
p2(2)=y1
p3(1)=x1
p3(2)=y2
p4(1)=x2
p4(2)=y2
p5(1)=x2
p5(2)=-y1
p6(1)=0.0d0
p6(2)=-y1

if(x_in(1) .ge. 0.0d0)then
  do i=1,2
   x(i) = x_in(i)
  enddo
else
  x(1)=-1.0d0*x_in(1)
  x(2)=x_in(2)
endif

! calculate the abs distance
  call dist_point_to_lined(p1,p2,x,dist1)
  call dist_point_to_lined(p2,p3,x,dist2)  
  call dist_point_to_lined(p3,p4,x,dist3)  
  call dist_point_to_lined(p4,p5,x,dist4)  
  call dist_point_to_lined(p5,p6,x,dist5)  
! find the sign
  if( x(1) .ge. 0.0d0 .and. x(1) .le. x2  .and. &
      x(2) .ge. -y1 .and. x(2) .le. y2)then
     ss = 1.0d0
  elseif( x(1) .ge. 0.0d0 .and. x(1) .le. x1  .and. &
      x(2) .ge. y2 .and. x(2) .le. y1)then
     ss = 1.0d0
  else
    ss = -1.0d0
  endif
 
  dist = ss*min(dist1,dist2,dist3,dist4,dist5)


end subroutine melting_ls

subroutine dist_point_to_lined(p1,p2,x,dist)
implicit none
! represent the line in parametric form,(v = f(s))
! v^x = x1 + (x2-x1)s
! v^y = y1 + (y2-y1)s
! v^z = z1 + (z2-z1)s
! 
! if the closest point on the line to point x  is outside p1 -- p2, 
! --------> return the distance from x either to p1 or p2 which is shorter.
! otherwise
! --------> return the distance from x to the cloest point

REAL_T,INTENT(in)      ::  p1(SDIM),p2(SDIM),x(SDIM)
REAL_T                 ::  dist


REAL_T                 :: diff10,diff21
REAL_T,allocatable     :: x10(:), x21(:)
INTEGER_T                      :: i
REAL_T                 :: s


dist = 0.0d0

allocate(x10(SDIM),x21(SDIM))
do i = 1,SDIM
 x10(i) = p1(i) - x(i)
 x21(i) = p2(i) - p1(i)
enddo

if (maxval(abs(x21)) .lt. 10d-8)then
 print *,"p1 and p2 are coincide with each other",p1,p2
 stop
endif

call l2norm(p1, x, diff10)
call l2norm(p2, p1, diff21)


s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)


if(s .gt. 1.0d0)then
 call l2norm(p2, x,dist)
elseif(s .lt. 0.0d0)then
 call l2norm(p1,x,dist)
else
 if(abs((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0) .lt. 1.0e-10)then
   dist= 0.0d0
 else
  dist = sqrt(((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0)/ &
        (diff21**2.0d0))
 endif
endif

deallocate(x10,x21) 

end subroutine dist_point_to_lined

subroutine dist_cuboid(dif1,dif2,dif3,dif4,dif5,dif6,dist)
! dif1-6 is the difference in each direction 
! dif1=x(1)-lox
! dif2=x(1)-hix
! dif3=x(2)-loy
! dif4=x(2)-hiy
! dif5=x(SDIM)-loz
! dif6=x(SDIM)-hiz
implicit none

REAL_T,intent(in) :: dif1,dif2,dif3,dif4,dif5,dif6
REAL_T            :: dist1,dist2,dist3,dist4,dist5,dist6
REAL_T            :: d1,d2,d3
REAL_T,intent(out):: dist
  dist1=dif1
  dist2=dif2
  dist3=dif3
  dist4=dif4
  dist5=dif5
  dist6=dif6

  d1=dist1*dist2
  d2=dist3*dist4
  d3=dist5*dist6

  if(d1.le.0.0d0)then
    if(d2.le.0.0d0)then  ! d1<0 d2<0 any d3
      if(d3.le.0.0d0)then ! d1<0 d2<0 d3<0
       dist=min(abs(dist5),abs(dist6),abs(dist1),abs(dist2),abs(dist3),abs(dist4))
      else   ! d1<0 d2<0 d3>0
       dist=min(abs(dist5),abs(dist6))
      endif
    else
      if(d3.le.0.0d0)then ! d1<0 d2>0 d3<0
       dist=min(abs(dist3),abs(dist4))
      else  ! d1<0 d2>0 d3>0
        dist=sqrt((min(abs(dist3),abs(dist4)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
      endif
    endif
  else
    if(d2.le.0.0d0)then
     if(d3.le.0.0d0)then  ! d1>0 d2<0 d3<0
       dist=min(abs(dist1),abs(dist2))
     else  ! d1>0 d2<0 d3>0
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
     endif
    else     
     if(d3.le.0.0d0)then  ! d1>0 d2>0 d3<0
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist3),abs(dist4)))**2.0d0)
     else   ! d1>0 d2>0 d3>0
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                (min(abs(dist3),abs(dist4)))**2.0d0+ &
                (min(abs(dist5),abs(dist6)))**2.0d0)
     endif  
    endif
  endif

end subroutine dist_cuboid

subroutine dist_rectangular(dif1,dif2,dif3,dif4,dist)
! dif1=x(1)-lox
! dif2=x(1)-hix
! dif3=x(2)-loy
! dif4=x(2)-hiy
implicit none

REAL_T,intent(in) :: dif1,dif2,dif3,dif4
REAL_T            :: dist1,dist2,dist3,dist4
REAL_T            :: d1,d2
REAL_T,intent(out):: dist

  dist1=dif1
  dist2=dif2
  dist3=dif3
  dist4=dif4

  d1=dist1*dist2
  d2=dist3*dist4

  if(d1.le.0.0d0)then
    if(d2.le.0.0d0)then  ! d1<0 d2<0 any d3
     dist=min(abs(dist1),abs(dist2),abs(dist3),abs(dist4))
    else
     dist=min(abs(dist3),abs(dist4)) 
    endif
  else
    if(d2.le.0.0d0)then
     dist=min(abs(dist1),abs(dist2)) 
    else     
     dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist3),abs(dist4)))**2.0d0)
    endif
  endif

end subroutine dist_rectangular

! ----------
! YANG LIU
! ----------
! From paper:  
!  Mu YT, Chen L, He YL, Kang QJ, Tao WQ. 
! Nucleate boiling performance evaluation of cavities at mesoscale level. 
! International Journal of Heat and Mass Transfer. 2017 Mar 31;106:708-19.
! 
! internal length unit ( mm )
! external length unit ( m )
! 
!  MATERIAL NUMBER
!  1:liquid   2: vapor   3: solid 

!  COORDINATE TYPE
!  1 3D cartisian
!  2 R-Z   

!  SHAPE 
!  1--> rectangular
!  2--> trapzoidal
!  3--> vshape
!  4--> semicircle
!  5--> spherical reentrant

subroutine cavity_distf_13(cavity_type, coord_type, x_in, dist)
use probcommon_module
implicit none

! liquid +  solid - 

INTEGER_T,INTENT(in)      :: coord_type
INTEGER_T, INTENT(in)     :: cavity_type
REAL_T,INTENT(in) :: x_in(SDIM)
REAL_T            :: x(SDIM)

INTEGER_T         :: dist_sign
REAL_T, INTENT(out) :: dist
REAL_T            :: dist1,dist2,dist3,dist4,dist5,dist6
REAL_T            :: x1(SDIM),x2(SDIM),x3(SDIM)
REAL_T            :: x4(SDIM),x5(SDIM)
REAL_T            :: trap,vshape,radius,lcinter,radius3d,bb
INTEGER_T i

REAL_T            :: pl1(SDIM),pl2(SDIM)
REAL_T            :: d1,d2,d3
REAL_T            :: val1,val2 
REAL_T            :: tp(SDIM)
REAL_T            :: v1(SDIM),v2(SDIM),v3(SDIM)
REAL_T            :: plg1(4,SDIM),plg2(4,SDIM)
REAL_T            :: plg3(3,SDIM),plg4(3,SDIM)
REAL_T            :: x1temp(SDIM)
REAL_T            :: x2temp(SDIM)
REAL_T            :: x3temp(SDIM)
REAL_T            :: ztest

INTEGER_T dir

INTEGER_T,parameter :: debugflag = 0
REAL_T,   parameter :: scale_factor=1.0d0
dist_sign = 1.0d0
dist = 0.0d0
dist1 = 0.0d0
dist2 = 0.0d0
dist3 = 0.0d0
do dir=1,SDIM
 x1(dir) = 0.0d0
 x2(dir) = 0.0d0
 x3(dir) = 0.0d0
 x4(dir) = 0.0d0
enddo
 
! convert from m to mm 
! April, 2018: convert from m to cm 
do i = 1,SDIM
  x(i) = x_in(i)*scale_factor
enddo
if (1.eq.0) then
 do i = 1,SDIM
  x(i) = x_in(i)
 enddo
endif

 if(debugflag .eq. 1)then
  print *,"case", cavity_type
  print *,"coord_type", coord_type
 endif


!  dimension sanity check
if(coord_type .eq. 1.and.cavity_type.ne.8.and.cavity_type.ne.12) then   ! 3D cartisian
 if(SDIM .ne. 3)then
  print *, "coord_type is not consistent with dimension"
  stop
 endif
elseif(coord_type.eq.2.and.cavity_type.ne.8.and. &
                cavity_type.ne.12) then ! R-Z axis symmetric
 if(SDIM .ne. 2) then
  print *,"coord_type is not consistent with dimension"
  stop
 endif
endif



if(cavity_type .eq. 1)then       ! 1--> rectangular
!  *                 * 
!  *       *********** 
!  *       *         * 
!  *       *         * 
!  *********         * 
!  
!*******************************************************************

 
 if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 if(x(1) .ge. 0.25d0) then   ! r >= 0.25
  dist1 = x(2) - 0.3d0        ! dist = z - 0.3   
  if(x(2) .le. 0.05d0)then
   dist2 = sqrt( (x(1)- 0.25d0)**2.0d0 + (x(2) - 0.05d0)**2.0d0)
   dist = -1.0d0*min(abs(dist1), abs(dist2))
  elseif(x(2) .lt. 0.3d0) then
   dist2 = 0.25d0 - x(1)
   dist = -1.0d0*min(abs(dist1),abs(dist2))
  else
   dist = dist1
  endif
 else  !   0 < r < 0.25
  dist1 = x(2) - 0.05
  if(x(2) .gt. 0.3d0)then
   dist2 = sqrt( (x(1)- 0.25d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0)
   dist = min(dist1, dist2)
  elseif(x(2) .gt. 0.05d0)then
   dist2 = 0.25 - x(1) 
   dist = min(dist1,dist2)
  else
   dist = dist1
  endif
 endif
 endif


 if(coord_type .eq. 1)then   ! 3D cartisian
 if(x(1) .ge. 0.25d0 .or. x(2) .ge. 0.25d0) then
  dist1 = x(SDIM) - 0.3d0 
  if(x(SDIM) .gt. 0.3d0)then         ! x(SDIM) > 0.3
   dist = dist1
  elseif(x(SDIM) .lt. 0.05d0) then   !    x(SDIM) < 0.05
    if(x(1) .ge. 0.25d0 .and. x(2) .le. 0.25d0)then
    x1(1) = 0.25d0
    x1(2) = 0.0d0
    x1(SDIM) = 0.05d0
    x2(1) = 0.25d0
    x2(2) = 0.25d0
    x2(SDIM) = 0.05d0
    call dist_point_to_line(x1,x2,x,dist2)
    dist = -1.0d0*min(abs(dist1),abs(dist2))
   elseif(x(1) .le. 0.25d0 .and. x(2) .ge. 0.25d0)then
    x1(1) = 0.0d0
    x1(2) = 0.25d0
    x1(SDIM) = 0.05d0
    x2(1) = 0.25d0
    x2(2) = 0.25d0
    x2(SDIM) = 0.05d0
    call dist_point_to_line(x1,x2,x,dist2)
    dist = -1.0d0*min(abs(dist1),abs(dist2))
   else
    x1(1) = 0.25d0
    x1(2) = 0.25d0
    x1(SDIM) = 0.05d0
    call l2norm(x,x1,dist2)
    dist = -1.0d0*min(abs(dist1),abs(dist2))    
   endif
  else   !  0.05 <  x(SDIM)  < 0.3
   if(x(1) .ge. 0.25d0 .and. x(2) .le. 0.25d0)then
    dist2 = 0.25d0 - x(1)
    dist = -1.0d0*min(abs(dist1),abs(dist2))
   elseif(x(1) .le. 0.25d0 .and. x(2) .ge. 0.25d0)then
    dist2 = 0.25d0 - x(2) 
    dist = -1.0d0*min(abs(dist1),abs(dist2))
   else
    x1(1) = 0.25d0
    x1(2) = 0.25d0
    x1(SDIM) = 0.3d0
    x2(1) = 0.25d0
    x2(2) = 0.25d0
    x2(SDIM) = 0.05d0
    call dist_point_to_line(x1,x2,x,dist2)
    dist = -1.0d0*min(abs(dist1),abs(dist2))    
   endif   
  endif
 else    ! x(1) .lt. 0.25d0 .and. x(2) .lt. 0.25d0
  dist1 = x(SDIM) - 0.05d0
  if(x(SDIM) .lt. 0.05d0 )then
   dist = dist1
  elseif(x(SDIM) .gt. 0.3d0)then
    x1(1) = 0.25d0
    x1(2) = 0.0d0
    x1(SDIM) = 0.3d0
    x2(1) = 0.25d0
    x2(2) = 0.25d0
    x2(SDIM) = 0.3d0
    call dist_point_to_line(x1,x2,x,dist2) 
    x3(1) = 0.0d0
    x3(2) = 0.25d0
    x3(SDIM) = 0.3d0
    x4(1) = 0.25d0
    x4(2) = 0.25d0
    x4(SDIM) = 0.3d0
    call dist_point_to_line(x3,x4,x,dist3) 
    if(dist2 .lt. 0.0d0 .or. dist3 .lt. 0.0d0) then
     print *,"dist2 and dist3 should be positive"
     stop 
    endif  
    dist = min(abs(dist1),abs(dist2),abs(dist3))
  else  !   0.05 <= x(SDIM) <= 0.3
    dist2 = 0.25d0 - x(1)   
    dist3 = 0.25d0 - x(2)
    if(dist2 .lt. 0.0d0 .or. dist3 .lt. 0.0d0) then
     print *,"dist2 and dist3 should be positive"
     stop 
    endif
    dist = min(abs(dist1),abs(dist2),abs(dist3))  
  endif
 endif
 endif   
!---------------------------------------------------------------



elseif(cavity_type .eq. 2)then    ! 2--> trapzoidal
!**************************************************************
!             *                 *  
!             *         x3      x4
!             *         *********
!             *       *
!             *     *
!             *****
!            x1   x2
!**************************************************************

if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 ! determine sign
 if(x(1) .ge. 0.25d0)then
  if(x(2) .ge. 0.3d0)then
   dist_sign = 1.0d0
  else
   dist_sign = -1.0d0
  endif 
 elseif(x(1) .le. 0.125d0)then
  if(x(2) .ge. 0.05d0)then
   dist_sign = 1.0d0
  else
   dist_sign = -1.0d0
  endif
 else
   trap = 2.0d0*(x(1)-0.125d0) + 0.05d0 - x(2)
   if(trap .gt. 0.0d0)then
    dist_sign = -1.0d0
   else
    dist_sign = +1.0d0
   endif
 endif
 ! determine distance
  x1(1) = 0.0d0
  x1(2) = 0.05d0 
  x2(1) = 0.125d0
  x2(2) = 0.05d0
  x3(1) = 0.25d0
  x3(2) = 0.3d0
  x4(1) = 0.6d0
  x4(2) = 0.3d0
  call dist_point_to_line(x1,x2,x,dist1)
  call dist_point_to_line(x2,x3,x,dist2) 
  call dist_point_to_line(x3,x4,x,dist3)  
  dist = dist_sign * min(dist1,dist2,dist3)

endif

if(coord_type .eq. 1)then   !  3D cartisian
 ! determine sign
 ! slant plane 1: x1 x2 x3
 ! slant plane 2: x2 x4 x5  
 x1(1) = 0.25d0 
 x1(2) = 0.0d0
 x1(SDIM) = 0.3d0
 x2(1) = 0.25d0
 x2(2) = 0.25d0
 x2(SDIM) = 0.3d0
 x3(1) = 0.125d0 
 x3(2) = 0.0d0
 x3(SDIM) = 0.05d0
 x4(1) = 0.0d0
 x4(2) = 0.25d0
 x4(SDIM) = 0.3d0
 x5(1) = 0.125d0
 x5(2) = 0.125d0
 x5(SDIM) = 0.05d0
 call make_3points_plane(x1,x2,x3,pl1,d1)
 call make_3points_plane(x2,x4,x5,pl2,d2)
 if(1 .eq. 0)then
 tp(1) = 0.0d0
 tp(2) = 0.0d0
 tp(SDIM) = 0.3d0
 call adjust_plane_sign(+1,tp,pl1,d1)
 call adjust_plane_sign(+1,tp,pl2,d2)
 endif
 call plane_distf(x,pl1,d1,val1)    ! val1 and val2 could be + or -
 call plane_distf(x,pl2,d2,val2)

 if(x(SDIM) .ge. 0.3d0)then
  dist_sign = +1.0d0
 elseif(x(SDIM) .ge. 0.05d0 .and. &
       val1 .ge. 0.0d0 .and. &
       val2 .ge. 0.0d0)then
  dist_sign = +1.0d0  
 else
  dist_sign = -1.0d0
 endif
 
 ! determine distance
 if(x(SDIM) .ge. 0.3d0 )then
  if(x(1) .ge. 0.25d0 .or. x(2) .ge. 0.25d0)then
   dist = x(SDIM) - 0.3d0
  else
   call dist_point_to_line(x1,x2,x,dist1)
   call dist_point_to_line(x2,x4,x,dist2)
   if(dist1 .lt. 0.0d0 .or. dist2 .lt. 0.0d0)then
    print *,"invalid dist1 or dist2"
    stop
   endif
   dist = min(dist1,dist2) 
  endif
 else
  plg1(1,1) = 0.25d0
  plg1(1,2) = 0.0d0
  plg1(1,SDIM)= 0.3d0
  plg1(2,1) = 0.25d0
  plg1(2,2) = 0.25d0
  plg1(2,SDIM)= 0.3d0
  plg1(3,1) = 0.125d0
  plg1(3,2) = 0.125d0
  plg1(3,SDIM)= 0.05d0
  plg1(4,1) = 0.125d0
  plg1(4,2) = 0.0d0
  plg1(4,SDIM)= 0.05d0
  ! x1 -> x2 -> x5 -> x3 
  ! from (0,0,0) its counterclockwise
  call find_p_plane_dist_wb(+1, 4, plg1, x, dist1)

  plg2(1,1) = 0.25d0
  plg2(1,2) = 0.25d0
  plg2(1,SDIM)= 0.3d0
  plg2(4,1) = 0.125d0
  plg2(4,2) = 0.125d0
  plg2(4,SDIM)= 0.05d0
  plg2(2,1) = 0.0d0
  plg2(2,2) = 0.25d0
  plg2(2,SDIM)= 0.3d0
  plg2(3,1) = 0.0d0
  plg2(3,2) = 0.125d0
  plg2(3,SDIM)= 0.05d0
  ! from (0,0,0) it's conterclockwise
  call find_p_plane_dist_wb(+1, 4, plg2, x, dist2) 
  
  if(x(1) .ge. 0.25d0 .or. x(2) .ge. 0.25d0)then
   dist3 = x(SDIM) - 0.3d0
   dist = dist_sign*min(abs(dist1),abs(dist2),abs(dist3)) 
  elseif(x(1) .le. 0.125d0 .and. x(2) .le. 0.125d0)then
   dist3 = x(SDIM) - 0.05d0
   dist = dist_sign*min(abs(dist1),abs(dist2),abs(dist3))
  else
   dist = dist_sign*min(abs(dist1),abs(dist2))   
  endif

 endif
endif
!------------------------------------------------------------

!  
! 
elseif(cavity_type .eq. 3) then   ! 3--> vshape
!**************************************************************
!             *            *  
!             *    x3      x4 
!             *    ********* 
!             *   *
!             * *
!             *
!            x1    
!**************************************************************


if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 ! determine sign
 if(x(1) .ge. 0.25d0)then
  if(x(2) .ge. 0.3d0)then
   dist_sign = 1.0d0
  else
   dist_sign = -1.0d0
  endif 
 else
  vshape = x(1)-x(2)+0.05d0
  if(vshape .gt. 0.0d0)then
   dist_sign = -1.0d0
  else
   dist_sign = 1.0d0
  endif
 endif
 ! determine distance
  x1(1) = 0.0d0
  x1(2) = 0.05d0 
  x3(1) = 0.25d0
  x3(2) = 0.3d0
  x4(1) = 0.6d0
  x4(2) = 0.3d0
  call dist_point_to_line(x1,x3,x,dist1) 
  call dist_point_to_line(x3,x4,x,dist3)  
  dist = dist_sign * min(dist1,dist3)

endif

if(coord_type .eq. 1)then   !  3D cartisian  
 ! determine sign
 if(x(SDIM) .ge. 0.3d0)then
  dist_sign = +1.0d0
 elseif(x(SDIM) .lt.0.05d0)then
  dist_sign = -1.0d0
 else
  !vshape is the radius of the circle section cut through the cone
  radius = x(SDIM) - 0.05d0
  vshape =  sqrt(x(1)**2.0d0+x(2)**2.d0)
  if( (vshape-radius) .ge. 0.0d0)then
   dist_sign = -1.0d0
  else
   dist_sign = +1.0d0
  endif
 endif

! determine distance
 vshape = sqrt(x(1)**2.0d0 + x(2)**2.0d0)

 if(x(SDIM) .ge. 0.3d0 )then
  if(vshape .ge. 0.25d0)then
   dist = x(SDIM) - 0.3d0
  elseif(vshape .ge. 0.0d0)then
   bb = -x(1)-x(2)+0.55d0
   if(bb .gt. 0.0d0)then
    dist = (0.25d0 - vshape + x(SDIM) - 0.3d0)/sqrt(2.0d0)
    if(dist .lt. 0.0d0)then
     print *,"invalid sign of dist"
     stop
    endif
   else
    dist = sqrt((0.25d0-vshape)**2.0d0+(x(SDIM)-0.3d0)**2.0d0)
   endif   
  else
   print *,"invalid vshape"
  endif
 else  ! x(SDIM) < 0.3
  v1(1) = 0.0d0
  v1(2) = 0.05d0
  v1(SDIM) = 0.0d0
  v2(1) = 0.25d0
  v2(2) = 0.3d0 
  v2(SDIM) = 0.0d0
  v3(1) = vshape
  v3(2) = x(SDIM)
  v3(SDIM) = 0.0d0
  call dist_point_to_line(v1,v2,v3,dist1)
  dist = dist_sign*dist1
  if(vshape .ge. 0.25d0)then
   dist2 =  0.3d0 - x(SDIM)
   dist = dist_sign*min(abs(dist1),abs(dist2))
  endif

 endif
endif

!------------------------------------------------------------

!  
! 
elseif(cavity_type .eq. 4) then   ! 4--> semicircle
!**************************************************************
!             *                 *  
!             *         x1      x2
!             *         ********* 
!             *        *
!             *       *
!             *    *
!             *   
!**************************************************************



if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 ! determine sign
 if(x(2) .ge. 0.3d0)then
   dist_sign = +1.0d0
 else
   radius = 0.25d0 - sqrt((x(1) - 0.0d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0)
   if(radius .ge. 0.0d0)then
    dist_sign = +1.0d0
   else
    dist_sign = -1.0d0
   endif  
 endif

 ! determine distance
  x1(1) = 0.25d0
  x1(2) = 0.3d0 
  x2(1) = 0.6d0
  x2(2) = 0.3d0
  if(x(2) .gt. 0.3d0)then
   call dist_point_to_line(x1,x2,x,dist1)
   dist = dist_sign*dist1
  else
   dist2 = 0.25d0 - sqrt((x(1))**2.0d0 + (x(2) - 0.3d0)**2.0d0)   
   call dist_point_to_line(x1,x2,x,dist3)  
   dist = dist_sign * min(abs(dist2),abs(dist3))
  endif
endif

if(coord_type .eq. 1)then   !  3D cartisian
 ! sign and distance at the same time
 radius = sqrt(x(1)**2.0d0 + x(2)**2.0d0) 
 if(x(SDIM) .ge. 0.3d00)then
  if(radius .ge. 0.25d0)then
   dist = x(SDIM) - 0.3d0
  else
   dist = sqrt((0.25d0 - radius)**2.0d0 + (x(SDIM)-0.3d0)**2.0d0)
  endif
 else
  radius3d = sqrt(x(1)**2.0d0 + x(2)**2.0d0 + (x(SDIM)-0.3d0)**2.0d0) 
  if(radius .le. 0.25d0)then
   dist = 0.25d0 - radius3d
  else 
   dist1 = 0.25d0 -radius3d
   dist2 = x(SDIM) - 0.3d0
   dist =  -1.0d0*min(abs(dist1),abs(dist2))
  endif 
 endif
endif


!------------------------------------------------------------

!  
! 
elseif(cavity_type .eq. 5) then   ! 5--> spherical reentrant
!**************************************************************
!             *  x1                          x2
!             *  *****************************
!             * * 
!             * * x3
!             *   * 
!             *     *
!             *       *
!             *        *         
!             *         *
!             *         * 
!             *       *
!             *     *
!             *   *
!             * *  
!**************************************************************



if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 ! determine sign
  lcinter = sqrt(0.25d0**2.0d0 - 0.15d0**2.0d0) + 0.3
 if(x(2) .ge. 0.6d0)then
   dist_sign = +1.0d0
 elseif(x(2) .lt. lcinter)then
   radius = 0.25d0 - sqrt((x(1) - 0.0d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0)
   if(radius .ge. 0.0d0)then
    dist_sign = +1.0d0
   else
    dist_sign = -1.0d0
   endif  
 else
   if(x(1) .le. 0.15d0)then
    dist_sign = +1.0d0
   else
    dist_sign = -1.0d0
   endif
 endif

 ! determine distance
  x1(1) = 0.15d0
  x1(2) = 0.6d0 
  x2(1) = 0.6d0
  x2(2) = 0.6d0
  x3(1) = 0.15d0
  x3(2) = lcinter

 if(x(2) .ge. 0.6d0)then
  call dist_point_to_line(x1,x2,x,dist1)
  dist = dist_sign*dist1 
 elseif(x(2) .lt. lcinter)then
   dist2 = 0.25d0 - sqrt((x(1) - 0.0d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0)
   call dist_point_to_line(x1,x2,x,dist3)  
   dist = dist_sign * min(abs(dist2),abs(dist3)) 
 else
   call dist_point_to_line(x1,x3,x,dist1)
   dist = dist_sign*dist1

   if(x(1) .ge. 0.15d0)then
    dist2 = 0.25d0 - sqrt((x(1) - 0.0d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0) 
    call dist_point_to_line(x1,x2,x,dist3)  
    dist = dist_sign * min(abs(dist2),abs(dist3),abs(dist1))
   endif 
 endif

endif

if(coord_type .eq. 1)then   !  3D cartisian
 ! determine sign
 lcinter = sqrt(0.25d0**2.0d0 - 0.15d0**2.0d0) + 0.3
 radius = sqrt((x(1))**2.0d0 + (x(2))**2.0d0)
 if(x(SDIM) .ge. 0.6d0)then
  if(radius .ge. 0.15d0)then
   dist = x(SDIM) - 0.6d0
  elseif(radius .ge. 0.0d0)then
   dist = sqrt((0.15d0 - radius)**2.0d0 + (x(SDIM)-0.6d0)**2.0d0)
  else
   print *,"invalid radius"
   stop
  endif 
 elseif(x(SDIM) .ge. lcinter)then
  dist = 0.15d0 - radius
  if(radius .gt. 0.15d0)then
   dist1 = x(SDIM) - 0.6d0
   radius3d = sqrt(x(1)**2.0d0 + x(2)**2.0d0 + (x(SDIM)-0.3d0)**2.0d0) 
   dist2 = 0.25d0 - radius3d
   dist  = -1.0d0*min(abs(dist),abs(dist1),abs(dist2))
  endif
 elseif((x(SDIM) .ge. 0.0d0).or.(1.eq.1)) then
  ztest=x(SDIM)
  if (ztest.lt.0.0d0) then
   ztest=0.0d0
  else if (ztest.ge.0.0d0) then
   ! do nothing
  else
   print *,"ztest corrupt"
   stop
  endif
  bb = (lcinter - 0.3d0)/0.15d0*radius + 0.3d0 - ztest
  radius3d = sqrt(x(1)**2.0d0 + x(2)**2.0d0 + (ztest-0.3d0)**2.0d0) 
  if(radius .lt. 0.15d0)then
   dist = 0.25d0 - radius3d
   if(bb .lt. 0.0)then
    dist1 = sqrt((radius-0.15d0)**2.0d0 + (ztest - lcinter)**2.0d0)
    dist = min(abs(dist),abs(dist1))
   endif
  else
   dist1 = 0.25d0 - radius3d
   dist2 = ztest - 0.6d0
   if(dist1 .ge. 0.0d0)then
    dist = dist1
   else
    dist = -1.0d0*min(abs(dist1),abs(dist2))
   endif
  endif
 else
  print *,"out of domain"
  stop
 endif
endif

elseif(cavity_type .eq. 6) then   ! 6--> 3D pyramid
 if(SDIM .ne. 3)then
  print *,"can't have 2d pyramid"
  stop
 endif
 plg3(1,1) = 0.25d0
 plg3(1,2) = 0.0d0
 plg3(1,SDIM)= 0.3d0
 plg3(2,1) = 0.25d0
 plg3(2,2) = 0.25d0
 plg3(2,SDIM)= 0.3d0
 plg3(3,1) = 0.0d0
 plg3(3,2) = 0.0d0
 plg3(3,SDIM)= 0.05d0
 
 call make_3points_plane(plg3(1,:),plg3(2,:),plg3(3,:),pl1,d1)
 call find_p_plane_dist_wb(+1, 3, plg3, x, dist1)

 plg4(1,1) = 0.0d0
 plg4(1,2) = 0.25d0
 plg4(1,SDIM)= 0.3d0
 plg4(2,1) = 0.0d0
 plg4(2,2) = 0.0d0
 plg4(2,SDIM)= 0.05d0
 plg4(3,1) = 0.25d0
 plg4(3,2) = 0.25d0
 plg4(3,SDIM)= 0.3d0

 do dir=1,SDIM
  x1temp(dir)=plg4(1,dir) 
  x2temp(dir)=plg4(2,dir) 
  x3temp(dir)=plg4(3,dir) 
 enddo
 call make_3points_plane(x1temp,x2temp,x3temp,pl2,d2)
 call find_p_plane_dist_wb(+1, 3, plg4, x, dist2) 
 
 
 if(x(SDIM) .ge. 0.3d0)then
  if(x(1) .ge. 0.25d0 .or. x(2) .ge. 0.25d0)then
   dist = x(SDIM) - 0.3d0
  else
   if(dist1 .lt. 0.0d0 .or. dist2 .lt. 0.0d0)then
    print *,"wrong sign of dist1 or dist2"
    stop
   endif
   dist = min(abs(dist1),abs(dist2))
  endif
 else   ! x(3) < 0.3d0
  if(x(1) .le. 0.25d0 .and. x(2) .le. 0.25d0)then
   call plane_distf(x,pl1,d1,val1)    
   call plane_distf(x,pl2,d2,val2)    
   if(val1 .ge. 0.0d0 .and. &
      val2 .ge. 0.0d0)then
    dist_sign  = +1.0d0
   else
    dist_sign = -1.0d0
   endif
   dist = dist_sign*min(abs(dist1),abs(dist2))
  else
   dist3 = x(SDIM) - 0.3d0
   dist = -1.0d0*min(abs(dist1),abs(dist2),abs(dist3))
  endif
 endif
else if (cavity_type.eq.7) then
 if(coord_type .eq. 2)then   !  R-Z  axis symmetric
  if(SDIM.ne.2)then
   print *,"R-Z, SDIM should be 2"
   stop
  endif
  dist=zblob-x(SDIM)  ! dist<0 in the solid (which is on top for this case)
 elseif(coord_type.eq.1)then   ! 3-D cartesian
  dist1=tan(radblob2)*x(1)+zblob-x(SDIM)
  dist=dsign(dist1*cos(radblob2),dist1)
 else
  print *,"coord_type error"        
  stop
 endif
else if (cavity_type.eq.8) then
 if(coord_type .eq. 2)then   !  R-Z  axis symmetric (2D)
!  if(SDIM.ne.2)then
!   print *,"R-Z, SDIM should be 2"
!   stop
!  endif
  dist=x(SDIM)-zblob  ! dist<0 in the solid (which is on bot for this case)
 elseif(coord_type.eq.1)then   ! 3-D cartesian
  dist1=x(SDIM)-(tan(radblob2)*x(1)+zblob)
  dist=dsign(dist1*cos(radblob2),dist1)
 else
  print *,"coord_type error"        
  stop
 endif
elseif(cavity_type.eq.9)then
  if(coord_type.ne.1)then
   print *,"cavity_type.eq.9 has to be on 3d"
   stop
  endif
  dist1=x(1)-5.0d0
  dist2=x(1)-1.0d0
  dist3=x(2)-5.0d0
  dist4=x(2)-1.0d0
  dist5=x(SDIM)-zblob2
  dist6=x(SDIM)-zblob
  d1=dist1*dist2
  d2=dist3*dist4
  d3=dist5*dist6

  if(d1.le.0.0d0)then
    if(d2.le.0.0d0)then
       dist=min(abs(dist5),abs(dist6))
    else
      if(d3.le.0.0d0)then
        dist=min(abs(dist3),abs(dist4))
      else
        dist=sqrt((min(abs(dist3),abs(dist4)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
      endif
    endif
  else
    if(d2.le.0.0d0)then
     if(d3.le.0.0d0)then
      dist=min(abs(dist1),abs(dist2)) 
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
     endif
    else
     if(d3.le.0.0d0)then
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist3),abs(dist4)))**2.0d0)
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                (min(abs(dist3),abs(dist4)))**2.0d0+ &
                (min(abs(dist5),abs(dist6)))**2.0d0)
     endif  
    endif
  endif

  if(x(SDIM).le.zblob.and.x(SDIM).ge.zblob2 .and.&
     x(1).le.5.0d0 .and. x(1).ge.1.0d0 .and. &     
     x(2).le.5.0d0 .and. x(2).ge.1.0d0)then
    dist=-1.0d0*dist
  else
    ! do nothing
  endif
elseif(cavity_type.eq.10)then
  if(coord_type.ne.1)then
   print *,"cavity_type.eq.10 has to be on 3d"
   stop
  endif
  dist1=x(1)-5.0e-3
  dist2=x(1)-1.0e-3
  dist3=x(2)-5.0e-3
  dist4=x(2)-1.0e-3
  dist5=x(SDIM)-xblob2
  dist6=x(SDIM)-zblob
  d1=dist1*dist2
  d2=dist3*dist4
  d3=dist5*dist6

  if(d1.le.0.0d0)then
    if(d2.le.0.0d0)then
       dist=min(abs(dist5),abs(dist6))
    else
      if(d3.le.0.0d0)then
        dist=min(abs(dist3),abs(dist4))
      else
        dist=sqrt((min(abs(dist3),abs(dist4)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
      endif
    endif
  else
    if(d2.le.0.0d0)then
     if(d3.le.0.0d0)then
      dist=min(abs(dist1),abs(dist2)) 
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
     endif
    else
     if(d3.le.0.0d0)then
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist3),abs(dist4)))**2.0d0)
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                (min(abs(dist3),abs(dist4)))**2.0d0+ &
                (min(abs(dist5),abs(dist6)))**2.0d0)
     endif  
    endif
  endif

  if(x(SDIM).le.zblob.and.x(SDIM).ge.xblob2 .and.&
     x(1).le.5.0e-3 .and. x(1).ge.1.0e-3 .and. &     
     x(2).le.5.0e-3 .and. x(2).ge.1.0e-3)then
    dist=-1.0d0*dist
  else
    ! do nothing
  endif
elseif(cavity_type.eq.11)then
  if(coord_type.ne.1)then
   print *,"cavity_type.eq.11 has to be on 3d"
   stop
  endif
  dist1=x(1)-problox
  dist2=x(1)-((probhix-problox)*0.5d0+problox)
  dist3=x(2)-probloy
  dist4=x(2)-probhiy
  dist5=x(SDIM)-probloz
  dist6=x(SDIM)-(probloz+(probhix-problox)*0.5d0)
  d1=dist1*dist2
  d2=dist3*dist4
  d3=dist5*dist6

  if(d1.le.0.0d0)then
    if(d2.le.0.0d0)then
       dist=min(abs(dist5),abs(dist6))
    else
      if(d3.le.0.0d0)then
        dist=min(abs(dist3),abs(dist4))
      else
        dist=sqrt((min(abs(dist3),abs(dist4)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
      endif
    endif
  else
    if(d2.le.0.0d0)then
     if(d3.le.0.0d0)then
      dist=min(abs(dist1),abs(dist2)) 
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
     endif
    else
     if(d3.le.0.0d0)then
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist3),abs(dist4)))**2.0d0)
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                (min(abs(dist3),abs(dist4)))**2.0d0+ &
                (min(abs(dist5),abs(dist6)))**2.0d0)
     endif  
    endif
  endif

  if(x(SDIM).le.probloz+(probhix-problox)*0.5d0 .and.&
     x(1).le.(probhix-problox)*0.5d0+problox)then
    dist=-1.0d0*dist
  else
    ! do nothing
  endif
elseif(cavity_type.eq.12)then
 if(SDIM.eq.3)then
  if(coord_type.ne.1)then
   print *,"cavity_type.eq.12 has to be on 3d"
   stop
  endif
  dist1=x(1)-problox
  dist2=x(1)-probhix
  dist3=x(2)-(probloy+0.105d0)
  dist4=x(2)-(probhiy-0.105d0)
  dist5=x(SDIM)-xblob2
  dist6=x(SDIM)-zblob
  d1=dist1*dist2
  d2=dist3*dist4
  d3=dist5*dist6

  call dist_cuboid(dist1,dist2,dist3,dist4,dist5,dist6,dist)

!  if(x(SDIM).le.zblob.and.x(SDIM).ge.xblob2 .and.&
!     x(1).le.probhix .and. x(1).ge.problox .and. &     
!     x(2).le.probhiy-0.105d0.and.x(2).ge.probloy+0.105d0)then
  if(d1.le.0.0d0.and.d2.le.0.0d0.and.d3.le.0.0d0)then
   dist=-1.0d0*dist
  else
    ! do nothing
  endif

 elseif(SDIM.eq.2)then
  dist3=x(1)-(problox+0.105d0)
  dist4=x(1)-(probhix-0.105d0)
  dist5=x(SDIM)-xblob2
  dist6=x(SDIM)-zblob
  d2=dist3*dist4
  d3=dist5*dist6

  call dist_rectangular(dist3,dist4,dist5,dist6,dist)  
  if(d2.le.0.0d0.and.d3.le.0.0d0)then
   dist=-1.0d0*dist
  else
    ! do nothing
  endif
  
 else
  print *,"SDIM invalid"
  stop
 endif


else
 print *, "invalid cavity type flag"
 stop
endif

! from cm to m
dist = dist/scale_factor

end subroutine cavity_distf_13


! --------------------------------------------------
! vapor +   liquid -

subroutine cavity_distf_12(coord_type, x_in, dist)
use probcommon_module
implicit none

INTEGER_T,INTENT(in) :: coord_type
REAL_T,INTENT(in) :: x_in(SDIM)
REAL_T            :: x(SDIM)

REAL_T, intent(out) :: dist
REAL_T            :: dist_temp
REAL_T            :: dist1,dist2,dist3
REAL_T            :: dist4,dist5,dist6,dist7
REAL_T            :: d1,d2,d3
REAL_T            :: center(SDIM)
REAL_T            :: film_thickness     ! thickness of the film
REAL_T, parameter :: scale_factor=1.0d0
INTEGER_T i

film_thickness=radblob4

! convert from m to cm 
do i = 1,SDIM
  x(i) = x_in(i)*scale_factor
enddo

if ((axis_dir.ge.1).and.(axis_dir.le.6)) then

 if(coord_type.eq.1.and.axis_dir.ne.12) then   ! 3D cartisian
  if(SDIM .ne. 3)then
   print *, "coord_type is not consistent with dimension"
   stop
  endif
 elseif(coord_type.eq.2.and.axis_dir.ne.12) then ! R-Z axis symmetric
  if(SDIM .ne. 2) then
   print *,"coord_type is not consistent with dimension"
   stop
  endif
 else
  print *,"invalid coord_type"
  stop
 endif

 if(coord_type .eq. 2)then 
   center(1) = 0.0d0
   center(2) = 0.05d0+film_thickness

   call l2norm(center, x , dist_temp)
   dist = 0.1d0 - dist_temp
   if(film_thickness .lt. 0.0d0)then
    print *,"film_thickness cannot be zero"
    stop
   elseif(film_thickness .lt. 1.0e-12)then
    ! do nothing
   else
    if(x(1) .gt. 0.1d0)then
     dist2 = 0.05d0+film_thickness - x(2)
     dist = sign(min(abs(dist2),abs(dist)), dist2)
    else
      ! do nothing
    endif
   endif
 endif

 if(coord_type .eq. 1)then 
   center(1) = 0.0d0
   center(2) = 0.0d0
   center(SDIM) = 0.05d0+film_thickness
   call l2norm(center, x , dist_temp)
   dist = 0.1d0 - dist_temp
   if(film_thickness .lt. 0.0d0)then
    print *,"film_thickness cannot be zero"
    stop
   elseif(film_thickness .lt. 1.0e-12)then
    ! do nothing
   else
    if(sqrt(x(1)**2.0d0 + x(2)**2.0d0) .gt. 0.1d0)then
     dist2 = 0.05d0+film_thickness - x(SDIM)
     dist = sign(min(abs(dist2),abs(dist)), dist2)
    else
      !  do nothing
    endif
   endif
 endif

else if (axis_dir.eq.7.or. &
         axis_dir.eq.8.or. &
         axis_dir.eq.9) then
 
 if(coord_type .eq. 2)then 

 if (SDIM.eq.2.and.axis_dir.eq.7) then
  ! do nothing
 else
  print *,"expecting 2D"
  stop
 endif
 center(1) = 0.0d0
 center(2) = 0.0d0
 center(SDIM) = zblob
 call l2norm(center, x , dist_temp)
! vapor +   liquid -
 dist=radblob-dist_temp
 elseif(coord_type .eq. 1 .and.axis_dir.eq.7)then
  if (SDIM.eq.3) then
   ! do nothing
  else
   print *,"expecting 3D"
   stop
  endif
  center(1) = xblob
  center(2) = yblob
  center(SDIM) =zblob+xblob*tan(radblob2)
  call l2norm(center,x,dist_temp)
  dist=radblob-dist_temp
 else
 print *,"coord_type invalid"
 stop
 endif
elseif(axis_dir.eq.10)then

  if(coord_type.ne.1)then
   print *,"cavity_type.eq.10 has to be on 3d"
   stop
  endif
  dist1=x(1)-(5.0e-3+zblob6)
  dist2=x(1)-(1.0e-3-zblob6)
  dist3=x(2)-(5.0e-3+zblob6)
  dist4=x(2)-(1.0e-3-zblob6)
  dist5=x(SDIM)-(xblob2-zblob6)
  dist6=x(SDIM)-(xblob2+1.0e-3)
  d1=dist1*dist2
  d2=dist3*dist4
  d3=dist5*dist6

  if(d1.le.0.0d0)then
    if(d2.le.0.0d0)then
       dist=min(abs(dist5),abs(dist6))
    else
      if(d3.le.0.0d0)then
        dist=min(abs(dist3),abs(dist4))
      else
        dist=sqrt((min(abs(dist3),abs(dist4)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
      endif
    endif
  else
    if(d2.le.0.0d0)then
     if(d3.le.0.0d0)then
      dist=min(abs(dist1),abs(dist2)) 
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
     endif
    else
     if(d3.le.0.0d0)then
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist3),abs(dist4)))**2.0d0)
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                (min(abs(dist3),abs(dist4)))**2.0d0+ &
                (min(abs(dist5),abs(dist6)))**2.0d0)
     endif  
    endif
  endif

  if(x(SDIM).le.xblob2+1.0e-3.and.x(SDIM).ge.xblob2-zblob6 .and.&
     x(1).le.5.0e-3+zblob6 .and. x(1).ge.1.0e-3-zblob6 .and. &     
     x(2).le.5.0e-3+zblob6 .and. x(2).ge.1.0e-3-zblob6)then
    ! do nothing
  else
    dist=-1.0d0*dist
    ! do nothing
  endif
elseif(axis_dir.eq.11)then

  if(coord_type.ne.1)then
   print *,"cavity_type.eq.10 has to be on 3d"
   stop
  endif
  dist1=x(1)-((probhix-problox)*0.5d0+problox+0.000375d0)
  dist2=x(1)-((probhix-problox)*0.5d0+problox)
  dist3=x(2)-(probloy)
  dist4=x(2)-(probhiy)
  dist5=x(SDIM)-(probloz+(probhix-problox)*0.5d0+0.000375d0)
  dist6=x(SDIM)-(probloz)
  d1=dist1*dist2
  d2=dist3*dist4
  d3=dist5*dist6

  if(d1.le.0.0d0)then
    if(d2.le.0.0d0)then
       dist=min(abs(dist5),abs(dist6))
    else
      if(d3.le.0.0d0)then
        dist=min(abs(dist3),abs(dist4))
      else
        dist=sqrt((min(abs(dist3),abs(dist4)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
      endif
    endif
  else
    if(d2.le.0.0d0)then
     if(d3.le.0.0d0)then
      dist=min(abs(dist1),abs(dist2)) 
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist5),abs(dist6)))**2.0d0)
     endif
    else
     if(d3.le.0.0d0)then
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                  (min(abs(dist3),abs(dist4)))**2.0d0)
     else
      dist=sqrt((min(abs(dist1),abs(dist2)))**2.0d0+ &
                (min(abs(dist3),abs(dist4)))**2.0d0+ &
                (min(abs(dist5),abs(dist6)))**2.0d0)
     endif  
    endif
  endif
  
  dist7=x(SDIM)-15e-3

   if(x(SDIM).le.probloz+0.00075d0 .and.&
      x(1).le.(probhix-problox)*0.5d0+problox+0.000375d0 .and. &
      x(1).ge.(probhix-problox)*0.5d0+problox)then      
     ! do nothing 
  else
    if(dist7.ge.0.0d0)then
      dist=dist7
    else
     dist=-1.0d0*min(dist,abs(dist7))
    ! do nothing
    endif
  endif
elseif(axis_dir.eq.12)then
 if(SDIM.eq.3)then
  if(coord_type.ne.1)then
   print *,"cavity_type.eq.11 has to be on 3d"
   stop
  endif
  dist1=x(1)-problox
  dist2=x(1)-probhix
  dist3=x(2)-(probloy+0.105d0)
  dist4=x(2)-(probhiy-0.105d0)
  dist5=x(SDIM)-(xblob2-zblob6)
  dist6=x(SDIM)-(xblob2+zblob6)
  d1=dist1*dist2
  d2=dist3*dist4
  d3=dist5*dist6

  call dist_cuboid(dist1,dist2,dist3,dist4,dist5,dist6,dist)

  dist7=x(SDIM)-0.3d0

!  if(x(SDIM).le.(xblob2+zblob6).and.x(SDIM).ge.(xblob2-zblob6) .and.&
!     x(1).le.probhix .and. x(1).ge.problox .and. &     
!     x(2).le.probhiy-0.1d0.and.x(2).ge.probloy+0.1d0)then
  if(d1.le.0.0d0.and.d2.le.0.0d0.and.d3.le.0.0d0)then
     ! do nothing  
  else
    if(dist7.ge.0.0d0)then
     dist=dist7
    else
     dist=-1.0d0*min(abs(dist7),dist)
    endif
  endif
 elseif(SDIM.eq.2)then
  dist3=x(1)-(problox+0.105d0)
  dist4=x(1)-(probhix-0.105d0)
  dist5=x(SDIM)-(xblob2-zblob6)
  dist6=x(SDIM)-(xblob2+zblob6)
  d2=dist3*dist4
  d3=dist5*dist6

  call dist_rectangular(dist3,dist4,dist5,dist6,dist)

  dist7=x(SDIM)-0.3d0

!  if(x(SDIM).le.(xblob2+zblob6).and.x(SDIM).ge.(xblob2-zblob6) .and.&
!     x(1).le.probhix .and. x(1).ge.problox .and. &     
!     x(2).le.probhiy-0.1d0.and.x(2).ge.probloy+0.1d0)then
  if(d2.le.0.0d0.and.d3.le.0.0d0)then
     ! do nothing  
  else
    if(dist7.ge.0.0d0)then
     dist=dist7
    else
     dist=-1.0d0*min(abs(dist7),dist)
    endif
  endif
     

 else
   print *,"SDIM invalid "
   stop
 endif

else
 print *,"axis_dir invalid"
 stop
endif


! from cm to m         
 dist = dist/scale_factor
! m      cm  
end subroutine cavity_distf_12


! dist>0 in the fluid 2d or 3d
subroutine CAVITY_soliddist(x,dist,im) 
use probcommon_module
use global_utility_module
implicit none
REAL_T, INTENT(in), dimension(SDIM) :: x !spatial coordinates
INTEGER_T, INTENT(in) :: im
REAL_T, INTENT(out) :: dist

if (num_materials.eq.3) then
 ! do nothing
else
 print *,"expecting num_materials.eq.3"
 stop
endif

if ((im.lt.1).or.(im.gt.num_materials)) then
 print *,"im invalid11"
 stop
endif

if (FSI_flag(im).eq.1) then ! prescribed solid (EUL)
 ! do nothing
else
 print *,"FSI_flag(im) invalid"
 stop
endif

dist=99999.0

! CAVITY_soliddist: dist>0 in fluid 2d or 3d
if (probtype.eq.710) then 
 call cavity_distf_13(axis_dir,levelrz+1,x,dist)
else
 print *,"expecting probtype.eq.710"
 stop
endif

end subroutine CAVITY_soliddist


! fluids tessellate the domain, solids are immersed. 
subroutine CAVITY_PHASE_CHANGE_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(out) :: LS(nmat)
INTEGER_T :: im
INTEGER_T :: im_solid_materialdist

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  im_solid_materialdist=im_solid_primary()

  if (probtype.eq.710) then

   do im=1,nmat
    if (FSI_flag(im).eq.1) then
     call CAVITY_soliddist(x,LS(im),nmat)  ! returns LS<0 in solid
     LS(im)=-LS(im)   ! now LS>0 in solid
    endif
   enddo

   if(axis_dir.eq.8.or.axis_dir.eq.9)then
!    print *,"axis_dir=8"
    LS(1)=1000.0d0
    LS(2)=-1000.0d0
   else
    ! vapor +   liquid -
    call cavity_distf_12(levelrz+1,x,LS(2))
    LS(1)=-LS(2)
   endif

   if (nmat.eq.3) then
    if (im_solid_materialdist.ne.3) then
     print *,"expecting im_solid_materialdist=3"
     stop
    endif
   else
    print *,"nmat invalid"
    stop
   endif

   if (is_rigid(nmat).ne.1) then
    print *,"expecting last material to be rigid"
    stop
   endif
  else
   print *,"expecting probtype.eq.710"
   stop
  endif

return
end subroutine CAVITY_PHASE_CHANGE_LS

subroutine CAVITY_PHASE_CHANGE_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: dx(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, INTENT(in) :: velsolid_flag

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

do dir=1,SDIM
 VEL(dir)=zero
enddo

if (probtype.eq.710) then

 do dir=1,SDIM
  VEL(dir)=zero
 enddo

else
 print *,"expecting probtype==710"
 stop
endif

return 
end subroutine CAVITY_PHASE_CHANGE_VEL


subroutine CAVITY_PHASE_CHANGE_PRES(x,t,LS,PRES,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: PRES
INTEGER_T :: gravity_dir
REAL_T :: gravity_dz

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
PRES=zero

call fort_derive_gravity_dir(gravity_vector,gravity_dir)

if (gravity_dir.eq.1) then
 gravity_dz=x(gravity_dir)-probhix
else if (gravity_dir.eq.2) then
 gravity_dz=x(gravity_dir)-probhiy
else if ((gravity_dir.eq.SDIM).and.(SDIM.eq.3)) then
 gravity_dz=x(gravity_dir)-probhiz
else
 print *,"gravity_dir invalid"
 stop
endif
PRES=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*gravity_dz

return 
end subroutine CAVITY_PHASE_CHANGE_PRES

subroutine CAVITY_PHASE_CHANGE_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
REAL_T   :: rtemp
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
if (probtype.eq.710) then
  do im=1,num_materials
   ibase=(im-1)*num_state_material
     ! density prescribed in the inputs file.
   STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im) 

   if (axis_dir.eq.8) then  ! axis_dir is cavity_type
    if (t.eq.zero) then
     if(x(SDIM).lt.zblob2)then   ! (zblob2,xblob2)heater
       STATE(ibase+ENUM_TEMPERATUREVAR+1)=xblob4
     elseif(x(SDIM).le.xblob2)then
       STATE(ibase+ENUM_TEMPERATUREVAR+1)=xblob4+zblob5
     elseif(x(SDIM).gt.xblob2.and.x(SDIM).le.zblob)then
      STATE(ibase+ENUM_TEMPERATUREVAR+1)=(xblob4+zblob5)-&
       (xblob4+zblob5-yblob4)/ &
       (zblob-xblob2)*(x(SDIM)-xblob2)
     elseif(x(SDIM).gt.zblob.and.x(SDIM).le.zblob+xblob5)then
      STATE(ibase+ENUM_TEMPERATUREVAR+1)=yblob4-&
       (yblob4-xblob4)/xblob5*(x(SDIM)-zblob)
     else
      STATE(ibase+ENUM_TEMPERATUREVAR+1)=xblob4
!      STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im)
     endif
    elseif(t.gt.zero)then
      ! do nothing 
    else
     print *,"t invalid"
     stop
    endif
   elseif(axis_dir.eq.9)then
    if(t.eq.zero)then
     if(x(SDIM).le.zblob.and.x(SDIM).ge.zblob2 .and.&
      x(1).le.5.0d0 .and. x(1).ge.1.0d0 .and. &     
      x(2).le.5.0d0 .and. x(2).ge.1.0d0)then
      
      if(x(SDIM).le.xblob2)then
       STATE(ibase+ENUM_TEMPERATUREVAR+1)=xblob4+zblob5
      elseif(x(SDIM).gt.xblob2.and.x(SDIM).le.zblob)then
       STATE(ibase+ENUM_TEMPERATUREVAR+1)=(xblob4+zblob5)-&
        zblob5/(zblob-xblob2)*(x(SDIM)-xblob2)
      endif
     else
       STATE(ibase+ENUM_TEMPERATUREVAR+1)=xblob4
     endif
    elseif(t.gt.zero)then
      ! do nothing 
    else
     print *,"t invalid"
     stop
    endif
   elseif(axis_dir.eq.10.or.axis_dir.eq.11)then
    if (t.eq.zero) then
     !initial temperature in inputs
     STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im) 
    else if (t.gt.zero) then
     STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im)
    else
     print *,"t invalid"
     stop
    endif
   elseif(axis_dir.eq.12)then
    if(SDIM.eq.3)then
     rtemp=0.2d0/(2.0*4.0*atan(1.0d0))- &
              sqrt((x(2)-0.2d0)**2.0d0+(x(SDIM)-0.15d0)**2.0d0)
    elseif(SDIM.eq.2)then
     rtemp=0.2d0/(2.0*4.0*atan(1.0d0))- &
             sqrt((x(1)-0.2d0)**2.0d0+(x(SDIM)-0.15d0)**2.0d0)
    else
     print *,"SDIM invalid"
     stop
    endif

    if (t.eq.zero) then
     !initial temperature in inputs
     if(im.ne.3)then
      if(abs(LS(3)).le.0.01d0)then
       STATE(ibase+ENUM_TEMPERATUREVAR+1)=yblob4-abs(LS(3))/0.01d0*(yblob4-xblob4)
      else
       STATE(ibase+ENUM_TEMPERATUREVAR+1)=xblob4 
      endif
     elseif(im.eq.3)then
      if(rtemp.ge.0.0d0)then
       STATE(ibase+ENUM_TEMPERATUREVAR+1)=xblob4+zblob5
      else
       if((0.2/(2.0d0*4.0d0*atan(1.0d0))-rtemp).le.0.05d0)then
        STATE(ibase+ENUM_TEMPERATUREVAR+1)=xblob4+zblob5- &
               abs(rtemp)/(0.05d0-0.2d0/(2.0*4.0*atan(1.0d0)))*(zblob5-(yblob4-xblob4)) 
       else
        STATE(ibase+ENUM_TEMPERATUREVAR+1)=yblob4
       endif  ! <=0.05
      endif ! rtemp
     else
      print *,"im invalid"
      stop
     endif
    else if (t.gt.zero) then
     ! do nothing
       !      STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
    else
     print *,"t invalid"
     stop
    endif
   else
    if (t.eq.zero) then
     !initial temperature in inputs
     STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im) 
    else if (t.gt.zero) then
     STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
    else
     print *,"t invalid"
     stop
    endif
   endif  ! cavity_type

   ! initial species in inputs?
  do n=1,num_species_var
   STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo

  ! bcflag: 0=called from initialize  1=called from bc
!    if (im.eq.1) then
!    call outside_temperature(t,x(1),x(2),x(SDIM),water_temp,im,bcflag)
!   STATE(ibase+ENUM_TEMPERATUREVAR+1)=water_temp  
!  endif

 enddo ! im=1..num_materials
else
 print *,"probtype invalid"
 stop
endif
 
return
end subroutine CAVITY_PHASE_CHANGE_STATE

 ! dir=1..sdim  side=1..2
subroutine CAVITY_PHASE_CHANGE_LS_BC(xwall,xghost,t,LS, &
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
 call CAVITY_PHASE_CHANGE_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CAVITY_PHASE_CHANGE_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine CAVITY_PHASE_CHANGE_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
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
if (probtype.eq.710) then
 ! do nothing
else
 print *,"expecting probtype==710"
 stop
endif
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call CAVITY_PHASE_CHANGE_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)

 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine CAVITY_PHASE_CHANGE_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine CAVITY_PHASE_CHANGE_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
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

 call CAVITY_PHASE_CHANGE_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CAVITY_PHASE_CHANGE_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine CAVITY_PHASE_CHANGE_STATE_BC(xwall,xghost,t,LS, &
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
 call CAVITY_PHASE_CHANGE_STATE(xghost,t,LS,local_STATE,local_bcflag, &
         nmat,num_state_material)
 ibase=(im-1)*num_state_material
 STATE=local_STATE(ibase+istate)
 im_crit=1
 call get_primary_material(LS,im_crit)
 ibase=(im_crit-1)*num_state_material
 STATE_merge=local_STATE(ibase+istate)
else
 print *,"istate invalid"
 stop
endif

return
end subroutine CAVITY_PHASE_CHANGE_STATE_BC


subroutine CAVITY_PHASE_CHANGE_velfreestream(problen,local_buffer)
use probcommon_module
IMPLICIT NONE
REAL_T, INTENT(inout) :: local_buffer(2*SDIM)
REAL_T, INTENT(in)    :: problen(SDIM)
REAL_T :: buf
INTEGER_T :: ibuf
INTEGER_T :: dirbc,side

if (probtype.eq.710) then
 dirbc=SDIM
 side=2
 buf=problen(SDIM)-yblob2
 ibuf=(side-1)*SDIM+dirbc
 if (local_buffer(ibuf).eq.zero) then
  local_buffer(ibuf)=buf
 endif
else
 print *,"expecting probtype==710"
 stop
endif

return
end subroutine CAVITY_PHASE_CHANGE_velfreestream

end module CAVITY_PHASE_CHANGE_module
