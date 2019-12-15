#undef BL_LANG_CC
#define BL_LANG_FORT

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
!------------------------------------------------------------
module mmat_FVM
use GeneralClass
use probcommon_module
use mof_routines_module
use MOF_pair_module

implicit none

integer,parameter :: shapeflag = 0   
! asteroid test case 0: asteroid 1:diamand
! 2: circle

integer probtype_mmat_FVM

contains


! Yang Liu distance function routines
!---------------------------------------------------------
subroutine dist_fns(imat,x,y,dist,probtype_in)

implicit none

integer,intent(in)               :: imat,probtype_in
real(kind=8), intent(in)         :: x,y
real(kind=8), intent(out)        :: dist
real(kind=8)                     :: dist1,dist2,dist3,dist4,dist5
real(kind=8)                     :: dist6
real(kind=8)                     :: dista,distb,distc,distd
real(kind=8)                     :: d1,d2,d3,d4,d5
real(kind=8)                     :: dd1,dd2,dd3,dd4
real(kind=8)                     :: xy(2)
real(kind=8)                     :: x1(2),x2(2),x3(2),x4(2),x5(2)
real(kind=8)                     :: x6(2)
real(kind=8)                     :: xxl(2),xxr(2)
integer                          :: i
real(kind=8)                     :: r1,r2,r3
real(kind=8)                     :: center(2),cc(2)

real(kind=8)                     :: x0,y0
real(kind=8)                     :: c1,c2,tt
integer                          :: pp1,pp2
real(kind=8)                     :: pcurve_crit(2)
integer                          :: cenflag,cccflag

real(kind=8)                     :: vt1(2),vt2(2),vt3(2),vt4(2)
real(kind=8)                     :: va(2),vb(2),vc(2),vd(2)
real(kind=8)                     :: vra(2),vrb(2),vrc(2),vrd(2)

integer                          :: ppcrit

real(kind=8)                     :: crossp(3)
real(kind=8)                     :: r_inner  ! radius of inner circle
real(kind=8)                     :: signtemp
real(kind=8)                     :: r_mat2
real(kind=8)                     :: material1_dist
real(kind=8)                     :: material2_dist
real(kind=8)                     :: material3_dist
real(kind=8)                     :: check_dist

if ((imat.ge.1).and.(imat.le.num_materials)) then
 ! do nothing
else
 print *,"imat invalid in dist_fns"
 stop
endif
 
if(probtype_in .eq. 19) then
 center(1) = 0.5d0
 center(2) = 0.5d0
  ! radcen and radeps defined in vof_cisl.F90
  ! radcen=0.25
  ! radeps=0.1
 if ((abs(radcen-0.25d0).le.1.0D-10).and. &
     (abs(radeps-0.1d0).le.1.0D-10)) then
  r1=radcen-radeps
  r2=radcen+radeps
 else
  print *,"expecting radcen=1/4, radeps=1/10"
  stop
 endif
 dist1= sqrt((x-center(1))**2.0d0+(y-center(2))**2.0d0)

 if (imat .eq. 1) then
   dist= r1 - dist1  ! inner material
 elseif(imat .eq. 2) then
   if(dist1 .le. r2 .and.  dist1 .ge. r1)then
      dist2 = r2 - dist1
      dist3 = dist1 - r1
      dist = min(dist2,dist3)
   else if (dist1.ge.r2) then
    dist=r2-dist1
   else if (dist1.le.r1) then
    dist=dist1-r1
   else
    print *,"dist1 bust"
    stop
   endif
 elseif(imat .eq. 3)then  
   dist = dist1 - r2  ! outer material
 else
   print *,"wrong number of mat"
   stop
 endif

elseif(probtype_in .eq. 14)then

 x0 = 2.0d0*x-1.0d0
 y0 = 2.0d0*y-1.0d0

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)
 
 tt = atan((y0-c2)/(x0-c1));  
 if((x0-c1) .ge. 0.0d0 .and. (y0-c2) .ge. 0.0d0)then
    ! do nothing
 elseif((x0-c1) .le. 0.0d0 .and. (y0-c2) .gt. 0.0d0)then
    tt = tt + Pi;
 elseif((x0-c1) .lt. 0.0d0 .and. (y0-c2) .lt. 0.0d0)then
    tt = tt +Pi;
 else
    tt = 2.0d0*Pi + tt;
 endif

 dist1 = -(sqrt((x0-c1)**2.0d0 + (y0-c2)**2.0d0) - &
         (0.5d0 + 0.2d0*sin(5.0d0*tt)))
 if(imat .eq. 1)then
  dist = dist1
 elseif(imat .eq. 2) then
  dist = -dist1
 else
  print *,"wrong imat flag in, 130"
  stop
 endif

elseif(probtype_in.eq.13) then

 x0 = 2.0d0*x-1.0d0
 y0 = 2.0d0*y-1.0d0

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)
 
 tt = atan((y0-c2)/(x0-c1));  
 if((x0-c1) .ge. 0.0d0 .and. (y0-c2) .ge. 0.0d0)then
    ! do nothing
 elseif((x0-c1) .le. 0.0d0 .and. (y0-c2) .gt. 0.0d0)then
    tt = tt + Pi;
 elseif((x0-c1) .lt. 0.0d0 .and. (y0-c2) .lt. 0.0d0)then
    tt = tt +Pi;
 else
    tt = 2.0d0*Pi + tt;
 endif

 dist1 = -(sqrt((x0-c1)**2.0d0 + (y0-c2)**2.0d0) - &
         (0.5d0 + 0.2d0*sin(5.0d0*tt)))
 dist2 = -(sqrt((x0-c1)**2.0d0 + (y0-c2)**2.0d0) - &
         (0.5d0 + pentaeps + 0.2d0*sin(5.0d0*tt)))

 if(imat .eq. 1)then
  dist = dist1  
 elseif(imat .eq. 3)then
  dist = -dist2
 elseif(imat .eq. 2)then
  if(dist1 .lt. 0.0d0 .and. dist2 .gt. 0.0d0)then
   dist = min(abs(dist1),abs(dist2))
  else
   dist = -min(abs(dist1),abs(dist2))
  endif
 else
  print *,"imat invalid for probtype_in.eq.13; imat=",imat
  stop
 endif

elseif(probtype_in.eq.15)then     ! asteroid 2 materials back
 xy(1)=x
 xy(2)=y

! center(1) = 0.02d0*sqrt(5.0d0)
! center(2) = 0.02d0*sqrt(5.0d0)
 center=0.0d0

 cc(1)=(center(1)+1.0d0)/2.0d0
 cc(2)=(center(2)+1.0d0)/2.0d0
!  print *,"xy",xy
 call rad_cal(xy,cc,tt)
 if((tt.ge.0.0d0).and.(tt.le.0.5d0*Pi))then
   pp1=1
   pp2=pcurve_num/4+1
 elseif(tt.le.Pi)then
   pp1=pcurve_num/4+1
   pp2=pcurve_num/2+1
 elseif(tt.le.1.5d0*Pi)then
   pp1=pcurve_num/2+1
   pp2=(pcurve_num/4)*3+1
 elseif(tt.le.2.0d0*Pi)then
   pp1=(pcurve_num/4)*3+1
   pp2=pcurve_num+1
 else
  print *,"invalid tt", tt
  stop
 endif 

 dist1=100000.0d0
! smrad=100000.0d0
 ppcrit=pp1

 if (pp1.le.pp2) then
  do i=pp1,pp2
   call l2normd(2,xy,pcurve_ls(:,i),dist2)
   if (dist2.lt.dist1)then
    pcurve_crit=pcurve_ls(:,i)
    dist1=dist2
    ppcrit=i
   endif
  enddo ! i=pp1,pp2
 else
  print *,"pp1 or pp2 invalid"
  stop
 endif
 

 !------------------------------------------------
  if(ppcrit.eq.pp2) then
   call l2normd(2,xy,pcurve_ls(:,ppcrit),dist5)
  else if ((ppcrit.ge.pp1).and.(ppcrit.lt.pp2)) then
   call dist_point_to_lined(2,pcurve_ls(:,ppcrit), &
                           pcurve_ls(:,ppcrit+1),xy,dist5)
  else
   print *,"ppcrit invalid"
   stop
  endif

 ! -------------------------------------------------
 check_dist=sqrt((x-pcurve_ls(1,ppcrit))**2+ &
                 (y-pcurve_ls(2,ppcrit))**2)
 if (check_dist.le.1.0D-14) then
  dist=0.0d0
  crossp=0.0d0
 else if (check_dist.ge.1.0D-14) then
  if(ppcrit.eq.pp2)then
   call YangLiu_cross_product(2,pcurve_ls(:,ppcrit-1),&
                     pcurve_ls(:,ppcrit),xy,crossp)
  else if ((ppcrit.ge.pp1).and.(ppcrit.lt.pp2)) then
   call YangLiu_cross_product(2,pcurve_ls(:,ppcrit), & 
                   pcurve_ls(:,ppcrit+1),xy,crossp) 
  else
   print *,"ppcrit invalid"
   stop
  endif

  if(crossp(3).eq.0.0d0)then
   dist=0.0d0
   crossp=0.0d0
!  print *,"YangLiu_cross_product is 0 check!"
!  print *,"xy",xy
!  print *, pcurve_ls(:,ppcrit-1)
!  print *,pcurve_ls(:,ppcrit)
!  print *,pcurve_ls(:,ppcrit+1)
!  print *,"crossp",crossp
!  stop
  endif

  dist=sign(min(dist1,dist5),crossp(3))

  if (imat.eq.1)then
   ! do nothing
  elseif (imat.eq.2) then
   dist=-dist
  else
   print *,"wrong num of materials for test 5"
   stop
  endif

 else
  print *,"check_dist invalid"
  stop
 endif

elseif(probtype_in .eq. 20)then
 r_inner=0.1d0
 xy(1)=x
 xy(2)=y
 cenflag=0    ! 0= axis symmetry aligned with grid
              ! 1= not aligned with grid 
 cccflag =0
 if(cenflag .eq. 1)then 
  center(1) = 0.02d0*sqrt(5.0d0)
  center(2) = 0.02d0*sqrt(5.0d0)
 elseif(cenflag .eq. 0)then
  center(1) = 0.0d0
  center(2) = 0.0d0
 else
  print *,"cenflag error"
  stop
 endif

 cc(1)=(center(1)+1.0d0)/2.0d0
 cc(2)=(center(2)+1.0d0)/2.0d0


   vra(1)=cc(1)+r_inner
   vra(2)=cc(2)
   vrb(1)=cc(1)
   vrb(2)=cc(2)+r_inner
   vrc(1)=cc(1)-r_inner
   vrc(2)=cc(2)
   vrd(1)=cc(1)
   vrd(2)=cc(2)-r_inner

   va(1)=cc(1)+0.4d0
   va(2)=cc(2)
   vb(1)=cc(1)
   vb(2)=cc(2)+0.4d0
   vc(1)=cc(1)-0.4d0
   vc(2)=cc(2)
   vd(1)=cc(1)
   vd(2)=cc(2)-0.4d0


 vt1(1)=0.6d0
 vt1(2)=0.5d0
 vt2(1)=0.5d0
 vt2(2)=0.6d0
 vt3(1)=0.4d0
 vt3(2)=0.5d0
 vt4(1)=0.5d0
 vt4(2)=0.4d0
 call dist_point_to_lined(2,vt1,vt2,xy,dd1)
 call dist_point_to_lined(2,vt2,vt3,xy,dd2) 
 call dist_point_to_lined(2,vt3,vt4,xy,dd3)
 call dist_point_to_lined(2,vt4,vt1,xy,dd4)

 dist6=min(dd1,dd2,dd3,dd4)

 d1=x+y-1.1d0
 d2=x+y-0.9d0
 d3=x-y-0.1d0
 d4=x-y+0.1d0

 if((d1.le.0.0d0).and.(d2.ge.0.0d0).and.(d3.le.0.0d0).and. &
    (d4.ge.0.0d0))then
  ! do nothing
 else
  dist6 = -dist6
 endif

 if(cenflag .eq. 0)then              ! center aligned with grid
  if(abs(xy(1) - cc(1)) .lt. 1.0e-12 &
     .and. abs(xy(2)-cc(2)) .lt. 1.0e-12)then
   cccflag=1    
  else
   call rad_cal(xy,cc,tt)
  endif 
 elseif(cenflag .eq. 1)then    ! center not aligned with grid
  call rad_cal(xy,cc,tt)
 endif

 if(cccflag .eq. 0)then

  if(tt .ge. 0.0d0 .and. tt .le. 0.5d0*Pi)then
   pp1=1
   pp2=pcurve_num/4+1
  elseif(tt .le. Pi)then
   pp1=pcurve_num/4+1
   pp2=pcurve_num/2+1
  elseif(tt .le. 1.5d0*Pi)then
   pp1= pcurve_num/2+1
   pp2=pcurve_num/4*3+1
  elseif(tt .le. 2.0d0*Pi)then
   pp1=pcurve_num/4*3+1
   pp2=pcurve_num+1
  else
   print *,"invalid tt", tt
   stop
  endif 


  dist1=100000.0d0
! smrad=100000.0d0
  ppcrit=pp1
  do i=pp1,pp2
   call l2normd(2,xy,pcurve_ls(:,i),dist2)
   if(dist2 .lt. dist1)then
    pcurve_crit=pcurve_ls(:,i)
    dist1=dist2
    ppcrit=i
   endif
!  if(abs(pcurve_rad(i)-tt) .lt. smrad)then
!   smrad=abs(pcurve_rad(i)-tt)
!   rad_crit= pcurve_ls(:,i)   
!  endif
  enddo

  if(pcurve_ls(1,ppcrit) .eq. x .and. pcurve_ls(2,ppcrit) .eq. y)then
   dist=0.0d0
   crossp=0.0d0
  else
   if(ppcrit .eq. pp2)then
    call YangLiu_cross_product(2,pcurve_ls(:,ppcrit-1),&
                   pcurve_ls(:,ppcrit),xy,crossp)
   else
    call YangLiu_cross_product(2,pcurve_ls(:,ppcrit), & 
                   pcurve_ls(:,ppcrit+1),xy,crossp) 
   endif

   if(crossp(3) .eq. 0.0d0)then
    dist=0.0d0
    crossp=0.0d0
!   print *,"YangLiu_cross_product is 0 check!"
!   print *,"xy",xy
!   print *, pcurve_ls(:,ppcrit-1)
!   print *,pcurve_ls(:,ppcrit)
!   print *,pcurve_ls(:,ppcrit+1)
!   print *,"crossp",crossp
!   stop
   endif
  endif

  if(ppcrit .eq. pp2)then
   ! do nothing
  else
   call dist_point_to_lined(2,pcurve_ls(:,ppcrit), &
                           pcurve_ls(:,ppcrit+1),xy,dist1)
!   write(13,*) "dist5",dist5,dist1
  endif

   call dist_point_to_lined(2,va,vra,xy,dist2)
   call dist_point_to_lined(2,vb,vrb,xy,dist3)
   call dist_point_to_lined(2,vc,vrc,xy,dist4)
   call dist_point_to_lined(2,vd,vrd,xy,dist5)

   call dist_point_to_lined(2,vra,vrb,xy,dista)
   call dist_point_to_lined(2,vrb,vrc,xy,distb)
   call dist_point_to_lined(2,vrc,vrd,xy,distc)
   call dist_point_to_lined(2,vrd,vra,xy,distd)

   if((dist2.lt.0.0d0).or.(dist3.lt.0.0d0).or. &
      (dist4.lt.0.0d0).or.(dist5.lt.0.0d0))then
    print *,"check dist2345"
    stop
   endif


!   first qudrant>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   if(tt .eq. 0.0d0)then
    if(imat .eq. 1)then
     if(crossp(3) .eq. 0.0d0)then
      dist=0.0d0
     else
      dist=sign(dist1,-crossp(3))
     endif
    elseif(imat .eq. 2)then
     if(xy(1) .lt. 0.5d0+r_inner)then
      dist=-dista
     elseif(xy(1) .gt. 0.9d0)then
      dist=-(xy(1)-0.9d0)
     else
      dist=0.0d0
     endif
    elseif(imat .eq. 3)then
     dist=-distb
    elseif(imat .eq. 4)then
     dist=-distc
    elseif(imat .eq. 5)then
     if(xy(1) .lt. 0.5d0+r_inner)then
      dist=-dista
     elseif(xy(1) .gt. 0.9d0)then
      dist=-(xy(1)-0.9d0)
     else
      dist=0.0d0
     endif
    elseif(imat .eq. 6)then
     dist=dist6
    else
      print *,"wrong material num"
      stop
    endif  

   elseif(tt .gt. 0.0d0 .and. tt .lt. 0.5d0*Pi)then

    if(imat .eq. 1)then    
     if(crossp(3) .eq. 0.0d0)then
      dist=0.0d0
     else
      dist=sign(dist1,-crossp(3))
     endif
    elseif(imat .eq. 2)then
     if(crossp(3) .ge. 0.0d0 .and. dist6 .le. 0.0d0)then
      dist = min(dist1,dist2,dist3,-dist6)
      if(dist1 .lt. 0.0d0 .or.  -dist6.lt.0.0d0)then
       print *,"check dist1,6"
       stop
      endif
     else
      dist= -min(dist1,-dist6)
     endif
    elseif(imat .eq. 3)then
     dist=-min(dist3,distb)
    elseif(imat .eq. 4)then
     dist=-distc
    elseif(imat .eq. 5)then
     dist=-dist2
    elseif(imat .eq. 6)then
     dist=dist6
    else
     print *,"wrong material num"
     stop
    endif

   elseif(tt .eq. 0.5d0*Pi)then 

    if(imat .eq. 1)then
     if(crossp(3) .eq. 0.0d0)then
      dist=0.0d0
     else
      dist=sign(dist1,-crossp(3))
     endif
    elseif(imat .eq. 2)then
     if(xy(2) .lt. 0.5d0+r_inner)then
      dist=-dista
     elseif(xy(2) .gt. 0.9d0)then
      dist=-(xy(2)-0.9d0)
     else
      dist=0.0d0
     endif
    elseif(imat .eq. 3)then
     if(xy(2) .lt. 0.5d0+r_inner)then
      dist=-distb
     elseif(xy(2) .gt. 0.9d0)then
      dist=-(xy(2)-0.9d0)
     else
      dist=0.0d0
     endif
    elseif(imat .eq. 4)then
     dist=-min(distc,dist4)
    elseif(imat .eq. 5)then
     dist=-min(dist2,distd)
    elseif(imat .eq. 6)then
     dist=dist6
    else
      print *,"wrong material num"
      stop
    endif  
    ! second quadrant >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   elseif(tt .lt. Pi)then

    if(imat .eq. 1)then
     if(crossp(3) .eq. 0.0d0)then
      dist=0.0d0
     else
      dist=sign(dist1,-crossp(3))
     endif
    elseif(imat .eq. 2)then
     dist=-min(dist3,dista)
    elseif(imat .eq. 3)then
     if(crossp(3) .gt. 0.0d0 .and. dist6 .lt. 0.0d0)then
      dist = min(dist1,dist3,dist4,-dist6)
      if(dist1 .lt. 0.0d0 .or.  -dist6.lt.0.0d0)then
       print *,"check dist1,6"
       stop
      endif
     else
      dist= - min(dist1,-dist6)
     endif
    elseif(imat .eq. 4)then
     dist=-min(dist4,distc)
    elseif(imat .eq. 5)then 
     dist=-distd
    elseif(imat .eq. 6)then
     dist=dist6
    else
     print *,"wrong material num"
     stop
    endif
   elseif(tt .eq. Pi)then
    if(imat .eq. 1)then
     if(crossp(3) .eq. 0.0d0)then
      dist=0.0d0
     else
      dist=sign(dist1,-crossp(3))
     endif
    elseif(imat .eq. 2)then
     dist=-min(dist3,dista)
    elseif(imat .eq. 3)then
     if(xy(1) .gt. 0.5d0-r_inner)then
      dist=-distb
     elseif(xy(1) .lt. 0.1d0)then
      dist=xy(1)-0.1d0
     else
      dist=0.0d0
     endif

    elseif(imat .eq. 4)then
     if(xy(1) .gt. 0.5d0-r_inner)then
      dist=-distc
     elseif(xy(1) .lt. 0.1d0)then
      dist=xy(1)-0.1d0
     else
      dist=0.0d0
     endif
    elseif(imat .eq. 5)then
     dist=-distd
    elseif(imat .eq. 6)then
     dist=dist6
    else
     print *,"wrong material num"
     stop
    endif 

   ! third quadrant
   elseif(tt .lt. 1.5d0*Pi)then

    if(imat .eq. 1)then
     if(crossp(3) .eq. 0.0d0)then
      dist=0.0d0
     else
      dist=sign(dist1,-crossp(3))
     endif
    elseif(imat .eq. 2)then
     dist=-dista
    elseif(imat .eq. 3)then
     dist=-min(dist4,distb)
 
    elseif(imat .eq. 4)then

     if(crossp(3) .gt. 0.0d0 .and. dist6 .lt. 0.0d0)then
      dist = min(dist1,dist5,dist4,-dist6)
      if(dist1 .lt. 0.0d0 .or.  -dist6.lt.0.0d0)then
       print *,"check dist1,6"
       stop
      endif
     else
      dist= - min(dist1,-dist6)
     endif
    elseif(imat .eq. 5)then
     dist=-min(dist5,distd)
    elseif(imat .eq. 6)then
     dist=dist6
    else
     print *,"wrong material num"
     stop
    endif
 
   elseif( tt .eq. 1.5d0*Pi) then

    if(imat .eq. 1)then
     if(crossp(3) .eq. 0.0d0)then
      dist=0.0d0
     else
      dist=sign(dist1,-crossp(3))
     endif
    elseif(imat .eq. 2)then
     dist = -min(dist2,dista)
    elseif(imat .eq. 3)then
     dist = -min(distb,dist4)

    elseif(imat .eq. 4)then

     if(xy(2) .gt. 0.5d0-r_inner)then
      dist=-distc
     elseif(xy(2) .lt. 0.1d0)then
      dist= xy(2)-0.1d0
     else
      dist=0.0d0
     endif

    elseif(imat .eq. 5)then
     if(xy(2) .gt. 0.5d0-r_inner)then
      dist=-distd
     elseif(xy(2) .lt. 0.1d0)then
      dist= xy(2)-0.1d0
     else
      dist=0.0d0
     endif
    elseif(imat .eq. 6)then
     dist=dist6
    else
      print *,"wrong material num"
      stop
    endif  

   elseif(tt .lt. 2.0d0*Pi)then


    if(imat .eq. 1)then
     if(crossp(3) .eq. 0.0d0)then
      dist=0.0d0
     else
      dist=sign(dist1,-crossp(3))
     endif
    elseif(imat .eq. 2)then
     dist=-min(dista,dist2)
 
    elseif(imat .eq. 3)then
     
     dist=-distb
 
    elseif(imat .eq. 4)then
      dist=-min(distc,dist5)

    elseif(imat .eq. 5)then  

     if(crossp(3) .gt. 0.0d0 .and. dist6 .lt. 0.0d0)then
      dist = min(dist1,dist2,dist5,-dist6)
      if(dist1 .lt. 0.0d0 .or.  -dist6.lt.0.0d0)then
       print *,"check dist1,6"
       stop
      endif
     else
      dist= -min(dist1,-dist6)
     endif
    elseif(imat .eq. 6)then
     dist=dist6
    else
     print *,"wrong material num"
     stop
    endif
   else
    print *,"invalid tt", tt
    stop
   endif  
 elseif(cccflag .eq. 1)then
  if(imat .eq. 1)then
   dist=-0.1d0*sqrt(2.0)
  elseif(imat .eq. 2 .or. imat .eq. 3 .or. imat .eq. 4 &
           .or. imat .eq. 5)then
   dist=-0.1d0/sqrt(2.0)
  else
   dist=+0.1d0/sqrt(2.0)
  endif
 else
  print *,"cccflag invald 353"
  stop
 endif

elseif(probtype_in .eq. 16)then      ! nucleate boiling set-up
 xy(1)=x
 xy(2)=y

  ! filament_test_type declared in vof_cisl.F90
  ! filament_test_type==0 for irregular material 2
 if (filament_test_type.eq.0) then
  ! dist1 
  x1(1)=0.0d0
  x1(2)=0.55d0
  x2(1)=0.21d0
  x2(2)=0.55d0
  x3(1)=0.79d0
  x3(2)=0.55d0
  x4(1)=1.0d0 
  x4(2)=0.55d0
  x5(1)=0.5d0
  x5(2)=0.55d0

  xxl(2)=(0.55d0**2.0+0.5d0**2.0-(0.29d0-thermal_delta)**2.0d0)/1.1d0
  xxr(2)=(0.55d0**2.0+0.5d0**2.0-(0.29d0-thermal_delta)**2.0d0)/1.1d0
  xxl(1)= 0.5d0-sqrt(0.5d0**2.0d0-xxl(2)**2.0)
  xxr(1)= 0.5d0+sqrt(0.5d0**2.0d0-xxr(2)**2.0)

  if(y .gt. 0.55d0)then
   call dist_point_to_line(x1,x2,xy,d1)
   call dist_point_to_line(x3,x4,xy,d2)
   dist1 = -min(d1,d2)
  else
   call dist_point_to_line(x1,x2,xy,d1)
   call dist_point_to_line(x3,x4,xy,d2)
   call l2normd(2,xy,x5, r1)
   d3=r1-0.29d0
   dist1=sign(min(abs(d3),abs(d1),abs(d2)),d3)
  endif 

! dist2
   x6(1)=0.5d0
   x6(2)=0.0d0
   call l2normd(2,xy,x6,r2)
   dist2= 0.5d0-r2
   
! dist3
   call l2normd(2,xy,x5,r3)
   dist3= (0.29d0-thermal_delta)-r3

  material1_dist=dist1

  call dist_point_to_arc(xy,xxr,xxl,x6,d4)
  call dist_point_to_arc(xy,xxl,xxr,x5,d5)

  if(dist2 .ge. 0.0d0 .and. dist3 .ge. 0.0d0 &
           .and. dist1 .le. 0.0d0)then
   material2_dist=min(d4,d5)
  else
   material2_dist=-min(d4,d5)
  endif

  if (material1_dist.ge.0.0d0) then
   material3_dist=-material1_dist
  else if (material2_dist.ge.0.0d0) then
   material3_dist=-material2_dist
  else
   material3_dist=min(-material1_dist,-material2_dist)
  endif 

  if(imat .eq. 1)then
   dist=material1_dist
  elseif(imat .eq. 2)then
   dist=material2_dist
  elseif(imat .eq. 3)then
   dist=material3_dist
  else
   print *,"probtype16 imat invalid"
   stop
  endif

 else if (filament_test_type.eq.1) then

  ! the signed distance function is symmetric with x=0.5d0
  ! find the reflection if x is bigger than 0.5d0
  if (x .gt. 0.5d0) then
   xy(1)=1.0d0-x
  elseif(x .le. 0.5d0)then 
   ! do nothing
  else
   print *, "x invalid", "x=",x  
   stop
  endif

  ! important points
  x1(1)=0.0d0
  x1(2)=0.55d0
  x2(1)=0.21d0
  x2(2)=0.55d0
  x3(1)=0.26d0
  x3(2)=0.5d0
  x4(1)=0.5d0
  x4(2)=0.26d0
  x5(1)=0.5d0
  x5(2)=0.5d0
  x6(1)=0.21d0
  x6(2)=0.5d0

  r_mat2=0.24d0-thermal_delta

  call l2normd(2,xy,x5,r1)
  dist1=r_mat2-r1
  
  call dist_point_to_line(x1,x2,xy,d1)
  call dist_point_to_arc(xy,x3,x2,x6,d2)
  call dist_point_to_arc(xy,x3,x4,x5,d3)

  if(xy(2) .gt. 0.55d0)then
   signtemp=-1.0d0
  elseif(xy(1) .ge. 0.21d0 .and. xy(1) .le. 0.26d0  &
         .and. xy(2) .ge. 0.5d0 .and. xy(2) .le. 0.55d0)then
   call l2normd(2,xy,x6,r2)
   dist2=0.05d0-r2
   signtemp=sign(1.0d0,dist2)
  elseif(xy(1) .ge. 0.26d0 .and. xy(1) .le. 0.5d0  &
         .and. xy(2) .ge. 0.5d0 .and. xy(2) .le. 0.55d0)then
   signtemp=-1.0d0
  elseif(xy(1) .ge. 0.26d0 .and. xy(1) .le. 0.5d0  &
         .and. xy(2) .ge. 0.26d0 .and. xy(2) .le. 0.5d0)then
   dist3=r1-0.24d0
   signtemp=sign(1.0d0,dist3)
  else
   signtemp=1.0d0 
  endif

  material2_dist=dist1
  material1_dist=signtemp*min(d1,d2,d3)

  if (material1_dist.ge.0.0d0) then
   material3_dist=-material1_dist
  else if (material2_dist.ge.0.0d0) then
   material3_dist=-material2_dist
  else
   material3_dist=min(-material1_dist,-material2_dist)
  endif 

  if(imat .eq. 2)then
   dist=material2_dist
  elseif(imat .eq. 1)then
   dist=material1_dist
  elseif(imat .eq. 3)then
   dist=material3_dist
  else
   print *,"imat invalid 1183"
   stop
  endif

 else
  print *,"filament_test_type invalid"
  stop
 endif

else
 print *,"probtype_in invalid"
 stop
endif


end subroutine dist_fns
!----------------------------------------------------
 subroutine dist_to_boundary(xy,dist)    ! distance from center to the boudary of comutational domain
 implicit none                           ! with the same rad

 real(kind=8),intent(in)   :: xy(2)
 real(kind=8)              :: cc(2)
 real(kind=8)              :: dist
 real(kind=8)              :: dd(2)
 real(kind=8)              :: theta

 cc=0.5d0
 if(cc(1) .eq. xy(1) .and. cc(2) .eq. xy(2))then
  print *,"warning, cc=xy,934"
  stop
 endif

 call rad_cal(xy,cc,theta) 
 dd=0.0d0

 if(theta .eq. 0.0d0 .or. theta .eq. Pi*0.5d0 &
     .or. theta .eq. Pi .or. theta .eq. Pi*1.5d0)then
   dist=0.5d0
 elseif(theta .eq. Pi*0.25d0 .or. theta .eq. Pi*0.75d0 &
      .or. theta .eq. Pi*1.25d0  .or. theta .eq. Pi*1.75d0)then
   dist=0.5d0*sqrt(2.0d0)
 elseif(theta .gt. Pi*0.25d0 .and. theta .lt. Pi*0.75d0)then
   dd(2)=1.0d0
   dd(1)=0.5d0/tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 elseif(theta .gt. Pi*0.75d0 .and. theta .lt. Pi*1.25d0)then
   dd(1)=0.0d0
   dd(2)=-0.5*tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 elseif(theta .gt. Pi*1.25d0 .and. theta .lt. Pi*1.75d0)then
   dd(2)=0.0d0
   dd(1)=-0.5d0/tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 elseif(theta .gt. 0.0d0 .and. theta .lt. Pi*0.25d0)then
   dd(1)=1.0d0
   dd(2)=0.5*tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 elseif(theta .gt. Pi*1.75d0 .and. theta .lt. Pi*2.0d0)then
   dd(1)=1.0d0
   dd(2)=0.5*tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 else
  print *,"check theta, 960", theta
  stop
 endif

 



 end subroutine
!------------------------------------------------- 
subroutine YangLiu_cross_product(sdim,x1,x2,cc,crossp)
implicit none

integer, intent(in)   :: sdim
real(kind=8),intent(in)  :: x1(sdim),x2(sdim),cc(sdim)
real(kind=8)           :: crossp(3)
real(kind=8)           :: u(sdim),v(sdim)
integer                :: i

 do i=1,sdim
  u(i)=x1(i)-cc(i)
  v(i)=x2(i)-cc(i)
 enddo

 if(sdim .eq. 2)then
   crossp(1)= 0.0d0
   crossp(2)= 0.0d0
   crossp(3)= u(1)*v(2)-u(2)*v(1)
 else
  print *,"check YangLiu_cross_product"
  stop
 endif

end subroutine YangLiu_cross_product
!-------------------------------------------------
subroutine getline_from_two_points(p1,p2,mx,my,c)
implicit none

real(kind=8),intent(in)    :: p1(2),p2(2)
real(kind=8),intent(out)   :: mx,my,c

if(abs(p1(1)-p2(1)) .lt. 1.0e-10 .and. abs(p1(2)-p2(2)) .lt. 1.0e-10 ) then
  print *,"p1 and p2 are coincide, 801"
  stop
elseif(abs(p1(1)-p2(1)) .lt. 1.0e-10)then
  mx=1.0d0
  my=0.0d0
  c=-p1(1)
elseif(abs(p1(2)-p2(2)) .lt. 1.0e-10)then
  mx=0.0d0
  my=1.0d0
  c=-p1(2)
else
  mx=(p1(2)-p2(2))/(p1(1)-p2(1))
  my=-1.0d0
  c=-mx*p1(1)+p1(2)
endif

end subroutine getline_from_two_points


! --------------------------------------------------
subroutine find_cloest_2d(sflag,dist,xin,xout)
implicit none

integer     ,intent(in)    :: sflag
real(kind=8),intent(in)    :: dist
real(kind=8),intent(in)    :: xin(2)
real(kind=8)               :: xout(2)
real(kind=8)               :: ss
real(kind=8)               :: xgrad(2)
real(kind=8)               :: xnorm
integer                    :: i

if(sflag .eq. 1)then
 ss = +1.0d0
elseif(sflag .eq. -1)then
 ss = -1.0d0
else
 print *,"invalid sign flag, 717"
 stop
endif

if(dist .eq. 0.0d0)then
 print *,"dist is 0, 723"
 stop
endif

xnorm = sqrt((xin(1)-0.5d0)**2.0d0 + (xin(2)-0.5d0)**2.0d0)
do i=1,2
 xgrad(i)=(xin(i)-0.5d0)/xnorm
 xout(i)= xin(i)+ss*xgrad(i)*abs(dist) 
enddo


end subroutine find_cloest_2d
!--------------------------------------------


subroutine dist_point_to_lined(sdim,p1,p2,x,dist)
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

integer,intent(in)           :: sdim
real(kind=8),intent(in)      ::  p1(sdim),p2(sdim),x(sdim)
real(kind=8)                 ::  dist


real(kind=8)                 :: diff10,diff21
real(kind=8),allocatable     :: x10(:), x21(:)
integer                      :: i
real(kind=8)                 :: s


dist = 0.0d0

allocate(x10(sdim),x21(sdim))
do i = 1,sdim
 x10(i) = p1(i) - x(i)
 x21(i) = p2(i) - p1(i)
enddo

if (maxval(abs(x21)) .lt. 10d-8)then
 print *,"p1 and p2 are coincide with each other",p1,p2
 stop
endif

call l2normd(sdim, p1, x, diff10)
call l2normd(sdim, p2, p1, diff21)

!print *,"diff10",diff10
!print *,"diff21",diff21

s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)

!write(13,*) "s=",s

if(s .gt. 1.0d0)then
 call l2normd(sdim, p2, x,dist)
elseif(s .lt. 0.0d0)then
 call l2normd(sdim,p1,x,dist)
else
! if(abs((diff10**2.0d0 )*(diff21**2.0d0) - & 
!        (dot_product(x10,x21))**2.0d0) .lt. 1.0e-10)then
!   dist= 0.0d0
! else
  dist = sqrt(((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0)/ &
        (diff21**2.0d0))
! endif
endif

deallocate(x10,x21) 

end subroutine dist_point_to_lined
! ---------------------------------------------------------------
subroutine dist_point_to_line_modify(sdim,p1,p2,x,dist,xc)
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

integer,intent(in)           :: sdim
real(kind=8),intent(in)      :: p1(sdim),p2(sdim),x(sdim)
real(kind=8)                 :: dist
real(kind=8)                 :: xc(sdim)


real(kind=8)                 :: diff10,diff21
real(kind=8),allocatable     :: x10(:), x21(:)
integer                      :: i
real(kind=8)                 :: s


dist = 0.0d0

allocate(x10(sdim),x21(sdim))
do i = 1,sdim
 x10(i) = p1(i) - x(i)
 x21(i) = p2(i) - p1(i)
enddo

if (maxval(abs(x21)) .lt. 10d-8)then
 print *,"p1 and p2 are coincide with each other",p1,p2
 stop
endif

call l2normd(sdim, p1, x, diff10)
call l2normd(sdim, p2, p1, diff21)

!print *,"diff10",diff10
!print *,"diff21",diff21

s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)

!print *,"s=",s

if(s .gt. 1.0d0)then
 call l2normd(sdim, p2, x,dist)
 xc=p2
elseif(s .lt. 0.0d0)then
 call l2normd(sdim,p1,x,dist)
 xc=p1
else
 if(abs((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0) .lt. 1.0e-10)then
   dist= 0.0d0
   xc=x
 else
  dist = sqrt(((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0)/ &
        (diff21**2.0d0))
  do i=1,sdim
   xc(i)=p1(i)+ s*(x(i)-p1(i))
  enddo
 endif
endif

deallocate(x10,x21) 

end subroutine dist_point_to_line_modify





subroutine starshape(xt)
implicit none

real(kind=8)         :: theta(pcurve_num)
real(kind=8),intent(out):: xt(2,pcurve_num)
integer          :: num
integer :: i

num = pcurve_num

do i = 1,num
 theta(i) = (i-1)*2.0d0*Pi/num
enddo

do i = 1,num
 xt(1,i) = 0.02d0*sqrt(5.0d0) + &
        (0.5d0 + 0.2d0*sin(5.0d0*theta(i)))*cos(theta(i))
 xt(2,i) = 0.02d0*sqrt(5.0d0) + &
        (0.5d0 + 0.2d0*sin(5.0d0*theta(i)))*sin(theta(i))
enddo


end subroutine starshape


subroutine starshape2(xt)
implicit none

real(kind=8)         :: theta(pcurve_num)
real(kind=8),intent(out):: xt(2,pcurve_num)
integer          :: num
integer :: i

num = pcurve_num

do i = 1,num
 theta(i) = (i-1)*2.0d0*Pi/num
enddo

do i = 1,num
 xt(1,i) = 0.02d0*sqrt(5.0d0) + &
        (0.5d0 -pentaeps + 0.2d0*sin(5.0d0*theta(i)))*cos(theta(i))
 xt(2,i) = 0.02d0*sqrt(5.0d0) + &
        (0.5d0 -pentaeps + 0.2d0*sin(5.0d0*theta(i)))*sin(theta(i))
enddo


end subroutine starshape2


subroutine asteroidshape(xt)
implicit none

real(kind=8)         :: theta(pcurve_num+1)
real(kind=8),intent(out):: xt(2,pcurve_num+1)
integer          :: num
integer :: i
integer :: cenflag
real(kind=8) :: cc(2),xvec(2)


num = pcurve_num
 cenflag= 0

do i = 1,num+1
 theta(i) = (i-1)*2.0d0*Pi/real(num,8)
enddo

if(cenflag .eq. 1)then
 cc(1) = 0.02d0*sqrt(5.0d0)
 cc(2) = 0.02d0*sqrt(5.0d0)
elseif(cenflag .eq. 0)then
 cc = 0.0d0
else
 print *,"wrong cenflag"
 stop
endif

if(shapeflag.eq.0)then

 do i = 1,num+1

  xt(1,i) = cc(1) + &
           0.6d0*cos(theta(i)) + 0.2d0*cos(3.0d0*theta(i))
  xt(2,i) = cc(2) + &
           0.6d0*sin(theta(i)) - 0.2d0*sin(3.0d0*theta(i))
  xvec(1)=xt(1,i)
  xvec(2)=xt(2,i)
  call rad_cal(xvec,cc,pcurve_rad(i))
 enddo ! i=1..num+1
elseif(shapeflag .eq. 1)then
  do i = 1,pcurve_num/4+1
   xt(1,i)= 0.5d0 -(i-1)*0.5d0/(real(pcurve_num,8)/4.0d0)
   xt(2,i)= 0.0d0 +(i-1)*0.5d0/(real(pcurve_num,8)/4.0d0)

   xt(1,pcurve_num/4+i)= 0.0d0 -(i-1)*0.5/(real(pcurve_num,8)/4.0)
   xt(2,pcurve_num/4+i)= 0.5d0 -(i-1)*0.5/(real(pcurve_num,8)/4.0)

   xt(1,pcurve_num/2+i)= -0.5d0 +(i-1)*0.5/(real(pcurve_num,8)/4.0)
   xt(2,pcurve_num/2+i)= 0.0d0 -(i-1)*0.5/(real(pcurve_num,8)/4.0)
 
   xt(1,pcurve_num/4*3+i)= 0.0d0 +(i-1)*0.5/(real(pcurve_num,8)/4.0)
   xt(2,pcurve_num/4*3+i)= -0.5d0 +(i-1)*0.5/(real(pcurve_num,8)/4.0)
  enddo 
elseif(shapeflag .eq. 2)then
 
 do i=1,num+1
  xt(1,i)= cc(1) + 0.5d0*cos(theta(i))
  xt(2,i)= cc(2) + 0.5d0*sin(theta(i))
 enddo

else
 print *,"wrong shapeflag in asteroidshape,xt"
 stop
endif

end subroutine asteroidshape




subroutine dist_point_to_arc(x,arc1,arc2,arcc,dist)
! for arc of a circle
!arc1 to arc2--arc start to arcend--counterclockwise
!arcc--arc center
implicit none

real(kind=8),intent(in)   :: x(2)
real(kind=8),intent(in)   :: arc1(2),arc2(2),arcc(2)
real(kind=8)              :: dist
real(kind=8)              :: tt(3)
real(kind=8)              :: r1,r2
integer                   :: i

dist=0.0d0
call rad_cal(arc1,arcc,tt(1))
call rad_cal(arc2,arcc,tt(2))
call rad_cal(x,arcc,tt(3))

do i=1,3
 if(tt(i) .gt. 2.0d0*Pi .or. tt(i) .lt. 0.0d0)then
  print *,"rad invalid 489"
  stop
 endif
enddo

if(tt(3) .lt. tt(1))then
 call l2normd(2,x,arc1,dist)
elseif(tt(3) .gt. tt(2))then
 call l2normd(2,x,arc2,dist)
else
 call l2normd(2,x,arcc,r1)
 call l2normd(2,arcc,arc1,r2)
 dist=abs(r1-r2)
endif

end subroutine dist_point_to_arc
!-------------------------------------------------------
subroutine rad_cal(x,c,theta)
! theta = [0,2Pi)
implicit none

real(kind=8),intent(in)    :: x(2),c(2)
real(kind=8)               :: theta             

 if(x(2) .eq. c(2) .and. x(1)-c(1) .gt. 0.0d0)then
   theta=0.0d0 
 elseif((x(1)-c(1)) .gt. 0.0d0 .and. (x(2)-c(2)) .gt. 0.0d0)then
  theta = atan((x(2)-c(2))/(x(1)-c(1))); 
 elseif(x(1) .eq. c(1) .and. x(2) .gt. c(2))then
   theta=0.5d0*Pi
 elseif((x(1)-c(1)) .lt. 0.0d0 .and. (x(2)-c(2)) .gt. 0.0d0)then
    theta = atan((x(2)-c(2))/(x(1)-c(1))); 
    theta = theta + Pi;
 elseif(x(2) .eq. c(2) .and. (x(1)-c(1)) .lt. 0.0d0)then
   theta=Pi
 elseif((x(1)-c(1)) .lt. 0.0d0 .and. (x(2)-c(2)) .lt. 0.0d0)then
    theta = atan((x(2)-c(2))/(x(1)-c(1))); 
    theta = theta +Pi;
 elseif(x(1) .eq. c(1) .and. x(2) .lt. c(2))then
   theta=1.5d0*Pi 
 elseif(x(1) .eq. c(1) .and. x(2) .eq. c(2))then
   theta=0.0d0
 else
    theta = atan((x(2)-c(2))/(x(1)-c(1))); 
    theta = 2.0d0*Pi + theta;
 endif

end subroutine rad_cal

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
integer,parameter            :: SDIM = 2
real(kind=8),intent(in)      ::  p1(SDIM),p2(SDIM),x(SDIM)
real(kind=8)                 ::  dist,dist1,dist2


real(kind=8)                 :: diff10,diff21
real(kind=8)                 :: x10(SDIM), x21(SDIM)
integer                      :: i
real(kind=8)                 :: s



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

!if (maxval(abs(x10)) .lt. 10d-8)then
! print *,"p1 and x are coincide with each other",p1,p2
! stop
!endif

call l2normd(2,p1, x, diff10)
call l2normd(2,p2, p1, diff21)


s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)

if(s .le. 1.0d0 .and. s .ge. 0.0)then

 do i = 1,SDIM
  dist = dist + (x10(i) + x21(i)*s)**2.0d0
 enddo
  dist = sqrt(dist)
else
 call l2normd(2,p1,x,dist1)
 call l2normd(2,p2,x,dist2)
 dist = min(dist1,dist2)
endif

end subroutine dist_point_to_line
!---------------------------------------------------------
subroutine l2normd(sdim,x1,x2, x1x2norm)
implicit none

integer,intent(in)       :: sdim
real(kind=8),intent(in)  :: x1(sDim),x2(SDIM)
real(kind=8)             :: x1x2norm

integer                  :: i
real(kind=8),allocatable :: diff(:)

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

end subroutine l2normd


subroutine dist_hypocycloid_back(xin,tt1,tt2,tt)
implicit none

real(kind=8),intent(in) :: xin(2)
real(kind=8)  :: tt1,tt2

real(kind=8)            :: x(2)
integer                 :: i
real(kind=8)            :: res
real(kind=8)            :: val1,val2,val3
real(kind=8)            :: tt


 do i=1,2
  x(i)=xin(i) 
 enddo
 print *,"xin",x

  call fund_hypocycloid(x,tt1,val1)
  call fund_hypocycloid(x,tt2,val2)
 if(abs(val1) .lt. 1.0e-8 .and. abs(val1).lt. abs(val2))then
  print *,"1",val1
  tt=tt1
 elseif(abs(val2) .lt. 1.0e-8 .and. abs(val2) .lt. abs(val1))then
  print *,"2",val2
  tt=tt2
 else 
 print *,"3"
 res=1.0e+5
 do while(res .gt. 1.0e-8)
  call fund_hypocycloid(x,tt1,val1)
  call fund_hypocycloid(x,tt2,val2)

  print *,"val1",tt1,val1
  print *,"val2",tt2,val2
  if(val1*val2 .ge. 0.0d0)then
   tt=0.5d0*(tt1+tt2)
   call fund_hypocycloid(x,tt,val3)  
   print *,"theta",tt1,tt2,tt
   print *,"val1 and val2 same sign 667",val1,val2,val3
   if(val3*val1 .lt. 0.0d0)then
    tt2=tt
   elseif(val3*val2 .lt. 0.0d0)then
    tt1=tt
   else
    print *,"check val1 and val2", val1,val2,val3
    stop
   endif
   res=val3
  else
   print *,"in bisection"
   tt=0.5d0*(tt1+tt2)
   call fund_hypocycloid(x,tt,val3)
!   print *,"val3",tt,val3
   if(val3*val1 .lt. 0.0d0)then
    tt2=tt
   elseif(val3*val2 .lt. 0.0d0)then
    tt1=tt
   else
    print *,"check val1 and val2", val1,val2,val3
    stop
   endif
   res=val3
  endif
 enddo
 endif


end subroutine dist_hypocycloid_back


subroutine dist_hypocycloid(xin,quadf,flag,tt)
implicit none

real(kind=8),intent(in) :: xin(2)
integer,intent(in)      :: quadf

real(kind=8)            :: x(2)
integer                 :: i
real(kind=8)            :: res
real(kind=8)            :: val1,val2,val3
real(kind=8)            :: tt
real(kind=8)            :: tt1,tt2
integer                 :: flag,step

real(kind=8),parameter  :: zdelta=1.0e-8

 flag = 0
 do i=1,2
  x(i)=xin(i) 
 enddo
! print *,"xin",x
 step=0
 if(quadf .eq. 1) then
  tt1= zdelta
  tt2= Pi/2.0d0-zdelta
 elseif(quadf .eq. 2)then
  tt1=Pi/2.0d0+zdelta 
  tt2=Pi-zdelta
 elseif(quadf .eq. 3)then
  tt1=Pi+zdelta 
  tt2=1.5d0*Pi-zdelta
 elseif(quadf .eq. 4)then
  tt1=1.5d0*Pi+zdelta 
  tt2=2.0d0*Pi-zdelta  
 else
  print *,"wrong quadrant flag"
  stop
 endif


 call fund_hypocycloid(x,tt1,val1)
 call fund_hypocycloid(x,tt2,val2)  
 if(val1 .eq. 0.0d0)then
  tt=tt1
 elseif(val2 .eq. 0.0d0)then
  tt=tt2
 elseif(val1*val2 .gt. 0.0d0)then
  flag=1
 elseif(val1*val2 .lt. 0.0d0)then
  res=1.0e+5
  step=0
  do while(res .gt. 1.0e-10 .and. step .lt. 1000)
  print *,"iter step",step 
  
  call fund_hypocycloid(x,tt1,val1)
  call fund_hypocycloid(x,tt2,val2)
   tt=0.5d0*(tt1+tt2)
   call fund_hypocycloid(x,tt,val3)
!   print *,"val3",tt,val3
   if(val3 .eq. 0.0d0)then
    !  do nothing
   elseif(val3*val1 .lt. 0.0d0)then
    tt2=tt
   elseif(val3*val2 .lt. 0.0d0)then
    tt1=tt
   else
    print *,"check val1 and val2", val1,val2,val3
    stop
   endif
   res=val3
   step=step+1 
  enddo

 else
  print *,"out of cases"
  stop
 endif  


! elseif(quadf .eq. 2) then
!  tt1=0.75d0*Pi
!  tt2=0.0d0
!  res = tt1
!  do while(res .gt. 1.0e-10 .and. step .lt. 1000)
!   step = step+1
!   if(step .gt. 100)then
!   print *,"quad",quadf,"step",step,"res",res
!   endif
!   call fund_hypocycloid(x,tt1,val1)
!   call fundd_hypocycloid(x,tt1,val2)
!   tt2=tt1- 0.5d0*val1/val2 
!   if(tt2 .gt. Pi .or. tt2 .le. 0.5d0*Pi)then
!    flag = 1
!    exit
!   endif
!   res=abs(tt1-tt2)
!   tt1=tt2
!  enddo
! elseif(quadf .eq. 3) then
!  tt1=1.25d0*Pi
!  tt2=0.0d0
!  res = tt1
!  do while(res .gt. 1.0e-10 .and. step .lt. 1000)
!   step=step+1
!   if(step .gt. 100)then
 !  print *,"quad",quadf,"step",step,"res",res
!   print *,"x",x
!   endif
!   call fund_hypocycloid(x,tt1,val1)
!   call fundd_hypocycloid(x,tt1,val2)
!   print *,"val1",val1,"val2",val2
!   tt2=tt1- 0.5d0*val1/val2 
!   print *,"tt1",tt1,"tt2",tt2

!   if(tt2 .gt. 1.5d0*Pi .or. tt2 .le. Pi)then
!    flag = 1
!    exit
!   endif
!   res=abs(tt1-tt2)
!   tt1=tt2
!  enddo

! elseif(quadf .eq. 4) then
!  tt1=1.75d0*Pi
!  tt2=0.0d0
!  res = tt1
!  do while(res .gt. 1.0e-10 .and. step .lt. 1000)
!   step=step+1
!   if(step .gt. 100)then
!   print *,"quad",quadf,"step",step,"res",res
!   endif
!   call fund_hypocycloid(x,tt1,val1)
!   call fundd_hypocycloid(x,tt1,val2)
!   tt2=tt1- 0.5d0*val1/val2 
!   if(tt2 .gt. 2.0d0*Pi .or. tt2 .le. 1.5d0*Pi)then
!    flag = 1
!    exit
!   endif
!   res=abs(tt1-tt2)
!   tt1=tt2
!  enddo



end subroutine dist_hypocycloid


subroutine fund_hypocycloid(x,theta,val)
! x in domain [0,1]
implicit none

real(kind=8),intent(in)  :: x(2)
real(kind=8),intent(in)  :: theta
real(kind=8)             :: val

real(kind=8)             :: c1,c2
real(kind=8)             :: xtheta,ytheta

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)

 xtheta=(0.6d0*cos(theta)+0.2d0*cos(3.0d0*theta) &
           + c1+1.0d0)/2.0d0 
 ytheta=(0.6d0*sin(theta)-0.2d0*sin(3.0d0*theta) &
           + c2+1.0d0)/2.0d0

 val=0.6d0*(xtheta-x(1))*(-sin(theta)-sin(3.0d0*theta))+ &
     0.6d0*(ytheta-x(2))*(cos(theta)-cos(3.0d0*theta))


end subroutine fund_hypocycloid


subroutine fundd_hypocycloid(x,theta,val)
implicit none

real(kind=8),intent(in)  :: x(2)
real(kind=8),intent(in)  :: theta
real(kind=8)             :: val

real(kind=8)             :: c1,c2
real(kind=8)             :: xtheta,ytheta

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)

 xtheta=(0.6d0*cos(theta)+0.2d0*cos(3.0d0*theta) &
           + c1+1.0d0)/2.0d0 
 ytheta=(0.6d0*sin(theta)-0.2d0*sin(3.0d0*theta) &
           + c2+1.0d0)/2.0d0

 val=0.5d0*(-0.6d0*sin(theta)-0.6d0*sin(3.0d0*theta))**2.0d0 + &
     (xtheta-x(1))*(-0.6d0*cos(theta)-1.8d0*cos(3.0d0*theta)) + &
     0.5d0*(0.6d0*cos(theta)-0.6d0*cos(3.0d0*theta))**2.0d0 + &
     (ytheta-x(2))*(-0.6d0*sin(theta)+1.8d0*sin(3.0d0*theta)) 


end subroutine




!---------------------------------------------------------
subroutine dist_concentric(imat,x,y,dist,probtype_in)
implicit none

integer,intent(in)               :: imat
integer,intent(in)               :: probtype_in
real(kind=8), intent(in)         :: x,y
real(kind=8), intent(out)        :: dist
real(kind=8)                     :: dist1,dist2,dist3
real(kind=8)                     :: r1,r2
real(kind=8)                     :: center(2)

if ((imat.ge.1).and.(imat.le.num_materials)) then
 ! do nothing
else
 print *,"imat invalid in dist_concentric"
 stop
endif

center(1) = 0.5d0
center(2) = 0.5d0

if(probtype_in.eq.1) then
 r1=radcen-radeps
 r2=radcen+radeps
 dist1= sqrt((x-center(1))**2.0d0+(y-center(2))**2.0d0)

 if (imat .eq. 1) then
   dist= r1 - dist1  ! inner material
 else if(imat .eq. 2) then
   if((dist1.le.r2).and.(dist1.ge.r1))then
      dist2 = r2 - dist1
      dist3 = dist1 - r1
      dist = min(dist2,dist3)
   else if (dist1.ge.r2) then
    dist=r2-dist1
   else if (dist1.le.r1) then
    dist=dist1-r1
   else
    print *,"dist1 bust"
    stop
   endif
 elseif(imat .eq. 3)then  
   dist = dist1 - r2  ! outer material
 else
   print *,"wrong number of mat"
   stop
 endif

else if(probtype_in.eq.0) then

  if(imat .eq. 1)then
    dist = 0.3d0 - y
  elseif(imat .eq. 2)then
    dist = y - 0.3d0
  else
   print *,"imat invalid imat=",imat
   print *,"x,y= ",x,y
   stop
  endif

else if(probtype_in.eq.2) then

  if(imat .eq. 1)then
    dist = 0.3d0 - x
  elseif(imat .eq. 2)then
    dist = x - 0.3d0
  else
   print *,"imat invalid imat=",imat
   print *,"x,y= ",x,y
   stop
  endif

else if (probtype_in.eq.3) then

 dist1= sqrt((x-center(1))**2.0d0+(y-center(2))**2.0d0)
 if (imat.eq.1) then
  dist=radblob-dist1
 else if (imat.eq.2) then
  dist=dist1-radblob
 else
  print *,"imat invalid"
  stop
 endif
else if (probtype_in.eq.4) then
 dist1=sqrt((x-center(1))**2.0d0+(y-center(2))**2.0d0)
 if (imat.eq.1) then
  dist=radblob-dist1
 else if (imat.eq.2) then
  dist=dist1-radblob
 else
  print *,"imat invalid"
  stop
 endif

else if (probtype_in.eq.5) then

 dist1=0.2d0-x 
 if (imat.eq.1) then
  dist=dist1
 else if (imat.eq.2) then
  dist=-dist1
 else
  print *,"imat invalid"
  stop
 endif

else if ((probtype_in.eq.19).or. &
         (probtype_in.eq.14).or. & 
         (probtype_in.eq.13).or. &
         (probtype_in.eq.15).or. &
         (probtype_in.eq.16).or. &
         (probtype_in.eq.400).or. &
         (probtype_in.eq.20)) then
 call dist_fns(imat,x,y,dist,probtype_in)

else
 print *,"probtype_in invalid6 ",probtype_in
 stop
endif


end subroutine dist_concentric

subroutine LS_sub(xsten,nhalf,dx,bfact,dist,nmat)
INTEGER_T, intent(in) :: bfact
INTEGER_T, intent(in) :: nhalf
REAL_T, intent(in) :: dx(2)
REAL_T, intent(in) :: xsten(-nhalf:nhalf,2)
INTEGER_T, intent(in) :: nmat
REAL_T, intent(out) :: dist(nmat)
INTEGER_T :: imat

 if (nmat.eq.num_materials) then
  do imat=1,nmat
   call dist_concentric(imat,xsten(0,1),xsten(0,2),dist(imat), &
          probtype_mmat_FVM)
  enddo
 else
  print *,"nmat invalid in LS_sub"
  stop
 endif

end subroutine LS_sub

!----------------------------------------------------------

!-----------------------------------------------------------------
subroutine update_volncen(nmat_in,dist,im,center,h,vol,cxtemp,cytemp)
implicit none

integer,intent(in)            :: nmat_in,im
real(kind=8),intent(in)            :: dist
real(kind=8),intent(in)       :: center(2)
real(kind=8),intent(in)       :: h
real(kind=8)                  :: vol(nmat_in)
real(kind=8)                  :: cxtemp(nmat_in),cytemp(nmat_in)

if(dist .gt. 0.0d0) then
  vol(im) = vol(im) + h*h
  cxtemp(im) = cxtemp(im) + h*h*center(1) 
  cytemp(im) = cytemp(im) + h*h*center(2)
endif

end subroutine update_volncen
!---------------------------------------------------------------------
subroutine renormalize_vf(nmat_in,vf)
implicit none

integer              :: nmat_in
real(kind=8)         :: vf(nmat_in)
real(kind=8)         :: vfsum
integer              :: i

vfsum = 0.0d0

do i = 1,nmat_in
  vfsum = vfsum + vf(i)
enddo
if (vfsum.le.0.0) then
 print *,"vfsum invalid"
 stop
endif
do  i =1,nmat_in
  vf(i) = vf(i)/vfsum
enddo


end subroutine
!-------------------------------------------------------
subroutine form_cell(vertex,cell)
implicit none

real(kind=8)             :: vertex(4,2)
type(polygon)            :: cell

if(cell%nodesnum .le. 2) then
  print *,"nodes num invalid in cell"
else
  
  cell%nodes(1)%val = vertex(1,:)
  cell%nodes(2)%val = vertex(2,:)
  cell%nodes(3)%val = vertex(4,:)
  cell%nodes(4)%val = vertex(3,:)  

  cell%SIDES(1)%PT(1)= cell%NODES(1)
  cell%SIDES(1)%PT(2)= cell%NODES(2)
  cell%SIDES(2)%PT(1)= cell%NODES(2)
  cell%SIDES(2)%PT(2)= cell%NODES(3)
  cell%SIDES(3)%PT(1)= cell%NODES(3)
  cell%SIDES(3)%PT(2)= cell%NODES(4)
  cell%SIDES(4)%PT(1)= cell%NODES(4)
  cell%SIDES(4)%PT(2)= cell%NODES(1)

endif


end subroutine form_cell
!-------------------------------------------------------
subroutine normalize_cutline(a,b,intercept)
implicit none

real(kind=8)              :: a, b, intercept
real(kind=8)              :: c

 c = sqrt(a*a + b*b)
 
 a = a/c
 b = b/c
 intercept = intercept/c

end subroutine normalize_cutline

!-------------------------------------------------------
subroutine check_sign(G,flag)
implicit none

real(kind=8)         :: g(4)
integer              :: flag
integer              :: i,j

flag = 0
do i = 1,4
  if(abs(g(i)) .lt. afrac_eps) then
     g(i)  = 0.0d0
  endif
enddo

do i = 1,4
   do j = 1,4
      if(g(i)*g(j) .lt. 0.0d0) then
        flag = 1
      endif
   enddo      
enddo     

end subroutine check_sign
!-------------------------------------------------------
subroutine AdaptQuad_sort(im,center,h,centers,dist,flag,probtype_in)
implicit none

integer                          :: flag
integer,intent(in)               :: im,probtype_in
real(kind=8),intent(in)          :: h
real(kind=8)                     :: dist
real(kind=8)                     :: center(2), centers(4,2)

  if (h.le.0.0) then
   print *,"h invalid"
   stop
  endif

  centers = 0.0d0
   flag = 0
   call dist_concentric(im,center(1),center(2),dist,probtype_in)
   if(abs(dist) .gt. h) then
      ! do nothing
      flag = 0
   else
     flag = 1  
     call AdaptQuad_sub(center, h, centers)  
   endif     


end subroutine AdaptQuad_sort


!--------------------------------------------------------
subroutine AdaptQuad_sub(center, h, centers)
implicit none

real(kind=8)             :: center(2)
real(kind=8)             :: centers(4,2)
real(kind=8),intent(in)  :: h

 if (h.le.0.0) then
  print *,"h invalid"
  stop
 endif
!    
! 3  4 
! 1  2
! 
 centers(1,1)  = center(1) - h*0.25d0
 centers(1,2)  = center(2) - h*0.25d0
 centers(2,1)  = center(1) + h*0.25d0
 centers(2,2)  = center(2) - h*0.25d0
 centers(3,1)  = center(1) - h*0.25d0
 centers(3,2)  = center(2) + h*0.25d0
 centers(4,1)  = center(1) + h*0.25d0
 centers(4,2)  = center(2) + h*0.25d0


end subroutine AdaptQuad_sub
!---------------------------------------------------------
subroutine find_vertex(center, h, vertex)
implicit none

real(kind=8)           :: center(2)
real(kind=8)           :: h
real(kind=8)           :: vertex(4,2)

 if (h.le.0.0) then
  print *,"h invalid"
  stop
 endif

!    
! 3  4 
! 1  2
! 
 vertex(1,1)  = center(1) - h*0.5d0
 vertex(1,2)  = center(2) - h*0.5d0
 vertex(2,1)  = center(1) + h*0.5d0
 vertex(2,2)  = center(2) - h*0.5d0
 vertex(3,1)  = center(1) - h*0.5d0
 vertex(3,2)  = center(2) + h*0.5d0
 vertex(4,1)  = center(1) + h*0.5d0
 vertex(4,2)  = center(2) + h*0.5d0


end subroutine find_vertex

! -----------------------------------------------------

subroutine slopecal(im,center,mx,my,probtype_in)
implicit none

integer, intent(in)                  :: probtype_in
integer                              :: im
real(kind=8),intent(in)              :: center(2)
real(kind=8)                         :: mx,my

real(kind=8)                     :: x,y,dist
real(kind=8)                     :: dist1
real(kind=8)                     :: r1,r2

 
 x = center(1)
 y = center(2)

 r1=radcen-radeps
 r2=radcen+radeps

 dist1= sqrt((x-0.5d0)**2.0d0+(y-0.5d0)**2.0d0)

! slope is grad phi/|grad phi|

if (probtype_in.eq.0) then
 mx=0.0
 if (im.eq.1) then
  my=-1.0
 else if (im.eq.2) then
  my=1.0
 else
  print *,"im invalid 1"
  stop
 endif

else if (probtype_in.eq.2) then
 my=0.0
 if (im.eq.1) then
  mx=-1.0
 else if (im.eq.2) then
  mx=1.0
 else
  print *,"im invalid 2"
  stop
 endif

else if (probtype_in.eq.3) then

 dist= radblob - dist1  ! inner distance function (material 1)
 mx = (x - 0.5d0)/(dist - radblob)
 my = (y - 0.5d0)/(dist - radblob)
 if (im.eq.1) then
  ! do nothing
 else if (im.eq.2) then
  mx=-mx
  my=-my
 else
  print *,"im invalid"
  stop
 endif

else if (probtype_in.eq.4) then
 ! normal points into the circle for material 1
 ! normal points out of the circle for material 2
 dist= radblob - dist1  ! inner distance function (material 1)
 mx = (x - 0.5d0)/(dist - radblob)
 my = (y - 0.5d0)/(dist - radblob)
 if (im.eq.1) then
  ! do nothing
 else if (im.eq.2) then
  mx=-mx
  my=-my
 else
  print *,"im invalid"
  stop
 endif

else if (probtype_in.eq.400) then
        ! call dist_concentric
else if (probtype_in.eq.5) then

 my=0.0d0
 mx=-1.0d0
 if (im.eq.1) then
  mx=-1.0d0
 else if (im.eq.2) then
  mx=1.0d0
 else
  print *,"im invalid"
  stop
 endif

else if (probtype_in.eq.1) then
 
 if (im.eq. 1) then ! inner material
   dist= r1 - dist1
   mx = (x - 0.5d0)/(dist - r1)
   my = (y - 0.5d0)/(dist - r1)
 elseif(im .eq. 2) then
   if (dist1.le.radcen) then
    mx = (x - 0.5d0)/dist1
    my = (y - 0.5d0)/dist1
   else if (dist1.ge.radcen) then
    mx = -(x - 0.5d0)/dist1
    my = -(y - 0.5d0)/dist1
   else
    print *,"dist1 is bad"
    stop
   endif 
 elseif(im .eq. 3)then
   dist = dist1 - r2  ! outer material
   mx = (x - 0.5d0)/(dist + r2)
   my = (y - 0.5d0)/(dist + r2)
 else
   print *,"wrong number of mat"
 endif
else 
 print *,"probtype_in invalid7 ",probtype_in
 stop
endif

end subroutine slopecal


!--------------------------------------------------------
!------------------------------------------------------
! Linear reconstruction normal and initial guess
!        ax+by+d=0
!  input:  G ,X, Y
!------------------------------------------------------
SUBROUTINE LS_LSF(method_flag,h,center,Gin,aout,bout,dout)
IMPLICIT NONE

! flag = 1  least square fit method   flag = 0, basic way 

integer,intent(in)                          :: method_flag
real(kind=8),intent(in)                     :: h
REAL(KIND=8)                                :: aa,bb,cc,dd,ee,ff,gg,hh
REAL(KIND=8)                                :: aout, bout, dout
INTEGER                                     :: i,j
REAL(KIND=8), EXTERNAL                      :: NORM_2d

REAL(KIND=8),INTENT(IN)                     :: CENTER(2)
REAL(KIND=8)                                :: VERTEX(4,2)
REAL(kind=8),intent(in)                     :: Gin(4)
real(kind=8)                                :: G(2,2)
real(kind=8)                                :: x(2),y(2)


G(1,1) = Gin(1) 
G(1,2) = Gin(2)
G(2,1) = Gin(3)
G(2,2) = Gin(4)

CALL find_vertex(center, h, vertex)

X(1) = VERTEX(1,1)
X(2) = VERTEX(4,1)
Y(1) = VERTEX(1,2)
Y(2) = VERTEX(4,2)


if (method_flag .eq. 1) then
 aa=0.0d0
 bb=0.0d0
 cc=0.0d0
 dd=0.0d0
 ee=0.0d0
 ff=0.0d0
 gg=0.0d0
 hh=0.0d0


 do i= 1, 2
  do j= 1, 2
     aa = aa + 2.0d0*((x(i)-center(1))**2.0d0)
     bb = bb + 2.0d0*(x(i)-center(1))*(y(j)-center(2))
     cc = cc + 2.0d0*(x(i)-center(1))
     dd = dd + 2.0d0*((y(j)-center(2))**2.0d0)
     ee = ee + 2.0d0*(y(j)-center(2))
     ff = ff + 2.0d0*(x(i)-center(1))*G(i,j)
     gg = gg + 2.0d0*(y(j)-center(2))*G(i,j)
     hh = hh + 2.0d0*G(i,j)
  enddo
 enddo

!!write(20,*) aa,bb,cc,dd,ee,ff,gg,hh

 if(bb .eq. 0.0d0  .and. ee .eq. 0.0d0) THEN
   bout = gg/dd
   aout = (ff-cc*hh)/(aa-cc*cc)
   dout = hh - cc*aout
 else
  dout = ((aa*hh-ff*cc)/(aa*ee-bb*cc)-(aa*gg-bb*ff)/(aa*dd-bb*bb))/ &
         &((aa-cc*cc)/(aa*ee-bb*cc)-(aa*ee-bb*cc)/(aa*dd-bb*bb))
  bout = (aa*gg-bb*ff)/(aa*dd-bb*bb)-((aa*ee-bb*cc)/(aa*dd-bb*bb))*dout
  aout = (ff-cc*dout-bb*bout)/aa
 endif


elseif(method_flag .eq. 2) then
 



endif


END SUBROUTINE LS_LSF
!----------------------------------------------------------------------


subroutine AdaptQuad_2d(iin,jin,nmat_in,dx,center,centroid,vf,probtype_in)
implicit none

integer,intent(in)            :: nmat_in,probtype_in
real(kind=8),intent(in)       :: dx(2)
real(kind=8),intent(in)       :: center(2)
integer     ,intent(in)       :: iin,jin

integer,parameter             :: ref_num = 5
integer                       :: im,i1,i2,i3,i4,i5,dir
integer                       :: ii
integer                       :: flag,sign_flag
real(kind=8)                  :: dist,h2,h3,h4,h5,h6
real(kind=8)                  :: x,y
real(kind=8)                  :: h
real(kind=8)                  :: center_hold(2)
real(kind=8)                  :: centers1(4,2),centers2(4,2)
real(kind=8)                  :: centers3(4,2),centers4(4,2),centers5(4,2)
real(kind=8)                  :: vertex(4,2)
real(kind=8)                  :: G(4)
real(kind=8)                  :: mx,my,c

real(kind=8)                  :: vf(nmat_in)
real(kind=8)                  :: centroid(nmat_in,2)
real(kind=8)                  :: cxtemp(nmat_in),cytemp(nmat_in)
real(kind=8)                  :: vol(nmat_in)
type(polygon)                 :: cell
type(polygon)                 :: POS_PLG,NEG_PLG
type(points)                  :: cen_temp
real(kind=8)                  :: vol_temp

! check
real(kind=8)                  :: g_center

type(points)                  :: seg(2)


x = center(1)
y = center(2)
h = dx(1)

if (h.le.0.0) then
 print *,"h invalid"
 stop
endif


 cell%nodesnum = 4 
 allocate(cell%nodes(4))
 vf = 0.0d0
 centroid = 0.0d0
 vol = 0.0d0
 cxtemp = 0.0d0
 cytemp = 0.0d0
 
do im = 1,nmat_in
   call dist_concentric(im,x,y,dist,probtype_in)

   if(abs(dist) .gt. h) then
     flag = 0
     call update_volncen(nmat_in,dist,im,center,h,vol,cxtemp,cytemp)
   else
     flag = 1   
     call AdaptQuad_sub(center, h, centers1)
   endif 
 
   if (flag.eq.0) then
    ! do nothing  
   else if(flag .eq. 1) then
    do i1 = 1,4 
     do dir=1,2
      center_hold(dir)=centers1(i1,dir)
     enddo
     h2=h/2.0d0
     call AdaptQuad_sort(im,center_hold,h2,centers2,dist,flag,probtype_in)
     if(flag .eq. 1) then
      do i2 = 1,4
       do dir=1,2
        center_hold(dir)=centers2(i2,dir)
       enddo
       h3=h/4.0d0
       call AdaptQuad_sort(im,center_hold,h3,centers3,dist,flag, &
         probtype_in)
       if(flag .eq. 1) then
        do i3 = 1,4
         do dir=1,2
          center_hold(dir)=centers3(i3,dir)
         enddo
         h4=h/8.0d0
         call AdaptQuad_sort(im,center_hold,h4,centers4,dist, &
              flag,probtype_in)
         if(flag .eq. 1) then
          do i4 = 1,4
           do dir=1,2
            center_hold(dir)=centers4(i4,dir)
           enddo
           h5=h/16.0d0
           call AdaptQuad_sort(im,center_hold,h5, centers5, &
             dist,flag,probtype_in)
           if(flag .eq. 1)then
            do i5 = 1,4
             do dir=1,2
              center_hold(dir)=centers5(i5,dir)
             enddo
             h6=h/32.0d0
             call find_vertex(center_hold, h6, vertex)

             do ii = 1,4
              call dist_concentric(im,vertex(ii,1),vertex(ii,2), &
                G(ii),probtype_in)
             enddo  ! ii
             call check_sign(G,sign_flag)

             if (1.eq.0) then
              print *,"after check_sign  sign_flag=",sign_flag
              do ii=1,4
               print *,"ii,vertex,G ",ii,vertex(ii,1),vertex(ii,2),G(ii)
              enddo
             endif

             if(sign_flag .eq. 1) then
              call slopecal(im,center_hold,mx,my,probtype_in)
               !============================================================        
               ! call LS_LSF(h/32.0d0,centers5(i5,:),G,mx,my,alphadump)
               ! ax + by +alpha = 0
              call dist_concentric(im,center_hold(1),center_hold(2), &
                 G_center,probtype_in)

  ! phi=mx x + my y + c
  ! mx x' + my y' + c = phi'
  ! c=phi'-mx x' -my y'
  !              normalize 
  !              call normalize_cutline(mx,my,alpha)
  !              c = -mx*center_hold(1)-my*center_hold(2)+alpha
              c = -mx*center_hold(1)-my*center_hold(2)+ G_center

              call form_cell(vertex,cell)
              call cell_line_inter(cell,mx,my,c,seg)

              if (1.eq.0) then
               print *,"x,y,h,im,mx,my,c ",x,y,h,im,mx,my,c
              endif

              CALL VOL_FRAC_CAL(cell,mx,my,c,POS_PLG,NEG_PLG)
                       
               ! call outputplg(pos_plg)
  
              CALL POLY_CENTROID(NEG_PLG,cen_temp)
              CALL poly_area(NEG_PLG,VOL_temp)                   
                 ! update centroid and volume
              vol(im) = vol(im) + vol_temp
              cxtemp(im) = cxtemp(im) + vol_temp*cen_temp%val(1) 
              cytemp(im) = cytemp(im) + vol_temp*cen_temp%val(2)
              CALL PLGDEL(POS_PLG)
              CALL PLGDEL(NEG_PLG)

              if (1.eq.0) then
               print *,"im,mx,my,c ",im,mx,my,c
               do ii=1,4
                print *,"ii,vertex,G ",ii,vertex(ii,1),vertex(ii,2),G(ii)
               enddo
               print *,"vfrac ",vol_temp/(h6*h6)
               print *,"cenx,ceny ",cen_temp%val(1),cen_temp%val(2)
              endif

             elseif(sign_flag .eq. 0) then 
              CALL dist_concentric(im,center_hold(1),center_hold(2), &
                 dist,probtype_in)
              CALL update_volncen(nmat_in,dist,im,center_hold, &
                 h/32.0d0,vol,cxtemp,cytemp)
             else
              print *,"sign_flag = ", sign_flag
              stop
             endif
            enddo  ! i5
           elseif(flag .eq. 0) then
            call update_volncen(nmat_in,dist,im,center_hold, &
              h/16.0d0,vol,cxtemp,cytemp)  
           else
            print *,"flag error in AdaptQuad_2d"
            stop
           endif
          enddo ! i4
         elseif(flag .eq. 0) then
          call update_volncen(nmat_in,dist,im,center_hold,h/8.0d0, &
            vol,cxtemp,cytemp)  
         else
          print *,"flag error in AdaptQuad_2d"
          stop
         endif
        enddo  ! i3
       elseif(flag .eq. 0) then
        call update_volncen(nmat_in,dist,im,center_hold,h/4.0d0, &
          vol,cxtemp,cytemp)  
       else
        print *,"flag error in AdaptQuad_2d" 
        stop
       endif
      enddo ! i2
     elseif(flag .eq. 0) then
      call update_volncen(nmat_in,dist,im,center_hold,h/2.0d0,vol, &
         cxtemp,cytemp)  
     else
       print *,"flag error in AdaptQuad_2d"    
     endif
    enddo ! i1

   else 
    print *,"flag invalid"
    stop
   endif

enddo ! im

if (h.le.0.0) then
 print *,"h invalid"
 stop
endif


do im = 1,nmat_in
   if(vol(im) .gt. afrac_eps)then
    centroid(im,1) = cxtemp(im)/vol(im)
    centroid(im,2) = cytemp(im)/vol(im)
   else
     vol(im) = 0.0d0
     centroid(im,:) = 0.0d0
   endif
   vf(im) = vol(im)/(h*h)

enddo


end subroutine AdaptQuad_2d
! -----------------------------------------------
!---------------------------------------------------------------------
subroutine vf_correct(iin,jin,nmat_in,vf,probtype_in)
implicit none

integer       ,intent(in)  :: nmat_in,iin,jin,probtype_in
real(kind=8)               :: vf(nmat_in)
real(kind=8)               :: vcheck
integer                    :: i

if ((probtype_in.eq.0).or. &
    (probtype_in.eq.2).or. &
    (probtype_in.eq.3).or. &
    (probtype_in.eq.4).or. &
    (probtype_in.eq.400).or. &
    (probtype_in.eq.5).or. &
    (probtype_in.eq.14).or. &
    (probtype_in.eq.15)) then

 if (nmat_in.ne.2) then
  print *,"nmat_in invalid"
  stop
 endif
 
 vcheck = vf(1) + vf(2)
 if(vf(1) .lt. 0.0d0 .or. vf(2) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+FACETOL_SANITY) then
  print *,"goes into vf_correct"
  stop
 endif

else if ((probtype_in.eq.1).or. &
         (probtype_in.eq.13).or. &
         (probtype_in.eq.19)) then
 if (nmat_in.ne.3) then
  print *,"nmat_in invalid"
  stop
 endif
 
 vcheck = vf(1) + vf(3)
 if(vf(1) .lt. 0.0d0 .or. vf(3) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+FACETOL_SANITY) then
  print *,"goes into vf_correct"
  if(vf(1) .gt. 1.0d0) then
    print *,"case1"
     vf(1) = 1.0d0
     vf(3) = 0.0d0
  elseif(vf(3) .gt. 1.0d0)then
    print *,"case2"
     vf(1) = 0.0d0
     vf(3) = 1.0d0
  else
    print *,"case3",iin,jin
    print *,"overshoot", vf(1) , vf(3)
    if(vf(1) .gt. vf(3))then
       vf(1) = 1.0d0
       vf(3) = 0.0d0
    elseif(vf(1) .lt. vf(3))then
       vf(1) = 0.0d0
       vf(3) = 1.0d0
    endif
  endif
 endif

elseif(probtype_in .eq. 16)then
 vcheck = vf(1) + vf(2)
 if(vf(1) .lt. 0.0d0 .or. vf(2) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+FACETOL_SANITY) then
  print *,"goes into vf_correct3",vcheck
  if(vf(1) .gt. 1.0d0) then
    print *,"case1"
     vf(1) = 1.0d0
     vf(2) = 0.0d0
  elseif(vf(2) .gt. 1.0d0)then
    print *,"case2"
     vf(1) = 0.0d0
     vf(2) = 1.0d0
  else
    print *,"case3",iin,jin
    print *,"overshoot", vf(1) , vf(2)
    if(vf(1) .gt. vf(2))then
       vf(1) = 1.0d0
       vf(2) = 0.0d0
    elseif(vf(1) .lt. vf(2))then
       vf(1) = 0.0d0
       vf(2) = 1.0d0
    endif
  endif
 endif
elseif(probtype_in .eq. 17)then
! vcheck = vf(1) + vf(2)+vf(3)+vf(4)+vf(5)
!  if(vcheck .gt. 1.0d0+FACETOL_SANITY) then
!  print *,"goes into vf_correct4",vcheck,vf
! endif
 do i=1,5
  if(vf(i) .lt. 0.0d0)then
   print *,"vf is negtive"
   stop
  endif
 enddo

 if(vf(1) .gt. 1.0d0)then
  vf(1)=1.0d0
  vf(2)=0.0d0
  vf(3)=0.0d0
  vf(4)=0.0d0
  vf(5)=0.0d0
 elseif(vf(2) .gt. 1.0d0)then
  vf(1)=0.0d0
  vf(2)=1.0d0
  vf(3)=0.0d0
  vf(4)=0.0d0
  vf(5)=0.0d0
 elseif(vf(3) .gt. 1.0d0)then
  vf(1)=0.0d0
  vf(2)=0.0d0
  vf(3)=1.0d0
  vf(4)=0.0d0
  vf(5)=0.0d0
 elseif(vf(4) .gt. 1.0d0)then
  vf(1)=0.0d0
  vf(2)=0.0d0
  vf(3)=0.0d0
  vf(4)=1.0d0
  vf(5)=0.0d0
 elseif(vf(5) .gt. 1.0d0)then
  vf(1)=0.0d0
  vf(2)=0.0d0
  vf(3)=0.0d0
  vf(4)=0.0d0
  vf(5)=1.0d0
 else
  ! do nothing
   
 endif
elseif(probtype_in .eq. 20)then
! vcheck = vf(1) + vf(2)+vf(3)+vf(4)+vf(5)
!  if(vcheck .gt. 1.0d0+FACETOL_SANITY) then
!  print *,"goes into vf_correct4",vcheck,vf
! endif

 if (nmat_in.ne.6) then
  print *,"nmat_in invalid"
  stop
 endif
 
 vcheck = vf(1) + vf(6)
 if(vf(1) .lt. 0.0d0 .or. vf(6) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+FACETOL_SANITY) then
  print *,"goes into vf_correct2",vcheck
  if(vf(1) .gt. 1.0d0) then
    print *,"case1"
     vf(1) = 1.0d0
     vf(6) = 0.0d0
  elseif(vf(6) .gt. 1.0d0)then
    print *,"case2"
     vf(1) = 0.0d0
     vf(6) = 1.0d0
  else
    print *,"case3",iin,jin
    print *,"overshoot", vf(1) , vf(3)
    if(vf(1) .gt. vf(6))then
       vf(1) = 1.0d0
       vf(6) = 0.0d0
    elseif(vf(1) .lt. vf(6))then
       vf(1) = 0.0d0
       vf(6) = 1.0d0
    endif
  endif
 endif


else
 print *,"probtype_in invalid10 ",probtype_in
 stop
endif

end subroutine vf_correct



!--------------------------------------------
SUBROUTINE init_vfncen(N,cell,nmat_in,dx,centroid,vf,probtype_in)
use stackvolume_module

IMPLICIT NONE

integer,intent(in)     :: N,nmat_in,probtype_in
type(polygon),intent(in)  :: cell(-1:N,-1:N)
real(kind=8),intent(in)  :: dx(2)
TYPE(POINTS),DIMENSION(-1:N,-1:N,nmat_in)  :: CENTROID
real(kind=8)                       :: center(2)
real(kind=8)                       :: cen_temp(nmat_in,2)
integer                            :: i,j,im,dir
real(kind=8)                       :: vf_temp(nmat_in)
real(kind=8)                       :: voldata(nmat_in,6)
real(kind=8)                       :: vofdark(nmat_in)
real(kind=8)                       :: voflight(nmat_in)
real(kind=8)                       :: cendark(nmat_in,2)
real(kind=8)                       :: cenlight(nmat_in,2)
real(kind=8)                       :: xsten(-3:3,2)
real(kind=8)                       :: vf(-1:N,-1:N,nmat_in)
integer                            :: stack_max_level
integer                            :: current_level
integer                            :: bfact
integer                            :: iloc,nhalf

 print *,"in init_vfncen "
 print *,"N,nmat_in ",N,nmat_in
 print *,"dx ",dx(1),dx(2) 

 probtype_mmat_FVM=probtype_in

 cen_temp = 0.0d0 ! (nmat_in,2)
 vf_temp = 0.0d0  ! (nmat_in)
 vf = 0.0d0
 do i = -1,N
  do j = -1,N
   do im= 1,nmat_in  
      centroid(i,j,im)%val = 0.0d0
   enddo
  enddo
 enddo

 do i = -1,N
   do j = -1,N
    do dir=1,2
     center(dir) = cell(i,j)%center%val(dir)
    enddo
     ! in: multimat_FVM.F90
    if (1.eq.0) then
     call AdaptQuad_2d(i,j,nmat_in,dx,center,cen_temp,vf_temp,probtype_in)
    else if (1.eq.1) then
     stack_max_level=6  ! was 4
     current_level=0
     nhalf=3
     do iloc=-nhalf,nhalf
      do dir=1,2
       xsten(iloc,dir)=center(dir)+0.5d0*dx(dir)*iloc
      enddo
     enddo 
     bfact=1
     call stackvolume_batch(xsten,nhalf,dx,bfact, &
       voldata,nmat_in,current_level,stack_max_level,LS_sub)
      ! centroid relative to the center (or centroid?) of the cell
     call extract_vof_cen_batch(voldata,vofdark,voflight, &
       cendark,cenlight,nmat_in)
     do im = 1, nmat_in
      vf_temp(im)=vofdark(im)
      do dir=1,2
       cen_temp(im,dir)=cendark(im,dir)+xsten(0,dir)
      enddo
     enddo
    else
     print *,"bust"
     stop
    endif

    call vf_correct(i,j,nmat_in,vf_temp,probtype_in)

if (probtype_in.eq.0) then
 ! do nothing
else if (probtype_in.eq.2) then
 ! do nothing
else if (probtype_in.eq.3) then
 ! do nothing
else if (probtype_in.eq.4) then
 ! do nothing
else if (probtype_in.eq.400) then
 ! do nothing
else if (probtype_in.eq.5) then
 ! do nothing
elseif(probtype_in .eq. 15)then   
 ! vf_temp(2)=1.0-vf_temp(1)

else if ((probtype_in.eq.1).or. &
         (probtype_in.eq.13).or. &
         (probtype_in.eq.19)) then
 vf_temp(2) = 1.0d0 - vf_temp(1) -vf_temp(3)

 if(abs(vf_temp(2)) .lt. FACETOL_DVOL)then
  vf_temp(2) = 0.0d0
 endif

 if(vf_temp(2) .lt. -FACETOL_SANITY) then
  print *,"vf_temp(2) is negative"
 endif

elseif(probtype_in .eq. 16)then      
 vf_temp(3) = 1.0d0 - vf_temp(1) -vf_temp(2)

 if(abs(vf_temp(3)) .lt. FACETOL_DVOL)then
  vf_temp(3) = 0.0d0
 endif
 if(vf_temp(3) .lt. -FACETOL_SANITY) then
  print *,"vf_temp(3) is negative"
 endif

elseif(probtype_in .eq. 17)then 

 if(center(1) .gt. 0.5d0 .and. center(2) .gt. 0.5d0)then
  vf_temp(2) = 1.0d0-vf_temp(1)
 elseif(center(1) .lt. 0.5d0 .and. center(2) .gt. 0.5d0)then
  vf_temp(3) = 1.0d0-vf_temp(1)  
 elseif(center(1) .lt. 0.5d0 .and. center(2) .lt. 0.5d0)then
  vf_temp(4) = 1.0d0-vf_temp(1)  
 elseif(center(1) .gt. 0.5d0 .and. center(2) .lt. 0.5d0)then
  vf_temp(5) = 1.0d0-vf_temp(1)  
 else
  print *,"center is not aligned with grid"
  stop
 endif

! do nothing
 if(1 .eq. 0)then
 if(vf_temp(1) .ne. 0.0d0 .and. vf_temp(2) .ne. 0.0d0 .and. &
    vf_temp(3) .ne. 0.0d0 .and. vf_temp(4) .eq. 0.0d0 .and. &
    vf_temp(5) .eq. 0.0d0)then
   vf_temp(3)=1.0d0-vf_temp(1)-vf_temp(2)
   do dir=1,2
    cen_temp(3,dir) =  &
     (center(dir) - vf_temp(1)*cen_temp(1,dir) - &
      vf_temp(2)*cen_temp(2,dir))/vf_temp(3)
   enddo
 elseif(vf_temp(1) .ne. 0.0d0 .and. vf_temp(2) .ne. 0.0d0 .and. &
        vf_temp(5) .ne. 0.0d0 .and. vf_temp(4) .eq. 0.0d0 .and. &
        vf_temp(3) .eq. 0.0d0)then
   vf_temp(5)=1.0d0-vf_temp(1)-vf_temp(2)
   do dir=1,2
    cen_temp(5,dir) =  &
     (center(dir) - vf_temp(1)*cen_temp(1,dir) - &
      vf_temp(2)*cen_temp(2,dir))/vf_temp(5)
   enddo
 elseif(vf_temp(1) .ne. 0.0d0 .and. vf_temp(4) .ne. 0.0d0 .and. &
        vf_temp(3) .ne. 0.0d0 .and. vf_temp(2) .eq. 0.0d0 .and. &
        vf_temp(5) .eq. 0.0d0)then
   vf_temp(4)=1.0d0-vf_temp(1)-vf_temp(3)
   do dir=1,2
    cen_temp(4,dir) =  &
     (center(dir) - vf_temp(1)*cen_temp(1,dir) - &
      vf_temp(3)*cen_temp(3,dir))/vf_temp(4)
   enddo
 elseif(vf_temp(1) .ne. 0.0d0 .and. vf_temp(4) .ne. 0.0d0 .and. &
        vf_temp(5) .ne. 0.0d0 .and. vf_temp(2) .eq. 0.0d0 .and. &
        vf_temp(3) .eq. 0.0d0)then
   vf_temp(4)=1.0d0-vf_temp(1)-vf_temp(5)
   do dir=1,2
    cen_temp(4,dir) =  &
     (center(dir) - vf_temp(1)*cen_temp(1,dir) - &
      vf_temp(5)*cen_temp(5,dir))/vf_temp(4)
   enddo
 else
   ! do nothing
 endif
 endif

elseif(probtype_in .eq. 20)then                                     

 if(center(1) .gt. 0.5d0 .and. center(2) .gt. 0.5d0)then
  vf_temp(2) = 1.0d0-vf_temp(1)-vf_temp(6)
 elseif(center(1) .lt. 0.5d0 .and. center(2) .gt. 0.5d0)then
  vf_temp(3) = 1.0d0-vf_temp(1)-vf_temp(6)  
 elseif(center(1) .lt. 0.5d0 .and. center(2) .lt. 0.5d0)then
  vf_temp(4) = 1.0d0-vf_temp(1)-vf_temp(6)  
 elseif(center(1) .gt. 0.5d0 .and. center(2) .lt. 0.5d0)then
  vf_temp(5) = 1.0d0-vf_temp(1)-vf_temp(6)  
 else
  print *,"center is not aligned with grid"
  stop
 endif

else
 print *,"probtype_in invalid8 ",probtype_in
 stop
endif

! sanity check
do im = 1,nmat_in

   if (abs(vf_temp(im)).le.FACETOL_DVOL) then
    vf_temp(im)=0.0
   else if (abs(vf_temp(im)-1.0).le.FACETOL_DVOL) then
    vf_temp(im)=1.0
   else if ((vf_temp(im).gt.0.0).and.(vf_temp(im).lt.1.0)) then
    ! do nothing
   else
    print *,im,vf_temp(im)
    print *, "Error in vf_temp 2" 
    stop
   endif

enddo ! im


if (probtype_in.eq.0) then
 ! do nothing
else if (probtype_in.eq.2) then
 ! do nothing
else if (probtype_in.eq.3) then
 ! do nothing
else if (probtype_in.eq.4) then
 ! do nothing
else if (probtype_in.eq.400) then
 ! do nothing
else if (probtype_in.eq.5) then
 ! do nothing
else if (probtype_in.eq.14) then
 ! do nothing
elseif(probtype_in .eq. 15) then
 !  do dir=1,2
 !   cen_temp(2,dir) =  &
 !   (center(dir) - vf_temp(1)*cen_temp(1,dir))/vf_temp(2)
 !  enddo  



else if ((probtype_in.eq.1).or. &
         (probtype_in.eq.13).or. &
         (probtype_in.eq.19)) then
 if(vf_temp(2) .gt. FACETOL_DVOL)then
  do dir=1,2
   cen_temp(2,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir) - vf_temp(3)*cen_temp(3,dir))/vf_temp(2)
  enddo
 endif

elseif(probtype_in.eq.16)then
 if(vf_temp(3) .gt. FACETOL_DVOL)then
  do dir=1,2
   cen_temp(3,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir) - vf_temp(2)*cen_temp(2,dir))/vf_temp(3)
  enddo
 endif

else if (probtype_in.eq.17) then
 if(vf_temp(2) .gt. FACETOL_DVOL .and. vf_temp(2) .lt. 1.0d0)then
  do dir=1,2
   cen_temp(2,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir))/vf_temp(2)
  enddo
 elseif(vf_temp(3) .gt. FACETOL_DVOL .and. vf_temp(3) .lt. 1.0d0)then
  do dir=1,2
   cen_temp(3,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir))/vf_temp(3)
  enddo
 elseif(vf_temp(4) .gt. FACETOL_DVOL .and. vf_temp(4) .lt. 1.0d0)then
  do dir=1,2
   cen_temp(4,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir))/vf_temp(4)
  enddo
 elseif(vf_temp(5) .gt. FACETOL_DVOL .and. vf_temp(5) .lt. 1.0d0)then
  do dir=1,2
   cen_temp(5,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir))/vf_temp(5)
  enddo 

 endif

else if (probtype_in.eq.20) then
 if(vf_temp(2) .gt. FACETOL_DVOL .and. vf_temp(2) .lt. 1.0d0)then
   do dir=1,2
   cen_temp(2,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir)-vf_temp(6)*cen_temp(6,dir))/vf_temp(2)
  enddo
 elseif(vf_temp(3) .gt. FACETOL_DVOL .and. vf_temp(3) .lt. 1.0d0)then
  do dir=1,2
   cen_temp(3,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir)-vf_temp(6)*cen_temp(6,dir))/vf_temp(3)
  enddo
 elseif(vf_temp(4) .gt. FACETOL_DVOL .and. vf_temp(4) .lt. 1.0d0)then
  do dir=1,2
   cen_temp(4,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir)-vf_temp(6)*cen_temp(6,dir))/vf_temp(4)
  enddo
 elseif(vf_temp(5) .gt. FACETOL_DVOL .and. vf_temp(5) .lt. 1.0d0)then
  do dir=1,2
   cen_temp(5,dir) =  &
    (center(dir) - vf_temp(1)*cen_temp(1,dir)-vf_temp(6)*cen_temp(6,dir))/vf_temp(5)
  enddo 

 endif

else
 print *,"probtype_in invalid9 ",probtype_in
 stop
endif

!call renormalize_vf_temp(nmat_in,vf_temp) 


    do im = 1, nmat_in
     do dir=1,2
      centroid(i,j,im)%val(dir) = cen_temp(im,dir)
     enddo
     vf(i,j,im) = vf_temp(im)

     if (1.eq.0) then
      print *,"i,j,im,vf,cenx,ceny ",i,j,im,vf_temp(im), &
        cen_temp(im,1),cen_temp(im,2)
     endif  
    enddo ! im

   enddo
 enddo

END SUBROUTINE init_vfncen
!---------------------------------------------------------------------------
!--------------------------------------------------------------
subroutine init_mofdata(N,sdim,dx,nmat_in,nten,CELL,vf,CENTROID,mofdata)
implicit none

integer,intent(in)               :: N,sdim,nmat_in
TYPE(POLYGON)                    :: CELL(-1:N,-1:N)
real(kind=8),dimension(-1:N,-1:N,nmat_in) :: vf
real(kind=8),intent(in)          :: dx(sdim)
TYPE(POINTS),DIMENSION(-1:N,-1:N,nmat_in) :: CENTROID

real(kind=8)                     :: mofdata(-1:N,-1:N,(sdim*2+3)*nmat_in)
real(kind=8)                     :: local_mofdata((sdim*2+3)*nmat_in)

integer                          :: i,j
integer                          :: im,dir 
integer                          :: ii,jj
integer                          :: bfact
integer                          :: mof_verbose
integer                          :: continuous_mof
integer                          :: caller_id
integer,intent(in)               :: nten
integer                          :: nhalf0
real(kind=8)                     :: xsten0(-3:3,sdim)
integer                          :: use_ls_data
REAL(kind=8)                     :: LS_stencil(D_DECL(-1:1,-1:1,-1:1),nmat_in)
real(kind=8)                     :: xtetlist_vof(sdim+1,sdim,2000)
real(kind=8)                     :: xtetlist_cen(sdim+1,sdim,2000)
real(kind=8)                     :: multi_centroidA(nmat_in,sdim)
integer                          :: ngeom_recon_in,vofcomp

print *,"in init_mofdata"
print *,"N,sdim,nmat_in,nten = ",N,sdim,nmat_in,nten
print *,"dx= ",dx(1),dx(2)

bfact = 1
nhalf0 = 3

mof_verbose = 0
LS_stencil = 0.0d0
xtetlist_vof = 0.0d0
xtetlist_cen = 0.0d0
ngeom_recon_in = 2*sdim+3

mofdata = 0.0d0

do  i = 0,N-1
 do j = 0,N-1

   do ii=-3,3
    do dir=1,sdim
     xsten0(ii,dir)=ii*0.5d0*dx(dir)+CELL(i,j)%center%val(dir)
    enddo
   enddo

   do im = 1,nmat_in
    if (abs(vf(i,j,im)).le.FACETOL_DVOL) then
     vf(i,j,im)=0.0
    else if (abs(vf(i,j,im)-1.0).le.FACETOL_DVOL) then
     vf(i,j,im)=1.0
    else if ((vf(i,j,im).gt.0.0).and.(vf(i,j,im).lt.1.0)) then
     ! do nothing
    else
     print *,"vf invalid vf=",vf(i,j,im)
     print *,"i,j,im,N,nmat_in,sdim,ngeom_recon_in ",i,j,im,N,nmat_in,sdim, &
      ngeom_recon_in
     stop
    endif

    vofcomp=(im-1)*ngeom_recon_in+1 
    mofdata(i,j,vofcomp) = vf(i,j,im)
    
    if(vf(i,j,im) .gt. 0.0 .and. vf(i,j,im) .lt. 1.0d0)then

     do dir = 1,sdim
      mofdata(i,j,dir+vofcomp) = &
       CENTROID(i,j,im)%val(dir)-CELL(i,j)%center%val(dir)
     enddo

    else if ((vf(i,j,im).eq.0.0).or.(vf(i,j,im).eq.1.0)) then
     ! do nothing
    else
     print *,"invalid volume fraction"
     stop   
    endif
    
    mofdata(i,j,3+vofcomp) = 0.0d0  
    mofdata(i,j,4+vofcomp) = 0.0d0  
    mofdata(i,j,5+vofcomp) = 0.0d0  
    mofdata(i,j,6+vofcomp) = 0.0d0  
   enddo ! im
   
!        print *, "i=", i, "j=", j
!        print *, mofdata(i,j,:)

        ! bfact = 1   xsten0 =     nhalf0 = 3
    continuous_mof=0  ! standard MOF
    caller_id=11
    use_ls_data=1
    do ii=-1,1
    do jj=-1,1
    do im = 1,nmat_in
     call dist_concentric(im,xsten0(2*ii,1),xsten0(2*jj,2), &
       LS_stencil(ii,jj,im),probtype_mmat_FVM)
    enddo
    enddo
    enddo
    do im=1,(sdim*2+3)*nmat_in
     local_mofdata(im)=mofdata(i,j,im)
    enddo
    call multimaterial_MOF( &
      bfact,dx,xsten0,nhalf0, &
      mof_verbose, &
      use_ls_data, & 
      LS_stencil, &
      xtetlist_vof, &
      xtetlist_cen, &
      1000, &  ! nlist_alloc=1000
      1000, &  ! nmax=1000
      local_mofdata, &
      multi_centroidA, &
      continuous_mof, &
      nmat_in, &
      sdim, &
      caller_id)
    do im=1,(sdim*2+3)*nmat_in
     mofdata(i,j,im)=local_mofdata(im)
    enddo

 enddo
enddo



end subroutine init_mofdata

!-------------------------------------------------------------------


subroutine get_filament_source(x_in,t_in,probtype_in,im,sdim,G_in)
IMPLICIT NONE

integer, intent(in) :: probtype_in
integer, intent(in) :: im
integer, intent(in) :: sdim
REAL*8, intent(in) :: x_in(sdim)
REAL*8, intent(in) :: t_in
REAL*8, intent(out) :: G_in
REAL*8 radius_in,theta_in
REAL*8 r1,r2
REAL*8 delx,dely
REAL*8 refc
REAL*8 mypi
real(kind=8) :: radial_slope

 mypi=4.0d0*atan(1.0d0)
 if (sdim.ne.2) then
  print *,"sdim invalid"
  stop
 endif
 if (t_in.ge.0.0d0) then
  ! do nothing
 else
  print *,"t invalid"
  stop
 endif

 if (probtype_in.eq.0) then
  if (im.eq.1) then
   G_in=0.0
  else if (im.eq.2) then
   G_in=0.0
  else
   print *,"im invalid 3"
   print *,"im=",im
   print *,"sdim=",sdim
   print *,"x,y,t ",x_in(1),x_in(2),t_in
   print *,"probtype_in ",probtype_in
   stop
  endif
 else if (probtype_in.eq.2) then
  if (im.eq.1) then
   G_in=0.0
  else if (im.eq.2) then
   G_in=0.0
  else
   print *,"im invalid 4"
   stop
  endif
 else if (probtype_in.eq.3) then
  if (im.eq.1) then
   G_in=0.0
  else if (im.eq.2) then
   G_in=0.0
  else
   print *,"im invalid 4"
   stop
  endif
 else if (probtype_in.eq.4) then
  if (im.eq.1) then
   G_in=0.0
  else if (im.eq.2) then
   G_in=0.0
  else
   print *,"im invalid 4"
   stop
  endif
 else if (probtype_in.eq.400) then
  if (im.eq.1) then
   G_in=0.0
  else if (im.eq.2) then
   G_in=0.0
  else
   print *,"im invalid 4"
   stop
  endif
 else if (probtype_in.eq.5) then
  if (im.eq.1) then
   G_in=0.0
  else if (im.eq.2) then
   G_in=0.0
  else
   print *,"im invalid 4"
   stop
  endif
 else if (probtype_in.eq.1) then
  if ((im.eq.1).or.(im.eq.3)) then
   G_in=0.0
  else if (im.eq.2) then

   r1=radcen-radeps
   r2=radcen+radeps
   if (r2.gt.r1) then
    ! do nothing
   else
    print *,"r2 or r1 invalid"
    stop
   endif

   delx=x_in(1)-0.5d0
   dely=x_in(2)-0.5d0

   radius_in = sqrt(delx**2.0d0 +dely**2.0d0)

    ! x=r cos(theta)
    ! y=r sin(theta)
   if (radius_in.le.radeps/1000.0) then
    theta_in=0.0
   else if ((delx.ge.0.0).and.(dely.ge.0.0)) then
    theta_in=acos(delx/radius_in)
   else if ((delx.le.0.0).and.(dely.ge.0.0)) then
    theta_in=acos(abs(delx)/radius_in)
    theta_in=mypi-theta_in
   else if ((delx.le.0.0).and.(dely.le.0.0)) then
    theta_in=acos(abs(delx)/radius_in)
    theta_in=mypi+theta_in
   else if ((delx.ge.0.0).and.(dely.le.0.0)) then
    theta_in=acos(delx/radius_in)
    theta_in=2.0d0*mypi-theta_in
   else
    print *,"delx or dely invalid"
    stop
   endif

   radial_slope=1.0d0/(r2-r1)
   if (radial_variation.eq.1) then
    ! do nothing
   else if (radial_variation.eq.0) then
    radial_slope=0.0d0
   else
    print *,"radial_variation invalid"
    stop
   endif

   ! T=2+slope*(r-radcen)+sin(theta) exp(-t/(rc^2))  alpha=1
   ! T_t - (T_rr + T_r/r + T_theta theta/r^2)=
   ! exp(-t/rc^2)(-sin(theta)/rc^2+sin(theta)/r^2)
   G_in=exp(-t_in/(radcen**2))*sin(theta_in)*(-1.0/(radcen**2)+ &
        1.0/(radius_in**2))-radial_slope/radius_in
  else
   print *,"im invalid 5"
   stop
  endif

 else if (probtype_in.eq.13) then

   G_in = 0.0d0

   ! declared in vof_cisl.F90 (set in main.F90)
   if (dirichlet_pentafoil.eq.1) then 

    if (im.eq.2) then

     if (1.eq.1) then

      ! T_{t}-(T_{xx}+T_{yy})
      refc= (0.02d0*sqrt(5.0d0)+1.0d0)/2.0d0 ! refc=0.0d0 default
      delx=x_in(1)-refc
      dely=x_in(2)-refc
      radius_in=sqrt(delx**2.0d0+dely**2.0d0)
      G_in=exp(-t_in)*(-(radius_in**2)-4.0d0)

     else if (1.eq.0) then

      delx=x_in(1)-0.5d0-0.02d0*sqrt(5.0)
      dely=x_in(2)-0.5d0-0.02d0*sqrt(5.0)  
      radius_in = sqrt(delx**2.0d0 +dely**2.0d0)

      ! x=r cos(theta)
      ! y=r sin(theta)
      if (radius_in.le.radeps/1000.0) then
       theta_in=0.0
      else if ((delx.ge.0.0).and.(dely.ge.0.0)) then
       theta_in=acos(delx/radius_in)
      else if ((delx.le.0.0).and.(dely.ge.0.0)) then
       theta_in=acos(abs(delx)/radius_in)
       theta_in=mypi-theta_in
      else if ((delx.le.0.0).and.(dely.le.0.0)) then
       theta_in=acos(abs(delx)/radius_in)
       theta_in=mypi+theta_in
      else if ((delx.ge.0.0).and.(dely.le.0.0)) then
       theta_in=acos(delx/radius_in)
       theta_in=2.0d0*mypi-theta_in
      else
       print *,"delx or dely invalid"
       stop
      endif
      ! T=2+sin(theta) exp(-t)  alpha=1
      ! T_t - (T_rr + T_r/r + T_theta theta/r^2)=
      ! exp(-t/rc^2)(-sin(theta)/rc^2+sin(theta)/r^2)
      G_in=exp(-t_in)*sin(theta_in)*(-1.0d0+1.0d0/(radius_in**2))

     else
      print *,"bust"
      stop
     endif

    endif
   else if (dirichlet_pentafoil.eq.0) then
    ! do nothing
   else
    print *,"dirichlet_pentafoil invalid"
    stop
   endif

 elseif(probtype_in.eq.14 )then
   G_in = (-4.0d0 -(x_in(1)**2.0d0+x_in(2)**2.0d0))*exp(-t_in)

 elseif(probtype_in.eq.15)then
   G_in = 0.0d0
 elseif(probtype_in.eq.16)then
   G_in = 0.0d0
 elseif(probtype_in.eq.17)then
   G_in = 0.0d0
 elseif(probtype_in.eq.19)then
   G_in = 0.0d0
 elseif(probtype_in.eq.20)then
   G_in = 0.0d0
 else
  print *,"probtype_in invalid11 ",probtype_in
  stop
 endif

return
end subroutine get_filament_source



subroutine set_polar_2d(sdim,N,M,kappa,tau,r,z,dr,dz,u,ptau,cycling_step)
implicit none


integer,intent(in)          :: sdim                 
integer,intent(in)          :: N      ! discretization in r direction
integer,intent(in)          :: M       ! discretization in theta direction
real(kind=8),intent(in)     :: tau
!real(kind=8),parameter      :: rlo=radcen-radeps
!real(kind=8),parameter      :: rhi=radcen+radeps

real(kind=8)                :: r(0:N)
real(kind=8)                :: z(0:M)

real(kind=8)                :: dr,dz
real(kind=8)                :: ptau
real(kind=8)                :: kappa

real(kind=8)                :: u(0:N,0:M)

integer                     :: i,j
integer                     :: cycling_step


if(sdim .ne. 2)then
 print *,"invalid dimension 2909"
 stop
endif
do i=0,N
do j=0,M
 u(i,j)=0.0d0
enddo
enddo

dr=(rhi-rlo)/N
dz=(2*Pi)/M

do i=0,N
 r(i)=rlo+i*dr
enddo
do i=0,M
 z(i)=i*dz
enddo

!ptau=0.5d0*kappa*min(((dr)**2.0d0), &
!       minval(r)*dr, &
!      (minval(r)**2.0d0)*((dz)**2.0d0))

ptau=0.5d0/(kappa*2.0d0*(1.0d0/(dr*dr)+1.0d0/(dz*dz*maxval(r)*maxval(r))))


if(ptau .gt. tau)then
 print *,"ptau is greater than tau"
 stop
endif

cycling_step = floor(tau/ptau) +1

ptau= tau/real(cycling_step,8)


!do i=0,N                                                              
! do j=0,M                                                             
!  u(i,j)=2.0d0                                                      
!  u(i,j) = 1.0d0 + 10.0d0*(r(i)-rlo)                                 
! enddo                                                               
!enddo                                                                

do i=0,N
 do j=0,M
  u(i,j)=BC_T1*(rhi-r(i))/(rhi-rlo) + &
         BC_T2*(r(i)-rlo)/(rhi-rlo) + &
         100.0d0*sin(z(j))*(r(i)-rlo)*(rhi-r(i))
 enddo
enddo
                                                                
!do j=0,M                                                             
! u(0,j)=T1                                                        
! u(N,j)=3.0d0                                                         
!enddo


end subroutine set_polar_2d

!----------------------------------------------------
subroutine polar_2d_heat(sdim,N,M,kappa,tau, r,z,dr,dz,u)
! polar coordinate 2-d heat equation solver.
! forward in time, central FD discretization
implicit none

integer,intent(in)          :: sdim                 
integer,intent(in)          :: N      ! discretization in r direction
integer,intent(in)          :: M       ! discretization in theta direction

real(kind=8),intent(in)     :: r(0:N)
real(kind=8),intent(in)     :: z(0:M)
real(kind=8)                :: u(0:N,0:M),u_new(0:N,0:M)
real(kind=8)                :: dr,dz
real(kind=8),intent(in)     :: tau
real(kind=8),intent(in)     :: kappa


integer                     :: i,j


 do i=1,N-1
  do j=1,M-1
   u_new(i,j)= (1.0d0+kappa*tau*(-2.0d0/(dr**2.0d0)- &
                       2.0d0/((r(i)**2.0d0)*(dz**2.0d0))))*u(i,j) &
                + kappa*tau*(  & 
                 (1.0d0/(dr**2.0d0)+1.0d0/(r(i)*2.0d0*dr))*u(i+1,j) &
                + (1.0d0/(dr**2.0d0)-1.0d0/(r(i)*2.0d0*dr))*u(i-1,j) &
                + 1.0d0/((r(i)**2.0d0)*(dz**2.0d0))*u(i,j+1) &
                + 1.0d0/((r(i)**2.0d0)*(dz**2.0d0))*u(i,j-1)) 
               
  enddo
 enddo
  
 do i=1,N-1
   u_new(i,0)= (1.0d0+kappa*tau*(-2.0d0/(dr**2.0d0)- &
                       2.0d0/((r(i)**2.0d0)*(dz**2.0d0))))*u(i,0) &
                + kappa*tau*( &
                + (1.0d0/(dr**2.0d0)+1.0d0/(r(i)*2.0d0*dr))*u(i+1,0) &
                + (1.0d0/(dr**2.0d0)-1.0d0/(r(i)*2.0d0*dr))*u(i-1,0) &
                + 1.0d0/((r(i)**2.0d0)*(dz**2.0d0))*u(i,1) &
                + 1.0d0/((r(i)**2.0d0)*(dz**2.0d0))*u(i,M-1)) 
 enddo

 do i=1,N-1
  u_new(i,M)=u_new(i,0)
 enddo
 do j=0,M
  u_new(0,j)=BC_T1
  u_new(N,j)=BC_T2
 enddo

! do j=0,M
!  write(91,*) u_new(0:N,j) 
! enddo


 do i=0,N
  do j=0,M
   u(i,j)=u_new(i,j)
  enddo
 enddo

!enddo



end subroutine polar_2d_heat




end module mmat_FVM

!-------------------------------------------------------

FUNCTION NORM_2d(X1,Y1,X2,Y2)
IMPLICIT NONE

REAL(KIND=8), INTENT(IN) :: X1,Y1,X2,Y2
REAL(KIND=8)             :: NORM_2d


NORM_2d = sqrt((X1-X2)**2.0d0 + (Y1-Y2)**2.0d0)

END FUNCTION norm_2d

!--------------------------------------------
!            function Uprescribe,Vprescribe
!---------------------------------------------
FUNCTION Uprescribe(local_probtype,iten,t,x,y)
USE GeneralClass
IMPLICIT NONE

integer, intent(in) :: local_probtype,iten
REAL(KIND=8),INTENT(IN):: t,x,y
REAL(KIND=8)           :: Uprescribe
REAL(KIND=8)           :: local_pi

local_pi=4.0d0*atan(1.0d0)
if (local_probtype.eq.0) then
        Uprescribe=0.0d0
else if (local_probtype.eq.1) then
        Uprescribe=0.0d0
else if (local_probtype.eq.2) then
        Uprescribe=0.0d0
else if (local_probtype.eq.3) then
        Uprescribe=0.0d0
else if (local_probtype.eq.4) then
        Uprescribe=0.0d0
else if (local_probtype.eq.400) then
        Uprescribe=0.0d0
else if (local_probtype.eq.5) then
        Uprescribe=0.0d0
else if (local_probtype.eq.13) then
        Uprescribe=0.0d0
else if (local_probtype.eq.14) then
        Uprescribe=0.0d0
else if (local_probtype.eq.15) then
        Uprescribe=0.0d0
else if (local_probtype.eq.16) then
        Uprescribe=0.0d0
else if (local_probtype.eq.17) then
        Uprescribe=0.0d0
else if (local_probtype.eq.19) then
        Uprescribe=0.0d0
else if (local_probtype.eq.20) then
        Uprescribe=0.0d0
else
        print *,"local_probtype invalid"
        stop
endif


!U=100.0d0*y*(1.0d0-y)*(local_pi/2.0d0-atan(x))
!U=2.0d0
!u=1.0d0
!if (nm .le. 0.25d0) THEN

!u= 2.0d0*local_pi*(2.0d0*y/100.0d0-1.0d0)*(1.0d0-((2*x/100.0d0-1.0d0)**2))
!u = 2.0d0*y -100.0d0

!   u = -2.0d0*(y-0.5d0)	
!else
!  u = 0.0d0
!endif

!u=-2.0d0*local_pi*y*(1.0d0-x**2.0d0)

END FUNCTION Uprescribe

!------------------------
FUNCTION Vprescribe(local_probtype,iten,t,x,y)
USE GeneralClass
IMPLICIT NONE

integer, intent(in) :: local_probtype,iten
REAL(KIND=8),INTENT(IN):: t,x,y
REAL(KIND=8)           :: Vprescribe
REAL(KIND=8)           :: local_pi

local_pi=4.0d0*atan(1.0d0)
if (local_probtype.eq.0) then
        Vprescribe=0.0d0
else if (local_probtype.eq.1) then
        Vprescribe=0.0d0
else if (local_probtype.eq.2) then
        Vprescribe=0.0d0
else if (local_probtype.eq.3) then
        Vprescribe=0.0d0
else if (local_probtype.eq.4) then
        Vprescribe=0.0d0
else if (local_probtype.eq.400) then
        Vprescribe=0.0d0
else if (local_probtype.eq.5) then
        Vprescribe=0.0d0
else if (local_probtype.eq.13) then
        Vprescribe=0.0d0
else if (local_probtype.eq.14) then
        Vprescribe=0.0d0
else if (local_probtype.eq.15) then
        Vprescribe=0.0d0
else if (local_probtype.eq.16) then
        Vprescribe=0.0d0
else if (local_probtype.eq.17) then
        Vprescribe=0.0d0
else if (local_probtype.eq.19) then
        Vprescribe=0.0d0
else if (local_probtype.eq.20) then
        Vprescribe=0.0d0
else
        print *,"local_probtype invalid"
        stop
endif

!V=atan(0.1d0*x*(1.0d0-x)*y*(1.0d0-y)*(1.0d0+t))
!V=1.0d0
!v=-3.0d0

!v= 2.0d0*(x-0.5d0)*exp(-((y-0.5d0)**2+(x-0.5d0)**2)/100.0d0)


!v = 2.0d0*local_pi*(2*x/100.0d0-1.0d0)*(((2.0d0*y/100.0d0-1.0d0)**2)-1.0d0)

!v = 100.0d0 - 2.0d0*x

!if(nm .le. 0.25d0)then
!   v= 2.0d0*(x-0.5d0)
!else
!   v=0
!endif

!v=-2.0d0*Pi*x*(y**2.0d0-1.0d0)

END FUNCTION Vprescribe


subroutine LSprescribe(local_probtype,local_nmat,local_sdim, &
                local_time,x,y,local_dx,local_LS)
use mmat_FVM, only: dist_concentric
IMPLICIT NONE

integer, intent(in) :: local_probtype,local_nmat,local_sdim
real(kind=8),intent(in) :: local_time,x,y
real(kind=8),intent(in) :: local_dx(local_sdim)
real(kind=8),intent(out) :: local_LS(local_nmat*(local_sdim+1))
integer imat,dir
real(kind=8) :: mag,xplus,xminus,yplus,yminus,LSplus,LSminus
real(kind=8) :: nrm(local_sdim)

 do imat=1,local_nmat
   ! dist_concentric also declared in multimat_FVM.F90
  call dist_concentric(imat,x,y,local_LS(imat),local_probtype)
  do dir=1,local_sdim
   xplus=x 
   xminus=x 
   yplus=y 
   yminus=y
   if (dir.eq.1) then
    xplus=x+0.5d0*local_dx(dir) 
    xminus=x-0.5d0*local_dx(dir) 
   else if (dir.eq.2) then
    yplus=y+0.5d0*local_dx(dir) 
    yminus=y-0.5d0*local_dx(dir) 
   else
    print *,"dir invalid"
    stop
   endif
   call dist_concentric(imat,xplus,yplus,LSplus,local_probtype)
   call dist_concentric(imat,xminus,yminus,LSminus,local_probtype)
   nrm(dir)=(LSplus-LSminus)/local_dx(dir)
  enddo ! dir=1..local_sdim
  mag=0.0d0
  do dir=1,local_sdim
   mag=mag+nrm(dir)**2
  enddo
  mag=sqrt(mag)
  do dir=1,local_sdim
   if (mag.eq.0.0d0) then
    nrm(dir)=0.0d0
   else if (mag.gt.0.0d0) then
    nrm(dir)=nrm(dir)/mag
   else
    print *,"mag invalid"
    stop
   endif
   local_LS(local_nmat+(imat-1)*local_sdim+dir)=nrm(dir)
  enddo !dir=1..local_sdim
 enddo ! imat=1..local_nmat
end subroutine LSprescribe

      subroutine polar_cart_interpolate(Np,Mp,upolar,pcenter,rlo,rhi,x_in, ux)
      use mmat_FVM, only : rad_cal

      implicit none
        
      integer,intent(in)          :: Np,Mp
      real(kind=8),intent(in),dimension(0:Np,0:Mp)  :: upolar
      real(kind=8),intent(in)          :: pcenter(2)
      real(kind=8),intent(in)          :: rlo,rhi
      real(kind=8),intent(in)          :: x_in(2)

      real(kind=8)                     :: ux
      real(kind=8)                     :: xr,xz
      integer                          :: rl,zl      
      real(kind=8)                     :: dr,dz
      real(kind=8)                     :: r_interp1,r_interp2
      real(kind=8)                     :: alpha,beta

      if ((rlo.gt.0.0d0).and.(rhi.gt.rlo)) then
       ! do nothing
      else
       print *,"rlo or rhi invalid"
       stop
      endif
      dr = (rhi-rlo)/real(Np,8)     ! 0..Np
      dz = 2.0d0*Pi/real(Mp,8)      ! 0..Mp

      xr = sqrt((x_in(1)-pcenter(1))**2.0d0 + (x_in(2)-pcenter(2))**2.0d0)
      if (xr.le.rlo) then
       xr=rlo
      else if (xr.ge.rhi) then
       xr=rhi
      else if ((xr.ge.rlo).and.(xr.le.rhi)) then
       ! do nothing
      else
       print *,"xr invalid"
       stop
      endif
      rl = floor((xr-rlo)/dr) 
      if (rl.eq.Np) then
       rl=Np-1
      endif
      
      call rad_cal(x_in,pcenter,xz)
      zl = floor(xz/dz)
      if (zl.lt.0) then
       zl=zl+Mp
      else if ((zl.ge.0).and.(zl.lt.Mp)) then
       ! do nothing
      else if (zl.eq.Mp) then
       zl=zl-Mp
      else
       print *,"zl invalid"
       stop
      endif
      if ((rl.ge.0).and.(rl.lt.Np).and. &
          (zl.ge.0).and.(zl.lt.Mp)) then
          
       alpha= (xr-(rlo+rl*dr))/dr
       beta=(xz-zl*dz)/dz

       r_interp1=upolar(rl,zl)*(1.0d0-alpha)+upolar(rl+1,zl)*alpha
       r_interp2=upolar(rl,zl+1)*(1.0d0-alpha)+upolar(rl+1,zl+1)*alpha
      
       ux=r_interp1*(1.0d0-beta)+r_interp2*beta

      else
       print *,"rl or zl invalid"
       stop
      endif
   
    
      end subroutine polar_cart_interpolate


    subroutine find_polar_cart_inter(Np,Mp,upolar,pcenter,rlo,rhi,x_in,  &
                    diflag, dux)
    use mmat_FVM, only : rad_cal
    implicit none 

     integer,intent(in)          :: Np,Mp
     real(kind=8),intent(in),dimension(0:Np,0:Mp)  :: upolar
     real(kind=8),intent(in)          :: pcenter(2)
     real(kind=8),intent(in)          :: rlo,rhi
     real(kind=8),intent(in)          :: x_in(2)

     real(kind=8)                     :: dux
     real(kind=8)                     :: xz
     integer                          :: zl      
     real(kind=8)                     :: dr,dz
     real(kind=8)                     :: du1,du2
     real(kind=8)                     :: beta
     
     integer                          :: diflag

      dr = (rhi-rlo)/real(Np,8)     
      dz = 2.0d0*Pi/real(Mp,8)

      call rad_cal(x_in,pcenter,xz)
      zl = floor(xz/dz)
      beta=(xz-zl*dz)/dz 

      if(diflag .eq. 1)then 
       du1= (upolar(1,zl)-upolar(0,zl))/dr
       du2=(upolar(1,zl+1)-upolar(0,zl+1))/dr
      elseif(diflag .eq. 2)then
       du1= (upolar(Np-1,zl)-upolar(Np,zl))/dr
       du2=(upolar(Np-1,zl+1)-upolar(Np,zl+1))/dr
      else
       print *,"check diflag",diflag
       stop
      endif 

      dux = du1*(1.0d0-beta)+du2*beta

    end subroutine find_polar_cart_inter


! GeneralClass is in: vof_cisl.F90
!------------------------------------- test_flag==2
function exact_temperature(x_vec,t,im,probtype_in,nmat_in,alpha)
Use GeneralClass
USE probmain_module 
USE supercooled_exact_sol
USE variable_temperature_drop
implicit none

integer,intent(in)         :: im,probtype_in,nmat_in
real(kind=8),intent(in)    :: alpha(nmat_in)
real(kind=8),intent(in)    :: t
real(kind=8),intent(in)    :: x_vec(3)
real(kind=8)               :: exact_temperature
real(kind=8)              :: radius,theta,r1,r2
real(kind=8)              :: TLO,THI,yI,yHI,a1,b1,a2,b2
real(kind=8)              :: mypi,delx,dely
real(kind=8)              :: base_time
real(kind=8)              :: exact_rad
real(kind=8)              :: refc
real(kind=8)              :: radial_slope

 mypi=4.0d0*atan(1.0d0)

 if (t.ge.0.0d0) then
  ! do nothing
 else
  print *,"t invalid"
  stop
 endif
 if ((im.ge.1).and.(im.le.nmat_in)) then
  ! do nothing
 else
  print *,"im out of range"
  stop
 endif

 if (probtype_in.eq.1) then

  if (nmat_in.eq.3) then
   ! do nothing
  else
   print *,"nmat_in invalid"
   stop
  endif

  if (im.eq.2) then
   r1=radcen-radeps
   r2=radcen+radeps
   if (r2.gt.r1) then
    ! do nothing
   else
    print *,"r2 or r1 invalid"
    stop
   endif
      
   delx=x_vec(1)-0.5d0
   dely=x_vec(2)-0.5d0
   radius = sqrt(delx**2.0d0 +dely**2.0d0)
 
       ! x=r cos(theta)
       ! y=r sin(theta)
   if (radius.le.radeps/1000.0) then
    theta=0.0
   else if ((delx.ge.0.0).and.(dely.ge.0.0)) then
    theta=acos(delx/radius)
   else if ((delx.le.0.0).and.(dely.ge.0.0)) then
    theta=acos(abs(delx)/radius)
    theta=mypi-theta
   else if ((delx.le.0.0).and.(dely.le.0.0)) then
    theta=acos(abs(delx)/radius)
    theta=mypi+theta
   else if ((delx.ge.0.0).and.(dely.le.0.0)) then
    theta=acos(delx/radius)
    theta=2.0d0*mypi-theta
   else
    print *,"delx or dely invalid"
    stop
   endif

   radial_slope=1.0d0/(r2-r1)
   if (radial_variation.eq.1) then
    ! do nothing
   else if (radial_variation.eq.0) then
    radial_slope=0.0d0
   else
    print *,"radial_variation invalid"
    stop
   endif

   exact_temperature=2.0d0+sin(theta)*exp(-t/(radcen**2))+ &
           radial_slope*(radius-radcen)

  else if ((im.eq.1).or.(im.eq.3)) then
   exact_temperature=0.0
  else
   print *,"im invalid 6"
   stop
  endif

 else if (probtype_in.eq.3) then
  if (nmat_in.ne.2) then
   print *,"nmat_in invalid"
   stop
  endif
  delx=x_vec(1)-0.5d0
  dely=x_vec(2)-0.5d0
  radius = sqrt(delx**2.0d0 +dely**2.0d0)
  call solidification_front_time_driver(base_time,radblob)
  call solidification_front_radius_driver(base_time+t,exact_rad)
  if (radius.le.exact_rad) then
   exact_temperature=saturation_temp(2) 
  else if (radius.ge.exact_rad) then
   call liquid_temperature_driver(radius,base_time+t,exact_temperature)
  else
   print *,"radius invalid"
   stop
  endif
 else if (probtype_in.eq.4) then
  if (nmat_in.ne.2) then
   print *,"nmat_in invalid"
   stop
  endif
  delx=x_vec(1)-0.5d0
  dely=x_vec(2)-0.5d0
  radius = sqrt(delx**2.0d0 +dely**2.0d0)
  call disk_interp_T(radius,2,exact_temperature)

 else if (probtype_in.eq.5) then
  if (nmat_in.ne.2) then
   print *,"nmat_in invalid"
   stop
  endif
  if (x_vec(1).le.0.2d0+t) then
   exact_temperature=273.0d0
  else if (x_vec(1).ge.0.2d0+t) then
   exact_temperature=272.0d0+exp(-(x_vec(1)-0.2d0-t))
  else
   print *,"x_vec bust"
   stop
  endif

 else if ((probtype_in.eq.0).or. &
          (probtype_in.eq.2)) then
  if (nmat_in.ne.2) then
   print *,"nmat_in invalid"
   stop
  endif
     ! T1=a1+b1 y
     ! T2=a2+b2 (y-yI)
     ! b1 * k1 = b2 * k2
     ! a1=TLO
     ! a2=a1+b1 * yI
     ! a2=TLO+b1 * yI
     ! a2+b2(yHI-yI)=THI
     ! a2+(b1 k1/k2)(yHI-yI)=THI
     ! TLO+b1 yI+(b1 k1/k2)(yHI-yI)=THI
     ! TLO k2 + b1 (yI k2 + k1 (yHI - yI))=THI k2
     ! b1(yI k2 + k1 (yHI-yI))=(THI-TLO)k2
  TLO=3.0d0
  THI=2.0d0
  yI=0.3d0
  yHI=1.0
  a1=TLO
  b1=(THI-TLO)*alpha(2)/(yI * alpha(2) +alpha(1)*(yHI-yI))
  b2=b1*alpha(1)/alpha(2)
  a2=a1+b1*yI 

  if (probtype_in.eq.0) then

   if (im.eq.1) then
    exact_temperature=a1+b1*x_vec(2)
   else if (im.eq.2) then
    exact_temperature=a2+b2*(x_vec(2)-yI)
   else
    print *,"im invalid 7"
    stop
   endif

  else if (probtype_in.eq.2) then

   if (im.eq.1) then
    exact_temperature=a1+b1*x_vec(1)
   else if (im.eq.2) then
    exact_temperature=a2+b2*(x_vec(1)-yI)
   else
    print *,"im invalid 8"
    stop
   endif
  else
   print *,"probtype_in invalid12 ",probtype_in
   stop
  endif
 elseif(probtype_in.eq.13)then

   ! declared in vof_cisl.F90
  if ((dirichlet_pentafoil.eq.1).and.(1.eq.0)) then 

   if (nmat_in.ne.3) then
    print *,"nmat_in invalid"
    stop
   endif
   if (im.eq.2) then
      
    delx=x_vec(1)-0.5d0-sqrt(5.0)*0.02d0
    dely=x_vec(2)-0.5d0-sqrt(5.0)*0.02d0
    radius = sqrt(delx**2.0d0 +dely**2.0d0)
 
!        x=r cos(theta)
!        y=r sin(theta)
    if (radius.le.radeps/1000.0) then
     theta=0.0
    else if ((delx.ge.0.0).and.(dely.ge.0.0)) then
     theta=acos(delx/radius)
    else if ((delx.le.0.0).and.(dely.ge.0.0)) then
     theta=acos(abs(delx)/radius)
     theta=mypi-theta
    else if ((delx.le.0.0).and.(dely.le.0.0)) then
     theta=acos(abs(delx)/radius)
     theta=mypi+theta
    else if ((delx.ge.0.0).and.(dely.le.0.0)) then
     theta=acos(delx/radius)
     theta=2.0d0*mypi-theta
    else
     print *,"delx or dely invalid"
     stop
    endif

    exact_temperature=2.0d0+sin(theta)*exp(-t)
   else if ((im.eq.1).or.(im.eq.3)) then
    exact_temperature=0.0
   else
    print *,"im invalid 6"
    stop
   endif

  else if ((dirichlet_pentafoil.eq.0).or.(1.eq.1)) then

   if (im.eq.2) then
    refc= (0.02d0*sqrt(5.0d0)+1.0d0)/2.0d0 ! refc=0.0d0 default
    exact_temperature=((x_vec(1)-refc)**2.0d0+(x_vec(2)-refc)**2.0d0)*exp(-t)
   else if ((im.eq.1).or.(im.eq.3)) then
    exact_temperature=0.0
   else
    print *,"im invalid"
    stop
   endif

  else
   print *,"dirichlet_pentafoil invalid"
   stop
  endif

 elseif(probtype_in .eq. 14)then
   exact_temperature = (x_vec(1)**2+x_vec(2)**2)*exp(-t)

 elseif(probtype_in .eq. 15)then
 
  exact_temperature=0.0d0

  ! refc= (1.0d0)/2.0d0
  ! exact_temperature = ((x-refc)**2.0d0 + (y-refc)**2.0d0)*exp(-t)

 elseif(probtype_in.eq.16)then

  call disk_interp_T(x_vec(2),2,exact_temperature)

 elseif(probtype_in.eq.17)then
  exact_temperature=0.0d0
 elseif(probtype_in.eq.20)then

  exact_temperature=0.0d0

 else if (probtype_in.eq.19) then
  if (im.eq.1) then
   exact_temperature=0.0d0
  else if (im.eq.3) then
   exact_temperature=0.0d0
  else if (im.eq.2) then
    ! Np,Mp,upolar,pcenter,rlo,rhi are defined in vof_cisl.F90
   call polar_cart_interpolate(Np,Mp,upolar,pcenter,rlo,rhi,x_vec, &
     exact_temperature)
  else
   print *,"im invalid in exact temp"
   stop
  endif
 else if (probtype_in.eq.400)then

  exact_temperature=0.0d0

 else
  print *,"probtype_in invalid13 ",probtype_in
  stop
 endif

end function exact_temperature

