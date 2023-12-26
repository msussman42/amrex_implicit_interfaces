#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "PROB_F.H"
#include "LEVEL_F.H"

#include "EXTRAP_COMP.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

#define maxnumline 10

module global_distance_module
use amrex_fort_module, only : amrex_real
use probcommon_module

implicit none

contains


! ----------
! CODY ESTEBE
! ----------

! from: nozzle2D.f90
! Phi>0 in the solid
!2d nozzle test case, lower half (Brusiani et al 2013)
!return nozzle signed distance function Phi (dist>0 in the substrate/solid)
subroutine nozzle2d(x_cm,y_cm,Phi) 
 implicit none
 !spatial coordinates, domain: x=[0:7000],y=[-300:2000] microns
 !  0<=x<=0.7cm
 ! -0.03 cm <= y <= 0.2 cm
 real(amrex_real), INTENT(in) :: x_cm, y_cm 
 real(amrex_real), INTENT(out) :: Phi !init to 0
 
 integer :: nozzletype, insideflag
 real(amrex_real) :: nl_width, nr_width, r, nl, nr !for nozzle entrance/exit
 real(amrex_real) :: l1, l2, a, bottom !for rounded corner and tangent nozzle line
 real(amrex_real) :: m, m2, phi1, phi2, vertDisplace !for inner slope of nozzle
 real(amrex_real) :: round_l, round_r !for nozzle corners
 real(amrex_real) :: x_um,y_um
 real(amrex_real) :: scaling_factor

 scaling_factor=radblob5

 x_um=x_cm*scaling_factor 
 y_um=(yblob/scaling_factor-y_cm)*scaling_factor

 !nozzle configuration
 !J:0, U:1, W:2
 nozzletype = NINT(radblob2)
 nl_width=radblob3
 nr_width=radblob4

 r=radblob6 !inlet radius of curvature

 !vertical height of nozzle entrance/exit
 !yblob should be the domain height divided by 2.
 nl = yblob-nl_width/2.0d0-r
 nr = yblob-nr_width/2.0d0

 if ((xblob.eq.0.0d0).and. &   !no channel on either side of core nozzle
     (nozzletype.eq.0).and. &!J nozzle
     (nl_width.eq.nr_width)) then ! nl_width = nr_width

  Phi=abs(y_um)-nl_width/2.0d0

 else if ((xblob.gt.0.0d0).or. &
          (nozzletype.eq.1).or. & ! U nozzle
          (nozzletype.eq.2).or. & ! W nozzle
          (nl_width.ne.nr_width)) then

  !line tangent to circle at point a,b (for nozzle inlet)
  !
  !                            .(l1,l2)
  !                          /
  !                        /
  !                      /
  !                    /
  !                  /
  !                /
  !              / _-----_
  !            /-          \
  !    (a,b).//              \
  !           |               |
  !          |                 |
  !          |       .(0,0)    |
  !          |                |
  !           \              /
  !             -          /
  !               -------
  
  l1=2.0d0*yblob-(xblob+r)
  l2=nr-(nl)
  a=(2*l1*(r**2)-sqrt(4*(l1**2)*(r**4)-4*((l1**2)+(l2**2))* &
    ((r**4)-(r**2)*(l2**2))))/(2*(l1**2+l2**2)) !-sqrt if slope pos and tangent to upper side of circle, else +sqrt
  if (l2.gt.zero) then
   bottom=(r**2-l1*a)/l2
  else if (l2.eq.zero) then
   bottom=zero
  else
   print *,"l2 invalid"
   stop
  endif
  a=a+(xblob+r)
  bottom=bottom+(nl)

  !        | \
  !        | 
  !        :  ...  
  !        | <vertDisplace
  !        |          \
  !    \\  |          /\
  !     \\ |      1/    \
  !      \\|phi1/        \
  !       \\ / phi2       \
  !        \\-------       \
  !         \\              \

  m = (nr-bottom)/(2.0d0*yblob-a) !slope nozzle
  m2 = -1.0/m
  phi1 = abs(atan(m2))
  phi2 = Pi/2.0-phi1
  vertDisplace = 1.0/cos(phi2)


  insideflag = 0
  if (((y_um.ge.((-1/m)*(x_um-a)+bottom)).and. &
      (x_um.le.2.0d0*yblob).and. &
      (y_um.le.(m*(x_um-a)+bottom))).or. &
     ((x_um.ge.xblob).and. &
      (x_um.le.2.0d0*yblob).and. &
      (y_um.le.nl))) then
   Phi = MIN(x_um-xblob, 2.0d0*yblob-x_um) !inside,column vert
   insideflag = 1
  endif
  if (((y_um.ge.((-1/m)*(x_um-a)+bottom)).and. &
       (x_um.le.2.0*yblob).and.  &
       (y_um.le.(m*(x_um-a)+bottom))) .or. &
      ((x_um.ge.xblob).and. &
       (x_um.le.2.0*yblob).and. &
       (y_um.le.nl))) then
   Phi = MIN(Phi, ((m*(x_um-a)+bottom)-y_um)/vertDisplace) !inside,column horiz
   insideflag = 1
  endif
  if ((x_um.ge.xblob).and. &
      (x_um.le.xblob+r).and. &
      (y_um.ge.nl).and. &
      (y_um.le.nl+r).and. &
      (y_um.le.((-1/m)*(x_um-a)+bottom))) then
   Phi = r-sqrt((x_um-(xblob+r))**2+(y_um-nl)**2) !inside, circle
   insideflag = 1
  endif
  if ((y_um.le.0).and. &
      (x_um.ge.xblob).and. &
      (x_um.le.2.0d0*yblob)) then
   round_l = sqrt((x_um-xblob)**2+(y_um-0)**2) !inside bottom left corner
   round_r = sqrt((x_um-2.0d0*yblob)**2+(y_um-0)**2) !inside bottom right corner
   Phi = MIN(round_l,round_r)
  endif
  if (y_um.le.0) then
   Phi = MAX(Phi,0-y_um) !inside,bottom
   insideflag = 1
  endif

  if (insideflag.eq.0) then
   !outside
   !left region,vert
   if ((y_um.gt.0).and. &
       (x_um.lt.xblob).and. &
       (y_um.lt.nl)) then
    Phi = x_um-xblob
   endif
   !right region,vert
   if ((y_um.gt.0).and. &
       (x_um.gt.2.0d0*yblob).and. &
       (y_um.lt.nr)) then
    Phi = 2.0d0*yblob-x_um
   endif
   !nozzle corners
   round_l=0
   round_r=0
   if ((x_um.le.2.0d0*yblob).and. &
       (y_um.gt.nl).and. &
       (y_um.le.((-1/m)*(x_um-2.0d0*yblob)+nr))) then
    round_l = r-sqrt((x_um-(xblob+r))**2+(y_um-nl)**2)
   endif
   if ((y_um.gt.((-1/m)*(x_um-2.0d0*yblob)+nr)).and. &
       (y_um.gt.nr)) then 
    round_r = 0-sqrt((x_um-2.0d0*yblob)**2+(y_um-nr)**2)
   endif
   Phi = MIN(Phi,round_l,round_r)  
   !nozzle
   if ((y_um.gt.(m*(x_um-a)+bottom)).and. &
       (x_um.le.2.0d0*yblob).and. &
       (y_um.gt.((-1/m)*(x_um-a)+bottom)).and. &
       (y_um.le.((-1/m)*(x_um-2.0d0*yblob)+nr))) then
    Phi = MAX(Phi, ((m*(x_um-a)+bottom)-y_um)/vertDisplace)
   endif
   !left region,horiz
   if ((y_um.gt.0).and.(x_um.lt.xblob)) then
    Phi = MAX(Phi,0-y_um)
   endif
   !right region,horiz
   if ((y_um.gt.0).and.(x_um.gt.2.0d0*yblob)) then
    Phi = MAX(Phi,0-y_um)
   endif
  endif

 else
  print *,"parameters invalid in nozzle2d"
  stop
 endif

 Phi=Phi/scaling_factor

end subroutine nozzle2d

! ------------------------
! END CODY ESTEBE
! ----------------------


        ! negative in the sphere
      subroutine spheredist(x,y,z,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,z,dist

      if (SDIM.eq.2) then
       if (abs(y-z).gt.EPS8) then
        print *,"expecting y=z"
        stop
       endif
       dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob
      else if (SDIM.eq.3) then
       dist=sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)-radblob
      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine spheredist

        ! dist>0 in the sphere
      subroutine stainless_steel_dist_rate( &
       x,y,z,time,dist,rate)
      IMPLICIT NONE

      real(amrex_real) x,y,z,time,dist,rate
      real(amrex_real) IMPACT_DIST,IMPACT_RATE,IMPACT_TIME
      real(amrex_real) T_DELAY,TPART1,UPART1,UPART2,ycen


      IMPACT_DIST=yblob-radblob
      if (IMPACT_DIST.le.zero) then
       print *,"IMPACT_DIST invalid"
       stop
      endif
      IMPACT_RATE=230.0D0
      IMPACT_TIME=IMPACT_DIST/IMPACT_RATE
      if (IMPACT_TIME.le.zero) then
       print *,"IMPACT_TIME invalid"
       stop
      endif
      T_DELAY=3.76688776419980306D-004
      TPART1=4.10014276627081750D-004
      UPART1=179.92175975809258D0
      UPART2=34.086095584802656D0
      if (time.le.IMPACT_TIME) then
       ycen=yblob-time*IMPACT_RATE
       rate=-IMPACT_RATE
      else if ((time.ge.IMPACT_TIME).and. &
               (time.le.IMPACT_TIME+T_DELAY)) then
       ycen=radblob
       rate=zero
      else if ((time.ge.IMPACT_TIME+T_DELAY).and. &
               (time.le.IMPACT_TIME+T_DELAY+TPART1)) then
       ycen=radblob+(time-IMPACT_TIME-T_DELAY)*UPART1
       rate=UPART1
      else if (time.ge.IMPACT_TIME+T_DELAY+TPART1) then
       ycen=radblob+TPART1*UPART1+ &
         (time-(IMPACT_TIME+T_DELAY+TPART1))*UPART2
       rate=UPART2
      endif
      if (ycen.lt.radblob) then
       print *,"ycen invalid"
       print *,"time=",time
       print *,"IMPACT_TIME=",IMPACT_TIME
       print *,"yblob=",yblob
       print *,"IMPACT_RATE=",IMPACT_RATE
       print *,"T_DELAY=",T_DELAY
       print *,"UPART1=",UPART1
       print *,"UPART2=",UPART2
       print *,"TPART1=",TPART1
       stop
      endif

      dist=radblob-sqrt((x-xblob)**2+(y-ycen)**2)

      return
      end subroutine stainless_steel_dist_rate

! positive in the solid
      subroutine materialdistsolid(x,y,z,dist,time,im)
      use global_utility_module
      use USERDEF_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      use TSPRAY_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x,y,z,time
      integer, INTENT(in) :: im
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real) tadv,xprime,yprime,zprime
      real(amrex_real) distz,disty,steel_rate
      real(amrex_real) xvec(SDIM)
      real(amrex_real) dist_array(num_materials)

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid10"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in materialdistsolid"
      else if (time.lt.zero) then
       print *,"time invalid in materialdistsolid"
       stop
      else
       print *,"time bust in materialdistsolid"
       stop
      endif

      if (is_rigid(im).ne.1) then
       print *,"is_rigid invalid GLOBALDIST.F90"
       stop
      endif

      if ((adv_dir.lt.1).or.(adv_dir.gt.2*SDIM+1)) then
       print *,"adv_dir invalid materialdistsolid (1)"
       stop
      endif

      xprime=x
      yprime=y
      zprime=z
      xvec(1)=x
      xvec(2)=y
      if (SDIM.eq.3) then
       xvec(SDIM)=z
      endif

      if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then

       if (is_in_probtype_list().eq.1) then

        call SUB_LS(xvec,time,dist_array,num_materials)
        dist=dist_array(im)

       else if (probtype.eq.402) then ! user defined thermal spray problem
        call TSPRAY_LS(xvec,time,dist_array)
        dist=dist_array(im)

       else if (probtype.eq.412) then ! user defined cavitation problem
        call CAV2Dstep_LS(xvec,time,dist_array)
        dist=dist_array(im)

       else if (probtype.eq.413) then ! zeyu's droplet impact problem
        call ZEYU_droplet_impact_LS(xvec,time,dist_array)
        dist=dist_array(im)

       else if (probtype.eq.311) then ! user defined
        call USERDEF_LS(xvec,time,dist_array)
        dist=dist_array(im)

        ! cavitation (in materialdistsolid)
       else if ((probtype.eq.46).and.(axis_dir.eq.10)) then
          ! dist>0 in the steel sphere
        call stainless_steel_dist_rate(xprime,yprime,zprime,time, &
         dist,steel_rate)

       else if ((probtype.eq.46).and.(axis_dir.eq.20)) then

          ! dist>0 in the substrate (solid) region.
          ! nozzle2d is declared in GLOBALDIST.F90
        call nozzle2d(xprime,yprime,dist)  

        ! materialdistsolid: flapping wing
       else if (probtype.eq.701) then
  
         ! dist>0 in the airfoil
        call naca_dist(x,y,z,time,dist)

        ! materialdistsolid above: user defined
        ! materialdistsolid below: cases that use soliddist.
       else  

        if (probtype.eq.32) then

         tadv=time

         if (adv_dir.eq.1) then
          if (advbot.ne.zero) then
           xprime=xprime-advbot*tadv
          else if ((xblob3.ne.zero).and.(xblob4.ne.zero)) then
           xprime=xprime-xblob3*(one-cos(two*Pi*tadv/xblob4))
          endif
         else if (adv_dir.eq.2) then
          if (advbot.ne.zero) then
           yprime=yprime-advbot*tadv
          else if ((xblob3.ne.zero).and.(xblob4.ne.zero)) then
           yprime=yprime-xblob3*(one-cos(two*Pi*tadv/xblob4))
          endif
         else if ((adv_dir.eq.SDIM).and.(SDIM.eq.3)) then
          if (advbot.ne.zero) then
           zprime=zprime-advbot*tadv
          else if ((xblob3.ne.zero).and.(xblob4.ne.zero)) then
           zprime=zprime-xblob3*(one-cos(two*Pi*tadv/xblob4))
          endif
         else
          print *,"adv_dir invalid probtype=32"
          stop
         endif
         ! sphere impact free surface - materialdistsolid
        else if ((probtype.eq.531).and.(SDIM.eq.2)) then 
         if (axis_dir.eq.0) then
          tadv=time
          yprime=yprime-tadv*advbot
         else if ((axis_dir.eq.1).or.(axis_dir.eq.2).or. &
                  (axis_dir.eq.3)) then
          ! do nothing
         else 
          print *,"axis_dir invalid probtype=531"
          stop
         endif
        else if (probtype.eq.bubbleInPackedColumn) then ! in materialdistsolid
         continue
        endif

        ! positive in fluid (now in materialdistsolid)
        call soliddist(xprime,yprime,zprime,dist,im)  
        ! now make dist>0 in the solid.
        dist=-dist
       endif ! cases in which soliddist is called.

      else if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
               (FSI_flag(im).eq.FSI_SHOELE_CTML)) then 

! dist>0 in the solid
! eventually as the program is running FSI_MF multifab will be
! copied to the fortran.
! closest value(s) on same processor are used.

       dist=-99999.0

      else
       print *,"FSI_flag invalid in materialdistsolid"
       print *,"im,FSI_flag(im) ",im,FSI_flag(im)
       stop
      endif

! microfluidics problem (in materialdistsolid)
! dist is positive in the solid
      if ((probtype.eq.5700).and.(SDIM.eq.3)) then

       if ((axis_dir.eq.0).or.(axis_dir.eq.1).or. &
           (axis_dir.eq.3)) then

        if (zprime.gt.half*(probloz+probhiz)) then
         distz=zprime-probhiz
        else
         distz=probloz-zprime
        endif
        if (distz.gt.dist) then
         dist=distz
        endif 
         ! comsol problem
        if (axis_dir.eq.1) then
         if (yprime.lt.probloy) then
          disty=probloy-yprime
          if (disty.gt.dist) then
           dist=disty
          endif
         endif
        endif
       else if (axis_dir.eq.4) then
        ! do nothing
       else if (axis_dir.eq.5) then
        ! do nothing
       else
        print *,"axis_dir invalid for 5700 probtype"
        stop
       endif

      endif

      return
      end subroutine materialdistsolid


      subroutine naca_velocity(x,y,z,time,vel)
      IMPLICIT NONE

      real(amrex_real) x,y,z,time
      real(amrex_real) vel(SDIM)
      integer dir
      real(amrex_real) amp

      if (probtype.eq.701) then
       amp=yblob
       do dir=1,SDIM
        vel(dir)=zero
       enddo
       vel(2)=-amp*two*Pi*sin(two*Pi*time)
      else
       print *,"probtype invalid"
       stop
      endif

      return
      end subroutine naca_velocity

       ! dist>0 in the airfoil
      subroutine naca_dist(x,y,z,time,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,z,time,dist
      real(amrex_real) xprime,yprime,thick,c,H,amp,offset

      if (probtype.eq.701) then
       amp=yblob
       offset=amp*(cos(two*Pi*time)-one)
       yprime=y-offset
       if (x.le.zero) then
        dist=-sqrt(x*x+yprime**2) 
       else if (x.ge.one) then
        dist=-sqrt((x-one)**2+yprime**2)
       else
        thick=12.0/100.0  ! naca0012
        c=one
        xprime=x/c
        H=5.0*thick*c*(0.2969*sqrt(xprime)-0.1260*xprime- &
         0.3516*(xprime**2)+0.2843*(xprime**3)-0.1015*(xprime**4))
        if (yprime.ge.zero) then
         dist=H-yprime
        else
         dist=H+yprime
        endif
       endif
      else
       print *,"probtype invalid"
       stop
      endif

      return
      end subroutine naca_dist

       ! dist>0 in the crystal
      subroutine crystal_distance(x,y,z,dist,dx,bfact,im_project)
      IMPLICIT NONE
   
      integer bfact
      integer nPoly
      parameter(nPoly=4)
      real(amrex_real) xp(nPoly),yp(nPoly)

      real(amrex_real) dx(SDIM)
      real(amrex_real) x,y,z,dist
      real(amrex_real) tmp,dist0,dist12,x1,x2,y1,y2,xproj,yproj
      integer i,j,isSolid,im_project
      real(amrex_real) xcen1,ycen1,zcen1
      real(amrex_real) xcen2,ycen2,zcen2

      print *,"crystal_distance needs to be updated"
      stop

      if (bfact.lt.1) then
       print *,"bfact invalid11"
       stop
      endif
      if ((im_project.lt.1).or.(im_project.gt.num_materials)) then
       print *,"im_project invalid"
       stop
      endif

      if (SDIM.eq.2) then
       if (abs(z-y).gt.EPS8) then
        print *,"expecting z=y if 2d"
        stop
       endif
      endif

       ! in: crystal_distance
      if ((probtype.eq.46).and.(SDIM.eq.2)) then 
       if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
        print *,"no rigid body expected"
        stop
       else if (axis_dir.eq.10) then
        if (im_project.eq.num_materials) then
         call spheredist(x,y,z,dist) ! dist<0 in the sphere
         dist=-dist
        else
         print *,"expecting im_project=num_materials; num_materials=",num_materials
         stop
        endif
       else
        print *,"axis_dir invalid"
        stop
       endif
       ! crystal_distance
      else if ((probtype.eq.36).and.(axis_dir.eq.0)) then
       if (im_project.eq.1) then 
        print *,"not expecting im_project=1"
        stop
       else if (im_project.eq.2) then
        call spheredist(x,y,z,dist) ! dist<0 in the sphere
        dist=-dist
       else
        print *,"im_project invalid"
        stop
       endif

       ! radblob2=distance between 2 ends=2.5mm
       ! radblob=radius of barbell = 250 microns=250E-6m=250E-4cm
       ! 250E-3 mm=.25mm
       ! radblob3=radius of connecting cylinder=.23mm
       ! radblob4=frequency=2.8 hertz = 2.8 cycles/s
       ! =2.8 * 180 degrees/s = 2.8 pi radians/s 
       ! T=1/2.8 seconds
       ! cycle is 45 degrees up, 45 down, 45 down 45 up
      else if ((probtype.eq.531).and.(axis_dir.eq.3).and.(SDIM.eq.2)) then
       xcen1=xblob-half*radblob2
       xcen2=xcen1+radblob2
       ycen1=zero
       ycen2=zero
       zcen1=yblob
       zcen2=yblob
#if (AMREX_SPACEDIM==2)
       call barbelldist(x,ycen1,y,xcen1,ycen1,zcen1, &
        xcen2,ycen2,zcen2,radblob,radblob,radblob3,dist)
       dist=-dist
#endif
      else if ((probtype.eq.531).and.(axis_dir.eq.2).and.(SDIM.eq.2)) then
       xp(1)=xblob-radblob*0.5
       yp(1)=yblob-radblob*0.5
       xp(2)=xblob-radblob*0.5
       yp(2)=yblob+radblob*0.5
       xp(3)=xblob+radblob*0.5
       yp(3)=yblob+radblob*0.5
       xp(4)=xblob+radblob*0.5
       yp(4)=yblob-radblob*0.5
       dist=1.0e9     ! hold a position
       do i=1,nPoly       ! find the distance to the polygon
        x1=xp(i)
        y1=yp(i)
        if (i.eq.nPoly) then
           x2=xp(1)
           y2=yp(1)
        else
           x2=xp(i+1)
           y2=yp(i+1)
        endif
        dist12=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
        tmp=((x-x1)*(x2-x1)+(y-y1)*(y2-y1))/dist12
        if (tmp.lt.0.0) then
         dist0=sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1))
        elseif (tmp.gt.1.0) then
         dist0=sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2))
        else
         xproj=x1+tmp*(x2-x1)
         yproj=y1+tmp*(y2-y1)
         dist0=sqrt((x-xproj)*(x-xproj)+ &
            (y-yproj)*(y-yproj))
        endif
        dist=min(dist,dist0)
       enddo

       isSolid= -1
       do i=1,nPoly
        if (i.eq.nPoly) then
           j=1
        else
           j=i+1
        endif
        if (((yp(i).ge.y) .neqv. (yp(j).ge.y)).and. &
         (x.le.(xp(j)-xp(i))*(y-yp(i))/(yp(j)-yp(i))+xp(i)))then
           isSolid=-isSolid
        endif
       enddo
 
       if (isSolid.eq.-1) then ! outside 
        dist=-dist
       endif
      else if ((probtype.eq.531).and.(axis_dir.eq.1).and.(SDIM.eq.2)) then
       dist=radblob-sqrt( (x-xblob)**2+(y-yblob)**2 )
      else
       print *,"invalid crystal_distance option"
       stop
      endif

      return
      end subroutine crystal_distance


      subroutine crystal_centroid(rigid_centroid,im_project)
      IMPLICIT NONE

      integer nPoly,im_project
      parameter(nPoly=4)

      real(amrex_real) rigid_centroid(SDIM)
      real(amrex_real) xp(nPoly),yp(nPoly)
      integer i,dir

      print *,"crystal_centroid needs to be updated"
      stop

      if ((im_project.lt.1).or.(im_project.gt.num_materials)) then
       print *,"im_project invalid"
       stop
      endif

       ! in: crystal centroid
      if ((probtype.eq.46).and.(SDIM.eq.2)) then 
       if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
        print *,"no rigid body expected"
        stop
       else if (axis_dir.eq.10) then
        rigid_centroid(1)=xblob
        rigid_centroid(2)=yblob
        if (SDIM.eq.3) then
         rigid_centroid(SDIM)=zblob
        endif
       else
        print *,"axis_dir invalid"
        stop
       endif 
       ! in: crystal centroid
      else if ((probtype.eq.36).and.(axis_dir.eq.0)) then

       rigid_centroid(1)=xblob
       rigid_centroid(2)=yblob
       if (SDIM.eq.3) then
        rigid_centroid(SDIM)=zblob
       endif

      else if ((probtype.eq.531).and.(axis_dir.eq.3).and.(SDIM.eq.2)) then

       rigid_centroid(1)=xblob
       rigid_centroid(2)=yblob

        ! falling crystal
      else if ((probtype.eq.531).and.(axis_dir.eq.2).and.(SDIM.eq.2)) then 
       xp(1)=xblob-radblob*0.5
       yp(1)=yblob-radblob*0.5
       xp(2)=xblob-radblob*0.5
       yp(2)=yblob+radblob*0.5
       xp(3)=xblob+radblob*0.5
       yp(3)=yblob+radblob*0.5
       xp(4)=xblob+radblob*0.5
       yp(4)=yblob-radblob*0.5

       do dir=1,SDIM
        rigid_centroid(dir)=zero
       enddo
       do i=1,nPoly
        dir=1
        rigid_centroid(dir)=rigid_centroid(dir)+xp(i)
        dir=2
        rigid_centroid(dir)=rigid_centroid(dir)+yp(i)
       enddo
       do dir=1,SDIM
        rigid_centroid(dir)=rigid_centroid(dir)/(nPoly+1.-1.0)
       enddo
 
      else if ((probtype.eq.531).and.(axis_dir.eq.1).and.(SDIM.eq.2)) then
       dir=1
       rigid_centroid(dir)=xblob
       dir=2
       rigid_centroid(dir)=yblob
      else
       print *,"invalid crystal_centroid option"
       stop
      endif

      return
      end subroutine crystal_centroid


       ! ycen1=zcen1=ycen2=zcen2=1/2
       ! phi<0 in the object
      subroutine barbelldist(x,y,z, &
         xcen1,ycen1,zcen1, &
         xcen2,ycen2,zcen2, &
         r1,r2,r3,phi)
      IMPLICIT NONE

      real(amrex_real) x,y,z
      real(amrex_real) xcen1,ycen1,zcen1
      real(amrex_real) xcen2,ycen2,zcen2
      real(amrex_real) r,r1,r2,r3,phi
      real(amrex_real) phi1,phi2,phi3

      phi1=sqrt( (x-xcen1)**2+(y-ycen1)**2+(z-zcen1)**2 ) - r1
      phi2=sqrt( (x-xcen2)**2+(y-ycen2)**2+(z-zcen2)**2 ) - r2
      phi3=sqrt( (y-ycen1)**2+(z-zcen1)**2 ) - r3

      if (x.lt.xcen1) then
       phi=phi1
      else if (x.ge.xcen2) then 
       phi=phi2
      else if ((xcen1.le.x).and.(x.le.half*(xcen1+xcen2))) then
       r=sqrt((y-ycen1)**2+(z-zcen1)**2)
       if ((r.gt.r3).and.(x.le.xcen1+r1)) then
        phi=phi1
       else if (phi1.le.0) then
        phi = -sqrt(phi1**2+phi3**2)
       else 
        phi = phi3
       endif
      else if ((half*(xcen1+xcen2).le.x).and.(x.le.xcen2)) then
       r=sqrt((y-ycen2)**2+(z-zcen2)**2)
       if ((r.gt.r3).and.(x.ge.xcen2-r2)) then
        phi=phi2
       else if (phi2.lt.0) then
        phi = -sqrt(phi2**2+phi3**2)
       else 
        phi = phi3
       endif
      endif

      return
      end subroutine barbelldist

       !jetting_plate_dist is called from soliddist.
       !soliddist is called from materialdistsolid.
       !dist<0 in the plate.
      subroutine jetting_plate_dist(x,y,z,dist)
      use global_utility_module
      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: x,y,z
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real) :: aspect,offset,distplate,hugedist

      hugedist=99999.0

       ! underwater explosion, growth, and collapse of explosion bubble
       ! near a flat plate.
      if (probtype.eq.42) then
       aspect=xblob2
!      offset=2.54d0
       offset=radblob2
       distplate=yblob2

       if (radblob2.eq.zero) then
        dist=hugedist
       else if (radblob2.gt.zero) then
          ! negative in the square
        call squaredist(x,y,-aspect,aspect,yblob+distplate, &
         yblob+distplate+offset,dist)
       else
        print *,"radblob2 invalid for bubble jetting problem"
        stop
       endif
      else
       print *,"probtype invalid in jetting_plate_dist"
       stop
      endif

      return
      end subroutine jetting_plate_dist

       ! probtype.eq.299
      subroutine INIT_LS_SOLID_MELT(x,y,z,time,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,z,time,dist,dist2
      integer nSphere
      real(amrex_real) xSphere
      real(amrex_real) ySphere
      real(amrex_real) delta_sphere,hugedist
      integer i
 
      hugedist=1.0D+10
 
      if (num_materials.ne.3) then
       print *,"num_materials invalid in INIT_LS_SOLID_MELT"
       stop
      endif

      if (SDIM.eq.2) then
       if (abs(z-y).gt.VOFTOL) then
        print *,"z=y required in INIT_LS_SOLID_MELT"
        print *,"x,y,z = ",x,y,z
        print *,"probtype,axis_dir ",probtype,axis_dir
        stop
       endif
      endif

      nSphere=3

      if (axis_dir.eq.0) then
       delta_sphere=two
      else if (axis_dir.eq.1) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.2) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.3) then
       delta_sphere=two
      else
       print *,"axis_dir invalid init ls solid melt"
       stop
      endif

      dist=-hugedist

      do i=1,nSphere
       xSphere=xblob+(i-one)*delta_sphere
       ySphere=yblob
       dist2=radblob-sqrt((x-xSphere)**2+(y-ySphere)**2) 
       dist=max(dist,dist2)
      enddo

      return
      end subroutine INIT_LS_SOLID_MELT

       ! probtype.eq.299
      subroutine INIT_LS_LIQUID_MELT(x,y,z,time,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,z,time,dist,dist2,dist3,dist4
      integer nSphere
      real(amrex_real) xSphere
      real(amrex_real) ySphere
      real(amrex_real) delta_sphere,hugedist
      integer i
 
      hugedist=1.0D+10
 
      if (num_materials.ne.3) then
       print *,"num_materials invalid in INIT_LS_LIQUID_MELT"
       stop
      endif
      if (radblob2.le.radblob) then
       print *,"radblob2 invalid"
       stop
      endif
      if (SDIM.eq.2) then
       if (abs(z-y).gt.VOFTOL) then
        print *,"z=y required in INIT_LS_LIQUID_MELT"
        print *,"x,y,z = ",x,y,z
        print *,"probtype,axis_dir ",probtype,axis_dir
        stop
       endif
      endif

      nSphere=3

      if (axis_dir.eq.0) then
       delta_sphere=two
      else if (axis_dir.eq.1) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.2) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.3) then
       delta_sphere=two
      else
       print *,"axis_dir invalid init ls liquid melt"
       stop
      endif

      dist=-hugedist

      do i=1,nSphere
       xSphere=xblob+(i-one)*delta_sphere
       ySphere=yblob
       dist2=radblob-sqrt((x-xSphere)**2+(y-ySphere)**2)  ! metal
       dist3=radblob2-sqrt((x-xSphere)**2+(y-ySphere)**2) ! air
       if (dist2.ge.zero) then ! in the metal
        dist4=-dist2
       else if (dist3.le.zero) then ! in the air
        dist4=dist3
       else if ((dist2.le.zero).and.(dist3.ge.zero)) then ! in the melt>0
        if (abs(dist2).lt.abs(dist3)) then
         dist4=abs(dist2)
        else
         dist4=abs(dist3)
        endif
       else
        print *,"dist bust"
        stop
       endif
       dist=max(dist,dist4) !find the shortest distance to the melt
      enddo ! i=1..nSphere

      return
      end subroutine INIT_LS_LIQUID_MELT


       ! probtype.eq.299
      subroutine INIT_LS_GAS_MELT(x,y,z,time,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,z,time,dist,dist2
      integer nSphere
      real(amrex_real) xSphere
      real(amrex_real) ySphere
      real(amrex_real) delta_sphere,hugedist
      integer i
 
      hugedist=1.0D+10
 
      if (num_materials.ne.3) then
       print *,"num_materials invalid in INIT_LS_GAS_MELT"
       stop
      endif

      if (SDIM.eq.2) then
       if (abs(z-y).gt.VOFTOL) then
        print *,"z=y required in INIT_LS_GAS_MELT"
        print *,"x,y,z = ",x,y,z
        print *,"probtype,axis_dir ",probtype,axis_dir
        stop
       endif
      endif

      nSphere=3

      if (axis_dir.eq.0) then
       delta_sphere=two
      else if (axis_dir.eq.1) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.2) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.3) then
       delta_sphere=two
      else
       print *,"axis_dir invalid init ls gas melt"
       stop
      endif

      dist=-hugedist
      do i=1,nSphere
       xSphere=xblob+(i-one)*delta_sphere
       ySphere=yblob
       dist2=radblob2-sqrt((x-xSphere)**2+(y-ySphere)**2) 
       dist=max(dist,dist2)
      enddo
      dist=-dist

      return
      end subroutine INIT_LS_GAS_MELT


       ! probtype==301
      subroutine INIT_LS_SOLID_AM(x,y,z,time,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,z,time,dist,dist2
      integer nSphere
      real(amrex_real) xSphere
      real(amrex_real) ySphere
      real(amrex_real) delta_sphere,hugedist
      integer i
 
      hugedist=1.0D+10
 
      if (num_materials.ne.3) then
       print *,"num_materials invalid in INIT_LS_SOLID_AM"
       stop
      endif

      if (SDIM.eq.2) then
       if (abs(z-y).gt.VOFTOL) then
        print *,"z=y required in INIT_LS_SOLID_AM"
        print *,"x,y,z = ",x,y,z
        print *,"probtype,axis_dir ",probtype,axis_dir
        stop
       endif
      endif

      nSphere=3

      if (axis_dir.eq.0) then
       delta_sphere=two
      else if (axis_dir.eq.1) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.2) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.3) then
       delta_sphere=two
      else
       print *,"axis_dir invalid init ls solid AM"
       stop
      endif

      dist=-hugedist

      do i=1,nSphere
       xSphere=xblob+(i-one)*delta_sphere
       ySphere=yblob
       dist2=radblob-sqrt((x-xSphere)**2+(y-ySphere)**2) 
       dist=max(dist,dist2)
      enddo

      return
      end subroutine INIT_LS_SOLID_AM

       ! probtype.eq.301
      subroutine INIT_LS_LIQUID_AM(x,y,z,time,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,z,time,dist
      integer nSphere
      real(amrex_real) delta_sphere,hugedist
 
      hugedist=1.0D+10
 
      if (num_materials.ne.3) then
       print *,"num_materials invalid in INIT_LS_LIQUID_AM"
       stop
      endif
      if (SDIM.eq.2) then
       if (abs(z-y).gt.VOFTOL) then
        print *,"z=y required in INIT_LS_LIQUID_AM"
        print *,"x,y,z = ",x,y,z
        print *,"probtype,axis_dir ",probtype,axis_dir
        stop
       endif
      endif

      nSphere=3

      if (axis_dir.eq.0) then
       delta_sphere=two
      else if (axis_dir.eq.1) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.2) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.3) then
       delta_sphere=two
      else
       print *,"axis_dir invalid init ls liquid melt"
       stop
      endif

      dist=-hugedist

      return
      end subroutine INIT_LS_LIQUID_AM


       ! probtype.eq.301
      subroutine INIT_LS_GAS_AM(x,y,z,time,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,z,time,dist,dist2
      integer nSphere
      real(amrex_real) xSphere
      real(amrex_real) ySphere
      real(amrex_real) delta_sphere,hugedist
      integer i
 
      hugedist=1.0D+10
 
      if (num_materials.ne.3) then
       print *,"num_materials invalid in INIT_LS_GAS_AM"
       stop
      endif

      if (SDIM.eq.2) then
       if (abs(z-y).gt.VOFTOL) then
        print *,"z=y required in INIT_LS_GAS_AM"
        print *,"x,y,z = ",x,y,z
        print *,"probtype,axis_dir ",probtype,axis_dir
        stop
       endif
      endif

      nSphere=3

      if (axis_dir.eq.0) then
       delta_sphere=two
      else if (axis_dir.eq.1) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.2) then
       delta_sphere=32.0e-6
      else if (axis_dir.eq.3) then
       delta_sphere=two
      else
       print *,"axis_dir invalid init ls gas AM"
       stop
      endif

      dist=-hugedist
      do i=1,nSphere
       xSphere=xblob+(i-one)*delta_sphere
       ySphere=yblob
       dist2=radblob-sqrt((x-xSphere)**2+(y-ySphere)**2) 
       dist=max(dist,dist2)
      enddo
      dist=-dist

      return
      end subroutine INIT_LS_GAS_AM

       ! soliddist is called by: subroutine materialdistsolid 
       ! dist<0 in the solid
       ! dist>0 in the fluid 
      subroutine soliddist(x,y,z,dist,im)
      use global_utility_module
      use bubbleControl_module

      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: x,y,z
      real(amrex_real), intent(out) :: dist
      real(amrex_real) dist1,temprad
      real(amrex_real) rr,xx,yy,zz,aspect,offset
      real(amrex_real) xlarge
      integer igeom
      real(amrex_real) costheta,sintheta
      real(amrex_real) xprime,yprime
      real(amrex_real) zprime
      real(amrex_real) radcross,dist2,radpt
      real(amrex_real) hugedist,dist3,dist4
      real(amrex_real) radx,radshrink
      real(amrex_real) pipexlo,pipexhi
      real(amrex_real) zmin,zmax
      real(amrex_real) angle_x,angle_y
      real(amrex_real), PARAMETER :: stub_zero=zero
      real(amrex_real), PARAMETER :: stub_five=five
      real(amrex_real), PARAMETER :: stub_one=one
      real(amrex_real), PARAMETER :: stub_eleven=11.0d0
      real(amrex_real) half_yblob
      real(amrex_real) half_zblob2
      real(amrex_real) factor_zblob
      real(amrex_real) diamblob
      integer i,j,iSphere

      if (num_materials.lt.1) then
       print *,"num_materials invalid in soliddist"
       stop
      endif

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid11"
       stop
      endif

      if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then
       ! do nothing
      else
       print *,"expecting FSI_flag(im)=FSI_PRESCRIBED_PROBF90 in soliddist"
       print *,"im,FSI_flag(im): ",im,FSI_flag(im)
       stop
      endif

      if (SDIM.eq.2) then
       if (abs(z-y).gt.VOFTOL) then
        print *,"z=y required in soliddist"
        print *,"x,y,z = ",x,y,z
        print *,"probtype,axis_dir ",probtype,axis_dir
        stop
       endif
      endif

      igeom=1

      hugedist=99999.0

      dist=hugedist
     
      if (probtype.eq.540) then  ! Rieber simulation no solid
       dist=99999.0
! airblast without nozzle: set soliddist
      else if (probtype.eq.529) then 
       dist=99999.0
      else if (((probtype.eq.531).or. &   ! soliddist
                (probtype.eq.390).or. &
                (probtype.eq.536).or. &
                (probtype.eq.537).or. &
                (probtype.eq.538).or. &
                (probtype.eq.541)).and. &
               (SDIM.eq.3)) then
       dist=99999.0

      else if (probtype.eq.102) then
! yblob3 = gas nozzle len
       if (yblob3.le.zero) then
        print *,"yblob3 invalid"
        stop
       endif
       if (radblob5.le.radblob2) then
        print *,"radblob5 invalid"
        stop
       endif
       if (radblob.gt.radblob2) then
        print *,"radblob2 invalid"
        stop
       endif
       if (y.gt.yblob3) then  ! expansion region
        radshrink=radblob7**2-radblob5**2
        if (radshrink.le.zero) then
         print *,"radshrink invalid"
         stop
        endif
        radshrink=sqrt(radshrink)
        radx=radblob4+(y-yblob3)*(radshrink-radblob4)/(probhiy-yblob3)
        dist=abs(x-(radx+half*radblob6))-half*radblob6
       else
        radx=radblob3+y*(radblob4-radblob3)/yblob3  ! radblob3=entry gas rad
        dist=abs(x-(radx+half*radblob6))-half*radblob6
         ! yblob=primary liquid nozzle len  yblob2=secondary len
        if (y.gt.yblob+yblob2) then 
          ! radblob2=secondary radius
         if (x.lt.radblob2) then
          dist1=sqrt( (x-radblob2)**2+(y-yblob-yblob2)**2 )
          ! radblob5=outer liquid nozzle radius
         else if (x.gt.radblob5) then
          dist1=sqrt( (x-radblob5)**2+(y-yblob-yblob2)**2 )
         else
          dist1=abs(y-yblob-yblob2)
         endif
        else if (y.gt.yblob) then
         dist1=abs(x-half*(radblob5+radblob2))-half*(radblob5-radblob2)
        else
         dist1=abs(x-half*(radblob5+radblob))-half*(radblob5-radblob)
        endif
        if (abs(dist1).lt.abs(dist)) then
         dist=dist1
        endif  
       endif
      else if (probtype.eq.42) then ! bubble jetting (soliddist)

       !dist<0 in the solid
       !soliddist is called by: subroutine materialdistsolid 
       call jetting_plate_dist(x,y,z,dist)

       ! soliddist
      else if (probtype.eq.bubbleInPackedColumn) then

       !open(11,file="sphere2d.dat",form="formatted")
       !read(11,*) nSphere
       nSphere=0
       do i=1,nSphere
          read(11,*) j, xSphere(i),ySphere(i),rSphere(i)
       enddo
       !close(11)

       dist=hugedist
       do iSphere=1,nSphere
          dist2=sqrt((x-xSphere(iSphere))**2+(y-ySphere(iSphere))**2) &
                     -rSphere(iSphere)
          dist=min(dist,dist2)
       enddo 
       ! in soliddist
      else if ((probtype.eq.46).and.(radblob2.gt.zero)) then
       print *,"the airgun problem needs to be redefined"
       stop
      else if ((probtype.eq.63).and.(SDIM.eq.3)) then
       call nozzlerad(z,radcross,stub_zero)
       radpt=sqrt(x**two+y**two)
       if (z.le.(two*xblob10)) then
        dist=radcross-radpt
       else if (radpt.ge.(half*xblob10)) then
        dist=z-two*xblob10
       else if (radpt.le.(half*xblob10)) then
        dist=sqrt((radcross-radpt)**two+(z-two*xblob10)**two)
       endif

      else if ((probtype.eq.63).and.(SDIM.eq.2)) then
        call nozzlerad(z,radcross,stub_zero)
  
        if (z.LE.2*xblob10) THEN
         dist=radcross-x
        elseif (x.LE.xblob10/2) THEN
         dist=dsqrt((x-xblob10/2)**two+(z-2*xblob10)**two)
        else 
         dist=z-2*xblob10
        endif
           
      else if (probtype.eq.64) then
       call nozzlerad(z,radcross,radblob)
  
       if (z.LE.(two*xblob10-radblob)) then
        dist=radcross-x
       elseif (x.LE.(xblob10/two+radblob)) then
        dist=dsqrt((x-(xblob10/two+radblob))**two+ &
        (z-(two*xblob10-radblob))**two)-radblob
       else 
        dist=z-two*xblob10
       endif

      else if ((probtype.eq.43).and.(zblob.gt.zero)) then ! soliddist
       half_yblob=half*yblob
       call squaredist(x,y,-xblob,xblob,half_yblob,half_yblob+zblob,dist) 
      else if ((probtype.eq.43).and.(radblob.gt.zero)) then
       offset=yblob/two+xblob/two
!      offset=1.1*yblob/two
       call squaredist(x,y,-radblob,radblob,-offset,offset,dist)
      else if (probtype.eq.48) then
       xlarge=1.0e+10
       call squaredist(x,y,xblob,xlarge,-yblob,yblob,dist)
        ! nozzle for bubble formation
      else if ((probtype.eq.25).and.(axis_dir.gt.0)) then ! soliddist
       if (zblob.le.zero) then
        print *,"zblob (nozzle ht) must be positive"
        stop
       endif
       aspect=radblob+radblob2
        ! dist<0 in solid
       factor_zblob=-1000.0d0*zblob
       call squaredist(x,y,radblob,aspect,factor_zblob,zblob,dist)
      else if (probtype.eq.62) then
       costheta=cos(xblob2)
       sintheta=sin(xblob2)
       xprime=costheta*(x-xblob)-sintheta*(y-yblob)
       yprime=sintheta*(x-xblob)+costheta*(y-yblob)
       half_zblob2=half*zblob2
       call squaredist(xprime,yprime,radblob,radblob2, &
         -half_zblob2,half_zblob2,dist) 
       call squaredist(xprime,yprime,-radblob2,-radblob, &
         -half_zblob2,half_zblob2,dist1) 
       if (dist1.lt.dist) then
        dist=dist1
       endif
      else if (probtype.eq.22) then
       call jetgeom(x,y,z,dist)
      else if (probtype.eq.35) then
       if (SDIM.eq.3) then
        rr=sqrt( (x-half*xblob)**2+(y-half*xblob)**2 )
        zz=z
       else
        rr=abs(x)
        zz=y
       endif
       aspect=Pi*30.0/180.0
       offset=zz*tan(aspect)
       dist=(radblob-offset)-rr
      else if ((probtype.eq.44).and.(axis_dir.eq.2)) then
        ! in 2d, y=z
       call damdist(x,z,dist,igeom)
      else if ((probtype.eq.30).or.(probtype.eq.32).or. &
               (probtype.eq.33).or.(probtype.eq.34).or. &
               (probtype.eq.21).or. &
               (probtype.eq.bubbleInPackedColumn)) then
       call selectgeom(x,y,z,dist)  ! in soliddist
      else if (probtype.eq.50) then
       print *,"this option obsolete"
       stop
      else if (probtype.eq.52) then
       print *,"obsolete"
       stop
        ! flapping wing, soliddist
      else if (probtype.eq.701) then 
       print *,"soliddist should not be called when probtype=701"
       stop

       ! 2D or 3D
      else if ((probtype.eq.538).or.(probtype.eq.541)) then
       dist=99999.0
      else if (probtype.eq.539) then  ! soliddist
       dist=99999.0
      else if (probtype.eq.53) then  ! dist>0 in fluid  soliddist
       if (axis_dir.eq.0) then ! no nozzle
        dist=9999.0
       else if ((axis_dir.eq.1).and.(SDIM.eq.2)) then
        dist1=sqrt( (x-xblob2)**2 + (y-yblob2)**2) -radblob2-radblob4
        dist2=sqrt( (x-xblob2)**2 + (y-yblob2)**2) -radblob2
        dist3=sqrt( (x-xblob2)**2 + (y-yblob2)**2) -radblob3
        dist4=sqrt( (x-xblob2)**2 + (y-yblob2)**2) -radblob3+radblob4

        if (x.ge.xblob2) then
         if (y.ge.radblob2+radblob4) then
          dist=sqrt( (x-xblob2)**2+(y-radblob2-radblob4)**2)
         else if (y.ge.radblob2) then
          dist=x-xblob2
         else if (y.ge.radblob3) then
          dist=sqrt( (x-xblob2)**2+(y-radblob3)**2)
         else if (y.ge.radblob3-radblob4) then
          dist=x-xblob2
         else
          dist=sqrt( (x-xblob2)**2+(y-radblob3+radblob4)**2)
         endif
        else if (dist1.ge.zero) then
         dist=dist1
        else if (dist2.ge.zero) then
         if (abs(dist1).lt.abs(dist2)) then
          dist=dist1
         else
          dist=-dist2
         endif
        else if (dist3.ge.zero) then
         if (abs(dist2).lt.abs(dist3)) then
          dist=-dist2
         else
          dist=dist3
         endif
        else if (dist4.ge.zero) then
         if (abs(dist3).lt.abs(dist4)) then
          dist=dist3
         else
          dist=-dist4
         endif
        else
         dist=-dist4
        endif

         ! pressure bc, cylindrical nozzle
         ! want dist>0 in fluid
       else if ((axis_dir.eq.2).and.(SDIM.eq.2)) then 
        zmin=-100.0*radblob
        zmax=3.0*radblob

! dist>0 outside square
        call squaredist(x,y,xblob-radblob,xblob+radblob,zmin,zmax,dist)
        diamblob=two*radblob
        call squaredist(x,y,xblob-diamblob,xblob+diamblob, &
         zmin,zmax,dist1)
        if (dist.le.zero) then
         dist=-dist
        else if (dist1.ge.zero) then
         dist=dist1
        else if (abs(dist).lt.abs(dist1)) then
         dist=-abs(dist)
        else
         dist=-abs(dist1)
        endif
       else if ((axis_dir.eq.1).or. &
                (axis_dir.eq.2)) then
        zmin=-100.0*radblob
        zmax=3.0*radblob

! dist>0 outside cylinder
        call cylinderdist(x,y,z,xblob,yblob,radblob,zmin,zmax,dist)
        diamblob=two*radblob
        call cylinderdist(x,y,z,xblob,yblob,diamblob,zmin,zmax,dist1)
        if (dist.le.zero) then
         dist=-dist
        else if (dist1.ge.zero) then
         dist=dist1
        else if (abs(dist).lt.abs(dist1)) then
         dist=-abs(dist)
        else
         dist=-abs(dist1)
        endif
       else if (axis_dir.eq.100) then  ! fan nozzle
        dist=99999.0
       else
        print *,"axis_dir invalid probtype=53"
        stop
       endif
         
      else if (probtype.eq.56) then
       print *,"obsolete"
       stop
      else if (probtype.eq.54) then
       print *,"obsolete"
       stop
        ! ship wave no solid prescribed here
      else if (probtype.eq.9) then  
       dist=99999.0
      else if ((probtype.eq.45).and.(axis_dir.eq.1)) then
       print *,"obsolete"
       stop
      else if ((probtype.eq.39).and.(radblob2.gt.zero).and. &
               (zblob2.gt.zero).and.(radblob2.lt.half*xblob)) then
       xx=abs(x)
       if (y.ge.zblob2) then
        dist=sqrt(xx**2+(y-zblob2)**2)
       else
        yy=radblob2*(zblob2-y)/zblob2
        dist=xx-yy
       endif
! natural convection in triangular enclosure (soliddist)
      else if (probtype.eq.81) then
       dist=y+yblob3-yblob*(one-x/xblob)
! sphere impact on flat surface (dist >0 in fluid)
! soliddist - falling sphere
      else if ((probtype.eq.531).and.(SDIM.eq.2)) then  
       if (axis_dir.eq.0) then
        dist=sqrt( (x-xblob)**2+(y-yblob)**2 ) -radblob
       else if (axis_dir.eq.1) then ! it should never come to this
        dist=99999.0
       else if (axis_dir.eq.2) then ! it should never come to this
        dist=99999.0
       else if (axis_dir.eq.3) then ! it should never come to this
        dist=99999.0
       else
        print *,"axis_dir invalid probtype=531"
        stop
       endif
       ! dist<0 in the solid
       ! dist>0 in the fluid 
       ! in: soliddist, turbulent cylindrical pipe
      else if ((probtype.eq.41).and.(axis_dir.eq.5)) then
       if (SDIM.eq.2) then
        dist=radblob-sqrt((y-yblob)**2)
       else if (SDIM.eq.3) then
        dist=radblob-sqrt((y-yblob)**2+(z-zblob)**2) 
       else
        print *,"dimension bust"
        stop
       endif 
      else if ((probtype.eq.41).and. &
               (axis_dir.eq.0).and.(SDIM.eq.3)) then
! x is free stream direction
       if (levelrz.eq.COORDSYS_CARTESIAN) then 
        dist=zblob2-sqrt(y**2+z**2)
       else
        print *,"levelrz invalid soliddist"
        stop
       endif

! pipe problem  soliddist: dist>0 fluid
      else if ((probtype.eq.41).and. &
               (SDIM.eq.2)) then  

         ! axis_dir=4 comparison with LSA
       if (axis_dir.eq.4) then
        dist=99999.0
       else
        pipexlo=problox
        pipexhi=probhix
        if ((axis_dir.eq.1).or.(axis_dir.eq.2)) then
         pipexlo=zero
         pipexhi=two*radblob3
        endif
 
        if (x.lt.pipexlo) then
         dist=pipexlo-x
        else if (x.gt.pipexhi) then
         dist=x-pipexhi
        else
         dist=-min( x-pipexlo, pipexhi-x )
        endif
        dist=-dist
       endif  ! axis_dir<> 4

! soliddist: dist>0 in fluid 2d or 3d
      else if (probtype.eq.59) then  ! inputs.block_ice_melt

       ! dist>0 in the substrate
       angle_x=radblob2
       angle_y=zero
       call ice_substrate_distance(x,y,z,angle_x,angle_y,dist)
       ! now make dist<0 in the substrate.
       dist=-dist

      else if ((probtype.eq.54).and.(SDIM.eq.3)) then
       zmin=zblob
       zmax=zblob+half*zblob
       temprad=(one+half)*radblob
       call cylinderdist(x,y,z,xblob,yblob,temprad,zmin,zmax,dist)
      else if ((probtype.eq.65).and. &
               (SDIM.eq.3)) then     
       costheta=cos(xblob2)
       sintheta=sin(xblob2)
       xprime=costheta*(x-xblob)-sintheta*(z-zblob)
       yprime=y-yblob
       zprime=sintheta*(x-xblob)+costheta*(z-zblob)
       half_zblob2=half*zblob2
       call tcylinderdist(xprime,yprime,zprime,stub_zero,stub_zero,radblob, &
        -half_zblob2,half_zblob2,dist)
       call tcylinderdist(xprime,yprime,zprime,stub_zero,stub_zero,radblob2, &
        -half_zblob2,half_zblob2,dist1)
       if (dist1.ge.zero) then
        dist=dist1
       else if ((dist.ge.zero).and.(dist1.le.zero)) then
        if (-dist1.lt.dist) then
         dist=dist1
        else
         dist=-dist
        endif
       else
         dist=-dist
       endif
!       lower cylinder first, center at (5,5,5)
       call cylinderdist(x,y,z,stub_five,stub_five,radblob4, &
               -stub_one,stub_eleven,dist1)
       dist1 = -dist1
       if (dist.gt.dist1) then
        dist = dist1
       endif

      else if ((probtype.eq.62).and.(SDIM.eq.3)) then
       costheta=cos(xblob2)
       sintheta=sin(xblob2)
       xprime=costheta*(x-xblob)-sintheta*(z-zblob)
       yprime=y-yblob
       zprime=sintheta*(x-xblob)+costheta*(z-zblob)
       half_zblob2=half*zblob2
       call cylinderdist(xprime,yprime,zprime,stub_zero,stub_zero,radblob, &
        -half_zblob2,half_zblob2,dist)
       call cylinderdist(xprime,yprime,zprime,stub_zero,stub_zero,radblob2, &
        -half_zblob2,half_zblob2,dist1)
       if (dist1.ge.zero) then
        dist=dist1
       else if ((dist.ge.zero).and.(dist1.le.zero)) then
        if (-dist1.lt.dist) then
         dist=dist1
        else
         dist=-dist
        endif
       else
         dist=-dist
       endif

      else if (probtype.eq.110) then ! pos distance in fluid
       call get_bottom_distIOWA(x,y,dist)
      else if (probtype.eq.5700) then  ! in soliddist
       if ((axis_dir.eq.0).or.(axis_dir.eq.1).or. &
           (axis_dir.eq.2).or.(axis_dir.eq.3)) then
        print *,"must have FSI_flag=FSI_PRESCRIBED_PROBF90"
        stop
       else if (axis_dir.eq.4) then
        call square_Tchannel_dist(x,y,z,dist)  ! dist>0 in solid
        dist=-dist
       else if (axis_dir.eq.5) then
        call squeeze_channel_dist(x,y,z,dist)  ! dist>0 in solid
        dist=-dist
       else
        print *,"axis_dir invalid probtype=5700"
        stop
       endif
      endif
   
      return
      end subroutine soliddist


      subroutine nozzlerad(zval,radcross,rounded)
      IMPLICIT NONE

      real(amrex_real) zval,radcross,rounded

      if (zval.le.xblob10) then
       radcross=xblob10
      else if (zval.le.(three*xblob10/two)) then
       radcross=-zval+two*xblob10
      else if (zval.le.two*xblob10-rounded) then
       radcross=xblob10/two
      else
       radcross=xblob10/two
      endif

      return
      end subroutine nozzlerad


      subroutine jetgeom(x,y,z,dist)
      IMPLICIT NONE
      real(amrex_real) x,y,z,dist
      real(amrex_real) NPT,HSB,NID,NOD,CHH,scaleCHH
      real(amrex_real) xx1(maxnumline),yy1(maxnumline)
      real(amrex_real) xx2(maxnumline),yy2(maxnumline)
      integer dd(maxnumline),lessflag(maxnumline),numline


      call initjetparms(HSB,NOD,NPT,NID,CHH,scaleCHH)

      xx1(1)=half*NOD
      xx2(1)=1.0D+10
      yy1(1)=HSB+NPT
      yy2(1)=HSB+NPT
      dd(1)=0
      lessflag(1)=1

      yy1(2)=HSB
      yy2(2)=HSB+NPT
      xx1(2)=half*NID
      xx2(2)=half*NOD
      dd(2)=1
      lessflag(2)=0

      xx1(3)=half*NID
      xx2(3)=half*scaleCHH
      yy1(3)=HSB
      yy2(3)=HSB
      dd(3)=0
      lessflag(3)=0
      if (xx1(3).gt.xx2(3)) then
       xx1(3)=half*scaleCHH
       xx2(3)=half*NID
       lessflag(3)=1
      endif
      
      yy1(4)=zero
      yy2(4)=HSB
      xx1(4)=half*scaleCHH
      xx2(4)=half*scaleCHH
      dd(4)=1
      lessflag(4)=0
    
      numline=4

      if ((axis_dir.ge.8).and.(axis_dir.le.10)) then
       print *,"this option disabled"
       stop
      else if ((axis_dir.eq.11).or.(axis_dir.eq.12).or. &
               (axis_dir.eq.13)) then
       call microfabgeom(x,y,dist)
      else
       call construct(x,y,xx1,yy1,xx2,yy2,dd,lessflag,numline,dist)
       if ((x.le.half*NOD).and.(x.le.half*NID)) then
        dist=abs(dist)
       endif
       if (y.ge.HSB+NPT) then
        dist=abs(dist)
       endif
       if ((x.ge.half*scaleCHH).and.(y.le.HSB)) then
        dist=-abs(dist)
       endif
       if ((x.ge.half*NOD).and.(x.ge.half*NID).and. &
           (y.le.HSB+NPT).and.(y.ge.HSB)) then
        dist=-abs(dist)
       endif
      endif

      return
      end subroutine jetgeom


      subroutine initjetparms(HSB,NOD,NPT,NID,CHH,scaleCHH)
      IMPLICIT NONE
      real(amrex_real) NPT,HSB,NID,NOD,CHH,scaleCHH


! get rid of floating exceptions !
      HSB=zero
      NOD=zero
      NPT=zero
      NID=zero
      CHH=zero
      scaleCHH=zero

      if (axis_dir.eq.0) then
       NPT=55.0
       NOD=23.5
       NID=41.0
       CHH=30.0
       scaleCHH=CHH*sqrt(four/Pi)
       HSB=50.0
      else if (axis_dir.eq.1) then
       NPT=30.0
       NOD=21.5
       NID=31.0
       CHH=30.0
       scaleCHH=CHH*sqrt(four/Pi)
       HSB=50.0
      else if (axis_dir.eq.2) then
       NPT=18.0
       NOD=21.0
       NID=26.7
       CHH=30.0
       scaleCHH=CHH*sqrt(four/Pi)
       HSB=50.0
      else if (axis_dir.eq.3) then
       NPT=18.0
       NOD=20.0
       NID=32.0
       CHH=30.0
       scaleCHH=CHH*sqrt(four/Pi)
       HSB=50.0
      else if (axis_dir.eq.4) then
       NPT=18.0
       NOD=28.0
       NID=19.0
       CHH=30.0
       scaleCHH=CHH*sqrt(four/Pi)
       HSB=50.0
      else if (axis_dir.eq.5) then
       NPT=18.0
       NOD=34.0
       NID=19.0
       CHH=30.0
       scaleCHH=CHH*sqrt(four/Pi)
       HSB=50.0
      else if ((axis_dir.eq.6).or.(axis_dir.eq.14)) then
       NPT=20.0
       HSB=30.0
       CHH=40.0
       scaleCHH=40.0
       NOD=25.0
       NID=40.0
      else if (axis_dir.gt.13) then
       print *,"axis_dir out of range in initjetparms!!!"
       stop
      endif

      return
      end subroutine initjetparms


      subroutine damdist(x,z,dist,igeom)
      use global_utility_module
      IMPLICIT NONE

      real(amrex_real) x,z,dist,y1,y2,xctr,yswap,xswap
      integer igeom
      real(amrex_real) x1,x2


      if (axis_dir.eq.4) then ! water if x<xblob y<yblob
       x1=-1.0D+10
       x2=xblob
       y1=-1.0D+10
       y2=yblob
       call squaredist(x,z,x1,x2,y1,y2,dist)
       dist=-dist
      else if (axis_dir.eq.3) then
       dist=-x
      else if ((axis_dir.eq.0).or.(axis_dir.eq.1).or.(axis_dir.eq.2)) then
       xswap=x
       if (axis_dir.eq.0) then
        xctr=half*xblob
        y1=yblob/four
        y2=three*yblob/four
       else if (axis_dir.eq.1) then
        xctr=xblob/four
        y1=-999999.0
        y2=half*yblob
       else if (axis_dir.eq.2) then
        if (igeom.eq.0) then
         xctr=xblob/two
         y1=two*yblob
         y2=four*yblob 
        else if (igeom.eq.1) then
         xctr=xblob/two
         y1=yblob
         y2=-999999.0
        else
         print *,"igeom invalid"
         stop
        endif
       else
         print *,"axis_dir out of range for dambreak problem"
         stop
       endif
   
       if (y1.gt.y2) then
        yswap=y1
        y1=y2
        y2=yswap
        xctr=xblob-xctr
        xswap=xblob-xswap
       endif
 
       if ((z.le.y1).or.((z.le.y2).and.(xswap.le.xctr))) then
        if (xswap.gt.xctr) then
         dist=z-y1
        else if (z.gt.y1) then
         dist=z-y2
         if (dist.lt.xswap-xctr) then
          dist=xswap-xctr
         endif        
        else
         dist=-sqrt((xswap-xctr)**2+(z-y1)**2)
         if (dist.lt.z-y2) then
          dist=z-y2
         endif
        endif
       else if (xswap.lt.xctr) then
        dist=z-y2
       else if (z.lt.y2) then
        dist=z-y1
        if (dist.gt.xswap-xctr) then
         dist=xswap-xctr
        endif  
       else
        dist=sqrt((xswap-xctr)**2+(z-y2)**2)
        if (dist.gt.z-y1) then
         dist=z-y1
        endif
       endif
       if (igeom.eq.0) then
        dist=-dist
       endif
      else
       print *,"axis_dir invalid damdist"
       stop
      endif

      return
      end subroutine damdist

      subroutine getslopeparms(xstart,ystart,zstart,xend,yend,zend)
      IMPLICIT NONE
      real(amrex_real) xstart,xend,ystart,yend,zstart,zend

      xstart=0.0
      xend=4.0
      ystart=1.0
      yend=1.0
      zstart=0.0
      zend=4.0

      return
      end subroutine getslopeparms


!       tapered cylinder
      subroutine tcylinderdist(x,y,z,xcen,ycen,rad,zmin,zmax,dist)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x,y,z,xcen,ycen,rad
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real), INTENT(in) :: zmin,zmax

      if (zmin.ge.zmax-EPS10) then
       print *,"invalid parameters ",zmin,zmax
       stop
      endif
      dist=sqrt((x-xcen)**2+(y-ycen)**2/9.0)-rad
      if (z.ge.zmax) then
       if (dist.le.zero) then
        dist=z-zmax
       else
        dist=sqrt(dist**2+(z-zmax)**2)
       endif
      else if (z.le.zmin) then
       if (dist.le.zero) then
        dist=zmin-z
       else
        dist=sqrt(dist**2+(zmin-z)**2)
       endif
      endif

      return
      end subroutine tcylinderdist


! dist is positive in the fluid
      subroutine get_bottom_distIOWA(x,y,dist)

      IMPLICIT NONE

      real(amrex_real) x,y,height,dist
      real(amrex_real) h,l


      if (probtype.ne.110) then
       print *,"probtype invalid get_bottom_distIOWA"
       stop
      endif

       ! inlet x/H=-52 outlet x/H=44
       ! inlet x=-52H=594.36 cm  outlet x=44H=502.92 cm
       ! zhi=5H=57.15 cm
       ! (x/l)=(x/(2.5 h))=2/5 (x/h)
       ! at inlet, (x/l)=-52 * 2/5 = -20.8
       ! at outlet, (x/l)=44 * 2/5 = 17.6
       ! 0<z/h<5
       ! dx_min=0.3mm  dz_min=0.2mm ?
       ! -1<x/l<1
      h=11.43  ! h=H  maximum bump height in centimeters
      l=2.5*h
      if (abs(x/l).le.one) then
       height=h*(one-two*((x/l)**2)+(x/l)**4)
       dist=y-height
      else if (x/l.lt.zero) then
       dist=sqrt((x+l)**2+y**2)
      else if (x/l.gt.zero) then
       dist=sqrt((x-l)**2+y**2)
      else
       print *,"bust"
       stop
      endif

      return
      end subroutine get_bottom_distIOWA

         ! dist>0 outside the channel
      subroutine square_Tchannel_dist(x,y,z,dist)
      IMPLICIT NONE
 
      real(amrex_real) x,y,z,zmid,dist,distz
      real(amrex_real) radchannel,xcen1,xcen2,xcen3,ychannel,ycen


      if ((probtype.eq.5700).and.(axis_dir.eq.4)) then
       if ((xblob2.lt.xblob3).and.(yblob2.lt.yblob3).and. &
           (zblob2.lt.zblob3)) then
       
        if (SDIM.eq.3) then
         zmid=z
        else if (SDIM.eq.2) then
         zmid=half*(zblob2+zblob3)
        else
         print *,"dimension bust"
         stop
        endif
 
        ! first find 2d distance in x-y plane

        if ((xblob4.le.zero).or.(x.le.xblob4)) then

         if (y.le.yblob2) then
          dist=yblob2-y
         else if (y.le.half*(yblob2+yblob3)) then
          dist=yblob2-y
         else if (y.le.yblob3) then
          if (x.le.xblob2) then
           dist=y-yblob3
          else if (x.ge.xblob3) then
           dist=y-yblob3
          else if (x.le.half*(xblob2+xblob3)) then
           dist=-sqrt( (y-yblob3)**2+(x-xblob2)**2 )
          else if (x.ge.half*(xblob2+xblob3)) then
           dist=-sqrt( (y-yblob3)**2+(x-xblob3)**2 )
          else
           print *,"x bust"
           stop
          endif
         else if ((x.ge.xblob2).and.(x.le.half*(xblob2+xblob3)).and. &
                  (y.ge.yblob3)) then
          dist=xblob2-x
         else if ((x.le.xblob3).and.(x.ge.half*(xblob2+xblob3)).and. &
                  (y.ge.yblob3)) then
          dist=x-xblob3
         else if ((x.le.xblob2).and.(y.ge.yblob3)) then
 
          if (y-yblob3.le.xblob2-x) then
           dist=y-yblob3
          else
           dist=xblob2-x
          endif
 
         else if ((x.ge.xblob3).and.(y.ge.yblob3)) then

          if (y-yblob3.le.x-xblob3) then
           dist=y-yblob3
          else
           dist=x-xblob3
          endif

         else
          print *,"parameter bust"
          stop
         endif 

          ! dist>0 outside channel
        else if (x.ge.xblob4) then
         radchannel=half*(yblob3-yblob2)
         xcen1=xblob4
         xcen2=xblob4+two*radblob4
         xcen3=xblob4+four*radblob4
         ychannel=yblob2+radchannel
         ycen=ychannel-radblob4
         if ((x.le.xcen2).and.(y.ge.ycen)) then
          dist=sqrt( (x-xcen1)**2+(y-ycen)**2 )
          if (dist.ge.radblob4) then
           dist=dist-radblob4-radchannel
          else
           dist=radblob4-radchannel-dist
          endif
         else if ((x.le.xcen3).and.(y.le.ycen)) then
          dist=sqrt( (x-xcen2)**2+(y-ycen)**2 )
          if (dist.ge.radblob4) then
           dist=dist-radblob4-radchannel
          else
           dist=radblob4-radchannel-dist
          endif
         else if ((x.le.xcen3).and.(y.ge.ycen)) then
          dist=sqrt( (x-xcen3)**2+(y-ycen)**2 )
          if (dist.ge.radblob4) then
           dist=dist-radblob4-radchannel
          else
           dist=radblob4-radchannel-dist
          endif
         else if (x.ge.xcen3) then
          if (y.le.ychannel) then
           dist=ychannel-y-radchannel
          else
           dist=y-ychannel-radchannel
          endif
         else
          print *,"channel dimension bust"
          stop
         endif

        else
         print *,"xblob4 invalid"
         stop
        endif

        if (zmid.le.zblob2) then
         distz=zblob2-zmid
        else if (zmid.ge.zblob3) then
         distz=zmid-zblob3
        else if (zmid.le.half*(zblob2+zblob3)) then
         distz=zblob2-zmid
        else if (zmid.ge.half*(zblob2+zblob3)) then
         distz=zmid-zblob3
        else
         print *,"bust for distz"
         stop
        endif
        
        if (distz.ge.zero) then
         if (dist.ge.zero) then
          dist=sqrt( dist**2 +distz**2 )
         else if (dist.le.zero) then
          dist=distz
         else
          print *,"dist invalid square_Tchannel_dist"
          stop
         endif
        else if (distz.le.zero) then
         if (dist.ge.zero) then
          ! do nothing
         else if (dist.le.zero) then
          if (abs(dist).le.abs(distz)) then
           ! do nothing
          else if (abs(dist).ge.abs(distz)) then
           dist=distz
          else
           print *,"dist invalid square_Tchannel_dist"
           stop
          endif
         else
          print *,"dist invalid square_Tchannel_dist"
          stop
         endif
        else
         print *,"distz invalid"
         stop
        endif
 
       else
        print *,"dimensions incorrect"
        stop
       endif
      else
       print *,"probtype or axis_dir incorrect"
       stop
      endif

      return
      end subroutine square_Tchannel_dist

      subroutine selectgeom(x,y,z,dist)
      use bubbleControl_module

      IMPLICIT NONE
      real(amrex_real) x,y,z,dist,dist1,dist2
      real(amrex_real) xstart,xend,ystart,yend,slope,intercept,xint,yint
      real(amrex_real) zstart,zend
      real(amrex_real) ylen,frontrad
      integer igeom,iSphere

      if (SDIM.eq.2) then
       if (abs(y-z).gt.VOFTOL) then
        print *,"y=z expected in 2D"
        stop
       endif
      endif


      igeom=1

      if ((probtype.eq.21).and.(SDIM.eq.2)) then
       ylen=yblob/5.0

       xstart=-xblob/two
       ystart=yblob/two-ylen/four
       xend=xstart+xblob*three/four
       yend=ystart+ylen
       frontrad=zero

!      xstart=xblob/four
!      ystart=yblob/two-ylen/four
!      xend=xstart+xblob/two
!      yend=ystart+ylen
!      frontrad=zero

       if ((x.le.xstart+frontrad).and.(x.ge.xstart).and. &
           (frontrad.gt.zero)) then
        ystart=ystart+(xstart+frontrad-x)*ylen/frontrad
       endif
       if ((x.ge.xend).and.(y.ge.ystart).and.(y.le.yend)) then
        dist=x-xend
       else if ((x.ge.xend).and.(y.le.ystart)) then
        dist=sqrt((x-xend)**2+(y-ystart)**2)
       else if ((x.ge.xend).and.(y.ge.yend)) then
        dist=sqrt((x-xend)**2+(y-yend)**2)
       else if ((y.ge.yend).and.(x.ge.xstart)) then
        dist=y-yend
       else if ((y.ge.yend).and.(x.le.xstart)) then
        dist=sqrt((x-xstart)**2+(y-yend)**2)
       else if ((y.le.ystart).and.(x.ge.xstart)) then
        dist=ystart-y
       else if (x.le.xstart) then
        dist=sqrt((x-xstart)**2+(y-yend)**2)
       else if ((xend-x.lt.x-xstart).and.(xend-x.lt.y-ystart).and. &
                (xend-x.lt.yend-y)) then
        dist=x-xend
       else if ((x-xstart.lt.y-ystart).and.(x-xstart.lt.yend-y)) then
        dist=xstart-x
       else if (y-ystart.lt.yend-y) then
        dist=ystart-y
       else 
        dist=y-yend
       endif

        ! selectgeom, dist<0 in the solid
      else if ((probtype.eq.32).and.(SDIM.eq.2)) then  
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        if (radblob.gt.zero) then
         dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
        else if (radblob.lt.zero) then
         dist=x-xblob
        else
         print *,"radblob cannot be zero for probtype=32"
         stop
        endif
       else
        print *,"levelrz invalid probtype 32"
        stop
       endif

       ! selectgeom
      else if ((probtype.eq.bubbleInPackedColumn).and.(SDIM.eq.2)) then
       ! there are nSize different size of solid particles in nSize zones
       !print*, " bubbleInPackedColumn"
       dist=1.0e9
       do iSphere=1,nSphere
          dist2=sqrt((x-xSphere(iSphere))**2+(y-ySphere(iSphere))**2) & 
                -rSphere(iSphere)
          dist=min(dist,dist2) 
       enddo
      else if ((probtype.eq.30).and.(SDIM.eq.2)) then
       dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
       dist1=x-xblob
       if ((dist.le.zero).and.(dist1.le.zero)) then
        if (dist1.gt.dist) then
         dist=dist1
        endif
       else if ((dist.ge.zero).and.(dist1.le.zero)) then
        dist=dist
       else if (abs(y-yblob).le.radblob) then
        dist=dist1
       else if (y.ge.yblob) then
        dist=sqrt((x-xblob)**2+(y-yblob-radblob)**2)
       else
        dist=sqrt((x-xblob)**2+(y-yblob+radblob)**2)
       endif
      else if ((probtype.eq.33).and.(SDIM.eq.2)) then
       call getslopeparms(xstart,ystart,zstart, &
        xend,yend,zend)

       slope=(yend-ystart)/(xend-xstart)
       intercept=ystart-slope*xstart
       xint=(x+slope*(y-intercept))/(one+slope*slope)  
       yint=slope*xint+intercept
       dist=sqrt( (x-xint)**2 + (y-yint)**2 )
       if (y.lt.slope*x+intercept) then
        dist=-dist
       endif
      else if ((probtype.eq.34).and.(SDIM.eq.2)) then
! xblob=yblob=zero,radblob=3/4 is used in 2d case
       dist=radblob-abs(x)
!      dist1=0.9-abs(y-1.0)
!      dist=min(dist,dist1)
      else
       print *,"probtype out of range in selectgeom"
      endif

      return
      end subroutine selectgeom


      subroutine microfabgeom(rr,z,dist)
      IMPLICIT NONE
      real(amrex_real) x,y,z,rr,dist


      real(amrex_real) NOD,NID,NPT,CHH,CHW,JLEN,incline,dist1,xdiff,ydiff

      call microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)

      if (SDIM.eq.2) then

      if ((rr.ge.half*NOD).and.(z.ge.JLEN)) then
       dist=z-JLEN
      else if (z.ge.JLEN) then
       dist=sqrt( (z-JLEN)**2 + (half*NOD-rr)**2 )
      else if ((z.ge.JLEN-NPT).and.(z.le.JLEN)) then
       incline=half*NOD+half*(NID-NOD)*(JLEN-z)/NPT
       dist=incline-rr
       if (rr.ge.incline) then
         if (dist.le.z-JLEN) then
          dist=z-JLEN
         endif
         if (dist.le.JLEN-NPT-z) then
          dist=JLEN-NPT-z
         endif
       endif
      else
       xdiff=half*CHH-rr

       if (xdiff.ge.zero) then
        dist=xdiff
        if (rr.ge.half*NID) then
         dist1=JLEN-NPT-z
        else
         dist1=sqrt( (JLEN-NPT-z)**2 + (rr-half*NID)**2 )
        endif 
        if (dist1.lt.dist) then
         dist=dist1
        endif
       else if (xdiff.le.zero) then
        dist=xdiff
       endif
      endif

! walls of domain coincide with nozzle here!!
      if (axis_dir.eq.13) then
       if ((rr.gt.NID*half).and.(z.lt.JLEN-NPT)) then
        dist=JLEN-NPT-z
       endif
      endif 

      else if (SDIM.eq.3) then

! we force the axis of symmetry at the origin for large geometry
      if (axis_dir.eq.13) then
       xblob=zero
       yblob=zero
      endif

      rr=sqrt( (x-xblob)**2 + (y-yblob)**2 )
      if ((rr.ge.half*NOD).and.(z.ge.JLEN)) then
       dist=z-JLEN
      else if (z.ge.JLEN) then
       dist=sqrt( (z-JLEN)**2 + (half*NOD-rr)**2 )
      else if ((z.ge.JLEN-NPT).and.(z.le.JLEN)) then
       incline=half*NOD+half*(NID-NOD)*(JLEN-z)/NPT
       dist=incline-rr
       if (rr.ge.incline) then
         if (dist.le.z-JLEN) then
          dist=z-JLEN
         endif
         if (dist.le.JLEN-NPT-z) then
          dist=JLEN-NPT-z
         endif
       endif
      else
       xdiff=half*CHH-abs(x-xblob)
       ydiff=half*CHW-abs(y-yblob)

       if ((xdiff.ge.zero).and.(ydiff.ge.zero)) then
        if (xdiff.lt.ydiff) then
         dist=xdiff
        else 
         dist=ydiff
        endif
        if (rr.ge.half*NID) then
         dist1=JLEN-NPT-z
        else
         dist1=sqrt( (JLEN-NPT-z)**2 + (rr-half*NID)**2 )
        endif 
        if (dist1.lt.dist) then
         dist=dist1
        endif
       else if ((xdiff.le.zero).and.(ydiff.ge.zero)) then
        dist=xdiff
       else if ((xdiff.ge.zero).and.(ydiff.le.zero)) then
        dist=ydiff
       else
        dist=-sqrt(xdiff**2 + ydiff**2)
       endif
      endif

      else
       print *,"dimension bust"
       stop
      endif


      return
      end subroutine microfabgeom


      subroutine construct(x,y,xx1,yy1,xx2,yy2,dd,lessflag,numline,dist)
      IMPLICIT NONE
      integer numline
      real(amrex_real) x,y
      real(amrex_real) xx1(maxnumline),yy1(maxnumline)
      real(amrex_real) xx2(maxnumline),yy2(maxnumline)
      integer dd(maxnumline)
      integer lessflag(maxnumline)
      real(amrex_real) dist,slope,intercept,localdist,xc,yc
      real(amrex_real) xmin,xmax,ymin,ymax,xtemp,ytemp
      real(amrex_real) xdiff,ydiff,linenorm,xdiff2,ydiff2

      integer iline,hitflag

      dist=1.0D+10
      do iline=1,numline
       hitflag=0

!      print *,"i,x1,y1,x2,y2,dd,ls ",iline,xx1(iline),yy1(iline),
!    &   xx2(iline),yy2(iline),dd(iline),lessflag(iline)
       ydiff=yy2(iline)-yy1(iline)
       xdiff=xx2(iline)-xx1(iline)
       linenorm=sqrt(xdiff**2 + ydiff**2)
       if (linenorm.gt.EPS10) then

       if (dd(iline).eq.0) then
        slope=ydiff/xdiff
        intercept=yy2(iline)-slope*xx2(iline)

        ytemp=slope*x+intercept
        call getminmax(ytemp,yy1(iline),yy2(iline),ymin,ymax)

        xc=(x+slope*(y-intercept))/(slope*slope+one)
        yc=slope*x+intercept
        if ((xc.ge.xx1(iline)).and.(xc.le.xx2(iline))) then
         localdist=sqrt( (x-xc)**2 + (y-yc)**2 )
         hitflag=1
        else if (xc.le.xx1(iline)) then
         xdiff2=x-xx1(iline)
         ydiff2=y-yy1(iline)
         localdist=sqrt( xdiff2**2 + ydiff2**2 )
         hitflag=1
        else 
         xdiff2=x-xx2(iline)
         ydiff2=y-yy2(iline)
         localdist=sqrt( xdiff2**2 + ydiff2**2 )
         hitflag=1
        endif
        if ( ((lessflag(iline).eq.1).and.(y.le.slope*x+intercept)).or. &
             ((lessflag(iline).eq.0).and.(y.ge.slope*x+intercept)) ) then
         localdist=-localdist
        endif
       else 
        slope=xdiff/ydiff
        intercept=xx2(iline)-slope*yy2(iline)

        xtemp=slope*y+intercept
        call getminmax(xtemp,xx1(iline),xx2(iline),xmin,xmax)

        yc=(y+slope*(x-intercept))/(slope*slope+one)
        xc=slope*yc+intercept
        if ((yc.ge.yy1(iline)).and.(yc.le.yy2(iline))) then
         localdist=sqrt( (x-xc)**2 + (y-yc)**2 )
         hitflag=1
        else if (yc.le.yy1(iline)) then
         xdiff2=x-xx1(iline)
         ydiff2=y-yy1(iline)
         localdist=sqrt( xdiff2**2 + ydiff2**2 )
         hitflag=1
        else 
         xdiff2=x-xx2(iline)
         ydiff2=y-yy2(iline)
         localdist=sqrt( xdiff2**2 + ydiff2**2 )
         hitflag=1
        endif
        if ( ((lessflag(iline).eq.1).and.(x.le.slope*y+intercept)).or. &
             ((lessflag(iline).eq.0).and.(x.ge.slope*y+intercept)) ) then
         localdist=-localdist
        endif
       endif
       endif

       if (hitflag.eq.1) then
       if (abs(localdist).lt.abs(dist)) then
        dist=localdist
       endif
       endif
      enddo

      return
      end subroutine construct


      subroutine microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)
      IMPLICIT NONE
      real(amrex_real) NOD,NID,NPT,CHH,CHW,JLEN


      NOD=32.0
      NID=52.0
      NPT=50.0
      CHH=74.0
      CHW=74.0
      JLEN=70.0

      if (axis_dir.eq.13) then
       NOD=30.0
! CHH=84 CHW=360 JLEN=8000
! choose CHH,CHW,JLEN so that ratio of surface area with velbc to volume
! is same as 3d
       CHH=168.0
       CHW=168.0
       JLEN=10913.5
      endif

      return
      end subroutine microfabparm


       ! squeezing channel: inlet width: 30 microns height: 10 microns
       ! outlet width: 30 microns outlet height: 10 microns
       ! xblob2=0.0  xblob3=0.003
       ! yblob2=-0.0015  yblob3=0.0015
       ! zblob2=-0.0005  zblob3=0.0005
       ! dist>0 outside the channel
      subroutine squeeze_channel_dist(x,y,z,dist)
      IMPLICIT NONE
 
      real(amrex_real) x,y,z,zmid,dist,distz


      if ((probtype.eq.5700).and.(axis_dir.eq.5)) then
       if ((xblob2.lt.xblob3).and.(yblob2.lt.yblob3).and. &
           (zblob2.lt.zblob3)) then
       
        if (SDIM.eq.3) then
         zmid=z
        else if (SDIM.eq.2) then
         zmid=half*(zblob2+zblob3)
        else
         print *,"dimension bust"
         stop
        endif
 
        ! first find 2d distance in x-y plane
        if (x.le.xblob2) then
         dist=xblob2-x
        else if (x.le.half*(xblob2+xblob3)) then
         dist=xblob2-x
        else if (x.le.xblob3) then
         if (y.ge.yblob3) then
          dist=x-xblob3
         else if (y.le.yblob2) then
          dist=x-xblob3
         else if (y.ge.zero) then
          dist=-sqrt((x-xblob3)**2+(y-yblob3)**2)
         else
          dist=-sqrt((x-xblob3)**2+(y-yblob2)**2)
         endif
        else if (y.ge.yblob3) then
         if (y-yblob3.le.x-xblob3) then
          dist=y-yblob3
         else
          dist=x-xblob3
         endif
        else if (y.le.yblob2) then
         if (yblob2-y.le.x-xblob3) then
          dist=yblob2-y
         else
          dist=x-xblob3
         endif
        else if (y.ge.zero) then
         dist=y-yblob3
        else
         dist=yblob2-y
        endif

        if (zmid.le.zblob2) then
         distz=zblob2-zmid
        else if (zmid.ge.zblob3) then
         distz=zmid-zblob3
        else if (zmid.le.half*(zblob2+zblob3)) then
         distz=zblob2-zmid
        else if (zmid.ge.half*(zblob2+zblob3)) then
         distz=zmid-zblob3
        else
         print *,"bust for distz"
         stop
        endif
 
        if (distz.ge.zero) then
         if (dist.ge.zero) then
          dist=sqrt( dist**2 +distz**2 )
         else if (dist.le.zero) then
          dist=distz
         else
          print *,"dist invalid squeeze_channel_dist"
          stop
         endif
        else if (distz.le.zero) then
         if (dist.ge.zero) then
          ! do nothing
         else if (dist.le.zero) then
          if (abs(dist).le.abs(distz)) then
           ! do nothing
          else if (abs(dist).ge.abs(distz)) then
           dist=distz
          else
           print *,"dist invalid squeeze_channel_dist"
           stop
          endif
         else
          print *,"dist invalid squeeze_channel_dist"
          stop
         endif
        else
         print *,"distz invalid"
         stop
        endif
 
       else
        print *,"dimensions incorrect"
        stop
       endif
      else
       print *,"probtype or axis_dir incorrect"
       stop
      endif

      return
      end subroutine squeeze_channel_dist


! lessflag=1 if y>mx+b => dist>0
! dd=0 if horizontal line y=mx+b, xx1<xx2
! dd=1 if vertical line x=my+b  , yy1<yy2

      subroutine getminmax(x1,x2,x3,xmin,xmax)
      IMPLICIT NONE
      real(amrex_real) x1,x2,x3,xmin,xmax

      if ((x1.ge.x2).and.(x1.ge.x3)) then
       xmax=x1
      else if ((x2.ge.x1).and.(x2.ge.x3)) then
       xmax=x2
      else
       xmax=x3
      endif
      if ((x1.le.x2).and.(x1.le.x3)) then
       xmin=x1
      else if ((x2.le.x1).and.(x2.le.x3)) then
       xmin=x2
      else
       xmin=x3
      endif
 
      return
      end subroutine getminmax

end module global_distance_module



