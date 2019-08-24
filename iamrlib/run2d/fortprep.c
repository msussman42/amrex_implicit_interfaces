#define _PROB_2D_F_  "%W% %G%"
#undef BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROB_AMR_F.H"
#include "PROB_F.H"
#include "LEVEL_F.H"
#include "COMPARE_F.H"
#include "ArrayLim.H"

#define SDIM 2
#define extendvar 0
#define maxnumline 10
#define MAXSTACK (499999)

#define IMAX_EB (132)
#define NCOMP_EB (4)

#define NOZZLEWIDTH (one/four)

      subroutine sloshing_dist(x,y,dist)
      IMPLICIT NONE
      REAL_T x,y,dist

#include "probdata.H"

      dist=zblob-y

      return
      end

c dist>0 inside of the square.
      subroutine sloshing_geom(x,y,dist,time)
      IMPLICIT NONE
      REAL_T x,y,dist,pos,time

#include "probdata.H"

      call sloshing_tank_pos(time,pos) 
      call squaredist(x,y,xblob-radblob+pos,xblob+radblob+pos,yblob-radblob,
     &  yblob+radblob,dist)
      dist=-dist

      return
      end

      subroutine sloshing_tank_pos(time,pos)
      IMPLICIT NONE
      REAL_T time,pos

#include "probdata.H"

      pos=-cos(vinletgas*time)/four

      return
      end

      subroutine sloshing_tank_vel(time,vel)
      IMPLICIT NONE
      REAL_T time,vel

#include "probdata.H"

      vel=vinletgas*sin(vinletgas*time)/four

      return
      end

      subroutine FORT_ICEINIT(xnew,xrad,ncomp,nballs)
      IMPLICIT NONE

      INTEGER_T nballs,ncomp

      REAL_T xnew(ncomp),xrad(ncomp)
      REAL_T xx,yy,zz,u,v,w
      REAL_T xicenew(3,nballs),xiceold(3,nballs)
      REAL_T radice(nballs)
      INTEGER_T i,ib,itot
      INTEGER_T dir

#include "probdata.H"

      if (3*nballs.ne.ncomp) then
       print *,"ncomp invalid"
       stop
      endif

      itot=1
      do i=1,nballs
       radice(i)=0.5
       do dir=1,3
        xrad(itot)=radice(i)
        itot=itot+1
       enddo
      enddo

      do i=1,nballs
       xx=zero
       yy=zero
       zz=zero

       xicenew(1,i)=xx
       xicenew(2,i)=yy
       xicenew(3,i)=zz
      enddo
      itot=1
      do i=1,nballs
       do dir=1,3
        xnew(itot)=xicenew(dir,i)
        itot=itot+1
       enddo
      enddo

      return
      end

c this only works for direction by direction stretching!

      subroutine icevel(xx,yy,zz,uu,vv,ww,
     & fdata,DIMS(fdata),flo,fhi,dt)
      IMPLICIT NONE

      REAL_T xx,yy,zz,uu,vv,ww,dt
      INTEGER_T flo(SDIM),fhi(SDIM)
      INTEGER_T DIMDEC(fdata)
      REAL_T fdata(DIMV(fdata),2*SDIM+1)

      INTEGER_T i,j,k,dir,ifound
      REAL_T xposlo,xposhi

      call checkbound(flo,fhi,DIMS(fdata),1,-1,1)

      ifound=0 
      j=flo(2)
      do i=flo(1),fhi(1)
       xposlo=half*(fdata(i-1,j,SDIM+2)+fdata(i,j,SDIM+2)) 
       xposhi=half*(fdata(i+1,j,SDIM+2)+fdata(i,j,SDIM+2)) 
       if ((xx.ge.xposlo).and.(xx.le.xposhi)) then
        ifound=1
        goto 10
       endif
      enddo
10    continue
      if (ifound.eq.1) then
       ifound=0
       do j=flo(2),fhi(2)
        xposlo=half*(fdata(i,j-1,SDIM+3)+fdata(i,j,SDIM+3))
        xposhi=half*(fdata(i,j+1,SDIM+3)+fdata(i,j,SDIM+3)) 
        if ((yy.ge.xposlo).and.(yy.le.xposhi)) then
         ifound=1
         goto 20
        endif
       enddo
20     continue
      endif

      uu=zero
      vv=zero
      ww=zero
      if (ifound.eq.1) then
       uu=fdata(i,j,1)*dt 
       vv=fdata(i,j,2)*dt 
      endif 

      return
      end

      subroutine FORT_ADVPART(xold,xnew,xrad,fdata,DIMS(fdata),
     & flo,fhi,dxfine,ncomp,nballs,dt)
      IMPLICIT NONE

      INTEGER_T ncomp,nballs
      INTEGER_T flo(SDIM),fhi(SDIM)
      INTEGER_T DIMDEC(fdata)
       
      REAL_T dxfine(SDIM),dt
      REAL_T xold(ncomp),xnew(ncomp),xrad(ncomp)
      REAL_T xiceold(3,nballs),xicenew(3,nballs)
      REAL_T velocity(3,nballs)
      REAL_T radice(nballs)
      REAL_T fdata(DIMV(fdata),2*SDIM+1)
      INTEGER_T ib,dir,itot

      call checkbound(flo,fhi,DIMS(fdata),1,-1,1)
      itot=1
      do ib=1,nballs
       radice(ib)=xrad(itot)
       do dir=1,3
        xiceold(dir,ib)=xold(itot)
        xicenew(dir,ib)=xnew(itot)
        itot=itot+1
       enddo
       call icevel(
     &  xiceold(1,ib),xiceold(2,ib),xiceold(3,ib),
     &  velocity(1,ib),velocity(2,ib),velocity(3,ib),
     &  fdata,DIMS(fdata),flo,fhi,dt)
      enddo
      
      itot=1
      do ib=1,nballs
       do dir=1,3
        xicenew(dir,ib)=xiceold(dir,ib)+velocity(dir,ib)
        xnew(itot)=xicenew(dir,ib)
        itot=itot+1
       enddo
      enddo

      return
      end

      subroutine FORT_PROBINIT ()

#include "probdata.H"

      call FORT_INPUTBOAT()

      return
      end


       subroutine FORT_OVERRIDE(rz_flag,
     &   ccvorterr,ccdenfact,ccvelfact,ccxblob,ccyblob,cczblob,
     &   ccradblob,
     &   ccxblob2,ccyblob2,cczblob2,ccradblob2,
     &   ccxblob3,ccyblob3,cczblob3,ccradblob3,
     &   ccxblob4,ccyblob4,cczblob4,ccradblob4,
     &   ccxblob5,ccyblob5,cczblob5,ccradblob5,
     &   ccxblob6,ccyblob6,cczblob6,ccradblob6,
     &   ccxblob7,ccyblob7,cczblob7,ccradblob7,
     &   ccxblob8,ccyblob8,cczblob8,ccradblob8,
     &   ccxblob9,ccyblob9,cczblob9,ccradblob9,
     &   ccxblob10,ccyblob10,cczblob10,ccradblob10,
     &   ccxactive,ccyactive,cczactive,
     &   ccractivex,
     &   ccractivey,
     &   ccractivez,
     &   cctwater,cctvapor,cctsatdef,ccpsatdef,
     &   cccpwater,cccpvapor,cccvwater,cccvvapor,
     &   ccgas_constant,
     &   ccgammawater,ccgammavapor,
     &   ccprobtype,ccadv_dir,ccadv_vel,ccaxis_dir,
     &   ccrgasinlet,cccontactangle,
     &   ccvinletgas,cctwall,ccwalltemp,cctcenter,ccsolidradius,
     &   ccdenair,ccdenwater,ccdenvapor,
     &   ccadvbot,ccviscburn,ccviscunburn,ccviscvapor,cchpressure,
     &   ccpcav,ccdencav,ccpchopp,ccdenchopp,ccsoundchopp,
     &   ccEchopp,ccbpres,ccbden,
     &   ccatait,ccbtait,ccrho0tait,ccdenatdep,ccrho0gas,ccajwl,
     &   ccbjwl,cccjwl,ccdjwl,ccejwl,ccr1jwl,ccr2jwl,ccsolvevapor,
     &   ccadiabatic,
     &   ccproblox,ccprobloy,ccprobloz,ccprobhix,ccprobhiy,ccprobhiz)

       IMPLICIT NONE
       INTEGER_T rz_flag
       INTEGER_T ccprobtype,ccadv_dir,ccaxis_dir,ccsolvevapor,ccadiabatic
       REAL_T ccvorterr,ccdenfact,ccvelfact,ccxblob,ccyblob,
     &   cczblob,ccradblob,
     &   ccxblob2,ccyblob2,cczblob2,ccradblob2,
     &   ccxblob3,ccyblob3,cczblob3,ccradblob3,
     &   ccxblob4,ccyblob4,cczblob4,ccradblob4,
     &   ccxblob5,ccyblob5,cczblob5,ccradblob5,
     &   ccxblob6,ccyblob6,cczblob6,ccradblob6,
     &   ccxblob7,ccyblob7,cczblob7,ccradblob7,
     &   ccxblob8,ccyblob8,cczblob8,ccradblob8,
     &   ccxblob9,ccyblob9,cczblob9,ccradblob9,
     &   ccxblob10,ccyblob10,cczblob10,ccradblob10,
     &   ccxactive,ccyactive,cczactive,
     &   ccractivex,
     &   ccractivey,
     &   ccractivez,
     &   cctwater,cctvapor,cctsatdef,ccpsatdef,
     &   cccpwater,cccpvapor,cccvwater,cccvvapor,
     &   ccgas_constant,
     &   ccgammawater,ccgammavapor,
     &   ccadv_vel,ccrgasinlet,cccontactangle,
     &   ccvinletgas,cctwall,ccwalltemp,
     &   cctcenter,ccsolidradius,
     &   ccdenair,ccdenwater,ccdenvapor,ccadvbot,
     &   ccviscburn,
     &   ccviscunburn,ccviscvapor,cchpressure,
     &   ccpcav,ccdencav,ccpchopp,ccdenchopp,ccsoundchopp,
     &   ccEchopp,ccbpres,ccbden,
     &   ccatait,ccbtait,ccrho0tait,ccdenatdep,ccrho0gas,ccajwl,
     &   ccbjwl,cccjwl,ccdjwl,ccejwl,ccr1jwl,ccr2jwl,
     &   ccproblox,ccprobloy,ccprobloz,ccprobhix,ccprobhiy,ccprobhiz


#include "probdata.H"

       character press_file*20
       INTEGER_T   error  

       problox=ccproblox
       probloy=ccprobloy
       probloz=ccprobloz
       probhix=ccprobhix
       probhiy=ccprobhiy
       probhiz=ccprobhiz

       levelrz=rz_flag
       vorterr=ccvorterr
       denfact=ccdenfact
       velfact=ccvelfact
       xblob=ccxblob
       yblob=ccyblob
       zblob=cczblob
       radblob=ccradblob

       xblob2=ccxblob2
       yblob2=ccyblob2
       zblob2=cczblob2
       radblob2=ccradblob2

       xblob3=ccxblob3
       yblob3=ccyblob3
       zblob3=cczblob3
       radblob3=ccradblob3

       xblob4=ccxblob4
       yblob4=ccyblob4
       zblob4=cczblob4
       radblob4=ccradblob4

       xblob5=ccxblob5
       yblob5=ccyblob5
       zblob5=cczblob5
       radblob5=ccradblob5

       xblob6=ccxblob6
       yblob6=ccyblob6
       zblob6=cczblob6
       radblob6=ccradblob6

       xblob7=ccxblob7
       yblob7=ccyblob7
       zblob7=cczblob7
       radblob7=ccradblob7

       xblob8=ccxblob8
       yblob8=ccyblob8
       zblob8=cczblob8
       radblob8=ccradblob8

       xblob9=ccxblob9
       yblob9=ccyblob9
       zblob9=cczblob9
       radblob9=ccradblob9

       xblob10=ccxblob10
       yblob10=ccyblob10
       zblob10=cczblob10
       radblob10=ccradblob10

       xactive=ccxactive
       yactive=ccyactive
       zactive=cczactive
       ractivex=ccractivex
       ractivey=ccractivey
       ractivez=ccractivez

       twater=cctwater
       tvapor=cctvapor
       tsatdef=cctsatdef
       psatdef=ccpsatdef
       cpwater=cccpwater
       cpvapor=cccpvapor
       cvwater=cccvwater
       cvvapor=cccvvapor
       gas_constant=ccgas_constant
       gammawater=ccgammawater
       gammavapor=ccgammavapor

       probtype=ccprobtype
       adv_dir=ccadv_dir
       adv_vel=ccadv_vel
       axis_dir=ccaxis_dir
       rgasinlet=ccrgasinlet
       contactangle=cccontactangle
       vinletgas=ccvinletgas
       twall=cctwall
       walltemp=ccwalltemp
       tcenter=cctcenter
       solidradius=ccsolidradius
       denair=ccdenair
       denwater=ccdenwater
       denvapor=ccdenvapor
       advbot=ccadvbot
       viscburn=ccviscburn
       viscunburn=ccviscunburn
       viscvapor=ccviscvapor
       hpressure=cchpressure

       pcav=ccpcav 
       dencav=ccdencav 
       pchopp=ccpchopp 
       denchopp=ccdenchopp
       soundchopp=ccsoundchopp
       Echopp=ccEchopp

       bpresglob=ccbpres
       bdenglob=ccbden

       atait=ccatait
       btait=ccbtait
       rho0tait=ccrho0tait
       denatdep=ccdenatdep
       rho0gas=ccrho0gas
       ajwl=ccajwl
       bjwl=ccbjwl
       cjwl=cccjwl
       djwl=ccdjwl
       ejwl=ccejwl
       r1jwl=ccr1jwl
       r2jwl=ccr2jwl

       solvevapor=ccsolvevapor
       adiabatic=ccadiabatic

       print *,"fort: problox,y,z,hix,y,z ",problox,probloy,probloz,
     &   probhix,probhiy,probhiz
       print *,"fort: rz ",levelrz
       print *,"fort: vorterr,denfact,velfact,xblob,yblob,zblob ",
     &   vorterr,denfact,velfact,xblob,yblob,zblob
       print *,"fort: radblob,probtype,adv_dir,adv_vel,axis_dir ",
     &   radblob,probtype,adv_dir,adv_vel,axis_dir
       print *,"fort: xblob2,yblob2,zblob2,radblob2",
     &   xblob2,yblob2,zblob2,radblob2
       print *,"fort: xblob3,yblob3,zblob3,radblob3",
     &   xblob3,yblob3,zblob3,radblob3
       print *,"fort: xblob4,yblob4,zblob4,radblob4",
     &   xblob4,yblob4,zblob4,radblob4
       print *,"fort: xblob5,yblob5,zblob5,radblob5",
     &   xblob5,yblob5,zblob5,radblob5
       print *,"fort: xblob6,yblob6,zblob6,radblob6",
     &   xblob6,yblob6,zblob6,radblob6
       print *,"fort: xblob7,yblob7,zblob7,radblob7",
     &   xblob7,yblob7,zblob7,radblob7
       print *,"fort: xblob8,yblob8,zblob8,radblob8",
     &   xblob8,yblob8,zblob8,radblob8
       print *,"fort: xblob9,yblob9,zblob9,radblob9",
     &   xblob9,yblob9,zblob9,radblob9
       print *,"fort: xblob10,yblob10,zblob10,radblob10",
     &   xblob10,yblob10,zblob10,radblob10

       print *,"fort: xactive,yactive,zactive,ractivexyz",
     &   xactive,yactive,zactive,ractivex,ractivey,ractivez
       print *,"fort: tw,tv,ts,ps ",twater,tvapor,tsatdef,psatdef
       print *,"fort: cpw,cpv,cvw,cvv,gammaw,gammav ",
     &   cpwater,cpvapor,cvwater,cvvapor,gammawater,gammavapor
       print *,"fort: gas_constant ",gas_constant
       print *,"fort:rgasinlet,vinletgas,twall,walltemp,tcenter,denair",
     &   rgasinlet,vinletgas,twall,walltemp,tcenter,denair
       print *,"fort:solidradius ",solidradius
       print *,"fort:contactangle ",contactangle
       print *,"fort: denwater,advbot,viscburn,viscunburn ",
     &   denwater,advbot,viscburn,viscunburn
       print *,"fort: denvapor,viscvapor,hpressure ",
     &   denvapor,viscvapor,hpressure
       print *,"fort: pcav,dencav,pchopp,denchopp,soundchopp",
     &   pcav,dencav,pchopp,denchopp,soundchopp
       print *,"fort: Echopp",Echopp
       print *,"fort:bpres,bden ",
     &   bpresglob,bdenglob
       print *,"fort:atait,btait,rho0tait,rho0gas,denatdep",atait,btait,
     &   rho0tait,rho0gas,denatdep
       print *,"fort:ajwl,bjwl,cjwl,djwl,ejwl,r1jwl,r2jwl ",ajwl,
     &   bjwl,cjwl,djwl,ejwl,r1jwl,r2jwl
       print *,"fort:solvevapor ",solvevapor
       print *,"fort:adiabatic ",adiabatic
     
       if ((probtype.eq.22).and.
     &    ((axis_dir.eq.11).or.(axis_dir.eq.12)) )  then
        print *,"this stuff no longer used"
        stop
       endif

       return
       end


      function CLS(phi,eps)
      REAL_T CLS,phi,eps,temp

      temp=phi/(two*eps)
      CLS=half*(sinh(temp)/cosh(temp)+one)

      return
      end 

      function myfact(n)
      INTEGER_T myfact,i,n

      if (n.eq.0) then
       myfact=1
      else
       myfact=1
       do i=2,n
        myfact=myfact*i 
       enddo
      endif

      return
      end

c lessflag=1 if y>mx+b => dist>0
c dd=0 if horizontal line y=mx+b, xx1<xx2
c dd=1 if vertical line x=my+b  , yy1<yy2

      subroutine getminmax(x1,x2,x3,xmin,xmax)
      IMPLICIT NONE
      REAL_T x1,x2,x3,xmin,xmax

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
      end

      subroutine inletpipedist(x,y,dist)
      IMPLICIT NONE

      REAL_T x,y,dist
#include "probdata.H"

      dist=x-xblob-radblob*sin(two*Pi*y/yblob)

      return
      end

      subroutine inletpipevel(vel,x,y,phase)
      IMPLICIT NONE

      REAL_T vel,x,y,dist
      REAL_T a1,a2,b1,b2,AA,LL,YY,coef1,coef2
      INTEGER_T phase

#include "probdata.H"

      call inletpipedist(x,y,dist)
      if (velfact.le.zero) then
       if (phase.eq.0) then
        vel=vinletgas
       else if (phase.eq.1) then
        vel=advbot
       else
        print *,"phase invalid"
        stop
       endif
      else
       if ((xblob5.ne.zero).or.(yblob5.ne.zero)) then
        print *,"inletpipevel must be modified for xblob5,yblob5<>0"
        print *,"or set velfact=0.0"
        stop
       endif
       if ((viscburn.le.zero).or.(viscunburn.le.zero)) then
        b1=zero
        b2=zero
        a1=zero
        a2=zero
       else
        AA=xblob-problox
        LL=probhix-problox
        a2=half*velfact/viscunburn
        b2=half*velfact/viscburn
        coef1=(advbot+b2*AA*AA-vinletgas)/(LL-AA)-
     &    a2*(LL-AA)
        coef2=-two*b2*AA-viscunburn*two*a2*(LL-AA)/viscburn
        a1=(coef1+coef2*AA/(LL-AA))/
     &     (one+AA*viscunburn/((LL-AA)*viscburn)) 
        b1=coef2-a1*viscunburn/viscburn
       endif
       if (phase.eq.0) then
        YY=dist+AA
        if (YY.gt.LL) then
         vel=vinletgas
        else if (YY.lt.AA) then
         vel=vinletgas+a1*(LL-AA)+a2*((LL-AA)**2)
        else
         vel=vinletgas+a1*(LL-YY)+a2*((LL-YY)**2)
        endif
       else if (phase.eq.1) then
        YY=AA+dist
        if (YY.gt.AA) then
         vel=advbot+b1*AA+b2*AA*AA
        else if (YY.lt.zero) then
         vel=advbot
        else
         vel=advbot+b1*YY+b2*YY*YY
        endif
       else
        print *,"phase invalid"
        stop
       endif
      endif

      return
      end


      subroutine rampvel(time,vel)
      IMPLICIT NONE

      REAL_T time,vel,tcutoff

#include "probdata.H"

      if ((probtype.eq.63).or.(probtype.eq.64)) then
       tcutoff=1.0
       if (time.gt.tcutoff) then
        vel=adv_vel
       else
        vel=time*adv_vel/tcutoff
       endif
      else
       vel=adv_vel
      endif

      return
      end

      subroutine construct(x,y,xx1,yy1,xx2,yy2,dd,lessflag,numline,dist)
      IMPLICIT NONE
      INTEGER_T numline
      REAL_T x,y
      REAL_T xx1(maxnumline),yy1(maxnumline)
      REAL_T xx2(maxnumline),yy2(maxnumline)
      INTEGER_T dd(maxnumline)
      INTEGER_T lessflag(maxnumline)
      REAL_T dist,slope,intercept,localdist,xc,yc
      REAL_T xmin,xmax,ymin,ymax,xtemp,ytemp
      REAL_T xdiff,ydiff,linenorm,xdiff2,ydiff2,costheta

      INTEGER_T iline,hitflag

      dist=1.0E+10
      do iline=1,numline
       hitflag=0

c      print *,"i,x1,y1,x2,y2,dd,ls ",iline,xx1(iline),yy1(iline),
c    &   xx2(iline),yy2(iline),dd(iline),lessflag(iline)
       ydiff=yy2(iline)-yy1(iline)
       xdiff=xx2(iline)-xx1(iline)
       linenorm=sqrt(xdiff**2 + ydiff**2)
       if (linenorm.gt.1.0E-10) then

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
        if ( ((lessflag(iline).eq.1).and.(y.le.slope*x+intercept)).or.
     &       ((lessflag(iline).eq.0).and.(y.ge.slope*x+intercept)) ) then
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
        if ( ((lessflag(iline).eq.1).and.(x.le.slope*y+intercept)).or.
     &       ((lessflag(iline).eq.0).and.(x.ge.slope*y+intercept)) ) then
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
      end

c time is in microseconds, LL is in microns,pressure is first in atmospheres
c and then converted since I scale
c the NS equations by
c L=10^-4 cm and T=10^-6 s and U=L/T=10^2 cm/s
c 1atm=1.013x10^6 dynes/cm^2=1.013x10^6 g/(s^2 cm)
c P is scaled by 1.013x10^2 
c LL is the length that the pressure will decay to zero.

c 0 < x < 2* xblob  0 < y < 2*yblob
c for now, xblob=yblob=18 (which used to be center of nozzle

c define Y_SHIFT 92.
#define Y_SHIFT 20

      subroutine distwall_r(x,y,z,dist,flag)
      IMPLICIT NONE
      REAL_T x,y,z,dist
      INTEGER_T flag
      REAL_T xr, yr, zr
      REAL_T nx, nz, dn
      nx = 82.
      nz = 11.
      dn = sqrt(nx*nx+nz*nz)
      nx = nx/dn
      nz = nz/dn
      xr =  nx*x - nz*z
      zr =  nz*x + nx*z 
      yr = y
      call distwall(xr,yr,zr,dist,flag)
      return
      end

      subroutine igorgeom(y,z,dist)
      IMPLICIT NONE
      REAL_T x,y,z,dist,JLEN

#include "probdata.H"

c change this for axis_dir.eq.8!!!
      x=36.0
      if (axis_dir.eq.8) then
       call distwall(x,y,z,dist,0)
      else if ((axis_dir.eq.9).or.(axis_dir.eq.10)) then
       JLEN=170.0
       call distwall(-z+JLEN,x-36.0,y-Y_SHIFT,dist,0)
      else
       print *,"invalid axis_dir in igorgeom"
       stop
      endif
 
      return
      end

      subroutine igordist(y,z,dist)
      IMPLICIT NONE
      REAL_T x,y,z,dist,JLEN

#include "probdata.H"

c change this for axis_dir.eq.8!!!
      x=36.0
      if (axis_dir.eq.8) then
       call distwall(x,y,z,dist,1)
      else if ((axis_dir.eq.9).or.(axis_dir.eq.10)) then
       JLEN=170.0
       call distwall(-z+JLEN,x-36.0,y-Y_SHIFT,dist,1)
      else
       print *,"invalid axis_dir in igordist"
       stop
      endif

      return
      end
 

      subroutine igordist_orig(y,z,dist)
      IMPLICIT NONE
      REAL_T x,y,z,dist,CHX,CHY,JLEN,distsquare
      INTEGER_T insquare

#include "probdata.H"

c change this for axis_dir.eq.8!!!
      x=36.0

       if (axis_dir.eq.8) then
        CHX=22.0
        CHY=13.0
        JLEN=70.0
       else if ((axis_dir.eq.9).or.(axis_dir.eq.10)) then
        CHX=21.5
        CHY=19.0
        JLEN=170.0
       else
        print *,"axis_dir out of range,igordist"
        stop
       endif

       if ( (abs(x-xblob).le.CHX/two).and.(abs(y-yblob).le.CHY/two) ) then
        insquare=1
       else
        insquare=0
       endif
       dist=JLEN-z
       if ( insquare.eq.0 ) then
        if (abs(x-xblob).le.CHX/two) then
         distsquare=abs(y-yblob)-CHY/two
        else if (abs(y-yblob).le.CHY/two) then
         distsquare=abs(x-xblob)-CHX/two
        else
         distsquare=sqrt( (CHX/two-abs(x-xblob))**2 +
     &                    (CHY/two-abs(y-yblob))**2 )
        endif  
        dist=dist-distsquare
       endif

      return
      end

      subroutine igorpressure(y,z,time,pressure,LL)
      IMPLICIT NONE
      REAL_T x,y,z,time,pressure,LL,JLEN,ALPHA,BETA
      REAL_T P0,PB,TC
      REAL_T zlow,zhigh,xradius,tscale

#include "probdata.H"

      x=36.0
      P0=73.5
      tscale=0.021
      PB=one
      ALPHA=one+time/0.021
      BETA=0.755
      TC=45.0
      pressure=P0*exp(-BETA*log(ALPHA))-PB
      P0=73.5
      ALPHA=0.36
      TC=13.3
      if (time.lt.ALPHA) then
       pressure=P0*(one-time/ALPHA)-PB
      else
       pressure=-PB
      endif

      if (time.lt.TC) then
c       pressure=pressure*exp(-time*time*time*time/(15.0*15.0*15.0*15.0))
      else
        pressure=zero
      endif
      pressure=pressure*1.013E+2

c     pressure=pressure*0.4

      if ((axis_dir.eq.8).or.(axis_dir.eq.9)) then
       
       if (z.gt.LL) then
        pressure=zero
       else 
        pressure=pressure*(LL-z)*(LL-z)/(LL*LL)
       endif
    
      else if (axis_dir.eq.10) then

       JLEN=170.0
       zlow=JLEN-105.0
       zhigh=JLEN-55.0
       xradius=14.0
       
       if ((z.gt.(JLEN+zhigh)/two).or.
     &     (y.gt.LL)) then
        pressure=zero
       else
        pressure=pressure*(LL-y)*(LL-y)/(LL*LL)
       endif
      else
       print *,"axis_dir out of range in pressure routine"
      endif 

      return 
      end
 
      subroutine igorforce(fy,fz,delta,y,z,time)
      IMPLICIT NONE
      REAL_T fx,fy,fz,delta,x,y,z,time,p2,p1,LL,deltasmall

#include "probdata.H"

      x=36.0
      if ((axis_dir.eq.8).or.(axis_dir.eq.9)) then
       LL=50.0
      else if (axis_dir.eq.10) then
c       LL=8.0
c       LL=12.
c       LL=6.0
       LL = Y_SHIFT-4.
      else
       print *,"axis_dir out of range in igorforce"
      endif
      
      deltasmall=delta/two
 
      if (axis_dir.eq.10) then
       call igorpressure(y+deltasmall,z,time,p2,LL)
       call igorpressure(y-deltasmall,z,time,p1,LL)
       fy=fy- (p2-p1)/(two*deltasmall)
      endif

      if ((axis_dir.eq.8).or.(axis_dir.eq.9)) then
       call igorpressure(y,z+deltasmall,time,p2,LL)
       call igorpressure(y,z-deltasmall,time,p1,LL)
       fz=fz- (p2-p1)/(two*deltasmall)
      endif

      return
      end
 
      subroutine microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)
      IMPLICIT NONE
      REAL_T NOD,NID,NPT,CHH,CHW,JLEN

#include "probdata.H" 

      NOD=32.0
      NID=52.0
      NPT=50.0
      CHH=74.0
      CHW=74.0
      JLEN=70.0

      if (axis_dir.eq.13) then
       NOD=30.0
c CHH=84 CHW=360 JLEN=8000
c choose CHH,CHW,JLEN so that ratio of surface area with velbc to volume
c is same as 3d
       CHH=168.0
       CHW=168.0
       JLEN=10913.5
      endif

      return
      end


c time is in microseconds, LL is in microns,PTERM is in atmospheres
c 1atm=1.013x10^6 dyne/cm^2

      subroutine microfabpressure(LL,PTERM,time)
      IMPLICIT NONE
      REAL_T x,y,z,rr,dist,LL,PTERM,time
      REAL_T realtime,realpress
      INTEGER_T error

#include "probdata.H" 

      REAL_T NOD,NID,NPT,CHH,CHW,JLEN,incline,dist1,xdiff,ydiff

      call microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)

      LL=half*(JLEN-NPT)
c simple pressure 2
c assume flat meniscus, so quiescent pressure is 1atm
      if (axis_dir.eq.13) then
       PTERM=zero
      else
       if (time.le.6.0) then
        PTERM=0.7-1.0
       else if (time.le.14.0) then
        PTERM=1.8-1.0
       else if (time.le.21.0) then
        PTERM=0.2-1.0
       else if (time.le.28.0) then
        PTERM=1.2-1.0
       else 
        PTERM=zero
       endif

       realtime=time*1.0E-6
       call pressure_bc(realpress,realtime,error)
       PTERM=realpress/1.00E+06 - one 
      endif
            
      return
      end

      subroutine microfabdist(rr,z,dist)
      IMPLICIT NONE
      REAL_T x,y,z,rr,dist

#include "probdata.H" 

      REAL_T NOD,NID,NPT,CHH,CHW,JLEN,incline,dist1,xdiff,ydiff

      call microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)

      incline=JLEN 
      if (rr.le.half*NOD) then
       if (axis_dir.eq.12) then
        incline=JLEN+sqrt(16.9**2-(half*NOD)**2)-sqrt(16.9**2-rr**2)
       endif
       dist=incline-z
      else
       if (rr.le.half*(NPT+NOD)) then
        dist=incline-z-(rr-half*NOD)
       else
        dist=incline-z-half*NPT
       endif
      endif

      return
      end
 
      subroutine microfabgeom(rr,z,dist)
      IMPLICIT NONE
      REAL_T x,y,z,rr,dist

#include "probdata.H" 

      REAL_T NOD,NID,NPT,CHH,CHW,JLEN,incline,dist1,xdiff,ydiff

      call microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)

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

c walls of domain coincide with nozzle here!!
      if (axis_dir.eq.13) then
       if ((rr.gt.NID*half).and.(z.lt.JLEN-NPT)) then
        dist=JLEN-NPT-z
       endif
      endif 

      return
      end

      subroutine pressure_bc( press, t    , error )
      IMPLICIT NONE

c     * Return the pressure boundary condition at time t provided by MicroFab
c     *  for the Okidata problems.

c     * CAUTION! Pressure is in dynes/ cm^2 and time is in seconds!
      
c     * This routine assumes that time is measured in seconds, the time values
c     * are equally spaced with spacing dt and that the pressure values are 
c     * known up to time t = 7.00E-05. For the MicroFab 
c     * test problems j= 0 ... = 70.  However,  to allow for future use of this
c     * type of BC, the size of the arrays time and pressbc may have to be 
c     * increased or passed as an argument to this subroutine.  Right now the
c     * number of data points in press_file is hard coded to be 70.

c     * Variables passed in ...

      INTEGER_T   error  
      
      real*8    press       , t
      
c     * The array pressbc(i,j).  The first array contains the time t, the
c     * second contains the value of the pressure on the inflow boundary at
c     * time t.

      real*8    dt          , sigma    

      real*8  time(0:100) ,  pressbc(0:100,1:3)
      INTEGER_T selectpress
      common / pressure_bcs / time      , pressbc   , dt, selectpress

c     INTEGER_T  1234567890, 1234567890, 1234567890, 1234567890, 1234567890

      INTEGER_T  j         
      
      error = 0
      j = int(t / dt)
      
      if (t .gt. 7.00E-05) then

c       * From 70 microseconds on, the pressure BC is 1 atm
        
        press = 1.00E+06
        
      else if ((time(j) .le. t) .and. (t .le. time(j+1))) then
        
c       * Lineraly interpolate between the given time values.
        
        sigma = (t - time(j)) /dt
        press = (1 - sigma) * pressbc(j,selectpress) + 
     &          sigma * pressbc(j+1,selectpress) 
        
      else
       error = -1
       print *,"error in pressure_bc, t= ",t
       stop
      end if

      return
      end

c
c     * This is a subroutine for setting inflow B.C.(Poiseuille flow)
c     * Nozzle radius (cm)   
c
      subroutine inflow_bc(system,x,phi,delx)
      IMPLICIT NONE
c
      INTEGER_T   system
      REAL_T    aveQ,  aveV,  radius,  x,  delx,  phi
      REAL_T Weber,Reynolds,Froude

#include "probdata.H" 

c       radius=0.085 
        radius=0.127

      go to (1,2,3,4,5,6,7,8,9,10,11), system
c
    1   aveQ = 1.3D-2
        go to 20
    2   aveQ = 4.8D-2
        go to 20
    3   aveQ = 2.0D-1
        go to 20
    4   aveQ = 5.0D-1
        go to 20
c   5   aveQ = one
    5   aveQ = 1.1
        go to 20
    6   aveQ = two
        go to 20
    7   aveQ = five
        go to 20
    8   aveQ = 6.8
        go to 20
    9   aveQ = 7.2
        go to 20
   10   aveQ = 15.0
        go to 20
   11   aveQ = 20.0
        go to 20

   20 continue
c
         aveV=aveQ/(Pi*radius**2)

c        Weber=(aveV**2)*radius*1.2238/66.0
c        Reynolds=1.2238*radius*aveV/1.26
         Weber=(aveV**2)*radius*0.996/51.1
         Reynolds=0.996*radius*aveV/0.00958

         Froude=(aveV**2)/(radius*980.0)
         if (1.eq.0) then
          print *,"Weber,Reynolds,Froude ",Weber,Reynolds,Froude
          print *,"1/Weber,1/Reynolds,1/Froude ",one/Weber,
     &       one/Reynolds,one/Froude
          stop
         endif
         if (radblob.ne.one) then
          print *,"dimensionless radius=1"
          stop
         endif
         if ((viscunburn.ne.one).or.(denwater.ne.one)) then
          print *,"dimensionless liquid parameters invalid"
          stop
         endif
         phi=0.0
         aveV=one
         radius=one
         if (abs(x).le.radius) then
          phi=2.0*aveV*(1.0-(abs(x)/radius)**2-((delx/radius)**2)/4.0)       
         endif
c
      return
      end

      subroutine vbc( velocity, t , yval,zval, error )
      IMPLICIT NONE
c velocity is in m/s (or microns/microseconds), t is in microseconds
c zval in microns

      INTEGER_T   error  
      
      real*8    velocity, t, yval,zval
      
      real*8    dt          , sigma    

      real*8  timehist(0:100),zpos(1:50),velbc(0:100,1:50),period,rigidwall
      INTEGER_T itime,ipos

      common / vel_bcs / timehist,zpos,velbc,period,rigidwall,itime,ipos

      INTEGER_T  i,j,istar,jstar 
      real*8 zdiff,tdiff,tlocal
      
      REAL_T NOD,NID,NPT,CHH,CHW,JLEN

      tlocal=t
      if (tlocal.lt.timehist(0)) then
       tlocal=timehist(0)
      endif

      error = 0
      call microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)
c     print *,"in vbc yval,jlen-npt-193,tlocal,period,t(it) ",yval,
c    &  JLEN-NPT-193.0,tlocal,period,timehist(itime)
      if ((yval.le.half).or.(yval.ge.JLEN-NPT-193.0)) then
       velocity=zero
      else

      if (tlocal.gt.period) then
       velocity=zero
      else
       if (tlocal.ge.timehist(itime)) then
        velocity=zero
       else
        do j=0,itime-1
         if ((timehist(j).le.tlocal).and.(timehist(j+1).ge.tlocal)) then
          jstar=j
         endif
        enddo
        if (BL_SPACEDIM.eq.2) then
         zval=half*zpos(ipos)
        endif
        if (zval.ge.zpos(ipos)) then
         velocity=zero
        else if (zval.le.zpos(1)) then
         velocity=zero
        else
         do i=1,ipos-1
          if ((zpos(i).le.zval).and.(zpos(i+1).ge.zval)) then
           istar=i
          endif
         enddo
         tdiff=(timehist(jstar+1)-tlocal)/(timehist(jstar+1)-timehist(jstar))
         zdiff=(zpos(istar+1)-zval)/(zpos(istar+1)-zpos(istar))
         velocity=zdiff*tdiff*velbc(jstar,istar)+
     &     zdiff*(one-tdiff)*velbc(jstar+1,istar)+
     &     (one-zdiff)*tdiff*velbc(jstar,istar+1)+
     &     (one-zdiff)*(one-tdiff)*velbc(jstar+1,istar+1)
        endif
       endif
      endif
c tlocal in range
      endif
c yval in range          
      return
      end

      subroutine readpress( press_file, error )
      IMPLICIT NONE

c     * Read in the pressure boundary conditions provided by MicroFab for the 
c     * Okidata problems and store in the array pressbc. 

c     * This routine assumes that time is measured in seconds, the time values
c     * are equally spaced with spacing dt and that the pressure values are 
c     * known up to time t = 7.00E-05. For the MicroFab 
c     * test problems j= 0 ... = 70.  However,  to allow for future use of this
c     * type of BC, the size of the arrays time and pressbc may have to be 
c     * increased or passed as an argument to this subroutine.  Right now the
c     * number of data points in press_file is hard coded to be 71.

c     * Variables passed in ...

      character press_file*20

      INTEGER_T   error  
      
      real*8    dt

c element 1-complicated 2-simple1 3-simple2

      real*8  time(0:100) ,  pressbc(0:100,1:3)

      INTEGER_T selectpress
      common / pressure_bcs / time      , pressbc   , dt, selectpress

c     INTEGER_T  1234567890, 1234567890, 1234567890, 1234567890, 1234567890

      INTEGER_T  j         
      
c     * Open the data file containing the pressure BCs.  The name of the file
c     * is passed as an argument to the subroutine.

c      open(7, file=press_file, form="formatted", status="old", err=900)

      selectpress=1
      print *,"selectpress = ",selectpress

      do j = 0, 70
         
c        read(7,100) time(j), pressbc(j,1), pressbc(j,2), pressbc(j,3)
        
      end do
 
c      close(7)
 
      do j=0,70
c      print *,"time,pressure ",time(j),pressbc(j,selectpress)
      enddo

      dt = time(1) - time(0)
      error = 0
      
      return
      
c 100   format(' ', e10.6, e9.6, e9.6, e9.6 )
c 900   write(6,910) press_file
c 910   format(" Can't open file = ",a," stopping ...")

      error = -1
      
      return
      end

      subroutine readvel( vel_file, error )
      IMPLICIT NONE

      character vel_file*20

      INTEGER_T   error  
c time is in microseconds, velocity in meter/s, zpos in microns 
    
      INTEGER_T itime,ipos 
      real*8  timehist(0:100),zpos(1:50),velbc(0:100,1:50),period,rigidwall

      common / vel_bcs / timehist,zpos,velbc,period,rigidwall,itime,ipos

      INTEGER_T  i,j         

      period=28.0
      rigidwall=180.0
      itime=11
      ipos=11
      
c      open(7, file=vel_file, form="formatted", status="old", err=901)

c      read(7,101) timehist(0),zpos(1),zpos(2),zpos(3),zpos(4),zpos(5),
c     &  zpos(6),zpos(7),zpos(8),zpos(9),zpos(10)
      zpos(11)=rigidwall

      do j = 0, itime
c       read(7,102) timehist(j),velbc(j,1),velbc(j,2),velbc(j,3),velbc(j,4),
c     &   velbc(j,5),velbc(j,6),velbc(j,7),velbc(j,8),velbc(j,9),velbc(j,10)
       velbc(j,11)=zero
      end do
c  101 format(f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2)  
c  102 format(f7.2,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3)  
      close(7)
 
      print *,"period,rigidwall ",period,rigidwall
      do j=1,ipos
       print *,"j,zpos ",j,zpos(j)
      enddo
      do j=0,itime
       print *,"j,timehist ",j,timehist(j)
      enddo
      do i=0,itime
      do j=1,ipos
       print *,"i,j,velbc ",i,j,velbc(i,j)
      enddo
      enddo
      
      error = 0
      
      return
      
c 901   write(6,911) vel_file
c 911   format(" Can't open file = ",a," stopping ...")

      error = -1
      
      return
      end


      subroutine localparam(Dbdry,grainbdry)
      IMPLICIT NONE
      REAL_T Dbdry,grainbdry

      Dbdry=zero
      grainbdry=zero

      return
      end

      subroutine initsplashparms(newradblob,ylevel,opening)
      IMPLICIT NONE
      REAL_T aval,bval,realvolume,cval,cvolume,newradblob,separation,
     &       opening,ylevel
      INTEGER_T i,j

#include "probdata.H"

      aval=radblob
      bval=three*radblob/two
      realvolume=four*Pi*radblob*radblob*radblob/three
      do i=1,100
       cval=(aval+bval)/two
       call partvolume(cval,contactangle,cvolume)
       if (cvolume.gt.realvolume) then
        bval=cval
       else
        aval=cval
       endif
      enddo
      newradblob=aval
      separation=newradblob*sin(contactangle-Pi/two)
      opening=newradblob*cos(contactangle-Pi/two)
      ylevel=yblob-separation

      return
      end

      subroutine initjetparms(HSB,NOD,NPT,NID,CHH,scaleCHH)
      IMPLICIT NONE
      REAL_T NPT,HSB,NID,NOD,CHH,scaleCHH

#include "probdata.H"

c get rid of floating exceptions !
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
      end

      subroutine rtdist(x,y,dist)
      IMPLICIT NONE
      REAL_T x,y,dist,rholevel,alpha,alpha2,beta,beta2
      REAL_T guess,yguess,xtest,ytest,rho1
      INTEGER_T i1,iter
      REAL_T fx,fxprime,tolerance

#include "probdata.H"

      dist=radblob*cos(xblob*Pi*x)-y
      rholevel=dist


      return
      end

      subroutine initstack()
      IMPLICIT NONE

      INTEGER_T istack(MAXSTACK)
      INTEGER_T istackptr

      common /stackstuff/ istack,istackptr

      istackptr=1

      return
      end

      subroutine stackempty(istatus)
      IMPLICIT NONE

      INTEGER_T istatus

      INTEGER_T istack(MAXSTACK)
      INTEGER_T istackptr

      common /stackstuff/ istack,istackptr

      if (istackptr.lt.1) then
       print *,"istackptr incorrect"
       stop
      endif

      if (istackptr.eq.1) then
       istatus=1
      else
       istatus=0
      endif

      return
      end

      subroutine stackfull(istatus)
      IMPLICIT NONE

      INTEGER_T istatus
      REAL_T upbound
      INTEGER_T istack(MAXSTACK)
      INTEGER_T istackptr

      common /stackstuff/ istack,istackptr

      upbound=MAXSTACK
      if (istackptr+100.ge.upbound) then
       istatus=1
       print *,"stack is almost full!!!"
      else
       istatus=0
      endif

      return
      end

      subroutine pushdata(ix,iy,iz,imult)
      IMPLICIT NONE

      INTEGER_T istack(MAXSTACK),imult,ix,iy,iz
      REAL_T upbound
      INTEGER_T istackptr

      common /stackstuff/ istack,istackptr

      upbound=MAXSTACK
      if (istackptr+4.ge.upbound) then
       print *,"stack overflow"
       stop
      else
       istack(istackptr)=ix
       istack(istackptr+1)=iy
       istack(istackptr+2)=iz
       istack(istackptr+3)=imult
       istackptr=istackptr+4
      endif
 
      return
      end
 
      subroutine popdata(ix,iy,iz,imult)
      IMPLICIT NONE

      INTEGER_T imult,istack(MAXSTACK),istackptr,ix,iy,iz

      common /stackstuff/ istack,istackptr

      if (istackptr-4.lt.1) then      
       print *,"stack underflow in popdata"   
       stop
      else
       imult=istack(istackptr-1)
       iz=istack(istackptr-2)
       iy=istack(istackptr-3)
       ix=istack(istackptr-4)
       istackptr=istackptr-4
      endif

      return
      end
 
      subroutine stackvolume(x,y,dxin,volfrac,ivapor)
      IMPLICIT NONE
      REAL_T x,y,z,dxin(SDIM),volfrac,lcenter,dxmax,hcutoff
      REAL_T xcenter,ycenter,zcenter,hx,hy,hz
      REAL_T ldata(3,3),rval,rlow,rhigh,volsum(32),dxsub(SDIM),newdx(SDIM)
      INTEGER_T rz_flag,i,j,k,isub,isubiter,dir,istatus,imult,ihalf
      INTEGER_T ix,iy,iz,ivapor
      REAL_T dx,dy,dz
     
#include "probdata.H"

      rz_flag=levelrz

      dxmax=dxin(1)
      do dir=2,SDIM
       if (dxmax.lt.dxin(dir)) then
        dxmax=dxin(dir)
       endif
      enddo

c 8 for production runs, 12 for testing curvature
      isub=8

      call initstack()      
      do i=1,isub+1
       volsum(i)=zero
      enddo
      call pushdata(0,0,0,0)

      do isubiter=1,MAXSTACK 
       call popdata(ix,iy,iz,imult)
       do dir=1,SDIM
        dxsub(dir)=dxin(dir)/REAL(2**imult)
       enddo
       xcenter=x+dxin(1)*REAL(ix)/(REAL(2**isub))+half*dxsub(1)-
     &   half*dxin(1)
       ycenter=y+dxin(2)*REAL(iy)/(REAL(2**isub))+half*dxsub(2)-
     &   half*dxin(2)
       rlow=xcenter-half*dxsub(1)
       rhigh=xcenter+half*dxsub(1)

       if (ivapor.eq.1) then 
        call vapordist(xcenter,ycenter,lcenter,dxsub(1))
       else
        call selectdist(xcenter,ycenter,lcenter,dxin(1))
       endif

       call stackfull(istatus)
       if ( (abs(lcenter).le.dxmax/REAL(2**imult)).and.(imult.lt.isub).and.
     &      (istatus.eq.0)) then
        do dir=1,SDIM
         newdx(dir)=dxsub(dir)*half
        enddo
        ihalf=2**(isub-imult-1)
        call pushdata(ix,iy,0,imult+1)
        call pushdata(ix+ihalf,iy,0,imult+1)
        call pushdata(ix+ihalf,iy+ihalf,0,imult+1)
        call pushdata(ix,iy+ihalf,0,imult+1)
       else 
        rval=one
        if ((rz_flag.eq.1).or.(rz_flag.eq.2)) then
         rval=xcenter
        endif
        if (lcenter.ge.dxmax/REAL(2**imult)) then
         volsum(imult+1)=volsum(imult+1)+rval
        else if (lcenter.le.-dxmax/REAL(2**imult)) then
         volsum(imult+1)=volsum(imult+1)
        else 
         do i=1,3
         do j=1,3
          if (ivapor.eq.1) then
           call vapordist(xcenter+(i-2)*dxsub(1),ycenter+(j-2)*dxsub(2),
     &                    ldata(i,j),dxsub(1))
          else
           call selectdist(xcenter+(i-2)*dxsub(1),ycenter+(j-2)*dxsub(2),
     &                     ldata(i,j),dxin(1))
          endif
         enddo 
         enddo 

         dx=dxsub(1)
         dy=dxsub(2)
         dz=dxsub(SDIM)

         call getvolume(ldata,rz_flag,rlow,rhigh,volfrac,dx,dy,dz)
         volsum(imult+1)=volsum(imult+1)+volfrac*rval
        endif
       endif
  
       call stackempty(istatus)
       if (istatus.eq.1) then 
        goto 100
       endif
      enddo
100   continue

      volfrac=zero
      do i=isub+1,1,-1
       volfrac=volfrac+volsum(i)/REAL(2**(2*(i-1)))
      enddo
      if ((rz_flag.eq.1).or.(rz_flag.eq.2)) then
       volfrac=volfrac/x
      endif

      if (isubiter.ge.MAXSTACK) then
       print *,"maximum iterations exceeded!!!"
      endif

      if ((isubiter.gt.1).and.(1.eq.0)) then
       print *,"x,y,vol,isubiter ",x,y,volfrac,isubiter
      endif

      return
      end

      subroutine rfile_dist( x,y,dist,rfilemax )
      IMPLICIT NONE

#include "probdata.H"

c     * j= 0 ... = rfilemax.

      REAL_T x,y,dist

      REAL_T zstatic(0:300),rstatic(0:300) 

      common / rfiledata / zstatic,rstatic

      INTEGER_T  j,distset,rfilemax 
      REAL_T deltaz,hval,rval,zval,z1,z2,r1,r2,rmid
      REAL_T distmin,rcoeff,zcoeff,ccoeff,dist1,zint,rint
      REAL_T dist2,dist3,rmin,rmax

      hval=zstatic(rfilemax)
      deltaz=zstatic(1)-zstatic(0)
      rval=x
      zval=yblob-y

      distset=0
      distmin=zero
      do j=1,rfilemax
       z1=zstatic(j-1)
       z2=zstatic(j)
       r1=rstatic(j-1)
       r2=rstatic(j)
       if (r1.lt.r2) then
        rmin=r1
        rmax=r2
       else
        rmin=r2
        rmax=r1
       endif
       if (abs(r2-r1).gt.abs(z2-z1)) then
        rcoeff=(z2-z1)/(r2-r1)
        zcoeff=-one
       else
        rcoeff=-one
        zcoeff=(r2-r1)/(z2-z1)
       endif
       dist=sqrt(rcoeff**2 + zcoeff**2)
       rcoeff=rcoeff/dist
       zcoeff=zcoeff/dist
       ccoeff=-rcoeff*r1-zcoeff*z1
       dist1=rcoeff*rval+zcoeff*zval+ccoeff
       rint=rval-dist1*rcoeff
       zint=zval-dist1*zcoeff
       if ((zint.lt.z1-1.0E-10).or.(zint.gt.z2+1.0E-10).or.
     &     (rint.lt.rmin-1.0E-10).or.(rint.gt.rmax+1.0E-10)) then
        dist2=sqrt( (rval-r1)**2 + (zval-z1)**2 )
        dist3=sqrt( (rval-r2)**2 + (zval-z2)**2 )
        if (dist2.le.dist3) then
         dist1=dist2
        else
         dist1=dist3
        endif
       endif
       if ((abs(dist1).lt.distmin).or.(distset.eq.0)) then
        distmin=abs(dist1)
        distset=1
       endif
      enddo
      
      if (zval.ge.hval) then
       dist=distmin  
      else if (zval.le.zero) then
       dist=rval-rstatic(0)
      else
       j=int( zval/deltaz )
       z1=j*deltaz
       r1=rstatic(j)
       r2=rstatic(j+1) 
       rmid=r1+(zval-z1)*(r2-r1)/deltaz 
       dist=rval-rmid
       if (dist.lt.zero) then
        dist=-distmin
       else
        dist=distmin
       endif
      endif
        
      return
      end

      subroutine readrfile( rfile, error, rfilemax )
      IMPLICIT NONE

#include "probdata.H"

c     *  0..rfilemax

c     * Variables passed in ...

      character rfile*20
      INTEGER_T   error,rfilemax 
      
      REAL_T  zstatic(0:300) ,  rstatic(0:300)

      common / rfiledata / zstatic,rstatic

      INTEGER_T  j         

      print *,"will assume rfilemax=",rfilemax
 
c     open(7, file=rfile, form="formatted", status="old", err=900)

      do j = 0,rfilemax 
         
c       read(7,100) zstatic(j),rstatic(j)
        
      end do
 
c     close(7)
 
      do j=0,rfilemax
       print *,"zstatic,rstatic ",zstatic(j),rstatic(j)
      enddo

      error = 0
      
      return
      
c 100  format(' ', e20.14, e21.14 )
c 900  write(6,910) rfile
c 910  format(" Can't open file = ",a," stopping ...")

      error = -1
      
      return
      end

      subroutine getslopeparms(xstart,ystart,xend,yend)
      IMPLICIT NONE
      REAL_T xstart,xend,ystart,yend

      xstart=0.0
      xend=4.0
      ystart=1.0
      yend=1.0

      return
      end

      subroutine selectgeom(x,y,dist)
      IMPLICIT NONE
      REAL_T x,y,z,dist,dist1
      REAL_T xstart,xend,ystart,yend,slope,intercept,xint,yint
      REAL_T ylen,frontrad,xcenter,ycenter
      INTEGER_T igeom

#include "probdata.H"

      igeom=1

      if (probtype.eq.21) then
       ylen=yblob/5.0

       xstart=-xblob/two
       ystart=yblob/two-ylen/four
       xend=xstart+xblob*three/four
       yend=ystart+ylen
       frontrad=zero

c      xstart=xblob/four
c      ystart=yblob/two-ylen/four
c      xend=xstart+xblob/two
c      yend=ystart+ylen
c      frontrad=zero

       if ((x.le.xstart+frontrad).and.(x.ge.xstart).and.
     &     (frontrad.gt.zero)) then
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
       else if ((xend-x.lt.x-xstart).and.(xend-x.lt.y-ystart).and.
     &          (xend-x.lt.yend-y)) then
        dist=x-xend
       else if ((x-xstart.lt.y-ystart).and.(x-xstart.lt.yend-y)) then
        dist=xstart-x
       else if (y-ystart.lt.yend-y) then
        dist=ystart-y
       else 
        dist=y-yend
       endif

      else if (probtype.eq.32) then
       dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
      else if (probtype.eq.30) then
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
      else if (probtype.eq.33) then
       call getslopeparms(xstart,ystart,xend,yend)

       slope=(yend-ystart)/(xend-xstart)
       intercept=ystart-slope*xstart
       xint=(x+slope*(y-intercept))/(one+slope*slope)  
       yint=slope*xint+intercept
       dist=sqrt( (x-xint)**2 + (y-yint)**2 )
       if (y.lt.slope*x+intercept) then
        dist=-dist
       endif
      else if (probtype.eq.34) then
c xblob=yblob=zero,radblob=3/4 is used in 2d case
       dist=radblob-abs(x)
c      dist1=0.9-abs(y-1.0)
c      dist=min(dist,dist1)
      else
       print *,"probtype out of range in selectgeom"
      endif

      return
      end
  
      subroutine contactparms(ymax,contactvar,contactleft,anglevar)
      IMPLICIT NONE
      REAL_T ymax,contactvar,contactleft,anglevar

#include "probdata.H"

      ymax=vinletgas
      contactvar=sqrt(radblob**2-(ymax-yblob)**2)
      contactleft=-one
      anglevar=acos(contactvar/radblob)
      anglevar=180.0*anglevar/Pi + 90.0

      return
      end

      subroutine damdist(x,y,dist,igeom)
      IMPLICIT NONE

      REAL_T x,y,dist,y1,y2,xctr,yswap,xswap
      INTEGER_T igeom

#include "probdata.H"

      if (axis_dir.eq.3) then
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
 
       if ((y.le.y1).or.((y.le.y2).and.(xswap.le.xctr))) then
        if (xswap.gt.xctr) then
         dist=y-y1
        else if (y.gt.y1) then
         dist=y-y2
         if (dist.lt.xswap-xctr) then
          dist=xswap-xctr
         endif        
        else
         dist=-sqrt((xswap-xctr)**2+(y-y1)**2)
         if (dist.lt.y-y2) then
          dist=y-y2
         endif
        endif
       else if (xswap.lt.xctr) then
        dist=y-y2
       else if (y.lt.y2) then
        dist=y-y1
        if (dist.gt.xswap-xctr) then
         dist=xswap-xctr
        endif  
       else
        dist=sqrt((xswap-xctr)**2+(y-y2)**2)
        if (dist.gt.y-y1) then
         dist=y-y1
        endif
       endif
       if (igeom.eq.0) then
        dist=-dist
       endif
      else
       print *,"axis_dir invalid"
       stop
      endif

      return
      end

      subroutine vapordist(x,y,dist,hx)
      IMPLICIT NONE
      REAL_T x,y,z,dist,hx

#include "probdata.H"

      REAL_T NPT,HSB,NID,NOD,CHH,scaleCHH,VRAD,dist1,dist2
      REAL_T xmin,xmax,ymin,ymax,zmin,zmax,zz,temprad
      REAL_T m,b,aspect
      REAL_T dist3,dist4,dist5,dist6,dist7,dist8,dist9,dist10
      INTEGER_T igeom
      REAL_T costheta,sintheta,xprime,yprime,zprime,delta
      REAL_T A(0:2,1:2),t1,t2,h1,h2,t3
      REAL_T hugedist

      hugedist=99999.0

      igeom=0

      dist=hugedist
      if (probtype.eq.102) then
       call sloshing_dist(x,y,dist)
      else if (probtype.eq.531) then
       dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob
      else if (probtype.eq.101) then
c dist<0 inside the square
       if (axis_dir.eq.0) then
        call squaredist(x,y,xblob-radblob,xblob+radblob,yblob-radblob,
     &   yblob+radblob,dist)
       else if (axis_dir.eq.1) then
        h1=1.10*radblob
        call squaredist(x,y,xblob-radblob-h1,xblob+radblob-h1,yblob-radblob,
     &   yblob+radblob,dist)
        call squaredist(x,y,xblob-radblob+h1,xblob+radblob+h1,yblob-radblob,
     &   yblob+radblob,dist2)
        if (dist2.lt.dist) then
         dist=dist2
        endif
       endif
      else if (probtype.eq.58) then
       dist=yblob+radblob*cos(two*Pi*x/xblob)-y
      else if ((probtype.eq.55).and.(axis_dir.eq.4)) then
       dist=yblob2-y
      else if ((probtype.eq.55).and.(axis_dir.eq.3)) then
       dist=yblob+radblob*cos(two*Pi*x/xblob)-y
       dist=-dist
      else if ((probtype.eq.55).and.(axis_dir.eq.6)) then
       dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob
       dist=-dist
      else if ((probtype.eq.55).and.(axis_dir.eq.5)) then
       if (radblob2.gt.zero) then
        dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob
       else
        dist=sqrt( (x-xblob)**2 + 
     &   (y-yblob-radblob*cos(contactangle))**2 )-radblob
       endif
      else if (probtype.eq.55) then
       dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob
c for 1d Welch: xblob=yblob=radblob=0  axis_dir.eq.2
       if ((axis_dir.eq.1).or.(axis_dir.eq.2)) then
        if (adv_dir.eq.1) then
         dist1=x-zblob
        else if (adv_dir.eq.2) then
         dist1=y-zblob
        else
         print *,"adv_dir invalid"
         stop
        endif
        if (dist1.lt.dist) then
         dist=dist1
        endif
       endif
      elseif ((probtype.eq.63).or.(probtype.eq.64)) then
        if (y.ge.two*xblob10) then
         dist=half*xblob10-x
        else
         dist=half*xblob10-x+(two*xblob10-y)
        endif
      else if ((probtype.eq.1).and.(axis_dir.eq.13)) then
       dist1=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
       dist2=hugedist
       dist3=hugedist
       dist4=hugedist
       dist5=hugedist
       dist6=hugedist
       dist7=hugedist
       dist8=hugedist
       dist9=hugedist
       dist10=hugedist
       if (radblob2.gt.zero) then
        dist2=sqrt((x-xblob2)**2+(y-yblob2)**2)-radblob2
        if (radblob3.gt.zero) then
         dist3=sqrt((x-xblob3)**2+(y-yblob3)**2)-radblob3
         if (radblob4.gt.zero) then
          dist4=sqrt((x-xblob4)**2+(y-yblob4)**2)-radblob4
          if (radblob5.gt.zero) then
           dist5=sqrt((x-xblob5)**2+(y-yblob5)**2)-radblob5
           if (radblob6.gt.zero) then
            dist6=sqrt((x-xblob6)**2+(y-yblob6)**2)-radblob6
            if (radblob7.gt.zero) then
             dist7=sqrt((x-xblob7)**2+(y-yblob7)**2)-radblob7
             if (radblob8.gt.zero) then
              dist8=sqrt((x-xblob8)**2+(y-yblob8)**2)-radblob8
              if (radblob9.gt.zero) then
               dist9=sqrt((x-xblob9)**2+(y-yblob9)**2)-radblob9
               if (radblob10.gt.zero) then
                dist10=sqrt((x-xblob10)**2+(y-yblob10)**2)-radblob10
               endif
              endif
             endif
            endif
           endif
          endif
         endif
        endif
       endif
       dist=min(dist1,dist2,dist3,dist4,dist5,dist6,dist7,dist8,dist9,
     %         dist10)
      else if ((probtype.eq.1).and.(axis_dir.eq.11)) then
       dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob
      else if ((probtype.eq.1).and.(axis_dir.eq.12)) then
       dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob
      else if ((probtype.eq.1).and.(axis_dir.eq.14)) then
       dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob
      else if ((probtype.eq.1).and.(axis_dir.eq.140)) then
       dist=max(-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob,
     &          -sqrt( (x-xblob)**2 + (y-yblob2)**2 )+radblob)
      else if ((probtype.eq.1).and.(axis_dir.eq.141)) then
       dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob
      else if (probtype.eq.5) then
       call ellipsedist(x,y,radblob,zblob,xblob,yblob,dist)
      else if (probtype.eq.13) then
       call legenddist(x,y,dist)
      else if ((probtype.eq.25).and.(axis_dir.eq.0)) then
       dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
       dist=-dist
      else if ((probtype.eq.25).and.(axis_dir.ge.1).and.
     &         (axis_dir.le.11)) then
       if (zblob.eq.zero) then
        aspect=hx
        call squaredist(x,y,-radblob,radblob,-hugedist,aspect,dist)
       else
        aspect=NOZZLEWIDTH*radblob
        if (aspect.lt.hx) then
         aspect=hx
        endif
        aspect=radblob+half*aspect
        call squaredist(x,y,-aspect,aspect,-zblob,zblob,dist)
       endif
      else if (probtype.eq.43) then
       xmin=-xblob-radblob
       xmax=xblob+radblob
       ymin=-yblob/two
       ymax=yblob/two
       if (zblob.gt.zero) then
        ymin=ymin-zblob/two
        ymax=ymax+zblob/two
       endif
       call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
      else if (probtype.eq.48) then
       dist1=y-yblob
       m=-one
       b=yblob-m*xblob
       dist2=y-(m*x+b)
       dist=dist1
       if (dist.lt.dist2) then
        dist=dist2
       endif
      else if ((probtype.eq.46).or.(probtype.eq.11)) then
       if ((probtype.eq.46).and.(radblob2.gt.zero)) then
        xmin=xblob-radblob3-radblob2+half*radblob6
        xmax=xblob+radblob3+radblob2-half*radblob6
        ymin=yblob-radblob4-radblob5+half*radblob6
        ymax=yblob+radblob4+radblob5-half*radblob6
        call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
       else
        temprad=radblob
        if (radblob.lt.hx) then
         temprad=hx
        endif
        xmin=-temprad
        xmax=temprad
        ymin=yblob-temprad
        ymax=yblob+temprad
        call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
       endif
      else if (probtype.eq.47) then
       if (radblob.ne.zero) then
         dist=radblob-abs(x-half*xblob)
       else
         print *,"radblob cannot equal zero!"
         stop
       endif
      else if (probtype.eq.35) then
       dist=half*yblob-y
      else if (probtype.eq.2) then
       dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
      else if ((probtype.eq.36).or.(probtype.eq.42)) then
       if (radblob.lt.hx) then
        print *,"radblob too small or base level too coarse"
        stop
       endif
       dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob

       if (probtype.eq.42) then
        dist1=zblob+yblob-y
        if (dist1.lt.dist) then
         dist=dist1
        endif
       endif
       if ((probtype.eq.36).and.(axis_dir.eq.100)) then
        dist=-dist
       endif
      else if (probtype.eq.9) then
       z=zero
       call geteblevel(x,y,z,1,dist)
      else if (probtype.eq.37) then
       dist=radblob-sqrt( (x-xblob)**2 + (y-yblob)**2 )
      else if ((probtype.eq.22).and.(axis_dir.eq.14)) then
       call initjetparms(HSB,NOD,NPT,NID,CHH,scaleCHH)
       VRAD=scaleCHH/eight
       dist=sqrt( (x-zero)**2 + (y-zero)**2 ) - VRAD  
      else if (probtype.eq.39) then
       dist=yblob+radblob*cos(two*Pi*x/xblob)-y
c Suvorov Formation of Taylor cone...
       if (radblob3.gt.zero) then
        dist=0.2*exp(-radblob3*x*x)-y+yblob3
       endif
      else if (probtype.eq.41) then
       if (axis_dir.eq.0) then
        dist=xblob+radblob*cos(two*Pi*y/yblob)-x
       else if ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &          (axis_dir.eq.3)) then
        call inletpipedist(x,y,dist)
       else
        print *,"axis_dir invalid"
        stop
       endif
      else if (probtype.eq.44) then
       call damdist(x,y,dist,igeom)
      else if (probtype.eq.45) then
       dist=yblob+1/(2*Pi)*(radblob*cos(2*Pi*(x-0.1)/xblob)+
     &    1/2*radblob**2*cos(4*Pi*(x-0.1)/xblob)+
     &    3/8*radblob**3*cos(6*Pi*(x-0.1)/xblob))-y
      else if (probtype.eq.21) then
       dist=yblob/two-y
      else if (probtype.eq.49) then
       dist=sqrt((x-xblob)**2+(y-yblob)**2) - radblob
       dist1= 1 - (y-yblob)
       if (dist1.lt.dist) then
        dist=dist1
       endif
      else if (probtype.eq.53) then
       ymin=zblob-1.0e+10
       ymax=zblob+two*hx
       xmin=xblob-radblob
       xmax=xblob+radblob
       call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
       dist=-dist 
      else if(probtype.eq.72) then
       xmin=xblob-radblob
       xmax=xblob+radblob
       ymin=zblob-1.0e+10
       ymax=zblob+two*hx
       call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
       dist=-dist
      else if (probtype.eq.54) then
       xmin=-xblob
       xmax=xblob
       ymin=-yblob
       ymax=yblob+half*yblob
       call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
       dist=-dist
      else if (probtype.eq.61) then
       dist=-(sqrt((x-xblob)**2+(y-yblob)**2) - radblob)
       if (axis_dir.eq.0) then
        dist1=-(3 + (y-yblob))
       else if (axis_dir.eq.1) then
        dist1=-(radblob + (y-yblob))
       else
        print *,"axis_dir invalid"
        stop
       endif
 
       if (dist1.gt.dist) then
        dist=dist1
       endif
      else if (probtype.eq.62) then
       dist=zblob3-y
       costheta=cos(xblob2)
       sintheta=sin(xblob2)
       xprime=costheta*(x-xblob)-sintheta*(y-yblob)
       yprime=sintheta*(x-xblob)+costheta*(y-yblob)
       delta=half*(radblob2-radblob)
       temprad=1.0E+10
       call squaredist(xprime,yprime,-radblob-delta,radblob+delta,
     &   -half*zblob2-delta,temprad,dist1)
       if (dist1.lt.dist) then
        dist=dist1
       endif 
      else if (probtype.eq.50) then
       xmin=-50.0
       xmax=half*twall
       ymin=-50.0
       ymax=50.0
       zmin=-50.0
       zmax=3.8
       zz=yblob
       call cubedist(xmin,xmax,ymin,ymax,zmin,zmax,
     &               x,zz,y,dist)
       dist=-dist

       xmin=-50.0
       xmax=50.0
       ymin=-50.0
       ymax=50.0
       zmin=-50.0
       zmax=0.2
       call cubedist(xmin,xmax,ymin,ymax,zmin,zmax,
     &               x,zz,y,dist1)
       dist1=-dist1

       if (dist1.gt.dist) then
        dist=dist1
       endif

      else if (probtype.eq.52) then
c      call handfree(x,yblob,y,dist)
       dist=zblob2-y
      else if (probtype.eq.56) then
       dist=zblob2-y
      else if (probtype.eq.66) then
       xprime=sqrt(three*radblob/(four*zblob*zblob*zblob))
       dist=zblob+radblob/(cosh(x*xprime)**2)-y
      else if (probtype.eq.38) then
       dist=yblob+radblob*cos(two*Pi*x/xblob)-y
      else if ((probtype.eq.26).and.(axis_dir.eq.1)) then
       dist=abs(y-half)-one/four
c natural convection in triangular enclosure
      else if (probtype.eq.81) then
       dist=yblob2-y
c rotating annulus
      else if (probtype.eq.82) then
       dist=hugedist
      endif

      return
      end

      subroutine jetgeom(x,y,dist)
      IMPLICIT NONE
      REAL_T x,y,z,dist
      REAL_T NPT,HSB,NID,NOD,CHH,scaleCHH
      REAL_T xx1(maxnumline),yy1(maxnumline)
      REAL_T xx2(maxnumline),yy2(maxnumline)
      INTEGER_T dd(maxnumline),lessflag(maxnumline),numline

#include "probdata.H"

      call initjetparms(HSB,NOD,NPT,NID,CHH,scaleCHH)

      xx1(1)=half*NOD
      xx2(1)=1.0E+10
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
       call igorgeom(x,y,dist)
      else if ((axis_dir.eq.11).or.(axis_dir.eq.12).or.
     &         (axis_dir.eq.13)) then
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
       if ((x.ge.half*NOD).and.(x.ge.half*NID).and.
     &     (y.le.HSB+NPT).and.(y.ge.HSB)) then
        dist=-abs(dist)
       endif
      endif

      return
      end

c negative on the inside of the triangle!
      subroutine triangledist(x,y,xlo,xhi,ylo,yhi,dist)
      IMPLICIT NONE

      REAL_T x,y,xlo,xhi,ylo,yhi,dist,dist1,dist2,dist3
      REAL_T m,b

      if ((xlo.ge.xhi-1.0E-10).or.(ylo.ge.yhi-1.0E-10)) then 
       print *,"invalid parameters ",xlo,xhi,ylo,yhi
       stop
      endif
      dist1=xlo-x
      dist2=ylo-y
      m=(yhi-ylo)/(xlo-xhi)
      b=ylo-m*xlo
      dist3=y-(m*x+b)
      dist=dist1
      if (dist2.gt.dist) then
       dist=dist2
      endif
      if (dist3.gt.dist) then
       dist=dist3
      endif
  
      return
      end

      subroutine airgunsolid(x,y,xcen,ycen,xhole,yhole,
     &                       height,width,gunthick,dist)
      IMPLICIT NONE

c        width/2
c ____    ____  
c  ___|  |__  |  
c  |        | | height/2
c _|        |_|          
c _          _  2yhole
c  |        | | 
c  |__    __| |
c  ___|  |____| gunthick  
c     2xhole
c Negative on the inside of the L-shaped solids!
c The subroutine below is intended to take advantage
c of axisymmetry so we only draw the right portion.
c
      REAL_T x,y,xcen,ycen,xhole,yhole
      REAL_T dist
      REAL_T height,width,gunthick
      REAL_T xlo1,xhi1,ylo1,yhi1
      REAL_T xlo2,xhi2,ylo2,yhi2
      REAL_T xlo3,xhi3,ylo3,yhi3
      REAL_T xlo4,xhi4,ylo4,yhi4
      
      if ( (xhole.lt.zero).or.(yhole.lt.zero).or.
     &    (height.le.zero).or.(width.le.zero).or.
     &    (gunthick.le.zero) ) then 
       print *,"These parameters must be positive ",xhole,
     &          yhole,height,width,gunthick
       stop
      endif
c bottom horizontal rectangle
      xlo1=xcen+xhole
      xhi1=xcen+xhole+half*width
      ylo1=ycen-half*height-yhole
      yhi1=ycen-half*height-yhole+gunthick
c top horizontal rectangle
      xlo2=xcen+xhole
      xhi2=xcen+xhole+half*width
      ylo2=ycen+half*height+yhole-gunthick
      yhi2=ycen+half*height+yhole
c bottom-right vertical rectangle
      xlo3=xcen+xhole+half*width-gunthick
      xhi3=xcen+xhole+half*width
      ylo3=ycen-half*height-yhole+gunthick
      yhi3=ycen-yhole
c top-right vertical rectangle
      xlo4=xcen+xhole+half*width-gunthick
      xhi4=xcen+xhole+half*width
      ylo4=ycen+yhole
      yhi4=ycen+half*height+yhole-gunthick
      
      if (y.le.yhi1) then
       call squaredist(x,y,xlo1,xhi1,ylo1,yhi1,dist)       
      endif
      if (y.gt.ylo2) then 
       call squaredist(x,y,xlo2,xhi2,ylo2,yhi2,dist)       
      endif
      if ((y.ge.ylo3).and.(y.le.(yhi3+yhole))) then
       call squaredist(x,y,xlo3,xhi3,ylo3,yhi3,dist)          
      endif
      if ((y.le.yhi4).and.(y.ge.(ylo4-yhole))) then
       call squaredist(x,y,xlo4,xhi4,ylo4,yhi4,dist)       
      endif
      
      return
      end

c negative on the inside of the square!
      subroutine squaredist(x,y,xlo,xhi,ylo,yhi,dist)
      IMPLICIT NONE

      REAL_T x,y,xlo,xhi,ylo,yhi,dist,dist1
      REAL_T xmid,ymid,xrad,yrad
 
      if ((xlo.ge.xhi-1.0E-10).or.(ylo.ge.yhi-1.0E-10)) then 
       print *,"invalid parameters ",xlo,xhi,ylo,yhi
       stop
      endif
      if ((x.le.xlo).and.(y.ge.ylo).and.(y.le.yhi)) then
       dist=xlo-x
      else if ((x.le.xlo).and.(y.ge.yhi)) then
       dist=sqrt( (x-xlo)**2 + (y-yhi)**2 )
      else if ((x.le.xlo).and.(y.le.ylo)) then
       dist=sqrt( (x-xlo)**2 + (y-ylo)**2 )
      else if ((x.ge.xhi).and.(y.ge.ylo).and.(y.le.yhi)) then
       dist=x-xhi
      else if ((x.ge.xhi).and.(y.ge.yhi)) then
       dist=sqrt( (x-xhi)**2 + (y-yhi)**2 )
      else if ((x.ge.xhi).and.(y.le.ylo)) then
       dist=sqrt( (x-xhi)**2 + (y-ylo)**2 )
      else if (y.ge.yhi) then
       dist=y-yhi
      else if (y.le.ylo) then
       dist=ylo-y
      else 
       xmid=half*(xlo+xhi)
       ymid=half*(ylo+yhi)
       xrad=half*(xhi-xlo)
       yrad=half*(yhi-ylo)

       dist=xrad-abs(x-xmid)
       dist1=yrad-abs(y-ymid)
       if ((dist.lt.zero).or.(dist1.lt.zero)) then
        print *,"dist,dist1 invalid",dist,dist1
        print *,"x,y,xlo,xhi,ylo,yhi ",x,y,xlo,xhi,ylo,yhi
        stop
       endif
       if (dist.lt.dist1) then
        dist=-dist
       else
        dist=-dist1
       endif
      endif

      return
      end

      subroutine nozzlerad(zval,radcross,rounded)
      IMPLICIT NONE
  
      REAL_T zval,radcross,rounded
#include "probdata.H"

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
      end
 
      subroutine soliddist(x,y,dist,hx,time)
      IMPLICIT NONE
      REAL_T x,y,z,dist,hx,time,dist1,temprad
      REAL_T rr,xx,yy,zz,aspect,offset,distplate
      REAL_T xlarge
      INTEGER_T igeom,onlypaddle
      REAL_T xsmall,ysmall,ylarge
      REAL_T costheta,sintheta,xprime,yprime
      REAL_T xmax,ymax,ydiag,yexit,radcross,dist2,A(0:2,1:2)
      REAL_T hugedist

#include "probdata.H"

      igeom=1

      hugedist=99999.0

      dist=hugedist

      if (probtype.eq.102) then
       call sloshing_geom(x,y,dist,time)
      else if (probtype.eq.42) then
       aspect=xblob2
c      offset=2.54d0
       offset=radblob2
       if (offset.lt.hx) then
        print *,"radblob2 (thickness plate) too small"
        print *,"or base grid too coarse"
        stop
       endif
       distplate=yblob2
       call squaredist(x,y,-aspect,aspect,yblob+distplate,
     &   yblob+distplate+offset,dist)
      else if ((probtype.eq.46).and.(radblob2.gt.zero)) then
       call airgunsolid(x,y,xblob,yblob,radblob3,radblob4,two*radblob5,
     &                  two*radblob2,radblob6,dist)
      
      elseif (probtype.eq.63) then
        call nozzlerad(y,radcross,zero)
  
        if (y.LE.2*xblob10) THEN
         dist=radcross-x
        elseif (x.LE.xblob10/2) THEN
         dist=dsqrt((x-xblob10/2)**two+(y-2*xblob10)**two)
        else 
         dist=y-2*xblob10
        endif
           
      elseif (probtype.eq.64) then
       call nozzlerad(y,radcross,radblob)
  
       if (y.LE.(two*xblob10-radblob)) then
        dist=radcross-x
       elseif (x.LE.(xblob10/two+radblob)) then
        dist=dsqrt((x-(xblob10/two+radblob))**two+
     &  (y-(two*xblob10-radblob))**two)-radblob
       else 
        dist=y-two*xblob10
       endif

      else if ((probtype.eq.43).and.(zblob.gt.zero)) then
       call squaredist(x,y,-xblob,xblob,yblob/two,yblob/two+zblob,dist) 
      else if ((probtype.eq.43).and.(radblob.gt.zero)) then
       offset=yblob/two+xblob/two
c      offset=1.1*yblob/two
       call squaredist(x,y,-radblob,radblob,-offset,offset,dist)
      else if (probtype.eq.48) then
       xlarge=1.0e+10
       call squaredist(x,y,xblob,xlarge,-yblob,yblob,dist)
      else if ((probtype.eq.25).and.(axis_dir.gt.0).and.(zblob.gt.0.0)) then
       aspect=NOZZLEWIDTH*radblob
       if (aspect.lt.hx) then
        aspect=hx
       endif
       aspect=aspect+radblob
    
       call squaredist(x,y,radblob,aspect,-zblob,zblob,dist)
      else if (probtype.eq.62) then
       costheta=cos(xblob2)
       sintheta=sin(xblob2)
       xprime=costheta*(x-xblob)-sintheta*(y-yblob)
       yprime=sintheta*(x-xblob)+costheta*(y-yblob)
       call squaredist(xprime,yprime,radblob,radblob2,
     &   -half*zblob2,half*zblob2,dist) 
       call squaredist(xprime,yprime,-radblob2,-radblob,
     &   -half*zblob2,half*zblob2,dist1) 
       if (dist1.lt.dist) then
        dist=dist1
       endif
      else if (probtype.eq.22) then
       call jetgeom(x,y,dist)
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
       call damdist(x,y,dist,igeom)
      else if ((probtype.eq.30).or.(probtype.eq.32).or.
     &         (probtype.eq.33).or.(probtype.eq.34).or.
     &         (probtype.eq.21)) then
       call selectgeom(x,y,dist)
      else if (probtype.eq.50) then
       zz=yblob
       onlypaddle=0
       call paddlegeom(x,zz,y,dist,onlypaddle)
      else if (probtype.eq.52) then
       call handgeom(x,yblob,y,dist,time)
      else if ((probtype.eq.55).and.(axis_dir.eq.5)) then
       dist=hugedist
      else if (probtype.eq.56) then
       call handgeom(x,yblob,y,dist,time)
      else if (probtype.eq.54) then
       print *,"obsolete"
       stop
      else if ( ((probtype.eq.9).and.(axis_dir.eq.1)).or.
     &          ((probtype.eq.45).and.(axis_dir.eq.1)) ) then
       if (probtype.eq.9) then
        zz=30.0
       else if (probtype.eq.45) then
        zz=0.5
       else
        print *,"probtype invalid"
        stop
       endif

       call rockgeom(x,zz,y,dist)
      else if ((probtype.eq.39).and.(radblob2.gt.zero).and.
     &         (zblob2.gt.zero).and.(radblob2.lt.half*xblob)) then
       xx=abs(x)
       if (y.ge.zblob2) then
        dist=sqrt(xx**2+(y-zblob2)**2)
       else
        yy=radblob2*(zblob2-y)/zblob2
        dist=xx-yy
       endif
c natural convection in triangular enclosure
      else if (probtype.eq.81) then
       dist=y+yblob3-yblob*(one-x/xblob)
c rotating annulus (cylindrical coordinates now)
      else if (probtype.eq.82) then
       if (1.eq.0) then
        dist=sqrt(x**2+y**2)-radblob
        dist1=radblob2-radblob-dist
        if (dist.gt.dist1) then
         dist=dist1
        endif 
       else
        dist=99999.0
       endif
c droplet/sphere impact
      else if (probtype.eq.531) then
       dist2=zero
       if (time.le.yblob3) then
        dist2=xblob3*time
       else
        dist2=xblob3*yblob3
       endif
       dist=sqrt((x-xblob2)**2+(y-yblob2-dist2)**2)-radblob2
      endif
   
      return
      end

#if (1==0)
      subroutine 2ddist(xmax,ymax,x,y,t1,t2,h1,h2,t3)
      IMPLICIT NONE

      REAL_T A(0:2,1:2),t1,t2,h1,h2,t3,dist2,xmax,ymax

c            c**** Section 3
               IF (y.GT.0.AND.y.LE.t1.AND.
     &        x.GE.h2.AND.x.LE.xmax) THEN
               dist=(h2-x)
c            c**** Section 4
            ELSEIF (y.GT.t1.AND.y.LE.x+t1-h2
     & .AND.x.GT.h2.AND.x.LE.xmax) THEN
               dist=-dsqrt((y-t1)**2.0d0+(x-h2)**2.0d0)
c            c**** Section 5
            ELSEIF (y.GE.x+t1-h2.AND.y.LE.
     & ((x-h2)/tan(1.1781))+t2.AND.
     & x.GE.-y+t2+h1.AND.x.LE.xmax) THEN
               A(1,1)=(t2-t1)**2.0d0
               A(1,2)=(t1-y)**2.0d0
               A(2,1)=(h1-h2)**2.0d0
               A(2,2)=(h1-x)**2.0d0
               A(0,1)=(h2-x)**2.0d0
               A(0,2)=(t2-y)**2.0d0
               dist=-dsqrt(A(0,1)+A(1,2)-
     & ((A(2,1)+A(1,1)+A(0,1)+A(1,2)-A(0,2)-A(2,2))**2.0d0)
     & /(4*(A(2,1)+A(1,1))))

c            c**** Section 6
            ELSEIF ((y.GE.((x-h2)/tan(1.1781))+t2).AND.
     & (y.LE.t3).AND.(x.GE.h1).AND.(x.LE.xmax)) THEN
               dist1=(h1-x)
               dist2=(y-t3)
               dist=Max(dist1,dist2)
c            c**** Section 7
            ELSEIF ((y.GT.t3).AND.(x.GT.h1).AND.
     & (x.LE.xmax).AND.(y.LE.ymax)) THEN
                dist=(y-t3)
c            c**** Section 8
            ELSEIF ((y.GE.t3).AND.(x.GE.0).AND.
     & (x.LE.h1).AND.(y.LE.ymax))  THEN
              dist=dsqrt((y-t3)**2.0d0+(x-h1)**2.0d0)
c            c**** Section 9
            ELSEIF ((y.GE.t2).AND.(y.LT.t3).AND.
     & (x.GE.0).AND.(x.LT.h1)) THEN
               dist=(h1-x)
c            c**** Section 10
            ELSEIF ((y.GE.(x+t2-h1)).AND.(y.LT.t2)
     & .AND.(x.GE.0)) THEN
               dist=dsqrt((y-t2)**2.0d0+(x-h1)**2.0d0)
c            c**** Section 11
            ELSEIF ((x.GE.0).AND.(((x.LE.(tan(1.1781)*
     & (y-t1)+h2)).AND.(y.LE.t1)).OR.((y.GE.t1).AND.
     & (x.GE.(y-t2+h1)).AND.(x.LE.(-y+t2+h1)))))THEN
               A(1,1)=(t2-t1)**2.0d0
               A(1,2)=(t1-y)**2.0d0
               A(2,1)=(h1-h2)**2.0d0
               A(2,2)=(h1-x)**2.0d0
               A(0,1)=(h2-x)**2.0d0
               A(0,2)=(t2-y)**2.0d0
               dist=dsqrt(A(0,1)+A(1,2)-((A(2,1)+A(1,1)+
     & A(0,1)+A(1,2)-A(0,2)-A(2,2))**2.0d0)/(4*(A(2,1)+A(1,1))))
c            c**** Section 12
            ELSEIF ((y.GT.0).AND.(y.LT.((x-h2)/tan(1.1781)
     & +t1)).AND.(x.GE.0).AND.(x.LT.h2)) THEN
               dist=(h2-x)
            endif         

      return
      end
#endif

      subroutine rockgeom(x,y,z,dist)
      IMPLICIT NONE

      REAL_T x,y,z,dist
      REAL_T xnew,ynew,znew
      REAL_T xf,yf,zf
      INTEGER_T i,j,k,dir
      REAL_T xx,yy,zz

      REAL_T  teblevel(IMAX_EB,IMAX_EB,IMAX_EB)
      INTEGER_T tidxmax(3)
      REAL_T  txmin(3),txmax(3),tdxboat(3)

      common / terrain_var / teblevel,txmin,txmax,tdxboat,tidxmax

 
#include "probdata.H"

      if ((probtype.eq.9).and.(axis_dir.eq.1)) then
       xnew=x-90.0
       ynew=y-30.0
       znew=z
      else if ((probtype.eq.45).and.(axis_dir.eq.1)) then
       xnew=x-0.5
       ynew=y-0.5
       znew=z
      else
       print *,"invalid in rockgeom"
       stop
      endif

      
      i=int ( (xnew-txmin(1))/tdxboat(1) ) + 1
      j=int ( (ynew-txmin(2))/tdxboat(2) ) + 1
      k=int ( (znew-txmin(3))/tdxboat(3) ) + 1
      xx=xnew
      yy=ynew
      zz=znew

      dir=1
      if (i.lt.1) then
       i=1
       xx=txmin(dir)
      endif
      if (i.ge.tidxmax(dir)) then
       i=tidxmax(dir)-1
       xx=txmax(dir)
      endif

      dir=2
      if (j.lt.1) then
       j=1
       yy=txmin(dir)
      endif
      if (j.ge.tidxmax(dir)) then
       j=tidxmax(dir)-1
       yy=txmax(dir)
      endif

      dir=3
      if (k.lt.1) then
       k=1
       zz=txmin(dir)
      endif
      if (k.ge.tidxmax(dir)) then
       k=tidxmax(dir)-1
       zz=txmax(dir)
      endif
      
      xf=(xx-(i-1)*tdxboat(1)-txmin(1))/tdxboat(1)
      yf=(yy-(j-1)*tdxboat(2)-txmin(2))/tdxboat(2)
      zf=(zz-(k-1)*tdxboat(3)-txmin(3))/tdxboat(3)
      dist=
     &  (one-xf)*(one-yf)*(one-zf)*teblevel(i,j,k)+
     &  (one-xf)*(one-yf)*(zf)*teblevel(i,j,k+1)+
     &  (one-xf)*(yf)*(one-zf)*teblevel(i,j+1,k)+
     &  (one-xf)*(yf)*(zf)*teblevel(i,j+1,k+1)+
     &  (xf)*(one-yf)*(one-zf)*teblevel(i+1,j,k)+
     &  (xf)*(one-yf)*(zf)*teblevel(i+1,j,k+1)+
     &  (xf)*(yf)*(one-zf)*teblevel(i+1,j+1,k)+
     &  (xf)*(yf)*(zf)*teblevel(i+1,j+1,k+1)
      return
      end

      subroutine vblobdist(x,y,dist)
      IMPLICIT NONE
      REAL_T x,y,z,dist,hx
      REAL_T xprime,yprime

#include "probdata.H"

      xprime=(x-xblob)*cos(contactangle)+(y-yblob)*sin(contactangle)
      yprime=-(x-xblob)*sin(contactangle)+(y-yblob)*cos(contactangle)
      yprime=yprime*radblob/zblob
      dist=sqrt(xprime**2+yprime**2)-radblob

      return
      end

c dist>0 inside of droplet
      subroutine legenddist(x,y,dist)
      IMPLICIT NONE
#include "probdata.H"

      REAL_T x,y,dist,mag,costheta,sintheta
      REAL_T cos2theta,sin2theta

      mag=sqrt((x-xblob)**2+(y-yblob)**2)
      if (mag.lt.1.0E-5) then
       dist=radblob+zblob
      else
       costheta=(y-yblob)/mag
       sintheta=(x-xblob)/mag
       cos2theta=costheta**2-sintheta**2
       sin2theta=two*costheta*sintheta
       if (levelrz.eq.0) then
        dist=radblob+zblob*(two*(costheta**2)-one)/two - mag
       else if (levelrz.eq.1) then 
        dist=radblob+zblob*(three*(costheta**2)-one)/two - mag
c       dist=radblob+zblob*(three*(cos2theta**2)-one)/two - mag
       else
        print *,"levelrz invalid"
        stop
       endif
      endif

      return
      end

c dist>0 inside of the ellipse
      subroutine ellipsedist(x,y,a,b,xc,yc,dist)
      IMPLICIT NONE

      REAL_T x,y,a,b,xc,yc,dist
      REAL_T xprime,yprime,xcritical,ycritical
      REAL_T factor

      xprime=x-xc
      yprime=y-yc
      if (xprime.lt.zero) then
       xprime=-xprime
      endif
      if (yprime.lt.zero) then
       yprime=-yprime
      endif
      if ((xprime.eq.zero).and.(yprime.eq.zero)) then
       dist=a
       if (dist.gt.b) then
        dist=b
       endif
      else if (xprime.gt.yprime) then
       factor=one/a**2 + (yprime/(xprime*b))**2
       xcritical=one/sqrt(factor)
       ycritical=yprime*xcritical/xprime
      else 
       factor=one/b**2 + (xprime/(yprime*a))**2
       ycritical=one/sqrt(factor)
       xcritical=xprime*ycritical/yprime
      endif
      dist=sqrt(xcritical**2+ycritical**2)-
     &     sqrt(xprime**2+yprime**2)

      return
      end
       
      subroutine selectdist(x,y,dist,hx)
      IMPLICIT NONE
      REAL_T x,y,z,dist,hx
      REAL_T  rho,dist1,dist2,rho1,rho2,hval
      INTEGER_T i,j,k,i1,j1,k1,ellipseflag,igeom
      REAL_T  deltax,deltay,xcritical,ycritical,factor
      REAL_T NPT,HSB,NID,NOD,CHH,scaleCHH
      REAL_T newradblob,opening,ylevel
      REAL_T xstart,xend,ystart,yend
      REAL_T height,theta,rval,zcenter
      REAL_T ymax
      REAL_T len,costheta,sintheta,xprime,yprime,localangle
      REAL_T xhit,yhit,slope,intercept,intercept1
      REAL_T sidelo,sidehi,xdiff,ydiff,aspect,hugedist

#include "probdata.H"

      igeom=0

      hugedist=99999.0
      dist=hugedist
      if ((probtype.eq.1).and.(axis_dir.eq.141)) then
       dist=-sqrt( (x-xblob)**2 + (y-yblob2)**2 )+radblob
      else if (probtype.eq.102) then
       dist=hugedist
      else if (probtype.eq.101) then
       dist=hugedist
      else if (probtype.eq.531) then
       dist=hugedist
      else if ((probtype.eq.1).and.(axis_dir.gt.10)) then
       dist=hugedist
      else if ((probtype.eq.1).or.(probtype.eq.11)) then
       dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob
      else if ((probtype.eq.63).or.(probtype.eq.64)) then
       dist=hugedist
      else if (probtype.eq.3) then
       dist=xblob+radblob*cos(two*Pi*y/yblob)-x
      else if (probtype.eq.4) then
       call rtdist(x,y,dist)
      else if (probtype.eq.7) then
       dist = radblob-sqrt((x-xblob)**2 + (y-yblob)**2)
       dist1 = yblob-two-y
       if (dist1.gt.dist) then
        dist=dist1
       endif
      else if (probtype.eq.8) then
       dist = -radblob+sqrt((x-xblob)**2 + (y-yblob)**2)
       dist1 = yblob+radblob-y+0.2
       if (dist1.lt.dist) then
        dist=dist1
       endif
      else if (probtype.eq.9) then
       dist=hugedist
      else if (probtype.eq.10) then
       print *,"option no longer used"
      else if (probtype.eq.12) then
       call legenddist(x,y,dist)
      else if (probtype.eq.13) then
       dist=hugedist
      else if (probtype.eq.14) then
       dist = sqrt((x-xblob)**2 + (y-yblob)**2)-radblob
      else if (probtype.eq.15) then
       print *,"fix me"
      else if (probtype.eq.16) then
        if (x.le.xblob+radblob) then
         dist=yblob-y-1.0E-3
        else
         dist=sqrt((x-xblob-radblob)*(x-xblob-radblob)+
     &             (yblob-y)*(yblob-y))
        endif
      else if (probtype.eq.17) then
       dist = radblob-sqrt((x-xblob)**2+(y-yblob-one)**2)
       dist1 = radblob-sqrt((x-xblob)**2+(y-yblob+one)**2)
       if (dist1.gt.dist) then
        dist=dist1
       endif
      else if (probtype.eq.18) then
       call initsplashparms(newradblob,ylevel,opening)
       call splashdist(x,y,newradblob,ylevel,opening,dist)
      else if (probtype.eq.19) then
       dist = radblob-sqrt((x-xblob)**2+(y-yblob)**2)
       if (contactangle.gt.two*radblob) then
        dist1=radblob*half-sqrt((x-xblob)**2+(y-yblob+contactangle)**2)
        if (dist1.gt.dist) then
         dist=dist1
        endif
       endif
      else if (probtype.eq.21) then
c      dist=yblob/two-y
       dist=hugedist
      else if (probtype.eq.22) then
       call initjetparms(HSB,NOD,NPT,NID,CHH,scaleCHH)
       call jettingdist(x,y,dist,HSB,NOD,NPT)
      else if (probtype.eq.23) then
       dist=yblob+radblob*cos(two*Pi*x/xblob)-y
      else if (probtype.eq.26) then
       dist=hugedist
      else if ((probtype.eq.24).or.(probtype.eq.27)) then
       dist=half-y
      else if (probtype.eq.25) then
       dist=hugedist
      else if (probtype.eq.28) then
       call zalesakdist(dist,x,y)
      else if (probtype.eq.29) then
       call deformdist(dist,x,y)
      else if ((probtype.eq.30).or.(probtype.eq.32)) then
c      dist=x
       dist=hugedist
      else if (probtype.eq.31) then
       dist=sqrt( (x-xblob)**2 + (y-yblob)**2 ) - radblob
      else if (probtype.eq.33) then
       call getslopeparms(xstart,ystart,xend,yend) 
       len=sqrt( (xend-xstart)**2 + (yend-ystart)**2 )
       costheta=(xend-xstart)/len
       sintheta=(yend-ystart)/len
       xprime=(x-xstart)*costheta+(y-ystart)*sintheta
       yprime=-(x-xstart)*sintheta+(y-ystart)*costheta
       xblob=half*(xstart+xend)
       yblob=zero
       radblob=one
       dist=sqrt( (xprime-xblob)**2 + (yprime-yblob)**2 )-radblob
      else if (probtype.eq.34) then
       localangle=contactangle
       localangle=two
       dist1=abs(x)
       if (dist1.gt.radblob) then
        height=1.0
       else if (localangle.gt.Pi/two+0.01) then
        theta=localangle-Pi/two
        rval=radblob/sin(theta)
        zcenter=1.0+rval*cos(theta)
        height=zcenter-sqrt(rval**2-dist1**2)
       else if (localangle.gt.Pi/two-0.01) then
        height=1.0
       else
        rval=radblob/cos(theta)
        zcenter=1.0-rval*sin(theta)
        height=zcenter+sqrt(rval**2-dist1**2)
       endif
       dist=height-y
      else if (probtype.eq.35) then
       dist=hugedist
      else if (probtype.eq.44) then
c      call damdist(x,y,dist,igeom)
       dist=hugedist
      else if ((probtype.eq.55).and.(axis_dir.eq.4)) then
       dist=hugedist
      else if (probtype.eq.46) then
        dist=zblob+yblob-y
      else if ((probtype.eq.36).or.(probtype.eq.37).or.
     &         (probtype.eq.39).or.(probtype.eq.41).or.
     &         (probtype.eq.42).or.(probtype.eq.43).or.
     &         (probtype.eq.44).or.(probtype.eq.45).or.
     &         (probtype.eq.11).or.
     &         (probtype.eq.49).or.(probtype.eq.53).or.
     &         (probtype.eq.72).or.
     &         (probtype.eq.61).or.(probtype.eq.5).or.
     &         (probtype.eq.47).or.(probtype.eq.48).or.
     &         (probtype.eq.54).or.(probtype.eq.55).or.
     &         (probtype.eq.62).or.(probtype.eq.58).or.
     &         (probtype.eq.81)) then
       dist=hugedist
      else if (probtype.eq.38) then
       dist=hugedist
      else if (probtype.eq.40) then
       call vblobdist(x,y,dist)
      else if ((probtype.eq.50).or.(probtype.eq.52).or.
     &         (probtype.eq.56).or.(probtype.eq.66)) then
       dist=hugedist
c rotating annulus
      else if (probtype.eq.82) then
       dist=99999.0
      else if (probtype.eq.2) then
       dist=99999.0
      else
       print *,"probtype out of range in selectdist: ",probtype
       stop
      endif

      return
      end


      subroutine FORT_INITUSOLID(nextra,usolid,DIMS(usolid),
     &  lo,hi,dxsub,time,xlo,xhi,
     &  xsolid,xpos,vrotation)
      IMPLICIT NONE
      INTEGER_T xsolid,xpos,nextra
      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T DIMDEC(usolid)
      REAL_T    time, dxsub(SDIM),vrotation
      REAL_T    xlo(SDIM), xhi(SDIM)
      REAL_T    usolid(DIMV(usolid),nextra)

      REAL_T  tstart,tfinish
      REAL_T timearr(0:IMAX_EB)
      INTEGER_T itimearr(0:IMAX_EB)
      INTEGER_T itimecount
      common / hand_var / tstart,tfinish,timearr,itimearr,itimecount

c     ::::: local variables
      INTEGER_T i, j,k,dir,i1,j1,k1,atfront1,atfront2,phase
      REAL_T  x, y,z
      REAL_T  dist1,dist2,xprime,yprime,zprime
      REAL_T  nn(3),dpdn,hx,hy,hz,vel,maxvel
      REAL_T  vofnew(3,3),vofold(3,3)

#include "probdata.H"

      call checkbound(lo,hi,DIMS(usolid),1,-1,1300)

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
       y = usolid(i,j,xpos+2)
       x = usolid(i,j,xpos+1)
       hx=half*(usolid(i+1,j,xpos+1)-usolid(i-1,j,xpos+1))
       hy=half*(usolid(i,j+1,xpos+2)-usolid(i,j-1,xpos+2))
       hz=hy

       do dir=1,SDIM
        usolid(i,j,xsolid+dir)=zero
       enddo
       if (probtype.eq.102) then
        call sloshing_tank_vel(time,usolid(i,j,xsolid+1))
c for pipe, free stream is y direction
       else if ((probtype.eq.41).and.(axis_dir.eq.2)) then
        usolid(i,j,xsolid+1)=zero
        phase=0
        call inletpipevel(usolid(i,j,xsolid+2),x,y,phase)
       else if ((probtype.eq.41).and.(axis_dir.eq.3)) then
        usolid(i,j,xsolid+1)=zero
        phase=1
        call inletpipevel(usolid(i,j,xsolid+2),x,y,phase)
       else if (probtype.eq.50) then
        print *,"obsolete"
        stop
       else if (probtype.eq.54) then
        usolid(i,j,xsolid+1)=zero
        usolid(i,j,xsolid+2)=advbot
       else if (probtype.eq.52) then
        call handvel(vel,time)
        usolid(i,j,xsolid+SDIM)=vel
       else if (probtype.eq.56) then
        if (tfinish-tstart.gt.1.0E-5) then
         maxvel=zero
         do i1=-1,1
         do j1=-1,1
          yprime = usolid(i+i1,j+j1,xpos+2)
          xprime = usolid(i+i1,j+j1,xpos+1)
          call handgeom(xprime,yblob,yprime,
     &           vofnew(i1+2,j1+2),tfinish)
          call handgeom(xprime,yblob,yprime,
     &           vofold(i1+2,j1+2),tstart)
          if (abs(vofnew(i1+2,j1+2)).gt.maxvel) then
           maxvel=abs(vofnew(i1+2,j1+2))
          endif
          if (abs(vofold(i1+2,j1+2)).gt.maxvel) then
           maxvel=abs(vofold(i1+2,j1+2))
          endif
         enddo 
         enddo 
         if (maxvel.eq.zero) then
          maxvel=one
         endif
         do i1=-1,1
         do j1=-1,1
          vofnew(i1+2,j1+2)=vofnew(i1+2,j1+2)/maxvel
          vofold(i1+2,j1+2)=vofold(i1+2,j1+2)/maxvel
         enddo
         enddo
         call fronterr(atfront1,vofold)
         call fronterr(atfront2,vofnew)
         if ((atfront1.eq.1).and.(atfront2.eq.1)) then
          nn(1)=(vofnew(3,2)-vofnew(1,2))/(two*hx)
          nn(2)=(vofnew(2,3)-vofnew(2,1))/(two*hy)
          nn(3)=zero
          dist1=sqrt(nn(1)**2+nn(2)**2+nn(3)**2)
          if (dist1.gt.1.0E-4) then
           maxvel=abs(adv_vel)
           dpdn=-(vofnew(2,2)-vofold(2,2))/
     &        ((tfinish-tstart)*dist1)
           if (dpdn.lt.-maxvel) then
            dpdn=-maxvel
           else if (dpdn.gt.maxvel) then
            dpdn=maxvel
           endif

           do dir=1,SDIM
            usolid(i,j,xsolid+dir)=usolid(i,j,xsolid+dir)+dpdn*nn(dir)/dist1
           enddo
          endif
         endif
c atfront
        endif  
c tfinish-tstart>0
       else if (probtype.eq.531) then
        if (time.le.yblob3) then
         usolid(i,j,xsolid+SDIM)=xblob3
        else
         usolid(i,j,xsolid+SDIM)=zero
        endif 
       endif
       if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
c       usolid(i,j,xsolid+1)=-vrotation*y
c       usolid(i,j,xsolid+2)=vrotation*x 
       endif    

      enddo
      enddo

      return
      end



      subroutine FORT_INITTEMPERATURE(xpos,temp,DIMS(temp),
     &  lo,hi,dx,time,xlo,xhi,phase,grav,tension,bpressure,
     &  bdensity)
      IMPLICIT NONE
      INTEGER_T    phase
      INTEGER_T    lo(SDIM),hi(SDIM)
      INTEGER_T    DIMDEC(temp)
      REAL_T     time, dx(SDIM),grav,tension
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T    temp(DIMV(temp))
      REAL_T    xpos(DIMV(temp),SDIM)
      REAL_T    density,bpressure
      REAL_T    bdensity
      REAL_T    dummyT

      INTEGER_T i, j,k
      REAL_T  x, y,z,rr
      REAL_T  hx, hy,hz,dist

#include "probdata.H"

      call checkbound(lo,hi,DIMS(temp),1,-1,1300)

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
       y = xpos(i,j,2)
       x = xpos(i,j,1)
       hx=half*(xpos(i+1,j,1)-xpos(i-1,j,1))
       hy=half*(xpos(i,j+1,2)-xpos(i,j-1,2))
       hz=hy

       temp(i,j)=zero
       if ((phase.eq.1).and.(solvevapor.eq.1)) then
        temp(i,j)=tvapor
       endif
       if ((phase.eq.0).and.(adiabatic.eq.0)) then
        temp(i,j)=twater
       endif
       if ((probtype.eq.55).and.(axis_dir.eq.4)) then 
        if (y.ge.yblob2) then
         temp(i,j)=twater
        else
         temp(i,j)=twater+(twater-tvapor)*(y-yblob2)/yblob2
        endif
       else if ((probtype.eq.55).and.(axis_dir.eq.6)) then
        if (phase.eq.0) then
         temp(i,j)=twater
        else if (phase.eq.1) then
         call vapordist(x,y,dist,hx)
         if (dist.ge.zero) then
          temp(i,j)=twater
         else
          temp(i,j)=tvapor
         endif
        else 
         print *,"phase invalid"
         stop
        endif
       else if ((probtype.eq.55).and.(axis_dir.eq.5)) then
        if (phase.eq.0) then
         temp(i,j)=twater
        else if (phase.eq.1) then
         temp(i,j)=tvapor
        else
         print *,"phase invalid"
         stop
        endif
       else if (probtype.eq.55) then
        if (phase.eq.0) then
         temp(i,j)=twater
        else if (phase.eq.1) then
         temp(i,j)=tvapor
c axis_dir.eq.2 1d Welch
         if (axis_dir.eq.2) then
          if (adv_dir.eq.1) then
           if (x.gt.zblob) then
            temp(i,j)=twater
           else
            temp(i,j)=tvapor+(twater-tvapor)*x/zblob
           endif
          else if (adv_dir.eq.2) then
           if (y.gt.zblob) then
            temp(i,j)=twater
           else
            temp(i,j)=tvapor+(twater-tvapor)*y/zblob
           endif
          else
           print *,"adv_dir invalid"
           stop
          endif
         else if (axis_dir.eq.3) then
          call vapordist(x,y,dist,hx)
          if (dist.ge.zero) then
           temp(i,j)=twater
          else
           temp(i,j)=tvapor+(twater-tvapor)*y/(y-dist)
          endif
         endif
        else 
         print *,"phase invalid"
         stop
        endif
       else if (probtype.eq.82) then
c       rr=sqrt(x*x+y*y)
        rr=x
        if (rr.le.radblob) then
         temp(i,j)=tsatdef
        else if (rr.lt.radblob2) then
         temp(i,j)=tsatdef+(walltemp-tsatdef)*(rr-radblob)/
     &               (radblob2-radblob)
        else
         temp(i,j)=walltemp
        endif
       endif

      enddo
      enddo

      return
      end



      subroutine FORT_INITDENSITY(xpos,density,DIMS(density),
     &  lo,hi,dx,time,xlo,xhi,phase,bdensity)
      IMPLICIT NONE
      INTEGER_T    phase
      INTEGER_T    lo(SDIM),hi(SDIM)
      INTEGER_T    DIMDEC(density)
      REAL_T     time, dx(SDIM),bdensity
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T    density(DIMV(density))
      REAL_T    xpos(DIMV(density),SDIM)

      INTEGER_T i, j,k
      REAL_T  x, y,z
      REAL_T  hx, hy,hz,taitden,hydropres,ptemp,dummyT

#include "probdata.H"

      call checkbound(lo,hi,DIMS(density),1,-1,1300)

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
       y = xpos(i,j,2)
       x = xpos(i,j,1)
       hx=half*(xpos(i+1,j,1)-xpos(i-1,j,1))
       hy=half*(xpos(i,j+1,2)-xpos(i,j-1,2))
       hz=hy

       dummyT=zero
    
       if (phase.eq.1) then
        density(i,j)=bdensity
       else if (phase.eq.0) then
        if ((probtype.eq.36).and.
     &      ((axis_dir.eq.4).or.(axis_dir.eq.2))) then
         ptemp=hydropres(y)
         density(i,j)=taitden(ptemp,dummyT)
c bubble jetting
        else if ((probtype.eq.42).and.(axis_dir.eq.1)) then
         density(i,j)=one+(denatdep-one)*(yblob+zblob-y)/zblob
c cavitation
        else if ((probtype.eq.46).or.(probtype.eq.11)) then
         density(i,j)=one+(denatdep-one)*(yblob+zblob-y)/zblob
        else if ((probtype.eq.43).or.(probtype.eq.48)) then
         ptemp=hydropres(y)
         density(i,j)=taitden(ptemp,dummyT)
        else if (probtype.eq.47) then
         density(i,j)=one
        else
         density(i,j)=denwater
        endif
       else
        print *,"phase invalid INITDENSITY"
        stop
       endif
      enddo
      enddo

      return
      end


      subroutine FORT_GETMAP(level,maxlevel,output,domlo,domhi,dir,
     &  problo,probhi)
      IMPLICIT NONE
 
      INTEGER_T level,maxlevel,domlo,domhi,dir
      REAL_T problo,probhi
      REAL_T output(domlo:domhi)
      REAL_T delta

      INTEGER_T i

      delta=(probhi-problo)/(domhi-domlo+one)

      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid"
       stop
      endif

      if (domlo.ne.0) then
       print *,"domlo invalid"
       stop
      endif

      do i=domlo,domhi
       output(i)=problo+(i+half)*delta
      enddo 

      return
      end
 
      subroutine FORT_LOADMAP(mf,DIMS(mf),lo,hi,
     &  ngrow,dir,time,level,maxlevel,problo,probhi,
     &  domlo,domhi)
      IMPLICIT NONE

      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T domlo(SDIM),domhi(SDIM)
      REAL_T problo(SDIM),probhi(SDIM)
      INTEGER_T DIMDEC(mf)
      INTEGER_T dir,ngrow,level,maxlevel
      REAL_T time
      REAL_T mf(DIMV(mf))

      INTEGER_T i,j,k
      REAL_T x,delta

      call checkbound(lo,hi,DIMS(mf),ngrow,-1,1001)
     
      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid in LOADMAP"
       stop
      endif
      do i=1,SDIM
       if (domlo(i).ne.0) then
        print *,"domlo invalid"
        stop
       endif
      enddo 

      delta=(probhi(dir+1)-problo(dir+1))/(domhi(dir+1)-domlo(dir+1)+one)
      do j = lo(2)-ngrow,hi(2)+ngrow
      do i = lo(1)-ngrow,hi(1)+ngrow
       if (dir.eq.0) then
        x=problo(dir+1)+(i+half)*delta
       else if (dir.eq.1) then
        x=problo(dir+1)+(j+half)*delta
       else 
        print *,"dir invalid"
        stop
       endif
       mf(i,j)=x
      enddo
      enddo

      return
      end

      subroutine FORT_FIXPARTPERIODIC(xpos,DIMS(xpos),
     &  mf,DIMS(mf),lo,hi,ngrow,scomp,ncomp,
     &  problen,problo,probhi,domlo,domhi,bc,numpart)
      IMPLICIT NONE

      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T domlo(SDIM),domhi(SDIM)
      INTEGER_T DIMDEC(xpos)
      INTEGER_T DIMDEC(mf)
      INTEGER_T scomp,ncomp,ngrow,numpart
      INTEGER_T bc(SDIM,2)
      REAL_T problen(SDIM)
      REAL_T problo(SDIM)
      REAL_T probhi(SDIM)
      REAL_T mf(DIMV(mf),ncomp)
      REAL_T xpos(DIMV(xpos),SDIM)

      INTEGER_T i,j,k,nn,dir,ibase

      if (scomp.ne.0) then
       print *,"scomp invalid"
       stop
      endif
      if (ncomp.ne.(SDIM+2)*numpart) then
       print *,"ncomp invalid"
       stop
      endif

      call checkbound(lo,hi,DIMS(mf),ngrow,-1,1001)
      call checkbound(lo,hi,DIMS(xpos),ngrow,-1,1001)
    
      do i=1,SDIM
       if (abs(probhi(i)-problo(i)-problen(i)).gt.1.0E-10) then
        print *,"something really bad with probhi,lo,len"
        stop
       endif
      enddo
 
      do nn=1,numpart
      do dir=0,SDIM-1
       ibase=(SDIM+2)*(nn-1)+dir+1

       if ((bc(1,1).eq.INT_DIR).and.(ARG_L1(mf).lt.domlo(1)).and.
     &    (dir.eq.0)) then
        do i=lo(1)-ngrow,domlo(1)-1
        do j=lo(2)-ngrow,hi(2)+ngrow
         mf(i,j,ibase)=mf(i,j,ibase)-problen(dir+1)
        enddo
        enddo
       endif

       if ((bc(1,2).eq.INT_DIR).and.(ARG_H1(mf).gt.domhi(1)).and.
     &     (dir.eq.0)) then
        do i=domhi(1)+1,hi(1)+ngrow
        do j=lo(2)-ngrow,hi(2)+ngrow
         mf(i,j,ibase)=mf(i,j,ibase)+problen(dir+1)
        enddo
        enddo
       endif

       if ((bc(2,1).eq.INT_DIR).and.(ARG_L2(mf).lt.domlo(2)).and.
     &     (dir.eq.1)) then
        do j=lo(2)-ngrow,domlo(2)-1
        do i=lo(1)-ngrow,hi(1)+ngrow
         mf(i,j,ibase)=mf(i,j,ibase)-problen(dir+1)
        enddo
        enddo
       endif

       if ((bc(2,2).eq.INT_DIR).and.(ARG_H2(mf).gt.domhi(2)).and.
     &     (dir.eq.1)) then
        do j=domhi(2)+1,hi(2)+ngrow
        do i=lo(1)-ngrow,hi(1)+ngrow
         mf(i,j,ibase)=mf(i,j,ibase)+problen(dir+1)
        enddo
        enddo
       endif

      enddo
      enddo

      return
      end

      subroutine FORT_INFLOWFIX(xpos,DIMS(xpos),mf,DIMS(mf),
     &  lo,hi,ngrow,state_idx,n,bc,vofgas,levelvapor,
     &  levelsolid,vofvapor,levelgas,xsolid,tempsolid,xvel,den,
     &  pres,temp,extraindex,time,domlo,domhi,vrotation)
      IMPLICIT NONE

      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T domlo(SDIM),domhi(SDIM)
      INTEGER_T DIMDEC(xpos)
      INTEGER_T DIMDEC(mf)
      INTEGER_T ngrow,state_idx,n,vofgas,levelvapor
      INTEGER_T levelsolid,vofvapor,levelgas,xsolid
      INTEGER_T tempsolid,xvel,den,pres,temp,extraindex
      INTEGER_T bc(SDIM,2)
      REAL_T time,dist,vrotation

      REAL_T xpos(DIMV(xpos),SDIM)
      REAL_T mf(DIMV(mf))

      INTEGER_T i,j,k,error
      REAL_T x,y,z,hx,hy,hz
      REAL_T x_vel,y_vel,z_vel
      REAL_T costheta,sintheta,xprime,yprime,zprime,rad
      REAL_T ydiff,ydiffcell
      REAL_T xzero,yzero,zzero

#include "probdata.H"

      call checkbound(lo,hi,DIMS(xpos),ngrow,-1,1001)
      call checkbound(lo,hi,DIMS(mf),ngrow,-1,1001)

      if (ngrow.le.0) then
       print *,"ngrow invalid"
       stop
      endif

      if ((state_idx.eq.0).or.(state_idx.eq.1)) then
c INFLOWFIX for horizontal velocity
       if (n.eq.xvel) then
        if (adv_dir.eq.1)then
         call rampvel(time,x_vel)
        else  
         x_vel = zero
        endif

        if ((bc(1,1).eq.EXT_DIR).and.(lo(1)-ngrow.lt.domlo(1))) then
         do i=lo(1)-ngrow,domlo(1)-1
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domlo(1),j,1)-xpos(domlo(1)-1,j,1)
          x=half*(xpos(domlo(1),j,1)+xpos(domlo(1)-1,j,1))
          y=xpos(i,j,2)
          if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
c          mf(i,j)=-vrotation*y
          endif
          if (probtype.eq.28) then
           call zalesakuu(mf(i,j),x,y)
          else if (probtype.eq.29) then
           call deformuu(mf(i,j),x,y,time,hx)
          else if (probtype.eq.31) then
           call circleuu(mf(i,j),x,y)
          else if ((probtype.eq.1).and.(axis_dir.eq.11)) then
           mf(i,j)=vinletgas*(y/yblob-one) 
          else if (probtype.eq.58) then
           mf(i,j)=zero
          endif
         enddo
         enddo 
        endif

c INFLOWFIX for horizontal velocity
        if ((bc(1,2).eq.EXT_DIR).and.(hi(1)+ngrow.gt.domhi(1))) then
         do i=domhi(1)+1,hi(1)+ngrow
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domhi(1)+1,j,1)-xpos(domhi(1),j,1)
          x=half*(xpos(domhi(1),j,1)+xpos(domhi(1)+1,j,1))
          y=xpos(i,j,2)
          if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
c          mf(i,j)=-vrotation*y
          endif
          if (probtype.eq.28) then
           call zalesakuu(mf(i,j),x,y)
          else if (probtype.eq.29) then
           call deformuu(mf(i,j),x,y,time,hx)
          else if (probtype.eq.31) then
           call circleuu(mf(i,j),x,y)
          else if ((probtype.eq.1).and.(axis_dir.eq.11)) then
           mf(i,j)=vinletgas*(y/yblob-one) 
          else if (probtype.eq.58) then
           mf(i,j)=zero
          else if ((probtype.eq.22).and.(axis_dir.eq.13)) then
           z=y
           y=zero
           call vbc(mf(i,j),time,z,y,error)
          endif
         enddo
         enddo 
        endif

c n.eq.xvel
        if ((bc(2,1).eq.EXT_DIR).and.(lo(2)-ngrow.lt.domlo(2))) then
         do j=lo(2)-ngrow,domlo(2)-1
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domlo(2),2)-xpos(i,domlo(2)-1,2)
          y=half*(xpos(i,domlo(2),2)+xpos(i,domlo(2)-1,2))
          x=xpos(i,j,1)
          if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
c          mf(i,j)=-vrotation*y
          endif
          if ((probtype.eq.1).and.(axis_dir.eq.13)) then
           mf(i,j)=zero
          else if ((probtype.eq.1).and.(axis_dir.ne.14)) then
           mf(i,j)=-vinletgas
          else if (probtype.eq.21) then
           mf(i,j)=adv_vel
          else if (probtype.eq.28) then
           call zalesakuu(mf(i,j),x,y)
          else if (probtype.eq.29) then
           call deformuu(mf(i,j),x,y,time,hy)
          else if (probtype.eq.31) then
           call circleuu(mf(i,j),x,y)
          else if (probtype.eq.58) then
           mf(i,j)=zero
          else if (probtype.eq.53) then
           mf(i,j)=zero
          else if (probtype.eq.10) then
           if ((x.ge.xblob-radblob).and.(x.le.xblob+radblob)) then
            mf(i,j)=zero
           endif
          endif
         enddo
         enddo 
        endif

        if ((bc(2,2).eq.EXT_DIR).and.(hi(2)+ngrow.gt.domhi(2))) then
         do j=domhi(2)+1,hi(2)+ngrow
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domhi(2)+1,2)-xpos(i,domhi(2),2)
          y=half*(xpos(i,domhi(2),2)+xpos(i,domhi(2)+1,2))
          x=xpos(i,j,1)
          if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
c          mf(i,j)=-vrotation*y
          endif
          if ((probtype.eq.1).and.(axis_dir.ne.14)) then
           mf(i,j)=vinletgas
          else if (probtype.eq.21) then
           mf(i,j)=adv_vel
          else if (probtype.eq.28) then
           call zalesakuu(mf(i,j),x,y)
          else if (probtype.eq.29) then
           call deformuu(mf(i,j),x,y,time,hy)
          else if (probtype.eq.31) then
           call circleuu(mf(i,j),x,y)
          else if (probtype.eq.58) then
           mf(i,j)=zero
          endif
         enddo
         enddo 
        endif
c INFLOWFIX for vertical velocity
       else if (n.eq.xvel+1) then
        if (adv_dir .eq. 2) then
         y_vel = adv_vel
        else
         y_vel = zero
        endif
        if (probtype.eq.53) then
         y_vel = zero
        endif

c INFLOWFIX for y-dir velocity
        if ((bc(1,1).eq.EXT_DIR).and.(lo(1)-ngrow.lt.domlo(1))) then
         do i=lo(1)-ngrow,domlo(1)-1
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domlo(1),j,1)-xpos(domlo(1)-1,j,1)
          x=half*(xpos(domlo(1),j,1)+xpos(domlo(1)-1,j,1))
          y=xpos(i,j,2)
          if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
           mf(i,j)=vrotation*x
          endif
          if (probtype.eq.28) then
           call zalesakvv(mf(i,j),x,y)
          else if (probtype.eq.29) then
           call deformvv(mf(i,j),x,y,time,hx)
          else if (probtype.eq.31) then
           call circlevv(mf(i,j),x,y)
c xlo wall gas pipe velocity
          else if ((probtype.eq.41).and.
     &             ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &              (axis_dir.eq.3))) then
           mf(i,j)=xblob5
          endif
         enddo
         enddo 
        endif

c INFLOWFIX for y-dir velocity
        if ((bc(1,2).eq.EXT_DIR).and.(hi(1)+ngrow.gt.domhi(1))) then
         do i=domhi(1)+1,hi(1)+ngrow
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domhi(1)+1,j,1)-xpos(domhi(1),j,1)
          x=half*(xpos(domhi(1),j,1)+xpos(domhi(1)+1,j,1))
          y=xpos(i,j,2)
          if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
           mf(i,j)=vrotation*x
          endif
          if (probtype.eq.28) then
           call zalesakvv(mf(i,j),x,y)
          else if (probtype.eq.29) then
           call deformvv(mf(i,j),x,y,time,hx)
          else if (probtype.eq.31) then
           call circlevv(mf(i,j),x,y)
          else if ((probtype.eq.36).and.(yblob10.ne.zero)) then
           mf(i,j)=yblob10
c wall liquid pipe velocity
          else if ((probtype.eq.41).and.
     &             ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &              (axis_dir.eq.3))) then
           mf(i,j)=yblob5
          endif
         enddo
         enddo 
        endif

c INFLOWFIX for y-dir velocity
c n.eq.xvel+1
        if ((bc(2,1).eq.EXT_DIR).and.(lo(2)-ngrow.lt.domlo(2))) then
         do j=lo(2)-ngrow,domlo(2)-1
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domlo(2),2)-xpos(i,domlo(2)-1,2)
          y=half*(xpos(i,domlo(2),2)+xpos(i,domlo(2)-1,2))
          x=xpos(i,j,1)
          if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
           mf(i,j)=vrotation*x
          endif
          if (probtype.eq.28) then
           call zalesakvv(mf(i,j),x,y)
          else if (probtype.eq.29) then
           call deformvv(mf(i,j),x,y,time,hy)
          else if (probtype.eq.31) then
           call circlevv(mf(i,j),x,y)
          else if ((probtype.eq.25).and.(axis_dir.gt.0)) then
           call inflow_bc(axis_dir,x,mf(i,j),hy)
          else if (probtype.eq.22) then
           mf(i,j)=y_vel*(xblob**2-x**2)*three/
     &          (two*xblob*xblob*xblob)
           if (x.gt.xblob) then
            mf(i,j)=zero
           endif
          else if (probtype.eq.53) then
           call get_orifice_velocity(x,mf(i,j))
          else if ((probtype.eq.63).or.(probtype.eq.64)) then
           mf(i,j)=y_vel
          else if (probtype.eq.10) then
           mf(i,j)=zero
           if ((x.ge.xblob-radblob).and.
     &         (x.le.xblob+radblob)) then
            mf(i,j)=advbot
           endif
          else if ((probtype.eq.36).and.(xblob10.gt.zero).and.
     &             (yblob10.ne.zero)) then
           mf(i,j)=x*yblob10/xblob10
          else if ((probtype.eq.41).and.
     &             ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &              (axis_dir.eq.3))) then
           yzero=zero
           call inletpipedist(x,yzero,dist)
           if (dist.lt.zero) then
c gas pipe velocity
            call inletpipevel(mf(i,j),x,y,1)
           else
c liquid pipe velocity
            call inletpipevel(mf(i,j),x,y,0)
           endif
          else if (probtype.eq.72) then
           call vapordist(x,y,dist,hx)
           if (dist.ge.zero) then
            mf(i,j)=advbot
           else
            mf(i,j)=adv_vel
           endif
          endif
         enddo
         enddo 
        endif

        if ((bc(2,2).eq.EXT_DIR).and.(hi(2)+ngrow.gt.domhi(2))) then
         do j=domhi(2)+1,hi(2)+ngrow
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domhi(2)+1,2)-xpos(i,domhi(2),2)
          y=half*(xpos(i,domhi(2),2)+xpos(i,domhi(2)+1,2))
          x=xpos(i,j,1)
          if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
           mf(i,j)=vrotation*x
          endif
          if (probtype.eq.28) then
           call zalesakvv(mf(i,j),x,y)
          else if (probtype.eq.29) then
           call deformvv(mf(i,j),x,y,time,hy)
          else if (probtype.eq.31) then
           call circlevv(mf(i,j),x,y)
          else if ((probtype.eq.16).or.
     &       ((probtype.eq.25).and.(axis_dir.eq.0)) ) then
           mf(i,j)=zero
           if ((x.ge.xblob-radblob).and.
     &         (x.le.xblob+radblob)) then
            mf(i,j)=-abs(advbot)
           endif
          else if (probtype.eq.62) then
           costheta=cos(xblob2)
           sintheta=sin(xblob2)
           xprime=costheta*(x-xblob)-sintheta*(y-yblob)
           yprime=sintheta*(x-xblob)+costheta*(y-yblob)
           if (xprime**2<radblob**2) then
            mf(i,j)=-radblob3
           endif
          else if ((probtype.eq.63).or.(probtype.eq.64)) then
           mf(i,j)=four*y_vel
          else if ((probtype.eq.36).and.(xblob10.gt.zero).and.
     &             (yblob10.ne.zero)) then
           mf(i,j)=x*yblob10/xblob10
          endif
         enddo
         enddo 
        endif
       endif
      else if (state_idx.eq.extraindex) then

c INFLOWFIX for old level set function
       if (n.eq.levelgas) then
        if ((bc(1,1).eq.EXT_DIR).and.(lo(1)-ngrow.lt.domlo(1))) then
         do i=lo(1)-ngrow,domlo(1)-1
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domlo(1),j,1)-xpos(domlo(1)-1,j,1)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          if (probtype.eq.6) then
           mf(i,j)=y-yblob
          endif
         enddo
         enddo 
        endif

        if ((bc(1,2).eq.EXT_DIR).and.(hi(1)+ngrow.gt.domhi(1))) then
         do i=domhi(1)+1,hi(1)+ngrow
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domhi(1)+1,j,1)-xpos(domhi(1),j,1)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
         enddo
         enddo 
        endif

        if ((bc(2,1).eq.EXT_DIR).and.(lo(2)-ngrow.lt.domlo(2))) then
         do j=lo(2)-ngrow,domlo(2)-1
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domlo(2),2)-xpos(i,domlo(2)-1,2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
         enddo
         enddo 
        endif

        if ((bc(2,2).eq.EXT_DIR).and.(hi(2)+ngrow.gt.domhi(2))) then
         do j=domhi(2)+1,hi(2)+ngrow
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domhi(2)+1,2)-xpos(i,domhi(2),2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
         enddo
         enddo 
        endif
c INFLOWFIX for old volume of fluid function
       else if (n.eq.vofgas) then

        if ((bc(1,1).eq.EXT_DIR).and.(lo(1)-ngrow.lt.domlo(1))) then
         do i=lo(1)-ngrow,domlo(1)-1
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domlo(1),j,1)-xpos(domlo(1)-1,j,1)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          if (probtype.eq.6) then
           mf(i,j)=zero
          endif
         enddo
         enddo 
        endif

        if ((bc(1,2).eq.EXT_DIR).and.(hi(1)+ngrow.gt.domhi(1))) then
         do i=domhi(1)+1,hi(1)+ngrow
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domhi(1)+1,j,1)-xpos(domhi(1),j,1)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
         enddo
         enddo 
        endif

        if ((bc(2,1).eq.EXT_DIR).and.(lo(2)-ngrow.lt.domlo(2))) then
         do j=lo(2)-ngrow,domlo(2)-1
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domlo(2),2)-xpos(i,domlo(2)-1,2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
         enddo
         enddo 
        endif

        if ((bc(2,2).eq.EXT_DIR).and.(hi(2)+ngrow.gt.domhi(2))) then
         do j=domhi(2)+1,hi(2)+ngrow
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domhi(2)+1,2)-xpos(i,domhi(2),2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
         enddo
         enddo 
        endif

c INFLOWFIX for new level set function
       else if (n.eq.levelvapor) then

        if ((bc(1,1).eq.EXT_DIR).and.(lo(1)-ngrow.lt.domlo(1))) then
         do i=lo(1)-ngrow,domlo(1)-1
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domlo(1),j,1)-xpos(domlo(1)-1,j,1)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          if (probtype.eq.58) then
           mf(i,j)=yblob2-y
          endif
         enddo
         enddo 
        endif

        if ((bc(1,2).eq.EXT_DIR).and.(hi(1)+ngrow.gt.domhi(1))) then
         do i=domhi(1)+1,hi(1)+ngrow
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domhi(1)+1,j,1)-xpos(domhi(1),j,1)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
         enddo
         enddo 
        endif

c INFLOWFIX for new level set function
        if ((bc(2,1).eq.EXT_DIR).and.(lo(2)-ngrow.lt.domlo(2))) then
         do j=lo(2)-ngrow,domlo(2)-1
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domlo(2),2)-xpos(i,domlo(2)-1,2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          if ((probtype.eq.25).and.(axis_dir.gt.0)) then
           mf(i,j)=abs(x)-radblob
          else if (probtype.eq.53) then
           mf(i,j)=radblob-abs(x-xblob)
          else if ((probtype.eq.41).and.
     &             ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &              (axis_dir.eq.3))) then
           yzero=zero
           call inletpipedist(x,yzero,dist)
           mf(i,j)=dist
          else if (probtype.eq.72) then
           mf(i,j)=radblob-abs(x-xblob)
          else if ((probtype.eq.55).and.(axis_dir.eq.3)) then
           if (mf(i,j).gt.y-radblob2) then
            mf(i,j)=y-radblob2
           endif
          endif
         enddo
         enddo 
        endif

c INFLOWFIX for new level set function (kluge)
        if ((bc(2,1).ne.INT_DIR).and.(lo(2)-ngrow.lt.domlo(2))) then
         do j=lo(2)-ngrow,domlo(2)-1
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domlo(2),2)-xpos(i,domlo(2)-1,2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          if ((probtype.eq.41).and.
     &        ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &         (axis_dir.eq.3))) then
           yzero=zero
           call inletpipedist(x,yzero,dist)
           mf(i,j)=dist
          endif
         enddo
         enddo 
        endif


c INFLOWFIX for new level set function
        if ((bc(2,2).eq.EXT_DIR).and.(hi(2)+ngrow.gt.domhi(2))) then
         do j=domhi(2)+1,hi(2)+ngrow
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domhi(2)+1,2)-xpos(i,domhi(2),2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          ydiff=y-half*(xpos(i,domhi(2),2)+xpos(i,domhi(2)+1,2))
          ydiffcell=y-xpos(i,domhi(2),2)
          if (probtype.eq.62) then
           mf(i,j)=zblob3-y
          else if ((probtype.eq.14).or.(probtype.eq.16).or.
     &      ((probtype.eq.25).and.(axis_dir.eq.0)) ) then
           if (probtype.eq.25) then
            mf(i,j)=radblob-x
           else if ((x.le.xblob+radblob).and.(probtype.eq.16)) then
            mf(i,j)=mf(i,domhi(2))-ydiffcell
           else
            mf(i,j)=two*mf(i,j-1)-mf(i,j-2)
           endif
          endif
         enddo
         enddo 
        endif

c INFLOWFIX for new volume of fluid function
       else if (n.eq.vofvapor) then

        if ((bc(1,1).eq.EXT_DIR).and.(lo(1)-ngrow.lt.domlo(1))) then
         do i=lo(1)-ngrow,domlo(1)-1
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domlo(1),j,1)-xpos(domlo(1)-1,j,1)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          if (probtype.eq.58) then
           if (y.ge.yblob2) then
            mf(i,j)=zero
           else if (y.le.yblob2) then
            mf(i,j)=one
           endif
          endif
         enddo
         enddo 
        endif

c INFLOWFIX for new volume of fluid function
        if ((bc(1,2).eq.EXT_DIR).and.(hi(1)+ngrow.gt.domhi(1))) then
         do i=domhi(1)+1,hi(1)+ngrow
         do j=lo(2)-ngrow,hi(2)+ngrow
          hx=xpos(domhi(1)+1,j,1)-xpos(domhi(1),j,1)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
         enddo
         enddo 
        endif

c INFLOWFIX for new volume of fluid function
        if ((bc(2,1).eq.EXT_DIR).and.(lo(2)-ngrow.lt.domlo(2))) then
         do j=lo(2)-ngrow,domlo(2)-1
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domlo(2),2)-xpos(i,domlo(2)-1,2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          if ((probtype.eq.25).and.(axis_dir.gt.0)) then
           if (abs(x).gt.radblob) then
            mf(i,j)=one
           else
            mf(i,j)=zero
           endif
          else if (probtype.eq.53) then
           if (abs(x-xblob).gt.radblob) then
            mf(i,j)=zero
           else 
            mf(i,j)=one
           endif
          else if ((probtype.eq.41).and.
     &             ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &              (axis_dir.eq.3))) then
           yzero=zero
           call inletpipedist(x,yzero,dist)
           if (dist.ge.zero) then
            mf(i,j)=one
           else
            mf(i,j)=zero
           endif
          else if (probtype.eq.72) then
           if (abs(x-xblob).gt.radblob) then
            mf(i,j)=zero
           else
            mf(i,j)=one
           endif
          else if ((probtype.eq.55).and.(axis_dir.eq.3)) then
           mf(i,j)=zero
          endif
         enddo
         enddo 
        endif


c INFLOWFIX for new volume of fluid function (kluge!)
        if ((bc(2,1).ne.INT_DIR).and.(lo(2)-ngrow.lt.domlo(2))) then
         do j=lo(2)-ngrow,domlo(2)-1
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domlo(2),2)-xpos(i,domlo(2)-1,2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          if ((probtype.eq.41).and.
     &        ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &         (axis_dir.eq.3))) then
           yzero=zero
           call inletpipedist(x,yzero,dist)
           if (dist.ge.zero) then
            mf(i,j)=one
           else
            mf(i,j)=zero
           endif
          endif
         enddo
         enddo 
        endif

c INFLOWFIX for new volume of fluid function
        if ((bc(2,2).eq.EXT_DIR).and.(hi(2)+ngrow.gt.domhi(2))) then
         do j=domhi(2)+1,hi(2)+ngrow
         do i=lo(1)-ngrow,hi(1)+ngrow
          hy=xpos(i,domhi(2)+1,2)-xpos(i,domhi(2),2)
          x=xpos(i,j,1)
          y=xpos(i,j,2)
          ydiff=y-half*(xpos(i,domhi(2),2)+xpos(i,domhi(2)+1,2))
          ydiffcell=y-xpos(i,domhi(2),2)
          if (probtype.eq.62) then
           mf(i,j)=zero
          else if ((probtype.eq.14).or.(probtype.eq.16).or.
     &      ((probtype.eq.25).and.(axis_dir.eq.0)) ) then
           if (probtype.eq.25) then
            mf(i,j)=one
           else if (probtype.eq.16) then
            mf(i,j)=zero
           endif
          endif
         enddo
         enddo 
        endif

       endif
      endif

      return
      end

      subroutine FORT_INITLEVELS(level,maxlevel,time,lo,hi,nscal,
     &	scal,DIMS(scal),dx,xlo,xhi,vofgas,levelvapor,levelsolid,
     &  vofvapor,levelgas,xpos,tsolid)
      IMPLICIT NONE
      INTEGER_T level,maxlevel,nscal
      INTEGER_T vofgas,levelvapor,levelsolid,vofvapor,levelgas,xpos 
      INTEGER_T tsolid
      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T DIMDEC(scal)
      REAL_T    time, dx(SDIM),dxmap(SDIM)
      REAL_T    xlo(SDIM), xhi(SDIM)
      REAL_T    scal(DIMV(scal),nscal)

c     ::::: local variables
      INTEGER_T i, j, k
      REAL_T  x, y, z
      REAL_T  hx, hy, hz

#include "probdata.H"

      call checkbound(lo,hi,DIMS(scal),1,-1,1300)
      hx = dx(1)
      hy = dx(2)

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
       y = scal(i,j,xpos+2)
       x = scal(i,j,xpos+1)
       dxmap(1)=half*(scal(i+1,j,xpos+1)-scal(i-1,j,xpos+1))
       dxmap(2)=half*(scal(i,j+1,xpos+2)-scal(i,j-1,xpos+2))

       call selectdist(x,y,scal(i,j,levelgas+1),hx)
       call stackvolume(x,y,dxmap,scal(i,j,vofgas+1),0)
       call vapordist(x,y,scal(i,j,levelvapor+1),dxmap(1))
       call stackvolume(x,y,dxmap,scal(i,j,vofvapor+1),1)
c      call tsolid_init(x,y,scal(i,j,tsolid+1))
      enddo
      enddo

      return
      end

      subroutine FORT_INITGEOM(nextra,lo,hi,dx,xlo,xhi,
     &	scal,DIMS(scal),time,levelsolid,xpos)
      IMPLICIT NONE
      INTEGER_T nextra,levelsolid,xpos
      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T DIMDEC(scal)
      REAL_T    xlo(SDIM), xhi(SDIM)
      REAL_T    dx(SDIM)
      REAL_T    scal(DIMV(scal),nextra)

      INTEGER_T i, j,k
      REAL_T  x, y,z,time
      REAL_T  hx, hy,hz

#include "probdata.H"

      call checkbound(lo,hi,DIMS(scal),1,-1,1300)

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
       y = scal(i,j,xpos+2)
       x = scal(i,j,xpos+1)
       hx=half*(scal(i+1,j,xpos+1)-scal(i-1,j,xpos+1))
       hy=half*(scal(i,j+1,xpos+2)-scal(i,j-1,xpos+2))
       hz=hy

       call soliddist(x,y,scal(i,j,levelsolid+1),hx,time)
      enddo
      enddo

      return
      end


      subroutine FORT_INITVELOCITY(level,time,lo,hi,
     &   xpos,vel,DIMS(vel),dx,xlo,xhi,phase,vrotation)
      IMPLICIT NONE
      INTEGER_T    level,phase
      INTEGER_T    lo(SDIM),hi(SDIM)
      INTEGER_T    DIMDEC(vel)
      REAL_T     time, dx(SDIM),vrotation
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     xpos(DIMV(vel),SDIM)
      REAL_T     vel(DIMV(vel),SDIM)

c     ::::: local variables
      INTEGER_T i,j,k
      REAL_T x,y,z
      REAL_T hx,hy,hz
      REAL_T x_vel,y_vel,z_vel,dist
      REAL_T xx_vel,yy_vel,zz_vel
      REAL_T xmax,ymax,zmax,contactvar,contactleft,anglevar
      REAL_T newradblob,ylevel,opening,ytop,radcross
      REAL_T s00,s01,s10,s11,xtemp

#include "probdata.H"

      call checkbound(lo,hi,DIMS(vel),1,-1,1300)

      if ((phase.ne.0).and.(phase.ne.1)) then
        print *,"phase invalid"
        stop
      endif

      print *,"levelrz=",levelrz
      if ((contactangle.lt.zero).or.(contactangle.gt.Pi+1.0E-2)) then
       print *,"contact angle (contactangle) out of range",contactangle
       stop
      endif
      print *,"contact angle (contactangle) is: ",contactangle

      if (adv_dir.eq.1) then
        xx_vel = adv_vel
        yy_vel = zero
        zz_vel = zero
      else if (adv_dir .eq. 2) then
        xx_vel = zero
        yy_vel = adv_vel
        zz_vel = zero
      else if ((adv_dir.eq.3).and.(SDIM.eq.3)) then
        xx_vel = zero
        yy_vel = zero
        zz_vel = adv_vel
      else if ((adv_dir.eq.3).and.(SDIM.eq.2)) then
        xx_vel = adv_vel
        yy_vel = adv_vel
        zz_vel = zero
      else if (probtype.ne.40) then
        print *,"adv_dir invalid: ",adv_dir
        stop
      endif

      if (adv_vel.ne.zero) then
        print *,"adv_dir,adv_vel ",adv_dir,adv_vel
      endif

c shear
      if (probtype.eq.1) then
        if (axis_dir.eq.0) then
         print *,"newtonian liquid"
        else if ((axis_dir.gt.0).and.(axis_dir.le.7)) then
         print *,"shear thinning liquid"
        else if (axis_dir.eq.11) then
         print *,"viscoelastic outer fluid"
        else if (axis_dir.eq.12) then
         print *,"viscoelastic drop"
        else if (axis_dir.eq.140) then
         print *,"droplets head on problem vinletgas=initial velocity"
        else if (axis_dir.eq.141) then
         print *,"diff. droplets head on problem vinletgas=initial velocity"
        else if (axis_dir.eq.13) then
         print *,"middle earth flow"
        else if (axis_dir.eq.14) then
         print *,"droplet collision problem vinletgas=initial velocity"
        else
         print *,"axis_dir invalid"
         stop
        endif
c bubble
      else if (probtype.eq.2) then
       if ((axis_dir.lt.0).or.(axis_dir.gt.7)) then
        print *,"axis_dir out of range in initbubble"
        stop
       else if (axis_dir.eq.0) then
        print *,"Newtonian liquid being computed...."
       else
        print *,"non-newtonian generalized cross carreau model liquid"
        print *,"axis_dir=",axis_dir
       endif
c capillary
      else if ((probtype.eq.3).or.(probtype.eq.41)) then
       if (axis_dir.eq.0) then
        print *,"RZ capillary problem"
        print *,"wavelen is yblob : ",yblob
        print *,"perturbation is radblob ",radblob
        print *,"base amplitude is xblob ",xblob
        print *,"r=xblob+radblob*cos(2pi y/yblob)"
        print *,"levelset < 0 in gas, levelset >0 in liquid"
        print *,"vfrac = 0 in gas, vfrac =1 in liquid"
       else if ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &          (axis_dir.eq.3)) then
        print *,"pipe flow problem"
        print *,"mean thickness of air region is xblob : ",xblob
        print *,"r=xblob+radblob*sin(2pi y/yblob)"
        print *,"initial inflow velocity of gas is advbot : ",advbot
        print *,"vertical velocity of gas at xlo : ",xblob5
        print *,"initial inflow velocity of liquid is vinletgas : ",vinletgas
        print *,"vertical velocity of water at xhi : ",yblob5
        print *,"set adv_dir=2"
        print *,"set adv_vel=0.0"
        print *,"levelset < 0 in gas, levelset >0 in liquid"
        if (axis_dir.eq.2) then
         print *,"outflow conditions at inlet for gas, inflow for liquid"
        endif
        if (axis_dir.eq.3) then
         print *,"outflow conditions at inlet for liquid, inflow for gas"
        endif
       else
        print *,"axis_dir invalid"
        stop
       endif
      else if (probtype.eq.4) then
       print *,"wavenumber is xblob : ",xblob
       print *,"y=radblob*cos(xblob*pi*x), xblob=2 for rt"
       print *,"xx_vel,yy_vel ",xx_vel,yy_vel
c gas burst
      else if (probtype.eq.8) then
       print *,"INITIALIZING RZ (axisym) GAS BURST PROBLEM "
      else if (probtype.eq.14) then
       call contactparms(ymax,contactvar,contactleft,anglevar)
       print *,"YMAX BETTER BE (vinletgas) ",vinletgas
       print *,"initial contactvar is: ",contactvar
       print *,"initial anglevar is: ",anglevar
      else if (probtype.eq.18) then
       print *, "initnosplash contactangle,advbot ",contactangle,advbot
       call initsplashparms(newradblob,ylevel,opening)
       print *,"newradblob,ylevel,opening",newradblob,ylevel,opening
c jetting 
      else if (probtype.eq.22) then
       print *,"jetting obselete"
c standing wave problem
      else if (probtype.eq.23) then
       print *,"standing wave problem (NOT r-z)"
       print *,"wavelen is xblob : ",xblob
       print *,"perturbation is radblob ",radblob
       print *,"base amplitude is yblob ",yblob
       print *,"y=yblob+radblob*cos(2pi x/xblob)"
       print *,"levelset < 0 in gas, levelset >0 in liquid"
       print *,"vfrac = 0 in gas, vfrac =1 in liquid"
c hanging
      else if (probtype.eq.25) then
       print *,"hanging drop problem or bubble column problem"
       print *,"axis_dir=1..11 if bubble column, axis_dir: ",axis_dir
       print *,"radius of orifice is radblob: ",radblob
       print *,"is axis_dir>0, zblob = height of column=",zblob
       print *,"otherwise zblob = radius preejected fluid=",zblob
       print *,"do not set zblob<0"
       print *,"advbot = rate water poured in =",advbot
       print *,"xblob should be 0, xblob=",xblob
       if (xblob.ne.0.0) then
        stop
       endif
       print *,"yblob=y value of inflow, yblob=",yblob
    
       if (axis_dir.eq.0) then 
        print *,"water density (denair) ",denair
        print *,"air density (denwater) ",denwater
        print *,"water viscosity (viscburn) ",viscburn
        print *,"air viscosity (viscunburn) ",viscunburn
       endif


       if ((axis_dir.gt.11).or.(axis_dir.lt.0)) then
        print *,"axis_dir out of range in inithanging"
       endif
       if ((axis_dir.gt.0).and.(zblob.lt.zero)) then
        print *,"zblob should be non-negative for bubble column problem"
        stop
       endif 
c shed
      else if ((probtype.eq.30).or.(probtype.eq.32).or.
     &         (probtype.eq.33).or.(probtype.eq.34) ) then
       print *,"probtype=30 means half circle, probtype=32 means full"
       print *,"probtype=33 means drop on a slope"
       print *,"probtype=34 means capillary tube"
       print *,"probtype=",probtype
c meniscus
      else if (probtype.eq.35) then
       print *,"radblob is NID/2 radblob= ",radblob
       print *,"yblob is NPT  yblob= ",yblob
       print *,"in 3d, xblob is domain base size xblob= ",xblob
      else if (probtype.eq.39) then
       print *,"standing wave problem (NOT r-z)"
       print *,"wavelen is xblob : ",xblob
       print *,"perturbation is radblob ",radblob
       print *,"base amplitude is yblob ",yblob
       print *,"y=yblob+radblob*cos(2pi x/xblob)"
       print *,"levelset < 0 in gas, levelset >0 in liquid"
       print *,"vfrac = 0 in gas, vfrac =1 in liquid"
      else if (probtype.eq.40) then
       if (adv_dir .eq. 1) then
         print *,"translation in x-direction with adv_vel=",adv_vel
       else if (adv_dir .eq. 2) then
         print *,"translation in y-direction with adv_vel=",adv_vel
       else if (adv_dir.eq.3) then
         print *,"translation in x and y-direction with adv_vel=",adv_vel
       else if (adv_dir.eq.4) then
         print *,"solid body rotation with adv_vel=",adv_vel
       else if (adv_dir.eq.5) then
         print *,"stretching with adv_vel=",adv_vel
       else
         write(6,*) "error: initvortpatch: adv_dir = ",adv_dir
         stop
       endif
c overturn
      else if (probtype.eq.45) then
       ytop=0.5
       print *,"using ytop=.5"
      endif
 
      do i = lo(1),hi(1)
      do j = lo(2),hi(2)
        x_vel=xx_vel
        y_vel=yy_vel
        z_vel=zz_vel

        y = xpos(i,j,2)
        x = xpos(i,j,1)
        hx=half*(xpos(i+1,j,1)-xpos(i-1,j,1))
        hy=half*(xpos(i,j+1,2)-xpos(i,j-1,2))
        hz=hy

        if (probtype.eq.1) then
         if (axis_dir.eq.11) then
          x_vel=vinletgas*(y/yblob-one)
         else if (axis_dir.eq.14) then
          dist=radblob-sqrt((x-xblob)**2+(y-yblob)**2)
          if (phase.eq.0) then
           x_vel=zero
           y_vel=vinletgas
          else if (phase.eq.1) then
           if (dist.ge.zero) then
            x_vel=zero
            y_vel=vinletgas
           endif
          else
           print *,"phase invalid"
           stop
          endif
         else if ((axis_dir.eq.140).or.(axis_dir.eq.141)) then
          dist=max(-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob,
     &             -sqrt( (x-xblob)**2 + (y-yblob2)**2 )+radblob)

          if (phase.eq.0) then
           x_vel=zero
           if( y > half*(yblob+yblob2))then
            y_vel=-abs(vinletgas)
           else
            y_vel=abs(vinletgas)
           endif
          else if (phase.eq.1) then
           if (dist.ge.zero) then
            x_vel=zero
            if( y > half*(yblob+yblob2))then
             y_vel=-abs(vinletgas)
            else
             y_vel=abs(vinletgas)
            endif
           endif
          else
           print *,"phase invalid"
           stop
          endif
         endif
        else if (probtype.eq.531) then
         dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob
         x_vel=zero
         y_vel=zero
         if (dist.ge.zero) then
          y_vel=vinletgas
         else
          y_vel=zero
         endif
        else if ((probtype.eq.3).or.(probtype.eq.41)) then
         if (axis_dir.eq.0) then
          if (phase.eq.0) then
          else if (phase.eq.1) then
           x_vel=zero
           y_vel=zero
          else
           print *,"phase invalid"
           stop
          endif
         else if ((axis_dir.eq.1).or.(axis_dir.eq.2).or.
     &            (axis_dir.eq.3)) then
          if (phase.eq.0) then
           x_vel=zero
           call inletpipevel(y_vel,x,y,phase)
          else if (phase.eq.1) then
           x_vel=zero
           call inletpipevel(y_vel,x,y,phase)
          else
           print *,"phase invalid"
           stop
          endif
         else
          print *,"axis_dir invalid"
          stop
         endif

c rotate
        else if (probtype.eq.5) then
         x_vel=zero
         y_vel=zero
c vortex patch
         if (phase.eq.0) then
c         x_vel=x_vel+adv_vel*(radblob**2)*(y-yblob)
c         y_vel=y_vel-adv_vel*(zblob**2)*(x-xblob)
          x_vel=x_vel+adv_vel*(y-yblob)
          y_vel=y_vel-adv_vel*(x-xblob)
         endif
c vortex stretching everywhere
         x_vel=x_vel+advbot*(x-xblob)
         y_vel=y_vel-advbot*(y-yblob)
c rotational vortex everywhere
         x_vel=x_vel+velfact*(y-yblob)
         y_vel=y_vel-velfact*(x-xblob)
c initbw
        else if (probtype.eq.9) then
         z=zero
         call geteblevel(x,y,z,2,x_vel)
         call geteblevel(x,y,z,3,y_vel)
c oilexpel
        else if (probtype.eq.16) then
         if ((y.gt.yblob).and.(x.le.xblob+radblob)) then
          y_vel = -abs(advbot)
         endif
        else if (probtype.eq.18) then
c dist<0 in air  dist>0 in water
c drop splash on pool of water
         call selectdist(x,y,dist,hx)
 
         x_vel = zero
         y_vel = zero
         if ((dist.gt.zero).and.(y.ge.ylevel)) then
          y_vel=-abs(advbot)
         endif
         call splashstream(x+hx/two,y+hy/two,newradblob,ylevel,
     &          opening,s11)
         call splashstream(x-hx/two,y+hy/two,newradblob,ylevel,
     &          opening,s01)
         call splashstream(x-hx/two,y-hy/two,newradblob,ylevel,
     &          opening,s00)
         call splashstream(x+hx/two,y-hy/two,newradblob,ylevel,
     &          opening,s10)
         x_vel=(s11+s01-s00-s10)/(two*hy)
         y_vel=-(s11+s10-s00-s01)/(two*hx)
        else if (probtype.eq.23) then
         call selectdist(x,y,dist,hx)
         if (dist.ge.zero) then
          x_vel=zero
          y_vel=zero
         endif
c validate
        else if (probtype.eq.24) then
         x_vel=-sin(Pi*x)*sin(Pi*x)*sin(two*Pi*y)
         y_vel=sin(Pi*y)*sin(Pi*y)*sin(two*Pi*x)
        else if (probtype.eq.25) then
         if (axis_dir.eq.0) then
          if ((y.gt.yblob).and.(x.le.xblob+radblob)) then
           y_vel = -abs(advbot)
          endif
         endif
c swirl
        else if (probtype.eq.26) then
         if (y.le.half) then
          x_vel=tanh( (y-one/four)*30.0 )
         else
          x_vel=tanh( (three/four-y)*30.0 )
         endif
         y_vel=0.05*sin(two*Pi*x)
        else if (probtype.eq.28) then
         call zalesakuu(x_vel,x,y)
         call zalesakvv(y_vel,x,y)
        else if (probtype.eq.29) then
         call deformuu(x_vel,x,y,zero,hy)
         call deformvv(y_vel,x,y,zero,hx)
c vbubble
        else if ((probtype.eq.36).or.(probtype.eq.37).or.
     &      (probtype.eq.42).or.(probtype.eq.46).or.
     &      (probtype.eq.11)) then
         if ((probtype.eq.36).and.
     &       ((axis_dir.eq.2).or.(axis_dir.eq.4))) then
          x_vel=zero
          y_vel=zero
         endif
         if ((probtype.eq.36).and.(xblob10.gt.zero).and.
     &       (yblob10.ne.zero)) then
          y_vel=x*yblob10/xblob10
         endif
c kh
        else if (probtype.eq.38) then
         call vapordist(x,y,dist,hx)
         if (dist.lt.zero) then
          y_vel=-x_vel*radblob*two*Pi*sin(two*Pi*x/xblob)/xblob
         else
          x_vel=zero
          y_vel=zero
         endif
c vstanding
        else if (probtype.eq.39) then
         if (phase.eq.0) then
          x_vel=zero
          y_vel=zero
         endif
c vortpatch
        else if (probtype.eq.40) then
         x_vel=zero
         y_vel=zero
c overturn
        else if (probtype.eq.45) then
         if (phase.eq.0) then
          x_vel=zblob*exp(-2*Pi*(ytop-y))*
     &        (1/(2*Pi)*(radblob*cos(2*Pi*(x-0.1)/xblob)+
     &        1/2*radblob**2*cos(4*Pi*(x-0.1)/xblob)+
     &        3/8*radblob**3*cos(6*Pi*(x-0.1)/xblob)))
          y_vel=zblob*exp(-2*Pi*(ytop-y))*
     &         (1/(2*Pi)*(radblob*sin(2*Pi*(x-0.1)/xblob)+
     &         1/2*radblob**2*sin(4*Pi*(x-0.1)/xblob)+
     &         3/8*radblob**3*sin(6*Pi*(x-0.1)/xblob)))
         else if (phase.eq.1) then
          x_vel=zero
          y_vel=zero
         else
          print *,"phase invalid"
          stop
         endif
c paddle
        else if (probtype.eq.50) then
         if (y.lt.zblob) then
          x_vel=zero
         endif
        else if (probtype.eq.58) then
         x_vel=zero
c jetbend
        else if (probtype.eq.53) then
         x_vel=zero
         call vapordist(x,y,dist,hx)
         if (dist.ge.zero) then
          call get_orifice_velocity(x,y_vel)
         else
          y_vel=zero
         endif
c 3D jet coaxial
        else if (probtype.eq.72) then
         call vapordist(x,y,dist,hx)
         if (dist.ge.zero) then
          x_vel=advbot
         else
          x_vel=adv_vel
         endif
c milkdrop
        else if (probtype.eq.61) then
         if (sqrt( (x-xblob)**2+(y-yblob)**2 ).le.radblob) then
          if (axis_dir.eq.1) then
           y_vel=-one
          endif
         endif
c nozzle
        else if ((probtype.eq.63).or.(probtype.eq.64)) then
         call nozzlerad(y,radcross,zero)
         if (x.gt.radcross) then
          y_vel=zero
         else
          y_vel=y_vel*(xblob10**2/radcross**2)
         endif
c pulse
        else if (probtype.eq.66) then
         xtemp=sqrt(three*radblob/(four*zblob*zblob*zblob))
         x_vel=sqrt(9.8*zblob)*(radblob/zblob)/(cosh(xtemp*x)**2)
         y_vel=sqrt(three*9.8*zblob)*((radblob/zblob)**(1.5))*
     &     (y/zblob)*tanh(xtemp*x)/(cosh(xtemp*x)**2)

        endif

        vel(i,j,1) = x_vel
        vel(i,j,2) = y_vel

        if ((vrotation.gt.zero).and.(axis_dir.eq.0)) then
c        vel(i,j,1)=-vrotation*y
         vel(i,j,2)=vrotation*x
        endif

      enddo
      enddo

      return
      end


      subroutine geteblevel(x,y,z,n,dist)
      IMPLICIT NONE
      INTEGER_T n
      REAL_T x,y,z,dist
      REAL_T xf,yf,zf

      REAL_T  eblevel(IMAX_EB,IMAX_EB,IMAX_EB,NCOMP_EB)
      INTEGER_T idxmax(3)
      REAL_T  xmin(3),xmax(3),dxboat(3),distmin,distmax

      common / eblevel_var / eblevel,xmin,xmax,dxboat,idxmax

      INTEGER_T i,j,k,dir
      REAL_T xx,yy,zz

      i=int ( (x-xmin(1))/dxboat(1) ) + 1
      j=int ( (y-xmin(2))/dxboat(2) ) + 1
      k=int ( (z-xmin(3))/dxboat(3) ) + 1
      xx=x
      yy=y
      zz=z

      dir=1
      if (i.lt.1) then
       i=1
       xx=xmin(dir)
      endif
      if (i.ge.idxmax(dir)) then
       i=idxmax(dir)-1
       xx=xmax(dir)
      endif

      dir=2
      if (j.lt.1) then
       j=1
       yy=xmin(dir)
      endif
      if (j.ge.idxmax(dir)) then
       j=idxmax(dir)-1
       yy=xmax(dir)
      endif

      dir=3
      if (k.lt.1) then
       k=1
       zz=xmin(dir)
      endif
      if (k.ge.idxmax(dir)) then
       k=idxmax(dir)-1
       zz=xmax(dir)
      endif
      
      xf=(xx-(i-1)*dxboat(1)-xmin(1))/dxboat(1)
      yf=(yy-(j-1)*dxboat(2)-xmin(2))/dxboat(2)
      zf=(zz-(k-1)*dxboat(3)-xmin(3))/dxboat(3)
      dist=
     &  (one-xf)*(one-yf)*(one-zf)*eblevel(i,j,k,n)+
     &  (one-xf)*(one-yf)*(zf)*eblevel(i,j,k+1,n)+
     &  (one-xf)*(yf)*(one-zf)*eblevel(i,j+1,k,n)+
     &  (one-xf)*(yf)*(zf)*eblevel(i,j+1,k+1,n)+
     &  (xf)*(one-yf)*(one-zf)*eblevel(i+1,j,k,n)+
     &  (xf)*(one-yf)*(zf)*eblevel(i+1,j,k+1,n)+
     &  (xf)*(yf)*(one-zf)*eblevel(i+1,j+1,k,n)+
     &  (xf)*(yf)*(zf)*eblevel(i+1,j+1,k+1,n)

      return
      end

c negative on the inside
      subroutine cubedist(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,dist)
      IMPLICIT NONE

      REAL_T xmin,xmax,ymin,ymax,zmin,zmax
      REAL_T x,y,z,dist
      REAL_T xcen,ycen,zcen,xrad,yrad,zrad
      REAL_T xdist,ydist,zdist

      xcen=half*(xmin+xmax)
      ycen=half*(ymin+ymax)
      zcen=half*(zmin+zmax)
      xrad=xmax-xcen
      yrad=ymax-ycen
      zrad=zmax-zcen

      xdist=abs(x-xcen)-xrad
      ydist=abs(y-ycen)-yrad
      zdist=abs(z-zcen)-zrad

      if ((xdist.le.zero).and.(ydist.le.zero).and.(zdist.le.zero)) then
       dist=xdist
       if (dist.lt.ydist) then
        dist=ydist
       endif
       if (dist.lt.zdist) then
        dist=zdist
       endif
      else
       if (xdist.lt.zero) then
        xdist=zero
       endif
       if (ydist.lt.zero) then
        ydist=zero
       endif
       if (zdist.lt.zero) then
        zdist=zero
       endif
       dist=sqrt(xdist**2+ydist**2+zdist**2)
      endif

      return
      end

c distance to star with center at origin

      subroutine stardist(x,y,z,radstar,radthick,dist)
      IMPLICIT NONE

      REAL_T x,y,z,dist
      REAL_T radstar,radthick,dist1,dist2

      call cubedist(-radstar,radstar,-radstar,radstar,
     &              -radthick,radthick,x,y,z,dist1)
      call cubedist(-radthick,radthick,-radstar,radstar,
     &              -radstar,radstar,x,y,z,dist2)

      if ((dist1.le.zero).and.(dist2.le.zero)) then
       dist=-sqrt(dist1**2+dist2**2)
      else if (dist1.le.zero) then
       dist=dist1
      else if (dist2.le.zero) then
       dist=dist2
      else
       dist=dist1
       if (dist.gt.dist2) then
        dist=dist2
       endif
      endif

      return
      end

c rotate is in radians, positive angle rotates clockwise.
      subroutine paddlegeom(x,y,z,dist,onlypaddle)
      IMPLICIT NONE

      REAL_T x,y,z,dist
      REAL_T xnew,ynew,znew,xrot,yrot,zrot
      REAL_T diststar,distcube,distcube1,distcube2,distcube3
      REAL_T xmin,xmax,ymin,ymax,zmin,zmax
      REAL_T xmin1,xmax1,ymin1,ymax1,zmin1,zmax1
      REAL_T xmin2,xmax2,ymin2,ymax2,zmin2,zmax2
      REAL_T xmin3,xmax3,ymin3,ymax3,zmin3,zmax3  
      INTEGER_T onlypaddle

#include "probdata.H"

      print *,"obsolete"
      stop

      return
      end

      subroutine FORT_GEOMINIT(time,dt)
      IMPLICIT NONE

      REAL_T time,dt
      character dwave*20
      INTEGER_T ewave,fwave

      character formtype*20
      INTEGER_T   error  

      REAL_T  tstart,tfinish
      REAL_T timearr(0:IMAX_EB)
      INTEGER_T itimearr(0:IMAX_EB)
      INTEGER_T itimecount
 
      REAL_T  eblevel(IMAX_EB,IMAX_EB,IMAX_EB,NCOMP_EB)
      INTEGER_T idxmax(3),dir
      REAL_T  xmin(3),xmax(3),dxboat(3)

      common / eblevel_var / eblevel,xmin,xmax,dxboat,idxmax

      common / hand_var / tstart,tfinish,timearr,itimearr,itimecount

      INTEGER_T i,j,k,i1,j1,k1
      REAL_T x,y,z,zlenwise,rcounter
   
      INTEGER_T iper,icrit
      REAL_T tper,tcrit
      INTEGER_T iread

#include "probdata.H"

      REAL_T    celvol  , dx      , dy      , dz
      REAL_T    dx2dy   , dy2dx   , dx2dz   , dz2dx   , dy2dz   , dz2dy
                                                      
      common / interface / celvol  , dx      , dy      , dz     ,
     &   dx2dy   , dy2dx   , dx2dz   , dz2dx   , dy2dz   , dz2dy

      INTEGER_T rz_flag
      REAL_T rval,ls(3,3,3),avg,cavg
      INTEGER_T ii,jj,kk,iter,maxiter,irb,i4two
                                                     
      rz_flag=0

      if ((probtype.eq.52).or.(probtype.eq.56)) then
       formtype="formatted"
       OPEN(unit=25,file="hcorresp.txt",access='sequential',
     &   form=formtype,status='old')
       READ(25,115) itimecount
       print *,"itimecount= ",itimecount
       if (itimecount.gt.IMAX_EB) then
        print *,"itimecount too big"
        stop
       endif

       do i1=1,itimecount
        READ(25,125) itimearr(i1),timearr(i1)
        print *,"index,time ",itimearr(i1),timearr(i1)
       enddo
       close(25)

       timearr(0)=zero
       itimearr(0)=itimearr(itimecount)

115    FORMAT(I4)
125    FORMAT(I4,f20.10)

       tper=time/timearr(itimecount)
       iper=INT(tper)
       tcrit=time-iper*timearr(itimecount)
       if (tcrit.gt.timearr(itimecount)+1.0E-8) then
        print *,"tcrit invalid"
        stop
       endif

       icrit=0
       do i1=1,itimecount
        if (tcrit.le.timearr(i1)+1.0E-8) then
         icrit=i1
         goto 135
        endif
       enddo
135    if (icrit.eq.0) then
        print *,"icrit invalid"
        stop
       endif

       if (dt.eq.zero) then
        tstart=-1.0
        tfinish=-1.0
       endif

       if ((tcrit.ge.tstart).and.(tcrit.le.tfinish)) then
        iread=0
       else
        iread=1
        tstart=timearr(icrit-1)
        tfinish=timearr(icrit)
        print *,"icrit,tcrit,tstart,tfinish ",icrit,tcrit,tstart,tfinish
       endif

       if (iread.eq.1) then

       do i1=1,2 
        if (i1.eq.1) then
         j=itimearr(icrit-1)
        else
         j=itimearr(icrit)
        endif
        if (j.lt.10) then
         if (probtype.eq.52) then
c         dwave = "h/hand00"//char(48+j)
          dwave = "h/hand000"//char(48+j)//".txt"
         else if (probtype.eq.56) then
          dwave = "h/swimmer000"//char(48+j)//".txt"
         endif
        else if (j.lt.100) then
         ewave = j/10
         if (probtype.eq.52) then
c         dwave = "h/hand0"//char(48+ewave)//char(48+j-ewave*10)
          dwave = "h/hand00"//char(48+ewave)//
     &            char(48+j-ewave*10)//".txt"
         else if (probtype.eq.56) then
          dwave = "h/swimmer00"//char(48+ewave)//
     &            char(48+j-ewave*10)//".txt"
         endif
        else
         ewave = j/100
         fwave = (j-ewave*100)/10
         if (probtype.eq.52) then
c         dwave = "h/hand"//char(48+ewave)//char(48+fwave)//
c    &             char(48+j-fwave*10-ewave*100)
          dwave = "h/hand0"//char(48+ewave)//char(48+fwave)//
     &             char(48+j-fwave*10-ewave*100)//".txt"
         else if (probtype.eq.56) then
          dwave = "h/swimmer0"//char(48+ewave)//char(48+fwave)//
     &             char(48+j-fwave*10-ewave*100)//".txt"
         endif
        endif

        formtype="formatted"
        print *,"opening ",dwave
        OPEN(unit=14,file=dwave,access='sequential',
     &       form=formtype,status='old')

        xmin(1)=-two*radblob
        xmax(1)=two*radblob
        xmin(2)=-two*radblob
        xmax(2)=two*radblob
        xmin(3)=-two*radblob
        xmax(3)=two*radblob
 
        idxmax(1)=128
        idxmax(2)=128
        idxmax(3)=128
        do dir=1,3
         dxboat(dir)=(xmax(dir)-xmin(dir))/REAL(idxmax(dir)-1)
        enddo 
        dx=dxboat(1)
        dy=dxboat(2)
        dz=dxboat(3)

        rcounter=zero
        do i=1,idxmax(1)
         do j=1,idxmax(2)
         do k=1,idxmax(3)
c         print *,"rcounter ",rcounter
          x=xmin(1)+REAL(i-1)*dxboat(1)
          y=xmin(2)+REAL(j-1)*dxboat(2)
          z=xmin(3)+REAL(k-1)*dxboat(3)
 
          rcounter=rcounter+one
          if (probtype.eq.52) then
           READ(14,100) eblevel(i,k,j,i1)
          else if (probtype.eq.56) then
           READ(14,100) eblevel(i,k,j,i1)
          endif
          eblevel(i,k,j,i1)=-radblob*eblevel(i,k,j,i1)
         enddo
         enddo
        enddo
        close(14)

        do i=1,idxmax(1)
        do j=1,idxmax(2)
        do k=1,idxmax(3)
         eblevel(i,k,j,i1+2)=zero
        enddo
        enddo
        enddo

100     FORMAT(f12.7)
       enddo
c i1

       endif
c iread.eq.1

      endif

      return
      end



      subroutine FORT_INPUTBOAT()
      IMPLICIT NONE

      character eblevel_file*20
      character formtype*20
      INTEGER_T   error  
     
      REAL_T  eblevel(IMAX_EB,IMAX_EB,IMAX_EB,NCOMP_EB)
      INTEGER_T idxmax(3),dir
      REAL_T  xmin(3),xmax(3),dxboat(3)

      REAL_T  teblevel(IMAX_EB,IMAX_EB,IMAX_EB)
      INTEGER_T tidxmax(3)
      REAL_T  txmin(3),txmax(3),tdxboat(3)

      REAL_T  lenscale,radthick,lenwidth

      common / eblevel_var / eblevel,xmin,xmax,dxboat,idxmax
      common / terrain_var / teblevel,txmin,txmax,tdxboat,tidxmax

      INTEGER_T count(IMAX_EB,IMAX_EB)
      INTEGER_T ebkwave(IMAX_EB,IMAX_EB)

      INTEGER_T factor,imax,index,k1vel,k2vel,kwave

      INTEGER_T i,j,k,i1,j1,k1
      REAL_T x,y,z,zlenwise
      REAL_T zkwave,x_vel,y_vel,z_vel
   
      REAL_T dist,f(199),g(199),s,m,t,fact
      REAL_T lpa,lpb
      REAL_T aa,bb,p(2),mm,a(2),b(2),pa(2),pb(2),ab(2)
      REAL_T length,temp(2),t1,t2,t3,t4
      COMPLEX*16 n,zed(199),II,d,fun(199)
      COMPLEX*16 avel(199),bvel(199),cvel(199)
      REAL_T fvel(199),gvel(199)
      REAL_T xx,yy,offset
      REAL_T time,dt

#include "probdata.H"

      time=zero
      dt=zero
      call GEOMINIT(time,dt)

      if (probtype.eq.50) then

       eblevel_file = "wheel.txt"
       formtype="formatted"
       if (axis_dir.eq.1) then
        print *,"opening ",eblevel_file
        OPEN(unit=13,file=eblevel_file,access='sequential',
     &   form=formtype,status='old')
       else if (axis_dir.ne.0) then
        print *,"axis_dir invalid"
        stop
       endif

       lenscale=radblob

       radthick=lenscale/8.0

       xmin(1)=-2.0*lenscale
       xmax(1)=2.0*lenscale
       xmin(2)=-2.0*lenscale
       xmax(2)=2.0*lenscale
       xmin(3)=-2.0*lenscale
       xmax(3)=2.0*lenscale
       idxmax(1)=128
       idxmax(2)=128
       idxmax(3)=128
       do dir=1,3
        dxboat(dir)=(xmax(dir)-xmin(dir))/REAL(idxmax(dir)-1)
       enddo 
       do i=1,idxmax(1)
       do j=1,idxmax(2)
       do k=1,idxmax(3)
        x=xmin(1)+REAL(i-1)*dxboat(1)
        y=xmin(2)+REAL(j-1)*dxboat(2)
        z=xmin(3)+REAL(k-1)*dxboat(3)
        if (axis_dir.eq.0) then
         call stardist(x,y,z,radblob,radthick,eblevel(i,j,k,1))
        else if (axis_dir.eq.1) then
         READ(13,100) eblevel(i,j,k,1)
         eblevel(i,j,k,1)=-lenscale*eblevel(i,j,k,1)
        else
         print *,"axis_dir invalid"
         stop
        endif

       enddo
       enddo
       enddo

100    FORMAT(f12.7)

       if (axis_dir.eq.1) then
        close(13)
       endif

      endif


      if ( ((probtype.eq.9).and.(axis_dir.eq.1)).or.
     &     ((probtype.eq.45).and.(axis_dir.eq.1)) ) then

       eblevel_file = "Terrain.txt"
       formtype="formatted"
       print *,"opening ",eblevel_file
       OPEN(unit=13,file=eblevel_file,access='sequential',
     &   form=formtype,status='old')

       if (probtype.eq.9) then
        lenscale=12.0
        lenwidth=100.0
       else if (probtype.eq.45) then
        lenscale=0.15
        lenwidth=0.5
       else
        print *,"probtype invalid"
        stop
       endif

       txmin(1)=-2.0*lenscale
       txmax(1)=2.0*lenscale
       txmin(2)=-2.0*lenscale
       txmax(2)=2.0*lenscale
       txmin(3)=-2.0*lenwidth
       txmax(3)=2.0*lenwidth

       tidxmax(1)=128
       tidxmax(2)=128
       tidxmax(3)=128
       do dir=1,3
        tdxboat(dir)=(txmax(dir)-txmin(dir))/REAL(tidxmax(dir)-1)
       enddo 
       do i=1,tidxmax(1)
       do j=1,tidxmax(2)
       do k=1,tidxmax(3)
        x=txmin(1)+REAL(i-1)*tdxboat(1)
        y=txmin(2)+REAL(j-1)*tdxboat(2)
        z=txmin(3)+REAL(k-1)*tdxboat(3)
        READ(13,131) teblevel(i,k,j)
        teblevel(i,k,j)=-lenscale*teblevel(i,k,j)
       enddo
       enddo
       enddo

131    FORMAT(f12.7)

       close(13)

      endif


      if (probtype.eq.9) then

       xmin(1)=-30.0
       xmax(1)=90.0
       xmin(2)=0.0
       xmax(2)=60.0
       xmin(3)=0.0
       xmax(3)=60.0

       idxmax(1)=128
       idxmax(2)=64
       idxmax(3)=4
       do dir=1,3
        dxboat(dir)=(xmax(dir)-xmin(dir))/REAL(idxmax(dir)-1)
       enddo 

       do k1=1,idxmax(3)
        print *,"initializing breaking wave slice,slicemax= ",k1,idxmax(3)

        zlenwise=xmin(3)+REAL(k1-1)*dxboat(3)
        t=zblob

        do i1=1,IMAX_EB
        do j1=1,IMAX_EB
         count(i1,j1)=0
        enddo
        enddo

        factor=6

        II=DCMPLX(cos(Pi/2),sin(Pi/2))
        call setupwave(zed,t)
        do i=1,199
	 f(i)=DREAL(zed(i))+30
	 g(i)=DIMAG(zed(i))+30	 
        enddo

        call setupwave(avel,t)
        call setupwave(bvel,t+0.01)
        do i = 1, 199
         cvel(i) = (bvel(i)-avel(i))/0.01
         fvel(i) = DREAL(cvel(i))
         gvel(i) = DIMAG(cvel(i))
        enddo

        do i = 1, 49
         fvel(50-i) = fvel(50)*(2*exp(i-50.0)-1)/4
         gvel(50-i) = gvel(50)*(50-i)/50
         gvel(i+150) = gvel(150)*exp(-(i-1.0)/10)
        enddo


        d=zed(150)-zed(50)
        aa=DIMAG(d)
        bb=DIMAG(zed(50))
        do i=1,199
         offset=REAL(i)*134.0/199.0
         xx=DREAL(zed(43))-offset
         yy=aa*exp(-(offset-90.0)**2/800.0) + bb
	 fun(i)=DCMPLX(xx,yy)
c	 print *,fun(i)
        enddo
      
        x=(radblob-one)/20.0
        k=1

        do i1=1,idxmax(1)
        do j1=1,idxmax(2)
         y=xmin(1)+REAL(i1-1)*dxboat(1)
         z=xmin(2)+REAL(j1-1)*dxboat(2)

 	 p(1)=y
	 p(2)=z
	 dist = 999999	

	 do i=1,199
	  a(1)=f(i)
	  a(2)=g(i)
	  b(1)=f(i+1)
	  b(2)=g(i+1)
     	  call distsub(p,a,b,mm)	
	  if (mm.lt.dist) then
	   dist = mm
	   k = i
	  endif
	 enddo
	
	 if (k.eq.199) then 
	  k=198
	 endif
	 n=zed(k+1)-zed(k)
	 s=(y-f(k))*DIMAG(n)-(z-g(k))*DREAL(n)
	 if (s.gt.0) then 
	  s = 1
	 else if (s.lt.0) then
	  s= -1
	 else
	  s= 0
	 endif
	 dist = -s*dist  
         eblevel(i1,j1,k1,1)=dist
         ebkwave(i1,j1)=k
        enddo
        enddo

        do i1=1,idxmax(1)
         index=1
         do j1=1,idxmax(2)
          y=xmin(1)+REAL(i1-1)*dxboat(1)
          z=xmin(2)+REAL(j1-1)*dxboat(2)

          dist=eblevel(i1,j1,k1,1)
          kwave=ebkwave(i1,j1)
	  y_vel = zero
	  z_vel = zero

          if (abs(dist).le.dxboat(1)) then
           y_vel = fvel(kwave)
	   z_vel = gvel(kwave)
  	   count(i1,1+index)=j1
           count(i1,1)=count(i1,1)+1
	   index=index+1
	  endif
          eblevel(i1,j1,k1,2) = y_vel/factor
          eblevel(i1,j1,k1,3) = z_vel/factor
         enddo
        enddo

        do i1=1,idxmax(1)
        do j1=1,idxmax(2)
         y=xmin(1)+REAL(i1-1)*dxboat(1)
         z=xmin(2)+REAL(j1-1)*dxboat(2)
      
         if (count(i1,1).eq.2) then
	  kwave=count(i1,2)
          zkwave = xmin(2)+REAL(kwave-1)*dxboat(2)
	  if (j1.le.kwave) then
           eblevel(i1,j1,k1,2)=eblevel(i1,kwave,k1,2)*
     &       (z-xmin(2))/(zkwave-xmin(2))
           eblevel(i1,j1,k1,3)=eblevel(i1,kwave,k1,3)*
     &       (z-xmin(2))/(zkwave-xmin(2))
	  else 
           eblevel(i1,j1,k1,2)=eblevel(i1,kwave+1,k1,2)*
     &       (xmax(2)-z)/(xmax(2)-zkwave)
           eblevel(i1,j1,k1,3)=eblevel(i1,kwave+1,k1,3)*
     &       (xmax(2)-z)/(xmax(2)-zkwave)
	  endif
	 else 
	  kwave=count(i1,2)
          zkwave = xmin(2)+REAL(kwave-1)*dxboat(2)
	  k1vel=count(i1,count(i1,1)+1)
	  k2vel=count(i1,count(i1,1)-1)		
          if (j1.lt.kwave) then
           eblevel(i1,j1,k1,2)=eblevel(i1,kwave,k1,2)*
     &       (z-xmin(2))/(zkwave-xmin(2))
           eblevel(i1,j1,k1,3)=eblevel(i1,kwave,k1,3)*
     &       (z-xmin(2))/(zkwave-xmin(2))
	  else if (j1.gt.k1vel) then
           eblevel(i1,j1,k1,2)=eblevel(i1,k1vel,k1,2)*
     &       (xmax(2)-z)/(xmax(2)-zkwave)
           eblevel(i1,j1,k1,3)=eblevel(i1,k1vel,k1,3)*
     &       (xmax(2)-z)/(xmax(2)-zkwave)
	  else if ((j1.ge.k2vel).and.(j1.lt.k1vel)) then
           eblevel(i1,j1,k1,2)=eblevel(i1,k2vel,k1,2)+
     &       (j1-k2vel)*
     &       (eblevel(i1,k1vel,k1,2)-eblevel(i1,k2vel,k1,2))/(k1vel-k2vel) 
           eblevel(i1,j1,k1,3)=eblevel(i1,k2vel,k1,3)+
     &       (j1-k2vel)*
     &       (eblevel(i1,k1vel,k1,3)-eblevel(i1,k2vel,k1,3))/(k1vel-k2vel) 
	  endif	
	 endif
	enddo
        enddo		

       enddo
c k1
      endif

      return
      end

      subroutine handoffset(ofs,time)
      IMPLICIT NONE

      REAL_T ofs,time

      ofs=zero
      if (time.le.0.4) then
       ofs=-time/8.0
      else if ((time.ge.0.4).and.(time.lt.0.8)) then
       ofs=-0.05
      else if ((time.ge.0.8).and.(time.lt.1.2)) then
       ofs=(time-0.9)/two
      else
       ofs=0.15
      endif

      return
      end
 
      subroutine handvel(vel,time)
      IMPLICIT NONE

      REAL_T vel,time

      vel=zero
      if (time.le.0.4) then
       vel=-one/eight
      else if ((time.ge.0.4).and.(time.lt.0.8)) then
       vel=zero
      else if ((time.ge.0.8).and.(time.lt.1.2)) then
       vel=half
      else
       vel=zero
      endif

      return
      end

      subroutine handgeom(x,y,z,dist,time)
      IMPLICIT NONE

      REAL_T x,y,z,dist,time
      REAL_T xnew,ynew,znew

      REAL_T  tstart,tfinish
      REAL_T timearr(0:IMAX_EB)
      INTEGER_T itimearr(0:IMAX_EB)
      INTEGER_T itimecount
      common / hand_var / tstart,tfinish,timearr,itimearr,itimecount
   
      REAL_T dist1,dist2
      INTEGER_T iper,icrit
      REAL_T tper,tcrit,theta,swapvar,ofs

#include "probdata.H"

      xnew=x-xblob
      ynew=y-yblob
      if ((probtype.eq.56).or.(probtype.eq.52)) then
       swapvar=xnew
       xnew=ynew
       ynew=swapvar
      endif

      if (probtype.eq.52) then
       call handoffset(ofs,time)
       znew=z-zblob-ofs
      else
       znew=z-zblob
      endif

      tper=time/timearr(itimecount)
      iper=INT(tper)
      tcrit=time-iper*timearr(itimecount)
      if (tcrit.gt.timearr(itimecount)+1.0E-8) then
       print *,"tcrit invalid"
       stop
      endif

      if (tcrit.le.tstart+1.0e-8) then
       theta=zero
      else if (tcrit.ge.tfinish-1.0e-8) then
       theta=one
      else
       theta=(tcrit-tstart)/(tfinish-tstart)
      endif

      call geteblevel(xnew,ynew,znew,1,dist1)
      call geteblevel(xnew,ynew,znew,2,dist2)

      dist=theta*dist2+(one-theta)*dist1

      return
      end

      subroutine handfree(x,y,z,dist)
      IMPLICIT NONE

      REAL_T x,y,z,dist,dist1
      REAL_T zthick
  
#include "probdata.H"

      zthick=9.0
      dist=sqrt( (x-twall)**2+(y-twall)**2 )-denfact
      dist1=zthick-z
      if (dist1.ge.zero) then
       if (dist.le.zero) then
        dist=-dist1
       else
        dist=sqrt(dist**2+dist1**2)
        dist=-dist
       endif
      else if (dist.ge.zero) then
       dist=-dist
      else
       dist=-dist
       dist1=-dist1
       if (dist.gt.dist1) then
        dist=dist1
       endif
      endif 

      return
      end


      subroutine deformuu(u,x,y,t,dx)
      IMPLICIT NONE
      REAL_T u,x,y,t,dx
      REAL_T aa,s1,s2,s3

#include "probdata.H"

      aa=cos(Pi*t/two)
      if (axis_dir.eq.0) then
       u=sin(four*Pi*(x+half))*sin(four*Pi*(y+half))*aa
      else if (axis_dir.eq.1) then
       s1=sin(Pi*y)
       s2=cos(Pi*y)
       s3=sin(Pi*x)
       u=-two*s1*s2*s3*s3
      else if (axis_dir.eq.2) then
       s1=sin(Pi*x)
       s2=sin(two*Pi*y)
       u=s1*s1*s2
       if (t.ge.one) then
        u=-u
       endif
      else
       print *,"axis_dir invalid"
       stop
      endif

      return
      end

      subroutine deformvv(v,x,y,t,dx)
      IMPLICIT NONE
      REAL_T v,x,y,t,dx
      REAL_T aa,s1,s2,s3

#include "probdata.H"

      aa=cos(Pi*t/two)
      if (axis_dir.eq.0) then
       v=cos(four*Pi*(x+half))*cos(four*Pi*(y+half))*aa
      else if (axis_dir.eq.1) then
       s1=sin(Pi*x)
       s2=cos(Pi*x)
       s3=sin(Pi*y)
       v=two*s1*s2*s3*s3
      else if (axis_dir.eq.2) then
       s1=sin(Pi*y)
       s2=sin(two*Pi*x)
       v=-s1*s1*s2
       if (t.ge.one) then
        v=-v
       endif
      else
       print *,"axis_dir invalid"
       stop
      endif

      return
      end

      subroutine deformdist(dist,x,y)
      IMPLICIT NONE
      REAL_T dist,x,y

      dist=sqrt( (x-half)**2 + (y-0.75)**2 )-0.15

      return
      end

      subroutine zalesakuu(u,x,y)
      IMPLICIT NONE
      REAL_T u,x,y

      u=-(Pi/314.0)*(y-50.0)

      return
      end

      subroutine zalesakvv(v,x,y)
      IMPLICIT NONE
      REAL_T v,x,y

      v=(Pi/314.0)*(x-50.0)

      return
      end

      subroutine zalesakdist(dist,x,y)
      IMPLICIT NONE
      REAL_T dist,x,y
      REAL_T dist1,dist2
    
      dist=sqrt((x-50.0)**2+(y-75.0)**2)-15.0
      if ((x.ge.47.5).and.(x.le.52.5)) then
       if (y.le.60.0) then
        if (x.lt.50.0) then
         dist=sqrt( (y-60.0)**2+(x-47.5)**2 )
        else
         dist=sqrt( (y-60.0)**2+(x-52.5)**2 )
        endif
       else if (y.le.85.0) then
        if (x.lt.50.0) then
         dist1=x-47.5
        else
         dist1=52.5-x
        endif
        dist2=85.0-y
        dist=min(dist1,dist2)
       else if ((y.le.90.0).and.(dist.le.zero)) then
        dist=max(dist,85.0-y)
       endif
      else if ((dist.lt.zero).and.(x.lt.47.5)) then
       if (y.le.85.0) then
        dist=max(dist,(x-47.5))
       else
        dist=max(dist,-sqrt( (x-47.5)**2+(y-85.0)**2 ) )
       endif
      else if ((dist.lt.zero).and.(x.gt.52.5)) then
       if (y.le.85.0) then
        dist=max(dist,(52.5-x))
       else
        dist=max(dist,-sqrt( (x-52.5)**2+(y-85.0)**2 ) )
       endif
      endif

      return
      end


      subroutine circleuu(u,x,y)
      IMPLICIT NONE
      REAL_T u,x,y

#include "probdata.H"

      u=zero
      if (adv_dir .eq. 1) then
         u = adv_vel
      else if (adv_dir .eq. 2) then
         u = zero
      else if (adv_dir.eq.3) then
         u = adv_vel
      endif

      return
      end

      subroutine circlevv(v,x,y)
      IMPLICIT NONE
      REAL_T v,x,y

#include "probdata.H"

      v=zero
      if (adv_dir .eq. 2) then
         v = adv_vel
      else if (adv_dir .eq. 1) then
         v = zero
      else if (adv_dir.eq.3) then
         v = adv_vel
      endif

      return
      end

      subroutine initcircle(level,time,lo,hi,
     &	 	            vel,DIMS(vel),
     &                      dx,xlo,xhi)
      IMPLICIT NONE

      INTEGER_T    level
      INTEGER_T    lo(SDIM), hi(SDIM)
      INTEGER_T    DIMDEC(vel)

      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     vel(DIMV(vel),SDIM)
      

c     ::::: local variables
      INTEGER_T i, j
      REAL_T  x, y, xtemp,ytemp
      REAL_T  hx, hy
      REAL_T  rho, dist, dist1,c, rhobub
      REAL_T  x_vel, y_vel
      REAL_T  rho1,rho2,hval
      INTEGER_T i1,j1,numsub,ellipseflag
      REAL_T xl1,xl2,yl1,yl2,xmin,ymin,height,lenarc,theta
      REAL_T x1,y1,x2,y2,area1,area2,area3,area4


#include "probdata.H"

      hx = dx(1)
      hy = dx(2)

      if (adv_dir .eq. 1) then
         x_vel = adv_vel
         y_vel = zero
      else if (adv_dir .eq. 2) then
         x_vel = zero
         y_vel = adv_vel
      else if (adv_dir.eq.3) then
         x_vel = adv_vel
         y_vel = adv_vel
      else
         write(6,*) "initcircle: adv_dir = ",adv_dir
         stop
      endif


         do j = lo(2), hi(2)
            y = xlo(2) + hy*(REAL(j-lo(2)) + half)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(REAL(i-lo(1)) + half)
               vel(i,j,1)=x_vel
               vel(i,j,2)=y_vel
	    enddo
         enddo


      return
      end

      subroutine jettingdist(x,y,dist,HSB,NOD,NPT)
      IMPLICIT NONE
      REAL_T x,y,dist,HSB,NOD,NPT,slope
      REAL_T xx1(maxnumline),yy1(maxnumline)
      REAL_T xx2(maxnumline),yy2(maxnumline)
      INTEGER_T dd(maxnumline),lessflag(maxnumline),numline

#include "probdata.H"

      xx1(1)=zero
      xx2(1)=half*NOD
      yy1(1)=HSB+NPT
      yy2(1)=HSB+NPT
      dd(1)=0
      lessflag(1)=0 

      slope = -0.4
      xx1(2)=half*NOD
      yy1(2)=HSB+NPT
      xx2(2)=1.0E+3
      yy2(2)=slope*(xx2(2)-xx1(2))+yy1(2)
      dd(2)=0
      lessflag(2)=0

      numline=2

      call construct(x,y,xx1,yy1,xx2,yy2,dd,lessflag,numline,dist)
      if (y.ge.HSB+NPT) then
       dist=-abs(dist)
      endif

      if ((axis_dir.ge.8).and.(axis_dir.le.10)) then
       call igordist(x,y,dist)
      else if ((axis_dir.eq.11).or.(axis_dir.eq.12).or.(axis_dir.eq.13)) then
       call microfabdist(x,y,dist)
      endif
    
      return
      end
 

      subroutine splashdist(xx,yy,radius,ylevel,opening,dist)
      IMPLICIT NONE
      REAL_T xx,yy,radius,ylevel,dist,dist1,dist2,opening

#include "probdata.H"

      dist=ylevel-yy
      dist1=sqrt((yy-ylevel)**2+(xx-opening)**2)
      dist2=radius-sqrt((xx-xblob)**2 + (yy-yblob)**2)

      if ((yy.le.ylevel).and.(abs(xx).ge.opening)) then
       dist=dist
      else if (yy.le.ylevel) then
       dist=dist1
      else if ((dist2.le.zero).and.(dist2.ge.dist)) then
       dist=dist2
      else if (dist2.le.zero) then
       dist=dist
      else if (yy.ge.yblob) then
       dist=dist2
      else if ( ((ylevel-yy)*(xx-xblob)/(yy-yblob)+xx).lt.opening) then
       dist=dist1
      else 
       dist=dist2
      endif
      return
      end

#define outerfactor (two)

      subroutine splashstream(xx,yy,radius,ylevel,opening,streamval)
      IMPLICIT NONE
      REAL_T xx,yy,radius,ylevel,dist,streamval,opening
      REAL_T dist1,dist2,dista,distb,distbb,xline,xrad

#include "probdata.H"

      call splashdist(xx,yy,radius,ylevel,opening,dist)

      dist1=sqrt((xx-xblob)**2+(yy-yblob)**2)
      distb=radius*(outerfactor+one)-dist1
      dista=dist1-radius

      if (yy.le.ylevel) then
       streamval=zero
      else if (dist.ge.zero) then
       streamval=abs(advbot)*xx
      else if (distb.le.zero) then
       streamval=zero
      else 
       xrad=radius*(xx-xblob)/dist1
       if (yy.le.yblob) then
        xline=(ylevel-yy)*(xx-xblob)/(yy-yblob) + xx
        distbb=sqrt((xline-xx)**2+(ylevel-yy)**2)
        if (distbb.lt.distb) then
         distb=distbb
        endif
       endif
       streamval=abs(advbot)*xrad*distb/(dista+distb)
      endif
      return 
      end


      subroutine partvolume(radius,contactangle,volume)
      IMPLICIT NONE
      REAL_T radius,contactangle,volume
      REAL_T hzero,temp1,r3,h3,r2,h2

      hzero=radius*sin(contactangle-Pi/two)
      r2=radius*radius
      r3=radius*r2
      h2=hzero*hzero
      h3=hzero*h2
      temp1=Pi*r2*(radius-hzero)-Pi*(r3 - h3)/three
      volume=four*Pi*r3/three - temp1
      return
      end
       

      subroutine length1(z,zout)
      IMPLICIT NONE
      REAL_T z(2),zout
      
      zout = sqrt(z(1)*z(1)+z(2)*z(2))
      return 

      end
      
      subroutine DOT1(a,b,c)
      IMPLICIT NONE
      REAL_T a(2),b(2),c
      
      c = a(1)*b(1)+a(2)*b(2)
      return 

      end    
      		  
      subroutine distsub(p,a,b,dist)
      IMPLICIT NONE
      REAL_T t,t1,t2,t3,t4,p(2),b(2),a(2)
      REAL_T dist,l,lpb,lpa,pb(2),pa(2),ab(2),temp(2)
      INTEGER_T p_inside_band,i
	
	do i=1,2
	  ab(i)=b(i)-a(i)
	  pb(i)=b(i)-p(i)
	  pa(i)=a(i)-p(i)
	enddo
c       one should now check if dot(pa,pb)/dot(pa,pa) is between 0 and 1
	call DOT1(pb,ab,t1) 
	call DOT1(ab,ab,t2)
	t = t1/t2	
c	t =DOT(pa,pb)/DOT(ab,ab)
	if ((t.gt.0).and.(t.lt.1))  then
	 call DOT1(pa,ab,t3)
	 do i=1,2
c        the second part of temp is the projection of p onto ab	  
	  temp(i) = t2*p(i)-t1*a(i)+t3*b(i)
c	  temp(i) =DOT(ab,ab)*p(i)-DOT(pb,ab)*a(i)+DOT(pa,ab)*b(i)
	 enddo
	 call length1(temp,t4)
	 dist = t4/t2
c	 dist = 1/DOT(ab,ab)*length(temp)
	else
	 call length1(pa,t1)
	 call length1(pb,t2)
         dist = min(t1,t2)
	endif	
	return 
      end
      
      
      subroutine setupwave(zed,t)
	IMPLICIT NONE
	REAL_T o(101),kk,t,ar(101),aa,bb
	INTEGER_T i,k
	COMPLEX*16 zed(199),A,B,A3,A4,A5,II
        COMPLEX*16 P0(101),P1(101),P2(101),Q0(101),Q1(101)
        COMPLEX*16 Q2(101),Q3(101),Q4(101),Q5(101), tang,d
#include "probdata.H"

        kk=20
	II=DCMPLX(cos(Pi/2),sin(Pi/2))
c       plunging breaker
        A=2.4*DCMPLX(cos(Pi/3),sin(Pi/3))
        B=2.4*DCMPLX(cos(Pi/3),sin(Pi/3))
        A3=DCMPLX(0,0.6)
        A4=DCMPLX(0,-0.6)
        A5=DCMPLX(-0.35,0)

	do i=1,101
	ar(i)=1
        o(i)=-2.2+(i-1)/(25.0)
	enddo
	do i=1,101
	P0(i) = t*ar(i)
	P1(i) = t*o(i)-sqrt(kk)*(1/(2*t)-1/(6*t**3))*II*ar(i)
	P2(i) = t*o(i)**2-kk*(1/(12*t**3)-7/(90*t**5)+1/(84*t**7))*ar(i)-
     &          sqrt(kk)*(1/t-1/(3*t**3))*II*o(i)

	Q0(i) = 1*ar(i)
	Q1(i) = o(i)-sqrt(kk)/2.*(1/(3*t**2)-1/(5*t**4))*II*ar(i)
	Q2(i) = o(i)**2-kk*(1/(60*t**4)-13/(630*t**6)+1/(180*t**8))*ar(i)-
     &          sqrt(kk)*(1/(3*t**2)-1/(5*t**4))*II*o(i)
	Q3(i) = o(i)**3-1.5*sqrt(kk)*(1/(3*t**2)-1/(5*t**4))*II*o(i)**2+
     &          3*kk*(1/(60*t**4)-13/(630*t**6)+1/(180*t**8))*o(i)+
     &          II*kk**1.5*(1/(840*t**6)-17/(7560*t**8))*ar(i)
	Q4(i) = o(i)**4-2*sqrt(kk)*(1/(3*t**2)-1/(5*t**4))*II*o(i)**3+
     &          6*kk*(1/(60*t**4)-13/(630*t**6)+1/(180*t**8))*o(i)**2+
     &          4*II*kk**1.5*(1/(840*t**6)-17/(7560*t**8))*o(i)+
     &          2*kk**2/(30240*t**8)*ar(i) 
	Q5(i) = o(i)**5-2.5*sqrt(kk)*(1/(3*t**2)-1/(5*t**4))*II*o(i)**4+
     &          10*kk*(1/(60*t**4)-13/(630*t**6)+1/(180*t**8))*o(i)**3+
     &          10*II*kk**1.5*(1/(840*t**6)-17/(7560*t**8))*o(i)**2+
     &          10*kk**2/(30240*t**8)*o(i)

	zed(i+49) = A*(P0(i)+P2(i)/2.)*II-B*(Q0(i)+Q2(i)/2.)+A3*Q3(i)+A4*Q4(i)+
     &          A5*Q5(i)-0.03*II*Q5(i)+8.1*t*ar(i)-t**1.5*II*ar(i)

	enddo


	tang = zed(50)-zed(51)
	tang = t*DREAL(tang)+II*DIMAG(tang)
	do i = 1,49
	zed(50-i) = zed(50-i+1)+tang
	tang = DREAL(tang)+II*DIMAG(tang)/i**(0.5)
	enddo


c       aa and bb are as in aa*exp(-t^2/c)+bb
	d=zed(150)-zed(50)
	aa=DIMAG(d)
	bb=DIMAG(zed(50))
        do i=1,49
	zed(i+150) = DREAL(zed(150))-i+II*(aa*exp(-REAL(i)**2/800)+bb)
        enddo
	return 
	end


      subroutine FORT_MASSEXPECT(mass,time)
      IMPLICIT NONE
      REAL_T mass,time,massadd
      REAL_T ymax,L,dropout,tsubtract

#include "probdata.H"

      if (probtype.eq.14) then
       ymax=vinletgas
       L=ymax-yblob
       dropout=two*Pi*radblob*radblob*radblob/three - Pi*L*radblob**2 +
     &    Pi*L*L*L/three
       mass=four*Pi*radblob*radblob*radblob/three - dropout
      else if (probtype.eq.16) then

       if (advbot.eq.zero) then
        mass=zblob
       else
        massadd=zero
        tsubtract=two*zblob*zblob*zblob/
     &               (three*radblob*radblob*advbot)
        if (zblob.lt.zero) then
         tsubtract=zero
         massadd=500.0
        endif
        mass=advbot*(tsubtract+time)*Pi*radblob**2+massadd
       endif
      else
       print *,"exact mass not known for this case"
      endif

      return
      end



c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the
c ::: distance from the interface,curv and vort. .  If fn
c ::: inbetween 0 and 1 then we tag as steep gradient.  If value less
c ::: than 1, then we tag as being inside or near the bubble.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: tag      <=  INTEGER_T tag array
c ::: lo,hi     => index extent of tag array
c ::: set       => INTEGER_T value to tag cell for refinement
c ::: clear     => INTEGER_T value to untag cell
c ::: vfrac     => derived array with flag for whethor near interface ...
c ::: ng        => number of ghost zones in phi array (should be 1)
c ::: nvar      => number of components in phi array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: level     => we look at this level in order to determine
c                   where to adapt at level+1
c ::: -----------------------------------------------------------

c 0 - not at interface
c 1 - at interface
c 2 - at interface and curvature "high"
c 3 - at interface and vorticity "high"
c 4 - at interface and curvature and vorticity "high"
      subroutine FORT_VFRACERROR (tag,DIMS(tag),set,clear,
     &      vfrac,DIMS(vfrac),xpos,DIMS(xpos),lo,hi,nvar,
     &      domlo,domhi,dx,xlo,problo,time,level,
     &      start_curvature_level,nblocks,xblocks,yblocks,zblocks,
     &      rxblocks,ryblocks,rzblocks,ncoarseblocks,
     &      xcoarseblocks,ycoarseblocks,zcoarseblocks,
     &      rxcoarseblocks,rycoarseblocks,rzcoarseblocks)
      IMPLICIT NONE

      INTEGER_T nblocks,ncoarseblocks,start_curvature_level
      REAL_T xblocks(10),yblocks(10),zblocks(10)
      REAL_T rxblocks(10),ryblocks(10),rzblocks(10)
      REAL_T xcoarseblocks(10),ycoarseblocks(10),zcoarseblocks(10)
      REAL_T rxcoarseblocks(10),rycoarseblocks(10),rzcoarseblocks(10)

      INTEGER_T   DIMDEC(tag)
      INTEGER_T   DIMDEC(vfrac)
      INTEGER_T   DIMDEC(xpos)
      INTEGER_T   nvar, set, clear, level
      INTEGER_T   lo(SDIM), hi(SDIM)
      INTEGER_T   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      INTEGER_T   tag(DIMV(tag))
      REAL_T    vfrac(DIMV(vfrac),nvar)
      REAL_T    xpos(DIMV(xpos),SDIM)

      REAL_T    x, y, ax, ay, aerr, rflag
      REAL_T    xdiff, ydiff, wallbot,walltop
      REAL_T    dist
      REAL_T    hx,hy,hz
      INTEGER_T   i, j,np
      INTEGER_T   tagflag

#include "probdata.H"

      if ((nblocks.lt.0).or.(ncoarseblocks.lt.0).or.
     &    (nblocks.ge.10).or.(ncoarseblocks.ge.10)) then
       print *,"nblocks or ncoarseblocks out of range"
       stop
      endif
      if (start_curvature_level<1) then
       print *,"start_curvature_level invalid"
       stop
      endif

      call checkbound(lo,hi,DIMS(xpos),1,-1,1400)
      call checkbound(lo,hi,DIMS(tag),0,-1,1400)
      call checkbound(lo,hi,DIMS(vfrac),0,-1,1400)

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
       y = xpos(i,j,2)
       x = xpos(i,j,1)
       hx=half*(xpos(i+1,j,1)-xpos(i-1,j,1))
       hy=half*(xpos(i,j+1,2)-xpos(i,j-1,2))
       hz=hy

       rflag=vfrac(i,j,1)
       tagflag=0

       if ((probtype.eq.25).and.(axis_dir.gt.0)) then
        if ( (abs(x-xblob).le.radblob).and.(y.le.zblob) ) then
         tagflag=1
        endif
        if ( (abs(x-xblob).le.radblob).and.(zblob.eq.zero).and.
     &       (y.le.hy) ) then
         tagflag=1
        endif

        if (rflag.gt.zero) then
         tagflag=1
        endif
       endif

       if ((level+1.ge.start_curvature_level).and.(tcenter.gt.zero)) then
        if ((rflag.eq.two).or.(rflag.eq.four)) then
         tagflag=1
        endif
       else
        if (rflag.gt.zero) then
         tagflag=1
        endif
       endif

       if (probtype.eq.53) then
        if ((x.gt.xblob+40.0*radblob).or.
     &      (y.gt.50.0*radblob)) then
         tagflag=0
        endif
        if ( (abs(x-xblob).le.two*hx).and.(y.le.two*hy)) then
         tagflag=1
        endif
       endif
      
       if (ractivex.gt.zero) then
        if ((abs(x-xactive).gt.ractivex).or.
     &      (abs(y-yactive).gt.ractivey)) then
         tagflag=0
        endif
       endif

       if (radblob10.gt.zero) then
       if ((probtype.ne.1).or.(axis_dir.ne.13)) then
        if ((abs(x-xblob10).le.radblob10).and.
     &      (abs(y-yblob10).le.radblob10)) then
         tagflag=1
        endif
       endif
       endif

       if (nblocks.gt.0) then
        do np=1,nblocks
         if ((abs(x-xblocks(np)).le.rxblocks(np)).and.
     &       (abs(y-yblocks(np)).le.ryblocks(np))) then
          tagflag=1
         endif
        enddo
       endif

       if (ncoarseblocks.gt.0) then
        do np=1,ncoarseblocks
         if ((abs(x-xcoarseblocks(np)).ge.rxcoarseblocks(np)).or.
     &       (abs(y-ycoarseblocks(np)).ge.rycoarseblocks(np))) then
          tagflag=0
         endif
        enddo
       endif

       if (tagflag.eq.1) then
        tag(i,j)=set
       endif
      enddo
      enddo

      return
      end



c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: rho      <=  density array
c ::: DIMS(rho) => index extent of rho array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------


c used for levelset functions
      subroutine FORT_DENFILL (rho,DIMS(rho),domlo,domhi,dx,
     &                         xlo,time,bc )
      IMPLICIT NONE

      INTEGER_T    DIMDEC(rho)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      INTEGER_T    bc(SDIM,2)
      INTEGER_T use_vfrac

      use_vfrac=10
      call denfillmain(rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc,use_vfrac)
      return
      end 


      subroutine FORT_VAPORFILL (rho,DIMS(rho),domlo,domhi,dx,
     &                         xlo,time,bc )
      IMPLICIT NONE

      INTEGER_T    DIMDEC(rho)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      INTEGER_T    bc(SDIM,2)
      INTEGER_T use_vfrac

      use_vfrac=12
      call denfillmain(rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc,use_vfrac)
      return
      end 

      subroutine FORT_VOFVAPFILL (rho,DIMS(rho),domlo,domhi,dx,
     &                         xlo,time,bc )
      IMPLICIT NONE

      INTEGER_T    DIMDEC(rho)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      INTEGER_T    bc(SDIM,2)
      INTEGER_T use_vfrac

      use_vfrac=101
      call denfillmain(rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc,use_vfrac)
      return
      end 

      subroutine FORT_SOLIDFILL (rho,DIMS(rho),domlo,domhi,dx,
     &                         xlo,time,bc )
      IMPLICIT NONE
      
      INTEGER_T    DIMDEC(rho)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      INTEGER_T    bc(SDIM,2)
      INTEGER_T use_vfrac

      use_vfrac=11
      call denfillmain(rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc,use_vfrac)
      return
      end

      subroutine FORT_VOFFILL (rho,DIMS(rho),domlo,domhi,dx,
     &                         xlo,time,bc )
      IMPLICIT NONE

      INTEGER_T    DIMDEC(rho)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      INTEGER_T    bc(SDIM,2)
      INTEGER_T use_vfrac

      use_vfrac=100
      call denfillmain(rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc,use_vfrac)
      return
      end 

c used for rho(phi)
      subroutine FORT_DENFILL3 (rho,DIMS(rho),domlo,domhi,dx,
     &                         xlo,time,bc )
      IMPLICIT NONE

      INTEGER_T    DIMDEC(rho)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      INTEGER_T    bc(SDIM,2)
      INTEGER_T use_vfrac

      use_vfrac=3
      call denfillmain(rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc,use_vfrac)
      return
      end 

      subroutine get_orifice_velocity(x,vel)
      IMPLICIT NONE
      REAL_T x,vel
        
#include "probdata.H"
          
      if ((x.ge.(xblob-radblob)).and.(x.le.(xblob+radblob))) then
       vel=advbot
      else
       vel=zero
      endif
                        
      return
      end

c use_vfrac=10 levelset
c use_vfrac=11 solid levelset
c use_vfrac=12 vapor levelset
c use_vfrac=100 VOF
c use_vfrac=101 vapor VOF
c use_vfrac=3  rho(phi) 

      subroutine denfillmain (rho,DIMS(rho),domlo,domhi,dx,
     &                         xlo,time,bc,use_vfrac )
      IMPLICIT NONE

      INTEGER_T    DIMDEC(rho),use_vfrac
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j,flagarray(2),flagend

#include "probdata.H"

      INTEGER_T lo(SDIM),hi(SDIM)
      REAL_T hx,hy,templevel,rhoy
      REAL_T swall,rad

      flagend=1

      hx=dx(1)
      hy=dx(2)

      lo(1)=ARG_L1(rho)
      lo(2)=ARG_L2(rho)
      hi(1)=ARG_H1(rho)
      hi(2)=ARG_H2(rho)
      if (lo(1).lt.domlo(1)) then
       lo(1)=domlo(1)
      endif
      if (hi(1).gt.domhi(1)) then
       hi(1)=domhi(1)
      endif
      if (lo(2).lt.domlo(2)) then
       lo(2)=domlo(2)
      endif
      if (hi(2).gt.domhi(2)) then
       hi(2)=domhi(2)
      endif

      call filcc(rho,DIMS(rho),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(rho).lt.domlo(1)) then
         do i = ARG_L1(rho), domlo(1)-1
         do j = ARG_L2(rho), ARG_H2(rho)
          rho(i,j)=rho(domlo(1),j)
	 enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(rho).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(rho)
         do j = ARG_L2(rho), ARG_H2(rho)
          rho(i,j)=rho(domhi(1),j)
	 enddo
	 enddo
      endif            


      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(rho).lt.domlo(2)) then
           do j = ARG_L2(rho), domlo(2)-1
           do i = ARG_L1(rho), ARG_H1(rho)
            rho(i,j)=rho(i,domlo(2))
	   enddo
	   enddo
      endif    

      if ( (bc(2,2).eq.EXT_DIR).and.(ARG_H2(rho).gt.domhi(2)) ) then
         do j = domhi(2)+1, ARG_H2(rho)
         do i = ARG_L1(rho), ARG_H1(rho)
          rho(i,j)=rho(i,domhi(2))
         enddo
         enddo
      endif
   
      return
      end


      subroutine FORT_TEMPFILL (adv,DIMS(adv),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(adv)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j

#include "probdata.H"

      call filcc(adv,DIMS(adv),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(adv).lt.domlo(1)) then
         do i = ARG_L1(adv), domlo(1)-1
            do j = ARG_L2(adv), ARG_H2(adv)
             adv(i,j) = twater
             if (probtype.eq.82) then
              adv(i,j) = tsatdef
             endif
	    enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(adv).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(adv)
            do j = ARG_L2(adv), ARG_H2(adv)
	     adv(i,j) = twater
             if (probtype.eq.82) then
              adv(i,j) = walltemp
             endif
	    enddo
	 enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(adv).lt.domlo(2)) then

         do j = ARG_L2(adv), domlo(2)-1
            do i = ARG_L1(adv), ARG_H1(adv)
             if ((probtype.eq.55).and.(axis_dir.eq.2)) then
              adv(i,j) = twater
             else
              adv(i,j) = walltemp
             endif
	    enddo
	 enddo

      endif            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do i = ARG_L1(adv), ARG_H1(adv)
	       adv(i,j) = twater
	    enddo
	 enddo
      endif            

      return
      end


      subroutine FORT_VTEMPFILL (adv,DIMS(adv),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(adv)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j

#include "probdata.H"

      call filcc(adv,DIMS(adv),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(adv).lt.domlo(1)) then
         do i = ARG_L1(adv), domlo(1)-1
            do j = ARG_L2(adv), ARG_H2(adv)
             adv(i,j) = tvapor
	    enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(adv).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(adv)
            do j = ARG_L2(adv), ARG_H2(adv)
	       adv(i,j) = tvapor
	    enddo
	 enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(adv).lt.domlo(2)) then

         do j = ARG_L2(adv), domlo(2)-1
            do i = ARG_L1(adv), ARG_H1(adv)
             adv(i,j) = walltemp
             if ((probtype.eq.55).and.(axis_dir.eq.0)) then
              adv(i,j)=tvapor
             endif
	    enddo
	 enddo

      endif            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do i = ARG_L1(adv), ARG_H1(adv)
	       adv(i,j) = tvapor
               if ((probtype.eq.55).and.
     &             ((axis_dir.eq.4).or.(axis_dir.eq.5))) then
                adv(i,j)=twater
               endif
	    enddo
	 enddo
      endif            

      return
      end


c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: adv      <=  advected quantity array
c ::: DIMS(adv) => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of adv array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(adv)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j

#include "probdata.H"

      call filcc(adv,DIMS(adv),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(adv).lt.domlo(1)) then
         do i = ARG_L1(adv), domlo(1)-1
            do j = ARG_L2(adv), ARG_H2(adv)
	       adv(i,j) = zero
	    enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(adv).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(adv)
            do j = ARG_L2(adv), ARG_H2(adv)
	       adv(i,j) = zero
	    enddo
	 enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(adv).lt.domlo(2)) then

         do j = ARG_L2(adv), domlo(2)-1
            do i = ARG_L1(adv), ARG_H1(adv)
	       adv(i,j) = one
	    enddo
	 enddo

      endif            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do i = ARG_L1(adv), ARG_H1(adv)
	       adv(i,j) = one
	    enddo
	 enddo
      endif            

      return
      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: u        <=  x velocity array
c ::: DIMS(u)   => index extent of u array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------


      subroutine FORT_TENSORFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(u)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j

#include "probdata.H"

      call filcc(u,DIMS(u),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(u).lt.domlo(1)) then
         do i = ARG_L1(u), domlo(1)-1
            do j = ARG_L2(u), ARG_H2(u)
	       u(i,j) = zero
	    enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(u).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(u)
            do j = ARG_L2(u), ARG_H2(u)
	       u(i,j) = zero
	    enddo
	 enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(u).lt.domlo(2)) then
         do j = ARG_L2(u), domlo(2)-1
            do i = ARG_L1(u), ARG_H1(u)
	       u(i,j) = zero
	    enddo
	 enddo
      endif            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(u).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(u)
            do i = ARG_L1(u), ARG_H1(u)
	       u(i,j) = zero
	    enddo
	 enddo
      endif            

      return
      end


      subroutine FORT_XPARTFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(u)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j

#include "probdata.H"

      call filcc(u,DIMS(u),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(u).lt.domlo(1)) then
         do i = ARG_L1(u), domlo(1)-1
            do j = ARG_L2(u), ARG_H2(u)
	       u(i,j) = zero
	    enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(u).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(u)
            do j = ARG_L2(u), ARG_H2(u)
	       u(i,j) = zero
	    enddo
	 enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(u).lt.domlo(2)) then
         do j = ARG_L2(u), domlo(2)-1
            do i = ARG_L1(u), ARG_H1(u)
	       u(i,j) = zero
	    enddo
	 enddo
      endif            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(u).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(u)
            do i = ARG_L1(u), ARG_H1(u)
	       u(i,j) = zero
	    enddo
	 enddo
      endif            

      return
      end



      subroutine FORT_XVELFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(u)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j,error
      REAL_T     x_vel

#include "probdata.H"

      if (adv_dir .eq. 1)then
         call rampvel(time,x_vel)
      else  
         x_vel = zero
      endif

      call filcc(u,DIMS(u),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(u).lt.domlo(1)) then
         do i = ARG_L1(u), domlo(1)-1
         do j = ARG_L2(u), ARG_H2(u)
	  u(i,j) = x_vel
	 enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(u).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(u)
         do j = ARG_L2(u), ARG_H2(u)
	  u(i,j) = x_vel
	 enddo
	 enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(u).lt.domlo(2)) then
         do j = ARG_L2(u), domlo(2)-1
         do i = ARG_L1(u), ARG_H1(u)
	  u(i,j) = x_vel
	 enddo
	 enddo
      endif            

      
      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(u).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(u)
         do i = ARG_L1(u), ARG_H1(u)
	  u(i,j) = x_vel
	 enddo
	 enddo
      endif            

      return
      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: v        <=  y velocity array
c ::: DIMS(v)  => index extent of v array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_YVELFILL (v,DIMS(v),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(v)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     v(DIMV(v))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j
      REAL_T     y_vel,hx,hy

#include "probdata.H"

      hx=dx(1)
      hy=dx(2)

      if (adv_dir .eq. 2) then
         y_vel = adv_vel
      else  
         y_vel = zero
      endif
      if (probtype.eq.53) then
       y_vel=zero
      endif

      call filcc(v,DIMS(v),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(v).lt.domlo(1)) then
         do i = ARG_L1(v), domlo(1)-1
         do j = ARG_L2(v),ARG_H2(v)
	  v(i,j) = y_vel
	 enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(v).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(v)
         do j = ARG_L2(v), ARG_H2(v)
          v(i,j) = y_vel
	 enddo
	 enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(v).lt.domlo(2)) then
         do j = ARG_L2(v), domlo(2)-1
         do i = ARG_L1(v), ARG_H1(v)
          v(i,j) = y_vel
	 enddo
	 enddo
      endif            


      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(v).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(v)
         do i = ARG_L1(v), ARG_H1(v)
	  v(i,j) = y_vel
         enddo
	 enddo
      endif            

      return
      end


c xlo=physical location of LL corner of u array
c xlo=problo+dx*(ARG_L-domlo)

      subroutine FORT_XPOSFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(u)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j
      REAL_T     xtest,ytest

#include "probdata.H"

      call filcc(u,DIMS(u),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(u).lt.domlo(1)) then
         do i = ARG_L1(u), domlo(1)-1
         do j = ARG_L2(u), ARG_H2(u)
          xtest=xlo(1)+dx(1)*(i-ARG_L1(u)+half)
          ytest=xlo(2)+dx(2)*(j-ARG_L2(u)+half)
	  u(i,j) = xtest
	 enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(u).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(u)
         do j = ARG_L2(u), ARG_H2(u)
          xtest=xlo(1)+dx(1)*(i-ARG_L1(u)+half)
          ytest=xlo(2)+dx(2)*(j-ARG_L2(u)+half)
	  u(i,j) = xtest
	 enddo
	 enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(u).lt.domlo(2)) then
         do j = ARG_L2(u), domlo(2)-1
         do i = ARG_L1(u), ARG_H1(u)
          xtest=xlo(1)+dx(1)*(i-ARG_L1(u)+half)
          ytest=xlo(2)+dx(2)*(j-ARG_L2(u)+half)
	  u(i,j) = xtest
	 enddo
	 enddo
      endif            

      
      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(u).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(u)
         do i = ARG_L1(u), ARG_H1(u)
          xtest=xlo(1)+dx(1)*(i-ARG_L1(u)+half)
          ytest=xlo(2)+dx(2)*(j-ARG_L2(u)+half)
	  u(i,j) = xtest
	 enddo
	 enddo
      endif            

      return
      end

      subroutine FORT_YPOSFILL (v,DIMS(v),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(v)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     v(DIMV(v))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j
      REAL_T  xtest,ytest

#include "probdata.H"

      call filcc(v,DIMS(v),domlo,domhi,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(v).lt.domlo(1)) then
         do i = ARG_L1(v), domlo(1)-1
         do j = ARG_L2(v),ARG_H2(v)
          xtest=xlo(1)+dx(1)*(i-ARG_L1(v)+half)
          ytest=xlo(2)+dx(2)*(j-ARG_L2(v)+half)
	  v(i,j) = ytest
	 enddo
	 enddo
      endif            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(v).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(v)
         do j = ARG_L2(v), ARG_H2(v)
          xtest=xlo(1)+dx(1)*(i-ARG_L1(v)+half)
          ytest=xlo(2)+dx(2)*(j-ARG_L2(v)+half)
	  v(i,j) = ytest
	 enddo
	 enddo
      endif            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(v).lt.domlo(2)) then
         do j = ARG_L2(v), domlo(2)-1
         do i = ARG_L1(v), ARG_H1(v)
          xtest=xlo(1)+dx(1)*(i-ARG_L1(v)+half)
          ytest=xlo(2)+dx(2)*(j-ARG_L2(v)+half)
	  v(i,j) = ytest
	 enddo
	 enddo
      endif            


      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(v).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(v)
         do i = ARG_L1(v), ARG_H1(v)
          xtest=xlo(1)+dx(1)*(i-ARG_L1(v)+half)
          ytest=xlo(2)+dx(2)*(j-ARG_L2(v)+half)
	  v(i,j) = ytest
         enddo
	 enddo
      endif            

      return
      end

      subroutine FORT_UMACFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(u)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      INTEGER_T    bc(SDIM,2)

#include "probdata.H"

      print *,"this routine obsolete"
      stop

      return
      end

      subroutine FORT_VMACFILL (v,DIMS(v),domlo,domhi,dx,xlo,time,bc)
      IMPLICIT NONE

      INTEGER_T    DIMDEC(v)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     v(DIMV(v))
      INTEGER_T    bc(SDIM,2)


#include "probdata.H"

      print *,"this routine obsolete"
      stop

      return
      end



c as opposed to the other fill routines, the exterior homogeneous
c dirichlet bc lives
c at the center of the GHOST CELL and not ON THE EDGE.

      subroutine FORT_CENPRESFILL (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)

      INTEGER_T    DIMDEC(p)
      INTEGER_T    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     p(DIMV(p))
      INTEGER_T    bc(SDIM,2)

      INTEGER_T    i, j
      INTEGER_T    ilo, ihi, jlo, jhi
      logical    fix_xlo, fix_xhi, fix_ylo, fix_yhi

      call filcc(p,DIMS(p),domlo,domhi,bc)

      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq.EXT_DIR)
      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq.EXT_DIR)
      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq.EXT_DIR)
      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq.EXT_DIR)

      ilo = max(ARG_L1(p),domlo(1))
      ihi = min(ARG_H1(p),domhi(1))
      jlo = max(ARG_L2(p),domlo(2))
      jhi = min(ARG_H2(p),domhi(2))

      if (fix_xlo) then
       do i = ARG_L1(p), domlo(1)-1
       do j = ARG_L2(p),ARG_H2(p)
        p(i,j) = zero
       enddo
       enddo
      endif
      if (fix_xhi) then
       do i = domhi(1)+1, ARG_H1(p)
       do j = ARG_L2(p),ARG_H2(p)
	p(i,j) = zero
       enddo
       enddo
      endif

      if (fix_ylo) then
       do j = ARG_L2(p), domlo(2)-1
       do i = ARG_L1(p), ARG_H1(p)
        p(i,j) = zero
       enddo
       enddo
      endif
      if (fix_yhi) then
       do j = domhi(2)+1, ARG_H2(p)
       do i = ARG_L1(p), ARG_H1(p)
        p(i,j) = zero
       enddo
       enddo
      endif

      return
      end

      subroutine FORT_FORCEVELOCITY(problo,probhi,
     &  vel,DIMS(vel),dir,xpos,DIMS(xpos),lo,hi,time)
      IMPLICIT NONE
      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T DIMDEC(vel)
      INTEGER_T DIMDEC(xpos)
      INTEGER_T dir
      REAL_T time
      REAL_T problo(SDIM),probhi(SDIM)
      REAL_T xpos(DIMV(xpos),SDIM)
      REAL_T vel(DIMV(vel))

#include "probdata.H"

      INTEGER_T i,j,k,ii,jj,kk
      REAL_T y

      call checkbound(lo,hi,DIMS(xpos),1,-1,2400)
      call checkbound(lo,hi,DIMS(vel),0,dir,2400)

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir out of range in FORCEVELOCITY"
       stop
      endif
      do i=lo(1),hi(1)+ii
      do j=lo(2),hi(2)+jj
       if(probtype.eq.72) then
        y=half*(xpos(i,j,2)+xpos(i,j-jj,2))
        if(y>yblob2) then
         if(dir.eq.0) then
          vel(i,j)=zero
         elseif (dir.eq.1) then
          vel(i,j)=adv_vel
         endif
        endif
       endif
c standard pipe problem
       if ((probtype.eq.41).and.(axis_dir.eq.1)) then
        y=half*(xpos(i,j,SDIM)+xpos(i,j-jj,SDIM))
        if (yblob2.gt.zero) then
         if (y.ge.probhi(SDIM)-yblob2) then
          if (dir.eq.0) then
           vel(i,j)=zero
          else if (dir.eq.1) then
           if (vel(i,j).lt.zero) then
            vel(i,j)=zero
           endif
          else
           print *,"dir invalid"
           stop
          endif
         endif
        endif
       endif 


      enddo
      enddo

      return
      end

