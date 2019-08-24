
      subroutine get_bottom_elevation(x,elevation)

      IMPLICIT NONE

      real*8 x,elevation
      real*8 h,l

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
      if (abs(x/l).le.1.0) then
       elevation=h*(1.0-2.0*((x/l)**2)+(x/l)**4)
      else if (x/l.lt.0.0) then
       elevation=0.0
      else if (x/l.gt.0.0) then
       elevation=0.0
      else
       print *,"bust"
       stop
      endif

      return
      end 

      subroutine get_left_velocityIOWA(t,u_left)

      IMPLICIT NONE
      integer k1,k2
      real*8 t,u_left

      real*8 inflow_time(3000)
      real*8 inflow_elevation(3000)
      real*8 inflow_velocity(3000)
      real*8 outflow_time(3000)
      real*8 outflow_elevation(3000)
      real*8 outflow_velocity(3000)
      integer inflow_count,outflow_count
      integer last_inflow_index,last_outflow_index
      common /fred_stern_data/ inflow_time,inflow_elevation, &
         inflow_velocity,outflow_time,outflow_elevation, &
         outflow_velocity,inflow_count,outflow_count, &
         last_inflow_index,last_outflow_index

      if (inflow_count.le.0) then
       print *,"inflow_count invalid"
       stop
      endif
     
      do while ((inflow_time(last_inflow_index).gt.t).and. &
                (last_inflow_index.gt.1)) 
       last_inflow_index=last_inflow_index-1
      enddo
      do while ((inflow_time(last_inflow_index).lt.t).and. &
                (last_inflow_index.lt.inflow_count))
       last_inflow_index=last_inflow_index+1
      enddo
      if (t.le.inflow_time(1)) then
       u_left=inflow_velocity(1)
      else if (t.ge.inflow_time(inflow_count)) then
       u_left=inflow_velocity(inflow_count)
      else 

       if (t.le.inflow_time(last_inflow_index)) then
        k1=last_inflow_index-1
       else
        k1=last_inflow_index
       endif
       k2=k1+1
       u_left=inflow_velocity(k1)+ &
        ( (inflow_velocity(k2)-inflow_velocity(k1))/ &
          (inflow_time(k2)-inflow_time(k1)) )*(t-inflow_time(k1))
      endif
      u_left=u_left*100.0  ! convert to cm/s

      return
      end 


      subroutine get_right_velocityIOWA(t,u_right)

      IMPLICIT NONE
      integer k1,k2
      real*8 t,u_right

      real*8 inflow_time(3000)
      real*8 inflow_elevation(3000)
      real*8 inflow_velocity(3000)
      real*8 outflow_time(3000)
      real*8 outflow_elevation(3000)
      real*8 outflow_velocity(3000)
      integer inflow_count,outflow_count
      integer last_inflow_index,last_outflow_index
      common /fred_stern_data/ inflow_time,inflow_elevation, &
         inflow_velocity,outflow_time,outflow_elevation, &
         outflow_velocity,inflow_count,outflow_count, &
         last_inflow_index,last_outflow_index

      if (outflow_count.le.0) then
       print *,"outflow_count invalid"
       stop
      endif
     
      do while ((outflow_time(last_outflow_index).gt.t).and. &
                (last_outflow_index.gt.1)) 
       last_outflow_index=last_outflow_index-1
      enddo
      do while ((outflow_time(last_outflow_index).lt.t).and. &
                (last_outflow_index.lt.outflow_count))
       last_outflow_index=last_outflow_index+1
      enddo
      if (t.le.outflow_time(1)) then
       u_right=outflow_velocity(1)
      else if (t.ge.outflow_time(outflow_count)) then
       u_right=outflow_velocity(outflow_count)
      else 

       if (t.le.outflow_time(last_outflow_index)) then
        k1=last_outflow_index-1
       else
        k1=last_outflow_index
       endif
       k2=k1+1
       u_right=outflow_velocity(k1)+ &
        ( (outflow_velocity(k2)-outflow_velocity(k1))/ &
          (outflow_time(k2)-outflow_time(k1)) )*(t-outflow_time(k1))
      endif
      u_right=u_right*100.0  ! convert to cm/s

      return
      end 


      subroutine get_left_elevationIOWA(t,elevation_left)

      IMPLICIT NONE
      integer k1,k2
      real*8 t,elevation_left

      real*8 inflow_time(3000)
      real*8 inflow_elevation(3000)
      real*8 inflow_velocity(3000)
      real*8 outflow_time(3000)
      real*8 outflow_elevation(3000)
      real*8 outflow_velocity(3000)
      integer inflow_count,outflow_count
      integer last_inflow_index,last_outflow_index
      common /fred_stern_data/ inflow_time,inflow_elevation, &
         inflow_velocity,outflow_time,outflow_elevation, &
         outflow_velocity,inflow_count,outflow_count, &
         last_inflow_index,last_outflow_index


      if (inflow_count.le.0) then
       print *,"inflow_count invalid"
       stop
      endif
     
      do while ((inflow_time(last_inflow_index).gt.t).and. &
                (last_inflow_index.gt.1)) 
       last_inflow_index=last_inflow_index-1
      enddo
      do while ((inflow_time(last_inflow_index).lt.t).and. &
                (last_inflow_index.lt.inflow_count))
       last_inflow_index=last_inflow_index+1
      enddo
      if (t.le.inflow_time(1)) then
       elevation_left=inflow_elevation(1)
      else if (t.ge.inflow_time(inflow_count)) then
       elevation_left=inflow_elevation(inflow_count)
      else 

       if (t.le.inflow_time(last_inflow_index)) then
        k1=last_inflow_index-1
       else
        k1=last_inflow_index
       endif
       k2=k1+1
       elevation_left=inflow_elevation(k1)+ &
        ( (inflow_elevation(k2)-inflow_elevation(k1))/ &
          (inflow_time(k2)-inflow_time(k1)) )*(t-inflow_time(k1))
      endif
      elevation_left=elevation_left*100.0  ! convert to cm

      return
      end 


      subroutine get_right_elevationIOWA(t,elevation_right)

      IMPLICIT NONE
      integer k1,k2
      real*8 t,elevation_right

      real*8 inflow_time(3000)
      real*8 inflow_elevation(3000)
      real*8 inflow_velocity(3000)
      real*8 outflow_time(3000)
      real*8 outflow_elevation(3000)
      real*8 outflow_velocity(3000)
      integer inflow_count,outflow_count
      integer last_inflow_index,last_outflow_index
      common /fred_stern_data/ inflow_time,inflow_elevation, &
         inflow_velocity,outflow_time,outflow_elevation, &
         outflow_velocity,inflow_count,outflow_count, &
         last_inflow_index,last_outflow_index


      if (outflow_count.le.0) then
       print *,"outflow_count invalid"
       stop
      endif
     
      do while ((outflow_time(last_outflow_index).gt.t).and. &
                (last_outflow_index.gt.1)) 
       last_outflow_index=last_outflow_index-1
      enddo
      do while ((outflow_time(last_outflow_index).lt.t).and. &
                (last_outflow_index.lt.outflow_count))
       last_outflow_index=last_outflow_index+1
      enddo
      if (t.le.outflow_time(1)) then
       elevation_right=outflow_elevation(1)
      else if (t.ge.outflow_time(outflow_count)) then
       elevation_right=outflow_elevation(outflow_count)
      else 

       if (t.le.outflow_time(last_outflow_index)) then
        k1=last_outflow_index-1
       else
        k1=last_outflow_index
       endif
       k2=k1+1
       elevation_right=outflow_elevation(k1)+ &
        ( (outflow_elevation(k2)-outflow_elevation(k1))/ &
          (outflow_time(k2)-outflow_time(k1)) )*(t-outflow_time(k1))
      endif
      elevation_right=elevation_right*100.0  ! convert to cm

      return
      end 

      subroutine init_data()

      character*12 namestr1
      character*13 namestr2
      integer i

      real*8 inflow_time(3000)
      real*8 inflow_elevation(3000)
      real*8 inflow_velocity(3000)
      real*8 outflow_time(3000)
      real*8 outflow_elevation(3000)
      real*8 outflow_velocity(3000)
      integer inflow_count,outflow_count
      integer last_inflow_index,last_outflow_index
      common /fred_stern_data/ inflow_time,inflow_elevation, &
         inflow_velocity,outflow_time,outflow_elevation, &
         outflow_velocity,inflow_count,outflow_count, &
         last_inflow_index,last_outflow_index

      print *,"opening InflowBC.dat and OutflowBC.dat"
      namestr1='InflowBC.dat' 
      namestr2='OutflowBC.dat' 

      open(unit=11,file=namestr1)
      read(11,*) inflow_count
      print *,"inflow_count= ",inflow_count
      do i=1,inflow_count
       read(11,*) inflow_time(i),inflow_velocity(i), &
         inflow_elevation(i)
      enddo
      close(11)

      open(unit=12,file=namestr2)
      read(12,*) outflow_count
      print *,"outflow_count= ",outflow_count
      do i=1,outflow_count
       read(12,*) outflow_time(i),outflow_velocity(i), &
         outflow_elevation(i)
      enddo
      close(12)
 
      return
      end


      subroutine doit(problo,probhi,ncell,dx,tstop)
      IMPLICIT NONE

      real*8 problo,probhi,dx,tstop
      integer ncell
      character*4 filename
      real*8 SS(-1:ncell)
      real*8 QQ(-1:ncell)
      real*8 DD(-1:ncell)  ! bottom topography
      real*8 xx(-1:ncell)
      real*8 SSflux(0:ncell)
      real*8 QQflux(0:ncell)
      real*8 start_elevation
      integer skip,nstep,i
      real*8 gravity,time
      real*8 elevation_right,elevation_left,u_right,u_left
      real*8 maxu,maxc,den,mom,uu,cc,dt,lambda
      real*8 denleft,denright,momleft,momright,pleft,pright
      real*8 minz,slope,maxslope
      integer icrit,icritslope

      skip=2000
      filename='wave' 
      print *,"output filename ",filename
      open(unit=11,file=filename)

      gravity=980.0
      start_elevation=22.862

      time=0.0
      nstep=0

      do i=-1,ncell
       xx(i)=problo+(i+0.5)*dx
       call get_bottom_elevation(xx(i),DD(i))
       QQ(i)=0.0
       SS(i)=start_elevation-DD(i)
      enddo

      do while (time.le.tstop-1.0E-10) 

       call get_right_elevationIOWA(time,elevation_right)
       call get_left_elevationIOWA(time,elevation_left)
       call get_right_velocityIOWA(time,u_right)
       call get_left_velocityIOWA(time,u_left)
       SS(-1)=elevation_left
       SS(ncell)=elevation_right
       QQ(-1)=elevation_left*u_left
       QQ(ncell)=elevation_right*u_right
 
       maxu=0.0
       maxc=0.0 
       do i=0,ncell-1
        den=SS(i)
        mom=QQ(i)
        if (den.le.0.0) then
         print *,"density must be positive"
         stop
        endif
        cc=sqrt(gravity*den)
        if (cc.gt.maxc) then
         maxc=cc
        endif
        uu=abs(mom/den)
        if (uu.gt.maxu) then
         maxu=abs(uu)
        endif
       enddo
       if (maxu+maxc.le.0.0) then
        print *,"must have acoustic waves"
        stop
       endif
       dt=0.8*dx/(maxu+maxc)
       if (time+dt.ge.tstop) then
        dt=tstop-time
       endif
       lambda=dt/dx

       do i=0,ncell
        denleft=SS(i-1)
        denright=SS(i)
        momleft=QQ(i-1)
        momright=QQ(i)
        pleft=0.5*gravity*(denleft**2)
        pright=0.5*gravity*(denright**2)
        SSflux(i)=0.5*(momleft+momright)-0.5*(denright-denleft)/lambda
        if (denright.le.0.0) then
         print *,"hydraulic section must be positive"
         stop
        endif
        if (denleft.le.0.0) then
         print *,"hydraulic section must be positive"
         stop
        endif
        QQflux(i)=0.5*(momleft**2/denleft+momright**2/denright+ &
          pleft+pright)-0.5*(momright-momleft)/lambda
       enddo
       do i=0,ncell-1
        SS(i)=SS(i)-lambda*(SSflux(i+1)-SSflux(i))
        QQ(i)=QQ(i)-lambda*(QQflux(i+1)-QQflux(i))- &
          lambda*SS(i)*gravity*0.5*(DD(i+1)-DD(i-1))
       enddo

       time=time+dt
       nstep=nstep+1

       if ((nstep-skip*(nstep/skip).eq.0).or. &
           (time.ge.tstop-1.0E-10)) then
        print *,"time=",time
        print *,"dt=",dt
        write(11,*) "  "
        icrit=0
        minz=SS(0)+DD(0)
        icritslope=0
        maxslope=0.0
        do i=0,ncell-1
         if (time.ge.tstop-1.0E-10) then
          write(11,*),time,xx(i),SS(i),SS(i)+DD(i),QQ(i)/SS(i)
         endif
         if (SS(i)+DD(i).lt.minz) then
          minz=SS(i)+DD(i)
          icrit=i
         endif
         if ((xx(i).gt.-100.0).and.(xx(i).lt.100.0)) then
          slope=abs((SS(i+1)+DD(i+1)-SS(i)-DD(i))/dx)
          if (slope.gt.maxslope) then
           maxslope=slope
           icritslope=i
          endif
         endif
        enddo
        print *,"icrit,x,minz ",icrit,xx(icrit),minz
        print *,"icritslope,x,maxslope ", &
         icritslope,xx(icritslope),maxslope
       endif
      enddo ! while time<tstop
 
      close(11)

      return
      end
 
      program main
      IMPLICIT NONE

      real*8 problo,probhi,tstop,dx
      integer ncell

      problo=-594.36
      probhi=502.92

       ! 0.1 mm  (0.2mm in Fred Stern's paper)
       ! 100000
      ncell=5000
      dx=(probhi-problo)/ncell
      tstop=11.0

      call init_data()
      call doit(problo,probhi,ncell,dx,tstop)

      return
      end

