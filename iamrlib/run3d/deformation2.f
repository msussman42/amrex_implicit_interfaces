      subroutine fsub(f,a,alpha,c,s,maxx,cenx)
      IMPLICIT NONE

      real*8 f,a,alpha,c,s,maxx,cenx
      real*8 xprime2,zprime2,x2
 
      xprime2=a*a/(1.0+s*s*alpha*alpha/((a**4)*c*c))
      zprime2=(1.0-xprime2/(a*a))*alpha*alpha/(a*a)
      x2=(sqrt(xprime2)*c-sqrt(zprime2)*s)**2
      f=sqrt(x2)-(maxx-cenx)

      return
      end

! grep "MAT=0 dir= 0 MIN INT=" run.out > minx
! grep "MAT=0 dir= 0 MAX INT=" run.out > maxx
! grep "MAT=0 dir= 1 MIN INT=" run.out > miny
! grep "MAT=0 dir= 1 MAX INT=" run.out > maxy
! grep "MAT=0 dir= 2 MIN INT=" run.out > minz
! grep "MAT=0 dir= 2 MAX INT=" run.out > maxz

      program main
      IMPLICIT NONE

      integer ntime,i,iter,j
      real*8 radblob,PI
      real*8 time,time2,minx,maxx,miny,maxy,minz,maxz
      real*8 cenx,ceny,cenz,V,angle,b,alpha,c,s
      real*8 alow,ahigh,flow,fhigh,amid,fmid,cc
      real*8 deformation,minerr,testerr
      real*8 zprime2,xprime2,z2

      print *,"this program reads 'minx','maxx','miny','maxy'"
      print *,"'minz','maxz'"
      print *,"and outputs 'deformation'"

      ntime=5090
      radblob=0.5
      print *,"ntime is the number of entries in the file"
      print *,"ntime=",ntime
      print *,"radblob is the initial radius of the sphere"
      print *,"radblob=",radblob
      if (radblob.le.0.0) then
       print *,"radblob invalid"
       stop
      endif
 
      open(unit=12, file='minx')
      open(unit=13, file='maxx')
      open(unit=14, file='miny')
      open(unit=15, file='maxy')
      open(unit=16, file='minz')
      open(unit=17, file='maxz')
      open(unit=18, file='deformation')

      PI=4.0*atan(1.0)

        ! x'=x cos(angle)+z sin(angle)
        ! z'=-x sin(angle) + z cos(angle)
        ! V=(4/3) pi abc
        ! c=alpha/a
        ! x'^2/a^2 + y'^2/b^2 + z'^2/c^2=1
      do i=1,ntime
       read(12,*) time,minx
       read(13,*) time2,maxx
       read(14,*) time2,miny
       read(15,*) time2,maxy
       read(16,*) time2,minz
       read(17,*) time2,maxz
       if (abs(time-time2).gt.1.0E-10) then
        print *,"time discrepancy"
        stop
       endif
       V=4.0*PI*(radblob**3)/3.0
       cenx=0.5*(minx+maxx)
       ceny=0.5*(miny+maxy)
       cenz=0.5*(minz+maxz)
       if (maxx-cenx.lt.1.0E-10) then
        print *,"maxx-minx too small"
        stop
       endif
       b=maxy-ceny
       alpha=V*3.0/(4.0*PI*b)


!      angle=atan((maxz-cenz)/(maxx-cenx))
       minerr=1.0E+10
       do j=1,99
        angle=j*0.5*PI/100.0
 
        c=cos(angle)
        s=sin(angle)

        alow=0.5*radblob
        ahigh=10.0*radblob
        call fsub(flow,alow,alpha,c,s,maxx,cenx)
        call fsub(fhigh,ahigh,alpha,c,s,maxx,cenx)
        if (flow*fhigh.le.0.0) then

         do iter=1,30
          if (flow*fhigh.gt.0.0) then
           print *,"flow,fhigh bust"
           stop
          endif
          amid=(alow+ahigh)/2.0
          call fsub(fmid,amid,alpha,c,s,maxx,cenx)
          if (flow*fmid.le.0.0) then
           ahigh=amid
           fhigh=fmid
          else
           alow=amid
           flow=fmid
          endif
         enddo 
         cc=alpha/amid

         zprime2=alpha*alpha/amid**2
         zprime2=zprime2/(1.0+(amid**4)*s*s/(alpha*alpha*c*c))
         xprime2=amid*amid*(1.0-amid*amid*zprime2/(alpha*alpha))
         z2=(sqrt(xprime2)*s+sqrt(zprime2)*c)**2
         testerr=abs( sqrt(z2)-(maxz-cenz))

         print *,"angle,err ",angle,testerr
 
         if (testerr.lt.minerr) then
          deformation=(amid-cc)/(amid+cc)       
          minerr=testerr
         endif

        endif

       enddo
       write(18,*) time,deformation
      enddo

      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)

      end

