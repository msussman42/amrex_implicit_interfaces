#define NCELL 128
#define PROBLO (0.d0)
#define PROBHI (2.d0)
#define probtype (0)
#define heatvisc (1.0e+5)
#define heatviscQ (1.0e+5)
#define cv (4.1855e+7)
#define den (1.0d0)
#define DTLEGACY (0.03125)
#define STOP_TIME (5.0)
#define maxstep (9999)
#define Pi (3.14159265358979d0)
#define basetemp (293.0)
#define walltemp (400.0)
#define clamped (0)
#define implicitsolve (0)

      subroutine init_temp(x,dx,temp)
      IMPLICIT NONE

      real*8 x,dx,temp

      temp=basetemp

      return
      end 

      subroutine temp_bc(xlo,dx,lo,hi,ngrow,temp)
      IMPLICIT NONE

      integer lo,hi,ngrow
      real*8 xlo,dx           
      real*8 temp(lo-ngrow:hi+ngrow)
      integer i,nper

      nper=hi-lo+1
      do i=lo-1,lo-ngrow,-1
       temp(i)=walltemp
      enddo
      do i=hi+1,hi+ngrow
       temp(i)=basetemp
      enddo

      return
      end

      subroutine output_temp(xlo,dx,lo,hi,temp)
      IMPLICIT NONE

      integer lo,hi
      real*8 xlo,dx
      real*8 temp(lo-1:hi+1)
      real*8 tempfab(lo-1:hi+1)
      integer i,isten,isub
      real*8 x_output,den_output,denprim
      character*4 filename

      filename='temp'
      print *,"filename ",filename
      open(unit=11,file=filename)

      do i=lo,hi
       tempfab(i)=temp(i)
      enddo
      call temp_bc(xlo,dx,lo,hi,1,tempfab)
      do i=lo,hi
       x_output=xlo+dx*(i-lo)+dx*0.5
       write(11,*) x_output,tempfab(i)
      enddo

      close(11)
300   FORMAT(I2)
400   FORMAT(A4,A2)

      return
      end
     
      subroutine advect(xlo,dx,lo,hi,temp,tempnew,dt)
      IMPLICIT NONE

      real*8 xlo,dx,dt
      integer lo,hi
      real*8 temp(lo-1:hi+1)
      real*8 tempnew(lo-1:hi+1)
      real*8 tempfab(lo-1:hi+1)
      real*8 rhs(lo-1:hi+1)
      real*8 diag(lo-1:hi+1)
      real*8 offdiag(lo-1:hi+1)
      integer i,igrid,isten,errormet,iter
      real*8 fluxplus,fluxminus,Q,error0,errorn
      

      do i=lo,hi
       tempfab(i)=temp(i)
      enddo
      call temp_bc(xlo,dx,lo,hi,1,tempfab)
      Q=-heatviscQ*(walltemp-basetemp)

      if (implicitsolve.eq.0) then

       do i=lo,hi
        fluxplus=heatvisc*(tempfab(i+1)-tempfab(i))/dx
        fluxminus=heatvisc*(tempfab(i)-tempfab(i-1))/dx
        if ((i.eq.lo).and.(clamped.eq.0)) then
         fluxminus=Q
        endif
       
        tempnew(i)=temp(i)+dt*(fluxplus-fluxminus)/ &
         (dx*den*cv)
       enddo  ! i

      else if (implicitsolve.eq.1) then
       do i=lo,hi
        rhs(i)=temp(i)
        diag(i)=1.0+2.0*dt*heatvisc/(dx*dx*den*cv)
        offdiag(i)=-dt*heatvisc/(dx*dx*den*cv)
        if (i.eq.lo) then
         rhs(i)=rhs(i)-dt*Q/(dx*den*cv)
         diag(i)=1.0+dt*heatvisc/(dx*dx*den*cv)
        endif
        if (i.eq.hi) then
         rhs(i)=rhs(i)+dt*heatvisc*basetemp/(dx*dx*den*cv)
        endif
        tempfab(i)=temp(i) 
       enddo
       errormet=0
       iter=0
       do while (errormet.eq.0) 
        do i=lo,hi
         if (i.eq.lo) then
          tempnew(i)=(rhs(i)-tempfab(i+1)*offdiag(i))/diag(i)
         else if (i.eq.hi) then
          tempnew(i)=(rhs(i)-tempfab(i-1)*offdiag(i))/diag(i)
         else
          tempnew(i)=(rhs(i)-tempfab(i+1)*offdiag(i)- &
             tempfab(i-1)*offdiag(i))/diag(i)
         endif
        enddo
        errorn=0.0
        do i=lo,hi
         errorn=errorn+(tempnew(i)-tempfab(i))**2
         tempfab(i)=tempnew(i)
        enddo
        errorn=sqrt(errorn)
        if (iter.eq.0) then
         error0=errorn
        else if (error0.eq.0.0) then
         errormet=1
        else if (errorn.le.1.0E-10*error0) then
         errormet=1
        endif
        print *,"error0,errorn,iter ",error0,errorn,iter
        iter=iter+1
       enddo
       
        
      else
       print *,"option invalid"
       stop
      endif

      return
      end

      program main
      IMPLICIT NONE
      real*8 dx,xlo,xhi
      integer lo,hi
      real*8 temp(-1:NCELL)
      real*8 tempnew(-1:NCELL)
      real*8 time,dt,x
      integer i,nstep


      lo=0
      hi=NCELL-1
      xlo=PROBLO
      dx=(PROBHI-PROBLO)/NCELL

      do i=lo,hi
       x=xlo+(i-lo+0.5d0)*dx
       call init_temp(x,dx,temp(i))
      enddo

      dt=DTLEGACY

      nstep=0
      do while ((time.le.STOP_TIME-1.0E-10).and. &
                (nstep.lt.maxstep))
       do i=lo,hi
        tempnew(i)=temp(i)
       enddo
       if (time+dt.ge.STOP_TIME) then
        dt=STOP_TIME-time
       endif
       call advect(xlo,dx,lo,hi,temp,tempnew,dt)
       do i=lo,hi
        temp(i)=tempnew(i)
       enddo
       time=time+dt
       nstep=nstep+1
      enddo

      call output_temp(xlo,dx,lo,hi,temp) 
      
      return
      end

