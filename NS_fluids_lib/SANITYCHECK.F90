#undef BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"

module CISL_SANITY_VARS_MODULE
IMPLICIT NONE

REAL_T :: time=zero
REAL_T :: dx=zero
INTEGER_T :: lo,hi
INTEGER_T :: nstep=0
! velocity,density,temperature,pressureMG,pressureEOS
! e=cv T=T
! p=(gamma-1)rho T
REAL_T, dimension(:,:), allocatable :: conserve
REAL_T, dimension(:,:), allocatable :: oldstate
REAL_T, dimension(:,:), allocatable :: newstate
REAL_T, dimension(:), allocatable :: velface
REAL_T, dimension(:), allocatable :: ustarface
REAL_T, dimension(:), allocatable :: denface
REAL_T, dimension(:), allocatable :: pface
REAL_T, dimension(:,:), allocatable :: massfaceside
REAL_T, dimension(:), allocatable :: lower
REAL_T, dimension(:), allocatable :: upper
REAL_T, dimension(:), allocatable :: diag
REAL_T, dimension(:), allocatable :: rhs
REAL_T, dimension(:), allocatable :: soln
REAL_T :: gamma=1.4
REAL_T :: omega=0.4
REAL_T :: problen=one
REAL_T :: DENL=one
REAL_T :: PRESL=one
REAL_T :: DENR=0.125
REAL_T :: PRESR=0.1
INTEGER_T :: nsolve

contains

end module CISL_SANITY_VARS_MODULE

module CISL_SANITY_MODULE
IMPLICIT NONE

contains

      subroutine compare_sanity(comparestate,scomp,ncomp,id)
      use CISL_SANITY_VARS_MODULE

      IMPLICIT NONE

      INTEGER_T scomp,ncomp,i,j,id
      REAL_T err,localerr
      REAL_T comparestate(lo-1:hi+1,5)

      err=zero

      call FLUSH(6)  ! unit=6 screen

      if (id.eq.1) then
       do i=lo,hi
        do j=scomp,scomp+ncomp-1
         localerr=abs(comparestate(i,j)-newstate(i,j))
         if (1.eq.0) then
          print *,"i,j,stateCLSMOF,stateSANITY,err ",i,j, &
           comparestate(i,j),newstate(i,j),localerr
         endif
         if (localerr.gt.err) then
          err=localerr
         endif
        enddo
       enddo
      else if (id.eq.7) then
       do i=lo,hi
        do j=scomp,scomp+ncomp-1
         localerr=abs(comparestate(i,j)-newstate(i,j))
         if (1.eq.0) then
          print *,"i,j,stateCLSMOF,stateSANITY,err ",i,j, &
           comparestate(i,j),newstate(i,j),localerr
         endif
         if (localerr.gt.err) then
          err=localerr
         endif
        enddo
       enddo
      else if (id.eq.8) then
       do i=lo,hi
        do j=scomp,scomp+ncomp-1
         localerr=abs(comparestate(i,j)-newstate(i,j))
         if (1.eq.0) then
          print *,"i,j,stateCLSMOF,stateSANITY,err ",i,j, &
           comparestate(i,j),newstate(i,j),localerr
         endif
         if (localerr.gt.err) then
          err=localerr
         endif
        enddo
       enddo
      else if (id.eq.2) then
       do i=lo,hi
        do j=scomp,scomp+ncomp-1
         localerr=abs(comparestate(i,j)-conserve(i,j))
         if (1.eq.0) then
          print *,"i,j,conserveCLSMOF,conserveSANITY,err ",i,j, &
           comparestate(i,j),conserve(i,j),localerr
         endif
         if (localerr.gt.err) then
          err=localerr
         endif
        enddo
       enddo
      else if (id.eq.3) then
       do i=lo,hi+1
        do j=scomp,scomp+ncomp-1
         localerr=abs(comparestate(i,j)-massfaceside(i,j))
         if (1.eq.0) then
          print *,"i,j,massfaceCLSMOF,massfaceSANITY,err ",i,j, &
           comparestate(i,j),massfaceside(i,j),localerr
         endif
         if (localerr.gt.err) then
          err=localerr
         endif
        enddo
       enddo
      else if (id.eq.4) then
       do i=lo,hi+1
        localerr=abs(comparestate(i,1)-denface(i))
        if (localerr.gt.err) then
         err=localerr
        endif
       enddo
      else if (id.eq.5) then
       do i=lo,hi+1
        localerr=abs(comparestate(i,1)-pface(i))
        if (localerr.gt.err) then
         err=localerr
        endif
       enddo
      else if (id.eq.6) then
       do i=lo,hi+1
        localerr=abs(comparestate(i,1)-velface(i))
        if (localerr.gt.err) then
         err=localerr
        endif
       enddo
      endif

      print *,"COMPARE SANITY ID=",id
      print *,"COMPARE SANITY ERR=",err
      call FLUSH(6)  ! unit=6 screen

      return
      end subroutine compare_sanity

      subroutine init_sanity(left_index,right_index)
      use CISL_SANITY_VARS_MODULE
      IMPLICIT NONE

      INTEGER_T i,j,left_index,right_index
      REAL_T x

      call FLUSH(6)  ! unit=6 screen
      lo=left_index
      hi=right_index
      dx=problen/(hi-lo+1)
      allocate(massfaceside(lo-1:hi+1,2))
      allocate(denface(lo-1:hi+1))
      allocate(oldstate(lo-1:hi+1,5))
      allocate(newstate(lo-1:hi+1,5))
      allocate(conserve(lo-1:hi+1,5))
      allocate(velface(lo-1:hi+1))
      allocate(pface(lo-1:hi+1))
      allocate(ustarface(lo-1:hi+1))
      allocate(lower(lo:hi))
      allocate(upper(lo:hi))
      allocate(diag(lo:hi))
      allocate(rhs(lo:hi))
      allocate(soln(lo:hi))

      do i=lo-1,hi+1
       x=(i+half)*dx
       velface(i)=zero
       if (x.lt.half) then
        oldstate(i,1)=zero
        oldstate(i,2)=DENL
        oldstate(i,3)=PRESL/(omega*DENL)
        oldstate(i,4)=PRESL
        oldstate(i,5)=PRESL
       else if (x.ge.half) then
        oldstate(i,1)=zero
        oldstate(i,2)=DENR
        oldstate(i,3)=PRESR/(omega*DENR)
        oldstate(i,4)=PRESR
        oldstate(i,5)=PRESR
       endif
       do j=1,5
        newstate(i,j)=oldstate(i,j)
       enddo

      enddo

      print *,"init_sanity: lo,hi ",lo,hi
      print *,"init_sanity: dx ",dx
      call FLUSH(6)  ! unit=6 screen

      return
      end subroutine init_sanity

      subroutine sanityBC()
      use CISL_SANITY_VARS_MODULE
      IMPLICIT NONE

      call FLUSH(6)  ! unit=6 screen
      print *,"in sanityBC()"

      oldstate(lo-1,1)=zero
      oldstate(lo-1,2)=DENL
      oldstate(lo-1,3)=PRESL/(omega*DENL)
      oldstate(lo-1,4)=PRESL
      oldstate(lo-1,5)=PRESL

      oldstate(hi+1,1)=zero
      oldstate(hi+1,2)=DENR
      oldstate(hi+1,3)=PRESR/(omega*DENR)
      oldstate(hi+1,4)=PRESR
      oldstate(hi+1,5)=PRESR


      newstate(lo-1,1)=zero
      newstate(lo-1,2)=DENL
      newstate(lo-1,3)=PRESL/(omega*DENL)
      newstate(lo-1,4)=PRESL
      newstate(lo-1,5)=PRESL

      newstate(hi+1,1)=zero
      newstate(hi+1,2)=DENR
      newstate(hi+1,3)=PRESR/(omega*DENR)
      newstate(hi+1,4)=PRESR
      newstate(hi+1,5)=PRESR

      call FLUSH(6)  ! unit=6 screen

      return
      end subroutine sanityBC

      subroutine CISL_sanity(dt,normdir,dir_counter,updateflag)
      use CISL_SANITY_VARS_MODULE
      IMPLICIT NONE

      REAL_T dt

      INTEGER_T i,j,normdir,dir_counter,updateflag
      REAL_T mass_bucket,mom_bucket,energy_bucket
      REAL_T xtargetlo,xtargethi
      REAL_T xdepartlo,xdeparthi
      REAL_T xdonatelo,xdonatehi
      REAL_T xintlo,xinthi

      call FLUSH(6)  ! unit=6 screen
      print *,"CISL_sanity: nstep= ",nstep
      print *,"CISL_sanity: time,dt ",time,dt
      call sanityBC()

      if (dir_counter.eq.0) then
       do i=lo-1,hi+1
        conserve(i,1)=oldstate(i,1)*oldstate(i,2)
        conserve(i,2)=oldstate(i,2)
        conserve(i,3)=(half*oldstate(i,1)**2+oldstate(i,3))*oldstate(i,2)
       enddo
      else if (dir_counter.eq.1) then
       do i=lo-1,hi+1
        conserve(i,1)=newstate(i,1)*newstate(i,2)
        conserve(i,2)=newstate(i,2)
        conserve(i,3)=(half*newstate(i,1)**2+newstate(i,3))*newstate(i,2)
       enddo
      else
       print *,"dir_counter invalid"
       stop
      endif

      if (updateflag.eq.1) then
       ! do nothing
      else if (updateflag.eq.2) then
       do i=lo,hi
        mass_bucket=zero
        mom_bucket=zero
        energy_bucket=zero
        do j=i-1,i+1
         xtargetlo=i*dx
         xtargethi=(i+1)*dx
         if (normdir.eq.0) then
          xdepartlo=xtargetlo-velface(i)*dt
          xdeparthi=xtargethi-velface(i+1)*dt
         else if (normdir.eq.1) then
          xdepartlo=xtargetlo
          xdeparthi=xtargethi
         else
          print *,"normdir invalid"
          stop
         endif
         xdonatelo=j*dx
         xdonatehi=(j+1)*dx
         xintlo=xdonatelo
         xinthi=xdonatehi
         if (xdepartlo.gt.xintlo) then
          xintlo=xdepartlo
         endif
         if (xdeparthi.lt.xinthi) then
          xinthi=xdeparthi
         endif
         if (xinthi.gt.xintlo) then
          mass_bucket=mass_bucket+(xinthi-xintlo)*conserve(j,2)
          mom_bucket=mom_bucket+(xinthi-xintlo)*conserve(j,1)
          energy_bucket=energy_bucket+(xinthi-xintlo)*conserve(j,3)
         endif
        enddo
        newstate(i,1)=mom_bucket/mass_bucket
        newstate(i,2)=mass_bucket/dx
        newstate(i,3)=energy_bucket/mass_bucket-half*(newstate(i,1)**2)  
       enddo
      else
       print *,"updateflag invalid"
       stop
      endif
        
      call FLUSH(6)  ! unit=6 screen
      return
      end subroutine CISL_sanity
      
      subroutine mass_face_weight()
      use CISL_SANITY_VARS_MODULE
      IMPLICIT NONE

      INTEGER_T i

      call FLUSH(6)  ! unit=6 screen
      print *,"SANITY mass_face_weight "

      call sanityBC()
      do i=lo,hi+1
       massfaceside(i,1)=half*dx*newstate(i-1,2)
       massfaceside(i,2)=half*dx*newstate(i,2)
       denface(i)=(massfaceside(i,1)+massfaceside(i,2))/dx
      enddo
      call FLUSH(6)  ! unit=6 screen

      return
      end subroutine mass_face_weight

        ! p/(rho c^2 dt^2)-div(gradp /rho)=peos/(rho c^2 dt^2)-div ustar/dt
        ! c^2=gamma p/rho
        ! p=omega rho e=omega rho T


      subroutine fast_solve(l,u,d,n,f,pressure_soln)
      use CISL_SANITY_VARS_MODULE
      IMPLICIT NONE

      integer n,i
      REAL_T l(n),u(n),d(n),f(n),pressure_soln(n)
      REAL_T ll(n),uu(n),dd(n),z(n)

      dd(1)=d(1)
      uu(1)=u(1)/dd(1)
      z(1)=f(1)/dd(1)
      do i=2,n-1
       ll(i)=l(i)
       dd(i)=d(i)-ll(i)*uu(i-1)
       uu(i)=u(i)/dd(i)
       z(i)=(f(i)-ll(i)*z(i-1))/dd(i)
      enddo
      ll(n)=l(n)
      dd(n)=d(n)-ll(n)*uu(n-1)
      z(n)=(f(n)-ll(n)*z(n-1))/dd(n)
      pressure_soln(n)=z(n)
      do i=n-1,1,-1
       pressure_soln(i)=z(i)-uu(i)*pressure_soln(i+1)
      enddo

      return
      end subroutine fast_solve


      subroutine CISL_projection(dt,sweep)
      use CISL_SANITY_VARS_MODULE
      IMPLICIT NONE

      REAL_T dt,peos,velbefore,velafter,x
      REAL_T c2
      INTEGER_T i,j,sweep

      call FLUSH(6)  ! unit=6 screen
      print *,"CISL_projection: dt,sweep ",dt,sweep
      call sanityBC()
      do i=lo,hi+1
       ustarface(i)=(massfaceside(i,1)*newstate(i-1,1)+ &
        massfaceside(i,2)*newstate(i,1))/ &
        (massfaceside(i,1)+massfaceside(i,2))
      enddo
      ustarface(lo)=zero
      ustarface(hi+1)=zero
      do i=lo,hi
       peos=omega*newstate(i,2)*newstate(i,3)
       c2=gamma*peos/newstate(i,2)
       rhs(i)=peos/(newstate(i,2)*c2*dt*dt)- &
        (ustarface(i+1)-ustarface(i))/(dx*dt)
       lower(i)=-one/(denface(i)*dx*dx)
       upper(i)=-one/(denface(i+1)*dx*dx)
       if (i.eq.lo) then
        lower(i)=zero
       endif
       if (i.eq.hi) then
        upper(i)=zero
       endif
       diag(i)=one/(newstate(i,2)*c2*dt*dt)-lower(i)-upper(i)
      enddo
      nsolve=hi-lo+1
      call fast_solve(lower,upper,diag,nsolve,rhs,soln)

      do i=lo,hi+1
       if (i.eq.lo) then
        pface(i)=soln(lo)
       else if (i.eq.hi+1) then
        pface(i)=soln(hi)
       else
        pface(i)=(massfaceside(i,1)*soln(i)+ &
           massfaceside(i,2)*soln(i-1))/ &
         (massfaceside(i,1)+massfaceside(i,2))
       endif
      enddo
      do i=lo,hi+1
       if (i.eq.lo) then
        velface(i)=zero
       else if (i.eq.hi+1) then 
        velface(i)=zero
       else
        velface(i)=ustarface(i)-dt*(soln(i)-soln(i-1))/ &
         (dx*denface(i))
       endif
      enddo
      do i=lo,hi
       velbefore=newstate(i,1)
       newstate(i,1)=newstate(i,1)-dt*(pface(i+1)-pface(i))/ &
         (newstate(i,2)*dx)
       velafter=newstate(i,1)
       newstate(i,3)=newstate(i,3)+half*(velbefore**2)- &
         half*(velafter**2)- &
         dt*(velface(i+1)*pface(i+1)-velface(i)*pface(i))/ &
         (dx*newstate(i,2))
      enddo
      call sanityBC()

      if (sweep.eq.2) then
       time=time+dt
       nstep=nstep+1
       do i=lo,hi
        do j=1,5
         oldstate(i,j)=newstate(i,j)
        enddo
       enddo 
       call sanityBC()
       if (time.ge.0.15-1.0E-3) then
        call FLUSH(6)  ! unit=6 screen
        print *,"SANITY OUTPUT x,u,rho,T"
        do i=lo,hi
         x=(i+half)*dx
         print *,x,newstate(i,1),newstate(i,2),newstate(i,3)
        enddo
        print *,"AFTER SANITY OUTPUT x,u,rho,T"
       endif
      else if (sweep.eq.1) then
       ! do nothing
      else
       print *,"sweep invalid"
       stop
      endif
      call FLUSH(6)  ! unit=6 screen

      return
      end subroutine CISL_projection


end module CISL_SANITY_MODULE

