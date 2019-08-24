! assume all arrays are dimensioned as:
!   data(is-1:ie+1,js-1:je+1)
! inputs: xpos (positions of cell centers)
!         levelnormal,levelfab (input level set function; levelnormal is
!                               a copy of levelfab)
!         levelflag (initially all zero)
!         maxdist=initially a negative number
!         mindist=initially a positive number
!         lo,hi (dimensions of "active region")
!         is,ie,js,je (dimensions of "active region" again)
!         narrow_band_size 
! outputs: levelfab, newfab (distance function from piecewise linear
!                            interface; newfab is a copy of levelfab)
!          levelflag (equals 1 if a distance was initialized)
!          maxdist,mindist (after exit of this routine, another routine
!                           should set newfab to be maxdist or mindist 
!                           in cells where levelflag=0) 
! assume that all cell centered data has one ghost cell.
      subroutine marching_triangles(newfab, &
        xpos, &
        levelnormal, &
        levelfab, &
        levelflag, &
        lo,hi,maxdist,mindist,narrow_band_size, &
        is,ie,js,je)
      IMPLICIT NONE

      integer is,ie,js,je
      integer narrow_band_size
      integer lo(2),hi(2)
      real*8 maxdist,mindist

      real*8  xpos(is-1:ie+1,js-1:je+1,2)
      real*8  newfab(is-1:ie+1,js-1:je+1)
      real*8  levelfab(is-1:ie+1,js-1:je+1)
      real*8  levelnormal(is-1:ie+1,js-1:je+1)
      real*8  levelflag(is-1:ie+1,js-1:je+1)
      real*8 lsnd(0:1,0:1)
      real*8 xnd(0:1,0:1,2)
      integer itriangle,iline,i1,j1,dir,i,j
      integer ngrow
      real*8 xtri(0:2,2) 
      real*8 lstri(0:2) 
      real*8 x1(2),x2(2),xtarget(2),xclose(2)
      real*8 aa,dist_line,mag,dot,tt,temp_dist,signdist

      ngrow=narrow_band_size
      if (ngrow.lt.3) then
       print *,"ngrow too small"
       stop
      endif

! assume all data has one ghost cell 
      do i=is-1,ie
      do j=js-1,je
       
       do i1=0,1
       do j1=0,1
        lsnd(i1,j1)=levelnormal(i+i1,j+j1)
        if (lsnd(i1,j1).eq.0.0) then
         lsnd(i1,j1)=1.0E-12
        endif
        do dir=1,2
         xnd(i1,j1,dir)=xpos(i+i1,j+j1,dir)
        enddo
       enddo
       enddo

       do itriangle=1,2
        iline=0
        if (itriangle.eq.1) then
         lstri(0)=lsnd(0,0)
         lstri(1)=lsnd(1,0)
         lstri(2)=lsnd(0,1)
         do dir=1,2
          xtri(0,dir)=xnd(0,0,dir)
          xtri(1,dir)=xnd(1,0,dir)
          xtri(2,dir)=xnd(0,1,dir)
         enddo
        else if (itriangle.eq.2) then
         lstri(0)=lsnd(1,1)
         lstri(1)=lsnd(1,0)
         lstri(2)=lsnd(0,1)
         do dir=1,2
          xtri(0,dir)=xnd(1,1,dir)
          xtri(1,dir)=xnd(1,0,dir)
          xtri(2,dir)=xnd(0,1,dir)
         enddo
        else
         print *,"itriangle invalid"
         stop
        endif

        if ((lstri(0)*lstri(1).le.0.0).and. &
            (lstri(1)*lstri(2).le.0.0)) then
          aa=abs(lstri(0))/(abs(lstri(0))+abs(lstri(1))+1.0E-12)
          do dir=1,2
           x1(dir)=(1.0-aa)*xtri(0,dir)+aa*xtri(1,dir)
          enddo
          aa=abs(lstri(1))/(abs(lstri(1))+abs(lstri(2))+1.0E-12)
          do dir=1,2
           x2(dir)=(1.0-aa)*xtri(1,dir)+aa*xtri(2,dir)
          enddo
          iline=1
        else if ((lstri(0)*lstri(2).le.0.0).and. &
                 (lstri(1)*lstri(2).le.0.0)) then
          aa=abs(lstri(0))/(abs(lstri(0))+abs(lstri(2))+1.0E-12)
          do dir=1,2
           x1(dir)=(1.0-aa)*xtri(0,dir)+aa*xtri(2,dir)
          enddo
          aa=abs(lstri(1))/(abs(lstri(1))+abs(lstri(2))+1.0E-12)
          do dir=1,2
           x2(dir)=(1.0-aa)*xtri(1,dir)+aa*xtri(2,dir)
          enddo
          iline=1
        else if ((lstri(0)*lstri(2).le.0.0).and. &
                 (lstri(0)*lstri(1).le.0.0)) then
          aa=abs(lstri(0))/(abs(lstri(0))+abs(lstri(2))+1.0E-12)
          do dir=1,2
           x1(dir)=(1.0-aa)*xtri(0,dir)+aa*xtri(2,dir)
          enddo
          aa=abs(lstri(0))/(abs(lstri(0))+abs(lstri(1))+1.0E-12)
          do dir=1,2
           x2(dir)=(1.0-aa)*xtri(0,dir)+aa*xtri(1,dir)
          enddo
          iline=1
        endif

        if (iline.eq.1) then
         do i1=-ngrow,ngrow
         do j1=-ngrow,ngrow
          if ((i1+i.ge.is).and.(i1+i.le.ie).and. &
              (j1+j.ge.js).and.(j1+j.le.je)) then
           do dir=1,2
            xtarget(dir)=xpos(i+i1,j+j1,dir)
           enddo
           dist_line=-1.0
! check critical points
! f(t)=(x(t)-x)^2 + (y(t)-y)^2
! x(t)=x1+t(x2-x1)
! y(t)=y1+t(y2-y1)
! f'(t)=2(x(t)-x)x'(t)+2(y(t)-y)y'(t)=0
! (x1-x+t(x2-x1))(x2-x1)+(y1-y+t(y2-y1))(y2-y1)=0
! t=( -(x1-x)(x2-x1)-(y1-y)(y2-y1) )/
!   ( (x2-x1)^2 + (y2-y1)^2 )
           mag=(x2(1)-x1(1))**2+(x2(2)-x1(2))**2
           dot=-(x1(1)-xtarget(1))*(x2(1)-x1(1))- &
                (x1(2)-xtarget(2))*(x2(2)-x1(2))
           if (mag.gt.0.0) then
            tt=dot/mag
            if ((tt.ge.0.0).and.(tt.le.1.0)) then
             dist_line=0.0
             do dir=1,2
              xclose(dir)=(1.0-tt)*x1(dir)+tt*x2(dir)
              dist_line=dist_line+(xclose(dir)-xtarget(dir))**2
             enddo
             dist_line=sqrt(dist_line)
            endif
           endif
! now check endpoints
           temp_dist=0.0
           do dir=1,2
            xclose(dir)=x1(dir)
            temp_dist=temp_dist+(xclose(dir)-xtarget(dir))**2
           enddo
           temp_dist=sqrt(temp_dist)
           if ((dist_line.lt.0.0).or.(temp_dist.lt.dist_line)) then
            dist_line=temp_dist
           endif

           temp_dist=0.0
           do dir=1,2
            xclose(dir)=x2(dir)
            temp_dist=temp_dist+(xclose(dir)-xtarget(dir))**2
           enddo
           temp_dist=sqrt(temp_dist)
           if ((dist_line.lt.0.0).or.(temp_dist.lt.dist_line)) then
            dist_line=temp_dist
           endif

           if (levelnormal(i+i1,j+j1).lt.0.0) then
            signdist=-1.0
           else
            signdist=1.0
           endif
           dist_line=dist_line*signdist

           if (levelflag(i+i1,j+j1).eq.0.0) then
            levelfab(i+i1,j+j1)=dist_line
           else
            if (abs(dist_line).lt. &
                abs(levelfab(i+i1,j+j1))) then
             levelfab(i+i1,j+j1)=dist_line
            endif
           endif
           levelflag(i+i1,j+j1)=1.0
          endif ! i+i1,j+j1 in the interior
         enddo
         enddo
        endif ! line was found
       enddo ! itriangle

      enddo
      enddo
 
      do i=is,ie
      do j=js,je
       newfab(i,j)=levelfab(i,j)
       if (levelflag(i,j).ne.0.0) then
        if (levelfab(i,j).lt.mindist) then
         mindist=levelfab(i,j)
        endif
        if (levelfab(i,j).gt.maxdist) then
         maxdist=levelfab(i,j)
        endif
       endif
      enddo
      enddo

      print *,"after marching triangles completed, maxdist,mindist ", &
       maxdist,mindist
         
      return
      end
