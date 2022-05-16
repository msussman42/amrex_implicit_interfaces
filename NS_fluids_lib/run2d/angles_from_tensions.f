      program main
      IMPLICIT NONE

      real*8 sigma12,sigma13,sigma23
      real*8 mypi
      real*8 s1,s2,s3,s3_new,asin3_1,asin3_2
      integer tol_met,iter
      integer select

      select=2

      if (select.eq.1) then
       sigma12=0.831d0   ! liquid/gas
       sigma23=0.1341d0  ! gas/ice
       sigma13=0.17407d0 ! liquid/ice
      else if (select.eq.2) then
       sigma12=sqrt(3.0d0)/2.0d0   ! liquid/gas
       sigma23=sigma12  ! gas/ice
       sigma13=sigma12  ! liquid/ice
      else
       print *,"select invalid"
       stop
      endif

      mypi=4.d0*atan(1.d0)

      s3=sin(120.0d0*mypi/180.0d0)
      tol_met=0
      iter=0
      do while (tol_met.eq.0) 
       asin3_1=asin(sigma23*s3/sigma12)
       asin3_2=asin(sigma13*s3/sigma12)

       s3_new=sin(2.0d0*mypi-asin3_1-asin3_2)
       if (s3_new.le.0.0d0) then
        s3_new=s3/2.0d0
       endif
       if (s3_new.gt.1.0d0) then
        s3_new=1.0d0
       endif
       if (abs(s3_new-s3).lt.1.0d-8) then
        tol_met=1
       endif   
       s3=s3_new
       iter=iter+1
       print *,"iter,s3,asin3_1,asin3_2 ",iter,s3,asin3_1,asin3_2
      enddo
 
      s1=sigma23*s3/sigma12
      s2=sigma13*s3/sigma12
      print *,"s1,s2,s3 ",s1,s2,s3
      print *,"theta 1 (degrees) ",asin(s1)*180.0d0/mypi 
      print *,"theta 2 (degrees) ",asin(s2)*180.0d0/mypi 
      print *,"theta 3 (degrees) ",asin(s3)*180.0d0/mypi 
      end

