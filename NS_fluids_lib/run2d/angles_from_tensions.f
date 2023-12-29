      program main
      IMPLICIT NONE

      real*8 sigma12,sigma13,sigma23
      real*8 mypi
      real*8 s1,s2,s3,s3_new,asin3_1,asin3_2
      real*8 s3_1,s3_2
      real*8 angle_test(0:3)
      real*8 s3_new_list(0:3)
      real*8 e_test(0:3)
      real*8 e_list(0:8)
      real*8 theta_list(0:8,0:2)
      real*8 ds
      integer irange_crit,irange,irange_start
      real*8 s3_range(0:2000)
      real*8 e_range(0:2000)
      integer icrit,itest
      integer i1,j1,k1
      integer select

      select=6

      if (select.eq.1) then
       sigma12=0.831d0   ! liquid/gas
       sigma23=0.8341d0  ! gas/ice
       sigma13=0.17407d0 ! liquid/ice
      else if (select.eq.2) then
       sigma12=sqrt(3.0d0)/2.0d0   ! liquid/gas
       sigma23=sigma12  ! gas/ice
       sigma13=sigma12  ! liquid/ice
      else if (select.eq.3) then
       sigma12=4.0d0/9.0d0
       sigma13=5.0d0/9.0d0
       sigma23=4.0d0/9.0d0
      else if (select.eq.4) then
       sigma12=0.831d0   ! liquid/gas
       sigma23=0.831d0  ! gas/ice
       sigma13=0.49977d0 ! liquid/ice
      else if (select.eq.5) then
       sigma12=0.072d0   ! liquid/gas
       sigma23=0.072d0  ! gas/ice
       sigma13=0.036d0 ! liquid/ice
      else if (select.eq.6) then
       sigma12=72.8d0   ! liquid/gas
       sigma23=29.1d0  ! gas/ice
       sigma13=29.1d0 ! liquid/ice
      else
       print *,"select invalid"
       stop
      endif

      mypi=4.d0*atan(1.d0)

      ds=1.0d0/1000.0d0

      irange_start=10
      irange_crit=irange_start
      do irange=irange_start,1000
       s3_range(irange)=irange*ds
       s3=s3_range(irange)

       s3_1=sigma23*s3/sigma12
       asin3_1=dasin(s3_1)
       s3_2=sigma13*s3/sigma12
       asin3_2=dasin(s3_2)

       angle_test(0)=2.0d0*mypi-asin3_1-asin3_2
       angle_test(1)=2.0d0*mypi-(mypi-asin3_1)-asin3_2
       angle_test(2)=2.0d0*mypi-(mypi-asin3_2)-asin3_1
       angle_test(3)=2.0d0*mypi-(mypi-asin3_2)-(mypi-asin3_1)

       do itest=0,3
        s3_new_list(itest)=sin(angle_test(itest))
        if ((angle_test(itest).gt.0.0d0).and.
     &      (angle_test(itest).lt.mypi)) then
         e_test(itest)=abs(s3_new_list(itest)-s3)
        else
         e_test(itest)=10.0d0
        endif
       enddo

       icrit=0
       do itest=1,3
        if (e_test(itest).le.e_test(icrit)) then
         icrit=itest
        endif
       enddo
       if ((s3_1.gt.1.0d0).or.(s3_2.gt.1.0d0)) then
        e_test(icrit)=10.0d0
       endif
       e_range(irange)=e_test(icrit)
       if (e_range(irange).lt.e_range(irange_crit)) then
        irange_crit=irange
       endif
      enddo
 
      s3=s3_range(irange_crit)
      s1=sigma23*s3/sigma12
      s2=sigma13*s3/sigma12
      print *,"irange_start=",irange_start
      print *,"irange_crit=",irange_crit
      print *,"s1,s2,s3 ",s1,s2,s3
      print *,"prior to enforcing 2 pi constraint:"
      theta_list(0,0)=dasin(s1)*180.0d0/mypi
      theta_list(0,1)=dasin(s2)*180.0d0/mypi
      theta_list(0,2)=dasin(s3)*180.0d0/mypi
      print *,"theta 1 (degrees) ",theta_list(0,0)
      print *,"theta 2 (degrees) ",theta_list(0,1)
      print *,"theta 3 (degrees) ",theta_list(0,2)
      print *,"after enforcing 2 pi constraint:"
      itest=1
      icrit=1
      do i1=0,1
      do j1=0,1
      do k1=0,1
       theta_list(itest,0)=theta_list(0,0) 
       theta_list(itest,1)=theta_list(0,1) 
       theta_list(itest,2)=theta_list(0,2) 
       if (i1.eq.1) then
        theta_list(itest,0)=180.0d0-theta_list(0,0) 
       endif
       if (j1.eq.1) then
        theta_list(itest,1)=180.0d0-theta_list(0,1) 
       endif
       if (k1.eq.1) then
        theta_list(itest,2)=180.0d0-theta_list(0,2) 
       endif
       e_list(itest)=abs(360.0d0-
     &  theta_list(itest,0)-
     &  theta_list(itest,1)-
     &  theta_list(itest,2))
       if (e_list(itest).lt.e_list(icrit)) then
        icrit=itest
       endif
       itest=itest+1 
      enddo
      enddo
      enddo
      print *,"theta 1 (degrees) ",theta_list(icrit,0)
      print *,"theta 2 (degrees) ",theta_list(icrit,1)
      print *,"theta 3 (degrees) ",theta_list(icrit,2)
      end

