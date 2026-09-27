      program main
      IMPLICIT NONE
      integer nplot,i,j,n_axis
      real*8 T_axis(6)
      real*8 mu_axis(6)
      real*8 temp_low
      real*8 temp_high
      real*8 temp_local
      real*8 visc

      n_axis=6
      T_axis(1)=240.0
      T_axis(2)=250.0
      T_axis(3)=260.0
      T_axis(4)=270.0
      T_axis(5)=280.0
      T_axis(6)=290.0
      mu_axis(1)=17.0
      mu_axis(2)=5.0
      mu_axis(3)=3.0
      mu_axis(4)=2.0
      mu_axis(5)=1.5
      mu_axis(6)=1.0
      nplot=100
      temp_low=240.0
      temp_high=290.0 

      do i=0,nplot
       temp_local=temp_low+(temp_high-temp_low)*i/nplot
       if (temp_local.le.T_axis(1)) then
        visc=mu_axis(1)
       else if (temp_local.ge.T_axis(n_axis)) then
        visc=mu_axis(n_axis)
       else
        do j=1,n_axis-1
         if ((temp_local.ge.T_axis(j)).and. &
             (temp_local.le.T_axis(j+1))) then
          visc=mu_axis(j)+(mu_axis(j+1)-mu_axis(j))* &
                (temp_local-T_axis(j))/(T_axis(j+1)-T_axis(j))
         endif
        enddo
       endif
       visc=visc/100.0
       print *,"temp,visc ",temp_local,visc
      enddo

      return
      end

