      program main
      IMPLICIT NONE
      integer nplot,i,j,n_axis
      real*8 temp_low
      real*8 temp_high
      real*8 temp_local
      real*8 cl

      temp_low=240.0
      temp_high=290.0 
      nplot=50

      do i=0,nplot
       temp_local=temp_low+(temp_high-temp_low)*i/nplot
       cl=-0.0011d0*(temp_local**3)+0.732d0*(temp_local**2)-163.0d0*temp_local+16582.0d0
       print *,"temp,cl ",temp_local,cl
      enddo

      return
      end

