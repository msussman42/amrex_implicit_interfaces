      PROGRAM time_average

      IMPLICIT NONE

      integer, PARAMETER :: NPLOT=1000
      Real*8, PARAMETER :: xstart=0.0d0
      Real*8, PARAMETER :: xend=1.2d0
      Real*8 :: xdata(0:NPLOT) 
      Real*8 :: ydata(0:NPLOT) 
      integer :: i

      print *,"xstart=",xstart
      print *,"xend=",xend
      print *,"NPLOT=",NPLOT
      print *,"redirect the output for gnuplot use"

      do i=0,NPLOT
       xdata(i)=xstart+i*(xend-xstart)/NPLOT
       ydata(i)=138.3d0*exp(-1.41d0*xdata(i))
       print *,xdata(i),"  ",ydata(i)
      enddo

      END PROGRAM
