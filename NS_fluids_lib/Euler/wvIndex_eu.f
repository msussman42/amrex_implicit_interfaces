      integer function wvIndex(i,dir)
      implicit none

c     Computes indexs of waves from "inside out", of xt diagram, 
c     in direction dir.
c
c     dir  wave direction parameter
c         = -1 for left  part of the solution to the RP
c         = +1 for right part of the solution to the RP
c
c     For a *N_WAVES* where N_WAVES is ODD. 
c 
c     i=1,(N_WAVES-1)/2
c
c------------------------------------------------------------------
c
c     Example for a *3* wave system (N_WAVES=3):
c
c     dir = -1 ( count waves from right to left ):
c        i       wvIndex   
c       ---      -------
c        1          1     left-facing wave (1-wave)
c
c     dir = +1 ( count waves from left to right ): 
c        i       wvIndex
c       ---      -------
c        1          3     right-facing wave (3-wave)

      include 'wave_constants.h'

c     Input.
      integer i,dir
c-----

      wvIndex = (N_WAVES+1)/2+dir*i

      end
