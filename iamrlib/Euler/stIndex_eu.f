      integer function stIndex(i,dir)
      implicit none

c     Computes indexs of states from "inside out", of xt diagram, 
c     in direction dir.
c
c     dir  wave direction parameter
c         = -1 for left  part of the solution to the RP
c         = +1 for right part of the solution to the RP
c
c     For a *N_STATES* where N_STATES is EVEN. 
c
c     i=1,N_STATES/2
c
c------------------------------------------------------------------
c
c     Example for a *4* state system (N_STATES=4):
c
c     dir = -1 ( count states from right to left ):
c        i       stIndex
c       ---      -------
c        1          2  middle state (behind) LEFT
c        2          1  initial state (ahead) LEFT
c
c     dir = +1 ( count states from left to right ): 
c        i       stIndex
c       ---      -------
c        1          3   middle state (behind) RIGHT       
c        2          4   initial state (ahead) RIGHT

      include 'wave_constants.h'

c     Input.
      integer i,dir
c-----

      stIndex = (N_STATES+1-dir)/2 + dir*i

      end
