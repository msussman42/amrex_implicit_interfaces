      subroutine outputPhasePoint(minCell,maxCell,uCell)
      implicit none

      include 'parameters.h'
      include 'wave_constants.h'
      include 'grid.h'

c     Input:
      integer minCell,maxCell
      real*8 uCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL,1:N_EXT_VARS)
c     Auxiliary:
      integer iCell

      end
