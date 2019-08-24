      subroutine outputSnapshot
     &     (minCell,maxCell,uCell,xCell)
      implicit none

      include 'parameters.h'
      include 'wave_constants.h'
      include 'grid.h'

c     Input:
      integer minCell,maxCell
      real*8 uCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL,1:N_EXT_VARS)
      real*8 xCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL)
c     Auxiliary:
      integer iCell

      end
