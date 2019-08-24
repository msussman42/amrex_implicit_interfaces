      function setDt(minCell,maxCell,uCell,dxCell,cfl)
      implicit none

      include 'parameters.h'
      include 'wave_constants.h'
      include 'grid.h'

c     Input:
      integer minCell,maxCell
      real*8 uCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL,1:N_EXT_VARS)
      real*8 dxCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL)
      real*8 cfl
c     Output:
      real*8 setDt
c     Auxiliary:
      integer iCell
      real*8 dt,dtCell,maxTime

      dt = t_infty
      do iCell = minCell,maxCell

         dtCell = maxTime
     &        (uCell(iCell,IVEL), uCell(iCell,ISND), dxCell(iCell))

         dt = dmin1( dtCell, dt )
      enddo

      setDt = cfl*dt

      end
