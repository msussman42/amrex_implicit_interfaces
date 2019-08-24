      subroutine checkDt(minCell,maxCell,uCell,dxCell,cfl,dt)
      implicit none

      include 'parameters.h'
      include 'wave_constants.h'
      include 'grid.h'

c     Input:
      integer minCell,maxCell
      real*8 uCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL,1:N_EXT_VARS)
      real*8 dxCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL)
      real*8 cfl
      real*8 dt
c     Auxiliary:
      integer iCell
      real*8 dtCell,maxTime

      do iCell = minCell,maxCell

         dtCell = maxTime
     &        ( uCell(iCell,IVEL), uCell(iCell,ISND), dxCell(iCell) )

         if (dtCell .lt. dt) then
            write(*,*) 'cell #',iCell
            write(*,*) 'dtCell =',dtCell, '   dt =',dt
            stop 'dtCell .lt. dt'
         endif
      enddo

      end
