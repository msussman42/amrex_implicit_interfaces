      subroutine setPrimExt(minCellG,maxCellG,cCell,uCell)
      implicit none

      include 'wave_constants.h'
      include 'grid.h'

c     Input:
      integer minCellG,maxCellG
      real*8 cCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL,1:N_EQUATIONS)
c     Output: 
      real*8 uCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL,1:N_EXT_VARS)
c     Auxiliary:
      integer iCell
      real*8 thermSt2snd

      do iCell = minCellG,maxCellG

         call stCons2stPrim(
     &        cCell(iCell,IRHO),cCell(iCell,IMOM),cCell(iCell,IENE),
     &        uCell(iCell,IRHO),uCell(iCell,IVEL),uCell(iCell,IPRE)) 

         uCell(iCell,ISND) = thermSt2snd(
     &        uCell(iCell,IRHO),uCell(iCell,IPRE))  

      enddo

      end
