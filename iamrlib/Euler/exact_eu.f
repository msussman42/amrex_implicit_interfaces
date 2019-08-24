      subroutine exact(time,minCell,maxCell,xNode,uExCell)
      implicit none

      include 'initial.h'
      include 'wave_constants.h'
      include 'grid.h'
      include 'exact.h'
c
c     Input:
      real*8 time
      integer minCell,maxCell
      real*8 xNode(I_FISRT_GHOST_NODE:I_LAST_GHOST_NODE)
c     Output:
      real*8 uExCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL,1:N_EQUATIONS)
c     Auxiliary:
      real*8 dxCurr
      integer iDisc,iState,iVar
      integer nStateCurr
      real*8 uCurr(0:N_INIT_STATE,1:N_EQUATIONS), xCurr(0:N_INIT_DISC)

      dxCurr = dxInit
      nStateCurr = nStateInit
      do iState = 1,nStateInit
         do iVar=1,N_EQUATIONS
            uCurr(iState,iVar) = cInit(iState,iVar)
         enddo
      enddo
      do iDisc = 1,nStateInit+1
         xCurr(iDisc) = xInit(iDisc) + (time-timeInit)*wInit(iDisc)
      enddo

      call initialStates(time,
     &     nStateCurr,uCurr,xCurr,dxCurr,xNode,minCell,maxCell,uExCell)

      end
