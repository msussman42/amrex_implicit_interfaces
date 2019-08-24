      subroutine setFluxes(minNode,maxNode,uCell,fNode)
      implicit none

      include 'parameters.h'
      include 'wave_constants.h'
      include 'grid.h'

c     Input:
      integer minNode,maxNode
      real*8 uCell(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL,1:N_EXT_VARS)
c     Input/Output: 
      real*8 fNode(I_FISRT_GHOST_NODE:I_LAST_GHOST_NODE,1:N_EQUATIONS)
c     Auxiliary:
      integer iNode,iCellL,iCellR
      real*8 DxDt
      real*8 rho_0(N_STATES),vel_0(N_STATES),pre_0(N_STATES)
      real*8 snd_0(N_STATES)
      real*8 w_0(N_EDGES)
      integer signif_0(N_WAVES),wave_type_0(N_WAVES)
      logical IsVacuum
      real*8 rho,vel,pre,snd

      do iNode = minNode, maxNode
         iCellL=leftCell(iNode)
         iCellR=rghtCell(iNode)
         call riemann(
     &        uCell(iCellL,IRHO),uCell(iCellL,IVEL),
     &        uCell(iCellL,IPRE),uCell(iCellL,ISND),
     &        uCell(iCellR,IRHO),uCell(iCellR,IVEL),
     &        uCell(iCellR,IPRE),uCell(iCellR,ISND),
     &        rho_0,vel_0,pre_0,snd_0,w_0,signif_0,wave_type_0,IsVacuum)
         DxDt=0.d0
         call slope2state(     
     &        rho_0,vel_0,pre_0,snd_0,w_0,signif_0,wave_type_0, DxDt,
     &        rho,vel,pre,snd)
         call st2flux(rho,vel,pre, DxDt,
     &        fNode(iNode,IRHO),fNode(iNode,IMOM),fNode(iNode,IENE))
      enddo

      end
