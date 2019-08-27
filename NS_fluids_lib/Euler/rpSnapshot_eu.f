      subroutine rpSnapshot(
     &     rho_0,vel_0,pre_0,snd_0,w_0,signif_0,wave_type_0,IsVacuum,
     &     xL,xM,xR,time)
      implicit none
      
      include 'wave_constants.h'

c     Input:
      real*8 rho_0(N_STATES),vel_0(N_STATES),pre_0(N_STATES)
      real*8 snd_0(N_STATES)
      real*8 w_0(N_EDGES)
      integer signif_0(N_WAVES),wave_type_0(N_WAVES)
      logical IsVacuum
      real*8 xL,xM,xR,time
c     Auxiliary:
      integer iWave,dir,iLedge,iRedge,iLstate,iRstate
      real*8 xLedge,xRedge


      end
