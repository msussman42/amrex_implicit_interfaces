      subroutine statesO(rhoA,velA,preA,sndA, velM,preM, IsVacuum,dir,
     &     rho_0,vel_0,pre_0,snd_0, w_0,signif_0,wave_type_0)
      implicit none

      include 'wave_constants.h'

c     Input:
      real*8 rhoA,velA,preA,sndA, velM,preM
      logical IsVacuum
      integer dir
c     Output:
      real*8 rho_0(N_STATES),vel_0(N_STATES),pre_0(N_STATES)
      real*8 snd_0(N_STATES)
      real*8 w_0(N_EDGES)
      integer signif_0(N_WAVES),wave_type_0(N_WAVES)    
c     Auxiliary:
      integer isA,isB,iw,ieA,ieB
      integer stIndex, wvIndex
      
c     states:
      isA = stIndex(2,dir)
      isB = stIndex(1,dir)
c     initial (ahead) state:
      call stPrimExt2stPrimExt(rhoA,velA,preA,sndA,
     &     rho_0(isA),vel_0(isA),pre_0(isA),snd_0(isA))
c     middle (behind) state:
      if (IsVacuum) then
         call setVacMidState(rhoA,velA,preA,sndA,dir, 
     &     rho_0(isB),vel_0(isB),pre_0(isB),snd_0(isB))
      else
         call setRegMidState(rhoA,velA,preA,sndA,velM,preM, 
     &     rho_0(isB),vel_0(isB),pre_0(isB),snd_0(isB))
      endif

c     wave 
      iw=wvIndex(1,dir)
      ieA=WAVE_TO_EDGE_A(iw)
      ieB=WAVE_TO_EDGE_B(iw)
      call wave(rho_0(isB),vel_0(isB),pre_0(isB),snd_0(isB),
     &     rhoA,velA,preA,sndA,dir,
     &     w_0(ieB),w_0(ieA),signif_0(iw),wave_type_0(iw))

      end
