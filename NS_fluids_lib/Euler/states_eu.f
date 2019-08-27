      subroutine states(
     &     rhoL,velL,preL,sndL, rhoR,velR,preR,sndR,
     &     velM,preM,IsVacuum,
     &     rho_0,vel_0,pre_0,snd_0, w_0,signif_0,wave_type_0)
      implicit none
c
c-----------------------------------------------------------------------
c                              INPUT 
c-----------------------------------------------------------------------
c
c
c       {Left,Middle,Right} STATE VARIABLES:
c       -----------------------------------
c    rho{L,R}  density
c    vel{L,M,R}  velocity
c    pre{L,M,R}  pressure
c    snd{L,R}  speed of sound
c-------------------------------------------------------------------------
c                            OUTPUT 
c-------------------------------------------------------------------------
c
c          STATE VARIABLES
c          ---------------
c    rho_0(1:N_STATES)   density
c    vel_0(1:N_STATES)   velocity
c    pre_0(1:N_STATES)   pressure
c    snd_0(1:N_STATES)   speed of sound
c     1 - Left Init. State;  2 - Left Middle State; 
c                            3 - Right Middle State; 4 - Right Init. State
c
c         WAVE STRUCTURE 
c         --------------
c    1 - left-facing wave; 2 - central wave; 3 - right-facing wave
c
c    wave_type_0(1:N_WAVES) -- type of the wave
c     = RFN rarefaction fan
c     = SFR shock front
c     = CD contact discontinuity
c
c    signif_0(1:N_WAVES) =1 if the wave is signigicant; =0 otherwise
c       THIS INFORMATION IS GIVEN FOR WAVE TRACKING PURPOSES: 
c       signif_0(i) = 1: i-th wave should be tracked
c       signif_0(i) = 0: i-th wave should not be tracked
c
c    EACH WAVE CAN BE EITHER A RAREFACTION FAN OR A SHARP FRONT.  
c    THE FOLLOWING CONVENTION IS USED: 
c    A SHARP FRONT IS DESCRIBED AS A FAN WHOSE EDGES COINCIDE.
c
c    w_0(1:N_EDGES)    velocities of the fan-edges:
c        w_0(1) = velocity of the outer  edge of the left-f fan
c        w_0(2) = velocity of the innner edge of the left-f fan
c        w_0(3) = velocity of the left   edge of the center fan
c        w_0(4) = velocity of the right  edge of the center fan
c        w_0(5) = velocity of the inner  edge of the right-f fan
c        w_0(6) = velocity of the outer  edge of the right-f fan
c
c
      include 'wave_constants.h'

c     Input:
      real*8 rhoL,velL,preL,sndL, rhoR,velR,preR,sndR, velM,preM
      logical IsVacuum
c     Output:
      real*8 rho_0(N_STATES),vel_0(N_STATES),pre_0(N_STATES)
      real*8 snd_0(N_STATES)
      real*8 w_0(N_EDGES)
      integer signif_0(N_WAVES),wave_type_0(N_WAVES)

c     left states and wave
      call statesO(rhoL,velL,preL,sndL, velM,preM, IsVacuum,DIR_L,
     &     rho_0,vel_0,pre_0,snd_0, w_0,signif_0,wave_type_0)

c     right states and wave
      call statesO(rhoR,velR,preR,sndR, velM,preM, IsVacuum,DIR_R,
     &     rho_0,vel_0,pre_0,snd_0, w_0,signif_0,wave_type_0)

c     center wave
      call waveCD(
     &     rho_0(2),vel_0(2),pre_0(2),snd_0(2),
     &     rho_0(3),vel_0(3),pre_0(3),snd_0(3),
     &     IsVacuum, 
     &     w_0(3),w_0(4),signif_0(2),wave_type_0(2))

      end
