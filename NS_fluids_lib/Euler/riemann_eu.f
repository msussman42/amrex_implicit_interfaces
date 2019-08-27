      subroutine riemann( 
     &     rho_l0,vel_l0,pre_l0,snd_l0, rho_r0,vel_r0,pre_r0,snd_r0,
     &     rho_0,vel_0,pre_0,snd_0, w_0,signif_0,wave_type_0,IsVacuum)
      implicit none

c
c-----------------------------------------------------------------------
c                              INPUT 
c-----------------------------------------------------------------------
c
c    LEFT STATE VARIABLES
c    --------------------
c  rho_l0    density
c  vel_l0    velocity
c  pre_l0    pressure
c  snd_l0    sound speed
c
c    RIGHT STATE VARIABLES
c    ---------------------
c  rho_r0    density
c  vel_r0    velocity
c  pre_r0    pressure
c  snd_r0    sound speed
c
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
c     1 - Left (Init.) State;  
c     2 - Left Middle State; 
c     3 - Right Middle State;
c     4 - Right (Init.) State;
c
c    NOTE: Left (initial) State  EQUALS Left  state of INPUT
c          Right (initial) State EQUALS Right state of INPUT
c
c         WAVE STRUCTURE 
c         --------------
c    1 - left-facing outer wave; 
c    2 - central wave
c    3 - right-facing outer wave
c
c    wave_type_0(1:N_WAVES) -- type of the wave
c     = RFN  rarefaction fan
c     = SFR  shock front
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
c--------------------------------------------------------------------------
c                             AUXILIARY 
c--------------------------------------------------------------------------
c
c    MIDDLE STATES VARIABLES returned by subroutine newton:
c    -----------------------
c    rho{L,R}    density
c    vel    velocity
c    pre    pressure
c
      include 'parameters.h'
      include 'wave_constants.h'

c     Input:
      real*8 rho_l0,vel_l0,pre_l0,snd_l0,rho_r0,vel_r0,pre_r0,snd_r0
c     Output:
      real*8 rho_0(N_STATES),vel_0(N_STATES),pre_0(N_STATES)
      real*8 snd_0(N_STATES)
      real*8 w_0(N_EDGES)
      integer signif_0(N_WAVES),wave_type_0(N_WAVES)
      logical IsVacuum
c     Auxiliary:
      real*8 vel,pre
      real*8 preWC_min,preWC_max
c
c

      if (rho_l0.le.0.0) then
       print *,"rho_l0 cannot be negative riemann_eu.f: ",rho_l0
       stop
      endif
      if (pre_l0.le.0.0) then
       print *,"pre_l0 cannot be negative riemann_eu.f: ",pre_l0
       stop
      endif
      if (snd_l0.le.0.0) then
       print *,"snd_l0 cannot be negative riemann_eu.f: ",snd_l0
       stop
      endif
      if (rho_r0.le.0.0) then
       print *,"rho_r0 cannot be negative riemann_eu.f: ",rho_r0
       stop
      endif
      if (pre_r0.le.0.0) then
       print *,"pre_r0 cannot be negative riemann_eu.f: ",pre_r0
       stop
      endif
      if (snd_r0.le.0.0) then
       print *,"snd_r0 cannot be negative riemann_eu.f: ",snd_r0
       stop
      endif

      call readEOSConst()

      call newton(
     &     rho_l0,vel_l0,pre_l0,snd_l0, rho_r0,vel_r0,pre_r0,snd_r0,
     &     vel,pre,IsVacuum)

      call states(
     &     rho_l0,vel_l0,pre_l0,snd_l0, rho_r0,vel_r0,pre_r0,snd_r0,
     &     vel,pre,IsVacuum,
     &     rho_0,vel_0,pre_0,snd_0, w_0,signif_0,wave_type_0)

      end
