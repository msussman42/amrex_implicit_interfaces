      subroutine waveCD(rhoL,velL,preL,sndL,
     &     rhoR,velR,preR,sndR,IsVacuum,
     &     wL,wR,signif,wave_type)
      implicit none
c
c-----------------------------------------------------------------------
c                              INPUT 
c-----------------------------------------------------------------------
c
c         LEFT STATE VARIABLES
c         ----------------------
c    rhoL  density
c    velL  velocity
c    preL  pressure
c    sndL  speed of sound
c
c         RIGHT STATE VARIABLES
c         ---------------------
c    rhoR  density
c    velR  velocity
c    preR  pressure
c    sndR  speed of sound
c
c-----------------------------------------------------------------------
c                            OUTPUT 
c-----------------------------------------------------------------------
c
c    wL  velocity of the left-edge of the CD
c    wR  velocity of the right-edge of the CD
c
c-----------------------------------------------------------------------
c
c    signif
c            = SIGNIFICANT if the wave is signigicant 
c            = INSIGNIFICANT otherwise
c       THIS INFORMATION IS GIVEN FOR WAVE TRACKING PURPOSES: 
c       signif = SIGNIFICANT: the wave should be tracked
c       signif = INSIGNIFICANT: the wave should not be tracked
c
c-----------------------------------------------------------------------
c
c    wave_type  -- type of the wave
c            = RFN  rarefaction fan
c            = SFR  shock front
c            = CD   contact disc.
c
      include 'wave_constants.h'

c     Input.
      real*8 rhoL,velL,preL,sndL,rhoR,velR,preR,sndR
      logical IsVacuum
c     Output.
      real*8 wL,wR
      integer signif, wave_type

      wL  = velL
      wR  = velR

      wave_type = CD

      if (IsVacuum) then
         signif = INSIGNIFICANT
      else
         signif = SIGNIFICANT
      endif

      end
