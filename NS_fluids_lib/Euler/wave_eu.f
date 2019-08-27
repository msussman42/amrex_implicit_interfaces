      subroutine wave(rhoB,velB,preB,sndB, rhoA,velA,preA,sndA, dir,
     &     wB,wA,signif,wave_type)
      implicit none
c
c-----------------------------------------------------------------------
c                              INPUT 
c-----------------------------------------------------------------------
c
c         BEHIND STATE VARIABLES
c         ----------------------
c    rhoB  density
c    velB  velocity
c    preB  pressure
c    sndB  speed of sound
c
c         AHEAD STATE VARIABLES
c         ---------------------
c    rhoA  density
c    velA  velocity
c    preA  pressure
c    sndA  speed of sound
c
c    dir  wave direction parameter
c        =-1 for left-facing wave
c        = 1 for right-facing wave
c
c-----------------------------------------------------------------------
c                            OUTPUT 
c-----------------------------------------------------------------------
c
c    wB velocity of the behind-edge of the fan
c    wA velocity of the ahead-edge of the fan
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
c            = CD   constant discontinuity
c
      include 'wave_constants.h'

c     Input.
      real*8 rhoB,velB,preB,sndB,rhoA,velA,preA,sndA
      integer dir
c     Output.
      real*8 wB,wA
      integer signif, wave_type
c     Auxiliary.
      real*8 primExtSt2W, veljump2vel

      if (preB .gt. preA) then
c        front
         wave_type = SFR
         wA = primExtSt2W(rhoB,velB,preB,sndB,rhoA,velA,preA,sndA,dir)
         wB = wA
      else
c        fan
         wave_type = RFN
         wA = veljump2vel(velA,sndA,dir)
         wB = veljump2vel(velB,sndB,dir)
      endif

      signif = SIGNIFICANT

      end
