      subroutine slope2state(     
     &     rho_0,vel_0,pre_0,snd_0,w_0,signif_0,wave_type_0, DxDt,
     &     rho,vel,pre,snd)
      implicit none
c
c-----------------------------------------------------------------------
c                              INPUT 
c-----------------------------------------------------------------------
c
c          STATE VARIABLES
c          ---------------
c    rho_0(1:N_STATES)   density
c    vel_0(1:N_STATES)   velocity
c    pre_0(1:N_STATES)   pressure
c    snd_0(1:N_STATES)   speed of sound
c    1 - Left (Init.) State;  2 - Left Middle State; 
c                3 - Right ddle State;  4  - Right (Init.) State
c
c    NOTE: Left (initial) State  EQUALS Left  state of INPUT
c          Right (initial) State EQUALS Right state of INPUT
c
c         WAVE STRUCTURE 
c         --------------
c    1 - left-facing wave; 2 - center wave; 3 - right-facing wave
c
c    wave_type_0(1:N_WAVES) -- type of the wave
c     = SFR  shock front
c     = RFN  rarefaction fan
c     = CD   contact discontinuity
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
c    w_0(1:N_EDGES) velocities of the fan-edges:
c        w_0(1) = velocity of the outer  edge of the left-f fan
c        w_0(2) = velocity of the innner edge of the left-f fan
c        w_0(3) = velocity of the left   edge of the center fan
c        w_0(4) = velocity of the right  edge of the center fan
c        w_0(5) = velocity of the inner  edge of the right-f fan
c        w_0(6) = velocity of the outer  edge of the right-f fan
c
c    DxDt --- slope of the ray on which the STATE is to be calculated
c
c-----------------------------------------------------------------------
c                              OUTPUT 
c-----------------------------------------------------------------------
c
c    STATE on the ray with the slope DxDt
c    ---------------------------------------------------------------
c  rho    density                         
c  vel    velocity
c  pre    pressure
c  snd    speed of sound
c
c
      include 'wave_constants.h'

c     Inoput:
      real*8 rho_0(N_STATES),vel_0(N_STATES),pre_0(N_STATES)
      real*8 snd_0(N_STATES)
      real*8 w_0(N_EDGES)
      integer signif_0(N_WAVES),wave_type_0(N_WAVES)
      real*8 DxDt
c     Output:
      real*8 rho,vel,pre,snd
c     Auxiliary:
      integer i,j,ka,k,dir

      do i=1,N_EDGES
         if (DxDt .lt. w_0(i)) goto 20
      enddo
 20   continue

c$$$      if (i .gt. 1) then
c$$$         write(*,*) 'w_0(',i-1,')=',w_0(i-1)
c$$$      endif
c$$$      write(*,*) 'DxDt =',DxDt
c$$$      write(*,*) 'w_0(',i,')=',w_0(i)
      
      if (IS_L_EDGE(i)) then
c        inside a state
c$$$         write(*,*) 'inside a state:'
         k = EDGE_TO_L_STATE(i)
c$$$         write(*,*) 'state k =',k
         call stPrimExt2stPrimExt(rho_0(k),vel_0(k),pre_0(k),snd_0(k),
     &        rho,vel,pre,snd) 
      else
c        inside a fan
c$$$         write(*,*) 'inside a fan:'
         j = EDGE_TO_L_WAVE(i)
c$$$         write(*,*) 'wave j =',j
         if ( wave_type_0(j) .ne. RFN) then
            write(*,*)
            write(*,*) 'i =',i,'  j =',j
            pause 'slope2state: it is NOT a FAN!'
         endif
         ka = WAVE_TO_STATE_A(j)
         dir = WAVE_TO_DIR(j)   
c$$$         write(*,*) 'state ka =',ka
c$$$         write(*,*) 'dir =',dir
         call insideFan(
     &        rho_0(ka),vel_0(ka),pre_0(ka),snd_0(ka),dir,AHEAD,DxDt, 
     &        rho,vel,pre,snd)
      endif
      
      end
