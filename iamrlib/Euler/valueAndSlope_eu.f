      subroutine valueAndSlope(
     &     pre,rhoC,velC,preC,sndC,impC,dir,centSide, 
     &     absDpreDvel,vel)
      implicit none
c
c     Works with the projections of wave-curves on the (pre,vel)-plane.
c     Wave curves are parametrized with pressure, pre.
c
c     Given pre (state W on the wave-curve is defined), the subroutine
c     computes velocity, vel, and abs(D(pressure)/D(velocity)), absDpreDvel.
c
c         state W lies ON a general wave-curve and is defined by pre;
c         the wave-curve is centered at state C.
c     
c     C is the state at which the wave-curve is centered.
c     With respect to the wave, 
c     state C could be either Ahead or Behind.
c             rhoC --- density at state C
c             velC --- velocity at state C 
c             preC --- pressure  at state C
c             sndC --- sound speed at state C
c             impC --- impedance at state C 
c
c     W is the state ON the wave-curve at which the VALUE is computed.
c     (W is the side of the wave opposite to C)
c     With respect to the wave, 
c     state W could be either Behind or Ahead.
c             pre --- pressure at state W
c
c     The state at which the wave-curve is centered is defined
c     by a triple {rhoC,velC,preC} (+ auxiliary impC).
c     The wave curve is parametrized with the pressure, pre. 
c 
c     Whether C is Ahead (Behind) and W is Behind (Ahead)
c     depends on the value of centSide
c             centSide =  -1 state C is Behind, state W is Ahead
c                         +1 state C  is Ahead, state W is Behind
c
c     dir  --- -1 1-wave-curve
c              +1 2-wave-curve
c 
c     Output: vel      --- value of the velocity on the wave-curve
c             absDpreDvel --- abs(D(pressure)/D(velocity)) on the wave-curve
c      
c     Auxiliary: mfl --- massflux at between states C and W
c                comp --- compression pre/preC, 
c                      where pre --- pressure at state W
c
c

      include 'parameters.h'

c     Input:
      real*8 pre,rhoC,velC,preC,sndC,impC
      integer dir,centSide
c     Output:
      real*8 absDpreDvel,vel
c     Auxiliary:
      real*8 veljump
      real*8 pre2veljumpExpan,pre2veljumpCompr
      external pre2veljumpExpan,pre2veljumpCompr
      real*8 pre2absDpreDvelExpan,pre2absDpreDvelCompr
      external pre2absDpreDvelExpan,pre2absDpreDvelCompr
      real*8 veljump2vel
      real*8 onWaveCurve
      external onWaveCurve
c----

      veljump     = onWaveCurve(pre,preC,sndC,centSide,
     &              pre2veljumpExpan,pre2veljumpCompr)
      absDpreDvel  = onWaveCurve(pre,preC,impC,centSide,
     &              pre2absDpreDvelExpan,pre2absDpreDvelCompr)
      vel = veljump2vel(velC,veljump,dir)

      end
