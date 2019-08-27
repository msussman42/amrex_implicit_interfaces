      function veljump2vel(velC,veljump,dir)
      implicit none
c
c     Computes the velocity at state W; 
c         state W lies ON a general wave-curve;
c         the wave-curve is centered at state C.
c     
c     C is the state at which the wave-curve is centered.
c     With respect to the wave, 
c     state C could be either Ahead or Behind.
c             velC --- velocity at state C 
c             veljump --- velocity jump from state C to state W
c                        ( = dir*(vel-velC) )
c
c     W is the state ON the wave-curve at which the VALUE is computed.
c     (W is the side of the wave opposite to C)
c     With respect to the wave, 
c     state W could be either Behind or Ahead.
c
c     The state at which the wave-curve is centered is defined
c     by a triple {rhoC,velC,preC} (+ auxiliary sndC).
c     The wave curve is parametrized with the pressure, pre. 
c 
c     dir  --- -1 1-wave-curve
c              +1 2-wave-curve
c 
c     Output: veljump2vel --- the velocity on the wave-curve
c

c     Input      
      real*8 velC,veljump
      integer dir
c     Output
      real*8 veljump2vel

      veljump2vel = velC + dir*veljump

      end
