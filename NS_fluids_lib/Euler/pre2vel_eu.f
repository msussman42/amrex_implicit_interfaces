      function pre2vel(pre,preC,velC,mflC,dir)
      implicit none
c
c     Computes the velocity at state W; 
c         state W lies ON a general wave-curve;
c         the wave-curve is centered at state C.
c     
c     C is the state at which the wave-curve is centered.
c     With respect to the wave, 
c     state C could be either Ahead or Behind.
c             preC --- pressure at state C
c             velC --- velocity at state C 
c             mflC --- massflux at state C
c
c     W is the state ON the wave-curve at which the VALUE is computed.
c     (W is the side of the wave opposite to C)
c     With respect to the wave, 
c     state W could be either Behind or Ahead.
c             pre --- pressure at state W
c
c     The state at which the wave-curve is centered is defined
c     by a triple {rhoC,velC,preC} (+ auxiliary impC and mflC).
c     The wave curve is parametrized with the pressure, pre. 
c 
c     dir  --- -1 1-wave-curve
c              +1 2-wave-curve
c 
c     Output: pre2vel --- the velocity on the wave-curve
c

c     Input      
      real*8 pre,preC,velC,mflC
      integer dir
c     Output
      real*8 pre2vel

      pre2vel = velC + dir*(pre-preC)/mflC

      end
