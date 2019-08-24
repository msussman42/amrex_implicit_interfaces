      function pre2densExpan(comp,rhoC)
      implicit none
c     Computes the density at state W;
c         state W lies ON the EXPANSION portion of a wave-curve;
c         the wave-curve is centered at state C.
c 
c     C is the state at which the wave-curve is centered.
c     (With respect to the wave, 
c     state C could be either Ahead or Behind.)
c
c     W is the state ON the wave-curve.
c     (W is the side of the wave opposite to C)
c     (With respect to the wave, 
c     state W could be either Behind or Ahead.)
c
c     The state at which the wave-curve is centered is defined
c     by a triple {rhoC,velC,preC} (+ auxiliary impC).
c     The wave curve is parametrized with compression, comp = pre/preC.
c
c------------------------------------------------------------------------
c------------------------------------------------------------------------ 
c 
c     Input: 
c             comp --- compression pre/preC, 
c                      where pre --- pressure at state W
c             rhoC --- density at state C
c
c 
c     Output: pre2densExpan --- density at state W         
c
      include 'mat_constants.h'

c     Input      
      real*8 comp,rhoC
c     Output
      real*8 pre2densExpan
c     Auxiliary
      real*8 compPower

      pre2densExpan = rhoC*( compPower(comp,gi) )

      end
