      function pre2relMassFluxCompr(comp)
      implicit none

c     Computes the RELATIVE mass flux from state C to state W;
c         state W lies ON the COMPRESSION portion of a wave-curve;
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
c     by a triple {rhoC,velC,preC} (+ auxiliary sndC).
c     The wave curve is parametrized with compression, comp = pre/preC.
c
c------------------------------------------------------------------------
c------------------------------------------------------------------------ 
c 
c     Input: 
c             comp --- compression pre/preC, 
c                      where pre --- pressure at state W
c 
c     Output: pre2relMassFluxCompr --- massflux from state C to state W 
c                         (relative to) over the impedance at state C       
c

      include 'mat_constants.h'

c     Input      
      real*8 comp
c     Output
      real*8 pre2relMassFluxCompr

      pre2relMassFluxCompr = dsqrt( 1.d0+gg2*(comp-1.d0) )

      end
