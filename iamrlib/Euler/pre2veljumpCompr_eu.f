      function pre2veljumpCompr(comp,sndC)
      implicit none
c
c     Computes the velocity jump (= dir*(velW-velC)) at state W;
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
c     by a pair {rhoC,velC,preC} (+ auxiliary sndC).
c     The wave curve is parametrized with compression, comp = pre/preC.
c
c------------------------------------------------------------------------
c------------------------------------------------------------------------ 
c 
c     Input: 
c             comp --- compression pre/preC, 
c                      where pre --- pressure at state W
c             sndC --- speed of sound at state C
c
c 
c     Output: pre2densCompr --- velocity jump at state W
c                                   (= dir*(velW-velC))
c
      include 'mat_constants.h'

c     Input      
      real*8 comp,sndC
c     Output
      real*8 pre2veljumpCompr
c     Auxiliary
      real*8 rmfl, pre2relMassFluxCompr

      rmfl = pre2relMassFluxCompr(comp)
      pre2veljumpCompr = sndC*( gi*(comp-1.d0)/rmfl )

      end
