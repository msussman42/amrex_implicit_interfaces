      function onWaveCurve(pre,preC,valC,centSide,
     &     onWaveCurveExpan,onWaveCurveCompr)
      implicit none
c
c     Computes a VALUE at state W; 
c         state W lies ON a general wave-curve;
c         the wave-curve is centered at state C.
c     
c
c     onWaveCurveExpan computes the VALUE on an EXPANSION portion of the curve 
c     onWaveCurveCompr computes the VALUE on a COMPRESSION portion of the curve
c
c     C is the state at which the wave-curve is centered.
c     With respect to the wave, 
c     state C could be either Ahead or Behind.
c             preC --- pressure  at state C
c             valC --- VALUE at state C 
c
c     W is the state ON the wave-curve at which the VALUE is computed.
c     (W is the side of the wave opposite to C)
c     With respect to the wave, 
c     state W could be either Behind or Ahead.
c             pre --- pressure at state W
c
c     The state at which the wave-curve is centered is defined
c     by a triple {rhoC,velC,preC}.
c     The wave curve is parametrized with the pressure, pre. 
c 
c     Whether C is Ahead (Behind) and W is Behind (Ahead)
c     depends on the value of centSide
c             centSide =  -1 state C is Behind, state W is Ahead
c                         +1 state C  is Ahead, state W is Behind
c
c     Output: onWaveCurve  --- value on the wave-curve
c
c     Auxiliary: comp --- compression pre/preC 
c
      
c     Input      
      real*8 pre,preC,valC
      integer centSide
      external onWaveCurveExpan,onWaveCurveCompr
      real*8 onWaveCurveExpan,onWaveCurveCompr
c     Output
      real*8 onWaveCurve
c     Auxiliary
      real*8 comp

      comp = pre/preC

      if ( (comp-1.d0)*centSide .le. 0.d0) then
         onWaveCurve = onWaveCurveExpan(comp,valC)
      else
         onWaveCurve = onWaveCurveCompr(comp,valC)
      endif

      end
