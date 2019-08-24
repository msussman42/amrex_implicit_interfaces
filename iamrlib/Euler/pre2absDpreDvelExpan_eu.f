      function pre2absDpreDvelExpan(comp,impC)
      implicit none
c
c     Computes the D(pressure)/D(velocity) at state W;
c         state W lies ON the EXPANSIONportion of a wave-curve;
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
c     by a triple {rhoC,velC,preC} (+ auxiliary impC and mflC).
c     The wave curve is parametrized with compression, comp = pre/preC.
c
c------------------------------------------------------------------------
c------------------------------------------------------------------------ 
c 
c     Input: 
c             comp --- compression pre/preC, 
c                      where pre --- pressure at state W
c             impC --- impedance at state C
c             mflC --- massflux at state C
c 
c     Output: rho2absDmomDrhoComp  --- 
c                    Abs( D(momentum)/D(density) )  at state W 
c
      include 'mat_constants.h'

c     Input      
      real*8 comp,impC
c     Output
      real*8 pre2absDpreDvelExpan
c     Auxiliary
      real*8 compPower

      pre2absDpreDvelExpan = impC*( compPower(comp,gg2) )

      end
