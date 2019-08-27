      function primExtSt2W(rho,vel,pre,snd,rhoC,velC,preC,sndC,dir)
      implicit none
c
c     Computes the Eulerian Speed of Discontinuity 
c                        connecting state C and state W;
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
c     Input:
c               STATE C
c               -----------------------
c            rhoC    density 
c            velC    velocity 
c            preC    pressure
c            sndC    speed of sound
c               STATE W
c               -----------------------
c            rho    density
c            vel    velocity
c            pre    pressure
c            snd    speed of sound
c
c            dir =+1 for right-facing wave, =-1 for left-facing wave
c     Output:
c            primExtSt2W --- Eulerian Speed of Discontinuity 
c
      include 'mat_constants.h'  
       
c     Input
      real*8 rho,vel,pre,snd,rhoC,velC,preC,sndC
      integer dir
c     Output
      real*8 primExtSt2W
c     Auxiliary
      real*8 comp, relSpDisc,pre2relMassFluxCompr, veljump2vel
c-----

      comp=pre/preC

      relSpDisc =  sndC*pre2relMassFluxCompr(comp)
      primExtSt2W = veljump2vel(velC,relSpDisc,dir)

      end


