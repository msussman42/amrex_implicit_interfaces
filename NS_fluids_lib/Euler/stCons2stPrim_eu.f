      subroutine stCons2stPrim(rhoCons,momCons,eneCons, rho,vel,pre)
      implicit none
c
c     INput: CONSERVATIVE STATE VARIABLES
c            rhoCons density
c            momCons momentum
c            eneCons energy
c
c     Output: PRIMITIVE STATE VARIABLES
c            rho    density 
c            vel    momentum
c            pre    pressure
c

      include 'mat_constants.h'

c     Input:
      real*8 rhoCons,momCons,eneCons
c     Output:
      real*8 rho,vel,pre

c------
      
      rho = rhoCons 
      vel = momCons/rhoCons
      pre = g1*( eneCons-0.5d0*(momCons**2)/rhoCons )

      end
