      subroutine stPrim2stCons(rho,vel,pre, rhoCons,momCons,eneCons)
      implicit none
c
c     Input: PRIMITIVE STATE VARIABLES
c            rho    density 
c            vel    momentum
c            pre    pressure
c
c     Output: CONSERVATIVE STATE VARIABLES
c            rhoCons density
c            momCons momentum
c            eneCons energy
c
      include 'mat_constants.h'

c     Input:
      real*8 rho,vel,pre
c     Output:
      real*8 rhoCons,momCons,eneCons
c------
      
      rhoCons = rho
      momCons = vel*rho
      eneCons = pre/g1+0.5d0*rho*(vel**2) 

      end
