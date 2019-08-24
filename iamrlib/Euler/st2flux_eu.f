      subroutine st2flux(rho,vel,pre, DxDt,  fRho,fMom,fEne)
      implicit none
c
c     Input: PRIMITIVE STATE VARIABLES
c            rho    density 
c            vel    momentum
c            pre    pressure
c
c     DxDt --- slope of the ray on which the FLUX is to be calculated
c
c     Output: FLUX
c            fRho   density flux
c            fMom   momentum flux 
c            fEne   energy flux 
c
c     Auxiliary: CONSERVATIVE STATE VARIABLES
c            rhoCons density
c            momCons momentum
c            eneCons energy
c

c     Input:
      real*8 rho,vel,pre,DxDt
c     Output:
      real*8 fRho,fMom,fEne
c     Auxiliary:
      real*8 rhoCons,momCons,eneCons
c------
      
      call stPrim2stCons(rho,vel,pre, rhoCons,momCons,eneCons)

      fRho = vel*rhoCons             -DxDt*rhoCons 
      fMom = vel*momCons+pre         -DxDt*momCons 
      fEne = vel*(eneCons+pre)       -DxDt*eneCons

      end
