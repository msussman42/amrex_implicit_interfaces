      function thermSt2snd(rho,pre)
c
c     Given a thermodynamic state  {rho,pre}
c     this subroutine computes the speed of sound thermSt2snd.
c
c     Input:
c         THERMODYNAIC STATE VARIABLES:
c       rho          --- density
c       pre          --- pressure
c     Output:
c       thermSt2snd  --- speed of sound    
c
c
      include 'mat_constants.h'
c
c     Input:
      real*8 rho,pre
c     Output:
      real*8 thermSt2snd
c
      thermSt2snd=dsqrt(pre/rho/gi)

      end
