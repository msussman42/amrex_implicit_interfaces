      function thermSt2invEntr(rho,pre)
c
c     Given a thermodynamic state  {rho,pre}
c     this subroutine computes the "invariant Entropy" thermSt2invEntr.
c
c     Input:
c         THERMODYNAIC STATE VARIABLES:
c       rho          --- density
c       pre          --- pressure
c     Output:
c       thermSt2invEntr  --- "invariant Entropy"  = pre/ rho**g / invEntrConst
c                             where invEntrConst = 1 or gamma-1
c
c
      include 'mat_constants.h'
c
c     Input:
      real*8 rho,pre
c     Output:
      real*8 thermSt2invEntr
c
      thermSt2invEntr=pre/rho**g/invEntrConst

      end
