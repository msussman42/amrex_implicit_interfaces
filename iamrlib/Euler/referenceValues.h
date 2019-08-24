c
c     Reference values for:
c
c     lo_0   = length
c     to_0   = time
c     rho_0 = density ( rho )
c     vo_0   = specific volumn
c     wo_0   = velocity
c     po_0   = pressure
c
c     uo_0(1:N_EQUATIONS) = conserved state variables
c
c     Used in non-dimensilization.
c     These are demensional quantities.

      double precision          lo_0,to_0,rho_0,vo_0,wo_0,po_0
      double precision          uo_0(1:N_EQUATIONS)
      common /NonDimConst_CM/   lo_0,to_0,rho_0,vo_0,wo_0,po_0,uo_0
