      subroutine setNonDimenConsts(time,timeF,nIS,uIS,xID)
      implicit none

c----------------------------------------------------------------------
c     Fill non dimensionalization constants.
c----------------------------------------------------------------------

      include 'mat_constants.h'
      include 'input_output.h'
      include 'initial.h'
      include 'wave_constants.h'
      include 'referenceValues.h'

c     Input:
      real*8 time,timeF
      integer nIS
      real*8 uIS(0:N_INIT_STATE,1:N_EXT_VARS), xID(0:N_INIT_DISC)
c-----

c     Primatives:
      lo_0   = 1.d0
      wo_0   = 1.d0
      rho_0  = 1.d0

c     Derived:
      po_0 = rho_0*wo_0**2
      to_0 = lo_0/wo_0
      vo_0 = 1.0d0/rho_0

c     For State Variables:
      uo_0(1) = rho_0
      uo_0(2) = wo_0
      uo_0(3) = po_0

      open(fid, file='scales.dat')
      write(fid,1000) rho_0,wo_0,lo_0,vo_0,po_0,to_0
      close(fid)

 1000 format('rho_0 =',g24.16E2,
     &     /,'wo_0 =',g24.16E2,/,'lo_0 =',g24.16E2,
     &     /,'vo_0 =',g24.16E2,/,'po_0 =',g24.16E2,
     &     /,'to_0 =',g24.16E2)

      end
