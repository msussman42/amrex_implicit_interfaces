      program main
      IMPLICIT NONE

c J/(kg K) 
      real*8 cpl,cpg,cps
c W/(m K) 
      real*8 kl,kg,ks
c kg/m^3
      real*8 rhol,rhog,rhos
c Pa s
      real*8 mul,mug
c N/m
      real*8 sigma0
c Kelvin
      real*8 TC,TM
      real*8 R
c J/kg
      real*8 Lh
      real*8 g
      real*8 pr,st,B0,tau_c

      cpl=801.0
      cpg=1230.0
      cps=801.0

      kl=64.0
      ks=32.0
      kg=0.101

      rhol=2583.0
      rhog=0.210
      rhos=2350.0

      mul=5.75D-4
      mug=6.1D-5

      sigma0=0.831

      TC=1422.0
      TM=1683.4

      R=4.73D-3 

      Lh=1.805D+6

      pr=cpl*mul/kl
      st=cpl*(TM-TC)/Lh

      g=9.8
      B0=rhol*g*R*R/sigma0
      tau_c=rhol*cpl*R*R/kl

      print *,"pr= ",pr
      print *,"st= ",st
      print *,"B0= ",B0
      print *,"tau_c= ",tau_c

      end

