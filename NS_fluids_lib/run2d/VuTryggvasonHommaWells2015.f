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

c     cpl=801.0
      cpl=910.0
c     cpg=1230.0
      cpg=910.0
c     cps=801.0
      cps=910.0

c     kl=64.0
      kl=152.3
c     ks=32.0
      ks=76.2
c     kg=0.101
      kg=1.52

      rhol=2583.0
c     rhog=0.210
      rhog=129.15
      rhos=2350.0

c     mul=5.75D-4
      mul=2.176D-3
c     mug=6.1D-5
      mug=1.08D-4

c     sigma0=0.831
      sigma0=1.74

c     TC=1422.0
      TC=1453.0
      TM=1683.4

c     R=4.73D-3 
c     R=6.0D-3 
c     R=5.1D-3 
      R=9.5D-3 

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

