      subroutine inv2st(inv1,inv2,inv3,inpOrder, rho,vel,pre,snd)
      implicit none
c
c     inv2 = invEntr =  
c            "invariant Entropy" = pre/ rho**g / invEntrConst
c                             where invEntrConst = 1 or gamma-1
c   
c     inpOrder = +1, inv1=invL, inv3=invR 
c       input ORDER: left  Riemann Invarient first then 
c                    right Riemann Invarient 
c     inpOrder = -1, inv1=invR, inv3=invL
c       input ORDER: right Riemann Invarient first then 
c                    left  Riemann Invarient 
c
c     invL = LEFT  going Riemann invarient (along dx/dt=vel-snd)
c          = el-vel
c     invR = RIGHT going Riemann invarient (along dx/dt=vel+snd)
c          = el+vel
c
      include 'mat_constants.h'

c     Input:
      real*8 inv1,inv2,inv3
      integer inpOrder
c     Output:
      real*8 rho,vel,pre,snd
c     Auxiliary:
      real*8 el, preOverRho, thermSt2snd
c-----

      el = 0.5d0*(inv3+inv1) 
      vel = 0.5d0*(inv3-inv1) * inpOrder

      preOverRho = rc5*el**2
      rho = ( preOverRho/(inv2*invEntrConst) )**rc6
      pre = preOverRho*rho

      snd = thermSt2snd(rho,pre)

      end
