      subroutine readEOSConst()
      implicit none

c     Reads Equation Of State (EOS) for the material. 

c     EOS is taken to be an EOS for Eulerian gas dynamics.

      include 'mat_constants.h'

      g=1.4d0

      g1=g-1.d0
      gg1=g1/2.d0/g
      gg2=(g+1.d0)/2.d0/g
      gi=1.d0/g
      gk=(g-1.d0)/(g+1.d0)
      gc=g*g1
      b1=g1/g
      b2=0.5d0/g
      rc1=2.d0/g1
      rc2=1.d0/gg1
      rc3=g*(g+1)/4.d0
      rc4=2.d0/gg2
      rc5=1.d0/g/rc1**2  
      rc6=1.d0/(g-1.d0)

      invEntrConst = 1.d0 ! or g-1.d0

      end

