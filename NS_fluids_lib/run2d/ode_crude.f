      program main
      IMPLICIT NONE

      real*8 reaction_rate
      real*8 dt
      integer nstep,iter
      real*8 YOLD,YNEW

      nstep=20000
      dt=0.0001d0
      reaction_rate=0.2d0
      YOLD=0.0d0
      YNEW=YOLD
      do iter=1,nstep
       YNEW=(YOLD+reaction_rate*dt)/(1.0d0+reaction_rate*dt)
       YOLD=YNEW
      enddo

      print *,"YNEW=",YNEW

      end

