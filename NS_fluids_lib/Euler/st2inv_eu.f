      function st2inv(vel,snd,dir)
      implicit none
c
c     invL = LEFT  going Riemann invarient (along dx/dt=vel-snd)
c          = -vel+el
c     invR = RIGHT going Riemann invarient (along dx/dt=vel+snd)
c          = +vel+el
c     where el = gg1*snd
c     dir = +1, inv=invR
c     dir = -1, inv=invL
c
      include 'mat_constants.h'

c     Input:
      real*8 vel,snd
      integer dir
c     Output:
      real*8 st2inv
c-----

      st2inv = rc1*snd + dir*vel 

      end
