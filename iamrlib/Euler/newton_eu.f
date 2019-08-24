      subroutine newton(r1,u1,p1,snd1,r2,u2,p2,snd2,
     &     uu,pp,IsVacuum)
      implicit none
c
c
      include 'parameters.h'
      include 'wave_constants.h'
c
      integer NITER
      parameter (NITER=1000)
      real*8 paccur
      parameter (paccur = pre_min)

c     Input:
      real*8 r1,u1,p1,snd1,r2,u2,p2,snd2
c     Output:
      real*8 uu,pp,ppsave
      logical IsVacuum
c     Auxiliary:
      real*8 c1,c2,z1,z2,uu1,uu2,dpp,dpp0
      logical wasVacuum
      integer i
c
c-------------------

      if (r1.le.0.0) then
       print *,"r1 cannot be negative newton_eu.f: ",r1
       stop
      endif
      if (p1.le.0.0) then
       print *,"p1 cannot be negative newton_eu.f: ",p1
       stop
      endif
      if (snd1.le.0.0) then
       print *,"snd1 cannot be negative newton_eu.f: ",snd1
       stop
      endif
      if (r2.le.0.0) then
       print *,"r2 cannot be negative newton_eu.f: ",r2
       stop
      endif
      if (p2.le.0.0) then
       print *,"p2 cannot be negative newton_eu.f: ",p2
       stop
      endif
      if (snd2.le.0.0) then
       print *,"snd2 cannot be negative newton_eu.f: ",snd2
       stop
      endif


      c1 = r1*snd1
      c2 = r2*snd2

      wasVacuum = .false.

      pp =(c2*p1+c1*p2-c1*c2*(u2-u1))/(c1+c2)

      do i=1,NITER
         if (pp .le. pre_min) then
            pp=pre_min
            if (wasVacuum) then
               IsVacuum = .true.
               return
            else
               wasVacuum=.true.
            endif
         else
            wasVacuum=.false.
         endif
c
         call valueAndSlope(pp,r1,u1,p1,snd1,c1,DIR_L,AHEAD,z1,uu1)
         call valueAndSlope(pp,r2,u2,p2,snd2,c2,DIR_R,AHEAD,z2,uu2)

         dpp=-z1*z2*(uu2-uu1)/(z1+z2)

         if (1.eq.0) then
          print *,"i= ",i
          print *,"pp,r1,u1,p1,snd1,c1,DIR_L,AHEAD,z1,uu1 ",
     &     pp,r1,u1,p1,snd1,c1,DIR_L,AHEAD,z1,uu1
          print *,"pp,r2,u2,p2,snd2,c2,DIR_R,AHEAD,z2,uu2 ",
     &     pp,r2,u2,p2,snd2,c2,DIR_R,AHEAD,z2,uu2
         endif

         if (dabs(dpp) .lt. paccur) then

            if (pp .le. pre_min) then
               IsVacuum = .true.
            else
               uu=(z1*uu1+z2*uu2)/(z1+z2)
               IsVacuum=.false.
            endif
            return
         endif

         ppsave=pp
         pp = pp + dpp
         if (i.eq.1) then
          dpp0=dpp
         endif

      enddo

      write(*,*)
      write(*,*) ' WARNING Number of iterations exceeded',NITER
      print *,"r1,u1,p1,snd1 ",r1,u1,p1,snd1
      print *,"r2,u2,p2,snd2 ",r2,u2,p2,snd2
      print *,"dpp0,dpp,paccur ",dpp0,dpp,paccur

      pp=ppsave
      if (pp .le. pre_min) then
       IsVacuum = .true.
      else
       uu=(z1*uu1+z2*uu2)/(z1+z2)
       IsVacuum=.false.
      endif
      return

1000  format(2g26.18E2)
      end

