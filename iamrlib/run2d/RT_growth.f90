      subroutine fy(y,alpha1,alpha2,Q,S,ff)
      IMPLICIT NONE

      Real*8 y,alpha1,alpha2,Q,S,ff

       ! page 444, Chandrasekhar, (121)
      ff=y**4+4.0*alpha1*alpha2*(y**3)+ &
       2.0*(1.0-6.0*alpha1*alpha2)*y*y- &
       4.0*(1.0-3.0*alpha1*alpha2)*y+ &
       (1.0-4.0*alpha1*alpha2)+ &
       Q*(alpha1-alpha2)+(Q**(1.0/3.0))*S

      return
      end


      PROGRAM RT_growth

      IMPLICIT NONE

      Real*8 :: mu1,mu2,rho1,rho2,tension,gravity
      Real*8 :: nu,pi,k,alpha1,alpha2
      Real*8 :: Q,S,kcrit,ylo,yhi,flo,fhi,ymid,fmid,n
      Real*8 :: atwood
      integer :: maxiter,numiter,option

      option=1

       ! page 441, "Hydrodynamic and Hydromagnetic stability"
       ! S. Chandrasekhar
       ! material 1 is the lower fluid
       ! material 2 is the upper fluid
       ! atwood=(rho2 - rho1)/(rho2 + rho1)
       ! (rho2+rho1)atwood = rho2-rho1
       ! rho1(1+atwood)=rho2(1-atwood)
       ! rho1=rho2(1-atwood)/(1+atwood)
      pi=4.0*atan(1.0)
      if (option.eq.1) then
       k=pi/4.0
       mu2=0.01
       rho2=1.0
       rho1=0.001225
       gravity=980.0
       tension=72.8
      else if (option.eq.2) then
       k=0.5
       gravity=1.0
       rho2=1.0
       mu2=1.0
       rho1=0.0
       tension=1.0
      else if (option.eq.3) then
       k=0.1
       gravity=1.0
       rho2=1.0
       mu2=1.0
       atwood=0.05
       rho1=rho2*(1.0-atwood)/(1.0+atwood)
       tension=1.0E-14
      else
       print *,"option invalid"
       stop
      endif

      atwood=(rho2-rho1)/(rho2+rho1)
      nu=mu2/rho2
      mu1=nu*rho1
      S=tension/((rho1+rho2)*((gravity*(nu**4))**(1.0/3.0)))
      alpha1=rho1/(rho1+rho2)
      alpha2=rho2/(rho1+rho2)
      Q=gravity/((k**3)*(nu**2))
      kcrit=sqrt((alpha2-alpha1)*(gravity**(2.0/3.0))/ &
       (S*(nu**(4.0/3.0))))
      print *,"atwood number ",atwood
      print *,"alpha1,alpha2,Q,S ",alpha1,alpha2,Q,S
      print *,"rho1,rho2,mu1,mu2,nu ",rho1,rho2,mu1,mu2,nu
      print *,"gravity,tension ",gravity,tension
      print *,"kcrit= ",kcrit
      print *,"k= (e^(ikx)) ",k
      if (k.ge.kcrit) then
       print *,"k must be less than kcrit"
       stop
      endif
      print *,"wave length is 2pi/k= ",2.0*pi/k
      ylo=1.0
      yhi=2.0
      call fy(ylo,alpha1,alpha2,Q,S,flo)
      if (flo.ge.0.0) then
       print *,"k too large (flo error)"
       stop
      endif
      call fy(yhi,alpha1,alpha2,Q,S,fhi)
      do while (fhi.le.0.0) 
       yhi=yhi*2.0
       call fy(yhi,alpha1,alpha2,Q,S,fhi)
      enddo
      print *,"ylo,flo (init) ",ylo,flo
      print *,"yhi,fhi (init) ",yhi,fhi

      numiter=0
      maxiter=100
      do while (numiter.le.maxiter)
       ymid=0.5*(ylo+yhi)
       call fy(ymid,alpha1,alpha2,Q,S,fmid)
       if (fmid.eq.0.0) then
        ylo=ymid
        yhi=ymid
       else if (fmid.gt.0.0) then
        yhi=ymid
       else
        ylo=ymid
       endif
       numiter=numiter+1
      enddo 
      print *,"ymid,fmid ",ymid,fmid

      n=(ymid*ymid-1.0)*k*k*nu

      print *,"the growth rate (e^(nt)) is n= ",n
      END PROGRAM
