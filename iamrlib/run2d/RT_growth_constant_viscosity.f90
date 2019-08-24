      subroutine f_of_n(k,n,alpha1,alpha2,nu1,nu2,grav,ff)
      IMPLICIT NONE

      Real*8 k,n,alpha1,alpha2,nu1,nu2,grav,ff
      Real*8 q1mk_over_n
      Real*8 q2mk_over_n
      Real*8 alpha_q_term

       ! page 443, Chandrasekhar, (113) (multiply through by n)
      q1mk_over_n=1.0d0/(nu1*(k+sqrt(k*k+n/nu1)))
      q2mk_over_n=1.0d0/(nu2*(k+sqrt(k*k+n/nu2)))
      alpha_q_term=alpha2*q1mk_over_n+alpha1*q2mk_over_n
      ff=grav*k*(alpha2-alpha1)*alpha_q_term-n*n*alpha_q_term- &
       4.0d0*n*k*alpha1*alpha2

      return
      end


      PROGRAM RT_growth

      IMPLICIT NONE

      Real*8 :: mu1,mu2,rho1,rho2,gravity
      Real*8 :: nu1,nu2,pi,k,alpha1,alpha2
      Real*8 :: kstart,kstop,delta_k
      Real*8 :: nlo,nhi,flo,fhi,nmid,fmid,n,nstokes
      Real*8 :: atwood
      integer :: maxiter,numiter,n_k,kloop,verbose

       ! page 441, "Hydrodynamic and Hydromagnetic stability"
       ! S. Chandrasekhar
       ! material 1 is the lower fluid
       ! material 2 is the upper fluid
       ! atwood=(rho2 - rho1)/(rho2 + rho1)
       ! (rho2+rho1)atwood = rho2-rho1
       ! rho1(1+atwood)=rho2(1-atwood)
       ! rho1=rho2(1-atwood)/(1+atwood)
       ! Gerya and Yuen
       ! L=500 km
       ! rho=3200 kg/m^3
       ! mu=10^21 Pa s
       ! g=9.8 m/s^2
       ! Re=rho L U/mu
       ! T=L/U
       ! u_t = -grad p/rho + mu (uxx+uyy)/rho + g yhat
       ! U^2/L u_t = ... + mu U/(rho L^2) + g yhat
       ! u_t = (1/Re) + Lg/U^2  yhat
       ! U^2=Lg
       ! Re=rho L^(3/2)g^(1/2)/mu=3.54E-9
       ! (1/Re)=2.82E+8
       ! T=L/U=(L/g)^(1/2)=225.9 seconds
       ! most dangerous mode: k=3.0E-7  n=4.3E-5
       ! wave length = 2 pi/k  k=2 pi/wave length
       ! if wave length=1/500 => k=500 (2 pi)=3141.6 and n=8.8E-15
       ! time to grow by factor of 2.7: nt=1  t=1/n dimensionless
       ! t=T/n  seconds=226/8.8E-15 sec=814 million years
       ! if wave length=1/5 => k=5 (2 pi)=31.4 and n=8.8E-13
       ! time to grow by factor of 2.7: nt=1  t=1/n dimensionless
       ! t=T/n  seconds=226/8.8E-13 sec=8 million years
       ! See also Turcotte and Schubert, GEODYNAMICS, 3rd edition,
       ! (6.158) in the limit as b->infinity. (page 291):
       ! t_a = 4 mu 2 pi/((rho1-rho2) g lambda)

      verbose=0

      pi=4.0d0*atan(1.0d0)
      kstart=0.0d0
!      kstop=2.0d0*pi
!      kstop=5.0E-7
!     kstop=3142.0
      kstop=31.4
      n_k=1000
      delta_k=(kstop-kstart)/n_k
      do kloop=0,n_k-1

       k=kstart+(kloop+0.5d0)*delta_k
       mu1=2.82E+8
       mu2=2.82E+8
       rho1=1.0d0
       rho2=33.0d0/32.0d0
       nu1=mu1/rho1
       nu2=mu2/rho2
       gravity=1.0d0

       atwood=(rho2-rho1)/(rho2+rho1)
       alpha1=rho1/(rho1+rho2)
       alpha2=rho2/(rho1+rho2)
       nlo=0.0
       nhi=2.0
       call f_of_n(k,nlo,alpha1,alpha2,nu1,nu2,gravity,flo)
       if (flo.le.0.0) then
        print *,"f(0) should be positive"
        stop
       endif
       call f_of_n(k,nhi,alpha1,alpha2,nu1,nu2,gravity,fhi)
       do while (fhi.ge.0.0) 
        nhi=nhi*2.0
        call f_of_n(k,nhi,alpha1,alpha2,nu1,nu2,gravity,fhi)
       enddo
       if (verbose.eq.1) then
        print *,"nlo,flo (init) ",nlo,flo
        print *,"nhi,fhi (init) ",nhi,fhi
       endif

       numiter=0
       maxiter=200
       do while (numiter.le.maxiter)
        nmid=0.5*(nlo+nhi)
        call f_of_n(k,nmid,alpha1,alpha2,nu1,nu2,gravity,fmid)
        if (fmid.eq.0.0) then
         nlo=nmid
         nhi=nmid
        else if (fmid.lt.0.0) then
         nhi=nmid
        else
         nlo=nmid
        endif
        numiter=numiter+1
       enddo
       if (verbose.eq.1) then
        print *,"k,nmid,fmid ",k,nmid,fmid
       endif

       n=nmid

       if (verbose.eq.1) then
        print *,"the growth rate (e^(nt)) is k,n= ",k,n
       endif
       nstokes=4.0d0*mu1*k/((rho2-rho1)*gravity)
       nstokes=1.0d0/nstokes
       print *,k,n,nstokes

      enddo !kloop

      END PROGRAM
