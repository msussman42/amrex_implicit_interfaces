      program main
      IMPLICIT NONE

      real*8 mu_L,mu_V,rho_L,rho_V
      real*8 F,Frho,ri,rj,rn,mu,rho,test_ratio
      real*8 max_mu_over_rho,Fcrit,Frho_crit
      integer i,j,N

      rho_L=1.0d0
      rho_V=0.001d0
      mu_L=0.01d0
      mu_V=0.0001d0
      N=100
      max_mu_over_rho=0.0d0
      Fcrit=0.0d0
      Frho_crit=0.0d0

      do i=0,N
       ri=i
       rn=N
       F=ri/rn
       mu=mu_L*mu_V/((1.0d0-F)*mu_L+F*mu_V)
       do j=0,N
        rj=j
        Frho=rj/rn
        if (F.lt.0.5d0) then
         if (Frho.gt.0.5d0+F) then
          Frho=0.5d0+F
         endif
        endif
        if (F.gt.0.5d0) then
         if (Frho.lt.F-0.50d0) then
          Frho=F-0.5d0
         endif
        endif
        rho=Frho*rho_L+(1.0d0-Frho)*rho_V
        test_ratio=mu/rho
        if (test_ratio.gt.max_mu_over_rho) then
         max_mu_over_rho=test_ratio
         Fcrit=F
         Frho_crit=Frho
        endif
       enddo
      enddo

      print *,"max_mu_over_rho: ",max_mu_over_rho
      print *,"Fcrit: ",Fcrit
      print *,"Frho_crit: ",Frho_crit
      print *,"mu_L/rho_L=",mu_L/rho_L
      print *,"mu_V/rho_V=",mu_V/rho_V

      end

