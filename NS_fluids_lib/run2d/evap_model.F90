      module evap_check
      IMPLICIT NONE
      CONTAINS
      
      subroutine Pgamma_Clausius_Clapyron(Pgamma,PSAT,Tgamma,TSAT,L,R,WV)
      IMPLICIT NONE

      real*8, intent(in) :: PSAT,Tgamma,TSAT,L,R,WV
      real*8, intent(out) :: Pgamma

      Pgamma=PSAT*exp(-(L*WV/R)*(1.0d0/Tgamma-1.0d0/TSAT))
      return
      end subroutine Pgamma_Clausius_Clapyron

      subroutine X_from_Tgamma(X,Tgamma,TSAT, &
        L,R,WV)
      IMPLICIT NONE

      real*8, intent(in) :: Tgamma,TSAT,L,R,WV
      real*8, intent(out) :: X

      X=exp(-(L*WV/R)*(1.0d0/Tgamma-1.0d0/TSAT))

      return
      end subroutine X_from_Tgamma


      subroutine Tgamma_from_TSAT_and_X(Tgamma,TSAT, &
       X,L,R,WV,Tgamma_min,Tgamma_max)
      IMPLICIT NONE

      real*8, intent(in) :: TSAT,X,L,R,WV,Tgamma_min,Tgamma_max
      real*8, intent(out) :: Tgamma

      if (X.eq.1.0d0) then
       Tgamma=TSAT
      else if ((X.gt.0.0d0).and.(X.lt.1.0d0)) then
       Tgamma=-log(X)*R/(L*WV)
       Tgamma=Tgamma+1.0d0/TSAT
       if (Tgamma.gt.0.0d0) then
        Tgamma=1.0d0/Tgamma
       else
        print *,"Tgamma invalid in Tgamma_from_TSAT_and_X"
        stop
       endif
      else
       print *,"X invalid in Tgamma_from_TSAT_and_X"
       stop
      endif
 
      return
      end subroutine Tgamma_from_TSAT_and_X


      subroutine massfrac_from_volfrac(X,Y,WA,WV)
      IMPLICIT NONE

      real*8, intent(in) :: X,WA,WV
      real*8, intent(out) :: Y

      Y=WV*X/(WA*(1.0d0-X)+WV*X)
      return
      end subroutine massfrac_from_volfrac

      end module evap_check

      program main
      use evap_check
      IMPLICIT NONE
      real*8 :: T_gamma_parm
      real*8 :: f_out
      real*8 :: X_gamma_loc,Y_gamma_loc
      real*8 :: local_R
      real*8 :: den_L,den_G,C_pG,k_G
      real*8 :: lambda,T_inf,T_sat,L_V,D_G,Y_inf,WV,WA,Le

      den_L = 0.7d0
      den_G = 0.001d0
      C_pG = 1.0d+7
      k_G = 0.052d+5
      lambda=k_G/(den_G*C_pG)
      T_inf = 700.0d0
      T_sat = 329.0d0
      L_V = 5.18d+9
      D_G = 0.52d0
      Y_inf=0.0d0
      WV=58.0d0
      WA=29.0d0
      Le=D_G*den_G*C_pG/k_G
      local_R=8.31446261815324d+7
      T_gamma_parm=294.94
      call X_from_Tgamma(X_gamma_loc,T_gamma_parm,T_sat,L_V, &
        local_R,WV)
      call massfrac_from_volfrac(X_gamma_loc,Y_gamma_loc,WA,WV)
      f_out=((Y_inf-1.0d0)/(Y_gamma_loc-1.0d0))**Le - 1.0d0 + &
       (T_gamma_parm-T_inf)*C_pG/L_V
      print *,"Y_gamma_loc,T_gamma_parm,f_out ", &
              Y_gamma_loc,T_gamma_parm,f_out
      print *,"den_L,den_G ",den_L,den_G
      print *,"C_pG,k_G,lambda ",C_pG,k_G,lambda
      print *,"T_inf,T_sat,L_V ",T_inf,T_sat,L_V
      print *,"D_G,Y_inf ",D_G,Y_inf
      print *,"WV,WA,Le ",WV,WA,Le
      print *,"local_R ",local_R

      return
      end

