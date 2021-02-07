      module evap_check
      IMPLICIT NONE

      real*8, PARAMETER :: den_L = 1.0d0
      real*8, PARAMETER :: den_G = 0.001d0
      real*8, PARAMETER :: C_pG = 1.0d+7
      real*8, PARAMETER :: k_G = 0.024d+5
      real*8, PARAMETER :: L_V = 2.26d+10
      real*8, PARAMETER :: D_G = 0.1d0
      real*8, PARAMETER :: WV_global = 18.02d0
      real*8, PARAMETER :: WA_global = 28.9d0
      real*8, PARAMETER :: R_global = 8.31446261815324d+7
      real*8, PARAMETER :: T_sat_global=373.15d0
      real*8 :: lambda
      real*8 :: Le
      real*8 :: T_inf_global
      real*8 :: Y_inf_global
     

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

      lambda=k_G/(den_G*C_pG)
      Le=D_G*den_G*C_pG/k_G

      T_inf_global = 300.5d0
      Y_inf_global=7.1d-3
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

