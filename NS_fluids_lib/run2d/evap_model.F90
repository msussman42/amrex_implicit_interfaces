
      subroutine Pgamma_Clausius_Clapyron(Pgamma,PSAT,Tgamma,TSAT,L,R,WV)
      IMPLICIT NONE

      REAL_T, intent(in) :: PSAT,Tgamma,TSAT,L,R,WV
      REAL_T, intent(out) :: Pgamma

      Pgamma=PSAT*exp(-(L*WV/R)*(one/Tgamma-one/TSAT))
      return
      end subroutine Pgamma_Clausius_Clapyron

      subroutine X_from_Tgamma(X,Tgamma,TSAT, &
        L,R,WV)
      IMPLICIT NONE

      REAL_T, intent(in) :: Tgamma,TSAT,L,R,WV
      REAL_T, intent(out) :: X

      X=exp(-(L*WV/R)*(one/Tgamma-one/TSAT))

      return
      end subroutine X_from_Tgamma


      subroutine Tgamma_from_TSAT_and_X(Tgamma,TSAT, &
       X,L,R,WV,Tgamma_min,Tgamma_max)
      IMPLICIT NONE

      REAL_T, intent(in) :: TSAT,X,L,R,WV,Tgamma_min,Tgamma_max
      REAL_T, intent(out) :: Tgamma

      if (X.eq.one) then
       Tgamma=TSAT
      else if ((X.gt.zero).and.(X.lt.one)) then
       Tgamma=-log(X)*R/(L*WV)
       Tgamma=Tgamma+one/TSAT
       if (Tgamma.gt.zero) then
        Tgamma=one/Tgamma
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

      REAL_T, intent(in) :: X,WA,WV
      REAL_T, intent(out) :: Y

      Y=WV*X/(WA*(one-X)+WV*X)
      return
      end subroutine massfrac_from_volfrac


      subroutine f_mdot(T_gamma_parm,f_out)
      REAL_T, intent(in) :: T_gamma_parm
      REAL_T, intent(out) :: f_out
      REAL_T :: X_gamma_loc,Y_gamma_loc

      call X_from_Tgamma(X_gamma_loc,T_gamma_parm,T_sat,L_V, &
        fort_R_Palmore_Desjardins,WV)
      call massfrac_from_volfrac(X_gamma_loc,Y_gamma_loc,WA,WV)

      f_out=((Y_inf-one)/(Y_gamma_loc-one))**Le - one - &
       (T_gamma_parm-T_inf)*C_pG/L_V
      return
      end subroutine f_mdot

      program main
      IMPLICIT NONE

      den_L = 0.7d0
      den_G = 0.001d0
      C_pG = 1.0d+7
      k_G = 0.052d+5
      lambda=k_G/(den_G*C_pG)
      T_inf = 700.0d0
      T_sat = 329.0d0
      L_V = 5.18d+9
D_G = fort_speciesviscconst(2)
Y_inf=fort_speciesconst(2)
WV=fort_species_molar_mass(1)  !num_species components
WA=fort_molar_mass(2)
! T_inf C_pG / L = 300.5 K * 1e+7 (erg/(g K)) / ((2.26e+10) erg/g)
!   = 0.133
! e.g. Le=0.1 cm^2/s * 0.001 g/cm^3 * 4.1855e+7 erg/(g K) /
!         (0.024e+5 g cm/(s^3 K)) = 1.74
! erg=g cm^2/s^2
! cm^2/s  * g/cm^3  * g cm^2 / s^2 * (1/(g K)) * (s^3 K)/(g cm) = dimensionless
!
Le=D_G*den_G*C_pG/k_G

      return
      end

