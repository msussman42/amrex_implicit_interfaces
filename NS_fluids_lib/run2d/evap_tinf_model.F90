      module evap_check
      IMPLICIT NONE

       ! probtype==0 => Borodulin test
       ! probtype==1 => Villegas et al test
      integer, PARAMETER :: probtype = 1

      integer :: find_TINF_from_TGAMMA
      real*8 :: radblob
      real*8 :: den_L
      real*8 :: den_G
      real*8 :: C_pG
      real*8 :: k_G
      real*8 :: L_V
      real*8 :: D_G
      real*8 :: WV_global
      real*8 :: WA_global
      real*8 :: R_global
      real*8 :: T_sat_global
      real*8 :: lambda
      real*8 :: Le
      real*8 :: T_inf_global
      real*8 :: Y_inf_global
      real*8 :: T_gamma
      real*8 :: X_gamma
      real*8 :: Y_gamma
      real*8 :: B_M
      real*8 :: D_not
      real*8 :: Sh
     

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


      subroutine drop_analytical_solution(time,x,D_gamma,T,Y,VEL, &
       VEL_I,LS_VAP)
      IMPLICIT NONE

      real*8, intent(in) :: time
      real*8, intent(in) :: x
      real*8, intent(out) :: D_gamma
      real*8, intent(out) :: T
      real*8, intent(out) :: Y
      real*8, intent(out) :: LS_VAP
      real*8, intent(out) :: VEL
      real*8, intent(out) :: VEL_I
      real*8 :: rr,mdot,vel_r,my_pi,expand_factor
     
      my_pi=4.0d0*atan(1.0d0)
 
      rr=x
      
      D_Gamma=D_not**2-8.0d0*den_G*D_G*log(1.0d0+B_M)*time/den_L
      if ((D_Gamma.le.D_not**2).and.(D_Gamma.ge.0.0d0)) then
       D_Gamma=sqrt(D_gamma)
      else
       print *,"D_gamma invalid"
       print *,"D_gamma= ",D_gamma
       print *,"D_not= ",D_not
       print *,"den_G= ",den_G
       print *,"D_G= ",D_G
       print *,"B_M= ",B_M
       print *,"den_L=",den_L
       print *,"time=",time
       stop
      endif
      mdot=my_pi*D_gamma*den_G*D_G*Sh*B_M
      
      LS_VAP=rr-0.5d0*D_gamma
     
      VEL_I=mdot/(my_pi*D_gamma*D_gamma*den_G)
      expand_factor=1.0d0/den_G-1.0d0/den_L
      VEL_I=-VEL_I/(expand_factor*den_L)
 
      if (LS_VAP.le.0.0d0) then
       VEL=0.0d0
       Y=Y_Gamma
       T=T_Gamma
      else if (LS_VAP.ge.0.0d0) then
       vel_r = mdot/(4.0*my_pi*rr*rr*den_G)
       VEL=vel_r
       Y=1.0d0+(Y_inf_global-1.0d0)*exp(-mdot/(4.0d0*my_pi*den_G*D_G*rr))
       if ((Y.ge.0.0d0).and.(Y.lt.1.0d0)) then
        ! do nothing
       else
        print *,"Y invalid"
        stop
       endif
       T=T_Gamma-L_V/C_pG+(T_inf_global-T_Gamma+L_V/C_pG)* &
               exp(-mdot*C_pG/(4.0d0*my_pi*k_G*rr))
       if (T.gt.0.0d0) then
        ! do nothing
       else
        print *,"T invalid"
        stop
       endif
      else
       print *,"rr invalid"
       stop
      endif
      
      return
      end subroutine drop_analytical_solution


      subroutine ff(T_gamma_parm,f_out)
      IMPLICIT NONE

      real*8, intent(in) :: T_gamma_parm
      real*8, intent(out) :: f_out
      real*8 :: X_gamma_loc
      real*8 :: Y_gamma_loc

      call X_from_Tgamma(X_gamma_loc,T_gamma_parm,T_sat_global,L_V, &
        R_global,WV_global)
      call massfrac_from_volfrac(X_gamma_loc,Y_gamma_loc, &
        WA_global,WV_global)
      f_out=((Y_inf_global-1.0d0)/(Y_gamma_loc-1.0d0))**Le - 1.0d0 + &
       (T_gamma_parm-T_inf_global)*C_pG/L_V

      end subroutine ff

      end module evap_check

      program main
      use evap_check
      IMPLICIT NONE
      real*8 :: f_out
      real*8 :: aa,bb,cc,fa,fb,fc
      integer :: iter
      real*8 TSTART,TSTOP,cur_time,dt
      real*8 cur_x,T,Y,VEL,VEL_I,LS,D_gamma
      integer nsteps,istep
      integer outer_iter,max_outer_iter


      if (probtype.eq.0) then ! Borodulin et al figure 8, row 3
!      find_TINF_from_TGAMMA=1
       find_TINF_from_TGAMMA=0
       radblob = 0.05d0
       cur_x=4.0d0*radblob
       den_L = 1.0d0
       den_G = 0.001d0
       C_pG = 1.0d+7
       k_G = 0.024d+5
!      L_V = 2.26d+10  
       L_V = 2.1d+10  
       D_G = 0.1d0
       WV_global = 18.02d0
       WA_global = 28.9d0
       R_global = 8.31446261815324d+7
       T_sat_global=373.15d0
       T_inf_global = 300.5d0
       Y_inf_global=7.1d-3
       T_gamma=300.5  
       cc=0.0d0
       TSTOP=1000.0d0
      else if (probtype.eq.1) then
       find_TINF_from_TGAMMA=0
       radblob = 0.005d0
       cur_x=4.0d0*radblob
       den_L = 0.7d0
       den_G = 0.001d0
       C_pG = 1.0d+7
       k_G = 0.052d+5
       L_V = 5.18D+9
       D_G = 0.52d0
       WV_global = 58.0d0
       WA_global = 29.0d0
       R_global = 8.31446261815324d+7
       T_sat_global=329.0d0
       T_inf_global = 700.0d0
       Y_inf_global=0.0d0
       T_gamma=300.5d0  ! placeholder
       cc=0.0d0
       TSTOP=0.028d0
      else
       print *,"probtype invalid"
       stop
      endif

      lambda=k_G/(den_G*C_pG)
      Le=D_G*den_G*C_pG/k_G

      outer_iter=0

      if (find_TINF_from_TGAMMA.eq.1) then
       max_outer_iter=100000
      else
       max_outer_iter=1
      endif

      do while ((cc.lt.T_gamma).and.(outer_iter.lt.max_outer_iter))

      outer_iter=outer_iter+1
      aa=100.0d0
      bb=T_sat_global
      call ff(aa,fa)
      call ff(bb,fb)
      if (fa*fb.lt.0.0d0) then
       do iter=1,100
        cc=0.5d0*(aa+bb)
        call ff(cc,fc)
        if (fa*fc.gt.0.0d0) then
         aa=cc
        else if (fb*fc.gt.0.0d0) then
         bb=cc
        else if (fa.eq.0.0d0) then
         bb=cc
        else if (fb.eq.0.0d0) then
         aa=cc
        else if (fc.eq.0.0d0) then
         aa=cc
        else
         print *,"fa or fb or fc bust"
         stop
        endif
       enddo ! iter=1..100
      else
       print *,"fa or fb bust"
       stop
      endif
      print *,"den_L,den_G ",den_L,den_G
      print *,"C_pG,k_G,lambda ",C_pG,k_G,lambda
      print *,"T_inf_global,T_sat_global,L_V ", &
       T_inf_global,T_sat_global,L_V
      print *,"cc= ",cc
      print *,"D_G,Y_inf_global ",D_G,Y_inf_global
      print *,"WV_global,WA_global,Le ",WV_global,WA_global,Le
      print *,"R_global ",R_global

      if (find_TINF_from_TGAMMA.eq.1) then
       T_inf_global=T_inf_global+0.01d0
      endif

      enddo

      T_gamma=cc

      call X_from_Tgamma(X_gamma,T_gamma,T_sat_global,L_V, &
        R_global,WV_global)
      call massfrac_from_volfrac(X_gamma,Y_gamma,WA_global,WV_global)

      B_M=(Y_gamma-Y_inf_global)/(1.0d0-Y_gamma)
      D_not=2.0d0*radblob
      if (B_M.gt.0.0d0) then
       Sh=2.0d0*log(1.0d0+B_M)/B_M
      else
       print *,"B_M must be positive"
       stop
      endif
      print *,"INIT_DROP_IN_SHEAR_MODULE T_gamma,Y_gamma ", &
        T_gamma,Y_gamma

      TSTART=0.0d0
      cur_time=TSTART
      nsteps=8000
      dt=(TSTOP-TSTART)/nsteps
      do istep=1,nsteps
       call drop_analytical_solution(cur_time,cur_x,D_gamma,T,Y, &
         VEL,VEL_I,LS)
       print *,cur_time," ",0.5d0*D_gamma," ",(D_gamma/(2.0d0*radblob))**2," ",T," ",Y, &
        " ",VEL," ",VEL_I
       cur_time=cur_time+dt
      enddo

      return
      end

