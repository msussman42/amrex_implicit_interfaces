      module evap_check
      IMPLICIT NONE

       ! probtype==0 => Borodulin test
       ! probtype==1 => Villegas et al test
      integer, PARAMETER :: probtype = 0

       ! evap_model==0 => Villegas/Palmore,Desjardins model
       ! evap_model==1 => Kassemi model  P_ref=P_gamma/X_gamma
       ! evap_model==2 => same as evap_model==0, except that initial
       ! volume fraction is zero.
      integer, PARAMETER :: evap_model=0

      integer, PARAMETER :: nsteps=1000
      integer, PARAMETER :: num_intervals=32

      integer :: find_TINF_from_TGAMMA
      real*8 :: radblob
      real*8 :: den_L
      real*8 :: den_G
      real*8 :: C_pG
      real*8 :: gamma_G
      real*8 :: accommodation_coefficient
      real*8 :: k_G
      real*8 :: alpha_G
      real*8 :: L_V
      real*8 :: D_G
      real*8 :: WV_global
      real*8 :: WA_global
      real*8 :: R_global
      real*8 :: T_sat_global
      real*8 :: e_sat_global
      real*8 :: P_sat_global
      real*8 :: Pgamma_init_global
      real*8 :: X_gamma_init_global
      real*8 :: e_gamma_global
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
     
      subroutine EOS_material(den,e,P)
      IMPLICIT NONE

      real*8, intent(in) :: den,e
      real*8, intent(out) :: P

      P=(gamma_G-1.0d0)*den*e

      return
      end subroutine EOS_material

      
      subroutine INTERNAL_material(den,T,e)
      IMPLICIT NONE

      real*8, intent(in) :: den,T
      real*8, intent(out) :: e
      real*8 :: C_vG

      C_vG=C_pG/gamma_G
      e=C_vG*T

      return
      end subroutine INTERNAL_material

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


      ! Kassemi, Kartuzova, Hylton
      ! Cryogenics 89(2018) 1-15, equation (6)
      subroutine MDOT_Kassemi(sigma,MolarMassFluid,R,Pgamma, &
       Pvapor_probe, &
       Tgamma,Tvapor_probe,MDOT)
       IMPLICIT NONE

      real*8, intent(in) :: sigma
      real*8, intent(in) :: MolarMassFluid
      real*8, intent(in) :: R
      real*8, intent(in) :: Pgamma
      real*8, intent(in) :: Pvapor_probe
      real*8, intent(in) :: Tgamma
      real*8, intent(in) :: Tvapor_probe
      real*8, intent(out) :: MDOT
      real*8 :: my_pi
     
      my_pi=4.0d0*atan(1.0d0)

      if ((sigma.ge.0.0d0).and.(sigma.lt.2.0d0).and. &
          (MolarMassFluid.gt.0.0d0).and. &
          (R.gt.0.0d0).and. &
          (Pgamma.gt.0.0d0).and. &
          (Pvapor_probe.gt.0.0d0).and. &
          (Tgamma.gt.0.0d0).and. &
          (Tvapor_probe.gt.0.0d0)) then
       MDOT=(2.0d0*sigma/(2.0d0-sigma))* &
         sqrt(MolarMassFluid/(2.0d0*my_pi*R))* &
         (Pgamma/sqrt(Tgamma)-Pvapor_probe/sqrt(Tvapor_probe))
      
      else
       print *,"parameter problems in MDOT_Kassemi"
       stop
      endif
      
      return
      end subroutine MDOT_Kassemi


      subroutine mdot_from_Y_probe(T_gamma,T_probe,Y_probe,dx,mdotY)
      IMPLICIT NONE

      real*8, intent(in) :: T_gamma,T_probe,Y_probe,dx
      real*8, intent(out) :: mdotY
      real*8 :: internal_energy
      real*8 :: Pvapor_probe
      real*8 :: Pgamma
      real*8 :: local_X
      real*8 :: local_Y

      if (evap_model.eq.1) then
       call Pgamma_Clausius_Clapyron( &
        Pgamma, &
        P_sat_global, &
        T_gamma, &
        T_sat_global, &
        L_V, &
        R_global, &
        WV_global)

       call INTERNAL_material(den_G,T_probe,internal_energy)
       call EOS_material(den_G,internal_energy,Pvapor_probe)
       call MDOT_Kassemi( &
              accommodation_coefficient, &
              WV_global, &
              R_global, &
              Pgamma, &
              Pvapor_probe, &
              T_gamma, &
              T_probe, &
              mdotY)
      else if ((evap_model.eq.0).or.(evap_model.eq.2)) then
       call X_from_Tgamma(local_X,T_gamma,T_sat_global, &
        L_V,R_global,WV_global)
       call massfrac_from_volfrac(local_X,local_Y,WA_global,WV_global)
       mdotY=den_G*D_G*(local_Y-Y_probe)/((1.0d0-local_Y)*dx)
      else
       print *,"evap_model invalid"
       stop
      endif

      end subroutine mdot_from_Y_probe


      subroutine mdot_from_T_probe(T_gamma,T_probe,dx,mdotT)
      IMPLICIT NONE

      real*8, intent(in) :: T_gamma,T_probe,dx
      real*8, intent(out) :: mdotT

      mdotT=k_G*(T_probe-T_gamma)/(L_V*dx)

      end subroutine mdot_from_T_probe

      subroutine mdot_diff_func(T_gamma,T_probe,Y_probe,dx,mdot_diff)
      IMPLICIT NONE

      real*8, intent(in) :: T_gamma,T_probe,Y_probe,dx
      real*8, intent(out) :: mdot_diff
      real*8 :: mdotT,mdotY

      call mdot_from_T_probe(T_gamma,T_probe,dx,mdotT)
      call mdot_from_Y_probe(T_gamma,T_probe,Y_probe,dx,mdotY)

      mdot_diff=mdotT-mdotY

      end subroutine mdot_diff_func


      subroutine tridiag_solve(l,u,d,n,f,soln)
      IMPLICIT NONE

      integer n,i
      real*8 l(n),u(n),d(n),f(n),soln(n)
      real*8 ll(n),uu(n),dd(n),z(n)

      dd(1)=d(1)
      uu(1)=u(1)/dd(1)
      z(1)=f(1)/dd(1)
      do i=2,n-1
       ll(i)=l(i)
       dd(i)=d(i)-ll(i)*uu(i-1)
       uu(i)=u(i)/dd(i)
       z(i)=(f(i)-ll(i)*z(i-1))/dd(i)
      enddo
      ll(n)=l(n)
      dd(n)=d(n)-ll(n)*uu(n-1)
      z(n)=(f(n)-ll(n)*z(n-1))/dd(n)
      soln(n)=z(n)
      do i=n-1,1,-1
       soln(i)=z(i)-uu(i)*soln(i+1)
      enddo

      return
      end subroutine tridiag_solve



      end module evap_check



      program main
      use evap_check
      IMPLICIT NONE
      real*8 :: f_out
      real*8 :: aa,bb,cc,fa,fb,fc
      integer :: iter
      real*8 TSTART,TSTOP,cur_time,dt
      real*8 cur_x,T,Y,VEL,VEL_I,LS
      real*8 D_gamma
      real*8 R_gamma
      real*8 R_gamma_NEW
      real*8 R_gamma_OLD
      real*8 probhi_R_domain
      real*8 dx
      real*8 dx_new
      integer istep
      integer outer_iter,max_outer_iter
      integer igrid,igrid_old
      real*8 xpos
      real*8 xpos_new
      real*8 xpos_old
      real*8 xpos_mh
      real*8 xpos_ph
      real*8 vel_mh
      real*8 vel_ph
      real*8 advect_plus
      real*8 advect_minus
      real*8 diffuse_plus
      real*8 diffuse_minus

      real*8, allocatable :: TNEW(:)
      real*8, allocatable :: TOLD(:)
      real*8, allocatable :: TOLD_grid(:)
      real*8, allocatable :: YNEW(:)
      real*8, allocatable :: YOLD(:)
      real*8, allocatable :: YOLD_grid(:)
      real*8, allocatable :: DT_CROSSING(:)
      real*8, allocatable :: RHS(:)
      real*8, allocatable :: DIAG(:)
      real*8, allocatable :: LDIAG(:)
      real*8, allocatable :: UDIAG(:)
      real*8, allocatable :: SOLN(:)
      real*8 :: T_gamma_a
      real*8 :: T_gamma_b
      real*8 :: T_gamma_c
      real*8 :: Y_gamma_c
      real*8 :: X_gamma_c
      real*8 :: mdot_diff_a
      real*8 :: mdot_diff_b
      real*8 :: mdot_diff_c
      real*8 :: mdotT
      real*8 :: burnvel
      real*8 :: e_grid,P_grid,X_grid


      if (probtype.eq.0) then ! Borodulin et al figure 8, row 3
!      find_TINF_from_TGAMMA=1
       find_TINF_from_TGAMMA=0
       radblob = 0.05d0  ! cm
       cur_x=4.0d0*radblob
       den_L = 1.0d0  ! g/cm^3
       den_G = 0.001d0 ! g/cm^3
       C_pG = 1.0d+7  ! erg/(g K)
       gamma_G = 1.4d0
       accommodation_coefficient=1.0d0
       k_G = 0.024d+5 ! erg/(cm s K)
!      L_V = 2.26d+10  
       L_V = 2.1d+10  ! erg/g
       D_G = 0.1d0  ! cm^2/s
       WV_global = 18.02d0  ! g/mol
       WA_global = 28.9d0   ! g/mol
       R_global = 8.31446261815324d+7  ! ergs/(mol K)
       T_sat_global=373.15d0  ! K
       T_inf_global = 300.5d0 ! K
       Y_inf_global=7.1d-3  ! dimensionless
       T_gamma=300.5   ! K
       cc=0.0d0
       TSTOP=1000.0d0  ! seconds
      else if (probtype.eq.1) then
       find_TINF_from_TGAMMA=0
       radblob = 0.005d0
       cur_x=4.0d0*radblob
       den_L = 0.7d0
       den_G = 0.001d0
       C_pG = 1.0d+7
       gamma_G = 1.4d0
       accommodation_coefficient=1.0d0
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

      X_gamma_init_global=X_gamma

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
      dt=(TSTOP-TSTART)/nsteps
      print *,"ANALYTICAL TIME R_Gamma ARATIO T Y  VEL VEL_I "
      
      do istep=1,nsteps
       call drop_analytical_solution(cur_time,cur_x,D_gamma,T,Y, &
         VEL,VEL_I,LS)
       print *,cur_time," ",0.5d0*D_gamma," ", &
        (D_gamma/(2.0d0*radblob))**2," ",T," ",Y, &
        " ",VEL," ",VEL_I
       cur_time=cur_time+dt
      enddo

      R_gamma_NEW=radblob
      R_gamma_OLD=radblob
     
      probhi_R_domain=8.0d0*radblob+R_gamma_NEW
      dx=(probhi_R_domain-R_gamma_NEW)/num_intervals
      cur_time=0.0d0

      allocate(TNEW(0:num_intervals))
      allocate(TOLD(0:num_intervals))
      allocate(TOLD_grid(0:num_intervals))
      allocate(YNEW(0:num_intervals))
      allocate(YOLD(0:num_intervals))
      allocate(YOLD_grid(0:num_intervals))
      allocate(DT_CROSSING(0:num_intervals))
      allocate(RHS(1:num_intervals))
      allocate(DIAG(1:num_intervals))
      allocate(LDIAG(1:num_intervals))
      allocate(UDIAG(1:num_intervals))
      allocate(SOLN(1:num_intervals))

      TNEW(0)=T_gamma
      TOLD(0)=T_gamma
      YNEW(0)=Y_gamma
      YOLD(0)=Y_gamma
      do igrid=1,num_intervals
       xpos=igrid*dx+R_gamma_new
       call drop_analytical_solution(cur_time,xpos,D_gamma,T,Y, &
         VEL,VEL_I,LS)
       TNEW(igrid)=T
       TOLD(igrid)=T
       YNEW(igrid)=Y
       YOLD(igrid)=Y
      enddo

      call INTERNAL_material(den_G,T_gamma,e_gamma_global)
      call EOS_material(den_G,e_gamma_global,Pgamma_init_global)
      P_sat_global=Pgamma_init_global/X_gamma_init_global

      if (evap_model.eq.0) then ! Villegas model
       ! do nothing
      else if (evap_model.eq.1) then ! Kassemi model
       ! do nothing
      else if (evap_model.eq.2) then ! same as Villegas, except X=P/P_ref
       do igrid=0,num_intervals
        call INTERNAL_material(den_G,TNEW(igrid),e_grid)
        call EOS_material(den_G,e_grid,P_grid)
!       X_grid=P_grid/P_sat_global
        X_grid=0.0d0
        call massfrac_from_volfrac(X_grid,YNEW(igrid), &
         WA_global,WV_global)
        YOLD(igrid)=YNEW(igrid)
       enddo
       Y_inf_global=YNEW(num_intervals)
      else
       print *,"evap_model invalid"
       stop
      endif

      istep=0
      print *,"probtype=",probtype
      print *,"evap_model=",evap_model
      print *,"nsteps=",nsteps
      print *,"num_intervals=",num_intervals
      print *,"radblob=",radblob
      print *,"probhi_R_domain=",probhi_R_domain
      print *,"STEP TIME ARAT TGAMMA YGAMMA "
      print *,istep," ",cur_time," ", &
         (R_gamma_NEW/radblob)**2," ", &
         TNEW(0)," ",YNEW(0)

      do istep=1,nsteps

       T_gamma_a=100.0d0
       T_gamma_b=T_sat_global
       call mdot_diff_func(T_gamma_a,TNEW(1),YNEW(1),dx,mdot_diff_a)
       call mdot_diff_func(T_gamma_b,TNEW(1),YNEW(1),dx,mdot_diff_b)

       if (mdot_diff_a*mdot_diff_b.lt.0.0d0) then
        do iter=1,100
         cc=0.5d0*(T_gamma_a+T_gamma_b)
         call mdot_diff_func(cc,TNEW(1),YNEW(1),dx,mdot_diff_c)
         if (mdot_diff_a*mdot_diff_c.gt.0.0d0) then
          T_gamma_a=cc
         else if (mdot_diff_b*mdot_diff_c.gt.0.0d0) then
          T_gamma_b=cc
         else if (mdot_diff_a.eq.0.0d0) then
          T_gamma_b=cc
         else if (mdot_diff_b.eq.0.0d0) then
          T_gamma_a=cc
         else if (mdot_diff_c.eq.0.0d0) then
          T_gamma_a=cc
         else
          print *,"mdot_diff_a or mdot_diff_b or mdot_diff_c bust"
          stop
         endif
        enddo ! iter=1..100
       else
        print *,"mdot_diff_a or mdot_diff_b bust"
        stop
       endif
       T_gamma_c=cc
       call X_from_Tgamma(X_gamma_c,T_gamma_c,T_sat_global, &
        L_V,R_global,WV_global)
       call massfrac_from_volfrac(X_gamma_c,Y_gamma_c, &
        WA_global,WV_global)

       call mdot_from_T_probe(T_gamma_c,TNEW(1),dx,mdotT)
       burnvel=mdotT/den_L
       R_gamma_NEW=R_gamma_OLD-dt*burnvel

       dx_new=(probhi_R_domain-R_gamma_NEW)/num_intervals
       igrid_old=0
       xpos_old=R_gamma_OLD
       do igrid=1,num_intervals-1
        xpos_new=igrid*dx_new+R_gamma_NEW
        if (xpos_new.lt.R_gamma_OLD) then
         TOLD_grid(igrid)=T_gamma_c
         YOLD_grid(igrid)=Y_gamma_c
         DT_CROSSING(igrid)=dt*(R_gamma_NEW-xpos_new)/ &
                               (R_gamma_NEW-R_gamma_OLD)
        else
         DT_CROSSING(igrid)=dt
         do while (xpos_new.gt.xpos_old+dx)
          igrid_old=igrid_old+1
          xpos_old=xpos_old+dx
         enddo
         TOLD_grid(igrid)=TOLD(igrid_old)+ &
          (TOLD(igrid_old+1)-TOLD(igrid_old))* &
          (xpos_new-xpos_old)/dx
         YOLD_grid(igrid)=YOLD(igrid_old)+ &
          (YOLD(igrid_old+1)-YOLD(igrid_old))* &
          (xpos_new-xpos_old)/dx
        endif
       enddo ! igrid=1..num_intervals-1

       DT_CROSSING(0)=0.0d0
       DT_CROSSING(num_intervals)=dt
       TOLD_grid(0)=T_gamma_c
       TOLD_grid(num_intervals)=T_inf_global
       YOLD_grid(0)=Y_gamma_c
       YOLD_grid(num_intervals)=Y_inf_global

       alpha_G=k_G/(den_G*C_pG)
       do igrid=1,num_intervals-1
        xpos_new=igrid*dx_new+R_gamma_NEW
        xpos_mh=xpos_new-0.5d0*dx_new
        xpos_ph=xpos_new+0.5d0*dx_new

        RHS(igrid)=(dx_new**2)* &
          (xpos_new**2)*TOLD_grid(igrid)/DT_CROSSING(igrid)
        DIAG(igrid)=(dx_new**2)* &
          (xpos_new**2)/DT_CROSSING(igrid)
        LDIAG(igrid)=0.0d0
        UDIAG(igrid)=0.0d0
        vel_mh=((R_gamma_NEW/xpos_mh)**2)*mdotT/den_G
        vel_ph=((R_gamma_NEW/xpos_ph)**2)*mdotT/den_G
        advect_plus=vel_ph*(xpos_ph**2)*dx_new
        advect_minus=vel_mh*(xpos_mh**2)*dx_new
        diffuse_plus=alpha_G*(xpos_ph**2)
        diffuse_minus=alpha_G*(xpos_mh**2)

        DIAG(igrid)=DIAG(igrid)+advect_plus
        DIAG(igrid)=DIAG(igrid)+diffuse_plus
        DIAG(igrid)=DIAG(igrid)+diffuse_minus
        if (igrid.eq.1) then
         RHS(igrid)=RHS(igrid)+advect_minus*T_gamma_c
         RHS(igrid)=RHS(igrid)+diffuse_minus*T_gamma_c
         LDIAG(igrid)=0.0d0
         UDIAG(igrid)=UDIAG(igrid)-diffuse_plus
        else if (igrid.eq.num_intervals-1) then
         RHS(igrid)=RHS(igrid)+diffuse_plus*T_inf_global
         UDIAG(igrid)=0.0d0
         LDIAG(igrid)=LDIAG(igrid)-diffuse_minus
         LDIAG(igrid)=LDIAG(igrid)-advect_minus
        else
         UDIAG(igrid)=UDIAG(igrid)-diffuse_plus
         LDIAG(igrid)=LDIAG(igrid)-diffuse_minus
         LDIAG(igrid)=LDIAG(igrid)-advect_minus
        endif

       enddo !igrid=1,num_intervals-1

       call tridiag_solve(LDIAG,UDIAG,DIAG,num_intervals-1,RHS,SOLN)

       TNEW(0)=T_gamma_c
       TNEW(num_intervals)=T_inf_global 
       do igrid=1,num_intervals-1
        TNEW(igrid)=SOLN(igrid)
       enddo

       do igrid=1,num_intervals-1
        xpos_new=igrid*dx_new+R_gamma_NEW
        xpos_mh=xpos_new-0.5d0*dx_new
        xpos_ph=xpos_new+0.5d0*dx_new

        RHS(igrid)=(dx_new**2)* &
          (xpos_new**2)*YOLD_grid(igrid)/DT_CROSSING(igrid)
        DIAG(igrid)=(dx_new**2)* &
          (xpos_new**2)/DT_CROSSING(igrid)
        LDIAG(igrid)=0.0d0
        UDIAG(igrid)=0.0d0
        vel_mh=((R_gamma_NEW/xpos_mh)**2)*mdotT/den_G
        vel_ph=((R_gamma_NEW/xpos_ph)**2)*mdotT/den_G
        advect_plus=vel_ph*(xpos_ph**2)*dx_new
        advect_minus=vel_mh*(xpos_mh**2)*dx_new
        diffuse_plus=D_G*(xpos_ph**2)
        diffuse_minus=D_G*(xpos_mh**2)

        DIAG(igrid)=DIAG(igrid)+advect_plus
        DIAG(igrid)=DIAG(igrid)+diffuse_plus
        DIAG(igrid)=DIAG(igrid)+diffuse_minus
        if (igrid.eq.1) then
         RHS(igrid)=RHS(igrid)+advect_minus*Y_gamma_c
         RHS(igrid)=RHS(igrid)+diffuse_minus*Y_gamma_c
         LDIAG(igrid)=0.0d0
         UDIAG(igrid)=UDIAG(igrid)-diffuse_plus
        else if (igrid.eq.num_intervals-1) then
         RHS(igrid)=RHS(igrid)+diffuse_plus*Y_inf_global
         UDIAG(igrid)=0.0d0
         LDIAG(igrid)=LDIAG(igrid)-diffuse_minus
         LDIAG(igrid)=LDIAG(igrid)-advect_minus
        else
         UDIAG(igrid)=UDIAG(igrid)-diffuse_plus
         LDIAG(igrid)=LDIAG(igrid)-diffuse_minus
         LDIAG(igrid)=LDIAG(igrid)-advect_minus
        endif

       enddo !igrid=1,num_intervals-1

       call tridiag_solve(LDIAG,UDIAG,DIAG,num_intervals-1,RHS,SOLN)

       YNEW(0)=Y_gamma_c
       YNEW(num_intervals)=Y_inf_global 
       do igrid=1,num_intervals-1
        YNEW(igrid)=SOLN(igrid)
       enddo

       do igrid=0,num_intervals
        TOLD(igrid)=TNEW(igrid)
        YOLD(igrid)=YNEW(igrid)
       enddo

       dx=dx_new
       R_gamma_OLD=R_gamma_NEW
       cur_time=cur_time+dt
       print *,istep," ",cur_time," ", &
         (R_gamma_NEW/radblob)**2," ", &
         TNEW(0)," ",YNEW(0)

      enddo ! istep=1..nsteps

      deallocate(TNEW)
      deallocate(TOLD)
      deallocate(TOLD_grid)
      deallocate(YNEW)
      deallocate(YOLD)
      deallocate(YOLD_grid)
      deallocate(DT_CROSSING)
      deallocate(RHS)
      deallocate(DIAG)
      deallocate(LDIAG)
      deallocate(UDIAG)
      deallocate(SOLN)

      return
      end

