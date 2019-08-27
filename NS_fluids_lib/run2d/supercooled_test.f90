!!$ Written by Mehdi Vahab
!!$ Last edit: 05/21/2015
!!$ Test program to use the supercooled_exact_sol module.
!!$ Check supercooled_exact_sol.f90 for more details



!!$ Written by Mehdi Vahab
!!$ Last edit: 05/21/2015
!!$ References: 
!!$
!!$ 
!!$ Functions and subroutines to evaluate the exact solution of 
!!$ solidification front in 2D cylindrical coordinate system. 
!!$ The region rd<R consists of solid at its melting point T1, 
!!$ and the region rd>R of supercooled liquid whose temperature
!!$  v2 -> V < T1 as rd -> inf.
!!$
!!$ The position of the surface of separation is given by
!!$
!!$  (1)    R = 2 lm*( k2*t)^0.5,                     
!!$
!!$ where lm is the root of
!!$
!!$  (2)    lm^2 * exp(-lm^2) * Ei(-lm**2) + St=0,    
!!$
!!$ where and k2 is thermal diffusivity of the liquid, 
!!$ Ei is exponential integral function 
!!$
!!$  (3)    Ei(x) = - int_{-x}^{\inf} exp(-t) / t dt,
!!$
!!$ and St is the Stefan number
!!$
!!$  (4)    St = c2*(T1-V)/L.                         
!!$
!!$ c2 is the specific heat capacity of the liquid and 
!!$ L is the latent heat.
!!$ The temperature in liquid at radius rd and time t is given by
!!$
!!$  (5)    v2(rd,t) = V + (L * St / c2)* Ei(-rd^2/(4*k2*t)) / Ei(-lm^2).
!!$
!!$ [1]. Conduction of Heat in Solids, H. S. Carslaw and J. C. Jaeger,
!!$        Oxford, Second Edition, 1959
!!$ [2]. Numerical Analysis, R. L. Burden, J. D. Faires, Brooks/Cole, 
!!$        Ninth Edition, 2011
!!$ [3]. http://en.wikipedia.org/wiki/Exponential_integral
 
module supercooled_exact_sol
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)

contains

!!$****************************************************

  REAL(kind=dp) function g(x)
  REAL(kind=dp),intent(in):: x

  g=exp(-x)/x

  end function g
  
  REAL(kind=dp) function Ei(x)
!!$    Evaluate and returns the exponantial integral function 
!!$    using the a series sum [3]:
!!$
!!$    y = euler_mascheroni_constant + ln|x| + sum_{k=1}^{inf} x^k/(k*k!)
!!$
!!$    input(s):
!!$    x     :    function argument, real number
!!$    output(s):
!!$    y     :    Ei(x), real number

    REAL(kind=dp),intent(in):: x
    REAL(kind=dp):: y, inc, xtok, kf
    REAL(kind=dp):: xlo,xhi,dx,a,b,slope
    REAL(kind=dp), parameter::em_const=0.5772156649015329_dp
    REAL(kind=dp), parameter::tol=1E-15_dp
    INTEGER:: k,N


    if (1.eq.1) then

     !! initiate result
     y = em_const + log(abs(x)) 
     !! counter
     k = 1
     !! k factorial
     kf = 1_dp
     !! x^k
     xtok = x
     !! sum incremental term
     inc = xtok/(k*kf);
     y = y + inc;

     !! loop, evaluate terms and add them up
     do
        !! write(*,"(ES30.15)") y
        k = k+1
        kf = kf * k
        xtok = xtok * x
        inc = xtok/(k*kf)
        y = y + inc
        if((abs(inc).lt.tol).or.(k.gt.100)) then
           exit
        end if
     end do

     !! check for convergence
     if(k.eq.101) then
        print *,"Ei(x) function did not comverge"
        print *,"x is", x
        stop
     end if

    else

     if (x.ge.0.0) then
      print *,"x must be negative in ei(x)"
      stop
     endif
     N=1000
     y=0.0
     if (abs(x).lt.1.0) then
      xhi=1.0
      xlo=abs(x)
      dx=(xhi-xlo)/N
      do k=0,N-1
       a=xlo+k*dx
       b=a+dx
       slope=(exp(-b)-exp(-a))/dx
       y=y+log(b/a)*(exp(-a)-slope*a)+slope*dx
      enddo 
     endif
     if (abs(x).lt.1.0) then
      xlo=1.0
     else
      xlo=abs(x)
     endif
     xhi=18.0*log(10.0)
     if (xlo.lt.xhi) then
      dx=(xhi-xlo)/N
      do k=0,N-1
       a=xlo+k*dx
       b=a+dx
       y=y+(dx/6.0)*(g(a)+4.0*g(a+0.5*dx)+g(b))
      enddo
     endif
     y=-y 
    endif

    Ei = y

  end function Ei

!!$****************************************************
  REAL(kind=dp) function f_lambda(lm,St)
!!$    Evaluates the function in Eq.(2)
!!$    input(s):
!!$    lm       :  \lambda, function argument, real number
!!$    St       :  Stefan number, function argument, real number
!!$    output(s):
!!$    f_lambda :  f(\lambda,St), function value, real number

    REAL(kind=dp),INTENT(in):: lm
    REAL(kind=dp),INTENT(in):: St
    
!    f_lambda = (lm**2.0) * exp(-lm**2.0) * Ei(-lm**2.0) + St;
    f_lambda = (lm**2.0) * exp(lm**2.0) * Ei(-lm**2.0) + St;

  end function f_lambda

!!$****************************************************
  subroutine find_lambda(lambda,St)
!!$    Find the smalest positive root of the function in Eq.(2)
!!$    input(s):
!!$    St       :  Stefan number, function argument, real number
!!$    output(s):
!!$    lm       :  \lambda, smallet positive root, real number


    REAL(kind=dp),INTENT(in):: St
    REAL(kind=dp),INTENT(out):: lambda

    REAL(kind=dp), parameter:: dx = 1E-3_dp     !! Marching step
    REAL(kind=dp), parameter:: tol = 1E-14_dp   !! tolerance of root finding 
    INTEGER, parameter:: max_itr = 50           !! max ietration number
    
    REAL(kind=dp) :: a, b, fa, fb, p, fp
    
    INTEGER :: k

    !! initialize the search region
    a = 1E-10_dp
    fa = f_lambda(a,St)
    b = a
    !! march forward to have a zero crossing in the search domain
    do
       b = b + dx
       fb = f_lambda(b,St)
       if((fa*fb).lt.0.0_dp) then
          exit
       end if
    end do

!!$    print *,"a= ", a
!!$    print *,"b= ", b

    !!b = 0.12_dp
    
    fb = f_lambda(b,St)

    if ((fa*fb).gt.0.0_dp) then
       print *,"Chosen range may not have a zero crossing!"
       stop
    end if

    
    !! Bisection root finding method based on [2]
    do k = 1, max_itr
!!$       print *,k
!!$       print *,p
       p = a + (b-a)/2.0_dp
       fp = f_lambda(p,St)

       if ((fp.eq.0.0_dp).or.((b-a)/2).lt.tol) then
          lambda = p
          return
       end if

       if((fp*fa).gt.0.0_dp) then
          a = p
          fa = fp
       else
          b = p
       end if
    end do
    
    
    print *,"Reached max iterations &
         but not converged to the requested tolerence"
    stop
  
  end subroutine find_lambda


  subroutine find_lambda2(lambda,St)

    REAL(kind=dp),INTENT(in):: St
    REAL(kind=dp),INTENT(out):: lambda

    REAL(kind=dp), parameter:: dx = 1E-3_dp     !! Marching step
    REAL(kind=dp), parameter:: tol = 1E-14_dp   !! tolerance of root finding 
    INTEGER, parameter:: max_itr = 50           !! max ietration number
    
    REAL(kind=dp) :: a, b, fa, fb, p, fp,lambda1
    
    INTEGER :: k

    call find_lambda(lambda1,St)

    !! initialize the search region
    a = lambda1+1E-10_dp
    fa = f_lambda(a,St)
    b = a
    !! march forward to have a zero crossing in the search domain
    do
       b = b + dx
       fb = f_lambda(b,St)
       if((fa*fb).lt.0.0_dp) then
          exit
       end if
    end do

!!$    print *,"a= ", a
!!$    print *,"b= ", b

    !!b = 0.12_dp
    
    fb = f_lambda(b,St)

    if ((fa*fb).gt.0.0_dp) then
       print *,"2nd: Chosen range may not have a zero crossing!"
       stop
    end if

    
    !! Bisection root finding method based on [2]
    do k = 1, max_itr
!!$       print *,k
!!$       print *,p
       p = a + (b-a)/2.0_dp
       fp = f_lambda(p,St)

       if ((fp.eq.0.0_dp).or.((b-a)/2).lt.tol) then
          lambda = p
          return
       end if

       if((fp*fa).gt.0.0_dp) then
          a = p
          fa = fp
       else
          b = p
       end if
    end do
    
    
    print *,"Reached max iterations &
         but not converged to the requested tolerence"
    stop
  
  end subroutine find_lambda2


!!$****************************************************
subroutine solidification_front_radius(lm, k2, t, R)
!!$    Evaluates the position of the solidification front 
!!$    based of Eq.(1)
!!$    input(s):
!!$    lm       :  \lambda, root of Eq.(2), real number
!!$    k2       :  thermal diffusivity of the liquid, real number
!!$    t        :  time, real number
!!$    output(s):
!!$    R        :  radius of solidification front, real number

  REAL(kind=dp),INTENT(in):: lm, k2, t
  REAL(kind=dp),INTENT(out):: R

  R = 2.0 * lm * sqrt(k2 * t)
end subroutine solidification_front_radius


subroutine solidification_front_time(lm, k2, t, R)

  REAL(kind=dp),INTENT(in):: lm, k2
  REAL(kind=dp),INTENT(in):: R
  REAL(kind=dp),INTENT(out):: t

  t=(R/(2.0*lm))**2.0/k2

end subroutine solidification_front_time



!!$****************************************************
subroutine liquid_temperature(lm, V, L, c2, St, rd, t, k2, v2)
!!$    Evaluates the temperature in liquid at given time and radius 
!!$    based of Eq.(5)
!!$    input(s):
!!$    lm       :  \lambda, root of Eq.(2), real number
!!$    V        :  temperature of liquid at infinity, real number
!!$    L        :  latent heat, real number
!!$    c2       :  specific heat capacity of the liquid , real number
!!$    St       :  Stefan number, real number
!!$    r        :  radius, function argument, real number
!!$    t        :  time, function argument, real number
!!$    k2       :  thermal diffusivity of the liquid, real number
!!$    output(s):
!!$    v2       :  temperature of liquid phase, real number
!!$
!!$    NOTE: This routine does not check that the given radius is
!!$          in the liquid phase. That is, be sure that rd > R, where 
!!$          R is the position of the solidification front!


  REAL(kind=dp),INTENT(in):: lm, V, L, c2, St, rd, t, k2
  REAL(kind=dp),INTENT(out):: v2


  v2 = V + (L * St / c2)* Ei(-rd*rd/(4.0*k2*t)) / Ei(-lm*lm)
end subroutine liquid_temperature
!!$****************************************************

end module supercooled_exact_sol

program supercooled_test
use supercooled_exact_sol
  IMPLICIT NONE

  REAL(kind=dp):: St, lm
  !! specific heat capacity of liquid water (erg . g^-1 . K^-1)
  REAL(kind=dp):: c2 = 4.1841E7_dp 
  !! thermal diffusivity of liquid water (cm^2 . s^-1)
  REAL(kind=dp):: k2 = 1.424440142444014e-03_dp
  !! latent heat of evaporation of fusion of water
  REAL(kind=dp):: L = 3.34E9_dp
  !! far-field temparature of supercooled liquid water
  REAL(kind=dp):: V = 263.15_dp
  !! melting point of water
  REAL(kind=dp):: T1 = 273.15_dp
  INTEGER::simple_parm=1
  !! time
  REAL(kind=dp):: t
  !! solidification front position at a given time
  REAL(kind=dp):: R
  !! radius
  REAL(kind=dp):: rd
  !! temperature in liquid in giver time and radius
  REAL(kind=dp):: v2
  !! time step size for the R vs t figure
  REAL(kind=dp):: dt

  REAL(kind=dp):: tm, max_rd, drd,radblob,time_start,dx,xhi,xlo
  INTEGER:: k,NPLOT
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! Usage examples !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (simple_parm.eq.1) then
   c2=1.0
   T1=273.0
   V=272.0 
   k2=1.0
   L=8.0
!   L=7.472
  endif

  !! find lambda for a given Stefan number
  St = c2*(T1-V)/L
  print *,"St=", st
  call find_lambda(lm, st)
  print *,"lambda=", lm

  radblob=0.05

  !! radius of the solidification front at a givern time
  t = 0.5_dp
  call solidification_front_radius(lm, k2, t, R)
  print *,"R=", R
  call solidification_front_time(lm, k2, time_start, radblob)

  !! temperature of liquid phase at a given time and radius
  t = 0.5_dp
  rd = radblob
  call liquid_temperature(lm, V, L, c2, st, rd, t, k2, v2)
  print *,"v2(t=",t,",r=",rd,")=", v2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  print *,"creating Ei_vs_x.dat"
  open (UNIT=10, FILE='Ei_vs_x.dat', STATUS='REPLACE', ACTION='WRITE')
  xhi=-0.1
  xlo=-5.0
  NPLOT=100
  dx=(xhi-xlo)/NPLOT
  do k = 0, NPLOT
   tm=xlo+dx*k
   write(10,*) tm, Ei(tm)
  end do
  close(UNIT=10)

  print *,"creating flm_vs_lm.dat"
  open (UNIT=10, FILE='flm_vs_lm.dat', STATUS='REPLACE', ACTION='WRITE')
  xhi=1.0
  xlo=1.0e-10
  NPLOT=100
  dx=(xhi-xlo)/NPLOT
  do k = 0, NPLOT
   tm=xlo+dx*k
   write(10,*) tm, f_lambda(tm,St)
  end do
  close(UNIT=10)
 
  !! store  R vs t data in file
  !! plot of solidification radius vs time
  NPLOT=500
  t = 0.2
  dt=(t/NPLOT)

  print *,"creating R_vs_t.dat"

  open (UNIT=10, FILE='R_vs_t.dat', STATUS='REPLACE', ACTION='WRITE')
  do k = 0, NPLOT
     tm = dt * k+time_start
     call solidification_front_radius(lm, k2, tm, R)
     write(10,*) tm-time_start, R
  end do
  close(UNIT=10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! plot of temperature vs radius profile at a given time 
  NPLOT=100
  t = 0.5_dp
  drd=t/NPLOT;
  max_rd = 0.1
  k = floor(max_rd/drd)

  !! radius of the solidification front at a givern time
  call solidification_front_radius(lm, k2, t, R)

  print *,"creating temp_vs_r.dat"

  open (UNIT=11, FILE='temp_vs_r.dat', STATUS='REPLACE', ACTION='WRITE')
  do k = 0, NPLOT
     rd = drd * k
     if(rd.le.R) then
        v2 = T1
     else
        call liquid_temperature(lm, V, L, c2, st, rd, t, k2, v2)
     end if
     write(11,*) rd, v2
  end do
  close(UNIT=11)
  
end program supercooled_test

