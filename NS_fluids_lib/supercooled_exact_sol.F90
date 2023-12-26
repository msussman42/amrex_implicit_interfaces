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
!!$  (2)    lm^2 * exp(lm^2) * Ei(-lm**2) + St=0,    
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
use amrex_fort_module, only : amrex_real
  IMPLICIT NONE

  REAL(amrex_real) :: supercooled_lm
  REAL(amrex_real) :: supercooled_temp_infinity
  REAL(amrex_real) :: supercooled_L
  REAL(amrex_real) :: supercooled_specific_heat
  REAL(amrex_real) :: supercooled_stefan_number
  REAL(amrex_real) :: supercooled_thermal_diff

contains

!!$****************************************************

  REAL(amrex_real) function Ei(x)
!!$    Evaluate and returns the exponantial integral function 
!!$    using the a series sum [3]:
!!$
!!$    y = euler_mascheroni_constant + ln|x| + sum_{k=1}^{inf} x^k/(k*k!)
!!$
!!$    input(s):
!!$    x     :    function argument, real number
!!$    output(s):
!!$    y     :    Ei(x), real number

    REAL(amrex_real),INTENT(in):: x
    REAL(amrex_real):: y, inc, xtok, kf
    REAL(amrex_real), parameter::em_const=0.5772156649015329D0
    REAL(amrex_real), parameter::tol=1.0D-7
    integer:: k


    !! initiate result
    y = em_const + log(abs(x)) 
    !! counter
    k = 1
    !! k factorial
    kf = 1.0D0
    !! x^k
    xtok = x
    !! sum incremental term
    inc = xtok/(k*kf)
    y = y + inc

    if (x.lt.-25.0D0) then
     y=0.0
    else 
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
    endif

    !! check for convergence
    if(k.eq.101) then
       print *,"Ei(x) function did not converge"
       print *,"x is", x
       print *,"returning Ei(x)=0"
       y=0.0
    end if

    Ei = y

  end function Ei

!!$****************************************************
  REAL(amrex_real) function f_lambda(lm,St)
!!$    Evaluates the function in Eq.(2)
!!$    input(s):
!!$    lm       :  \lambda, function argument, real number
!!$    St       :  Stefan number, function argument, real number
!!$    output(s):
!!$    f_lambda :  f(\lambda,St), function value, real number

    REAL(amrex_real),INTENT(in):: lm
    REAL(amrex_real),INTENT(in):: St
    
    f_lambda = lm**2 * exp(lm**2) * Ei(-lm**2) + St

  end function f_lambda
!!$****************************************************
  subroutine find_lambda(lambda,St)
!!$    Find the smalest positive root of the function in Eq.(2)
!!$    input(s):
!!$    St       :  Stefan number, function argument, real number
!!$    output(s):
!!$    lm       :  \lambda, smallet positive root, real number


    REAL(amrex_real),INTENT(in):: St
    REAL(amrex_real),INTENT(out):: lambda

    REAL(amrex_real), parameter:: dx = 1.0D-2     !! Marching step
    REAL(amrex_real), parameter:: tol = 1.0D-7   !! tolerance of root finding 
    integer, parameter:: max_itr = 50           !! max iteration number
    integer, parameter:: max_bracketing= 1000  !! max iteration number
    
    REAL(amrex_real) :: a, b, fa, fb, p, fp
    
    integer :: k
    integer :: nbracket

    !! initialize the search region
    a = 1.0D-7
    fa = f_lambda(a,St)
    b = a
    !! march forward to have a zero crossing in the search domain
    nbracket=0
    do
       b = b + dx
       fb = f_lambda(b,St)
       if((fa*fb).lt.0.0D0) then
          exit
       end if
       nbracket=nbracket+1
       if (nbracket.gt.max_bracketing) then
        exit
       endif
    end do

!!$    print *,"a= ", a
!!$    print *,"b= ", b

    !!b = 0.12D0
    
    fb = f_lambda(b,St)

    if ((fa*fb).gt.0.0D0) then
       print *,"Chosen range may not have a zero crossing!"
       print *,"find_lambda failed, setting lambda=0"
       lambda=0.0d0
       return
    end if

    
    !! Bisection root finding method based on [2]
    do k = 1, max_itr
!!$       print *,k
!!$       print *,p
       p = a + (b-a)/2.0D0
       fp = f_lambda(p,St)

       if ((fp.eq.0.0D0).or.((b-a)/2).lt.tol) then
          lambda = p
          return
       end if

       if((fp*fa).gt.0.0D0) then
          a = p
          fa = fp
       else
          b = p
       end if
    end do
    
    
    print *,"Reached max iterations &
         but not converged to the requested tolerence"
    print *,"find_lambda failed, setting lambda=0"
    lambda=0.0d0
  
  end subroutine find_lambda
!!$****************************************************


subroutine solidification_front_speed_driver(t, speed)
!!$    t        :  time, real number
!!$    output(s):
!!$    R        :  radius of solidification front, real number
  REAL(amrex_real),INTENT(in):: t
  REAL(amrex_real),INTENT(out):: speed

  if (supercooled_lm.ge.0.0d0) then

   if (supercooled_thermal_diff.ge.0.0d0) then

    if (t.gt.0.0d0) then
     speed=supercooled_lm*supercooled_thermal_diff/ &
       sqrt(supercooled_thermal_diff  * t)
    else
     print *,"t invalid"
     stop
    endif
   else
    print *,"supercooled_thermal_diff invalid"
    print *,"supercooled_thermal_diff=",supercooled_thermal_diff
    stop
   endif
  else
   print *,"supercooled_lm invalid"
   stop
  endif

end subroutine solidification_front_speed_driver

subroutine solidification_front_radius_driver(t, R)
!!$    t        :  time, real number
!!$    output(s):
!!$    R        :  radius of solidification front, real number
  REAL(amrex_real),INTENT(in):: t
  REAL(amrex_real),INTENT(out):: R

  if (supercooled_lm.ge.0.0d0) then

   if (supercooled_thermal_diff.ge.0.0d0) then

    if (t.ge.0.0d0) then
     R = 2.0d0 * supercooled_lm * sqrt(supercooled_thermal_diff  * t)
    else
     print *,"t invalid"
     stop
    endif
   else
    print *,"supercooled_thermal_diff invalid"
    print *,"supercooled_thermal_diff(2)=",supercooled_thermal_diff
    stop
   endif
  else
   print *,"supercooled_lm invalid"
   stop
  endif

end subroutine solidification_front_radius_driver

subroutine solidification_front_radius(lm, k2, t, R)
!!$    Evaluates the position of the solidification front 
!!$    based of Eq.(1)
!!$    input(s):
!!$    lm       :  \lambda, root of Eq.(2), real number
!!$    k2       :  thermal diffusivity of the liquid, real number
!!$    t        :  time, real number
!!$    output(s):
!!$    R        :  radius of solidification front, real number

  REAL(amrex_real),INTENT(in):: lm, k2, t
  REAL(amrex_real),INTENT(out):: R

  supercooled_lm=lm
  supercooled_thermal_diff=k2
  call solidification_front_radius_driver(t,R)

end subroutine solidification_front_radius


subroutine solidification_front_time_driver(t, R)

  REAL(amrex_real),INTENT(in):: R
  REAL(amrex_real),INTENT(out):: t

  if (supercooled_lm.gt.0.0d0) then
   if (supercooled_thermal_diff.gt.0.0d0) then
    t=(R/(2.0d0*supercooled_lm))**2/supercooled_thermal_diff
   else if (supercooled_thermal_diff.eq.0.0d0) then
    print *,"supercooled_thermal_diff cannot be 0.0"
    stop
   else
    print *,"supercooled_thermal_diff invalid"
    stop
   endif
  else 
   print *,"supercooled_lm invalid setting t=0.0"
   t=0.0d0
  endif

end subroutine solidification_front_time_driver

subroutine solidification_front_time(lm, k2, t, R)

  REAL(amrex_real),INTENT(in):: lm, k2
  REAL(amrex_real),INTENT(in):: R
  REAL(amrex_real),INTENT(out):: t

  supercooled_lm=lm
  supercooled_thermal_diff=k2

  call solidification_front_time_driver(t,R)

end subroutine solidification_front_time


subroutine liquid_temperature_driver(rd, t, v2)
!!$    rd       :  radius, function argument, real number
!!$    t        :  time, function argument, real number
!!$    output(s):
!!$    v2       :  temperature of liquid phase, real number
!!$

  REAL(amrex_real),INTENT(in):: rd, t
  REAL(amrex_real),INTENT(out):: v2
  REAL(amrex_real) :: Ei_parm

  if (supercooled_L.ge.0.0d0) then

   if (supercooled_temp_infinity.gt.0.0d0) then

    if (t.gt.0.0d0) then

     if (supercooled_specific_heat.gt.0.0d0) then

      if (supercooled_stefan_number.gt.0.0d0) then

       if (supercooled_thermal_diff.gt.0.0d0) then

        if (supercooled_lm.ge.0.0d0) then  
         Ei_parm=-rd*rd/(4.0d0*supercooled_thermal_diff*t)

         v2 = supercooled_temp_infinity +  &
          (supercooled_L * supercooled_stefan_number / &
           supercooled_specific_heat)* &
           Ei(Ei_parm) / &
           Ei(-(supercooled_lm**2))

         if (v2.le.0.0d0) then
          print *,"supercooled_lm ",supercooled_lm
          print *,"supercooled_L ",supercooled_L
          print *,"supercooled_stefan_number ",supercooled_stefan_number
          print *,"supercooled_specific_heat ",supercooled_specific_heat
          print *,"supercooled_thermal_diff ",supercooled_thermal_diff
          print *,"rd ",rd
          print *,"t ",t
          Ei_parm=-rd*rd/(4.0d0*supercooled_thermal_diff*t)
          print *,"Ei numerator ", &
            Ei(Ei_parm)
          print *,"Ei denominator ", &
            Ei(-(supercooled_lm**2))
          stop
         endif
        else
         print *,"supercooled_lm invalid"
         stop
        endif
       else
        print *,"supercooled_thermal_diff invalid"
        stop
       endif
      else
       print *,"supercooled_stefan_number invalid"
       stop
      endif
     else
      print *,"supercooled_specific_heat invalid"
      stop
     endif
    else
     print *,"t invalid"
     stop
    endif
   else
    print *,"supercooled_temp_infinity invalid"
    stop
   endif

  else
   print *,"supercooled_L invalid"
   stop
  endif

end subroutine liquid_temperature_driver

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
!!$    rd       :  radius, function argument, real number
!!$    t        :  time, function argument, real number
!!$    k2       :  thermal diffusivity of the liquid, real number
!!$    output(s):
!!$    v2       :  temperature of liquid phase, real number
!!$
!!$    NOTE: This routine does not check that the given radius is
!!$          in the liquid phase. That is, be sure that rd > R, where 
!!$          R is the position of the solidification front!


  REAL(amrex_real),INTENT(in):: lm, V, L, c2, St, rd, t, k2
  REAL(amrex_real),INTENT(out):: v2

  supercooled_lm=lm
  supercooled_temp_infinity=V
  supercooled_L=L
  supercooled_specific_heat=c2
  supercooled_stefan_number=St
  supercooled_thermal_diff=k2

  call liquid_temperature_driver(rd,t,v2)
  
end subroutine liquid_temperature
!!$****************************************************

end module supercooled_exact_sol
