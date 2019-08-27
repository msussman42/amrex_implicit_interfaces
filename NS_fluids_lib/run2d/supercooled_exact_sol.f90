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
    REAL(kind=dp), parameter::em_const=0.5772156649015329_dp
    REAL(kind=dp), parameter::tol=1E-15_dp
    INTEGER:: k


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
    
    f_lambda = lm**2 * exp(-lm**2) * Ei(-lm**2) + St;

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

    REAL(kind=dp), parameter:: dx = 1E-2_dp     !! Marching step
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

  R = 2 * lm * sqrt(k2 * t)
end subroutine solidification_front_radius
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


  v2 = V + (L * St / c2)* Ei(-rd*rd/(4*k2*t)) / Ei(-lm*lm)
end subroutine liquid_temperature
!!$****************************************************

end module supercooled_exact_sol
