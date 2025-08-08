subroutine mod_bessel_first_kind(x,n,I_n)
IMPLICIT NONE

real*8, intent(in) :: x
real*8, intent(out) :: I_n
integer, intent(in) :: n
real*8 m_fact,mpn_fact,xpower
integer i

m_fact=1.0d0
mpn_fact=1.0d0
do i=1,n
 mpn_fact=mpn_fact*i
enddo
xpower=(0.5d0*x)**(n)
I_n=0.0d0
do i=0,100
 if (i.gt.0) then
  m_fact=m_fact*i
  mpn_fact=mpn_fact*(i+n)
  xpower=xpower*x*x/4.0d0
 endif
 I_n=I_n+xpower/(m_fact*mpn_fact)
enddo

return
end subroutine mod_bessel_first_kind

! c=f(k)  k=wave number  0< k r_0 < 1
! page 5860 section 6.5 Popinet 2009  
! u_t = N(u) + F   (1)
! N(u) is a nonlinear differential operator.  After discretization,
! N(u) is some complication nonlinear function of u where
! "u" is a vector in R^{d} and N is a vector function also in R^{d}.
! given a base state in which N(u_base)+F=0, we perturb the base state
! u_perturb=u_base + du e^{beta t} and plug back into (1).
! beta du e^{beta t} = grad N(u_base) du e^{beta t} + O(du^2)
! define A=grad N(u_base)  A is a "Jacobian matrix" with dimensions
! d \Times d.  A_{ij}=\partial N_{i}/\partial u_{j}.
! beta du = A du + O(du^2)
! Objective: find the eigenvector "du" with the largest eigenvalue |Re beta|
! advantage of the power method: no need to explicitly find "A"
! Disadvantages of the "non intrusive Linear Stability Analysis" method:
! 1. N(u) might not have a continuous derivative or second derivative.
! 2. max |\lambda(A)| might have a real part that is less than or equal to
! zero.
! growth_rate^2 = I_1(k*r0)*k*r0*(1-k*r0)/(I_0(k*r0))
! anecdotedly, k r0 = .7 is the critical value.
! wave length = 2 pi/(.7/r0)=2*pi*r0/.7
program rayleigh_capillary_growth
real*8 density,r0,sigma,k,x,I_0,I_1,c,c2
integer i

density=1.0d0
r0=1.0d0
sigma=1.0d0

k=0.7
x=k*r0
call mod_bessel_first_kind(x,0,I_0)
call mod_bessel_first_kind(x,1,I_1)
c2=sigma*I_1*x*(1.0d0-x*x)/(density*(r0**3)*I_0)
c=sqrt(c2)
print *,"specific data at a single point: k,c ",k,c

do i=0,100
 k=(i/100.0d0)/r0
 x=k*r0
 call mod_bessel_first_kind(x,0,I_0)
 call mod_bessel_first_kind(x,1,I_1)
 c2=sigma*I_1*x*(1.0d0-x*x)/(density*(r0**3)*I_0)
 c=sqrt(c2)
 print *,"k,c ",k,c
enddo

end program
