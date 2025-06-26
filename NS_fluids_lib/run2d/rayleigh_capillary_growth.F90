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

 
program rayleigh_capillary_growth
real*8 density,r0,sigma,k,x,I_0,I_1,c,c2
integer i

density=1.0d0
r0=1.0d0
sigma=1.0d0
do i=0,100
 k=i/100.0d0
 x=k*r0
 call mod_bessel_first_kind(x,0,I_0)
 call mod_bessel_first_kind(x,1,I_1)
 c2=sigma*I_1*x*(1.0d0-x*x)/(density*(r0**3)*I_0)
 c=sqrt(c2)
 print *,"k,c ",k,c
enddo

end program
