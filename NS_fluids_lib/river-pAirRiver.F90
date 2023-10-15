#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

MODULE River
use amrex_fort_module, only : amrex_real

!Parameter functions for two phase incompressible
implicit none
CONTAINS

! option=0,1,2, ...
subroutine RiverVelocity(x,y,z,vel,option,probloz,probhiz)
        real(amrex_real), INTENT(in) :: x,y,z,probloz,probhiz
        integer, INTENT(in) :: option
        real(amrex_real),INTENT(out) :: vel(SDIM)
        real(amrex_real) h,theta,velfreestream
        integer dir

 call RiverHeight(x,y,h,option)

 if (SDIM.eq.2) then
  if (abs(z-y).gt.1.0E-8) then
   print *,"expecting z=y"
   stop
  endif
 endif

 velfreestream=1.0
 do dir=1,SDIM
  vel(dir)=0.0
 enddo

 if (z.le.h) then
  ! do nothing, velocity initially 0 in the liquid.
 else if (z.ge.probhiz) then
  vel(1)=velfreestream
 else 
  theta=velfreestream*(z-h)/(probhiz-h)
  if ((theta.lt.0.0).or.(theta.gt.1.0)) then
   print *,"theta out of range"
   stop
  endif
  vel(1)=theta 
 endif 

end subroutine RiverVelocity

! option is axis_dir in the inputs file.
! option=0,1,2,...
SUBROUTINE RiverHeight(x,y,h,option)
        integer,INTENT(in) :: option
        real(amrex_real),INTENT(in) :: x,y
        real(amrex_real),INTENT(out) :: h
        real(amrex_real) :: a

        a = 0.05

        if (option == 0) then 
         h = 0.0     !(constant height interface)
        else if (option == 1) then 
         if (SDIM.eq.3) then
                h = a*sin(x + y) !sinsoidic height interface with amplitude a
         else if (SDIM.eq.2) then
                h = a*sin(x) !sinsoidic height interface with amplitude a
         else
          print *,"dimension bust"
          stop
         endif
        end if 

END SUBROUTINE RiverHeight


! option=0,1,2, ...
subroutine RiverPressure(x,y,z,t,p,denair,denwater,option)
use probcommon_module
use global_utility_module

real(amrex_real), INTENT(in)  :: x,y,z,t,denair,denwater
integer, INTENT(in)  :: option
real(amrex_real), INTENT(out) :: p
real(amrex_real) h
integer :: gravity_dir

 call fort_derive_gravity_dir(gravity_vector,gravity_dir)

 if (SDIM.eq.2) then
  if (abs(z-y).gt.1.0E-8) then
   print *,"expecting z=y"
   stop
  endif
 endif
 call RiverHeight(x,y,h,option)
 if (z.ge.h) then
  p=denair*abs(gravity_vector(gravity_dir))*(h-z)
  if (x.lt.0.0) then
   p=p+abs(gravity_vector(gravity_dir))
  else if (x.gt.4.0) then
   ! do nothing
  else
   p=p+(abs(gravity_vector(gravity_dir))*((4.0-x)/4.0))
  endif
 else
  p=denwater*abs(gravity_vector(gravity_dir))*(h-z)
  if (x.lt.0.0) then
   p=p+abs(gravity_vector(gravity_dir))
  else if (x.gt.4.0) then
   ! do nothing
  else
   p=p+(abs(gravity_vector(gravity_dir))*((4.0-x)/4.0))
  endif
 endif
   
end subroutine RiverPressure


END MODULE River
