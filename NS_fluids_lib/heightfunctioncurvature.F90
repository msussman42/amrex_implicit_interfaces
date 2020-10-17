#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

module height_method_module

contains

subroutine get_curvature_heightf(nmat,vf_in,delta_x,kappa)
implicit none

INTEGER_T, intent(in)  :: nmat
REAL_T, dimension(nmat,-3:3,-3:3), intent(in) :: vf_in
REAL_T,intent(in)  :: delta_x
!REAL_T,intent(in)  :: phi_center  ! signed distant function in the center cell
INTEGER_T :: i,imat
REAL_T       :: FL1,F1,FR1,FL2,F2,FR2
REAL_T       :: diff1,diff2
REAL_T       :: hprime,hdoubleprime
REAL_T       :: kappa(nmat)     ! curvature output

FL1=0.0d0
F1=0.0d0
FR1=0.0d0
FL2=0.0d0
F2=0.0d0
FR2=0.0d0

if(nmat.ne.2)then
 print *,"nmat has to be 2"
 stop
endif

do imat=1,2
!  horizontal
! sum up the volume fraction of the chosen material in each column
  do i=-3,3
   FL1=FL1+vf_in(imat,-1,i)
   F1=F1+vf_in(imat,0,i)
   FR1=FR1+vf_in(imat,1,i)
   diff1=abs(FL1-FR1)
  enddo
! vertical 
! sum up the volume fraction of the chosen material in each row
  do i=-3,3
   FL2=FL2+vf_in(imat,i,-1)
   F2=F2+vf_in(imat,i,0)
   FR2=FR2+vf_in(imat,i,1) 
   diff2=abs(FL2-FR2)
  enddo

! determine the target cell interface orientation from center cell interface normal 
  if(diff1.lt.diff2)then
   ! more horizontal
   hprime=(FR1-FL1)/(2.0d0*delta_x)
   hdoubleprime=(FR1+FL1-2.0d0*F1)/(delta_x*delta_x)
  else   
   ! more vertical 
   hprime=(FR2-FL2)/(2.0d0*delta_x)
   hdoubleprime=(FR2+FL2-2.0d0*F2)/(delta_x*delta_x)
  endif

  kappa(imat)=hdoubleprime/((1.0d0+hprime**2.0d0))**1.5d0

enddo

! approximate the curvature of at the interface ( which has the same x or y coordinate with cell center) in the center cell 



end subroutine get_curvature_heightf

end module height_method_module
