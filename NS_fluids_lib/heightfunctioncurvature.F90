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

subroutine count_sign_change(LSCOL,SIGN_COUNT)
IMPLICIT NONE

REAL_T, intent(in) :: LSCOL(-3:3)
INTEGER_T, intent(out) :: SIGN_COUNT
INTEGER_T :: i
REAL_T :: LS1,LS2

SIGN_COUNT=0
do i=-3,2
 LS1=LSCOL(i)  ! F-1/2
 LS2=LSCOL(i+1)  ! F-1/2
 if ((LS1.ge.0.0d0).and.(LS2.ge.0.0d0)) then
  ! do nothing
 else if ((LS1.lt.0.0d0).and.(LS2.lt.0.0d0)) then
  ! do nothing
 else if ((LS1.ge.0.0d0).and.(LS2.lt.0.0d0)) then
  SIGN_COUNT=SIGN_COUNT+1
 else if ((LS1.lt.0.0d0).and.(LS2.ge.0.0d0)) then
  SIGN_COUNT=SIGN_COUNT+1
 else
  print *,"LS1 or LS2 NaN"
  stop
 endif
enddo

return
end subroutine count_sign_change

subroutine get_curvature_heightf(nmat,ls_in,vf_in,delta_x,kappa,curv_valid)
implicit none

INTEGER_T, intent(in)  :: nmat
INTEGER_T, intent(inout)  :: curv_valid ! 1=success   0=failure
REAL_T, dimension(nmat,-3:3,-3:3), intent(in) :: vf_in
REAL_T, dimension(nmat,-3:3,-3:3), intent(in) :: ls_in
REAL_T,intent(in)  :: delta_x
!REAL_T,intent(in)  :: phi_center  ! signed distant function in the center cell
INTEGER_T :: i,j,imat
REAL_T :: FVERTICAL(-1:1)
INTEGER_T :: SIGN_VERTICAL(-1:1)
REAL_T :: FHORIZONTAL(-1:1)
INTEGER_T :: SIGN_HORIZONTAL(-1:1)
REAL_T :: LSCOL(-3:3)
REAL_T       :: diff1,diff2
REAL_T       :: hprime,hdoubleprime
REAL_T, intent(out) :: kappa(nmat+1)     ! curvature output
REAL_T :: kappa_scale

if(nmat.ne.2)then
 print *,"nmat has to be 2"
 stop
endif

if (curv_valid.eq.1) then
 ! do nothing
else
 print *,"curv_valid invalid"
 stop
endif

do imat=1,2

  do i=-1,1
   FVERTICAL(i)=0.0d0  ! vertical orientation
   FHORIZONTAL(i)=0.0d0  ! horizontal orientation
  enddo

!  horizontal
! sum up the volume fraction of the chosen material in each column
  do j=-1,1
  do i=-3,3
   FHORIZONTAL(j)=FHORIZONTAL(j)+vf_in(imat,j,i)  ! column sum
   LSCOL(i)=ls_in(imat,j,i)
  enddo
  call count_sign_change(LSCOL,SIGN_HORIZONTAL(j))
  enddo
  diff1=abs(FHORIZONTAL(1)-FHORIZONTAL(-1))
! vertical 
! sum up the volume fraction of the chosen material in each row
  do j=-1,1
  do i=-3,3
   FVERTICAL(j)=FVERTICAL(j)+vf_in(imat,i,j)
   LSCOL(i)=ls_in(imat,i,j)
  enddo
  call count_sign_change(LSCOL,SIGN_VERTICAL(j))
  enddo
  diff2=abs(FVERTICAL(1)-FVERTICAL(-1))

! determine the target cell interface orientation from center cell interface normal 
  if(diff1.lt.diff2)then
   ! more horizontal
   hprime=(FHORIZONTAL(1)-FHORIZONTAL(-1))/(2.0d0*delta_x)
   hdoubleprime=(FHORIZONTAL(1)+FHORIZONTAL(-1)-2.0d0*FHORIZONTAL(0))/ &
    (delta_x*delta_x)

   if ((SIGN_HORIZONTAL(-1).eq.1).and. &
       (SIGN_HORIZONTAL(0).eq.1).and. &
       (SIGN_HORIZONTAL(1).eq.1)) then
    ! do nothing
   else
    curv_valid=0
   endif

  else if (diff1.ge.diff2) then
   ! more vertical 
   hprime=(FVERTICAL(1)-FVERTICAL(-1))/(2.0d0*delta_x)
   hdoubleprime=(FVERTICAL(1)+FVERTICAL(-1)-2.0d0*FVERTICAL(0))/ &
    (delta_x*delta_x)

   if ((SIGN_VERTICAL(-1).eq.1).and. &
       (SIGN_VERTICAL(0).eq.1).and. &
       (SIGN_VERTICAL(1).eq.1)) then
    ! do nothing
   else
    curv_valid=0
   endif
  else
   print *,"diff1 or diff2 NaN"
   stop
  endif
   ! assume delta_x=delta_y
  hprime=hprime*delta_x
  hdoubleprime=hdoubleprime*delta_x

  if (curv_valid.eq.1) then
   kappa(imat)=hdoubleprime/((1.0d0+hprime**2.0d0)**1.5d0)
  else if (curv_valid.eq.0) then
   ! do nothing
  else
   print *,"curv_valid invalid"
   stop
  endif

enddo ! imat=1,2

if (curv_valid.eq.0) then
 kappa(1)=0.0d0
 kappa(2)=0.0d0
else if (curv_valid.eq.1) then
 kappa_scale=max(abs(kappa(1)),abs(kappa(2)))
 kappa_scale=max(kappa_scale,1.0d0)
 if (abs(kappa(1)+kappa(2)).ge.1.0D-12*kappa_scale) then
  print *,"kappa(1) and kappa(2) failed sanity check"
  print *,"kappa(1)=",kappa(1)
  print *,"kappa(2)=",kappa(2)
  print *,"kappa_scale=",kappa_scale
  stop
 endif
endif

! approximate the curvature of at the interface ( which has the same x or y coordinate with cell center) in the center cell 



end subroutine get_curvature_heightf

end module height_method_module
