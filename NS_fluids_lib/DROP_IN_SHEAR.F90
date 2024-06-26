#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

#include "AMReX_ArrayLim.H"
#include "EXTRAP_COMP.H"


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

! probtype==424
module DROP_IN_SHEAR_module
use amrex_fort_module, only : amrex_real

implicit none                   

real(amrex_real), PARAMETER :: Y_inf_default = 7.1D-3

real(amrex_real) :: DEF_VAPOR_GAMMA
real(amrex_real) :: den_G,C_pG,k_G,lambda,T_inf,T_sat,L_V,D_G,Y_inf
real(amrex_real) :: den_L
real(amrex_real) :: WV,WA,Le
real(amrex_real) :: D_not,B_M,Sh
real(amrex_real) :: T_gamma,X_gamma,Y_gamma
real(amrex_real) :: T_gamma_liquid
real(amrex_real) :: dr_max_level

contains

subroutine f_mdot(T_gamma_parm,f_out)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: T_gamma_parm
real(amrex_real), INTENT(out) :: f_out
real(amrex_real) :: X_gamma_loc,Y_gamma_loc

call X_from_Tgamma(X_gamma_loc,T_gamma_parm,T_sat,L_V, &
 fort_R_Palmore_Desjardins,WV)
call massfrac_from_volfrac(X_gamma_loc,Y_gamma_loc,WA,WV)

if ((Y_gamma_loc.ge.zero).and.(Y_gamma_loc.lt.one)) then
 if (L_V.gt.zero) then
  if (C_pG.gt.zero) then
   if ((Y_inf.ge.zero).and.(Y_inf.lt.one)) then
    if ((T_inf.gt.zero).and.(T_inf.lt.10.0d0*T_sat)) then
     f_out=((Y_inf-one)/(Y_gamma_loc-one))**Le - one + &
      (T_gamma_parm-T_inf)*C_pG/L_V
    else
     print *,"T_inf invalid"
     stop
    endif
   else
    print *,"Y_inf invalid"
    stop
   endif
  else
   print *,"C_pG invalid"
   stop
  endif
 else
  print *,"L_V invalid"
  stop
 endif
else
 print *,"Y_gamma invalid"
 stop
endif

return
end subroutine f_mdot

  ! do any initial preparation needed
subroutine INIT_DROP_IN_SHEAR_MODULE()
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real) :: a,b,c,f_a,f_b,f_c
real(amrex_real) :: a_min,f_a_min
integer :: iter
integer :: ilev

if (fort_n_cell(1).gt.0) then
 dr_max_level=(probhix-problox)/fort_n_cell(1)
 do ilev=1,fort_max_level
  dr_max_level=0.5d0*dr_max_level
 enddo
else
 print *,"fort_n_cell(1) invalid"
 stop
endif
 
DEF_VAPOR_GAMMA =  1.666666667D0

! ergs/(mol kelvin) is the default for NavierStokes::R_Palmore_Desjardins

den_L = fort_denconst(1)
den_G = fort_denconst(2)
C_pG = fort_stiffCP(2)
k_G = fort_heatviscconst(2)
lambda=k_G/(den_G*C_pG)
T_gamma_liquid=fort_tempconst(1)
T_inf = fort_tempconst(2)
T_sat = fort_saturation_temp(1)
L_V = get_user_latent_heat(1,room_temperature,1)
D_G = fort_speciesviscconst(2)
Y_inf=fort_speciesconst(2)
if (Y_inf.eq.1.0d0) then
 Y_inf=radblob2
 if (Y_inf.eq.0.0d0) then
  Y_inf=Y_inf_default
 else if ((Y_inf.gt.0.0d0).and.(Y_inf.lt.1.0d0)) then
  ! do nothing
 else
  print *,"radblob2 invalid"
  stop
 endif
else if ((Y_inf.ge.0.0d0).and.(Y_inf.lt.1.0d0)) then
 ! do nothing
else
 print *,"Y_inf invalid"
 stop
endif
WV=fort_species_molar_mass(1)  !num_species components
WA=fort_molar_mass(2)
! T_inf C_pG / L = T_inf (K) * C_pG (erg/(g K)) / (L (erg/g))
!   = dimensionless
! e.g. Le=D_G cm^2/s * den_G g/cm^3 * C_pG erg/(g K) / 
!         (k_G g cm/(s^3 K)) = dimensionless
! erg=g cm^2/s^2
! cm^2/s  * g/cm^3  * g cm^2 / s^2 * (1/(g K)) * (s^3 K)/(g cm) = dimensionless
!
Le=D_G*den_G*C_pG/k_G

print *,"Parameters for analytical solution for evaporating droplet: "
print *,"den_L,den_G,C_pG,k_G,lambda ",den_L,den_G,C_pG,k_G,lambda
print *,"T_inf,T_sat,L_V,D_G,Y_inf,WV,WA,Le ",T_inf,T_sat,L_V,D_G,Y_inf,WV,WA,Le
a=EVAP_BISECTION_TOL
b=T_sat*(one-EVAP_BISECTION_TOL)
call f_mdot(a,f_a)
call f_mdot(b,f_b)
if (f_b.gt.zero) then
 ! do nothing
else
 print *,"f_b invalid"
 stop
endif
a_min=a
f_a_min=f_a
do while ((f_a.gt.zero).and.(a.lt.b))
 print *,"searching for bracketing interval,a,b,f_a,f_b ", &
   a,b,f_a,f_b
 a=a+0.001d0*T_sat
 if (a.lt.b) then
  call f_mdot(a,f_a)
  if (abs(f_a).lt.abs(f_a_min)) then
   a_min=a
   f_a_min=f_a
  endif
 else
  print *,"a cannot exceed b: a,f_a,b,f_b ",a,f_a,b,f_b
  print *,"a_min,f_a_min ",a_min,f_a_min
  stop
 endif
enddo
if (f_a*f_b.le.zero) then
 do iter=1,100
  c=0.5d0*(a+b)
  call f_mdot(c,f_c)
  if (f_a*f_c.le.zero) then
   b=c
   f_b=f_c
  else if (f_b*f_c.le.zero) then
   a=c
   f_a=f_c
  else
   print *,"f_c became corrupt"
   stop
  endif
 enddo ! iter=1..100
else
 print *,"f_a and f_b must have different signs"
 stop
endif

T_gamma=c  
call X_from_Tgamma(X_gamma,T_gamma,T_sat,L_V, &
 fort_R_Palmore_Desjardins,WV)
call massfrac_from_volfrac(X_gamma,Y_gamma,WA,WV)

B_M=(Y_gamma-Y_inf)/(one-Y_gamma)
D_not=two*radblob
if (B_M.gt.zero) then
 Sh=two*log(one+B_M)/B_M
else
 print *,"B_M must be positive"
 stop
endif

print *,"INIT_DROP_IN_SHEAR_MODULE T_gamma,T_gamma_liquid,Y_gamma ", &
        T_gamma,T_gamma_liquid,Y_gamma
print *,"dr_max_level=",dr_max_level

return
end subroutine INIT_DROP_IN_SHEAR_MODULE

! dir=velocity component
subroutine DROP_IN_SHEAR_CFL_HELPER(time,dir,uu,dx)
use probcommon_module
implicit none
integer, INTENT(in) :: dir
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) :: uu
real(amrex_real), INTENT(in) :: dx(SDIM)

real(amrex_real) utest


if ((dir.lt.0).or.(dir.ge.SDIM)) then
 print *,"dir invalid"
 stop
endif

if (dir.eq.0) then
 ! do nothing
else if (dir.eq.1) then
 ! do nothing
else if ((dir.eq.2).and.(SDIM.eq.3)) then
 ! do nothing
else
 print *,"dir invalid DROP_IN_SHEAR_CFL_HELPER"
 stop
endif

! shear flow is in the x direction
if (probtype.eq.424) then
 if (dir.eq.1) then
  utest=abs(vinletgas)
  uu=max(abs(uu),abs(utest))
 endif
else
 print *,"unexpected probtype"
 stop
endif

return
end subroutine DROP_IN_SHEAR_CFL_HELPER

! dist>0 in the fluid
subroutine DROP_IN_SHEAR_soliddist(x,dist,im) 
use probcommon_module
use global_utility_module
implicit none
real(amrex_real), INTENT(in), dimension(SDIM) :: x !spatial coordinates
integer, INTENT(in) :: im
real(amrex_real), INTENT(out) :: dist

if (num_materials.lt.1) then
 print *,"num_materials invalid in soliddist"
 stop
endif

if ((im.lt.1).or.(im.gt.num_materials)) then
 print *,"im invalid11"
 stop
endif

if (probtype.eq.424) then
 print *,"no embedded solids"
 stop
else
 print *,"expecting probtype.eq.424"
 stop
endif

end subroutine DROP_IN_SHEAR_soliddist


 ! fluids tessellate the domain, solids are immersed. 
subroutine DROP_IN_SHEAR_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(out) :: LS(nmat)
integer :: im
integer :: im_solid_materialdist
real(amrex_real) :: initial_time

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  im_solid_materialdist=im_solid_primary()
  initial_time=zero

  if (probtype.eq.424) then

   do im=1,nmat
     ! prescribed rigid solid (PROB.F90)
    if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then
     call DROP_IN_SHEAR_soliddist(x,LS(im),nmat)  ! returns LS<0 in solid
     LS(im)=-LS(im)   ! now LS>0 in solid
    endif
   enddo

   if (axis_dir.eq.0) then  ! drop in shear flow
    if (SDIM.eq.2) then
     LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)   ! liquid
    else if (SDIM.eq.3) then
     LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+ &
             (x(SDIM)-zblob)**2)   ! liquid
    else
     print *,"dimension bust"
     stop
    endif
    LS(2)=-LS(1) ! ambient gas

    ! above: drop in shear flow
    ! below: gas flow over a flat liquid interface
   else if (axis_dir.eq.1) then 
    LS(1)=radblob-x(SDIM)
    LS(2)=-LS(1) ! ambient gas
   else
    print *,"axis_dir invalid"
    stop
   endif

  else
   print *,"expecting probtype.eq.424"
   stop
  endif

return
end subroutine DROP_IN_SHEAR_LS

! initial velocity is some kind of shear flow
subroutine DROP_IN_SHEAR_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: VEL(SDIM)
real(amrex_real) :: pres_analytical
integer dir
integer, INTENT(in) :: velsolid_flag
real(amrex_real) :: vert_lo,vert_hi
real(amrex_real) :: D_gamma,T_analytical,Y_analytical,LS_analytical

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

if ((velsolid_flag.eq.0).or. &
    (velsolid_flag.eq.1)) then
 ! do nothing
else 
 print *,"velsolid_flag invalid"
 stop
endif

do dir=1,SDIM
 if (dx(dir).gt.zero) then
  ! do nothing
 else
  print *,"dx invalid"
  stop
 endif
enddo

do dir=1,SDIM
 VEL(dir)=zero
enddo

if (probtype.eq.424) then

 do dir=1,SDIM
  VEL(dir)=zero
 enddo

 if (SDIM.eq.2) then
  vert_lo=probloy
  vert_hi=probhiy
 else if (SDIM.eq.3) then
  vert_lo=probloz
  vert_hi=probhiz
 else
  print *,"dimension bust"
  stop
 endif 
 if (vert_hi-vert_lo.gt.zero) then
  ! do nothing
 else
  print *,"vert_hi-vert_lo invalid"
  stop
 endif

 if (axis_dir.eq.0) then ! drop in shear
  VEL(1)=-vinletgas+two*vinletgas*(x(SDIM)-vert_lo)/(vert_hi-vert_lo)
  if (vinletgas.eq.zero) then
   if (1.eq.1) then
    ! do nothing
   else if (1.eq.0) then
    call drop_analytical_solution(t,x,D_gamma,T_analytical, &
      Y_analytical,VEL,LS_analytical,pres_analytical)
   endif
  else if (vinletgas.ne.zero) then
   ! do nothing
  else
   print *,"vinletgas invalid"
   stop
  endif
 else if (axis_dir.eq.1) then ! liquid layer in shear
  if (x(SDIM).le.radblob) then
   VEL(1)=zero ! liquid (periodic boundary conditions xlo and xhi)
  else if (x(SDIM).ge.radblob) then
   VEL(1)=vinletgas*(x(SDIM)-radblob)/(vert_hi-radblob) ! gas
  else
   print *,"x(SDIM) invalid"
   stop
  endif
 else
  print *,"axis_dir invalid"
  stop
 endif
else
 print *,"expecting probtype==424"
 stop
endif

return 
end subroutine DROP_IN_SHEAR_VEL


! this routine used as a default when
! pressure boundary conditions are prescribed.
! For the case when only top wall is 
! "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine DROP_IN_SHEAR_PRES(x,t,LS,PRES,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: PRES
real(amrex_real) :: D_gamma,T_analytical,Y_analytical,LS_analytical
real(amrex_real) :: pres_analytical
real(amrex_real) :: VEL(SDIM)

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if (probtype.eq.424) then

 PRES=zero

 if (axis_dir.eq.0) then
  if (vinletgas.eq.zero) then
   call drop_analytical_solution(t,x,D_gamma,T_analytical, &
      Y_analytical,VEL,LS_analytical,pres_analytical)
   PRES=pres_analytical
  endif
 else if (axis_dir.eq.1) then
  ! do nothing
 else
  print *,"axis_dir invalid"
  stop
 endif

else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif

return 
end subroutine DROP_IN_SHEAR_PRES



subroutine DROP_IN_SHEAR_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
integer, INTENT(in) :: nmat
integer, INTENT(in) :: nstate_mat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: STATE(nmat*nstate_mat)
integer im,ibase,n
real(amrex_real) :: D_gamma,T_analytical,Y_analytical,LS_analytical
real(amrex_real) :: pres_analytical
real(amrex_real) :: VEL(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if (nstate_mat.eq.num_state_material) then
 ! do nothing
else
 print *,"nstate_mat invalid"
 stop
endif

if (probtype.eq.424) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im) 
  if (t.eq.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im) 
  else if (t.gt.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif

   ! initial species in inputs?
  do n=1,num_species_var
   STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo

 enddo ! im=1..num_materials

 if (axis_dir.eq.0) then
   !stagnant
  if (vinletgas.eq.zero) then
   call drop_analytical_solution(t,x,D_gamma,T_analytical, &
      Y_analytical,VEL,LS_analytical,pres_analytical)
   do im=1,num_materials
    ibase=(im-1)*num_state_material
    STATE(ibase+ENUM_TEMPERATUREVAR+1)=T_analytical
    STATE(ibase+ENUM_SPECIESVAR+1)=Y_analytical
   enddo
  endif
 else if (axis_dir.eq.1) then
  ! do nothing
 else
  print *,"axis_dir invalid"
  stop
 endif
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine DROP_IN_SHEAR_STATE

 ! dir=1..sdim  side=1..2
subroutine DROP_IN_SHEAR_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(inout) :: LS(nmat)
real(amrex_real), INTENT(in) :: LS_in(nmat)
integer, INTENT(in) :: dir,side
real(amrex_real), INTENT(in) ::  dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call DROP_IN_SHEAR_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine DROP_IN_SHEAR_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine DROP_IN_SHEAR_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(inout) :: VEL
real(amrex_real), INTENT(in) :: VEL_in
integer, INTENT(in) :: veldir,dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) local_VEL(SDIM)
integer velsolid_flag

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if (probtype.eq.424) then
 ! do nothing
else
 print *,"expecting probtype==424"
 stop
endif
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call DROP_IN_SHEAR_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine DROP_IN_SHEAR_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine DROP_IN_SHEAR_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(inout) :: PRES
real(amrex_real), INTENT(in) :: PRES_in
integer, INTENT(in) :: dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) :: rr

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

 call DROP_IN_SHEAR_PRES(xghost,t,LS,PRES,nmat)

 rr=xghost(1)**2+xghost(2)**2
 if (SDIM.eq.3) then
  rr=rr+xghost(SDIM)**2
 endif
 rr=sqrt(rr)
 if (1.eq.0) then
  print *,"x,y,r,t,pres ",xghost(1),xghost(2),rr,t,PRES
 endif

return
end subroutine DROP_IN_SHEAR_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine DROP_IN_SHEAR_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real) :: local_STATE(nmat*num_state_material)
real(amrex_real), INTENT(inout) :: STATE
real(amrex_real), INTENT(inout) :: STATE_merge
real(amrex_real), INTENT(in) :: STATE_in
integer, INTENT(in) :: dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)
integer, INTENT(in) :: istate,im
integer ibase,im_crit
integer local_bcflag

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
local_bcflag=1

if ((istate.ge.1).and. &
    (istate.le.num_state_material).and. &
    (im.ge.1).and. &
    (im.le.num_materials)) then
 call DROP_IN_SHEAR_STATE(xghost,t,LS,local_STATE,local_bcflag, &
         nmat,num_state_material)
 ibase=(im-1)*num_state_material
 STATE=local_STATE(ibase+istate)
 call get_primary_material(LS,im_crit)
 ibase=(im_crit-1)*num_state_material
 STATE_merge=local_STATE(ibase+istate)

 if (probtype.eq.424) then
  ! do nothing
 else
  print *,"expecting probtype == 424"
  stop
 endif
else
 print *,"istate invalid"
 stop
endif

return
end subroutine DROP_IN_SHEAR_STATE_BC

subroutine DROP_IN_SHEAR_HEATSOURCE(im,VFRAC,time,x, &
     xsten,nhalf,temp, &
     heat_source,den,CV,dt,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
integer, INTENT(in) :: im
real(amrex_real), INTENT(in) :: VFRAC(nmat)
real(amrex_real), INTENT(in) :: time
integer, INTENT(in) :: nhalf
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
real(amrex_real), INTENT(in) :: temp(nmat)
real(amrex_real), INTENT(in) :: den(nmat)
real(amrex_real), INTENT(in) :: CV(nmat)
real(amrex_real), INTENT(in) :: dt
real(amrex_real), INTENT(out) :: heat_source

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
 
return
end subroutine DROP_IN_SHEAR_HEATSOURCE


subroutine DROP_IN_SHEAR_EB_heat_source(time,dt,xsten,nhalf, &
      heat_flux,heat_dir,heat_side)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nhalf
real(amrex_real), dimension(-nhalf:nhalf,SDIM), INTENT(in) :: xsten
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(in) :: dt
real(amrex_real), INTENT(out) :: heat_flux
integer, INTENT(out) :: heat_dir
integer, INTENT(out) :: heat_side

if (time.lt.zero) then
 print *,"time invalid"
 stop
endif
if (dt.le.zero) then
 print *,"dt invalid"
 stop
endif

if (probtype.eq.424) then

 heat_flux=zero
 heat_dir=0
 heat_side=0

else
 print *,"expecting probtype==424"
 stop
endif

return
end subroutine DROP_IN_SHEAR_EB_heat_source

  ! only called at faces with an adjoining solid cell and
  ! an adjoining fluid cell.
subroutine DROP_IN_SHEAR_microcell_heat_coeff(heatcoeff,dx,veldir)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: dx(SDIM)
integer, INTENT(in) :: veldir
real(amrex_real), INTENT(inout) :: heatcoeff

if (probtype.eq.424) then

 ! do nothing

else
 print *,"expecting probtype==424"
 stop
endif

return
end subroutine DROP_IN_SHEAR_microcell_heat_coeff

subroutine DROP_IN_SHEAR_velfreestream(problen,local_buffer)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(inout) :: local_buffer(2*SDIM)
real(amrex_real), INTENT(in)    :: problen(SDIM)

if (probtype.eq.424) then
 ! do nothing
else
 print *,"expecting probtype==424"
 stop
endif

return
end subroutine DROP_IN_SHEAR_velfreestream


subroutine DROP_IN_SHEAR_nucleation(nucleate_in,xsten,nhalf,make_seed)
use probcommon_module_types
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: nhalf
real(amrex_real), dimension(-nhalf:nhalf,SDIM), INTENT(in) :: xsten
integer, INTENT(inout) :: make_seed
type(nucleation_parm_type_input), INTENT(in) :: nucleate_in
real(amrex_real) :: LL

LL=nucleate_in%LL

if (probtype.eq.424) then
 ! do nothing
else
 print *,"expecting probtype==424"
 stop
endif

return
end subroutine DROP_IN_SHEAR_nucleation

subroutine drop_analytical_solution(time,x,D_gamma,T,Y,VEL,LS_VAP,PRES)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(out) :: D_gamma
real(amrex_real), INTENT(out) :: T
real(amrex_real), INTENT(out) :: Y
real(amrex_real), INTENT(out) :: LS_VAP
real(amrex_real), INTENT(out) :: PRES
real(amrex_real), INTENT(out) :: VEL(SDIM)
real(amrex_real) :: rr,mdot,vel_r
real(amrex_real) :: VELCOEFF

if (SDIM.eq.2) then
 rr=sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)
else if (SDIM.eq.3) then
 rr=sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)
else
 print *,"dimension bust"
 stop
endif

if ((den_G.lt.fort_denconst(1)).and. &
    (den_G.gt.zero)) then
 ! do nothing
else
 print *,"den_G (fort_denconst(2)) invalid"
 stop
endif
if (num_species_var.eq.1) then
 ! do nothing
else
 print *,"num_species_var invalid"
 stop
endif

if (B_M.gt.zero) then
 ! do nothing
else
 print *,"expecting B_M>0"
 stop
endif

D_Gamma=D_not**2-eight*den_G*D_G*log(one+B_M)*time/den_L
if ((D_Gamma.le.D_not**2).and.(D_Gamma.ge.zero)) then
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
mdot=Pi*D_gamma*den_G*D_G*Sh*B_M
if (mdot.gt.zero) then
 ! do nothing
else
 print *,"expecting mdot>0"
 stop
endif

LS_VAP=rr-half*D_gamma

if (LS_VAP.le.zero) then
 VEL(1)=zero
 VEL(2)=zero
 VEL(SDIM)=zero
 Y=Y_gamma
 T=T_gamma_liquid
 PRES=zero
else if (LS_VAP.ge.zero) then
 VELCOEFF = mdot/(4.0d0*Pi*den_G)
 vel_r = VELCOEFF/(rr*rr)
  ! for pressure:
  ! div( u^2 ) = -p_r/den_G
  ! u=M/r^2  u^2=M^2/r^4  div u^2=(1/r^2) (M^2/r^2)_r=
  ! -2/r^5  *  M^2
  ! p=-int(-2/r^5  *  M^2)*den_G=int(2/r^5  *  M^2)*den_G= 
  ! den_G(C-(1/2) M^2/r^4)
 PRES=half*den_G*(VELCOEFF**2)*(one/(radblob**4)-one/(rr**4))

 VEL(1)=vel_r*(x(1)-xblob)/rr
 VEL(2)=vel_r*(x(2)-yblob)/rr
 if (SDIM.eq.3) then
  VEL(SDIM)=vel_r*(x(SDIM)-zblob)/rr
 endif
 Y=one+(Y_inf-one)*exp(-mdot/(four*Pi*den_G*D_G*rr))
 if ((Y.ge.zero).and.(Y.lt.one)) then
  ! do nothing
 else
  print *,"Y invalid"
  stop
 endif
 T=T_gamma-L_V/C_pG+(T_inf-T_gamma+L_V/C_pG)* &
         exp(-mdot*C_pG/(four*Pi*k_G*rr))
 if (T.gt.zero) then
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


subroutine DROP_IN_SHEAR_ASSIMILATE( &
     assimilate_in,assimilate_out, &
     i,j,k,cell_flag)
use probcommon_module
IMPLICIT NONE

type(assimilate_parm_type), INTENT(in) :: assimilate_in
type(assimilate_out_parm_type), INTENT(inout) :: assimilate_out
integer, INTENT(in) :: i,j,k,cell_flag

integer :: nstate,nstate_test
real(amrex_real) :: rr,r_exact,tcrit
real(amrex_real) :: xcrit(SDIM)
integer :: dir
integer :: im
integer :: ibase
real(amrex_real) T_exact,Y_exact,LS_VAP_exact
real(amrex_real) VEL_exact(SDIM)
real(amrex_real) PRES_exact
real(amrex_real) :: D_gamma
integer :: ok_to_overwrite_vel

ok_to_overwrite_vel=0

nstate=assimilate_in%nstate

nstate_test=STATE_NCOMP
if (nstate.eq.nstate_test) then
 ! do nothing
else
 print *,"nstate invalid"
 print *,"nstate=",nstate
 print *,"nstate_test=",nstate_test
 stop
endif

if ((num_materials.eq.2).and. &
    (num_state_material.eq.3).and. & ! density, temperature, species
    (probtype.eq.424)) then
 if (axis_dir.eq.0) then
  if (vinletgas.eq.zero) then
   do dir=1,SDIM
    xcrit(dir)=assimilate_in%xsten(0,dir)
   enddo
   if (SDIM.eq.2) then
    rr=sqrt((xcrit(1)-xblob)**2+(xcrit(2)-yblob)**2)
   else if (SDIM.eq.3) then
    rr=sqrt((xcrit(1)-xblob)**2+(xcrit(2)-yblob)**2+(xcrit(SDIM)-zblob)**2)
   else
    print *,"dimension bust"
    stop
   endif
 
   r_exact=(probhix-(xblob+radblob))/5.0d0
   tcrit=assimilate_in%cur_time  ! cur_time_slab

   if (rr.ge.probhix-xblob-r_exact) then
    call drop_analytical_solution(tcrit,xcrit,D_gamma,T_exact,Y_exact, &
     VEL_exact,LS_VAP_exact,PRES_exact)

    if (cell_flag.eq.0) then ! MAC GRID X
     if (ok_to_overwrite_vel.eq.1) then
      assimilate_out%macx(D_DECL(i,j,k))=VEL_exact(1)
     endif
    else if (cell_flag.eq.1) then ! MAC GRID Y
     if (ok_to_overwrite_vel.eq.1) then
      assimilate_out%macy(D_DECL(i,j,k))=VEL_exact(2)
     endif
    else if ((cell_flag.eq.2).and.(SDIM.eq.3)) then ! MAC GRID Z
     if (ok_to_overwrite_vel.eq.1) then
      assimilate_out%macz(D_DECL(i,j,k))=VEL_exact(SDIM)
     endif
    else if (cell_flag.eq.-1) then
     if (ok_to_overwrite_vel.eq.1) then
      do dir=1,SDIM
       assimilate_out%state(D_DECL(i,j,k),dir)=VEL_exact(dir)
      enddo
     endif
     do im=1,num_materials
      ibase=STATECOMP_STATES+(im-1)*num_state_material
      assimilate_out%state(D_DECL(i,j,k),ibase+ENUM_TEMPERATUREVAR+1)=T_exact
      assimilate_out%state(D_DECL(i,j,k),ibase+ENUM_SPECIESVAR+1)=Y_exact
     enddo
    else 
     print *,"cell_flag invalid"
     stop
    endif
   else if (rr.le.probhix-xblob-r_exact) then
    ! do nothing
   else
    print *,"rr bust"
    stop
   endif
  else if (vinletgas.ne.zero) then
   ! do nothing
  else
   print *,"vinletgas bust"
   stop
  endif
 else if (axis_dir.eq.1) then
  ! do nothing
 else
  print *,"axis_dir invalid"
  stop
 endif
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine DROP_IN_SHEAR_ASSIMILATE


! This routine is called from FORT_SUMMASS
! set ns.ncomp_sum_int_user=
subroutine DROP_IN_SHEAR_SUMINT(GRID_DATA_IN,increment_out1, &
                increment_out2,nsum1,nsum2,isweep)
use probcommon_module_types
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nsum1,nsum2,isweep
type(user_defined_sum_int_type), INTENT(in) :: GRID_DATA_IN
real(amrex_real), INTENT(inout) :: increment_out1(nsum1)
real(amrex_real), INTENT(inout) :: increment_out2(nsum2)
integer :: i,j,k,dir
real(amrex_real) :: xlocal(SDIM)
real(amrex_real) :: D_gamma,T_analytical,Y_analytical,LS_VAP_analytical
real(amrex_real) :: vel_analytical(SDIM)
real(amrex_real) :: pres_analytical
real(amrex_real) :: vel_compute(SDIM)
real(amrex_real) :: LS_compute
real(amrex_real) :: TEMPERATURE_compute
real(amrex_real) :: Y_compute
integer :: im_crit,tcomp
real(amrex_real) :: interface_thick_rad

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid

if (nsum1.gt.0) then
if (isweep.eq.0) then
 do dir=1,nsum1
  increment_out1(dir)=zero
 enddo
 if (axis_dir.eq.0) then
  if (vinletgas.eq.zero) then
   if (nsum1.eq.3+SDIM) then
    do dir=1,SDIM
     xlocal(dir)=GRID_DATA_IN%xsten(0,dir)
    enddo
    call drop_analytical_solution(GRID_DATA_IN%time, &
     xlocal,D_gamma,T_analytical,Y_analytical, &
     vel_analytical,LS_VAP_analytical,pres_analytical)
    LS_compute=GRID_DATA_IN%lsfab(D_DECL(i,j,k),2)
    im_crit=2
    tcomp=(im_crit-1)*num_state_material+ENUM_TEMPERATUREVAR+1
    TEMPERATURE_compute=GRID_DATA_IN%den(D_DECL(i,j,k),tcomp)
    Y_compute=GRID_DATA_IN%den(D_DECL(i,j,k),tcomp+1)
    do dir=1,SDIM
     vel_compute(dir)=GRID_DATA_IN%vel(D_DECL(i,j,k),dir)
    enddo
    interface_thick_rad=two*GRID_DATA_IN%dx(1)
 
    if (abs(LS_VAP_analytical).lt.interface_thick_rad) then
     increment_out1(1)= &
       GRID_DATA_IN%volgrid*abs(LS_compute-LS_VAP_analytical)/ &
       (two*interface_thick_rad)
    else
     increment_out1(1)=zero
    endif
    if (LS_VAP_analytical.gt.zero) then
     increment_out1(2)=GRID_DATA_IN%volgrid* &
        abs(TEMPERATURE_compute-T_analytical)
     increment_out1(3)=GRID_DATA_IN%volgrid* &
          abs(Y_compute-Y_analytical)
     do dir=1,SDIM
      increment_out1(3+dir)=GRID_DATA_IN%volgrid* &
          abs(vel_compute(dir)-vel_analytical(dir))
     enddo
    else if (LS_VAP_analytical.le.zero) then
     ! do nothing
    else
     print *,"LS_VAP_analytical invalid"
     stop
    endif
   else
    print *,"nsum1 invalid"
    stop
   endif
  else 
   print *,"expecting vinletgas==0.0"
   stop
  endif
 else
  print *,"expecting axis_dir==0"
  stop
 endif
else if (isweep.eq.1) then
 ! do nothing
else 
 print *,"isweep invalid"
 stop
endif

else if (nsum1.eq.0) then
 ! do nothing
else
 print *,"nsum1 invalid"
 print *,"nsum1 ",nsum1
 stop
endif

end subroutine DROP_IN_SHEAR_SUMINT

end module DROP_IN_SHEAR_module
