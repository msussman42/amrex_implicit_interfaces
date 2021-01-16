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

! probtype==424
module DROP_IN_SHEAR_module

implicit none                   

REAL_T :: DEF_VAPOR_GAMMA
REAL_T :: LOCAL_R_Palmore_Desjardins
REAL_T :: den_G,C_pG,k_G,lambda,T_inf,T_sat,L_V,D_G,Y_inf
REAL_T :: den_L
REAL_T :: WV,WA,Le
REAL_T :: D_not,B_M,Sh
REAL_T :: T_gamma,X_gamma,Y_gamma

contains

subroutine f_mdot(T_gamma_parm,f_out)
use global_utility_module
IMPLICIT NONE

REAL_T, intent(in) :: T_gamma_parm
REAL_T, intent(out) :: f_out
REAL_T :: X_gamma_loc,Y_gamma_loc

call X_from_Tgamma(X_gamma_loc,T_gamma_parm,T_sat,L_V, &
 LOCAL_R_Palmore_Desjardins,WV)
call massfrac_from_volfrac(X_gamma_loc,Y_gamma_loc,WA,WV)

if ((Y_gamma_loc.ge.zero).and.(Y_gamma_loc.lt.one)) then
 if (L_V.gt.zero) then
  if (C_pG.gt.zero) then
   if ((Y_inf.ge.zero).and.(Y_inf.lt.one)) then
    if ((T_inf.gt.zero).and.(T_inf.lt.T_sat)) then
     f_out=((Y_inf-one)/(Y_gamma_loc-one))**Le - one - &
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

REAL_T :: a,b,c,f_a,f_b,f_c
INTEGER_T :: iter

DEF_VAPOR_GAMMA =  1.666666667D0

! ergs/(mol kelvin)  (same as default for NavierStokes::R_Palmore_Desjardins)
LOCAL_R_Palmore_Desjardins=8.31446261815324D+7

den_L = fort_denconst(1)
den_G = fort_denconst(2)
C_pG = fort_stiffCP(2)
k_G = fort_heatviscconst(2)
lambda=k_G/(den_G*C_pG)
T_inf = fort_tempconst(2)
T_sat = fort_saturation_temp(1)
L_V = fort_latent_heat(1)
D_G = fort_speciesviscconst(2)
Y_inf=fort_speciesconst(2)
WV=fort_species_molar_mass(1)  !num_species components
WA=fort_molar_mass(2)
! T_inf C_pG / L = 300.5 K * 1e+7 (erg/(g K)) / ((2.26e+10) erg/g)
!   = 0.133
! e.g. Le=0.1 cm^2/s * 0.001 g/cm^3 * 4.1855e+7 erg/(g K) / 
!         (0.024e+5 g cm/(s^3 K)) = 1.74
! erg=g cm^2/s^2
! cm^2/s  * g/cm^3  * g cm^2 / s^2 * (1/(g K)) * (s^3 K)/(g cm) = dimensionless
!
Le=D_G*den_G*C_pG/k_G

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
do while ((f_a.gt.zero).and.(a.lt.b))
 a=a+0.001*T_sat
 call f_mdot(a,f_a)
enddo
if (f_a*f_b.le.zero) then
 do iter=1,100
  c=0.5d0*(a+b)
  call f_mdot(c,f_c)
  if (f_a*f_c.le.zero) then
   b=c
   f_b=f_c
  else if (f_b*f_c.le.zero) then
   a=c;
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
 LOCAL_R_Palmore_Desjardins,WV)
call massfrac_from_volfrac(X_gamma,Y_gamma,WA,WV)

B_M=(Y_gamma-Y_inf)/(one-Y_gamma)
D_not=two*radblob
if (B_M.gt.zero) then
 Sh=two*log(one+B_M)/B_M
else
 print *,"B_M must be positive"
 stop
endif

return
end subroutine INIT_DROP_IN_SHEAR_MODULE

! dir=velocity component
subroutine DROP_IN_SHEAR_CFL_HELPER(time,dir,uu,dx)
use probcommon_module
implicit none
INTEGER_T, intent(in) :: dir
REAL_T, intent(in) :: time
REAL_T, intent(inout) :: uu
REAL_T, intent(in) :: dx(SDIM)

REAL_T utest


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
REAL_T, intent(in), dimension(SDIM) :: x !spatial coordinates
INTEGER_T, intent(in) :: im
REAL_T, intent(out) :: dist
INTEGER_T :: nmat

nmat=num_materials
if (nmat.lt.1) then
 print *,"nmat invalid in soliddist"
 stop
endif

if ((im.lt.1).or.(im.gt.nmat)) then
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

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(out) :: LS(nmat)
INTEGER_T :: im
INTEGER_T :: im_solid_materialdist
REAL_T :: initial_time

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  im_solid_materialdist=im_solid_primary()
  initial_time=zero

   ! TODO (supermesh prototype code first)
   !       (i) quadratic least squares for heat flux and mass fraction
   !           flux
   !       (ii) twin material for capturing thermal layer. 
   ! for evaporation and boiling JCP paper:
   ! AMR, MOF, Space-time spectral
   ! 1. evaporating drop in shear
   ! 2. forced convective boiling from heated cylinder/sphere
   ! 3. i) planar boiling/evaporation front
   !    ii) expansing vapor bubble in superheated liquid.
   !    iii) nucleate boiling in either quiescent or flowing liquid.
   !    iv) pool boiling in quiescent or flowing liquid 
   ! 4. verification of space-time spectral accuracy for a problem
   !    using the Boussinesq approximation. (no phase change)
   ! for freezing and melting JCP paper?
   ! 1. static freezing test
   ! 2. ice cube melting on top of water on top of substrate.
   ! 3. thermal spray - convective heat transfer melting and solidification.
   ! 4. same as (3) for the evaporation and boiling JCP paper.
  if (probtype.eq.424) then

   do im=1,nmat
    if (FSI_flag(im).eq.1) then
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

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, intent(in) :: velsolid_flag
REAL_T :: vert_lo,vert_hi
REAL_T :: D_gamma,T_analytical,Y_analytical,LS_analytical

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
   call drop_analytical_solution(t,x,D_gamma,T_analytical, &
      Y_analytical,VEL,LS_analytical)
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
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: PRES

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
PRES=zero

return 
end subroutine DROP_IN_SHEAR_PRES



subroutine DROP_IN_SHEAR_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, intent(in) :: nmat
INTEGER_T, intent(in) :: nstate_mat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: STATE(nmat*nstate_mat)
INTEGER_T im,ibase,n
REAL_T :: D_gamma,T_analytical,Y_analytical,LS_analytical
REAL_T :: VEL(SDIM)

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
  STATE(ibase+1)=fort_denconst(im) ! density prescribed in the inputs file.
  if (t.eq.zero) then
   STATE(ibase+2)=fort_initial_temperature(im) !initial temperature in inputs
  else if (t.gt.zero) then
   STATE(ibase+2)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif

   ! initial species in inputs?
  do n=1,num_species_var
   STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo

 enddo ! im=1..num_materials

 if (axis_dir.eq.0) then
  if (vinletgas.eq.zero) then
   call drop_analytical_solution(t,x,D_gamma,T_analytical, &
      Y_analytical,VEL,LS_analytical)
   do im=1,num_materials
    ibase=(im-1)*num_state_material
    STATE(ibase+2)=T_analytical
    STATE(ibase+3)=Y_analytical
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

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(inout) :: LS(nmat)
REAL_T, intent(in) :: LS_in(nmat)
INTEGER_T, intent(in) :: dir,side
REAL_T, intent(in) ::  dx(SDIM)

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

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(inout) :: VEL
REAL_T, intent(in) :: VEL_in
INTEGER_T, intent(in) :: veldir,dir,side
REAL_T, intent(in) :: dx(SDIM)
REAL_T local_VEL(SDIM)
INTEGER_T velsolid_flag

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

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(inout) :: PRES
REAL_T, intent(in) :: PRES_in
INTEGER_T, intent(in) :: dir,side
REAL_T, intent(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

 call DROP_IN_SHEAR_PRES(xghost,t,LS,PRES,nmat)


return
end subroutine DROP_IN_SHEAR_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine DROP_IN_SHEAR_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T :: local_STATE(nmat*num_state_material)
REAL_T, intent(inout) :: STATE
REAL_T, intent(inout) :: STATE_merge
REAL_T, intent(in) :: STATE_in
INTEGER_T, intent(in) :: dir,side
REAL_T, intent(in) :: dx(SDIM)
INTEGER_T, intent(in) :: istate,im
INTEGER_T ibase,im_crit,im_loop
INTEGER_T local_bcflag

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
 im_crit=1
 do im_loop=2,num_materials
  if (LS(im_loop).gt.LS(im_crit)) then
   im_crit=im_loop
  endif
 enddo
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

INTEGER_T, intent(in) :: nmat
INTEGER_T, intent(in) :: im
REAL_T, intent(in) :: VFRAC(nmat)
REAL_T, intent(in) :: time
INTEGER_T, intent(in) :: nhalf
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
REAL_T, intent(in) :: temp(nmat)
REAL_T, intent(in) :: den(nmat)
REAL_T, intent(in) :: CV(nmat)
REAL_T, intent(in) :: dt
REAL_T, intent(out) :: heat_source

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

INTEGER_T, intent(in) :: nhalf
REAL_T, dimension(-nhalf:nhalf,SDIM), intent(in) :: xsten
REAL_T, intent(in) :: time
REAL_T, intent(in) :: dt
REAL_T, intent(out) :: heat_flux
INTEGER_T, intent(out) :: heat_dir
INTEGER_T, intent(out) :: heat_side

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

REAL_T, intent(in) :: dx(SDIM)
INTEGER_T, intent(in) :: veldir
REAL_T, intent(inout) :: heatcoeff

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
REAL_T, intent(inout) :: local_buffer(2*SDIM)
REAL_T, intent(in)    :: problen(SDIM)

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
INTEGER_T, intent(in) :: nhalf
REAL_T, dimension(-nhalf:nhalf,SDIM), intent(in) :: xsten
INTEGER_T, intent(inout) :: make_seed
type(nucleation_parm_type_input), intent(in) :: nucleate_in
REAL_T :: LL

LL=nucleate_in%LL

if (probtype.eq.424) then
 ! do nothing
else
 print *,"expecting probtype==424"
 stop
endif

return
end subroutine DROP_IN_SHEAR_nucleation

subroutine drop_analytical_solution(time,x,D_gamma,T,Y,VEL,LS_VAP)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: time
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(out) :: D_gamma
REAL_T, intent(out) :: T
REAL_T, intent(out) :: Y
REAL_T, intent(out) :: LS_VAP
REAL_T, intent(out) :: VEL(SDIM)
REAL_T :: rr,mdot,vel_r

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

D_Gamma=D_not**2-four*den_G*D_G*log(one+B_M)*t/den_L
if (D_Gamma.gt.D_not**2) then
 D_Gamma=sqrt(D_gamma)
else
 print *,"D_gamma invalid"
 stop
endif
mdot=Pi*D_gamma*den_G*D_G*Sh*B_M

LS_VAP=rr-half*D_gamma

if (LS_VAP.le.zero) then
 VEL(1)=zero
 VEL(2)=zero
 VEL(SDIM)=zero
 Y=Y_Gamma
 T=T_Gamma
else if (LS_VAP.ge.zero) then
 vel_r = mdot/(4.0*Pi*rr*rr*den_G)
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
 T=T_Gamma+L_V/C_pG+(T_inf-T_Gamma-L_V/C_pG)* &
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

! This routine is called from FORT_SUMMASS
! set ns.ncomp_sum_int_user=
subroutine DROP_IN_SHEAR_SUMINT(GRID_DATA_IN,increment_out,nsum)
use probcommon_module_types
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nsum
type(user_defined_sum_int_type), intent(in) :: GRID_DATA_IN
REAL_T, intent(out) :: increment_out(nsum)
INTEGER_T :: i,j,k,dir

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid

if (nsum.gt.0) then
 do dir=1,nsum
  increment_out(dir)=zero
 enddo
 print *,"nothing here yet"
 stop
else if (nsum.eq.0) then
 ! do nothing
else
 print *,"nsum invalid"
 print *,"nsum ",nsum
 stop
endif

end subroutine DROP_IN_SHEAR_SUMINT

end module DROP_IN_SHEAR_module
