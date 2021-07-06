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

! probtype==2002 (see run2d/inputs.SIMPLE_PALMORE_DESJARDINS)
module SIMPLE_PALMORE_DESJARDINS_module

implicit none                   

REAL_T :: DEF_VAPOR_GAMMA
REAL_T :: l_verification

contains

  ! do any initial preparation needed
subroutine INIT_SIMPLE_PALMORE_DESJARDINS_MODULE()
use probcommon_module
IMPLICIT NONE
REAL_T dummy_x,dummy_t,dummy_TorY,dummy_LS
REAL_T t_physical_init
INTEGER_T use_T
REAL_T C_pG,k_G,den_G,lambda

  DEF_VAPOR_GAMMA =  1.666666667D0

  call SIMPLE_PALMORE_DESJARDINS_GetDiffusionLayer(l_verification)

  use_t=1
  dummy_x=zero
  dummy_t=zero
  call SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC( &
   dummy_x,dummy_t,use_T,dummy_TorY, &
   dummy_LS,t_physical_init)

  C_pG = fort_stiffCP(2)
  k_G = fort_heatviscconst(2)
  den_G = fort_denconst(2)
  lambda=k_G/(den_G*C_pG)

  print *,"expected interface location is: 2 l sqrt(lambda*(tphys+t))-xblob2"
  print *,"l=",l_verification
  print *,"lambda=",lambda
  print *,"tphys=",t_physical_init
  print *,"xblob2=",xblob2

return
end subroutine INIT_SIMPLE_PALMORE_DESJARDINS_MODULE

! gas on the left (material 2)
! liquid on the right (material 1)
! xblob=initial interface location
subroutine SIMPLE_PALMORE_DESJARDINS_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(out) :: LS(nmat)
REAL_T :: TEMPERATURE_analytical
REAL_T :: LS_analytical
REAL_T :: t_physical_init
INTEGER_T :: use_T

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

if ((num_materials.eq.2).and.(probtype.eq.2002)) then

 use_T=1
 call SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC( &
   x(1),t,use_T,TEMPERATURE_analytical, &
   LS_analytical,t_physical_init)

 LS(1)=LS_analytical  ! liquid
 LS(2)=-LS(1)      ! gas

else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine SIMPLE_PALMORE_DESJARDINS_LS

! initial velocity is zero
subroutine SIMPLE_PALMORE_DESJARDINS_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, intent(in) :: velsolid_flag

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

return 
end subroutine SIMPLE_PALMORE_DESJARDINS_VEL


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine SIMPLE_PALMORE_DESJARDINS_PRES(x,t,LS,PRES,nmat)
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
end subroutine SIMPLE_PALMORE_DESJARDINS_PRES


subroutine SIMPLE_PALMORE_DESJARDINS_DiffusionLayer(l,f) 
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE

 REAL_T, intent(in) :: l !diffusion layer value
 REAL_T, intent(out) :: f
 
 REAL_T :: T_inf, T_gamma, L_V, C_pG, erf_result, T_sat
 REAL_T :: k_G, den_G, D_G
 REAL_T :: Y_gamma,Y_G,WV,WA,R,X_gamma,Y_gamma_test,X_gamma_test
 REAL_T :: T_gamma_min
 REAL_T :: T_gamma_max
 REAL_T :: T_gamma_test
 REAL_T :: lambda
 INTEGER_T :: JINT
 
 T_inf = fort_tempconst(2)
 T_gamma = fort_tempconst(1)
 T_sat = fort_saturation_temp(1)
 L_V = fort_latent_heat(1)
 C_pG = fort_stiffCP(2)
 k_G = fort_heatviscconst(2)
 den_G = fort_denconst(2)
 D_G = fort_speciesviscconst(2)
 Y_gamma=fort_speciesconst(1)
 Y_G=fort_speciesconst(2)
 WV=fort_species_molar_mass(1)
 WA=fort_molar_mass(2)
 R=fort_R_Palmore_Desjardins
 call volfrac_from_massfrac(X_gamma,Y_gamma,WA,WV)
 call massfrac_from_volfrac(X_gamma,Y_gamma_test,WA,WV)
 if (abs(Y_gamma-Y_gamma_test).le.1.0D-8) then
  ! do nothing
 else
  print *,"Y_gamma_test invalid"
  stop
 endif

 T_gamma_min=0.0d0
 T_gamma_max=1.0D+20

 call Tgamma_from_TSAT_and_X(T_gamma_test,T_sat,X_gamma,L_V,R,WV, &
   T_gamma_min,T_gamma_max)
 if (abs(T_gamma-T_gamma_test).le.1.0D-8) then
  ! do nothing
 else
  print *,"T_gamma_test invalid1"
  print *,"T_gamma= ",T_gamma
  print *,"T_gamma_test= ",T_gamma_test
  print *,"T_sat= ",T_sat
  print *,"X_gamma= ",X_gamma
  print *,"L_V= ",L_V
  print *,"R= ",R
  print *,"WV= ",WV
  print *,"T_gamma_min ",T_gamma_min
  print *,"T_gamma_max ",T_gamma_max
  stop
 endif

 call X_from_Tgamma(X_gamma_test,T_gamma,T_sat,L_V,R,WV) 
 if (abs(X_gamma-X_gamma_test).le.1.0D-10) then
  ! do nothing
 else
  print *,"X_gamma_test invalid 2"
  print *,"X_gamma= ",X_gamma
  print *,"X_gamma_test= ",X_gamma_test
  print *,"T_gamma= ",T_gamma
  print *,"T_sat= ",T_sat
  print *,"L_V= ",L_V
  print *,"R= ",R
  print *,"WV= ",WV
  stop
 endif

 lambda=k_G/(den_G*C_pG)

   ! required that D_G=lambda
 if ((T_inf.gt.T_gamma).and. &
     (T_gamma.le.T_sat).and. &
     (l.ge.zero).and. &
     (Y_G.ge.zero).and. &
     (Y_G.le.one).and. &
     (abs(lambda-D_G).lt.1.0D-8)) then
  ! do nothing
 else
  print *,"T_inf, T_gamma, T_sat, l, or D_G invalid"
  print *,"T_inf=",T_inf
  print *,"T_gamma=",T_gamma
  print *,"T_sat=",T_sat
  print *,"l=",l
  print *,"Y_G=",Y_G
  print *,"lambda=",lambda
  print *,"D_G=",D_G
  stop
 endif
 JINT=0 ! JINT=2 => exp(l^2)erf(l) but xneg<l<xmax
 call calerf(l,erf_result,JINT)

 f=l*EXP(l**2)*erf_result - C_pG*(T_inf-T_gamma)/(sqrt(Pi)*L_V)

end subroutine SIMPLE_PALMORE_DESJARDINS_DiffusionLayer

subroutine SIMPLE_PALMORE_DESJARDINS_GetDiffusionLayer(l)
 !bisection method to find the diffusion layer value
 IMPLICIT NONE
 
 REAL_T, intent(out) :: l
 REAL_T :: a, b, c, fa, fb, fc
 INTEGER_T :: iter
 
 !endpoints
 a = 0.0d0
 b = 5.0d0
 call SIMPLE_PALMORE_DESJARDINS_DiffusionLayer(a,fa)
 call SIMPLE_PALMORE_DESJARDINS_DiffusionLayer(b,fb)
 if (fa*fb.le.0.0d0) then
  ! do nothing
 else
  print *,"a,b is not a bracketing interval"
  stop
 endif

 iter = 1
 do while (iter.LT.100)
  c = (a+b)/2.0d0
  call SIMPLE_PALMORE_DESJARDINS_DiffusionLayer(c,fc)
  iter = iter + 1
  call SIMPLE_PALMORE_DESJARDINS_DiffusionLayer(a,fa)
  if (fc*fa .GE. zero ) then
   a = c
  else if (fc*fa.LE.zero) then
   b = c
  else
   print *,"bracketing interval property has been lost"
   stop
  endif
 enddo ! iter.LT.100
 l=c
end subroutine SIMPLE_PALMORE_DESJARDINS_GetDiffusionLayer

! lambda=k/(rho C_p)
! xblob is the interface at the start of the computation.
! xblob2 is the physical position corresponding to x_compute=0
! x_0 is zero  
! xblob2>x_0
! xblob+xblob2 is the physical location of the interface at t_compute=0
subroutine SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC( &
  x, t, use_T, TorY, LS_exact,t_physical_init) 
 use global_utility_module
 use probcommon_module
 !returns either temperature or mass frac
 ! for mass_frac: T_inf = Y_inf, T_gamma = Y_gamma, l_Y = l, lambda = D
 IMPLICIT NONE
 
 REAL_T, intent(in) :: x, t
 INTEGER_T, intent(in) :: use_T
 REAL_T :: T_inf, T_gamma, T_sat, lambda
 REAL_T :: k_G, den_G, D_G
 REAL_T :: L_V,C_pG
 REAL_T, intent(out) :: TorY,LS_exact
 REAL_T, intent(out) :: t_physical_init
 INTEGER_T :: JINT
 REAL_T erf_result_x 
 REAL_T erf_result_l 
 REAL_T arg_x
 REAL_T x_gamma_physical
 REAL_T x_gamma_domain
 REAL_T X_gamma_test
 REAL_T Y_gamma_test
 REAL_T X_gamma
 REAL_T Y_gamma
 REAL_T Y_inf
 REAL_T Y_inf_test
 REAL_T WV,WA,R
 REAL_T :: T_gamma_min
 REAL_T :: T_gamma_max
 REAL_T :: T_gamma_test

 T_inf = fort_tempconst(2)
 T_gamma = fort_tempconst(1)
 T_sat = fort_saturation_temp(1)
 L_V = fort_latent_heat(1)
 C_pG = fort_stiffCP(2)
 k_G = fort_heatviscconst(2)
 den_G = fort_denconst(2)
 D_G = fort_speciesviscconst(2)
 Y_gamma=fort_speciesconst(1)  
 Y_inf=fort_speciesconst(2)

 WV=fort_species_molar_mass(1)
 WA=fort_molar_mass(2)
 R=fort_R_Palmore_Desjardins
 call volfrac_from_massfrac(X_gamma,Y_gamma,WA,WV)
 call massfrac_from_volfrac(X_gamma,Y_gamma_test,WA,WV)
 if (abs(Y_gamma-Y_gamma_test).le.1.0D-8) then
  ! do nothing
 else
  print *,"Y_gamma_test invalid"
  stop
 endif

 T_gamma_min=0.0d0
 T_gamma_max=1.0D+20

 call Tgamma_from_TSAT_and_X(T_gamma_test,T_sat,X_gamma,L_V,R,WV, &
   T_gamma_min,T_gamma_max)
 if (abs(T_gamma-T_gamma_test).le.1.0D-8) then
  ! do nothing
 else
  print *,"T_gamma_test invalid2"
  print *,"T_gamma= ",T_gamma
  print *,"T_gamma_test= ",T_gamma_test
  stop
 endif

 call X_from_Tgamma(X_gamma_test,T_gamma,T_sat,L_V,R,WV) 
 if (abs(X_gamma-X_gamma_test).le.1.0D-10) then
  ! do nothing
 else
  print *,"X_gamma_test invalid 1"
  print *,"X_gamma= ",X_gamma
  print *,"X_gamma_test= ",X_gamma_test
  print *,"T_gamma= ",T_gamma
  print *,"T_sat= ",T_sat
  print *,"L_V= ",L_V
  print *,"R= ",R
  print *,"WV= ",WV
  stop
 endif

 if ((xblob2.gt.zero).and. &
     (k_G.gt.zero).and. &
     (den_G.gt.zero).and. &
     (C_pG.gt.zero).and. &
     (t.ge.zero)) then
  ! do nothing
 else
  print *,"xblob2, k_G, den_G, C_pG, or t invalid"
  stop
 endif
 lambda=k_G/(den_G*C_pG)

 if (Y_gamma.eq.1.0d0) then
  if (Y_inf.eq.1.0d0) then
   T_gamma_test=T_sat
   Y_inf_test=1.0d0
  else
   print *,"Y_inf invalid"
   stop
  endif
 else if ((Y_gamma.ge.0.0d0).and. &
          (Y_gamma.lt.1.0d0)) then
  T_gamma_test=T_inf+(L_V/C_pG)*(Y_gamma-Y_inf)/(Y_gamma-1.0d0)
  Y_inf_test=Y_gamma-(T_gamma-T_inf)*(Y_gamma-1.0d0)*(C_pG/L_V)
 else
  print *,"Y_gamma invalid"
  stop
 endif

   ! required that D_G=lambda
 if ((T_inf.gt.T_gamma).and. &
     (Y_gamma.ge.Y_inf).and. &
     (T_gamma.le.T_sat).and. &
     (l_verification.ge.zero).and. &
     (abs(lambda-D_G).lt.1.0D-8).and. &
     (abs(T_gamma_test-T_gamma).lt.1.0D-8).and. &
     (abs(Y_inf_test-Y_inf).lt.1.0D-8)) then
  ! do nothing
 else
  print *,"T_inf, Y_inf, T_sat, T_gamma, l, D_G, or TY_eqn invalid"
  print *,"T_inf=",T_inf
  print *,"T_gamma=",T_gamma
  print *,"T_sat=",T_sat
  print *,"l_verification=",l_verification
  print *,"Y_gamma=",Y_gamma
  print *,"Y_inf=",Y_inf
  print *,"lambda=",lambda
  print *,"D_G=",D_G
  print *,"T_gamma_test=",T_gamma_test
  print *,"Y_inf_test=",Y_inf_test
  stop
 endif

  ! xgamma=x0+ 2 l sqrt(lambda t)
  ! (xgamma-x0)/(2l) = sqrt(lambda t)
  ! [(xgamma-x0)/(2l)]^2 = lambda t
  ! t=[(xgamma-x0)/(2l)]^2/lambda

 t_physical_init=((xblob+xblob2)/(2.0d0*l_verification))**2/lambda
 x_gamma_physical=2.0d0*l_verification*sqrt(lambda*(t_physical_init+t))

 x_gamma_domain=x_gamma_physical-xblob2

 LS_exact=x-x_gamma_domain

 if (x_gamma_physical.gt.xblob+xblob2-1.0D-8) then
  ! do nothing
 else
  print *,"x_gamma_physical invalid"
  stop
 endif

 JINT=0 ! JINT=0 => erf(l)  JINT=2 => exp(l^2)erf(l) but xneg<l<xmax
 call calerf(l_verification,erf_result_l,JINT)

 if (x.gt.x_gamma_domain) then
  if (use_T.eq.1) then
   TorY=T_gamma
  else if (use_T.eq.0) then
   TorY=Y_gamma
  else
   print *,"use_T invalid"
   stop
  endif
 else if (x.le.x_gamma_domain) then
  arg_x=(x+xblob2)/(2.0d0*SQRT(lambda*(t_physical_init+t)))
  call calerf(arg_x,erf_result_x,JINT)

  if (use_T.eq.1) then
   TorY = T_inf +(T_gamma-T_inf)*erf_result_x/erf_result_l
  else if (use_T.eq.0) then
   TorY = Y_inf +(Y_gamma-Y_inf)*erf_result_x/erf_result_l
  else
   print *,"use_T invalid"
   stop
  endif
 else
  print *,"x is corrupt"
  stop
 endif
 
end subroutine SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC

subroutine SIMPLE_PALMORE_DESJARDINS_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
INTEGER_T im,ibase,use_T
REAL_T LS_exact
REAL_T t_physical_init

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
if ((num_materials.eq.2).and. &
    (num_state_material.eq.3).and. & ! density, temperature, species
    (probtype.eq.2002)) then
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

   ! TEMPERATURE
  use_T=1
  call SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC( &
   x(1),t,use_T,STATE(ibase+2),LS_exact,t_physical_init)
   ! MASS FRACTION
  use_T=0
  call SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC( &
   x(1),t,use_T,STATE(ibase+3),LS_exact,t_physical_init)

  if (im.eq.1) then ! water
   state(ibase+2)=fort_tempconst(im)
   state(ibase+3)=fort_speciesconst(im)
  endif

 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine SIMPLE_PALMORE_DESJARDINS_STATE



subroutine SIMPLE_PALMORE_DESJARDINS_ASSIMILATE( &
     assimilate_in,assimilate_out, &
     i,j,k,cell_flag)
use probcommon_module
IMPLICIT NONE

type(assimilate_parm_type), intent(in) :: assimilate_in
type(assimilate_out_parm_type), intent(inout) :: assimilate_out
INTEGER_T, intent(in) :: i,j,k,cell_flag

INTEGER_T :: nmat,nstate,nstate_test
REAL_T :: x_exact,xcrit,tcrit
INTEGER_T :: use_T
INTEGER_T :: dir
INTEGER_T :: im
INTEGER_T :: ibase
REAL_T local_temp,local_massfrac,LS_exact,t_physical_init

nmat=assimilate_in%nmat
nstate=assimilate_in%nstate

nstate_test=(SDIM+1)+nmat*(num_state_material+ngeom_raw)+1
if (nstate.eq.nstate_test) then
 ! do nothing
else
 print *,"nstate invalid"
 print *,"nstate=",nstate
 print *,"nstate_test=",nstate_test
 stop
endif

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((num_materials.eq.2).and. &
    (num_state_material.eq.3).and. & ! density, temperature, species
    (probtype.eq.2002)) then
 x_exact=probhix/8.0d0
 xcrit=assimilate_in%xsten(0,1)
 tcrit=assimilate_in%time  ! cur_time_slab
 if ((xcrit.le.x_exact).or. &
     (xcrit.ge.probhix-x_exact)) then
  if (cell_flag.eq.0) then ! MAC GRID X
   if (xcrit.le.x_exact) then
    assimilate_out%macx(D_DECL(i,j,k))=0.0d0
   endif
  else if (cell_flag.eq.1) then ! MAC GRID Y
   if (xcrit.le.x_exact) then
    assimilate_out%macy(D_DECL(i,j,k))=0.0d0
   endif
  else if ((cell_flag.eq.2).and.(SDIM.eq.3)) then ! MAC GRID Z
   if (xcrit.le.x_exact) then
    assimilate_out%macz(D_DECL(i,j,k))=0.0d0
   endif
  else if (cell_flag.eq.-1) then
    ! TEMPERATURE
    use_T=1
    call SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC( &
      xcrit,tcrit,use_T,local_temp,LS_exact,t_physical_init)
     ! MASS FRACTION
    use_T=0
    call SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC( &
      xcrit,tcrit,use_T,local_massfrac,LS_exact,t_physical_init)
    if (xcrit.le.x_exact) then
     do dir=1,SDIM
      assimilate_out%state(D_DECL(i,j,k),dir)=0.0d0
     enddo
    endif

    if (xcrit.le.x_exact) then
     do im=1,num_materials
      ibase=SDIM+1+(im-1)*num_state_material
      assimilate_out%state(D_DECL(i,j,k),ibase+2)=local_temp
      assimilate_out%state(D_DECL(i,j,k),ibase+3)=local_massfrac
     enddo
    endif
  else 
   print *,"cell_flag invalid"
   stop
  endif
 else if ((xcrit.ge.x_exact).and. &
          (xcrit.le.probhix-x_exact)) then
  ! do nothing
 else
  print *,"xcrit invalid"
  stop
 endif
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine SIMPLE_PALMORE_DESJARDINS_ASSIMILATE



 ! dir=1..sdim  side=1..2
subroutine SIMPLE_PALMORE_DESJARDINS_LS_BC(xwall,xghost,t,LS, &
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
 call SIMPLE_PALMORE_DESJARDINS_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine SIMPLE_PALMORE_DESJARDINS_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine SIMPLE_PALMORE_DESJARDINS_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
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
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call SIMPLE_PALMORE_DESJARDINS_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine SIMPLE_PALMORE_DESJARDINS_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine SIMPLE_PALMORE_DESJARDINS_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
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
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call SIMPLE_PALMORE_DESJARDINS_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine SIMPLE_PALMORE_DESJARDINS_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine SIMPLE_PALMORE_DESJARDINS_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
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
 call SIMPLE_PALMORE_DESJARDINS_STATE(xghost,t,LS,local_STATE,local_bcflag, &
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
else
 print *,"istate invalid"
 stop
endif

return
end subroutine SIMPLE_PALMORE_DESJARDINS_STATE_BC

subroutine SIMPLE_PALMORE_DESJARDINS_HEATSOURCE(im,VFRAC,time,x, &
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
 
 ! set a hot temperature in the liquid
if ((num_materials.eq.2).and.(probtype.eq.2002)) then

 heat_source=zero

else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine SIMPLE_PALMORE_DESJARDINS_HEATSOURCE

! This routine is called from FORT_SUMMASS
subroutine SIMPLE_PALMORE_DESJARDINS_SUMINT(GRID_DATA_IN,increment_out1, &
                increment_out2,nsum1,nsum2,isweep)
use probcommon_module_types
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nsum1,nsum2,isweep
type(user_defined_sum_int_type), intent(in) :: GRID_DATA_IN
REAL_T, intent(inout) :: increment_out1(nsum1)
REAL_T, intent(inout) :: increment_out2(nsum2)
INTEGER_T :: i,j,k
INTEGER_T :: dir
INTEGER_T :: im_crit
INTEGER_T :: tcomp
INTEGER_T :: use_T
REAL_T :: xlocal(SDIM)
REAL_T :: VOF_analytical
REAL_T :: VOF_compute
REAL_T :: x_analytical,x_left,x_right
REAL_T :: LS_analytical
REAL_T :: LS_compute
REAL_T :: TEMPERATURE_analytical
REAL_T :: Y_analytical
REAL_T :: TEMPERATURE_compute
REAL_T :: Y_compute
REAL_T :: interface_thick_rad
REAL_T :: t_physical_init

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid

if (nsum1.eq.6) then

if (isweep.eq.0) then

 do dir=1,SDIM
  xlocal(dir)=GRID_DATA_IN%xsten(0,dir)
 enddo
 x_left=GRID_DATA_IN%xsten(-1,1)
 x_right=GRID_DATA_IN%xsten(1,1)

 use_T=1
 call SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC( &
   xlocal(1),GRID_DATA_IN%time,use_T,TEMPERATURE_analytical, &
   LS_analytical,t_physical_init)
 use_T=0
 call SIMPLE_PALMORE_DESJARDINS_TEMPorMASSFRAC( &
   xlocal(1),GRID_DATA_IN%time,use_T,Y_analytical, &
   LS_analytical,t_physical_init)

 x_analytical=xlocal(1)-LS_analytical

 interface_thick_rad=two*GRID_DATA_IN%dx(1)

  ! slopes: 1..nmat*ngeom_recon
 VOF_compute=GRID_DATA_IN%slopes(D_DECL(i,j,k),1)
 if (x_right.gt.x_left) then
  if (x_analytical.le.x_left) then
   VOF_analytical=one
  else if (x_analytical.ge.x_right) then
   VOF_analytical=zero
  else if ((x_analytical.ge.x_left).and. &
           (x_analytical.le.x_right)) then
   VOF_analytical=(x_right-x_analytical)/(x_right-x_left)
  else
   print *,"x_analytical invalid"
   stop
  endif
 else
  print *,"x_right or x_left invalid"
  stop
 endif
 if ((VOF_analytical.ge.zero).and. &
     (VOF_analytical.le.one)) then
  ! symmetric difference error is measured too in FORT_SUMMASS.
  increment_out1(1)=GRID_DATA_IN%volgrid*abs(VOF_compute-VOF_analytical)
 else
  print *,"VOF_analytical invalid"
  stop
 endif
  
 LS_compute=GRID_DATA_IN%lsfab(D_DECL(i,j,k),1)
 if (abs(LS_analytical).lt.interface_thick_rad) then
  increment_out1(2)=GRID_DATA_IN%volgrid*abs(LS_compute-LS_analytical)/ &
    (two*interface_thick_rad)
 else
  increment_out1(2)=zero
 endif
 increment_out1(3)=zero
 increment_out1(4)=zero
 increment_out1(5)=zero
 increment_out1(6)=zero
 if (VOF_analytical.le.VOFTOL) then
  im_crit=2
  tcomp=(im_crit-1)*num_state_material+2
  TEMPERATURE_compute=GRID_DATA_IN%den(D_DECL(i,j,k),tcomp)
  increment_out1(3)=GRID_DATA_IN%volgrid* &
          abs(TEMPERATURE_compute-TEMPERATURE_analytical)
  Y_compute=GRID_DATA_IN%den(D_DECL(i,j,k),tcomp+1)
  increment_out1(4)=GRID_DATA_IN%volgrid* &
          abs(Y_compute-Y_analytical)
 else if (VOF_analytical.ge.one-VOFTOL) then
  im_crit=1
  tcomp=(im_crit-1)*num_state_material+2
  TEMPERATURE_compute=GRID_DATA_IN%den(D_DECL(i,j,k),tcomp)
  increment_out1(5)=GRID_DATA_IN%volgrid* &
          abs(TEMPERATURE_compute-TEMPERATURE_analytical)
  Y_compute=GRID_DATA_IN%den(D_DECL(i,j,k),tcomp+1)
  increment_out1(6)=GRID_DATA_IN%volgrid* &
          abs(Y_compute-Y_analytical)
 else if ((VOF_analytical.ge.VOFTOL).and. &
          (VOF_analytical.le.one-VOFTOL)) then
  ! do nothing
 else
  print *,"VOF_analytical invalid"
  stop
 endif

else if (isweep.eq.1) then
  ! do nothing
else
  print *,"isweep invalid"
  stop
endif

else
 print *,"nsum1 invalid"
 print *,"nsum1 ",nsum1
 stop
endif

end subroutine SIMPLE_PALMORE_DESJARDINS_SUMINT



end module SIMPLE_PALMORE_DESJARDINS_module
