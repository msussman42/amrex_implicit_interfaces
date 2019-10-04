
subroutine H_epsilon(x,epsilon,H_epsilon_output)

  use amrex_fort_module, only : amrex_real
  implicit none
  real(amrex_real), intent(in) :: x
  real(amrex_real), intent(in) :: epsilon
  real(amrex_real), intent(out) :: H_epsilon_output

  H_epsilon_output=1.0d0/(1.0d0+exp(-x/epsilon))

  return
  end subroutine H_epsilon

  subroutine custom_flux(flux,uleft,uright,Delta,a,flux_type)
  use amrex_fort_module, only : amrex_real
  implicit none
  
  real(amrex_real), intent(out) :: flux
  real(amrex_real), intent(in) :: uleft
  real(amrex_real), intent(in) :: uright
  real(amrex_real), intent(in) :: Delta
  real(amrex_real), intent(in) :: a
  real(amrex_real) :: eps1,eps2,eps3,eps4

   ! u_t = div(F(u)
   ! flux_type==0  F(u)=grad u
   ! flux_type==1  F(u)=-a_vector u
   ! flux_type==2  F(u)=-a_vector u^2/2
  integer, intent(in) :: flux_type

  real(amrex_real) :: du,sleft,sright,H1,H2,H3,H4

  eps1=1.0D-3
  eps2=eps1
  eps3=eps2
  eps4=eps3

  if (flux_type.eq.0) then
   flux=(uright-uleft)/Delta
  else if (flux_type.eq.1) then
   call H_epsilon(a,eps1,H1)
   flux=-( H1*uleft*a + (1.0d0-H1)*uright*a )
  else if (flux_type.eq.2) then
   du=a*(uleft-uright)
   call H_epsilon(du,eps1,H1)
   sleft=H1*0.5d0*a*(uleft+uright)+(1.0d0-H1)*a*uleft
   call H_epsilon(du,eps2,H2)
   sright=H2*0.5d0*a*(uleft+uright)+(1.0d0-H2)*a*uright
   call H_epsilon(sleft,eps3,H3)
   call H_epsilon(sright,eps4,H4)
   flux=-( H3*a*uleft*uleft*0.5d0+ &
           (1.0d0-H3)*( (1.0d0-H4)*a*uright*uright*0.5d0 ) )
  else
   print *,"flux_type invalid"
   stop
  endif

  return
  end subroutine custom_flux
   
subroutine compute_flux (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi, &
                         dx,bc_vector,bc_value, &
                         dirichlet_condition, &
                         neumann_condition, &
                         periodic_condition, &
                         a_vector, &
                         flux_type) bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  integer, intent(in) :: bc_vector(4)
  real(amrex_real), intent(in) :: bc_value(4)
   ! bc_vector can be one of three values:
   ! "dirichlet_condition", "neumann_condition", "periodic_condition"
   ! bc_vector(1) => xlo wall  (dir=0 side=0)
   ! bc_vector(2) => xhi wall  (dir=0 side=1)
   ! bc_vector(3) => ylo wall  (dir=1 side=0)
   ! bc_vector(4) => yhi wall  (dir=1 side=1)
  integer, intent(in) :: dirichlet_condition
  integer, intent(in) :: neumann_condition
  integer, intent(in) :: periodic_condition
  real(amrex_real), intent(in) :: a_vector(2)
   ! u_t = div F(u)
   ! flux_type==0  F(u)=grad u
   ! flux_type==1  F(u)=-a_vector u
   ! flux_type==2  F(u)=-a_vector u^2/2
  integer, intent(in) :: flux_type

  ! local variables
  integer i,j
  real(amrex_real) :: phi_left,phi_right,Delta

  ! x-fluxes
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)+1
     phi_left=phi(i-1,j)
     phi_right=phi(i,j)
     Delta=dx(1)

     if (i.eq.lo(1)) then
      if (bc_vector(1).eq.dirichlet_condition) then
       Delta=0.5d0*dx(1)
       phi_left=bc_value(1)
      else if  (bc_vector(1).eq.neumann_condition) then
       Delta=0.5d0*dx(1)
        ! (phi_right-phi_left)/Delta=-bc_value(1)
        ! phi_left=Delta*bc_value(1)+phi_right
       phi_left=Delta*bc_value(1)+phi_right
      else if (bc_vector(1).eq.periodic_condition) then
       ! do nothing
      else
       print *,"bc_vector(1) invalid"
       stop
      endif
     endif

     if (i.eq.hi(1)+1) then
      if (bc_vector(2).eq.dirichlet_condition) then
       Delta=0.5d0*dx(1)
       phi_right=bc_value(2)
      else if  (bc_vector(2).eq.neumann_condition) then
       Delta=0.5d0*dx(1)
        ! (phi_right-phi_left)/Delta=bc_value(2)
        ! phi_right=Delta*bc_value(2)+phi_left
       phi_right=Delta*bc_value(2)+phi_left
      else if (bc_vector(2).eq.periodic_condition) then
       ! do nothing
      else
       print *,"bc_vector(2) invalid"
       stop
      endif
     endif

     call custom_flux(fluxx(i,j),phi_left,phi_right,Delta, &
       a_vector(1),flux_type)

  end do
  end do

  ! y-fluxes
  do j = lo(2), hi(2)+1
  do i = lo(1), hi(1)

     phi_left=phi(i,j-1)
     phi_right=phi(i,j)
     Delta=dx(2)

     if (j.eq.lo(2)) then
      if (bc_vector(3).eq.dirichlet_condition) then
       Delta=0.5d0*dx(2)
       phi_left=bc_value(3)
      else if  (bc_vector(3).eq.neumann_condition) then
       Delta=0.5d0*dx(2)
        ! (phi_right-phi_left)/Delta=-bc_value(3)
        ! phi_left=Delta*bc_value(1)+phi_right
       phi_left=Delta*bc_value(3)+phi_right
      else if (bc_vector(3).eq.periodic_condition) then
       ! do nothing
      else
       print *,"bc_vector(3) invalid"
       stop
      endif
     endif

     if (i.eq.hi(2)+1) then
      if (bc_vector(4).eq.dirichlet_condition) then
       Delta=0.5d0*dx(2)
       phi_right=bc_value(4)
      else if  (bc_vector(4).eq.neumann_condition) then
       Delta=0.5d0*dx(2)
        ! (phi_right-phi_left)/Delta=bc_value(4)
        ! phi_right=Delta*bc_value(4)+phi_left
       phi_right=Delta*bc_value(4)+phi_left
      else if (bc_vector(4).eq.periodic_condition) then
       ! do nothing
      else
       print *,"bc_vector(4) invalid"
       stop
      endif
     endif

     call custom_flux(fluxy(i,j),phi_left,phi_right,Delta, &
       a_vector(2),flux_type)

  end do
  end do

end subroutine compute_flux


subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, &
                       dx, dt) bind(C, name="update_phi")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2))
  real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), intent(in)    :: dt

  ! local variables
  integer i,j
  real(amrex_real) :: dtdx(2)

  dtdx = dt/dx

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)

     phinew(i,j) = phiold(i,j) &
          + dtdx(1) * (fluxx(i+1,j  ) - fluxx(i,j)) &
          + dtdx(2) * (fluxy(i  ,j+1) - fluxy(i,j))

  end do
  end do

end subroutine update_phi
