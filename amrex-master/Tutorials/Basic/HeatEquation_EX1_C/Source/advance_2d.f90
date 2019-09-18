
subroutine compute_flux (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi, &
                         dx,bc_vector,bc_value, &
                         dirichlet_condition, &
                         neumann_condition, &
                         periodic_condition) bind(C, name="compute_flux")

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

  ! local variables
  integer i,j

  ! x-fluxes
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)+1
     fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
     if ((bc_vector(1).eq.dirichlet_condition).and. &
         (i.eq.lo(1))) then
      fluxx(i,j)= (phi(i,j)-bc_value(1))/(0.5d0*dx(1))
     endif
     if ((bc_vector(2).eq.dirichlet_condition).and. &
         (i.eq.hi(1)+1)) then
      fluxx(i,j)= (bc_value(2)-phi(i-1,j))/(0.5d0*dx(1))
     endif
  end do
  end do

  ! y-fluxes
  do j = lo(2), hi(2)+1
  do i = lo(1), hi(1)
     fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
     if ((bc_vector(3).eq.dirichlet_condition).and. &
         (j.eq.lo(2))) then
      fluxy(i,j)= (phi(i,j)-bc_value(3))/(0.5d0*dx(2))
     endif
     if ((bc_vector(4).eq.dirichlet_condition).and. &
         (j.eq.hi(2)+1)) then
      fluxy(i,j)= (bc_value(4)-phi(i,j-1))/(0.5d0*dx(2))
     endif
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
