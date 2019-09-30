subroutine init_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi) bind(C, name="init_phi")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 

  integer          :: i,j
  double precision :: x,y,r2

  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

!        r2 = ((x-0.25d0)**2 + (y-0.25d0)**2) / 0.01d0
!         r2 = (((x-0.25d0)/4.0d0)**2 + (y-0.25d0)**2)/0.01d0
         r2 = x
        phi(i,j) = 1.d0 + exp(-r2)

     end do
  end do

end subroutine init_phi


subroutine compute_dt(lo, hi, phi, philo, phihi, a_vector, flux_type, dt, dx, prob_lo, prob_hi) bind(C, name="compute_dt")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
  real(amrex_real), intent(in   ) :: a_vector(2)
  real(amrex_real), intent(inout) :: dt
  integer, intent(in) :: flux_type
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 

  integer          :: i,j
  double precision :: x,y,r2
  double precision :: local_dt
  double precision :: a_magnitude
  double precision :: u

  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

        u=phi(i,j)

        if (flux_type.eq.0) then ! F(u)=-grad u
         local_dt=0.9d0*dx(1)*dx(1)/(2.0d0*2)
        else if (flux_type.eq.1) then ! F(u)= a_vector u

         a_magnitude=sqrt(a_vector(1)**2+a_vector(2)**2)
         if (a_magnitude>0.0d0) then
          local_dt=0.9d0*dx(1)/a_magnitude
         else
          print *,"a_magnitude must be positive"
          stop
         endif
         
        else if (flux_type.eq.2) then ! F(u)= a_vector u^2/2

         a_magnitude=sqrt(a_vector(1)**2+a_vector(2)**2)*abs(u)
         if (a_magnitude>0.0d0) then
          local_dt=0.9d0*dx(1)/a_magnitude
         else
          local_dt=1.0D+20
         endif
         
        else
         print *,"flux_type invalid"
         stop
        endif

        if (local_dt.lt.dt) then
         dt=local_dt
        endif

     end do
  end do

end subroutine compute_dt

