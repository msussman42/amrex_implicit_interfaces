subroutine init_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi) bind(C, name="init_phi")

  use amrex_fort_module, only : amrex_real, amrex_spacedim
  implicit none

  integer, intent(in) :: lo(3), hi(3), philo(3), phihi(3)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
  real(amrex_real), intent(in   ) :: dx(3) 
  real(amrex_real), intent(in   ) :: prob_lo(3) 
  real(amrex_real), intent(in   ) :: prob_hi(3) 

  integer          :: i,j,k
  double precision :: x,y,z

  do k = lo(3), hi(3)
     z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
     do j = lo(2), hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i = lo(1), hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

            if (AMREX_SPACEDIM.eq.2) then
				phi(i,j,k) = 0.15d0 - sqrt((x-0.5d0)**2 + (y-0.75d0)**2) ! TODO: 2d not handled correctly
			else 
				phi(i,j,k) = 0.15d0 - sqrt((x-0.35d0)**2 + (y-0.35d0)**2 + (z-0.35d0)**2)
                !phi(i,j,k) = 0.15d0 - sqrt((x-0.5d0)**2 + (y-0.75d0)**2)
			endif

        end do
     end do
  end do
  


end subroutine init_phi

subroutine init_vel(lo, hi, domlo, domhi, &
					vel_u, ulo, uhi, vel_v, vlo, vhi, vel_w, wlo, whi, &
					dx, prob_lo, prob_hi, time) bind(C, name="init_vel")

  use amrex_fort_module, only : amrex_real, amrex_spacedim
  implicit none

  integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3)
  integer, intent(in) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)
  real(amrex_real), intent(inout) :: vel_u( ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3))
  real(amrex_real), intent(inout) :: vel_v( vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3))
  real(amrex_real), intent(inout) :: vel_w( wlo(1):whi(1), wlo(2):whi(2), wlo(3):whi(3))
  real(amrex_real), intent(in	) :: dx(3)
  real(amrex_real), intent(in   ) :: prob_lo(3) 
  real(amrex_real), intent(in   ) :: prob_hi(3) 
  real(amrex_real), intent(in) :: time

  integer          :: i,j,k
  double precision :: x,y,z, PI

  PI=4.D0*DATAN(1.D0)
  
  ! u-vel
  do k = lo(3), hi(3)
	  z = prob_lo(3) + (dble(k)+0.0d0) * dx(3)
	  do j = lo(2), hi(2)
		  y = prob_lo(2) + (dble(j)+0.0d0) * dx(2)
		  do i = lo(1), hi(1)+1
			x = prob_lo(1) + (dble(i)+0.0d0) * dx(1)
			
			if (AMREX_SPACEDIM.eq.2) then
				vel_u(i,j,k) = sin(PI*x)**2 * sin(2*PI*y)
			else 
				vel_u(i,j,k) = 2*sin(PI*x)**2 * sin(2*PI*y) * sin(2*PI*z)
                !vel_u(i,j,k) = sin(PI*x)**2 * sin(2*PI*y)
                if (time .gt. 1.0) then
                    vel_u(i,j,k) = -2*sin(PI*x)**2 * sin(2*PI*y) * sin(2*PI*z)
                endif
			endif
			
			
		  end do
	  end do
  end do
  
  ! v-vel
  do k = lo(3), hi(3)
	  z = prob_lo(3) + (dble(k)+0.0d0) * dx(3)
	  do j = lo(2), hi(2)+1
		  y = prob_lo(2) + (dble(j)+0.0d0) * dx(2)
		  do i = lo(1), hi(1)
			x = prob_lo(1) + (dble(i)+0.0d0) * dx(1)
			
			if (AMREX_SPACEDIM.eq.2) then
				vel_v(i,j,k) = -sin(PI*y)**2 * sin(2*PI*x)
			else 
				vel_v(i,j,k) = -sin(PI*y)**2 * sin(2*PI*x) * sin(2*PI*z)
                !vel_v(i,j,k) = -sin(PI*y)**2 * sin(2*PI*x)
                if (time .gt. 1.0) then
                    vel_v(i,j,k) = sin(PI*y)**2 * sin(2*PI*x) * sin(2*PI*z)
                endif
			endif
			
			
		  end do
	  end do
  end do
  
  ! w-vel
  do k = lo(3), hi(3)+1
	  z = prob_lo(3) + (dble(k)+0.0d0) * dx(3)
	  do j = lo(2), hi(2)
		  y = prob_lo(2) + (dble(j)+0.0d0) * dx(2)
		  do i = lo(1), hi(1)
			x = prob_lo(1) + (dble(i)+0.0d0) * dx(1)
			
			if (AMREX_SPACEDIM.eq.2) then
				!no w-vel
			else 
				vel_w(i,j,k) = -sin(PI*z)**2 * sin(2*PI*x) * sin(2*PI*y)
                !vel_w(i,j,k) = 0
                if (time .gt. 1.0) then
                    vel_w(i,j,k) = sin(PI*z)**2 * sin(2*PI*x) * sin(2*PI*y)
                endif
			endif
			
			
		  end do
	  end do
  end do
  

end subroutine init_vel