!cell sorting (each cell contains particle indices of which particles lie within it), taken from AMReX "CellSortedParticles tutorial"
module cell_sorted_particle_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int
  
  implicit none
  
  private
  
  public remove_particle_from_cell
  
  contains
  
    subroutine remove_particle_from_cell(cell_parts, cell_np, new_np, i)
    
      use iso_c_binding, only: c_int
    
      implicit none
    
      integer(c_int), intent(in   ) :: cell_np
      integer(c_int), intent(inout) :: cell_parts(cell_np)
      integer(c_int), intent(inout) :: new_np
      integer(c_int), intent(in   ) :: i 

      cell_parts(i) = cell_parts(new_np)
      new_np = new_np - 1
        
    end subroutine remove_particle_from_cell
  
end module cell_sorted_particle_module


!subroutine move_particles(particles, np, lo, hi, &
!     cell_part_ids, cell_part_cnt, clo, chi, plo, dx, dt) &
!     bind(c,name="move_particles")
!  
!  use amrex_fort_module, only: amrex_real
!  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
!  use cell_sorted_particle_module, only: particle_t, remove_particle_from_cell
!  
!  implicit none
!
!  type(particle_t), intent(inout), target :: particles(np)
!  integer(c_int), intent(in) :: np
!  integer(c_int), intent(in) :: lo(3), hi(3)
!  integer(c_int), intent(in) :: clo(3), chi(3)
!  type(c_ptr), intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
!  integer(c_int), intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
!  real(amrex_real), intent(in) :: plo(3)
!  real(amrex_real), intent(in) :: dx(3)
!  real(amrex_real), intent(in) :: dt
!  
!  integer :: i, j, k, p, cell_np, new_np
!  integer :: cell(3)
!  integer(c_int), pointer :: cell_parts(:)
!  type(particle_t), pointer :: part
!  real(amrex_real) inv_dx(3)
!
!  inv_dx = 1.d0/dx
!  
!  do k = lo(3), hi(3)
!     do j = lo(2), hi(2)
!        do i = lo(1), hi(1)
!           cell_np = cell_part_cnt(i,j,k)
!           call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])
!
!           new_np = cell_np
!           p = 1
!           do while (p <= new_np)
!              part => particles(cell_parts(p))
!              
!              ! move the particle in a straight line
!              part%pos = part%pos + dt*part%vel
!
!              ! if it has changed cells, remove from vector.
!              ! otherwise continue
!              cell = floor((part%pos - plo)*inv_dx)              
!              if ((cell(1) /= i) .or. (cell(2) /= j) .or. (cell(3) /= k)) then
!                 part%sorted = 0
!                 call remove_particle_from_cell(cell_parts, cell_np, new_np, p)  
!              else
!                 p = p + 1
!              end if
!           end do
!
!           cell_part_cnt(i,j,k) = new_np
!           
!        end do
!     end do
!  end do
!  
!end subroutine move_particles



subroutine minimizeError(particles, np, &
                         lo, hi, &
                         phi, philo, phihi, &
                         cell_part_ids, cell_part_cnt, clo, chi, &
                         mi_i, mi_j, mi_k, maskSize, & 
                         p_lo, dx, &
                         polyOrder, numrows, numcols, w, P, A_inv) bind(C, name="minimizeError")

 use amrex_fort_module, only : amrex_real
 use iso_c_binding, only: c_ptr, c_int, c_f_pointer
 use IntfPtc_module, only : particle_t
 implicit none
 !interface for to get interpolating poly coeffs
 
 INTERFACE
  SUBROUTINE getPolyInterpCoeffs(polyOrder, numrows, numcols, dx, w, P, A_inv, &
                               phi, philo, phihi, indexVec, a_coeffs) BIND(C, name="C_getPolyInterpCoeffs")
   use amrex_fort_module, only : amrex_real
   IMPLICIT NONE
   integer, intent(in) :: polyOrder, numrows, numcols
   real(amrex_real), dimension(3), intent(in) :: dx
   integer, intent(in) ::  philo(3), phihi(3)
   real(amrex_real), intent(in) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
   real(amrex_real), dimension(numrows,numcols), intent(in) :: P !polynomial matrix
   real(amrex_real), dimension(numrows), intent(in) :: w !stencil weights
   real(amrex_real), dimension(numcols,numcols), intent(in) :: A_inv 
   real(amrex_real), dimension(numcols), intent(out) :: a_coeffs !coefficients of poly. approx. to phi
   integer, dimension(3), intent(in) :: indexVec
  END SUBROUTINE getPolyInterpCoeffs
 END INTERFACE
 
 type(particle_t), intent(in), target :: particles(np)
 integer, intent(in) :: np
 integer(c_int), intent(in) :: lo(3), hi(3)
 integer(c_int), intent(in) :: clo(3), chi(3)
 type(c_ptr), intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
 integer(c_int), intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
       
 integer, intent(in) :: philo(3), phihi(3)
 real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 integer, intent(in) :: maskSize
 integer, dimension(maskSize), intent(in) :: mi_i, mi_j, mi_k
 !real(amrex_real), dimension(maskSize), intent(in) :: a_coeffs_mask0,a_coeffs_mask1,a_coeffs_mask2,a_coeffs_mask3
 real(amrex_real), intent(in) :: dx(3), p_lo(3)
 
 integer, intent(in) :: polyOrder, numrows, numcols
 real(amrex_real), dimension(numrows,numcols), intent(in) :: P !polynomial matrix
 real(amrex_real), dimension(numrows), intent(in) :: w !stencil weights
 real(amrex_real), dimension(numcols,numcols), intent(in) :: A_inv 
 real(amrex_real), dimension(numcols) :: a_coeffs !coefficients of poly. approx. to phi
 
 integer :: i, j, k, cell_np, new_np, m, countptc
 integer :: r,t,s
 integer :: index_a
 integer :: cell(3)
 integer(c_int), pointer :: cell_parts(:)
 type(particle_t), pointer :: part
 real(amrex_real) inv_dx(3)
 real(amrex_real) :: lambda, weight_p, sum_weight_p, error_cell
 
 integer :: p_iter
 real(amrex_real) :: cellCenter_x, cellCenter_y, cellCenter_z, phi_tilde
 real(amrex_real) :: PtcPos(3)
 integer, dimension(3) :: currentIndex
 
  inv_dx = 1.d0/dx
  
  !countptc = 0;
  !!can improve by looping only over cells containing particles (use mask)
  !do k = lo(3), hi(3)
  !do j = lo(2), hi(2)
  !do i = lo(1), hi(1)
  !         cell_np = cell_part_cnt(i,j,k)
  !         call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])
  !         
  !         !if (cell_np .GT. 0) then
  !         ! print '(A23,I2,A2,I2,A2,I2,A4,3I3)', "num particles in cell (", i,", ", j,", ",k, ") : ", cell_np
  !         !endif
  !         countptc = countptc+cell_np
  !         
  !end do !end do i
  !end do !end do j
  !end do !end do k
  !print *, "sum of particles in cells : ", countptc
  
  countptc = 0
  do m = 1, maskSize !iterating over outer narrowband (can make inner, but easier do w/ outer since init. in same func as interface ptc's)
   i = mi_i(m)
   j = mi_j(m)
   k = mi_k(m)
   
   !!put in a check so that stencil never goes outside of narrowband (interface particles should be in center of narrowband)
   
   cell_np = cell_part_cnt(i,j,k)
   call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])
   !print '(A23,I2,A2,I2,A2,I2,A4,3I3)', "num particles in cell (", i,", ", j,", ",k, ") : ", cell_np
   countptc = countptc+cell_np
   
   cellCenter_x = p_lo(1) + (i+0.5)*dx(1)
   cellCenter_y = p_lo(2) + (j+0.5)*dx(2)
   cellCenter_z = p_lo(3) + (k+0.5)*dx(3)
   currentIndex = (/i, j, k /)
   call getPolyInterpCoeffs(polyOrder, numrows, numcols, dx, w, P, A_inv, phi,philo,phihi, currentIndex, a_coeffs)
    
   p_iter = 1
   error_cell = 0.0
   sum_weight_p = 0.0
   lambda = 0.0
   do while (p_iter <= cell_np) !iterate through particles in cell
    part => particles(cell_parts(p_iter))

    PtcPos(1) = part%pos(1)
    PtcPos(2) = part%pos(2)
    PtcPos(3) = part%pos(3)
   
    index_a = 1
    phi_tilde = 0.0
    do r = 0, polyOrder
    do t = 0, polyOrder
    do s = 0, polyOrder
     if (s+t+r .GT. polyOrder) then
      CYCLE
     else
      !print *, 'str ', s+t+r
      phi_tilde = phi_tilde + a_coeffs(index_a)*(PtcPos(1)-cellCenter_x)**s *(PtcPos(2)-cellCenter_y)**t *(PtcPos(3)-cellCenter_z)**r  
      index_a = index_a+1
     end if
    end do !end do s
    end do !end do t
    end do !end do r
    
    !print *, "phidiff : ", phi_tilde-phi(i,j,k)
    
    !find \lambda which minimizes the error
    !weight_p = sqrt((PtcPos(1)-cellCenter_x)**2 +(PtcPos(2)-cellCenter_y)**2 +(PtcPos(3)-cellCenter_z)**2)
    weight_p = 1.0/((PtcPos(1)-cellCenter_x)**2 +(PtcPos(2)-cellCenter_y)**2 +(PtcPos(3)-cellCenter_z)**2 +10**-8)
    !print*, weight_p
    !print '(A23,I2,A2,I2,A2,I2,A4,F2.6)', "weights in cell (", i,", ", j,", ",k, ") : ", weight_p
    sum_weight_p = sum_weight_p + weight_p
    error_cell = error_cell + weight_p*(phi_tilde-0)
    
    
    p_iter = p_iter + 1
   end do !end do while
   
   
   if (cell_np .GT. 0) then
    lambda = -error_cell/sum_weight_p
    !print *, "lambda : ", lambda
    phi(i,j,k) = phi(i,j,k)+lambda
   endif
   
   !print *, "m : ", m, "   mSize : ", maskSize
  end do ! end do m mask
  
  if (countptc .NE. np) then
   print *, "warning: lost particles"
   print *, "total number of particles (in block) : ", np
   print *, "sum of particles in cells (in block) : ", countptc
   ! stop
  endif
  
end subroutine minimizeError
