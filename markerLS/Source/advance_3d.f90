!module mask_particles
! use amrex_fort_module, only: amrex_real, amrex_particle_real
! use iso_c_binding ,    only: c_int
! implicit none
! 
! public particle_t
! type, bind(C)  :: particle_t
!  real(amrex_particle_real) :: pos(3)
!  integer(c_int)   :: id
!  integer(c_int)   :: cpu
!  integer(c_int)   :: i
!  integer(c_int)   :: j
!  integer(c_int)   :: k
! end type particle_t
!
!end module mask_particles

function Delta_op(a, b)
 use amrex_fort_module, only : amrex_real
 implicit none

 real(amrex_real), intent(in) :: a, b 
 real(amrex_real) :: Delta_op !output

 !Delta operator
 !Delta_plus(Phi_k) = Phi_kp1 - Phi_k
 !Delta_minus(Phi_k) = Phi_k - Phi_km1
 Delta_op = a-b
 RETURN
end function Delta_op

function Delta_minusplus(a, b, c)
 use amrex_fort_module, only : amrex_real
 implicit none

 real(amrex_real), intent(in) :: a, b, c
 real(amrex_real) :: Delta_minusplus !output

 !Delta_minus(Delta_plus(b))
 Delta_minusplus = a-2*b+c
 RETURN
end function Delta_minusplus

function Phi_WENO(a, b, c, d) result(phiWENO)
 use amrex_fort_module, only : amrex_real
 use constants
 implicit none
 
 real(amrex_real), intent(in) :: a, b, c, d
 
 real(amrex_real) phiWENO !output
 real(amrex_real) IS_0, IS_1, IS_2
 real(amrex_real) alpha_0, alpha_1, alpha_2
 real(amrex_real) omega_0, omega_2
 
 IS_0 = 13*(a-b)**2 + 3*(a-3*b)**2
 IS_1 = 13*(b-c)**2 + 3*(b+c)**2
 IS_2 = 13*(c-d)**2 + 3*(3*c-d)**2
	
 alpha_0 = 1/((M_EPSILON+IS_0)**2)
 alpha_1 = 6/((M_EPSILON+IS_1)**2)
 alpha_2 = 3/((M_EPSILON+IS_2)**2)

 omega_0 = alpha_0/(alpha_0+alpha_1+alpha_2)
 omega_2 = alpha_2/(alpha_0+alpha_1+alpha_2)
	
 phiWENO = (1.0/3.0)*omega_0*(a-2*b+c)+(1.0/6.0)*(omega_2-0.5)*(b-2*c+d)

end function Phi_WENO

subroutine dPhi_dx(dx,Phi_im3,Phi_im2,Phi_im1,Phi_i,Phi_ip1,Phi_ip2,Phi_ip3, &
			       dPhi_dx_L,dPhi_dx_R) !output (u_minus, u_plus, ...)
 !returns left and right biased derivatives of Phi
 use amrex_fort_module, only : amrex_real
 implicit none

 real(amrex_real), intent(in) :: dx
 real(amrex_real), intent(in) :: Phi_im3,Phi_im2,Phi_im1,Phi_i,Phi_ip1,Phi_ip2,Phi_ip3
 real(amrex_real), intent(out) :: dPhi_dx_L,dPhi_dx_R
 !dPhi_dx_L: stencil spans Phi_im3 to Phi_ip2
 !dPhi_dx_R: stencil spans Phi_im2 to Phi_ip3
 
 ! local variables
 real(amrex_real) :: Dp_Phi_im2, Dp_Phi_im1, Dp_Phi_i, Dp_Phi_ip1
 real(amrex_real) :: a_L, b_L, c_L, d_L, a_R, b_R, c_R, d_R
 real(amrex_real) :: Delta_op, Delta_minusplus, Phi_WENO !functions

 !Dp : Delta_plus, second var subtracted from first
 Dp_Phi_im2 = Delta_op(Phi_im1, Phi_im2)
 Dp_Phi_im1 = Delta_op(Phi_i, Phi_im1)
 Dp_Phi_i   = Delta_op(Phi_ip1, Phi_i)
 Dp_Phi_ip1 = Delta_op(Phi_ip2, Phi_ip1)
 
 !Delta_minusplus operator is applied to second var passed
 a_L = Delta_minusplus(Phi_im1,Phi_im2,Phi_im3)/dx
 b_L = Delta_minusplus(Phi_i,Phi_im1,Phi_im2)/dx
 c_L = Delta_minusplus(Phi_ip1,Phi_i,Phi_im1)/dx
 d_L = Delta_minusplus(Phi_ip2,Phi_ip1,Phi_i)/dx
	
 a_R = Delta_minusplus(Phi_ip3,Phi_ip2,Phi_ip1)/dx
 b_R = Delta_minusplus(Phi_ip2,Phi_ip1,Phi_i)/dx
 c_R = Delta_minusplus(Phi_ip1,Phi_i,Phi_im1)/dx
 d_R = Delta_minusplus(Phi_i,Phi_im1,Phi_im2)/dx

 dPhi_dx_L = (-Dp_Phi_im2+7*Dp_Phi_im1+7*Dp_Phi_i-Dp_Phi_ip1)/(12*dx) - Phi_WENO(a_L,b_L,c_L,d_L)
 dPhi_dx_R = (-Dp_Phi_im2+7*Dp_Phi_im1+7*Dp_Phi_i-Dp_Phi_ip1)/(12*dx) + Phi_WENO(a_R,b_R,c_R,d_R)

end subroutine dPhi_dx

!!fix H_func, alpha_op - pass vars instead of hardcode
function H_adv(xspeed, yspeed, zspeed, u, v, w)
 use amrex_fort_module, only : amrex_real
 implicit none
 
 real(amrex_real) :: H_adv !output
 real(amrex_real), intent(in) :: xspeed, yspeed, zspeed, u, v, w
 !phi_t + H(x,t,phi,phi_x) = 0
 !H = vel dot grad(phi) // H = a*phi_x 
 H_adv = xspeed*u + yspeed*v + zspeed*w
 RETURN
end function H_adv

function alpha_op(speed, u_plus, u_minus)
 use amrex_fort_module, only : amrex_real
 implicit none
 
 real(amrex_real) :: alpha_op !output
 real(amrex_real), intent(in) :: speed, u_plus, u_minus
	
 real(amrex_real) :: a
 
 a = speed
	
 alpha_op = abs(a)
 RETURN
end function alpha_op

	
function LLF_flux(xspeed, yspeed, zspeed, u_minus, u_plus, v_minus, v_plus, w_minus, w_plus)
 use amrex_fort_module, only : amrex_real
 implicit none

 real(amrex_real) :: LLF_flux !output
 real(amrex_real) :: H_adv, alpha_op !functions
 
 real(amrex_real), intent(in) :: xspeed, yspeed, zspeed
 real(amrex_real), intent(in) :: u_minus, u_plus, v_minus, v_plus, w_minus, w_plus
 ! u_minus, u_plus: right going, left going characteristics
	
 !dPhi_dx(dx, Phi_im3, Phi_im2, Phi_im1, Phi_i, Phi_ip1, Phi_ip2, Phi_ip3, u_minus, u_plus) !use weno scheme to get u_minus, u_plus
 LLF_flux = H_adv(xspeed,yspeed,zspeed,(u_plus+u_minus)/2,(v_plus+v_minus)/2,(w_plus+w_minus)/2) &
            - alpha_op(xspeed,u_plus,u_minus)*(u_plus-u_minus)/2 &
            - alpha_op(yspeed,v_plus,v_minus)*(v_plus-v_minus)/2 &
            - alpha_op(zspeed,w_plus,w_minus)*(w_plus-w_minus)/2
	
 RETURN
end function LLF_flux


subroutine TimeDeriv_WENO(dphidt,phi,philo,phihi, &
                     xVel,ulo,uhi,yVel,vlo,vhi,zVel,wlo,whi, &
                     dx, mask,masklo,maskhi,mi_i,mi_j,mi_k,mask_size)

 use amrex_fort_module, only : amrex_real
 implicit none

 real(amrex_real), dimension(mask_size), intent(out) :: dphidt
 integer, intent(in) ::  philo(3), phihi(3), masklo(3), maskhi(3), ulo(3),uhi(3), vlo(3),vhi(3), wlo(3),whi(3)
 real(amrex_real), intent(in) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 real(amrex_real), intent(in) :: xVel( ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3))
 real(amrex_real), intent(in) :: yVel( vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3))
 real(amrex_real), intent(in) :: zVel( wlo(1):whi(1), wlo(2):whi(2), wlo(3):whi(3))
 integer, intent(in) :: mask(masklo(1):maskhi(1),masklo(2):maskhi(2),masklo(3):maskhi(3))
 real(amrex_real), intent(in) :: dx(3)
 integer, intent(in) :: mask_size
 integer, dimension(mask_size), intent(in) :: mi_i, mi_j, mi_k
 
 !local variables
 real(amrex_real) :: LLF_flux !function
 !real(amrex_real) :: Phi_im3,Phi_im2,Phi_im1,Phi_i,Phi_ip1,Phi_ip2,Phi_ip3
 !real(amrex_real) :: Phi_jm3,Phi_jm2,Phi_jm1,Phi_j,Phi_jp1,Phi_jp2,Phi_jp3
 !real(amrex_real) :: Phi_km3,Phi_km2,Phi_km1,Phi_k,Phi_kp1,Phi_kp2,Phi_kp3
 !real(amrex_real) :: u_minus, u_plus, v_minus, v_plus, w_minus, w_plus
 real(amrex_real) :: xspeed, yspeed, zspeed
 real(amrex_real) :: u_minus, u_plus, v_minus, v_plus, w_minus, w_plus
 integer :: i,j,k, m
 
 do m = 1, mask_size
  i = mi_i(m)
  j = mi_j(m)
  k = mi_k(m)
  !dphidx solve --> return u_minus,u_plus
  !dPhi_dx(dx(0), Phi_im3, Phi_im2, Phi_im1, Phi_i, Phi_ip1, Phi_ip2, Phi_ip3, u_minus, u_plus)
  call dPhi_dx(dx(1), phi(i-3,j,k), phi(i-2,j,k), phi(i-1,j,k), phi(i,j,k), phi(i+1,j,k), phi(i+2,j,k), phi(i+3,j,k), u_minus, u_plus)
  !dphidy solve --> return v_minus,v_plus
  !dPhi_dx(dx(1), Phi_jm3, Phi_jm2, Phi_jm1, Phi_j, Phi_jp1, Phi_jp2, Phi_jp3, v_minus, v_plus)
  call dPhi_dx(dx(2), phi(i,j-3,k), phi(i,j-2,k), phi(i,j-1,k), phi(i,j,k), phi(i,j+1,k), phi(i,j+2,k), phi(i,j+3,k), v_minus, v_plus)
  !dphidz solve --> return w_minus,w_plus
  !dPhi_dx(dx(1), Phi_km3, Phi_km2, Phi_km1, Phi_k, Phi_kp1, Phi_kp2, Phi_kp3, w_minus, w_plus)
  call dPhi_dx(dx(3), phi(i,j,k-3), phi(i,j,k-2), phi(i,j,k-1), phi(i,j,k), phi(i,j,k+1), phi(i,j,k+2), phi(i,j,k+3), w_minus, w_plus)
  
  !!!interp, change this
  xspeed = (xVel(i+1,j,k)+xVel(i,j,k))/2.0
  yspeed = (yVel(i,j+1,k)+yVel(i,j,k))/2.0
  zspeed = (zVel(i,j,k+1)+zVel(i,j,k))/2.0
  
  dphidt(m) = -LLF_flux(xspeed,yspeed,zspeed, u_minus,u_plus,v_minus,v_plus,w_minus,w_plus)
 end do !m=1,NBmask_size
end subroutine TimeDeriv_WENO

subroutine RK3_Williamson (lo, hi, phi_old, phi_oldlo, phi_oldhi, phi, philo, phihi, &
                       mask, masklo, maskhi, &
                       vel_u, ulo, uhi, vel_v, vlo, vhi, vel_w, wlo, whi, &
                       dx, dt, time, &
                       mi_i, mi_j, mi_k, &
                       mask_size, iter, G, Glo, Ghi) bind(C, name="RK3_Williamson")

 use amrex_fort_module, only : amrex_real
 implicit none

 integer, intent(in) :: lo(3), hi(3), philo(3), phihi(3), masklo(3), maskhi(3)
 integer, intent(in) :: phi_oldlo(3), phi_oldhi(3)
 integer, intent(in) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)
 real(amrex_real), intent(inout) :: phi_old(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 integer, intent(in)    :: mask(masklo(1):maskhi(1),masklo(2):maskhi(2),masklo(3):maskhi(3))
 real(amrex_real), intent(in)    :: vel_u( ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3))
 real(amrex_real), intent(in)    :: vel_v( vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3))
 real(amrex_real), intent(in)    :: vel_w( wlo(1):whi(1), wlo(2):whi(2), wlo(3):whi(3))
 real(amrex_real), intent(in)    :: dx(3)
 real(amrex_real), intent(in)    :: dt, time
 integer, intent(in) :: mask_size
 integer, dimension(mask_size), intent(in) :: mi_i, mi_j, mi_k
 integer, intent(in) :: iter
 
 
 !local variables
  
 !call RK3_Williamson(phinew,pnlo,pnhi,mask,masklo,maskhi,dx,dt,time,dphidt,mi_i,mi_j,mi_k)

 real(amrex_real),parameter,dimension(3) :: a_coeffs = (/ 0.0, -5.0/9.0, -153.0/128.0 /)
 real(amrex_real),parameter,dimension(3) :: b_coeffs = (/ 0.0, 1.0/3.0, 3.0/4.0 /)
 real(amrex_real),parameter,dimension(3) :: g_coeffs = (/ 1.0/3.0, 15.0/16.0, 8.0/15.0 /)
 
 integer, intent(in) :: Glo(3), Ghi(3)
 real(amrex_real), intent(inout) :: G(Glo(1):Ghi(1),Glo(2):Ghi(2),Glo(3):Ghi(3))
 !real(amrex_real), dimension(mask_size) :: G
 real(amrex_real), dimension(mask_size) :: dphidt!(SIZE(mi_i),SIZE(mi_j),SIZE(mi_k))
 real(amrex_real) t 
 !integer iter, m !iterators
 integer :: m
 integer i,j,k
 
 !do iter = 1, 3 !rk intervals
 ! t = time + b_coeffs(iter)*dt
 ! dphi_dt = TimeDeriv()!!TODO
 ! do k = lo(3), hi(3)
 ! do j = lo(2), hi(2)
 ! do i = lo(1), hi(1)
 !  G(i,j,k) = a_coeffs(iter)*G(i,j,k) + dphi_dt(i,j,k)
 !  phi(i,j,k) = phi(i,j,k) + g_coeffs(iter)*dt*G(i,j,k)
 ! end do
 ! end do
 ! end do
 !end do !end iter
 
 !if (iter .EQ. 1) then
 !   G = 0 !init array G to zero on first RK step
 !end if
 
 !print *, 'Iter ', iter
 !G(m) --> G(i,j,k) temp??
  
 ! using mask
 !do iter = 1, 3 !rk intervals
  t = time + b_coeffs(iter)*dt
  call TimeDeriv_WENO(dphidt,phi,philo,phihi,vel_u,ulo,uhi,vel_v,vlo,vhi,vel_w,wlo,whi,&
                 dx,mask,masklo,maskhi,mi_i,mi_j,mi_k,mask_size) !return dphidt
  do m = 1, mask_size
   i = mi_i(m)
   j = mi_j(m)
   k = mi_k(m)
   !if (mask(mi_i(m),mi_j(m),mi_k(m)) .EQ. 2) then
    !G(m) = a_coeffs(iter)*G(m) + dphidt(m)
    !phi(mi_i(m),mi_j(m),mi_k(m)) = phi(mi_i(m),mi_j(m),mi_k(m)) + g_coeffs(iter)*dt*G(m)
    
    G(i,j,k) = a_coeffs(iter)*G(i,j,k) + dphidt(m)
    phi(i,j,k) = phi(i,j,k) + g_coeffs(iter)*dt*G(i,j,k)
   !end if
  end do
 !end do

end subroutine RK3_Williamson



subroutine attractToIntf(particles_data, num_Ptc, &
                    cell_part_ids, cell_part_cnt, clo, chi, &
                    phi, philo, phihi, &
                    n_grid, nlo, nhi, &
                    vel_u, ulo, uhi, vel_v, vlo, vhi, vel_w, wlo, whi, &
                    maskSize, &
                    mi_i, mi_j, mi_k, dx, p_lo, &
                    polyOrder, numrows, numcols, w, P, A_inv) bind(C, name="attractToIntf")
                    
 use amrex_fort_module, only : amrex_real
 use iso_c_binding, only: c_ptr, c_int, c_f_pointer
 use IntfPtc_module, only : particle_t
 implicit none
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
   real(amrex_real), dimension(numcols), intent(in) :: a_coeffs !coefficients of poly. approx. to phi
   integer, dimension(3), intent(in) :: indexVec
  END SUBROUTINE getPolyInterpCoeffs
 END INTERFACE
 
 integer, intent(in) :: num_Ptc, maskSize
 type(particle_t), intent(inout), target :: particles_data(num_Ptc)      
 integer, intent(in) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)
 integer, intent(in) :: philo(3), phihi(3)
 real(amrex_real), intent(in) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 integer, intent(in) :: nlo(3), nhi(3)
 real(amrex_real), intent(in) :: n_grid(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),1:3)
 real(amrex_real), intent(in) :: vel_u( ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3))
 real(amrex_real), intent(in) :: vel_v( vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3))
 real(amrex_real), intent(in) :: vel_w( wlo(1):whi(1), wlo(2):whi(2), wlo(3):whi(3))
 integer, dimension(maskSize), intent(in) :: mi_i, mi_j, mi_k
 real(amrex_real), intent(in) :: dx(3), p_lo(3)
 integer, intent(in) :: polyOrder, numrows, numcols
 real(amrex_real), dimension(numrows,numcols), intent(in) :: P !polynomial matrix
 real(amrex_real), dimension(numrows), intent(in) :: w !stencil weights
 real(amrex_real), dimension(numcols,numcols), intent(in) :: A_inv 
 
 real(amrex_real), dimension(numcols):: a_coeffs_x, a_coeffs_y, a_coeffs_z
 integer :: i,j,k, m,step
 real(amrex_real) :: phi_x, phi_y, phi_z
 real(amrex_real) :: nx, ny, nz
 real(amrex_real) :: lambda
 !real(amrex_real) :: PtcPos_old(3)
 integer, dimension(3) :: currentIndex
 integer :: s,t,r
 real(amrex_real) :: cellCenter_x, cellCenter_y, cellCenter_z
 integer :: index_a
 real(amrex_real) :: n_gridx(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3))
 real(amrex_real) :: n_gridy(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3))
 real(amrex_real) :: n_gridz(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3))
 
 integer(c_int), intent(in) :: clo(3), chi(3)
 type(c_ptr), intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
 integer(c_int), intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
 integer :: cell_np
 integer(c_int), pointer :: cell_parts(:)
 type(particle_t), pointer :: part
 integer :: p_iter
 real(amrex_real) :: PtcPos(3)
 
 do m = 1, maskSize !num_Ptc !loop over cells, then go over particles within the cell
  !cell index containing particle
  !i = NINT((particles_data(m)%pos(1)-p_lo(1))/dx(1)-0.5)
  !j = NINT((particles_data(m)%pos(2)-p_lo(2))/dx(2)-0.5)
  !k = NINT((particles_data(m)%pos(3)-p_lo(3))/dx(3)-0.5)
  i = mi_i(m)
  j = mi_j(m)
  k = mi_k(m)
  lambda = 1.0
  
  cell_np = cell_part_cnt(i,j,k)
  if (cell_np .GT. 0) then
   call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])
  
   !do step = 1, 1 !?don't need since know exact dist from interface
    !!find normal dir (LS)
    !phi_x = (phi(i+1,j,k)-phi(i-1,j,k))/(2*dx(1))
    !phi_y = (phi(i,j+1,k)-phi(i,j-1,k))/(2*dx(2))
    !phi_z = (phi(i,j,k+1)-phi(i,j,k-1))/(2*dx(3))
    ! 
    !n_x = phi_x/abs(phi_x)
    !n_y = phi_y/abs(phi_y)
    !n_z = phi_z/abs(phi_z)
   cellCenter_x = p_lo(1) + (i+0.5)*dx(1)
   cellCenter_y = p_lo(2) + (j+0.5)*dx(2)
   cellCenter_z = p_lo(3) + (k+0.5)*dx(3)
   currentIndex = (/i, j, k /)
   
   n_gridx(:,:,:) = n_grid(:,:,:,1)
   n_gridy(:,:,:) = n_grid(:,:,:,2)
   n_gridz(:,:,:) = n_grid(:,:,:,3)
   
   !print*, "ngridx: ", n_gridx(i,j,k),"ngrid1: ", n_grid(i,j,k,1)
   
   !!get polynomial interpolation of grid normal
   !!add start component, numcomponent
   !!!only need to do this once per cell, not per particle
   call getPolyInterpCoeffs(polyOrder, numrows, numcols, dx, w, P, A_inv, n_gridx,nlo,nhi, currentIndex, a_coeffs_x)
   call getPolyInterpCoeffs(polyOrder, numrows, numcols, dx, w, P, A_inv, n_gridy,nlo,nhi, currentIndex, a_coeffs_y)
   call getPolyInterpCoeffs(polyOrder, numrows, numcols, dx, w, P, A_inv, n_gridz,nlo,nhi, currentIndex, a_coeffs_z)
   !print *, "a_coeffs_x: ", a_coeffs_x
   !print *, "a_coeffs_y: ", a_coeffs_y
   
   !PtcPos_old(1) = particles_data(m)%pos(1)
   !PtcPos_old(2) = particles_data(m)%pos(2)
   !PtcPos_old(3) = particles_data(m)%pos(3)
   p_iter = 1
   do while (p_iter <= cell_np) !iterate through particles in cell
    part => particles_data(cell_parts(p_iter))
    PtcPos(1) = part%pos(1)
    PtcPos(2) = part%pos(2)
    PtcPos(3) = part%pos(3)
   
   
    nx=0.0; ny=0.0; nz=0.0
    index_a = 1
    do r = 0, polyOrder
    do t = 0, polyOrder
    do s = 0, polyOrder
     if (s+t+r .GT. polyOrder) then
      CYCLE
     else
      !print *, 'str ', s+t+r
      nx = nx + a_coeffs_x(index_a)*(PtcPos(1)-cellCenter_x)**s *(PtcPos(2)-cellCenter_y)**t *(PtcPos(3)-cellCenter_z)**r 
      ny = ny + a_coeffs_y(index_a)*(PtcPos(1)-cellCenter_x)**s *(PtcPos(2)-cellCenter_y)**t *(PtcPos(3)-cellCenter_z)**r 
      nz = nz + a_coeffs_z(index_a)*(PtcPos(1)-cellCenter_x)**s *(PtcPos(2)-cellCenter_y)**t *(PtcPos(3)-cellCenter_z)**r 
      index_a = index_a+1
     end if
    end do !end do s
    end do !end do t
    end do !end do r
   
    !print *, 'nx ', nx, "n_grid(i,j,k,1)", n_grid(i,j,k,1)
   
    !x_new = x_p + lambda*(phi_goal-phi_p)*n(x_p)
    !particles_data(m)%pos(1) = PtcPos_old(1) + lambda*(0-particles_data(p_iter)%phi(1))*nx
    !particles_data(m)%pos(2) = PtcPos_old(2) + lambda*(0-particles_data(p_iter)%phi(1))*ny
    !particles_data(m)%pos(3) = PtcPos_old(3) + lambda*(0-particles_data(p_iter)%phi(1))*nz
    part%pos(1) = PtcPos(1) + lambda*(0-part%phi(1))*nx!n_grid(i,j,k,1)!nx
    part%pos(2) = PtcPos(2) + lambda*(0-part%phi(1))*ny!n_grid(i,j,k,2)!ny
    part%pos(3) = PtcPos(3) + lambda*(0-part%phi(1))*nz!n_grid(i,j,k,3)!nz
   
    !!TODO: change phi_cell to phi_p -- need to interp to phi to particle (DONE) -- replace least squares w/ tricubic interp
   
    !lambda = lambda/2.0
   
    p_iter = p_iter + 1
   end do !end do while p_iter
  end if !end if cell_np>0
 end do !end do m mask
 
 
end subroutine attractToIntf


subroutine RK3_Williamson_Ptc(particles_data, &
                       vel_u, ulo, uhi, vel_v, vlo, vhi, vel_w, wlo, whi, &
                       dt, time, &
                       num_Ptc, &
                       pi_i, pi_j, pi_k, dx, p_lo, &
                       iter) bind(C, name="RK3_Williamson_Ptc") !, G_Ptc, Glo, Ghi
                      
 use amrex_fort_module, only : amrex_real
 use IntfPtc_module, only : particle_t
 implicit none
 
 integer, intent(in) :: num_Ptc
 type(particle_t), intent(inout) :: particles_data(num_Ptc)      
 integer, intent(in) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)
 real(amrex_real), intent(in)    :: vel_u( ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3))
 real(amrex_real), intent(in)    :: vel_v( vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3))
 real(amrex_real), intent(in)    :: vel_w( wlo(1):whi(1), wlo(2):whi(2), wlo(3):whi(3)) 
 real(amrex_real), intent(in)    :: dt, time 
 integer, intent(in) :: iter
 !integer, intent(in) :: Glo(3), Ghi(3)
 !real(amrex_real), intent(inout) :: G_Ptc(Glo(1):Ghi(1),Glo(2):Ghi(2),Glo(3):Ghi(3),1:3)
 integer, dimension(num_Ptc), intent(in) :: pi_i, pi_j, pi_k
 real(amrex_real), intent(in) :: dx(3), p_lo(3)
 
 !local variables
 real(amrex_real),parameter,dimension(3) :: a_coeffs = (/ 0.0, -5.0/9.0, -153.0/128.0 /)
 real(amrex_real),parameter,dimension(3) :: b_coeffs = (/ 0.0, 1.0/3.0, 3.0/4.0 /)
 real(amrex_real),parameter,dimension(3) :: g_coeffs = (/ 1.0/3.0, 15.0/16.0, 8.0/15.0 /)
 real(amrex_real) t 
 integer m, dir
 integer i,j,k
 real(amrex_real) vel
 integer PtcCell_i, PtcCell_j, PtcCell_k
 real(amrex_real) velPosL_i, velPosL_j, velPosL_k
 real(amrex_real) xd, yd, zd
 real(amrex_real) ptc_u, ptc_v, ptc_w
 
 
 t = time + b_coeffs(iter)*dt
 
 do m = 1, num_Ptc
 
  if (iter .EQ. 1) then
   do dir = 1,3
    particles_data(m)%G_Ptc(dir) = 0.0
   end do
  end if
  !!position of left index of cell vel.
  !velPosL_i = p_lo(1) + (pi_i(m)+0.0)*dx(1)
  !velPosL_j = p_lo(2) + (pi_j(m)+0.0)*dx(2)
  !velPosL_k = p_lo(3) + (pi_k(m)+0.0)*dx(3)
  
  !cell index containing particle
  PtcCell_i = NINT((particles_data(m)%pos(1)-p_lo(1))/dx(1)-0.5)
  PtcCell_j = NINT((particles_data(m)%pos(2)-p_lo(2))/dx(2)-0.5)
  PtcCell_k = NINT((particles_data(m)%pos(3)-p_lo(3))/dx(3)-0.5)
  
  !position of left index of cell vel.
  velPosL_i = p_lo(1) + (PtcCell_i+0.0)*dx(1)
  velPosL_j = p_lo(2) + (PtcCell_j+0.0)*dx(2)
  velPosL_k = p_lo(3) + (PtcCell_k+0.0)*dx(3)
  
  !!i = pi_i(m)
  !!j = pi_j(m)
  !!k = pi_k(m)
  i = PtcCell_i
  j = PtcCell_j
  k = PtcCell_k
   
  !linear interp. velocity
  xd = (particles_data(m)%pos(1)-velPosL_i)/dx(1)
  yd = (particles_data(m)%pos(2)-velPosL_j)/dx(2)
  zd = (particles_data(m)%pos(3)-velPosL_k)/dx(3)
  ptc_u = vel_u(i,j,k)*(1-xd) + vel_u(i+1,j,k)*xd
  ptc_v = vel_v(i,j,k)*(1-yd) + vel_v(i,j+1,k)*yd
  ptc_w = vel_w(i,j,k)*(1-zd) + vel_w(i,j,k+1)*zd
   
  do dir = 1, 3
   if (dir .EQ. 1) then
    !vel = vel_u(i,j,k)
    vel = ptc_u
   elseif (dir .EQ. 2) then
    !vel = vel_v(i,j,k)
    vel = ptc_v
   else
    !vel = vel_w(i,j,k)
    vel = ptc_w
   endif
  
   !G_Ptc(i,j,k,dir) = a_coeffs(iter)*G_Ptc(i,j,k,dir) + vel
   !particles_data(m)%pos(dir) = particles_data(m)%pos(dir) + g_coeffs(iter)*dt*G_Ptc(i,j,k,dir)
   
   particles_data(m)%G_Ptc(dir) = a_coeffs(iter)*particles_data(m)%G_Ptc(dir) + vel
   particles_data(m)%pos(dir) = particles_data(m)%pos(dir) + g_coeffs(iter)*dt*particles_data(m)%G_Ptc(dir)
  end do
 end do
 !!!TODO: Need to update particle cell index


end subroutine RK3_Williamson_Ptc




subroutine return_PlusMinus(var, varplus, varminus)
 use amrex_fort_module, only : amrex_real
 implicit none
 
 real(amrex_real), intent(in) :: var
 real(amrex_real), intent(out) :: varplus, varminus
    
 if (var>0) then 
  varplus = var
  varminus = 0
 else 
  varplus = 0
  varminus = var
 endif

end subroutine return_PlusMinus

function Redist(s, a, b, c, d, e, f)
 use amrex_fort_module, only : amrex_real
 implicit none
 
 real(amrex_real), intent(in) :: s
 real(amrex_real), intent(in) :: a, b, c, d, e, f
 real(amrex_real) :: Redist !output
 !real(amrex_real) H
 real(amrex_real) splus,sminus
 real(amrex_real) aplus, bplus, cplus, dplus, eplus, fplus, &
                  aminus, bminus, cminus, dminus, eminus, fminus
    
 call return_PlusMinus(a,aplus,aminus)
 call return_PlusMinus(b,bplus,bminus)
 call return_PlusMinus(c,cplus,cminus)
 call return_PlusMinus(d,dplus,dminus)
 call return_PlusMinus(e,eplus,eminus)
 call return_PlusMinus(f,fplus,fminus)
    
 call return_PlusMinus(s,splus,sminus)
	
 !Redist = splus*(sqrt(aplus**2+bminus**2+cplus**2+dminus**2+eplus**2+fminus**2)-1.0) &
 !         +sminus*(sqrt(aminus**2+bplus**2+cminus**2+dplus**2+eminus**2+fplus**2)-1.0)
 Redist = splus*(sqrt(MAX(aplus**2,bminus**2)+MAX(cplus**2,dminus**2)+MAX(eplus**2,fminus**2))-1.0) &
          +sminus*(sqrt(MAX(aminus**2,bplus**2)+MAX(cminus**2,dplus**2)+MAX(eminus**2,fplus**2))-1.0)
 RETURN
end function Redist

subroutine D_pm_phi(dx, d_im1, d_i, d_ip1, a, b)
 use amrex_fort_module, only : amrex_real
 implicit none

 real(amrex_real), intent(in) :: dx, d_im1, d_i, d_ip1
 real(amrex_real), intent(out) :: a, b
 
 a = (d_i - d_im1)/dx !D^(-)
 b = (d_ip1 - d_i)/dx !D^(+)
!!replace w/ high order
end subroutine D_pm_phi

subroutine TimeDeriv_Redist(dDdt, d_n, dlo, dhi, &
                     dx, mask, masklo, maskhi, mi_i, mi_j, mi_k, mask_size)

 use amrex_fort_module, only : amrex_real
 implicit none

 real(amrex_real), dimension(mask_size), intent(out) :: dDdt
 integer, intent(in) :: dlo(3), dhi(3), masklo(3), maskhi(3)
 real(amrex_real), intent(in) :: d_n(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
 integer, intent(in) :: mask(masklo(1):maskhi(1),masklo(2):maskhi(2),masklo(3):maskhi(3))
 real(amrex_real), intent(in) :: dx(3)
 integer, intent(in) :: mask_size
 integer, dimension(mask_size), intent(in) :: mi_i, mi_j, mi_k
 
 !local variables
 real(amrex_real) :: s, a, b, c, d, e, f
 real(amrex_real) :: d_ijk
 real(amrex_real) :: d_im1,d_ip1
 real(amrex_real) :: d_jm1,d_jp1
 real(amrex_real) :: d_km1,d_kp1
 real(amrex_real) :: Redist
 integer :: i,j,k, m
 
 
 
 
 do m = 1, mask_size
  i = mi_i(m)
  j = mi_j(m)
  k = mi_k(m)
  
  d_ijk = d_n(i,j,k)
  d_im1 = d_n(i-1,j,k)
  d_ip1 = d_n(i+1,j,k)
  d_jm1 = d_n(i,j-1,k)
  d_jp1 = d_n(i,j+1,k)
  d_km1 = d_n(i,j,k-1)
  d_kp1 = d_n(i,j,k+1)
  
  !!return D_x^(-) and D_x^(+)  --> a, b
  !call D_pm_phi(dx(1), d_n(i-1,j,k), d_n(i,j,k), d_n(i+1,j,k), a, b)
  !!return D_y^(-) and D_y^(+)  --> c, d
  !call D_pm_phi(dx(2), d_n(i,j-1,k), d_n(i,j,k), d_n(i,j+1,k), c, d)
  !!return D_z^(-) and D_z^(+)  --> e, f
  !call D_pm_phi(dx(3), d_n(i,j,k-1), d_n(i,j,k), d_n(i,j,k+1), e, f)
  

  
  !return D_x^(-) and D_x^(+)  --> a, b
  call D_pm_phi(dx(1), d_im1, d_ijk, d_ip1, a, b)
  !return D_y^(-) and D_y^(+)  --> c, d
  call D_pm_phi(dx(2), d_jm1, d_ijk, d_jp1, c, d)
  !return D_z^(-) and D_z^(+)  --> e, f
  call D_pm_phi(dx(3), d_km1, d_ijk, d_kp1, e, f)
  
    if (mask(i-1,j,k) .EQ. 0) then
      !a=b
      a=0
    endif
    if (mask(i+1,j,k) .EQ. 0) then
      !b=a
      b=0
    endif
    if (mask(i,j-1,k) .EQ. 0) then
      !c=d
      c=0
    endif
    if (mask(i,j+1,k) .EQ. 0) then
      !d=c
      d=0
    endif
    if (mask(i,j,k-1) .EQ. 0) then
      !e=f
      e=0
    endif
    if (mask(i,j,k+1) .EQ. 0) then
      !f=e
      f=0
    endif
  
  s = d_n(i,j,k)/sqrt(d_n(i,j,k)**2 + dx(1)**2) !!temporary, fix this
  
  dDdt(m) = Redist(s, a, b, c, d, e, f)
  
 end do !m=1,NBmask_size
end subroutine TimeDeriv_Redist

subroutine Redistance(lo, hi, phi, philo, phihi, &
                      mask, masklo, maskhi, &
                      dx, dtau, tau, &
                      mi_i, mi_j, mi_k, &
                      mask_size) bind(C, name="Redistance")
                       
 use amrex_fort_module, only : amrex_real
 implicit none

 integer, intent(in) ::  lo(3), hi(3)
 integer, intent(in) ::  philo(3), phihi(3), masklo(3), maskhi(3)
 real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 integer, intent(in) :: mask(masklo(1):maskhi(1),masklo(2):maskhi(2),masklo(3):maskhi(3))
 real(amrex_real), intent(in) :: dx(3)
 integer, intent(in) :: mask_size
 integer, dimension(mask_size), intent(in) :: mi_i, mi_j, mi_k
 
 real(amrex_real), intent(in) :: dtau
 real(amrex_real), intent(inout) :: tau 
    
 real(amrex_real) :: d_n  (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 real(amrex_real), dimension(mask_size) :: d_np1
 real(amrex_real), dimension(mask_size) :: dDdt
 integer i,j,k, m
 
 real(amrex_real) maxdiff
    
 !!this way does not work w/ tiling!! ghost values not copied   
 !do m = 1, mask_size !iterate through mask
 ! i = mi_i(m) 
 ! j = mi_j(m)
 ! k = mi_k(m)
 ! d_n(i,j,k) = phi(i,j,k)
 ! !if (mask(i,j,k).EQ.0)then
 ! !  d_n(i,j,k)=10
 ! !endif
 ! !if (mask(i,j,k).EQ.1)then
 ! !  d_n(i,j,k)=20
 ! !endif
 ! !if (mask(i,j,k).EQ.2)then
 ! !  d_n(i,j,k)=30
 ! !endif
 ! !if (mask(i,j,k).EQ.3)then
 ! !  d_n(i,j,k)=40
 ! !endif
 !end do
	
 do k = philo(3),phihi(3)
 do j = philo(2),phihi(2)
 do i = philo(1),phihi(1)
  d_n(i,j,k) = phi(i,j,k)
 enddo
 enddo
 enddo
	
 !tau = 0.0 !pseudotime
 !dtau = dx(1)/3 !redistancing time step restriction to be monotone: (dtau/dx)|s_ij|<=1/2
 !!print *, 'dtau ', dtau
 !iter_redist = 0 !should only take one or two iterations within tube
 !do while (iter_redist < 3)
  tau = tau + dtau
  call TimeDeriv_Redist(dDdt, d_n, philo, phihi, dx, mask, masklo, maskhi, mi_i, mi_j, mi_k, mask_size) !return dphidt
  do m = 1, mask_size !iterate through mask
   i = mi_i(m) 
   j = mi_j(m)
   k = mi_k(m)
   d_np1(m) = d_n(i,j,k)-(dtau)*dDdt(m)
  end do
 
  maxdiff = 0
  do m = 1, mask_size !iterate through mask
   !!!add difference between d_n and d_np1 to check if converged
   i = mi_i(m) 
   j = mi_j(m)
   k = mi_k(m)
   if( abs(d_np1(m)-d_n(i,j,k)) .GT. maxdiff) then
    maxdiff = abs(d_np1(m)-d_n(i,j,k)) 
   endif
   d_n(i,j,k) = d_np1(m)
  end do
  
!  print *, 'redist_diff ', maxdiff
 
 ! iter_redist = iter_redist + 1
 !end do !end while
 
 !!only needed at end of redistancing, doing every iter for debugging/ ease of passing from cpp routine
 do m = 1, mask_size !iterate through mask
  i = mi_i(m) 
  j = mi_j(m)
  k = mi_k(m)
  phi(i,j,k) = d_np1(m)
 end do
	
end subroutine Redistance







subroutine GridNormals (n_grid, nlo, nhi, &
                        phi, philo, phihi, &
                        mi_i, mi_j, mi_k, mask_size, &
                        dx) bind(C, name="GridNormals")
                    
 use amrex_fort_module, only : amrex_real
 implicit none

 integer, intent(in) :: nlo(3), nhi(3), philo(3), phihi(3)
 real(amrex_real), intent(out) :: n_grid(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),1:3)
 real(amrex_real), intent(in) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 integer, intent(in) :: mask_size
 integer, dimension(mask_size), intent(in) :: mi_i, mi_j, mi_k
 real(amrex_real), intent(in) :: dx(3)
 
 integer :: i,j,k, m
 real(amrex_real) :: phi_x, phi_y, phi_z
 real(amrex_real) :: norm_grad_phi
 
 do m = 1, mask_size !iterate through mask
  i = mi_i(m) 
  j = mi_j(m)
  k = mi_k(m)
  !find normal dir (LS)
  phi_x = (phi(i+1,j,k)-phi(i-1,j,k))/(2*dx(1))
  phi_y = (phi(i,j+1,k)-phi(i,j-1,k))/(2*dx(2))
  phi_z = (phi(i,j,k+1)-phi(i,j,k-1))/(2*dx(3))
  norm_grad_phi = sqrt(phi_x**2+phi_y**2+phi_z**2)
  
  n_grid(i,j,k,1) = phi_x/norm_grad_phi
  n_grid(i,j,k,2) = phi_y/norm_grad_phi
  n_grid(i,j,k,3) = phi_z/norm_grad_phi
 end do

end subroutine GridNormals



subroutine compute_flux (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
						 vel_u, ulo, uhi, vel_v, vlo, vhi, vel_w, wlo, whi, &
                         dx) bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(3), hi(3), domlo(3), domhi(3)
  integer philo(3), phihi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
  integer ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2), fxlo(3): fxhi(3))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2), fylo(3): fyhi(3))
  real(amrex_real), intent(inout) :: fluxz( fzlo(1): fzhi(1), fzlo(2): fzhi(2), fzlo(3): fzhi(3))
  real(amrex_real), intent(in   ) :: vel_u( ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3))
  real(amrex_real), intent(in   ) :: vel_v( vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3))
  real(amrex_real), intent(in   ) :: vel_w( wlo(1):whi(1), wlo(2):whi(2), wlo(3):whi(3))
  real(amrex_real), intent(in)    :: dx(3)
  
  ! local variables
  integer i,j,k

  
  do k = lo(3), hi(3)+1
  do j = lo(2), hi(2)+1
  do i = lo(1), hi(1)+1
  if ( (i .LE. hi(1)).AND.(j .LE. hi(2)).AND.(k .LE. hi(3)) ) then 
     fluxx(i,j,k) = 0.5*( phi(i,j,k) + phi(i-1,j,k) ) * vel_u(i,j,k)
     fluxy(i,j,k) = 0.5*( phi(i,j,k) + phi(i,j-1,k) ) * vel_v(i,j,k)
     fluxz(i,j,k) = 0.5*( phi(i,j,k) + phi(i,j,k-1) ) * vel_w(i,j,k)
    elseif ( (i .EQ. hi(1)+1).AND.(j .LE. hi(2)).AND.(k .LE. hi(3))) then
     fluxx(i,j,k) = 0.5*( phi(i,j,k) + phi(i-1,j,k) ) * vel_u(i,j,k)
    elseif ( (i .LE. hi(1)).AND.(j .EQ. hi(2)+1).AND.(k .LE. hi(3))) then
     fluxy(i,j,k) = 0.5*( phi(i,j,k) + phi(i,j-1,k) ) * vel_v(i,j,k)
    elseif ( (i .LE. hi(1)).AND.(j .LE. hi(2)).AND.(k .EQ. hi(3)+1)) then
     fluxz(i,j,k) = 0.5*( phi(i,j,k) + phi(i,j,k-1) ) * vel_w(i,j,k)
    endif
  end do
  end do
  end do


end subroutine compute_flux

subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
                       dx, dt) bind(C, name="update_phi")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(3), hi(3), polo(3), pohi(3), pnlo(3), pnhi(3), &
       fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
  real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2),polo(3):pohi(3))
  real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3))
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
  real(amrex_real), intent(in   ) :: fluxz (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
  real(amrex_real), intent(in)    :: dx(3)
  real(amrex_real), intent(in)    :: dt

  ! local variables
  integer i,j,k
  real(amrex_real) :: dtdx(3)

  dtdx = dt/dx

  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)

     phinew(i,j,k) = phiold(i,j,k) &
          - dtdx(1) * (fluxx(i+1,j  ,k  ) - fluxx(i,j,k)) &
          - dtdx(2) * (fluxy(i  ,j+1,k  ) - fluxy(i,j,k)) &
          - dtdx(3) * (fluxz(i  ,j  ,k+1) - fluxz(i,j,k))

  end do
  end do
  end do

end subroutine update_phi



!!!crossing time

subroutine crossingTime(lo, hi, dx, &
                        phi, philo, phihi, &
                        phi_old, phi_oldlo, phi_oldhi, &
                        dist, distlo, disthi, &
                        mask, masklo, maskhi, &
                        n_grid, nlo, nhi, &
                        mi_i, mi_j, mi_k, maskSize, &
                        tau, dtau, &
                        iter_redist, max_iter) bind(C, name="crossingTime")
                        
            
 use amrex_fort_module, only : amrex_real
 implicit none

 integer, intent(in) ::  lo(3), hi(3)
 integer, intent(in) ::  philo(3), phihi(3), masklo(3), maskhi(3), nlo(3), nhi(3)
 real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 real(amrex_real), intent(in) :: mask(masklo(1):maskhi(1),masklo(2):maskhi(2),masklo(3):maskhi(3))
 real(amrex_real), intent(in) :: n_grid(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),1:3)
 real(amrex_real), intent(in) :: dx(3)
 integer, intent(in) :: maskSize
 integer, dimension(maskSize), intent(in) :: mi_i, mi_j, mi_k
 integer, intent(in) :: iter_redist, max_iter
 
 real(amrex_real), intent(inout) :: tau
 real(amrex_real), intent(in) :: dtau 
 !real(amrex_real) :: dist(maskSize) !to set phi equal to distance
 integer, intent(in) ::  phi_oldlo(3), phi_oldhi(3), distlo(3), disthi(3)
 real(amrex_real), intent(inout) :: phi_old(phi_oldlo(1):phi_oldhi(1),phi_oldlo(2):phi_oldhi(2),phi_oldlo(3):phi_oldhi(3))
 real(amrex_real), intent(inout) :: dist(distlo(1):disthi(1),distlo(2):disthi(2),distlo(3):disthi(3))
 real(amrex_real) :: phi_new(maskSize)!, phi_old(maskSize) !phi values for crossing time approach
 !real(amrex_real) :: dist
 real(amrex_real) :: grad_phi(3)
 real(amrex_real) :: maxdiff
 logical :: flag_distfound
 
 integer i,j,k, m, n
 real(amrex_real) norm_grad_phi, normal
 
 real(amrex_real) phi_ijk, phi_ip1, phi_im1, phi_jp1, phi_jm1, phi_kp1, phi_km1
 real(amrex_real) a, b, c, d, e, f
 real(amrex_real) aplus, bplus, cplus, dplus, eplus, fplus, &
                  aminus, bminus, cminus, dminus, eminus, fminus
 
 !dtau = dx(1)/2
 
 !tau = 0.0
 
 n = iter_redist
 if (n .EQ. 1) then 
  do m = 1, maskSize
   i = mi_i(m); j = mi_j(m); k = mi_k(m)
   dist(i,j,k) = 999942 !phi(i,j,k) !init dist = phi_0 !instead of flag, init dist to impossible int/flag (don't want to overwrite if already found)
   !phi_old(i,j,k) = phi(i,j,k) !done on cpp side 
  enddo 
 endif
 
 !do n = 1, 24 !crossing time iterations !should be no more than max width of narrow band !!fix this
  do m = 1, maskSize !iterate through mask
   i = mi_i(m); j = mi_j(m); k = mi_k(m)
   
   phi_ijk = phi_old(i,j,k)
   phi_ip1 = phi_old(i+1,j,k)
   phi_im1 = phi_old(i-1,j,k)
   phi_jp1 = phi_old(i,j+1,k)
   phi_jm1 = phi_old(i,j-1,k)
   phi_kp1 = phi_old(i,j,k+1)
   phi_km1 = phi_old(i,j,k-1)
   
   !return D_x^(-) and D_x^(+)  --> a, b
   call D_pm_phi(dx(1), phi_im1, phi_ijk, phi_ip1, a, b)
   !return D_y^(-) and D_y^(+)  --> c, d
   call D_pm_phi(dx(2), phi_jm1, phi_ijk, phi_jp1, c, d)
   !return D_z^(-) and D_z^(+)  --> e, f
   call D_pm_phi(dx(3), phi_km1, phi_ijk, phi_kp1, e, f)
   
!   !print *, 'a= ', a
!   if (mask(i-1,j,k) .EQ. 0 .OR.  mask(i-1,j,k) .EQ. 1) then
!      a=0
!    endif
!    if (mask(i+1,j,k) .EQ. 0 .OR.  mask(i+1,j,k) .EQ. 1) then
!      b=0
!    endif
!    if (mask(i,j-1,k) .EQ. 0 .OR.  mask(i,j-1,k) .EQ. 1) then
!      c=0
!    endif
!    if (mask(i,j+1,k) .EQ. 0 .OR.  mask(i,j+1,k) .EQ. 1) then
!      d=0
!    endif
!    if (mask(i,j,k-1) .EQ. 0 .OR.  mask(i,j,k-1) .EQ. 1) then
!      e=0
!    endif
!    if (mask(i,j,k+1) .EQ. 0 .OR.  mask(i,j,k+1) .EQ. 1) then
!      f=0
!    endif
!    
!   call return_PlusMinus(a,aplus,aminus)
!   call return_PlusMinus(b,bplus,bminus)
!   call return_PlusMinus(c,cplus,cminus)
!   call return_PlusMinus(d,dplus,dminus)
!   call return_PlusMinus(e,eplus,eminus)
!   call return_PlusMinus(f,fplus,fminus)
!   
!   norm_grad_phi = sqrt(aplus**2+bminus**2+cplus**2+dminus**2+eplus**2+fminus**2 +aminus**2+bplus**2+cminus**2+dplus**2+eminus**2+fplus**2)
!   
!   !norm_grad_phi = sqrt(MAX(aplus**2,bminus**2)+MAX(cplus**2,dminus**2)+MAX(eplus**2,fminus**2)) &
!   !                 +sqrt(MAX(aminus**2,bplus**2)+MAX(cminus**2,dplus**2)+MAX(eminus**2,fplus**2))
   
   flag_distfound = .FALSE.
   !print *, 'norm_grad_phi= ', norm_grad_phi
   
   grad_phi(1) = (phi_old(i+1,j,k)-phi_old(i-1,j,k))/(2*dx(1))
   grad_phi(2) = (phi_old(i,j+1,k)-phi_old(i,j-1,k))/(2*dx(2))
   grad_phi(3) = (phi_old(i,j,k+1)-phi_old(i,j,k-1))/(2*dx(3))
   
   if ((phi_ijk.GE.phi_ip1 .AND. n.GT.21) .OR. (phi_ijk.LT.phi_ip1 .AND. n.LE.21))then
    grad_phi(1) = (phi_old(i+1,j,k)-phi_old(i,j,k))/dx(1)
   elseif ((phi_ijk.GE.phi_im1 .AND. n.GT.21) .OR. (phi_ijk.LT.phi_im1 .AND. n.LE.21)) then
    grad_phi(1) = (phi_old(i,j,k)-phi_old(i-1,j,k))/dx(1)
   endif
   
   if ((phi_ijk.GE.phi_jp1 .AND. n.GT.21) .OR. (phi_ijk.LT.phi_jp1 .AND. n.LE.21))then
    grad_phi(2) = (phi_old(i,j+1,k)-phi_old(i,j,k))/dx(1)
   elseif ((phi_ijk.GE.phi_jm1 .AND. n.GT.21) .OR. (phi_ijk.LT.phi_jm1 .AND. n.LE.21)) then
    grad_phi(2) = (phi_old(i,j,k)-phi_old(i,j-1,k))/dx(1)
   endif
   
   if ((phi_ijk.GE.phi_kp1 .AND. n.GT.21) .OR. (phi_ijk.LT.phi_kp1 .AND. n.LE.21))then
    grad_phi(3) = (phi_old(i,j,k+1)-phi_old(i,j,k))/dx(1)
   elseif ((phi_ijk.GE.phi_km1 .AND. n.GT.21) .OR. (phi_ijk.LT.phi_km1 .AND. n.LE.21)) then
    grad_phi(3) = (phi_old(i,j,k)-phi_old(i,j,k-1))/dx(1)
   endif
  ! if (( mask(i,j-1,k) .EQ. 0  .OR.  mask(i,j-1,k) .EQ. 1 ) .OR. &
  !     (phi_ijk.GT.phi_jp1 .AND. phi_ijk.GT.0) .OR. (phi_ijk.LT.phi_jp1 .AND. phi_ijk.LT.0)) then
  !  grad_phi(2) = (phi_old(i,j+1,k)-phi_old(i,j,k))/dx(2)
  ! elseif ((mask(i,j+1,k) .EQ. 0  .OR.  mask(i,j+1,k) .EQ. 1 ) .OR. &
  !     (phi_ijk.GT.phi_jm1 .AND. phi_ijk.GT.0) .OR. (phi_ijk.LT.phi_jm1 .AND. phi_ijk.LT.0)) then 
  !  grad_phi(2) = (phi_old(i,j,k)-phi_old(i,j-1,k))/dx(2)
  ! endif
  ! if (( mask(i,j,k-1) .EQ. 0  .OR.  mask(i,j,k-1) .EQ. 1 ) .OR. &
  !     (phi_ijk.GT.phi_kp1 .AND. phi_ijk.GT.0) .OR. (phi_ijk.LT.phi_kp1 .AND. phi_ijk.LT.0)) then
  !  grad_phi(3) = (phi_old(i,j,k+1)-phi_old(i,j,k))/dx(3)
  ! elseif ((mask(i,j,k+1) .EQ. 0  .OR.  mask(i,j,k+1) .EQ. 1 ) .OR. &
  !     (phi_ijk.GT.phi_km1 .AND. phi_ijk.GT.0) .OR. (phi_ijk.LT.phi_km1 .AND. phi_ijk.LT.0)) then
  !  grad_phi(3) = (phi_old(i,j,k)-phi_old(i,j,k-1))/dx(3)
  ! endif
   norm_grad_phi = sqrt(grad_phi(1)**2+grad_phi(2)**2+grad_phi(3)**2)
   !if ((n .EQ. 1) .AND. (n_grid(i,j,k,1) .NE. grad_phi(1)/norm_grad_phi))then
   ! print *, "n_grid(i,j,k)", n_grid(i,j,k,:), "grad_phi, ", grad_phi/norm_grad_phi
   !endif
   
   !normal = -grad_phi/|grad_phi|
   !note : normal*grad_phi = |grad_phi|
   
   if (n.LE.21) then !phi(i,j,k) .LT. 0.0 .AND.  !phi<0
    phi_new(m) = phi_old(i,j,k)+norm_grad_phi*dtau !-DOT_PRODUCT(-grad_phi/norm_grad_phi,grad_phi)*dtau !+norm_grad_phi*dtau
    !print *, "phi_old: ", phi_old(i,j,k), " phi_new: ", phi_new(m)
    if (phi_new(m) .GE. 0.0 .AND. phi_old(i,j,k) .LT. 0.0 .AND. dist(i,j,k) .EQ. 999942) then
     dist(i,j,k) = -tau - dtau*phi_old(i,j,k)/(phi_old(i,j,k)-phi_new(m))
     flag_distfound = .TRUE.
     !print *, "phi_0 : ", phi(i,j,k), " , dist: ", dist
    ! print *, 'phi < 0, dist found, n= ', n
    endif
   elseif (n.GT.21 .AND. 1==1) then !phi(i,j,k) .GT. 0.0 .AND.   !phi>0
    phi_new(m) = phi_old(i,j,k)-norm_grad_phi*dtau !+DOT_PRODUCT(-grad_phi/norm_grad_phi,grad_phi)*dtau !-norm_grad_phi*dtau
    !print *, "phi_old: ", phi_old(i,j,k), " phi_new: ", phi_new(m)
    if (phi_new(m) .LE. 0.0 .AND. phi_old(i,j,k) .GT. 0.0 .AND. dist(i,j,k) .EQ. 999942) then
     dist(i,j,k) = tau + dtau*phi_old(i,j,k)/(phi_old(i,j,k)-phi_new(m))
     flag_distfound = .TRUE.
     !print *, "phi_0 : ", phi(i,j,k), " , dist: ", dist(i,j,k)
    ! print *, 'phi > 0, dist found, n= ', n
    endif
   else
   ! phi_new(m) = phi_old(i,j,k)
   endif ! end if phi<0, phi>0
   !print *, "next "
   
  enddo ! end m
  
  !update phi_old !do not combine w/ previous loop (can't update phi_old on the fly)
  do m = 1, maskSize
   i = mi_i(m); j = mi_j(m); k = mi_k(m)
   phi_old(i,j,k) = phi_new(m)
  enddo ! end m
  tau = tau + dtau
 !enddo !end n !redistancing iterations
  
 if (iter_redist .EQ. max_iter) then 
  maxdiff = 0.0
  do m = 1, maskSize
   i = mi_i(m); j = mi_j(m); k = mi_k(m)
   if( abs(dist(i,j,k)-phi(i,j,k)) .GT. maxdiff) then
    maxdiff = abs(dist(i,j,k)-phi(i,j,k)) 
   endif
   
   if (dist(i,j,k).EQ.999942)then
    dist(i,j,k) = phi(i,j,k)
   endif
   
   !print *, "phi: ", phi(i,j,k), " dist: ", dist(m) !add check to see that sign of update stays the same
   if ( (phi(i,j,k) .GT. 0.0 .AND. dist(i,j,k) .LT. 0.0) .OR. (phi(i,j,k) .LT. 0.0 .AND. dist(i,j,k) .GT. 0.0) ) then 
    print *, "phi changed signs, phi: ", phi(i,j,k), " dist: ", dist(i,j,k)
   endif
   phi(i,j,k) = dist(i,j,k)
  enddo
  print *, 'redist_diff ', maxdiff
 endif

end subroutine crossingTime
 
 !do m = 1, maskSize !iterate through mask
 ! i = mi_i(m) 
 ! j = mi_j(m)
 ! k = mi_k(m)
 ! phi_0 = phi(i,j,k)
 ! phi_old = phi_0
 ! dist = 0.0
 ! 
 ! gradphi(1) = (phi(i+1,j,k)-phi(i-1,j,k))/(2*dx(1))
 ! gradphi(2) = (phi(i,j+1,k)-phi(i,j-1,k))/(2*dx(2))
 ! gradphi(3) = (phi(i,j,k+1)-phi(i,j,k-1))/(2*dx(3))
 ! norm_grad_phi = sqrt(gradphi(1)**2+gradphi(2)**2+gradphi(3)**2)
 ! 
 ! !note : n*gradphi = |gradphi|
 ! 
 ! tau = 0.0
 ! if (phi_0 .LT. 0.0) then
 !  do n = 1, 720000 !!change this --max need half width of NB //6dx for 5th order WENO, 4dx for 3rd order WENO, use mask
 !   phi_new = phi_old+norm_grad_phi*dtau!+DOT_PRODUCT(n_grid(i,j,k,:),gradphi)*dtau
 !   print *, "phi_old: ", phi_old, " phi_new: ", phi_new
 !   !print *, "n_grid(i,j,k)", n_grid(i,j,k,:), "gradphi, ", gradphi/norm_grad_phi
 !   if (phi_new .GE. 0.0 .AND. phi_old .LT. 0.0) then
 !    dist = tau - dtau*phi_old/(phi_old-phi_new)
 !    !print *, 'phi < 0, dist found, n= ', n
 !    EXIT
 !   endif
 !   !if (n .GE. 720000) then
 !   ! print *, "uh oh "
 !   ! !print *, "n_grid(i,j,k)", n_grid(i,j,k,:), "gradphi, ", gradphi/norm_grad_phi
 !   !endif
 !   phi_old = phi_new
 !   tau = tau - dtau
 !  enddo 
 ! elseif (phi_0 .GT. 0.0) then
 !  do n = 1, 720000 !!change this --max need half width of NB
 !   phi_new = phi_old-norm_grad_phi*dtau!-DOT_PRODUCT(n_grid(i,j,k,:),gradphi)*dtau
 !   !print * , "phi " , phi_new
 !   if (phi_new .LE. 0.0 .AND. phi_old .GT. 0.0) then
 !    dist = tau + dtau*phi_old/(phi_old-phi_new)
 !    !print *, 'phi > 0, dist found, n= ', n
 !    EXIT
 !   endif
 !   !if (n .GE. 240) then
 !   ! print *, "uh oh"
 !   !endif
 !   phi_old = phi_new
 !   tau = tau + dtau
 !  enddo 
 ! endif
 ! !print *, "next "
 ! 
 ! phi_update(m) = dist !update phi
 !enddo !end m maskSize
 