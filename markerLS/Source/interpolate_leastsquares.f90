!!TODO: define global consts in a module -> eg SDIM, PI, Tol, zero one half, etc
subroutine polyMatrixDim(numrows, numcols, polyOrder) bind(C, name="polyMatrixDim")
 !polyMatrixDim finds the total number of rows and columns in polynomial matrix/stencil P for memory allocation
 use amrex_fort_module, only : amrex_spacedim
 implicit none
 
 integer, intent(in) :: polyOrder
 integer, intent(out) :: numrows, numcols
 integer :: numcolsUpdate, colmodify

 numrows = (2*CEILING(polyOrder/2.0)+3)**AMREX_SPACEDIM !as many rows as support cells (2 extra to avoid illconditioned)
 numcols = (polyOrder+1)*(polyOrder+2)/2.0 !(polyOrder+1)^SDIM is full size, modify to reduce terms in summation: s+t+r<=p
 print *, "polyOrder : ", polyOrder
 print *,"2drows: ", numrows, "2dcols: ",numcols
 if (AMREX_SPACEDIM.eq.3) then
  numcolsUpdate = numcols
  do colmodify = (polyOrder+1),2,-1
   numcolsUpdate = numcolsUpdate-colmodify
   numcols = numcols + numcolsUpdate
  enddo
 endif
 
 print *,"rows: ", numrows, "cols: ",numcols
end subroutine polyMatrixDim


subroutine init_Pmatrix(polyOrder, numrows, numcols, dx, w_poly, P) bind(C, name="init_Pmatrix")
 !initialize the polynomial matrix P
 use amrex_fort_module, only : amrex_real, amrex_spacedim
 implicit none

 integer, intent(in) :: polyOrder
 integer, intent(in) :: numrows, numcols
 real(amrex_real), dimension(3), intent(in) :: dx
 
 real(amrex_real), dimension(numrows,numcols), intent(out) :: P !polynomial matrix
 real(amrex_real), dimension(numrows), intent(out) :: w_poly !stencil weights
 
 integer :: row, col, shift
 integer :: r, s, t, supportcellsX, supportcellsY, supportcellsZ
 real(amrex_real) :: sumweights
 
 integer :: i,j

 sumweights = 0.0
 row = 1
 col = 1
 shift = CEILING(polyOrder/2.0) !polyOrder+1 support cells in xyz dirs
 do supportcellsZ = -(shift+1),shift+1 !the '+1' gives polyOrder+3 to avoid ill-conditioned
  do supportcellsY = -(shift+1),shift+1
   do supportcellsX = -(shift+1),shift+1
    if (AMREX_SPACEDIM.eq.2) then
     w_poly(row) = 22.0**((-abs(supportcellsX)+shift) + (-abs(supportcellsY)+shift))
     sumweights = sumweights + w_poly(row)
    elseif (AMREX_SPACEDIM.eq.3) then
     w_poly(row) = 22.0**((-abs(supportcellsX)+shift) + (-abs(supportcellsY)+shift) + (-abs(supportcellsZ)+shift))
     sumweights = sumweights + w_poly(row)
    endif
    do r = 0,polyOrder !! can reduce iterations, but not priority since P created only once
     do t = 0,polyOrder
      do s = 0,polyOrder
       if (s+t+r > polyOrder) then
        CYCLE
       else
        !print*, "row: ", row, "col: ", col
        P(row,col) = (supportcellsX*dx(1))**s * (supportcellsY*dx(2))**t * (supportcellsZ*dx(3))**r
        col = col+1
       endif
      enddo !end s
     enddo !end t
     if (AMREX_SPACEDIM.eq.2) then
      EXIT !r=0 for 2d 
     endif
    enddo !end r
    row = row+1
    col = 1
   enddo !end supportcellsX
  enddo !end supportcellsY
  if (AMREX_SPACEDIM.eq.2) then
   EXIT !skip Z if 2d
  endif
 enddo !end supportcellsZ
 
 print*, "w_poly"
 do i = 1, numrows
  !w_poly(i) = w_poly(i)/sumweights
  write(*,"(10E14.5)") (w_poly(i))
 enddo
 
 print*, "P"
 do i = 1, numrows
  write(*,"(10E14.5)") (P(i,j), j=1,numcols)
 enddo
end subroutine init_Pmatrix

subroutine init_LeastSquaresMatrix(numrows, numcols, w, P, A_inv) bind(C, name="init_LeastSquaresMatrix")
 !store matrix A^-1 for solving least squares problem Ax=b for the polynomial matrix
 !!currently using Gaussian Elimination to get inverse -> switch to Householder QR
 !minimization problem sum_{i=1}^{num_supportcells} w_i*( P(x_i)-phi_i )^2
 use amrex_fort_module, only : amrex_real
 implicit none
 
 integer, intent(in) :: numrows, numcols
 real(amrex_real), dimension(numrows,numcols), intent(in) :: P !polynomial matrix
 real(amrex_real), dimension(numrows), intent(in) :: w !stencil weights
 
 real(amrex_real), dimension(numcols,numcols), intent(out) :: A_inv
 
 real(amrex_real), dimension(numcols,2*numcols) :: A_ge !augmented matrix A_ge = [A I] for gaussian elimination
 real(amrex_real) :: sumval
 integer :: r,s !row, col index of A
 integer :: ii !supporting cells index
 integer :: Prow,Pcol !pivot row, col index
 integer :: m,n !augmented matrix dim
 real(amrex_real) :: swap_temp !temporary variable for swapping rows
 integer :: i,j !row, col index augmented matrix
 integer :: i_max !index of Mmaxval
 real(amrex_real) :: Mmaxval !matrix max value
 real(amrex_real) :: multval, rowdiv !value to multiply/divide row by
 
 !A_ge = augmented matrix [A I]
 do r = 1,numcols
  do s = 1,numcols
   !init A
   sumval = 0.0
   do ii=1,numrows
    sumval = sumval + w(ii)*P(ii,r)*P(ii,s)
   enddo
   A_ge(r,s) = sumval
   !init I
   if (s.eq.r) then
    A_ge(r,s+numcols) = 1.0
   else
    A_ge(r,s+numcols) = 0.0
   endif
  enddo
 enddo
 
 !obtaining inverse by gaussian elimination
 !A_ge = augmented matrix [A I] -> [I A^-1]
 !refer to pseudocode https://en.wikipedia.org/wiki/Gaussian_elimination
 Prow = 1 !init pivot row, col
 Pcol = 1
 m = numcols !num matrix rows
 n = 2*numcols !num matrix cols
 swap_temp = 0.0
 do while (Prow.le.m .and. Pcol.le.n)
  !find the k-th pivot (using largest abs value helps with numerical stability)
  i_max = 0
  Mmaxval = 0.0
  do i = Prow,m
   if ( Mmaxval .lt. abs(A_ge(i, Pcol)) ) then
    Mmaxval = abs(A_ge(i, Pcol))
    i_max = i
   endif
  enddo
  if (Mmaxval.eq.0) then
   !no pivot in this column, pass to next column 
   Pcol = Pcol+1
  else
   !swap rows: Prow <--> i_max
   do j=1,n
    swap_temp=A_ge(Prow,j)
    A_ge(Prow,j)=A_ge(i_max,j)
    A_ge(i_max,j)=swap_temp
   enddo
   !do for all rows below pivot: 
   do i = Prow+1,m
    multval = A_ge(i,Pcol) / A_ge(Prow,Pcol)
    !fill with zeros the lower part of pivot column: 
    A_ge(i,Pcol) = 0
    !do for all remaining elements in current row: 
    do j = Pcol+1,n
     A_ge(i,j) = A_ge(i,j) - A_ge(Prow,j)*multval
    enddo !end j
   enddo !end i
   !increase pivot row and column index
   Prow = Prow+1
   Pcol = Pcol+1
  endif
 enddo !end do while
 
 !currently only row echelon, sweep in reverse to obtain RRE
 Prow = numcols
 Pcol = numcols
 do while (Prow.ge.1 .and. Pcol.ge.1)
  !do for all rows above pivot: 
  if (A_ge(Prow,Pcol).eq.0) then
   !no pivot in this column, pass to next column 
   Pcol = Pcol-1
  else
   !divide out current row by pivot
   rowdiv = A_ge(Prow,Pcol)
   A_ge(Prow,Pcol)=1
   do j = Pcol+1,n
    A_ge(Prow,j)=A_ge(Prow,j)/rowdiv
   enddo
   do i = Prow-1,1,-1
    multval = A_ge(i,Pcol)
    !fill with zeros the upper part of pivot column: 
    A_ge(i,Pcol) = 0
    !do for all remaining elements in current row: 
    if (abs (multval) .gt. 1.0D-8) then !!replace with TOL?
     do j = Pcol+1,n
      A_ge(i,j) = A_ge(i,j) - A_ge(Prow,j)*multval
     enddo !end j
    endif
   enddo !end i
   !decrease pivot row and column 
   Prow = Prow-1
   Pcol = Pcol-1
  endif
 enddo !end do while
 
 !extract A^-1 from the augmented matrix
 do j=1,numcols
  do i=1,numcols
   A_inv(i,j) = A_ge(i,j+numcols)
  enddo
 enddo

 print*, "A_inv"
 do i = 1, numcols
  write(*,"(10E14.5)") (A_inv(i,j), j=1,numcols)
 enddo
end subroutine init_LeastSquaresMatrix


subroutine getPolyInterpCoeffs(polyOrder, numrows, numcols, dx, w, P, A_inv, &
                               phi, philo, phihi, indexVec, a_coeffs) bind(C, name="C_getPolyInterpCoeffs")
 !return the coefficients for the polynomial interpolation of phi
 use amrex_fort_module, only : amrex_real, amrex_spacedim
 implicit none
 
 integer, intent(in) :: polyOrder, numrows, numcols
 real(amrex_real), dimension(3), intent(in) :: dx
 integer, intent(in) ::  philo(3), phihi(3)
 real(amrex_real), intent(in) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
 real(amrex_real), dimension(numrows,numcols), intent(in) :: P !polynomial matrix
 real(amrex_real), dimension(numrows), intent(in) :: w !stencil weights
 real(amrex_real), dimension(numcols,numcols), intent(in) :: A_inv 
 real(amrex_real), dimension(numcols), intent(out) :: a_coeffs !coefficients of poly. approx. to phi
 integer, dimension(3), intent(in) :: indexVec

 !local variables
 real(amrex_real), dimension(numrows) :: phivec !values of phi in stencil
 real(amrex_real), dimension(numcols) :: b !values of phi in stencil
 integer row, shift, i,j,k, r, ii
 integer xcoord, ycoord, zcoord
 real(amrex_real) sumval
 
 i = indexVec(1); j = indexVec(2); k = indexVec(3)
 
 !init phivector
 row = 1
 shift = CEILING(polyOrder/2.0) !polyOrder+1 support cells in xyz dirs
 do zcoord = -(shift+1), shift+1
  do ycoord = -(shift+1), shift+1
   do xcoord = -(shift+1), shift+1
    if (AMREX_SPACEDIM.EQ.2) then
     phivec(row) = phi(i+xcoord, j+ycoord, k)
    elseif (AMREX_SPACEDIM.EQ.3) then
     phivec(row) = phi(i+xcoord, j+ycoord, k+zcoord)
    end if
    row = row+1
   end do!end xcoord
  end do!end ycoord
  if (AMREX_SPACEDIM.EQ.2) then
   EXIT !skip Z if 2d
  end if
 end do!end zcoord
 
 !init b vector
 do r = 1, numcols
  sumval = 0
  do ii = 1, numrows
   sumval = sumval + w(ii)*phivec(ii)*P(ii,r)
  end do
  b(r) = sumval
 end do
                
 !!a_coeffs = A_inv*b !solve P*a=phivec for a
 !find a_coeffs
 do i = 1, numcols
  a_coeffs(i) = 0.0
  do j = 1, numcols
   a_coeffs(i) = a_coeffs(i) + A_inv(i,j)*b(j)
  end do !end do j
 end do !end do i
 
 !!test
 !!!a_coeffs(1) = phi(i,j,k)
 !!print*, "phi_ijk", phi(i,j,k), "a_coeffs", a_coeffs(1)
 
 !print*, "a_coeffs"
 !do i = 1, numcols
 ! write(*,"(10E14.5)") (a_coeffs(i))
 !enddo
 
 !print*, "phivec"
 !do i = 1, numrows
 ! write(*,"(10E14.5)") (phivec(i))
 !enddo
 
 !print*, "Phi_ijk - a_coeffs(1): ", phi(i,i,k)-a_coeffs(1)

end subroutine getPolyInterpCoeffs


!function AlmostEqual(a, b) result(almostEq) !IN: a,b   OUT: bool
! use amrex_fort_module, only : amrex_real
! implicit none
! 
! double precision :: TOL_EPSILON
! 
! real(amrex_real), intent(in) :: a, b
! logical :: almostEq
! 
! TOL_EPSILON=10E-8
!
! if (a.eq.0 .or. b.eq.0) then
!  if (abs(a-b) .le. 2.0*TOL_EPSILON) then
!   almostEq = .true.
!  else
!   almostEq = .false.
!  endif
! else
!  if (abs(a-b).le.TOL_EPSILON*abs(a) .and. abs(a-b).le.TOL_EPSILON*abs(b)) then
!   almostEq = .true.
!  else
!   almostEq = .false.
!  endif
! endif
!end function AlmostEqual



