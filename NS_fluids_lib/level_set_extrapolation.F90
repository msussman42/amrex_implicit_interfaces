  MODULE ZEYU_LS_extrapolation
  IMPLICIT NONE

  contains

! calculate Givens rotation matrix
!    _     _ T _  _     _ _
!   |  c  s | |  a |   | r |
!                    = 
!   |_-s  c_| |_ b_|   |_0_|
!
  subroutine givens(a, b, c, s)
    implicit none

    double precision, intent(in ) :: a, b
    double precision, intent(out) :: c, s

    if (abs(b) .le. 1.d-16) then
       c = 1.d0
       s = 0.d0
    else
       if (abs(b) .gt. abs(a)) then
          s = 1.d0 / sqrt(1.d0 + (a / b)**2)
          c = s * (-a / b)
       else
          c = 1.d0 / sqrt(1.d0 + (b / a)**2)
          s = c * (-b / a)
       endif
    endif

  end subroutine givens

! QR factorization for least squares problems
! A: m x n, x: sn, b: sm
! To solve min||Ax-b||^2 for x:
! 1. A = QR;
! 2. d=QTb;
! 3. solve Rx=d.
  subroutine least_squares_QR(A, x, b, m, n)
    implicit none

    integer,          intent(in ) :: m, n
    double precision, intent(in ) :: A(m, n), b(m)
    double precision, intent(out) :: x(n)

    integer i, j, k, it
    double precision R(m, n+1), y, z, c, s, temp1, temp2
    double precision ATA(n, n), residual_verify, ATAx(n), ATb(n)

    do i = 1, m
       do j = 1, n + 1
          if (j .eq. n + 1) then
             R(i, j) = b(i)
          else
             R(i, j) = A(i, j)
          endif
       enddo
    enddo

    if (1 .eq. 0) then
       print *, "original R = "
       do i = 1, m
          print "(5e13.5)", (R(i, j), j =1, n+1)
       enddo
    endif

    do it = 1, n
       do i = it + 1, m
          y = R(it, it)
          z = R(i, it)
          if (abs(z) .lt. 1.d-16) then
             cycle
          else
             call givens(y, z, c, s)
             do j = it, n + 1
                temp1 = c * R(it, j) + (-s) * R(i, j)
                temp2 = s * R(it, j) + c * R(i, j)
                R(it, j) = temp1
                R(i, j) = temp2
             enddo
          endif
       enddo
    enddo

    if (1 .eq. 0) then
       print *, "R = "
       do i = 1, m
          print "(5e13.5)", (R(i, j), j =1, n+1)
       enddo
    endif

    if (abs(R(n, n)) .le. 1.d-16) then
       print *, "Warning! R(", n, ",", n, &
                ") is close to zero in QR fatorization!"
       R(n, n) = 1.d-16
    endif
    x(n) = R(n, n+1) / R(n, n)
    do i = 1, n-1
       x(n-i) = R(n-i, n+1)
       do j = n-i+1, n
          x(n-i) = x(n-i) - R(n-i, j) * x(j)
       enddo
       if (abs(R(n-i, n-i)) .le. 1.d-16) then
          print *, "Warning! R(", n-i, ",", n-i, &
                   ") is close to zero in QR fatorization!"
          R(n-i, n-i) = 1.d-16
       endif
       x(n-i) = x(n-i) / R(n-i, n-i)
    enddo

! sanity check
    do i = 1, n
       do j = 1, n
          ATA(i, j) = 0.d0
          do k = 1, m
             ATA(i, j) = ATA(i, j) + A(k, i) * A(k, j)
          enddo
       enddo
    enddo

    residual_verify = 0.d0
    do i = 1, n
       ATAx(i) = 0.d0
       ATb(i) = 0.d0
       do j = 1, n
          ATAx(i) = ATAx(i) + ATA(i, j) * x(j)
          ATb(i) = ATb(i) + A(j, i) * b(j)
       enddo
       do j = n+1, m
          ATb(i) = ATb(i) + A(j, i) * b(j)
       enddo
       residual_verify = residual_verify + (ATAx(i) - ATb(i))**2
    enddo

    residual_verify = sqrt(residual_verify) / n
    if (residual_verify .le. 1.d-8) then
       ! do nothing
    else
       print *, "Error! ||ATAx - ATb|| = ", residual_verify
       stop
    endif
! end sanity check

  end subroutine least_squares_QR

! if 2D, nk should be 1, pos_z should be 0
! inputs:
!   nmat = number of materials
!   pos_xyz(ni,nj,nk,dim)
!   ls(ni,nj,nk,nmat)
!   weights(ni,nj,nk)
!   is_fluid(nmat)   =  1 if fluid and needs to be extrapolated, 0 otherwise
!    is_fluid(im)=1 if material "im" is a fluid material (i.e. gets extrap)
!    is_fluid(im)=0 if material "im" is not a fluid (i.e. does not get extrap)
! output:
!   ls_extrap(im=1...nmat) = extrapolated distance for materials in which 
!                            is_fluid(im).eq.1
  subroutine level_set_extrapolation(pos_xyz, ls, weights, is_fluid, &
                                     ls_extrap, ni, nj, nk, nmat, dim)
    implicit none

    integer,          intent(in ) :: ni, nj, nk, nmat, dim
    integer,          intent(in ) :: is_fluid(nmat)
    double precision, intent(in ) :: pos_xyz(ni, nj, nk, dim), &
                                     ls(ni, nj, nk, nmat), &
                                     weights(ni, nj, nk)
    double precision, intent(out) :: ls_extrap(nmat)

    integer i, j, k, m, n, im
    double precision A(ni*nj*nk,dim+1), b(ni*nj*nk), pos0(dim), var(dim+1)
    !ls_extrap = a . (x - x0) + b
    !find a(3) and b while minimize SUM(w_ij*(phi_ij-phi_fluid_ij)^2)
    !calculate phi_fluid at x0
    !var(dim+1) = {a(dim), b}

    if (((ni+1)/2)*2.eq.ni+1) then
     ! do nothing
    else
     print *,"ni invalid"
     stop
    endif
    if (((nj+1)/2)*2.eq.nj+1) then
     ! do nothing
    else
     print *,"nj invalid"
     stop
    endif
    if (((nk+1)/2)*2.eq.nk+1) then
     ! do nothing
    else
     print *,"nk invalid"
     stop
    endif
    if (ni.eq.nj) then
     ! do nothing
    else 
     print *,"ni or nj invalid"
     stop
    endif
    if (dim.eq.3) then
     if (ni.eq.nk) then
      ! do nothing
     else 
      print *,"ni or nk invalid"
      stop
     endif
    else if (dim.eq.2) then
     if (nk.eq.1) then
      ! do nothing
     else
      print *,"nk invalid"
      stop
     endif
    else
     print *,"dim invalid"
     stop
    endif

    pos0(1) = pos_xyz((ni+1)/2, (nj+1)/2, (nk+1)/2, 1)
    pos0(2) = pos_xyz((ni+1)/2, (nj+1)/2, (nk+1)/2, 2)
    if (dim .eq. 3) then
       pos0(dim) = pos_xyz((ni+1)/2, (nj+1)/2, (nk+1)/2, 3)
    endif
    m = ni * nj * nk
    n = dim + 1
    do i = 1, m
       do j = 1, n
          A(i, j) = 0.d0
       enddo
       b(i) = 0.d0
    enddo
    do i = 1, n
       var(i) = 0.d0
    enddo
    do im = 1, nmat
       ls_extrap(im) = 0.d0
    enddo

    do im = 1, nmat
       if (is_fluid(im) .eq. 1) then
          do i = 1, ni
             do j = 1, nj
                do k = 1, nk
                   if (weights(i,j,k) .lt. 0.d0) then
                      print *, "Error! Weight factor less than zero!"
                      stop
                   endif
                   A((i-1)*nj*nk+(j-1)*nk+k,1) = weights(i,j,k) &
                                         * (pos_xyz(i,j,k,1)-pos0(1))
                   A((i-1)*nj*nk+(j-1)*nk+k,2) = weights(i,j,k) &
                                         * (pos_xyz(i,j,k,2)-pos0(2))
                   if (dim .eq. 3) then
                      A((i-1)*nj*nk+(j-1)*nk+k,dim) = weights(i,j,k) &
                                          * (pos_xyz(i,j,k,dim)-pos0(dim))
                   endif
                   A((i-1)*nj*nk+(j-1)*nk+k, dim+1) = 1.d0
                   b((i-1)*nj*nk+(j-1)*nk+k) = sqrt(weights(i,j,k)) &
                                             * ls(i,j,k,im)
                enddo
             enddo       
          enddo

          call least_squares_QR(A, var, b, m, n)

          ls_extrap(im) = var(dim + 1)
       else!is_fluid(im) = 0, soild
          ls_extrap(im) = ls((ni+1)/2, (nj+1)/2, (nk+1)/2, im)
       endif
    enddo ! im=1..nmat

  end subroutine level_set_extrapolation

  END MODULE ZEYU_LS_extrapolation


