! --------------------------------------------------------------------
!  Flow Physics and Computation Lab
!   
!
!  Matrix solver for FSI modules in
!  VICAR3D, a viscous, Cartesian, 3D flow solver.
!
!  This is a contineously developing project.
!
!  Starting Developers:
!  Kourosh Shoele
!
!
!  Final Filename: mgmresVER3.F90
!  Latest Modification: June, 01 2016 
!  by Kourosh Shoele
! --------------------------------------------------------------------



subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input, real*8 X(N), the vector to be multiplied by A'.
!
!    Output, real*8 W(N), the value of A'*X.
!
  implicit none

  integer  n
  integer  nz_num

  real*8 a(nz_num)
  integer  i
  integer  ia(n+1)
  integer  ja(nz_num)
  integer  k
  integer  k1
  integer  k2
  real*8 w(n)
  real*8 x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(ja(k1:k2)) = w(ja(k1:k2)) + a(k1:k2) * x(i)
  end do

  return
end
subroutine atx_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! ATX_ST computes A'*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input, real*8 X(N), the vector to be multiplied by A'.
!
!    Output, real*8 W(N), the value of A'*X.
!
  implicit none

  integer  n
  integer  nz_num

  real*8 a(nz_num)
  integer  i
  integer  ia(nz_num)
  integer  j
  integer  ja(nz_num)
  integer  k
  real*8 w(n)
  real*8 x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(j) = w(j) + a(k) * x(i)
  end do

  return
end
subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_CR computes A*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input, real*8 X(N), the vector to be multiplied by A.
!
!    Output, real*8 W(N), the value of A*X.
!
  implicit none

  integer  n
  integer  nz_num

  real*8 a(nz_num)
  integer  i
  integer  ia(n+1)
  integer  ja(nz_num)
  integer  k
  integer  k1
  integer  k2
  real*8 w(n)
  real*8 x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

  return
end
subroutine ax_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_ST computes A*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input, real*8 X(N), the vector to be multiplied by A.
!
!    Output, real*8 W(N), the value of A*X.
!
  implicit none

  integer  n
  integer  nz_num

  real*8 a(nz_num)
  integer  i
  integer  ia(nz_num)
  integer  j
  integer  ja(nz_num)
  integer  k
  real*8 w(n)
  real*8 x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(i) = w(i) + a(k) * x(j)
  end do

  return
end
subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

!*****************************************************************************80
!
!! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!    The array UA can be used to locate the diagonal elements of the matrix.
!
!    It is assumed that every row of the matrix includes a diagonal element,
!    and that the elements of each row have been ascending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!    On output, the order of the entries of JA may have changed because of
!    the sorting.
!
!    Output, integer  UA(N), the index of the diagonal element
!    of each row.
!
  implicit none

  integer  n
  integer  nz_num

  integer  i
  integer  ia(n+1)
  integer  k
  integer  ja(nz_num)
  integer  ua(n)

  ua(1:n) = -1

  do i = 1, n
    do k = ia(i), ia(i+1) - 1
      if ( ja(k) == i ) then
        ua(i) = k
      end if
    end do
  end do

  return
end
subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )

!*****************************************************************************80
!
!! ILU_CR computes the incomplete LU factorization of a matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input, integer  UA(N), the index of the diagonal element
!    of each row.
!
!    Output, real*8 L(NZ_NUM), the ILU factorization of A.
!
  implicit none

  integer  n
  integer  nz_num

  real*8 a(nz_num)
  integer  i
  integer  ia(n+1)
  integer  iw(n)
  integer  j
  integer  ja(nz_num)
  integer  jj
  integer  jrow
  integer  jw
  integer  k
  real*8 l(nz_num)
  real*8 tl
  integer  ua(n)

!
!  Copy A.
!
  l(1:nz_num) = a(1:nz_num)

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i+1) - 1
      iw(ja(k)) = k
    end do

    do j = ia(i), ia(i+1) - 1
      jrow = ja(j)
      if ( i <= jrow ) then
        exit
      end if
      tl = l(j) * l(ua(jrow))
      l(j) = tl
      do jj = ua(jrow) + 1, ia(jrow+1) - 1
        jw = iw(ja(jj))
        if ( jw /= -1 ) then
          l(jw) = l(jw) - tl * l(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow /= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( l(j) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    l(j) = 1.0D+00 / l(j)

  end do

  l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

  return
end
subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

!*****************************************************************************80
!
!! LUS_CR applies the incomplete LU preconditioner.
!
!  Discussion:
!
!    The linear system M * Z = R is solved for Z.  M is the incomplete
!    LU preconditioner matrix, and R is a vector supplied by the user.
!    So essentially, we're solving L * U * Z = R.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 L(NZ_NUM), the matrix values.
!
!    Input, integer  UA(N), the index of the diagonal element
!    of each row.
!
!    Input, real*8 R(N), the right hand side.
!
!    Output, real*8 Z(N), the solution of the system M * Z = R.
!
  implicit none

  integer  n
  integer  nz_num

  integer  i
  integer  ia(n+1)
  integer  j
  integer  ja(nz_num)
  real*8 l(nz_num)
  real*8 r(n)
  integer  ua(n)
  real*8 w(n)
  real*8 z(n)
!
!  Copy R in.
!
  w(1:n) = r(1:n)
!
!  Solve L * w = w where L is unit lower triangular.
!
  do i = 2, n
    do j = ia(i), ua(i) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
  end do
!
!  Solve U * w = w, where U is upper triangular.
!
  do i = n, 1, -1
    do j = ua(i) + 1, ia(i+1) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
    w(i) = w(i) / l(ua(i))
  end do
!
!  Copy Z out.
!
  z(1:n) = w(1:n)

  return
end
subroutine mgmres_st ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, tol_abs, &
  tol_rel )

!*****************************************************************************80
!
!! MGMRES_ST applies restarted GMRES to a sparse triplet matrix.
!
!  Discussion:
!
!    The linear system A*X=B is solved iteratively.
!
!    The matrix A is assumed to be stored in sparse triplet form.  Only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the linear system.
!
!    Input, integer  NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer  IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input/output, real*8 X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real*8 RHS(N), the right hand side of the linear system.
!
!    Input, integer  ITR_MAX, the maximum number of (outer)
!    iterations to take.
!
!    Input, integer  MR, the maximum number of (inner) iterations
!    to take.  0 < MR <= N.
!
!    Input, real*8 TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real*8 TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer  mr
  integer  n
  integer  nz_num

  real*8 a(nz_num)
  real*8 av
  real*8 c(1:mr)
  real*8, parameter :: delta = 1.0D-03
  real*8 g(1:mr+1)
  real*8 h(1:mr+1,1:mr)
  real*8 htmp
  integer  i
  integer  ia(nz_num)
  integer  itr
  integer  itr_max
  integer  itr_used
  integer  j
  integer  ja(nz_num)
  integer  k
  integer  k_copy
  real*8 mu
  real*8 r(1:n)
  real*8 rho
  real*8 rho_tol
  real*8 rhs(1:n)
  real*8 s(1:mr)
  real*8 tol_abs
  real*8 tol_rel
  real*8 v(1:n,1:mr+1)
  logical, parameter :: verbose = .true.
  real*8 x(1:n)
  real*8 y(1:mr+1)

  itr_used = 0

  if ( n < mr ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
    write ( *, '(a)' ) '  N < MR.'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  MR = ', mr
    stop
  end if

  do itr = 1, itr_max

    call ax_st ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    rho = sqrt ( dot_product ( r(1:n), r(1:n) ) )

    if ( verbose ) then
      write ( *, '(a,i8,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_st ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - h(j,k) * v(1:n,j)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( av + delta * h(k+1,k) == av ) then

        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do

        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then

        y(1:k+1) = h(1:k+1,k)

        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y(1:k+1) )
        end do

        h(1:k+1,k) = y(1:k+1)

      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )
      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g(1:k+1) )
      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i8,a,g14.6)' ) '  K =   ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)'       ) ' '
    write ( *, '(a)'       ) 'MGMRES_ST:'
    write ( *, '(a,i8)'    ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end
subroutine mult_givens ( c, s, k, g )

!*****************************************************************************80
!
!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
!
!  Discussion:
!
!    In order to make it easier to compare this code with the Original C,
!    the vector indexing is 0-based.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, real*8 C, S, the cosine and sine of a Givens
!    rotation.
!
!    Input, integer  K, indicates the location of the first
!    vector entry.
!
!    Input/output, real*8 G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
!
  implicit none

  integer  k

  real*8 c
  real*8 g(1:k+1)
  real*8 g1
  real*8 g2
  real*8 s

  g1 = c * g(k) - s * g(k+1)
  g2 = s * g(k) + c * g(k+1)

  g(k)   = g1
  g(k+1) = g2

  return
end


!! Solver disabled for compiling with GNU
!subroutine direct_solverPardiso ( n, nz_num, ia, ja, a, x, rhs)
!!      include 'mkl_pardiso.f77'
!!.. Internal solver memory pointer for 64-bit architectures
!!.. INTEGER*8 pt(64)
!!.. Internal solver memory pointer for 32-bit architectures
!!.. INTEGER*4 pt(64)
!!.. This is OK in both cases
!      INTEGER*8 pt(64)
!!.. All other variables
!      INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
!      INTEGER iparm(64)
!      INTEGER idum(1)
!      REAL*8  ddum(1)
!
!      integer  nz_num
!      integer  ia(n+1)
!      integer  ja(nz_num)
!      real*8 a(nz_num)
!      real*8 x(n)
!      real*8 rhs(n)
!      integer  i
!
!
!!.. Fill all arrays containing matrix data.
!!..
!!.. Set up PARDISO control parameter
!!..
!      do i = 1, 64
!         iparm(i) = 0
!      end do
!      iparm(1) = 1 ! no solver default
!      iparm(2) = 2 ! fill-in reordering from METIS
!      iparm(3) = 1 ! numbers of processors
!      iparm(4) = 0 ! no iterative-direct algorithm
!      iparm(5) = 0 ! no user fill-in reducing permutation
!      iparm(6) = 0 ! =0 solution on the first n compoments of x
!      iparm(7) = 0 ! not in use
!      iparm(8) = 9 ! numbers of iterative refinement steps
!      iparm(9) = 0 ! not in use
!      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
!      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
!      iparm(12) = 0 ! not in use
!      iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
!      iparm(14) = 0 ! Output: number of perturbed pivots
!      iparm(15) = 0 ! not in use
!      iparm(16) = 0 ! not in use
!      iparm(17) = 0 ! not in use
!      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
!      iparm(19) = -1 ! Output: Mflops for LU factorization
!      iparm(20) = 0 ! Output: Numbers of CG Iterations
!      error = 0 ! initialize error flag
!      msglvl = 1 ! print statistical information
!      mtype = 11 ! real unsymmetric
!!.. Initiliaze the internal solver memory pointer. This is only
!! necessary for the FIRST call of the PARDISO solver.
!      do i = 1, 64
!         pt(i) = 0
!      end do
!!.. Reordering and Symbolic Factorization, This step also allocates
!! all memory that is necessary for the factorization
!
!      phase = 11 ! only reordering and symbolic factorization
!      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,idum, nrhs, iparm, msglvl, ddum, ddum, error)
!      WRITE(*,*) 'Reordering completed ... '
!      IF (error .NE. 0) THEN
!         WRITE(*,*) 'The following ERROR was detected: ', error
!         STOP 1
!      END IF
!      WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
!      WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
!!.. Factorization.
!      phase = 22 ! only factorization
!      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,idum, nrhs, iparm, msglvl, ddum, ddum, error)
!      WRITE(*,*) 'Factorization completed ... '
!      IF (error .NE. 0) THEN
!         WRITE(*,*) 'The following ERROR was detected: ', error
!         STOP 1
!      ENDIF
!!.. Back substitution and iterative refinement
!      iparm(8) = 2 ! max numbers of iterative refinement steps
!      phase = 33 ! only factorization
!      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,idum, nrhs, iparm, msglvl, RHS, x, error)
!      WRITE(*,*) 'Solve completed ... '
!      WRITE(*,*) 'The solution of the system is '
!!.. Termination and release of memory
!      phase = -1 ! release internal memory
!      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,idum, nrhs, iparm, msglvl, ddum, ddum, error)
!return
!end subroutine




subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
  tol_abs, tol_rel,Prec_flag,verbose_flag)

!*****************************************************************************80
!
!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the linear system.
!
!    Input, integer  NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input/output, real*8 X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real*8 RHS(N), the right hand side of the linear system.
!
!    Input, integer  ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer  MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real*8 TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real*8 TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer  mr
  integer  n
  integer  nz_num
  integer  Prec_flag

  real*8 a(nz_num)
  real*8 av
  real*8 c(mr+1)
  real*8, parameter :: delta = 1.0D-03
  real*8 g(mr+1)
  real*8 h(mr+1,mr)
  real*8 htmp
  integer  i
  integer  ia(n+1)
  integer  itr
  integer  itr_max
  integer  itr_used
  integer  j
  integer  ja(nz_num)
  integer  k
  integer  k_copy
  real*8 mu
  real*8,ALLOCATABLE :: l(:)
  real*8 r(n)
  real*8 rho
  real*8 rho_tol
  real*8 rhs(n)
  real*8 s(mr+1)
  real*8 tol_abs
  real*8 tol_rel
  integer  ua(n)
  real*8 v(n,mr+1)
  logical verbose_flag,verbose,verbose0
  real*8 x(n)
  real*8 y(mr+1)
  verbose = .false.
  verbose0 = verbose_flag
  itr_used = 0

  if(Prec_flag .eq. 2) allocate(l(nz_num))


  call rearrange_cr ( n, nz_num, ia, ja, a )
  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
  if(Prec_flag .eq. 2) call ilu_cr ( n, nz_num, ia, ja, a, ua, l )
  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

  do itr = 1, itr_max

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    if(Prec_flag .eq. 1)then
    do i=1,n
         r(i)=r(i)/a(ua(i))
    enddo
    else if(Prec_flag .eq. 2) then
       call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )
    endif

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if
    !print*, 'rho=======', rho 
    if (rho .le. 1.0d-24) return

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

    if(Prec_flag .eq. 1)then
       do i=1,n
         v(i,k+1)=v(i,k+1)/a(ua(i))
      enddo
    else if(Prec_flag .eq. 2) then
      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )
    endif

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if (verbose0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  if(Prec_flag .eq. 2) deallocate(l)
  return
end

subroutine pmgmres_ilu2_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
  tol_abs, tol_rel,Prec_flag,ilucal_flag,verbose_flag)

!*****************************************************************************80
!
!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the linear system.
!
!    Input, integer  NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input/output, real*8 X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real*8 RHS(N), the right hand side of the linear system.
!
!    Input, integer  ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer  MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real*8 TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real*8 TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer  mr
  integer  n
  integer  nz_num
  integer  Prec_flag

  real*8 a(nz_num)
  real*8 av
  real*8 c(mr+1)
  real*8, parameter :: delta = 1.0D-03
  real*8 g(mr+1)
  real*8 h(mr+1,mr)
  real*8 htmp
  integer  i
  integer  ia(n+1)
  integer  itr
  integer  itr_max
  integer  itr_used
  integer  j
  integer  ja(nz_num)
  integer  k
  integer  k_copy
  real*8 mu
  real*8,save, ALLOCATABLE :: l(:)
  real*8 r(n)
  real*8 rho
  real*8 rho_tol
  real*8 rhs(n)
  real*8 s(mr+1)
  real*8 tol_abs
  real*8 tol_rel
  integer  ua(n)
  real*8 v(n,mr+1)
  logical verbose_flag,verbose,verbose0,ilucal_flag
  real*8 x(n)
  real*8 y(mr+1)
  verbose = .false.
  verbose0 = verbose_flag
  itr_used = 0

  if(Prec_flag .eq. 2  .and. ilucal_flag) then
     if(allocated(l))  then  
     deallocate(l)
     endif
     allocate(l(nz_num))
     l=0.0
  endif


  call rearrange_cr ( n, nz_num, ia, ja, a )
  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
  if((Prec_flag .eq. 2) .and. ilucal_flag) call ilu_cr ( n, nz_num, ia, ja, a, ua, l )
  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

  do itr = 1, itr_max

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    if(Prec_flag .eq. 1)then
    do i=1,n
         r(i)=r(i)/a(ua(i))
    enddo
    else if(Prec_flag .eq. 2) then
       call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )
    endif

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if
    !print*, 'rho=======', rho 
    if (rho .le. 1.0d-24) return

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

    if(Prec_flag .eq. 1)then
       do i=1,n
         v(i,k+1)=v(i,k+1)/a(ua(i))
      enddo
    else if(Prec_flag .eq. 2) then
      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )
    endif

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if (verbose0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

!  if(Prec_flag .eq. 2) deallocate(l)
  return
end

subroutine pmgmres_ilu2_crRowScale ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
  tol_abs, tol_rel,Prec_flag,ilucal_flag,verbose_flag)

!*****************************************************************************80
!
!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the linear system.
!
!    Input, integer  NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input/output, real*8 X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real*8 RHS(N), the right hand side of the linear system.
!
!    Input, integer  ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer  MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real*8 TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real*8 TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer  mr
  integer  n
  integer  nz_num
  integer  Prec_flag

  real*8 a(nz_num)
  real*8 av
  real*8 c(mr+1)
  real*8, parameter :: delta = 1.0D-03
  real*8 g(mr+1)
  real*8 h(mr+1,mr)
  real*8 htmp
  integer  i
  integer  ia(n+1)
  integer  itr
  integer  itr_max
  integer  itr_used
  integer  j
  integer  ja(nz_num)
  integer  k
  integer  k_copy
  real*8 mu
  real*8,save, ALLOCATABLE :: l(:)
  real*8 r(n),diag(n)
  real*8 rho
  real*8 rho_tol
  real*8 rhs(n)
  real*8 s(mr+1)
  real*8 tol_abs

  real*8 tol_rel
  integer  ua(n)
  real*8 v(n,mr+1)
  logical verbose_flag,verbose,verbose0,ilucal_flag
  real*8 x(n)
  real*8 y(mr+1)
  integer ierr,k1,k2,ii,nrm
  real*8 scal
  verbose = .false.
  verbose0 = verbose_flag
  itr_used = 0

  if(Prec_flag .eq. 2  .and. ilucal_flag) then
     if(allocated(l))  then  
     deallocate(l)
     endif
     allocate(l(nz_num))
     l=0.0
  endif


  call rearrange_cr ( n, nz_num, ia, ja, a )

!     
!     normalize each row 
!     

  nrm=2
  do ii = 1, n

    scal = 0.0D+00
    k1 = ia(ii)
    k2 = ia(ii+1) - 1

    if ( nrm .eq. 0 ) then
      do k = k1, k2
        scal = max ( scal, abs ( a(k) ) )
      end do
    else if ( nrm .eq. 1 ) then
      do k = k1, k2
        scal = scal + abs ( a(k) )
      end do
    else
      do k = k1, k2
        scal = scal + a(k)**2
      end do
    end if

    if ( nrm .eq. 2 ) then
      scal = sqrt ( scal )
    end if

    diag(ii) = scal

  end do

  ierr = 0
  do  j=1, n
         if (diag(j) .eq. 0.0d0) then
            ierr = j 
            
            stop "diag(j) =0"
         else
            diag(j) = 1.0d0/diag(j)
         endif
  enddo
  do ii=1,n
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         scal = diag(ii) 
         do  k=k1, k2
           a(k) = a(k)*scal
         enddo
         rhs(ii)=rhs(ii)*scal
  enddo    


  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
  if((Prec_flag .eq. 2) .and. ilucal_flag) call ilu_cr ( n, nz_num, ia, ja, a, ua, l )
  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

  do itr = 1, itr_max

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    if(Prec_flag .eq. 1)then
    do i=1,n
         r(i)=r(i)/a(ua(i))
    enddo
    else if(Prec_flag .eq. 2) then
       call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )
    endif

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if
    !print*, 'rho=======', rho 
    if (rho .le. 1.0d-24) return

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

    if(Prec_flag .eq. 1)then
       do i=1,n
         v(i,k+1)=v(i,k+1)/a(ua(i))
      enddo
    else if(Prec_flag .eq. 2) then
      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )
    endif

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if (verbose0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

!  if(Prec_flag .eq. 2) deallocate(l)
  return
end

subroutine pmgmres_ilu_crMRHS ( n, nz_num, ia, ja, a, &
 x1,x2,x3, rhs1,rhs2,rhs3, &
 itr_max, mr, &
  tol_abs, tol_rel,Prec_flag,verbose_flag)

!*****************************************************************************80
!
!! PMGMRES_ILU_CRMRHS applies the preconditioned restarted GMRES algorithm for multiple right hand sides.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the linear system.
!
!    Input, integer  NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input/output, real*8 X1(N),X2(N),X3(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real*8 RHS1(N),RHS2(N),RHS3(N) the right hand side of the linear system.
!
!    Input, integer  ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer  MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real*8 TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real*8 TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer  mr
  integer  n
  integer  nz_num
  integer  Prec_flag

  real*8 a(nz_num)
  real*8 av
  real*8 c(mr+1)
  real*8, parameter :: delta = 1.0D-03
  real*8 g(mr+1)
  real*8 h(mr+1,mr)
  real*8 htmp
  integer  i
  integer  ia(n+1)
  integer  itr
  integer  itr_max
  integer  itr_used
  integer  itr_usedwrite(3)
  integer  j
  integer  ja(nz_num)
  integer  k
  integer  k_copy
  integer ivec
  real*8 mu
  real*8,ALLOCATABLE :: l(:)
  real*8 r(n)
  real*8 rho
  real*8 rho_used(3)
  real*8 rho_tol
  real*8 rhs1(n),rhs2(n),rhs3(n),rhs(n)
  real*8 s(mr+1)
  real*8 tol_abs
  real*8 tol_rel
  integer  ua(n)
  real*8 v(n,mr+1)
  logical verbose_flag,verbose,verbose0
  real*8 x1(n),x2(n),x3(n),x(n)
  real*8 y(mr+1)
  verbose = .false.
  verbose0 = verbose_flag
  itr_used = 0

  if(Prec_flag .eq. 2) allocate(l(nz_num))


  call rearrange_cr ( n, nz_num, ia, ja, a )
  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
  if(Prec_flag .eq. 2) call ilu_cr ( n, nz_num, ia, ja, a, ua, l )
  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

 do ivec=1,3
  itr_used = 0
  if (ivec==1) then
   x(1:n)  =x1(1:n)
   rhs(1:n)=rhs1(1:n) 
  elseif(ivec==2) then
   x(1:n)  =x2(1:n)
   rhs(1:n)=rhs2(1:n) 
  else
   x(1:n)  =x3(1:n)
   rhs(1:n)=rhs3(1:n) 
  endif
  do itr = 1, itr_max

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    if(Prec_flag .eq. 1)then
    do i=1,n
         r(i)=r(i)/a(ua(i))
    enddo
    else if(Prec_flag .eq. 2) then
       call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )
    endif

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,i4,a,g14.6)' ) '  VEC = ', ivec,'  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if
    !print*, 'rho=======', rho 
    if (rho .le. 1.0d-24) then 
       exit
    endif

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

    if(Prec_flag .eq. 1)then
       do i=1,n
         v(i,k+1)=v(i,k+1)/a(ua(i))
      enddo
    else if(Prec_flag .eq. 2) then
      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )
    endif

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,i4,a,g14.6)' ) '  VEC = ', ivec,'  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
       exit
    end if

  end do

  itr_usedwrite(ivec)=itr_used 
  rho_used(ivec)=rho

  if (ivec==1) then
   x1(1:n)  =x(1:n)
  elseif(ivec==2) then
   x2(1:n)  =x(1:n)
  else
   x3(1:n)  =x(1:n)
  endif
  enddo

  if (verbose0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CRMRHS:'
    write ( *, '(a,i4,a,i4,a,i4)' ) '  Iterations = ', itr_usedwrite(1),' ', itr_usedwrite(2),' ', itr_usedwrite(3)
    write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) '  Final residual = ', rho_used(1),' ', rho_used(2),' ', rho_used(3)

  end if

  if(Prec_flag .eq. 2) deallocate(l)
  return
end

subroutine pmgmres_ilu2_crBigRHS ( n, nz_num, ia, ja, a, &
 x1,x2,x3, rhs1,rhs2,rhs3, &
 itr_max, mr, &
 tol_abs, tol_rel,Prec_flag,ilucal_flag,verbose_flag)

!*****************************************************************************80
!
!! PMGMRES_ILU2_CRBigRHS applies the preconditioned restarted GMRES algorithm 
!                       for a combination of the right hand side.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the linear system.
!
!    Input, integer  NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input/output, real*8 X1(N),X2(N),X3(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real*8 RHS1(N),RHS2(N),RHS3(N) the right hand side of the linear system.
!
!    Input, integer  ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer  MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real*8 TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real*8 TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer  mr
  integer  n, n3
  integer  nz_num
  integer  Prec_flag

  real*8 a(nz_num)
  real*8 av
  real*8 c(mr+1)
  real*8, parameter :: delta = 1.0D-03
  real*8 g(mr+1)
  real*8 h(mr+1,mr)
  real*8 htmp
  integer  i
  integer  ia(3*n+1)
  integer  itr
  integer  itr_max
  integer  itr_used
  integer  itr_usedwrite(3)
  integer  j
  integer  ja(nz_num)
  integer  k
  integer  k_copy
  integer ivec
  real*8 mu
  real*8,save, ALLOCATABLE :: l(:)
  real*8 r(3*n)
  real*8 rho
  real*8 rho_used(3)
  real*8 rho_tol
  real*8 rhs1(n),rhs2(n),rhs3(n),rhs(3*n)
  real*8 s(mr+1)
  real*8 tol_abs
  real*8 tol_rel
  integer  ua(3*n)
  real*8 v(3*n,mr+1)
  logical verbose_flag,verbose,verbose0,ilucal_flag
  real*8 x1(n),x2(n),x3(n),x(3*n)
  real*8 y(mr+1)
  verbose = .false.
  verbose0 = verbose_flag
  itr_used = 0
  n3=3*n
  if(Prec_flag .eq. 2  .and. ilucal_flag) then
     if(allocated(l))  then  
     deallocate(l)
     endif
     allocate(l(nz_num))
     l=0.0
  endif

  call rearrange_cr ( n3, nz_num, ia, ja, a )
  call diagonal_pointer_cr ( n3, nz_num, ia, ja, ua )
  if((Prec_flag .eq. 2) .and. ilucal_flag) call ilu_cr ( n3, nz_num, ia, ja, a, ua, l )
  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if


   itr_used = 0
     x(1:n)      =x1(1:n)
   rhs(1:n)      =rhs1(1:n) 
     x(n+1:2*n)  =x2(1:n)
   rhs(n+1:2*n)  =rhs2(1:n) 
     x(2*n+1:3*n)=x3(1:n)
   rhs(2*n+1:3*n)=rhs3(1:n) 

  ivec=1

  do itr = 1, itr_max

    call ax_cr ( n3, nz_num, ia, ja, a, x, r )

    r(1:n3) = rhs(1:n3) - r(1:n3)

    if(Prec_flag .eq. 1)then
    do i=1,n3
         r(i)=r(i)/a(ua(i))
    enddo
    else if(Prec_flag .eq. 2) then
       call lus_cr ( n3, nz_num, ia, ja, l, ua, r, r )
    endif

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,i4,a,g14.6)' ) '  VEC = ', ivec,'  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if
    !print*, 'rho=======', rho 
    if (rho .le. 1.0d-24) then 
       exit
    endif

    v(1:n3,1) = r(1:n3) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n3, nz_num, ia, ja, a, v(1:n3,k), v(1:n3,k+1) ) 

    if(Prec_flag .eq. 1)then
       do i=1,n3
         v(i,k+1)=v(i,k+1)/a(ua(i))
      enddo
    else if(Prec_flag .eq. 2) then
      call lus_cr ( n3, nz_num, ia, ja, l, ua, v(1:n3,k+1), v(1:n3,k+1))
    endif

      av = sqrt ( dot_product ( v(1:n3,k+1), v(1:n3,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n3,k+1), v(1:n3,j) )
        v(1:n3,k+1) = v(1:n3,k+1) - v(1:n3,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n3,k+1), v(1:n3,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n3,k+1), v(1:n3,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n3,k+1) = v(1:n3,k+1) - htmp * v(1:n3,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n3,k+1), v(1:n3,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n3,k+1) = v(1:n3,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,i4,a,g14.6)' ) '  VEC = ', ivec,'  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n3
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
       exit
    end if

  end do

  itr_usedwrite(ivec)=itr_used 
  rho_used(ivec)=rho

   x1(1:n)  =x(1:n)
   x2(1:n)  =x(n+1:2*n)
   x3(1:n)  =x(2*n+1:3*n)

  if (verbose0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CRBigRHS:'
    write ( *, '(a,i4,a,i4,a,i4)' ) '  Iterations = ', itr_usedwrite(1)
    write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) '  Final residual = ', rho_used(1)

  end if

!  if(Prec_flag .eq. 2) deallocate(l)
  return
end


subroutine pmgmres_ilu_crBigRHS ( n, nz_num, ia, ja, a, &
 x1,x2,x3, rhs1,rhs2,rhs3, &
 itr_max, mr, &
 tol_abs, tol_rel,Prec_flag,verbose_flag)

!*****************************************************************************80
!
!! PMGMRES_ILU_CRBigRHS applies the preconditioned restarted GMRES algorithm 
!                       for a combination of the right hand side.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the linear system.
!
!    Input, integer  NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real*8 A(NZ_NUM), the matrix values.
!
!    Input/output, real*8 X1(N),X2(N),X3(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real*8 RHS1(N),RHS2(N),RHS3(N) the right hand side of the linear system.
!
!    Input, integer  ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer  MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real*8 TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real*8 TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer  mr
  integer  n, n3
  integer  nz_num
  integer  Prec_flag

  real*8 a(nz_num)
  real*8 av
  real*8 c(mr+1)
  real*8, parameter :: delta = 1.0D-03
  real*8 g(mr+1)
  real*8 h(mr+1,mr)
  real*8 htmp
  integer  i
  integer  ia(3*n+1)
  integer  itr
  integer  itr_max
  integer  itr_used
  integer  itr_usedwrite(3)
  integer  j
  integer  ja(nz_num)
  integer  k
  integer  k_copy
  integer ivec
  real*8 mu
  real*8,ALLOCATABLE :: l(:)
  real*8 r(3*n)
  real*8 rho
  real*8 rho_used(3)
  real*8 rho_tol
  real*8 rhs1(n),rhs2(n),rhs3(n),rhs(3*n)
  real*8 s(mr+1)
  real*8 tol_abs
  real*8 tol_rel
  integer  ua(3*n)
  real*8 v(3*n,mr+1)
  logical verbose_flag,verbose,verbose0
  real*8 x1(n),x2(n),x3(n),x(3*n)
  real*8 y(mr+1)
  verbose = .false.
  verbose0 = verbose_flag
  itr_used = 0
  n3=3*n
  if(Prec_flag .eq. 2) allocate(l(nz_num))


  call rearrange_cr ( n3, nz_num, ia, ja, a )
  call diagonal_pointer_cr ( n3, nz_num, ia, ja, ua )
  if(Prec_flag .eq. 2) call ilu_cr ( n3, nz_num, ia, ja, a, ua, l )
  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if


   itr_used = 0
     x(1:n)      =x1(1:n)
   rhs(1:n)      =rhs1(1:n) 
     x(n+1:2*n)  =x2(1:n)
   rhs(n+1:2*n)  =rhs2(1:n) 
     x(2*n+1:3*n)=x3(1:n)
   rhs(2*n+1:3*n)=rhs3(1:n) 

  ivec=1

  do itr = 1, itr_max

    call ax_cr ( n3, nz_num, ia, ja, a, x, r )

    r(1:n3) = rhs(1:n3) - r(1:n3)

    if(Prec_flag .eq. 1)then
    do i=1,n3
         r(i)=r(i)/a(ua(i))
    enddo
    else if(Prec_flag .eq. 2) then
       call lus_cr ( n3, nz_num, ia, ja, l, ua, r, r )
    endif

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,i4,a,g14.6)' ) '  VEC = ', ivec,'  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if
    !print*, 'rho=======', rho 
    if (rho .le. 1.0d-24) then 
       exit
    endif

    v(1:n3,1) = r(1:n3) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n3, nz_num, ia, ja, a, v(1:n3,k), v(1:n3,k+1) ) 

    if(Prec_flag .eq. 1)then
       do i=1,n3
         v(i,k+1)=v(i,k+1)/a(ua(i))
      enddo
    else if(Prec_flag .eq. 2) then
      call lus_cr ( n3, nz_num, ia, ja, l, ua, v(1:n3,k+1), v(1:n3,k+1))
    endif

      av = sqrt ( dot_product ( v(1:n3,k+1), v(1:n3,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n3,k+1), v(1:n3,j) )
        v(1:n3,k+1) = v(1:n3,k+1) - v(1:n3,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n3,k+1), v(1:n3,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n3,k+1), v(1:n3,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n3,k+1) = v(1:n3,k+1) - htmp * v(1:n3,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n3,k+1), v(1:n3,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n3,k+1) = v(1:n3,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,i4,a,g14.6)' ) '  VEC = ', ivec,'  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n3
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
       exit
    end if

  end do

  itr_usedwrite(ivec)=itr_used 
  rho_used(ivec)=rho

   x1(1:n)  =x(1:n)
   x2(1:n)  =x(n+1:2*n)
   x3(1:n)  =x(2*n+1:3*n)

  if (verbose0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CRBigRHS:'
    write ( *, '(a,i4,a,i4,a,i4)' ) '  Iterations = ', itr_usedwrite(1)
    write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) '  Final residual = ', rho_used(1)

  end if

  if(Prec_flag .eq. 2) deallocate(l)
  return
end



subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real*8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer  N, the number of entries in the vector.
!
!    Input/output, integer  SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real*8 R(N), the vector of pseudorandom values.
!
  implicit none

  integer  n

  integer  i
  integer  k
  integer  seed
  real*8 r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine rearrange_cr ( n, nz_num, ia, ja, a )

!*****************************************************************************80
!
!! REARRANGE_CR sorts a sparse compressed row matrix.
!
!  Discussion:
!
!    This routine guarantees that the entries in the CR matrix
!    are properly sorted.
!
!    After the sorting, the entries of the matrix are rearranged in such
!    a way that the entries of each column are listed in ascending order
!    of their column values.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), the compressed row indices.
!
!    Input/output, integer  JA(NZ_NUM), the column indices.
!    On output, these may have been rearranged by the sorting.
!
!    Input/output, real*8 A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
  implicit none

  integer  n
  integer  nz_num

  real*8 a(nz_num)
  integer  i
  integer  ia(n+1)
  integer  i4temp
  integer  ja(nz_num)
  integer  k
  integer  l
  real*8 r8temp

  do i = 1, n

    do k = ia(i), ia(i+1) - 2
      do l = k + 1, ia(i+1) - 1

        if ( ja(l) < ja(k) ) then
          i4temp = ja(l)
          ja(l)  = ja(k)
          ja(k)  = i4temp

          r8temp = a(l)
          a(l)   = a(k)
          a(k)   = r8temp
        end if

      end do
    end do

  end do

  return
end

! SUSSMAN: modified the name to avoid conflict
subroutine timestamp_Burkardt( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer  d
  integer  h
  integer  m
  integer  mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer  n
  integer  s
  integer  values(8)
  integer  y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
