MODULE Inverters

	use TreeCodeGlobal
	implicit none

	! used to restrict variables to prevent accidental calling outside
	private      ! everything defaults to 'private' except those expressly labeled otherwise.
	public			:: LEGS, GAUSS, gaussianElimination, LinSolv4	! subroutines
	real(KIND=r8), parameter :: DEFAULT_SMALLEST_PIVOT=1e-10

contains 

!**GAUSS*****************************************************************
! Subroutine to find solution of a linear system of N equations in N   *
! unknowns using Gaussian elimination, provided a unique solution      *
! exists.  The coefficients and constants of the linear system are     *
! stored in the matrix LIN, which has N rows and N+1 columns.			*
! If the system is singular, SINGUL is returned as true, and the       *
! solution X is undefined.  Local identifiers used are:                *                        
!     I,J,K  : subscripts                                              *
!     MULT   : multiplier used to eliminate an unknown                 *
!     ABSPIV : absolute value of pivot element                         *
!     PIVROW : row containing pivot element                            *
!     EPSIL  : a small positive real value ("almost zero")             *
!     TEMP   : used to interchange rows of matrix                      *
!
! Accepts: Two-dimensional array LIM, integers LIMROW				   *
! Returns: One-dimensional array X and logical value SINGUL            *
!************************************************************************

	SUBROUTINE GAUSS(A, LIMROW, X, B)

		! incoming variables
		integer, intent(in)								:: LIMROW
		real(KIND=r8), dimension(LIMROW, LIMROW), intent(in) :: A
		real(KIND=r8), dimension(LIMROW), intent(in)	:: B

		! outsourced variables
		real(KIND=r8), dimension(LIMROW),intent(inout) 	:: X

		integer											:: PIVROW, I, J, K
		real(KIND=r8), parameter						:: EPSIL=1e-15
		real(KIND=r8), dimension(LIMROW, LIMROW+1)		:: LIN
		real(KIND=r8)									:: TEMP, MULT, ABSPIV

		DO I = 1, LIMROW
			DO J = 1, LIMROW
				LIN(I,J) = A(I,J)
			END DO
			LIN(I,LIMROW+1) = B(I)
		END DO

		ILoop: DO I = 1, LIMROW ! Locate pivot element
			ABSPIV = abs(LIN(I,I))
			PIVROW = I
			DO K = I + 1, LIMROW
				IF (ABS(LIN(K,I)) > ABSPIV) THEN
					ABSPIV = ABS(LIN(K,I))
					PIVROW = K
				ENDIF
			END DO

	! Check if matrix is (nearly) singular
			IF (ABSPIV < EPSIL) THEN
				stop 'nearly singular matrix. Stopping'
			ENDIF

	! It isn't, so interchange rows PIVROW and I if necessary
			IF (PIVROW /= I) THEN
				DO J = 1, LIMROW + 1
					TEMP = LIN(I,J)
					LIN(I,J) = LIN(PIVROW,J)
					LIN(PIVROW,J) = TEMP
				END DO
			ENDIF

	! Eliminate Ith unknown from equations I + 1, ..., N
			JLoop: DO J = I + 1, LIMROW
				MULT = -LIN(J,I) / LIN(I,I)
				KLoop: DO K = I, LIMROW + 1
					LIN(J,K) = LIN(J,K) +  MULT * LIN(I,K)
				END DO KLoop
			END DO JLoop
		END DO ILoop

! Find the solutions by back substitution
		X(LIMROW) = LIN(LIMROW, LIMROW + 1) / LIN(LIMROW,LIMROW)
		DO J = LIMROW - 1, 1, -1
			X(J) = LIN(J, LIMROW + 1)
			DO K = J + 1, LIMROW
				X(J) = X(J) - LIN(J,K) * X(K)
			END DO
			X(J) = X(J) / LIN(J,J)
		END DO

	END SUBROUTINE GAUSS

! ***************
! SUBROUTINE GMRes
! ***************
! DESCRIPTION: an O(N^2 log N) Matrix solver process for iteratively solving Ax=b
!
! INPUT: A,b
!
! OUTPUT: x
!
! CALLING PROGRAM: MAT_Greens
!
! By:           Jerry Emhoff
! Rewriten By:  Andrew J. Christlieb
! Date:         8/26/03
! Converted to Fortran 90/95 By: Anton VanderWyst
!
!Contact: Andrew J. Christlieb, Ph. D.
!         Assistant Professor
!         Department of Mathematics
!         University of Michigan
!         2470 East Hall
!         Ann Arbor, MI 48109.
!         
!         Office: 4851 East Hall
!         E-mail: christli@umich.edu
!         Tel:    (734) 763-5725
!------------------------------------------------------
! ** to be worked on later

! Updated 10/24/2001.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.3   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!
! http://www.physics.unlv.edu/~pang/comp3/code43.f90
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! An example of solving linear equation set A(N,N)*X(N) = B(N)
! with the partial-pivoting Gaussian elimination scheme. 

	SUBROUTINE LEGS (A,N,B,X,INDX)
!
! Subroutine to solve the equation A(N,N)*X(N) = B(N) with the
! partial-pivoting Gaussian elimination scheme.
! Copyright (c) Tao Pang 2001.
!
		IMPLICIT NONE
		INTEGER, INTENT (IN) :: N
		INTEGER :: I,J
		INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
		REAL*8, INTENT (INOUT), DIMENSION (N,N) :: A
		REAL*8, INTENT (INOUT), DIMENSION (N) :: B
		REAL*8, INTENT (INOUT), DIMENSION (N) :: X
!
		CALL ELGS (A,N,INDX)

		DO I = 1, N-1
			DO J = I+1, N
				B(INDX(J)) = B(INDX(J))-A(INDX(J),I)*B(INDX(I))
			END DO
		END DO
!
		X(N) = B(INDX(N))/A(INDX(N),N)
		DO I = N-1, 1, -1
			X(I) = B(INDX(I))
			DO J = I+1, N
				X(I) = X(I)-A(INDX(I),J)*X(J)
			END DO
			X(I) =  X(I)/A(INDX(I),I)
		END DO
!
	END SUBROUTINE LEGS
!
	SUBROUTINE ELGS (A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
!
		IMPLICIT NONE
		INTEGER, INTENT (IN) :: N
		INTEGER :: I,J,K,ITMP
		INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
		REAL*8 :: C1,PI,PI1,PJ
		REAL*8, INTENT (INOUT), DIMENSION (N,N) :: A
		REAL*8, DIMENSION (N) :: C
!
! Initialize the index
!
		DO I = 1, N
			INDX(I) = I
		END DO
!
! Find the rescaling factors, one from each row
!
		DO I = 1, N
			C1= 0.0
			DO J = 1, N
				C1 = DMAX1(C1,ABS(A(I,J)))
			END DO
			C(I) = C1
		END DO
!
! Search the pivoting (largest) element from each column
!
		DO J = 1, N-1
			PI1 = 0.0
			DO I = J, N
				PI = ABS(A(INDX(I),J))/C(INDX(I))
				IF (PI.GT.PI1) THEN
					PI1 = PI
					K   = I
				ENDIF
			END DO
!
! Interchange the rows via INDX(N) to record pivoting order
!
			ITMP    = INDX(J)
			INDX(J) = INDX(K)
			INDX(K) = ITMP
			DO I = J+1, N
				PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
				A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!

				DO K = J+1, N
					A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
				END DO
			END DO
		END DO
!
	END SUBROUTINE ELGS

   ! Use Gaussian elimination to calculate the solution to the linear 
   ! system, A x = b.  No partial pivoting is done.  If the threshold 
   ! argument is present, it is used as the smallest allowable pivot 
   ! encountered in the computation; otherwise, DEFAULT_SMALLEST_PIVOT, 
   ! defined in this module, is used as the default threshold.  The status
   ! of the computation is a logical returned by the function indicating
   ! the existence of a unique solution (.true.), or the nonexistence of
   ! a unique solution or threshold passed (.false.).

   ! Note that this is an inappropriate method for some linear systems.
   ! In particular, the linear system, M x = b, where M = 10e-12 I, will 
   ! cause this routine to fail due to the presence of small pivots.  
   ! However, this system is perfectly conditioned, with solution x = b.
   ! http://www.scd.ucar.edu/tcg/consweb/Fortran90/F90Tutorial/node26.html

	SUBROUTINE gaussianElimination( A, b, x, YesSolve, threshold )
		implicit none
		logical, intent(out) 							:: YesSolve
		real(KIND=r8), dimension( :, : ), intent( in ) 	:: A   	! Assume the shape of A.
		real(KIND=r8), dimension( : ), intent( in ) 	:: b    ! Assume the shape of b.
		real(KIND=r8), dimension( : ), intent( out ) 	:: x    ! Assume the shape of x.

      ! The optional attribute specifies that the indicated argument
      ! is not required to be present in a call to the function.  The
      ! presence of optional arguments, such as threshold, may be checked
      ! using the intrinsic logical function, present (see below).

		real(KIND=r8), optional, intent( in ) :: threshold

		integer			:: i, j   ! Local index variables.
		integer 		:: N      ! Order of the linear system.
		real(KIND=r8) 	:: m         ! Multiplier.
		real(KIND=r8) 	:: smallestPivot = DEFAULT_SMALLEST_PIVOT

      ! Pointers to the appropriate rows of the matrix during the elmination.
		real(KIND=r8), dimension( : ), pointer :: pivotRow
		real(KIND=r8), dimension( : ), pointer :: currentRow

      ! Copies of the input arguments.  These copies are modified during
      ! the computation.
      ! The target attribute is used to indicate that the specified 
      ! variable may be the target of a pointer.  Rows of ACopy are targets
      ! of pivotRow and currentRow, defined above.

		real(KIND=r8), dimension( size( b ), size( b ) ), target :: ACopy
		real(KIND=r8), dimension( size( b ) ) :: bCopy

      ! Status of the computation.  The return value of the function.
		logical successful

      ! Change the smallestPivot if the threshold argument was included.
		if ( present( threshold ) ) smallestPivot = abs( threshold )

      ! Setup the order of the system by using the intrinsic function size.
      ! size returns the number of elements in the specified dimension of
      ! an array or the total number of elements if the dimension is not
      ! specified.  Also assume that a unique solution exists initially.

		N = size( b )   
		ACopy = A
		bCopy = b
		successful = .true.

      ! Begin the Gaussian elimination algorithm.
      ! Note the use of array sections in the following loops.  These 
      ! eliminate the need for many do loops that are common in Fortran 
      ! 77 code.
      ! Pointers are also used below and enhance the readability of the
      ! elimination process.

      ! Begin with the first row.
		i = 1

      ! Reduce the system to upper triangular.
		do while ( ( successful ) .and. ( i <= N-1 ) )

         ! The following statement is called pointer assignment and uses
         ! the pointer assignment operator `=>'.  This causes pivotRow 
         ! to be an alias for the ith row of ACopy.  Note that this does
         ! not cause any movement of data.

         ! Assign the pivot row.
			pivotRow => ACopy( i, : )

			if (1==i) then
				print *, 'testing ACopy'
			endif
         ! Verify that the current pivot is not smaller than smallestPivot.
			successful = abs( pivotRow( i ) ) >= smallestPivot

			if ( successful ) then

	    ! Eliminate the entries in the pivot column below the pivot row.

				do j = i+1, N
	       ! Assign the current row.
					currentRow => ACopy( j, : )

               ! Calculate the multiplier.
					m = currentRow( i ) / pivotRow( i ) 

               ! Perform the elimination step on currentRow and right 
               ! hand side, bCopy.
					currentRow = m * pivotRow - currentRow
					bCopy( j ) = m * bCopy( i ) - bCopy( j )
				end do

			else
				print *, i, pivotRow( i )
			end if

         ! Move to the next row.
			i = i + 1

		end do

      ! Check the last pivot.
		pivotRow => ACopy( N, : )
		if ( successful ) successful = abs( pivotRow( N ) ) >= smallestPivot

		if ( successful ) then
			do i = N, 2, -1   ! Backward substitution.

            ! Determine the ith unknown, x( i ).
				x( i ) = bCopy( i ) / ACopy( i, i )

            ! Substitute the now known value of x( i ), reducing the order of 
            ! the system by 1.
				bCopy = bCopy - x( i ) * ACopy( :, i )

			end do
		end if

      ! Determine the value of x( 1 ) as a special case.
		if ( successful ) x( 1 ) = bCopy( 1 ) / ACopy( 1, 1 )

      ! Prepare the return value of the function.
		YesSolve = successful

	end SUBROUTINE gaussianElimination


   ! The LU decomposition of a matrix may be represented in a compact form
   ! existing in a single matrix, M,  if the assignments M=L and M=U are 
   ! done (in that order).  The diagonal entries in L are assumed to be 
   ! unity so that no storage space is necessary.  Instead, the diagonal
   ! of M is used to hold the diagonal entries of U.  This is a common 
   ! method of storing the LU decomposition of a matrix.

   ! The algorithm belows makes an additional assumption concerning the 
   ! pivots or diagonal elements of U.  Computation terminates if one of
   ! these pivots is smaller than the given or default threshold.  In this
   ! case, the LU decomposition is not formed.  Note that this algorithm 
   ! successfully terminates if such an LU can be computed.  In this case
   ! the coefficient matrix, A, is nonsingular.  (No attempt for recovery,
   ! such as permutation of rows, is done.)

   ! Compute the LU decomposition of A, storing the result in LU so that
   ! A is not overwritten.  If the threshold argument is present, it is used 
   ! as the smallest allowable pivot encountered in the computation;
   ! otherwise, DEFAULT_SMALLEST_PIVOT, defined in this module, is used as
   ! the default threshold during the computation.  The status of the
   ! computation is a logical returned by the function indicating the
   ! success (.true.) or failure (.false.) of the factorization
   ! After the computation, LU will contain the multipliers below the main
   ! diagonal (L) and the result after elimination on and above the main
   ! diagonal (U), so that A = L * U.

	function LUFactor ( A, LU, threshold ) 
		implicit none
		logical LUFactor
		real(KIND=r8), dimension( :, : ), intent( in ) :: A
		real(KIND=r8), dimension( :, : ), intent( out ) :: LU 
		real(KIND=r8), optional, intent( in ) :: threshold

		integer k, i
		integer N
		logical successful   ! Status of the computation.
		real(KIND=r8) :: smallestPivot = DEFAULT_SMALLEST_PIVOT

      ! Reassign the smallestPivot, set the order of the system, and 
      ! copy A into LU as it will be written to during the factorization.

		if ( present( threshold ) ) smallestPivot = abs( threshold )
		N = size( A, 1 )   
		LU = A

      ! Begin the LU factorization algorithm.
      ! The status of the computation is initially successful.
		successful = .true.

		k = 1   ! Begin with the first column.
		do while ( ( successful ) .and. ( k <= N-1 ) )

         ! Verify that the kth pivot is not smaller than smallestPivot.
			successful = abs( LU( k, k ) ) >= smallestPivot

			if ( successful ) then
            ! Calculate the multipliers (L) for the current column.
				LU( k+1:N, k ) = LU( k+1:N, k ) / LU( k, k )

            ! Perform elimination on the upper portion of the matrix (U). 
				do i = k+1, N
					LU( i, k+1:N ) = LU( i, k+1:N ) - LU( i, k ) * LU( k, k+1:N )
				enddo

				k = k + 1   ! Move to the next column.
			end if

		enddo

      ! Prepare the return value of the function.
		LUFactor = successful

	end function LUFactor


   ! Let A = L*U where LU represents the LU decomposition of A stored in the
   ! format produced by LUFactor, A, L, U in R**(NxN).
   ! Solve the linear system, A x = b, using the LU decomposition of A stored
   ! in LU.  Since LU is the LU decomposition of A, A is nonsingular.  
   ! Consequently, the columns of A constitute a basis for R**N.   So, there 
   ! must exist a unique solution to the linear system A x = b.
   ! LUSolve returns the solution to this linear system.

	function LUSolve( LU, b ) result( x )
		implicit none
		real(KIND=r8), dimension( :, : ), intent( in ) :: LU
		real(KIND=r8), dimension( : ), intent( in ) :: b
		real(KIND=r8), dimension( size( b ) ) :: x

		integer k
		integer N
		real(KIND=r8), dimension( size( b ) ) :: bCopy

      ! Determine the order of the system and store a copy of b in bCopy
      ! as it is written during the computation.
		N = size( b )
		bCopy = b

      ! Assume LU is in the form of LU and solve the system in two steps.
      ! First, using forward elmination to solve L y = b, then using
      ! backward elmination to solve U x = y.  In both cases, the right
      ! hand side is overwritten with the solution as it is computed.

      ! Forward elimination.  Store the solution into the right hand side.
		do k = 1, N-1
			bCopy( k+1:N ) = bCopy( k+1:N ) - bCopy( k ) * LU( k+1:N, k )
		end do

      ! Backward elimination.  Store the solution into the right hand side.
		do k = N, 2, -1
			bCopy( k ) = bcopy( k ) / LU( k, k )
			bCopy( 1:k-1 ) = bCopy( 1:k-1 ) - bCopy( k ) * LU( 1:k-1, k )
		end do

      ! Solve for the 1st unknown as a special case.
		bCopy( 1 ) = bCopy( 1 ) / LU( 1, 1 )

      ! Assign a return value for the function via its result variable, x.
		x = bCopy

	end function LUSolve


   ! Output A in Matlab format, using name in the Matlab assignment statement.
	subroutine printMatrix( A, name )
		implicit none
		real(KIND=r8), dimension( :, : ) :: A   ! Assume the shape of A.
		character name  ! Name for use in assignment, ie, name = ......

		integer n, m, i, j

		n = size( A, 1 )
		m = size( A, 2 )

		write( *, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      ! Output the matrix, except for the last row, which needs no `;'.
		do i = 1, n-1

         ! Output current row.
			do j = 1, m-1
				write( *, fmt="(f10.6,a2)", advance = "no" ) A( i, j ), ', '
			end do 

         ! Output last element in row and end current row.
			write( *, fmt="(f10.6,a1)" ) A( i, m ), ';'

		end do 

      ! Output the last row.
		do j = 1, m-1
			write( *, fmt="(f10.6,a2)", advance = "no" ) A( i, j ), ', '
		end do 

      ! Output last element in row and end.
		write( *, fmt="(f10.6,a1)" ) A( i, m ), ']'

	end subroutine printMatrix


   ! Output b in Matlab format, using name in the Matlab assignment statement.
	subroutine printVector( b, name )
		implicit none
		real(KIND=r8), dimension( : ) :: b   ! Assume the shape of b.
		character name   ! Name for use in assignment, ie, name = ......

		integer n, i

		n = size( b )

		write( *, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

		do i = 1, n-1
			write( *, fmt = "(f10.6,a2)", advance = "no" ) b( i ), ', '
		end do

		write( *, fmt = "(f10.6,a2)" ) b( n ), ']'''

	end subroutine printVector

	!**************************************************************************
! Features of LinSolv4:
! LinSolv4 like Linear_Solver but uses FULL (simultaneous row & column) pivoting
!  Postulate:
!  1. Changing columns (i.e. changing the order of variables) does not affect the determinant.
!  Pivoting:-
!     1. For the submatrix led by the element at position (j,j), find the largest element 
!        as the pivot. Record its position in prow, pcol
!     2. Swap rows and cols to make the pivot at position (j,j)

	!COMPOSITE AND SMART STRUCTURES GROUP (CSSG)
	!School of Aerospace, Mechanical and Mechatronic Engineering,
	!Building J07 University of Sydney,
	!Sydney, NSW, Australia, 2006.
	!Tel: 61-2-9351 2338 Fax: 61-2-9351 4841

	subroutine LinSolv4(a,b,row,soln,ans,detmt)
		implicit none
		INTEGER, INTENT(IN) :: row,soln
		REAL(KIND=r8) , INTENT(IN), TARGET :: a(row,row), b(row,soln)
		REAL(KIND=r8) , INTENT(OUT) :: ans(row,soln)
		REAL(KIND=r8) , INTENT(OUT) :: detmt

		INTEGER i,j, p          	! need to specify p with dimension one to
                                 ! use MAXLOC intrinsic function
		INTEGER :: memoryvec(row), iw,ix,iy,iz, itmp, prow, pcol
		REAL(KIND=r8) , POINTER  :: swap_ij(:), c(:,:)
		REAL(KIND=r8)  :: tmp
		REAL(KIND=r8), ALLOCATABLE :: ctmp(:,:)

		ALLOCATE (c(row,row+soln), swap_ij(row+soln))
		ALLOCATE (ctmp(row,row+soln))
		detmt = 1.0d0

! ERROR Check #1
! Checking for rows of zeros (Theorem 1)
		do i = 1, row
			tmp = 0.0d0
			do j = 1, row 
				tmp = tmp + dabs( a(i,j) )
			end do
			if (tmp.lt.smallepsil*smallepsil) then
				PRINT *, 'Error the sum of row ',i,' is less than ', smallepsil*smallepsil
				PRINT *, 'High possibility that matrix is singular!!!'
				call EXIT(10)
			end if
		end do

! Constructing ONE single augmented matrix (actually c is a pointer)
		do i = 1, row
			c(i, 1:row) = a(i,:)
			c(i, row+1 : row+soln) = b(i,:)
		end do

! Initializing memory vector for the index/position for solution 
! Must be used if there is column swapping
		do i = 1, row
			memoryvec(i) = i
		end do

! Begin Row-Reduction Operations
		do j = 1, row                 ! for the j-th column
			WRITE(60,*)'Current Matrix: Original matrix ',j,'-th operation'
			do p = 1,row
				WRITE(60,500) (c(p,i), i = 1,row+soln)
			end do
			WRITE(60,510)

!  PIVOTING - of the submatrix for the j-th operation
!  Search for Pivot from whole submatrix led by (j,j)
			tmp = 0.0d0;      prow = j; 		pcol = j
			do i = j, row
				do p = j, row
					if ( tmp .lt. dabs(c(i,p)) ) then         ! tmp < c(i,p)
						tmp = dabs(c(i,p))
						prow = i;   pcol = p
					end if
				end do
			end do         ! i loop

! Swap columns NB swapping column must be accompanied by swapping row.
! Swap row/col  pcol <->  j
			if (pcol.ne.j) then
				WRITE(60,*) 'col change from - to ',pcol,j
				do iw = 1, row
					if (iw.eq.pcol) then
						ix=j
					elseif (iw.eq.j) then
						ix=pcol
					else
						ix=iw
					ENDIF
					do iy = 1, row
						if (iy.eq.pcol) then
							iz=j
						elseif (iy.eq.j) then
							iz=pcol
						else
							iz=iy
						ENDIF
						ctmp(ix,iz) = c(iw,iy)     ! Building up temporary matrix to represent
					end do     ! end loop iy      ! the matrix where the row/col is changed imm<->i
					do iy = row+1 , row+soln
						ctmp(ix,iy) = c(iw,iy)
					end do
				end do    	  ! end for iw
				itmp = memoryvec(j)             	! Keep track of swapped variables/solutions
				memoryvec(j) = memoryvec(pcol)
				memoryvec(pcol) = itmp
! Actual swapping of row/col imm<->j
				c(:,:) = ctmp(:,:)
! Row position of pivot might also change after changing columns/rows
				if (prow.eq.j) then
					prow = pcol
				ELSEIF (prow.eq.pcol) then
					prow = j
				end if
			end if               ! end column swap
                           ! ASSUMING swapping col/row do NOT change determinant



! Swap rows
			if (prow .ne. j) then
				swap_ij(:) = c(j,:);    c(j,:) = c(prow,:);   c(prow,:) = swap_ij(:)
				detmt = -1.0d0 * detmt
			end if
! Determinant change signs if rows are swapped (Theorem 3b.)


! ERROR Check #2
! If after PIVOTING the diagonal element(now having the largest value)
! is still very small then Matrix is singular
! provided c(j,j) is not the last diagonal element
			if ( (dabs(c(j,j)) < 1.0d-12*smallepsil).AND.(j.ne.row)) then
				PRINT *, 'ERROR: Matrix is Singular. Found at ',j,'-th row operation'
				PRINT *, 'diagonal element is ', c(j,j)
				call exit(10)
			end if

! ERROR Check #3
! If at the j-th row, all other elements are zero and element(j,j) is very small
! then might be singular or inconsistent matrix
			tmp = 0.0d0
			do i = 1, j-1
				tmp = tmp + dabs(c(j,i))
			end do
			do i = j+1, row
				tmp = tmp + dabs(c(j,i))
			end do                 
			if ( (tmp.lt.smallepsil) .AND. (dabs(c(j,j)).lt.1.0d-12*smallepsil) ) then
				PRINT *, 'ERROR: Matrix is Singular/Inconsistent. Found at ',j,'-th row operation'
				PRINT *, 'Diagonal element is too small ', c(j,j)
				PRINT *, 'And all other elements in this row are almost zero'
				call exit(10)
			end if

! Divide the j-th row by leading diagonal
			tmp = c(j,j)
			c(j,j) = 1.0d0
			DO i = j+1, row+soln
				c(j,i) = c(j,i) / tmp
			END do  


! Finding Determinant as the factors of the diagoanals (Theorem 3a.)
			detmt = detmt * tmp

! Subtract multiple of j-th row from all rows (except j-th row itself)
! This leaves the j-th column with only one "1" in the diagonal position.
			do i = 1, row                                      ! 1 0 ~ ~
				if (i .ne. j) then                              ! 0 1 ~ ~
					tmp = c(i,j)                                 ! 0 0 1 ~
					c(i,j) = 0.0d0                               ! 0 0 ~ ~
					write (60,*) i, tmp
					do p = j+1, row+soln
						c(i,p) = c(i,p) - tmp * c(j,p)
					end do
				end if
			end do
		end do      ! j loop

		ans(memoryvec(:),:) = c(:,row+1:row+soln)

		deallocate (c, ctmp, swap_ij)

		WRITE(60,*)'Linear Solver Subroutine: Original matrix'
		do j = 1, row
			WRITE(60,500) (a(j,i), i = 1,row)
		end do
		WRITE(60,510)

		WRITE(60,*)'Linear Solver Subroutine: Solution matrix'
		do j = 1, row
			WRITE(60,500) (b(j,i), i = 1,soln)
		end do
		WRITE(60,510)

		WRITE(60,*)'Result from Linear Solver Subroutine'
		do j = 1, row
			WRITE(60,500) (ans(j,i), i = 1,soln)
		end do
		WRITE(60,510)
		WRITE(60,*)'determinant is ',detmt

		500  FORMAT(300g15.4)
		510  FORMAT(//)

	END subroutine LinSolv4
end MODULE Inverters
