      program main
      IMPLICIT NONE
* cd amrex_implicit_interfaces/NS_fluids_lib/BLASLAPACK/LINEAR
* make
* copy "liblinear.a" into the same director as svd_experiment.f
* gfortran svd_experiment.f -L. -llinear
* ./a.out
      integer M,N
      parameter (M=6,N=5)
      integer LDA
      integer LDU
      integer LDVT
      parameter (LDA=M,LDU=M,LDVT=N)
      integer LWMAX
      parameter (LWMAX=1000)
      integer INFO,LWORK
      integer i,j,k

      DOUBLE PRECISION A( LDA, N ), U( LDU, M ), VT( LDVT, N ), S( N ),
     $                 WORK( LWMAX )
      DOUBLE PRECISION USIGMA(LDA,N)
      DOUBLE PRECISION ATEST(LDA,N)
      DOUBLE PRECISION ACOPY(LDA,N)
      DOUBLE PRECISION frobenius_norm
      DATA             A/
     $  8.79, 6.11,-9.15, 9.57,-3.49, 9.84,
     $  9.93, 6.91,-7.93, 1.64, 4.02, 0.15,
     $  9.83, 5.04, 4.86, 8.83, 9.80,-8.99,
     $  5.45,-0.27, 4.85, 0.74,10.00,-6.02,
     $  3.16, 7.98, 3.01, 5.80, 4.27,-5.31
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         DGESVD
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DGESVD Example Program Results'

      CALL PRINT_MATRIX( 'A',
     $                   M, N, A, LDA )

      do i=1,M
       do j=1,N
        ACOPY(i,j)=A(i,j)
       enddo
      enddo
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
     $             WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
     $             WORK, LWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF
*
*     Print singular values.
*
      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
*
*     Print left singular vectors.
*
      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)',
     $                   M, N, U, LDU )
*
*     Print right singular vectors.
*
      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)',
     $                   N, N, VT, LDVT )

      do i=1,M
       do j=1,N
        USIGMA(i,j)=0.0d0
        do k=1,M
         if (k.eq.j) then
          USIGMA(i,j)=USIGMA(i,j)+U(i,k)*S(k)
         endif
        enddo
       enddo
      enddo
      do i=1,M
       do j=1,N
        ATEST(i,j)=0.0d0
        do k=1,N
         ATEST(i,j)=ATEST(i,j)+USIGMA(i,k)*VT(k,j)
        enddo
       enddo
      enddo
      CALL PRINT_MATRIX( 'ATEST',
     $                   M, N, ATEST, LDA )

      frobenius_norm=0.0d0
      do i=1,M
       do j=1,N
        frobenius_norm=frobenius_norm+(ACOPY(i,j)-ATEST(i,j))**2
       enddo
      enddo
      frobenius_norm=sqrt(frobenius_norm)
      print *,"frobenius_norm= ",frobenius_norm
      STOP
      END
*
*     End of DGESVD Example.
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
      END

