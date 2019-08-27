MODULE GMRes_Tree
  !v2	21 Nov	implementing pointers
  !v3	24 Nov	removing pointers, changing subroutine minimize
  !v4	29 Nov	changing to Jerry's newest, removing 'minimize'
  !v5	06 Dec	v4 working! Changing to single row allocatable arrays
  !				  to ease stack handling
  !v6   27 Feb  v5 still working. Changing to treecode for all i=1,n and j=1,n
  !				  matrix*vector computations. Aim to minimize mat.vec multiplications
  !v7   08 Mar  v6 working. Adding restarts. Revamping Mlk, Tlk calls
  !v8	25 Mar	v7 working for non-tree. Changing Tlk calls again

	use TreeCode, only			: TreeCode_Matrix
	use PotentialTree, only		: ComputeAn
	use TreeOperations, only	: ClusterMoment, TaylorCoeff
	use Inverters, only			: LEGS, gaussianElimination
	use TreeCodeGlobal

	implicit none

	! used to restrict variables to prevent accidental calling outside
	private      ! everything defaults to 'private' except those expressly labeled otherwise.
	public	:: gmresIter

contains 

	SUBROUTINE gmresIter(x, A, b, n, tol)
		! incoming variables
		integer, intent(in)							:: n
		real(kind=r8), intent(in)					:: tol
		real(kind=r8), dimension(n),intent(in) 		:: b
		real(kind=r8), dimension(n,n), intent(in) 	:: A

		! outgoing variable
		real(kind=r8), dimension(n),intent(inout)	:: x

  		! regular variables
		integer							:: i, m, k, flag, outer, inner, imin, &
			jmin, initer, outiter
		integer, dimension(n)			:: indx
		integer, dimension(2)			:: iter
		real(kind=r8)					:: n2b, tolb, relres, normr, normrmin, alpha
		real(kind=r8), dimension(n)		:: r, u, v, w, xmin, y, ytemp, mattemp, &
			additive
		real(kind=r8), dimension(2,2)	:: Jhat
		real(kind=r8), dimension(n,n)	:: U1, Psol, R1, J1, Jtemp		

		continue
  		!**********
  		! main program commands begin
  		!**********
  		! defining parameters in main program
		outer=1; inner=n

		! Check for all zero right hand side vector => all zero solution
		n2b = norm(b)					! Norm of rhs vector, b
		if (n2b == 0) then            	! if rhs vector is all zeros
			print *, 'the answer is x={ 0 }'
			stop
		endif

		!Set up for the method
		flag = 1
		xmin = x                        ! Iterate which has minimal residual so far
		tolb = tol * n2b;               ! Relative tolerance

		mattemp = MATMUL(A,x)
		do i=1,n
			r(i) = b(i) - mattemp(i)
		enddo
		normr = norm(r)  				! Norm of residual

		if (normr <= tolb) then         ! Initial guess is good enough
			flag = 0
			relres = normr / n2b
			stop
		endif

		normrmin = normr            	! Norm of residual from xmin

		do outiter = 1, outer
			u = r + sign(1._r8,r(1))*normr*unit_vec(1,n)
			u = u / norm(u)
			U1(:,:)=0
			Psol(:,:)= 0
			R1(:,:)= 0
			U1(:,1) = u

			w = r - 2.* u * DOT_PRODUCT(u,r)

			do initer= 1,inner
				v = unit_vec(initer,n) - 2*u*u(initer)
				do k = (initer-1),1,-1
					v = v - 2*U1(:,k) * DOT_PRODUCT(U1(:,k),v)
				enddo

				! give P1*P2*P3...Pm*em
				Psol(:,initer) = v
				v = MATMUL(A, v)

				do k = 1,initer
					relres = DOT_PRODUCT(U1(:,k),v)
					v = v - 2*U1(:,k) * DOT_PRODUCT(U1(:,k),v)
				enddo

				! gives Pm*Pm-1*...P1*A*P1*P2*..Pm*em
				! determine Pm+1
				if (.NOT. (initer==SIZE(v) .OR. all( v(initer+1:n)==0) )) then
					u(:) = 0
					alpha = -sign(1._r8,v(initer+1)) * norm(v(initer+1:n))
					u(initer+1:n) = v(initer+1:n) - alpha*unit_vec(1,n-initer)
					u = u / norm(u)
					U1(:,initer+1) = u

					! apply Pm+1 to v
					v = v - 2*U1(:,k) * DOT_PRODUCT(U1(:,k),v)
				endif

				if (1==initer) then
					J1(:,:) = 0.0; do m = 1, n; J1(m,m) = 1.0; enddo
				else
					J1 = MATMUL(Jtemp,J1)
					v = MATMUL(J1,v)
				endif

				! find Given's rotation Jm
				if (.NOT. ( initer==SIZE(v) )) then
					if ( 0==v(initer+1) ) then
						Jtemp(:,:) = 0.0; do m = 1, n; Jtemp(m,m) = 1.0; enddo
					else
						Jtemp = 0.0; do m = 1, n; Jtemp(m,m) = 1.0; enddo
						CALL planerot(Jhat,v(initer:initer+1), v(initer:initer+1) )
						Jtemp(initer:initer+1,initer:initer+1) = Jhat
						w(initer:initer+1) = matmul(Jhat,w(initer:initer+1))
					endif
				endif

				R1(:,initer) = v

				if ( initer<inner ) then
					normr = abs(w(initer+1))
				endif

				if ( normr <= normrmin ) then
					normrmin = normr
					imin = outiter
					jmin = initer
				endif

				if ( normr < tolb ) then
					flag = 0
					iter = [outiter, initer]
					exit
				endif
			enddo        ! ends innner loop

			!y = R1(1:jmin,:) \ w(1:jmin)
			!SUBROUTINE LEGS (A,N,B,X,INDX)			
			CALL LEGS(R1(1:jmin,:), jmin, w(1:jmin), ytemp, indx)
			do m=1,jmin !(? unsure about jmin...)
				y(m)=ytemp(indx(m))
			enddo

			!call gaussianElimination( R1(1:jmin,:), w(1:jmin), y, GaussYes )
			!if (.NOT.(GaussYes)) then
			!stop '  ** Gauss elimination failed! **'
			!endif

			additive = matmul(Psol,y)
			x = x + additive
			xmin = x

			r = b - matmul(A,x)
			normr = norm(r)

			if ( normr <= normrmin ) then
				xmin = x
				normrmin = normr
			endif

			if ( norm(additive) < 1e-12 ) then
				exit		! no change in outer iterate
			endif

			if ( normr < tolb ) then
				flag = 0
				iter = [outiter, initer]
				exit
			endif

			print *, 'normr=', normr
		enddo        		! ends outer loop

		! returned solution is that with minimum residual
		if ( flag == 0 ) then
			relres = normr / n2b
		else
			x = xmin
			iter = [imin, jmin]
			relres = normrmin / n2b
		endif

		write(*, fmt=10), 'GMRes took ', iter(1)+iter(2)-1, ' iterations.'
		write(*, fmt=11), 'GMRes residual is ', relres
		goto 999

		!*****************
		! FORMAT STATEMENTS
		10 FORMAT (a, i5, a)
		11 FORMAT (a, e8.3)
	999 END SUBROUTINE gmresIter

	FUNCTION unit_vec(k,n)
		integer, intent(in)	:: k,n
		real(kind=r8), dimension(n) :: unit_vec

		continue

		unit_vec(:) = 0.; unit_vec(k) = 1.
	END FUNCTION unit_vec

	FUNCTION allnoZero(v,sizev)
		integer, intent(in)		:: sizev
		real(kind=r8), dimension(sizev), intent(in) :: v
		integer 				:: allnoZero

		integer					:: i

		continue

		allnoZero=1
		do i=1,sizev
			if (0==v(i)) then
				allnoZero=0
			endif
		enddo
	END FUNCTION allnoZero

	FUNCTION norm(x)
		real(kind=r8), dimension(:), intent(in) :: x
		real(kind=r8) 					:: norm

		integer					:: i, sizex

		continue

		sizex=size(x)
		norm=0.
		do i=1,sizex
			norm=norm+1.*x(i)*x(i)
		enddo
		norm=SQRT(norm)
	END FUNCTION norm

	SUBROUTINE DROTG(da,db,c,s)
	! http://www.netlib.org/blas/drotg.f
	! construct givens plane rotation.
	! jack dongarra, linpack, 3/11/78.

		real(kind=r8), intent(inout)	:: da, db, c, s

		real(kind=r8)					:: roe, scale, r, z

		continue 

		if (abs(da) > abs(db)) then
			roe=da
		else 
			roe=db
		endif

		scale = abs(da) + abs(db)

		if ( scale /= 0) then
			r = scale*SQRT( (da/scale)**2 + (db/scale)**2 )
			r = sign(1.0_r8,roe)*r
			s = db/r
			z = 1.0_r8

			if ( abs(da) > abs(db) ) then
				z=s
			endif

			if ( abs(db) >= abs(da) .AND. ( c /= 0. ) ) then
				z=1.0_r8/c
			endif
		else
			s = 0.0_r8
			r = 0.0_r8
			z = 0.0_r8
		endif

		da=r
		db=z
	END SUBROUTINE DROTG

	SUBROUTINE planerot(G, yPlane, x)
		real(kind=r8), dimension(2), intent(in)	:: x

		real(kind=r8), dimension(2), intent(out)	:: yPlane
		real(kind=r8), dimension(2,2), intent(out) :: G

		real(kind=r8)	:: rPlane

		continue

		!PLANEROT Givens plane rotation.
		![G,Y] = PLANEROT(X), where X is a 2-component column vector,
		!returns a 2-by-2 orthogonal matrix G so that Y = G*X has Y(2) = 0.

		if (x(2) /= 0 ) then
			rPlane = norm(x)
			G(1,:)= [x(1)/rPlane, x(2)/rPlane]
			G(2,:)= [-x(2)/rPlane, x(1)/rPlane]
		!G = [x'; -x(2) x(1)]/r;
			yPlane = [rPlane, 0._r8]
		else
			G(1,:) = [1._r8, 0._r8]
			G(2,:) = [0._r8, 1._r8]
		! G = eye(2,class(x));
		endif
	END SUBROUTINE planerot
end MODULE GMRes_Tree
