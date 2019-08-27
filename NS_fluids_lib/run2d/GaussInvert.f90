MODULE GaussInvert

	implicit none

	private
	public :: Simple_Gauss, Pivot_Gauss

contains

	SUBROUTINE Simple_Gauss(Ainv, Grid)		! Gauss elimination - not max pivot.
		real*8	:: 	Grid(:,:)			! The matrix to be reduced.
		real*8	:: 	Ainv(size(Grid,1))			! Returns the solution vector.
		real*8	:: 	G(size(Grid,1),size(Grid,2))		! G is a local work array.
		logical	::	Not_pivot_row(size(Grid,1),size(Grid,2		))	! Pivot row mask.

		integer :: N, L
		continue

		N = size(Grid)				
		G = Grid						! Work on G, not Grid.

		do L=1,N						! G(L,L) is next pivot element.

			if (abs(G(L,L))< 1E-6) then
				stop 'zero encountered in pivot'			
			endif

			G(L,:) = G(L,:)/G(L,L)		

			Not_pivot_row = .true.; Not_pivot_row(L,:) = .false.		
			where ( Not_pivot_row )			&			
			G = G-G(:,spread(L,1,N+1))*G(spread(L,1,N),:)			
		end do							! Repeat for all pivots.

		Ainv = G(:,N+1)					
	END SUBROUTINE Simple_Gauss						

	function  Pivot_Gauss(Grid)			! Gauss elimination, max pivot.
		real*8	:: 	Grid(:,:)			! The matrix to be reduced.
		real*8	:: 	Pivot_Gauss(size(Grid,1)		)	! Returns the solution vector.
		real*8	:: 	G(size(Grid,1),size(Grid,2))			! G is a local work array.
		integer	:: 	P(size(Grid,1),2)			! P is array of pivots.
		logical	::	Not_pivot_row(size(Grid,1),size(Grid,2))				! Mask current pivot row only.
		logical	::	Not_pivot_rows_or_cols(size(Grid,1),size(Grid,1))		! Mask out all
								! previous pivots
								! rows and columns.

		integer :: N, L
		continue

		N = size(Grid)				
		G = Grid							! Work on G, not Grid.

		do L=1,N							! L is next pivot number.

			Not_pivot_rows_or_cols = .true.						
			Not_pivot_rows_or_cols(P(1:L-1,1),:) = .false.			 		 	
			Not_pivot_rows_or_cols(:,P(1:L-1,2)) = .false.			 		 	
			P(L,:) = maxloc(abs(G(:,1:N)),mask=Not_pivot_rows_or_cols)			
			if (abs(G(P(L,1),P(L,2))).lt.1E-4)  stop 'ill-conditioned matrix'				

			G(P(L,1),:) = G(P(L,1),:)/G(P(L,1),P(L,2))						

			Not_pivot_row = .true.; Not_pivot_row(P(L,1),:) = .false.			 			
			where ( Not_pivot_row )	&						
			G = G-G(:,spread(P(L,2),1,N+1))*G(spread(P(L,1),1,N),:)	
		enddo							! Repeat for all pivots.

		Pivot_Gauss(P(:,2)) = G(P(:,1),N+1)
	end function Pivot_Gauss	
END	MODULE GaussInvert
