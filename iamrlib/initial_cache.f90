      subroutine init_cache(order_r,typ)
      IMPLICIT NONE

      INTEGER_T order_r
      REAL_T yleft(0:order_r)
      REAL_T yright(0:order_r)
      REAL_T y(0:order_r)
      REAL_T yGL(0:order_r)
      REAL_T y_extend(0:order_r+1)
      REAL_T yGL_extend(0:order_r+1)
      REAL_T yLT(0:order_r+1)
      REAL_T yRT(0:order_r+1)
      REAL_T yLRT(0:order_r+1)

      REAL_T, dimension(:,:), allocatable :: deriv_matrix
      REAL_T, dimension(:), allocatable :: tempx,tempw

      INTEGER_T i,j,typ,i1,j1

      if (order_r.lt.1) then
       print *,"order_r invalid"
       stop
      endif

      allocate(cache_gauss(1:order_r+1,0:order_r+1))
      allocate(cache_gauss_lobatto(1:order_r+1,0:order_r+1))
      allocate(cache_gauss_w(1:order_r+1,0:order_r+1))
      allocate(cache_gauss_lobatto_w(1:order_r+1,0:order_r+1))

      allocate(cache_wMATGL(1:order_r,0:order_r,0:order_r))
      allocate(cache_wMAT(1:order_r,0:order_r,0:order_r))
      allocate(cache_w_left(1:order_r,0:order_r,0:order_r))
      allocate(cache_w_right(1:order_r,0:order_r,0:order_r))

      allocate(cache_wMAT_extend(1:order_r,0:order_r+1,0:order_r+1))
      allocate(cache_wRT(1:order_r,0:order_r+1,0:order_r+1))
      allocate(cache_wLT(1:order_r,0:order_r+1,0:order_r+1))
      allocate(cache_wLRT(1:order_r,0:order_r+1,0:order_r+1))



      do i=1,order_r+1
      
allocate(tempx(0:i))
allocate(tempw(0:i))
      
       if (typ.eq.0) then
        call LegendreGaussLobattoNodesAndWeights(i,tempx,tempw) 
       else if (typ.eq.1) then
        call ClenshawGaussLobattoNodesAndWeights(i,tempx,tempw) 
       else
        print *,"typ invalid"
        stop
       endif

       do j=0,i
        cache_gauss_lobatto(i,j)=tempx(j)
        cache_gauss_lobatto_w(i,j)=tempw(j)
       enddo

       if (typ.eq.0) then
        call LegendreGaussNodesAndWeights(i-1,tempx(0:i-1),tempw(0:i-1)) 
       else if (typ.eq.1) then
        call ClenshawGaussNodesAndWeights(i-1,tempx(0:i-1),tempw(0:i-1)) 
       else
        print *,"typ invalid"
        stop
       endif

       do j=0,i-1
        cache_gauss(i,j)=tempx(j)
        cache_gauss_w(i,j)=tempw(j)
       enddo

deallocate(tempx,tempw)

      enddo ! i=1.. order_r+1

      do i=1,order_r

       do i1=0,i-1
        y(i1)=cache_gauss(i,i1)
       enddo
       do i1=0,i
        yGL(i1)=cache_gauss_lobatto(i,i1)
       enddo
       yright(0)=-one
       yleft(i)=one
       do i1=1,i
        yright(i1)=y(i1-1)
        yleft(i1-1)=y(i1-1)
       enddo

       y_extend(0)=-one-abs(y(0)+one)
       do i1=1,i
        y_extend(i1)=y(i1-1)
       enddo
       y_extend(i+1)=one+abs(one-y(i-1))
       do i1=0,i+1
        yGL_extend(i1)=cache_gauss_lobatto(i+1,i1)
       enddo
       yLT(0)=-one
       yLT(i+1)=y_extend(i+1)
       yRT(0)=y_extend(0)
       yRT(i+1)=one
       yLRT(0)=-one
       yLRT(i+1)=one
       do i1=1,i
        yLT(i1)=y(i1-1)
        yRT(i1)=y(i1-1)
        yLRT(i1)=y(i1-1)
       enddo ! i1

       allocate(deriv_matrix(0:i+1,0:i+1))
       call polyinterp_Dmatrix(i+1,yGL_extend,deriv_matrix)
       do i1=0,i+1
       do j1=0,i+1
        cache_wMAT_extend(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i+1,0:i+1))
       call polyinterp_Dmatrix(i+1,yRT,deriv_matrix)
       do i1=0,i+1
       do j1=0,i+1
        cache_wRT(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i+1,0:i+1))
       call polyinterp_Dmatrix(i+1,yLT,deriv_matrix)
       do i1=0,i+1
       do j1=0,i+1
        cache_wLT(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i+1,0:i+1))
       call polyinterp_Dmatrix(i+1,yLRT,deriv_matrix)
       do i1=0,i+1
       do j1=0,i+1
        cache_wLRT(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)


       allocate(deriv_matrix(0:i,0:i))
       call polyinterp_Dmatrix(i,yGL,deriv_matrix)
       do i1=0,i
       do j1=0,i
        cache_wMATGL(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i-1,0:i-1))
       call polyinterp_Dmatrix(i-1,y,deriv_matrix)
       do i1=0,i-1
       do j1=0,i-1
        cache_wMAT(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i,0:i))
       call polyinterp_Dmatrix(i,yright,deriv_matrix)
       do i1=0,i
       do j1=0,i
        cache_w_right(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i,0:i))
       call polyinterp_Dmatrix(i,yleft,deriv_matrix)
       do i1=0,i
       do j1=0,i
        cache_w_left(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

      enddo ! i

      return
      end subroutine init_cache