!=============
! input:
! center stencil:
! data(0:bfact-1)   gauss points
! 
! bctype(1) left neumann, dirichlet, interior
! bctype(2) right neumann, dirichlet, interior
! bcvalue(1) or bcvalue(2)
! 
! if interior boundary, 
!then the location is at the neighboring Gauss point.
! if exterior boundary, 
!then the location is at the element boundary.
! 
! left and right interior boundary case:
! x | x x     x x | x
! 
! left not interior:
!   | x x     x x | x
!   x
! 
! interpolate to:
!   |x x       x x|  GL points: 0..bfact+1
!   x             x
! 
! gradient:
!   | x    x     x | GL points: 0..bfact
!   x              x
!=========================

! weakflag=0
 subroutine lineGRAD_new(bctype,bcvalue,source,DP,bfact,dx)
 IMPLICIT NONE

 REAL_T dx
 INTEGER_T bctype(1:2)
 REAL_T bcvalue(2)
 INTEGER_T bfact
 REAL_T source(0:bfact-1) 
 REAL_T DP(0:bfact)
 
 REAL_T y(0:bfact-1)
 REAL_T yGL(0:bfact)
 REAL_T wMAT(0:bfact-1,0:bfact-1)
 REAL_T wMATGL(0:bfact,0:bfact)
 REAL_T y_extend(0:bfact+1)
 REAL_T yGL_extend(0:bfact+1)
 REAL_T wMAT_extend(0:bfact+1,0:bfact+1)
 REAL_T yleftbc(0:bfact+1)
 REAL_T yrightbc(0:bfact+1)
 REAL_T w_leftbc(0:bfact+1,0:bfact+1)
 REAL_T w_rightbc(0:bfact+1,0:bfact+1)
 INTEGER_T i,j,i1,j1
 REAL_T dx_element
 REAL_T PLINE(0:bfact+1)
 REAL_T PLINE2(0:bfact+1)
 REAL_T pts(0:bfact+1)
 REAL_T w_pts(0:bfact+1,0:bfact+1)
 REAL_T sum,f

  do i1=0,bfact-1
   y(i1)=cache_gauss(bfact,i1)
  enddo
  do i1=0,bfact
   yGL(i1)=cache_gauss_lobatto(bfact,i1)
  enddo
  do j1=0,bfact-1
   do i1=0,bfact-1
     wMAT(i1,j1)=cache_wMAT(bfact,i1,j1)
   enddo
  enddo
  do j1=0,bfact
   do i1=0,bfact
    wMATGL(i1,j1)=cache_wMATGL(bfact,i1,j1)
   enddo
  enddo

! precomputed
!   do j1=0,bfact
!    do i1=0,bfact
!     w_right(i1,j1)=cache_w_right(bfact,i1,j1)
!    enddo
!   enddo
!   do j1=0,bfact
!    do i1=0,bfact
!     w_left(i1,j1)=cache_w_left(bfact,i1,j1)
!    enddo
!   enddo

!------for interior
  y_extend(0)=-one-abs(y(0)+one)
  do i=1,bfact
   y_extend(i)=y(i-1)
  enddo
  y_extend(bfact+1)=one+abs(one-y(bfact-1))
  do i1=0,bfact+1
   yGL_extend(i1)=cache_gauss_lobatto(bfact+1,i1)
  enddo

!------left or right boundary
  yleftbc(0)=-one
  yrightbc(bfact+1)=one
  do i1=1,bfact
   yleftbc(i1)=y(i1-1)
   yrightbc(i1)=y(i1)
  enddo
  yleftbc(bfact+1)=one+abs(one-y(bfact-1))
  yrightbc(0)=-one-abs(y(0)+one)

  do i1=1,bfact+1
   pts(i1)=yrightbc(i1)
  enddo
  pts(0)=-one

! precomputed??
  call polyinterp_Dmatrix(bfact+1,yGL_extend,wMAT_extend)
  call polyinterp_Dmatrix(bfact+1,yrightbc,w_rightbc)
  call polyinterp_Dmatrix(bfact+1,yleftbc,w_leftbc)
  call polyinterp_Dmatrix(bfact+1,pts,w_pts)

!------------
  PLINE(0)=bcvalue(1)
  PLINE(bfact+1)=bcvalue(2) 
  do i1=1,bfact
   PLINE(i1)=source(i1-1)
  enddo


  i=1
  j=2
  if (bctype(i)==0) then
   !left interior

   if (bctype(j)==0) then
    !right interior
    call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
     y_extend,yGL_extend)
    dx_element=dx*bfact
    call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
     wMAT_extend,yGL_extend,yGL,dx_element)

    elseif (bctype(j)==1) then
     !right Dirichlet
     call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
      yrightbc,yGL_extend)
     dx_element=dx*bfact
     call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
      wMAT_extend,yGL_extend,yGL,dx_element)

    elseif (bctype(j)==2) then
     !right Neumann
     !yrightbc
     sum=bcvalue(1)*w_rightbc(0,bfact+1) 
     do i1=1,bfact
     sum=sum+source(i1-1)*w_rightbc(i1,bfact+1)
     enddo
     sum=bcvalue(2)*(dx_element/two)-sum
     PLINE(bfact+1)=sum/w_rightbc(bfact+1,bfact+1)
     call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
      yrightbc,yGL_extend)
     dx_element=dx*bfact
     call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
      wMAT_extend,yGL_extend,yGL,dx_element)
      
    else
      print*,'Not valid bc type for right(interior)'
    endif


  elseif (bctype(i)==1) then
   !left Dirichlet

   if (bctype(j)==0) then
    !right interior
    call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
     yleftbc,yGL_extend)
    dx_element=dx*bfact
    call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
     wMAT_extend,yGL_extend,yGL,dx_element)

    elseif (bctype(j)==1) then
     !right Dirichlet
     call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
      pts,yGL_extend)
     dx_element=dx*bfact
     call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
      wMAT_extend,yGL_extend,yGL,dx_element)
      
    elseif (bctype(j)==2) then
     !right Neumann
     !pts
     sum=bcvalue(1)*w_pts(0,bfact+1) 
     do i1=1,bfact
     sum=sum+source(i1-1)*w_pts(i1,bfact+1)
     enddo
     sum=bcvalue(2)*(dx_element/two)-sum
     PLINE(bfact+1)=sum/w_pts(bfact+1,bfact+1)
     call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
      pts,yGL_extend)
     dx_element=dx*bfact
     call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
      wMAT_extend,yGL_extend,yGL,dx_element)


    else
      print*,'Not valid bc type for right(Dirichlet)'
    endif


  elseif (bctype(i)==2) then
   !left Neumann 

   if (bctype(j)==0) then
    !right interior
    !yleftbc
     sum=bcvalue(2)*w_leftbc(bfact+1,0) 
     do i1=1,bfact
     sum=sum+source(i1-1)*w_leftbc(i1,0)
     enddo
     sum=bcvalue(1)*(dx_element/two)-sum
     PLINE(0)=sum/w_leftbc(0,0)
     call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
      yleftbc,yGL_extend)
     dx_element=dx*bfact
     call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
      wMAT_extend,yGL_extend,yGL,dx_element)

    elseif (bctype(j)==1) then
     !right Dirichlet
     !pts
     sum=bcvalue(2)*w_pts(bfact+1,0) 
     do i1=1,bfact
     sum=sum+source(i1-1)*w_pts(i1,0)
     enddo
     sum=bcvalue(1)*(dx_element/two)-sum
     PLINE(0)=sum/w_pts(0,0)
     call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
      pts,yGL_extend)
     dx_element=dx*bfact
     call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
      wMAT_extend,yGL_extend,yGL,dx_element)

    elseif (bctype(j)==2) then
     !right Neumann
     !pts
     sum=0.0D0 
     do i1=1,bfact
     sum=sum+source(i1-1)*w_pts(i1,0)
     enddo
     sum=sum*(dx_element/two)-sum
     f=0.0D0 
     do i1=1,bfact
     f=f+source(i1-1)*w_pts(i1,bfact+1)
     enddo
     f=f*(dx_element/two)-f
     PLINE(bfact+1)=(sum*w_pts(0,bfact+1)- &
       f*w_pts(0,0))/(w_pts(bfact+1,0)*w_pts(0,bfact+1)- &
       w_pts(0,0)*w_pts(bfact+1,bfact+1))
     PLINE(0)=sum-(PLINE(bfact+1)*w_pts(bfact+1,0))/(w_pts(0,0))

     call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
      pts,yGL_extend)
     dx_element=dx*bfact
     call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
      wMAT_extend,yGL_extend,yGL,dx_element)

    else
      print*,'Not valid bc type for right(Neumann)'
    endif


  else
    print*,'Not valid bc type for left'
  endif


  return
 end subroutine lineGRAD_new
