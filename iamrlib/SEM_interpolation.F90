!============================
! ghost points of fine grid in the coarse grid
!===========================

 subroutine Interpolation_c2f(SDIM,bfact_c,gridtype, &
    lo,hi,dx_c,xfine,xlo,fcoarse,fxfine)
 IMPLICIT NONE
 INTEGER_T SDIM,bfact_c
 INTEGER_T gridtype ! 0:ggg; 1:lgg; 2:glg; 3:ggl
 INTEGER_T lo(1:SDIM),hi(1:SDIM)
 REAL_T dx_c
 REAL_T xfine(1:SDIM)
 REAL_T xlo(1:SDIM)
 REAL_T fcoarse(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
 REAL_T fxfine
 
 INTEGER_T i1,i,j,k
 INTEGER_T index(1:SDIM)
 INTEGER_T ie(1:SDIM)
 REAL_T y(0:bfact_c-1)
 REAL_T yGL(0:bfact_c)
 REAL_T, dimension(:),allocatable :: temp,f
 REAL_T, dimension(:),allocatable :: ypoints,bwG
 REAL_T, dimension(:,:),allocatable :: lg,g
 REAL_T sum
 
 
 if (bfact_c.lt.4) then
  print "Not valid order on the coarse grid"
  stop
 endif
 
 if (SDIM.lt.2) then
  print "Not valid dimension size"
  stop
 endif
 
 if ((SDIM-gridtype).lt.0) then
  print*,'Not valid gridtype'
 endif

 do i1=0,bfact_c-1
  y(i1)=cache_gauss(bfact_c,i1)
 enddo
 do i1=0,bfact_c
  yGL(i1)=cache_gauss_lobatto(bfact_c,i1)
 enddo
 
 do i=1,SDIM
  ie(i)=floor((xfine(i)-xlo(i))/(bfact_c*dx_c))
  index(i)=ie(i)*bfact_c
 enddo
 
  allocate(ypoints(0:bfact_c-1),bwG(0:bfact_c-1))
  allocate(lg(0:bfact_c-1,1:SDIM))
  
  do j=1,SDIM
   do i=0,bfact_c-1
    ypoints(i)=xlo(j)+ie*bfact_c*dx_c+ &
       (y(i)+1.0d0)*(bfact_c*dx_c)*0.5D0
   enddo
    allocate(temp(0:bfact_c-1))
    call BarycentricWeights(bfact_c-1,ypoints,bwG)
    call LagrangeInterpolatingPolynomial(bfact_c-1, &
       xfine(j),ypoints,bwG,temp)
    do i=0,bfact_c-1
     lg(i,j)=temp(i)
    enddo
    deallocate(temp)
  enddo

  deallocate(bwG,ypoints)
  
 if (gridtype==0) then
 ! ggg
  if (SDIM==2) then
   allocate(f(0:bfact_c-1))
   f=0.0d0
   do j=0,bfact_c-1
    do i=0,bfact_c-1
     f(j)=f(j)+fcoarse(index(1)+i,index(2)+j,lo(3))*lg(i,1)
    enddo
   enddo
   sum=0.0d0
   do i=0,bfact_c-1
    sum=sum+f(i)*lg(i,2)
   enddo
   deallocate(f)
  elseif (SDIM==3) then
   allocate(g(0:bfact_c-1,0:bfact_c-1),f(0:bfact_c-1))
   g=0.0d0
   do k=0,bfact_c-1
    do j=0,bfact_c-1
     do i=0,bfact_c-1
      g(j,k)=g(j,k)+fcoarse(index(1)+i,index(2)+j,index(3)+k)*lg(i,1)
     enddo
    enddo
   enddo
   f=0.0d0
   do j=0,bfact_c-1
    do i=0,bfact_c-1
     f(j)=f(j)+g(i,j)*lg(i,2)
    enddo
   enddo
   sum=0
   do i=0,bfact_c-1
    sum=sum+f(i)*lg(i,3)
   enddo
   deallocate(f,g)
  endif
  
  
 elseif
 
  allocate(ypoints(0:bfact_c),bwG(0:bfact_c))
  allocate(temp(0:bfact_c))
  
   do i=0,bfact_c
    ypoints(i)=xlo(gridtype)+ie*bfact_c*dx_c+ &
       (yGL(i)+1.0d0)*(bfact_c*dx_c)*0.5D0
   enddo
    call BarycentricWeights(bfact_c,ypoints,bwG)
    call LagrangeInterpolatingPolynomial(bfact_c, &
       xfine(gridtype),ypoints,bwG,temp)


  if (SDIM==2) then
  
   if (gridtype==1) then
   allocate(f(0:bfact_c-1))
   f=0.0d0
   do j=0,bfact_c-1
    do i=0,bfact_c
     f(j)=f(j)+fcoarse(index(1)+i,index(2)+j,lo(3))*temp(i)
    enddo
   enddo
   sum=0.0d0
   do i=0,bfact_c-1
    sum=sum+f(i)*lg(i,2)
   enddo
   deallocate(f)
   
   elseif (gridtype==2) then
   allocate(f(0:bfact_c))
   f=0.0d0
   do j=0,bfact_c
    do i=0,bfact_c-1
     f(j)=f(j)+fcoarse(index(1)+i,index(2)+j,lo(3))*lg(i,1)
    enddo
   enddo
   sum=0.0d0
   do i=0,bfact_c
    sum=sum+f(i)*temp(i)
   enddo
   deallocate(f)
   endif
   
  elseif (SDIM==3) then
  
  if (gridtype==1) then
   allocate(g(0:bfact_c-1,0:bfact_c-1),f(0:bfact_c-1))
   g=0.0d0
   do k=0,bfact_c-1
    do j=0,bfact_c-1
     do i=0,bfact_c
      g(j,k)=g(j,k)+fcoarse(index(1)+i,index(2)+j,index(3)+k)*temp(i)
     enddo
    enddo
   enddo
   f=0.0d0
   do j=0,bfact_c-1
    do i=0,bfact_c-1
     f(j)=f(j)+g(i,j)*lg(i,2)
    enddo
   enddo
   sum=0
   do i=0,bfact_c-1
    sum=sum+f(i)*lg(i,3)
   enddo
   deallocate(f,g)
   
  elseif (gridtype==2) then
   allocate(g(0:bfact_c,0:bfact_c-1),f(0:bfact_c-1))
   g=0.0d0
   do k=0,bfact_c-1
    do j=0,bfact_c
     do i=0,bfact_c-1
      g(j,k)=g(j,k)+fcoarse(index(1)+i,index(2)+j,index(3)+k)*lg(i,1)
     enddo
    enddo
   enddo
   f=0.0d0
   do j=0,bfact_c-1
    do i=0,bfact_c
     f(j)=f(j)+g(i,j)*temp(i)
    enddo
   enddo
   sum=0
   do i=0,bfact_c-1
    sum=sum+f(i)*lg(i,3)
   enddo
   deallocate(f,g)
   
  elseif (gridtype==3) then
   allocate(g(0:bfact_c-1,0:bfact_c),f(0:bfact_c))
   g=0.0d0
   do k=0,bfact_c
    do j=0,bfact_c-1
     do i=0,bfact_c-1
      g(j,k)=g(j,k)+fcoarse(index(1)+i,index(2)+j,index(3)+k)*lg(i,1)
     enddo
    enddo
   enddo
   f=0.0d0
   do j=0,bfact_c
    do i=0,bfact_c-1
     f(j)=f(j)+g(i,j)*lg(i,2)
    enddo
   enddo
   sum=0
   do i=0,bfact_c
    sum=sum+f(i)*temp(i)
   enddo
   deallocate(f,g)
   
  endif  !! gridtype

  endif  !! SDIM


 deallocate(ypoints,bwG,temp)
 
 
 endif  !gridtype

  deallocate(lg)

 fxfine=sum

 return
 end subroutine Interpolation_c2f
