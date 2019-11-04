#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>


 subroutine gsrb3d(is,ie,js,je,ks,ke, &
   bxleft,bxright,byleft,byright,bzleft,bzright, &
   diag,diagsing,soln,rhs,mask,nsmooth)

 IMPLICIT NONE

 INTEGER_T :: i,j,k,is,ie,js,je,ks,ke,nsmooth,n
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bxleft
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bxright
 REAL_T, dimension(is:ie,js:je,ks:ke) :: byleft
 REAL_T, dimension(is:ie,js:je,ks:ke) :: byright
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bzleft
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bzright
 REAL_T, dimension(is:ie,js:je,ks:ke) :: diag
 REAL_T, dimension(is:ie,js:je,ks:ke) :: diagsing
 REAL_T, dimension(is:ie,js:je,ks:ke) :: soln
 REAL_T, dimension(is:ie,js:je,ks:ke) :: rhs
 REAL_T, dimension(is:ie,js:je,ks:ke) :: mask
 REAL_T, allocatable, dimension(:,:,:) :: ax
 REAL_T, allocatable, dimension(:,:,:) :: solnsave
 REAL_T, allocatable, dimension(:,:,:) :: rhssave
 REAL_T, allocatable, dimension(:,:,:) :: redsoln
 REAL_T, allocatable, dimension(:,:,:) :: blacksoln

 allocate(ax(is:ie,js:je,ks:ke))
 allocate(solnsave(is:ie,js:je,ks:ke))
 allocate(rhssave(is:ie,js:je,ks:ke))
 allocate(redsoln(is:ie,js:je,ks:ke))
 allocate(blacksoln(is:ie,js:je,ks:ke))

 do i=is+1,ie-1
 do j=js+1,je-1
 do k=ks+1,ke-1
  rhssave(i,j,k)=rhs(i,j,k)+ &
    bxleft(i-1,j,k)*soln(i-1,j,k)+bxright(i+1,j,k)*soln(i+1,j,k)+ &
    byleft(i,j-1,k)*soln(i,j-1,k)+byright(i,j+1,k)*soln(i,j+1,k)+ &
    bzleft(i,j,k-1)*soln(i,j,k-1)+bzright(i,j,k+1)*soln(i,j,k+1)- &
    diagsing(i,j,k)*soln(i,j,k)
 enddo
 enddo
 enddo

!$acc region 

 do i=is,ie
 do j=js,je
 do k=ks,ke
  solnsave(i,j,k)=0.0
 enddo
 enddo
 enddo

  do n=1,nsmooth

   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    ax(i,j,k)=rhssave(i,j,k)+ &
     bxleft(i-1,j,k)*solnsave(i-1,j,k)+bxright(i+1,j,k)*solnsave(i+1,j,k)+ &
     byleft(i,j-1,k)*solnsave(i,j-1,k)+byright(i,j+1,k)*solnsave(i,j+1,k)+ &
     bzleft(i,j,k-1)*solnsave(i,j,k-1)+bzright(i,j,k+1)*solnsave(i,j,k+1)- &
     diagsing(i,j,k)*solnsave(i,j,k)
   enddo
   enddo
   enddo
   do i=is,ie
   do j=js,je
   do k=ks,ke
    redsoln(i,j,k)=0.0
    blacksoln(i,j,k)=0.0
   enddo
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    redsoln(i,j,k)=ax(i,j,k)/diag(i,j,k)
   enddo
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    blacksoln(i,j,k)=(ax(i,j,k)+ &
     bxleft(i-1,j,k)*redsoln(i-1,j,k)+bxright(i+1,j,k)*redsoln(i+1,j,k)+ &
     byleft(i,j-1,k)*redsoln(i,j-1,k)+byright(i,j+1,k)*redsoln(i,j+1,k)+ &
     bzleft(i,j,k-1)*redsoln(i,j,k-1)+bzright(i,j,k+1)*redsoln(i,j,k+1))/ &
     diag(i,j,k)
   enddo
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    redsoln(i,j,k)=(ax(i,j,k)+ &
     bxleft(i-1,j,k)*blacksoln(i-1,j,k)+bxright(i+1,j,k)*blacksoln(i+1,j,k)+ &
     byleft(i,j-1,k)*blacksoln(i,j-1,k)+byright(i,j+1,k)*blacksoln(i,j+1,k)+ &
     bzleft(i,j,k-1)*blacksoln(i,j,k-1)+bzright(i,j,k+1)*blacksoln(i,j,k+1))/ &
     diag(i,j,k)
   enddo
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    solnsave(i,j,k)=solnsave(i,j,k)+mask(i,j,k)*redsoln(i,j,k)+ &
     (1.0-mask(i,j,k))*blacksoln(i,j,k)
   enddo
   enddo
   enddo

  enddo

!$acc end region

  do i=is+1,ie-1
  do j=js+1,je-1
  do k=ks+1,ke-1
    soln(i,j,k)=soln(i,j,k)+solnsave(i,j,k)
  enddo
  enddo
  enddo

  deallocate(solnsave)
  deallocate(rhssave)
  deallocate(redsoln)
  deallocate(blacksoln)
  deallocate(ax)

  return
  end


 subroutine icrb3d(is,ie,js,je,ks,ke, &
   bxleft,bxright,byleft,byright,bzleft,bzright, &
   diag,diagsing,soln,rhs,icdiag,mask,nsmooth)

 IMPLICIT NONE

 INTEGER_T :: i,j,k,is,ie,js,je,ks,ke,nsmooth,n
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bxleft
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bxright
 REAL_T, dimension(is:ie,js:je,ks:ke) :: byleft
 REAL_T, dimension(is:ie,js:je,ks:ke) :: byright
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bzleft
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bzright
 REAL_T, dimension(is:ie,js:je,ks:ke) :: diag
 REAL_T, dimension(is:ie,js:je,ks:ke) :: diagsing
 REAL_T, dimension(is:ie,js:je,ks:ke) :: soln
 REAL_T, dimension(is:ie,js:je,ks:ke) :: rhs
 REAL_T, dimension(is:ie,js:je,ks:ke) :: icdiag
 REAL_T, dimension(is:ie,js:je,ks:ke) :: mask
 REAL_T, allocatable, dimension(:,:,:) :: ax
 REAL_T, allocatable, dimension(:,:,:) :: solnsave
 REAL_T, allocatable, dimension(:,:,:) :: rhssave
 REAL_T, allocatable, dimension(:,:,:) :: redsoln
 REAL_T, allocatable, dimension(:,:,:) :: blacksoln

 allocate(ax(is:ie,js:je,ks:ke))
 allocate(solnsave(is:ie,js:je,ks:ke))
 allocate(rhssave(is:ie,js:je,ks:ke))
 allocate(redsoln(is:ie,js:je,ks:ke))
 allocate(blacksoln(is:ie,js:je,ks:ke))

 do i=is+1,ie-1
 do j=js+1,je-1
 do k=ks+1,ke-1
  rhssave(i,j,k)=rhs(i,j,k)+ &
    bxleft(i-1,j,k)*soln(i-1,j,k)+bxright(i+1,j,k)*soln(i+1,j,k)+ &
    byleft(i,j-1,k)*soln(i,j-1,k)+byright(i,j+1,k)*soln(i,j+1,k)+ &
    bzleft(i,j,k-1)*soln(i,j,k-1)+bzright(i,j,k+1)*soln(i,j,k+1)- &
    diagsing(i,j,k)*soln(i,j,k)
 enddo
 enddo
 enddo

!$acc region 

 do i=is,ie
 do j=js,je
 do k=ks,ke
  solnsave(i,j,k)=0.0
 enddo
 enddo
 enddo

  do n=1,nsmooth

   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    ax(i,j,k)=rhssave(i,j,k)+ &
     bxleft(i-1,j,k)*solnsave(i-1,j,k)+bxright(i+1,j,k)*solnsave(i+1,j,k)+ &
     byleft(i,j-1,k)*solnsave(i,j-1,k)+byright(i,j+1,k)*solnsave(i,j+1,k)+ &
     bzleft(i,j,k-1)*solnsave(i,j,k-1)+bzright(i,j,k+1)*solnsave(i,j,k+1)- &
     diagsing(i,j,k)*solnsave(i,j,k)
   enddo
   enddo
   enddo
   do i=is,ie
   do j=js,je
   do k=ks,ke
    redsoln(i,j,k)=0.0
    blacksoln(i,j,k)=0.0
   enddo
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    redsoln(i,j,k)=ax(i,j,k)/diag(i,j,k)
   enddo
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    blacksoln(i,j,k)=(ax(i,j,k)+ &
     bxleft(i-1,j,k)*redsoln(i-1,j,k)+bxright(i+1,j,k)*redsoln(i+1,j,k)+ &
     byleft(i,j-1,k)*redsoln(i,j-1,k)+byright(i,j+1,k)*redsoln(i,j+1,k)+ &
     bzleft(i,j,k-1)*redsoln(i,j,k-1)+bzright(i,j,k+1)*redsoln(i,j,k+1))/ &
     icdiag(i,j,k)
   enddo
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    redsoln(i,j,k)=(ax(i,j,k)+ &
     bxleft(i-1,j,k)*blacksoln(i-1,j,k)+bxright(i+1,j,k)*blacksoln(i+1,j,k)+ &
     byleft(i,j-1,k)*blacksoln(i,j-1,k)+byright(i,j+1,k)*blacksoln(i,j+1,k)+ &
     bzleft(i,j,k-1)*blacksoln(i,j,k-1)+bzright(i,j,k+1)*blacksoln(i,j,k+1))/ &
     diag(i,j,k)
   enddo
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    solnsave(i,j,k)=solnsave(i,j,k)+mask(i,j,k)*redsoln(i,j,k)+ &
     (1.0-mask(i,j,k))*blacksoln(i,j,k)
   enddo
   enddo
   enddo

  enddo

!$acc end region

  do i=is+1,ie-1
  do j=js+1,je-1
  do k=ks+1,ke-1
    soln(i,j,k)=soln(i,j,k)+solnsave(i,j,k)
  enddo
  enddo
  enddo

  deallocate(solnsave)
  deallocate(rhssave)
  deallocate(redsoln)
  deallocate(blacksoln)
  deallocate(ax)

  return
  end


 subroutine icc3d(is,ie,js,je,ks,ke, &
  bxleft,bxright,byleft,byright,bzleft,bzright,diag,diagsing, &
  soln,rhs,icdiag,icbx,icby,icbz,nsmooth)

 IMPLICIT NONE

 INTEGER_T :: i,j,k,is,ie,js,je,ks,ke,nsmooth,n
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bxleft
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bxright
 REAL_T, dimension(is:ie,js:je,ks:ke) :: byleft
 REAL_T, dimension(is:ie,js:je,ks:ke) :: byright
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bzleft
 REAL_T, dimension(is:ie,js:je,ks:ke) :: bzright
 REAL_T, dimension(is:ie,js:je,ks:ke) :: diag
 REAL_T, dimension(is:ie,js:je,ks:ke) :: diagsing
 REAL_T, dimension(is:ie,js:je,ks:ke) :: soln
 REAL_T, dimension(is:ie,js:je,ks:ke) :: rhs
 REAL_T, dimension(is:ie,js:je,ks:ke) :: icdiag
 REAL_T, dimension(is:ie,js:je,ks:ke) :: icbx
 REAL_T, dimension(is:ie,js:je,ks:ke) :: icby
 REAL_T, dimension(is:ie,js:je,ks:ke) :: icbz
 REAL_T, allocatable, dimension(:,:,:) :: ax
 REAL_T, allocatable, dimension(:,:,:) :: solnsave
 REAL_T, allocatable, dimension(:,:,:) :: rhssave
 REAL_T, allocatable, dimension(:,:,:) :: redsoln
 REAL_T, allocatable, dimension(:,:,:) :: blacksoln
 REAL_T YY,XX

 allocate(ax(is:ie,js:je,ks:ke))
 allocate(solnsave(is:ie,js:je,ks:ke))
 allocate(rhssave(is:ie,js:je,ks:ke))
 allocate(redsoln(is:ie,js:je,ks:ke))
 allocate(blacksoln(is:ie,js:je,ks:ke))

 do i=is+1,ie-1
 do j=js+1,je-1
 do k=ks+1,ke-1
  rhssave(i,j,k)=rhs(i,j,k)+ &
    bxleft(i-1,j,k)*soln(i-1,j,k)+bxright(i+1,j,k)*soln(i+1,j,k)+ &
    byleft(i,j-1,k)*soln(i,j-1,k)+byright(i,j+1,k)*soln(i,j+1,k)+ &
    bzleft(i,j,k-1)*soln(i,j,k-1)+bzright(i,j,k+1)*soln(i,j,k+1)- &
    diagsing(i,j,k)*soln(i,j,k)
 enddo
 enddo
 enddo

 do i=is,ie
 do j=js,je
 do k=ks,ke
  solnsave(i,j,k)=0.0
 enddo
 enddo
 enddo

  do n=1,nsmooth

   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    ax(i,j,k)=rhssave(i,j,k)+ &
     bxleft(i-1,j,k)*solnsave(i-1,j,k)+bxright(i+1,j,k)*solnsave(i+1,j,k)+ &
     byleft(i,j-1,k)*solnsave(i,j-1,k)+byright(i,j+1,k)*solnsave(i,j+1,k)+ &
     bzleft(i,j,k-1)*solnsave(i,j,k-1)+bzright(i,j,k+1)*solnsave(i,j,k+1)- &
     diagsing(i,j,k)*solnsave(i,j,k)
   enddo
   enddo
   enddo

   do k=ks+1,ke-1
   do j=js+1,je-1
   do i=is+1,ie-1
    YY=ax(i,j,k)
    if (i.gt.is+1) then
     YY=YY-icbx(i,j,k)*redsoln(i-1,j,k)
    endif
    if (j.gt.js+1) then
     YY=YY-icby(i,j,k)*redsoln(i,j-1,k)
    endif
    if (k.gt.ks+1) then
     YY=YY-icbz(i,j,k)*redsoln(i,j,k-1)
    endif
    redsoln(i,j,k)=YY
   enddo
   enddo
   enddo
   do k=ks+1,ke-1
   do j=js+1,je-1
   do i=is+1,ie-1
    blacksoln(i,j,k)=redsoln(i,j,k)/icdiag(i,j,k)
   enddo 
   enddo 
   enddo 
   do k=ke-1,ks+1,-1
   do j=je-1,js+1,-1
   do i=ie-1,is+1,-1
    XX=blacksoln(i,j,k)
    if (i.lt.ie-1) then
     XX=XX-icbx(i+1,j,k)*redsoln(i+1,j,k)
    endif
    if (j.lt.je-1) then
     XX=XX-icby(i,j+1,k)*redsoln(i,j+1,k)
    endif
    if (k.lt.ke-1) then
     XX=XX-icbz(i,j,k+1)*redsoln(i,j,k+1)
    endif
    redsoln(i,j,k)=XX
   enddo
   enddo
   enddo

   do i=is+1,ie-1
   do j=js+1,je-1
   do k=ks+1,ke-1
    solnsave(i,j,k)=solnsave(i,j,k)+redsoln(i,j,k)
   enddo
   enddo
   enddo

  enddo

  do i=is+1,ie-1
  do j=js+1,je-1
  do k=ks+1,ke-1
    soln(i,j,k)=soln(i,j,k)+solnsave(i,j,k)
  enddo
  enddo
  enddo

  deallocate(solnsave)
  deallocate(rhssave)
  deallocate(redsoln)
  deallocate(blacksoln)
  deallocate(ax)

  return
  end



 subroutine gsrb2d(is,ie,js,je, &
   bxleft,bxright,byleft,byright, &
   diag,diagsing,soln,rhs,mask,nsmooth)

 IMPLICIT NONE

 INTEGER_T :: i,j,is,ie,js,je,nsmooth,n
 REAL_T, dimension(is:ie,js:je) :: bxleft
 REAL_T, dimension(is:ie,js:je) :: bxright
 REAL_T, dimension(is:ie,js:je) :: byleft
 REAL_T, dimension(is:ie,js:je) :: byright
 REAL_T, dimension(is:ie,js:je) :: diag
 REAL_T, dimension(is:ie,js:je) :: diagsing
 REAL_T, dimension(is:ie,js:je) :: soln
 REAL_T, dimension(is:ie,js:je) :: rhs
 REAL_T, dimension(is:ie,js:je) :: mask
 REAL_T, allocatable, dimension(:,:) :: ax
 REAL_T, allocatable, dimension(:,:) :: solnsave
 REAL_T, allocatable, dimension(:,:) :: rhssave
 REAL_T, allocatable, dimension(:,:) :: redsoln
 REAL_T, allocatable, dimension(:,:) :: blacksoln

 allocate(ax(is:ie,js:je))
 allocate(solnsave(is:ie,js:je))
 allocate(rhssave(is:ie,js:je))
 allocate(redsoln(is:ie,js:je))
 allocate(blacksoln(is:ie,js:je))

 do i=is+1,ie-1
 do j=js+1,je-1
  rhssave(i,j)=rhs(i,j)+ &
    bxleft(i-1,j)*soln(i-1,j)+bxright(i+1,j)*soln(i+1,j)+ &
    byleft(i,j-1)*soln(i,j-1)+byright(i,j+1)*soln(i,j+1)- &
    diagsing(i,j)*soln(i,j)
 enddo
 enddo

!$acc region 

 do i=is,ie
 do j=js,je
  solnsave(i,j)=0.0
 enddo
 enddo

  do n=1,nsmooth

   do i=is+1,ie-1
   do j=js+1,je-1
    ax(i,j)=rhssave(i,j)+ &
     bxleft(i-1,j)*solnsave(i-1,j)+bxright(i+1,j)*solnsave(i+1,j)+ &
     byleft(i,j-1)*solnsave(i,j-1)+byright(i,j+1)*solnsave(i,j+1)- &
     diagsing(i,j)*solnsave(i,j)
   enddo
   enddo
   do i=is,ie
   do j=js,je
    redsoln(i,j)=0.0
    blacksoln(i,j)=0.0
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
    redsoln(i,j)=ax(i,j)/diag(i,j)
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
    blacksoln(i,j)=(ax(i,j)+ &
     bxleft(i-1,j)*redsoln(i-1,j)+bxright(i+1,j)*redsoln(i+1,j)+ &
     byleft(i,j-1)*redsoln(i,j-1)+byright(i,j+1)*redsoln(i,j+1))/ &
     diag(i,j)
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
    redsoln(i,j)=(ax(i,j)+ &
     bxleft(i-1,j)*blacksoln(i-1,j)+bxright(i+1,j)*blacksoln(i+1,j)+ &
     byleft(i,j-1)*blacksoln(i,j-1)+byright(i,j+1)*blacksoln(i,j+1))/ &
     diag(i,j)
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
    solnsave(i,j)=solnsave(i,j)+mask(i,j)*redsoln(i,j)+ &
     (1.0-mask(i,j))*blacksoln(i,j)
   enddo
   enddo

  enddo

!$acc end region

  do i=is+1,ie-1
  do j=js+1,je-1
    soln(i,j)=soln(i,j)+solnsave(i,j)
  enddo
  enddo

  deallocate(solnsave)
  deallocate(rhssave)
  deallocate(redsoln)
  deallocate(blacksoln)
  deallocate(ax)

  return
  end


 subroutine icrb2d(is,ie,js,je, &
   bxleft,bxright,byleft,byright, &
   diag,diagsing,soln,rhs,icdiag,mask,nsmooth)

 IMPLICIT NONE

 INTEGER_T :: i,j,is,ie,js,je,nsmooth,n
 REAL_T, dimension(is:ie,js:je) :: bxleft
 REAL_T, dimension(is:ie,js:je) :: bxright
 REAL_T, dimension(is:ie,js:je) :: byleft
 REAL_T, dimension(is:ie,js:je) :: byright
 REAL_T, dimension(is:ie,js:je) :: diag
 REAL_T, dimension(is:ie,js:je) :: diagsing
 REAL_T, dimension(is:ie,js:je) :: soln
 REAL_T, dimension(is:ie,js:je) :: rhs
 REAL_T, dimension(is:ie,js:je) :: icdiag
 REAL_T, dimension(is:ie,js:je) :: mask
 REAL_T, allocatable, dimension(:,:) :: ax
 REAL_T, allocatable, dimension(:,:) :: solnsave
 REAL_T, allocatable, dimension(:,:) :: rhssave
 REAL_T, allocatable, dimension(:,:) :: redsoln
 REAL_T, allocatable, dimension(:,:) :: blacksoln

 allocate(ax(is:ie,js:je))
 allocate(rhssave(is:ie,js:je))
 allocate(solnsave(is:ie,js:je))
 allocate(redsoln(is:ie,js:je))
 allocate(blacksoln(is:ie,js:je))

 do i=is+1,ie-1
 do j=js+1,je-1
  rhssave(i,j)=rhs(i,j)+ &
    bxleft(i-1,j)*soln(i-1,j)+bxright(i+1,j)*soln(i+1,j)+ &
    byleft(i,j-1)*soln(i,j-1)+byright(i,j+1)*soln(i,j+1)- &
    diagsing(i,j)*soln(i,j)
 enddo
 enddo

!$acc region 

 do i=is,ie
 do j=js,je
  solnsave(i,j)=0.0
 enddo
 enddo

  do n=1,nsmooth

   do i=is+1,ie-1
   do j=js+1,je-1
    ax(i,j)=rhssave(i,j)+ &
     bxleft(i-1,j)*solnsave(i-1,j)+bxright(i+1,j)*solnsave(i+1,j)+ &
     byleft(i,j-1)*solnsave(i,j-1)+byright(i,j+1)*solnsave(i,j+1)- &
     diagsing(i,j)*solnsave(i,j)
   enddo
   enddo
   do i=is,ie
   do j=js,je
    redsoln(i,j)=0.0
    blacksoln(i,j)=0.0
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
    redsoln(i,j)=ax(i,j)/diag(i,j)
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
    blacksoln(i,j)=(ax(i,j)+ &
     bxleft(i-1,j)*redsoln(i-1,j)+bxright(i+1,j)*redsoln(i+1,j)+ &
     byleft(i,j-1)*redsoln(i,j-1)+byright(i,j+1)*redsoln(i,j+1))/ &
     icdiag(i,j)
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
    redsoln(i,j)=(ax(i,j)+ &
     bxleft(i-1,j)*blacksoln(i-1,j)+bxright(i+1,j)*blacksoln(i+1,j)+ &
     byleft(i,j-1)*blacksoln(i,j-1)+byright(i,j+1)*blacksoln(i,j+1))/ &
     diag(i,j)
   enddo
   enddo
   do i=is+1,ie-1
   do j=js+1,je-1
    solnsave(i,j)=solnsave(i,j)+mask(i,j)*redsoln(i,j)+ &
     (1.0-mask(i,j))*blacksoln(i,j)
   enddo
   enddo

  enddo

!$acc end region

  do i=is+1,ie-1
  do j=js+1,je-1
    soln(i,j)=soln(i,j)+solnsave(i,j)
  enddo
  enddo

  deallocate(solnsave)
  deallocate(rhssave)
  deallocate(redsoln)
  deallocate(blacksoln)
  deallocate(ax)

  return
  end


 subroutine icc2d(is,ie,js,je, &
  bxleft,bxright,byleft,byright,diag,diagsing, &
  soln,rhs,icdiag,icbx,icby,nsmooth)

 IMPLICIT NONE

 INTEGER_T :: i,j,is,ie,js,je,nsmooth,n
 REAL_T, dimension(is:ie,js:je) :: bxleft
 REAL_T, dimension(is:ie,js:je) :: bxright
 REAL_T, dimension(is:ie,js:je) :: byleft
 REAL_T, dimension(is:ie,js:je) :: byright
 REAL_T, dimension(is:ie,js:je) :: diag
 REAL_T, dimension(is:ie,js:je) :: diagsing
 REAL_T, dimension(is:ie,js:je) :: soln
 REAL_T, dimension(is:ie,js:je) :: rhs
 REAL_T, dimension(is:ie,js:je) :: icdiag
 REAL_T, dimension(is:ie,js:je) :: icbx
 REAL_T, dimension(is:ie,js:je) :: icby
 REAL_T, allocatable, dimension(:,:) :: ax
 REAL_T, allocatable, dimension(:,:) :: solnsave
 REAL_T, allocatable, dimension(:,:) :: rhssave
 REAL_T, allocatable, dimension(:,:) :: redsoln
 REAL_T, allocatable, dimension(:,:) :: blacksoln
 REAL_T YY,XX

 allocate(ax(is:ie,js:je))
 allocate(solnsave(is:ie,js:je))
 allocate(rhssave(is:ie,js:je))
 allocate(redsoln(is:ie,js:je))
 allocate(blacksoln(is:ie,js:je))

 do i=is+1,ie-1
 do j=js+1,je-1
  rhssave(i,j)=rhs(i,j)+ &
    bxleft(i-1,j)*soln(i-1,j)+bxright(i+1,j)*soln(i+1,j)+ &
    byleft(i,j-1)*soln(i,j-1)+byright(i,j+1)*soln(i,j+1)- &
    diagsing(i,j)*soln(i,j)
 enddo
 enddo

 do i=is,ie
 do j=js,je
  solnsave(i,j)=0.0
 enddo
 enddo

  do n=1,nsmooth

   do i=is+1,ie-1
   do j=js+1,je-1
    ax(i,j)=rhssave(i,j)+ &
     bxleft(i-1,j)*solnsave(i-1,j)+bxright(i+1,j)*solnsave(i+1,j)+ &
     byleft(i,j-1)*solnsave(i,j-1)+byright(i,j+1)*solnsave(i,j+1)- &
     diagsing(i,j)*solnsave(i,j)
   enddo
   enddo

   do j=js+1,je-1
   do i=is+1,ie-1
    YY=ax(i,j)
    if (i.gt.is+1) then
     YY=YY-icbx(i,j)*redsoln(i-1,j)
    endif
    if (j.gt.js+1) then
     YY=YY-icby(i,j)*redsoln(i,j-1)
    endif
    redsoln(i,j)=YY
   enddo
   enddo
   do j=js+1,je-1
   do i=is+1,ie-1
    blacksoln(i,j)=redsoln(i,j)/icdiag(i,j)
   enddo 
   enddo 
   do j=je-1,js+1,-1
   do i=ie-1,is+1,-1
    XX=blacksoln(i,j)
    if (i.lt.ie-1) then
     XX=XX-icbx(i+1,j)*redsoln(i+1,j)
    endif
    if (j.lt.je-1) then
     XX=XX-icby(i,j+1)*redsoln(i,j+1)
    endif
    redsoln(i,j)=XX
   enddo
   enddo

   do i=is+1,ie-1
   do j=js+1,je-1
    solnsave(i,j)=solnsave(i,j)+redsoln(i,j)
   enddo
   enddo

  enddo

  do i=is+1,ie-1
  do j=js+1,je-1
    soln(i,j)=soln(i,j)+solnsave(i,j)
  enddo
  enddo

  deallocate(solnsave)
  deallocate(rhssave)
  deallocate(redsoln)
  deallocate(blacksoln)
  deallocate(ax)

  return
  end




