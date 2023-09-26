c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine cal_tension(dt2,free_bc,free_bcb
     &                      ,tdx0_1,tdx0_2,tdx0_3
     &                      ,tdxnm1_1,tdxnm1_2,tdxnm1_3,ibdy,ibdyg)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c solve for tension inside inextensible fibers
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      integer i,j,ibdy,flaginex,jq,ibdyg,i2g
      real temp1,temp2,temp3,temp3dt,dt2,rhs(nr_ibm_fib,ns_ibm_fib)
      real temp4_1_ibm(nr_ibm_fib,0:ns_ibm_fib)
     &    ,temp4_2_ibm(nr_ibm_fib,0:ns_ibm_fib)
     &    ,temp4_3_ibm(nr_ibm_fib,0:ns_ibm_fib)
     &      ,a(ns_ibm_fib)
     &      ,b(ns_ibm_fib),c(ns_ibm_fib)
     &      ,tdx0_1(nr_ibm_fib)
     &      ,tdx0_2(nr_ibm_fib)
     &      ,tdx0_3(nr_ibm_fib)
     &      ,tdxnm1_1(nr_ibm_fib)
     &      ,tdxnm1_2(nr_ibm_fib)
     &      ,tdxnm1_3(nr_ibm_fib)
      logical free_bc,free_bcb,implicit_tension

      implicit_tension=.true.

      temp3dt = dt2
!     here dt acually is delta_t^2 in input to subroutine
c here {dt^2} bring to t^{n+1/2} to have a good matrix shape
      do i=ibdy,ibdy
      i2g=ibdyg 
      flaginex=iflaginext(i2g)
         do j=1,ns_ibm_r_fib(i)-1
c adding (dx.dx)^{n-1} to rhs  rhs <-- \frac{(1+(dx.dx)^{n-1})}
            rhs(i,j)=(density_coef(i2g)+density_coefp(i2g))
     &           *(1.0d0+fs_1_ibmo_fib(i,j))/2.0d0
!     +density_coefp(i)/2.0d0
c computing (dx.dx)^{n}
! num_divu_outer increased ....?
! perhaps increase "maxiter"
! 
            if (implicit_tension.eqv..true.) then
             temp1=(gx_ibm_massive_fib(i,j+1)-gx_ibm_massive_fib(i,j))/
     &          dsf_ibm_fib(i,j)
             temp2=(gy_ibm_massive_fib(i,j+1)-gy_ibm_massive_fib(i,j))/
     &          dsf_ibm_fib(i,j) 
            else if (implicit_tension.eqv..false.) then
             temp1=(gx_ibm_fib(i,j+1)-gx_ibm_fib(i,j))/dsf_ibm_fib(i,j)
             temp2=(gy_ibm_fib(i,j+1)-gy_ibm_fib(i,j))/dsf_ibm_fib(i,j) 
            else
             print *,"implicit_tension corrupt logical"
             stop
            endif

            temp3=0.0
            if(ndim .gt. 2) 
     &      temp3=(gz_ibm_fib(i,j+1)-gz_ibm_fib(i,j))/dsf_ibm_fib(i,j) 
c saving computed term for using later time
            fs_1_ibm_fib(i,j)=temp1**2.+temp2**2.+temp3**2.
            if(flaginex .eq. 0) then
            fs_2_ibm_fib(i,j)=ks_ibm_r(i2g)
     &       *(fs_1_ibm_fib(i,j)**0.5-1.0)*temp3dt
            goto 900
            endif
c rhs <-- rhs - \frac{ (dx.dx)^{n} }{dt^2}
            rhs(i,j)=rhs(i,j)-(density_coef(i2g)+density_coefp(i2g))
     &      *fs_1_ibm_fib(i,j)
c rhs <-- rhs + cs_ibm(i)/2 \frac{ 1 - (dx.dx)^{n} }{dt} 
           rhs(i,j)=rhs(i,j)
     &       +cs_ibm(i2g)/2.d0*(1.d0-fs_1_ibmo_fib(i,j))
     &      *temp3dt**0.5
c rhs <---rhs - \(dv.dv)^{n}
             temp1=(vibm1_fib(i,j+1)-vibm1_fib(i,j))/dsf_ibm_fib(i,j)
             temp2=(vibm2_fib(i,j+1)-vibm2_fib(i,j))/dsf_ibm_fib(i,j)
             if(ndim .gt. 2) 
     &       temp3=(vibm3_fib(i,j+1)-vibm3_fib(i,j))/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-(density_coef(i2g)+density_coefp(i2g))
     &      *(temp1**2.+temp2**2.+temp3**2.)*temp3dt
c computing (dx)^{*}
            temp4_1_ibm(i,j)=
     &          (gx_ibm_massive_fib(i,j+1)-gx_ibm_massive_fib(i,j))
     &       /dsf_ibm_fib(i,j)
            temp4_2_ibm(i,j)=
     &          (gy_ibm_massive_fib(i,j+1)-gy_ibm_massive_fib(i,j))
     &       /dsf_ibm_fib(i,j)   
           if(ndim .gt. 2) then
            temp4_3_ibm(i,j)=
     &          (gz_ibm_massive_fib(i,j+1)-gz_ibm_massive_fib(i,j))
     &       /dsf_ibm_fib(i,j)
            else
            temp4_3_ibm(i,j)=0.d0
            endif
900     continue
        end do
        if(flaginex .ne. 0) then
   
        do j=1,ns_ibm_r_fib(i)-1
         
            if ((j .eq. 1).and. .not.(free_bc)) then
            temp2=0.0d0
            temp1=fb_1_ibm_fib(i,j+1) + f_impuls1_fib(i,j+1)
     &           +finer_1_ibm_fib(i,j+1)
     &           -f_link1_fib(i,j+1)*ds_ibm_fib(i,2)
     &              -fibm1_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_1_ibm(i,j)*temp2*temp3dt
            temp2=0.0d0
            temp1=fb_2_ibm_fib(i,j+1) + f_impuls2_fib(i,j+1)
     &           +finer_2_ibm_fib(i,j+1)
     &           -f_link2_fib(i,j+1)*ds_ibm_fib(i,2)
     &              -fibm2_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_2_ibm(i,j)*temp2*temp3dt
             if(ndim .gt. 2) then
               temp2=0.0d0
               temp1=fb_3_ibm_fib(i,j+1) + f_impuls3_fib(i,j+1)
     &           +finer_3_ibm_fib(i,j+1)
     &           -f_link3_fib(i,j+1)*ds_ibm_fib(i,2)
     &              -fibm3_fib(i,j+1)
               temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
               rhs(i,j)=rhs(i,j)-temp4_3_ibm(i,j)*temp2*temp3dt
             endif
            else if ((j .eq. 1).and. free_bc) then
            temp2=fb_1_ibm_fib(i,j)  + f_impuls1_fib(i,j)
     &        +finer_1_ibm_fib(i,j)
     &        -f_link1_fib(i,j)*ds_ibm_fib(i,1)-fibm1_fib(i,j)
            temp1=fb_1_ibm_fib(i,j+1)  + f_impuls1_fib(i,j+1)
     &           +finer_1_ibm_fib(i,j+1)
     &          -f_link1_fib(i,j+1)*ds_ibm_fib(i,2)
     &              -fibm1_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_1_ibm(i,j)*temp2*temp3dt
     &            -temp4_1_ibm(i,j)
     &           *tdx0_1(i)*temp3dt*2./dsf_ibm_fib(i,j)**2

            temp2=fb_2_ibm_fib(i,j) + f_impuls2_fib(i,j)
     &        +finer_2_ibm_fib(i,j)
     &         -f_link2_fib(i,j)*ds_ibm_fib(i,1)-fibm2_fib(i,j)
            temp1=fb_2_ibm_fib(i,j+1)  + f_impuls2_fib(i,j+1)
     &           +finer_2_ibm_fib(i,j+1)
     &          -f_link2_fib(i,j+1)*ds_ibm_fib(i,2)
     &              -fibm2_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_2_ibm(i,j)*temp2*temp3dt
     &            -temp4_2_ibm(i,j)
     &            *tdx0_2(i)*temp3dt*2./dsf_ibm_fib(i,j)**2
            if(ndim .gt. 2) then
            temp2=fb_3_ibm_fib(i,j) + f_impuls3_fib(i,j)
     &        +finer_3_ibm_fib(i,j)
     &       - f_link3_fib(i,j)*ds_ibm_fib(i,1)-fibm3_fib(i,j)
            temp1=fb_3_ibm_fib(i,j+1)  + f_impuls3_fib(i,j+1)
     &        +finer_3_ibm_fib(i,j+1)-f_link3_fib(i,j+1)*ds_ibm_fib(i,2)
     &             -fibm3_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_3_ibm(i,j)*temp2*temp3dt
     &            -temp4_3_ibm(i,j)
     &            *tdx0_3(i)*temp3dt*2./dsf_ibm_fib(i,j)**2
            endif
        elseif ((j .eq. ns_ibm_r_fib(i)-1).and. .not.(free_bcb)) then
            temp2=fb_1_ibm_fib(i,j)  + f_impuls1_fib(i,j)
     &        +finer_1_ibm_fib(i,j)
     &        -f_link1_fib(i,j)*ds_ibm_fib(i,j)-fibm1_fib(i,j)
            temp1=0.0d0
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_1_ibm(i,j)*temp2*temp3dt
            temp2=fb_2_ibm_fib(i,j) + f_impuls2_fib(i,j)
     &        +finer_2_ibm_fib(i,j)
     &         -f_link2_fib(i,j)*ds_ibm_fib(i,j)-fibm2_fib(i,j)
            temp1=0.0d0
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_2_ibm(i,j)*temp2*temp3dt
            if(ndim .gt. 2) then
            temp2=fb_3_ibm_fib(i,j) + f_impuls3_fib(i,j)
     &        +finer_3_ibm_fib(i,j)
     &       - f_link3_fib(i,j)*ds_ibm_fib(i,j)-fibm3_fib(i,j)
            temp1=0.0d0
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)

            rhs(i,j)=rhs(i,j)-temp4_3_ibm(i,j)*temp2*temp3dt
            endif
         elseif ((j .eq. ns_ibm_r_fib(i)-1).and. free_bcb) then
!     check
            temp2=fb_1_ibm_fib(i,j) + f_impuls1_fib(i,j)
     &      +finer_1_ibm_fib(i,j)
     &      -f_link1_fib(i,j)*ds_ibm_fib(i,j)-fibm1_fib(i,j)
            temp1=fb_1_ibm_fib(i,j+1)  + f_impuls1_fib(i,j+1)
     &      +finer_1_ibm_fib(i,j+1)-f_link1_fib(i,j+1)*ds_ibm_fib(i,j)
     &              -fibm1_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)
     &              -temp4_1_ibm(i,j)*temp2*temp3dt
     &          -temp4_1_ibm(i,j)
     &          *tdxnm1_1(i)*temp3dt*2./dsf_ibm_fib(i,j)**2
            temp2=fb_2_ibm_fib(i,j) + f_impuls2_fib(i,j)
     &       +finer_2_ibm_fib(i,j)
     &       -f_link2_fib(i,j)*ds_ibm_fib(i,j)-fibm2_fib(i,j)
            temp1=fb_2_ibm_fib(i,j+1) + f_impuls2_fib(i,j+1)
     &      +finer_2_ibm_fib(i,j+1)-f_link2_fib(i,j+1)*ds_ibm_fib(i,j)
     &              -fibm2_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)
     &              -temp4_2_ibm(i,j)*temp2*temp3dt
     &          -temp4_2_ibm(i,j)
     &          *tdxnm1_2(i)*temp3dt*2./dsf_ibm_fib(i,j)**2
            if(ndim .gt. 2) then
            temp2=fb_3_ibm_fib(i,j) + f_impuls3_fib(i,j)
     &       +finer_3_ibm_fib(i,j)
     &       -f_link3_fib(i,j)*ds_ibm_fib(i,j)-fibm3_fib(i,j)
            temp1=fb_3_ibm_fib(i,j+1) + f_impuls3_fib(i,j+1)
     &      +finer_3_ibm_fib(i,j+1)-f_link3_fib(i,j+1)*ds_ibm_fib(i,j)
     &              -fibm3_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)
     &              -temp4_3_ibm(i,j)*temp2*temp3dt
     &          -temp4_3_ibm(i,j)
     &          *tdxnm1_3(i)*temp3dt*2./dsf_ibm_fib(i,j)**2
            endif
            else
            temp2=fb_1_ibm_fib(i,j) + f_impuls1_fib(i,j)
     &           +finer_1_ibm_fib(i,j)-f_link1_fib(i,j)-fibm1_fib(i,j)
            temp1=fb_1_ibm_fib(i,j+1) + f_impuls1_fib(i,j+1)
     &           +finer_1_ibm_fib(i,j+1)
     &           -f_link1_fib(i,j+1)-fibm1_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_1_ibm(i,j)*temp2*temp3dt
            temp2=fb_2_ibm_fib(i,j) + f_impuls2_fib(i,j)
     &           +finer_2_ibm_fib(i,j)-f_link2_fib(i,j)-fibm2_fib(i,j)
            temp1=fb_2_ibm_fib(i,j+1) + f_impuls2_fib(i,j+1)
     &           +finer_2_ibm_fib(i,j+1)
     &           -f_link2_fib(i,j+1)-fibm2_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_2_ibm(i,j)*temp2*temp3dt
            if(ndim .gt. 2) then
            temp2=fb_3_ibm_fib(i,j) + f_impuls3_fib(i,j)
     &           +finer_3_ibm_fib(i,j)-f_link3_fib(i,j)-fibm3_fib(i,j)
            temp1=fb_3_ibm_fib(i,j+1) + f_impuls3_fib(i,j+1)
     &           +finer_3_ibm_fib(i,j+1)
     &           -f_link3_fib(i,j+1)-fibm3_fib(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm_fib(i,j)
            rhs(i,j)=rhs(i,j)-temp4_3_ibm(i,j)*temp2*temp3dt
            endif
            end if
         end do
       do j=1,ns_ibm_r_fib(i)-1
         temp2=1.0d0/dsf_ibm_fib(i,j)/ds_ibm_fib(i,j)
         temp3=1.0d0/dsf_ibm_fib(i,j)/ds_ibm_fib(i,j+1)
         if ((j .eq. 1) .and. .not.(free_bc)) then
          c(j)=0.0d0
          a(j)=-temp3*( temp4_1_ibm(i,j)**2.0
     &          +       temp4_2_ibm(i,j)**2.0
     &          +       temp4_3_ibm(i,j)**2.0 )
          b(j)=temp3*( temp4_1_ibm(i,j+1)*temp4_1_ibm(i,j)
     &          +      temp4_2_ibm(i,j+1)*temp4_2_ibm(i,j)
     &          +      temp4_3_ibm(i,j+1)*temp4_3_ibm(i,j) )
         elseif ((j .eq. 1) .and. free_bc) then
          temp1=2.0d0/dsf_ibm_fib(i,j)/dsf_ibm_fib(i,j)
          c(j)=0.0d0
          a(j)=-(temp3+1.0d0*temp1)*( temp4_1_ibm(i,j)**2.0
     &          +                temp4_2_ibm(i,j)**2.0 
     &          +                temp4_3_ibm(i,j)**2.0 )
          b(j)=temp3*( temp4_1_ibm(i,j+1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j+1)*temp4_2_ibm(i,j) 
     &          +         temp4_3_ibm(i,j+1)*temp4_3_ibm(i,j) )
         elseif ((j .eq. ns_ibm_r_fib(i)-1).and. free_bcb) then
          temp1=2.0d0/dsf_ibm_fib(i,j)/dsf_ibm_fib(i,j)
          c(j)=temp2*( temp4_1_ibm(i,j-1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j-1)*temp4_2_ibm(i,j) 
     &          +         temp4_3_ibm(i,j-1)*temp4_3_ibm(i,j) )
          a(j)=-(temp2+1.0d0*temp1)*( temp4_1_ibm(i,j)**2.0
     &          +                temp4_2_ibm(i,j)**2.0 
     &          +                temp4_3_ibm(i,j)**2.0 )
          b(j)=0.0d0
        elseif ((j .eq. ns_ibm_r_fib(i)-1).and. .not.(free_bcb)) then
          c(j)=temp2*( temp4_1_ibm(i,j-1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j-1)*temp4_2_ibm(i,j) 
     &          +         temp4_3_ibm(i,j-1)*temp4_3_ibm(i,j) )
          a(j)=-(temp2)*( temp4_1_ibm(i,j)**2.0
     &          +                temp4_2_ibm(i,j)**2.0 
     &          +                temp4_3_ibm(i,j)**2.0 )
          b(j)=0.0d0
         else 
          c(j)=temp2*( temp4_1_ibm(i,j-1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j-1)*temp4_2_ibm(i,j) 
     &          +         temp4_3_ibm(i,j-1)*temp4_3_ibm(i,j) )
          a(j)=-(temp2+temp3)*( temp4_1_ibm(i,j)**2.0
     &          +                temp4_2_ibm(i,j)**2.0 
     &          +                temp4_3_ibm(i,j)**2.0 )
          b(j)=temp3*( temp4_1_ibm(i,j+1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j+1)*temp4_2_ibm(i,j) 
     &          +         temp4_3_ibm(i,j+1)*temp4_3_ibm(i,j) )
         end if
       end do
       call thomas(ns_ibm_r_fib(i)-1,a,b,c,rhs(i,1:ns_ibm_r_fib(i)-1)
     &          ,fs_2_ibm_fib(i,1:ns_ibm_r_fib(i)-1),ns_ibm_r_fib(i)-1)
       endif
      enddo
      return
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
