c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine eval_n_dn_ddn(v,w, 
     *                           s,dsdv,dsdw,ddsdvdv,ddsdwdw,ddsdvdw, 
     *                           n, phi, l, iv, jtype
     *                           ,jtypemax,maxordern)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c evaluate shape functions 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       implicit none
       integer n , jtype
       real v, w
       integer jtypemax,maxordern
       real  l(jtypemax,maxordern)
       real   iv(jtypemax,maxordern,maxordern)
       real  phi(jtypemax,0:2,maxordern,12)
       integer npow, k, i, ii, j
       real pow2, coef, pl
       real nn(maxordern), s(maxordern), m(12)
       real dnndv(maxordern), dsdv(maxordern), dmdv(12,2) 
       real dnndw(maxordern), dsdw (maxordern)
       real ddnndvdv(maxordern), ddsdvdv(maxordern), ddmdv2(12,3) 
       real ddnndwdw(maxordern), ddsdwdw(maxordern)
       real ddnndvdw(maxordern), ddsdvdw(maxordern)
 
       if (n .ne. 0) then

       npow = int(1-log(v+w+1.0d-12)/log(2.0))
       pow2= 2.0**(npow-1)
       v = v*pow2
       w = w*pow2
       if ( v .gt. 0.5 ) then
           k=0
           v=2.0*v-1.0
           w=2.0*w
       else if ( w .gt. 0.5 ) then
           k=2
           v=2.0*v
           w=2.0*w-1.0
       else
           k=1
           v=1.0-2.0*v
           w=1.0-2.0*w
       end if
       if( k .eq. 1) then
          coef=-1.0
       else
          coef=1.0
       end if
       call shapen_dn_ddn(v,w,m,dmdv,ddmdv2)
       do i=1,n+6
          ii=i-1
          pl=l(jtype,i)**(npow-1)
          nn(i)=0.0
          dnndv(i)=0.0
          dnndw(i)=0.0
          ddnndvdv(i)=0.0
          ddnndwdw(i)=0.0
          ddnndvdw(i)=0.0
          do j=1,12
                 nn(i)=nn(i)+phi(jtype,k,i,j)*m(j)
                 dnndv(i)=dnndv(i)+phi(jtype,k,i,j)*dmdv(j,1)
                 dnndw(i)=dnndw(i)+phi(jtype,k,i,j)*dmdv(j,2)
                 ddnndvdv(i)=ddnndvdv(i)+phi(jtype,k,i,j)*ddmdv2(j,1)
                 ddnndwdw(i)=ddnndwdw(i)+phi(jtype,k,i,j)*ddmdv2(j,2)
                 ddnndvdw(i)=ddnndvdw(i)+phi(jtype,k,i,j)*ddmdv2(j,3)
          end do
	  nn(i)=pl*nn(i)	
	  dnndv(i)=pl*dnndv(i)*(coef*2.0)**npow	
	  dnndw(i)=pl*dnndw(i)*(coef*2.0)**npow	
	  ddnndvdv(i)=pl*ddnndvdv(i)*(coef*2.0)**(2*npow)
	  ddnndwdw(i)=pl*ddnndwdw(i)*(coef*2.0)**(2*npow)
	  ddnndvdw(i)=pl*ddnndvdw(i)*(coef*2.0)**(2*npow)
       end do
       do i=1,n+6
          s(i)=0.0
          dsdv(i)=0.0
          dsdw(i)=0.0
          ddsdvdv(i)=0.0
          ddsdwdw(i)=0.0
          ddsdvdw(i)=0.0
           do j=1,n+6
              s(i)   =s(i)   + nn(j)   *iv(jtype,j,i)       
              dsdv(i)=dsdv(i)+ dnndv(j)*iv(jtype,j,i)       
              dsdw(i)=dsdw(i)+ dnndw(j)*iv(jtype,j,i)       
              ddsdvdv(i)=ddsdvdv(i)+ ddnndvdv(j)*iv(jtype,j,i)       
              ddsdwdw(i)=ddsdwdw(i)+ ddnndwdw(j)*iv(jtype,j,i)       
              ddsdvdw(i)=ddsdvdw(i)+ ddnndvdw(j)*iv(jtype,j,i)       
           end do
        end do  
        
        else
          call shapen_dn_ddn(v,w,m,dmdv,ddmdv2)
          do j=1,12
                 s(j)=m(j)
                 dsdv(j)= dmdv(j,1)
                 dsdw(j)= dmdv(j,2)
                 ddsdvdv(j)= ddmdv2(j,1)
                 ddsdwdw(j)= ddmdv2(j,2)
                 ddsdvdw(j)= ddmdv2(j,3)
          end do
         end if


        return
        end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
