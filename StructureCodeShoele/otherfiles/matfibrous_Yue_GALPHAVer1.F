c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine  matfibrous_Yue_GALPHA(
     &         ntermiMAXrow,ntermiMAX
     &        ,iacsrfibershell,jacsrfibershell,aacsrfibershell
     &        ,nq_ibmin,ns_ibmin
     &        ,BoundaryFibShellin
     &        ,Sig11  
     &        ,Sig12 
     &        ,Sig21 
     &        ,Sig22           
     &        ,Siga_d 
     &        ,gam11  
     &        ,gam12
     &        ,gam21
     &        ,gam22
     &        ,gama_d 
     &        ,dh1  
     &        ,ds1  
     &        ,dg1 
     &        ,dk1    
     &        ,dh2  
     &        ,ds2  
     &        ,dg2 
     &        ,dk2
     &        ,SigB
     &        ,SigB_d
     &        ,e4coefint 
     &        ,iflagsimulation    
     &        ,GX   
     &        ,dGxdt
     &        ,nq_ibmBin,ns_ibmBin
     &        ,aamat_in
     &        ,ForceOut
     &        ,MaterialTypeIBM
     &        ,alpha_coefK
     &        ,scal_K
     &        ,FlagForceCalculation 
     &        )
c-----------------------------------------------------------------------
      implicit none
      integer nq_ibmin,ns_ibmin,ntermiMAXrow,ntermiMAX
     &       ,nq_ibmBin,ns_ibmBin
      integer iacsrfibershell(ntermiMAXrow+1)
     &       ,jacsrfibershell(ntermiMAX)
      real*8  aacsrfibershell(ntermiMAX)
     &        ,ForceOut(nq_ibmin*ns_ibmin*3)
     &        ,aamat_in(nq_ibmin*ns_ibmin,3,3),aamat_intmp(3,3)
     &        ,alpha_coefK,scal_K
     &        ,GX(1:3,1-nq_ibmBin:nq_ibmin+nq_ibmBin
     &               ,1-ns_ibmBin:ns_ibmin+ns_ibmBin) 
     &        ,dGXdt(1:3,1-nq_ibmBin:nq_ibmin+nq_ibmBin
     &               ,1-ns_ibmBin:ns_ibmin+ns_ibmBin) 
      integer BoundaryFibShellin(2,2),MaterialTypeIBM
      real*8   Sig11(nq_ibmin,ns_ibmin)
     &        ,Sig12(nq_ibmin,ns_ibmin) 
     &        ,Sig21(nq_ibmin,ns_ibmin) 
     &        ,Sig22(nq_ibmin,ns_ibmin) 
     &        ,Siga_d(nq_ibmin,ns_ibmin,5)              
     &        ,gam11(nq_ibmin,ns_ibmin)  
     &        ,gam12(nq_ibmin,ns_ibmin)
     &        ,gam21(nq_ibmin,ns_ibmin)
     &        ,gam22(nq_ibmin,ns_ibmin)
     &        ,gama_d(nq_ibmin,ns_ibmin,3)
     &        ,dh1(nq_ibmin,ns_ibmin)
     &        ,ds1(nq_ibmin,ns_ibmin)  
     &        ,dg1(nq_ibmin,ns_ibmin) 
     &        ,dk1(nq_ibmin,ns_ibmin)    
     &        ,dh2(nq_ibmin,ns_ibmin)  
     &        ,ds2(nq_ibmin,ns_ibmin)  
     &        ,dg2(nq_ibmin,ns_ibmin) 
     &        ,dk2(nq_ibmin,ns_ibmin)

      real*8 SigB(4,nq_ibmin+ns_ibmin,3)
     &      ,SigB_d(4,nq_ibmin+ns_ibmin,3)
     &      ,e4coefint,coefX00

      real*8 DX_i_j(3),DX_ip_j(3),DX_im_j(3)
     &      ,DY_i_j(3),DY_ip_j(3),DY_im_j(3)
     &      ,DX_i_jp(3),DX_ip_jp(3),DX_im_jp(3)
     &      ,DY_i_jp(3),DY_ip_jp(3),DY_im_jp(3)
     &      ,DX_i_jm(3),DX_ip_jm(3),DX_im_jm(3)
     &      ,DY_i_jm(3),DY_ip_jm(3),DY_im_jm(3)
     &      ,DXX_i_j(3),DYY_i_j(3),DXY_i_j(3),DYX_i_j(3)
     
      integer i,j,im1,im2,im3,im4,jm1,jm2,jm3,jm4
     &  ,ihm1,ihm2,ihm3,ihm4,jhm1,jhm2,jhm3,jhm4,i4
      integer ip1,ip2,ip3,ip4,jp1,jp2,jp3,jp4
     &  ,ihp1,ihp2,ihp3,ihp4,jhp1,jhp2,jhp3,jhp4,k
      real*8 X(-4:4,-4:4),coef,X00,Y(3,3,-4:4,-4:4),temp_ibm3
     &  ,t(220),e5c,e6e,P(0:2,-1:1),Q(0:2,-1:1),A0(10)
      integer iAdd(-4:4,-4:4,1:3),iAddp(-4:4,-4:4,1:3)
     &  ,iiAdd(-4:4),jiAdd(-4:4),i2,j2,ntemi,noj,iA0(10,2)
     &  ,ist,ien,jst,jen,iflagsimulation,nsnq,i3,j3,nojtmp
      integer FlagForceCalculation(nq_ibmin,ns_ibmin)
      integer i1s(2),i2s(2),iu,iv,iw,ju,jv,jw

      integer formulationFlaga,formulationFlagb
      parameter(formulationFlagb=0)
!! ================== Inside i,j  2<i<nx-1; 2<j<ny-1

!!      print*, 'BoundaryFibShellin(1,1),BoundaryFibShellin(1,2)'
!!      print*, BoundaryFibShellin(1,1),BoundaryFibShellin(1,2)
!!      print*, 'BoundaryFibShellin(2,1),BoundaryFibShellin(2,2)'
!!      print*, BoundaryFibShellin(2,1),BoundaryFibShellin(2,2)

      ForceOut=0.0
      if((e4coefint .gt. 50.0) .and. (e4coefint .lt. 100.0)) 
     &   e4coefint=e4coefint-50.0

      e5c=2.0-e4coefint 
      e6e=e4coefint
      nsnq=ns_ibmin*nq_ibmin  !total grid points
!      dt=sqrt(dt2)
      iacsrfibershell=0
      jacsrfibershell=0
      aacsrfibershell=0.0

      ntemi=0;noj=0
      do i=1,nq_ibmin
 
      do j=1,ns_ibmin

          k=(i-1)*ns_ibmin
          aamat_intmp=aamat_in(k+j,1:3,1:3)

         ist=1
         ien=nq_ibmin

         jst=1
         jen=ns_ibmin


         if(BoundaryFibShellin(1,1) .eq. 3) ist=1
         if(BoundaryFibShellin(1,2) .eq. 3) ien=nq_ibmin
            


         if(BoundaryFibShellin(2,1) .eq. 3) jst=1
         if(BoundaryFibShellin(2,2) .eq. 3) jen=ns_ibmin

         FlagForceCalculation(i,j)=0

         im1=mod(i-1+nq_ibmin,nq_ibmin)
         if(im1 .eq. 0) im1=nq_ibmin
         im2=mod(i-2+nq_ibmin,nq_ibmin)
         if(im2 .eq. 0) im2=nq_ibmin
         im3=mod(i-3+nq_ibmin,nq_ibmin)
         if(im3 .eq. 0) im3=nq_ibmin
         im4=mod(i-4+nq_ibmin,nq_ibmin)
         if(im4 .eq. 0) im4=nq_ibmin

         ip1=mod(i+1+nq_ibmin,nq_ibmin)
         if(ip1 .eq. 0) ip1=nq_ibmin
         ip2=mod(i+2+nq_ibmin,nq_ibmin)
         if(ip2 .eq. 0) ip2=nq_ibmin
         ip3=mod(i+3+nq_ibmin,nq_ibmin)
         if(ip3 .eq. 0) ip3=nq_ibmin
         ip4=mod(i+4+nq_ibmin,nq_ibmin)
         if(ip4 .eq. 0) ip4=nq_ibmin

         ihm1=im1
         ihm2=im2
         ihm3=im3
         ihm4=im4

         ihp1=ip1
         ihp2=ip2
         ihp3=ip3
         ihp4=ip4

         jm1=mod(j-1+ns_ibmin,ns_ibmin)
         if(jm1 .eq. 0) jm1=ns_ibmin
         jm2=mod(j-2+ns_ibmin,ns_ibmin)
         if(jm2 .eq. 0) jm2=ns_ibmin
         jm3=mod(j-3+ns_ibmin,ns_ibmin)
         if(jm3 .eq. 0) jm3=ns_ibmin
         jm4=mod(j-4+ns_ibmin,ns_ibmin)
         if(jm4 .eq. 0) jm4=ns_ibmin

         jp1=mod(j+1+ns_ibmin,ns_ibmin)
         if(jp1 .eq. 0) jp1=ns_ibmin
         jp2=mod(j+2+ns_ibmin,ns_ibmin)
         if(jp2 .eq. 0) jp2=ns_ibmin
         jp3=mod(j+3+ns_ibmin,ns_ibmin)
         if(jp3 .eq. 0) jp3=ns_ibmin
         jp4=mod(j+4+ns_ibmin,ns_ibmin)
         if(jp4 .eq. 0) jp4=ns_ibmin

         jhm1=jm1
         jhm2=jm2
         jhm3=jm3
         jhm4=jm4
      
         jhp1=jp1
         jhp2=jp2
         jhp3=jp3
         jhp4=jp4

         
         iiAdd(-4:4)=(/im4,im3,im2,im1,i,ip1,ip2,ip3,ip4/)
         jiAdd(-4:4)=(/jm4,jm3,jm2,jm1,j,jp1,jp2,jp3,jp4/)

        do i2=-4,4
         do j2=-4,4
             do j3=1,3
                    iAdd(i2,j2,j3)
     &                 =(iiAdd(i2)-1)*ns_ibmin*3+(jiAdd(j2)-1)*3+j3
             enddo
           enddo  
         enddo        

!         pause
         X=0.d0
         X00=0.0

       if((j .le. 2).and. (BoundaryFibShellin(1,1) .eq. 0))then
          X(-4:4,-4:4)=0.0;X00=1.0d0;coefX00=1.0d0
          FlagForceCalculation(i,j)=1

       elseif((j .ge. ns_ibmin-1)
     &       .and. (BoundaryFibShellin(1,2) .eq. 0))then
          X(-4:4,-4:4)=0.0;X00=1.0d0;coefX00=1.0d0
          FlagForceCalculation(i,j)=1

       elseif((i .le. 2).and. (BoundaryFibShellin(2,1) .eq. 0))then
          X(-4:4,-4:4)=0.0;X00=1.0d0;coefX00=1.0d0
          FlagForceCalculation(i,j)=1

       elseif((i .ge. nq_ibmin-1)
     &       .and. (BoundaryFibShellin(2,2) .eq. 0))then
          X(-4:4,-4:4)=0.0;X00=1.0d0;coefX00=1.0d0
          FlagForceCalculation(i,j)=1


       elseif((j .le. 1).and. (BoundaryFibShellin(1,1) .eq. 2))then
          X(-4:4,-4:4)=0.0;X00=1.0d0;coefX00=1.0d0
          FlagForceCalculation(i,j)=1

       elseif((j .ge. ns_ibmin)
     &       .and. (BoundaryFibShellin(1,2) .eq. 2))then
          X(-4:4,-4:4)=0.0;X00=1.0d0;coefX00=1.0d0
          FlagForceCalculation(i,j)=1

       elseif((i .le. 1).and. (BoundaryFibShellin(2,1) .eq. 2))then
          X(-4:4,-4:4)=0.0;X00=1.0d0;coefX00=1.0d0
          FlagForceCalculation(i,j)=1

       elseif((i .ge. nq_ibmin)
     &       .and. (BoundaryFibShellin(2,2) .eq. 2))then
          X(-4:4,-4:4)=0.0;X00=1.0d0;coefX00=1.0d0
          FlagForceCalculation(i,j)=1

!! ================== BNDY 1,1
       elseif(((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   .and.
     &    ((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   )then
!= ===============  1,1 VER2
!CCCCCCCC+X-xCCCCCCCC

!= ===============  1,1 VER2
!CCCCCCCC+X-xCCCCCCCC

      X(0,0)=((1.5/(dh2(i,j)*(1.5*ds1(i,j) - 0.5*ds1(ip1,j))) + 1.5/(d
     &h1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))))*(6.0*gam21(i,j) - 4.0*
     &e4coefint*gam21(i,j)))/(dh1(i,j)*dh2(i,j)) - (2.0*sigB(1,j,2))/(
     &dh2(i,j)*ds2(i,j)) - (1.5*sigB(3,i,3))/(dh1(i,j)*(1.5*ds2(i,j) -
     & 0.5*ds2(i,jp1))) - (1.5*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,j) - 
     &0.5*ds1(ip1,j))) - (2.0*sigB(3,i,1))/(dh1(i,j)*ds1(i,j)) + gam11
     &(ip1,j)/(dh1(i,j)*dh1(ip1,j)**2*ds1(i,j)) + gam22(i,jp1)/(dh2(i,
     &j)*dh2(i,jp1)**2*ds2(i,j)) - (1.0*e4coefint*gam11(i,j)*((1.5*(2.
     &0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*d
     &s2(i,jp1)) - 2.0))/(dh1(i,j)**2*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))*
     &*2) - (1.0*e4coefint*gam22(i,j)*((1.5*(2.0*ds1(i,j) - 3.0*ds1(ip
     &1,j) + ds1(ip2,j)))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - 2.0))/(dh2
     &(i,j)**2*(1.5*ds1(i,j) - 0.5*ds1(ip1,j))**2)

      X(0,+1)=(2.0*sigB(1,j,2))/(dh2(i,j)*ds2(i,j)) + (2.0*sigB(3,i,3)
     &)/(dh1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))) - (1.5*sigB(1,j,3))
     &/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) - (gam22(i,jp1)*
     &((1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1) + 1.0/(dh2(i,jp2)*ds
     &2(i,jp1))))/(dh2(i,j)*dh2(i,jp1)) - (1.0*gam22(i,j)*((1.5*(2.0*d
     &s2(i,jp1) - 3.0*ds2(i,jp2) + ds2(i,jp3)))/(1.5*ds2(i,jp1) - 0.5*
     &ds2(i,jp2)) - 2.0))/(dh2(i,j)**2*(1.5*ds2(i,jp1) - 0.5*ds2(i,jp2
     &))**2) - ((1.5/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) + 
     &2.0/(dh1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))))*(6.0*gam21(i,j) 
     &- 4.0*e4coefint*gam21(i,j)))/(dh1(i,j)*dh2(i,j)) - (1.0*gam21(i,
     &jp1)*((1.5*(2.0*ds1(i,jp1) - 3.0*ds1(ip1,jp1) + ds1(ip2,jp1)))/(
     &1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1)) - 2.0))/(dh2(i,j)*dh2(i,jp1)*
     &(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))**2) + (1.0*e4coefint*gam11(i
     &,j)*((2.0*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2
     &(i,j) - 0.5*ds2(i,jp1)) - 5.0))/(dh1(i,j)**2*(1.5*ds2(i,j) - 0.5
     &*ds2(i,jp1))**2)

      X(0,+2)=(gam22(i,jp1)*((1.0*(1.0/ds2(i,jp1) + 1.0/ds2(i,jp2)))/d
     &h2(i,jp2) + 1.0/(dh2(i,jp1)*ds2(i,jp1))))/(dh2(i,j)*dh2(i,jp1)) 
     &- (0.5*sigB(3,i,3))/(dh1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))) +
     & (1.0*gam22(i,j)*((2.0*(2.0*ds2(i,jp1) - 3.0*ds2(i,jp2) + ds2(i,
     &jp3)))/(1.5*ds2(i,jp1) - 0.5*ds2(i,jp2)) - 5.0))/(dh2(i,j)**2*(1
     &.5*ds2(i,jp1) - 0.5*ds2(i,jp2))**2) + (0.5*(6.0*gam21(i,j) - 4.0
     &*e4coefint*gam21(i,j)))/(dh1(i,j)**2*dh2(i,j)*(1.5*ds2(i,j) - 0.
     &5*ds2(i,jp1))) + (1.0*gam21(i,jp1)*((1.5*(2.0*ds1(i,jp2) - 3.0*d
     &s1(ip1,jp2) + ds1(ip2,jp2)))/(1.5*ds1(i,jp2) - 0.5*ds1(ip1,jp2))
     & - 2.0))/(dh2(i,j)*dh2(i,jp1)*(1.5*ds1(i,jp2) - 0.5*ds1(ip1,jp2)
     &)**2) - (1.0*e4coefint*gam11(i,j)*((0.5*(2.0*ds2(i,j) - 3.0*ds2(
     &i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 4.0))/(d
     &h1(i,j)**2*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**2)

      X(0,+3)=- (1.0*gam22(i,j)*((0.5*(2.0*ds2(i,jp1) - 3.0*ds2(i,jp2)
     & + ds2(i,jp3)))/(1.5*ds2(i,jp1) - 0.5*ds2(i,jp2)) - 4.0))/(dh2(i
     &,j)**2*(1.5*ds2(i,jp1) - 0.5*ds2(i,jp2))**2) - (1.0*e4coefint*ga
     &m11(i,j))/(dh1(i,j)**2*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**2) - (1.
     &0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*dh2(i,jp2)*ds2(i,jp2))

      X(0,+4)=-(1.0*gam22(i,j))/(dh2(i,j)**2*(1.5*ds2(i,jp1) - 0.5*ds2
     &(i,jp2))**2)

      X(+1,0)=(2.0*sigB(3,i,1))/(dh1(i,j)*ds1(i,j)) + (2.0*sigB(1,j,3)
     &)/(dh2(i,j)*(1.5*ds1(i,j) - 0.5*ds1(ip1,j))) - (1.5*sigB(3,i,3))
     &/(dh1(i,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (gam11(ip1,j)*
     &((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + 1.0/(dh1(ip2,j)*ds
     &1(ip1,j))))/(dh1(i,j)*dh1(ip1,j)) - (1.0*gam11(i,j)*((1.5*(2.0*d
     &s1(ip1,j) - 3.0*ds1(ip2,j) + ds1(ip3,j)))/(1.5*ds1(ip1,j) - 0.5*
     &ds1(ip2,j)) - 2.0))/(dh1(i,j)**2*(1.5*ds1(ip1,j) - 0.5*ds1(ip2,j
     &))**2) - ((2.0/(dh2(i,j)*(1.5*ds1(i,j) - 0.5*ds1(ip1,j))) + 1.5/
     &(dh1(i,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))))*(6.0*gam21(i,j) 
     &- 4.0*e4coefint*gam21(i,j)))/(dh1(i,j)*dh2(i,j)) - (1.0*gam21(ip
     &1,j)*((1.5*(2.0*ds2(ip1,j) - 3.0*ds2(ip1,jp1) + ds2(ip1,jp2)))/(
     &1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1)) - 2.0))/(dh1(i,j)*dh1(ip1,j)*
     &(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))**2) + (1.0*e4coefint*gam22(i
     &,j)*((2.0*(2.0*ds1(i,j) - 3.0*ds1(ip1,j) + ds1(ip2,j)))/(1.5*ds1
     &(i,j) - 0.5*ds1(ip1,j)) - 5.0))/(dh2(i,j)**2*(1.5*ds1(i,j) - 0.5
     &*ds1(ip1,j))**2)

      X(+1,+1)=(2.0*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1(i
     &p1,jp1))) + (2.0*sigB(3,i,3))/(dh1(i,j)*(1.5*ds2(ip1,j) - 0.5*ds
     &2(ip1,jp1))) + ((6.0*gam21(i,j) - 4.0*e4coefint*gam21(i,j))*(2.0
     &/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) + 2.0/(dh1(i,j)*
     &(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1)))))/(dh1(i,j)*dh2(i,j)) + (1.
     &0*gam21(i,jp1)*((2.0*(2.0*ds1(i,jp1) - 3.0*ds1(ip1,jp1) + ds1(ip
     &2,jp1)))/(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1)) - 5.0))/(dh2(i,j)*d
     &h2(i,jp1)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))**2) + (1.0*gam21(i
     &p1,j)*((2.0*(2.0*ds2(ip1,j) - 3.0*ds2(ip1,jp1) + ds2(ip1,jp2)))/
     &(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1)) - 5.0))/(dh1(i,j)*dh1(ip1,j)
     &*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))**2)

      X(+1,+2)=- (0.5*sigB(3,i,3))/(dh1(i,j)*(1.5*ds2(ip1,j) - 0.5*ds2
     &(ip1,jp1))) - (0.5*(6.0*gam21(i,j) - 4.0*e4coefint*gam21(i,j)))/
     &(dh1(i,j)**2*dh2(i,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (1.
     &0*gam21(ip1,j)*((0.5*(2.0*ds2(ip1,j) - 3.0*ds2(ip1,jp1) + ds2(ip
     &1,jp2)))/(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1)) - 4.0))/(dh1(i,j)*d
     &h1(ip1,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))**2) - (1.0*gam21(i
     &,jp1)*((2.0*(2.0*ds1(i,jp2) - 3.0*ds1(ip1,jp2) + ds1(ip2,jp2)))/
     &(1.5*ds1(i,jp2) - 0.5*ds1(ip1,jp2)) - 5.0))/(dh2(i,j)*dh2(i,jp1)
     &*(1.5*ds1(i,jp2) - 0.5*ds1(ip1,jp2))**2)

      X(+1,+3)=-(1.0*gam21(ip1,j))/(dh1(i,j)*dh1(ip1,j)*(1.5*ds2(ip1,j
     &) - 0.5*ds2(ip1,jp1))**2)

      X(+2,0)=(gam11(ip1,j)*((1.0*(1.0/ds1(ip1,j) + 1.0/ds1(ip2,j)))/d
     &h1(ip2,j) + 1.0/(dh1(ip1,j)*ds1(ip1,j))))/(dh1(i,j)*dh1(ip1,j)) 
     &- (0.5*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,j) - 0.5*ds1(ip1,j))) +
     & (1.0*gam11(i,j)*((2.0*(2.0*ds1(ip1,j) - 3.0*ds1(ip2,j) + ds1(ip
     &3,j)))/(1.5*ds1(ip1,j) - 0.5*ds1(ip2,j)) - 5.0))/(dh1(i,j)**2*(1
     &.5*ds1(ip1,j) - 0.5*ds1(ip2,j))**2) + (0.5*(6.0*gam21(i,j) - 4.0
     &*e4coefint*gam21(i,j)))/(dh1(i,j)*dh2(i,j)**2*(1.5*ds1(i,j) - 0.
     &5*ds1(ip1,j))) + (1.0*gam21(ip1,j)*((1.5*(2.0*ds2(ip2,j) - 3.0*d
     &s2(ip2,jp1) + ds2(ip2,jp2)))/(1.5*ds2(ip2,j) - 0.5*ds2(ip2,jp1))
     & - 2.0))/(dh1(i,j)*dh1(ip1,j)*(1.5*ds2(ip2,j) - 0.5*ds2(ip2,jp1)
     &)**2) - (1.0*e4coefint*gam22(i,j)*((0.5*(2.0*ds1(i,j) - 3.0*ds1(
     &ip1,j) + ds1(ip2,j)))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - 4.0))/(d
     &h2(i,j)**2*(1.5*ds1(i,j) - 0.5*ds1(ip1,j))**2)

      X(+2,+1)=- (0.5*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1
     &(ip1,jp1))) - (0.5*(6.0*gam21(i,j) - 4.0*e4coefint*gam21(i,j)))/
     &(dh1(i,j)*dh2(i,j)**2*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) - (1.
     &0*gam21(i,jp1)*((0.5*(2.0*ds1(i,jp1) - 3.0*ds1(ip1,jp1) + ds1(ip
     &2,jp1)))/(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1)) - 4.0))/(dh2(i,j)*d
     &h2(i,jp1)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))**2) - (1.0*gam21(i
     &p1,j)*((2.0*(2.0*ds2(ip2,j) - 3.0*ds2(ip2,jp1) + ds2(ip2,jp2)))/
     &(1.5*ds2(ip2,j) - 0.5*ds2(ip2,jp1)) - 5.0))/(dh1(i,j)*dh1(ip1,j)
     &*(1.5*ds2(ip2,j) - 0.5*ds2(ip2,jp1))**2)

      X(+2,+2)=(1.0*gam21(i,jp1)*((0.5*(2.0*ds1(i,jp2) - 3.0*ds1(ip1,j
     &p2) + ds1(ip2,jp2)))/(1.5*ds1(i,jp2) - 0.5*ds1(ip1,jp2)) - 4.0))
     &/(dh2(i,j)*dh2(i,jp1)*(1.5*ds1(i,jp2) - 0.5*ds1(ip1,jp2))**2) + 
     &(1.0*gam21(ip1,j)*((0.5*(2.0*ds2(ip2,j) - 3.0*ds2(ip2,jp1) + ds2
     &(ip2,jp2)))/(1.5*ds2(ip2,j) - 0.5*ds2(ip2,jp1)) - 4.0))/(dh1(i,j
     &)*dh1(ip1,j)*(1.5*ds2(ip2,j) - 0.5*ds2(ip2,jp1))**2)

      X(+2,+3)=(1.0*gam21(ip1,j))/(dh1(i,j)*dh1(ip1,j)*(1.5*ds2(ip2,j)
     & - 0.5*ds2(ip2,jp1))**2)

      X(+3,0)=- (1.0*gam11(i,j)*((0.5*(2.0*ds1(ip1,j) - 3.0*ds1(ip2,j)
     & + ds1(ip3,j)))/(1.5*ds1(ip1,j) - 0.5*ds1(ip2,j)) - 4.0))/(dh1(i
     &,j)**2*(1.5*ds1(ip1,j) - 0.5*ds1(ip2,j))**2) - (1.0*e4coefint*ga
     &m22(i,j))/(dh2(i,j)**2*(1.5*ds1(i,j) - 0.5*ds1(ip1,j))**2) - (1.
     &0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh1(ip2,j)*ds1(ip2,j))

      X(+3,+1)=-(1.0*gam21(i,jp1))/(dh2(i,j)*dh2(i,jp1)*(1.5*ds1(i,jp1
     &) - 0.5*ds1(ip1,jp1))**2)

      X(+3,+2)=(1.0*gam21(i,jp1))/(dh2(i,j)*dh2(i,jp1)*(1.5*ds1(i,jp2)
     & - 0.5*ds1(ip1,jp2))**2)

      X(+4,0)=-(1.0*gam11(i,j))/(dh1(i,j)**2*(1.5*ds1(ip1,j) - 0.5*ds1
     &(ip2,j))**2)


      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY Nx,1
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   )then
      coefX00=1.0d0
 
      X(-2,0)=
     &0.375/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(1.5*ds2(ihm2,j) 
     &- 0.5*ds2(ihm2,jhp1))) + 0.375/((0.5*ds1(ihm2,j) - 1.5*ds1(i
     &hm1,j))*(1.5*ds2(i,j) - 0.5*ds2(i,jhp1)))

      X(-2,+1)=
     &- 0.5/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(1.5*ds2(ihm2,j) 
     &- 0.5*ds2(ihm2,jhp1))) - 0.5/((0.5*ds1(ihm2,jhp1) - 1.5*ds1(
     &ihm1,jhp1))*(1.5*ds2(i,j) - 0.5*ds2(i,jhp1)))

      X(-2,+2)=
     &0.125/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(1.5*ds2(ihm2,j) 
     &- 0.5*ds2(ihm2,jhp1))) + 0.125/((0.5*ds1(ihm2,jhp2) - 1.5*ds
     &1(ihm1,jhp2))*(1.5*ds2(i,j) - 0.5*ds2(i,jhp1)))

      X(-1,0)=
     &- 1.5/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(1.5*ds2(ihm1,j) 
     &- 0.5*ds2(ihm1,jhp1))) - 1.5/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm
     &1,j))*(1.5*ds2(i,j) - 0.5*ds2(i,jhp1)))

      X(-1,+1)=
     &2.0/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(1.5*ds2(ihm1,j) - 
     &0.5*ds2(ihm1,jhp1))) + 2.0/((0.5*ds1(ihm2,jhp1) - 1.5*ds1(ih
     &m1,jhp1))*(1.5*ds2(i,j) - 0.5*ds2(i,jhp1)))

      X(-1,+2)=
     &- 0.5/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(1.5*ds2(ihm1,j) 
     &- 0.5*ds2(ihm1,jhp1))) - 0.5/((0.5*ds1(ihm2,jhp2) - 1.5*ds1(
     &ihm1,jhp2))*(1.5*ds2(i,j) - 0.5*ds2(i,jhp1)))

      X(0,0)=
     &2.25/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(1.5*ds2(i,j) - 0.
     &5*ds2(i,jhp1)))

      X(0,+1)=
     &- 1.5/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(1.5*ds2(i,j) - 0
     &.5*ds2(i,jhp1))) - 1.5/((0.5*ds1(ihm2,jhp1) - 1.5*ds1(ihm1,j
     &hp1))*(1.5*ds2(i,j) - 0.5*ds2(i,jhp1)))

      X(0,+2)=
     &0.375/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(1.5*ds2(i,j) - 0
     &.5*ds2(i,jhp1))) + 0.375/((0.5*ds1(ihm2,jhp2) - 1.5*ds1(ihm1
     &,jhp2))*(1.5*ds2(i,j) - 0.5*ds2(i,jhp1)))

!! ================== BNDY 1,Ny
       elseif(((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
      coefX00=1.0d0
 
      X(0,-2)=
     &0.375/((1.5*ds1(i,jhm2) - 0.5*ds1(ihp1,jhm2))*(0.5*ds2(i,jhm
     &2) - 1.5*ds2(i,jhm1))) + 0.375/((1.5*ds1(i,j) - 0.5*ds1(ihp1
     &,j))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(0,-1)=
     &- 1.5/((1.5*ds1(i,jhm1) - 0.5*ds1(ihp1,jhm1))*(0.5*ds2(i,jhm
     &2) - 1.5*ds2(i,jhm1))) - 1.5/((1.5*ds1(i,j) - 0.5*ds1(ihp1,j
     &))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(0,0)=
     &2.25/((1.5*ds1(i,j) - 0.5*ds1(ihp1,j))*(0.5*ds2(i,jhm2) - 1.
     &5*ds2(i,jhm1)))

      X(+1,-2)=
     &- 0.5/((1.5*ds1(i,jhm2) - 0.5*ds1(ihp1,jhm2))*(0.5*ds2(i,jhm
     &2) - 1.5*ds2(i,jhm1))) - 0.5/((1.5*ds1(i,j) - 0.5*ds1(ihp1,j
     &))*(0.5*ds2(ihp1,jhm2) - 1.5*ds2(ihp1,jhm1)))

      X(+1,-1)=
     &2.0/((1.5*ds1(i,jhm1) - 0.5*ds1(ihp1,jhm1))*(0.5*ds2(i,jhm2)
     & - 1.5*ds2(i,jhm1))) + 2.0/((1.5*ds1(i,j) - 0.5*ds1(ihp1,j))
     &*(0.5*ds2(ihp1,jhm2) - 1.5*ds2(ihp1,jhm1)))

      X(+1,0)=
     &- 1.5/((1.5*ds1(i,j) - 0.5*ds1(ihp1,j))*(0.5*ds2(i,jhm2) - 1
     &.5*ds2(i,jhm1))) - 1.5/((1.5*ds1(i,j) - 0.5*ds1(ihp1,j))*(0.
     &5*ds2(ihp1,jhm2) - 1.5*ds2(ihp1,jhm1)))

      X(+2,-2)=
     &0.125/((1.5*ds1(i,jhm2) - 0.5*ds1(ihp1,jhm2))*(0.5*ds2(i,jhm
     &2) - 1.5*ds2(i,jhm1))) + 0.125/((1.5*ds1(i,j) - 0.5*ds1(ihp1
     &,j))*(0.5*ds2(ihp2,jhm2) - 1.5*ds2(ihp2,jhm1)))

      X(+2,-1)=
     &- 0.5/((1.5*ds1(i,jhm1) - 0.5*ds1(ihp1,jhm1))*(0.5*ds2(i,jhm
     &2) - 1.5*ds2(i,jhm1))) - 0.5/((1.5*ds1(i,j) - 0.5*ds1(ihp1,j
     &))*(0.5*ds2(ihp2,jhm2) - 1.5*ds2(ihp2,jhm1)))

      X(+2,0)=
     &0.375/((1.5*ds1(i,j) - 0.5*ds1(ihp1,j))*(0.5*ds2(i,jhm2) - 1
     &.5*ds2(i,jhm1))) + 0.375/((1.5*ds1(i,j) - 0.5*ds1(ihp1,j))*(
     &0.5*ds2(ihp2,jhm2) - 1.5*ds2(ihp2,jhm1)))

!! ================== BNDY Nx,Ny
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
      coefX00=1.0d0
      X(-2,-2)=
     &0.125/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(0.5*ds2(ihm2,jhm
     &2) - 1.5*ds2(ihm2,jhm1))) + 0.125/((0.5*ds1(ihm2,jhm2) - 1.5
     &*ds1(ihm1,jhm2))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(-2,-1)=
     &- 0.5/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(0.5*ds2(ihm2,jhm
     &2) - 1.5*ds2(ihm2,jhm1))) - 0.5/((0.5*ds1(ihm2,jhm1) - 1.5*d
     &s1(ihm1,jhm1))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(-2,0)=
     &0.375/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(0.5*ds2(ihm2,jhm
     &2) - 1.5*ds2(ihm2,jhm1))) + 0.375/((0.5*ds1(ihm2,j) - 1.5*ds
     &1(ihm1,j))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(-1,-2)=
     &- 0.5/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(0.5*ds2(ihm1,jhm
     &2) - 1.5*ds2(ihm1,jhm1))) - 0.5/((0.5*ds1(ihm2,jhm2) - 1.5*d
     &s1(ihm1,jhm2))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(-1,-1)=
     &2.0/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(0.5*ds2(ihm1,jhm2)
     & - 1.5*ds2(ihm1,jhm1))) + 2.0/((0.5*ds1(ihm2,jhm1) - 1.5*ds1
     &(ihm1,jhm1))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(-1,0)=
     &- 1.5/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(0.5*ds2(ihm1,jhm
     &2) - 1.5*ds2(ihm1,jhm1))) - 1.5/((0.5*ds1(ihm2,j) - 1.5*ds1(
     &ihm1,j))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(0,-2)=
     &0.375/((0.5*ds1(ihm2,jhm2) - 1.5*ds1(ihm1,jhm2))*(0.5*ds2(i,
     &jhm2) - 1.5*ds2(i,jhm1))) + 0.375/((0.5*ds1(ihm2,j) - 1.5*ds
     &1(ihm1,j))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(0,-1)=
     &- 1.5/((0.5*ds1(ihm2,jhm1) - 1.5*ds1(ihm1,jhm1))*(0.5*ds2(i,
     &jhm2) - 1.5*ds2(i,jhm1))) - 1.5/((0.5*ds1(ihm2,j) - 1.5*ds1(
     &ihm1,j))*(0.5*ds2(i,jhm2) - 1.5*ds2(i,jhm1)))

      X(0,0)=
     &2.25/((0.5*ds1(ihm2,j) - 1.5*ds1(ihm1,j))*(0.5*ds2(i,jhm2) -
     & 1.5*ds2(i,jhm1)))
!! ================== BNDY 2,2
       elseif (  (
     &         (  (j.eq. 2) 
     &            .and. (    (BoundaryFibShellin(1,1) .eq. 1)
     &                   .or.(BoundaryFibShellin(1,1) .eq. 2)
     &                   ) 
     &         ) 
     &        .and.
     &        (   (i.eq. 2) 
     &            .and. (    (BoundaryFibShellin(2,1) .eq. 1)
     &                   .or.(BoundaryFibShellin(2,1) .eq. 2)
     &                   ) 
     &         ) 
     &       ) )  then
!= ===============  2,2 VER2


!CCCCCCCC+X-xCCCCCCCC

      X(-1,-1)=(0.25*sig12(im1,jm1))/(dh1(i,j)*ds2(im1,jm1)) + (0.25*s
     &ig21(im1,jm1))/(dh2(i,j)*ds1(im1,jm1)) - (1.0*gam12(im1,jm1))/(d
     &g1(im1,jm1)*dh1(i,j)*dk2(im1,j)*ds2(im1,jm1)) - (1.0*gam21(im1,j
     &m1))/(dg2(im1,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(im1,jm1)) + (e4coefin
     &t*gam11(im1,j))/(dh1(i,j)*dh2(im1,j)*ds1(im1,j)*ds2(im1,jm1)) + 
     &(e4coefint*gam22(i,jm1))/(dh1(i,jm1)*dh2(i,j)*ds1(im1,jm1)*ds2(i
     &,jm1))

      X(-1,0)=(1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*dk1(i,jm1)*ds1(
     &im1,j)) + (1.0*gam21(im1,j))/(dg2(im1,j)*dk1(i,j)*ds1(im1,j))))/
     &dh2(i,j) + ((gam11(i,j)/(dh1(i,j)*ds1(im1,j)) - (e4coefint*gam11
     &(im1,j)*(1.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/dh2(im1,j))/ds1(im1
     &,j) + gam11(i,j)/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) - ((0.
     &25*sig12(im1,jm1))/ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(im1,j)
     &)/dh1(i,j) + ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.25*sig21(im1
     &,j))/ds1(im1,j))/dh2(i,j) + sig11(im1,j)/(dh1(i,j)*ds1(im1,j)) +
     & (1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(im1,jm1)) + (1.0*g
     &am12(im1,j))/(dg1(im1,j)*ds2(im1,j))))/(dh1(i,j)*dk2(im1,j))

      X(-1,+1)=(1.0*e4coefint*gam11(im1,j))/(dh1(i,j)*dh2(im1,j)*ds1(i
     &m1,j)*ds2(im1,j)) - (0.25*sig21(im1,j))/(dh2(i,j)*ds1(im1,jp1)) 
     &- (1.0*gam12(im1,j))/(dg1(im1,j)*dh1(i,j)*dk2(im1,j)*ds2(im1,j))
     & - (1.0*gam21(im1,j))/(dg2(im1,j)*dh2(i,j)*dk1(i,j)*ds1(im1,jp1)
     &) - (0.25*sig12(im1,j))/(dh1(i,j)*ds2(im1,j))

      X(0,-1)=(1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*dk2(im1,j)*ds2(
     &i,jm1)) + (1.0*gam12(i,jm1))/(dg1(i,jm1)*dk2(i,j)*ds2(i,jm1))))/
     &dh1(i,j) + ((gam22(i,j)/(dh2(i,j)*ds2(i,jm1)) - (e4coefint*gam22
     &(i,jm1)*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/dh1(i,jm1))/ds2(i,j
     &m1) + gam22(i,j)/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + ((0.
     &25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(i,jm1))/ds2(i,jm1))/
     &dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,jm1) - (0.25*sig21(i,j
     &m1))/ds1(i,jm1))/dh2(i,j) + sig22(i,jm1)/(dh2(i,j)*ds2(i,jm1)) +
     & (1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*ds1(im1,jm1)) + (1.0*g
     &am21(i,jm1))/(dg2(i,jm1)*ds1(i,jm1))))/(dh2(i,j)*dk1(i,jm1))

      X(0,0)=- (1.0*(((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(i,jm1)) +
     & (1.0*gam12(im1,j))/(dg1(im1,j)*ds2(i,j)))/dk2(im1,j) + (1.0*((1
     &.0*gam12(i,jm1))/(dg1(i,jm1)*ds2(i,jm1)) + (1.0*gam12(i,j))/(dg1
     &(i,j)*ds2(i,j))))/dk2(i,j)))/dh1(i,j) - (1.0*(((1.0*gam21(im1,jm
     &1))/(dg2(im1,jm1)*ds1(im1,j)) + (1.0*gam21(i,jm1))/(dg2(i,jm1)*d
     &s1(i,j)))/dk1(i,jm1) + (1.0*((1.0*gam21(im1,j))/(dg2(im1,j)*ds1(
     &im1,j)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j))))/dk1(i,j)))/dh2(
     &i,j) - ((1.0*sig11(im1,j))/ds1(im1,j) + (1.0*sig11(i,j))/ds1(i,j
     &))/dh1(i,j) - ((1.0*sig22(i,jm1))/ds2(i,jm1) + (1.0*sig22(i,j))/
     &ds2(i,j))/dh2(i,j) - (((1.0*gam11(ip1,j))/(dh1(ip1,j)*ds1(i,j)) 
     &+ (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,j))/ds1(i,j
     &) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds1(i
     &m1,j)))/dh1(i,j) - (((1.0*gam22(i,jp1))/(dh2(i,jp1)*ds2(i,j)) + 
     &(gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j))/ds2(i,j) 
     &+ (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds2(i,j
     &m1)))/dh2(i,j) - ((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12
     &(im1,j))/ds2(i,j) - (0.25*sig12(i,jm1))/ds2(i,jm1) + (0.25*sig12
     &(i,j))/ds2(i,j))/dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,j) - 
     &(0.25*sig21(im1,j))/ds1(im1,j) - (0.25*sig21(i,jm1))/ds1(i,j) + 
     &(0.25*sig21(i,j))/ds1(i,j))/dh2(i,j)

      X(0,+1)=(1.0*((1.0*gam12(im1,j))/(dg1(im1,j)*dk2(im1,j)*ds2(i,j)
     &) + (1.0*gam12(i,j))/(dg1(i,j)*dk2(i,j)*ds2(i,j))))/dh1(i,j) - (
     &(0.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,j))/ds2(i,j))/dh1(i
     &,j) + ((0.25*sig21(im1,j))/ds1(im1,jp1) - (0.25*sig21(i,j))/ds1(
     &i,jp1))/dh2(i,j) + (((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0
     &*gam22(i,jp1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1))/ds2(i
     &,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) 
     &+ (1.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*((1.0*gam21(im1,j)
     &)/(dg2(im1,j)*ds1(im1,jp1)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j
     &p1))))/(dh2(i,j)*dk1(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(+1,-1)=(1.0*e4coefint*gam22(i,jm1))/(dh1(i,jm1)*dh2(i,j)*ds1(i
     &,jm1)*ds2(i,jm1)) - (0.25*sig21(i,jm1))/(dh2(i,j)*ds1(i,jm1)) - 
     &(1.0*gam12(i,jm1))/(dg1(i,jm1)*dh1(i,j)*dk2(i,j)*ds2(ip1,jm1)) -
     & (1.0*gam21(i,jm1))/(dg2(i,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(i,jm1)) 
     &- (0.25*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1))

      X(+1,0)=(1.0*((1.0*gam21(i,jm1))/(dg2(i,jm1)*dk1(i,jm1)*ds1(i,j)
     &) + (1.0*gam21(i,j))/(dg2(i,j)*dk1(i,j)*ds1(i,j))))/dh2(i,j) + (
     &(0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j))/ds2(ip1,j))
     &/dh1(i,j) - ((0.25*sig21(i,jm1))/ds1(i,j) - (0.25*sig21(i,j))/ds
     &1(i,j))/dh2(i,j) + (((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0
     &*gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/dh1(ip1,j))/ds1(i
     &,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) 
     &+ (1.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*((1.0*gam12(i,jm1)
     &)/(dg1(i,jm1)*ds2(ip1,jm1)) + (1.0*gam12(i,j))/(dg1(i,j)*ds2(ip1
     &,j))))/(dh1(i,j)*dk2(i,j))

      X(+1,+1)=(0.25*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (0.25*sig21(i
     &,j))/(dh2(i,j)*ds1(i,jp1)) - (1.0*gam12(i,j))/(dg1(i,j)*dh1(i,j)
     &*dk2(i,j)*ds2(ip1,j)) - (1.0*gam21(i,j))/(dg2(i,j)*dh2(i,j)*dk1(
     &i,j)*ds1(i,jp1))

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY Nx-1,2
           elseif(  (
     &         (  (j.eq. 2) 
     &            .and. (    (BoundaryFibShellin(1,1) .eq. 1)
     &                   .or.(BoundaryFibShellin(1,1) .eq. 2)
     &                   ) 
     &         ) 
     &        .and.
     &        (   (i.eq. nq_ibmin-1) 
     &            .and. (    (BoundaryFibShellin(2,2) .eq. 1)
     &                   .or.(BoundaryFibShellin(2,2) .eq. 2)
     &                   ) 
     &         ) 
     &       ) )  then

      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,-1)=(0.25*sig12(im1,jm1))/(dh1(i,j)*ds2(im1,jm1)) + (0.25*s
     &ig21(im1,jm1))/(dh2(i,j)*ds1(im1,jm1)) - (1.0*gam12(im1,jm1))/(d
     &g1(im1,jm1)*dh1(i,j)*dk2(im1,j)*ds2(im1,jm1)) - (1.0*gam21(im1,j
     &m1))/(dg2(im1,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(im1,jm1)) + (e4coefin
     &t*gam22(i,jm1))/(dh1(i,jm1)*dh2(i,j)*ds1(im1,jm1)*ds2(i,jm1))

      X(-1,0)=(1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*dk1(i,jm1)*ds1(
     &im1,j)) + (1.0*gam21(im1,j))/(dg2(im1,j)*dk1(i,j)*ds1(im1,j))))/
     &dh2(i,j) + ((1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)) + (gam1
     &1(im1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(im1,j)))/ds1(im1
     &,j) + gam11(i,j)/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) - ((0.
     &25*sig12(im1,jm1))/ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(im1,j)
     &)/dh1(i,j) + ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.25*sig21(im1
     &,j))/ds1(im1,j))/dh2(i,j) + sig11(im1,j)/(dh1(i,j)*ds1(im1,j)) +
     & (1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(im1,jm1)) + (1.0*g
     &am12(im1,j))/(dg1(im1,j)*ds2(im1,j))))/(dh1(i,j)*dk2(im1,j))

      X(-1,+1)=- (0.25*sig12(im1,j))/(dh1(i,j)*ds2(im1,j)) - (0.25*sig
     &21(im1,j))/(dh2(i,j)*ds1(im1,jp1)) - (1.0*gam12(im1,j))/(dg1(im1
     &,j)*dh1(i,j)*dk2(im1,j)*ds2(im1,j)) - (1.0*gam21(im1,j))/(dg2(im
     &1,j)*dh2(i,j)*dk1(i,j)*ds1(im1,jp1))

      X(0,-1)=(1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*dk2(im1,j)*ds2(
     &i,jm1)) + (1.0*gam12(i,jm1))/(dg1(i,jm1)*dk2(i,j)*ds2(i,jm1))))/
     &dh1(i,j) + ((gam22(i,j)/(dh2(i,j)*ds2(i,jm1)) - (e4coefint*gam22
     &(i,jm1)*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/dh1(i,jm1))/ds2(i,j
     &m1) + gam22(i,j)/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + ((0.
     &25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(i,jm1))/ds2(i,jm1))/
     &dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,jm1) - (0.25*sig21(i,j
     &m1))/ds1(i,jm1))/dh2(i,j) + sig22(i,jm1)/(dh2(i,j)*ds2(i,jm1)) +
     & (1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*ds1(im1,jm1)) + (1.0*g
     &am21(i,jm1))/(dg2(i,jm1)*ds1(i,jm1))))/(dh2(i,j)*dk1(i,jm1))

      X(0,0)=- (1.0*(((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(i,jm1)) +
     & (1.0*gam12(im1,j))/(dg1(im1,j)*ds2(i,j)))/dk2(im1,j) + (1.0*((1
     &.0*gam12(i,jm1))/(dg1(i,jm1)*ds2(i,jm1)) + (1.0*gam12(i,j))/(dg1
     &(i,j)*ds2(i,j))))/dk2(i,j)))/dh1(i,j) - (1.0*(((1.0*gam21(im1,jm
     &1))/(dg2(im1,jm1)*ds1(im1,j)) + (1.0*gam21(i,jm1))/(dg2(i,jm1)*d
     &s1(i,j)))/dk1(i,jm1) + (1.0*((1.0*gam21(im1,j))/(dg2(im1,j)*ds1(
     &im1,j)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j))))/dk1(i,j)))/dh2(
     &i,j) - ((1.0*sig11(im1,j))/ds1(im1,j) + (1.0*sig11(i,j))/ds1(i,j
     &))/dh1(i,j) - ((1.0*sig22(i,jm1))/ds2(i,jm1) + (1.0*sig22(i,j))/
     &ds2(i,j))/dh2(i,j) - ((1.0*((1.0*gam11(im1,j))/(dh1(im1,j)*ds1(i
     &m1,j)) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,
     &j)))/ds1(im1,j) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(
     &dh1(i,j)*ds1(i,j)))/dh1(i,j) - (((1.0*gam22(i,jp1))/(dh2(i,jp1)*
     &ds2(i,j)) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j
     &))/ds2(i,j) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(
     &i,j)*ds2(i,jm1)))/dh2(i,j) - ((0.25*sig12(im1,jm1))/ds2(i,jm1) -
     & (0.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,jm1))/ds2(i,jm1) +
     & (0.25*sig12(i,j))/ds2(i,j))/dh1(i,j) - ((0.25*sig21(im1,jm1))/d
     &s1(im1,j) - (0.25*sig21(im1,j))/ds1(im1,j) - (0.25*sig21(i,jm1))
     &/ds1(i,j) + (0.25*sig21(i,j))/ds1(i,j))/dh2(i,j)

      X(0,+1)=(1.0*((1.0*gam12(im1,j))/(dg1(im1,j)*dk2(im1,j)*ds2(i,j)
     &) + (1.0*gam12(i,j))/(dg1(i,j)*dk2(i,j)*ds2(i,j))))/dh1(i,j) - (
     &(0.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,j))/ds2(i,j))/dh1(i
     &,j) + ((0.25*sig21(im1,j))/ds1(im1,jp1) - (0.25*sig21(i,j))/ds1(
     &i,jp1))/dh2(i,j) + (((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0
     &*gam22(i,jp1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1))/ds2(i
     &,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) 
     &+ (1.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*((1.0*gam21(im1,j)
     &)/(dg2(im1,j)*ds1(im1,jp1)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j
     &p1))))/(dh2(i,j)*dk1(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(+1,-1)=(e4coefint*gam11(ip1,j))/(dh1(i,j)*dh2(ip1,j)*ds1(i,j)*
     &ds2(ip1,jm1)) - (0.25*sig21(i,jm1))/(dh2(i,j)*ds1(i,jm1)) - (1.0
     &*gam12(i,jm1))/(dg1(i,jm1)*dh1(i,j)*dk2(i,j)*ds2(ip1,jm1)) - (1.
     &0*gam21(i,jm1))/(dg2(i,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(i,jm1)) - (0
     &.25*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1)) + (1.0*e4coefint*gam22
     &(i,jm1))/(dh1(i,jm1)*dh2(i,j)*ds1(i,jm1)*ds2(i,jm1))

      X(+1,0)=(1.0*((1.0*gam21(i,jm1))/(dg2(i,jm1)*dk1(i,jm1)*ds1(i,j)
     &) + (1.0*gam21(i,j))/(dg2(i,j)*dk1(i,j)*ds1(i,j))))/dh2(i,j) + (
     &(0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j))/ds2(ip1,j))
     &/dh1(i,j) - ((0.25*sig21(i,jm1))/ds1(i,j) - (0.25*sig21(i,j))/ds
     &1(i,j))/dh2(i,j) + (((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) - (e4c
     &oefint*gam11(ip1,j)*(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/dh2(ip1
     &,j))/ds1(i,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))
     &/dh1(i,j) + (1.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*((1.0*ga
     &m12(i,jm1))/(dg1(i,jm1)*ds2(ip1,jm1)) + (1.0*gam12(i,j))/(dg1(i,
     &j)*ds2(ip1,j))))/(dh1(i,j)*dk2(i,j))

      X(+1,+1)=(0.25*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (0.25*sig21(i
     &,j))/(dh2(i,j)*ds1(i,jp1)) - (1.0*gam12(i,j))/(dg1(i,j)*dh1(i,j)
     &*dk2(i,j)*ds2(ip1,j)) - (1.0*gam21(i,j))/(dg2(i,j)*dh2(i,j)*dk1(
     &i,j)*ds1(i,jp1)) + (1.0*e4coefint*gam11(ip1,j))/(dh1(i,j)*dh2(ip
     &1,j)*ds1(i,j)*ds2(ip1,j))
      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY 2,Ny-1
           elseif(  (
     &         (  (j.eq. ns_ibmin-1) 
     &            .and. (    (BoundaryFibShellin(1,2) .eq. 1)
     &                   .or.(BoundaryFibShellin(1,2) .eq. 2)
     &                   ) 
     &         ) 
     &        .and.
     &        (   (i.eq. 2) 
     &            .and. (    (BoundaryFibShellin(2,1) .eq. 1)
     &                   .or.(BoundaryFibShellin(2,1) .eq. 2)
     &                   ) 
     &         ) 
     &       ) ) then

      X(-1,-1)=(0.25*sig12(im1,jm1))/(dh1(i,j)*ds2(im1,jm1)) + (0.25*s
     &ig21(im1,jm1))/(dh2(i,j)*ds1(im1,jm1)) - (1.0*gam12(im1,jm1))/(d
     &g1(im1,jm1)*dh1(i,j)*dk2(im1,j)*ds2(im1,jm1)) - (1.0*gam21(im1,j
     &m1))/(dg2(im1,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(im1,jm1)) + (e4coefin
     &t*gam11(im1,j))/(dh1(i,j)*dh2(im1,j)*ds1(im1,j)*ds2(im1,jm1))

      X(-1,0)=(1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*dk1(i,jm1)*ds1(
     &im1,j)) + (1.0*gam21(im1,j))/(dg2(im1,j)*dk1(i,j)*ds1(im1,j))))/
     &dh2(i,j) + ((gam11(i,j)/(dh1(i,j)*ds1(im1,j)) - (e4coefint*gam11
     &(im1,j)*(1.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/dh2(im1,j))/ds1(im1
     &,j) + gam11(i,j)/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) - ((0.
     &25*sig12(im1,jm1))/ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(im1,j)
     &)/dh1(i,j) + ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.25*sig21(im1
     &,j))/ds1(im1,j))/dh2(i,j) + sig11(im1,j)/(dh1(i,j)*ds1(im1,j)) +
     & (1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(im1,jm1)) + (1.0*g
     &am12(im1,j))/(dg1(im1,j)*ds2(im1,j))))/(dh1(i,j)*dk2(im1,j))

      X(-1,+1)=(1.0*e4coefint*gam11(im1,j))/(dh1(i,j)*dh2(im1,j)*ds1(i
     &m1,j)*ds2(im1,j)) - (0.25*sig21(im1,j))/(dh2(i,j)*ds1(im1,jp1)) 
     &- (1.0*gam12(im1,j))/(dg1(im1,j)*dh1(i,j)*dk2(im1,j)*ds2(im1,j))
     & - (1.0*gam21(im1,j))/(dg2(im1,j)*dh2(i,j)*dk1(i,j)*ds1(im1,jp1)
     &) - (0.25*sig12(im1,j))/(dh1(i,j)*ds2(im1,j)) + (e4coefint*gam22
     &(i,jp1))/(dh1(i,jp1)*dh2(i,j)*ds1(im1,jp1)*ds2(i,j))

      X(0,-2)=-(1.0*gam22(i,jm1))/(dh2(i,jm1)*dh2(i,j)*ds2(i,jm2)*ds2(
     &i,jm1))

      X(0,-1)=(1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*dk2(im1,j)*ds2(
     &i,jm1)) + (1.0*gam12(i,jm1))/(dg1(i,jm1)*dk2(i,j)*ds2(i,jm1))))/
     &dh1(i,j) + ((1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)) + (gam2
     &2(i,jm1)*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i,jm1)))/ds2(i,j
     &m1) + gam22(i,j)/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + ((0.
     &25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(i,jm1))/ds2(i,jm1))/
     &dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,jm1) - (0.25*sig21(i,j
     &m1))/ds1(i,jm1))/dh2(i,j) + sig22(i,jm1)/(dh2(i,j)*ds2(i,jm1)) +
     & (1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*ds1(im1,jm1)) + (1.0*g
     &am21(i,jm1))/(dg2(i,jm1)*ds1(i,jm1))))/(dh2(i,j)*dk1(i,jm1))

      X(0,0)=- (1.0*(((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(i,jm1)) +
     & (1.0*gam12(im1,j))/(dg1(im1,j)*ds2(i,j)))/dk2(im1,j) + (1.0*((1
     &.0*gam12(i,jm1))/(dg1(i,jm1)*ds2(i,jm1)) + (1.0*gam12(i,j))/(dg1
     &(i,j)*ds2(i,j))))/dk2(i,j)))/dh1(i,j) - (1.0*(((1.0*gam21(im1,jm
     &1))/(dg2(im1,jm1)*ds1(im1,j)) + (1.0*gam21(i,jm1))/(dg2(i,jm1)*d
     &s1(i,j)))/dk1(i,jm1) + (1.0*((1.0*gam21(im1,j))/(dg2(im1,j)*ds1(
     &im1,j)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j))))/dk1(i,j)))/dh2(
     &i,j) - ((1.0*sig11(im1,j))/ds1(im1,j) + (1.0*sig11(i,j))/ds1(i,j
     &))/dh1(i,j) - ((1.0*sig22(i,jm1))/ds2(i,jm1) + (1.0*sig22(i,j))/
     &ds2(i,j))/dh2(i,j) - ((1.0*((1.0*gam22(i,jm1))/(dh2(i,jm1)*ds2(i
     &,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,
     &j)))/ds2(i,jm1) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(
     &dh2(i,j)*ds2(i,j)))/dh2(i,j) - (((1.0*gam11(ip1,j))/(dh1(ip1,j)*
     &ds1(i,j)) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,j
     &))/ds1(i,j) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(
     &i,j)*ds1(im1,j)))/dh1(i,j) - ((0.25*sig12(im1,jm1))/ds2(i,jm1) -
     & (0.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,jm1))/ds2(i,jm1) +
     & (0.25*sig12(i,j))/ds2(i,j))/dh1(i,j) - ((0.25*sig21(im1,jm1))/d
     &s1(im1,j) - (0.25*sig21(im1,j))/ds1(im1,j) - (0.25*sig21(i,jm1))
     &/ds1(i,j) + (0.25*sig21(i,j))/ds1(i,j))/dh2(i,j)

      X(0,+1)=(1.0*((1.0*gam12(im1,j))/(dg1(im1,j)*dk2(im1,j)*ds2(i,j)
     &) + (1.0*gam12(i,j))/(dg1(i,j)*dk2(i,j)*ds2(i,j))))/dh1(i,j) - (
     &(0.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,j))/ds2(i,j))/dh1(i
     &,j) + ((0.25*sig21(im1,j))/ds1(im1,jp1) - (0.25*sig21(i,j))/ds1(
     &i,jp1))/dh2(i,j) + (((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) - (e4c
     &oefint*gam22(i,jp1)*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/dh1(i,j
     &p1))/ds2(i,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))
     &/dh2(i,j) + (1.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*((1.0*ga
     &m21(im1,j))/(dg2(im1,j)*ds1(im1,jp1)) + (1.0*gam21(i,j))/(dg2(i,
     &j)*ds1(i,jp1))))/(dh2(i,j)*dk1(i,j))

      X(+1,-1)=- (0.25*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1)) - (0.25*s
     &ig21(i,jm1))/(dh2(i,j)*ds1(i,jm1)) - (1.0*gam12(i,jm1))/(dg1(i,j
     &m1)*dh1(i,j)*dk2(i,j)*ds2(ip1,jm1)) - (1.0*gam21(i,jm1))/(dg2(i,
     &jm1)*dh2(i,j)*dk1(i,jm1)*ds1(i,jm1))

      X(+1,0)=(1.0*((1.0*gam21(i,jm1))/(dg2(i,jm1)*dk1(i,jm1)*ds1(i,j)
     &) + (1.0*gam21(i,j))/(dg2(i,j)*dk1(i,j)*ds1(i,j))))/dh2(i,j) + (
     &(0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j))/ds2(ip1,j))
     &/dh1(i,j) - ((0.25*sig21(i,jm1))/ds1(i,j) - (0.25*sig21(i,j))/ds
     &1(i,j))/dh2(i,j) + (((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0
     &*gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/dh1(ip1,j))/ds1(i
     &,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) 
     &+ (1.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*((1.0*gam12(i,jm1)
     &)/(dg1(i,jm1)*ds2(ip1,jm1)) + (1.0*gam12(i,j))/(dg1(i,j)*ds2(ip1
     &,j))))/(dh1(i,j)*dk2(i,j))

      X(+1,+1)=(0.25*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (0.25*sig21(i
     &,j))/(dh2(i,j)*ds1(i,jp1)) - (1.0*gam12(i,j))/(dg1(i,j)*dh1(i,j)
     &*dk2(i,j)*ds2(ip1,j)) - (1.0*gam21(i,j))/(dg2(i,j)*dh2(i,j)*dk1(
     &i,j)*ds1(i,jp1)) + (1.0*e4coefint*gam22(i,jp1))/(dh1(i,jp1)*dh2(
     &i,j)*ds1(i,jp1)*ds2(i,j))

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY Nx-1,Ny-1
           elseif(  (
     &         (  (j.eq. ns_ibmin-1) 
     &            .and. (    (BoundaryFibShellin(1,2) .eq. 1)
     &                   .or.(BoundaryFibShellin(1,2) .eq. 2)
     &                   ) 
     &         ) 
     &        .and.
     &        (   (i.eq. nq_ibmin-1) 
     &            .and. (    (BoundaryFibShellin(2,2) .eq. 1)
     &                   .or.(BoundaryFibShellin(2,2) .eq. 2)
     &                   ) 
     &         ) 
     &       ) ) then

      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,-1)=(0.25*sig12(im1,jm1))/(dh1(i,j)*ds2(im1,jm1)) + (0.25*s
     &ig21(im1,jm1))/(dh2(i,j)*ds1(im1,jm1)) - (1.0*gam12(im1,jm1))/(d
     &g1(im1,jm1)*dh1(i,j)*dk2(im1,j)*ds2(im1,jm1)) - (1.0*gam21(im1,j
     &m1))/(dg2(im1,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(im1,jm1))

      X(-1,0)=(1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*dk1(i,jm1)*ds1(
     &im1,j)) + (1.0*gam21(im1,j))/(dg2(im1,j)*dk1(i,j)*ds1(im1,j))))/
     &dh2(i,j) + ((1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)) + (gam1
     &1(im1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(im1,j)))/ds1(im1
     &,j) + gam11(i,j)/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) - ((0.
     &25*sig12(im1,jm1))/ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(im1,j)
     &)/dh1(i,j) + ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.25*sig21(im1
     &,j))/ds1(im1,j))/dh2(i,j) + sig11(im1,j)/(dh1(i,j)*ds1(im1,j)) +
     & (1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(im1,jm1)) + (1.0*g
     &am12(im1,j))/(dg1(im1,j)*ds2(im1,j))))/(dh1(i,j)*dk2(im1,j))

      X(-1,+1)=(e4coefint*gam22(i,jp1))/(dh1(i,jp1)*dh2(i,j)*ds1(im1,j
     &p1)*ds2(i,j)) - (0.25*sig21(im1,j))/(dh2(i,j)*ds1(im1,jp1)) - (1
     &.0*gam12(im1,j))/(dg1(im1,j)*dh1(i,j)*dk2(im1,j)*ds2(im1,j)) - (
     &1.0*gam21(im1,j))/(dg2(im1,j)*dh2(i,j)*dk1(i,j)*ds1(im1,jp1)) - 
     &(0.25*sig12(im1,j))/(dh1(i,j)*ds2(im1,j))

      X(0,-2)=-(1.0*gam22(i,jm1))/(dh2(i,jm1)*dh2(i,j)*ds2(i,jm2)*ds2(
     &i,jm1))

      X(0,-1)=(1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*dk2(im1,j)*ds2(
     &i,jm1)) + (1.0*gam12(i,jm1))/(dg1(i,jm1)*dk2(i,j)*ds2(i,jm1))))/
     &dh1(i,j) + ((1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)) + (gam2
     &2(i,jm1)*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i,jm1)))/ds2(i,j
     &m1) + gam22(i,j)/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + ((0.
     &25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(i,jm1))/ds2(i,jm1))/
     &dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,jm1) - (0.25*sig21(i,j
     &m1))/ds1(i,jm1))/dh2(i,j) + sig22(i,jm1)/(dh2(i,j)*ds2(i,jm1)) +
     & (1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*ds1(im1,jm1)) + (1.0*g
     &am21(i,jm1))/(dg2(i,jm1)*ds1(i,jm1))))/(dh2(i,j)*dk1(i,jm1))

      X(0,0)=- (1.0*(((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(i,jm1)) +
     & (1.0*gam12(im1,j))/(dg1(im1,j)*ds2(i,j)))/dk2(im1,j) + (1.0*((1
     &.0*gam12(i,jm1))/(dg1(i,jm1)*ds2(i,jm1)) + (1.0*gam12(i,j))/(dg1
     &(i,j)*ds2(i,j))))/dk2(i,j)))/dh1(i,j) - (1.0*(((1.0*gam21(im1,jm
     &1))/(dg2(im1,jm1)*ds1(im1,j)) + (1.0*gam21(i,jm1))/(dg2(i,jm1)*d
     &s1(i,j)))/dk1(i,jm1) + (1.0*((1.0*gam21(im1,j))/(dg2(im1,j)*ds1(
     &im1,j)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j))))/dk1(i,j)))/dh2(
     &i,j) - ((1.0*sig11(im1,j))/ds1(im1,j) + (1.0*sig11(i,j))/ds1(i,j
     &))/dh1(i,j) - ((1.0*sig22(i,jm1))/ds2(i,jm1) + (1.0*sig22(i,j))/
     &ds2(i,j))/dh2(i,j) - ((1.0*((1.0*gam11(im1,j))/(dh1(im1,j)*ds1(i
     &m1,j)) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,
     &j)))/ds1(im1,j) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(
     &dh1(i,j)*ds1(i,j)))/dh1(i,j) - ((1.0*((1.0*gam22(i,jm1))/(dh2(i,
     &jm1)*ds2(i,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j
     &)))/dh2(i,j)))/ds2(i,jm1) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds
     &2(i,j)))/(dh2(i,j)*ds2(i,j)))/dh2(i,j) - ((0.25*sig12(im1,jm1))/
     &ds2(i,jm1) - (0.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,jm1))/
     &ds2(i,jm1) + (0.25*sig12(i,j))/ds2(i,j))/dh1(i,j) - ((0.25*sig21
     &(im1,jm1))/ds1(im1,j) - (0.25*sig21(im1,j))/ds1(im1,j) - (0.25*s
     &ig21(i,jm1))/ds1(i,j) + (0.25*sig21(i,j))/ds1(i,j))/dh2(i,j)

      X(0,+1)=(1.0*((1.0*gam12(im1,j))/(dg1(im1,j)*dk2(im1,j)*ds2(i,j)
     &) + (1.0*gam12(i,j))/(dg1(i,j)*dk2(i,j)*ds2(i,j))))/dh1(i,j) - (
     &(0.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,j))/ds2(i,j))/dh1(i
     &,j) + ((0.25*sig21(im1,j))/ds1(im1,jp1) - (0.25*sig21(i,j))/ds1(
     &i,jp1))/dh2(i,j) + (((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) - (e4c
     &oefint*gam22(i,jp1)*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/dh1(i,j
     &p1))/ds2(i,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))
     &/dh2(i,j) + (1.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*((1.0*ga
     &m21(im1,j))/(dg2(im1,j)*ds1(im1,jp1)) + (1.0*gam21(i,j))/(dg2(i,
     &j)*ds1(i,jp1))))/(dh2(i,j)*dk1(i,j))

      X(+1,-1)=(e4coefint*gam11(ip1,j))/(dh1(i,j)*dh2(ip1,j)*ds1(i,j)*
     &ds2(ip1,jm1)) - (0.25*sig21(i,jm1))/(dh2(i,j)*ds1(i,jm1)) - (1.0
     &*gam12(i,jm1))/(dg1(i,jm1)*dh1(i,j)*dk2(i,j)*ds2(ip1,jm1)) - (1.
     &0*gam21(i,jm1))/(dg2(i,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(i,jm1)) - (0
     &.25*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1))

      X(+1,0)=(1.0*((1.0*gam21(i,jm1))/(dg2(i,jm1)*dk1(i,jm1)*ds1(i,j)
     &) + (1.0*gam21(i,j))/(dg2(i,j)*dk1(i,j)*ds1(i,j))))/dh2(i,j) + (
     &(0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j))/ds2(ip1,j))
     &/dh1(i,j) - ((0.25*sig21(i,jm1))/ds1(i,j) - (0.25*sig21(i,j))/ds
     &1(i,j))/dh2(i,j) + (((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) - (e4c
     &oefint*gam11(ip1,j)*(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/dh2(ip1
     &,j))/ds1(i,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))
     &/dh1(i,j) + (1.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*((1.0*ga
     &m12(i,jm1))/(dg1(i,jm1)*ds2(ip1,jm1)) + (1.0*gam12(i,j))/(dg1(i,
     &j)*ds2(ip1,j))))/(dh1(i,j)*dk2(i,j))

      X(+1,+1)=(0.25*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (0.25*sig21(i
     &,j))/(dh2(i,j)*ds1(i,jp1)) - (1.0*gam12(i,j))/(dg1(i,j)*dh1(i,j)
     &*dk2(i,j)*ds2(ip1,j)) - (1.0*gam21(i,j))/(dg2(i,j)*dh2(i,j)*dk1(
     &i,j)*ds1(i,jp1)) + (1.0*e4coefint*gam11(ip1,j))/(dh1(i,j)*dh2(ip
     &1,j)*ds1(i,j)*ds2(ip1,j)) + (1.0*e4coefint*gam22(i,jp1))/(dh1(i,
     &jp1)*dh2(i,j)*ds1(i,jp1)*ds2(i,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY 1,2  
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .eq. 2)) then

!= ===============  1,2 VER2
!CCCCCCCC+X-xCCCCCCCC

      X(0,-1)=((gam22(i,j)/(dh2(i,j)*ds2(i,jm1)) - (1.0*e4coefint*gam2
     &2(i,jm1)*((1.5*(2.0*ds1(i,jm1) - 3.0*ds1(ip1,jm1) + ds1(ip2,jm1)
     &))/(1.5*ds1(i,jm1) - 0.5*ds1(ip1,jm1)) - 2.0))/(1.5*ds1(i,jm1) -
     & 0.5*ds1(ip1,jm1))**2)/ds2(i,jm1) + gam22(i,j)/(dh2(i,j)*ds2(i,j
     &m1)*ds2(i,j)))/dh2(i,j) - (0.5*sig12(i,jm1))/(dh1(i,j)*ds2(i,jm1
     &)) + sigB(1,jm1,2)/(dh2(i,j)*ds2(i,jm1)) + (0.75*sigB(1,jm1,3))/
     &(dh2(i,jm1)*(1.5*ds1(i,jm1) - 0.5*ds1(ip1,jm1))) - (1.0*(e4coefi
     &nt - 2.0))/(dh1(i,j)**2*dh2(i,j)*ds2(i,jm1)) + (2.0*e6e*gam21(i,
     &j))/(dh1(i,j)*dh2(i,j)*ds2(i,jm1))

      X(0,0)=(2.0*((0.25*sig12(i,jm1))/ds2(i,jm1) - (0.25*sig12(i,j))/
     &ds2(i,j)))/dh1(i,j) - ((1.0*sigB(1,jm1,2))/ds2(i,jm1) + (1.0*sig
     &B(1,j,2))/ds2(i,j))/dh2(i,j) + ((0.75*sigB(1,jm1,3))/(1.5*ds1(i,
     &j) - 0.5*ds1(ip1,j)) - (0.75*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds
     &1(ip1,j)))/dh2(i,jm1) - (((1.0*gam22(i,jp1))/(dh2(i,jp1)*ds2(i,j
     &)) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j))/ds2(
     &i,j) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds
     &2(i,jm1)))/dh2(i,j) - (2.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.
     &0*(e4coefint - 2.0)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh1(i,j)**
     &2*dh2(i,j)) - (1.0*gam11(i,j))/(dh1(i,j)**2*dh1(ip1,j)*ds1(i,j))
     & + gam11(ip1,j)/(dh1(i,j)*dh1(ip1,j)**2*ds1(i,j)) - (2.0*e6e*gam
     &21(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh1(i,j)*dh2(i,j))

      X(0,+1)=(((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*gam22(i,jp
     &1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1))/ds2(i,j) + (1.0*
     &gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + (0.5*sig1
     &2(i,j))/(dh1(i,j)*ds2(i,j)) + (1.0*sigB(1,j,2))/(dh2(i,j)*ds2(i,
     &j)) - (0.75*sigB(1,j,3))/(dh2(i,jm1)*(1.5*ds1(i,jp1) - 0.5*ds1(i
     &p1,jp1))) - (1.0*(e4coefint - 2.0))/(dh1(i,j)**2*dh2(i,j)*ds2(i,
     &j)) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*dh2(i,j)*ds2(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(+1,-1)=gam21(ip1,j)/(dh1(i,j)*dh1(ip1,j)*dh2(ip1,j)*ds2(ip1,jm
     &1)) - (1.0*sigB(1,jm1,3))/(dh2(i,jm1)*(1.5*ds1(i,jm1) - 0.5*ds1(
     &ip1,jm1))) - (2.0*gam11(i,j))/(dh1(i,j)**2*dh2(ip1,j)*ds2(ip1,jm
     &1)) - (2.0*e6e*gam21(i,j))/(dh1(i,j)*dh2(ip1,j)*ds2(ip1,jm1)) - 
     &(0.5*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1)) + (1.0*e4coefint*gam2
     &2(i,jm1)*((2.0*(2.0*ds1(i,jm1) - 3.0*ds1(ip1,jm1) + ds1(ip2,jm1)
     &))/(1.5*ds1(i,jm1) - 0.5*ds1(ip1,jm1)) - 5.0))/(dh2(i,j)*ds2(i,j
     &m1)*(1.5*ds1(i,jm1) - 0.5*ds1(ip1,jm1))**2)

      X(+1,0)=(2.0*((0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j
     &))/ds2(ip1,j)))/dh1(i,j) - ((1.0*sigB(1,jm1,3))/(1.5*ds1(i,j) - 
     &0.5*ds1(ip1,j)) - (1.0*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,
     &j)))/dh2(i,jm1) + (2.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*ga
     &m11(i,j)*((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + (2.0*(1.0
     &/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/dh2(ip1,j)))/dh1(i,j)**2 - (gam
     &11(ip1,j)*((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + 1.0/(dh1
     &(ip2,j)*ds1(ip1,j))))/(dh1(i,j)*dh1(ip1,j)) - (gam21(ip1,j)*(1.0
     &/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/(dh1(i,j)*dh1(ip1,j)*dh2(ip1,j)
     &) + (2.0*e6e*gam21(i,j)*(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/(dh
     &1(i,j)*dh2(ip1,j))

      X(+1,+1)=(0.5*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (1.0*sigB(1,j,
     &3))/(dh2(i,jm1)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) - (2.0*gam1
     &1(i,j))/(dh1(i,j)**2*dh2(ip1,j)*ds2(ip1,j)) - (2.0*e6e*gam21(i,j
     &))/(dh1(i,j)*dh2(ip1,j)*ds2(ip1,j)) + (1.0*gam21(ip1,j))/(dh1(i,
     &j)*dh1(ip1,j)*dh2(ip1,j)*ds2(ip1,j))

      X(+2,-1)=(0.25*sigB(1,jm1,3))/(dh2(i,jm1)*(1.5*ds1(i,jm1) - 0.5*
     &ds1(ip1,jm1))) - (1.0*gam21(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh2(ip2
     &,j)*ds2(ip2,jm1)) - (1.0*e4coefint*gam22(i,jm1)*((0.5*(2.0*ds1(i
     &,jm1) - 3.0*ds1(ip1,jm1) + ds1(ip2,jm1)))/(1.5*ds1(i,jm1) - 0.5*
     &ds1(ip1,jm1)) - 4.0))/(dh2(i,j)*ds2(i,jm1)*(1.5*ds1(i,jm1) - 0.5
     &*ds1(ip1,jm1))**2)

      X(+2,0)=((0.25*sigB(1,jm1,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - 
     &(0.25*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)))/dh2(i,jm1) -
     & (1.0*gam11(i,j))/(dh1(i,j)**2*dh1(ip1,j)*ds1(ip1,j)) + (gam11(i
     &p1,j)*((1.0*(1.0/ds1(ip1,j) + 1.0/ds1(ip2,j)))/dh1(ip2,j) + 1.0/
     &(dh1(ip1,j)*ds1(ip1,j))))/(dh1(i,j)*dh1(ip1,j)) + (1.0*gam21(ip1
     &,j)*(1.0/ds2(ip2,jm1) + 1.0/ds2(ip2,j)))/(dh1(i,j)*dh1(ip1,j)*dh
     &2(ip2,j))

      X(+2,+1)=- (0.25*sigB(1,j,3))/(dh2(i,jm1)*(1.5*ds1(i,jp1) - 0.5*
     &ds1(ip1,jp1))) - (1.0*gam21(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh2(ip2
     &,j)*ds2(ip2,j))

      X(+3,-1)=-(1.0*e4coefint*gam22(i,jm1))/(dh2(i,j)*ds2(i,jm1)*(1.5
     &*ds1(i,jm1) - 0.5*ds1(ip1,jm1))**2)

      X(+3,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh1(ip2,j)*ds1(
     &ip2,j))


      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY 1,N-1  
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .eq. ns_ibmin-1)) then

!= ===============  1,N-1 VER2
!CCCCCCCC+X-xCCCCCCCC

      X(0,-2)=-(1.0*gam22(i,jm1))/(dh2(i,jm1)*dh2(i,j)*ds2(i,jm2)*ds2(
     &i,jm1))

      X(0,-1)=(1.0*(((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)) + (gam22(i
     &,jm1)*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i,jm1))/ds2(i,jm1) 
     &+ (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) - (
     &0.5*sig12(i,jm1))/(dh1(i,j)*ds2(i,jm1)) + sigB(1,jm1,2)/(dh2(i,j
     &)*ds2(i,jm1)) + (0.75*sigB(1,jm1,3))/(dh2(i,jm1)*(1.5*ds1(i,jm1)
     & - 0.5*ds1(ip1,jm1))) - (1.0*(e4coefint - 2.0))/(dh1(i,j)**2*dh2
     &(i,j)*ds2(i,jm1)) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*dh2(i,j)*ds2(
     &i,jm1))

      X(0,0)=(2.0*((0.25*sig12(i,jm1))/ds2(i,jm1) - (0.25*sig12(i,j))/
     &ds2(i,j)))/dh1(i,j) - ((1.0*sigB(1,jm1,2))/ds2(i,jm1) + (1.0*sig
     &B(1,j,2))/ds2(i,j))/dh2(i,j) - (1.0*(((1.0*gam22(i,jm1))/(dh2(i,
     &jm1)*ds2(i,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j
     &)))/dh2(i,j))/ds2(i,jm1) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0
     &/ds2(i,j)))/(dh2(i,j)*ds2(i,j))))/dh2(i,j) + ((0.75*sigB(1,jm1,3
     &))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - (0.75*sigB(1,j,3))/(1.5*ds1
     &(i,j) - 0.5*ds1(ip1,j)))/dh2(i,jm1) - (2.0*sig11(i,j))/(dh1(i,j)
     &*ds1(i,j)) + (1.0*(e4coefint - 2.0)*(1.0/ds2(i,jm1) + 1.0/ds2(i,
     &j)))/(dh1(i,j)**2*dh2(i,j)) - (1.0*gam11(i,j))/(dh1(i,j)**2*dh1(
     &ip1,j)*ds1(i,j)) + gam11(ip1,j)/(dh1(i,j)*dh1(ip1,j)**2*ds1(i,j)
     &) - (2.0*e6e*gam21(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh1(i,
     &j)*dh2(i,j))

      X(0,+1)=(1.0*((1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) - (1.0*
     &e4coefint*gam22(i,jp1)*((1.5*(2.0*ds1(i,jp1) - 3.0*ds1(ip1,jp1) 
     &+ ds1(ip2,jp1)))/(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1)) - 2.0))/(1.
     &5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))**2))/ds2(i,j) + (1.0*gam22(i,j)
     &)/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) + (0.5*sig12(i,j))/(
     &dh1(i,j)*ds2(i,j)) + (1.0*sigB(1,j,2))/(dh2(i,j)*ds2(i,j)) - (0.
     &75*sigB(1,j,3))/(dh2(i,jm1)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1)))
     & - (1.0*(e4coefint - 2.0))/(dh1(i,j)**2*dh2(i,j)*ds2(i,j)) + (2.
     &0*e6e*gam21(i,j))/(dh1(i,j)*dh2(i,j)*ds2(i,j))

      X(+1,-1)=gam21(ip1,j)/(dh1(i,j)*dh1(ip1,j)*dh2(ip1,j)*ds2(ip1,jm
     &1)) - (1.0*sigB(1,jm1,3))/(dh2(i,jm1)*(1.5*ds1(i,jm1) - 0.5*ds1(
     &ip1,jm1))) - (2.0*gam11(i,j))/(dh1(i,j)**2*dh2(ip1,j)*ds2(ip1,jm
     &1)) - (2.0*e6e*gam21(i,j))/(dh1(i,j)*dh2(ip1,j)*ds2(ip1,jm1)) - 
     &(0.5*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1))

      X(+1,0)=(2.0*((0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j
     &))/ds2(ip1,j)))/dh1(i,j) - ((1.0*sigB(1,jm1,3))/(1.5*ds1(i,j) - 
     &0.5*ds1(ip1,j)) - (1.0*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,
     &j)))/dh2(i,jm1) + (2.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*ga
     &m11(i,j)*((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + (2.0*(1.0
     &/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/dh2(ip1,j)))/dh1(i,j)**2 - (gam
     &11(ip1,j)*((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + 1.0/(dh1
     &(ip2,j)*ds1(ip1,j))))/(dh1(i,j)*dh1(ip1,j)) - (gam21(ip1,j)*(1.0
     &/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/(dh1(i,j)*dh1(ip1,j)*dh2(ip1,j)
     &) + (2.0*e6e*gam21(i,j)*(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/(dh
     &1(i,j)*dh2(ip1,j))

      X(+1,+1)=(0.5*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (1.0*sigB(1,j,
     &3))/(dh2(i,jm1)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) - (2.0*gam1
     &1(i,j))/(dh1(i,j)**2*dh2(ip1,j)*ds2(ip1,j)) - (2.0*e6e*gam21(i,j
     &))/(dh1(i,j)*dh2(ip1,j)*ds2(ip1,j)) + (1.0*gam21(ip1,j))/(dh1(i,
     &j)*dh1(ip1,j)*dh2(ip1,j)*ds2(ip1,j)) + (1.0*e4coefint*gam22(i,jp
     &1)*((2.0*(2.0*ds1(i,jp1) - 3.0*ds1(ip1,jp1) + ds1(ip2,jp1)))/(1.
     &5*ds1(i,jp1) - 0.5*ds1(ip1,jp1)) - 5.0))/(dh2(i,j)*ds2(i,j)*(1.5
     &*ds1(i,jp1) - 0.5*ds1(ip1,jp1))**2)

      X(+2,-1)=(0.25*sigB(1,jm1,3))/(dh2(i,jm1)*(1.5*ds1(i,jm1) - 0.5*
     &ds1(ip1,jm1))) - (1.0*gam21(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh2(ip2
     &,j)*ds2(ip2,jm1))

      X(+2,0)=((0.25*sigB(1,jm1,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - 
     &(0.25*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)))/dh2(i,jm1) -
     & (1.0*gam11(i,j))/(dh1(i,j)**2*dh1(ip1,j)*ds1(ip1,j)) + (gam11(i
     &p1,j)*((1.0*(1.0/ds1(ip1,j) + 1.0/ds1(ip2,j)))/dh1(ip2,j) + 1.0/
     &(dh1(ip1,j)*ds1(ip1,j))))/(dh1(i,j)*dh1(ip1,j)) + (1.0*gam21(ip1
     &,j)*(1.0/ds2(ip2,jm1) + 1.0/ds2(ip2,j)))/(dh1(i,j)*dh1(ip1,j)*dh
     &2(ip2,j))

      X(+2,+1)=- (0.25*sigB(1,j,3))/(dh2(i,jm1)*(1.5*ds1(i,jp1) - 0.5*
     &ds1(ip1,jp1))) - (1.0*gam21(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh2(ip2
     &,j)*ds2(ip2,j)) - (1.0*e4coefint*gam22(i,jp1)*((0.5*(2.0*ds1(i,j
     &p1) - 3.0*ds1(ip1,jp1) + ds1(ip2,jp1)))/(1.5*ds1(i,jp1) - 0.5*ds
     &1(ip1,jp1)) - 4.0))/(dh2(i,j)*ds2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1
     &(ip1,jp1))**2)

      X(+3,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh1(ip2,j)*ds1(
     &ip2,j))

      X(+3,+1)=-(1.0*e4coefint*gam22(i,jp1))/(dh2(i,j)*ds2(i,j)*(1.5*d
     &s1(i,jp1) - 0.5*ds1(ip1,jp1))**2)



      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0





!! ================== BNDY 2,1  
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .eq. 2)) then

!= ===============  2,1 VER2
!CCCCCCCC+X-xCCCCCCCC

      X(-1,0)=(gam11(i,j)/(dh1(i,j)*ds1(im1,j)) - (1.0*e4coefint*gam11
     &(im1,j)*((1.5*(2.0*ds2(im1,j) - 3.0*ds2(im1,jp1) + ds2(im1,jp2))
     &)/(1.5*ds2(im1,j) - 0.5*ds2(im1,jp1)) - 2.0))/(1.5*ds2(im1,j) - 
     &0.5*ds2(im1,jp1))**2)/(dh1(i,j)*ds1(im1,j)) + sigB(3,im1,1)/(dh1
     &(i,j)*ds1(im1,j)) - (0.5*sig12(im1,j))/(dh2(i,j)*ds1(im1,j)) + (
     &0.75*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.5*ds2(im1,jp
     &1))) - (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(i,j)*dh2(i,j)**2*
     &ds1(im1,j)) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*dh2(i,j)**2*ds1(im1
     &,j))

      X(-1,+1)=(2.0*gam21(i,jp1))/(dh1(i,jp1)*dh2(i,j)*dh2(i,jp1)*ds1(
     &im1,jp1)) - (1.0*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.
     &5*ds2(im1,jp1))) - (2.0*gam22(i,j))/(dh1(i,jp1)*dh2(i,j)**2*ds1(
     &im1,jp1)) - (2.0*e6e*gam21(i,j))/(dh1(i,jp1)*dh2(i,j)**2*ds1(im1
     &,jp1)) - (0.5*sig12(im1,j))/(dh2(i,j)*ds1(im1,jp1)) + (1.0*e4coe
     &fint*gam11(im1,j)*((2.0*(2.0*ds2(im1,j) - 3.0*ds2(im1,jp1) + ds2
     &(im1,jp2)))/(1.5*ds2(im1,j) - 0.5*ds2(im1,jp1)) - 5.0))/(dh1(i,j
     &)*ds1(im1,j)*(1.5*ds2(im1,j) - 0.5*ds2(im1,jp1))**2)

      X(-1,+2)=(0.25*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.5*
     &ds2(im1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*dh2(i,j)*dh2(i,j
     &p1)*ds1(im1,jp2)) - (1.0*e4coefint*gam11(im1,j)*((0.5*(2.0*ds2(i
     &m1,j) - 3.0*ds2(im1,jp1) + ds2(im1,jp2)))/(1.5*ds2(im1,j) - 0.5*
     &ds2(im1,jp1)) - 4.0))/(dh1(i,j)*ds1(im1,j)*(1.5*ds2(im1,j) - 0.5
     &*ds2(im1,jp1))**2)

      X(-1,+3)=-(1.0*e4coefint*gam11(im1,j))/(dh1(i,j)*ds1(im1,j)*(1.5
     &*ds2(im1,j) - 0.5*ds2(im1,jp1))**2)

      X(0,0)=(2.0*((0.25*sig12(im1,j))/ds1(im1,j) - (0.25*sig12(i,j))/
     &ds1(i,j)))/dh2(i,j) - ((1.0*sigB(3,im1,1))/ds1(im1,j) + (1.0*sig
     &B(3,i,1))/ds1(i,j))/dh1(i,j) - ((1.0*(gam11(ip1,j)/(dh1(ip1,j)*d
     &s1(i,j)) - (1.0*e4coefint*gam11(i,j)*((1.5*(2.0*ds2(i,j) - 3.0*d
     &s2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 2.0))
     &/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**2))/ds1(i,j) + (gam11(i,j)*(1.
     &0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds1(im1,j)))/dh1(i,j) + 
     &((0.75*sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - (0.75*si
     &gB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) - (2.0*si
     &g22(i,j))/(dh2(i,j)*ds2(i,j)) - (1.0*gam22(i,j))/(dh2(i,j)**2*dh
     &2(i,jp1)*ds2(i,j)) + gam22(i,jp1)/(dh2(i,j)*dh2(i,jp1)**2*ds2(i,
     &j)) + (1.0*gam22(i,j)*(e4coefint - 2.0)*(1.0/ds1(im1,j) + 1.0/ds
     &1(i,j)))/(dh1(i,j)*dh2(i,j)**2) - (2.0*e6e*gam21(i,j)*(1.0/ds1(i
     &m1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*dh2(i,j)**2)

      X(0,+1)=(2.0*((0.25*sig12(im1,j))/ds1(im1,jp1) - (0.25*sig12(i,j
     &))/ds1(i,jp1)))/dh2(i,j) - ((1.0*sigB(3,im1,3))/(1.5*ds2(i,j) - 
     &0.5*ds2(i,jp1)) - (1.0*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1)))/dh1(im1,j) + (2.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*ga
     &m22(i,j)*((2.0*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/dh1(i,jp1) +
     & (1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1)))/dh2(i,j)**2 - (gam
     &22(i,jp1)*((1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1) + 1.0/(dh2
     &(i,jp2)*ds2(i,jp1))))/(dh2(i,j)*dh2(i,jp1)) - (2.0*gam21(i,jp1)*
     &(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,jp1)*dh2(i,j)*dh2(i,
     &jp1)) + (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))
     &/(dh1(i,jp1)*dh2(i,j)**2) - (1.0*e4coefint*gam11(i,j)*((2.0*(2.0
     &*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds
     &2(i,jp1)) - 5.0))/(dh1(i,j)*ds1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,j
     &p1))**2)

      X(0,+2)=((0.25*sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 
     &(0.25*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) -
     & (1.0*gam22(i,j))/(dh2(i,j)**2*dh2(i,jp1)*ds2(i,jp1)) + (gam22(i
     &,jp1)*((1.0*(1.0/ds2(i,jp1) + 1.0/ds2(i,jp2)))/dh2(i,jp2) + 1.0/
     &(dh2(i,jp1)*ds2(i,jp1))))/(dh2(i,j)*dh2(i,jp1)) + (2.0*gam21(i,j
     &p1)*(1.0/ds1(im1,jp2) + 1.0/ds1(i,jp2)))/(dh1(i,jp2)*dh2(i,j)*dh
     &2(i,jp1)) + (1.0*e4coefint*gam11(i,j)*((0.5*(2.0*ds2(i,j) - 3.0*
     &ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 4.0)
     &)/(dh1(i,j)*ds1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**2)

      X(0,+3)=(1.0*e4coefint*gam11(i,j))/(dh1(i,j)*ds1(i,j)*(1.5*ds2(i
     &,j) - 0.5*ds2(i,jp1))**2) - (1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,j
     &p1)*dh2(i,jp2)*ds2(i,jp2))

      X(+1,0)=((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)) + (1.0*
     &gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/(dh1(ip1,j)*ds1(i,
     &j)))/dh1(i,j) + (1.0*sigB(3,i,1))/(dh1(i,j)*ds1(i,j)) + (0.5*sig
     &12(i,j))/(dh2(i,j)*ds1(i,j)) - (0.75*sigB(3,i,3))/(dh1(im1,j)*(1
     &.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (1.0*gam22(i,j)*(e4coefint 
     &- 2.0))/(dh1(i,j)*dh2(i,j)**2*ds1(i,j)) + (2.0*e6e*gam21(i,j))/(
     &dh1(i,j)*dh2(i,j)**2*ds1(i,j))

      X(+1,+1)=(0.5*sig12(i,j))/(dh2(i,j)*ds1(i,jp1)) + (1.0*sigB(3,i,
     &3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (2.0*gam2
     &2(i,j))/(dh1(i,jp1)*dh2(i,j)**2*ds1(i,jp1)) - (2.0*e6e*gam21(i,j
     &))/(dh1(i,jp1)*dh2(i,j)**2*ds1(i,jp1)) + (2.0*gam21(i,jp1))/(dh1
     &(i,jp1)*dh2(i,j)*dh2(i,jp1)*ds1(i,jp1))

      X(+1,+2)=- (0.25*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*
     &ds2(ip1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*dh2(i,j)*dh2(i,j
     &p1)*ds1(i,jp2))

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))


      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY N-1,1  
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .eq. nq_ibmin-1)) then
!= ===============  N-1,1 VER2
!CCCCCCCC+X-xCCCCCCCC

      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,0)=(1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)) + 
     &(gam11(im1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/(dh1(im1,j)*ds1
     &(im1,j))))/dh1(i,j) + sigB(3,im1,1)/(dh1(i,j)*ds1(im1,j)) - (0.5
     &*sig12(im1,j))/(dh2(i,j)*ds1(im1,j)) + (0.75*sigB(3,im1,3))/(dh1
     &(im1,j)*(1.5*ds2(im1,j) - 0.5*ds2(im1,jp1))) - (1.0*gam22(i,j)*(
     &e4coefint - 2.0))/(dh1(i,j)*dh2(i,j)**2*ds1(im1,j)) + (2.0*e6e*g
     &am21(i,j))/(dh1(i,j)*dh2(i,j)**2*ds1(im1,j))

      X(-1,+1)=(2.0*gam21(i,jp1))/(dh1(i,jp1)*dh2(i,j)*dh2(i,jp1)*ds1(
     &im1,jp1)) - (1.0*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.
     &5*ds2(im1,jp1))) - (2.0*gam22(i,j))/(dh1(i,jp1)*dh2(i,j)**2*ds1(
     &im1,jp1)) - (2.0*e6e*gam21(i,j))/(dh1(i,jp1)*dh2(i,j)**2*ds1(im1
     &,jp1)) - (0.5*sig12(im1,j))/(dh2(i,j)*ds1(im1,jp1))

      X(-1,+2)=(0.25*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.5*
     &ds2(im1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*dh2(i,j)*dh2(i,j
     &p1)*ds1(im1,jp2))

      X(0,0)=(2.0*((0.25*sig12(im1,j))/ds1(im1,j) - (0.25*sig12(i,j))/
     &ds1(i,j)))/dh2(i,j) - ((1.0*sigB(3,im1,1))/ds1(im1,j) + (1.0*sig
     &B(3,i,1))/ds1(i,j))/dh1(i,j) + ((0.75*sigB(3,im1,3))/(1.5*ds2(i,
     &j) - 0.5*ds2(i,jp1)) - (0.75*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds
     &2(i,jp1)))/dh1(im1,j) - (1.0*(((1.0*gam11(im1,j))/(dh1(im1,j)*ds
     &1(im1,j)) - (1.0*e4coefint*gam11(i,j)*((1.5*(2.0*ds2(i,j) - 3.0*
     &ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 2.0)
     &)/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**2)/ds1(im1,j) + (1.0*gam11(i,
     &j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds1(i,j))))/dh1(i,
     &j) - (2.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) - (1.0*gam22(i,j))/(dh
     &2(i,j)**2*dh2(i,jp1)*ds2(i,j)) + gam22(i,jp1)/(dh2(i,j)*dh2(i,jp
     &1)**2*ds2(i,j)) + (1.0*gam22(i,j)*(e4coefint - 2.0)*(1.0/ds1(im1
     &,j) + 1.0/ds1(i,j)))/(dh1(i,j)*dh2(i,j)**2) - (2.0*e6e*gam21(i,j
     &)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*dh2(i,j)**2)

      X(0,+1)=(2.0*((0.25*sig12(im1,j))/ds1(im1,jp1) - (0.25*sig12(i,j
     &))/ds1(i,jp1)))/dh2(i,j) - ((1.0*sigB(3,im1,3))/(1.5*ds2(i,j) - 
     &0.5*ds2(i,jp1)) - (1.0*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1)))/dh1(im1,j) + (2.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*ga
     &m22(i,j)*((2.0*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/dh1(i,jp1) +
     & (1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1)))/dh2(i,j)**2 - (gam
     &22(i,jp1)*((1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1) + 1.0/(dh2
     &(i,jp2)*ds2(i,jp1))))/(dh2(i,j)*dh2(i,jp1)) - (2.0*gam21(i,jp1)*
     &(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,jp1)*dh2(i,j)*dh2(i,
     &jp1)) + (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))
     &/(dh1(i,jp1)*dh2(i,j)**2) - (1.0*e4coefint*gam11(i,j)*((2.0*(2.0
     &*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds
     &2(i,jp1)) - 5.0))/(dh1(i,j)*ds1(im1,j)*(1.5*ds2(i,j) - 0.5*ds2(i
     &,jp1))**2)

      X(0,+2)=((0.25*sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 
     &(0.25*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) -
     & (1.0*gam22(i,j))/(dh2(i,j)**2*dh2(i,jp1)*ds2(i,jp1)) + (gam22(i
     &,jp1)*((1.0*(1.0/ds2(i,jp1) + 1.0/ds2(i,jp2)))/dh2(i,jp2) + 1.0/
     &(dh2(i,jp1)*ds2(i,jp1))))/(dh2(i,j)*dh2(i,jp1)) + (2.0*gam21(i,j
     &p1)*(1.0/ds1(im1,jp2) + 1.0/ds1(i,jp2)))/(dh1(i,jp2)*dh2(i,j)*dh
     &2(i,jp1)) + (1.0*e4coefint*gam11(i,j)*((0.5*(2.0*ds2(i,j) - 3.0*
     &ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 4.0)
     &)/(dh1(i,j)*ds1(im1,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**2)

      X(0,+3)=(1.0*e4coefint*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*(1.5*ds2
     &(i,j) - 0.5*ds2(i,jp1))**2) - (1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i
     &,jp1)*dh2(i,jp2)*ds2(i,jp2))

      X(+1,0)=(1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) - (1.0*e4coef
     &int*gam11(ip1,j)*((1.5*(2.0*ds2(ip1,j) - 3.0*ds2(ip1,jp1) + ds2(
     &ip1,jp2)))/(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1)) - 2.0))/(1.5*ds2(
     &ip1,j) - 0.5*ds2(ip1,jp1))**2))/(dh1(i,j)*ds1(i,j)) + (1.0*sigB(
     &3,i,1))/(dh1(i,j)*ds1(i,j)) + (0.5*sig12(i,j))/(dh2(i,j)*ds1(i,j
     &)) - (0.75*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip
     &1,jp1))) - (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(i,j)*dh2(i,j)
     &**2*ds1(i,j)) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*dh2(i,j)**2*ds1(i
     &,j))

      X(+1,+1)=(0.5*sig12(i,j))/(dh2(i,j)*ds1(i,jp1)) + (1.0*sigB(3,i,
     &3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (2.0*gam2
     &2(i,j))/(dh1(i,jp1)*dh2(i,j)**2*ds1(i,jp1)) - (2.0*e6e*gam21(i,j
     &))/(dh1(i,jp1)*dh2(i,j)**2*ds1(i,jp1)) + (2.0*gam21(i,jp1))/(dh1
     &(i,jp1)*dh2(i,j)*dh2(i,jp1)*ds1(i,jp1)) + (1.0*e4coefint*gam11(i
     &p1,j)*((2.0*(2.0*ds2(ip1,j) - 3.0*ds2(ip1,jp1) + ds2(ip1,jp2)))/
     &(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1)) - 5.0))/(dh1(i,j)*ds1(i,j)*(
     &1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))**2)

      X(+1,+2)=- (0.25*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*
     &ds2(ip1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*dh2(i,j)*dh2(i,j
     &p1)*ds1(i,jp2)) - (1.0*e4coefint*gam11(ip1,j)*((0.5*(2.0*ds2(ip1
     &,j) - 3.0*ds2(ip1,jp1) + ds2(ip1,jp2)))/(1.5*ds2(ip1,j) - 0.5*ds
     &2(ip1,jp1)) - 4.0))/(dh1(i,j)*ds1(i,j)*(1.5*ds2(ip1,j) - 0.5*ds2
     &(ip1,jp1))**2)

      X(+1,+3)=-(1.0*e4coefint*gam11(ip1,j))/(dh1(i,j)*ds1(i,j)*(1.5*d
     &s2(ip1,j) - 0.5*ds2(ip1,jp1))**2)



      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY 1,j  2<=j<=Ny-1
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .gt. jst) .and. (j.lt. jen)) then
!= ===============  1,j VER2
!CCCCCCCC+X-xCCCCCCCC

      X(0,-2)=-(1.0*gam22(i,jm1))/(dh2(i,jm1)*dh2(i,j)*ds2(i,jm2)*ds2(
     &i,jm1))

      X(0,-1)=(1.0*(((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)) + (gam22(i
     &,jm1)*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i,jm1))/ds2(i,jm1) 
     &+ (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) - (
     &0.5*sig12(i,jm1))/(dh1(i,j)*ds2(i,jm1)) + sigB(1,jm1,2)/(dh2(i,j
     &)*ds2(i,jm1)) + (0.75*sigB(1,jm1,3))/(dh2(i,jm1)*(1.5*ds1(i,jm1)
     & - 0.5*ds1(ip1,jm1))) - (1.0*gam11(i,j)*(e4coefint - 2.0))/(dh1(
     &i,j)**2*dh2(i,j)*ds2(i,jm1)) + (2.0*e6e*gam21(i,j))/(dh1(i,j)**2
     &*dh2(i,j)*ds2(i,jm1))

      X(0,0)=(2.0*((0.25*sig12(i,jm1))/ds2(i,jm1) - (0.25*sig12(i,j))/
     &ds2(i,j)))/dh1(i,j) - ((1.0*sigB(1,jm1,2))/ds2(i,jm1) + (1.0*sig
     &B(1,j,2))/ds2(i,j))/dh2(i,j) - (1.0*(((1.0*gam22(i,jm1))/(dh2(i,
     &jm1)*ds2(i,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j
     &)))/dh2(i,j))/ds2(i,jm1) + (1.0*((1.0*gam22(i,jp1))/(dh2(i,jp1)*
     &ds2(i,j)) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j
     &)))/ds2(i,j)))/dh2(i,j) + ((0.75*sigB(1,jm1,3))/(1.5*ds1(i,j) - 
     &0.5*ds1(ip1,j)) - (0.75*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1
     &,j)))/dh2(i,jm1) - (2.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) - (1.0*g
     &am11(i,j))/(dh1(i,j)**2*dh1(ip1,j)*ds1(i,j)) + gam11(ip1,j)/(dh1
     &(i,j)*dh1(ip1,j)**2*ds1(i,j)) + (1.0*gam11(i,j)*(e4coefint - 2.0
     &)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh1(i,j)**2*dh2(i,j)) - (2.0
     &*e6e*gam21(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh1(i,j)**2*dh
     &2(i,j))

      X(0,+1)=(1.0*((1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*
     &gam22(i,jp1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1)))/ds2(i
     &,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j)
     & + (0.5*sig12(i,j))/(dh1(i,j)*ds2(i,j)) + (1.0*sigB(1,j,2))/(dh2
     &(i,j)*ds2(i,j)) - (0.75*sigB(1,j,3))/(dh2(i,jm1)*(1.5*ds1(i,jp1)
     & - 0.5*ds1(ip1,jp1))) - (1.0*gam11(i,j)*(e4coefint - 2.0))/(dh1(
     &i,j)**2*dh2(i,j)*ds2(i,j)) + (2.0*e6e*gam21(i,j))/(dh1(i,j)**2*d
     &h2(i,j)*ds2(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(+1,-1)=(2.0*gam21(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh2(ip1,j)*ds2(
     &ip1,jm1)) - (1.0*sigB(1,jm1,3))/(dh2(i,jm1)*(1.5*ds1(i,jm1) - 0.
     &5*ds1(ip1,jm1))) - (2.0*gam11(i,j))/(dh1(i,j)**2*dh2(ip1,j)*ds2(
     &ip1,jm1)) - (2.0*e6e*gam21(i,j))/(dh1(i,j)**2*dh2(ip1,j)*ds2(ip1
     &,jm1)) - (0.5*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1))

      X(+1,0)=(2.0*((0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j
     &))/ds2(ip1,j)))/dh1(i,j) - ((1.0*sigB(1,jm1,3))/(1.5*ds1(i,j) - 
     &0.5*ds1(ip1,j)) - (1.0*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,
     &j)))/dh2(i,jm1) + (2.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*ga
     &m11(i,j)*((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + (2.0*(1.0
     &/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/dh2(ip1,j)))/dh1(i,j)**2 - (gam
     &11(ip1,j)*((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + 1.0/(dh1
     &(ip2,j)*ds1(ip1,j))))/(dh1(i,j)*dh1(ip1,j)) - (2.0*gam21(ip1,j)*
     &(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/(dh1(i,j)*dh1(ip1,j)*dh2(ip
     &1,j)) + (2.0*e6e*gam21(i,j)*(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))
     &/(dh1(i,j)**2*dh2(ip1,j))

      X(+1,+1)=(0.5*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (1.0*sigB(1,j,
     &3))/(dh2(i,jm1)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) - (2.0*gam1
     &1(i,j))/(dh1(i,j)**2*dh2(ip1,j)*ds2(ip1,j)) - (2.0*e6e*gam21(i,j
     &))/(dh1(i,j)**2*dh2(ip1,j)*ds2(ip1,j)) + (2.0*gam21(ip1,j))/(dh1
     &(i,j)*dh1(ip1,j)*dh2(ip1,j)*ds2(ip1,j))

      X(+2,-1)=(0.25*sigB(1,jm1,3))/(dh2(i,jm1)*(1.5*ds1(i,jm1) - 0.5*
     &ds1(ip1,jm1))) - (2.0*gam21(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh2(ip2
     &,j)*ds2(ip2,jm1))

      X(+2,0)=((0.25*sigB(1,jm1,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - 
     &(0.25*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)))/dh2(i,jm1) -
     & (1.0*gam11(i,j))/(dh1(i,j)**2*dh1(ip1,j)*ds1(ip1,j)) + (gam11(i
     &p1,j)*((1.0*(1.0/ds1(ip1,j) + 1.0/ds1(ip2,j)))/dh1(ip2,j) + 1.0/
     &(dh1(ip1,j)*ds1(ip1,j))))/(dh1(i,j)*dh1(ip1,j)) + (2.0*gam21(ip1
     &,j)*(1.0/ds2(ip2,jm1) + 1.0/ds2(ip2,j)))/(dh1(i,j)*dh1(ip1,j)*dh
     &2(ip2,j))

      X(+2,+1)=- (0.25*sigB(1,j,3))/(dh2(i,jm1)*(1.5*ds1(i,jp1) - 0.5*
     &ds1(ip1,jp1))) - (2.0*gam21(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh2(ip2
     &,j)*ds2(ip2,j))

      X(+3,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*dh1(ip2,j)*ds1(
     &ip2,j))


      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0

!! ================== BNDY i,1  1<=i<=Nx-1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i.gt. ist) .and. (i.lt. ien)) then

!= ===============  1,i VER2
!CCCCCCCC+X-xCCCCCCCC

      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,0)=(1.0*(((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)) + (gam11(i
     &m1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(im1,j))/ds1(im1,j) 
     &+ (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j) + s
     &igB(3,im1,1)/(dh1(i,j)*ds1(im1,j)) - (0.5*sig12(im1,j))/(dh2(i,j
     &)*ds1(im1,j)) + (0.75*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j)
     & - 0.5*ds2(im1,jp1))) - (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(
     &i,j)*dh2(i,j)**2*ds1(im1,j)) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*dh
     &2(i,j)**2*ds1(im1,j))

      X(-1,+1)=(2.0*gam21(i,jp1))/(dh1(i,jp1)*dh2(i,j)*dh2(i,jp1)*ds1(
     &im1,jp1)) - (1.0*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.
     &5*ds2(im1,jp1))) - (2.0*gam22(i,j))/(dh1(i,jp1)*dh2(i,j)**2*ds1(
     &im1,jp1)) - (2.0*e6e*gam21(i,j))/(dh1(i,jp1)*dh2(i,j)**2*ds1(im1
     &,jp1)) - (0.5*sig12(im1,j))/(dh2(i,j)*ds1(im1,jp1))

      X(-1,+2)=(0.25*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.5*
     &ds2(im1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*dh2(i,j)*dh2(i,j
     &p1)*ds1(im1,jp2))

      X(0,0)=(2.0*((0.25*sig12(im1,j))/ds1(im1,j) - (0.25*sig12(i,j))/
     &ds1(i,j)))/dh2(i,j) - ((1.0*sigB(3,im1,1))/ds1(im1,j) + (1.0*sig
     &B(3,i,1))/ds1(i,j))/dh1(i,j) - (1.0*(((1.0*gam11(im1,j))/(dh1(im
     &1,j)*ds1(im1,j)) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j
     &)))/dh1(i,j))/ds1(im1,j) + (1.0*((1.0*gam11(ip1,j))/(dh1(ip1,j)*
     &ds1(i,j)) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,j
     &)))/ds1(i,j)))/dh1(i,j) + ((0.75*sigB(3,im1,3))/(1.5*ds2(i,j) - 
     &0.5*ds2(i,jp1)) - (0.75*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,j
     &p1)))/dh1(im1,j) - (2.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) - (1.0*g
     &am22(i,j))/(dh2(i,j)**2*dh2(i,jp1)*ds2(i,j)) + gam22(i,jp1)/(dh2
     &(i,j)*dh2(i,jp1)**2*ds2(i,j)) + (1.0*gam22(i,j)*(e4coefint - 2.0
     &)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*dh2(i,j)**2) - (2.0
     &*e6e*gam21(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*dh2(i
     &,j)**2)

      X(0,+1)=(2.0*((0.25*sig12(im1,j))/ds1(im1,jp1) - (0.25*sig12(i,j
     &))/ds1(i,jp1)))/dh2(i,j) - ((1.0*sigB(3,im1,3))/(1.5*ds2(i,j) - 
     &0.5*ds2(i,jp1)) - (1.0*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1)))/dh1(im1,j) + (2.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*ga
     &m22(i,j)*((2.0*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/dh1(i,jp1) +
     & (1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1)))/dh2(i,j)**2 - (gam
     &22(i,jp1)*((1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1) + 1.0/(dh2
     &(i,jp2)*ds2(i,jp1))))/(dh2(i,j)*dh2(i,jp1)) - (2.0*gam21(i,jp1)*
     &(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,jp1)*dh2(i,j)*dh2(i,
     &jp1)) + (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))
     &/(dh1(i,jp1)*dh2(i,j)**2)

      X(0,+2)=((0.25*sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 
     &(0.25*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) -
     & (1.0*gam22(i,j))/(dh2(i,j)**2*dh2(i,jp1)*ds2(i,jp1)) + (gam22(i
     &,jp1)*((1.0*(1.0/ds2(i,jp1) + 1.0/ds2(i,jp2)))/dh2(i,jp2) + 1.0/
     &(dh2(i,jp1)*ds2(i,jp1))))/(dh2(i,j)*dh2(i,jp1)) + (2.0*gam21(i,j
     &p1)*(1.0/ds1(im1,jp2) + 1.0/ds1(i,jp2)))/(dh1(i,jp2)*dh2(i,j)*dh
     &2(i,jp1))

      X(0,+3)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*dh2(i,jp2)*ds2(
     &i,jp2))

      X(+1,0)=(1.0*((1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*
     &gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/dh1(ip1,j)))/ds1(i
     &,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j)
     & + (1.0*sigB(3,i,1))/(dh1(i,j)*ds1(i,j)) + (0.5*sig12(i,j))/(dh2
     &(i,j)*ds1(i,j)) - (0.75*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j)
     & - 0.5*ds2(ip1,jp1))) - (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(
     &i,j)*dh2(i,j)**2*ds1(i,j)) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*dh2(
     &i,j)**2*ds1(i,j))

      X(+1,+1)=(0.5*sig12(i,j))/(dh2(i,j)*ds1(i,jp1)) + (1.0*sigB(3,i,
     &3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (2.0*gam2
     &2(i,j))/(dh1(i,jp1)*dh2(i,j)**2*ds1(i,jp1)) - (2.0*e6e*gam21(i,j
     &))/(dh1(i,jp1)*dh2(i,j)**2*ds1(i,jp1)) + (2.0*gam21(i,jp1))/(dh1
     &(i,jp1)*dh2(i,j)*dh2(i,jp1)*ds1(i,jp1))

      X(+1,+2)=- (0.25*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*
     &ds2(ip1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*dh2(i,j)*dh2(i,j
     &p1)*ds1(i,jp2))

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))


      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY N,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin)   
     &        .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .gt. jst) .and. (j.lt. jen)) then
      coefX00=1.0d0
       if(formulationFlagb .eq. 1) then

      X(-4,0)=(1.5*gam11(i,j))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**3

      X(-3,0)=-(gam11(i,j)*((1.0*(1.5*ds1(im2,j) - 4.5*ds1(im1,j))*(ds
     &1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1,j)))/(0.5*ds1(im2,j) - 1
     &.5*ds1(im1,j))**2 + 7.0))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**3

      X(-2,-1)=(0.5*e5c*gam21(i,j))/(dh2(im2,j)*ds2(im2,jm1)*(0.5*ds1(
     &im2,j) - 1.5*ds1(im1,j)))

      X(-2,0)=(gam11(i,j)*((((0.5*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*d
     &s1(im1,j)))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) + 4.0)*(1.5*ds1(im
     &2,j) - 4.5*ds1(im1,j))*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im
     &1,j)))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**2 - (0.5*(1.5*ds1(im4,
     &j) - 5.5*ds1(im3,j) + 6.5*ds1(im2,j) - 2.5*ds1(im1,j)))/(0.5*ds1
     &(im2,j) - 1.5*ds1(im1,j)) + 12.0))/(0.5*ds1(im2,j) - 1.5*ds1(im1
     &,j))**3 - (0.5*sigB(2,j,1))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) - 
     &(0.5*e5c*gam21(i,j)*(1.0/ds2(im2,jm1) + 1.0/ds2(im2,j)))/(dh2(im
     &2,j)*(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))



      X(-2,+1)=(0.5*e5c*gam21(i,j))/(dh2(im2,j)*ds2(im2,j)*(0.5*ds1(im
     &2,j) - 1.5*ds1(im1,j)))

      X(-1,-1)=-(2.0*e5c*gam21(i,j))/(dh2(im1,j)*ds2(im1,jm1)*(0.5*ds1
     &(im2,j) - 1.5*ds1(im1,j)))

      X(-1,0)=(2.0*sigB(2,j,1))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) - (g
     &am11(i,j)*((((2.0*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1,j))
     &)/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) + 5.0)*(1.5*ds1(im2,j) - 4.5
     &*ds1(im1,j))*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1,j)))/(0.
     &5*ds1(im2,j) - 1.5*ds1(im1,j))**2 - (2.0*(1.5*ds1(im4,j) - 5.5*d
     &s1(im3,j) + 6.5*ds1(im2,j) - 2.5*ds1(im1,j)))/(0.5*ds1(im2,j) - 
     &1.5*ds1(im1,j)) + 9.0))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**3 + (
     &2.0*e5c*gam21(i,j)*(1.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/(dh2(im1
     &,j)*(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))

      X(-1,+1)=-(2.0*e5c*gam21(i,j))/(dh2(im1,j)*ds2(im1,j)*(0.5*ds1(i
     &m2,j) - 1.5*ds1(im1,j)))

      X(0,-1)=(1.5*e5c*gam21(i,j))/(dh2(i,j)*ds2(i,jm1)*(0.5*ds1(im2,j
     &) - 1.5*ds1(im1,j))) - (0.5*sigB(2,j,3))/ds2(i,jm1)

      X(0,0)=1.0*sigB(2,j,3)*(0.5/ds2(i,jm1) - 0.5/ds2(i,j)) - (1.5*si
     &gB(2,j,1))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) + (gam11(i,j)*((((1
     &.5*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1,j)))/(0.5*ds1(im2,
     &j) - 1.5*ds1(im1,j)) + 2.0)*(1.5*ds1(im2,j) - 4.5*ds1(im1,j))*(d
     &s1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1,j)))/(0.5*ds1(im2,j) - 
     &1.5*ds1(im1,j))**2 - (1.5*(1.5*ds1(im4,j) - 5.5*ds1(im3,j) + 6.5
     &*ds1(im2,j) - 2.5*ds1(im1,j)))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))
     & + 2.5))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**3 - (1.5*e5c*gam21(i
     &,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*(0.5*ds1(im2,j) -
     & 1.5*ds1(im1,j)))

      X(0,+1)=(0.5*sigB(2,j,3))/ds2(i,j) + (1.5*e5c*gam21(i,j))/(dh2(i
     &,j)*ds2(i,j)*(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))

      else

      X(-4,0)=(1.5*gam11(i,j))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**3

      X(-3,0)=-(gam11(i,j)*((1.0*(1.5*ds1(im2,j) - 4.5*ds1(im1,j))*(ds
     &1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1,j)))/(0.5*ds1(im2,j) - 1
     &.5*ds1(im1,j))**2 + 7.0))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**3

      X(-2,-1)=(0.5*e5c*gam21(i,j))/(dh2(i,j)*ds2(i,jm1)*(0.5*ds1(im2,
     &jm1) - 1.5*ds1(im1,jm1)))

      X(-2,0)=(gam11(i,j)*((((0.5*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*d
     &s1(im1,j)))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) + 4.0)*(1.5*ds1(im
     &2,j) - 4.5*ds1(im1,j))*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im
     &1,j)))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**2 - (0.5*(1.5*ds1(im4,
     &j) - 5.5*ds1(im3,j) + 6.5*ds1(im2,j) - 2.5*ds1(im1,j)))/(0.5*ds1
     &(im2,j) - 1.5*ds1(im1,j)) + 12.0))/(0.5*ds1(im2,j) - 1.5*ds1(im1
     &,j))**3 - (e5c*gam21(i,j)*(0.5/(ds2(i,jm1)*(0.5*ds1(im2,j) - 1.5
     &*ds1(im1,j))) + 0.5/(ds2(i,j)*(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))
     &))/dh2(i,j) - (0.5*sigB(2,j,1))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)
     &)

      X(-2,+1)=(0.5*e5c*gam21(i,j))/(dh2(i,j)*ds2(i,j)*(0.5*ds1(im2,jp
     &1) - 1.5*ds1(im1,jp1)))

      X(-1,-1)=-(2.0*e5c*gam21(i,j))/(dh2(i,j)*ds2(i,jm1)*(0.5*ds1(im2
     &,jm1) - 1.5*ds1(im1,jm1)))

      X(-1,0)=(2.0*sigB(2,j,1))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) + (e
     &5c*gam21(i,j)*(2.0/(ds2(i,jm1)*(0.5*ds1(im2,j) - 1.5*ds1(im1,j))
     &) + 2.0/(ds2(i,j)*(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))))/dh2(i,j) 
     &- (gam11(i,j)*((((2.0*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1
     &,j)))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) + 5.0)*(1.5*ds1(im2,j) -
     & 4.5*ds1(im1,j))*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1,j)))
     &/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**2 - (2.0*(1.5*ds1(im4,j) - 5
     &.5*ds1(im3,j) + 6.5*ds1(im2,j) - 2.5*ds1(im1,j)))/(0.5*ds1(im2,j
     &) - 1.5*ds1(im1,j)) + 9.0))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**3

      X(-1,+1)=-(2.0*e5c*gam21(i,j))/(dh2(i,j)*ds2(i,j)*(0.5*ds1(im2,j
     &p1) - 1.5*ds1(im1,jp1)))

      X(0,-1)=(1.5*e5c*gam21(i,j))/(dh2(i,j)*ds2(i,jm1)*(0.5*ds1(im2,j
     &m1) - 1.5*ds1(im1,jm1))) - (0.5*sigB(2,j,3))/ds2(i,jm1)

      X(0,0)=1.0*sigB(2,j,3)*(0.5/ds2(i,jm1) - 0.5/ds2(i,j)) - (1.5*si
     &gB(2,j,1))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) - (e5c*gam21(i,j)*(
     &1.5/(ds2(i,jm1)*(0.5*ds1(im2,j) - 1.5*ds1(im1,j))) + 1.5/(ds2(i,
     &j)*(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))))/dh2(i,j) + (gam11(i,j)*(
     &(((1.5*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1,j)))/(0.5*ds1(
     &im2,j) - 1.5*ds1(im1,j)) + 2.0)*(1.5*ds1(im2,j) - 4.5*ds1(im1,j)
     &)*(ds1(im3,j) - 3.0*ds1(im2,j) + 2.0*ds1(im1,j)))/(0.5*ds1(im2,j
     &) - 1.5*ds1(im1,j))**2 - (1.5*(1.5*ds1(im4,j) - 5.5*ds1(im3,j) +
     & 6.5*ds1(im2,j) - 2.5*ds1(im1,j)))/(0.5*ds1(im2,j) - 1.5*ds1(im1
     &,j)) + 2.5))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j))**3

      X(0,+1)=(0.5*sigB(2,j,3))/ds2(i,j) + (1.5*e5c*gam21(i,j))/(dh2(i
     &,j)*ds2(i,j)*(0.5*ds1(im2,jp1) - 1.5*ds1(im1,jp1)))

      endif
      X(-4:4,-4:4)=-X(-4:4,-4:4)
!! ================== BNDY i,1 2<=i<=Nx-1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .gt. ist) .and. (i.lt. ien)) then
      coefX00=1.0d0
       if(formulationFlagb .eq. 1) then

      X(-1,0)=(1.5*e5c*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*(1.5*ds2(i,j) 
     &- 0.5*ds2(i,jp1))) - (0.5*sigB(3,i,3))/ds1(im1,j)

      X(-1,+1)=-(2.0*e5c*gam21(i,j))/(dh1(i,jp1)*ds1(im1,jp1)*(1.5*ds2
     &(i,j) - 0.5*ds2(i,jp1)))

      X(-1,+2)=(0.5*e5c*gam21(i,j))/(dh1(i,jp2)*ds1(im1,jp2)*(1.5*ds2(
     &i,j) - 0.5*ds2(i,jp1)))

      X(0,0)=1.0*sigB(3,i,3)*(0.5/ds1(im1,j) - 0.5/ds1(i,j)) - (1.5*si
     &gB(3,i,2))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) + (gam22(i,j)*((((1.5
     &*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0
     &.5*ds2(i,jp1)) - 2.0)*(4.5*ds2(i,j) - 1.5*ds2(i,jp1))*(2.0*ds2(i
     &,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1))**2 - (1.5*(2.5*ds2(i,j) - 6.5*ds2(i,jp1) + 5.5*ds2(i,jp2) - 
     &1.5*ds2(i,jp3)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) + 2.5))/(1.5*ds
     &2(i,j) - 0.5*ds2(i,jp1))**3 - (1.5*e5c*gam21(i,j)*(1.0/ds1(im1,j
     &) + 1.0/ds1(i,j)))/(dh1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))

      X(0,+1)=(2.0*sigB(3,i,2))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - (gam
     &22(i,j)*((((2.0*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1
     &.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 5.0)*(4.5*ds2(i,j) - 1.5*ds2(i,j
     &p1))*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j)
     & - 0.5*ds2(i,jp1))**2 - (2.0*(2.5*ds2(i,j) - 6.5*ds2(i,jp1) + 5.
     &5*ds2(i,jp2) - 1.5*ds2(i,jp3)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) 
     &+ 9.0))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**3 + (2.0*e5c*gam21(i,j)
     &*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,jp1)*(1.5*ds2(i,j) 
     &- 0.5*ds2(i,jp1)))

      X(0,+2)=(gam22(i,j)*((((0.5*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2
     &(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 4.0)*(4.5*ds2(i,j) -
     & 1.5*ds2(i,jp1))*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(
     &1.5*ds2(i,j) - 0.5*ds2(i,jp1))**2 - (0.5*(2.5*ds2(i,j) - 6.5*ds2
     &(i,jp1) + 5.5*ds2(i,jp2) - 1.5*ds2(i,jp3)))/(1.5*ds2(i,j) - 0.5*
     &ds2(i,jp1)) + 12.0))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**3 - (0.5*s
     &igB(3,i,2))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - (0.5*e5c*gam21(i,j
     &)*(1.0/ds1(im1,jp2) + 1.0/ds1(i,jp2)))/(dh1(i,jp2)*(1.5*ds2(i,j)
     & - 0.5*ds2(i,jp1)))

      X(0,+3)=(gam22(i,j)*(((4.5*ds2(i,j) - 1.5*ds2(i,jp1))*(2.0*ds2(i
     &,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1))**2 - 7.0))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**3

      X(0,+4)=(1.5*gam22(i,j))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**3

      X(+1,0)=(0.5*sigB(3,i,3))/ds1(i,j) + (1.5*e5c*gam21(i,j))/(dh1(i
     &,j)*ds1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))

      X(+1,+1)=-(2.0*e5c*gam21(i,j))/(dh1(i,jp1)*ds1(i,jp1)*(1.5*ds2(i
     &,j) - 0.5*ds2(i,jp1)))

      X(+1,+2)=(0.5*e5c*gam21(i,j))/(dh1(i,jp2)*ds1(i,jp2)*(1.5*ds2(i,
     &j) - 0.5*ds2(i,jp1)))

      else


      X(-1,0)=(1.5*e5c*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*(1.5*ds2(im1,j
     &) - 0.5*ds2(im1,jp1))) - (0.5*sigB(3,i,3))/ds1(im1,j)

      X(-1,+1)=-(2.0*e5c*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*(1.5*ds2(im1
     &,j) - 0.5*ds2(im1,jp1)))

      X(-1,+2)=(0.5*e5c*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*(1.5*ds2(im1,
     &j) - 0.5*ds2(im1,jp1)))

      X(0,0)=1.0*sigB(3,i,3)*(0.5/ds1(im1,j) - 0.5/ds1(i,j)) - (1.5*si
     &gB(3,i,2))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - (e5c*gam21(i,j)*(1.
     &5/(ds1(im1,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))) + 1.5/(ds1(i,j)*(
     &1.5*ds2(i,j) - 0.5*ds2(i,jp1)))))/dh1(i,j) + (gam22(i,j)*((((1.5
     &*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0
     &.5*ds2(i,jp1)) - 2.0)*(4.5*ds2(i,j) - 1.5*ds2(i,jp1))*(2.0*ds2(i
     &,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1))**2 - (1.5*(2.5*ds2(i,j) - 6.5*ds2(i,jp1) + 5.5*ds2(i,jp2) - 
     &1.5*ds2(i,jp3)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) + 2.5))/(1.5*ds
     &2(i,j) - 0.5*ds2(i,jp1))**3

      X(0,+1)=(2.0*sigB(3,i,2))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) + (e5c
     &*gam21(i,j)*(2.0/(ds1(im1,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))) + 
     &2.0/(ds1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))))/dh1(i,j) - (gam
     &22(i,j)*((((2.0*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1
     &.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 5.0)*(4.5*ds2(i,j) - 1.5*ds2(i,j
     &p1))*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j)
     & - 0.5*ds2(i,jp1))**2 - (2.0*(2.5*ds2(i,j) - 6.5*ds2(i,jp1) + 5.
     &5*ds2(i,jp2) - 1.5*ds2(i,jp3)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) 
     &+ 9.0))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**3

      X(0,+2)=(gam22(i,j)*((((0.5*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2
     &(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 4.0)*(4.5*ds2(i,j) -
     & 1.5*ds2(i,jp1))*(2.0*ds2(i,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(
     &1.5*ds2(i,j) - 0.5*ds2(i,jp1))**2 - (0.5*(2.5*ds2(i,j) - 6.5*ds2
     &(i,jp1) + 5.5*ds2(i,jp2) - 1.5*ds2(i,jp3)))/(1.5*ds2(i,j) - 0.5*
     &ds2(i,jp1)) + 12.0))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**3 - (e5c*g
     &am21(i,j)*(0.5/(ds1(im1,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1))) + 0.
     &5/(ds1(i,j)*(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))))/dh1(i,j) - (0.5*s
     &igB(3,i,2))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))

      X(0,+3)=(gam22(i,j)*(((4.5*ds2(i,j) - 1.5*ds2(i,jp1))*(2.0*ds2(i
     &,j) - 3.0*ds2(i,jp1) + ds2(i,jp2)))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1))**2 - 7.0))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**3

      X(0,+4)=(1.5*gam22(i,j))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1))**3

      X(+1,0)=(0.5*sigB(3,i,3))/ds1(i,j) + (1.5*e5c*gam21(i,j))/(dh1(i
     &,j)*ds1(i,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1)))

      X(+1,+1)=-(2.0*e5c*gam21(i,j))/(dh1(i,j)*ds1(i,j)*(1.5*ds2(ip1,j
     &) - 0.5*ds2(ip1,jp1)))

      X(+1,+2)=(0.5*e5c*gam21(i,j))/(dh1(i,j)*ds1(i,j)*(1.5*ds2(ip1,j)
     & - 0.5*ds2(ip1,jp1)))
      endif

      X(-4:4,-4:4)=-X(-4:4,-4:4)
!! ================== BNDY i,N 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin)   
     &        .and. (BoundaryFibShellin(1,2) .eq. 1)  
     &    .and. (i .gt. ist) .and. (i.lt. ien)) then
      coefX00=1.0d0
       if(formulationFlagb .eq. 1) then
      X(-1,-2)=(0.5*e5c*gam21(i,j))/(dh1(i,jm2)*ds1(im1,jm2)*(0.5*ds2(
     &i,jm2) - 1.5*ds2(i,jm1)))

      X(-1,-1)=-(2.0*e5c*gam21(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*(0.5*ds2
     &(i,jm2) - 1.5*ds2(i,jm1)))

      X(-1,0)=(1.5*e5c*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*(0.5*ds2(i,jm2
     &) - 1.5*ds2(i,jm1))) - (0.5*sigB(4,i,3))/ds1(im1,j)

      X(0,-4)=(1.5*gam22(i,j))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**3

      X(0,-3)=-(gam22(i,j)*((1.0*(1.5*ds2(i,jm2) - 4.5*ds2(i,jm1))*(ds
     &2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,jm1)))/(0.5*ds2(i,jm2) - 1
     &.5*ds2(i,jm1))**2 + 7.0))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**3

      X(0,-2)=(gam22(i,j)*((((0.5*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*d
     &s2(i,jm1)))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) + 4.0)*(1.5*ds2(i,
     &jm2) - 4.5*ds2(i,jm1))*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,
     &jm1)))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**2 - (0.5*(1.5*ds2(i,jm
     &4) - 5.5*ds2(i,jm3) + 6.5*ds2(i,jm2) - 2.5*ds2(i,jm1)))/(0.5*ds2
     &(i,jm2) - 1.5*ds2(i,jm1)) + 12.0))/(0.5*ds2(i,jm2) - 1.5*ds2(i,j
     &m1))**3 - (0.5*sigB(4,i,2))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) - 
     &(0.5*e5c*gam21(i,j)*(1.0/ds1(im1,jm2) + 1.0/ds1(i,jm2)))/(dh1(i,
     &jm2)*(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))

      X(0,-1)=(2.0*sigB(4,i,2))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) - (g
     &am22(i,j)*((((2.0*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,jm1))
     &)/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) + 5.0)*(1.5*ds2(i,jm2) - 4.5
     &*ds2(i,jm1))*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,jm1)))/(0.
     &5*ds2(i,jm2) - 1.5*ds2(i,jm1))**2 - (2.0*(1.5*ds2(i,jm4) - 5.5*d
     &s2(i,jm3) + 6.5*ds2(i,jm2) - 2.5*ds2(i,jm1)))/(0.5*ds2(i,jm2) - 
     &1.5*ds2(i,jm1)) + 9.0))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**3 + (
     &2.0*e5c*gam21(i,j)*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/(dh1(i,j
     &m1)*(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))

      X(0,0)=1.0*sigB(4,i,3)*(0.5/ds1(im1,j) - 0.5/ds1(i,j)) - (1.5*si
     &gB(4,i,2))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) + (gam22(i,j)*((((1
     &.5*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,jm1)))/(0.5*ds2(i,jm
     &2) - 1.5*ds2(i,jm1)) + 2.0)*(1.5*ds2(i,jm2) - 4.5*ds2(i,jm1))*(d
     &s2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,jm1)))/(0.5*ds2(i,jm2) - 
     &1.5*ds2(i,jm1))**2 - (1.5*(1.5*ds2(i,jm4) - 5.5*ds2(i,jm3) + 6.5
     &*ds2(i,jm2) - 2.5*ds2(i,jm1)))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))
     & + 2.5))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**3 - (1.5*e5c*gam21(i
     &,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*(0.5*ds2(i,jm2) -
     & 1.5*ds2(i,jm1)))

      X(+1,-2)=(0.5*e5c*gam21(i,j))/(dh1(i,jm2)*ds1(i,jm2)*(0.5*ds2(i,
     &jm2) - 1.5*ds2(i,jm1)))

      X(+1,-1)=-(2.0*e5c*gam21(i,j))/(dh1(i,jm1)*ds1(i,jm1)*(0.5*ds2(i
     &,jm2) - 1.5*ds2(i,jm1)))

      X(+1,0)=(0.5*sigB(4,i,3))/ds1(i,j) + (1.5*e5c*gam21(i,j))/(dh1(i
     &,j)*ds1(i,j)*(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))

      else

      X(-1,-2)=(0.5*e5c*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*(0.5*ds2(im1,
     &jm2) - 1.5*ds2(im1,jm1)))

      X(-1,-1)=-(2.0*e5c*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*(0.5*ds2(im1
     &,jm2) - 1.5*ds2(im1,jm1)))

      X(-1,0)=(1.5*e5c*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*(0.5*ds2(im1,j
     &m2) - 1.5*ds2(im1,jm1))) - (0.5*sigB(4,i,3))/ds1(im1,j)

      X(0,-4)=(1.5*gam22(i,j))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**3

      X(0,-3)=-(gam22(i,j)*((1.0*(1.5*ds2(i,jm2) - 4.5*ds2(i,jm1))*(ds
     &2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,jm1)))/(0.5*ds2(i,jm2) - 1
     &.5*ds2(i,jm1))**2 + 7.0))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**3

      X(0,-2)=(gam22(i,j)*((((0.5*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*d
     &s2(i,jm1)))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) + 4.0)*(1.5*ds2(i,
     &jm2) - 4.5*ds2(i,jm1))*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,
     &jm1)))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**2 - (0.5*(1.5*ds2(i,jm
     &4) - 5.5*ds2(i,jm3) + 6.5*ds2(i,jm2) - 2.5*ds2(i,jm1)))/(0.5*ds2
     &(i,jm2) - 1.5*ds2(i,jm1)) + 12.0))/(0.5*ds2(i,jm2) - 1.5*ds2(i,j
     &m1))**3 - (e5c*gam21(i,j)*(0.5/(ds1(im1,j)*(0.5*ds2(i,jm2) - 1.5
     &*ds2(i,jm1))) + 0.5/(ds1(i,j)*(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))
     &))/dh1(i,j) - (0.5*sigB(4,i,2))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)
     &)

      X(0,-1)=(2.0*sigB(4,i,2))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) + (e
     &5c*gam21(i,j)*(2.0/(ds1(im1,j)*(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))
     &) + 2.0/(ds1(i,j)*(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))))/dh1(i,j) 
     &- (gam22(i,j)*((((2.0*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,j
     &m1)))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) + 5.0)*(1.5*ds2(i,jm2) -
     & 4.5*ds2(i,jm1))*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,jm1)))
     &/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**2 - (2.0*(1.5*ds2(i,jm4) - 5
     &.5*ds2(i,jm3) + 6.5*ds2(i,jm2) - 2.5*ds2(i,jm1)))/(0.5*ds2(i,jm2
     &) - 1.5*ds2(i,jm1)) + 9.0))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**3

      X(0,0)=1.0*sigB(4,i,3)*(0.5/ds1(im1,j) - 0.5/ds1(i,j)) - (1.5*si
     &gB(4,i,2))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) - (e5c*gam21(i,j)*(
     &1.5/(ds1(im1,j)*(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))) + 1.5/(ds1(i,
     &j)*(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))))/dh1(i,j) + (gam22(i,j)*(
     &(((1.5*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,jm1)))/(0.5*ds2(
     &i,jm2) - 1.5*ds2(i,jm1)) + 2.0)*(1.5*ds2(i,jm2) - 4.5*ds2(i,jm1)
     &)*(ds2(i,jm3) - 3.0*ds2(i,jm2) + 2.0*ds2(i,jm1)))/(0.5*ds2(i,jm2
     &) - 1.5*ds2(i,jm1))**2 - (1.5*(1.5*ds2(i,jm4) - 5.5*ds2(i,jm3) +
     & 6.5*ds2(i,jm2) - 2.5*ds2(i,jm1)))/(0.5*ds2(i,jm2) - 1.5*ds2(i,j
     &m1)) + 2.5))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1))**3

      X(+1,-2)=(0.5*e5c*gam21(i,j))/(dh1(i,j)*ds1(i,j)*(0.5*ds2(ip1,jm
     &2) - 1.5*ds2(ip1,jm1)))

      X(+1,-1)=-(2.0*e5c*gam21(i,j))/(dh1(i,j)*ds1(i,j)*(0.5*ds2(ip1,j
     &m2) - 1.5*ds2(ip1,jm1)))

      X(+1,0)=(0.5*sigB(4,i,3))/ds1(i,j) + (1.5*e5c*gam21(i,j))/(dh1(i
     &,j)*ds1(i,j)*(0.5*ds2(ip1,jm2) - 1.5*ds2(ip1,jm1)))

      endif
      X(-4:4,-4:4)=-X(-4:4,-4:4)
!! ================== BNDY 2,j  2<=j<=Ny-1
       elseif((i.eq. 2) .and.
     &    ((BoundaryFibShellin(2,1) .eq. 1) 
     &     .or. (BoundaryFibShellin(2,1) .eq. 2) )
     &    .and. (j .ge. jst) .and. (j.le. jen)) then

      X(-1,-1)=(0.25*sig12(im1,jm1))/(dh1(i,j)*ds2(im1,jm1)) + (0.25*s
     &ig21(im1,jm1))/(dh2(i,j)*ds1(im1,jm1)) - (1.0*gam12(im1,jm1))/(d
     &g1(im1,jm1)*dh1(i,j)*dk2(im1,j)*ds2(im1,jm1)) - (1.0*gam21(im1,j
     &m1))/(dg2(im1,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(im1,jm1)) + (e4coefin
     &t*gam11(im1,j))/(dh1(i,j)*dh2(im1,j)*ds1(im1,j)*ds2(im1,jm1))

      X(-1,0)=(1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*dk1(i,jm1)*ds1(
     &im1,j)) + (1.0*gam21(im1,j))/(dg2(im1,j)*dk1(i,j)*ds1(im1,j))))/
     &dh2(i,j) + ((gam11(i,j)/(dh1(i,j)*ds1(im1,j)) - (e4coefint*gam11
     &(im1,j)*(1.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/dh2(im1,j))/ds1(im1
     &,j) + gam11(i,j)/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) - ((0.
     &25*sig12(im1,jm1))/ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(im1,j)
     &)/dh1(i,j) + ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.25*sig21(im1
     &,j))/ds1(im1,j))/dh2(i,j) + sig11(im1,j)/(dh1(i,j)*ds1(im1,j)) +
     & (1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(im1,jm1)) + (1.0*g
     &am12(im1,j))/(dg1(im1,j)*ds2(im1,j))))/(dh1(i,j)*dk2(im1,j))

      X(-1,+1)=(1.0*e4coefint*gam11(im1,j))/(dh1(i,j)*dh2(im1,j)*ds1(i
     &m1,j)*ds2(im1,j)) - (0.25*sig21(im1,j))/(dh2(i,j)*ds1(im1,jp1)) 
     &- (1.0*gam12(im1,j))/(dg1(im1,j)*dh1(i,j)*dk2(im1,j)*ds2(im1,j))
     & - (1.0*gam21(im1,j))/(dg2(im1,j)*dh2(i,j)*dk1(i,j)*ds1(im1,jp1)
     &) - (0.25*sig12(im1,j))/(dh1(i,j)*ds2(im1,j))

      X(0,-2)=-(1.0*gam22(i,jm1))/(dh2(i,jm1)*dh2(i,j)*ds2(i,jm2)*ds2(
     &i,jm1))

      X(0,-1)=(1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*dk2(im1,j)*ds2(
     &i,jm1)) + (1.0*gam12(i,jm1))/(dg1(i,jm1)*dk2(i,j)*ds2(i,jm1))))/
     &dh1(i,j) + (1.0*(((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)) + (gam2
     &2(i,jm1)*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i,jm1))/ds2(i,jm
     &1) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) 
     &+ ((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(i,jm1))/ds2(i,
     &jm1))/dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,jm1) - (0.25*sig
     &21(i,jm1))/ds1(i,jm1))/dh2(i,j) + sig22(i,jm1)/(dh2(i,j)*ds2(i,j
     &m1)) + (1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*ds1(im1,jm1)) + 
     &(1.0*gam21(i,jm1))/(dg2(i,jm1)*ds1(i,jm1))))/(dh2(i,j)*dk1(i,jm1
     &))

      X(0,0)=- (1.0*(((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(i,jm1)) +
     & (1.0*gam12(im1,j))/(dg1(im1,j)*ds2(i,j)))/dk2(im1,j) + (1.0*((1
     &.0*gam12(i,jm1))/(dg1(i,jm1)*ds2(i,jm1)) + (1.0*gam12(i,j))/(dg1
     &(i,j)*ds2(i,j))))/dk2(i,j)))/dh1(i,j) - (1.0*(((1.0*gam21(im1,jm
     &1))/(dg2(im1,jm1)*ds1(im1,j)) + (1.0*gam21(i,jm1))/(dg2(i,jm1)*d
     &s1(i,j)))/dk1(i,jm1) + (1.0*((1.0*gam21(im1,j))/(dg2(im1,j)*ds1(
     &im1,j)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j))))/dk1(i,j)))/dh2(
     &i,j) - ((1.0*sig11(im1,j))/ds1(im1,j) + (1.0*sig11(i,j))/ds1(i,j
     &))/dh1(i,j) - ((1.0*sig22(i,jm1))/ds2(i,jm1) + (1.0*sig22(i,j))/
     &ds2(i,j))/dh2(i,j) - (1.0*(((1.0*gam22(i,jm1))/(dh2(i,jm1)*ds2(i
     &,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,
     &j))/ds2(i,jm1) + (1.0*((1.0*gam22(i,jp1))/(dh2(i,jp1)*ds2(i,j)) 
     &+ (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j)))/ds2(i,
     &j)))/dh2(i,j) - (((1.0*gam11(ip1,j))/(dh1(ip1,j)*ds1(i,j)) + (ga
     &m11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,j))/ds1(i,j) + (
     &gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds1(im1,j)
     &))/dh1(i,j) - ((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(im
     &1,j))/ds2(i,j) - (0.25*sig12(i,jm1))/ds2(i,jm1) + (0.25*sig12(i,
     &j))/ds2(i,j))/dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.
     &25*sig21(im1,j))/ds1(im1,j) - (0.25*sig21(i,jm1))/ds1(i,j) + (0.
     &25*sig21(i,j))/ds1(i,j))/dh2(i,j)

      X(0,+1)=(1.0*((1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*
     &gam22(i,jp1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1)))/ds2(i
     &,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j)
     & + (1.0*((1.0*gam12(im1,j))/(dg1(im1,j)*dk2(im1,j)*ds2(i,j)) + (
     &1.0*gam12(i,j))/(dg1(i,j)*dk2(i,j)*ds2(i,j))))/dh1(i,j) - ((0.25
     &*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,j))/ds2(i,j))/dh1(i,j) +
     & ((0.25*sig21(im1,j))/ds1(im1,jp1) - (0.25*sig21(i,j))/ds1(i,jp1
     &))/dh2(i,j) + (1.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*((1.0*
     &gam21(im1,j))/(dg2(im1,j)*ds1(im1,jp1)) + (1.0*gam21(i,j))/(dg2(
     &i,j)*ds1(i,jp1))))/(dh2(i,j)*dk1(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(+1,-1)=- (0.25*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1)) - (0.25*s
     &ig21(i,jm1))/(dh2(i,j)*ds1(i,jm1)) - (1.0*gam12(i,jm1))/(dg1(i,j
     &m1)*dh1(i,j)*dk2(i,j)*ds2(ip1,jm1)) - (1.0*gam21(i,jm1))/(dg2(i,
     &jm1)*dh2(i,j)*dk1(i,jm1)*ds1(i,jm1))

      X(+1,0)=(1.0*((1.0*gam21(i,jm1))/(dg2(i,jm1)*dk1(i,jm1)*ds1(i,j)
     &) + (1.0*gam21(i,j))/(dg2(i,j)*dk1(i,j)*ds1(i,j))))/dh2(i,j) + (
     &(0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j))/ds2(ip1,j))
     &/dh1(i,j) - ((0.25*sig21(i,jm1))/ds1(i,j) - (0.25*sig21(i,j))/ds
     &1(i,j))/dh2(i,j) + (((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0
     &*gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/dh1(ip1,j))/ds1(i
     &,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) 
     &+ (1.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*((1.0*gam12(i,jm1)
     &)/(dg1(i,jm1)*ds2(ip1,jm1)) + (1.0*gam12(i,j))/(dg1(i,j)*ds2(ip1
     &,j))))/(dh1(i,j)*dk2(i,j))

      X(+1,+1)=(0.25*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (0.25*sig21(i
     &,j))/(dh2(i,j)*ds1(i,jp1)) - (1.0*gam12(i,j))/(dg1(i,j)*dh1(i,j)
     &*dk2(i,j)*ds2(ip1,j)) - (1.0*gam21(i,j))/(dg2(i,j)*dh2(i,j)*dk1(
     &i,j)*ds1(i,jp1))

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY N-1,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin-1) .and.
     &    ((BoundaryFibShellin(2,2) .eq. 1) 
     &     .or. (BoundaryFibShellin(2,2) .eq. 2) )  
     &    .and. (j .gt. jst) .and. (j.lt. jen)) then
      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,-1)=(0.25*sig12(im1,jm1))/(dh1(i,j)*ds2(im1,jm1)) + (0.25*s
     &ig21(im1,jm1))/(dh2(i,j)*ds1(im1,jm1)) - (1.0*gam12(im1,jm1))/(d
     &g1(im1,jm1)*dh1(i,j)*dk2(im1,j)*ds2(im1,jm1)) - (1.0*gam21(im1,j
     &m1))/(dg2(im1,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(im1,jm1))

      X(-1,0)=(1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*dk1(i,jm1)*ds1(
     &im1,j)) + (1.0*gam21(im1,j))/(dg2(im1,j)*dk1(i,j)*ds1(im1,j))))/
     &dh2(i,j) + ((1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)) + (gam1
     &1(im1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(im1,j)))/ds1(im1
     &,j) + gam11(i,j)/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) - ((0.
     &25*sig12(im1,jm1))/ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(im1,j)
     &)/dh1(i,j) + ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.25*sig21(im1
     &,j))/ds1(im1,j))/dh2(i,j) + sig11(im1,j)/(dh1(i,j)*ds1(im1,j)) +
     & (1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(im1,jm1)) + (1.0*g
     &am12(im1,j))/(dg1(im1,j)*ds2(im1,j))))/(dh1(i,j)*dk2(im1,j))

      X(-1,+1)=- (0.25*sig12(im1,j))/(dh1(i,j)*ds2(im1,j)) - (0.25*sig
     &21(im1,j))/(dh2(i,j)*ds1(im1,jp1)) - (1.0*gam12(im1,j))/(dg1(im1
     &,j)*dh1(i,j)*dk2(im1,j)*ds2(im1,j)) - (1.0*gam21(im1,j))/(dg2(im
     &1,j)*dh2(i,j)*dk1(i,j)*ds1(im1,jp1))

      X(0,-2)=-(1.0*gam22(i,jm1))/(dh2(i,jm1)*dh2(i,j)*ds2(i,jm2)*ds2(
     &i,jm1))

      X(0,-1)=(1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*dk2(im1,j)*ds2(
     &i,jm1)) + (1.0*gam12(i,jm1))/(dg1(i,jm1)*dk2(i,j)*ds2(i,jm1))))/
     &dh1(i,j) + (1.0*(((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)) + (gam2
     &2(i,jm1)*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i,jm1))/ds2(i,jm
     &1) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) 
     &+ ((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(i,jm1))/ds2(i,
     &jm1))/dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,jm1) - (0.25*sig
     &21(i,jm1))/ds1(i,jm1))/dh2(i,j) + sig22(i,jm1)/(dh2(i,j)*ds2(i,j
     &m1)) + (1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*ds1(im1,jm1)) + 
     &(1.0*gam21(i,jm1))/(dg2(i,jm1)*ds1(i,jm1))))/(dh2(i,j)*dk1(i,jm1
     &))

      X(0,0)=- (1.0*(((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(i,jm1)) +
     & (1.0*gam12(im1,j))/(dg1(im1,j)*ds2(i,j)))/dk2(im1,j) + (1.0*((1
     &.0*gam12(i,jm1))/(dg1(i,jm1)*ds2(i,jm1)) + (1.0*gam12(i,j))/(dg1
     &(i,j)*ds2(i,j))))/dk2(i,j)))/dh1(i,j) - (1.0*(((1.0*gam21(im1,jm
     &1))/(dg2(im1,jm1)*ds1(im1,j)) + (1.0*gam21(i,jm1))/(dg2(i,jm1)*d
     &s1(i,j)))/dk1(i,jm1) + (1.0*((1.0*gam21(im1,j))/(dg2(im1,j)*ds1(
     &im1,j)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j))))/dk1(i,j)))/dh2(
     &i,j) - ((1.0*sig11(im1,j))/ds1(im1,j) + (1.0*sig11(i,j))/ds1(i,j
     &))/dh1(i,j) - ((1.0*sig22(i,jm1))/ds2(i,jm1) + (1.0*sig22(i,j))/
     &ds2(i,j))/dh2(i,j) - ((1.0*((1.0*gam11(im1,j))/(dh1(im1,j)*ds1(i
     &m1,j)) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,
     &j)))/ds1(im1,j) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(
     &dh1(i,j)*ds1(i,j)))/dh1(i,j) - (1.0*(((1.0*gam22(i,jm1))/(dh2(i,
     &jm1)*ds2(i,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j
     &)))/dh2(i,j))/ds2(i,jm1) + (1.0*((1.0*gam22(i,jp1))/(dh2(i,jp1)*
     &ds2(i,j)) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j
     &)))/ds2(i,j)))/dh2(i,j) - ((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0
     &.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,jm1))/ds2(i,jm1) + (0
     &.25*sig12(i,j))/ds2(i,j))/dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(
     &im1,j) - (0.25*sig21(im1,j))/ds1(im1,j) - (0.25*sig21(i,jm1))/ds
     &1(i,j) + (0.25*sig21(i,j))/ds1(i,j))/dh2(i,j)

      X(0,+1)=(1.0*((1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*
     &gam22(i,jp1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1)))/ds2(i
     &,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j)
     & + (1.0*((1.0*gam12(im1,j))/(dg1(im1,j)*dk2(im1,j)*ds2(i,j)) + (
     &1.0*gam12(i,j))/(dg1(i,j)*dk2(i,j)*ds2(i,j))))/dh1(i,j) - ((0.25
     &*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,j))/ds2(i,j))/dh1(i,j) +
     & ((0.25*sig21(im1,j))/ds1(im1,jp1) - (0.25*sig21(i,j))/ds1(i,jp1
     &))/dh2(i,j) + (1.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*((1.0*
     &gam21(im1,j))/(dg2(im1,j)*ds1(im1,jp1)) + (1.0*gam21(i,j))/(dg2(
     &i,j)*ds1(i,jp1))))/(dh2(i,j)*dk1(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(+1,-1)=(e4coefint*gam11(ip1,j))/(dh1(i,j)*dh2(ip1,j)*ds1(i,j)*
     &ds2(ip1,jm1)) - (0.25*sig21(i,jm1))/(dh2(i,j)*ds1(i,jm1)) - (1.0
     &*gam12(i,jm1))/(dg1(i,jm1)*dh1(i,j)*dk2(i,j)*ds2(ip1,jm1)) - (1.
     &0*gam21(i,jm1))/(dg2(i,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(i,jm1)) - (0
     &.25*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1))

      X(+1,0)=(1.0*((1.0*gam21(i,jm1))/(dg2(i,jm1)*dk1(i,jm1)*ds1(i,j)
     &) + (1.0*gam21(i,j))/(dg2(i,j)*dk1(i,j)*ds1(i,j))))/dh2(i,j) + (
     &(0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j))/ds2(ip1,j))
     &/dh1(i,j) - ((0.25*sig21(i,jm1))/ds1(i,j) - (0.25*sig21(i,j))/ds
     &1(i,j))/dh2(i,j) + (((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) - (e4c
     &oefint*gam11(ip1,j)*(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/dh2(ip1
     &,j))/ds1(i,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))
     &/dh1(i,j) + (1.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*((1.0*ga
     &m12(i,jm1))/(dg1(i,jm1)*ds2(ip1,jm1)) + (1.0*gam12(i,j))/(dg1(i,
     &j)*ds2(ip1,j))))/(dh1(i,j)*dk2(i,j))

      X(+1,+1)=(0.25*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (0.25*sig21(i
     &,j))/(dh2(i,j)*ds1(i,jp1)) - (1.0*gam12(i,j))/(dg1(i,j)*dh1(i,j)
     &*dk2(i,j)*ds2(ip1,j)) - (1.0*gam21(i,j))/(dg2(i,j)*dh2(i,j)*dk1(
     &i,j)*ds1(i,jp1)) + (1.0*e4coefint*gam11(ip1,j))/(dh1(i,j)*dh2(ip
     &1,j)*ds1(i,j)*ds2(ip1,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY i,2 2<=i<=Nx-1
       elseif((j.eq. 2) .and.  
     &    ((BoundaryFibShellin(1,1) .eq. 1) 
     &     .or. (BoundaryFibShellin(1,1) .eq. 2) )
     &    .and. (i .gt. ist) .and. (i.lt. ien)) then
      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,-1)=(0.25*sig12(im1,jm1))/(dh1(i,j)*ds2(im1,jm1)) + (0.25*s
     &ig21(im1,jm1))/(dh2(i,j)*ds1(im1,jm1)) - (1.0*gam12(im1,jm1))/(d
     &g1(im1,jm1)*dh1(i,j)*dk2(im1,j)*ds2(im1,jm1)) - (1.0*gam21(im1,j
     &m1))/(dg2(im1,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(im1,jm1)) + (e4coefin
     &t*gam22(i,jm1))/(dh1(i,jm1)*dh2(i,j)*ds1(im1,jm1)*ds2(i,jm1))

      X(-1,0)=(1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*dk1(i,jm1)*ds1(
     &im1,j)) + (1.0*gam21(im1,j))/(dg2(im1,j)*dk1(i,j)*ds1(im1,j))))/
     &dh2(i,j) + (1.0*(((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)) + (gam1
     &1(im1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(im1,j))/ds1(im1,
     &j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j) 
     &- ((0.25*sig12(im1,jm1))/ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(
     &im1,j))/dh1(i,j) + ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.25*sig
     &21(im1,j))/ds1(im1,j))/dh2(i,j) + sig11(im1,j)/(dh1(i,j)*ds1(im1
     &,j)) + (1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(im1,jm1)) + 
     &(1.0*gam12(im1,j))/(dg1(im1,j)*ds2(im1,j))))/(dh1(i,j)*dk2(im1,j
     &))

      X(-1,+1)=- (0.25*sig12(im1,j))/(dh1(i,j)*ds2(im1,j)) - (0.25*sig
     &21(im1,j))/(dh2(i,j)*ds1(im1,jp1)) - (1.0*gam12(im1,j))/(dg1(im1
     &,j)*dh1(i,j)*dk2(im1,j)*ds2(im1,j)) - (1.0*gam21(im1,j))/(dg2(im
     &1,j)*dh2(i,j)*dk1(i,j)*ds1(im1,jp1))

      X(0,-1)=(1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*dk2(im1,j)*ds2(
     &i,jm1)) + (1.0*gam12(i,jm1))/(dg1(i,jm1)*dk2(i,j)*ds2(i,jm1))))/
     &dh1(i,j) + ((gam22(i,j)/(dh2(i,j)*ds2(i,jm1)) - (e4coefint*gam22
     &(i,jm1)*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/dh1(i,jm1))/ds2(i,j
     &m1) + gam22(i,j)/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + ((0.
     &25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(i,jm1))/ds2(i,jm1))/
     &dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,jm1) - (0.25*sig21(i,j
     &m1))/ds1(i,jm1))/dh2(i,j) + sig22(i,jm1)/(dh2(i,j)*ds2(i,jm1)) +
     & (1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*ds1(im1,jm1)) + (1.0*g
     &am21(i,jm1))/(dg2(i,jm1)*ds1(i,jm1))))/(dh2(i,j)*dk1(i,jm1))

      X(0,0)=- (1.0*(((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(i,jm1)) +
     & (1.0*gam12(im1,j))/(dg1(im1,j)*ds2(i,j)))/dk2(im1,j) + (1.0*((1
     &.0*gam12(i,jm1))/(dg1(i,jm1)*ds2(i,jm1)) + (1.0*gam12(i,j))/(dg1
     &(i,j)*ds2(i,j))))/dk2(i,j)))/dh1(i,j) - (1.0*(((1.0*gam21(im1,jm
     &1))/(dg2(im1,jm1)*ds1(im1,j)) + (1.0*gam21(i,jm1))/(dg2(i,jm1)*d
     &s1(i,j)))/dk1(i,jm1) + (1.0*((1.0*gam21(im1,j))/(dg2(im1,j)*ds1(
     &im1,j)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j))))/dk1(i,j)))/dh2(
     &i,j) - ((1.0*sig11(im1,j))/ds1(im1,j) + (1.0*sig11(i,j))/ds1(i,j
     &))/dh1(i,j) - ((1.0*sig22(i,jm1))/ds2(i,jm1) + (1.0*sig22(i,j))/
     &ds2(i,j))/dh2(i,j) - (1.0*(((1.0*gam11(im1,j))/(dh1(im1,j)*ds1(i
     &m1,j)) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,
     &j))/ds1(im1,j) + (1.0*((1.0*gam11(ip1,j))/(dh1(ip1,j)*ds1(i,j)) 
     &+ (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,j)))/ds1(i,
     &j)))/dh1(i,j) - (((1.0*gam22(i,jp1))/(dh2(i,jp1)*ds2(i,j)) + (ga
     &m22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j))/ds2(i,j) + (
     &gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds2(i,jm1)
     &))/dh2(i,j) - ((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(im
     &1,j))/ds2(i,j) - (0.25*sig12(i,jm1))/ds2(i,jm1) + (0.25*sig12(i,
     &j))/ds2(i,j))/dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.
     &25*sig21(im1,j))/ds1(im1,j) - (0.25*sig21(i,jm1))/ds1(i,j) + (0.
     &25*sig21(i,j))/ds1(i,j))/dh2(i,j)

      X(0,+1)=(1.0*((1.0*gam12(im1,j))/(dg1(im1,j)*dk2(im1,j)*ds2(i,j)
     &) + (1.0*gam12(i,j))/(dg1(i,j)*dk2(i,j)*ds2(i,j))))/dh1(i,j) - (
     &(0.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,j))/ds2(i,j))/dh1(i
     &,j) + ((0.25*sig21(im1,j))/ds1(im1,jp1) - (0.25*sig21(i,j))/ds1(
     &i,jp1))/dh2(i,j) + (((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0
     &*gam22(i,jp1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1))/ds2(i
     &,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) 
     &+ (1.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*((1.0*gam21(im1,j)
     &)/(dg2(im1,j)*ds1(im1,jp1)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j
     &p1))))/(dh2(i,j)*dk1(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(+1,-1)=(1.0*e4coefint*gam22(i,jm1))/(dh1(i,jm1)*dh2(i,j)*ds1(i
     &,jm1)*ds2(i,jm1)) - (0.25*sig21(i,jm1))/(dh2(i,j)*ds1(i,jm1)) - 
     &(1.0*gam12(i,jm1))/(dg1(i,jm1)*dh1(i,j)*dk2(i,j)*ds2(ip1,jm1)) -
     & (1.0*gam21(i,jm1))/(dg2(i,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(i,jm1)) 
     &- (0.25*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1))

      X(+1,0)=(1.0*((1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*
     &gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/dh1(ip1,j)))/ds1(i
     &,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j)
     & + (1.0*((1.0*gam21(i,jm1))/(dg2(i,jm1)*dk1(i,jm1)*ds1(i,j)) + (
     &1.0*gam21(i,j))/(dg2(i,j)*dk1(i,j)*ds1(i,j))))/dh2(i,j) + ((0.25
     &*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j))/ds2(ip1,j))/dh1(
     &i,j) - ((0.25*sig21(i,jm1))/ds1(i,j) - (0.25*sig21(i,j))/ds1(i,j
     &))/dh2(i,j) + (1.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*((1.0*
     &gam12(i,jm1))/(dg1(i,jm1)*ds2(ip1,jm1)) + (1.0*gam12(i,j))/(dg1(
     &i,j)*ds2(ip1,j))))/(dh1(i,j)*dk2(i,j))

      X(+1,+1)=(0.25*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (0.25*sig21(i
     &,j))/(dh2(i,j)*ds1(i,jp1)) - (1.0*gam12(i,j))/(dg1(i,j)*dh1(i,j)
     &*dk2(i,j)*ds2(ip1,j)) - (1.0*gam21(i,j))/(dg2(i,j)*dh2(i,j)*dk1(
     &i,j)*ds1(i,jp1))

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0

!! ================== BNDY i,N-1 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin-1)  .and.
     &    ((BoundaryFibShellin(1,2) .eq. 1) 
     &     .or. (BoundaryFibShellin(1,2) .eq. 2) )
     &    .and. (i .gt. ist) .and. (i.lt. ien)) then
      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,-1)=(0.25*sig12(im1,jm1))/(dh1(i,j)*ds2(im1,jm1)) + (0.25*s
     &ig21(im1,jm1))/(dh2(i,j)*ds1(im1,jm1)) - (1.0*gam12(im1,jm1))/(d
     &g1(im1,jm1)*dh1(i,j)*dk2(im1,j)*ds2(im1,jm1)) - (1.0*gam21(im1,j
     &m1))/(dg2(im1,jm1)*dh2(i,j)*dk1(i,jm1)*ds1(im1,jm1))

      X(-1,0)=(1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*dk1(i,jm1)*ds1(
     &im1,j)) + (1.0*gam21(im1,j))/(dg2(im1,j)*dk1(i,j)*ds1(im1,j))))/
     &dh2(i,j) + (1.0*(((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)) + (gam1
     &1(im1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(im1,j))/ds1(im1,
     &j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j) 
     &- ((0.25*sig12(im1,jm1))/ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(
     &im1,j))/dh1(i,j) + ((0.25*sig21(im1,jm1))/ds1(im1,j) - (0.25*sig
     &21(im1,j))/ds1(im1,j))/dh2(i,j) + sig11(im1,j)/(dh1(i,j)*ds1(im1
     &,j)) + (1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(im1,jm1)) + 
     &(1.0*gam12(im1,j))/(dg1(im1,j)*ds2(im1,j))))/(dh1(i,j)*dk2(im1,j
     &))

      X(-1,+1)=(e4coefint*gam22(i,jp1))/(dh1(i,jp1)*dh2(i,j)*ds1(im1,j
     &p1)*ds2(i,j)) - (0.25*sig21(im1,j))/(dh2(i,j)*ds1(im1,jp1)) - (1
     &.0*gam12(im1,j))/(dg1(im1,j)*dh1(i,j)*dk2(im1,j)*ds2(im1,j)) - (
     &1.0*gam21(im1,j))/(dg2(im1,j)*dh2(i,j)*dk1(i,j)*ds1(im1,jp1)) - 
     &(0.25*sig12(im1,j))/(dh1(i,j)*ds2(im1,j))

      X(0,-2)=-(1.0*gam22(i,jm1))/(dh2(i,jm1)*dh2(i,j)*ds2(i,jm2)*ds2(
     &i,jm1))

      X(0,-1)=(1.0*((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*dk2(im1,j)*ds2(
     &i,jm1)) + (1.0*gam12(i,jm1))/(dg1(i,jm1)*dk2(i,j)*ds2(i,jm1))))/
     &dh1(i,j) + ((1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)) + (gam2
     &2(i,jm1)*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i,jm1)))/ds2(i,j
     &m1) + gam22(i,j)/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + ((0.
     &25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(i,jm1))/ds2(i,jm1))/
     &dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(im1,jm1) - (0.25*sig21(i,j
     &m1))/ds1(i,jm1))/dh2(i,j) + sig22(i,jm1)/(dh2(i,j)*ds2(i,jm1)) +
     & (1.0*((1.0*gam21(im1,jm1))/(dg2(im1,jm1)*ds1(im1,jm1)) + (1.0*g
     &am21(i,jm1))/(dg2(i,jm1)*ds1(i,jm1))))/(dh2(i,j)*dk1(i,jm1))

      X(0,0)=- (1.0*(((1.0*gam12(im1,jm1))/(dg1(im1,jm1)*ds2(i,jm1)) +
     & (1.0*gam12(im1,j))/(dg1(im1,j)*ds2(i,j)))/dk2(im1,j) + (1.0*((1
     &.0*gam12(i,jm1))/(dg1(i,jm1)*ds2(i,jm1)) + (1.0*gam12(i,j))/(dg1
     &(i,j)*ds2(i,j))))/dk2(i,j)))/dh1(i,j) - (1.0*(((1.0*gam21(im1,jm
     &1))/(dg2(im1,jm1)*ds1(im1,j)) + (1.0*gam21(i,jm1))/(dg2(i,jm1)*d
     &s1(i,j)))/dk1(i,jm1) + (1.0*((1.0*gam21(im1,j))/(dg2(im1,j)*ds1(
     &im1,j)) + (1.0*gam21(i,j))/(dg2(i,j)*ds1(i,j))))/dk1(i,j)))/dh2(
     &i,j) - ((1.0*sig11(im1,j))/ds1(im1,j) + (1.0*sig11(i,j))/ds1(i,j
     &))/dh1(i,j) - ((1.0*sig22(i,jm1))/ds2(i,jm1) + (1.0*sig22(i,j))/
     &ds2(i,j))/dh2(i,j) - ((1.0*((1.0*gam22(i,jm1))/(dh2(i,jm1)*ds2(i
     &,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,
     &j)))/ds2(i,jm1) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(
     &dh2(i,j)*ds2(i,j)))/dh2(i,j) - (1.0*(((1.0*gam11(im1,j))/(dh1(im
     &1,j)*ds1(im1,j)) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j
     &)))/dh1(i,j))/ds1(im1,j) + (1.0*((1.0*gam11(ip1,j))/(dh1(ip1,j)*
     &ds1(i,j)) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,j
     &)))/ds1(i,j)))/dh1(i,j) - ((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0
     &.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,jm1))/ds2(i,jm1) + (0
     &.25*sig12(i,j))/ds2(i,j))/dh1(i,j) - ((0.25*sig21(im1,jm1))/ds1(
     &im1,j) - (0.25*sig21(im1,j))/ds1(im1,j) - (0.25*sig21(i,jm1))/ds
     &1(i,j) + (0.25*sig21(i,j))/ds1(i,j))/dh2(i,j)

      X(0,+1)=(1.0*((1.0*gam12(im1,j))/(dg1(im1,j)*dk2(im1,j)*ds2(i,j)
     &) + (1.0*gam12(i,j))/(dg1(i,j)*dk2(i,j)*ds2(i,j))))/dh1(i,j) - (
     &(0.25*sig12(im1,j))/ds2(i,j) - (0.25*sig12(i,j))/ds2(i,j))/dh1(i
     &,j) + ((0.25*sig21(im1,j))/ds1(im1,jp1) - (0.25*sig21(i,j))/ds1(
     &i,jp1))/dh2(i,j) + (((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) - (e4c
     &oefint*gam22(i,jp1)*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/dh1(i,j
     &p1))/ds2(i,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))
     &/dh2(i,j) + (1.0*sig22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*((1.0*ga
     &m21(im1,j))/(dg2(im1,j)*ds1(im1,jp1)) + (1.0*gam21(i,j))/(dg2(i,
     &j)*ds1(i,jp1))))/(dh2(i,j)*dk1(i,j))

      X(+1,-1)=- (0.25*sig12(i,jm1))/(dh1(i,j)*ds2(ip1,jm1)) - (0.25*s
     &ig21(i,jm1))/(dh2(i,j)*ds1(i,jm1)) - (1.0*gam12(i,jm1))/(dg1(i,j
     &m1)*dh1(i,j)*dk2(i,j)*ds2(ip1,jm1)) - (1.0*gam21(i,jm1))/(dg2(i,
     &jm1)*dh2(i,j)*dk1(i,jm1)*ds1(i,jm1))

      X(+1,0)=(1.0*((1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*
     &gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/dh1(ip1,j)))/ds1(i
     &,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j)
     & + (1.0*((1.0*gam21(i,jm1))/(dg2(i,jm1)*dk1(i,jm1)*ds1(i,j)) + (
     &1.0*gam21(i,j))/(dg2(i,j)*dk1(i,j)*ds1(i,j))))/dh2(i,j) + ((0.25
     &*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j))/ds2(ip1,j))/dh1(
     &i,j) - ((0.25*sig21(i,jm1))/ds1(i,j) - (0.25*sig21(i,j))/ds1(i,j
     &))/dh2(i,j) + (1.0*sig11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*((1.0*
     &gam12(i,jm1))/(dg1(i,jm1)*ds2(ip1,jm1)) + (1.0*gam12(i,j))/(dg1(
     &i,j)*ds2(ip1,j))))/(dh1(i,j)*dk2(i,j))

      X(+1,+1)=(0.25*sig12(i,j))/(dh1(i,j)*ds2(ip1,j)) + (0.25*sig21(i
     &,j))/(dh2(i,j)*ds1(i,jp1)) - (1.0*gam12(i,j))/(dg1(i,j)*dh1(i,j)
     &*dk2(i,j)*ds2(ip1,j)) - (1.0*gam21(i,j))/(dg2(i,j)*dh2(i,j)*dk1(
     &i,j)*ds1(i,jp1)) + (1.0*e4coefint*gam22(i,jp1))/(dh1(i,jp1)*dh2(
     &i,j)*ds1(i,jp1)*ds2(i,j))

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0

!other points
       else
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
       X(-2,0)=
     & -gam11(im1,j)/ (dh1(im1,j)* dh1(i,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))

       X(-1,-1)=
     & sig12(ihm1,jhm1)/ (4.0d0* dh1(i,j)* ds2(im1,jhm1)) + 
     & sig21(ihm1,jhm1)/ (4.0d0* dh2(i,j)* ds1(ihm1,jm1)) - 
     & gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* dh1(i,j)* dk2(ihm1,j)* 
     & ds2(im1,jhm1)) - gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)* ds1(ihm1,jm1))

       X(-1,0)=
     & (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* dk1(i,jhm1)* 
     & ds1(ihm1,j)) + gam21(ihm1,j)/ (dg2(ihm1,j)* dk1(i,j)* 
     & ds1(ihm1,j)))/ dh2(i,j) - (sig12(ihm1,jhm1)/ (4.0d0* 
     & ds2(im1,jhm1)) - sig12(ihm1,j)/ (4.0d0* ds2(im1,j)))/ 
     & dh1(i,j) + (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,j)) - 
     & sig21(ihm1,j)/ (4.0d0* ds1(ihm1,j)))/ dh2(i,j) + 
     & ((gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)) + (gam11(im1,j)* 
     & (1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j)))/ dh1(im1,j))/ 
     & ds1(ihm1,j) + gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)* 
     & ds1(i,j)))/ dh1(i,j) + sig11(ihm1,j)/ (dh1(i,j)* 
     & ds1(ihm1,j)) + (gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* 
     & ds2(im1,jhm1)) + gam12(ihm1,j)/ (dg1(ihm1,j)* ds2(im1,j)))/ 
     & (dh1(i,j)* dk2(ihm1,j))

       X(-1,+1)=
     & - sig12(ihm1,j)/ (4.0d0* dh1(i,j)* ds2(im1,j)) - 
     &   sig21(ihm1,j)/ (4.0d0* dh2(i,j)* ds1(ihm1,jp1)) - 
     & gam12(ihm1,j)/ (dg1(ihm1,j)* dh1(i,j)* dk2(ihm1,j)* 
     & ds2(im1,j)) - gam21(ihm1,j)/ (dg2(ihm1,j)* dh2(i,j)* 
     & dk1(i,j)* ds1(ihm1,jp1))

       X(0,-2)=
     & -gam22(i,jm1)/ (dh2(i,jm1)* dh2(i,j)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-1)=
     & (gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* dk2(ihm1,j)* 
     & ds2(i,jhm1)) + gam12(i,jhm1)/ (dg1(i,jhm1)* dk2(i,j)* 
     & ds2(i,jhm1)))/ dh1(i,j) + (sig12(ihm1,jhm1)/ (4.0d0* 
     & ds2(i,jhm1)) - sig12(i,jhm1)/ (4.0d0* ds2(i,jhm1)))/ 
     & dh1(i,j) - (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,jm1)) - 
     & sig21(i,jhm1)/ (4.0d0* ds1(i,jm1)))/ dh2(i,j) + 
     & ((gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)) + (gam22(i,jm1)* 
     & (1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1)))/ dh2(i,jm1))/ 
     & ds2(i,jhm1) + gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* 
     & ds2(i,j)))/ dh2(i,j) + sig22(i,jhm1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & + (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,jm1)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,jm1)))/ (dh2(i,j)* 
     & dk1(i,jhm1))

       X(0,0)=
     & - ((gam11(im1,j)/ (dh1(im1,j)* ds1(ihm1,j)) + (gam11(i,j)* 
     &   (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ dh1(i,j))/ 
     & ds1(ihm1,j) + (gam11(ip1,j)/ (dh1(ip1,j)* ds1(i,j)) + 
     & (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ 
     & dh1(i,j))/ ds1(i,j))/ dh1(i,j) - ((gam22(i,jm1)/ 
     & (dh2(i,jm1)* ds2(i,jhm1)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) 
     & + 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,jhm1) + (gam22(i,jp1)/ 
     & (dh2(i,jp1)* ds2(i,j)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 
     & 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,j))/ dh2(i,j) - 
     & (sig11(ihm1,j)/ ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ 
     & dh1(i,j) - (sig22(i,jhm1)/ ds2(i,jhm1) + sig22(i,j)/ 
     & ds2(i,j))/ dh2(i,j) - ((gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* 
     & ds2(i,jhm1)) + gam12(ihm1,j)/ (dg1(ihm1,j)* ds2(i,j)))/ 
     & dk2(ihm1,j) + (gam12(i,jhm1)/ (dg1(i,jhm1)* ds2(i,jhm1)) + 

     & gam12(i,j)/ (dg1(i,j)* ds2(i,j)))/ dk2(i,j))/ dh1(i,j) - 
     & ((gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,j)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,j)))/ dk1(i,jhm1) + 
     & (gam21(ihm1,j)/ (dg2(ihm1,j)* ds1(ihm1,j)) + gam21(i,j)/ 
     & (dg2(i,j)* ds1(i,j)))/ dk1(i,j))/ dh2(i,j) - 
     & (sig12(ihm1,jhm1)/ (4.0d0* ds2(i,jhm1)) - sig12(ihm1,j)/ 
     & (4.0d0* ds2(i,j)) - sig12(i,jhm1)/ (4.0d0* ds2(i,jhm1)) + 
     & sig12(i,j)/ (4.0d0* ds2(i,j)))/ dh1(i,j) - 
     & (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,j)) - sig21(ihm1,j)/ 
     & (4.0d0* ds1(ihm1,j)) - sig21(i,jhm1)/ (4.0d0* ds1(i,j)) + 
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j)

       X(0,+1)=
     & (gam12(ihm1,j)/ (dg1(ihm1,j)* dk2(ihm1,j)* ds2(i,j)) + 
     & gam12(i,j)/ (dg1(i,j)* dk2(i,j)* ds2(i,j)))/ dh1(i,j) - 
     & (sig12(ihm1,j)/ (4.0d0* ds2(i,j)) - sig12(i,j)/ (4.0d0* 

     & ds2(i,j)))/ dh1(i,j) + (sig21(ihm1,j)/ (4.0d0* 
     & ds1(ihm1,jp1)) - sig21(i,j)/ (4.0d0* ds1(i,jp1)))/ dh2(i,j) 
     & + ((gam22(i,j)/ (dh2(i,j)* ds2(i,j)) + (gam22(i,jp1)* 
     & (1.0d0/ ds2(i,j) + 1.0d0/ ds2(i,jhp1)))/ dh2(i,jp1))/ 
     & ds2(i,j) + gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ 
     & dh2(i,j) + sig22(i,j)/ (dh2(i,j)* ds2(i,j)) + 
     & (gam21(ihm1,j)/ (dg2(ihm1,j)* ds1(ihm1,jp1)) + gam21(i,j)/ 
     & (dg2(i,j)* ds1(i,jp1)))/ (dh2(i,j)* dk1(i,j))

       X(0,+2)=
     & -gam22(i,jp1)/ (dh2(i,j)* dh2(i,jp1)* ds2(i,j)* ds2(i,jhp1))

       X(+1,-1)=
     & - sig12(i,jhm1)/ (4.0d0* dh1(i,j)* ds2(ip1,jhm1)) - 
     &   sig21(i,jhm1)/ (4.0d0* dh2(i,j)* ds1(i,jm1)) - 
     & gam12(i,jhm1)/ (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)* 
     & ds2(ip1,jhm1)) - gam21(i,jhm1)/ (dg2(i,jhm1)* dh2(i,j)* 
     & dk1(i,jhm1)* ds1(i,jm1))

       X(+1,0)=
     & (gam21(i,jhm1)/ (dg2(i,jhm1)* dk1(i,jhm1)* ds1(i,j)) + 
     & gam21(i,j)/ (dg2(i,j)* dk1(i,j)* ds1(i,j)))/ dh2(i,j) + 
     & (sig12(i,jhm1)/ (4.0d0* ds2(ip1,jhm1)) - sig12(i,j)/ (4.0d0* 
     & ds2(ip1,j)))/ dh1(i,j) - (sig21(i,jhm1)/ (4.0d0* ds1(i,j)) - 
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) + ((gam11(i,j)/ 
     & (dh1(i,j)* ds1(i,j)) + (gam11(ip1,j)* (1.0d0/ ds1(i,j) + 
     & 1.0d0/ ds1(ihp1,j)))/ dh1(ip1,j))/ ds1(i,j) + gam11(i,j)/ 
     & (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + sig11(i,j)/ 
     & (dh1(i,j)* ds1(i,j)) + (gam12(i,jhm1)/ (dg1(i,jhm1)* 
     & ds2(ip1,jhm1)) + gam12(i,j)/ (dg1(i,j)* ds2(ip1,j)))/ 
     & (dh1(i,j)* dk2(i,j))

       X(+1,+1)=
     & sig12(i,j)/ (4.0d0* dh1(i,j)* ds2(ip1,j)) + sig21(i,j)/ 
     & (4.0d0* dh2(i,j)* ds1(i,jp1)) - gam12(i,j)/ (dg1(i,j)* 
     & dh1(i,j)* dk2(i,j)* ds2(ip1,j)) - gam21(i,j)/ (dg2(i,j)* 
     & dh2(i,j)* dk1(i,j)* ds1(i,jp1))

       X(+2,0)=
     & -gam11(ip1,j)/ (dh1(i,j)* dh1(ip1,j)* ds1(i,j)* ds1(ihp1,j))
       X(-4:4,-4:4)=-X(-4:4,-4:4)
       endif

        do i2=-4,4
          do j2=-4,4
            if(abs(X(i2,j2)).gt. 0.d0) then

            do k=1,3
                nojtmp=(k-1)*nsnq+(i-1)*ns_ibmin+j
!              if((i2.eq.0).and.(j2.eq.0)) then
                  ForceOut(nojtmp)=ForceOut(nojtmp)+X(i2,j2)*
     &            (GX(k,i+i2,j+j2)+dGXdt(k,i+i2,j+j2))
!              else
!                  ForceOut(nojtmp)=ForceOut(nojtmp)+X(i2,j2)*
!     &            GX(k,i+i2,j+j2)
!              endif
              enddo
              !X(i2,j2)=X(i2,j2)  !*alpha_coefK
            endif

          enddo
        enddo

       X(0,0)=X(0,0)+X00*scal_K

       Y=0.0

!       goto 1313

!===nonlinear parts

       if((j .le. 2).and. (BoundaryFibShellin(1,1) .eq. 0))then
          !Y(1:3,1:3,-4:4,-4:4)=0.0
              do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
              enddo
       elseif((j .ge. ns_ibmin-1)
     &       .and. (BoundaryFibShellin(1,2) .eq. 0))then
          !Y(1:3,1:3,-4:4,-4:4)=0.0
              do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
              enddo
       elseif((i .le. 2).and. (BoundaryFibShellin(2,1) .eq. 0))then
          !Y(1:3,1:3,-4:4,-4:4)=0.0
              do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
              enddo
       elseif((i .ge. nq_ibmin-1)
     &       .and. (BoundaryFibShellin(2,2) .eq. 0))then
          !Y(1:3,1:3,-4:4,-4:4)=0.0
              do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
              enddo
       elseif((j .le. 1).and. (BoundaryFibShellin(1,1) .eq. 2))then
          !Y(1:3,1:3,-4:4,-4:4)=0.0
              do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
              enddo
       elseif((j .ge. ns_ibmin)
     &       .and. (BoundaryFibShellin(1,2) .eq. 2))then
         ! Y(1:3,1:3,-4:4,-4:4)=0.0
              do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
              enddo
       elseif((i .le. 1).and. (BoundaryFibShellin(2,1) .eq. 2))then
          !Y(1:3,1:3,-4:4,-4:4)=0.0
              do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
              enddo
       elseif((i .ge. nq_ibmin)
     &       .and. (BoundaryFibShellin(2,2) .eq. 2))then
          !Y(1:3,1:3,-4:4,-4:4)=0.0
              do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
              enddo
!! ================== BNDY 1,1
       elseif(((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   .and.
     &    ((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   )then

       do i3=1,3
        P(0:2,-1:1)=GX(i3,i:ip2,j:jp2)
        do j3=1,3
          Q(0:2,-1:1)=GX(j3,i:ip2,j:jp2)
          call KnonlinearBNDY11(A0(1:8),iA0(1:8,1:2)
     &         ,P,Q
     &         ,ds1(i:ip1,j:jp1),ds2(i:ip1,j:jp1)
     &         ,dh1(i:ip1,j:jp1),dh2(i:ip1,j:jp1)
     &         ,SigB_d(3,i,1)
     &         ,SigB_d(3,i,2)
     &         ,SigB_d(1,j,2)
     &         ,SigB_d(1,j,1)
     &         ,SigB_d(3,i,3)
     &         ,SigB_d(1,j,3)
     &         )
           do i4=1,8
                Y(i3,j3,iA0(i4,1),iA0(i4,2))=-A0(i4)
           enddo
         enddo
        enddo
        do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
        enddo
!! ================== BNDY Nx,1
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   )then
      t(2) = ds2(im2,j)*3.0D0
      t(3) = ds2(i,j)*3.0D0
      t(4) = ds1(im1,j)*(-3.0D0)
      t(5) = ds1(im2,j)+t(4)
      t(6) = 1.0D0/t(5)
      t(34) = ds2(im2,jp1)*(-1.0D0)
      t(7) = -t(34)
      t(8) = t(2)-t(7)
      t(9) = 1.0D0/t(8)
      t(13) = ds2(i,jp1)*(-1.0D0)
      t(10) = -t(13)
      t(11) = t(3)-t(10)
      t(12) = 1.0D0/t(11)
      t(14) = ds2(im1,j)*3.0D0
      t(15) = ds2(i,j)*1.5D0
      t(26) = ds2(i,jp1)*(-5.0D-1)
      t(16) = -t(26)
      t(17) = t(15)-t(16)
      t(18) = 1.0D0/t(17)
      t(19) = ds1(im1,jp1)*(-3.0D0)
      t(20) = ds1(im2,jp1)+t(19)
      t(21) = 1.0D0/t(20)
      t(22) = ds2(im1,j)*1.5D0
      t(44) = ds2(im1,jp1)*(-5.0D-1)
      t(23) = -t(44)
      t(24) = t(22)-t(23)
      t(25) = 1.0D0/t(24)
      t(27) = ds1(im1,jp2)*(-3.0D0)
      t(28) = ds1(im2,jp2)+t(27)
      t(29) = 1.0D0/t(28)
      t(30) = t(15)+t(26)
      t(31) = 1.0D0/t(30)
      t(32) = t(3)+t(13)
      t(33) = 1.0D0/t(32)
      t(35) = ds2(im2,j)*1.5D0
      t(55) = ds2(im2,jp1)*(-5.0D-1)
      t(36) = -t(55)
      t(37) = t(35)-t(36)
      t(38) = 1.0D0/t(37)
      t(39) = t(6)*t(38)*(-1.0D0)
      t(40) = t(2)+t(34)
      t(41) = 1.0D0/t(40)
      t(42) = ds2(im1,jp1)*(-1.0D0)
      t(43) = t(3)+t(13)+t(14)+t(42)
      t(45) = t(22)+t(44)
      t(46) = 1.0D0/t(45)
      t(47) = t(29)*t(31)*(-1.0D0)
      t(48) = t(6)*t(31)*4.5D0
      t(49) = ds1(im2,j)+ds1(im2,jp1)+t(4)+t(19)
      t(50) = t(6)*t(21)*t(33)*t(49)*(-6.0D0)
      t(51) = ds1(im2,j)+ds1(im2,jp2)+t(4)+t(27)
      t(52) = t(6)*t(29)*t(33)*t(51)*1.5D0
      t(53) = t(2)+t(3)+t(13)+t(34)
      t(54) = t(6)*t(33)*t(41)*t(53)*1.5D0
      t(56) = t(21)*t(31)*(-1.0D0)
      t(57) = t(6)*t(41)*(1.0D0/2.0D0)
      t(58) = t(29)*t(33)*(1.0D0/2.0D0)
      t(59) = t(57)+t(58)
      t(60) = t(14)+t(42)
      t(61) = 1.0D0/t(60)
      t(62) = t(6)*t(33)*t(43)*t(61)*(-6.0D0)
      t(63) = t(6)*t(46)*4.0D0
      t(64) = t(21)*t(31)*4.0D0
      t(65) = t(63)+t(64)
      t(66) = t(6)*t(46)*(-1.0D0)
      t(67) = t(47)+t(66)
      Y(1,1,-2,0) = t(6)*t(9)*t(12)*(t(2)+t(3)-t(7)-t(10))*1.5D0
      Y(1,1,-2,1) = t(39)-t(18)*t(21)*1.0D0
      Y(1,1,-2,2) = t(6)*t(9)*(1.0D0/2.0D0)+t(12)*t(29)*(1.0D0/2.0
     &D0)
      Y(1,1,-1,0) = (t(6)*t(33)*t(43)*6.0D0)/(ds2(im1,jp1)*1.0D0-t
     &(14))
      Y(1,1,-1,1) = t(6)*t(25)*4.0D0+t(18)*t(21)*4.0D0
      Y(1,1,-1,2) = t(47)-t(6)*t(25)*1.0D0
      Y(1,1,0,0) = t(48)
      Y(1,1,0,1) = t(50)
      Y(1,1,0,2) = t(52)
      Y(2,2,-2,0) = t(54)
      Y(2,2,-2,1) = t(39)+t(56)
      Y(2,2,-2,2) = t(59)
      Y(2,2,-1,0) = t(62)
      Y(2,2,-1,1) = t(65)
      Y(2,2,-1,2) = t(67)
      Y(2,2,0,0) = t(48)
      Y(2,2,0,1) = t(50)
      Y(2,2,0,2) = t(52)
      Y(3,3,-2,0) = t(54)
      Y(3,3,-2,1) = t(56)-(t(6)*1.0D0)/(t(35)+t(55))
      Y(3,3,-2,2) = t(59)
      Y(3,3,-1,0) = t(62)
      Y(3,3,-1,1) = t(65)
      Y(3,3,-1,2) = t(67)
      Y(3,3,0,0) = t(48)
      Y(3,3,0,1) = t(50)
      Y(3,3,0,2) = t(52)
              !Y(1:3,1:3,-4:4,-4:4)=alpha_coefK*Y(1:3,1:3,-4:4,-4:4)

!! ================== BNDY 1,Ny
       elseif(((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
      t(2) = ds1(i,jm2)*3.0D0
      t(3) = ds1(i,j)*3.0D0
      t(11) = ds1(ip1,j)*(-1.0D0)
      t(4) = -t(11)
      t(5) = t(3)-t(4)
      t(6) = 1.0D0/t(5)
      t(7) = ds2(i,jm1)*(-3.0D0)
      t(8) = ds2(i,jm2)+t(7)
      t(9) = 1.0D0/t(8)
      t(10) = ds1(i,jm1)*3.0D0
      t(12) = ds1(i,j)*1.5D0
      t(16) = ds1(ip1,j)*(-5.0D-1)
      t(13) = -t(16)
      t(14) = t(12)-t(13)
      t(15) = 1.0D0/t(14)
      t(17) = ds2(ip1,jm1)*(-3.0D0)
      t(18) = ds2(ip1,jm2)+t(17)

      t(19) = 1.0D0/t(18)
      t(34) = ds1(ip1,jm2)*(-1.0D0)
      t(20) = -t(34)
      t(21) = t(2)-t(20)
      t(22) = 1.0D0/t(21)
      t(23) = t(3)+t(11)
      t(24) = 1.0D0/t(23)
      t(25) = ds1(i,jm1)*1.5D0
      t(42) = ds1(ip1,jm1)*(-5.0D-1)
      t(26) = -t(42)
      t(27) = t(25)-t(26)
      t(28) = 1.0D0/t(27)
      t(29) = t(12)+t(16)
      t(30) = 1.0D0/t(29)
      t(31) = ds2(ip2,jm1)*(-3.0D0)
      t(32) = ds2(ip2,jm2)+t(31)
      t(33) = 1.0D0/t(32)
      t(35) = ds1(ip1,jm1)*(-1.0D0)
      t(36) = t(3)+t(10)+t(11)+t(35)
      t(37) = ds1(i,jm2)*1.5D0
      t(60) = ds1(ip1,jm2)*(-5.0D-1)
      t(38) = -t(60)
      t(39) = t(37)-t(38)
      t(40) = 1.0D0/t(39)
      t(41) = t(9)*t(40)*(-1.0D0)
      t(43) = t(19)*t(30)*4.0D0
      t(44) = ds2(i,jm2)+ds2(ip1,jm2)+t(7)+t(17)
      t(45) = t(9)*t(19)*t(24)*t(44)*(-6.0D0)
      t(46) = t(2)+t(34)
      t(47) = 1.0D0/t(46)
      t(48) = t(24)*t(33)*(1.0D0/2.0D0)
      t(49) = t(25)+t(42)
      t(50) = 1.0D0/t(49)
      t(51) = t(30)*t(33)*(-1.0D0)
      t(52) = ds2(i,jm2)+ds2(ip2,jm2)+t(7)+t(31)
      t(53) = t(9)*t(24)*t(33)*t(52)*1.5D0
      t(54) = t(2)+t(3)+t(11)+t(34)
      t(55) = t(9)*t(24)*t(47)*t(54)*1.5D0
      t(56) = t(10)+t(35)
      t(57) = 1.0D0/t(56)
      t(58) = t(9)*t(24)*t(36)*t(57)*(-6.0D0)
      t(59) = t(9)*t(30)*4.5D0
      t(61) = t(19)*t(30)*(-1.0D0)
      t(62) = t(9)*t(50)*4.0D0
      t(63) = t(43)+t(62)
      t(64) = t(9)*t(47)*(1.0D0/2.0D0)
      t(65) = t(48)+t(64)
      t(66) = t(9)*t(50)*(-1.0D0)
      t(67) = t(51)+t(66)
      Y(1,1,0,-2) = t(6)*t(9)*t(22)*(t(2)+t(3)-t(4)-t(20))*1.5D0
      Y(1,1,0,-1) = (t(6)*t(9)*t(36)*6.0D0)/(ds1(ip1,jm1)*1.0D0-t(
     &10))
      Y(1,1,0,0) = t(9)*t(15)*4.5D0
      Y(1,1,1,-2) = t(41)-t(15)*t(19)*1.0D0
      Y(1,1,1,-1) = t(43)+t(9)*t(28)*4.0D0
      Y(1,1,1,0) = t(45)
      Y(1,1,2,-2) = t(48)+t(9)*t(22)*(1.0D0/2.0D0)
      Y(1,1,2,-1) = t(51)-t(9)*t(28)*1.0D0
      Y(1,1,2,0) = t(53)
      Y(2,2,0,-2) = t(55)
      Y(2,2,0,-1) = t(58)
      Y(2,2,0,0) = t(59)
      Y(2,2,1,-2) = t(41)+t(61)
      Y(2,2,1,-1) = t(63)
      Y(2,2,1,0) = t(45)
      Y(2,2,2,-2) = t(65)
      Y(2,2,2,-1) = t(67)
      Y(2,2,2,0) = t(53)
      Y(3,3,0,-2) = t(55)
      Y(3,3,0,-1) = t(58)
      Y(3,3,0,0) = t(59)
      Y(3,3,1,-2) = t(61)-(t(9)*1.0D0)/(t(37)+t(60))
      Y(3,3,1,-1) = t(63)
      Y(3,3,1,0) = t(45)
      Y(3,3,2,-2) = t(65)
      Y(3,3,2,-1) = t(67)
      Y(3,3,2,0) = t(53)
              !Y(1:3,1:3,-4:4,-4:4)=alpha_coefK*Y(1:3,1:3,-4:4,-4:4)

!! ================== BNDY Nx,Ny
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
      t(2) = ds1(im1,j)*(-3.0D0)
      t(3) = ds1(im2,j)+t(2)
      t(4) = 1.0D0/t(3)
      t(5) = ds2(im2,jm1)*(-3.0D0)
      t(6) = ds2(im2,jm2)+t(5)
      t(7) = 1.0D0/t(6)
      t(8) = ds2(i,jm1)*(-3.0D0)
      t(9) = ds2(i,jm2)+t(8)
      t(10) = 1.0D0/t(9)
      t(11) = ds1(im1,jm2)*(-3.0D0)
      t(12) = ds1(im2,jm2)+t(11)
      t(13) = 1.0D0/t(12)
      t(14) = ds2(im1,jm1)*(-3.0D0)
      t(15) = ds2(im1,jm2)+t(14)
      t(16) = 1.0D0/t(15)
      t(17) = ds1(im1,jm1)*(-3.0D0)
      t(18) = ds1(im2,jm1)+t(17)
      t(19) = 1.0D0/t(18)
      t(20) = t(4)*t(7)*(1.0D0/2.0D0)
      t(21) = t(10)*t(13)*(1.0D0/2.0D0)
      t(22) = t(20)+t(21)
      t(23) = t(4)*t(7)*(-2.0D0)
      t(24) = t(10)*t(19)*(-2.0D0)
      t(25) = t(23)+t(24)
      t(26) = ds2(im2,jm2)+ds2(i,jm2)+t(5)+t(8)
      t(27) = t(4)*t(7)*t(10)*t(26)*1.5D0
      t(28) = t(4)*t(16)*(-2.0D0)
      t(29) = t(10)*t(13)*(-2.0D0)
      t(30) = t(28)+t(29)
      t(31) = t(4)*t(16)*8.0D0
      t(32) = t(10)*t(19)*8.0D0
      t(33) = t(31)+t(32)
      t(34) = ds2(im1,jm2)+ds2(i,jm2)+t(8)+t(14)
      t(35) = t(4)*t(10)*t(16)*t(34)*(-6.0D0)
      t(36) = ds1(im2,jm2)+ds1(im2,j)+t(2)+t(11)
      t(37) = t(4)*t(10)*t(13)*t(36)*1.5D0
      t(38) = ds1(im2,jm1)+ds1(im2,j)+t(2)+t(17)
      t(39) = t(4)*t(10)*t(19)*t(38)*(-6.0D0)
      t(40) = t(4)*t(10)*9.0D0
      Y(1,1,-2,-2) = t(22)
      Y(1,1,-2,-1) = t(25)
      Y(1,1,-2,0) = t(27)
      Y(1,1,-1,-2) = t(30)
      Y(1,1,-1,-1) = t(33)
      Y(1,1,-1,0) = t(35)
      Y(1,1,0,-2) = t(37)
      Y(1,1,0,-1) = t(39)
      Y(1,1,0,0) = t(40)
      Y(2,2,-2,-2) = t(22)
      Y(2,2,-2,-1) = t(25)
      Y(2,2,-2,0) = t(27)
      Y(2,2,-1,-2) = t(30)
      Y(2,2,-1,-1) = t(33)
      Y(2,2,-1,0) = t(35)
      Y(2,2,0,-2) = t(37)
      Y(2,2,0,-1) = t(39)
      Y(2,2,0,0) = t(40)
      Y(3,3,-2,-2) = t(22)
      Y(3,3,-2,-1) = t(25)
      Y(3,3,-2,0) = t(27)
      Y(3,3,-1,-2) = t(30)
      Y(3,3,-1,-1) = t(33)
      Y(3,3,-1,0) = t(35)
      Y(3,3,0,-2) = t(37)
      Y(3,3,0,-1) = t(39)
      Y(3,3,0,0) = t(40)
             ! Y(1:3,1:3,-4:4,-4:4)=alpha_coefK*Y(1:3,1:3,-4:4,-4:4)
!! ================== BNDY 1,j  2<=j<=Ny-1
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .gt. jst) .and. (j.lt. jen)) then

       do i3=1,3
        P(0:2,-1:1)=GX(i3,i:ip2,jm1:jp1)
        do j3=1,3
          Q(0:2,-1:1)=GX(j3,i:ip2,jm1:jp1)
          call KnonlinearBNDY1j(A0(1:9),iA0(1:9,1:2)
     &         ,P,Q
     &         ,ds1(i:ip1,jm1:jp1),ds2(i:ip1,jm1:jp1)
     &         ,dh1(i:ip1,jm1:jp1),dh2(i:ip1,jm1:jp1)
     &         ,SigB_d(1,jm1:jp1,3)
     &         ,Siga_d(i,jm1:jp1,5)
     &         ,SigB_d(1,jm1:jp1,2),SigB_d(1,jm1:jp1,1)
     &         ,Siga_d(i,j,1),Siga_d(i,j,2) )
           do i4=1,9
                Y(i3,j3,iA0(i4,1),iA0(i4,2))=-A0(i4)
           enddo
         enddo
        enddo
        do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
        enddo
!! ================== BNDY i,1  2<=i<=Nx-1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .gt. ist) .and. (i.lt. ien)) then

       do i3=1,3
        P(0:2,-1:1)=GX(i3,im1:ip1,j:jp2)
        do j3=1,3
          Q(0:2,-1:1)=GX(j3,im1:ip1,j:jp2)
          call KnonlinearBNDY1i(A0(1:9),iA0(1:9,1:2)
     &         ,P,Q
     &         ,ds1(im1:ip1,j:jp1),ds2(im1:ip1,j:jp1)
     &         ,dh1(im1:ip1,j:jp1),dh2(im1:ip1,j:jp1)
     &         ,SigB_d(3,im1:ip1,3)
     &         ,Siga_d(im1:ip1,j,5)
     &         ,SigB_d(3,im1:ip1,1),SigB_d(3,im1:ip1,2)
     &         ,Siga_d(i,j,2),Siga_d(i,j,1) )
           do i4=1,9
                Y(i3,j3,iA0(i4,1),iA0(i4,2))=-A0(i4)
           enddo
         enddo
        enddo
        do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)

        enddo

!! ================== BNDY N,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin)   
     &        .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .ge. jst) .and. (j.le. jen)) then

      if(formulationFlagb .eq. 1) then
      t(2) = ds1(im1,j)*(-3.0D0)
      t(3) = ds1(im2,j)+t(2)
      t(4) = 1.0D0/t(3)
      t(5) = 1.0D0/t(3)**3
      t(15) = GX(1,im1,j)*(-4.0D0)
      t(16) = GX(1,i,j)*3.0D0
      t(6) = t(15)+t(16)+GX(1,im2,j)
      t(20) = ds2(i,jm1)*GX(1,i,j)
      t(21) = ds2(i,j)*GX(1,i,jm1)
      t(22) = ds2(i,jm1)*GX(1,i,jp1)*(-1.0D0)
      t(23) = ds2(i,j)*GX(1,i,j)*(-1.0D0)
      t(7) = t(20)+t(21)+t(22)+t(23)
      t(8) = 1.0D0/dh2(im2,j)
      t(9) = 1.0D0/ds2(im2,jm1)
      t(10) = 1.0D0/ds2(im2,j)
      t(11) = 1.0D0/t(3)**5
      t(12) = ds1(im3,j)**2
      t(13) = ds1(im2,j)**2
      t(14) = ds1(im1,j)**2
      t(17) = t(6)**2
      t(18) = 1.0D0/ds2(i,jm1)**2
      t(19) = 1.0D0/ds2(i,j)**2
      t(24) = t(7)**2
      t(25) = 1.0D0/dh2(im1,j)
      t(26) = 1.0D0/ds2(im1,jm1)
      t(27) = 1.0D0/ds2(im1,j)
      t(28) = 1.0D0/ds2(i,jm1)
      t(29) = ds1(im2,j)*5.0D-1
      t(41) = ds1(im1,j)*(-1.5D0)
      t(30) = -t(41)
      t(31) = t(29)-t(30)
      t(32) = 1.0D0/t(31)
      t(33) = GX(1,i,j)*(-1.0D0)
      t(34) = t(33)+GX(1,i,jm1)
      t(35) = t(28)*t(34)*5.0D-1
      t(36) = 1.0D0/ds2(i,j)
      t(37) = GX(1,i,jp1)*(-1.0D0)
      t(38) = t(37)+GX(1,i,j)
      t(39) = t(36)*t(38)*5.0D-1
      t(40) = t(35)+t(39)
      t(42) = GX(1,im2,j)*5.0D-1
      t(43) = GX(1,im1,j)*(-2.0D0)
      t(44) = GX(1,i,j)*1.5D0
      t(45) = t(42)+t(43)+t(44)
      t(46) = 1.0D0/dh2(i,j)
      t(47) = ds2(i,j)*(-1.0D0)
      t(48) = ds2(i,jm1)+t(47)
      t(49) = GX(1,im2,j)*(1.0D0/2.0D0)
      t(50) = t(43)+t(44)+t(49)
      t(51) = t(29)+t(41)
      t(52) = 1.0D0/t(51)
      t(53) = 1.0D0/t(51)**3
      t(54) = GX(2,im2,j)*5.0D-1
      t(55) = GX(2,im1,j)*(-2.0D0)
      t(56) = GX(2,i,j)*1.5D0
      t(57) = t(54)+t(55)+t(56)
      t(58) = GX(2,i,j)*(-1.0D0)
      t(59) = t(58)+GX(2,i,jm1)
      t(60) = t(28)*t(59)*5.0D-1
      t(61) = GX(2,i,jp1)*(-1.0D0)
      t(62) = t(61)+GX(2,i,j)
      t(63) = t(36)*t(62)*5.0D-1
      t(64) = t(60)+t(63)
      t(65) = t(28)*5.0D-1
      t(79) = t(36)*(-5.0D-1)
      t(66) = -t(79)
      t(67) = t(65)-t(66)
      t(68) = GX(3,im2,j)*5.0D-1
      t(69) = GX(3,im1,j)*(-2.0D0)
      t(70) = GX(3,i,j)*1.5D0
      t(71) = t(68)+t(69)+t(70)
      t(72) = GX(3,i,j)*(-1.0D0)
      t(73) = t(72)+GX(3,i,jm1)
      t(74) = t(28)*t(73)*5.0D-1
      t(75) = GX(3,i,jp1)*(-1.0D0)
      t(76) = t(75)+GX(3,i,j)
      t(77) = t(36)*t(76)*5.0D-1
      t(78) = t(74)+t(77)
      t(80) = t(65)+t(79)
      t(81) = SigB_d(2,j,1)*t(45)*t(53)*t(57)*(-5.0D-1)
      t(82) = SigB_d(2,j,3)*t(40)*t(52)*t(64)*(-5.0D-1)
      t(83) = t(81)+t(82)
      t(84) = SigB_d(2,j,1)*t(45)*t(53)*t(57)*2.0D0
      t(85) = SigB_d(2,j,3)*t(40)*t(52)*t(64)*2.0D0
      t(86) = t(84)+t(85)
      t(87) = SigB_d(2,j,1)*t(45)*t(53)*t(57)*(-1.5D0)
      t(88) = SigB_d(2,j,3)*t(40)*t(52)*t(64)*(-1.5D0)
      t(89) = gam11(i,j)*t(5)*1.2D1
      t(90) = 1.0D0/t(3)**4
      t(91) = ds1(im3,j)*(-6.0D0)
      t(92) = ds1(im2,j)*1.1D1
      t(93) = ds1(im1,j)*9.0D0
      t(94) = t(91)+t(92)+t(93)
      t(95) = gam11(i,j)*t(90)*t(94)*8.0D0
      t(96) = e5c*gam12(i,j)*t(4)*t(8)*t(9)
      t(97) = sigB(2,j,1)*t(4)*(-1.0D0)
      t(98) = ds1(im4,j)*ds1(im2,j)*3.0D0
      t(99) = ds1(im4,j)*ds1(im1,j)*(-9.0D0)
      t(100) = ds1(im3,j)*ds1(im2,j)*1.3D1
      t(101) = ds1(im3,j)*ds1(im1,j)*1.29D2
      t(102) = ds1(im2,j)*ds1(im1,j)*(-2.84D2)
      t(103) = t(12)*(-1.2D1)
      t(104) = t(13)*2.5D1
      t(105) = t(14)*3.9D1
      t(106) = t(98)+t(99)+t(100)+t(101)+t(102)+t(103)+t(104)+t(10
     &5)
      t(107) = gam11(i,j)*t(11)*t(106)*(-4.0D0)
      t(116) = GX(2,im1,j)*(-4.0D0)
      t(117) = GX(2,i,j)*3.0D0
      t(108) = t(116)+t(117)+GX(2,im2,j)
      t(119) = ds2(i,jm1)*GX(2,i,j)
      t(120) = ds2(i,j)*GX(2,i,jm1)
      t(121) = ds2(i,jm1)*GX(2,i,jp1)*(-1.0D0)
      t(122) = ds2(i,j)*GX(2,i,j)*(-1.0D0)
      t(109) = t(119)+t(120)+t(121)+t(122)
      t(110) = ds2(im2,jm1)+ds2(im2,j)
      t(111) = e5c*gam12(i,j)*t(4)*t(8)*t(9)*t(10)*t(110)*(-1.0D0)
      t(112) = e5c*gam12(i,j)*t(4)*t(8)*t(10)
      t(113) = e5c*gam12(i,j)*t(4)*t(25)*t(26)*(-4.0D0)
      t(114) = sigB(2,j,1)*t(4)*4.0D0
      t(163) = ds1(im4,j)*ds1(im2,j)*6.0D0
      t(164) = ds1(im4,j)*ds1(im1,j)*(-1.8D1)
      t(165) = ds1(im3,j)*ds1(im2,j)*9.2D1
      t(166) = ds1(im3,j)*ds1(im1,j)*6.0D1
      t(167) = ds1(im2,j)*ds1(im1,j)*(-7.6D1)
      t(168) = t(12)*(-2.4D1)
      t(169) = t(13)*(-1.09D2)
      t(170) = t(14)*3.3D1
      t(171) = t(163)+t(164)+t(165)+t(166)+t(167)+t(168)+t(169)+t(
     &170)
      t(115) = gam11(i,j)*t(11)*t(171)*8.0D0
      t(118) = t(108)**2
      t(123) = t(109)**2
      t(124) = ds2(im1,jm1)+ds2(im1,j)
      t(125) = e5c*gam12(i,j)*t(4)*t(25)*t(26)*t(27)*t(124)*4.0D0
      t(126) = e5c*gam12(i,j)*t(4)*t(25)*t(27)*(-4.0D0)
      t(127) = sigB(2,j,3)*t(28)*(-5.0D-1)
      t(128) = sigB(2,j,1)*t(4)*(-3.0D0)
      t(129) = ds1(im4,j)*ds1(im2,j)*9.0D0
      t(130) = ds1(im4,j)*ds1(im1,j)*(-2.7D1)
      t(131) = ds1(im3,j)*ds1(im2,j)*1.59D2
      t(132) = ds1(im3,j)*ds1(im1,j)*2.7D1
      t(133) = ds1(im2,j)*ds1(im1,j)*6.6D1
      t(134) = t(12)*(-3.6D1)
      t(135) = t(13)*(-2.18D2)
      t(136) = t(129)+t(130)+t(131)+t(132)+t(133)+t(134)+t(135)
      t(137) = gam11(i,j)*t(11)*t(136)*(-4.0D0)
      t(138) = sigB(2,j,3)*t(28)*t(36)*t(48)*(-5.0D-1)
      t(139) = ds2(i,jm1)+ds2(i,j)
      t(140) = e5c*gam12(i,j)*t(4)*t(28)*t(36)*t(46)*t(139)*(-3.0D
     &0)
      t(141) = GX(2,im2,j)*(1.0D0/2.0D0)
      t(142) = t(55)+t(56)+t(141)
      t(143) = sigB(2,j,3)*t(36)*5.0D-1
      t(144) = e5c*gam12(i,j)*t(36)*t(46)*t(52)*1.5D0
      t(145) = SigB_d(2,j,1)*t(45)*t(53)*t(71)*(-5.0D-1)
      t(146) = SigB_d(2,j,3)*t(40)*t(52)*t(78)*(-5.0D-1)
      t(147) = t(145)+t(146)
      t(148) = SigB_d(2,j,1)*t(45)*t(53)*t(71)*2.0D0
      t(149) = SigB_d(2,j,3)*t(40)*t(52)*t(78)*2.0D0
      t(150) = t(148)+t(149)
      t(151) = SigB_d(2,j,1)*t(45)*t(53)*t(71)*(-1.5D0)
      t(152) = SigB_d(2,j,3)*t(40)*t(52)*t(78)*(-1.5D0)
      t(153) = SigB_d(2,j,1)*t(53)*t(57)*t(71)*(-5.0D-1)
      t(154) = SigB_d(2,j,3)*t(52)*t(64)*t(78)*(-5.0D-1)
      t(155) = t(153)+t(154)
      t(156) = SigB_d(2,j,1)*t(53)*t(57)*t(71)*2.0D0
      t(157) = SigB_d(2,j,3)*t(52)*t(64)*t(78)*2.0D0
      t(158) = t(156)+t(157)
      t(159) = SigB_d(2,j,1)*t(53)*t(57)*t(71)*(-1.5D0)
      t(160) = SigB_d(2,j,3)*t(52)*t(64)*t(78)*(-1.5D0)
      t(172) = GX(3,im1,j)*(-4.0D0)
      t(173) = GX(3,i,j)*3.0D0
      t(161) = t(172)+t(173)+GX(3,im2,j)
      t(175) = ds2(i,jm1)*GX(3,i,j)
      t(176) = ds2(i,j)*GX(3,i,jm1)
      t(177) = ds2(i,jm1)*GX(3,i,jp1)*(-1.0D0)
      t(178) = ds2(i,j)*GX(3,i,j)*(-1.0D0)
      t(162) = t(175)+t(176)+t(177)+t(178)
      t(174) = t(161)**2
      t(179) = t(162)**2
      t(180) = e5c*gam12(i,j)*t(28)*t(46)*t(52)*1.5D0
      t(181) = GX(3,im2,j)*(1.0D0/2.0D0)
      t(182) = t(69)+t(70)+t(181)
      Y(1,1,-4,0) = t(89)
      Y(1,1,-3,0) = t(95)
      Y(1,1,-2,-1) = t(96)
      Y(1,1,-2,0) = t(97)+t(107)+t(111)-SigB_d(2,j,1)*t(5)*t(17)*1
     &.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(24)*2.5D-1
      Y(1,1,-2,1) = t(112)
      Y(1,1,-1,-1) = t(113)
      Y(1,1,-1,0) = t(114)+t(115)+t(125)+SigB_d(2,j,1)*t(5)*t(17)*
     &4.0D0+SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(24)*1.0D0
      Y(1,1,-1,1) = t(126)
      Y(1,1,0,-1) = t(127)+e5c*gam12(i,j)*t(28)*t(32)*t(46)*1.5D0-
     &SigB_d(2,j,3)*t(28)*t(32)*t(40)*t(45)*5.0D-1-SigB_d(2,j,2)*t(28)*t
     &(40)*t(45)*t(52)*5.0D-1
      Y(1,1,0,0) = t(128)+t(137)+t(138)+t(140)-SigB_d(2,j,1)*t(5)*
     &t(17)*3.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(24)*7.5D-1-SigB_d(2,j
     &,3)*t(4)*t(7)*t(18)*t(19)*t(48)*t(50)*5.0D-1-SigB_d(2,j,2)*t(4)*t(
     &7)*t(18)*t(19)*t(48)*t(50)*5.0D-1
      Y(1,1,0,1) = t(143)+t(144)+SigB_d(2,j,3)*t(36)*t(40)*t(45)*t
     &(52)*5.0D-1+SigB_d(2,j,2)*t(36)*t(40)*t(45)*t(52)*5.0D-1
      Y(1,2,-2,0) = t(83)
      Y(1,2,-1,0) = t(86)
      Y(1,2,0,-1) = SigB_d(2,j,3)*t(28)*t(40)*t(52)*t(57)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(45)*t(52)*t(64)*5.0D-1
      Y(1,2,0,0) = t(87)+t(88)+SigB_d(2,j,3)*t(40)*t(52)*t(57)*t(6
     &7)+SigB_d(2,j,2)*t(45)*t(52)*t(64)*t(67)
      Y(1,2,0,1) = SigB_d(2,j,3)*t(36)*t(40)*t(52)*t(57)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(45)*t(52)*t(64)*5.0D-1
      Y(1,3,-2,0) = t(147)
      Y(1,3,-1,0) = t(150)
      Y(1,3,0,-1) = SigB_d(2,j,3)*t(28)*t(40)*t(52)*t(71)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(45)*t(52)*t(78)*5.0D-1
      Y(1,3,0,0) = t(151)+t(152)+SigB_d(2,j,3)*t(40)*t(52)*t(71)*t
     &(80)+SigB_d(2,j,2)*t(45)*t(52)*t(78)*t(80)
      Y(1,3,0,1) = SigB_d(2,j,3)*t(36)*t(40)*t(52)*t(71)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(45)*t(52)*t(78)*5.0D-1
      Y(2,1,-2,0) = t(83)
      Y(2,1,-1,0) = t(86)
      Y(2,1,0,-1) = SigB_d(2,j,3)*t(28)*t(45)*t(52)*t(64)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(40)*t(52)*t(57)*5.0D-1
      Y(2,1,0,0) = t(87)+t(88)+SigB_d(2,j,3)*t(45)*t(52)*t(64)*t(8
     &0)+SigB_d(2,j,2)*t(40)*t(52)*t(57)*t(80)
      Y(2,1,0,1) = SigB_d(2,j,3)*t(36)*t(45)*t(52)*t(64)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(40)*t(52)*t(57)*5.0D-1
      Y(2,2,-4,0) = t(89)
      Y(2,2,-3,0) = t(95)
      Y(2,2,-2,-1) = t(96)
      Y(2,2,-2,0) = t(97)+t(107)+t(111)-SigB_d(2,j,1)*t(5)*t(118)*
     &1.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(123)*2.5D-1
      Y(2,2,-2,1) = t(112)
      Y(2,2,-1,-1) = t(113)
      Y(2,2,-1,0) = t(114)+t(115)+t(125)+SigB_d(2,j,1)*t(5)*t(118)
     &*4.0D0+SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(123)*1.0D0
      Y(2,2,-1,1) = t(126)
      Y(2,2,0,-1) = t(127)+t(180)-SigB_d(2,j,3)*t(28)*t(52)*t(57)*
     &t(64)*5.0D-1-SigB_d(2,j,2)*t(28)*t(52)*t(57)*t(64)*5.0D-1
      Y(2,2,0,0) = t(128)+t(137)+t(138)+t(140)-SigB_d(2,j,1)*t(5)*
     &t(118)*3.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(123)*7.5D-1-SigB_d(2
     &,j,3)*t(4)*t(18)*t(19)*t(48)*t(109)*t(142)*5.0D-1-SigB_d(2,j,2)*t(
     &4)*t(18)*t(19)*t(48)*t(109)*t(142)*5.0D-1
      Y(2,2,0,1) = t(143)+t(144)+SigB_d(2,j,3)*t(36)*t(52)*t(57)*t
     &(64)*5.0D-1+SigB_d(2,j,2)*t(36)*t(52)*t(57)*t(64)*5.0D-1
      Y(2,3,-2,0) = t(155)
      Y(2,3,-1,0) = t(158)
      Y(2,3,0,-1) = SigB_d(2,j,3)*t(28)*t(52)*t(64)*t(71)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(52)*t(57)*t(78)*5.0D-1
      Y(2,3,0,0) = t(159)+t(160)+SigB_d(2,j,3)*t(52)*t(64)*t(71)*t
     &(80)+SigB_d(2,j,2)*t(52)*t(57)*t(78)*t(80)
      Y(2,3,0,1) = SigB_d(2,j,3)*t(36)*t(52)*t(64)*t(71)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(52)*t(57)*t(78)*5.0D-1
      Y(3,1,-2,0) = t(147)
      Y(3,1,-1,0) = t(150)
      Y(3,1,0,-1) = SigB_d(2,j,3)*t(28)*t(45)*t(52)*t(78)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(40)*t(52)*t(71)*5.0D-1
      Y(3,1,0,0) = t(151)+t(152)+SigB_d(2,j,3)*t(45)*t(52)*t(78)*t
     &(80)+SigB_d(2,j,2)*t(40)*t(52)*t(71)*t(80)
      Y(3,1,0,1) = SigB_d(2,j,3)*t(36)*t(45)*t(52)*t(78)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(40)*t(52)*t(71)*5.0D-1
      Y(3,2,-2,0) = t(155)
      Y(3,2,-1,0) = t(158)
      Y(3,2,0,-1) = SigB_d(2,j,3)*t(28)*t(52)*t(57)*t(78)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(52)*t(64)*t(71)*5.0D-1
      Y(3,2,0,0) = t(159)+t(160)+SigB_d(2,j,3)*t(52)*t(57)*t(78)*t
     &(80)+SigB_d(2,j,2)*t(52)*t(64)*t(71)*t(80)
      Y(3,2,0,1) = SigB_d(2,j,3)*t(36)*t(52)*t(57)*t(78)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(52)*t(64)*t(71)*5.0D-1
      Y(3,3,-4,0) = t(89)
      Y(3,3,-3,0) = t(95)
      Y(3,3,-2,-1) = t(96)
      Y(3,3,-2,0) = t(97)+t(107)+t(111)-SigB_d(2,j,1)*t(5)*t(174)*
     &1.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(179)*2.5D-1
      Y(3,3,-2,1) = t(112)
      Y(3,3,-1,-1) = t(113)
      Y(3,3,-1,0) = t(114)+t(115)+t(125)+SigB_d(2,j,1)*t(5)*t(174)
     &*4.0D0+SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(179)*1.0D0
      Y(3,3,-1,1) = t(126)
      Y(3,3,0,-1) = t(127)+t(180)-SigB_d(2,j,3)*t(28)*t(52)*t(71)*
     &t(78)*5.0D-1-SigB_d(2,j,2)*t(28)*t(52)*t(71)*t(78)*5.0D-1
      Y(3,3,0,0) = t(128)+t(137)+t(138)+t(140)-SigB_d(2,j,1)*t(5)*
     &t(174)*3.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(179)*7.5D-1-SigB_d(2
     &,j,3)*t(4)*t(18)*t(19)*t(48)*t(162)*t(182)*5.0D-1-SigB_d(2,j,2)*t(
     &4)*t(18)*t(19)*t(48)*t(162)*t(182)*5.0D-1
      Y(3,3,0,1) = t(143)+t(144)+SigB_d(2,j,3)*t(36)*t(52)*t(71)*t
     &(78)*5.0D-1+SigB_d(2,j,2)*t(36)*t(52)*t(71)*t(78)*5.0D-1

      else
      t(2) = ds1(im1,j)*(-3.0D0)
      t(3) = ds1(im2,j)+t(2)
      t(4) = 1.0D0/t(3)
      t(5) = 1.0D0/t(3)**3
      t(15) = GX(1,im1,j)*(-4.0D0)
      t(16) = GX(1,i,j)*3.0D0
      t(6) = t(15)+t(16)+GX(1,im2,j)
      t(20) = ds2(i,jm1)*GX(1,i,j)
      t(21) = ds2(i,j)*GX(1,i,jm1)
      t(22) = ds2(i,jm1)*GX(1,i,jp1)*(-1.0D0)
      t(23) = ds2(i,j)*GX(1,i,j)*(-1.0D0)
      t(7) = t(20)+t(21)+t(22)+t(23)
      t(8) = 1.0D0/dh2(im2,j)
      t(9) = 1.0D0/ds2(im2,jm1)
      t(10) = 1.0D0/ds2(im2,j)
      t(11) = 1.0D0/t(3)**5
      t(12) = ds1(im3,j)**2
      t(13) = ds1(im2,j)**2
      t(14) = ds1(im1,j)**2
      t(17) = t(6)**2
      t(18) = 1.0D0/ds2(i,jm1)**2
      t(19) = 1.0D0/ds2(i,j)**2
      t(24) = t(7)**2
      t(25) = 1.0D0/dh2(im1,j)
      t(26) = 1.0D0/ds2(im1,jm1)
      t(27) = 1.0D0/ds2(im1,j)
      t(28) = 1.0D0/ds2(i,jm1)
      t(29) = ds1(im2,j)*5.0D-1
      t(41) = ds1(im1,j)*(-1.5D0)
      t(30) = -t(41)
      t(31) = t(29)-t(30)
      t(32) = 1.0D0/t(31)
      t(33) = GX(1,i,j)*(-1.0D0)
      t(34) = t(33)+GX(1,i,jm1)
      t(35) = t(28)*t(34)*5.0D-1
      t(36) = 1.0D0/ds2(i,j)
      t(37) = GX(1,i,jp1)*(-1.0D0)
      t(38) = t(37)+GX(1,i,j)
      t(39) = t(36)*t(38)*5.0D-1
      t(40) = t(35)+t(39)
      t(42) = GX(1,im2,j)*5.0D-1
      t(43) = GX(1,im1,j)*(-2.0D0)
      t(44) = GX(1,i,j)*1.5D0
      t(45) = t(42)+t(43)+t(44)
      t(46) = 1.0D0/dh2(i,j)
      t(47) = ds2(i,j)*(-1.0D0)
      t(48) = ds2(i,jm1)+t(47)
      t(49) = GX(1,im2,j)*(1.0D0/2.0D0)
      t(50) = t(43)+t(44)+t(49)
      t(51) = t(29)+t(41)
      t(52) = 1.0D0/t(51)
      t(53) = 1.0D0/t(51)**3
      t(54) = GX(2,im2,j)*5.0D-1
      t(55) = GX(2,im1,j)*(-2.0D0)
      t(56) = GX(2,i,j)*1.5D0
      t(57) = t(54)+t(55)+t(56)
      t(58) = GX(2,i,j)*(-1.0D0)
      t(59) = t(58)+GX(2,i,jm1)
      t(60) = t(28)*t(59)*5.0D-1
      t(61) = GX(2,i,jp1)*(-1.0D0)
      t(62) = t(61)+GX(2,i,j)
      t(63) = t(36)*t(62)*5.0D-1
      t(64) = t(60)+t(63)
      t(65) = t(28)*5.0D-1
      t(79) = t(36)*(-5.0D-1)
      t(66) = -t(79)
      t(67) = t(65)-t(66)
      t(68) = GX(3,im2,j)*5.0D-1
      t(69) = GX(3,im1,j)*(-2.0D0)
      t(70) = GX(3,i,j)*1.5D0
      t(71) = t(68)+t(69)+t(70)
      t(72) = GX(3,i,j)*(-1.0D0)
      t(73) = t(72)+GX(3,i,jm1)
      t(74) = t(28)*t(73)*5.0D-1
      t(75) = GX(3,i,jp1)*(-1.0D0)
      t(76) = t(75)+GX(3,i,j)
      t(77) = t(36)*t(76)*5.0D-1
      t(78) = t(74)+t(77)
      t(80) = t(65)+t(79)
      t(81) = SigB_d(2,j,1)*t(45)*t(53)*t(57)*(-5.0D-1)
      t(82) = SigB_d(2,j,3)*t(40)*t(52)*t(64)*(-5.0D-1)
      t(83) = t(81)+t(82)
      t(84) = SigB_d(2,j,1)*t(45)*t(53)*t(57)*2.0D0
      t(85) = SigB_d(2,j,3)*t(40)*t(52)*t(64)*2.0D0
      t(86) = t(84)+t(85)
      t(87) = SigB_d(2,j,1)*t(45)*t(53)*t(57)*(-1.5D0)
      t(88) = SigB_d(2,j,3)*t(40)*t(52)*t(64)*(-1.5D0)
      t(89) = gam11(i,j)*t(5)*1.2D1
      t(90) = 1.0D0/t(3)**4
      t(91) = ds1(im3,j)*(-6.0D0)
      t(92) = ds1(im2,j)*1.1D1
      t(93) = ds1(im1,j)*9.0D0
      t(94) = t(91)+t(92)+t(93)
      t(95) = gam11(i,j)*t(90)*t(94)*8.0D0
      t(96) = e5c*gam12(i,j)*t(4)*t(8)*t(9)
      t(97) = sigB(2,j,1)*t(4)*(-1.0D0)
      t(98) = ds1(im4,j)*ds1(im2,j)*3.0D0
      t(99) = ds1(im4,j)*ds1(im1,j)*(-9.0D0)
      t(100) = ds1(im3,j)*ds1(im2,j)*1.3D1
      t(101) = ds1(im3,j)*ds1(im1,j)*1.29D2
      t(102) = ds1(im2,j)*ds1(im1,j)*(-2.84D2)
      t(103) = t(12)*(-1.2D1)
      t(104) = t(13)*2.5D1
      t(105) = t(14)*3.9D1
      t(106) = t(98)+t(99)+t(100)+t(101)+t(102)+t(103)+t(104)+t(10
     &5)
      t(107) = gam11(i,j)*t(11)*t(106)*(-4.0D0)
      t(116) = GX(2,im1,j)*(-4.0D0)
      t(117) = GX(2,i,j)*3.0D0
      t(108) = t(116)+t(117)+GX(2,im2,j)
      t(119) = ds2(i,jm1)*GX(2,i,j)
      t(120) = ds2(i,j)*GX(2,i,jm1)
      t(121) = ds2(i,jm1)*GX(2,i,jp1)*(-1.0D0)
      t(122) = ds2(i,j)*GX(2,i,j)*(-1.0D0)
      t(109) = t(119)+t(120)+t(121)+t(122)
      t(110) = ds2(im2,jm1)+ds2(im2,j)
      t(111) = e5c*gam12(i,j)*t(4)*t(8)*t(9)*t(10)*t(110)*(-1.0D0)
      t(112) = e5c*gam12(i,j)*t(4)*t(8)*t(10)
      t(113) = e5c*gam12(i,j)*t(4)*t(25)*t(26)*(-4.0D0)
      t(114) = sigB(2,j,1)*t(4)*4.0D0
      t(163) = ds1(im4,j)*ds1(im2,j)*6.0D0
      t(164) = ds1(im4,j)*ds1(im1,j)*(-1.8D1)
      t(165) = ds1(im3,j)*ds1(im2,j)*9.2D1
      t(166) = ds1(im3,j)*ds1(im1,j)*6.0D1
      t(167) = ds1(im2,j)*ds1(im1,j)*(-7.6D1)
      t(168) = t(12)*(-2.4D1)
      t(169) = t(13)*(-1.09D2)
      t(170) = t(14)*3.3D1
      t(171) = t(163)+t(164)+t(165)+t(166)+t(167)+t(168)+t(169)+t(
     &170)
      t(115) = gam11(i,j)*t(11)*t(171)*8.0D0
      t(118) = t(108)**2
      t(123) = t(109)**2
      t(124) = ds2(im1,jm1)+ds2(im1,j)
      t(125) = e5c*gam12(i,j)*t(4)*t(25)*t(26)*t(27)*t(124)*4.0D0
      t(126) = e5c*gam12(i,j)*t(4)*t(25)*t(27)*(-4.0D0)
      t(127) = sigB(2,j,3)*t(28)*(-5.0D-1)
      t(128) = sigB(2,j,1)*t(4)*(-3.0D0)
      t(129) = ds1(im4,j)*ds1(im2,j)*9.0D0
      t(130) = ds1(im4,j)*ds1(im1,j)*(-2.7D1)
      t(131) = ds1(im3,j)*ds1(im2,j)*1.59D2
      t(132) = ds1(im3,j)*ds1(im1,j)*2.7D1
      t(133) = ds1(im2,j)*ds1(im1,j)*6.6D1
      t(134) = t(12)*(-3.6D1)
      t(135) = t(13)*(-2.18D2)
      t(136) = t(129)+t(130)+t(131)+t(132)+t(133)+t(134)+t(135)
      t(137) = gam11(i,j)*t(11)*t(136)*(-4.0D0)
      t(138) = sigB(2,j,3)*t(28)*t(36)*t(48)*(-5.0D-1)
      t(139) = ds2(i,jm1)+ds2(i,j)
      t(140) = e5c*gam12(i,j)*t(4)*t(28)*t(36)*t(46)*t(139)*(-3.0D
     &0)
      t(141) = GX(2,im2,j)*(1.0D0/2.0D0)
      t(142) = t(55)+t(56)+t(141)
      t(143) = sigB(2,j,3)*t(36)*5.0D-1
      t(144) = e5c*gam12(i,j)*t(36)*t(46)*t(52)*1.5D0
      t(145) = SigB_d(2,j,1)*t(45)*t(53)*t(71)*(-5.0D-1)
      t(146) = SigB_d(2,j,3)*t(40)*t(52)*t(78)*(-5.0D-1)
      t(147) = t(145)+t(146)
      t(148) = SigB_d(2,j,1)*t(45)*t(53)*t(71)*2.0D0
      t(149) = SigB_d(2,j,3)*t(40)*t(52)*t(78)*2.0D0
      t(150) = t(148)+t(149)
      t(151) = SigB_d(2,j,1)*t(45)*t(53)*t(71)*(-1.5D0)
      t(152) = SigB_d(2,j,3)*t(40)*t(52)*t(78)*(-1.5D0)
      t(153) = SigB_d(2,j,1)*t(53)*t(57)*t(71)*(-5.0D-1)
      t(154) = SigB_d(2,j,3)*t(52)*t(64)*t(78)*(-5.0D-1)
      t(155) = t(153)+t(154)
      t(156) = SigB_d(2,j,1)*t(53)*t(57)*t(71)*2.0D0
      t(157) = SigB_d(2,j,3)*t(52)*t(64)*t(78)*2.0D0
      t(158) = t(156)+t(157)
      t(159) = SigB_d(2,j,1)*t(53)*t(57)*t(71)*(-1.5D0)
      t(160) = SigB_d(2,j,3)*t(52)*t(64)*t(78)*(-1.5D0)
      t(172) = GX(3,im1,j)*(-4.0D0)
      t(173) = GX(3,i,j)*3.0D0
      t(161) = t(172)+t(173)+GX(3,im2,j)
      t(175) = ds2(i,jm1)*GX(3,i,j)
      t(176) = ds2(i,j)*GX(3,i,jm1)
      t(177) = ds2(i,jm1)*GX(3,i,jp1)*(-1.0D0)
      t(178) = ds2(i,j)*GX(3,i,j)*(-1.0D0)
      t(162) = t(175)+t(176)+t(177)+t(178)
      t(174) = t(161)**2
      t(179) = t(162)**2
      t(180) = e5c*gam12(i,j)*t(28)*t(46)*t(52)*1.5D0
      t(181) = GX(3,im2,j)*(1.0D0/2.0D0)
      t(182) = t(69)+t(70)+t(181)
      Y(1,1,-4,0) = t(89)
      Y(1,1,-3,0) = t(95)
      Y(1,1,-2,-1) = t(96)
      Y(1,1,-2,0) = t(97)+t(107)+t(111)-SigB_d(2,j,1)*t(5)*t(17)*1
     &.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(24)*2.5D-1
      Y(1,1,-2,1) = t(112)
      Y(1,1,-1,-1) = t(113)
      Y(1,1,-1,0) = t(114)+t(115)+t(125)+SigB_d(2,j,1)*t(5)*t(17)*
     &4.0D0+SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(24)*1.0D0
      Y(1,1,-1,1) = t(126)
      Y(1,1,0,-1) = t(127)+e5c*gam12(i,j)*t(28)*t(32)*t(46)*1.5D0-
     &SigB_d(2,j,3)*t(28)*t(32)*t(40)*t(45)*5.0D-1-SigB_d(2,j,2)*t(28)*t
     &(40)*t(45)*t(52)*5.0D-1
      Y(1,1,0,0) = t(128)+t(137)+t(138)+t(140)-SigB_d(2,j,1)*t(5)*
     &t(17)*3.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(24)*7.5D-1-SigB_d(2,j
     &,3)*t(4)*t(7)*t(18)*t(19)*t(48)*t(50)*5.0D-1-SigB_d(2,j,2)*t(4)*t(
     &7)*t(18)*t(19)*t(48)*t(50)*5.0D-1
      Y(1,1,0,1) = t(143)+t(144)+SigB_d(2,j,3)*t(36)*t(40)*t(45)*t
     &(52)*5.0D-1+SigB_d(2,j,2)*t(36)*t(40)*t(45)*t(52)*5.0D-1
      Y(1,2,-2,0) = t(83)
      Y(1,2,-1,0) = t(86)
      Y(1,2,0,-1) = SigB_d(2,j,3)*t(28)*t(40)*t(52)*t(57)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(45)*t(52)*t(64)*5.0D-1
      Y(1,2,0,0) = t(87)+t(88)+SigB_d(2,j,3)*t(40)*t(52)*t(57)*t(6
     &7)+SigB_d(2,j,2)*t(45)*t(52)*t(64)*t(67)
      Y(1,2,0,1) = SigB_d(2,j,3)*t(36)*t(40)*t(52)*t(57)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(45)*t(52)*t(64)*5.0D-1
      Y(1,3,-2,0) = t(147)
      Y(1,3,-1,0) = t(150)
      Y(1,3,0,-1) = SigB_d(2,j,3)*t(28)*t(40)*t(52)*t(71)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(45)*t(52)*t(78)*5.0D-1
      Y(1,3,0,0) = t(151)+t(152)+SigB_d(2,j,3)*t(40)*t(52)*t(71)*t
     &(80)+SigB_d(2,j,2)*t(45)*t(52)*t(78)*t(80)
      Y(1,3,0,1) = SigB_d(2,j,3)*t(36)*t(40)*t(52)*t(71)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(45)*t(52)*t(78)*5.0D-1
      Y(2,1,-2,0) = t(83)
      Y(2,1,-1,0) = t(86)
      Y(2,1,0,-1) = SigB_d(2,j,3)*t(28)*t(45)*t(52)*t(64)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(40)*t(52)*t(57)*5.0D-1
      Y(2,1,0,0) = t(87)+t(88)+SigB_d(2,j,3)*t(45)*t(52)*t(64)*t(8
     &0)+SigB_d(2,j,2)*t(40)*t(52)*t(57)*t(80)
      Y(2,1,0,1) = SigB_d(2,j,3)*t(36)*t(45)*t(52)*t(64)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(40)*t(52)*t(57)*5.0D-1
      Y(2,2,-4,0) = t(89)
      Y(2,2,-3,0) = t(95)
      Y(2,2,-2,-1) = t(96)
      Y(2,2,-2,0) = t(97)+t(107)+t(111)-SigB_d(2,j,1)*t(5)*t(118)*
     &1.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(123)*2.5D-1
      Y(2,2,-2,1) = t(112)
      Y(2,2,-1,-1) = t(113)
      Y(2,2,-1,0) = t(114)+t(115)+t(125)+SigB_d(2,j,1)*t(5)*t(118)
     &*4.0D0+SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(123)*1.0D0
      Y(2,2,-1,1) = t(126)
      Y(2,2,0,-1) = t(127)+t(180)-SigB_d(2,j,3)*t(28)*t(52)*t(57)*
     &t(64)*5.0D-1-SigB_d(2,j,2)*t(28)*t(52)*t(57)*t(64)*5.0D-1
      Y(2,2,0,0) = t(128)+t(137)+t(138)+t(140)-SigB_d(2,j,1)*t(5)*
     &t(118)*3.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(123)*7.5D-1-SigB_d(2
     &,j,3)*t(4)*t(18)*t(19)*t(48)*t(109)*t(142)*5.0D-1-SigB_d(2,j,2)*t(
     &4)*t(18)*t(19)*t(48)*t(109)*t(142)*5.0D-1
      Y(2,2,0,1) = t(143)+t(144)+SigB_d(2,j,3)*t(36)*t(52)*t(57)*t
     &(64)*5.0D-1+SigB_d(2,j,2)*t(36)*t(52)*t(57)*t(64)*5.0D-1
      Y(2,3,-2,0) = t(155)
      Y(2,3,-1,0) = t(158)
      Y(2,3,0,-1) = SigB_d(2,j,3)*t(28)*t(52)*t(64)*t(71)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(52)*t(57)*t(78)*5.0D-1
      Y(2,3,0,0) = t(159)+t(160)+SigB_d(2,j,3)*t(52)*t(64)*t(71)*t
     &(80)+SigB_d(2,j,2)*t(52)*t(57)*t(78)*t(80)
      Y(2,3,0,1) = SigB_d(2,j,3)*t(36)*t(52)*t(64)*t(71)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(52)*t(57)*t(78)*5.0D-1
      Y(3,1,-2,0) = t(147)
      Y(3,1,-1,0) = t(150)
      Y(3,1,0,-1) = SigB_d(2,j,3)*t(28)*t(45)*t(52)*t(78)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(40)*t(52)*t(71)*5.0D-1
      Y(3,1,0,0) = t(151)+t(152)+SigB_d(2,j,3)*t(45)*t(52)*t(78)*t
     &(80)+SigB_d(2,j,2)*t(40)*t(52)*t(71)*t(80)
      Y(3,1,0,1) = SigB_d(2,j,3)*t(36)*t(45)*t(52)*t(78)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(40)*t(52)*t(71)*5.0D-1
      Y(3,2,-2,0) = t(155)
      Y(3,2,-1,0) = t(158)
      Y(3,2,0,-1) = SigB_d(2,j,3)*t(28)*t(52)*t(57)*t(78)*(-5.0D-1
     &)-SigB_d(2,j,2)*t(28)*t(52)*t(64)*t(71)*5.0D-1
      Y(3,2,0,0) = t(159)+t(160)+SigB_d(2,j,3)*t(52)*t(57)*t(78)*t
     &(80)+SigB_d(2,j,2)*t(52)*t(64)*t(71)*t(80)
      Y(3,2,0,1) = SigB_d(2,j,3)*t(36)*t(52)*t(57)*t(78)*5.0D-1+Si
     &gB_d(2,j,2)*t(36)*t(52)*t(64)*t(71)*5.0D-1
      Y(3,3,-4,0) = t(89)
      Y(3,3,-3,0) = t(95)
      Y(3,3,-2,-1) = t(96)
      Y(3,3,-2,0) = t(97)+t(107)+t(111)-SigB_d(2,j,1)*t(5)*t(174)*
     &1.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(179)*2.5D-1
      Y(3,3,-2,1) = t(112)
      Y(3,3,-1,-1) = t(113)
      Y(3,3,-1,0) = t(114)+t(115)+t(125)+SigB_d(2,j,1)*t(5)*t(174)
     &*4.0D0+SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(179)*1.0D0
      Y(3,3,-1,1) = t(126)
      Y(3,3,0,-1) = t(127)+t(180)-SigB_d(2,j,3)*t(28)*t(52)*t(71)*
     &t(78)*5.0D-1-SigB_d(2,j,2)*t(28)*t(52)*t(71)*t(78)*5.0D-1
      Y(3,3,0,0) = t(128)+t(137)+t(138)+t(140)-SigB_d(2,j,1)*t(5)*
     &t(174)*3.0D0-SigB_d(2,j,3)*t(4)*t(18)*t(19)*t(179)*7.5D-1-SigB_d(2
     &,j,3)*t(4)*t(18)*t(19)*t(48)*t(162)*t(182)*5.0D-1-SigB_d(2,j,2)*t(
     &4)*t(18)*t(19)*t(48)*t(162)*t(182)*5.0D-1
      Y(3,3,0,1) = t(143)+t(144)+SigB_d(2,j,3)*t(36)*t(52)*t(71)*t
     &(78)*5.0D-1+SigB_d(2,j,2)*t(36)*t(52)*t(71)*t(78)*5.0D-1

      endif
       Y(1:3,1:3,-4:4,-4:4)=-Y(1:3,1:3,-4:4,-4:4)

!! ================== BNDY i,1 2<=i<=Nx-1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .ge. ist) .and. (i.le. ien)) then

      if(formulationFlagb .eq. 1) then
      t(2) = 1.0D0/ds1(im1,j)
      t(3) = ds2(i,j)*1.5D0
      t(15) = ds2(i,jp1)*(-5.0D-1)
      t(4) = -t(15)
      t(5) = t(3)-t(4)
      t(6) = 1.0D0/t(5)
      t(7) = GX(1,i,j)*(-1.0D0)
      t(8) = t(7)+GX(1,im1,j)
      t(9) = t(2)*t(8)*5.0D-1
      t(10) = 1.0D0/ds1(i,j)
      t(11) = GX(1,ip1,j)*(-1.0D0)
      t(12) = t(11)+GX(1,i,j)
      t(13) = t(10)*t(12)*5.0D-1
      t(14) = t(9)+t(13)
      t(16) = GX(1,i,j)*1.5D0
      t(17) = GX(1,i,jp1)*(-2.0D0)
      t(18) = GX(1,i,jp2)*5.0D-1
      t(19) = t(16)+t(17)+t(18)
      t(20) = t(3)+t(15)
      t(21) = 1.0D0/t(20)
      t(43) = GX(1,i,j)*3.0D0
      t(44) = GX(1,i,jp1)*(-4.0D0)
      t(22) = t(43)+t(44)+GX(1,i,jp2)
      t(23) = ds2(i,j)*3.0D0
      t(28) = ds2(i,jp1)*(-1.0D0)
      t(24) = -t(28)
      t(25) = t(23)-t(24)
      t(33) = ds1(im1,j)*GX(1,i,j)
      t(34) = ds1(i,j)*GX(1,im1,j)
      t(35) = ds1(im1,j)*GX(1,ip1,j)*(-1.0D0)
      t(36) = ds1(i,j)*GX(1,i,j)*(-1.0D0)
      t(26) = t(33)+t(34)+t(35)+t(36)
      t(27) = 1.0D0/dh1(i,j)
      t(29) = 1.0D0/ds1(im1,j)**2
      t(30) = 1.0D0/ds1(i,j)**2
      t(31) = ds1(i,j)*(-1.0D0)
      t(32) = ds1(im1,j)+t(31)
      t(37) = GX(1,i,jp2)*(1.0D0/2.0D0)
      t(38) = t(16)+t(17)+t(37)
      t(39) = t(23)+t(28)
      t(40) = ds2(i,jp1)**2
      t(41) = ds2(i,jp2)**2
      t(42) = 1.0D0/t(20)**3
      t(45) = t(22)**2
      t(46) = t(26)**2
      t(47) = 1.0D0/dh1(i,jp1)
      t(48) = 1.0D0/ds1(im1,jp1)
      t(49) = 1.0D0/t(39)
      t(50) = ds2(i,j)**2
      t(51) = 1.0D0/dh1(i,jp2)
      t(52) = 1.0D0/ds1(im1,jp2)
      t(53) = 1.0D0/ds1(i,jp1)
      t(54) = 1.0D0/ds1(i,jp2)
      t(55) = GX(2,i,j)*1.5D0
      t(56) = GX(2,i,jp1)*(-2.0D0)
      t(57) = GX(2,i,jp2)*5.0D-1
      t(58) = t(55)+t(56)+t(57)
      t(59) = GX(2,i,j)*(-1.0D0)
      t(60) = t(59)+GX(2,im1,j)
      t(61) = t(2)*t(60)*5.0D-1
      t(62) = GX(2,ip1,j)*(-1.0D0)
      t(63) = t(62)+GX(2,i,j)
      t(64) = t(10)*t(63)*5.0D-1
      t(65) = t(61)+t(64)
      t(66) = t(2)*5.0D-1
      t(80) = t(10)*(-5.0D-1)
      t(67) = -t(80)
      t(68) = t(66)-t(67)
      t(69) = GX(3,i,j)*1.5D0
      t(70) = GX(3,i,jp1)*(-2.0D0)
      t(71) = GX(3,i,jp2)*5.0D-1
      t(72) = t(69)+t(70)+t(71)
      t(73) = GX(3,i,j)*(-1.0D0)
      t(74) = t(73)+GX(3,im1,j)
      t(75) = t(2)*t(74)*5.0D-1
      t(76) = GX(3,ip1,j)*(-1.0D0)
      t(77) = t(76)+GX(3,i,j)
      t(78) = t(10)*t(77)*5.0D-1
      t(79) = t(75)+t(78)
      t(81) = t(66)+t(80)
      t(82) = SigB_d(3,i,2)*t(19)*t(42)*t(58)*(-1.5D0)
      t(83) = SigB_d(3,i,3)*t(14)*t(21)*t(65)*(-1.5D0)
      t(84) = SigB_d(3,i,2)*t(19)*t(42)*t(58)*2.0D0
      t(85) = SigB_d(3,i,3)*t(14)*t(21)*t(65)*2.0D0
      t(86) = t(84)+t(85)
      t(87) = SigB_d(3,i,2)*t(19)*t(42)*t(58)*(-5.0D-1)
      t(88) = SigB_d(3,i,3)*t(14)*t(21)*t(65)*(-5.0D-1)
      t(89) = t(87)+t(88)
      t(90) = sigB(3,i,3)*t(2)*(-5.0D-1)
      t(91) = e5c*gam12(i,j)*t(21)*t(47)*t(48)*(-2.0D0)
      t(92) = sigB(3,i,2)*t(21)*(-1.5D0)
      t(124) = GX(2,i,j)*3.0D0
      t(125) = GX(2,i,jp1)*(-4.0D0)
      t(93) = t(124)+t(125)+GX(2,i,jp2)
      t(94) = 1.0D0/t(39)**5
      t(95) = ds2(i,j)*ds2(i,jp1)*6.6D1
      t(96) = ds2(i,j)*ds2(i,jp2)*2.7D1
      t(97) = ds2(i,j)*ds2(i,jp3)*(-2.7D1)
      t(98) = ds2(i,jp1)*ds2(i,jp2)*1.59D2
      t(99) = ds2(i,jp1)*ds2(i,jp3)*9.0D0
      t(100) = t(40)*(-2.18D2)
      t(101) = t(41)*(-3.6D1)
      t(102) = t(95)+t(96)+t(97)+t(98)+t(99)+t(100)+t(101)
      t(103) = sigB(3,i,3)*t(2)*t(10)*t(32)*(-5.0D-1)
      t(107) = ds1(im1,j)*GX(2,i,j)
      t(108) = ds1(i,j)*GX(2,im1,j)
      t(109) = ds1(im1,j)*GX(2,ip1,j)*(-1.0D0)
      t(110) = ds1(i,j)*GX(2,i,j)*(-1.0D0)
      t(104) = t(107)+t(108)+t(109)+t(110)
      t(105) = ds1(im1,j)+ds1(i,j)
      t(106) = e5c*gam12(i,j)*t(2)*t(10)*t(27)*t(49)*t(105)*(-3.0D
     &0)
      t(111) = GX(2,i,jp2)*(1.0D0/2.0D0)
      t(112) = t(55)+t(56)+t(111)
      t(113) = sigB(3,i,2)*t(21)*2.0D0
      t(114) = ds2(i,j)*ds2(i,jp1)*7.6D1
      t(115) = ds2(i,j)*ds2(i,jp2)*(-6.0D1)
      t(116) = ds2(i,j)*ds2(i,jp3)*1.8D1
      t(117) = ds2(i,jp1)*ds2(i,jp2)*(-9.2D1)
      t(118) = ds2(i,jp1)*ds2(i,jp3)*(-6.0D0)
      t(119) = t(50)*(-3.3D1)
      t(120) = t(40)*1.09D2
      t(121) = t(41)*2.4D1
      t(122) = t(114)+t(115)+t(116)+t(117)+t(118)+t(119)+t(120)+t(
     &121)
      t(123) = gam11(i,j)*t(94)*t(122)*(-8.0D0)
      t(126) = t(93)**2
      t(127) = t(104)**2
      t(128) = ds1(im1,jp1)+ds1(i,jp1)
      t(129) = e5c*gam12(i,j)*t(47)*t(48)*t(49)*t(53)*t(128)*4.0D0
      t(130) = sigB(3,i,2)*t(21)*(-5.0D-1)
      t(131) = 1.0D0/t(39)**2
      t(132) = ds2(i,j)*ds2(i,jp1)*(-2.84D2)
      t(133) = ds2(i,j)*ds2(i,jp2)*1.29D2
      t(134) = ds2(i,j)*ds2(i,jp3)*(-9.0D0)
      t(135) = ds2(i,jp1)*ds2(i,jp2)*1.3D1
      t(136) = ds2(i,jp1)*ds2(i,jp3)*3.0D0
      t(137) = t(50)*3.9D1
      t(138) = t(40)*2.5D1
      t(139) = t(41)*(-1.2D1)
      t(140) = t(132)+t(133)+t(134)+t(135)+t(136)+t(137)+t(138)+t(
     &139)
      t(141) = gam11(i,j)*t(42)*t(131)*t(140)*(-5.0D-1)
      t(142) = ds1(im1,jp2)+ds1(i,jp2)
      t(143) = e5c*gam12(i,j)*t(49)*t(51)*t(52)*t(54)*t(142)*(-1.0
     &D0)
      t(144) = ds2(i,j)*9.0D0
      t(145) = ds2(i,jp1)*1.1D1
      t(146) = ds2(i,jp2)*(-6.0D0)
      t(147) = t(144)+t(145)+t(146)
      t(148) = gam11(i,j)*t(42)*t(49)*t(147)*(-1.0D0)
      t(149) = gam11(i,j)*t(42)*1.5D0
      t(150) = sigB(3,i,3)*t(10)*5.0D-1
      t(151) = e5c*gam12(i,j)*t(10)*t(21)*t(27)*1.5D0
      t(152) = e5c*gam12(i,j)*t(21)*t(47)*t(53)*(-2.0D0)
      t(153) = e5c*gam12(i,j)*t(49)*t(51)*t(54)
      t(154) = SigB_d(3,i,2)*t(19)*t(42)*t(72)*(-1.5D0)
      t(155) = SigB_d(3,i,3)*t(14)*t(21)*t(79)*(-1.5D0)
      t(156) = SigB_d(3,i,2)*t(19)*t(42)*t(72)*2.0D0
      t(157) = SigB_d(3,i,3)*t(14)*t(21)*t(79)*2.0D0
      t(158) = t(156)+t(157)
      t(159) = SigB_d(3,i,2)*t(19)*t(42)*t(72)*(-5.0D-1)
      t(160) = SigB_d(3,i,3)*t(14)*t(21)*t(79)*(-5.0D-1)
      t(161) = t(159)+t(160)
      t(162) = SigB_d(3,i,2)*t(42)*t(58)*t(72)*(-1.5D0)
      t(163) = SigB_d(3,i,3)*t(21)*t(65)*t(79)*(-1.5D0)
      t(164) = SigB_d(3,i,2)*t(42)*t(58)*t(72)*2.0D0
      t(165) = SigB_d(3,i,3)*t(21)*t(65)*t(79)*2.0D0
      t(166) = t(164)+t(165)
      t(167) = SigB_d(3,i,2)*t(42)*t(58)*t(72)*(-5.0D-1)
      t(168) = SigB_d(3,i,3)*t(21)*t(65)*t(79)*(-5.0D-1)
      t(169) = t(167)+t(168)
      t(170) = e5c*gam12(i,j)*t(2)*t(21)*t(27)*1.5D0
      t(171) = e5c*gam12(i,j)*t(49)*t(51)*t(52)
      t(181) = GX(3,i,j)*3.0D0
      t(182) = GX(3,i,jp1)*(-4.0D0)
      t(172) = t(181)+t(182)+GX(3,i,jp2)
      t(173) = gam11(i,j)*t(94)*t(102)*(-4.0D0)
      t(175) = ds1(im1,j)*GX(3,i,j)
      t(176) = ds1(i,j)*GX(3,im1,j)
      t(177) = ds1(im1,j)*GX(3,ip1,j)*(-1.0D0)
      t(178) = ds1(i,j)*GX(3,i,j)*(-1.0D0)
      t(174) = t(175)+t(176)+t(177)+t(178)
      t(179) = GX(3,i,jp2)*(1.0D0/2.0D0)
      t(180) = t(69)+t(70)+t(179)
      t(183) = t(172)**2
      t(184) = t(174)**2
      Y(1,1,-1,0) = t(90)+e5c*gam12(i,j)*t(2)*t(6)*t(27)*1.5D0-Sig
     &B_d(3,i,3)*t(2)*t(6)*t(14)*t(19)*5.0D-1-SigB_d(3,i,1)*t(2)*t(14)*t
     &(19)*t(21)*5.0D-1
      Y(1,1,-1,1) = t(91)
      Y(1,1,-1,2) = (e5c*gam12(i,j)*t(51)*t(52))/t(25)
      Y(1,1,0,0) = t(92)+t(103)+t(106)-gam11(i,j)*1.0D0/t(25)**5*t
     &(102)*4.0D0-SigB_d(3,i,2)*t(42)*t(45)*3.75D-1-SigB_d(3,i,3)*t(21)*
     &t(29)*t(30)*t(46)*3.75D-1-SigB_d(3,i,3)*t(21)*t(26)*t(29)*t(30)*t(
     &32)*t(38)*2.5D-1-SigB_d(3,i,1)*t(21)*t(26)*t(29)*t(30)*t(32)*t(38)
     &*2.5D-1
      Y(1,1,0,1) = t(113)+t(123)+t(129)+SigB_d(3,i,2)*t(42)*t(45)*
     &5.0D-1+SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(46)*5.0D-1
      Y(1,1,0,2) = t(130)+t(141)+t(143)-SigB_d(3,i,2)*t(42)*t(45)*
     &1.25D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(46)*1.25D-1
      Y(1,1,0,3) = t(148)
      Y(1,1,0,4) = t(149)
      Y(1,1,1,0) = t(150)+t(151)+SigB_d(3,i,3)*t(10)*t(14)*t(19)*t
     &(21)*5.0D-1+SigB_d(3,i,1)*t(10)*t(14)*t(19)*t(21)*5.0D-1
      Y(1,1,1,1) = t(152)
      Y(1,1,1,2) = t(153)
      Y(1,2,-1,0) = SigB_d(3,i,3)*t(2)*t(14)*t(21)*t(58)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(19)*t(21)*t(65)*5.0D-1
      Y(1,2,0,0) = t(82)+t(83)+SigB_d(3,i,3)*t(14)*t(21)*t(58)*t(6
     &8)+SigB_d(3,i,1)*t(19)*t(21)*t(65)*t(68)
      Y(1,2,0,1) = t(86)
      Y(1,2,0,2) = t(89)
      Y(1,2,1,0) = SigB_d(3,i,3)*t(10)*t(14)*t(21)*t(58)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(19)*t(21)*t(65)*5.0D-1
      Y(1,3,-1,0) = SigB_d(3,i,3)*t(2)*t(14)*t(21)*t(72)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(19)*t(21)*t(79)*5.0D-1
      Y(1,3,0,0) = t(154)+t(155)+SigB_d(3,i,3)*t(14)*t(21)*t(72)*t
     &(81)+SigB_d(3,i,1)*t(19)*t(21)*t(79)*t(81)
      Y(1,3,0,1) = t(158)
      Y(1,3,0,2) = t(161)
      Y(1,3,1,0) = SigB_d(3,i,3)*t(10)*t(14)*t(21)*t(72)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(19)*t(21)*t(79)*5.0D-1
      Y(2,1,-1,0) = SigB_d(3,i,3)*t(2)*t(19)*t(21)*t(65)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(14)*t(21)*t(58)*5.0D-1
      Y(2,1,0,0) = t(82)+t(83)+SigB_d(3,i,3)*t(19)*t(21)*t(65)*t(8
     &1)+SigB_d(3,i,1)*t(14)*t(21)*t(58)*t(81)
      Y(2,1,0,1) = t(86)
      Y(2,1,0,2) = t(89)
      Y(2,1,1,0) = SigB_d(3,i,3)*t(10)*t(19)*t(21)*t(65)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(14)*t(21)*t(58)*5.0D-1
      Y(2,2,-1,0) = t(90)+t(170)-SigB_d(3,i,3)*t(2)*t(21)*t(58)*t(
     &65)*5.0D-1-SigB_d(3,i,1)*t(2)*t(21)*t(58)*t(65)*5.0D-1
      Y(2,2,-1,1) = t(91)
      Y(2,2,-1,2) = t(171)
      Y(2,2,0,0) = t(92)+t(103)+t(106)+t(173)-SigB_d(3,i,2)*t(42)*
     &t(126)*3.75D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(127)*3.75D-1-SigB
     &_d(3,i,3)*t(21)*t(29)*t(30)*t(32)*t(104)*t(112)*2.5D-1-SigB_d(3,i,
     &1)*t(21)*t(29)*t(30)*t(32)*t(104)*t(112)*2.5D-1
      Y(2,2,0,1) = t(113)+t(123)+t(129)+SigB_d(3,i,2)*t(42)*t(126)
     &*5.0D-1+SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(127)*5.0D-1
      Y(2,2,0,2) = t(130)+t(141)+t(143)-SigB_d(3,i,2)*t(42)*t(126)
     &*1.25D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(127)*1.25D-1
      Y(2,2,0,3) = t(148)
      Y(2,2,0,4) = t(149)
      Y(2,2,1,0) = t(150)+t(151)+SigB_d(3,i,3)*t(10)*t(21)*t(58)*t
     &(65)*5.0D-1+SigB_d(3,i,1)*t(10)*t(21)*t(58)*t(65)*5.0D-1
      Y(2,2,1,1) = t(152)
      Y(2,2,1,2) = t(153)
      Y(2,3,-1,0) = SigB_d(3,i,3)*t(2)*t(21)*t(65)*t(72)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(21)*t(58)*t(79)*5.0D-1
      Y(2,3,0,0) = t(162)+t(163)+SigB_d(3,i,3)*t(21)*t(65)*t(72)*t
     &(81)+SigB_d(3,i,1)*t(21)*t(58)*t(79)*t(81)
      Y(2,3,0,1) = t(166)
      Y(2,3,0,2) = t(169)
      Y(2,3,1,0) = SigB_d(3,i,3)*t(10)*t(21)*t(65)*t(72)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(21)*t(58)*t(79)*5.0D-1
      Y(3,1,-1,0) = SigB_d(3,i,3)*t(2)*t(19)*t(21)*t(79)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(14)*t(21)*t(72)*5.0D-1
      Y(3,1,0,0) = t(154)+t(155)+SigB_d(3,i,3)*t(19)*t(21)*t(79)*t
     &(81)+SigB_d(3,i,1)*t(14)*t(21)*t(72)*t(81)
      Y(3,1,0,1) = t(158)
      Y(3,1,0,2) = t(161)
      Y(3,1,1,0) = SigB_d(3,i,3)*t(10)*t(19)*t(21)*t(79)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(14)*t(21)*t(72)*5.0D-1
      Y(3,2,-1,0) = SigB_d(3,i,3)*t(2)*t(21)*t(58)*t(79)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(21)*t(65)*t(72)*5.0D-1
      Y(3,2,0,0) = t(162)+t(163)+SigB_d(3,i,3)*t(21)*t(58)*t(79)*t
     &(81)+SigB_d(3,i,1)*t(21)*t(65)*t(72)*t(81)
      Y(3,2,0,1) = t(166)
      Y(3,2,0,2) = t(169)
      Y(3,2,1,0) = SigB_d(3,i,3)*t(10)*t(21)*t(58)*t(79)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(21)*t(65)*t(72)*5.0D-1
      Y(3,3,-1,0) = t(90)+t(170)-SigB_d(3,i,3)*t(2)*t(21)*t(72)*t(
     &79)*5.0D-1-SigB_d(3,i,1)*t(2)*t(21)*t(72)*t(79)*5.0D-1
      Y(3,3,-1,1) = t(91)
      Y(3,3,-1,2) = t(171)
      Y(3,3,0,0) = t(92)+t(103)+t(106)+t(173)-SigB_d(3,i,2)*t(42)*
     &t(183)*3.75D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(184)*3.75D-1-SigB
     &_d(3,i,3)*t(21)*t(29)*t(30)*t(32)*t(174)*t(180)*2.5D-1-SigB_d(3,i,
     &1)*t(21)*t(29)*t(30)*t(32)*t(174)*t(180)*2.5D-1
      Y(3,3,0,1) = t(113)+t(123)+t(129)+SigB_d(3,i,2)*t(42)*t(183)
     &*5.0D-1+SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(184)*5.0D-1
      Y(3,3,0,2) = t(130)+t(141)+t(143)-SigB_d(3,i,2)*t(42)*t(183)
     &*1.25D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(184)*1.25D-1
      Y(3,3,0,3) = t(148)
      Y(3,3,0,4) = t(149)
      Y(3,3,1,0) = t(150)+t(151)+SigB_d(3,i,3)*t(10)*t(21)*t(72)*t
     &(79)*5.0D-1+SigB_d(3,i,1)*t(10)*t(21)*t(72)*t(79)*5.0D-1
      Y(3,3,1,1) = t(152)
      Y(3,3,1,2) = t(153)

      else
      t(2) = 1.0D0/ds1(im1,j)
      t(3) = ds2(i,j)*1.5D0
      t(15) = ds2(i,jp1)*(-5.0D-1)
      t(4) = -t(15)
      t(5) = t(3)-t(4)
      t(6) = 1.0D0/t(5)
      t(7) = GX(1,i,j)*(-1.0D0)
      t(8) = t(7)+GX(1,im1,j)
      t(9) = t(2)*t(8)*5.0D-1
      t(10) = 1.0D0/ds1(i,j)
      t(11) = GX(1,ip1,j)*(-1.0D0)
      t(12) = t(11)+GX(1,i,j)
      t(13) = t(10)*t(12)*5.0D-1
      t(14) = t(9)+t(13)
      t(16) = GX(1,i,j)*1.5D0
      t(17) = GX(1,i,jp1)*(-2.0D0)
      t(18) = GX(1,i,jp2)*5.0D-1
      t(19) = t(16)+t(17)+t(18)
      t(20) = t(3)+t(15)
      t(21) = 1.0D0/t(20)
      t(43) = GX(1,i,j)*3.0D0
      t(44) = GX(1,i,jp1)*(-4.0D0)
      t(22) = t(43)+t(44)+GX(1,i,jp2)
      t(23) = ds2(i,j)*3.0D0
      t(28) = ds2(i,jp1)*(-1.0D0)
      t(24) = -t(28)
      t(25) = t(23)-t(24)
      t(33) = ds1(im1,j)*GX(1,i,j)
      t(34) = ds1(i,j)*GX(1,im1,j)
      t(35) = ds1(im1,j)*GX(1,ip1,j)*(-1.0D0)
      t(36) = ds1(i,j)*GX(1,i,j)*(-1.0D0)
      t(26) = t(33)+t(34)+t(35)+t(36)
      t(27) = 1.0D0/dh1(i,j)
      t(29) = 1.0D0/ds1(im1,j)**2
      t(30) = 1.0D0/ds1(i,j)**2
      t(31) = ds1(i,j)*(-1.0D0)
      t(32) = ds1(im1,j)+t(31)
      t(37) = GX(1,i,jp2)*(1.0D0/2.0D0)
      t(38) = t(16)+t(17)+t(37)
      t(39) = t(23)+t(28)
      t(40) = ds2(i,jp1)**2
      t(41) = ds2(i,jp2)**2
      t(42) = 1.0D0/t(20)**3
      t(45) = t(22)**2
      t(46) = t(26)**2
      t(47) = 1.0D0/dh1(i,jp1)
      t(48) = 1.0D0/ds1(im1,jp1)
      t(49) = 1.0D0/t(39)
      t(50) = ds2(i,j)**2
      t(51) = 1.0D0/dh1(i,jp2)
      t(52) = 1.0D0/ds1(im1,jp2)
      t(53) = 1.0D0/ds1(i,jp1)
      t(54) = 1.0D0/ds1(i,jp2)
      t(55) = GX(2,i,j)*1.5D0
      t(56) = GX(2,i,jp1)*(-2.0D0)
      t(57) = GX(2,i,jp2)*5.0D-1
      t(58) = t(55)+t(56)+t(57)
      t(59) = GX(2,i,j)*(-1.0D0)
      t(60) = t(59)+GX(2,im1,j)
      t(61) = t(2)*t(60)*5.0D-1
      t(62) = GX(2,ip1,j)*(-1.0D0)
      t(63) = t(62)+GX(2,i,j)
      t(64) = t(10)*t(63)*5.0D-1
      t(65) = t(61)+t(64)
      t(66) = t(2)*5.0D-1
      t(80) = t(10)*(-5.0D-1)
      t(67) = -t(80)
      t(68) = t(66)-t(67)
      t(69) = GX(3,i,j)*1.5D0
      t(70) = GX(3,i,jp1)*(-2.0D0)
      t(71) = GX(3,i,jp2)*5.0D-1
      t(72) = t(69)+t(70)+t(71)
      t(73) = GX(3,i,j)*(-1.0D0)
      t(74) = t(73)+GX(3,im1,j)
      t(75) = t(2)*t(74)*5.0D-1
      t(76) = GX(3,ip1,j)*(-1.0D0)
      t(77) = t(76)+GX(3,i,j)
      t(78) = t(10)*t(77)*5.0D-1
      t(79) = t(75)+t(78)
      t(81) = t(66)+t(80)
      t(82) = SigB_d(3,i,2)*t(19)*t(42)*t(58)*(-1.5D0)
      t(83) = SigB_d(3,i,3)*t(14)*t(21)*t(65)*(-1.5D0)
      t(84) = SigB_d(3,i,2)*t(19)*t(42)*t(58)*2.0D0
      t(85) = SigB_d(3,i,3)*t(14)*t(21)*t(65)*2.0D0
      t(86) = t(84)+t(85)
      t(87) = SigB_d(3,i,2)*t(19)*t(42)*t(58)*(-5.0D-1)
      t(88) = SigB_d(3,i,3)*t(14)*t(21)*t(65)*(-5.0D-1)
      t(89) = t(87)+t(88)
      t(90) = sigB(3,i,3)*t(2)*(-5.0D-1)
      t(91) = e5c*gam12(i,j)*t(21)*t(47)*t(48)*(-2.0D0)
      t(92) = sigB(3,i,2)*t(21)*(-1.5D0)
      t(124) = GX(2,i,j)*3.0D0
      t(125) = GX(2,i,jp1)*(-4.0D0)
      t(93) = t(124)+t(125)+GX(2,i,jp2)
      t(94) = 1.0D0/t(39)**5
      t(95) = ds2(i,j)*ds2(i,jp1)*6.6D1
      t(96) = ds2(i,j)*ds2(i,jp2)*2.7D1
      t(97) = ds2(i,j)*ds2(i,jp3)*(-2.7D1)
      t(98) = ds2(i,jp1)*ds2(i,jp2)*1.59D2
      t(99) = ds2(i,jp1)*ds2(i,jp3)*9.0D0
      t(100) = t(40)*(-2.18D2)
      t(101) = t(41)*(-3.6D1)
      t(102) = t(95)+t(96)+t(97)+t(98)+t(99)+t(100)+t(101)
      t(103) = sigB(3,i,3)*t(2)*t(10)*t(32)*(-5.0D-1)
      t(107) = ds1(im1,j)*GX(2,i,j)
      t(108) = ds1(i,j)*GX(2,im1,j)
      t(109) = ds1(im1,j)*GX(2,ip1,j)*(-1.0D0)
      t(110) = ds1(i,j)*GX(2,i,j)*(-1.0D0)
      t(104) = t(107)+t(108)+t(109)+t(110)
      t(105) = ds1(im1,j)+ds1(i,j)
      t(106) = e5c*gam12(i,j)*t(2)*t(10)*t(27)*t(49)*t(105)*(-3.0D
     &0)
      t(111) = GX(2,i,jp2)*(1.0D0/2.0D0)
      t(112) = t(55)+t(56)+t(111)
      t(113) = sigB(3,i,2)*t(21)*2.0D0
      t(114) = ds2(i,j)*ds2(i,jp1)*7.6D1
      t(115) = ds2(i,j)*ds2(i,jp2)*(-6.0D1)
      t(116) = ds2(i,j)*ds2(i,jp3)*1.8D1
      t(117) = ds2(i,jp1)*ds2(i,jp2)*(-9.2D1)
      t(118) = ds2(i,jp1)*ds2(i,jp3)*(-6.0D0)
      t(119) = t(50)*(-3.3D1)
      t(120) = t(40)*1.09D2
      t(121) = t(41)*2.4D1
      t(122) = t(114)+t(115)+t(116)+t(117)+t(118)+t(119)+t(120)+t(
     &121)
      t(123) = gam11(i,j)*t(94)*t(122)*(-8.0D0)
      t(126) = t(93)**2
      t(127) = t(104)**2
      t(128) = ds1(im1,jp1)+ds1(i,jp1)
      t(129) = e5c*gam12(i,j)*t(47)*t(48)*t(49)*t(53)*t(128)*4.0D0
      t(130) = sigB(3,i,2)*t(21)*(-5.0D-1)
      t(131) = 1.0D0/t(39)**2
      t(132) = ds2(i,j)*ds2(i,jp1)*(-2.84D2)
      t(133) = ds2(i,j)*ds2(i,jp2)*1.29D2
      t(134) = ds2(i,j)*ds2(i,jp3)*(-9.0D0)
      t(135) = ds2(i,jp1)*ds2(i,jp2)*1.3D1
      t(136) = ds2(i,jp1)*ds2(i,jp3)*3.0D0
      t(137) = t(50)*3.9D1
      t(138) = t(40)*2.5D1
      t(139) = t(41)*(-1.2D1)
      t(140) = t(132)+t(133)+t(134)+t(135)+t(136)+t(137)+t(138)+t(
     &139)
      t(141) = gam11(i,j)*t(42)*t(131)*t(140)*(-5.0D-1)
      t(142) = ds1(im1,jp2)+ds1(i,jp2)
      t(143) = e5c*gam12(i,j)*t(49)*t(51)*t(52)*t(54)*t(142)*(-1.0
     &D0)
      t(144) = ds2(i,j)*9.0D0
      t(145) = ds2(i,jp1)*1.1D1
      t(146) = ds2(i,jp2)*(-6.0D0)
      t(147) = t(144)+t(145)+t(146)
      t(148) = gam11(i,j)*t(42)*t(49)*t(147)*(-1.0D0)
      t(149) = gam11(i,j)*t(42)*1.5D0
      t(150) = sigB(3,i,3)*t(10)*5.0D-1
      t(151) = e5c*gam12(i,j)*t(10)*t(21)*t(27)*1.5D0
      t(152) = e5c*gam12(i,j)*t(21)*t(47)*t(53)*(-2.0D0)
      t(153) = e5c*gam12(i,j)*t(49)*t(51)*t(54)
      t(154) = SigB_d(3,i,2)*t(19)*t(42)*t(72)*(-1.5D0)
      t(155) = SigB_d(3,i,3)*t(14)*t(21)*t(79)*(-1.5D0)
      t(156) = SigB_d(3,i,2)*t(19)*t(42)*t(72)*2.0D0
      t(157) = SigB_d(3,i,3)*t(14)*t(21)*t(79)*2.0D0
      t(158) = t(156)+t(157)
      t(159) = SigB_d(3,i,2)*t(19)*t(42)*t(72)*(-5.0D-1)
      t(160) = SigB_d(3,i,3)*t(14)*t(21)*t(79)*(-5.0D-1)
      t(161) = t(159)+t(160)
      t(162) = SigB_d(3,i,2)*t(42)*t(58)*t(72)*(-1.5D0)
      t(163) = SigB_d(3,i,3)*t(21)*t(65)*t(79)*(-1.5D0)
      t(164) = SigB_d(3,i,2)*t(42)*t(58)*t(72)*2.0D0
      t(165) = SigB_d(3,i,3)*t(21)*t(65)*t(79)*2.0D0
      t(166) = t(164)+t(165)
      t(167) = SigB_d(3,i,2)*t(42)*t(58)*t(72)*(-5.0D-1)
      t(168) = SigB_d(3,i,3)*t(21)*t(65)*t(79)*(-5.0D-1)
      t(169) = t(167)+t(168)
      t(170) = e5c*gam12(i,j)*t(2)*t(21)*t(27)*1.5D0
      t(171) = e5c*gam12(i,j)*t(49)*t(51)*t(52)
      t(181) = GX(3,i,j)*3.0D0
      t(182) = GX(3,i,jp1)*(-4.0D0)
      t(172) = t(181)+t(182)+GX(3,i,jp2)
      t(173) = gam11(i,j)*t(94)*t(102)*(-4.0D0)
      t(175) = ds1(im1,j)*GX(3,i,j)
      t(176) = ds1(i,j)*GX(3,im1,j)
      t(177) = ds1(im1,j)*GX(3,ip1,j)*(-1.0D0)
      t(178) = ds1(i,j)*GX(3,i,j)*(-1.0D0)
      t(174) = t(175)+t(176)+t(177)+t(178)
      t(179) = GX(3,i,jp2)*(1.0D0/2.0D0)
      t(180) = t(69)+t(70)+t(179)
      t(183) = t(172)**2
      t(184) = t(174)**2
      Y(1,1,-1,0) = t(90)+e5c*gam12(i,j)*t(2)*t(6)*t(27)*1.5D0-Sig
     &B_d(3,i,3)*t(2)*t(6)*t(14)*t(19)*5.0D-1-SigB_d(3,i,1)*t(2)*t(14)*t
     &(19)*t(21)*5.0D-1
      Y(1,1,-1,1) = t(91)
      Y(1,1,-1,2) = (e5c*gam12(i,j)*t(51)*t(52))/t(25)
      Y(1,1,0,0) = t(92)+t(103)+t(106)-gam11(i,j)*1.0D0/t(25)**5*t
     &(102)*4.0D0-SigB_d(3,i,2)*t(42)*t(45)*3.75D-1-SigB_d(3,i,3)*t(21)*
     &t(29)*t(30)*t(46)*3.75D-1-SigB_d(3,i,3)*t(21)*t(26)*t(29)*t(30)*t(
     &32)*t(38)*2.5D-1-SigB_d(3,i,1)*t(21)*t(26)*t(29)*t(30)*t(32)*t(38)
     &*2.5D-1
      Y(1,1,0,1) = t(113)+t(123)+t(129)+SigB_d(3,i,2)*t(42)*t(45)*
     &5.0D-1+SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(46)*5.0D-1
      Y(1,1,0,2) = t(130)+t(141)+t(143)-SigB_d(3,i,2)*t(42)*t(45)*
     &1.25D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(46)*1.25D-1
      Y(1,1,0,3) = t(148)
      Y(1,1,0,4) = t(149)
      Y(1,1,1,0) = t(150)+t(151)+SigB_d(3,i,3)*t(10)*t(14)*t(19)*t
     &(21)*5.0D-1+SigB_d(3,i,1)*t(10)*t(14)*t(19)*t(21)*5.0D-1
      Y(1,1,1,1) = t(152)
      Y(1,1,1,2) = t(153)
      Y(1,2,-1,0) = SigB_d(3,i,3)*t(2)*t(14)*t(21)*t(58)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(19)*t(21)*t(65)*5.0D-1
      Y(1,2,0,0) = t(82)+t(83)+SigB_d(3,i,3)*t(14)*t(21)*t(58)*t(6
     &8)+SigB_d(3,i,1)*t(19)*t(21)*t(65)*t(68)
      Y(1,2,0,1) = t(86)
      Y(1,2,0,2) = t(89)
      Y(1,2,1,0) = SigB_d(3,i,3)*t(10)*t(14)*t(21)*t(58)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(19)*t(21)*t(65)*5.0D-1
      Y(1,3,-1,0) = SigB_d(3,i,3)*t(2)*t(14)*t(21)*t(72)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(19)*t(21)*t(79)*5.0D-1
      Y(1,3,0,0) = t(154)+t(155)+SigB_d(3,i,3)*t(14)*t(21)*t(72)*t
     &(81)+SigB_d(3,i,1)*t(19)*t(21)*t(79)*t(81)
      Y(1,3,0,1) = t(158)
      Y(1,3,0,2) = t(161)
      Y(1,3,1,0) = SigB_d(3,i,3)*t(10)*t(14)*t(21)*t(72)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(19)*t(21)*t(79)*5.0D-1
      Y(2,1,-1,0) = SigB_d(3,i,3)*t(2)*t(19)*t(21)*t(65)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(14)*t(21)*t(58)*5.0D-1
      Y(2,1,0,0) = t(82)+t(83)+SigB_d(3,i,3)*t(19)*t(21)*t(65)*t(8
     &1)+SigB_d(3,i,1)*t(14)*t(21)*t(58)*t(81)
      Y(2,1,0,1) = t(86)
      Y(2,1,0,2) = t(89)
      Y(2,1,1,0) = SigB_d(3,i,3)*t(10)*t(19)*t(21)*t(65)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(14)*t(21)*t(58)*5.0D-1
      Y(2,2,-1,0) = t(90)+t(170)-SigB_d(3,i,3)*t(2)*t(21)*t(58)*t(
     &65)*5.0D-1-SigB_d(3,i,1)*t(2)*t(21)*t(58)*t(65)*5.0D-1
      Y(2,2,-1,1) = t(91)
      Y(2,2,-1,2) = t(171)
      Y(2,2,0,0) = t(92)+t(103)+t(106)+t(173)-SigB_d(3,i,2)*t(42)*
     &t(126)*3.75D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(127)*3.75D-1-SigB
     &_d(3,i,3)*t(21)*t(29)*t(30)*t(32)*t(104)*t(112)*2.5D-1-SigB_d(3,i,
     &1)*t(21)*t(29)*t(30)*t(32)*t(104)*t(112)*2.5D-1
      Y(2,2,0,1) = t(113)+t(123)+t(129)+SigB_d(3,i,2)*t(42)*t(126)
     &*5.0D-1+SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(127)*5.0D-1
      Y(2,2,0,2) = t(130)+t(141)+t(143)-SigB_d(3,i,2)*t(42)*t(126)
     &*1.25D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(127)*1.25D-1
      Y(2,2,0,3) = t(148)
      Y(2,2,0,4) = t(149)
      Y(2,2,1,0) = t(150)+t(151)+SigB_d(3,i,3)*t(10)*t(21)*t(58)*t
     &(65)*5.0D-1+SigB_d(3,i,1)*t(10)*t(21)*t(58)*t(65)*5.0D-1
      Y(2,2,1,1) = t(152)
      Y(2,2,1,2) = t(153)
      Y(2,3,-1,0) = SigB_d(3,i,3)*t(2)*t(21)*t(65)*t(72)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(21)*t(58)*t(79)*5.0D-1
      Y(2,3,0,0) = t(162)+t(163)+SigB_d(3,i,3)*t(21)*t(65)*t(72)*t
     &(81)+SigB_d(3,i,1)*t(21)*t(58)*t(79)*t(81)
      Y(2,3,0,1) = t(166)
      Y(2,3,0,2) = t(169)
      Y(2,3,1,0) = SigB_d(3,i,3)*t(10)*t(21)*t(65)*t(72)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(21)*t(58)*t(79)*5.0D-1
      Y(3,1,-1,0) = SigB_d(3,i,3)*t(2)*t(19)*t(21)*t(79)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(14)*t(21)*t(72)*5.0D-1
      Y(3,1,0,0) = t(154)+t(155)+SigB_d(3,i,3)*t(19)*t(21)*t(79)*t
     &(81)+SigB_d(3,i,1)*t(14)*t(21)*t(72)*t(81)
      Y(3,1,0,1) = t(158)
      Y(3,1,0,2) = t(161)
      Y(3,1,1,0) = SigB_d(3,i,3)*t(10)*t(19)*t(21)*t(79)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(14)*t(21)*t(72)*5.0D-1
      Y(3,2,-1,0) = SigB_d(3,i,3)*t(2)*t(21)*t(58)*t(79)*(-5.0D-1)
     &-SigB_d(3,i,1)*t(2)*t(21)*t(65)*t(72)*5.0D-1
      Y(3,2,0,0) = t(162)+t(163)+SigB_d(3,i,3)*t(21)*t(58)*t(79)*t
     &(81)+SigB_d(3,i,1)*t(21)*t(65)*t(72)*t(81)
      Y(3,2,0,1) = t(166)
      Y(3,2,0,2) = t(169)
      Y(3,2,1,0) = SigB_d(3,i,3)*t(10)*t(21)*t(58)*t(79)*5.0D-1+Si
     &gB_d(3,i,1)*t(10)*t(21)*t(65)*t(72)*5.0D-1
      Y(3,3,-1,0) = t(90)+t(170)-SigB_d(3,i,3)*t(2)*t(21)*t(72)*t(
     &79)*5.0D-1-SigB_d(3,i,1)*t(2)*t(21)*t(72)*t(79)*5.0D-1
      Y(3,3,-1,1) = t(91)
      Y(3,3,-1,2) = t(171)
      Y(3,3,0,0) = t(92)+t(103)+t(106)+t(173)-SigB_d(3,i,2)*t(42)*
     &t(183)*3.75D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(184)*3.75D-1-SigB
     &_d(3,i,3)*t(21)*t(29)*t(30)*t(32)*t(174)*t(180)*2.5D-1-SigB_d(3,i,
     &1)*t(21)*t(29)*t(30)*t(32)*t(174)*t(180)*2.5D-1
      Y(3,3,0,1) = t(113)+t(123)+t(129)+SigB_d(3,i,2)*t(42)*t(183)
     &*5.0D-1+SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(184)*5.0D-1
      Y(3,3,0,2) = t(130)+t(141)+t(143)-SigB_d(3,i,2)*t(42)*t(183)
     &*1.25D-1-SigB_d(3,i,3)*t(21)*t(29)*t(30)*t(184)*1.25D-1
      Y(3,3,0,3) = t(148)
      Y(3,3,0,4) = t(149)
      Y(3,3,1,0) = t(150)+t(151)+SigB_d(3,i,3)*t(10)*t(21)*t(72)*t
     &(79)*5.0D-1+SigB_d(3,i,1)*t(10)*t(21)*t(72)*t(79)*5.0D-1
      Y(3,3,1,1) = t(152)
      Y(3,3,1,2) = t(153)
      endif
              Y(1:3,1:3,-4:4,-4:4)=-Y(1:3,1:3,-4:4,-4:4)

!! ================== BNDY i,N 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin)   
     &        .and. (BoundaryFibShellin(1,2) .eq. 1)  
     &    .and. (i .ge. ist) .and. (i.le. ien)) then
      if(formulationFlagb .eq. 1) then
      t(2) = ds2(i,jm1)*(-3.0D0)
      t(3) = ds2(i,jm2)+t(2)
      t(4) = 1.0D0/t(3)
      t(5) = 1.0D0/ds1(im1,j)
      t(6) = ds2(i,jm2)*5.0D-1
      t(18) = ds2(i,jm1)*(-1.5D0)
      t(7) = -t(18)
      t(8) = t(6)-t(7)
      t(9) = 1.0D0/t(8)
      t(10) = GX(1,i,j)*(-1.0D0)
      t(11) = t(10)+GX(1,im1,j)
      t(12) = t(5)*t(11)*5.0D-1
      t(13) = 1.0D0/ds1(i,j)
      t(14) = GX(1,ip1,j)*(-1.0D0)
      t(15) = t(14)+GX(1,i,j)
      t(16) = t(13)*t(15)*5.0D-1
      t(17) = t(12)+t(16)
      t(19) = GX(1,i,jm2)*5.0D-1
      t(20) = GX(1,i,jm1)*(-2.0D0)
      t(21) = GX(1,i,j)*1.5D0
      t(22) = t(19)+t(20)+t(21)
      t(23) = 1.0D0/t(3)**3
      t(32) = GX(1,i,jm1)*(-4.0D0)
      t(33) = GX(1,i,j)*3.0D0
      t(24) = t(32)+t(33)+GX(1,i,jm2)
      t(37) = ds1(im1,j)*GX(1,i,j)
      t(38) = ds1(i,j)*GX(1,im1,j)
      t(39) = ds1(im1,j)*GX(1,ip1,j)*(-1.0D0)
      t(40) = ds1(i,j)*GX(1,i,j)*(-1.0D0)
      t(25) = t(37)+t(38)+t(39)+t(40)
      t(26) = 1.0D0/dh1(i,jm2)
      t(27) = 1.0D0/ds1(im1,jm2)
      t(28) = 1.0D0/t(3)**5
      t(29) = ds2(i,jm3)**2
      t(30) = ds2(i,jm2)**2
      t(31) = ds2(i,jm1)**2
      t(34) = t(24)**2
      t(35) = 1.0D0/ds1(im1,j)**2
      t(36) = 1.0D0/ds1(i,j)**2
      t(41) = t(25)**2
      t(42) = 1.0D0/dh1(i,jm1)
      t(43) = 1.0D0/ds1(im1,jm1)
      t(44) = 1.0D0/dh1(i,j)
      t(45) = ds1(i,j)*(-1.0D0)
      t(46) = ds1(im1,j)+t(45)
      t(47) = GX(1,i,jm2)*(1.0D0/2.0D0)
      t(48) = t(20)+t(21)+t(47)
      t(49) = 1.0D0/ds1(i,jm2)
      t(50) = 1.0D0/ds1(i,jm1)
      t(51) = t(6)+t(18)
      t(52) = 1.0D0/t(51)
      t(53) = GX(2,i,jm2)*5.0D-1
      t(54) = GX(2,i,jm1)*(-2.0D0)
      t(55) = GX(2,i,j)*1.5D0
      t(56) = t(53)+t(54)+t(55)
      t(57) = GX(2,i,j)*(-1.0D0)
      t(58) = t(57)+GX(2,im1,j)
      t(59) = t(5)*t(58)*5.0D-1
      t(60) = GX(2,ip1,j)*(-1.0D0)
      t(61) = t(60)+GX(2,i,j)
      t(62) = t(13)*t(61)*5.0D-1
      t(63) = t(59)+t(62)
      t(64) = 1.0D0/t(51)**3
      t(65) = t(5)*5.0D-1
      t(79) = t(13)*(-5.0D-1)
      t(66) = -t(79)
      t(67) = t(65)-t(66)
      t(68) = GX(3,i,jm2)*5.0D-1
      t(69) = GX(3,i,jm1)*(-2.0D0)
      t(70) = GX(3,i,j)*1.5D0
      t(71) = t(68)+t(69)+t(70)
      t(72) = GX(3,i,j)*(-1.0D0)
      t(73) = t(72)+GX(3,im1,j)
      t(74) = t(5)*t(73)*5.0D-1
      t(75) = GX(3,ip1,j)*(-1.0D0)
      t(76) = t(75)+GX(3,i,j)
      t(77) = t(13)*t(76)*5.0D-1
      t(78) = t(74)+t(77)
      t(80) = t(65)+t(79)
      t(81) = SigB_d(4,i,2)*t(22)*t(56)*t(64)*(-5.0D-1)
      t(82) = SigB_d(4,i,3)*t(17)*t(52)*t(63)*(-5.0D-1)
      t(83) = t(81)+t(82)
      t(84) = SigB_d(4,i,2)*t(22)*t(56)*t(64)*2.0D0
      t(85) = SigB_d(4,i,3)*t(17)*t(52)*t(63)*2.0D0
      t(86) = t(84)+t(85)
      t(87) = SigB_d(4,i,2)*t(22)*t(56)*t(64)*(-1.5D0)
      t(88) = SigB_d(4,i,3)*t(17)*t(52)*t(63)*(-1.5D0)
      t(89) = e5c*gam12(i,j)*t(4)*t(26)*t(27)
      t(90) = e5c*gam12(i,j)*t(4)*t(42)*t(43)*(-4.0D0)
      t(91) = sigB(4,i,3)*t(5)*(-5.0D-1)
      t(92) = gam11(i,j)*t(23)*1.2D1
      t(93) = 1.0D0/t(3)**4
      t(94) = ds2(i,jm3)*(-6.0D0)
      t(95) = ds2(i,jm2)*1.1D1
      t(96) = ds2(i,jm1)*9.0D0
      t(97) = t(94)+t(95)+t(96)
      t(98) = gam11(i,j)*t(93)*t(97)*8.0D0
      t(99) = sigB(4,i,2)*t(4)*(-1.0D0)
      t(100) = ds2(i,jm4)*ds2(i,jm2)*3.0D0
      t(101) = ds2(i,jm4)*ds2(i,jm1)*(-9.0D0)
      t(102) = ds2(i,jm3)*ds2(i,jm2)*1.3D1
      t(103) = ds2(i,jm3)*ds2(i,jm1)*1.29D2
      t(104) = ds2(i,jm2)*ds2(i,jm1)*(-2.84D2)
      t(105) = t(29)*(-1.2D1)
      t(106) = t(30)*2.5D1
      t(107) = t(31)*3.9D1
      t(108) = t(100)+t(101)+t(102)+t(103)+t(104)+t(105)+t(106)+t(
     &107)
      t(109) = gam11(i,j)*t(28)*t(108)*(-4.0D0)
      t(116) = GX(2,i,jm1)*(-4.0D0)
      t(117) = GX(2,i,j)*3.0D0
      t(110) = t(116)+t(117)+GX(2,i,jm2)
      t(119) = ds1(im1,j)*GX(2,i,j)
      t(120) = ds1(i,j)*GX(2,im1,j)
      t(121) = ds1(im1,j)*GX(2,ip1,j)*(-1.0D0)
      t(122) = ds1(i,j)*GX(2,i,j)*(-1.0D0)
      t(111) = t(119)+t(120)+t(121)+t(122)
      t(112) = ds1(im1,jm2)+ds1(i,jm2)
      t(113) = e5c*gam12(i,j)*t(4)*t(26)*t(27)*t(49)*t(112)*(-1.0D
     &0)
      t(114) = sigB(4,i,2)*t(4)*4.0D0
      t(164) = ds2(i,jm4)*ds2(i,jm2)*6.0D0
      t(165) = ds2(i,jm4)*ds2(i,jm1)*(-1.8D1)
      t(166) = ds2(i,jm3)*ds2(i,jm2)*9.2D1
      t(167) = ds2(i,jm3)*ds2(i,jm1)*6.0D1
      t(168) = ds2(i,jm2)*ds2(i,jm1)*(-7.6D1)
      t(169) = t(29)*(-2.4D1)
      t(170) = t(30)*(-1.09D2)
      t(171) = t(31)*3.3D1
      t(172) = t(164)+t(165)+t(166)+t(167)+t(168)+t(169)+t(170)+t(
     &171)
      t(115) = gam11(i,j)*t(28)*t(172)*8.0D0
      t(118) = t(110)**2
      t(123) = t(111)**2
      t(124) = ds1(im1,jm1)+ds1(i,jm1)
      t(125) = e5c*gam12(i,j)*t(4)*t(42)*t(43)*t(50)*t(124)*4.0D0
      t(126) = sigB(4,i,2)*t(4)*(-3.0D0)
      t(127) = ds2(i,jm4)*ds2(i,jm2)*9.0D0
      t(128) = ds2(i,jm4)*ds2(i,jm1)*(-2.7D1)
      t(129) = ds2(i,jm3)*ds2(i,jm2)*1.59D2
      t(130) = ds2(i,jm3)*ds2(i,jm1)*2.7D1
      t(131) = ds2(i,jm2)*ds2(i,jm1)*6.6D1
      t(132) = t(29)*(-3.6D1)
      t(133) = t(30)*(-2.18D2)
      t(134) = t(127)+t(128)+t(129)+t(130)+t(131)+t(132)+t(133)
      t(135) = gam11(i,j)*t(28)*t(134)*(-4.0D0)
      t(136) = sigB(4,i,3)*t(5)*t(13)*t(46)*(-5.0D-1)
      t(137) = ds1(im1,j)+ds1(i,j)
      t(138) = e5c*gam12(i,j)*t(4)*t(5)*t(13)*t(44)*t(137)*(-3.0D0
     &)
      t(139) = GX(2,i,jm2)*(1.0D0/2.0D0)
      t(140) = t(54)+t(55)+t(139)
      t(141) = e5c*gam12(i,j)*t(4)*t(26)*t(49)
      t(142) = e5c*gam12(i,j)*t(4)*t(42)*t(50)*(-4.0D0)
      t(143) = sigB(4,i,3)*t(13)*5.0D-1
      t(144) = e5c*gam12(i,j)*t(13)*t(44)*t(52)*1.5D0
      t(145) = SigB_d(4,i,2)*t(22)*t(64)*t(71)*(-5.0D-1)
      t(146) = SigB_d(4,i,3)*t(17)*t(52)*t(78)*(-5.0D-1)
      t(147) = t(145)+t(146)
      t(148) = SigB_d(4,i,2)*t(22)*t(64)*t(71)*2.0D0
      t(149) = SigB_d(4,i,3)*t(17)*t(52)*t(78)*2.0D0
      t(150) = t(148)+t(149)
      t(151) = SigB_d(4,i,2)*t(22)*t(64)*t(71)*(-1.5D0)
      t(152) = SigB_d(4,i,3)*t(17)*t(52)*t(78)*(-1.5D0)
      t(153) = SigB_d(4,i,2)*t(56)*t(64)*t(71)*(-5.0D-1)
      t(154) = SigB_d(4,i,3)*t(52)*t(63)*t(78)*(-5.0D-1)
      t(155) = t(153)+t(154)
      t(156) = SigB_d(4,i,2)*t(56)*t(64)*t(71)*2.0D0
      t(157) = SigB_d(4,i,3)*t(52)*t(63)*t(78)*2.0D0
      t(158) = t(156)+t(157)
      t(159) = SigB_d(4,i,2)*t(56)*t(64)*t(71)*(-1.5D0)
      t(160) = SigB_d(4,i,3)*t(52)*t(63)*t(78)*(-1.5D0)
      t(161) = e5c*gam12(i,j)*t(5)*t(44)*t(52)*1.5D0
      t(173) = GX(3,i,jm1)*(-4.0D0)
      t(174) = GX(3,i,j)*3.0D0
      t(162) = t(173)+t(174)+GX(3,i,jm2)
      t(176) = ds1(im1,j)*GX(3,i,j)
      t(177) = ds1(i,j)*GX(3,im1,j)
      t(178) = ds1(im1,j)*GX(3,ip1,j)*(-1.0D0)
      t(179) = ds1(i,j)*GX(3,i,j)*(-1.0D0)
      t(163) = t(176)+t(177)+t(178)+t(179)
      t(175) = t(162)**2
      t(180) = t(163)**2
      t(181) = GX(3,i,jm2)*(1.0D0/2.0D0)
      t(182) = t(69)+t(70)+t(181)
      Y(1,1,-1,-2) = t(89)
      Y(1,1,-1,-1) = t(90)
      Y(1,1,-1,0) = t(91)+e5c*gam12(i,j)*t(5)*t(9)*t(44)*1.5D0-Sig
     &B_d(4,i,3)*t(5)*t(9)*t(17)*t(22)*5.0D-1-SigB_d(4,i,1)*t(5)*t(17)*t
     &(22)*t(52)*5.0D-1
      Y(1,1,0,-4) = t(92)
      Y(1,1,0,-3) = t(98)
      Y(1,1,0,-2) = t(99)+t(109)+t(113)-SigB_d(4,i,2)*t(23)*t(34)*
     &1.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(41)*2.5D-1
      Y(1,1,0,-1) = t(114)+t(115)+t(125)+SigB_d(4,i,2)*t(23)*t(34)
     &*4.0D0+SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(41)*1.0D0
      Y(1,1,0,0) = t(126)+t(135)+t(136)+t(138)-SigB_d(4,i,2)*t(23)
     &*t(34)*3.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(41)*7.5D-1-SigB_d(4,
     &i,3)*t(4)*t(25)*t(35)*t(36)*t(46)*t(48)*5.0D-1-SigB_d(4,i,1)*t(4)*
     &t(25)*t(35)*t(36)*t(46)*t(48)*5.0D-1
      Y(1,1,1,-2) = t(141)
      Y(1,1,1,-1) = t(142)
      Y(1,1,1,0) = t(143)+t(144)+SigB_d(4,i,3)*t(13)*t(17)*t(22)*t
     &(52)*5.0D-1+SigB_d(4,i,1)*t(13)*t(17)*t(22)*t(52)*5.0D-1
      Y(1,2,-1,0) = SigB_d(4,i,3)*t(5)*t(17)*t(52)*t(56)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(22)*t(52)*t(63)*5.0D-1
      Y(1,2,0,-2) = t(83)
      Y(1,2,0,-1) = t(86)
      Y(1,2,0,0) = t(87)+t(88)+SigB_d(4,i,3)*t(17)*t(52)*t(56)*t(6
     &7)+SigB_d(4,i,1)*t(22)*t(52)*t(63)*t(67)
      Y(1,2,1,0) = SigB_d(4,i,3)*t(13)*t(17)*t(52)*t(56)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(22)*t(52)*t(63)*5.0D-1
      Y(1,3,-1,0) = SigB_d(4,i,3)*t(5)*t(17)*t(52)*t(71)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(22)*t(52)*t(78)*5.0D-1
      Y(1,3,0,-2) = t(147)
      Y(1,3,0,-1) = t(150)
      Y(1,3,0,0) = t(151)+t(152)+SigB_d(4,i,3)*t(17)*t(52)*t(71)*t
     &(80)+SigB_d(4,i,1)*t(22)*t(52)*t(78)*t(80)
      Y(1,3,1,0) = SigB_d(4,i,3)*t(13)*t(17)*t(52)*t(71)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(22)*t(52)*t(78)*5.0D-1
      Y(2,1,-1,0) = SigB_d(4,i,3)*t(5)*t(22)*t(52)*t(63)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(17)*t(52)*t(56)*5.0D-1
      Y(2,1,0,-2) = t(83)
      Y(2,1,0,-1) = t(86)
      Y(2,1,0,0) = t(87)+t(88)+SigB_d(4,i,3)*t(22)*t(52)*t(63)*t(8
     &0)+SigB_d(4,i,1)*t(17)*t(52)*t(56)*t(80)
      Y(2,1,1,0) = SigB_d(4,i,3)*t(13)*t(22)*t(52)*t(63)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(17)*t(52)*t(56)*5.0D-1
      Y(2,2,-1,-2) = t(89)
      Y(2,2,-1,-1) = t(90)
      Y(2,2,-1,0) = t(91)+t(161)-SigB_d(4,i,3)*t(5)*t(52)*t(56)*t(
     &63)*5.0D-1-SigB_d(4,i,1)*t(5)*t(52)*t(56)*t(63)*5.0D-1
      Y(2,2,0,-4) = t(92)
      Y(2,2,0,-3) = t(98)
      Y(2,2,0,-2) = t(99)+t(109)+t(113)-SigB_d(4,i,2)*t(23)*t(118)
     &*1.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(123)*2.5D-1
      Y(2,2,0,-1) = t(114)+t(115)+t(125)+SigB_d(4,i,2)*t(23)*t(118
     &)*4.0D0+SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(123)*1.0D0
      Y(2,2,0,0) = t(126)+t(135)+t(136)+t(138)-SigB_d(4,i,2)*t(23)
     &*t(118)*3.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(123)*7.5D-1-SigB_d(
     &4,i,3)*t(4)*t(35)*t(36)*t(46)*t(111)*t(140)*5.0D-1-SigB_d(4,i,1)*t
     &(4)*t(35)*t(36)*t(46)*t(111)*t(140)*5.0D-1
      Y(2,2,1,-2) = t(141)
      Y(2,2,1,-1) = t(142)
      Y(2,2,1,0) = t(143)+t(144)+SigB_d(4,i,3)*t(13)*t(52)*t(56)*t
     &(63)*5.0D-1+SigB_d(4,i,1)*t(13)*t(52)*t(56)*t(63)*5.0D-1
      Y(2,3,-1,0) = SigB_d(4,i,3)*t(5)*t(52)*t(63)*t(71)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(52)*t(56)*t(78)*5.0D-1
      Y(2,3,0,-2) = t(155)
      Y(2,3,0,-1) = t(158)
      Y(2,3,0,0) = t(159)+t(160)+SigB_d(4,i,3)*t(52)*t(63)*t(71)*t
     &(80)+SigB_d(4,i,1)*t(52)*t(56)*t(78)*t(80)
      Y(2,3,1,0) = SigB_d(4,i,3)*t(13)*t(52)*t(63)*t(71)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(52)*t(56)*t(78)*5.0D-1
      Y(3,1,-1,0) = SigB_d(4,i,3)*t(5)*t(22)*t(52)*t(78)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(17)*t(52)*t(71)*5.0D-1
      Y(3,1,0,-2) = t(147)
      Y(3,1,0,-1) = t(150)
      Y(3,1,0,0) = t(151)+t(152)+SigB_d(4,i,3)*t(22)*t(52)*t(78)*t
     &(80)+SigB_d(4,i,1)*t(17)*t(52)*t(71)*t(80)
      Y(3,1,1,0) = SigB_d(4,i,3)*t(13)*t(22)*t(52)*t(78)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(17)*t(52)*t(71)*5.0D-1
      Y(3,2,-1,0) = SigB_d(4,i,3)*t(5)*t(52)*t(56)*t(78)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(52)*t(63)*t(71)*5.0D-1
      Y(3,2,0,-2) = t(155)
      Y(3,2,0,-1) = t(158)
      Y(3,2,0,0) = t(159)+t(160)+SigB_d(4,i,3)*t(52)*t(56)*t(78)*t
     &(80)+SigB_d(4,i,1)*t(52)*t(63)*t(71)*t(80)
      Y(3,2,1,0) = SigB_d(4,i,3)*t(13)*t(52)*t(56)*t(78)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(52)*t(63)*t(71)*5.0D-1
      Y(3,3,-1,-2) = t(89)
      Y(3,3,-1,-1) = t(90)
      Y(3,3,-1,0) = t(91)+t(161)-SigB_d(4,i,3)*t(5)*t(52)*t(71)*t(
     &78)*5.0D-1-SigB_d(4,i,1)*t(5)*t(52)*t(71)*t(78)*5.0D-1
      Y(3,3,0,-4) = t(92)
      Y(3,3,0,-3) = t(98)
      Y(3,3,0,-2) = t(99)+t(109)+t(113)-SigB_d(4,i,2)*t(23)*t(175)
     &*1.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(180)*2.5D-1
      Y(3,3,0,-1) = t(114)+t(115)+t(125)+SigB_d(4,i,2)*t(23)*t(175
     &)*4.0D0+SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(180)*1.0D0
      Y(3,3,0,0) = t(126)+t(135)+t(136)+t(138)-SigB_d(4,i,2)*t(23)
     &*t(175)*3.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(180)*7.5D-1-SigB_d(
     &4,i,3)*t(4)*t(35)*t(36)*t(46)*t(163)*t(182)*5.0D-1-SigB_d(4,i,1)*t
     &(4)*t(35)*t(36)*t(46)*t(163)*t(182)*5.0D-1
      Y(3,3,1,-2) = t(141)
      Y(3,3,1,-1) = t(142)
      Y(3,3,1,0) = t(143)+t(144)+SigB_d(4,i,3)*t(13)*t(52)*t(71)*t
     &(78)*5.0D-1+SigB_d(4,i,1)*t(13)*t(52)*t(71)*t(78)*5.0D-1
       else
      t(2) = ds2(i,jm1)*(-3.0D0)
      t(3) = ds2(i,jm2)+t(2)
      t(4) = 1.0D0/t(3)
      t(5) = 1.0D0/ds1(im1,j)
      t(6) = ds2(i,jm2)*5.0D-1
      t(18) = ds2(i,jm1)*(-1.5D0)
      t(7) = -t(18)
      t(8) = t(6)-t(7)
      t(9) = 1.0D0/t(8)
      t(10) = GX(1,i,j)*(-1.0D0)
      t(11) = t(10)+GX(1,im1,j)
      t(12) = t(5)*t(11)*5.0D-1
      t(13) = 1.0D0/ds1(i,j)
      t(14) = GX(1,ip1,j)*(-1.0D0)
      t(15) = t(14)+GX(1,i,j)
      t(16) = t(13)*t(15)*5.0D-1
      t(17) = t(12)+t(16)
      t(19) = GX(1,i,jm2)*5.0D-1
      t(20) = GX(1,i,jm1)*(-2.0D0)
      t(21) = GX(1,i,j)*1.5D0
      t(22) = t(19)+t(20)+t(21)
      t(23) = 1.0D0/t(3)**3
      t(32) = GX(1,i,jm1)*(-4.0D0)
      t(33) = GX(1,i,j)*3.0D0
      t(24) = t(32)+t(33)+GX(1,i,jm2)
      t(37) = ds1(im1,j)*GX(1,i,j)
      t(38) = ds1(i,j)*GX(1,im1,j)
      t(39) = ds1(im1,j)*GX(1,ip1,j)*(-1.0D0)
      t(40) = ds1(i,j)*GX(1,i,j)*(-1.0D0)
      t(25) = t(37)+t(38)+t(39)+t(40)
      t(26) = 1.0D0/dh1(i,jm2)
      t(27) = 1.0D0/ds1(im1,jm2)
      t(28) = 1.0D0/t(3)**5
      t(29) = ds2(i,jm3)**2
      t(30) = ds2(i,jm2)**2
      t(31) = ds2(i,jm1)**2
      t(34) = t(24)**2
      t(35) = 1.0D0/ds1(im1,j)**2
      t(36) = 1.0D0/ds1(i,j)**2
      t(41) = t(25)**2
      t(42) = 1.0D0/dh1(i,jm1)
      t(43) = 1.0D0/ds1(im1,jm1)
      t(44) = 1.0D0/dh1(i,j)
      t(45) = ds1(i,j)*(-1.0D0)
      t(46) = ds1(im1,j)+t(45)
      t(47) = GX(1,i,jm2)*(1.0D0/2.0D0)
      t(48) = t(20)+t(21)+t(47)
      t(49) = 1.0D0/ds1(i,jm2)
      t(50) = 1.0D0/ds1(i,jm1)
      t(51) = t(6)+t(18)
      t(52) = 1.0D0/t(51)
      t(53) = GX(2,i,jm2)*5.0D-1
      t(54) = GX(2,i,jm1)*(-2.0D0)
      t(55) = GX(2,i,j)*1.5D0
      t(56) = t(53)+t(54)+t(55)
      t(57) = GX(2,i,j)*(-1.0D0)
      t(58) = t(57)+GX(2,im1,j)
      t(59) = t(5)*t(58)*5.0D-1
      t(60) = GX(2,ip1,j)*(-1.0D0)
      t(61) = t(60)+GX(2,i,j)
      t(62) = t(13)*t(61)*5.0D-1
      t(63) = t(59)+t(62)
      t(64) = 1.0D0/t(51)**3
      t(65) = t(5)*5.0D-1
      t(79) = t(13)*(-5.0D-1)
      t(66) = -t(79)
      t(67) = t(65)-t(66)
      t(68) = GX(3,i,jm2)*5.0D-1
      t(69) = GX(3,i,jm1)*(-2.0D0)
      t(70) = GX(3,i,j)*1.5D0
      t(71) = t(68)+t(69)+t(70)
      t(72) = GX(3,i,j)*(-1.0D0)
      t(73) = t(72)+GX(3,im1,j)
      t(74) = t(5)*t(73)*5.0D-1
      t(75) = GX(3,ip1,j)*(-1.0D0)
      t(76) = t(75)+GX(3,i,j)
      t(77) = t(13)*t(76)*5.0D-1
      t(78) = t(74)+t(77)
      t(80) = t(65)+t(79)
      t(81) = SigB_d(4,i,2)*t(22)*t(56)*t(64)*(-5.0D-1)

      t(82) = SigB_d(4,i,3)*t(17)*t(52)*t(63)*(-5.0D-1)
      t(83) = t(81)+t(82)
      t(84) = SigB_d(4,i,2)*t(22)*t(56)*t(64)*2.0D0
      t(85) = SigB_d(4,i,3)*t(17)*t(52)*t(63)*2.0D0
      t(86) = t(84)+t(85)
      t(87) = SigB_d(4,i,2)*t(22)*t(56)*t(64)*(-1.5D0)
      t(88) = SigB_d(4,i,3)*t(17)*t(52)*t(63)*(-1.5D0)
      t(89) = e5c*gam12(i,j)*t(4)*t(26)*t(27)
      t(90) = e5c*gam12(i,j)*t(4)*t(42)*t(43)*(-4.0D0)
      t(91) = sigB(4,i,3)*t(5)*(-5.0D-1)
      t(92) = gam11(i,j)*t(23)*1.2D1
      t(93) = 1.0D0/t(3)**4
      t(94) = ds2(i,jm3)*(-6.0D0)
      t(95) = ds2(i,jm2)*1.1D1
      t(96) = ds2(i,jm1)*9.0D0
      t(97) = t(94)+t(95)+t(96)
      t(98) = gam11(i,j)*t(93)*t(97)*8.0D0
      t(99) = sigB(4,i,2)*t(4)*(-1.0D0)
      t(100) = ds2(i,jm4)*ds2(i,jm2)*3.0D0
      t(101) = ds2(i,jm4)*ds2(i,jm1)*(-9.0D0)
      t(102) = ds2(i,jm3)*ds2(i,jm2)*1.3D1
      t(103) = ds2(i,jm3)*ds2(i,jm1)*1.29D2
      t(104) = ds2(i,jm2)*ds2(i,jm1)*(-2.84D2)
      t(105) = t(29)*(-1.2D1)
      t(106) = t(30)*2.5D1
      t(107) = t(31)*3.9D1
      t(108) = t(100)+t(101)+t(102)+t(103)+t(104)+t(105)+t(106)+t(
     &107)
      t(109) = gam11(i,j)*t(28)*t(108)*(-4.0D0)
      t(116) = GX(2,i,jm1)*(-4.0D0)
      t(117) = GX(2,i,j)*3.0D0
      t(110) = t(116)+t(117)+GX(2,i,jm2)
      t(119) = ds1(im1,j)*GX(2,i,j)
      t(120) = ds1(i,j)*GX(2,im1,j)
      t(121) = ds1(im1,j)*GX(2,ip1,j)*(-1.0D0)
      t(122) = ds1(i,j)*GX(2,i,j)*(-1.0D0)
      t(111) = t(119)+t(120)+t(121)+t(122)
      t(112) = ds1(im1,jm2)+ds1(i,jm2)
      t(113) = e5c*gam12(i,j)*t(4)*t(26)*t(27)*t(49)*t(112)*(-1.0D
     &0)
      t(114) = sigB(4,i,2)*t(4)*4.0D0
      t(164) = ds2(i,jm4)*ds2(i,jm2)*6.0D0
      t(165) = ds2(i,jm4)*ds2(i,jm1)*(-1.8D1)
      t(166) = ds2(i,jm3)*ds2(i,jm2)*9.2D1
      t(167) = ds2(i,jm3)*ds2(i,jm1)*6.0D1
      t(168) = ds2(i,jm2)*ds2(i,jm1)*(-7.6D1)
      t(169) = t(29)*(-2.4D1)
      t(170) = t(30)*(-1.09D2)
      t(171) = t(31)*3.3D1
      t(172) = t(164)+t(165)+t(166)+t(167)+t(168)+t(169)+t(170)+t(
     &171)
      t(115) = gam11(i,j)*t(28)*t(172)*8.0D0
      t(118) = t(110)**2
      t(123) = t(111)**2
      t(124) = ds1(im1,jm1)+ds1(i,jm1)
      t(125) = e5c*gam12(i,j)*t(4)*t(42)*t(43)*t(50)*t(124)*4.0D0
      t(126) = sigB(4,i,2)*t(4)*(-3.0D0)
      t(127) = ds2(i,jm4)*ds2(i,jm2)*9.0D0
      t(128) = ds2(i,jm4)*ds2(i,jm1)*(-2.7D1)
      t(129) = ds2(i,jm3)*ds2(i,jm2)*1.59D2
      t(130) = ds2(i,jm3)*ds2(i,jm1)*2.7D1
      t(131) = ds2(i,jm2)*ds2(i,jm1)*6.6D1
      t(132) = t(29)*(-3.6D1)
      t(133) = t(30)*(-2.18D2)
      t(134) = t(127)+t(128)+t(129)+t(130)+t(131)+t(132)+t(133)
      t(135) = gam11(i,j)*t(28)*t(134)*(-4.0D0)
      t(136) = sigB(4,i,3)*t(5)*t(13)*t(46)*(-5.0D-1)
      t(137) = ds1(im1,j)+ds1(i,j)
      t(138) = e5c*gam12(i,j)*t(4)*t(5)*t(13)*t(44)*t(137)*(-3.0D0
     &)
      t(139) = GX(2,i,jm2)*(1.0D0/2.0D0)
      t(140) = t(54)+t(55)+t(139)
      t(141) = e5c*gam12(i,j)*t(4)*t(26)*t(49)
      t(142) = e5c*gam12(i,j)*t(4)*t(42)*t(50)*(-4.0D0)
      t(143) = sigB(4,i,3)*t(13)*5.0D-1
      t(144) = e5c*gam12(i,j)*t(13)*t(44)*t(52)*1.5D0
      t(145) = SigB_d(4,i,2)*t(22)*t(64)*t(71)*(-5.0D-1)
      t(146) = SigB_d(4,i,3)*t(17)*t(52)*t(78)*(-5.0D-1)
      t(147) = t(145)+t(146)
      t(148) = SigB_d(4,i,2)*t(22)*t(64)*t(71)*2.0D0
      t(149) = SigB_d(4,i,3)*t(17)*t(52)*t(78)*2.0D0
      t(150) = t(148)+t(149)
      t(151) = SigB_d(4,i,2)*t(22)*t(64)*t(71)*(-1.5D0)
      t(152) = SigB_d(4,i,3)*t(17)*t(52)*t(78)*(-1.5D0)
      t(153) = SigB_d(4,i,2)*t(56)*t(64)*t(71)*(-5.0D-1)
      t(154) = SigB_d(4,i,3)*t(52)*t(63)*t(78)*(-5.0D-1)
      t(155) = t(153)+t(154)
      t(156) = SigB_d(4,i,2)*t(56)*t(64)*t(71)*2.0D0
      t(157) = SigB_d(4,i,3)*t(52)*t(63)*t(78)*2.0D0
      t(158) = t(156)+t(157)
      t(159) = SigB_d(4,i,2)*t(56)*t(64)*t(71)*(-1.5D0)
      t(160) = SigB_d(4,i,3)*t(52)*t(63)*t(78)*(-1.5D0)
      t(161) = e5c*gam12(i,j)*t(5)*t(44)*t(52)*1.5D0
      t(173) = GX(3,i,jm1)*(-4.0D0)
      t(174) = GX(3,i,j)*3.0D0
      t(162) = t(173)+t(174)+GX(3,i,jm2)
      t(176) = ds1(im1,j)*GX(3,i,j)
      t(177) = ds1(i,j)*GX(3,im1,j)
      t(178) = ds1(im1,j)*GX(3,ip1,j)*(-1.0D0)
      t(179) = ds1(i,j)*GX(3,i,j)*(-1.0D0)
      t(163) = t(176)+t(177)+t(178)+t(179)
      t(175) = t(162)**2
      t(180) = t(163)**2
      t(181) = GX(3,i,jm2)*(1.0D0/2.0D0)
      t(182) = t(69)+t(70)+t(181)
      Y(1,1,-1,-2) = t(89)
      Y(1,1,-1,-1) = t(90)
      Y(1,1,-1,0) = t(91)+e5c*gam12(i,j)*t(5)*t(9)*t(44)*1.5D0-Sig
     &B_d(4,i,3)*t(5)*t(9)*t(17)*t(22)*5.0D-1-SigB_d(4,i,1)*t(5)*t(17)*t
     &(22)*t(52)*5.0D-1
      Y(1,1,0,-4) = t(92)
      Y(1,1,0,-3) = t(98)
      Y(1,1,0,-2) = t(99)+t(109)+t(113)-SigB_d(4,i,2)*t(23)*t(34)*
     &1.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(41)*2.5D-1
      Y(1,1,0,-1) = t(114)+t(115)+t(125)+SigB_d(4,i,2)*t(23)*t(34)
     &*4.0D0+SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(41)*1.0D0
      Y(1,1,0,0) = t(126)+t(135)+t(136)+t(138)-SigB_d(4,i,2)*t(23)
     &*t(34)*3.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(41)*7.5D-1-SigB_d(4,
     &i,3)*t(4)*t(25)*t(35)*t(36)*t(46)*t(48)*5.0D-1-SigB_d(4,i,1)*t(4)*
     &t(25)*t(35)*t(36)*t(46)*t(48)*5.0D-1
      Y(1,1,1,-2) = t(141)
      Y(1,1,1,-1) = t(142)
      Y(1,1,1,0) = t(143)+t(144)+SigB_d(4,i,3)*t(13)*t(17)*t(22)*t
     &(52)*5.0D-1+SigB_d(4,i,1)*t(13)*t(17)*t(22)*t(52)*5.0D-1
      Y(1,2,-1,0) = SigB_d(4,i,3)*t(5)*t(17)*t(52)*t(56)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(22)*t(52)*t(63)*5.0D-1
      Y(1,2,0,-2) = t(83)
      Y(1,2,0,-1) = t(86)
      Y(1,2,0,0) = t(87)+t(88)+SigB_d(4,i,3)*t(17)*t(52)*t(56)*t(6
     &7)+SigB_d(4,i,1)*t(22)*t(52)*t(63)*t(67)
      Y(1,2,1,0) = SigB_d(4,i,3)*t(13)*t(17)*t(52)*t(56)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(22)*t(52)*t(63)*5.0D-1
      Y(1,3,-1,0) = SigB_d(4,i,3)*t(5)*t(17)*t(52)*t(71)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(22)*t(52)*t(78)*5.0D-1
      Y(1,3,0,-2) = t(147)
      Y(1,3,0,-1) = t(150)
      Y(1,3,0,0) = t(151)+t(152)+SigB_d(4,i,3)*t(17)*t(52)*t(71)*t
     &(80)+SigB_d(4,i,1)*t(22)*t(52)*t(78)*t(80)
      Y(1,3,1,0) = SigB_d(4,i,3)*t(13)*t(17)*t(52)*t(71)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(22)*t(52)*t(78)*5.0D-1
      Y(2,1,-1,0) = SigB_d(4,i,3)*t(5)*t(22)*t(52)*t(63)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(17)*t(52)*t(56)*5.0D-1
      Y(2,1,0,-2) = t(83)
      Y(2,1,0,-1) = t(86)
      Y(2,1,0,0) = t(87)+t(88)+SigB_d(4,i,3)*t(22)*t(52)*t(63)*t(8
     &0)+SigB_d(4,i,1)*t(17)*t(52)*t(56)*t(80)
      Y(2,1,1,0) = SigB_d(4,i,3)*t(13)*t(22)*t(52)*t(63)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(17)*t(52)*t(56)*5.0D-1
      Y(2,2,-1,-2) = t(89)
      Y(2,2,-1,-1) = t(90)
      Y(2,2,-1,0) = t(91)+t(161)-SigB_d(4,i,3)*t(5)*t(52)*t(56)*t(
     &63)*5.0D-1-SigB_d(4,i,1)*t(5)*t(52)*t(56)*t(63)*5.0D-1
      Y(2,2,0,-4) = t(92)
      Y(2,2,0,-3) = t(98)
      Y(2,2,0,-2) = t(99)+t(109)+t(113)-SigB_d(4,i,2)*t(23)*t(118)
     &*1.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(123)*2.5D-1
      Y(2,2,0,-1) = t(114)+t(115)+t(125)+SigB_d(4,i,2)*t(23)*t(118
     &)*4.0D0+SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(123)*1.0D0
      Y(2,2,0,0) = t(126)+t(135)+t(136)+t(138)-SigB_d(4,i,2)*t(23)
     &*t(118)*3.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(123)*7.5D-1-SigB_d(
     &4,i,3)*t(4)*t(35)*t(36)*t(46)*t(111)*t(140)*5.0D-1-SigB_d(4,i,1)*t
     &(4)*t(35)*t(36)*t(46)*t(111)*t(140)*5.0D-1
      Y(2,2,1,-2) = t(141)
      Y(2,2,1,-1) = t(142)
      Y(2,2,1,0) = t(143)+t(144)+SigB_d(4,i,3)*t(13)*t(52)*t(56)*t
     &(63)*5.0D-1+SigB_d(4,i,1)*t(13)*t(52)*t(56)*t(63)*5.0D-1
      Y(2,3,-1,0) = SigB_d(4,i,3)*t(5)*t(52)*t(63)*t(71)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(52)*t(56)*t(78)*5.0D-1
      Y(2,3,0,-2) = t(155)
      Y(2,3,0,-1) = t(158)
      Y(2,3,0,0) = t(159)+t(160)+SigB_d(4,i,3)*t(52)*t(63)*t(71)*t
     &(80)+SigB_d(4,i,1)*t(52)*t(56)*t(78)*t(80)
      Y(2,3,1,0) = SigB_d(4,i,3)*t(13)*t(52)*t(63)*t(71)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(52)*t(56)*t(78)*5.0D-1
      Y(3,1,-1,0) = SigB_d(4,i,3)*t(5)*t(22)*t(52)*t(78)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(17)*t(52)*t(71)*5.0D-1
      Y(3,1,0,-2) = t(147)
      Y(3,1,0,-1) = t(150)
      Y(3,1,0,0) = t(151)+t(152)+SigB_d(4,i,3)*t(22)*t(52)*t(78)*t
     &(80)+SigB_d(4,i,1)*t(17)*t(52)*t(71)*t(80)
      Y(3,1,1,0) = SigB_d(4,i,3)*t(13)*t(22)*t(52)*t(78)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(17)*t(52)*t(71)*5.0D-1
      Y(3,2,-1,0) = SigB_d(4,i,3)*t(5)*t(52)*t(56)*t(78)*(-5.0D-1)
     &-SigB_d(4,i,1)*t(5)*t(52)*t(63)*t(71)*5.0D-1
      Y(3,2,0,-2) = t(155)
      Y(3,2,0,-1) = t(158)
      Y(3,2,0,0) = t(159)+t(160)+SigB_d(4,i,3)*t(52)*t(56)*t(78)*t
     &(80)+SigB_d(4,i,1)*t(52)*t(63)*t(71)*t(80)
      Y(3,2,1,0) = SigB_d(4,i,3)*t(13)*t(52)*t(56)*t(78)*5.0D-1+Si
     &gB_d(4,i,1)*t(13)*t(52)*t(63)*t(71)*5.0D-1
      Y(3,3,-1,-2) = t(89)
      Y(3,3,-1,-1) = t(90)
      Y(3,3,-1,0) = t(91)+t(161)-SigB_d(4,i,3)*t(5)*t(52)*t(71)*t(
     &78)*5.0D-1-SigB_d(4,i,1)*t(5)*t(52)*t(71)*t(78)*5.0D-1
      Y(3,3,0,-4) = t(92)
      Y(3,3,0,-3) = t(98)
      Y(3,3,0,-2) = t(99)+t(109)+t(113)-SigB_d(4,i,2)*t(23)*t(175)
     &*1.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(180)*2.5D-1
      Y(3,3,0,-1) = t(114)+t(115)+t(125)+SigB_d(4,i,2)*t(23)*t(175
     &)*4.0D0+SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(180)*1.0D0
      Y(3,3,0,0) = t(126)+t(135)+t(136)+t(138)-SigB_d(4,i,2)*t(23)
     &*t(175)*3.0D0-SigB_d(4,i,3)*t(4)*t(35)*t(36)*t(180)*7.5D-1-SigB_d(
     &4,i,3)*t(4)*t(35)*t(36)*t(46)*t(163)*t(182)*5.0D-1-SigB_d(4,i,1)*t
     &(4)*t(35)*t(36)*t(46)*t(163)*t(182)*5.0D-1
      Y(3,3,1,-2) = t(141)
      Y(3,3,1,-1) = t(142)
      Y(3,3,1,0) = t(143)+t(144)+SigB_d(4,i,3)*t(13)*t(52)*t(71)*t
     &(78)*5.0D-1+SigB_d(4,i,1)*t(13)*t(52)*t(71)*t(78)*5.0D-1
       endif
              Y(1:3,1:3,-4:4,-4:4)=-Y(1:3,1:3,-4:4,-4:4)
       else
!! ================== BNDY 2,2
!! ================== BNDY Nx-1,2
!! ================== BNDY 2,Ny-1
!! ================== BNDY Nx-1,Ny-1

!! ================== BNDY 2,j  2<=j<=Ny-1
!! ================== BNDY N-1,j  2<=j<=Ny-1
!! ================== BNDY i,2 2<=i<=Nx-1
!! ================== BNDY i,N-1 2<=i<=Nx-1
!! ================== Middle Points
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &            *(GX(i2,ihm1,j) - GX(i2,i,j))/(dh1(i,j)*ds1(ihm1,j))   
              enddo;enddo 

              DY_i_jp(1:3)=
     &         (GX(1:3,ihm1,jhp1)-GX(1:3,ihm1,j))/ds2(ihm1,j)

              DY_i_jm(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm1,jhm1))/ds2(ihm1,jhm1)                      

              DY_ip_j(1:3)=
     &        (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j)                    

              DY_ip_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1)                

              DY_im_j(1:3)=0.25*(
     &         DY_i_jp(1:3)+DY_i_jm(1:3)+DY_ip_j(1:3)+DY_ip_jm(1:3) )  

              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,-1)
     &                =Y(i2,j2,-1:0,-1)
     &                +Siga_d(ihm1,j,2)*DY_im_j(j2)*0.25
     &                *(/-1.0/ds2(ihm1,jhm1),-1.0/ds2(i,jhm1)/)
     &            *(GX(i2,ihm1,j) - GX(i2,i,j))/(dh1(i,j)*ds1(ihm1,j))   
              enddo;enddo 

              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,2)*DY_im_j(j2)*0.25
     &                *(/-1.0/ds2(ihm1,j)+1.0/ds2(ihm1,jhm1)
     &                  ,-1.0/ds2(i,   j)+1.0/ds2(i,   jhm1)/)
     &            *(GX(i2,ihm1,j) - GX(i2,i,j))/(dh1(i,j)*ds1(ihm1,j))   
              enddo;enddo 

              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,1)
     &                =Y(i2,j2,-1:0,1)
     &                +Siga_d(ihm1,j,2)*DY_im_j(j2)*0.25
     &                *(/1.0/ds2(ihm1,j),1.0/ds2(i,j)/)
     &            *(GX(i2,ihm1,j) - GX(i2,i,j))/(dh1(i,j)*ds1(ihm1,j))   
              enddo;enddo 



              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,0:1,0)
     &                =Y(i2,j2,0:1,0)
     &                +Siga_d(i,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(i,j),1.0/ds1(i,j)/)
     &            *-(GX(i2,i,j) - GX(i2,ihp1,j))/(dh1(i,j)*ds1(i,j))   
              enddo;enddo

              DY_i_jp(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j)

              DY_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1)                      

              DY_ip_j(1:3)=
     &        (GX(1:3,ihp1,jhp1)-GX(1:3,ihp1,j))/ds2(ihp1,j)                    

              DY_ip_jm(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,ihp1,jhm1))/ds2(ihp1,jhm1)                

              DY_im_j(1:3)=0.25*(
     &         DY_i_jp(1:3)+DY_i_jm(1:3)+DY_ip_j(1:3)+DY_ip_jm(1:3) )  

              do i2=1,3;do j2=1,3
                 Y(i2,j2,0:1,-1)
     &                =Y(i2,j2,0:1,-1)
     &                +Siga_d(i,j,2)*DY_im_j(j2)*0.25
     &                *(/-1.0/ds2(i,jhm1),-1.0/ds2(ihp1,jhm1)/)
     &            *-(GX(i2,i,j) - GX(i2,ihp1,j))/(dh1(i,j)*ds1(i,j))   
              enddo;enddo 

              do i2=1,3;do j2=1,3
                 Y(i2,j2,0:1,0)
     &                =Y(i2,j2,0:1,0)
     &                +Siga_d(i,j,2)*DY_im_j(j2)*0.25
     &                *(/-1.0/ds2(i,j)+1.0/ds2(i,jhm1)
     &                  ,-1.0/ds2(ihp1,j)+1.0/ds2(ihp1,jhm1)/)
     &            *-(GX(i2,i,j) - GX(i2,ihp1,j))/(dh1(i,j)*ds1(i,j))    
              enddo;enddo 

              do i2=1,3;do j2=1,3
                 Y(i2,j2,0:1,1)
     &                =Y(i2,j2,0:1,1)
     &                +Siga_d(i,j,2)*DY_im_j(j2)*0.25
     &                *(/1.0/ds2(i,j),1.0/ds2(ihp1,j)/)
     &            *-(GX(i2,i,j) - GX(i2,ihp1,j))/(dh1(i,j)*ds1(i,j))   
              enddo;enddo 

!sig12 -1,-1

              DY_ip_j(1:3)=
     &        (GX(1:3,i,j)-GX(1:3,i,jhm1))
     &                     /ds2(i,jhm1) 
              DY_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm1,jhm1))/ds2(ihm1,jhm1) 
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))

              DX_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))
     &                     /ds1(ihm1,j) 
              DX_i_jm(1:3)=
     &         (GX(1:3,i,jhm1)-GX(1:3,ihm1,jhm1))
     &                     /ds1(ihm1,jhm1) 
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(ihm1,jhm1,5)*(
     &                ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))
     &               /(4.0*ds2(im1,jhm1)) + (GX(i2,i,jhm1) - GX(i2,i,j))
     &               /(4.0*ds2(i,jhm1)))/dh1(i,j)
     &               +((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))
     &               /(4.0*ds1(ihm1,jm1)) + (GX(i2,ihm1,j) - GX(i2,i,j))
     &               /(4.0*ds1(ihm1,j)))/dh2(i,j)
     &                )
              Y(i2,j2,-1:0,0)=
     &                Y(i2,j2,-1:0,0)
     &                +DY_i_j(j2)
     &                *0.5*(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &                *temp_ibm3 

              Y(i2,j2,-1:0,-1)=
     &                Y(i2,j2,-1:0,-1)
     &                +DY_i_j(j2)
     &                *0.5*(/-1.0/ds1(ihm1,jhm1),1.0/ds1(ihm1,jhm1)/)
     &                *temp_ibm3 

              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +DX_i_j(j2)
     &                *0.5*(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
     &                *temp_ibm3 

              Y(i2,j2,-1,-1:0)=
     &                Y(i2,j2,-1,-1:0)
     &                +DX_i_j(j2)
     &                *0.5*(/-1.0/ds2(ihm1,jhm1),1.0/ds2(ihm1,jhm1)/)
     &                *temp_ibm3 
              enddo;enddo 

!sig12 -1,0

              DY_ip_j(1:3)=
     &        (GX(1:3,i,jhp1)-GX(1:3,i,j))
     &                     /ds2(i,j) 
              DY_im_j(1:3)=
     &         (GX(1:3,ihm1,jhp1)-GX(1:3,ihm1,j))/ds2(ihm1,j) 
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))

              DX_i_jp(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,ihm1,jhp1))
     &                     /ds1(ihm1,jhp1) 
              DX_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))
     &                     /ds1(ihm1,j) 
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(ihm1,j,5)*(
     &                ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))
     &               /(4.0*ds2(im1,j)) 
     &               + (GX(i2,i,j) - GX(i2,i,jhp1))
     &               /(4.0*ds2(i,j)))/dh1(i,j)
     &               -((GX(i2,ihm1,j) - GX(i2,i,j))
     &               /(4.0*ds1(ihm1,j)) 
     &               + (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))
     &               /(4.0*ds1(ihm1,jp1)))/dh2(i,j)
     &                )
              Y(i2,j2,-1:0,1)=
     &                Y(i2,j2,-1:0,1)
     &                +DY_i_j(j2)
     &                *0.5*(/-1.0/ds1(ihm1,jhp1),1.0/ds1(ihm1,jhp1)/)
     &                *temp_ibm3 

              Y(i2,j2,-1:0,0)=
     &                Y(i2,j2,-1:0,0)
     &                +DY_i_j(j2)

     &                *0.5*(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &                *temp_ibm3 

              Y(i2,j2,0,0:1)=
     &                Y(i2,j2,0,0:1)
     &                +DX_i_j(j2)
     &                *0.5*(/-1.0/ds2(i,j),1.0/ds2(i,j)/)
     &                *temp_ibm3 

              Y(i2,j2,-1,0:1)=
     &                Y(i2,j2,-1,0:1)
     &                +DX_i_j(j2)
     &                *0.5*(/-1.0/ds2(ihm1,j),1.0/ds2(ihm1,j)/)
     &                *temp_ibm3 
              enddo;enddo 

!sig12 0,-1

              DY_ip_j(1:3)=
     &        (GX(1:3,ihp1,j)-GX(1:3,ihp1,jhm1))
     &                     /ds2(ihp1,jhm1) 
              DY_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))

              DX_i_jp(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))
     &                     /ds1(i,j) 
              DX_i_jm(1:3)=
     &         (GX(1:3,ihp1,jhm1)-GX(1:3,i,jhm1))
     &                     /ds1(i,jhm1) 
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(i,jhm1,5)*(
     &                -((GX(i2,i,jhm1) - GX(i2,i,j))
     &               /(4.0*ds2(i,jhm1)) 
     &               + (GX(i2,ihp1,jhm1)- GX(i2,ihp1,j))
     &               /(4.0*ds2(ip1,jhm1)))/dh1(i,j)
     &               +((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))
     &               /(4.0*ds1(i,jm1)) 
     &               + (GX(i2,i,j) - GX(i2,ihp1,j))
     &               /(4.0*ds1(i,j)))/dh2(i,j)
     &                )
              Y(i2,j2,0:1,0)=
     &                Y(i2,j2,0:1,0)
     &                +DY_i_j(j2)
     &                *0.5*(/-1.0/ds1(i,j),1.0/ds1(i,j)/)
     &                *temp_ibm3 

              Y(i2,j2,0:1,-1)=
     &                Y(i2,j2,0:1,-1)
     &                +DY_i_j(j2)
     &                *0.5*(/-1.0/ds1(i,jhm1),1.0/ds1(i,jhm1)/)
     &                *temp_ibm3 

              Y(i2,j2,1,-1:0)=
     &                Y(i2,j2,1,-1:0)
     &                +DX_i_j(j2)
     &                *0.5*(/-1.0/ds2(ihp1,jhm1),1.0/ds2(ihp1,jhm1)/)
     &                *temp_ibm3 

              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +DX_i_j(j2)
     &                *0.5*(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
     &                *temp_ibm3 
              enddo;enddo  
!sig12 0,0

              DY_ip_j(1:3)=
     &        (GX(1:3,ihp1,jhp1)-GX(1:3,ihp1,j))
     &                     /ds2(ihp1,j) 
              DY_im_j(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))

              DX_i_jp(1:3)=
     &         (GX(1:3,ihp1,jhp1)-GX(1:3,i,jhp1))
     &                     /ds1(i,jhp1) 
              DX_i_jm(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))
     &                     /ds1(i,j) 
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(i,j,5)*(
     &               -((GX(i2,i,j) - GX(i2,i,jhp1))
     &               /(4.0*ds2(i,j)) 
     &               + (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))
     &               /(4.0*ds2(ip1,j)))/dh1(i,j)
     &               -((GX(i2,i,j) - GX(i2,ihp1,j))
     &               /(4.0*ds1(i,j)) 
     &               + (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))
     &               /(4.0*ds1(i,jp1)))/dh2(i,j)
     &                )
              Y(i2,j2,0:1,1)=
     &                Y(i2,j2,0:1,1)
     &                +DY_i_j(j2)
     &                *0.5*(/-1.0/ds1(i,jhp1),1.0/ds1(i,jhp1)/)
     &                *temp_ibm3 

              Y(i2,j2,0:1,0)=
     &                Y(i2,j2,0:1,0)
     &                +DY_i_j(j2)
     &                *0.5*(/-1.0/ds1(i,j),1.0/ds1(i,j)/)
     &                *temp_ibm3 

              Y(i2,j2,1,0:1)=
     &                Y(i2,j2,1,0:1)
     &                +DX_i_j(j2)
     &                *0.5*(/-1.0/ds2(ihp1,j),1.0/ds2(ihp1,j)/)
     &                *temp_ibm3 

              Y(i2,j2,0,0:1)=
     &                Y(i2,j2,0,0:1)
     &                +DX_i_j(j2)
     &                *0.5*(/-1.0/ds2(i,j),1.0/ds2(i,j)/)
     &                *temp_ibm3 
              enddo;enddo  
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2


              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2

              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
     &             *(GX(i2,i,jhm1) - GX(i2,i,j))/(dh2(i,j)*ds2(i,jhm1))  
              enddo;enddo 

              DY_i_jp(1:3)=
     &         (GX(1:3,ihp1,jhm1)-GX(1:3,i,jhm1))/ds1(i,jhm1)

              DY_i_jm(1:3)=
     &         (GX(1:3,i,jhm1)-GX(1:3,ihm1,jhm1))/ds1(ihm1,jhm1)                      

              DY_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)

              DY_ip_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)                      

              DY_im_j(1:3)=0.25*(
     &         DY_i_jp(1:3)+DY_i_jm(1:3)+DY_ip_j(1:3)+DY_ip_jm(1:3) ) 

              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1,-1:0)
     &                =Y(i2,j2,-1,-1:0)
     &                +Siga_d(ihm1,j,3)*DY_im_j(j2)*0.25
     &                *(/-1.0/ds1(ihm1,jhm1),-1.0/ds1(ihm1,j)/)
     &             *(GX(i2,i,jhm1) - GX(i2,i,j))/(dh2(i,j)*ds2(i,jhm1))    
              enddo;enddo 

              do i2=1,3;do j2=1,3
                 Y(i2,j2,0,-1:0)
     &                =Y(i2,j2,0,-1:0)
     &                +Siga_d(ihm1,j,3)*DY_im_j(j2)*0.25
     &               *(/-1.0/ds1(i,jhm1)+1.0/ds1(ihm1,jhm1)
     &                 ,-1.0/ds1(i,   j)+1.0/ds1(ihm1,   j)/)
     &             *(GX(i2,i,jhm1) - GX(i2,i,j))/(dh2(i,j)*ds2(i,jhm1))    
              enddo;enddo 

              do i2=1,3;do j2=1,3
                 Y(i2,j2,1,-1:0)
     &                =Y(i2,j2,1,-1:0)
     &                +Siga_d(ihm1,j,3)*DY_im_j(j2)*0.25
     &                *(/1.0/ds1(i,jhm1),1.0/ds1(i,j)/)
     &             *(GX(i2,i,jhm1) - GX(i2,i,j))/(dh2(i,j)*ds2(i,jhm1))    
              enddo;enddo 


              DY_i_jp(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2

              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)=
     &                Y(i2,j2,0,0:1)
     &                +Siga_d(i,j,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,j ),1.0/ds2(i,j )/)
     &             *-(GX(i2,i,j) - GX(i2,i,jhp1))/(dh2(i,j)*ds2(i,j))  
              enddo;enddo 

              DY_i_jp(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)

              DY_i_jm(1:3)=

     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)                      

              DY_ip_j(1:3)=
     &         (GX(1:3,ihp1,jhp1)-GX(1:3,i,jhp1))/ds1(i,jhp1)

              DY_ip_jm(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,ihm1,jhp1))/ds1(ihm1,jhp1)                      

              DY_im_j(1:3)=0.25*(
     &         DY_i_jp(1:3)+DY_i_jm(1:3)+DY_ip_j(1:3)+DY_ip_jm(1:3) ) 

              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1,0:1)
     &                =Y(i2,j2,-1,0:1)
     &                +Siga_d(ihm1,j,3)*DY_im_j(j2)*0.25
     &                *(/-1.0/ds1(ihm1,j),-1.0/ds1(ihm1,jhp1)/)
     &             *-(GX(i2,i,j) - GX(i2,i,jhp1))/(dh2(i,j)*ds2(i,j))  
              enddo;enddo 
 
              do i2=1,3;do j2=1,3
                 Y(i2,j2,0,0:1)
     &                =Y(i2,j2,0,0:1)
     &                +Siga_d(ihm1,j,3)*DY_im_j(j2)*0.25
     &               *(/-1.0/ds1(i,j)+1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(i,jhp1)+1.0/ds1(ihm1,jhp1)/)
     &             *-(GX(i2,i,j) - GX(i2,i,jhp1))/(dh2(i,j)*ds2(i,j))  
              enddo;enddo 

              do i2=1,3;do j2=1,3
                 Y(i2,j2,1,0:1)
     &                =Y(i2,j2,1,0:1)
     &                +Siga_d(ihm1,j,3)*DY_im_j(j2)*0.25
     &                *(/1.0/ds1(i,j),1.0/ds1(i,jhp1)/)
     &             *-(GX(i2,i,j) - GX(i2,i,jhp1))/(dh2(i,j)*ds2(i,j))  
              enddo;enddo 
              Y(1:3,1:3,-4:4,-4:4)=-Y(1:3,1:3,-4:4,-4:4)
              do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
              enddo
              
       endif


      if(MaterialTypeIBM .eq. 1) then
        write(*,*) 'not yet implemneted'
      endif

!===TEST=================================== 
!1313    continue
!              do i3=1,3
!                 Y(i3,i3,-4:4,-4:4)=
!     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
!              enddo
!======================================      
        do i3=1,3
         do j3=1,3
            Y(i3,j3,0,0)=Y(i3,j3,0,0)+(1.0-coefX00)*
     &           aamat_intmp(i3,j3)/alpha_coefK
         enddo
        enddo
      
        do i3=1,3
         noj=noj+1
         iacsrfibershell(noj)=ntemi+1
         do j3=1,3
          do i2=-4,4
           do j2=-4,4
            if(abs(Y(i3,j3,i2,j2)).gt. 1.d-16) then
              ntemi=ntemi+1
              jacsrfibershell(ntemi)=iAdd(i2,j2,j3)
              aacsrfibershell(ntemi)=alpha_coefK*Y(i3,j3,i2,j2)
            endif
           enddo
          enddo
          enddo !j3
         enddo !i3
      enddo !j
      enddo !i
      noj=noj+1
      iacsrfibershell(noj)=ntemi+1 
      return 
      end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine KnonlinearBNDY1j(A0,iA0,P,Q,ds1,ds2,dh1,dh2,
     &   dBsig21,dsig12,dBsig22A,dBsig22B,dsig11A,dsig11B)
      implicit none
      integer iA0(9,2)
      real*8 t(150),A0(9)
      real*8 ds1(0:1,-1:1),ds2(0:1,-1:1),dh2(0:1,-1:1),dh1(0:1,-1:1)
     &  ,P(0:2,-1:1),Q(0:2,-1:1)
     &  ,dBsig21(-1:1),dsig12(-1:1)
     &  ,dBsig22A(-1:1),dBsig22B(-1:1)
     &  ,dsig11A,dsig11B
      t(2) = 1.0D0/ds2(0,-1)
      t(3) = ds1(0,-1)*1.5D0
      t(7) = ds1(+1,-1)*(-5.0D-1)
      t(4) = -t(7)
      t(5) = t(3)-t(4)
      t(6) = 1.0D0/t(5)
      t(8) = ds1(0,0)*1.5D0
      t(31) = ds1(+1,0)*(-5.0D-1)
      t(9) = -t(31)
      t(10) = t(8)-t(9)
      t(11) = 1.0D0/t(10)
      t(12) = 1.0D0/dh1(0,0)
      t(13) = P(0,0)*(-1.0D0)
      t(14) = P(0,-1)+t(13)
      t(15) = P(+1,0)*(-1.0D0)
      t(16) = 1.0D0/ds1(0,0)
      t(17) = 1.0D0/ds1(0,-1)
      t(18) = t(2)*t(14)*2.5D-1
      t(19) = 1.0D0/ds2(+1,-1)
      t(20) = P(+1,-1)+t(15)
      t(21) = t(19)*t(20)*2.5D-1
      t(22) = t(18)+t(21)
      t(23) = Q(+1,0)*(-1.0D0)
      t(24) = 1.0D0/dh2(0,-1)
      t(25) = t(3)+t(7)
      t(26) = 1.0D0/t(25)
      t(27) = P(0,-1)*1.5D0
      t(28) = P(+1,-1)*(-2.0D0)
      t(29) = P(+2,-1)*5.0D-1
      t(30) = t(27)+t(28)+t(29)
      t(32) = P(0,0)*1.5D0
      t(33) = P(+1,0)*(-2.0D0)
      t(34) = P(+2,0)*5.0D-1
      t(35) = t(32)+t(33)+t(34)
      t(36) = Q(0,0)*(-1.0D0)
      t(37) = Q(0,-1)+t(36)
      t(38) = P(0,0)+t(15)
      t(39) = 1.0D0/dh2(0,0)
      t(40) = t(8)+t(31)
      t(41) = 1.0D0/t(40)
      t(42) = 1.0D0/ds2(0,0)
      t(43) = t(26)*t(30)*5.0D-1
      t(44) = t(35)*t(41)*5.0D-1
      t(45) = t(43)+t(44)
      t(46) = Q(0,-1)*1.5D0
      t(47) = Q(+1,-1)*(-2.0D0)
      t(48) = Q(+2,-1)*5.0D-1
      t(49) = t(46)+t(47)+t(48)
      t(50) = t(26)*t(49)*5.0D-1
      t(51) = Q(0,0)*1.5D0
      t(52) = Q(+1,0)*(-2.0D0)
      t(53) = Q(+2,0)*5.0D-1
      t(54) = t(51)+t(52)+t(53)
      t(55) = t(41)*t(54)*5.0D-1
      t(56) = ds1(0,+1)*1.5D0
      t(79) = ds1(+1,+1)*(-5.0D-1)
      t(57) = -t(79)
      t(58) = t(56)-t(57)
      t(59) = 1.0D0/t(58)
      t(60) = Q(+1,-1)*(-1.0D0)
      t(61) = Q(0,-1)+t(60)
      t(62) = t(17)*t(61)
      t(63) = Q(0,0)+t(23)
      t(64) = t(16)*t(63)
      t(65) = t(62)+t(64)
      t(66) = P(0,+1)*(-1.0D0)
      t(67) = P(0,0)+t(66)
      t(68) = t(2)*t(37)
      t(69) = Q(+1,-1)+t(23)
      t(70) = t(19)*t(69)
      t(71) = t(68)+t(70)
      t(72) = t(42)*t(67)*2.5D-1
      t(73) = 1.0D0/ds2(+1,0)
      t(74) = P(+1,+1)*(-1.0D0)
      t(75) = P(+1,0)+t(74)
      t(76) = t(73)*t(75)*2.5D-1
      t(77) = t(72)+t(76)
      t(78) = Q(+1,+1)*(-1.0D0)
      t(80) = P(0,+1)*1.5D0
      t(81) = P(+1,+1)*(-2.0D0)
      t(82) = P(+2,+1)*5.0D-1
      t(83) = t(80)+t(81)+t(82)
      t(84) = Q(0,+1)*(-1.0D0)
      t(85) = Q(0,0)+t(84)
      t(86) = t(56)+t(79)
      t(87) = 1.0D0/t(86)
      t(88) = t(83)*t(87)*5.0D-1
      t(89) = t(44)+t(88)
      t(90) = Q(0,+1)*1.5D0
      t(91) = Q(+1,+1)*(-2.0D0)
      t(92) = Q(+2,+1)*5.0D-1
      t(93) = t(90)+t(91)+t(92)
      t(94) = 1.0D0/ds1(0,+1)
      t(95) = Q(0,+1)+t(78)
      t(96) = t(94)*t(95)
      t(97) = t(64)+t(96)
      t(98) = t(42)*t(85)
      t(99) = Q(+1,0)+t(78)
      t(100) = t(73)*t(99)
      t(101) = t(98)+t(100)
      A0(1) = t(2)*t(14)*t(39)*(dBsig22A(-1)*t(2)+dBsig22B(-
     &1)*t(6)*7.5D-1)*(-1.0D0)+dBsig21(-1)*t(2)*t(24)*(t(50)+t(11)*t
     &(54)*5.0D-1)*(t(6)*t(30)*5.0D-1+t(11)*t(35)*5.0D-1)+dsig11B*4.0
     &*t(2)*t(12)*t(16)*t(38)*1.25D-1-dsig12(-1)*4.0*t(2)*t(12)*t(22)*t(
     &65)*1.25D-1-dsig12(-1)*4.0*t(12)*t(17)*t(22)*t(71)*1.25D-1+dBsig21
     &(-1)*t(2)*t(24)*t(26)*t(37)*t(45)*7.5D-1
      iA0(1,1)=0; iA0(1,2)=-1
      A0(2) = t(12)*t(16)*t(38)*(dsig11A*4.0*t(16)*1.0D0-dsig11
     &B*4.0*(t(2)*2.5D-1-t(42)*2.5D-1))*5.0D-1+t(2)*t(14)*t(39)*(dBsi
     &g22A(-1)*t(2)*1.0D0-dBsig22B(-1)*t(41)*7.5D-1)*1.0D0+t(39)*t
     &(42)*t(67)*(dBsig22A(0)*t(42)+dBsig22B(0)*t(41)*7.5D-1)-dBsi
     &g21(0)*t(24)*t(42)*(t(44)+t(59)*t(83)*5.0D-1)*(t(55)+t(59)*t(9
     &3)*5.0D-1)*1.0D0-dBsig21(-1)*t(2)*t(24)*t(45)*(t(50)+t(55))*1.
     &0D0+dsig12(-1)*4.0*t(2)*t(12)*t(22)*t(65)*1.25D-1
     &-dsig12(-1)*4.0*t(1
     &2)*t(16)*t(22)*t(71)*1.25D-1-dsig12(0)*4.0*t(12)*t(16)*t(77)*t(101
     &)*1.25D-1-dsig12(0)*4.0*t(12)*t(42)*t(77)*t(97)*1.25D-1+dBsig21(
     &-1)*t(2)*t(24)*t(37)*t(41)*t(45)*7.5D-1-dBsig21(0)*t(24)*t(41)
     &*t(42)*t(85)*t(89)*7.5D-1
      iA0(2,1)=0; iA0(2,2)=0
      A0(3) = -t(39)*t(42)*t(67)*(dBsig22A(0)*t(42)*1.0D0-dBsi
     &g22B(0)*t(87)*7.5D-1)+dBsig21(0)*t(24)*t(42)*t(89)*(t(55)+t(
     &87)*t(93)*5.0D-1)-dsig11B*4.0*t(12)*t(16)*t(38)*t(42)*1.25D-1+d
     &sig12(0)*4.0*t(12)*t(42)*t(77)*t(97)*1.25D-1
     &-dsig12(0)*4.0*t(12)*t(7
     &7)*t(94)*t(101)*1.25D-1-dBsig21(0)*t(24)*t(42)*t(85)*t(87)*t(8
     &9)*7.5D-1
      iA0(3,1)=0; iA0(3,2)= 1
      A0(4) = dBsig22B(-1)*t(2)*t(14)*t(26)*t(39)+dsig11B*4.0
     &*t(12)*t(16)*t(19)*t(38)*1.25D-1-dsig12(-1)*4.0*t(12)*t(19)*t(22)*
     &t(65)*1.25D-1+dsig12(-1)*4.0*t(12)*t(17)*t(22)*t(71)*1.25D-1-dBsig
     &21(-1)*t(2)*t(24)*t(26)*t(37)*t(45)*1.0D0
      iA0(4,1)=1; iA0(4,2)=-1
      A0(5) = t(12)*t(16)*t(38)*(dsig11A*4.0*t(16)+dsig11B*4.0
     &*(t(19)*2.5D-1-t(73)*2.5D-1))*(-5.0D-1)+dBsig22B(-1)*t(2)*t(14
     &)*t(39)*t(41)-dBsig22B(0)*t(39)*t(41)*t(42)*t(67)*1.0D0+dsig12
     &(-1)*4.0*t(12)*t(19)*t(22)*t(65)*1.25D-1+dsig12(-1)
     &*4.0*t(12)*t(16)*
     &t(22)*t(71)*1.25D-1+dsig12(0)*4.0*t(12)*t(16)*t(77)*t(101)*1.25D-1
     &-dsig12(0)*4.0*t(12)*t(73)*t(77)*t(97)*1.25D-1-dBsig21(-1)*t(2)*
     &t(24)*t(37)*t(41)*t(45)*1.0D0+dBsig21(0)*t(24)*t(41)*t(42)*t(8
     &5)*t(89)
      iA0(5,1)=1; iA0(5,2)=0
      A0(6) = dBsig22B(0)*t(39)*t(42)*t(67)*t(87)*(-1.0D0)-dsi
     &g11B*4.0*t(12)*t(16)*t(38)*t(73)*1.25D-1+dsig12(0)*4.0*t(12)*t(73
     &)*t(77)*t(97)*1.25D-1+dsig12(0)*4.0*t(12)*t(77)*t(94)*t(101)*1.25D
     &-1+dBsig21(0)*t(24)*t(42)*t(85)*t(87)*t(89)
      iA0(6,1)=1; iA0(6,2)=1
      A0(7) = dBsig22B(-1)*t(2)*t(14)*t(26)*t(39)*(-2.5D-1)+dB
     &sig21(-1)*t(2)*t(24)*t(26)*t(37)*t(45)*2.5D-1
      iA0(7,1)=2; iA0(7,2)=-1
      A0(8) = dBsig22B(-1)*t(2)*t(14)*t(39)*t(41)*(-2.5D-1)+dB
     &sig22B(0)*t(39)*t(41)*t(42)*t(67)*2.5D-1+dBsig21(-1)*t(2)*t(
     &24)*t(37)*t(41)*t(45)*2.5D-1-dBsig21(0)*t(24)*t(41)*t(42)*t(85
     &)*t(89)*2.5D-1
      iA0(8,1)=2; iA0(8,2)=0
      A0(9) = dBsig22B(0)*t(39)*t(42)*t(67)*t(87)*2.5D-1-dBsig
     &21(0)*t(24)*t(42)*t(85)*t(87)*t(89)*2.5D-1
      iA0(9,1)=2; iA0(9,2)= 1

      return
      end    
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      subroutine KnonlinearBNDY1i(A0,iA0,P,Q,ds1,ds2,dh1,dh2,
     &   dBsig21,dsig12,dBsig11A,dBsig11B,dsig22A,dsig22B)
      implicit none
      integer iA0(9,2)
      real*8 t(150),A0(9)
      real*8 ds1(-1:1,0:1),ds2(-1:1,0:1),dh2(-1:1,0:1),dh1(-1:1,0:1)
     &  ,P(-1:1,0:2),Q(-1:1,0:2)
     &  ,dBsig21(-1:1),dsig12(-1:1)
     &  ,dBsig11A(-1:1),dBsig11B(-1:1)
     &  ,dsig22A,dsig22B

      t(2) = 1.0D0/ds1(-1,0)
      t(3) = ds2(-1,0)*1.5D0
      t(7) = ds2(-1,+1)*(-5.0D-1)
      t(4) = -t(7)
      t(5) = t(3)-t(4)
      t(6) = 1.0D0/t(5)
      t(8) = ds2(0,0)*1.5D0
      t(31) = ds2(0,+1)*(-5.0D-1)
      t(9) = -t(31)
      t(10) = t(8)-t(9)
      t(11) = 1.0D0/t(10)
      t(12) = P(0,0)*(-1.0D0)
      t(13) = P(-1,0)+t(12)
      t(14) = 1.0D0/ds1(-1,+1)
      t(15) = 1.0D0/dh2(0,0)
      t(16) = t(2)*t(13)*2.5D-1
      t(17) = P(0,+1)*(-1.0D0)
      t(18) = P(-1,+1)+t(17)
      t(19) = t(14)*t(18)*2.5D-1
      t(20) = t(16)+t(19)
      t(21) = 1.0D0/ds2(-1,0)
      t(22) = Q(0,+1)*(-1.0D0)
      t(23) = 1.0D0/ds2(0,0)
      t(24) = 1.0D0/dh1(-1,0)
      t(25) = t(3)+t(7)
      t(26) = 1.0D0/t(25)
      t(27) = P(-1,0)*1.5D0
      t(28) = P(-1,+1)*(-2.0D0)
      t(29) = P(-1,+2)*5.0D-1
      t(30) = t(27)+t(28)+t(29)
      t(32) = P(0,0)*1.5D0
      t(33) = P(0,+1)*(-2.0D0)
      t(34) = P(0,+2)*5.0D-1
      t(35) = t(32)+t(33)+t(34)
      t(36) = Q(0,0)*(-1.0D0)
      t(37) = Q(-1,0)+t(36)
      t(38) = t(2)*t(37)
      t(39) = Q(-1,+1)+t(22)
      t(40) = t(14)*t(39)
      t(41) = t(38)+t(40)
      t(42) = Q(-1,+1)*(-1.0D0)
      t(43) = Q(-1,0)+t(42)
      t(44) = t(21)*t(43)
      t(45) = Q(0,0)+t(22)
      t(46) = t(23)*t(45)
      t(47) = t(44)+t(46)
      t(48) = P(0,0)+t(17)
      t(49) = 1.0D0/dh1(0,0)
      t(50) = t(26)*t(30)*5.0D-1
      t(51) = t(8)+t(31)
      t(52) = 1.0D0/t(51)
      t(53) = t(35)*t(52)*5.0D-1
      t(54) = t(50)+t(53)
      t(55) = 1.0D0/ds1(0,0)
      t(56) = Q(-1,0)*1.5D0
      t(57) = Q(-1,+1)*(-2.0D0)
      t(58) = Q(-1,+2)*5.0D-1
      t(59) = t(56)+t(57)+t(58)
      t(60) = t(26)*t(59)*5.0D-1
      t(61) = Q(0,0)*1.5D0
      t(62) = Q(0,+1)*(-2.0D0)
      t(63) = Q(0,+2)*5.0D-1
      t(64) = t(61)+t(62)+t(63)
      t(65) = t(52)*t(64)*5.0D-1
      t(66) = ds2(+1,0)*1.5D0
      t(79) = ds2(+1,+1)*(-5.0D-1)
      t(67) = -t(79)
      t(68) = t(66)-t(67)
      t(69) = 1.0D0/t(68)
      t(70) = P(+1,0)*(-1.0D0)
      t(71) = P(0,0)+t(70)
      t(72) = 1.0D0/ds1(0,+1)
      t(73) = t(55)*t(71)*2.5D-1
      t(74) = P(+1,+1)*(-1.0D0)
      t(75) = P(0,+1)+t(74)
      t(76) = t(72)*t(75)*2.5D-1
      t(77) = t(73)+t(76)
      t(78) = Q(+1,+1)*(-1.0D0)
      t(80) = P(+1,0)*1.5D0
      t(81) = P(+1,+1)*(-2.0D0)
      t(82) = P(+1,+2)*5.0D-1
      t(83) = t(80)+t(81)+t(82)
      t(84) = Q(+1,0)*(-1.0D0)
      t(85) = Q(0,0)+t(84)
      t(86) = t(55)*t(85)
      t(87) = Q(0,+1)+t(78)
      t(88) = t(72)*t(87)
      t(89) = t(86)+t(88)
      t(90) = 1.0D0/ds2(+1,0)
      t(91) = Q(+1,0)+t(78)
      t(92) = t(90)*t(91)
      t(93) = t(46)+t(92)
      t(94) = t(66)+t(79)
      t(95) = 1.0D0/t(94)
      t(96) = t(83)*t(95)*5.0D-1
      t(97) = t(53)+t(96)
      t(98) = Q(+1,0)*1.5D0
      t(99) = Q(+1,+1)*(-2.0D0)
      t(100) = Q(+1,+2)*5.0D-1
      t(101) = t(98)+t(99)+t(100)
      A0(1) = t(2)*t(13)*t(49)*(dBsig11A(-1)*t(2)+dBsig11B(-1)*t
     &(6)*7.5D-1)*(-1.0D0)+dBsig21(-1)*t(2)*t(24)*(t(60)+t(11)*t(64)*5
     &.0D-1)*(t(6)*t(30)*5.0D-1+t(11)*t(35)*5.0D-1)
     &-dsig12(-1)*4.0*t(2)*t(
     &15)*t(20)*t(47)*1.25D-1
     &-dsig12(-1)*4.0*t(15)*t(20)*t(21)*t(41)*1.25D
     &-1+dsig22B*4.0*t(2)*t(15)*t(23)*t(48)*1.25D-1
     &+dBsig21(-1)*t(2)*t(24)
     &*t(26)*t(37)*t(54)*7.5D-1
      A0(2) = dBsig11B(-1)*t(2)*t(13)*t(26)*t(49)-dsig12(-1)*4.0*t(
     &14)*t(15)*t(20)*t(47)*1.25D-1
     &+dsig12(-1)*4.0*t(15)*t(20)*t(21)*t(41)
     &*1.25D-1+dsig22B*4.0*t(14)*t(15)*t(23)*t(48)*1.25D-1
     &-dBsig21(-1)*t(2
     &)*t(24)*t(26)*t(37)*t(54)*1.0D0
      A0(3) = dBsig11B(-1)*t(2)*t(13)*t(26)*t(49)*(-2.5D-1)+dBsi
     &g21(-1)*t(2)*t(24)*t(26)*t(37)*t(54)*2.5D-1
      A0(4) = t(2)*t(13)*t(49)*(dBsig11A(-1)*t(2)*1.0D0-dBsig11B
     &(-1)*t(52)*7.5D-1)*1.0D0
     &+t(15)*t(23)*t(48)*(dsig22A*4.0*t(23)*1.0D0-
     &dsig22B*4.0*(t(2)*2.5D-1-t(55)*2.5D-1))*5.0D-1
     &+t(49)*t(55)*t(71)*(dB
     &sig11B(0)*t(52)*7.5D-1+dBsig11A(0)*t(55))-dBsig21(0)*t(24)*t(55)
     &*(t(53)+t(69)*t(83)*5.0D-1)*(t(65)+t(69)*t(101)*5.0D-1)*1.0D0-dB
     &sig21(-1)*4.0*t(2)*t(24)*t(54)*(t(60)+t(65))*1.0D0
     &+dsig12(-1)*4.0*t(2)*t
     &(15)*t(20)*t(47)*1.25D-1
     &-dsig12(-1)*4.0*t(15)*t(20)*t(23)*t(41)*1.25
     &D-1-dsig12(0)*4.0*t(15)*t(23)*t(77)*t(89)*1.25D-1
     &-dsig12(0)*4.0*t(15)*t(
     &55)*t(77)*t(93)*1.25D-1+dBsig21(-1)*t(2)*t(24)*t(37)*t(52)*t(54)
     &*7.5D-1-dBsig21(0)*t(24)*t(52)*t(55)*t(85)*t(97)*7.5D-1
      A0(5) = t(15)*t(23)*t(48)*(dsig22A*4.0*t(23)+dsig22B*4.0*(t(14)*2.
     &5D-1-t(72)*2.5D-1))*(-5.0D-1)+dBsig11B(-1)*t(2)*t(13)*t(49)*t(52
     &)-dBsig11B(0)*t(49)*t(52)*t(55)*t(71)*1.0D0+dsig12(-1)*t(14)*t(1
     &5)*t(20)*t(47)*1.25D-1
     &+dsig12(-1)*4.0*t(15)*t(20)*t(23)*t(41)*1.25D-
     &1+dsig12(0)*4.0*t(15)*t(23)*t(77)*t(89)*1.25D-1
     &-dsig12(0)*4.0*t(15)*t(72
     &)*t(77)*t(93)*1.25D-1-dBsig21(-1)*t(2)*t(24)*t(37)*t(52)*t(54)*1
     &.0D0+dBsig21(0)*t(24)*t(52)*t(55)*t(85)*t(97)
      A0(6) = dBsig11B(-1)*t(2)*t(13)*t(49)*t(52)*(-2.5D-1)+dBsi
     &g11B(0)*t(49)*t(52)*t(55)*t(71)*2.5D-1+dBsig21(-1)*t(2)*t(24)*t(
     &37)*t(52)*t(54)*2.5D-1-dBsig21(0)*t(24)*t(52)*t(55)*t(85)*t(97)*
     &2.5D-1
      A0(7) = -t(49)*t(55)*t(71)*(dBsig11A(0)*t(55)*1.0D0-dBsig1
     &1B(0)*t(95)*7.5D-1)+dBsig21(0)*t(24)*t(55)*t(97)*(t(65)+t(95)*t(
     &101)*5.0D-1)+dsig12(0)*4.0*t(15)*t(55)*t(77)*t(93)*1.25D-1
     &-dsig12(0)*4.0
     &*t(15)*t(77)*t(89)*t(90)*1.25D-1
     &-dsig22B*4.0*t(15)*t(23)*t(48)*t(55)
     &*1.25D-1-dBsig21(0)*t(24)*t(55)*t(85)*t(95)*t(97)*7.5D-1
      A0(8) = dBsig11B(0)*t(49)*t(55)*t(71)*t(95)*(-1.0D0)+dsig1
     &2(0)*4.0*t(15)*t(72)*t(77)*t(93)*1.25D-1
     &+dsig12(0)*4.0*t(15)*t(77)*t(89)
     &*t(90)*1.25D-1-dsig22B*4.0*t(15)*t(23)*t(48)*t(72)*1.25D-1
     &+dBsig21(0
     &)*t(24)*t(55)*t(85)*t(95)*t(97)
      A0(9) = dBsig11B(0)*t(49)*t(55)*t(71)*t(95)*2.5D-1-dBsig21
     &(0)*t(24)*t(55)*t(85)*t(95)*t(97)*2.5D-1
      iA0(1,1)=-1; iA0(1,2)=0
      iA0(2,1)=-1; iA0(2,2)=1
      iA0(3,1)=-1; iA0(3,2)=2
      iA0(4,1)=0; iA0(4,2)=0
      iA0(5,1)=0; iA0(5,2)=1
      iA0(6,1)=0; iA0(6,2)=2
      iA0(7,1)=1; iA0(7,2)=0
      iA0(8,1)=1; iA0(8,2)=1
      iA0(9,1)=1; iA0(9,2)=2


      return
      end
  

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      subroutine KnonlinearBNDY11(A0,iA0,P,Q,ds1,ds2,dh1,dh2,
     &   dBsig11A,dBsig11B,dBsig22A,dBsig22B,dBsig12F1,dBsig12F2)
      implicit none
      integer iA0(8,2)
      real*8 t(90),A0(8)
      real*8 ds1(0:1,0:1),ds2(0:1,0:1),dh2(0:1,0:1),dh1(0:1,0:1)
     &  ,P(0:2,0:2),Q(0:2,0:2)
     &  ,dBsig11A,dBsig11B,dBsig22A,dBsig22B,dBsig12F1,dBsig12F2

      t(2) = 1.0D0/ds1(0,0)
      t(3) = 1.0D0/ds2(0,0)
      t(4) = ds1(0,0)*1.5D0
      t(8) = ds1(+1,0)*(-5.0D-1)
      t(5) = -t(8)
      t(6) = t(4)-t(5)
      t(7) = 1.0D0/t(6)
      t(9) = ds1(0,+1)*1.5D0
      t(32) = ds1(+1,+1)*(-5.0D-1)
      t(10) = -t(32)
      t(11) = t(9)-t(10)
      t(12) = 1.0D0/t(11)
      t(13) = ds2(0,0)*1.5D0
      t(19) = ds2(0,+1)*(-5.0D-1)
      t(14) = -t(19)
      t(15) = t(13)-t(14)
      t(16) = 1.0D0/t(15)
      t(17) = P(0,0)*1.5D0
      t(18) = P(+1,+1)*(-2.0D0)
      t(20) = Q(0,0)*1.5D0
      t(21) = ds2(+1,0)*1.5D0
      t(42) = ds2(+1,+1)*(-5.0D-1)
      t(22) = -t(42)
      t(23) = t(21)-t(22)
      t(24) = 1.0D0/t(23)
      t(25) = Q(+1,+1)*(-2.0D0)
      t(26) = 1.0D0/dh2(0,0)
      t(27) = t(4)+t(8)
      t(28) = 1.0D0/t(27)
      t(29) = P(+1,0)*(-2.0D0)
      t(30) = P(+2,0)*5.0D-1
      t(31) = t(17)+t(29)+t(30)
      t(33) = P(0,+1)*1.5D0
      t(34) = P(+2,+1)*5.0D-1
      t(35) = t(18)+t(33)+t(34)
      t(36) = 1.0D0/dh1(0,0)
      t(37) = t(13)+t(19)
      t(38) = 1.0D0/t(37)
      t(39) = P(0,+1)*(-2.0D0)
      t(40) = P(0,+2)*5.0D-1
      t(41) = t(17)+t(39)+t(40)
      t(43) = P(+1,0)*1.5D0
      t(44) = P(+1,+2)*5.0D-1
      t(45) = t(18)+t(43)+t(44)
      t(46) = 1.0D0/dh2(0,0)
      t(47) = t(9)+t(32)
      t(48) = 1.0D0/t(47)
      t(49) = P(0,+1)*(-1.0D0)
      t(50) = P(0,0)+t(49)
      t(51) = t(28)*t(31)*5.0D-1
      t(52) = t(35)*t(48)*5.0D-1
      t(53) = t(51)+t(52)
      t(54) = Q(+1,0)*(-2.0D0)
      t(55) = Q(+2,0)*5.0D-1
      t(56) = t(20)+t(54)+t(55)
      t(57) = t(28)*t(56)*5.0D-1
      t(58) = Q(0,+1)*1.5D0
      t(59) = Q(+2,+1)*5.0D-1
      t(60) = t(25)+t(58)+t(59)
      t(61) = 1.0D0/dh1(0,0)
      t(62) = P(+1,0)*(-1.0D0)
      t(63) = P(0,0)+t(62)
      t(64) = Q(0,+1)*(-1.0D0)
      t(65) = Q(0,0)+t(64)
      t(66) = t(38)*t(41)*5.0D-1
      t(67) = t(21)+t(42)
      t(68) = 1.0D0/t(67)
      t(69) = t(45)*t(68)*5.0D-1
      t(70) = t(66)+t(69)
      t(71) = Q(+1,0)*(-1.0D0)
      t(72) = Q(0,0)+t(71)
      t(73) = Q(0,+1)*(-2.0D0)
      t(74) = Q(0,+2)*5.0D-1
      t(75) = t(20)+t(73)+t(74)
      t(76) = t(38)*t(75)*5.0D-1
      t(77) = Q(+1,0)*1.5D0
      t(78) = Q(+1,+2)*5.0D-1
      t(79) = t(25)+t(77)+t(78)
      A0(1) = t(2)*t(61)*t(63)*(dBsig11A*t(2)+dBsig11B*t(16)*7
     &.5D-1)*2.0D0+t(3)*t(46)*t(50)*(dBsig22A*t(3)+dBsig22A*t(7)*7.5D-
     &1)*2.0D0-dBsig12F1*t(2)*t(36)*(t(76)+t(24)*t(79)*5.0D-1)*(t(16)*
     &t(41)*5.0D-1+t(24)*t(45)*5.0D-1)*2.0D0-dBsig12F2*t(3)*t(26)*(t(5
     &7)+t(12)*t(60)*5.0D-1)*(t(7)*t(31)*5.0D-1+t(12)*t(35)*5.0D-1)*2.
     &0D0-dBsig12F1*t(2)*t(36)*t(38)*t(70)*t(72)*1.5D0-dBsig12F2*t(3)*
     &t(26)*t(28)*t(53)*t(65)*1.5D0

      A0(2) = t(3)*t(46)*t(50)*(dBsig22A*t(3)*1.0D0-dBsig22A*t
     &(48)*7.5D-1)*(-2.0D0)+dBsig12F2*t(3)*t(26)*t(53)*(t(57)+t(48)*t(
     &60)*5.0D-1)*2.0D0-dBsig11B*t(2)*t(38)*t(61)*t(63)*2.0D0+dBsig12F
     &1*t(2)*t(36)*t(38)*t(70)*t(72)*2.0D0-dBsig12F2*t(3)*t(26)*t(48)*
     &t(53)*t(65)*1.5D0

      A0(3) = dBsig11B*t(2)*t(38)*t(61)*t(63)*5.0D-1-dBsig12F1
     &*t(2)*t(36)*t(38)*t(70)*t(72)*5.0D-1

      A0(4) = t(2)*t(61)*t(63)*(dBsig11A*t(2)*1.0D0-dBsig11B*t
     &(68)*7.5D-1)*(-2.0D0)+dBsig12F1*t(2)*t(36)*t(70)*(t(76)+t(68)*t(
     &79)*5.0D-1)*2.0D0-dBsig22A*t(3)*t(28)*t(46)*t(50)*2.0D0-dBsig12F
     &1*t(2)*t(36)*t(68)*t(70)*t(72)*1.5D0+dBsig12F2*t(3)*t(26)*t(28)*
     &t(53)*t(65)*2.0D0

      A0(5) = dBsig11B*t(2)*t(61)*t(63)*t(68)*(-2.0D0)-dBsig22
     &A*t(3)*t(46)*t(48)*t(50)*2.0D0+dBsig12F1*t(2)*t(36)*t(68)*t(70)*
     &t(72)*2.0D0+dBsig12F2*t(3)*t(26)*t(48)*t(53)*t(65)*2.0D0

      A0(6) = dBsig11B*t(2)*t(61)*t(63)*t(68)*5.0D-1-dBsig12F1
     &*t(2)*t(36)*t(68)*t(70)*t(72)*5.0D-1

      A0(7) = dBsig22A*t(3)*t(28)*t(46)*t(50)*5.0D-1-dBsig12F2
     &*t(3)*t(26)*t(28)*t(53)*t(65)*5.0D-1

      A0(8) = dBsig22A*t(3)*t(46)*t(48)*t(50)*5.0D-1-dBsig12F2
     &*t(3)*t(26)*t(48)*t(53)*t(65)*5.0D-1
      iA0(1,1)=0; iA0(1,2)=0
      iA0(2,1)=0; iA0(2,2)=1
      iA0(3,1)=0; iA0(3,2)=2
      iA0(4,1)=1; iA0(4,2)=0
      iA0(5,1)=1; iA0(5,2)=1
      iA0(6,1)=1; iA0(6,2)=2
      iA0(7,1)=2; iA0(7,2)=0
      iA0(8,1)=2; iA0(8,2)=1
      return
      end

 

