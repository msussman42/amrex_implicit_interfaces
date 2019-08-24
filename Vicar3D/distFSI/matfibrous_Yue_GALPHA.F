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
     &  ,t(220),e5c,e6e,P(-2:2,-2:2),Q(-2:2,-2:2),A0(10)
      integer iAdd(-4:4,-4:4,1:3)
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
       
       X(0,0)=1.0/ds2(i,j)**2*ds1(i,j)**2
       X(1,1)=X(0,0)
       X(1,0)=-X(0,0)
       X(0,1)=-X(0,0)   
       FlagForceCalculation(i,j)=0
       coefX00=1.0d0   
!! ================== BNDY Nx,1
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   )then
       X(0,0)=1.0/ds2(i-1,j)**2*ds1(i-1,j)**2
       X(-1,1)=X(0,0)
       X(-1,0)=-X(0,0)
       X(0,1)=-X(0,0)   
       FlagForceCalculation(i,j)=0
       coefX00=1.0d0   
!! ================== BNDY 1,Ny
       elseif(((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
       X(0,0)=1.0/ds2(i,j-1)**2*ds1(i,j-1)**2
       X(1,-1)=X(0,0)
       X(1,0)=-X(0,0)
       X(0,-1)=-X(0,0)   
       FlagForceCalculation(i,j)=0
       coefX00=1.0d0 

!! ================== BNDY Nx,Ny
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
       X(0,0)=1.0/ds2(i-1,j-1)**2*ds1(i-1,j-1)**2
       X(-1,-1)=X(0,0)
       X(-1,0)=-X(0,0)
       X(0,-1)=-X(0,0)   
       FlagForceCalculation(i,j)=0
       coefX00=1.0d0 
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

      X(0,-1)=(gam22(i,j)/(dh2(i,j)*ds2(i,jm1)**2) + gam22(i,j)/(dh2(i
     &,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + sigB(1,jm1,2)/(dh2(i,j)*ds2
     &(i,jm1)) - (0.5*sig12(i,jm1))/(ds1(i,j)*ds2(i,jm1)) + (0.75*sigB
     &(1,jm1,3))/(dh2(i,j)*(1.5*ds1(i,jm1) - 0.5*ds1(ip1,jm1))) - (1.0
     &*gam11(i,j)*(e4coefint - 2.0))/(dh2(i,j)*ds1(i,j)**2*ds2(i,jm1))
     & + (2.0*e6e*gam21(i,j))/(dh2(i,j)*ds1(i,j)**2*ds2(i,jm1))

      X(0,0)=(2.0*((0.25*sig12(i,jm1))/ds2(i,jm1) - (0.25*sig12(i,j))/
     &ds2(i,j)))/ds1(i,j) - ((1.0*sigB(1,jm1,2))/ds2(i,jm1) + (1.0*sig
     &B(1,j,2))/ds2(i,j))/dh2(i,j) - (2.0*sig11(i,j))/ds1(i,j)**2 + ((
     &0.75*sigB(1,jm1,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - (0.75*sigB
     &(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)))/dh2(i,j) - (((1.0*gam2
     &2(i,jp1))/(dh2(i,jp1)*ds2(i,j)) + (gam22(i,j)*(1.0/ds2(i,jm1) + 
     &1.0/ds2(i,j)))/dh2(i,j))/ds2(i,j) + (gam22(i,j)*(1.0/ds2(i,jm1) 
     &+ 1.0/ds2(i,j)))/(dh2(i,j)*ds2(i,jm1)))/dh2(i,j) - (1.0*gam11(i,
     &j))/(dh1(ip1,j)*ds1(i,j)**3) + gam11(ip1,j)/(dh1(ip1,j)*ds1(i,j)
     &**2*ds1(ip1,j)) + (1.0*gam11(i,j)*(e4coefint - 2.0)*(1.0/ds2(i,j
     &m1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds1(i,j)**2) - (2.0*e6e*gam21(i,j
     &)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds1(i,j)**2)

      X(0,+1)=(((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*gam22(i,jp
     &1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1))/ds2(i,j) + (1.0*
     &gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + (1.0*sigB
     &(1,j,2))/(dh2(i,j)*ds2(i,j)) + (0.5*sig12(i,j))/(ds1(i,j)*ds2(i,
     &j)) - (0.75*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1
     &,jp1))) - (1.0*gam11(i,j)*(e4coefint - 2.0))/(dh2(i,j)*ds1(i,j)*
     &*2*ds2(i,j)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*ds1(i,j)**2*ds2(i,
     &j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(+1,-1)=(2.0*gam21(ip1,j))/(dh2(ip1,j)*ds1(i,j)*ds1(ip1,j)*ds2(
     &ip1,jm1)) - (1.0*sigB(1,jm1,3))/(dh2(i,j)*(1.5*ds1(i,jm1) - 0.5*
     &ds1(ip1,jm1))) - (2.0*gam11(i,j))/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip
     &1,jm1)) - (2.0*e6e*gam21(i,j))/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip1,j
     &m1)) - (0.5*sig12(i,jm1))/(ds1(i,j)*ds2(ip1,jm1))

      X(+1,0)=(2.0*((0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j
     &))/ds2(ip1,j)))/ds1(i,j) - ((1.0*sigB(1,jm1,3))/(1.5*ds1(i,j) - 
     &0.5*ds1(ip1,j)) - (1.0*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,
     &j)))/dh2(i,j) + (2.0*sig11(i,j))/ds1(i,j)**2 + (1.0*gam11(i,j)*(
     &(1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + (2.0*(1.0/ds2(ip1,j
     &m1) + 1.0/ds2(ip1,j)))/dh2(ip1,j)))/ds1(i,j)**2 - (gam11(ip1,j)*
     &((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + 1.0/(dh1(ip2,j)*ds
     &1(ip1,j))))/(ds1(i,j)*ds1(ip1,j)) - (2.0*gam21(ip1,j)*(1.0/ds2(i
     &p1,jm1) + 1.0/ds2(ip1,j)))/(dh2(ip1,j)*ds1(i,j)*ds1(ip1,j)) + (2
     &.0*e6e*gam21(i,j)*(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/(dh2(ip1,
     &j)*ds1(i,j)**2)

      X(+1,+1)=(0.5*sig12(i,j))/(ds1(i,j)*ds2(ip1,j)) + (1.0*sigB(1,j,
     &3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) - (2.0*gam11(
     &i,j))/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip1,j)) - (2.0*e6e*gam21(i,j))
     &/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip1,j)) + (2.0*gam21(ip1,j))/(dh2(i
     &p1,j)*ds1(i,j)*ds1(ip1,j)*ds2(ip1,j))

      X(+2,-1)=(0.25*sigB(1,jm1,3))/(dh2(i,j)*(1.5*ds1(i,jm1) - 0.5*ds
     &1(ip1,jm1))) - (2.0*gam21(ip1,j))/(dh2(ip2,j)*ds1(i,j)*ds1(ip1,j
     &)*ds2(ip2,jm1))

      X(+2,0)=((0.25*sigB(1,jm1,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - 
     &(0.25*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)))/dh2(i,j) - (
     &1.0*gam11(i,j))/(dh1(ip1,j)*ds1(i,j)**2*ds1(ip1,j)) + (gam11(ip1
     &,j)*((1.0*(1.0/ds1(ip1,j) + 1.0/ds1(ip2,j)))/dh1(ip2,j) + 1.0/(d
     &h1(ip1,j)*ds1(ip1,j))))/(ds1(i,j)*ds1(ip1,j)) + (2.0*gam21(ip1,j
     &)*(1.0/ds2(ip2,jm1) + 1.0/ds2(ip2,j)))/(dh2(ip2,j)*ds1(i,j)*ds1(
     &ip1,j))

      X(+2,+1)=- (0.25*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds
     &1(ip1,jp1))) - (2.0*gam21(ip1,j))/(dh2(ip2,j)*ds1(i,j)*ds1(ip1,j
     &)*ds2(ip2,j))

      X(+3,0)=-(1.0*gam11(ip1,j))/(dh1(ip2,j)*ds1(i,j)*ds1(ip1,j)*ds1(
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
     &+ (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) + s
     &igB(1,jm1,2)/(dh2(i,j)*ds2(i,jm1)) - (0.5*sig12(i,jm1))/(ds1(i,j
     &)*ds2(i,jm1)) + (0.75*sigB(1,jm1,3))/(dh2(i,j)*(1.5*ds1(i,jm1) -
     & 0.5*ds1(ip1,jm1))) - (1.0*gam11(i,j)*(e4coefint - 2.0))/(dh2(i,
     &j)*ds1(i,j)**2*ds2(i,jm1)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*ds1(
     &i,j)**2*ds2(i,jm1))

      X(0,0)=(2.0*((0.25*sig12(i,jm1))/ds2(i,jm1) - (0.25*sig12(i,j))/
     &ds2(i,j)))/ds1(i,j) - ((1.0*sigB(1,jm1,2))/ds2(i,jm1) + (1.0*sig
     &B(1,j,2))/ds2(i,j))/dh2(i,j) - (1.0*(((1.0*gam22(i,jm1))/(dh2(i,
     &jm1)*ds2(i,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j
     &)))/dh2(i,j))/ds2(i,jm1) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0
     &/ds2(i,j)))/(dh2(i,j)*ds2(i,j))))/dh2(i,j) - (2.0*sig11(i,j))/ds
     &1(i,j)**2 + ((0.75*sigB(1,jm1,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)
     &) - (0.75*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)))/dh2(i,j)
     & - (1.0*gam11(i,j))/(dh1(ip1,j)*ds1(i,j)**3) + gam11(ip1,j)/(dh1
     &(ip1,j)*ds1(i,j)**2*ds1(ip1,j)) + (1.0*gam11(i,j)*(e4coefint - 2
     &.0)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds1(i,j)**2) - (2
     &.0*e6e*gam21(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds1
     &(i,j)**2)

      X(0,+1)=(1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)**2) + (1.0*gam
     &22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) + (1.0*sigB(1
     &,j,2))/(dh2(i,j)*ds2(i,j)) + (0.5*sig12(i,j))/(ds1(i,j)*ds2(i,j)
     &) - (0.75*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,j
     &p1))) - (1.0*gam11(i,j)*(e4coefint - 2.0))/(dh2(i,j)*ds1(i,j)**2
     &*ds2(i,j)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*ds1(i,j)**2*ds2(i,j)
     &)

      X(+1,-1)=(2.0*gam21(ip1,j))/(dh2(ip1,j)*ds1(i,j)*ds1(ip1,j)*ds2(
     &ip1,jm1)) - (1.0*sigB(1,jm1,3))/(dh2(i,j)*(1.5*ds1(i,jm1) - 0.5*
     &ds1(ip1,jm1))) - (2.0*gam11(i,j))/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip
     &1,jm1)) - (2.0*e6e*gam21(i,j))/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip1,j
     &m1)) - (0.5*sig12(i,jm1))/(ds1(i,j)*ds2(ip1,jm1))

      X(+1,0)=(2.0*((0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j
     &))/ds2(ip1,j)))/ds1(i,j) - ((1.0*sigB(1,jm1,3))/(1.5*ds1(i,j) - 
     &0.5*ds1(ip1,j)) - (1.0*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,
     &j)))/dh2(i,j) + (2.0*sig11(i,j))/ds1(i,j)**2 + (1.0*gam11(i,j)*(
     &(1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + (2.0*(1.0/ds2(ip1,j
     &m1) + 1.0/ds2(ip1,j)))/dh2(ip1,j)))/ds1(i,j)**2 - (gam11(ip1,j)*
     &((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + 1.0/(dh1(ip2,j)*ds
     &1(ip1,j))))/(ds1(i,j)*ds1(ip1,j)) - (2.0*gam21(ip1,j)*(1.0/ds2(i
     &p1,jm1) + 1.0/ds2(ip1,j)))/(dh2(ip1,j)*ds1(i,j)*ds1(ip1,j)) + (2
     &.0*e6e*gam21(i,j)*(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/(dh2(ip1,
     &j)*ds1(i,j)**2)

      X(+1,+1)=(0.5*sig12(i,j))/(ds1(i,j)*ds2(ip1,j)) + (1.0*sigB(1,j,
     &3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) - (2.0*gam11(
     &i,j))/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip1,j)) - (2.0*e6e*gam21(i,j))
     &/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip1,j)) + (2.0*gam21(ip1,j))/(dh2(i
     &p1,j)*ds1(i,j)*ds1(ip1,j)*ds2(ip1,j))

      X(+2,-1)=(0.25*sigB(1,jm1,3))/(dh2(i,j)*(1.5*ds1(i,jm1) - 0.5*ds
     &1(ip1,jm1))) - (2.0*gam21(ip1,j))/(dh2(ip2,j)*ds1(i,j)*ds1(ip1,j
     &)*ds2(ip2,jm1))

      X(+2,0)=((0.25*sigB(1,jm1,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - 
     &(0.25*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)))/dh2(i,j) - (
     &1.0*gam11(i,j))/(dh1(ip1,j)*ds1(i,j)**2*ds1(ip1,j)) + (gam11(ip1
     &,j)*((1.0*(1.0/ds1(ip1,j) + 1.0/ds1(ip2,j)))/dh1(ip2,j) + 1.0/(d
     &h1(ip1,j)*ds1(ip1,j))))/(ds1(i,j)*ds1(ip1,j)) + (2.0*gam21(ip1,j
     &)*(1.0/ds2(ip2,jm1) + 1.0/ds2(ip2,j)))/(dh2(ip2,j)*ds1(i,j)*ds1(
     &ip1,j))

      X(+2,+1)=- (0.25*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds
     &1(ip1,jp1))) - (2.0*gam21(ip1,j))/(dh2(ip2,j)*ds1(i,j)*ds1(ip1,j
     &)*ds2(ip2,j))

      X(+3,0)=-(1.0*gam11(ip1,j))/(dh1(ip2,j)*ds1(i,j)*ds1(ip1,j)*ds1(
     &ip2,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY N,2  
       elseif((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .eq. 2)) then
!= ===============  N,2 VER2
!CCCCCCCC+X-xCCCCCCCC
      X(-3,0)=-(1.0*gam11(im1,j))/(dh1(im2,j)*ds1(im3,j)*ds1(im2,j)*ds
     &1(im1,j))

      X(-2,-1)=(0.25*sigB(2,jm1,3))/(dh2(i,j)*(0.5*ds1(im2,jm1) - 1.5*
     &ds1(im1,jm1))) - (2.0*gam21(im1,j))/(dh2(im2,j)*ds1(im2,j)*ds1(i
     &m1,j)*ds2(im2,jm1))

      X(-2,0)=((0.25*sigB(2,jm1,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) 
     &- (0.25*sigB(2,j,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))/dh2(i,j)
     & + gam11(i,j)/(dh1(im1,j)*ds1(im2,j)*ds1(im1,j)**2) + (1.0*gam11
     &(im1,j)*((1.0/ds1(im3,j) + 1.0/ds1(im2,j))/dh1(im2,j) + 1.0/(dh1
     &(im1,j)*ds1(im2,j))))/(ds1(im2,j)*ds1(im1,j)) + (2.0*gam21(im1,j
     &)*(1.0/ds2(im2,jm1) + 1.0/ds2(im2,j)))/(dh2(im2,j)*ds1(im2,j)*ds
     &1(im1,j))

      X(-2,+1)=- (0.25*sigB(2,j,3))/(dh2(i,j)*(0.5*ds1(im2,jp1) - 1.5*
     &ds1(im1,jp1))) - (2.0*gam21(im1,j))/(dh2(im2,j)*ds1(im2,j)*ds1(i
     &m1,j)*ds2(im2,j))

      X(-1,-1)=(0.5*sig12(im1,jm1))/(ds1(im1,j)*ds2(im1,jm1)) - (1.0*s
     &igB(2,jm1,3))/(dh2(i,j)*(0.5*ds1(im2,jm1) - 1.5*ds1(im1,jm1))) +
     & (2.0*gam11(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,jm1)) - (2.0
     &*e6e*gam21(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,jm1)) + (2.0*
     &gam21(im1,j))/(dh2(im1,j)*ds1(im2,j)*ds1(im1,j)*ds2(im1,jm1))

      X(-1,0)=(2.0*sig11(im1,j))/ds1(im1,j)**2 - ((1.0*sigB(2,jm1,3))/
     &(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) - (1.0*sigB(2,j,3))/(0.5*ds1(i
     &m2,j) - 1.5*ds1(im1,j)))/dh2(i,j) - (2.0*((0.25*sig12(im1,jm1))/
     &ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(im1,j)))/ds1(im1,j) - (ga
     &m11(i,j)*((1.0/ds1(im2,j) + 1.0/ds1(im1,j))/dh1(im1,j) + (2.0*(1
     &.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/dh2(im1,j)))/ds1(im1,j)**2 - 
     &(1.0*gam11(im1,j)*((1.0*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(i
     &m1,j) + 1.0/(dh1(im2,j)*ds1(im2,j))))/(ds1(im2,j)*ds1(im1,j)) - 
     &(2.0*gam21(im1,j)*(1.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/(dh2(im1,
     &j)*ds1(im2,j)*ds1(im1,j)) + (2.0*e6e*gam21(i,j)*(1.0/ds2(im1,jm1
     &) + 1.0/ds2(im1,j)))/(dh2(im1,j)*ds1(im1,j)**2)

      X(-1,+1)=(1.0*sigB(2,j,3))/(dh2(i,j)*(0.5*ds1(im2,jp1) - 1.5*ds1
     &(im1,jp1))) - (0.5*sig12(im1,j))/(ds1(im1,j)*ds2(im1,j)) + (2.0*
     &gam11(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,j)) - (2.0*e6e*gam
     &21(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,j)) + (2.0*gam21(im1,
     &j))/(dh2(im1,j)*ds1(im2,j)*ds1(im1,j)*ds2(im1,j))

      X(0,-1)=(gam22(i,j)/(dh2(i,j)*ds2(i,jm1)**2) + gam22(i,j)/(dh2(i
     &,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + sigB(2,jm1,2)/(dh2(i,j)*ds2
     &(i,jm1)) + (0.5*sig12(im1,jm1))/(ds1(im1,j)*ds2(i,jm1)) + (0.75*
     &sigB(2,jm1,3))/(dh2(i,j)*(0.5*ds1(im2,jm1) - 1.5*ds1(im1,jm1))) 
     &+ (gam11(i,j)*(e4coefint - 2.0))/(dh2(i,j)*ds1(im1,j)**2*ds2(i,j
     &m1)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*ds1(im1,j)**2*ds2(i,jm1))

      X(0,0)=((0.75*sigB(2,jm1,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) -
     & (0.75*sigB(2,j,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))/dh2(i,j) 
     &- (2.0*((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(im1,j))/d
     &s2(i,j)))/ds1(im1,j) - (2.0*sig11(im1,j))/ds1(im1,j)**2 - ((1.0*
     &sigB(2,jm1,2))/ds2(i,jm1) + (1.0*sigB(2,j,2))/ds2(i,j))/dh2(i,j)
     & - (((1.0*gam22(i,jp1))/(dh2(i,jp1)*ds2(i,j)) + (gam22(i,j)*(1.0
     &/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j))/ds2(i,j) + (gam22(i,j)*(1
     &.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds2(i,jm1)))/dh2(i,j) +
     & (1.0*gam11(i,j))/(dh1(im1,j)*ds1(im1,j)**3) + (1.0*gam11(im1,j)
     &)/(dh1(im1,j)*ds1(im2,j)*ds1(im1,j)**2) - (gam11(i,j)*(e4coefint
     & - 2.0)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds1(im1,j)**2
     &) - (2.0*e6e*gam21(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,
     &j)*ds1(im1,j)**2)

      X(0,+1)=(((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*gam22(i,jp
     &1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1))/ds2(i,j) + (1.0*
     &gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j)))/dh2(i,j) + (1.0*sigB
     &(2,j,2))/(dh2(i,j)*ds2(i,j)) - (0.5*sig12(im1,j))/(ds1(im1,j)*ds
     &2(i,j)) - (0.75*sigB(2,j,3))/(dh2(i,j)*(0.5*ds1(im2,jp1) - 1.5*d
     &s1(im1,jp1))) + (1.0*gam11(i,j)*(e4coefint - 2.0))/(dh2(i,j)*ds1
     &(im1,j)**2*ds2(i,j)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*ds1(im1,j)
     &**2*ds2(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
      FlagForceCalculation(i,j)=1
      coefX00=0.0d0
!! ================== BNDY N,N-1  
       elseif((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .eq. ns_ibmin-1)) then
!= ===============  N,N-1 VER2
!CCCCCCCC+X-xCCCCCCCC

      X(-3,0)=-(1.0*gam11(im1,j))/(dh1(im2,j)*ds1(im3,j)*ds1(im2,j)*ds
     &1(im1,j))

      X(-2,-1)=(0.25*sigB(2,jm1,3))/(dh2(i,j)*(0.5*ds1(im2,jm1) - 1.5*
     &ds1(im1,jm1))) - (2.0*gam21(im1,j))/(dh2(im2,j)*ds1(im2,j)*ds1(i
     &m1,j)*ds2(im2,jm1))

      X(-2,0)=((0.25*sigB(2,jm1,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) 
     &- (0.25*sigB(2,j,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))/dh2(i,j)
     & + gam11(i,j)/(dh1(im1,j)*ds1(im2,j)*ds1(im1,j)**2) + (1.0*gam11
     &(im1,j)*((1.0/ds1(im3,j) + 1.0/ds1(im2,j))/dh1(im2,j) + 1.0/(dh1
     &(im1,j)*ds1(im2,j))))/(ds1(im2,j)*ds1(im1,j)) + (2.0*gam21(im1,j
     &)*(1.0/ds2(im2,jm1) + 1.0/ds2(im2,j)))/(dh2(im2,j)*ds1(im2,j)*ds
     &1(im1,j))

      X(-2,+1)=- (0.25*sigB(2,j,3))/(dh2(i,j)*(0.5*ds1(im2,jp1) - 1.5*
     &ds1(im1,jp1))) - (2.0*gam21(im1,j))/(dh2(im2,j)*ds1(im2,j)*ds1(i
     &m1,j)*ds2(im2,j))

      X(-1,-1)=(0.5*sig12(im1,jm1))/(ds1(im1,j)*ds2(im1,jm1)) - (1.0*s
     &igB(2,jm1,3))/(dh2(i,j)*(0.5*ds1(im2,jm1) - 1.5*ds1(im1,jm1))) +
     & (2.0*gam11(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,jm1)) - (2.0
     &*e6e*gam21(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,jm1)) + (2.0*
     &gam21(im1,j))/(dh2(im1,j)*ds1(im2,j)*ds1(im1,j)*ds2(im1,jm1))

      X(-1,0)=(2.0*sig11(im1,j))/ds1(im1,j)**2 - ((1.0*sigB(2,jm1,3))/
     &(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) - (1.0*sigB(2,j,3))/(0.5*ds1(i
     &m2,j) - 1.5*ds1(im1,j)))/dh2(i,j) - (2.0*((0.25*sig12(im1,jm1))/
     &ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(im1,j)))/ds1(im1,j) - (ga
     &m11(i,j)*((1.0/ds1(im2,j) + 1.0/ds1(im1,j))/dh1(im1,j) + (2.0*(1
     &.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/dh2(im1,j)))/ds1(im1,j)**2 - 
     &(1.0*gam11(im1,j)*((1.0*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(i
     &m1,j) + 1.0/(dh1(im2,j)*ds1(im2,j))))/(ds1(im2,j)*ds1(im1,j)) - 
     &(2.0*gam21(im1,j)*(1.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/(dh2(im1,
     &j)*ds1(im2,j)*ds1(im1,j)) + (2.0*e6e*gam21(i,j)*(1.0/ds2(im1,jm1
     &) + 1.0/ds2(im1,j)))/(dh2(im1,j)*ds1(im1,j)**2)

      X(-1,+1)=(1.0*sigB(2,j,3))/(dh2(i,j)*(0.5*ds1(im2,jp1) - 1.5*ds1
     &(im1,jp1))) - (0.5*sig12(im1,j))/(ds1(im1,j)*ds2(im1,j)) + (2.0*
     &gam11(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,j)) - (2.0*e6e*gam
     &21(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,j)) + (2.0*gam21(im1,
     &j))/(dh2(im1,j)*ds1(im2,j)*ds1(im1,j)*ds2(im1,j))

      X(0,-2)=-(1.0*gam22(i,jm1))/(dh2(i,jm1)*dh2(i,j)*ds2(i,jm2)*ds2(
     &i,jm1))

      X(0,-1)=(1.0*(((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)) + (gam22(i
     &,jm1)*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i,jm1))/ds2(i,jm1) 
     &+ (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) + s
     &igB(2,jm1,2)/(dh2(i,j)*ds2(i,jm1)) + (0.5*sig12(im1,jm1))/(ds1(i
     &m1,j)*ds2(i,jm1)) + (0.75*sigB(2,jm1,3))/(dh2(i,j)*(0.5*ds1(im2,
     &jm1) - 1.5*ds1(im1,jm1))) + (gam11(i,j)*(e4coefint - 2.0))/(dh2(
     &i,j)*ds1(im1,j)**2*ds2(i,jm1)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*
     &ds1(im1,j)**2*ds2(i,jm1))

      X(0,0)=((0.75*sigB(2,jm1,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) -
     & (0.75*sigB(2,j,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))/dh2(i,j) 
     &- (2.0*((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(im1,j))/d
     &s2(i,j)))/ds1(im1,j) - (1.0*(((1.0*gam22(i,jm1))/(dh2(i,jm1)*ds2
     &(i,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(
     &i,j))/ds2(i,jm1) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j
     &)))/(dh2(i,j)*ds2(i,j))))/dh2(i,j) - (2.0*sig11(im1,j))/ds1(im1,
     &j)**2 - ((1.0*sigB(2,jm1,2))/ds2(i,jm1) + (1.0*sigB(2,j,2))/ds2(
     &i,j))/dh2(i,j) + (1.0*gam11(i,j))/(dh1(im1,j)*ds1(im1,j)**3) + (
     &1.0*gam11(im1,j))/(dh1(im1,j)*ds1(im2,j)*ds1(im1,j)**2) - (gam11
     &(i,j)*(e4coefint - 2.0)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,
     &j)*ds1(im1,j)**2) - (2.0*e6e*gam21(i,j)*(1.0/ds2(i,jm1) + 1.0/ds
     &2(i,j)))/(dh2(i,j)*ds1(im1,j)**2)

      X(0,+1)=(1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)**2) + (1.0*gam
     &22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) + (1.0*sigB(2
     &,j,2))/(dh2(i,j)*ds2(i,j)) - (0.5*sig12(im1,j))/(ds1(im1,j)*ds2(
     &i,j)) - (0.75*sigB(2,j,3))/(dh2(i,j)*(0.5*ds1(im2,jp1) - 1.5*ds1
     &(im1,jp1))) + (1.0*gam11(i,j)*(e4coefint - 2.0))/(dh2(i,j)*ds1(i
     &m1,j)**2*ds2(i,j)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*ds1(im1,j)**
     &2*ds2(i,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
      FlagForceCalculation(i,j)=1
      coefX00=0.0d0

!! ================== BNDY 2,N  
       elseif((i.eq. 2) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (j .eq. ns_ibmin)) then
!= ===============  2,N VER2
!CCCCCCCC+X-xCCCCCCCC

      X(-1,-2)=(0.25*sigB(4,im1,3))/(dh1(i,j)*(0.5*ds2(im1,jm2) - 1.5*
     &ds2(im1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(im1,jm2)*ds2
     &(i,jm2)*ds2(i,jm1))

      X(-1,-1)=(0.5*sig12(im1,jm1))/(ds1(im1,jm1)*ds2(i,jm1)) - (1.0*s
     &igB(4,im1,3))/(dh1(i,j)*(0.5*ds2(im1,jm2) - 1.5*ds2(im1,jm1))) +
     & (2.0*gam22(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) - (2.0
     &*e6e*gam21(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) + (2.0*
     &gam21(i,jm1))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm2)*ds2(i,jm1))

      X(-1,0)=(gam11(i,j)/(dh1(i,j)*ds1(im1,j)**2) + gam11(i,j)/(dh1(i
     &,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) + sigB(4,im1,1)/(dh1(i,j)*ds1
     &(im1,j)) + (0.5*sig12(im1,jm1))/(ds1(im1,j)*ds2(i,jm1)) + (0.75*
     &sigB(4,im1,3))/(dh1(i,j)*(0.5*ds2(im1,jm2) - 1.5*ds2(im1,jm1))) 
     &+ (gam22(i,j)*(e4coefint - 2.0))/(dh1(i,j)*ds1(im1,j)*ds2(i,jm1)
     &**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*ds2(i,jm1)**2)

      X(0,-3)=-(1.0*gam22(i,jm1))/(dh2(i,jm2)*ds2(i,jm3)*ds2(i,jm2)*ds
     &2(i,jm1))

      X(0,-2)=((0.25*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) 
     &- (0.25*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(i,j)
     & + gam22(i,j)/(dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) + (1.0*gam22
     &(i,jm1)*((1.0/ds2(i,jm3) + 1.0/ds2(i,jm2))/dh2(i,jm2) + 1.0/(dh2
     &(i,jm1)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) + (2.0*gam21(i,jm1
     &)*(1.0/ds1(im1,jm2) + 1.0/ds1(i,jm2)))/(dh1(i,jm2)*ds2(i,jm2)*ds
     &2(i,jm1))

      X(0,-1)=(2.0*sig22(i,jm1))/ds2(i,jm1)**2 - ((1.0*sigB(4,im1,3))/
     &(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) - (1.0*sigB(4,i,3))/(0.5*ds2(i
     &,jm2) - 1.5*ds2(i,jm1)))/dh1(i,j) - (2.0*((0.25*sig12(im1,jm1))/
     &ds1(im1,jm1) - (0.25*sig12(i,jm1))/ds1(i,jm1)))/ds2(i,jm1) - (ga
     &m22(i,j)*((2.0*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/dh1(i,jm1) +
     & (1.0/ds2(i,jm2) + 1.0/ds2(i,jm1))/dh2(i,jm1)))/ds2(i,jm1)**2 - 
     &(1.0*gam22(i,jm1)*((1.0*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i
     &,jm1) + 1.0/(dh2(i,jm2)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) - 
     &(2.0*gam21(i,jm1)*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/(dh1(i,jm
     &1)*ds2(i,jm2)*ds2(i,jm1)) + (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,jm1
     &) + 1.0/ds1(i,jm1)))/(dh1(i,jm1)*ds2(i,jm1)**2)

      X(0,0)=((0.75*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) -
     & (0.75*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(i,j) 
     &- (2.0*((0.25*sig12(im1,jm1))/ds1(im1,j) - (0.25*sig12(i,jm1))/d
     &s1(i,j)))/ds2(i,jm1) - (2.0*sig22(i,jm1))/ds2(i,jm1)**2 - ((1.0*
     &sigB(4,im1,1))/ds1(im1,j) + (1.0*sigB(4,i,1))/ds1(i,j))/dh1(i,j)
     & - (((1.0*gam11(ip1,j))/(dh1(ip1,j)*ds1(i,j)) + (gam11(i,j)*(1.0
     &/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,j))/ds1(i,j) + (gam11(i,j)*(1
     &.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds1(im1,j)))/dh1(i,j) +
     & (1.0*gam22(i,j))/(dh2(i,jm1)*ds2(i,jm1)**3) + (1.0*gam22(i,jm1)
     &)/(dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) - (gam22(i,j)*(e4coefint
     & - 2.0)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds2(i,jm1)**2
     &) - (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,
     &j)*ds2(i,jm1)**2)

      X(+1,-2)=- (0.25*sigB(4,i,3))/(dh1(i,j)*(0.5*ds2(ip1,jm2) - 1.5*
     &ds2(ip1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(i,jm2)*ds2(i
     &,jm2)*ds2(i,jm1))

      X(+1,-1)=(1.0*sigB(4,i,3))/(dh1(i,j)*(0.5*ds2(ip1,jm2) - 1.5*ds2
     &(ip1,jm1))) - (0.5*sig12(i,jm1))/(ds1(i,jm1)*ds2(i,jm1)) + (2.0*
     &gam22(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) - (2.0*e6e*gam
     &21(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) + (2.0*gam21(i,jm
     &1))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm2)*ds2(i,jm1))

      X(+1,0)=(((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*gam11(ip1,
     &j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/dh1(ip1,j))/ds1(i,j) + (1.0*
     &gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)))/dh1(i,j) + (1.0*sigB
     &(4,i,1))/(dh1(i,j)*ds1(i,j)) - (0.5*sig12(i,jm1))/(ds1(i,j)*ds2(
     &i,jm1)) - (0.75*sigB(4,i,3))/(dh1(i,j)*(0.5*ds2(ip1,jm2) - 1.5*d
     &s2(ip1,jm1))) + (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(i,j)*ds1
     &(i,j)*ds2(i,jm1)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*ds1(i,j)*d
     &s2(i,jm1)**2)

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
      FlagForceCalculation(i,j)=1
      coefX00=0.0d0
!! ================== BNDY N-1,N  
       elseif((i.eq. nq_ibmin-1) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (j .eq. ns_ibmin)) then
!= ===============  N-1,N VER2
!CCCCCCCC+X-xCCCCCCCC

      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,-2)=(0.25*sigB(4,im1,3))/(dh1(i,j)*(0.5*ds2(im1,jm2) - 1.5*
     &ds2(im1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(im1,jm2)*ds2
     &(i,jm2)*ds2(i,jm1))

      X(-1,-1)=(0.5*sig12(im1,jm1))/(ds1(im1,jm1)*ds2(i,jm1)) - (1.0*s
     &igB(4,im1,3))/(dh1(i,j)*(0.5*ds2(im1,jm2) - 1.5*ds2(im1,jm1))) +
     & (2.0*gam22(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) - (2.0
     &*e6e*gam21(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) + (2.0*
     &gam21(i,jm1))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm2)*ds2(i,jm1))

      X(-1,0)=(1.0*(((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)) + (gam11(i
     &m1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(im1,j))/ds1(im1,j) 
     &+ (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j) + s
     &igB(4,im1,1)/(dh1(i,j)*ds1(im1,j)) + (0.5*sig12(im1,jm1))/(ds1(i
     &m1,j)*ds2(i,jm1)) + (0.75*sigB(4,im1,3))/(dh1(i,j)*(0.5*ds2(im1,
     &jm2) - 1.5*ds2(im1,jm1))) + (gam22(i,j)*(e4coefint - 2.0))/(dh1(
     &i,j)*ds1(im1,j)*ds2(i,jm1)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*
     &ds1(im1,j)*ds2(i,jm1)**2)

      X(0,-3)=-(1.0*gam22(i,jm1))/(dh2(i,jm2)*ds2(i,jm3)*ds2(i,jm2)*ds
     &2(i,jm1))

      X(0,-2)=((0.25*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) 
     &- (0.25*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(i,j)
     & + gam22(i,j)/(dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) + (1.0*gam22
     &(i,jm1)*((1.0/ds2(i,jm3) + 1.0/ds2(i,jm2))/dh2(i,jm2) + 1.0/(dh2
     &(i,jm1)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) + (2.0*gam21(i,jm1
     &)*(1.0/ds1(im1,jm2) + 1.0/ds1(i,jm2)))/(dh1(i,jm2)*ds2(i,jm2)*ds
     &2(i,jm1))

      X(0,-1)=(2.0*sig22(i,jm1))/ds2(i,jm1)**2 - ((1.0*sigB(4,im1,3))/
     &(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) - (1.0*sigB(4,i,3))/(0.5*ds2(i
     &,jm2) - 1.5*ds2(i,jm1)))/dh1(i,j) - (2.0*((0.25*sig12(im1,jm1))/
     &ds1(im1,jm1) - (0.25*sig12(i,jm1))/ds1(i,jm1)))/ds2(i,jm1) - (ga
     &m22(i,j)*((2.0*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/dh1(i,jm1) +
     & (1.0/ds2(i,jm2) + 1.0/ds2(i,jm1))/dh2(i,jm1)))/ds2(i,jm1)**2 - 
     &(1.0*gam22(i,jm1)*((1.0*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i
     &,jm1) + 1.0/(dh2(i,jm2)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) - 
     &(2.0*gam21(i,jm1)*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/(dh1(i,jm
     &1)*ds2(i,jm2)*ds2(i,jm1)) + (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,jm1
     &) + 1.0/ds1(i,jm1)))/(dh1(i,jm1)*ds2(i,jm1)**2)

      X(0,0)=((0.75*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) -
     & (0.75*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(i,j) 
     &- (2.0*((0.25*sig12(im1,jm1))/ds1(im1,j) - (0.25*sig12(i,jm1))/d
     &s1(i,j)))/ds2(i,jm1) - (1.0*(((1.0*gam11(im1,j))/(dh1(im1,j)*ds1
     &(im1,j)) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(
     &i,j))/ds1(im1,j) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j
     &)))/(dh1(i,j)*ds1(i,j))))/dh1(i,j) - (2.0*sig22(i,jm1))/ds2(i,jm
     &1)**2 - ((1.0*sigB(4,im1,1))/ds1(im1,j) + (1.0*sigB(4,i,1))/ds1(
     &i,j))/dh1(i,j) + (1.0*gam22(i,j))/(dh2(i,jm1)*ds2(i,jm1)**3) + (
     &1.0*gam22(i,jm1))/(dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) - (gam22
     &(i,j)*(e4coefint - 2.0)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,
     &j)*ds2(i,jm1)**2) - (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,j) + 1.0/ds
     &1(i,j)))/(dh1(i,j)*ds2(i,jm1)**2)

      X(+1,-2)=- (0.25*sigB(4,i,3))/(dh1(i,j)*(0.5*ds2(ip1,jm2) - 1.5*
     &ds2(ip1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(i,jm2)*ds2(i
     &,jm2)*ds2(i,jm1))

      X(+1,-1)=(1.0*sigB(4,i,3))/(dh1(i,j)*(0.5*ds2(ip1,jm2) - 1.5*ds2
     &(ip1,jm1))) - (0.5*sig12(i,jm1))/(ds1(i,jm1)*ds2(i,jm1)) + (2.0*
     &gam22(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) - (2.0*e6e*gam
     &21(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) + (2.0*gam21(i,jm
     &1))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm2)*ds2(i,jm1))

      X(+1,0)=(1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)**2) + (1.0*gam
     &11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j) + (1.0*sigB(4
     &,i,1))/(dh1(i,j)*ds1(i,j)) - (0.5*sig12(i,jm1))/(ds1(i,j)*ds2(i,
     &jm1)) - (0.75*sigB(4,i,3))/(dh1(i,j)*(0.5*ds2(ip1,jm2) - 1.5*ds2
     &(ip1,jm1))) + (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(i,j)*ds1(i
     &,j)*ds2(i,jm1)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*ds1(i,j)*ds2
     &(i,jm1)**2)

      X(-4:4,-4:4)=-X(-4:4,-4:4)
      FlagForceCalculation(i,j)=1
      coefX00=0.0d0
!! ================== BNDY 2,1  
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .eq. 2)) then
!= ===============  2,1 VER2
!CCCCCCCC+X-xCCCCCCCC

      X(-1,0)=gam11(i,j)/(dh1(i,j)**2*ds1(im1,j)**2) + sigB(3,im1,1)/(
     &dh1(i,j)*ds1(im1,j)) - (0.5*sig12(im1,j))/(ds1(im1,j)*ds2(i,j)) 
     &+ (0.75*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.5*ds2(im1
     &,jp1))) - (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(i,j)*ds1(im1,j
     &)*ds2(i,j)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*ds2(i
     &,j)**2)

      X(-1,+1)=(2.0*gam21(i,jp1))/(dh1(i,jp1)*ds1(im1,jp1)*ds2(i,j)*ds
     &2(i,jp1)) - (1.0*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.
     &5*ds2(im1,jp1))) - (2.0*gam22(i,j))/(dh1(i,jp1)*ds1(im1,jp1)*ds2
     &(i,j)**2) - (2.0*e6e*gam21(i,j))/(dh1(i,jp1)*ds1(im1,jp1)*ds2(i,
     &j)**2) - (0.5*sig12(im1,j))/(ds1(im1,jp1)*ds2(i,j))

      X(-1,+2)=(0.25*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.5*
     &ds2(im1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*ds1(im1,jp2)*ds2
     &(i,j)*ds2(i,jp1))

      X(0,0)=(2.0*((0.25*sig12(im1,j))/ds1(im1,j) - (0.25*sig12(i,j))/
     &ds1(i,j)))/ds2(i,j) - ((1.0*sigB(3,im1,1))/ds1(im1,j) + (1.0*sig
     &B(3,i,1))/ds1(i,j))/dh1(i,j) - ((1.0*gam11(ip1,j))/(dh1(ip1,j)*d
     &s1(i,j)**2) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(
     &i,j)*ds1(im1,j)))/dh1(i,j) - (2.0*sig22(i,j))/ds2(i,j)**2 + ((0.
     &75*sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - (0.75*sigB(3
     &,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) - (1.0*gam22(
     &i,j))/(dh2(i,jp1)*ds2(i,j)**3) + gam22(i,jp1)/(dh2(i,jp1)*ds2(i,
     &j)**2*ds2(i,jp1)) + (1.0*gam22(i,j)*(e4coefint - 2.0)*(1.0/ds1(i
     &m1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds2(i,j)**2) - (2.0*e6e*gam21(i
     &,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds2(i,j)**2)

      X(0,+1)=(2.0*((0.25*sig12(im1,j))/ds1(im1,jp1) - (0.25*sig12(i,j
     &))/ds1(i,jp1)))/ds2(i,j) - ((1.0*sigB(3,im1,3))/(1.5*ds2(i,j) - 
     &0.5*ds2(i,jp1)) - (1.0*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1)))/dh1(im1,j) + (2.0*sig22(i,j))/ds2(i,j)**2 + (1.0*gam22(i,j)
     &*((2.0*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/dh1(i,jp1) + (1.0/ds
     &2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1)))/ds2(i,j)**2 - (gam22(i,jp1
     &)*((1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1) + 1.0/(dh2(i,jp2)*
     &ds2(i,jp1))))/(ds2(i,j)*ds2(i,jp1)) - (2.0*gam21(i,jp1)*(1.0/ds1
     &(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,jp1)*ds2(i,j)*ds2(i,jp1)) + 
     &(2.0*e6e*gam21(i,j)*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,
     &jp1)*ds2(i,j)**2)

      X(0,+2)=((0.25*sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 
     &(0.25*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) -
     & (1.0*gam22(i,j))/(dh2(i,jp1)*ds2(i,j)**2*ds2(i,jp1)) + (gam22(i
     &,jp1)*((1.0*(1.0/ds2(i,jp1) + 1.0/ds2(i,jp2)))/dh2(i,jp2) + 1.0/
     &(dh2(i,jp1)*ds2(i,jp1))))/(ds2(i,j)*ds2(i,jp1)) + (2.0*gam21(i,j
     &p1)*(1.0/ds1(im1,jp2) + 1.0/ds1(i,jp2)))/(dh1(i,jp2)*ds2(i,j)*ds
     &2(i,jp1))

      X(0,+3)=-(1.0*gam22(i,jp1))/(dh2(i,jp2)*ds2(i,j)*ds2(i,jp1)*ds2(
     &i,jp2))

      X(+1,0)=((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)) + (1.0*
     &gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/(dh1(ip1,j)*ds1(i,
     &j)))/dh1(i,j) + (1.0*sigB(3,i,1))/(dh1(i,j)*ds1(i,j)) + (0.5*sig
     &12(i,j))/(ds1(i,j)*ds2(i,j)) - (0.75*sigB(3,i,3))/(dh1(im1,j)*(1
     &.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (1.0*gam22(i,j)*(e4coefint 
     &- 2.0))/(dh1(i,j)*ds1(i,j)*ds2(i,j)**2) + (2.0*e6e*gam21(i,j))/(
     &dh1(i,j)*ds1(i,j)*ds2(i,j)**2)

      X(+1,+1)=(0.5*sig12(i,j))/(ds1(i,jp1)*ds2(i,j)) + (1.0*sigB(3,i,
     &3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (2.0*gam2
     &2(i,j))/(dh1(i,jp1)*ds1(i,jp1)*ds2(i,j)**2) - (2.0*e6e*gam21(i,j
     &))/(dh1(i,jp1)*ds1(i,jp1)*ds2(i,j)**2) + (2.0*gam21(i,jp1))/(dh1
     &(i,jp1)*ds1(i,jp1)*ds2(i,j)*ds2(i,jp1))

      X(+1,+2)=- (0.25*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*
     &ds2(ip1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*ds1(i,jp2)*ds2(i
     &,j)*ds2(i,jp1))

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
     &*sig12(im1,j))/(ds1(im1,j)*ds2(i,j)) + (0.75*sigB(3,im1,3))/(dh1
     &(im1,j)*(1.5*ds2(im1,j) - 0.5*ds2(im1,jp1))) - (1.0*gam22(i,j)*(
     &e4coefint - 2.0))/(dh1(i,j)*ds1(im1,j)*ds2(i,j)**2) + (2.0*e6e*g
     &am21(i,j))/(dh1(i,j)*ds1(im1,j)*ds2(i,j)**2)

      X(-1,+1)=(2.0*gam21(i,jp1))/(dh1(i,jp1)*ds1(im1,jp1)*ds2(i,j)*ds
     &2(i,jp1)) - (1.0*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.
     &5*ds2(im1,jp1))) - (2.0*gam22(i,j))/(dh1(i,jp1)*ds1(im1,jp1)*ds2
     &(i,j)**2) - (2.0*e6e*gam21(i,j))/(dh1(i,jp1)*ds1(im1,jp1)*ds2(i,
     &j)**2) - (0.5*sig12(im1,j))/(ds1(im1,jp1)*ds2(i,j))

      X(-1,+2)=(0.25*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.5*
     &ds2(im1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*ds1(im1,jp2)*ds2
     &(i,j)*ds2(i,jp1))

      X(0,0)=(2.0*((0.25*sig12(im1,j))/ds1(im1,j) - (0.25*sig12(i,j))/
     &ds1(i,j)))/ds2(i,j) - ((1.0*sigB(3,im1,1))/ds1(im1,j) + (1.0*sig
     &B(3,i,1))/ds1(i,j))/dh1(i,j) - (2.0*sig22(i,j))/ds2(i,j)**2 + ((
     &0.75*sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - (0.75*sigB
     &(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) - (1.0*((1.
     &0*gam11(im1,j))/(dh1(im1,j)*ds1(im1,j)**2) + (1.0*gam11(i,j)*(1.
     &0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds1(i,j))))/dh1(i,j) - (
     &1.0*gam22(i,j))/(dh2(i,jp1)*ds2(i,j)**3) + gam22(i,jp1)/(dh2(i,j
     &p1)*ds2(i,j)**2*ds2(i,jp1)) + (1.0*gam22(i,j)*(e4coefint - 2.0)*
     &(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds2(i,j)**2) - (2.0*e
     &6e*gam21(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds2(i,j
     &)**2)

      X(0,+1)=(2.0*((0.25*sig12(im1,j))/ds1(im1,jp1) - (0.25*sig12(i,j
     &))/ds1(i,jp1)))/ds2(i,j) - ((1.0*sigB(3,im1,3))/(1.5*ds2(i,j) - 
     &0.5*ds2(i,jp1)) - (1.0*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1)))/dh1(im1,j) + (2.0*sig22(i,j))/ds2(i,j)**2 + (1.0*gam22(i,j)
     &*((2.0*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/dh1(i,jp1) + (1.0/ds
     &2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1)))/ds2(i,j)**2 - (gam22(i,jp1
     &)*((1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1) + 1.0/(dh2(i,jp2)*
     &ds2(i,jp1))))/(ds2(i,j)*ds2(i,jp1)) - (2.0*gam21(i,jp1)*(1.0/ds1
     &(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,jp1)*ds2(i,j)*ds2(i,jp1)) + 
     &(2.0*e6e*gam21(i,j)*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,
     &jp1)*ds2(i,j)**2)

      X(0,+2)=((0.25*sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 
     &(0.25*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) -
     & (1.0*gam22(i,j))/(dh2(i,jp1)*ds2(i,j)**2*ds2(i,jp1)) + (gam22(i
     &,jp1)*((1.0*(1.0/ds2(i,jp1) + 1.0/ds2(i,jp2)))/dh2(i,jp2) + 1.0/
     &(dh2(i,jp1)*ds2(i,jp1))))/(ds2(i,j)*ds2(i,jp1)) + (2.0*gam21(i,j
     &p1)*(1.0/ds1(im1,jp2) + 1.0/ds1(i,jp2)))/(dh1(i,jp2)*ds2(i,j)*ds
     &2(i,jp1))

      X(0,+3)=-(1.0*gam22(i,jp1))/(dh2(i,jp2)*ds2(i,j)*ds2(i,jp1)*ds2(
     &i,jp2))

      X(+1,0)=(1.0*gam11(i,j))/(dh1(i,j)**2*ds1(i,j)**2) + (1.0*sigB(3
     &,i,1))/(dh1(i,j)*ds1(i,j)) + (0.5*sig12(i,j))/(ds1(i,j)*ds2(i,j)
     &) - (0.75*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1
     &,jp1))) - (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(i,j)*ds1(i,j)*
     &ds2(i,j)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*ds1(i,j)*ds2(i,j)*
     &*2)

      X(+1,+1)=(0.5*sig12(i,j))/(ds1(i,jp1)*ds2(i,j)) + (1.0*sigB(3,i,
     &3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (2.0*gam2
     &2(i,j))/(dh1(i,jp1)*ds1(i,jp1)*ds2(i,j)**2) - (2.0*e6e*gam21(i,j
     &))/(dh1(i,jp1)*ds1(i,jp1)*ds2(i,j)**2) + (2.0*gam21(i,jp1))/(dh1
     &(i,jp1)*ds1(i,jp1)*ds2(i,j)*ds2(i,jp1))

      X(+1,+2)=- (0.25*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*
     &ds2(ip1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*ds1(i,jp2)*ds2(i
     &,j)*ds2(i,jp1))

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
     &+ (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) + s
     &igB(1,jm1,2)/(dh2(i,j)*ds2(i,jm1)) - (0.5*sig12(i,jm1))/(ds1(i,j
     &)*ds2(i,jm1)) + (0.75*sigB(1,jm1,3))/(dh2(i,j)*(1.5*ds1(i,jm1) -
     & 0.5*ds1(ip1,jm1))) - (1.0*gam11(i,j)*(e4coefint - 2.0))/(dh2(i,
     &j)*ds1(i,j)**2*ds2(i,jm1)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*ds1(
     &i,j)**2*ds2(i,jm1))

      X(0,0)=(2.0*((0.25*sig12(i,jm1))/ds2(i,jm1) - (0.25*sig12(i,j))/
     &ds2(i,j)))/ds1(i,j) - ((1.0*sigB(1,jm1,2))/ds2(i,jm1) + (1.0*sig
     &B(1,j,2))/ds2(i,j))/dh2(i,j) - (1.0*(((1.0*gam22(i,jm1))/(dh2(i,
     &jm1)*ds2(i,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j
     &)))/dh2(i,j))/ds2(i,jm1) + (1.0*((1.0*gam22(i,jp1))/(dh2(i,jp1)*
     &ds2(i,j)) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j
     &)))/ds2(i,j)))/dh2(i,j) - (2.0*sig11(i,j))/ds1(i,j)**2 + ((0.75*
     &sigB(1,jm1,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - (0.75*sigB(1,j,
     &3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)))/dh2(i,j) - (1.0*gam11(i,j))
     &/(dh1(ip1,j)*ds1(i,j)**3) + gam11(ip1,j)/(dh1(ip1,j)*ds1(i,j)**2
     &*ds1(ip1,j)) + (1.0*gam11(i,j)*(e4coefint - 2.0)*(1.0/ds2(i,jm1)
     & + 1.0/ds2(i,j)))/(dh2(i,j)*ds1(i,j)**2) - (2.0*e6e*gam21(i,j)*(
     &1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds1(i,j)**2)

      X(0,+1)=(1.0*((1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*
     &gam22(i,jp1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1)))/ds2(i
     &,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j)
     & + (1.0*sigB(1,j,2))/(dh2(i,j)*ds2(i,j)) + (0.5*sig12(i,j))/(ds1
     &(i,j)*ds2(i,j)) - (0.75*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,jp1) -
     & 0.5*ds1(ip1,jp1))) - (1.0*gam11(i,j)*(e4coefint - 2.0))/(dh2(i,
     &j)*ds1(i,j)**2*ds2(i,j)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*ds1(i,
     &j)**2*ds2(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(+1,-1)=(2.0*gam21(ip1,j))/(dh2(ip1,j)*ds1(i,j)*ds1(ip1,j)*ds2(
     &ip1,jm1)) - (1.0*sigB(1,jm1,3))/(dh2(i,j)*(1.5*ds1(i,jm1) - 0.5*
     &ds1(ip1,jm1))) - (2.0*gam11(i,j))/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip
     &1,jm1)) - (2.0*e6e*gam21(i,j))/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip1,j
     &m1)) - (0.5*sig12(i,jm1))/(ds1(i,j)*ds2(ip1,jm1))

      X(+1,0)=(2.0*((0.25*sig12(i,jm1))/ds2(ip1,jm1) - (0.25*sig12(i,j
     &))/ds2(ip1,j)))/ds1(i,j) - ((1.0*sigB(1,jm1,3))/(1.5*ds1(i,j) - 
     &0.5*ds1(ip1,j)) - (1.0*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,
     &j)))/dh2(i,j) + (2.0*sig11(i,j))/ds1(i,j)**2 + (1.0*gam11(i,j)*(
     &(1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + (2.0*(1.0/ds2(ip1,j
     &m1) + 1.0/ds2(ip1,j)))/dh2(ip1,j)))/ds1(i,j)**2 - (gam11(ip1,j)*
     &((1.0/ds1(i,j) + 1.0/ds1(ip1,j))/dh1(ip1,j) + 1.0/(dh1(ip2,j)*ds
     &1(ip1,j))))/(ds1(i,j)*ds1(ip1,j)) - (2.0*gam21(ip1,j)*(1.0/ds2(i
     &p1,jm1) + 1.0/ds2(ip1,j)))/(dh2(ip1,j)*ds1(i,j)*ds1(ip1,j)) + (2
     &.0*e6e*gam21(i,j)*(1.0/ds2(ip1,jm1) + 1.0/ds2(ip1,j)))/(dh2(ip1,
     &j)*ds1(i,j)**2)

      X(+1,+1)=(0.5*sig12(i,j))/(ds1(i,j)*ds2(ip1,j)) + (1.0*sigB(1,j,
     &3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds1(ip1,jp1))) - (2.0*gam11(
     &i,j))/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip1,j)) - (2.0*e6e*gam21(i,j))
     &/(dh2(ip1,j)*ds1(i,j)**2*ds2(ip1,j)) + (2.0*gam21(ip1,j))/(dh2(i
     &p1,j)*ds1(i,j)*ds1(ip1,j)*ds2(ip1,j))

      X(+2,-1)=(0.25*sigB(1,jm1,3))/(dh2(i,j)*(1.5*ds1(i,jm1) - 0.5*ds
     &1(ip1,jm1))) - (2.0*gam21(ip1,j))/(dh2(ip2,j)*ds1(i,j)*ds1(ip1,j
     &)*ds2(ip2,jm1))

      X(+2,0)=((0.25*sigB(1,jm1,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)) - 
     &(0.25*sigB(1,j,3))/(1.5*ds1(i,j) - 0.5*ds1(ip1,j)))/dh2(i,j) - (
     &1.0*gam11(i,j))/(dh1(ip1,j)*ds1(i,j)**2*ds1(ip1,j)) + (gam11(ip1
     &,j)*((1.0*(1.0/ds1(ip1,j) + 1.0/ds1(ip2,j)))/dh1(ip2,j) + 1.0/(d
     &h1(ip1,j)*ds1(ip1,j))))/(ds1(i,j)*ds1(ip1,j)) + (2.0*gam21(ip1,j
     &)*(1.0/ds2(ip2,jm1) + 1.0/ds2(ip2,j)))/(dh2(ip2,j)*ds1(i,j)*ds1(
     &ip1,j))

      X(+2,+1)=- (0.25*sigB(1,j,3))/(dh2(i,j)*(1.5*ds1(i,jp1) - 0.5*ds
     &1(ip1,jp1))) - (2.0*gam21(ip1,j))/(dh2(ip2,j)*ds1(i,j)*ds1(ip1,j
     &)*ds2(ip2,j))

      X(+3,0)=-(1.0*gam11(ip1,j))/(dh1(ip2,j)*ds1(i,j)*ds1(ip1,j)*ds1(
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
     &igB(3,im1,1)/(dh1(i,j)*ds1(im1,j)) - (0.5*sig12(im1,j))/(ds1(im1
     &,j)*ds2(i,j)) + (0.75*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j)
     & - 0.5*ds2(im1,jp1))) - (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(
     &i,j)*ds1(im1,j)*ds2(i,j)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*ds
     &1(im1,j)*ds2(i,j)**2)

      X(-1,+1)=(2.0*gam21(i,jp1))/(dh1(i,jp1)*ds1(im1,jp1)*ds2(i,j)*ds
     &2(i,jp1)) - (1.0*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.
     &5*ds2(im1,jp1))) - (2.0*gam22(i,j))/(dh1(i,jp1)*ds1(im1,jp1)*ds2
     &(i,j)**2) - (2.0*e6e*gam21(i,j))/(dh1(i,jp1)*ds1(im1,jp1)*ds2(i,
     &j)**2) - (0.5*sig12(im1,j))/(ds1(im1,jp1)*ds2(i,j))

      X(-1,+2)=(0.25*sigB(3,im1,3))/(dh1(im1,j)*(1.5*ds2(im1,j) - 0.5*
     &ds2(im1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*ds1(im1,jp2)*ds2
     &(i,j)*ds2(i,jp1))

      X(0,0)=(2.0*((0.25*sig12(im1,j))/ds1(im1,j) - (0.25*sig12(i,j))/
     &ds1(i,j)))/ds2(i,j) - ((1.0*sigB(3,im1,1))/ds1(im1,j) + (1.0*sig
     &B(3,i,1))/ds1(i,j))/dh1(i,j) - (1.0*(((1.0*gam11(im1,j))/(dh1(im
     &1,j)*ds1(im1,j)) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j
     &)))/dh1(i,j))/ds1(im1,j) + (1.0*((1.0*gam11(ip1,j))/(dh1(ip1,j)*
     &ds1(i,j)) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,j
     &)))/ds1(i,j)))/dh1(i,j) - (2.0*sig22(i,j))/ds2(i,j)**2 + ((0.75*
     &sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - (0.75*sigB(3,i,
     &3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) - (1.0*gam22(i,j
     &))/(dh2(i,jp1)*ds2(i,j)**3) + gam22(i,jp1)/(dh2(i,jp1)*ds2(i,j)*
     &*2*ds2(i,jp1)) + (1.0*gam22(i,j)*(e4coefint - 2.0)*(1.0/ds1(im1,
     &j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds2(i,j)**2) - (2.0*e6e*gam21(i,j)
     &*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds2(i,j)**2)

      X(0,+1)=(2.0*((0.25*sig12(im1,j))/ds1(im1,jp1) - (0.25*sig12(i,j
     &))/ds1(i,jp1)))/ds2(i,j) - ((1.0*sigB(3,im1,3))/(1.5*ds2(i,j) - 
     &0.5*ds2(i,jp1)) - (1.0*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp
     &1)))/dh1(im1,j) + (2.0*sig22(i,j))/ds2(i,j)**2 + (1.0*gam22(i,j)
     &*((2.0*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/dh1(i,jp1) + (1.0/ds
     &2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1)))/ds2(i,j)**2 - (gam22(i,jp1
     &)*((1.0/ds2(i,j) + 1.0/ds2(i,jp1))/dh2(i,jp1) + 1.0/(dh2(i,jp2)*
     &ds2(i,jp1))))/(ds2(i,j)*ds2(i,jp1)) - (2.0*gam21(i,jp1)*(1.0/ds1
     &(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,jp1)*ds2(i,j)*ds2(i,jp1)) + 
     &(2.0*e6e*gam21(i,j)*(1.0/ds1(im1,jp1) + 1.0/ds1(i,jp1)))/(dh1(i,
     &jp1)*ds2(i,j)**2)

      X(0,+2)=((0.25*sigB(3,im1,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)) - 
     &(0.25*sigB(3,i,3))/(1.5*ds2(i,j) - 0.5*ds2(i,jp1)))/dh1(im1,j) -
     & (1.0*gam22(i,j))/(dh2(i,jp1)*ds2(i,j)**2*ds2(i,jp1)) + (gam22(i
     &,jp1)*((1.0*(1.0/ds2(i,jp1) + 1.0/ds2(i,jp2)))/dh2(i,jp2) + 1.0/
     &(dh2(i,jp1)*ds2(i,jp1))))/(ds2(i,j)*ds2(i,jp1)) + (2.0*gam21(i,j
     &p1)*(1.0/ds1(im1,jp2) + 1.0/ds1(i,jp2)))/(dh1(i,jp2)*ds2(i,j)*ds
     &2(i,jp1))

      X(0,+3)=-(1.0*gam22(i,jp1))/(dh2(i,jp2)*ds2(i,j)*ds2(i,jp1)*ds2(
     &i,jp2))

      X(+1,0)=(1.0*((1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*
     &gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/dh1(ip1,j)))/ds1(i
     &,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j)
     & + (1.0*sigB(3,i,1))/(dh1(i,j)*ds1(i,j)) + (0.5*sig12(i,j))/(ds1
     &(i,j)*ds2(i,j)) - (0.75*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j)
     & - 0.5*ds2(ip1,jp1))) - (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(
     &i,j)*ds1(i,j)*ds2(i,j)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*ds1(
     &i,j)*ds2(i,j)**2)

      X(+1,+1)=(0.5*sig12(i,j))/(ds1(i,jp1)*ds2(i,j)) + (1.0*sigB(3,i,
     &3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*ds2(ip1,jp1))) - (2.0*gam2
     &2(i,j))/(dh1(i,jp1)*ds1(i,jp1)*ds2(i,j)**2) - (2.0*e6e*gam21(i,j
     &))/(dh1(i,jp1)*ds1(i,jp1)*ds2(i,j)**2) + (2.0*gam21(i,jp1))/(dh1
     &(i,jp1)*ds1(i,jp1)*ds2(i,j)*ds2(i,jp1))

      X(+1,+2)=- (0.25*sigB(3,i,3))/(dh1(im1,j)*(1.5*ds2(ip1,j) - 0.5*
     &ds2(ip1,jp1))) - (2.0*gam21(i,jp1))/(dh1(i,jp2)*ds1(i,jp2)*ds2(i
     &,j)*ds2(i,jp1))

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
       FlagForceCalculation(i,j)=1
       coefX00=0.0d0
!! ================== BNDY N,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin)   
     &        .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .gt. jst) .and. (j.lt. jen)) then
!= ===============  N,j VER2
!CCCCCCCC+X-xCCCCCCCC 
      X(-3,0)=-(1.0*gam11(im1,j))/(dh1(im2,j)*ds1(im3,j)*ds1(im2,j)*ds
     &1(im1,j))

      X(-2,-1)=(0.25*sigB(2,jm1,3))/(dh2(i,j)*(0.5*ds1(im2,jm1) - 1.5*
     &ds1(im1,jm1))) - (2.0*gam21(im1,j))/(dh2(im2,j)*ds1(im2,j)*ds1(i
     &m1,j)*ds2(im2,jm1))

      X(-2,0)=((0.25*sigB(2,jm1,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) 
     &- (0.25*sigB(2,j,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))/dh2(i,j)
     & + gam11(i,j)/(dh1(im1,j)*ds1(im2,j)*ds1(im1,j)**2) + (1.0*gam11
     &(im1,j)*((1.0/ds1(im3,j) + 1.0/ds1(im2,j))/dh1(im2,j) + 1.0/(dh1
     &(im1,j)*ds1(im2,j))))/(ds1(im2,j)*ds1(im1,j)) + (2.0*gam21(im1,j
     &)*(1.0/ds2(im2,jm1) + 1.0/ds2(im2,j)))/(dh2(im2,j)*ds1(im2,j)*ds
     &1(im1,j))

      X(-2,+1)=- (0.25*sigB(2,j,3))/(dh2(i,j)*(0.5*ds1(im2,jp1) - 1.5*
     &ds1(im1,jp1))) - (2.0*gam21(im1,j))/(dh2(im2,j)*ds1(im2,j)*ds1(i
     &m1,j)*ds2(im2,j))

      X(-1,-1)=(0.5*sig12(im1,jm1))/(ds1(im1,j)*ds2(im1,jm1)) - (1.0*s
     &igB(2,jm1,3))/(dh2(i,j)*(0.5*ds1(im2,jm1) - 1.5*ds1(im1,jm1))) +
     & (2.0*gam11(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,jm1)) - (2.0
     &*e6e*gam21(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,jm1)) + (2.0*
     &gam21(im1,j))/(dh2(im1,j)*ds1(im2,j)*ds1(im1,j)*ds2(im1,jm1))

      X(-1,0)=(2.0*sig11(im1,j))/ds1(im1,j)**2 - ((1.0*sigB(2,jm1,3))/
     &(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) - (1.0*sigB(2,j,3))/(0.5*ds1(i
     &m2,j) - 1.5*ds1(im1,j)))/dh2(i,j) - (2.0*((0.25*sig12(im1,jm1))/
     &ds2(im1,jm1) - (0.25*sig12(im1,j))/ds2(im1,j)))/ds1(im1,j) - (ga
     &m11(i,j)*((1.0/ds1(im2,j) + 1.0/ds1(im1,j))/dh1(im1,j) + (2.0*(1
     &.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/dh2(im1,j)))/ds1(im1,j)**2 - 
     &(1.0*gam11(im1,j)*((1.0*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(i
     &m1,j) + 1.0/(dh1(im2,j)*ds1(im2,j))))/(ds1(im2,j)*ds1(im1,j)) - 
     &(2.0*gam21(im1,j)*(1.0/ds2(im1,jm1) + 1.0/ds2(im1,j)))/(dh2(im1,
     &j)*ds1(im2,j)*ds1(im1,j)) + (2.0*e6e*gam21(i,j)*(1.0/ds2(im1,jm1
     &) + 1.0/ds2(im1,j)))/(dh2(im1,j)*ds1(im1,j)**2)

      X(-1,+1)=(1.0*sigB(2,j,3))/(dh2(i,j)*(0.5*ds1(im2,jp1) - 1.5*ds1
     &(im1,jp1))) - (0.5*sig12(im1,j))/(ds1(im1,j)*ds2(im1,j)) + (2.0*
     &gam11(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,j)) - (2.0*e6e*gam
     &21(i,j))/(dh2(im1,j)*ds1(im1,j)**2*ds2(im1,j)) + (2.0*gam21(im1,
     &j))/(dh2(im1,j)*ds1(im2,j)*ds1(im1,j)*ds2(im1,j))

      X(0,-2)=-(1.0*gam22(i,jm1))/(dh2(i,jm1)*dh2(i,j)*ds2(i,jm2)*ds2(
     &i,jm1))

      X(0,-1)=(1.0*(((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)) + (gam22(i
     &,jm1)*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i,jm1))/ds2(i,jm1) 
     &+ (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j) + s
     &igB(2,jm1,2)/(dh2(i,j)*ds2(i,jm1)) + (0.5*sig12(im1,jm1))/(ds1(i
     &m1,j)*ds2(i,jm1)) + (0.75*sigB(2,jm1,3))/(dh2(i,j)*(0.5*ds1(im2,
     &jm1) - 1.5*ds1(im1,jm1))) + (gam11(i,j)*(e4coefint - 2.0))/(dh2(
     &i,j)*ds1(im1,j)**2*ds2(i,jm1)) + (2.0*e6e*gam21(i,j))/(dh2(i,j)*
     &ds1(im1,j)**2*ds2(i,jm1))

      X(0,0)=((0.75*sigB(2,jm1,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)) -
     & (0.75*sigB(2,j,3))/(0.5*ds1(im2,j) - 1.5*ds1(im1,j)))/dh2(i,j) 
     &- (2.0*((0.25*sig12(im1,jm1))/ds2(i,jm1) - (0.25*sig12(im1,j))/d
     &s2(i,j)))/ds1(im1,j) - (1.0*(((1.0*gam22(i,jm1))/(dh2(i,jm1)*ds2
     &(i,jm1)) + (1.0*gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(
     &i,j))/ds2(i,jm1) + (1.0*((1.0*gam22(i,jp1))/(dh2(i,jp1)*ds2(i,j)
     &) + (gam22(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/dh2(i,j)))/ds2(
     &i,j)))/dh2(i,j) - (2.0*sig11(im1,j))/ds1(im1,j)**2 - ((1.0*sigB(
     &2,jm1,2))/ds2(i,jm1) + (1.0*sigB(2,j,2))/ds2(i,j))/dh2(i,j) + (1
     &.0*gam11(i,j))/(dh1(im1,j)*ds1(im1,j)**3) + (1.0*gam11(im1,j))/(
     &dh1(im1,j)*ds1(im2,j)*ds1(im1,j)**2) - (gam11(i,j)*(e4coefint - 
     &2.0)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*ds1(im1,j)**2) -
     & (2.0*e6e*gam21(i,j)*(1.0/ds2(i,jm1) + 1.0/ds2(i,j)))/(dh2(i,j)*
     &ds1(im1,j)**2)

      X(0,+1)=(1.0*((1.0*((1.0*gam22(i,j))/(dh2(i,j)*ds2(i,j)) + (1.0*
     &gam22(i,jp1)*(1.0/ds2(i,j) + 1.0/ds2(i,jp1)))/dh2(i,jp1)))/ds2(i
     &,j) + (1.0*gam22(i,j))/(dh2(i,j)*ds2(i,jm1)*ds2(i,j))))/dh2(i,j)
     & + (1.0*sigB(2,j,2))/(dh2(i,j)*ds2(i,j)) - (0.5*sig12(im1,j))/(d
     &s1(im1,j)*ds2(i,j)) - (0.75*sigB(2,j,3))/(dh2(i,j)*(0.5*ds1(im2,
     &jp1) - 1.5*ds1(im1,jp1))) + (1.0*gam11(i,j)*(e4coefint - 2.0))/(
     &dh2(i,j)*ds1(im1,j)**2*ds2(i,j)) + (2.0*e6e*gam21(i,j))/(dh2(i,j
     &)*ds1(im1,j)**2*ds2(i,j))

      X(0,+2)=-(1.0*gam22(i,jp1))/(dh2(i,j)*dh2(i,jp1)*ds2(i,j)*ds2(i,
     &jp1))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
      FlagForceCalculation(i,j)=1
      coefX00=0.0d0
!! ================== BNDY i,N 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin)   
     &        .and. (BoundaryFibShellin(1,2) .eq. 1)  
     &    .and. (i .gt. ist) .and. (i.lt. ien)) then
!= ===============  i,N VER2
!CCCCCCCC+X-xCCCCCCCC
      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,-2)=(0.25*sigB(4,im1,3))/(dh1(i,j)*(0.5*ds2(im1,jm2) - 1.5*
     &ds2(im1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(im1,jm2)*ds2
     &(i,jm2)*ds2(i,jm1))

      X(-1,-1)=(0.5*sig12(im1,jm1))/(ds1(im1,jm1)*ds2(i,jm1)) - (1.0*s
     &igB(4,im1,3))/(dh1(i,j)*(0.5*ds2(im1,jm2) - 1.5*ds2(im1,jm1))) +
     & (2.0*gam22(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) - (2.0
     &*e6e*gam21(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) + (2.0*
     &gam21(i,jm1))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm2)*ds2(i,jm1))

      X(-1,0)=(1.0*(((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)) + (gam11(i
     &m1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/dh1(im1,j))/ds1(im1,j) 
     &+ (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j) + s
     &igB(4,im1,1)/(dh1(i,j)*ds1(im1,j)) + (0.5*sig12(im1,jm1))/(ds1(i
     &m1,j)*ds2(i,jm1)) + (0.75*sigB(4,im1,3))/(dh1(i,j)*(0.5*ds2(im1,
     &jm2) - 1.5*ds2(im1,jm1))) + (gam22(i,j)*(e4coefint - 2.0))/(dh1(
     &i,j)*ds1(im1,j)*ds2(i,jm1)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*
     &ds1(im1,j)*ds2(i,jm1)**2)

      X(0,-3)=-(1.0*gam22(i,jm1))/(dh2(i,jm2)*ds2(i,jm3)*ds2(i,jm2)*ds
     &2(i,jm1))

      X(0,-2)=((0.25*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) 
     &- (0.25*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(i,j)
     & + gam22(i,j)/(dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) + (1.0*gam22
     &(i,jm1)*((1.0/ds2(i,jm3) + 1.0/ds2(i,jm2))/dh2(i,jm2) + 1.0/(dh2
     &(i,jm1)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) + (2.0*gam21(i,jm1
     &)*(1.0/ds1(im1,jm2) + 1.0/ds1(i,jm2)))/(dh1(i,jm2)*ds2(i,jm2)*ds
     &2(i,jm1))

      X(0,-1)=(2.0*sig22(i,jm1))/ds2(i,jm1)**2 - ((1.0*sigB(4,im1,3))/
     &(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) - (1.0*sigB(4,i,3))/(0.5*ds2(i
     &,jm2) - 1.5*ds2(i,jm1)))/dh1(i,j) - (2.0*((0.25*sig12(im1,jm1))/
     &ds1(im1,jm1) - (0.25*sig12(i,jm1))/ds1(i,jm1)))/ds2(i,jm1) - (ga
     &m22(i,j)*((2.0*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/dh1(i,jm1) +
     & (1.0/ds2(i,jm2) + 1.0/ds2(i,jm1))/dh2(i,jm1)))/ds2(i,jm1)**2 - 
     &(1.0*gam22(i,jm1)*((1.0*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2(i
     &,jm1) + 1.0/(dh2(i,jm2)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) - 
     &(2.0*gam21(i,jm1)*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/(dh1(i,jm
     &1)*ds2(i,jm2)*ds2(i,jm1)) + (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,jm1
     &) + 1.0/ds1(i,jm1)))/(dh1(i,jm1)*ds2(i,jm1)**2)

      X(0,0)=((0.75*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) -
     & (0.75*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(i,j) 
     &- (2.0*((0.25*sig12(im1,jm1))/ds1(im1,j) - (0.25*sig12(i,jm1))/d
     &s1(i,j)))/ds2(i,jm1) - (1.0*(((1.0*gam11(im1,j))/(dh1(im1,j)*ds1
     &(im1,j)) + (1.0*gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(
     &i,j))/ds1(im1,j) + (1.0*((1.0*gam11(ip1,j))/(dh1(ip1,j)*ds1(i,j)
     &) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/dh1(i,j)))/ds1(
     &i,j)))/dh1(i,j) - (2.0*sig22(i,jm1))/ds2(i,jm1)**2 - ((1.0*sigB(
     &4,im1,1))/ds1(im1,j) + (1.0*sigB(4,i,1))/ds1(i,j))/dh1(i,j) + (1
     &.0*gam22(i,j))/(dh2(i,jm1)*ds2(i,jm1)**3) + (1.0*gam22(i,jm1))/(
     &dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) - (gam22(i,j)*(e4coefint - 
     &2.0)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds2(i,jm1)**2) -
     & (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*
     &ds2(i,jm1)**2)

      X(+1,-2)=- (0.25*sigB(4,i,3))/(dh1(i,j)*(0.5*ds2(ip1,jm2) - 1.5*
     &ds2(ip1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(i,jm2)*ds2(i
     &,jm2)*ds2(i,jm1))

      X(+1,-1)=(1.0*sigB(4,i,3))/(dh1(i,j)*(0.5*ds2(ip1,jm2) - 1.5*ds2
     &(ip1,jm1))) - (0.5*sig12(i,jm1))/(ds1(i,jm1)*ds2(i,jm1)) + (2.0*
     &gam22(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) - (2.0*e6e*gam
     &21(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) + (2.0*gam21(i,jm
     &1))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm2)*ds2(i,jm1))


      X(+1,0)=(1.0*((1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(i,j)) + (1.0*
     &gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/dh1(ip1,j)))/ds1(i
     &,j) + (1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j))))/dh1(i,j)
     & + (1.0*sigB(4,i,1))/(dh1(i,j)*ds1(i,j)) - (0.5*sig12(i,jm1))/(d
     &s1(i,j)*ds2(i,jm1)) - (0.75*sigB(4,i,3))/(dh1(i,j)*(0.5*ds2(ip1,
     &jm2) - 1.5*ds2(ip1,jm1))) + (1.0*gam22(i,j)*(e4coefint - 2.0))/(
     &dh1(i,j)*ds1(i,j)*ds2(i,jm1)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j
     &)*ds1(i,j)*ds2(i,jm1)**2)

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))

      X(-4:4,-4:4)=-X(-4:4,-4:4)
      FlagForceCalculation(i,j)=1
      coefX00=0.0d0
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
!     X(i2,j2)=X(i2,j2)  !*alpha_coefK
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
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
        enddo
!! ================== BNDY Nx,1
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   )then
        do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
        enddo

!! ================== BNDY 1,Ny
       elseif(((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
        do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
        enddo

!! ================== BNDY Nx,Ny
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
        do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
        enddo
!! ================== BNDY 1,j  2<=j<=Ny-1
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .gt. jst) .and. (j.lt. jen)) then

       do i3=1,3
        P(0:2,-1:1)=GX(i3,i:ip2,jm1:jp1)
        do j3=1,3
          Q(0:2,-1:1)=GX(j3,i:ip2,jm1:jp1)
          call KnonlinearBNDY1j(A0(1:9),iA0(1:9,1:2)
     &         ,P(0:2,-1:1),Q(0:2,-1:1)
     &         ,ds1(i:ip1,jm1:jp1),ds2(i:ip1,jm1:jp1)
     &         ,dh1(i:ip1,jm1:jp1),dh2(i:ip1,jm1:jp1)
     &         ,SigB_d(1,jm1:j,3)
     &         ,Siga_d(i,jm1:j,5)
     &         ,SigB_d(1,jm1:j,2),SigB_d(1,jm1:j,1)
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
        P(-1:1,0:2)=GX(i3,im1:ip1,j:jp2)
        do j3=1,3
          Q(-1:1,0:2)=GX(j3,im1:ip1,j:jp2)
          call KnonlinearBNDY1i(A0(1:9),iA0(1:9,1:2)
     &         ,P(-1:1,0:2),Q(-1:1,0:2)
     &         ,ds1(im1:ip1,j:jp1),ds2(im1:ip1,j:jp1)
     &         ,dh1(im1:ip1,j:jp1),dh2(im1:ip1,j:jp1)
     &         ,SigB_d(3,im1:i,3)
     &         ,Siga_d(im1:i,j,5)
     &         ,SigB_d(3,im1:i,1),SigB_d(3,im1:i,2)
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
     &    .and. (j .gt. jst) .and. (j.lt. jen)) then

       do i3=1,3
        P(-2:0,-1:1)=GX(i3,im2:i,jm1:jp1)
        do j3=1,3
          Q(-2:0,-1:1)=GX(j3,im2:i,jm1:jp1)
          call KnonlinearBNDYNj(A0(1:9),iA0(1:9,1:2)
     &         ,P(-2:0,-1:1),Q(-2:0,-1:1)
     &         ,ds1(im2:im1,jm1:jp1),ds2(im1:i,jm1:jp1)
     &         ,dh1(im1:i,jm1:jp1),dh2(im1:i,jm1:jp1)
     &         ,SigB_d(2,jm1:j,3)
     &         ,Siga_d(i-1,jm1:j,5)
     &         ,SigB_d(2,jm1:j,2),SigB_d(2,jm1:j,1)
     &         ,Siga_d(i-1,j,1),Siga_d(i-1,j,2) )
           !A0=0.0
           do i4=1,9
                Y(i3,j3,iA0(i4,1),iA0(i4,2))=-A0(i4)
           enddo
         enddo
        enddo
        do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)
        enddo

!! ================== BNDY i,N 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin)   
     &        .and. (BoundaryFibShellin(1,2) .eq. 1)  
     &    .and. (i .gt. ist) .and. (i.lt. ien)) then

       do i3=1,3
        P(-1:1,-2:0)=GX(i3,im1:ip1,jm2:j)
        do j3=1,3
          Q(-1:1,-2:0)=GX(j3,im1:ip1,jm2:j)
          call KnonlinearBNDYNi(A0(1:9),iA0(1:9,1:2)
     &         ,P(-1:1,-2:0),Q(-1:1,-2:0)
     &         ,ds1(im1:ip1,jm1:j),ds2(im1:ip1,jm2:jm1)
     &         ,dh1(im1:ip1,jm1:j),dh2(im1:ip1,jm1:j)
     &         ,SigB_d(4,im1:i,3)
     &         ,Siga_d(im1:i,j-1,5)
     &         ,SigB_d(4,im1:i,1),SigB_d(4,im1:i,2)
     &         ,Siga_d(i,j-1,2),Siga_d(i,j-1,1) )
           !A0=0.0
           do i4=1,9
                Y(i3,j3,iA0(i4,1),iA0(i4,2))=-A0(i4)
           enddo
         enddo
        enddo
        do i3=1,3
                 Y(i3,i3,-4:4,-4:4)=
     &           Y(i3,i3,-4:4,-4:4)+X(-4:4,-4:4)

        enddo

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
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)
!     x+1/2,y
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
     &             (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)
!     x+1/2,y
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
     &             (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1)
!     x,y+1/2


              DY_i_jp(1:3)=
     &             (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1)
!     x,y+1/2

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
     &             (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j)
!     x,y+1/2

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
     &  ,dBsig21(-1:0),dsig12(-1:0)
     &  ,dBsig22A(-1:0),dBsig22B(-1:0)
     &  ,dsig11A,dsig11B
      t(2) = 1.0D0/ds2(0,-1)
      t(3) = 1.0D0/dh2(0,0)
      t(4) = ds1(0,-1)*1.5D0
      t(8) = ds1(+1,-1)*(-5.0D-1)
      t(5) = -t(8)
      t(6) = t(4)-t(5)
      t(7) = 1.0D0/t(6)
      t(9) = ds1(0,0)*1.5D0
      t(30) = ds1(+1,0)*(-5.0D-1)
      t(10) = -t(30)
      t(11) = t(9)-t(10)
      t(12) = 1.0D0/t(11)
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
      t(24) = t(4)+t(8)
      t(25) = 1.0D0/t(24)
      t(26) = P(0,-1)*1.5D0
      t(27) = P(+1,-1)*(-2.0D0)
      t(28) = P(+2,-1)*5.0D-1
      t(29) = t(26)+t(27)+t(28)
      t(31) = P(0,0)*1.5D0
      t(32) = P(+1,0)*(-2.0D0)
      t(33) = P(+2,0)*5.0D-1
      t(34) = t(31)+t(32)+t(33)
      t(35) = Q(0,0)*(-1.0D0)
      t(36) = Q(0,-1)+t(35)
      t(37) = 1.0D0/ds1(0,0)**2
      t(38) = P(0,0)+t(15)
      t(39) = 1.0D0/ds2(0,0)
      t(40) = t(9)+t(30)
      t(41) = 1.0D0/t(40)
      t(42) = t(2)*t(36)
      t(43) = Q(+1,-1)+t(23)
      t(44) = t(19)*t(43)
      t(45) = t(42)+t(44)
      t(46) = P(0,+1)*(-1.0D0)
      t(47) = P(0,0)+t(46)
      t(48) = 1.0D0/ds2(+1,0)
      t(49) = t(25)*t(29)*5.0D-1
      t(50) = t(34)*t(41)*5.0D-1
      t(51) = t(49)+t(50)
      t(52) = Q(0,-1)*1.5D0
      t(53) = Q(+1,-1)*(-2.0D0)
      t(54) = Q(+2,-1)*5.0D-1
      t(55) = t(52)+t(53)+t(54)
      t(56) = t(25)*t(55)*5.0D-1
      t(57) = Q(0,0)*1.5D0
      t(58) = Q(+1,0)*(-2.0D0)
      t(59) = Q(+2,0)*5.0D-1
      t(60) = t(57)+t(58)+t(59)
      t(61) = t(41)*t(60)*5.0D-1
      t(62) = ds1(0,+1)*1.5D0
      t(78) = ds1(+1,+1)*(-5.0D-1)
      t(63) = -t(78)
      t(64) = t(62)-t(63)
      t(65) = 1.0D0/t(64)
      t(66) = Q(+1,-1)*(-1.0D0)
      t(67) = Q(0,-1)+t(66)
      t(68) = t(17)*t(67)
      t(69) = Q(0,0)+t(23)
      t(70) = t(16)*t(69)
      t(71) = t(68)+t(70)
      t(72) = t(39)*t(47)*2.5D-1
      t(73) = P(+1,+1)*(-1.0D0)
      t(74) = P(+1,0)+t(73)
      t(75) = t(48)*t(74)*2.5D-1
      t(76) = t(72)+t(75)
      t(77) = Q(+1,+1)*(-1.0D0)
      t(79) = P(0,+1)*1.5D0
      t(80) = P(+1,+1)*(-2.0D0)
      t(81) = P(+2,+1)*5.0D-1
      t(82) = t(79)+t(80)+t(81)
      t(83) = Q(0,+1)*(-1.0D0)
      t(84) = Q(0,0)+t(83)
      t(85) = t(62)+t(78)
      t(86) = 1.0D0/t(85)
      t(87) = t(82)*t(86)*5.0D-1
      t(88) = t(50)+t(87)
      t(89) = Q(0,+1)*1.5D0
      t(90) = Q(+1,+1)*(-2.0D0)
      t(91) = Q(+2,+1)*5.0D-1
      t(92) = t(89)+t(90)+t(91)
      t(93) = 1.0D0/ds1(0,+1)
      t(94) = Q(0,+1)+t(77)
      t(95) = t(93)*t(94)
      t(96) = t(70)+t(95)
      t(97) = t(39)*t(84)
      t(98) = Q(+1,0)+t(77)
      t(99) = t(48)*t(98)
      t(100) = t(97)+t(99)

      A0(1) = dsig11B*t(2)*t(37)*t(38)*5.0D-1-t(2)*t(3)*t(14)*(d
     &Bsig22A(-1)*t(2)+dBsig22B(-1)*t(7)*7.5D-1)*1.0D0+dBsig21(-1)*t(2
     &)*t(3)*(t(56)+t(12)*t(60)*5.0D-1)*(t(7)*t(29)*5.0D-1+t(12)*t(34)
     &*5.0D-1)-dsig12(-1)*t(16)*t(17)*t(22)*t(45)*5.0D-1-dsig12(-1)*t(
     &2)*t(16)*t(22)*t(71)*5.0D-1+dBsig21(-1)*t(2)*t(3)*t(25)*t(36)*t(
     &51)*7.5D-1
      A0(2) = t(37)*t(38)*(dsig11A*t(16)*1.0D0-dsig11B*(t(2)*2.5
     &D-1-t(39)*2.5D-1))*2.0D0-dsig12(-1)*t(22)*t(37)*t(45)*5.0D-1-dsi
     &g12(0)*t(37)*t(76)*t(100)*5.0D-1+t(2)*t(3)*t(14)*(dBsig22A(-1)*t
     &(2)*1.0D0-dBsig22B(-1)*t(41)*7.5D-1)*1.0D0+t(3)*t(39)*t(47)*(dBs
     &ig22A(0)*t(39)+dBsig22B(0)*t(41)*7.5D-1)-dBsig21(0)*t(3)*t(39)*(
     &t(50)+t(65)*t(82)*5.0D-1)*(t(61)+t(65)*t(92)*5.0D-1)*1.0D0-dBsig
     &21(-1)*t(2)*t(3)*t(51)*(t(56)+t(61))*1.0D0+dsig12(-1)*t(2)*t(16)
     &*t(22)*t(71)*5.0D-1-dsig12(0)*t(16)*t(39)*t(76)*t(96)*5.0D-1+dBs
     &ig21(-1)*t(2)*t(3)*t(36)*t(41)*t(51)*7.5D-1-dBsig21(0)*t(3)*t(39
     &)*t(41)*t(84)*t(88)*7.5D-1
      A0(3) = dsig11B*t(37)*t(38)*t(39)*(-5.0D-1)-t(3)*t(39)*t(4
     &7)*(dBsig22A(0)*t(39)*1.0D0-dBsig22B(0)*t(86)*7.5D-1)+dBsig21(0)
     &*t(3)*t(39)*t(88)*(t(61)+t(86)*t(92)*5.0D-1)+dsig12(0)*t(16)*t(3
     &9)*t(76)*t(96)*5.0D-1-dsig12(0)*t(16)*t(76)*t(93)*t(100)*5.0D-1-
     &dBsig21(0)*t(3)*t(39)*t(84)*t(86)*t(88)*7.5D-1
      A0(4) = dsig11B*t(19)*t(37)*t(38)*5.0D-1+dBsig22B(-1)*t(2)
     &*t(3)*t(14)*t(25)+dsig12(-1)*t(16)*t(17)*t(22)*t(45)*5.0D-1-dsig
     &12(-1)*t(16)*t(19)*t(22)*t(71)*5.0D-1-dBsig21(-1)*t(2)*t(3)*t(25
     &)*t(36)*t(51)*1.0D0
      A0(5) = t(37)*t(38)*(dsig11A*t(16)+dsig11B*(t(19)*2.5D-1-t
     &(48)*2.5D-1))*(-2.0D0)+dsig12(-1)*t(22)*t(37)*t(45)*5.0D-1+dsig1
     &2(0)*t(37)*t(76)*t(100)*5.0D-1+dBsig22B(-1)*t(2)*t(3)*t(14)*t(41
     &)-dBsig22B(0)*t(3)*t(39)*t(41)*t(47)*1.0D0+dsig12(-1)*t(16)*t(19
     &)*t(22)*t(71)*5.0D-1-dsig12(0)*t(16)*t(48)*t(76)*t(96)*5.0D-1-dB
     &sig21(-1)*t(2)*t(3)*t(36)*t(41)*t(51)*1.0D0+dBsig21(0)*t(3)*t(39
     &)*t(41)*t(84)*t(88)
      A0(6) = dsig11B*t(37)*t(38)*t(48)*(-5.0D-1)-dBsig22B(0)*t(
     &3)*t(39)*t(47)*t(86)*1.0D0+dsig12(0)*t(16)*t(48)*t(76)*t(96)*5.0
     &D-1+dsig12(0)*t(16)*t(76)*t(93)*t(100)*5.0D-1+dBsig21(0)*t(3)*t(
     &39)*t(84)*t(86)*t(88)
      A0(7) = dBsig22B(-1)*t(2)*t(3)*t(14)*t(25)*(-2.5D-1)+dBsig
     &21(-1)*t(2)*t(3)*t(25)*t(36)*t(51)*2.5D-1
      A0(8) = dBsig22B(-1)*t(2)*t(3)*t(14)*t(41)*(-2.5D-1)+dBsig
     &22B(0)*t(3)*t(39)*t(41)*t(47)*2.5D-1+dBsig21(-1)*t(2)*t(3)*t(36)
     &*t(41)*t(51)*2.5D-1-dBsig21(0)*t(3)*t(39)*t(41)*t(84)*t(88)*2.5D
     &-1
      A0(9) = dBsig22B(0)*t(3)*t(39)*t(47)*t(86)*2.5D-1-dBsig21(
     &0)*t(3)*t(39)*t(84)*t(86)*t(88)*2.5D-1


      iA0(1,1)=0; iA0(1,2)=-1
      iA0(2,1)=0; iA0(2,2)=0
      iA0(3,1)=0; iA0(3,2)=1
      iA0(4,1)=1; iA0(4,2)=0
      iA0(5,1)=1; iA0(5,2)=-1
      iA0(6,1)=1; iA0(6,2)=1
      iA0(7,1)=2; iA0(7,2)=-1
      iA0(8,1)=2; iA0(8,2)=0
      iA0(9,1)=2; iA0(9,2)=1
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
     &  ,dBsig21(-1:0),dsig12(-1:0)
     &  ,dBsig11A(-1:0),dBsig11B(-1:0)
     &  ,dsig22A,dsig22B

      t(2) = 1.0D0/ds1(-1,0)
      t(3) = ds2(-1,0)*1.5D0
      t(7) = ds2(-1,+1)*(-5.0D-1)
      t(4) = -t(7)
      t(5) = t(3)-t(4)
      t(6) = 1.0D0/t(5)
      t(8) = ds2(0,0)*1.5D0
      t(30) = ds2(0,+1)*(-5.0D-1)
      t(9) = -t(30)
      t(10) = t(8)-t(9)
      t(11) = 1.0D0/t(10)
      t(12) = P(0,0)*(-1.0D0)
      t(13) = P(-1,0)+t(12)
      t(14) = P(0,+1)*(-1.0D0)
      t(15) = 1.0D0/ds1(-1,+1)
      t(16) = 1.0D0/ds2(0,0)
      t(17) = t(2)*t(13)*2.5D-1
      t(18) = P(-1,+1)+t(14)
      t(19) = t(15)*t(18)*2.5D-1
      t(20) = t(17)+t(19)
      t(21) = 1.0D0/ds2(-1,0)
      t(22) = Q(0,+1)*(-1.0D0)
      t(23) = 1.0D0/dh1(-1,0)
      t(24) = t(3)+t(7)
      t(25) = 1.0D0/t(24)
      t(26) = P(-1,0)*1.5D0
      t(27) = P(-1,+1)*(-2.0D0)
      t(28) = P(-1,+2)*5.0D-1
      t(29) = t(26)+t(27)+t(28)
      t(31) = P(0,0)*1.5D0
      t(32) = P(0,+1)*(-2.0D0)
      t(33) = P(0,+2)*5.0D-1
      t(34) = t(31)+t(32)+t(33)
      t(35) = Q(0,0)*(-1.0D0)
      t(36) = Q(-1,0)+t(35)
      t(37) = 1.0D0/ds2(0,0)**2
      t(38) = P(0,0)+t(14)
      t(39) = t(2)*t(36)
      t(40) = Q(-1,+1)+t(22)
      t(41) = t(15)*t(40)
      t(42) = t(39)+t(41)
      t(43) = Q(-1,+1)*(-1.0D0)
      t(44) = Q(-1,0)+t(43)
      t(45) = t(21)*t(44)
      t(46) = Q(0,0)+t(22)
      t(47) = t(16)*t(46)
      t(48) = t(45)+t(47)
      t(49) = 1.0D0/dh1(0,0)
      t(50) = t(25)*t(29)*5.0D-1
      t(51) = t(8)+t(30)
      t(52) = 1.0D0/t(51)
      t(53) = t(34)*t(52)*5.0D-1
      t(54) = t(50)+t(53)
      t(55) = 1.0D0/ds1(0,0)
      t(56) = P(+1,0)*(-1.0D0)
      t(57) = P(0,0)+t(56)
      t(58) = 1.0D0/ds1(0,+1)
      t(59) = Q(-1,0)*1.5D0
      t(60) = Q(-1,+1)*(-2.0D0)
      t(61) = Q(-1,+2)*5.0D-1
      t(62) = t(59)+t(60)+t(61)
      t(63) = t(25)*t(62)*5.0D-1
      t(64) = Q(0,0)*1.5D0
      t(65) = Q(0,+1)*(-2.0D0)
      t(66) = Q(0,+2)*5.0D-1
      t(67) = t(64)+t(65)+t(66)
      t(68) = t(52)*t(67)*5.0D-1
      t(69) = ds2(+1,0)*1.5D0
      t(79) = ds2(+1,+1)*(-5.0D-1)
      t(70) = -t(79)
      t(71) = t(69)-t(70)
      t(72) = 1.0D0/t(71)
      t(73) = t(55)*t(57)*2.5D-1
      t(74) = P(+1,+1)*(-1.0D0)
      t(75) = P(0,+1)+t(74)
      t(76) = t(58)*t(75)*2.5D-1
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
      t(88) = t(58)*t(87)
      t(89) = t(86)+t(88)
      t(90) = 1.0D0/ds2(+1,0)
      t(91) = Q(+1,0)+t(78)
      t(92) = t(90)*t(91)
      t(93) = t(47)+t(92)
      t(94) = t(69)+t(79)
      t(95) = 1.0D0/t(94)
      t(96) = t(83)*t(95)*5.0D-1
      t(97) = t(53)+t(96)
      t(98) = Q(+1,0)*1.5D0
      t(99) = Q(+1,+1)*(-2.0D0)
      t(100) = Q(+1,+2)*5.0D-1
      t(101) = t(98)+t(99)+t(100)

      A0(1) = dsig22B*t(2)*t(37)*t(38)*5.0D-1-t(2)*t(13)*t(49)*(
     &dBsig11A(-1)*t(2)+dBsig11B(-1)*t(6)*7.5D-1)*1.0D0+dBsig21(-1)*t(
     &2)*t(23)*(t(63)+t(11)*t(67)*5.0D-1)*(t(6)*t(29)*5.0D-1+t(11)*t(3
     &4)*5.0D-1)-dsig12(-1)*t(2)*t(16)*t(20)*t(48)*5.0D-1-dsig12(-1)*t
     &(16)*t(20)*t(21)*t(42)*5.0D-1+dBsig21(-1)*t(2)*t(23)*t(25)*t(36)
     &*t(54)*7.5D-1
      A0(2) = dsig22B*t(15)*t(37)*t(38)*5.0D-1+dBsig11B(-1)*t(2)
     &*t(13)*t(25)*t(49)-dsig12(-1)*t(15)*t(16)*t(20)*t(48)*5.0D-1+dsi
     &g12(-1)*t(16)*t(20)*t(21)*t(42)*5.0D-1-dBsig21(-1)*t(2)*t(23)*t(
     &25)*t(36)*t(54)*1.0D0
      A0(3) = dBsig11B(-1)*t(2)*t(13)*t(25)*t(49)*(-2.5D-1)+dBsi
     &g21(-1)*t(2)*t(23)*t(25)*t(36)*t(54)*2.5D-1
      A0(4) = t(37)*t(38)*(dsig22A*t(16)*1.0D0-dsig22B*(t(2)*2.5
     &D-1-t(55)*2.5D-1))*2.0D0-dsig12(-1)*t(20)*t(37)*t(42)*5.0D-1-dsi
     &g12(0)*t(37)*t(77)*t(89)*5.0D-1+t(2)*t(13)*t(49)*(dBsig11A(-1)*t
     &(2)*1.0D0-dBsig11B(-1)*t(52)*7.5D-1)*1.0D0+t(49)*t(55)*t(57)*(dB
     &sig11B(0)*t(52)*7.5D-1+dBsig11A(0)*t(55))-dBsig21(0)*t(23)*t(55)
     &*(t(53)+t(72)*t(83)*5.0D-1)*(t(68)+t(72)*t(101)*5.0D-1)*1.0D0-dB
     &sig21(-1)*t(2)*t(23)*t(54)*(t(63)+t(68))*1.0D0+dsig12(-1)*t(2)*t
     &(16)*t(20)*t(48)*5.0D-1-dsig12(0)*t(16)*t(55)*t(77)*t(93)*5.0D-1
     &+dBsig21(-1)*t(2)*t(23)*t(36)*t(52)*t(54)*7.5D-1-dBsig21(0)*t(23
     &)*t(52)*t(55)*t(85)*t(97)*7.5D-1
      A0(5) = t(37)*t(38)*(dsig22A*t(16)+dsig22B*(t(15)*2.5D-1-t
     &(58)*2.5D-1))*(-2.0D0)+dsig12(-1)*t(20)*t(37)*t(42)*5.0D-1+dsig1
     &2(0)*t(37)*t(77)*t(89)*5.0D-1+dBsig11B(-1)*t(2)*t(13)*t(49)*t(52
     &)-dBsig11B(0)*t(49)*t(52)*t(55)*t(57)*1.0D0+dsig12(-1)*t(15)*t(1
     &6)*t(20)*t(48)*5.0D-1-dsig12(0)*t(16)*t(58)*t(77)*t(93)*5.0D-1-d
     &Bsig21(-1)*t(2)*t(23)*t(36)*t(52)*t(54)*1.0D0+dBsig21(0)*t(23)*t
     &(52)*t(55)*t(85)*t(97)
      A0(6) = dBsig11B(-1)*t(2)*t(13)*t(49)*t(52)*(-2.5D-1)+dBsi
     &g11B(0)*t(49)*t(52)*t(55)*t(57)*2.5D-1+dBsig21(-1)*t(2)*t(23)*t(
     &36)*t(52)*t(54)*2.5D-1-dBsig21(0)*t(23)*t(52)*t(55)*t(85)*t(97)*
     &2.5D-1
      A0(7) = dsig22B*t(37)*t(38)*t(55)*(-5.0D-1)-t(49)*t(55)*t(
     &57)*(dBsig11A(0)*t(55)*1.0D0-dBsig11B(0)*t(95)*7.5D-1)+dBsig21(0
     &)*t(23)*t(55)*t(97)*(t(68)+t(95)*t(101)*5.0D-1)+dsig12(0)*t(16)*
     &t(55)*t(77)*t(93)*5.0D-1-dsig12(0)*t(16)*t(77)*t(89)*t(90)*5.0D-
     &1-dBsig21(0)*t(23)*t(55)*t(85)*t(95)*t(97)*7.5D-1
      A0(8) = dsig22B*t(37)*t(38)*t(58)*(-5.0D-1)-dBsig11B(0)*t(
     &49)*t(55)*t(57)*t(95)*1.0D0+dsig12(0)*t(16)*t(58)*t(77)*t(93)*5.
     &0D-1+dsig12(0)*t(16)*t(77)*t(89)*t(90)*5.0D-1+dBsig21(0)*t(23)*t
     &(55)*t(85)*t(95)*t(97)
      A0(9) = dBsig11B(0)*t(49)*t(55)*t(57)*t(95)*2.5D-1-dBsig21
     &(0)*t(23)*t(55)*t(85)*t(95)*t(97)*2.5D-1


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
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine KnonlinearBNDYNj(A0,iA0,P,Q,ds1,ds2,dh1,dh2,
     &   dBsig21,dsig12,dBsig22A,dBsig22B,dsig11A,dsig11B)
      implicit none
      integer iA0(9,2)
      real*8 t(150),A0(9)
      real*8 ds1(-2:-1,-1:1),ds2(-1:0,-1:1)
     &  ,dh2(-1:0,-1:1),dh1(-1:0,-1:1)
     &  ,P(-2:0,-1:1),Q(-2:0,-1:1)
     &  ,dBsig21(-1:0),dsig12(-1:0)
     &  ,dBsig22A(-1:0),dBsig22B(-1:0)
     &  ,dsig11A,dsig11B
      t(2) = 1.0D0/dh2(0,0)
      t(3) = 1.0D0/ds2(0,-1)
      t(4) = ds1(-2,-1)*5.0D-1
      t(8) = ds1(-1,-1)*(-1.5D0)
      t(5) = -t(8)
      t(6) = t(4)-t(5)
      t(7) = 1.0D0/t(6)
      t(9) = ds1(-2,0)*5.0D-1
      t(15) = ds1(-1,0)*(-1.5D0)
      t(10) = -t(15)
      t(11) = t(9)-t(10)
      t(12) = 1.0D0/t(11)
      t(13) = P(0,0)*(-1.0D0)
      t(14) = P(0,-1)+t(13)
      t(16) = t(4)+t(8)
      t(17) = 1.0D0/t(16)
      t(18) = P(-2,-1)*5.0D-1
      t(19) = P(-1,-1)*(-2.0D0)
      t(20) = P(0,-1)*1.5D0
      t(21) = t(18)+t(19)+t(20)
      t(22) = t(9)+t(15)
      t(23) = 1.0D0/t(22)
      t(24) = P(-2,0)*5.0D-1
      t(25) = P(-1,0)*(-2.0D0)
      t(26) = P(0,0)*1.5D0
      t(27) = t(24)+t(25)+t(26)
      t(28) = Q(0,0)*(-1.0D0)
      t(29) = Q(0,-1)+t(28)
      t(30) = 1.0D0/ds2(0,0)
      t(31) = t(23)*t(27)*5.0D-1
      t(32) = ds1(-2,+1)*5.0D-1
      t(38) = ds1(-1,+1)*(-1.5D0)
      t(33) = -t(38)
      t(34) = t(32)-t(33)
      t(35) = 1.0D0/t(34)
      t(36) = P(0,+1)*(-1.0D0)
      t(37) = P(0,0)+t(36)
      t(39) = P(-2,+1)*5.0D-1
      t(40) = P(-1,+1)*(-2.0D0)
      t(41) = P(0,+1)*1.5D0
      t(42) = t(39)+t(40)+t(41)
      t(43) = t(32)+t(38)
      t(44) = 1.0D0/t(43)
      t(45) = Q(0,+1)*(-1.0D0)
      t(46) = Q(0,0)+t(45)
      t(47) = 1.0D0/ds2(-1,-1)
      t(48) = 1.0D0/ds1(-1,0)
      t(49) = 1.0D0/ds1(-1,-1)
      t(50) = P(-1,0)*(-1.0D0)
      t(51) = P(-1,-1)+t(50)
      t(52) = t(47)*t(51)*2.5D-1
      t(53) = t(3)*t(14)*2.5D-1
      t(54) = t(52)+t(53)
      t(55) = t(17)*t(21)*5.0D-1
      t(56) = t(31)+t(55)
      t(57) = 1.0D0/ds1(-1,0)**2
      t(58) = P(-1,0)+t(13)
      t(59) = Q(-1,0)*(-1.0D0)
      t(60) = Q(-1,-1)+t(59)
      t(61) = t(47)*t(60)
      t(62) = t(3)*t(29)
      t(63) = t(61)+t(62)
      t(64) = 1.0D0/ds2(-1,0)
      t(65) = Q(0,-1)*(-1.0D0)
      t(66) = Q(-1,-1)+t(65)
      t(67) = t(49)*t(66)
      t(68) = Q(-1,0)+t(28)
      t(69) = t(48)*t(68)
      t(70) = t(67)+t(69)
      t(71) = P(-1,+1)*(-1.0D0)
      t(72) = P(-1,0)+t(71)
      t(73) = t(64)*t(72)*2.5D-1
      t(74) = t(30)*t(37)*2.5D-1
      t(75) = t(73)+t(74)
      t(76) = t(42)*t(44)*5.0D-1
      t(77) = t(31)+t(76)
      t(78) = 1.0D0/ds1(-1,+1)
      t(79) = Q(-1,+1)+t(45)
      t(80) = t(78)*t(79)
      t(81) = t(69)+t(80)
      t(82) = Q(-1,+1)*(-1.0D0)
      t(83) = Q(-1,0)+t(82)
      t(84) = t(64)*t(83)
      t(85) = t(30)*t(46)
      t(86) = t(84)+t(85)
      t(87) = Q(-2,-1)*5.0D-1
      t(88) = Q(-1,-1)*(-2.0D0)
      t(89) = Q(0,-1)*1.5D0
      t(90) = t(87)+t(88)+t(89)
      t(91) = t(17)*t(90)*5.0D-1
      t(92) = Q(-2,0)*5.0D-1
      t(93) = Q(-1,0)*(-2.0D0)
      t(94) = Q(0,0)*1.5D0
      t(95) = t(92)+t(93)+t(94)
      t(96) = t(23)*t(95)*5.0D-1
      t(97) = t(91)+t(96)
      t(98) = Q(-2,+1)*5.0D-1
      t(99) = Q(-1,+1)*(-2.0D0)
      t(100) = Q(0,+1)*1.5D0
      t(101) = t(98)+t(99)+t(100)
      t(102) = t(44)*t(101)*5.0D-1
      t(103) = t(96)+t(102)

      A0(1) = dBsig22B(-1)*t(2)*t(3)*t(7)*t(14)*(-2.5D-1)+dBsig2
     &1(-1)*t(2)*t(3)*t(17)*t(29)*(t(7)*t(21)*5.0D-1+t(12)*t(27)*5.0D-
     &1)*2.5D-1
      A0(2) = dBsig22B(-1)*t(2)*t(3)*t(12)*t(14)*(-2.5D-1)+dBsig
     &22B(0)*t(2)*t(23)*t(30)*t(37)*2.5D-1+dBsig21(-1)*t(2)*t(3)*t(23)
     &*t(29)*t(56)*2.5D-1-dBsig21(0)*t(2)*t(23)*t(30)*t(46)*(t(31)+t(3
     &5)*t(42)*5.0D-1)*2.5D-1
      A0(3) = dBsig22B(0)*t(2)*t(30)*t(35)*t(37)*2.5D-1-dBsig21(
     &0)*t(2)*t(30)*t(44)*t(46)*t(77)*2.5D-1
      A0(4) = dsig11B*t(47)*t(57)*t(58)*(-5.0D-1)+dBsig22B(-1)*t
     &(2)*t(3)*t(14)*t(17)+dsig12(-1)*t(48)*t(49)*t(54)*t(63)*5.0D-1+d
     &sig12(-1)*t(47)*t(48)*t(54)*t(70)*5.0D-1-dBsig21(-1)*t(2)*t(3)*t
     &(17)*t(29)*t(56)*1.0D0
      A0(5) = t(57)*t(58)*(dsig11A*t(48)*1.0D0-dsig11B*(t(47)*2.
     &5D-1-t(64)*2.5D-1))*(-2.0D0)+dsig12(-1)*t(54)*t(57)*t(63)*5.0D-1
     &+dsig12(0)*t(57)*t(75)*t(86)*5.0D-1+dBsig22B(-1)*t(2)*t(3)*t(14)
     &*t(23)-dBsig22B(0)*t(2)*t(23)*t(30)*t(37)*1.0D0-dsig12(-1)*t(47)
     &*t(48)*t(54)*t(70)*5.0D-1+dsig12(0)*t(48)*t(64)*t(75)*t(81)*5.0D
     &-1-dBsig21(-1)*t(2)*t(3)*t(23)*t(29)*t(56)*1.0D0+dBsig21(0)*t(2)
     &*t(23)*t(30)*t(46)*t(77)
      A0(6) = dsig11B*t(57)*t(58)*t(64)*5.0D-1-dBsig22B(0)*t(2)*
     &t(30)*t(37)*t(44)*1.0D0-dsig12(0)*t(48)*t(64)*t(75)*t(81)*5.0D-1
     &+dsig12(0)*t(48)*t(75)*t(78)*t(86)*5.0D-1+dBsig21(0)*t(2)*t(30)*
     &t(44)*t(46)*t(77)
      A0(7) = dsig11B*t(3)*t(57)*t(58)*(-5.0D-1)-t(2)*t(3)*t(14)
     &*(dBsig22A(-1)*t(3)+dBsig22B(-1)*t(17)*7.5D-1)*1.0D0+dBsig21(-1)
     &*t(2)*t(3)*t(56)*t(97)+dsig12(-1)*t(3)*t(48)*t(54)*t(70)*5.0D-1-
     &dsig12(-1)*t(48)*t(49)*t(54)*t(63)*5.0D-1+dBsig21(-1)*t(2)*t(3)*
     &t(17)*t(29)*t(56)*7.5D-1
      A0(8) = t(57)*t(58)*(dsig11A*t(48)+dsig11B*(t(3)*2.5D-1-t(
     &30)*2.5D-1))*2.0D0-dsig12(-1)*t(54)*t(57)*t(63)*5.0D-1-dsig12(0)
     &*t(57)*t(75)*t(86)*5.0D-1+t(2)*t(3)*t(14)*(dBsig22A(-1)*t(3)*1.0
     &D0-dBsig22B(-1)*t(23)*7.5D-1)*1.0D0+t(2)*t(30)*t(37)*(dBsig22B(0
     &)*t(23)*7.5D-1+dBsig22A(0)*t(30))-dBsig21(-1)*t(2)*t(3)*t(56)*t(
     &97)*1.0D0-dBsig21(0)*t(2)*t(30)*t(77)*t(103)*1.0D0-dsig12(-1)*t(
     &3)*t(48)*t(54)*t(70)*5.0D-1+dsig12(0)*t(30)*t(48)*t(75)*t(81)*5.
     &0D-1+dBsig21(-1)*t(2)*t(3)*t(23)*t(29)*t(56)*7.5D-1-dBsig21(0)*t
     &(2)*t(23)*t(30)*t(46)*t(77)*7.5D-1
      A0(9) = dsig11B*t(30)*t(57)*t(58)*5.0D-1-t(2)*t(30)*t(37)*
     &(dBsig22A(0)*t(30)*1.0D0-dBsig22B(0)*t(44)*7.5D-1)+dBsig21(0)*t(
     &2)*t(30)*t(77)*t(103)-dsig12(0)*t(30)*t(48)*t(75)*t(81)*5.0D-1-d
     &sig12(0)*t(48)*t(75)*t(78)*t(86)*5.0D-1-dBsig21(0)*t(2)*t(30)*t(
     &44)*t(46)*t(77)*7.5D-1



      iA0(1,1)=-2; iA0(1,2)=-1  !3,4
      iA0(2,1)=-2; iA0(2,2)=0   !3,5
      iA0(3,1)=-2; iA0(3,2)=1   !3,6
      iA0(4,1)=-1; iA0(4,2)=-1  !4,4
      iA0(5,1)=-1; iA0(5,2)=0  !4,5
      iA0(6,1)=-1; iA0(6,2)=1   !4,6
      iA0(7,1)=0; iA0(7,2)=-1  !5,4
      iA0(8,1)=0; iA0(8,2)=0   !5,5
      iA0(9,1)=0; iA0(9,2)=1   !5,6
      return
      end    
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine KnonlinearBNDYNi(A0,iA0,P,Q,ds1,ds2,dh1,dh2,
     &   dBsig21,dsig12,dBsig11A,dBsig11B,dsig22A,dsig22B)
      implicit none
      integer iA0(9,2)
      real*8 t(150),A0(9)
      real*8 ds1(-1:1,-1:0),ds2(-1:1,-2:-1)
     &  ,dh2(-1:1,-1:0),dh1(-1:1,-1:0)
     &  ,P(-1:1,-2:0),Q(-1:1,-2:0)
     &  ,dBsig21(-1:0),dsig12(-1:0)
     &  ,dBsig11A(-1:0),dBsig11B(-1:0)
     &  ,dsig22A,dsig22B

      t(2) = 1.0D0/dh1(0,0)
      t(3) = 1.0D0/ds1(-1,0)
      t(4) = ds2(-1,-2)*5.0D-1
      t(8) = ds2(-1,-1)*(-1.5D0)
      t(5) = -t(8)
      t(6) = t(4)-t(5)
      t(7) = 1.0D0/t(6)
      t(9) = P(0,0)*(-1.0D0)
      t(10) = 1.0D0/ds1(-1,-1)
      t(11) = P(-1,0)+t(9)
      t(12) = Q(0,0)*(-1.0D0)
      t(13) = Q(-1,0)+t(12)
      t(14) = 1.0D0/ds2(0,-1)
      t(15) = P(0,-1)*(-1.0D0)
      t(16) = P(-1,-1)+t(15)
      t(17) = t(10)*t(16)*2.5D-1
      t(18) = t(3)*t(11)*2.5D-1
      t(19) = t(17)+t(18)
      t(20) = 1.0D0/ds2(-1,-1)
      t(21) = t(4)+t(8)
      t(22) = 1.0D0/t(21)
      t(23) = P(-1,-2)*5.0D-1
      t(24) = P(-1,-1)*(-2.0D0)
      t(25) = P(-1,0)*1.5D0
      t(26) = t(23)+t(24)+t(25)
      t(27) = ds2(0,-2)*5.0D-1
      t(39) = ds2(0,-1)*(-1.5D0)
      t(28) = -t(39)
      t(29) = t(27)-t(28)
      t(30) = 1.0D0/t(29)
      t(31) = P(0,-2)*5.0D-1
      t(32) = P(0,-1)*(-2.0D0)
      t(33) = P(0,0)*1.5D0
      t(34) = t(31)+t(32)+t(33)
      t(35) = t(30)*t(34)*5.0D-1
      t(36) = 1.0D0/ds2(0,-1)**2
      t(37) = P(0,-1)+t(9)
      t(38) = t(22)*t(26)*5.0D-1
      t(40) = t(27)+t(39)
      t(41) = 1.0D0/t(40)
      t(42) = Q(0,-1)*(-1.0D0)
      t(43) = Q(-1,-1)+t(42)
      t(44) = t(10)*t(43)
      t(45) = t(3)*t(13)
      t(46) = t(44)+t(45)
      t(47) = Q(-1,0)*(-1.0D0)
      t(48) = Q(-1,-1)+t(47)
      t(49) = t(20)*t(48)
      t(50) = Q(0,-1)+t(12)
      t(51) = t(14)*t(50)
      t(52) = t(49)+t(51)
      t(53) = t(34)*t(41)*5.0D-1
      t(54) = t(38)+t(53)
      t(55) = 1.0D0/ds1(0,0)
      t(56) = 1.0D0/ds1(0,-1)
      t(57) = P(+1,0)*(-1.0D0)
      t(58) = P(0,0)+t(57)
      t(59) = Q(+1,0)*(-1.0D0)
      t(60) = Q(0,0)+t(59)
      t(61) = P(+1,-1)*(-1.0D0)
      t(62) = P(0,-1)+t(61)
      t(63) = t(56)*t(62)*2.5D-1
      t(64) = t(55)*t(58)*2.5D-1
      t(65) = t(63)+t(64)
      t(66) = ds2(+1,-2)*5.0D-1
      t(92) = ds2(+1,-1)*(-1.5D0)
      t(67) = -t(92)
      t(68) = t(66)-t(67)
      t(69) = 1.0D0/t(68)
      t(70) = P(+1,-2)*5.0D-1
      t(71) = P(+1,-1)*(-2.0D0)
      t(72) = P(+1,0)*1.5D0
      t(73) = t(70)+t(71)+t(72)
      t(74) = t(69)*t(73)*5.0D-1
      t(75) = t(53)+t(74)
      t(76) = Q(+1,-1)*(-1.0D0)
      t(77) = Q(0,-1)+t(76)
      t(78) = t(56)*t(77)
      t(79) = t(55)*t(60)
      t(80) = t(78)+t(79)
      t(81) = Q(-1,-2)*5.0D-1
      t(82) = Q(-1,-1)*(-2.0D0)
      t(83) = Q(-1,0)*1.5D0
      t(84) = t(81)+t(82)+t(83)
      t(85) = t(22)*t(84)*5.0D-1
      t(86) = Q(0,-2)*5.0D-1
      t(87) = Q(0,-1)*(-2.0D0)
      t(88) = Q(0,0)*1.5D0
      t(89) = t(86)+t(87)+t(88)
      t(90) = t(41)*t(89)*5.0D-1
      t(91) = t(85)+t(90)
      t(93) = t(66)+t(92)
      t(94) = 1.0D0/t(93)
      t(95) = 1.0D0/ds2(+1,-1)
      t(96) = Q(+1,-1)+t(59)
      t(97) = t(95)*t(96)
      t(98) = t(51)+t(97)
      t(99) = t(73)*t(94)*5.0D-1
      t(100) = t(53)+t(99)
      t(101) = Q(+1,-2)*5.0D-1
      t(102) = Q(+1,-1)*(-2.0D0)
      t(103) = Q(+1,0)*1.5D0
      t(104) = t(101)+t(102)+t(103)
      t(105) = t(94)*t(104)*5.0D-1
      t(106) = t(90)+t(105)

      A0(1) = dBsig11B(-1)*t(2)*t(3)*t(7)*t(11)*(-2.5D-1)+dB
     &sig21(-1)*t(2)*t(3)*t(13)*t(22)*(t(35)+t(7)*t(26)*5.0D-1)*2.5D-1
      A0(2) = dsig22B*t(10)*t(36)*t(37)*(-5.0D-1)+dBsig11B(-1
     &)*t(2)*t(3)*t(11)*t(22)+dsig12(-1)*t(10)*t(14)*t(19)*t(52)*5.0D
     &-1+dsig12(-1)*t(14)*t(19)*t(20)*t(46)*5.0D-1-dBsig21(-1)*t(2)*t(
     &3)*t(13)*t(22)*(t(35)+t(38))*1.0D0
      A0(3) = dsig22B*t(3)*t(36)*t(37)*(-5.0D-1)-t(2)*t(3)*t(1
     &1)*(dBsig11A(-1)*t(3)+dBsig11B(-1)*t(22)*7.5D-1)*1.0D0+dBsig21
     &(-1)*t(2)*t(3)*t(54)*t(91)+dsig12(-1)*t(3)*t(14)*t(19)*t(52)*5.0
     &D-1-dsig12(-1)*t(14)*t(19)*t(20)*t(46)*5.0D-1+dBsig21(-1)*t(2)*t
     &(3)*t(13)*t(22)*t(54)*7.5D-1
      A0(4) = dBsig11B(-1)*t(2)*t(3)*t(11)*t(41)*(-2.5D-1)+d
     &Bsig11B(0)*t(2)*t(41)*t(55)*t(58)*2.5D-1+dBsig21(-1)*t(2)*t(3)
     &*t(13)*t(41)*t(54)*2.5D-1-dBsig21(0)*t(2)*t(41)*t(55)*t(60)*t(75
     &)*2.5D-1
      A0(5) = t(36)*t(37)*(dsig22A*t(14)*1.0D0-dsig22B*(t(10)*
     &2.5D-1-t(56)*2.5D-1))*(-2.0D0)+dsig12(-1)*t(19)*t(36)*t(46)*5.0D
     &-1+dsig12(0)*t(36)*t(65)*t(80)*5.0D-1+dBsig11B(-1)*t(2)*t(3)*t
     &(11)*t(41)-dBsig11B(0)*t(2)*t(41)*t(55)*t(58)*1.0D0-dsig12(-1)
     &*t(10)*t(14)*t(19)*t(52)*5.0D-1+dsig12(0)*t(14)*t(56)*t(65)*t(98
     &)*5.0D-1-dBsig21(-1)*t(2)*t(3)*t(13)*t(41)*t(54)*1.0D0+dBsig21(0
     &)*t(2)*t(41)*t(55)*t(60)*t(75)
      A0(6) = t(36)*t(37)*(dsig22A*t(14)+dsig22B*(t(3)*2.5D-1-
     &t(55)*2.5D-1))*2.0D0-dsig12(-1)*t(19)*t(36)*t(46)*5.0D-1-dsig12(
     &0)*t(36)*t(65)*t(80)*5.0D-1+t(2)*t(3)*t(11)*(dBsig11A(-1)*t(3)*1
     &.0D0-dBsig11B(-1)*t(41)*7.5D-1)*1.0D0+t(2)*t(55)*t(58)*(dBsig1
     &1B(0)*t(41)*7.5D-1+dBsig11A(0)*t(55))-dBsig21(-1)*t(2)*t(3)*t(
     &54)*t(91)*1.0D0-dBsig21(0)*t(2)*t(55)*t(100)*t(106)*1.0D0-dsig12
     &(-1)*t(3)*t(14)*t(19)*t(52)*5.0D-1+dsig12(0)*t(14)*t(55)*t(65)*t
     &(98)*5.0D-1+dBsig21(-1)*t(2)*t(3)*t(13)*t(41)*t(54)*7.5D-1-dBsig
     &21(0)*t(2)*t(41)*t(55)*t(60)*t(100)*7.5D-1
      A0(7) = dBsig11B(0)*t(2)*t(55)*t(58)*t(94)*2.5D-1-dBsi
     &g21(0)*t(2)*t(55)*t(60)*t(94)*t(100)*2.5D-1
      A0(8) = dsig22B*t(36)*t(37)*t(56)*5.0D-1-dBsig11B(0)*t
     &(2)*t(55)*t(58)*t(94)*1.0D0-dsig12(0)*t(14)*t(56)*t(65)*t(98)*5.
     &0D-1+dsig12(0)*t(14)*t(65)*t(80)*t(95)*5.0D-1+dBsig21(0)*t(2)*t(
     &55)*t(60)*t(94)*t(100)
      A0(9) = dsig22B*t(36)*t(37)*t(55)*5.0D-1-t(2)*t(55)*t(58
     &)*(dBsig11A(0)*t(55)*1.0D0-dBsig11B(0)*t(94)*7.5D-1)+dBsig21(0
     &)*t(2)*t(55)*t(100)*t(106)-dsig12(0)*t(14)*t(55)*t(65)*t(98)*5.0
     &D-1-dsig12(0)*t(14)*t(65)*t(80)*t(95)*5.0D-1-dBsig21(0)*t(2)*t(5
     &5)*t(60)*t(94)*t(100)*7.5D-1

      iA0(1,1)=-1; iA0(1,2)=-2  !4,3
      iA0(2,1)=-1; iA0(2,2)=-1  !4,4  
      iA0(3,1)=-1; iA0(3,2)=0  !4,5
      iA0(4,1)=0; iA0(4,2)=-2   !5,3
      iA0(5,1)=0; iA0(5,2)=-1   !5,4
      iA0(6,1)=0; iA0(6,2)=0   !5,5
      iA0(7,1)=1; iA0(7,2)=-2   !6,3
      iA0(8,1)=1; iA0(8,2)=-1   !6,4
      iA0(9,1)=1; iA0(9,2)=0   !6,5
      return
      end
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

!-S--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|
!-S--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|
!-S--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|
!-S--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|
!-S--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|
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
      A0=.0
      return
      end
!-E--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|
!-E--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|
!-E--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|
!-E--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|
!-E--*|OBSOLETE--.---------.---------.---------.---------.---------.-|-------|





 

