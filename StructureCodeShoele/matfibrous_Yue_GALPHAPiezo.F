c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine  matfibrous_Yue_GALPHAPiezo(
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
     &        ,Piezo_Formulation_Flag
     &        ,al
     &        ,Vpiezo
     &        ,dXnadt
     &        ,d2Xnadt2
     &        ,aamatPiezo_in,n_aamatPiezo_fsh
     &        ,VecPiezo
     &        ,ImplicitElecMech
     &        ,FlagForceCalculation 
     &        )
c-----------------------------------------------------------------------
      implicit none
      integer nq_ibmin,ns_ibmin,ntermiMAXrow,ntermiMAX
     &       ,nq_ibmBin,ns_ibmBin
      integer iacsrfibershell(ntermiMAXrow+1)
     &       ,jacsrfibershell(ntermiMAX)
      real*8  aacsrfibershell(ntermiMAX)
     &        ,ForceOut(ntermiMAXrow+1)
     &        ,aamat_in(nq_ibmin*ns_ibmin,3,3),aamat_intmp(3,3)
     &        ,alpha_coefK(5),scal_K
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
      real*8 X(-4:4,-4:4),coef,X00,Y(4,4,-4:4,-4:4),temp_ibm3
     &  ,t(220),e5c,e6e,P(-2:2,-2:2),Q(-2:2,-2:2)	
      real*8 YpiezoM(3,-3:3,-3:3)
      integer iAdd(-4:4,-4:4,1:4)
     &  ,iiAdd(-4:4),jiAdd(-4:4),i2,j2,ntemi,noj,iA0(10,2)
     &  ,ist,ien,jst,jen,iflagsimulation,nsnq,i3,j3,nojtmp
      integer FlagForceCalculation(nq_ibmin,ns_ibmin)
      integer i1s(2),i2s(2),iu,iv,iw,ju,jv,jw

!! piezo parts
      integer nYmax,Piezo_Formulation_Flag,TAG0(350,2)
     &  ,nA0max,ixyz,iTAG,MAPTAG0ij(4),iminG,imaxG,jminG,jmaxG

      logical Electromechanical,ImplicitElecMech
      parameter(Electromechanical=.true.)

      integer n_aamatPiezo_fsh
      real*8   al(nq_ibmin,ns_ibmin)
     &        ,Vpiezo(nq_ibmin,ns_ibmin)
     &        ,A0(350)
     &        ,dXnadt(1:3,nq_ibmin,ns_ibmin)
     &        ,d2Xnadt2(1:3,nq_ibmin,ns_ibmin)
     &        ,GdX(1:3,-1:1,-1:1)
     &        ,GddX(1:3,-1:1,-1:1)
     &        ,ddKapa
     &        ,BB(81)
     &        ,aamatPiezo_in(n_aamatPiezo_fsh)
     &        ,VecPiezo(1:3,1:4,nq_ibmin,ns_ibmin)
!!
      real*8 Ypiezo,aacsrfibershellpiezo(nq_ibmin*ns_ibmin*3)
     &      ,aamatPiezo_intmp
      integer nYmax_piezo
!! ================== Inside i,j  2<i<nx-1; 2<j<ny-1
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
      aamatPiezo_intmp=0.0

      ntemi=0;noj=0
      do i=1,nq_ibmin
 
      do j=1,ns_ibmin

        k=(i-1)*ns_ibmin
        A0=0.0;BB=0.0
         TAG0=0;MAPTAG0ij=0

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

         if(n_aamatPiezo_fsh .eq. 1) then
          do i2=-4,4
          do j2=-4,4
             do j3=1,3
                    iAdd(i2,j2,j3)
     &                 =(iiAdd(i2)-1)*ns_ibmin*3+(jiAdd(j2)-1)*3+j3
             enddo
           enddo  
           enddo     
         else   
          do i2=-4,4
          do j2=-4,4
             do j3=1,4
                    iAdd(i2,j2,j3)
     &                 =(iiAdd(i2)-1)*ns_ibmin*4+(jiAdd(j2)-1)*4+j3
             enddo
           enddo  
           enddo
         endif 

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
        do j3=1,3
         do i3=1,3
            Y(i3,j3,0,0)=Y(i3,j3,0,0)+(1.0-coefX00)*
     &           aamat_intmp(i3,j3)/alpha_coefK(1)
            Y(i3,j3,-4:4,-4:4)=alpha_coefK(1)*Y(i3,j3,-4:4,-4:4)
         enddo
        enddo
      nYmax=3


         im1=i-1!+nq_ibmin,nq_ibmin)
        ! if(im1 .eq. 0) im1=nq_ibmin
         im2=i-2!+nq_ibmin,nq_ibmin)
        ! if(im2 .eq. 0) im2=nq_ibmin
         im3=i-3!+nq_ibmin,nq_ibmin)
        ! if(im3 .eq. 0) im3=nq_ibmin
         ip1=i+1!+nq_ibmin,nq_ibmin)
        ! if(ip1 .eq. 0) ip1=nq_ibmin
         ip2=i+2!+nq_ibmin,nq_ibmin)
        ! if(ip2 .eq. 0) ip2=nq_ibmin
         ip3=i+3!+nq_ibmin,nq_ibmin)
        ! if(ip3 .eq. 0) ip3=nq_ibmin

         ihm1=im1
         ihm2=im2
         ihm3=im3

         ihp1=ip1
         ihp2=ip2
         ihp3=ip3


         jm1=j-1!+ns_ibmin,ns_ibmin)
       !  if(jm1 .eq. 0) jm1=ns_ibmin
         jm2=j-2!+ns_ibmin,ns_ibmin)
       !  if(jm2 .eq. 0) jm2=ns_ibmin
         jm3=j-3!+ns_ibmin,ns_ibmin)
       !  if(jm3 .eq. 0) jm3=ns_ibmin
         jp1=j+1!+ns_ibmin,ns_ibmin)
       !  if(jp1 .eq. 0) jp1=ns_ibmin
         jp2=j+2!+ns_ibmin,ns_ibmin)
       !  if(jp2 .eq. 0) jp2=ns_ibmin
         jp3=j+3!+ns_ibmin,ns_ibmin)

         jhm1=jm1
         jhm2=jm2
         jhm3=jm3
      
         jhp1=jp1
         jhp2=jp2
         jhp3=jp3

      A0=0.0
      if(Electromechanical .and. ImplicitElecMech) then
       if (
     &         (ihm2.gt. 0).and.(ihp2.lt.nq_ibmin)
     &         .and.
     &         (jhm2.gt. 0).and.(jhp2.lt.ns_ibmin)
     &        ) then 

        if(Piezo_Formulation_Flag .eq. 0) then
         call PiezoforceFSimple(
     &       ds1(ihm2:ihp2,jhm2:jhp2)
     &      ,ds2(ihm2:ihp2,jhm2:jhp2)
     &      ,dg1(ihm2:ihp2,jhm2:jhp2)
     &      ,dg2(ihm2:ihp2,jhm2:jhp2)
     &      ,dh1(ihm2:ihp2,jhm2:jhp2)
     &      ,dh2(ihm2:ihp2,jhm2:jhp2)
     &      ,GX(1,ihm2:ihp2,jhm2:jhp2)
     &      ,GX(2,ihm2:ihp2,jhm2:jhp2)
     &      ,GX(3,ihm2:ihp2,jhm2:jhp2)
     &      ,al(ihm2:ihp2,jhm2:jhp2)
     &      ,Vpiezo(ihm2:ihp2,jhm2:jhp2)
     &      ,A0(1:200),TAG0(1:200,1:2))   
             nA0max=177  
        elseif(Piezo_Formulation_Flag .eq. 1) then
           call PiezoforceF(
     &       ds1(ihm2:ihp2,jhm2:jhp2)
     &      ,ds2(ihm2:ihp2,jhm2:jhp2)
     &      ,dg1(ihm2:ihp2,jhm2:jhp2)
     &      ,dg2(ihm2:ihp2,jhm2:jhp2)
     &      ,dh1(ihm2:ihp2,jhm2:jhp2)
     &      ,dh2(ihm2:ihp2,jhm2:jhp2)
     &      ,GX(1,ihm2:ihp2,jhm2:jhp2)
     &      ,GX(2,ihm2:ihp2,jhm2:jhp2)
     &      ,GX(3,ihm2:ihp2,jhm2:jhp2)
     &      ,al(ihm2:ihp2,jhm2:jhp2)
     &      ,Vpiezo(ihm2:ihp2,jhm2:jhp2)
     &      ,VecPiezo(1:3,1:4,ihm1:ihp1,jhm1:jhp1)
     &      ,A0(1:350),TAG0(1:350,1:2))      
             nA0max=345
        endif

        if(n_aamatPiezo_fsh .gt. 1) then
           do j3=1,nA0max
               ixyz=TAG0(j3,1)
               iTAG=TAG0(j3,2)
               if(iTAG .le. 100) then
                 call TAGMAPCross(iTAG,MAPTAG0ij) 
                 Y(ixyz,MAPTAG0ij(3),MAPTAG0ij(1),MAPTAG0ij(2))=
     &              Y(ixyz,MAPTAG0ij(3),MAPTAG0ij(1),MAPTAG0ij(2))
     &             +alpha_coefK(2)*A0(j3)
               endif
            enddo
            nYmax=4
            Y(nYmax,nYmax,0,0)=Y(nYmax,nYmax,0,0)+
     &          aamatPiezo_in((i-1)*ns_ibmin+j)  !kmatrix 
            nYmax_piezo=0
        else
           do j3=1,nA0max
               ixyz=TAG0(j3,1)
               iTAG=TAG0(j3,2)
               if(iTAG .le. 100) then
                 call TAGMAPCross(iTAG,MAPTAG0ij) 
                 if(MAPTAG0ij(3) .eq. 4) then
                 Y(ixyz,MAPTAG0ij(3),0,0)=
     &              Y(ixyz,MAPTAG0ij(3),0,0)
     &             +alpha_coefK(2)*A0(j3)
                 else
                 Y(ixyz,MAPTAG0ij(3),MAPTAG0ij(1),MAPTAG0ij(2))=
     &              Y(ixyz,MAPTAG0ij(3),MAPTAG0ij(1),MAPTAG0ij(2))
     &             +alpha_coefK(2)*A0(j3)
                 endif
               endif
            enddo
            nYmax=3
            if((i.eq. 1) .and. (j.eq.1))Ypiezo=aamatPiezo_in(1) !kmatrix 
            nYmax_piezo=1
        endif



              iminG=-1;imaxG=1
              if(i.eq. 1) iminG=0
              if(i.eq. nq_ibmin) imaxG=0

              jminG=-1;jmaxG=1
              if(j.eq. 1) jminG=0
              if(j.eq. ns_ibmin) jmaxG=0
              GdX=0.0;GddX=0.0

              GdX(1:3,iminG:imaxG,jminG:jmaxG)
     &          =dXnadt(1:3
     &                 ,max(1,i-1):min(nq_ibmin,i+1)
     &                 ,max(1,j-1):min(ns_ibmin,j+1))
              GddX(1:3,iminG:imaxG,jminG:jmaxG)
     &          =d2Xnadt2(1:3
     &                 ,max(1,i-1):min(nq_ibmin,i+1)
     &                 ,max(1,j-1):min(ns_ibmin,j+1))

        call ddKapaCalculation(
     &       ds1(ihm1:ihp1,jhm1:jhp1)
     &      ,ds2(ihm1:ihp1,jhm1:jhp1)
     &      ,dg1(ihm1:ihp1,jhm1:jhp1)
     &      ,dg2(ihm1:ihp1,jhm1:jhp1)
     &      ,dh1(ihm1:ihp1,jhm1:jhp1)
     &      ,dh2(ihm1:ihp1,jhm1:jhp1)
     &      ,GX(1,ihm1:ihp1,jhm1:jhp1)
     &      ,GX(2,ihm1:ihp1,jhm1:jhp1)
     &      ,GX(3,ihm1:ihp1,jhm1:jhp1)
     &      ,GdX(1,-1:1,-1:1)
     &      ,GdX(2,-1:1,-1:1)
     &      ,GdX(3,-1:1,-1:1)
     &      ,GddX(1,-1:1,-1:1)
     &      ,GddX(2,-1:1,-1:1)
     &      ,GddX(3,-1:1,-1:1)
     &      ,Piezo_Formulation_Flag,.true.
     &      ,ddKapa 
     &      ,BB(81) 
     &      )
       if(n_aamatPiezo_fsh .gt. 1) then
        aamatPiezo_in((i-1)*ns_ibmin+j)=ddKapa ! output
        YpiezoM=0.0
        do j3=1,81
           call TAGMAPddKapa(j3,MAPTAG0ij)
             Y(nYmax,MAPTAG0ij(2),MAPTAG0ij(3),MAPTAG0ij(4))=
     &              Y(nYmax,MAPTAG0ij(2),MAPTAG0ij(3),MAPTAG0ij(4))
     &             +alpha_coefK(2+MAPTAG0ij(1))*BB(j3)
        enddo  
       else
        aamatPiezo_intmp=aamatPiezo_intmp+ddKapa*ds1(i,j)*ds2(i,j)
        do j3=1,81
           call TAGMAPddKapa(j3,MAPTAG0ij)
             YpiezoM(MAPTAG0ij(2),MAPTAG0ij(3),MAPTAG0ij(4))=
     &              alpha_coefK(2+MAPTAG0ij(1))*BB(j3)*ds1(i,j)*ds2(i,j)
        enddo  
       endif

      else
        YpiezoM=0.0
         aamatPiezo_in((i-1)*ns_ibmin+j)=0.0
      endif      
                
      endif !(Electromechanical .and. ImplicitElecMech) then 



      do i3=1,nYmax
      noj=noj+1
      iacsrfibershell(noj)=ntemi+1
      do j3=1,nYmax
       do i2=-4,4
        do j2=-4,4
            if(abs(Y(i3,j3,i2,j2)).gt. 0.d0) then
              ntemi=ntemi+1
              jacsrfibershell(ntemi)=iAdd(i2,j2,j3)
              aacsrfibershell(ntemi)=Y(i3,j3,i2,j2)
            endif
        enddo
       enddo
       enddo !j3
       do j3=1,nYmax_piezo
       do i2=0,0
        do j2=0,0
            if(abs(Y(i3,3+j3,i2,j2)).gt. 0.d0) then
              ntemi=ntemi+1
              jacsrfibershell(ntemi)=nq_ibmin*ns_ibmin*3+1
              aacsrfibershell(ntemi)=Y(i3,3+j3,i2,j2)
            endif
        enddo
       enddo
       enddo !j3
      enddo !i3

      do i3=1,nYmax_piezo
      do j3=1,nYmax
       do i2=-3,3
        do j2=-3,3
            if(abs(YpiezoM(j3,i2,j2)).gt. 0.d0) then
              aacsrfibershellpiezo(iAdd(i2,j2,j3))=
     &        aacsrfibershellpiezo(iAdd(i2,j2,j3))
     &        +YpiezoM(j3,i2,j2)
            endif
        enddo
       enddo
       enddo !j3
      enddo !i3

      enddo !j
      enddo !i


      do i3=1,nYmax_piezo
      iacsrfibershell(noj)=ntemi+1
      noj=noj+1
      do i2=1,nq_ibmin*ns_ibmin*3
            temp_ibm3=aacsrfibershellpiezo(i2)
            if(temp_ibm3.gt. 0.d0) then  
              ntemi=ntemi+1
              jacsrfibershell(ntemi)=i2
              aacsrfibershell(ntemi)=temp_ibm3
            endif
      enddo
        
       ntemi=ntemi+1
       jacsrfibershell(ntemi)=1+nq_ibmin*ns_ibmin*3
       aacsrfibershell(ntemi)=Ypiezo
       aamatPiezo_in(1)=aamatPiezo_intmp
       enddo !i3

      noj=noj+1
      iacsrfibershell(noj)=ntemi+1 
      if(ntemi .ge. ntermiMAX) 
     &      stop "Larger matrix space for ja aa in piezo"
      return 
      end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
