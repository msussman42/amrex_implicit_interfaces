c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine matfibrous_Huang_GALPHAPiezo(
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


      real*8 DX_i_j(3),DX_ip_j(3),DX_im_j(3)
     &      ,DY_i_j(3),DY_ip_j(3),DY_im_j(3)
     &      ,DX_i_jp(3),DX_ip_jp(3),DX_im_jp(3)
     &      ,DY_i_jp(3),DY_ip_jp(3),DY_im_jp(3)
     &      ,DX_i_jm(3),DX_ip_jm(3),DX_im_jm(3)
     &      ,DY_i_jm(3),DY_ip_jm(3),DY_im_jm(3)
     &      ,DXX_i_j(3),DYY_i_j(3),DXY_i_j(3),DYX_i_j(3)
     
      integer i,j,im1,im2,im3,jm1,jm2,jm3,ihm1,ihm2,ihm3,jhm1,jhm2,jhm3
      integer ip1,ip2,ip3,jp1,jp2,jp3,ihp1,ihp2,ihp3,jhp1,jhp2,jhp3,k
      real*8 X(-3:3,-3:3),coef,X00,Y(4,4,-3:3,-3:3),temp_ibm3
      real*8 YpiezoM(3,-3:3,-3:3)
      integer iAdd(-3:3,-3:3,1:4)
     &  ,iiAdd(-3:3),jiAdd(-3:3),i2,j2,ntemi,noj
     &  ,ist,ien,jst,jen,iflagsimulation,nsnq,i3,j3,nojtmp

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

         im1=mod(i-1+nq_ibmin,nq_ibmin)
         if(im1 .eq. 0) im1=nq_ibmin
         im2=mod(i-2+nq_ibmin,nq_ibmin)
         if(im2 .eq. 0) im2=nq_ibmin
         im3=mod(i-3+nq_ibmin,nq_ibmin)
         if(im3 .eq. 0) im3=nq_ibmin
         ip1=mod(i+1+nq_ibmin,nq_ibmin)
         if(ip1 .eq. 0) ip1=nq_ibmin
         ip2=mod(i+2+nq_ibmin,nq_ibmin)
         if(ip2 .eq. 0) ip2=nq_ibmin
         ip3=mod(i+3+nq_ibmin,nq_ibmin)
         if(ip3 .eq. 0) ip3=nq_ibmin

         ihm1=im1
         ihm2=im2
         ihm3=im3

         ihp1=ip1
         ihp2=ip2
         ihp3=ip3


         jm1=mod(j-1+ns_ibmin,ns_ibmin)
         if(jm1 .eq. 0) jm1=ns_ibmin
         jm2=mod(j-2+ns_ibmin,ns_ibmin)
         if(jm2 .eq. 0) jm2=ns_ibmin
         jm3=mod(j-3+ns_ibmin,ns_ibmin)
         if(jm3 .eq. 0) jm3=ns_ibmin
         jp1=mod(j+1+ns_ibmin,ns_ibmin)
         if(jp1 .eq. 0) jp1=ns_ibmin
         jp2=mod(j+2+ns_ibmin,ns_ibmin)
         if(jp2 .eq. 0) jp2=ns_ibmin
         jp3=mod(j+3+ns_ibmin,ns_ibmin)
         if(jp3 .eq. 0) jp3=ns_ibmin

         jhm1=jm1
         jhm2=jm2
         jhm3=jm3
      
         jhp1=jp1
         jhp2=jp2
         jhp3=jp3

         
         iiAdd(-3:3)=(/im3,im2,im1,i,ip1,ip2,ip3/)
         jiAdd(-3:3)=(/jm3,jm2,jm1,j,jp1,jp2,jp3/)

         if(n_aamatPiezo_fsh .eq. 1) then
            do i2=-3,3
              do j2=-3,3
                do j3=1,3
                    iAdd(i2,j2,j3)
     &                 =(iiAdd(i2)-1)*ns_ibmin*3+(jiAdd(j2)-1)*3+j3
                enddo
              enddo  
            enddo      
         else
            do i2=-3,3
              do j2=-3,3
                do j3=1,4
                    iAdd(i2,j2,j3)
     &                 =(iiAdd(i2)-1)*ns_ibmin*4+(jiAdd(j2)-1)*4+j3
                enddo
              enddo  
            enddo   
         endif

       X=0.d0
       X00=0.0
       if((j .le. 2).and. (BoundaryFibShellin(1,1) .eq. 0))then
          X(-3:3,-3:3)=0.0;X00=1.0d0
       elseif((j .ge. ns_ibmin-1)
     &       .and. (BoundaryFibShellin(1,2) .eq. 0))then
          X(-3:3,-3:3)=0.0;X00=1.0d0
       elseif((i .le. 2).and. (BoundaryFibShellin(2,1) .eq. 0))then
          X(-3:3,-3:3)=0.0;X00=1.0d0
       elseif((i .ge. nq_ibmin-1)
     &       .and. (BoundaryFibShellin(2,2) .eq. 0))then
          X(-3:3,-3:3)=0.0;X00=1.0d0

       elseif((j .le. 1).and. (BoundaryFibShellin(1,1) .eq. 2))then
          X(-3:3,-3:3)=0.0;X00=1.0d0
       elseif((j .ge. ns_ibmin)
     &       .and. (BoundaryFibShellin(1,2) .eq. 2))then
          X(-3:3,-3:3)=0.0;X00=1.0d0
       elseif((i .le. 1).and. (BoundaryFibShellin(2,1) .eq. 2))then
          X(-3:3,-3:3)=0.0;X00=1.0d0
       elseif((i .ge. nq_ibmin)
     &       .and. (BoundaryFibShellin(2,2) .eq. 2))then
          X(-3:3,-3:3)=0.0;X00=1.0d0
!! ================== BNDY 1,1
       elseif(((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   .and.
     &    ((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   )then
 
       X(0,0)=
     & gam11(ip1,j)/ (dh1(ip1,j)* ds1(i,j)**2.0d0* ds1(ihp1,j)) - 
     & (2.0d0* sig22(i,j))/ ds2(i,j)**2 - (2.0d0* sig11(i,j))/ 
     & ds1(i,j)**2 + gam22(i,jp1)/ (dh2(i,jp1)* ds2(i,j)**2.0d0* 
     & ds2(i,jhp1)) - (4.0d0* gam12(i,j))/ (dg1(i,j)* dg2(i,j)* 
     & ds1(i,j)* ds2(i,j)) - (4.0d0* gam21(i,j))/ (dg1(i,j)* 
     & dg2(i,j)* ds1(i,j)* ds2(i,j))

       X(0,+1)=
     & (2.0d0* sig22(i,j))/ ds2(i,j)**2 - (gam22(i,jp1)* ((1.0d0/ 
     & ds2(i,j) + 1.0d0/ ds2(i,jhp1))/ dh2(i,jp1) + 1.0d0/ 
     & (dh2(i,jp2)* ds2(i,jhp1))))/ (ds2(i,j)* ds2(i,jhp1)) + 
     & (4.0d0* gam12(i,j))/ (dg1(i,j)* dg2(i,j)* ds1(i,j)* 
     & ds2(i,j)) + (4.0d0* gam21(i,j))/ (dg1(i,j)* dg2(i,j)* 
     & ds1(i,jp1)* ds2(i,j))

       X(0,+2)=
     & (gam22(i,jp1)* ((1.0d0/ ds2(i,jhp1) + 1.0d0/ ds2(i,jhp2))/ 
     & dh2(i,jp2) + 1.0d0/ (dh2(i,jp1)* ds2(i,jhp1))))/ (ds2(i,j)* 
     & ds2(i,jhp1))

       X(0,+3)=
     & -gam22(i,jp1)/ (dh2(i,jp2)* ds2(i,j)* ds2(i,jhp1)* 
     & ds2(i,jhp2))

       X(+1,0)=
     & (2.0d0* sig11(i,j))/ ds1(i,j)**2 - (gam11(ip1,j)* ((1.0d0/ 
     & ds1(i,j) + 1.0d0/ ds1(ihp1,j))/ dh1(ip1,j) + 1.0d0/ 
     & (dh1(ip2,j)* ds1(ihp1,j))))/ (ds1(i,j)* ds1(ihp1,j)) + 
     & (4.0d0* gam12(i,j))/ (dg1(i,j)* dg2(i,j)* ds1(i,j)* 
     & ds2(ip1,j)) + (4.0d0* gam21(i,j))/ (dg1(i,j)* dg2(i,j)* 
     & ds1(i,j)* ds2(i,j))

       X(+1,+1)=
     & - (4.0d0* gam12(i,j))/ (dg1(i,j)* dg2(i,j)* ds1(i,j)* 
     &   ds2(ip1,j)) - (4.0d0* gam21(i,j))/ (dg1(i,j)* dg2(i,j)* 
     & ds1(i,jp1)* ds2(i,j))

       X(+2,0)=
     & (gam11(ip1,j)* ((1.0d0/ ds1(ihp1,j) + 1.0d0/ ds1(ihp2,j))/ 
     & dh1(ip2,j) + 1.0d0/ (dh1(ip1,j)* ds1(ihp1,j))))/ (ds1(i,j)* 
     & ds1(ihp1,j))

       X(+3,0)=
     & -gam11(ip1,j)/ (dh1(ip2,j)* ds1(i,j)* ds1(ihp1,j)* 
     & ds1(ihp2,j))


!! ================== BNDY Nx,1
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   )then
 
       X(-3,0)=
     & -gam11(im1,j)/ (dh1(im2,j)* ds1(ihm3,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))

       X(-2,0)=
     & (gam11(im1,j)* ((1.0d0/ ds1(ihm3,j) + 1.0d0/ ds1(ihm2,j))/ 
     & dh1(im2,j) + 1.0d0/ (dh1(im1,j)* ds1(ihm2,j))))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j))

       X(-1,0)=
     & (2.0d0* sig11(ihm1,j))/ ds1(ihm1,j)**2 - (gam11(im1,j)* 
     & ((1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j))/ dh1(im1,j) + 
     & 1.0d0/ (dh1(im2,j)* ds1(ihm2,j))))/ (ds1(ihm2,j)* 
     & ds1(ihm1,j)) + (4.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* ds1(ihm1,j)* ds2(im1,j)) + (4.0d0* 
     & gam21(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* ds1(ihm1,j)* 
     & ds2(i,j))

       X(-1,+1)=
     & - (4.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* 
     &   ds1(ihm1,j)* ds2(im1,j)) - (4.0d0* gam21(ihm1,j))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* ds1(ihm1,jp1)* ds2(i,j))

       X(0,0)=
     & gam11(im1,j)/ (dh1(im1,j)* ds1(ihm2,j)* ds1(ihm1,j)**2) - 
     & (2.0d0* sig22(i,j))/ ds2(i,j)**2 - (2.0d0* sig11(ihm1,j))/ 
     & ds1(ihm1,j)**2 + gam22(i,jp1)/ (dh2(i,jp1)* ds2(i,j)**2.0d0* 
     & ds2(i,jhp1)) - (4.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* ds1(ihm1,j)* ds2(i,j)) - (4.0d0* 
     & gam21(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* ds1(ihm1,j)* 
     & ds2(i,j))

       X(0,+1)=
     & (2.0d0* sig22(i,j))/ ds2(i,j)**2 - (gam22(i,jp1)* ((1.0d0/ 
     & ds2(i,j) + 1.0d0/ ds2(i,jhp1))/ dh2(i,jp1) + 1.0d0/ 
     & (dh2(i,jp2)* ds2(i,jhp1))))/ (ds2(i,j)* ds2(i,jhp1)) + 
     & (4.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* 
     & ds1(ihm1,j)* ds2(i,j)) + (4.0d0* gam21(ihm1,j))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* ds1(ihm1,jp1)* ds2(i,j))

       X(0,+2)=
     & (gam22(i,jp1)* ((1.0d0/ ds2(i,jhp1) + 1.0d0/ ds2(i,jhp2))/ 
     & dh2(i,jp2) + 1.0d0/ (dh2(i,jp1)* ds2(i,jhp1))))/ (ds2(i,j)* 
     & ds2(i,jhp1))

       X(0,+3)=
     & -gam22(i,jp1)/ (dh2(i,jp2)* ds2(i,j)* ds2(i,jhp1)* 
     & ds2(i,jhp2))

!! ================== BNDY 1,Ny
       elseif(((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
 
       X(0,-3)=
     & -gam22(i,jm1)/ (dh2(i,jm2)* ds2(i,jhm3)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-2)=
     & (gam22(i,jm1)* ((1.0d0/ ds2(i,jhm3) + 1.0d0/ ds2(i,jhm2))/ 
     & dh2(i,jm2) + 1.0d0/ (dh2(i,jm1)* ds2(i,jhm2))))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1))

       X(0,-1)=
     & (2.0d0* sig22(i,j-1))/ ds2(i,jhm1)**2 - (gam22(i,jm1)* 
     & ((1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1))/ dh2(i,jm1) + 
     & 1.0d0/ (dh2(i,jm2)* ds2(i,jhm2))))/ (ds2(i,jhm2)* 
     & ds2(i,jhm1)) + (4.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds1(i,j)* ds2(i,jhm1)) + (4.0d0* 
     & gam21(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,jm1)* 
     & ds2(i,jhm1))

       X(0,0)=
     & gam11(ip1,j)/ (dh1(ip1,j)* ds1(i,j)**2.0d0* ds1(ihp1,j)) - 
     & (2.0d0* sig22(i,j-1))/ ds2(i,jhm1)**2 - (2.0d0* sig11(i,j))/ 
     & ds1(i,j)**2 + gam22(i,jm1)/ (dh2(i,jm1)* ds2(i,jhm2)* 
     & ds2(i,jhm1)**2) - (4.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds1(i,j)* ds2(i,jhm1)) - (4.0d0* 
     & gam21(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,j)* 
     & ds2(i,jhm1))

       X(+1,-1)=
     & - (4.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* 
     &   ds1(i,j)* ds2(ip1,jhm1)) - (4.0d0* gam21(i,jhm1))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,jm1)* ds2(i,jhm1))

       X(+1,0)=
     & (2.0d0* sig11(i,j))/ ds1(i,j)**2 - (gam11(ip1,j)* ((1.0d0/ 
     & ds1(i,j) + 1.0d0/ ds1(ihp1,j))/ dh1(ip1,j) + 1.0d0/ 
     & (dh1(ip2,j)* ds1(ihp1,j))))/ (ds1(i,j)* ds1(ihp1,j)) + 
     & (4.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,j)* 
     & ds2(ip1,jhm1)) + (4.0d0* gam21(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds1(i,j)* ds2(i,jhm1))

       X(+2,0)=
     & (gam11(ip1,j)* ((1.0d0/ ds1(ihp1,j) + 1.0d0/ ds1(ihp2,j))/ 
     & dh1(ip2,j) + 1.0d0/ (dh1(ip1,j)* ds1(ihp1,j))))/ (ds1(i,j)* 
     & ds1(ihp1,j))

       X(+3,0)=
     & -gam11(ip1,j)/ (dh1(ip2,j)* ds1(i,j)* ds1(ihp1,j)* 
     & ds1(ihp2,j))


!! ================== BNDY Nx,Ny
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
 
       X(-3,0)=
     & -gam11(im1,j)/ (dh1(im2,j)* ds1(ihm3,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))

       X(-2,0)=
     & (gam11(im1,j)* ((1.0d0/ ds1(ihm3,j) + 1.0d0/ ds1(ihm2,j))/ 
     & dh1(im2,j) + 1.0d0/ (dh1(im1,j)* ds1(ihm2,j))))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j))

       X(-1,-1)=
     & - (4.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* 
     &   dg2(ihm1,jhm1)* ds1(ihm1,j)* ds2(im1,jhm1)) - (4.0d0* 
     & gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,jm1)* ds2(i,jhm1))

       X(-1,0)=
     & (2.0d0* sig11(ihm1,j))/ ds1(ihm1,j)**2 - (gam11(im1,j)* 
     & ((1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j))/ dh1(im1,j) + 
     & 1.0d0/ (dh1(im2,j)* ds1(ihm2,j))))/ (ds1(ihm2,j)* 
     & ds1(ihm1,j)) + (4.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* 
     & dg2(ihm1,jhm1)* ds1(ihm1,j)* ds2(im1,jhm1)) + (4.0d0* 
     & gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,j)* ds2(i,jhm1))

       X(0,-3)=
     & -gam22(i,jm1)/ (dh2(i,jm2)* ds2(i,jhm3)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-2)=
     & (gam22(i,jm1)* ((1.0d0/ ds2(i,jhm3) + 1.0d0/ ds2(i,jhm2))/ 
     & dh2(i,jm2) + 1.0d0/ (dh2(i,jm1)* ds2(i,jhm2))))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1))

       X(0,-1)=
     & (2.0d0* sig22(i,j-1))/ ds2(i,jhm1)**2 - (gam22(i,jm1)* 
     & ((1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1))/ dh2(i,jm1) + 
     & 1.0d0/ (dh2(i,jm2)* ds2(i,jhm2))))/ (ds2(i,jhm2)* 
     & ds2(i,jhm1)) + (4.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* 
     & dg2(ihm1,jhm1)* ds1(ihm1,j)* ds2(i,jhm1)) + (4.0d0* 
     & gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,jm1)* ds2(i,jhm1))

       X(0,0)=
     & gam11(im1,j)/ (dh1(im1,j)* ds1(ihm2,j)* ds1(ihm1,j)**2) - 
     & (2.0d0* sig22(i,j-1))/ ds2(i,jhm1)**2 - (2.0d0* 
     & sig11(ihm1,j))/ ds1(ihm1,j)**2 + gam22(i,jm1)/ (dh2(i,jm1)* 
     & ds2(i,jhm2)* ds2(i,jhm1)**2) - (4.0d0* gam12(ihm1,jhm1))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* ds1(ihm1,j)* ds2(i,jhm1)) - 
     & (4.0d0* gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,j)* ds2(i,jhm1))

!! ================== BNDY 1,2
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .eq. 2)) then
       X(0,-1)=
     & (gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)**2) + gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ 
     & (dh2(i,j)* ds2(i,jhm1)) - sig12(i,jhm1)/ (2.0d0* ds1(i,j)* 
     & ds2(i,jhm1)) + (2.0d0* gam21(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* dk2(i,j)* ds1(i,jm1)) + (2.0d0* gam12(i,jhm1))/ 
     & (dg1(i,jhm1)* dk2(i,j)* ds1(i,j)* ds2(i,jhm1))

       X(0,0)=
     & (sig12(i,jhm1)/ (2.0d0* ds2(i,jhm1)) - sig12(i,j)/ (2.0d0* 
     & ds2(i,j)))/ ds1(i,j) - ((2.0d0* gam21(i,jhm1))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,j)) + (2.0d0* gam21(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds1(i,j)))/ dk2(i,j) - (sig22(i,j-1)/ 
     & ds2(i,jhm1) + sig22(i,j)/ ds2(i,j))/ dh2(i,j) - (2.0d0* 
     & sig11(i,j))/ ds1(i,j)**2 - ((gam22(i,jp1)/ (dh2(i,jp1)* 
     & ds2(i,j)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 1.0d0/ 
     & ds2(i,j)))/ dh2(i,j))/ ds2(i,j) + (gam22(i,j)* (1.0d0/ 
     & ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* ds2(i,jhm1)))/ 
     & dh2(i,j) - (2.0d0* (gam12(i,jhm1)/ (dg1(i,jhm1)* 
     & ds2(i,jhm1)) + gam12(i,j)/ (dg1(i,j)* ds2(i,j))))/ 
     & (dk2(i,j)* ds1(i,j)) + gam11(ip1,j)/ (dh1(ip1,j)* 
     & ds1(i,j)**2.0d0* ds1(ihp1,j))

       X(0,+1)=
     & ((gam22(i,j)/ (dh2(i,j)* ds2(i,j)) + (gam22(i,jp1)* (1.0d0/ 
     & ds2(i,j) + 1.0d0/ ds2(i,jhp1)))/ dh2(i,jp1))/ ds2(i,j) + 
     & gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + 
     & sig22(i,j)/ (dh2(i,j)* ds2(i,j)) + sig12(i,j)/ (2.0d0* 
     & ds1(i,j)* ds2(i,j)) + (2.0d0* gam21(i,j))/ (dg1(i,j)* 
     & dg2(i,j)* dk2(i,j)* ds1(i,jp1)) + (2.0d0* gam12(i,j))/ 
     & (dg1(i,j)* dk2(i,j)* ds1(i,j)* ds2(i,j))

       X(0,+2)=
     & -gam22(i,jp1)/ (dh2(i,j)* dh2(i,jp1)* ds2(i,j)* ds2(i,jhp1))

       X(+1,-1)=
     & - sig12(i,jhm1)/ (2.0d0* ds1(i,j)* ds2(ip1,jhm1)) - (2.0d0* 
     &   gam21(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* dk2(i,j)* 
     & ds1(i,jm1)) - (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dk2(i,j)* ds1(i,j)* ds2(ip1,jhm1))

       X(+1,0)=
     & (2.0d0* sig11(i,j))/ ds1(i,j)**2 + ((2.0d0* gam21(i,jhm1))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,j)) + (2.0d0* gam21(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds1(i,j)))/ dk2(i,j) + (sig12(i,jhm1)/ 
     & (2.0d0* ds2(ip1,jhm1)) - sig12(i,j)/ (2.0d0* ds2(ip1,j)))/ 
     & ds1(i,j) + (2.0d0* (gam12(i,jhm1)/ (dg1(i,jhm1)* 
     & ds2(ip1,jhm1)) + gam12(i,j)/ (dg1(i,j)* ds2(ip1,j))))/ 
     & (dk2(i,j)* ds1(i,j)) - (gam11(ip1,j)* ((1.0d0/ ds1(i,j) + 
     & 1.0d0/ ds1(ihp1,j))/ dh1(ip1,j) + 1.0d0/ (dh1(ip2,j)* 
     & ds1(ihp1,j))))/ (ds1(i,j)* ds1(ihp1,j))

       X(+1,+1)=
     & sig12(i,j)/ (2.0d0* ds1(i,j)* ds2(ip1,j)) - (2.0d0* 
     & gam21(i,j))/ (dg1(i,j)* dg2(i,j)* dk2(i,j)* ds1(i,jp1)) - 
     & (2.0d0* gam12(i,j))/ (dg1(i,j)* dk2(i,j)* ds1(i,j)* 
     & ds2(ip1,j))

       X(+2,0)=
     & (gam11(ip1,j)* ((1.0d0/ ds1(ihp1,j) + 1.0d0/ ds1(ihp2,j))/ 
     & dh1(ip2,j) + 1.0d0/ (dh1(ip1,j)* ds1(ihp1,j))))/ (ds1(i,j)* 
     & ds1(ihp1,j))

       X(+3,0)=
     & -gam11(ip1,j)/ (dh1(ip2,j)* ds1(i,j)* ds1(ihp1,j)* 
     & ds1(ihp2,j))
!! ================== BNDY 1,Ny-1
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .eq. ns_ibmin-1)) then
       X(0,-2)=
     & -gam22(i,jm1)/ (dh2(i,jm1)* dh2(i,j)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-1)=
     & ((gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)) + (gam22(i,jm1)* 
     & (1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1)))/ dh2(i,jm1))/ 
     & ds2(i,jhm1) + gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* 
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & - sig12(i,jhm1)/ (2.0d0* ds1(i,j)* ds2(i,jhm1)) + (2.0d0* 
     & gam21(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* dk2(i,j)* 
     & ds1(i,jm1)) + (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dk2(i,j)* ds1(i,j)* ds2(i,jhm1))

       X(0,0)=
     & (sig12(i,jhm1)/ (2.0d0* ds2(i,jhm1)) - sig12(i,j)/ (2.0d0* 
     & ds2(i,j)))/ ds1(i,j) - ((2.0d0* gam21(i,jhm1))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,j)) + (2.0d0* gam21(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds1(i,j)))/ dk2(i,j) - (sig22(i,j-1)/ 
     & ds2(i,jhm1) + sig22(i,j)/ ds2(i,j))/ dh2(i,j) - (2.0d0* 
     & sig11(i,j))/ ds1(i,j)**2 - ((gam22(i,jm1)/ (dh2(i,jm1)* 
     & ds2(i,jhm1)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 1.0d0/ 
     & ds2(i,j)))/ dh2(i,j))/ ds2(i,jhm1) + (gam22(i,j)* (1.0d0/ 
     & ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* ds2(i,j)))/ 
     & dh2(i,j) - (2.0d0* (gam12(i,jhm1)/ (dg1(i,jhm1)* 
     & ds2(i,jhm1)) + gam12(i,j)/ (dg1(i,j)* ds2(i,j))))/ 
     & (dk2(i,j)* ds1(i,j)) + gam11(ip1,j)/ (dh1(ip1,j)* 
     & ds1(i,j)**2.0d0* ds1(ihp1,j))

       X(0,+1)=
     & (gam22(i,j)/ (dh2(i,j)* ds2(i,j)**2) + gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j)/ 
     & (dh2(i,j)* ds2(i,j)) + sig12(i,j)/ (2.0d0* ds1(i,j)* 
     & ds2(i,j)) + (2.0d0* gam21(i,j))/ (dg1(i,j)* dg2(i,j)* 
     & dk2(i,j)* ds1(i,jp1)) + (2.0d0* gam12(i,j))/ (dg1(i,j)* 
     & dk2(i,j)* ds1(i,j)* ds2(i,j))

       X(+1,-1)=
     & - sig12(i,jhm1)/ (2.0d0* ds1(i,j)* ds2(ip1,jhm1)) - (2.0d0* 
     &   gam21(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* dk2(i,j)* 
     & ds1(i,jm1)) - (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dk2(i,j)* ds1(i,j)* ds2(ip1,jhm1))

       X(+1,0)=
     & (2.0d0* sig11(i,j))/ ds1(i,j)**2 + ((2.0d0* gam21(i,jhm1))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,j)) + (2.0d0* gam21(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds1(i,j)))/ dk2(i,j) + (sig12(i,jhm1)/ 
     & (2.0d0* ds2(ip1,jhm1)) - sig12(i,j)/ (2.0d0* ds2(ip1,j)))/ 
     & ds1(i,j) + (2.0d0* (gam12(i,jhm1)/ (dg1(i,jhm1)* 
     & ds2(ip1,jhm1)) + gam12(i,j)/ (dg1(i,j)* ds2(ip1,j))))/ 
     & (dk2(i,j)* ds1(i,j)) - (gam11(ip1,j)* ((1.0d0/ ds1(i,j) + 
     & 1.0d0/ ds1(ihp1,j))/ dh1(ip1,j) + 1.0d0/ (dh1(ip2,j)* 
     & ds1(ihp1,j))))/ (ds1(i,j)* ds1(ihp1,j))

       X(+1,+1)=
     & sig12(i,j)/ (2.0d0* ds1(i,j)* ds2(ip1,j)) - (2.0d0* 
     & gam21(i,j))/ (dg1(i,j)* dg2(i,j)* dk2(i,j)* ds1(i,jp1)) - 
     & (2.0d0* gam12(i,j))/ (dg1(i,j)* dk2(i,j)* ds1(i,j)* 
     & ds2(ip1,j))

       X(+2,0)=
     & (gam11(ip1,j)* ((1.0d0/ ds1(ihp1,j) + 1.0d0/ ds1(ihp2,j))/ 
     & dh1(ip2,j) + 1.0d0/ (dh1(ip1,j)* ds1(ihp1,j))))/ (ds1(i,j)* 
     & ds1(ihp1,j))

       X(+3,0)=
     & -gam11(ip1,j)/ (dh1(ip2,j)* ds1(i,j)* ds1(ihp1,j)* 
     & ds1(ihp2,j))

!! ================== BNDY Nx,2
       elseif((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .eq. 2)) then

       X(-3,0)=
     & -gam11(im1,j)/ (dh1(im2,j)* ds1(ihm3,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))

       X(-2,0)=
     & (gam11(im1,j)* ((1.0d0/ ds1(ihm3,j) + 1.0d0/ ds1(ihm2,j))/ 
     & dh1(im2,j) + 1.0d0/ (dh1(im1,j)* ds1(ihm2,j))))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j))

       X(-1,-1)=
     & sig12(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)* ds2(im1,jhm1)) - 
     & (2.0d0* gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & dk2(ihm1,j)* ds1(ihm1,jm1)) - (2.0d0* gam12(ihm1,jhm1))/ 
     & (dg1(ihm1,jhm1)* dk2(ihm1,j)* ds1(ihm1,j)* ds2(im1,jhm1))

       X(-1,0)=
     & (2.0d0* sig11(ihm1,j))/ ds1(ihm1,j)**2 + ((2.0d0* 
     & gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,j)) + (2.0d0* gam21(ihm1,j))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* ds1(ihm1,j)))/ dk2(ihm1,j) - (sig12(ihm1,jhm1)/ 
     & (2.0d0* ds2(im1,jhm1)) - sig12(ihm1,j)/ (2.0d0* 
     & ds2(im1,j)))/ ds1(ihm1,j) + (2.0d0* (gam12(ihm1,jhm1)/ 
     & (dg1(ihm1,jhm1)* ds2(im1,jhm1)) + gam12(ihm1,j)/ 
     & (dg1(ihm1,j)* ds2(im1,j))))/ (dk2(ihm1,j)* ds1(ihm1,j)) - 
     & (gam11(im1,j)* ((1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j))/ 
     & dh1(im1,j) + 1.0d0/ (dh1(im2,j)* ds1(ihm2,j))))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j))

       X(-1,+1)=
     & - sig12(ihm1,j)/ (2.0d0* ds1(ihm1,j)* ds2(im1,j)) - (2.0d0* 
     &   gam21(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* dk2(ihm1,j)* 
     & ds1(ihm1,jp1)) - (2.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* 
     & dk2(ihm1,j)* ds1(ihm1,j)* ds2(im1,j))

       X(0,-1)=
     & (gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)**2) + gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ 
     & (dh2(i,j)* ds2(i,jhm1)) + sig12(ihm1,jhm1)/ (2.0d0* 
     & ds1(ihm1,j)* ds2(i,jhm1)) + (2.0d0* gam21(ihm1,jhm1))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* dk2(ihm1,j)* ds1(ihm1,jm1)) 
     & + (2.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dk2(ihm1,j)* 
     & ds1(ihm1,j)* ds2(i,jhm1))

       X(0,0)=
     & gam11(im1,j)/ (dh1(im1,j)* ds1(ihm2,j)* ds1(ihm1,j)**2) - 
     & ((2.0d0* gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,j)) + (2.0d0* gam21(ihm1,j))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* ds1(ihm1,j)))/ dk2(ihm1,j) - (sig22(i,j-1)/ 
     & ds2(i,jhm1) + sig22(i,j)/ ds2(i,j))/ dh2(i,j) - 
     & (sig12(ihm1,jhm1)/ (2.0d0* ds2(i,jhm1)) - sig12(ihm1,j)/ 
     & (2.0d0* ds2(i,j)))/ ds1(ihm1,j) - ((gam22(i,jp1)/ 
     & (dh2(i,jp1)* ds2(i,j)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 
     & 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,j) + (gam22(i,j)* 
     & (1.0d0/ ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* 
     & ds2(i,jhm1)))/ dh2(i,j) - (2.0d0* (gam12(ihm1,jhm1)/ 
     & (dg1(ihm1,jhm1)* ds2(i,jhm1)) + gam12(ihm1,j)/ (dg1(ihm1,j)* 
     & ds2(i,j))))/ (dk2(ihm1,j)* ds1(ihm1,j)) - (2.0d0* 
     & sig11(ihm1,j))/ ds1(ihm1,j)**2

       X(0,+1)=
     & ((gam22(i,j)/ (dh2(i,j)* ds2(i,j)) + (gam22(i,jp1)* (1.0d0/ 
     & ds2(i,j) + 1.0d0/ ds2(i,jhp1)))/ dh2(i,jp1))/ ds2(i,j) + 
     & gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + 
     & sig22(i,j)/ (dh2(i,j)* ds2(i,j)) - sig12(ihm1,j)/ (2.0d0* 
     & ds1(ihm1,j)* ds2(i,j)) + (2.0d0* gam21(ihm1,j))/ 

     & (dg1(ihm1,j)* dg2(ihm1,j)* dk2(ihm1,j)* ds1(ihm1,jp1)) + 
     & (2.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* dk2(ihm1,j)* 
     & ds1(ihm1,j)* ds2(i,j))

       X(0,+2)=
     & -gam22(i,jp1)/ (dh2(i,j)* dh2(i,jp1)* ds2(i,j)* ds2(i,jhp1))

!! ================== BNDY Nx,Ny-1
       elseif((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .eq. ns_ibmin-1)) then
       X(-3,0)=
     & -gam11(im1,j)/ (dh1(im2,j)* ds1(ihm3,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))

       X(-2,0)=
     & (gam11(im1,j)* ((1.0d0/ ds1(ihm3,j) + 1.0d0/ ds1(ihm2,j))/ 
     & dh1(im2,j) + 1.0d0/ (dh1(im1,j)* ds1(ihm2,j))))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j))

       X(-1,-1)=
     & sig12(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)* ds2(im1,jhm1)) - 
     & (2.0d0* gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & dk2(ihm1,j)* ds1(ihm1,jm1)) - (2.0d0* gam12(ihm1,jhm1))/ 
     & (dg1(ihm1,jhm1)* dk2(ihm1,j)* ds1(ihm1,j)* ds2(im1,jhm1))


       X(-1,0)=
     & (2.0d0* sig11(ihm1,j))/ ds1(ihm1,j)**2 + ((2.0d0* 
     & gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,j)) + (2.0d0* gam21(ihm1,j))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* ds1(ihm1,j)))/ dk2(ihm1,j) - (sig12(ihm1,jhm1)/ 
     & (2.0d0* ds2(im1,jhm1)) - sig12(ihm1,j)/ (2.0d0* 
     & ds2(im1,j)))/ ds1(ihm1,j) + (2.0d0* (gam12(ihm1,jhm1)/ 
     & (dg1(ihm1,jhm1)* ds2(im1,jhm1)) + gam12(ihm1,j)/ 
     & (dg1(ihm1,j)* ds2(im1,j))))/ (dk2(ihm1,j)* ds1(ihm1,j)) - 
     & (gam11(im1,j)* ((1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j))/ 
     & dh1(im1,j) + 1.0d0/ (dh1(im2,j)* ds1(ihm2,j))))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j))

       X(-1,+1)=
     & - sig12(ihm1,j)/ (2.0d0* ds1(ihm1,j)* ds2(im1,j)) - (2.0d0* 
     &   gam21(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* dk2(ihm1,j)* 
     & ds1(ihm1,jp1)) - (2.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* 
     & dk2(ihm1,j)* ds1(ihm1,j)* ds2(im1,j))

       X(0,-2)=
     & -gam22(i,jm1)/ (dh2(i,jm1)* dh2(i,j)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-1)=
     & ((gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)) + (gam22(i,jm1)* 
     & (1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1)))/ dh2(i,jm1))/ 
     & ds2(i,jhm1) + gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* 
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & + sig12(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)* ds2(i,jhm1)) + 
     & (2.0d0* gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & dk2(ihm1,j)* ds1(ihm1,jm1)) + (2.0d0* gam12(ihm1,jhm1))/ 
     & (dg1(ihm1,jhm1)* dk2(ihm1,j)* ds1(ihm1,j)* ds2(i,jhm1))

       X(0,0)=
     & gam11(im1,j)/ (dh1(im1,j)* ds1(ihm2,j)* ds1(ihm1,j)**2) - 
     & ((2.0d0* gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,j)) + (2.0d0* gam21(ihm1,j))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* ds1(ihm1,j)))/ dk2(ihm1,j) - (sig22(i,j-1)/ 
     & ds2(i,jhm1) + sig22(i,j)/ ds2(i,j))/ dh2(i,j) - 

     & (sig12(ihm1,jhm1)/ (2.0d0* ds2(i,jhm1)) - sig12(ihm1,j)/ 
     & (2.0d0* ds2(i,j)))/ ds1(ihm1,j) - ((gam22(i,jm1)/ 
     & (dh2(i,jm1)* ds2(i,jhm1)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) 
     & + 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,jhm1) + (gam22(i,j)* 
     & (1.0d0/ ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* 
     & ds2(i,j)))/ dh2(i,j) - (2.0d0* (gam12(ihm1,jhm1)/ 
     & (dg1(ihm1,jhm1)* ds2(i,jhm1)) + gam12(ihm1,j)/ (dg1(ihm1,j)* 
     & ds2(i,j))))/ (dk2(ihm1,j)* ds1(ihm1,j)) - (2.0d0* 
     & sig11(ihm1,j))/ ds1(ihm1,j)**2

       X(0,+1)=
     & (gam22(i,j)/ (dh2(i,j)* ds2(i,j)**2) + gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j)/ 
     & (dh2(i,j)* ds2(i,j)) - sig12(ihm1,j)/ (2.0d0* ds1(ihm1,j)* 
     & ds2(i,j)) + (2.0d0* gam21(ihm1,j))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* dk2(ihm1,j)* ds1(ihm1,jp1)) + (2.0d0* 
     & gam12(ihm1,j))/ (dg1(ihm1,j)* dk2(ihm1,j)* ds1(ihm1,j)* 
     & ds2(i,j))

!! ================== BNDY 2,1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .eq. 2)) then
       X(-1,0)=
     & (gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)**2) + gam11(i,j)/ 
     & (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + 
     & sig11(ihm1,j)/ (dh1(i,j)* ds1(ihm1,j)) - sig21(ihm1,j)/ 
     & (2.0d0* ds1(ihm1,j)* ds2(i,j)) + (2.0d0* gam12(ihm1,j))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* dk1(i,j)* ds2(im1,j)) + (2.0d0* 
     & gam21(ihm1,j))/ (dg2(ihm1,j)* dk1(i,j)* ds1(ihm1,j)* 
     & ds2(i,j))

       X(-1,+1)=
     & - sig21(ihm1,j)/ (2.0d0* ds1(ihm1,jp1)* ds2(i,j)) - (2.0d0* 
     &   gam12(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* dk1(i,j)* 
     & ds2(im1,j)) - (2.0d0* gam21(ihm1,j))/ (dg2(ihm1,j)* 
     & dk1(i,j)* ds1(ihm1,jp1)* ds2(i,j))

       X(0,0)=
     & (sig21(ihm1,j)/ (2.0d0* ds1(ihm1,j)) - sig21(i,j)/ (2.0d0* 
     & ds1(i,j)))/ ds2(i,j) - ((2.0d0* gam12(ihm1,j))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* ds2(i,j)) + (2.0d0* gam12(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds2(i,j)))/ dk1(i,j) - (sig11(ihm1,j)/ 
     & ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ dh1(i,j) - (2.0d0* 
     & sig22(i,j))/ ds2(i,j)**2 - ((gam11(ip1,j)/ (dh1(ip1,j)* 
     & ds1(i,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 1.0d0/ 
     & ds1(i,j)))/ dh1(i,j))/ ds1(i,j) + (gam11(i,j)* (1.0d0/ 
     & ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* ds1(ihm1,j)))/ 
     & dh1(i,j) - (2.0d0* (gam21(ihm1,j)/ (dg2(ihm1,j)* 
     & ds1(ihm1,j)) + gam21(i,j)/ (dg2(i,j)* ds1(i,j))))/ 
     & (dk1(i,j)* ds2(i,j)) + gam22(i,jp1)/ (dh2(i,jp1)* 
     & ds2(i,j)**2.0d0* ds2(i,jhp1))

       X(0,+1)=
     & (2.0d0* sig22(i,j))/ ds2(i,j)**2 + ((2.0d0* gam12(ihm1,j))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* ds2(i,j)) + (2.0d0* gam12(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds2(i,j)))/ dk1(i,j) + (sig21(ihm1,j)/ 
     & (2.0d0* ds1(ihm1,jp1)) - sig21(i,j)/ (2.0d0* ds1(i,jp1)))/ 
     & ds2(i,j) + (2.0d0* (gam21(ihm1,j)/ (dg2(ihm1,j)* 
     & ds1(ihm1,jp1)) + gam21(i,j)/ (dg2(i,j)* ds1(i,jp1))))/ 
     & (dk1(i,j)* ds2(i,j)) - (gam22(i,jp1)* ((1.0d0/ ds2(i,j) + 
     & 1.0d0/ ds2(i,jhp1))/ dh2(i,jp1) + 1.0d0/ (dh2(i,jp2)* 
     & ds2(i,jhp1))))/ (ds2(i,j)* ds2(i,jhp1))

       X(0,+2)=
     & (gam22(i,jp1)* ((1.0d0/ ds2(i,jhp1) + 1.0d0/ ds2(i,jhp2))/ 
     & dh2(i,jp2) + 1.0d0/ (dh2(i,jp1)* ds2(i,jhp1))))/ (ds2(i,j)* 
     & ds2(i,jhp1))

       X(0,+3)=
     & -gam22(i,jp1)/ (dh2(i,jp2)* ds2(i,j)* ds2(i,jhp1)* 
     & ds2(i,jhp2))

       X(+1,0)=
     & ((gam11(i,j)/ (dh1(i,j)* ds1(i,j)) + (gam11(ip1,j)* (1.0d0/ 
     & ds1(i,j) + 1.0d0/ ds1(ihp1,j)))/ dh1(ip1,j))/ ds1(i,j) + 
     & gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + 
     & sig11(i,j)/ (dh1(i,j)* ds1(i,j)) + sig21(i,j)/ (2.0d0* 
     & ds1(i,j)* ds2(i,j)) + (2.0d0* gam12(i,j))/ (dg1(i,j)* 
     & dg2(i,j)* dk1(i,j)* ds2(ip1,j)) + (2.0d0* gam21(i,j))/ 
     & (dg2(i,j)* dk1(i,j)* ds1(i,j)* ds2(i,j))

       X(+1,+1)=
     & sig21(i,j)/ (2.0d0* ds1(i,jp1)* ds2(i,j)) - (2.0d0* 
     & gam12(i,j))/ (dg1(i,j)* dg2(i,j)* dk1(i,j)* ds2(ip1,j)) - 
     & (2.0d0* gam21(i,j))/ (dg2(i,j)* dk1(i,j)* ds1(i,jp1)* 
     & ds2(i,j))

       X(+2,0)=
     & -gam11(ip1,j)/ (dh1(i,j)* dh1(ip1,j)* ds1(i,j)* ds1(ihp1,j))

!! ================== BNDY Nx-1,1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .eq. nq_ibmin-1)) then
       X(-2,0)=
     & -gam11(im1,j)/ (dh1(im1,j)* dh1(i,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))

       X(-1,0)=
     & ((gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)) + (gam11(im1,j)* 
     & (1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j)))/ dh1(im1,j))/ 
     & ds1(ihm1,j) + gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)* 
     & ds1(i,j)))/ dh1(i,j) + sig11(ihm1,j)/ (dh1(i,j)* 
     & ds1(ihm1,j)) - sig21(ihm1,j)/ (2.0d0* ds1(ihm1,j)* ds2(i,j)) 
     & + (2.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* 
     & dk1(i,j)* ds2(im1,j)) + (2.0d0* gam21(ihm1,j))/ 
     & (dg2(ihm1,j)* dk1(i,j)* ds1(ihm1,j)* ds2(i,j))

       X(-1,+1)=
     & - sig21(ihm1,j)/ (2.0d0* ds1(ihm1,jp1)* ds2(i,j)) - (2.0d0* 
     &   gam12(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* dk1(i,j)* 
     & ds2(im1,j)) - (2.0d0* gam21(ihm1,j))/ (dg2(ihm1,j)* 
     & dk1(i,j)* ds1(ihm1,jp1)* ds2(i,j))

       X(0,0)=
     & (sig21(ihm1,j)/ (2.0d0* ds1(ihm1,j)) - sig21(i,j)/ (2.0d0* 
     & ds1(i,j)))/ ds2(i,j) - ((2.0d0* gam12(ihm1,j))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* ds2(i,j)) + (2.0d0* gam12(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds2(i,j)))/ dk1(i,j) - (sig11(ihm1,j)/ 
     & ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ dh1(i,j) - (2.0d0* 
     & sig22(i,j))/ ds2(i,j)**2 - ((gam11(im1,j)/ (dh1(im1,j)* 
     & ds1(ihm1,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 1.0d0/ 
     & ds1(i,j)))/ dh1(i,j))/ ds1(ihm1,j) + (gam11(i,j)* (1.0d0/ 
     & ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* ds1(i,j)))/ 
     & dh1(i,j) - (2.0d0* (gam21(ihm1,j)/ (dg2(ihm1,j)* 
     & ds1(ihm1,j)) + gam21(i,j)/ (dg2(i,j)* ds1(i,j))))/ 
     & (dk1(i,j)* ds2(i,j)) + gam22(i,jp1)/ (dh2(i,jp1)* 
     & ds2(i,j)**2.0d0* ds2(i,jhp1))

       X(0,+1)=
     & (2.0d0* sig22(i,j))/ ds2(i,j)**2 + ((2.0d0* gam12(ihm1,j))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* ds2(i,j)) + (2.0d0* gam12(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds2(i,j)))/ dk1(i,j) + (sig21(ihm1,j)/ 
     & (2.0d0* ds1(ihm1,jp1)) - sig21(i,j)/ (2.0d0* ds1(i,jp1)))/ 
     & ds2(i,j) + (2.0d0* (gam21(ihm1,j)/ (dg2(ihm1,j)* 
     & ds1(ihm1,jp1)) + gam21(i,j)/ (dg2(i,j)* ds1(i,jp1))))/ 
     & (dk1(i,j)* ds2(i,j)) - (gam22(i,jp1)* ((1.0d0/ ds2(i,j) + 
     & 1.0d0/ ds2(i,jhp1))/ dh2(i,jp1) + 1.0d0/ (dh2(i,jp2)* 
     & ds2(i,jhp1))))/ (ds2(i,j)* ds2(i,jhp1))

       X(0,+2)=
     & (gam22(i,jp1)* ((1.0d0/ ds2(i,jhp1) + 1.0d0/ ds2(i,jhp2))/ 
     & dh2(i,jp2) + 1.0d0/ (dh2(i,jp1)* ds2(i,jhp1))))/ (ds2(i,j)* 
     & ds2(i,jhp1))

       X(0,+3)=
     & -gam22(i,jp1)/ (dh2(i,jp2)* ds2(i,j)* ds2(i,jhp1)* 
     & ds2(i,jhp2))

       X(+1,0)=
     & (gam11(i,j)/ (dh1(i,j)* ds1(i,j)**2) + gam11(i,j)/ 
     & (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + sig11(i,j)/ 
     & (dh1(i,j)* ds1(i,j)) + sig21(i,j)/ (2.0d0* ds1(i,j)* 
     & ds2(i,j)) + (2.0d0* gam12(i,j))/ (dg1(i,j)* dg2(i,j)* 
     & dk1(i,j)* ds2(ip1,j)) + (2.0d0* gam21(i,j))/ (dg2(i,j)* 
     & dk1(i,j)* ds1(i,j)* ds2(i,j))

       X(+1,+1)=
     & sig21(i,j)/ (2.0d0* ds1(i,jp1)* ds2(i,j)) - (2.0d0* 
     & gam12(i,j))/ (dg1(i,j)* dg2(i,j)* dk1(i,j)* ds2(ip1,j)) - 
     & (2.0d0* gam21(i,j))/ (dg2(i,j)* dk1(i,j)* ds1(i,jp1)* 
     & ds2(i,j))

!! ================== BNDY 2,Ny
       elseif((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (i .eq. 2)) then
       X(-1,-1)=
     & sig21(ihm1,jhm1)/ (2.0d0* ds1(ihm1,jm1)* ds2(i,jhm1)) - 
     & (2.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & dk1(i,jhm1)* ds2(im1,jhm1)) - (2.0d0* gam21(ihm1,jhm1))/ 
     & (dg2(ihm1,jhm1)* dk1(i,jhm1)* ds1(ihm1,jm1)* ds2(i,jhm1))

       X(-1,0)=
     & (gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)**2) + gam11(i,j)/ 
     & (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + 
     & sig11(ihm1,j)/ (dh1(i,j)* ds1(ihm1,j)) + sig21(ihm1,jhm1)/ 
     & (2.0d0* ds1(ihm1,j)* ds2(i,jhm1)) + (2.0d0* 
     & gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & dk1(i,jhm1)* ds2(im1,jhm1)) + (2.0d0* gam21(ihm1,jhm1))/ 
     & (dg2(ihm1,jhm1)* dk1(i,jhm1)* ds1(ihm1,j)* ds2(i,jhm1))

       X(0,-3)=
     & -gam22(i,jm1)/ (dh2(i,jm2)* ds2(i,jhm3)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-2)=
     & (gam22(i,jm1)* ((1.0d0/ ds2(i,jhm3) + 1.0d0/ ds2(i,jhm2))/ 
     & dh2(i,jm2) + 1.0d0/ (dh2(i,jm1)* ds2(i,jhm2))))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1))

       X(0,-1)=
     & (2.0d0* sig22(i,j-1))/ ds2(i,jhm1)**2 + ((2.0d0* 
     & gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds2(i,jhm1)) + (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds2(i,jhm1)))/ dk1(i,jhm1) - (sig21(ihm1,jhm1)/ 
     & (2.0d0* ds1(ihm1,jm1)) - sig21(i,jhm1)/ (2.0d0* 
     & ds1(i,jm1)))/ ds2(i,jhm1) + (2.0d0* (gam21(ihm1,jhm1)/ 
     & (dg2(ihm1,jhm1)* ds1(ihm1,jm1)) + gam21(i,jhm1)/ 
     & (dg2(i,jhm1)* ds1(i,jm1))))/ (dk1(i,jhm1)* ds2(i,jhm1)) - 
     & (gam22(i,jm1)* ((1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1))/ 
     & dh2(i,jm1) + 1.0d0/ (dh2(i,jm2)* ds2(i,jhm2))))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1))

       X(0,0)=
     & gam22(i,jm1)/ (dh2(i,jm1)* ds2(i,jhm2)* ds2(i,jhm1)**2) - 
     & ((2.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds2(i,jhm1)) + (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds2(i,jhm1)))/ dk1(i,jhm1) - (sig11(ihm1,j)/ 
     & ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ dh1(i,j) - 
     & (sig21(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)) - sig21(i,jhm1)/ 
     & (2.0d0* ds1(i,j)))/ ds2(i,jhm1) - ((gam11(ip1,j)/ 
     & (dh1(ip1,j)* ds1(i,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 
     & 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(i,j) + (gam11(i,j)* 
     & (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* 
     & ds1(ihm1,j)))/ dh1(i,j) - (2.0d0* (gam21(ihm1,jhm1)/ 
     & (dg2(ihm1,jhm1)* ds1(ihm1,j)) + gam21(i,jhm1)/ (dg2(i,jhm1)* 
     & ds1(i,j))))/ (dk1(i,jhm1)* ds2(i,jhm1)) - (2.0d0* 
     & sig22(i,j-1))/ ds2(i,jhm1)**2

       X(+1,-1)=
     & - sig21(i,jhm1)/ (2.0d0* ds1(i,jm1)* ds2(i,jhm1)) - (2.0d0* 
     &   gam12(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* dk1(i,jhm1)* 
     & ds2(ip1,jhm1)) - (2.0d0* gam21(i,jhm1))/ (dg2(i,jhm1)* 
     & dk1(i,jhm1)* ds1(i,jm1)* ds2(i,jhm1))

       X(+1,0)=
     & ((gam11(i,j)/ (dh1(i,j)* ds1(i,j)) + (gam11(ip1,j)* (1.0d0/ 
     & ds1(i,j) + 1.0d0/ ds1(ihp1,j)))/ dh1(ip1,j))/ ds1(i,j) + 
     & gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + 
     & sig11(i,j)/ (dh1(i,j)* ds1(i,j)) - sig21(i,jhm1)/ (2.0d0* 
     & ds1(i,j)* ds2(i,jhm1)) + (2.0d0* gam12(i,jhm1))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* dk1(i,jhm1)* ds2(ip1,jhm1)) + 
     & (2.0d0* gam21(i,jhm1))/ (dg2(i,jhm1)* dk1(i,jhm1)* ds1(i,j)* 
     & ds2(i,jhm1))

       X(+2,0)=
     & -gam11(ip1,j)/ (dh1(i,j)* dh1(ip1,j)* ds1(i,j)* ds1(ihp1,j))

!! ================== BNDY Nx-1,Ny
       elseif((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (i .eq. nq_ibmin-1)) then

       X(-2,0)=
     & -gam11(im1,j)/ (dh1(im1,j)* dh1(i,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))

       X(-1,-1)=
     & sig21(ihm1,jhm1)/ (2.0d0* ds1(ihm1,jm1)* ds2(i,jhm1)) - 
     & (2.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & dk1(i,jhm1)* ds2(im1,jhm1)) - (2.0d0* gam21(ihm1,jhm1))/ 
     & (dg2(ihm1,jhm1)* dk1(i,jhm1)* ds1(ihm1,jm1)* ds2(i,jhm1))

       X(-1,0)=
     & ((gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)) + (gam11(im1,j)* 
     & (1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j)))/ dh1(im1,j))/ 
     & ds1(ihm1,j) + gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)* 
     & ds1(i,j)))/ dh1(i,j) + sig11(ihm1,j)/ (dh1(i,j)* 
     & ds1(ihm1,j)) + sig21(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)* 
     & ds2(i,jhm1)) + (2.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* 
     & dg2(ihm1,jhm1)* dk1(i,jhm1)* ds2(im1,jhm1)) + (2.0d0* 
     & gam21(ihm1,jhm1))/ (dg2(ihm1,jhm1)* dk1(i,jhm1)* 
     & ds1(ihm1,j)* ds2(i,jhm1))

       X(0,-3)=
     & -gam22(i,jm1)/ (dh2(i,jm2)* ds2(i,jhm3)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-2)=
     & (gam22(i,jm1)* ((1.0d0/ ds2(i,jhm3) + 1.0d0/ ds2(i,jhm2))/ 
     & dh2(i,jm2) + 1.0d0/ (dh2(i,jm1)* ds2(i,jhm2))))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1))

       X(0,-1)=
     & (2.0d0* sig22(i,j-1))/ ds2(i,jhm1)**2 + ((2.0d0* 
     & gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds2(i,jhm1)) + (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds2(i,jhm1)))/ dk1(i,jhm1) - (sig21(ihm1,jhm1)/ 
     & (2.0d0* ds1(ihm1,jm1)) - sig21(i,jhm1)/ (2.0d0* 
     & ds1(i,jm1)))/ ds2(i,jhm1) + (2.0d0* (gam21(ihm1,jhm1)/ 
     & (dg2(ihm1,jhm1)* ds1(ihm1,jm1)) + gam21(i,jhm1)/ 
     & (dg2(i,jhm1)* ds1(i,jm1))))/ (dk1(i,jhm1)* ds2(i,jhm1)) - 
     & (gam22(i,jm1)* ((1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1))/ 
     & dh2(i,jm1) + 1.0d0/ (dh2(i,jm2)* ds2(i,jhm2))))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1))

       X(0,0)=
     & gam22(i,jm1)/ (dh2(i,jm1)* ds2(i,jhm2)* ds2(i,jhm1)**2) - 
     & ((2.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds2(i,jhm1)) + (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds2(i,jhm1)))/ dk1(i,jhm1) - (sig11(ihm1,j)/ 
     & ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ dh1(i,j) - 
     & (sig21(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)) - sig21(i,jhm1)/ 
     & (2.0d0* ds1(i,j)))/ ds2(i,jhm1) - ((gam11(im1,j)/ 
     & (dh1(im1,j)* ds1(ihm1,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) 
     & + 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(ihm1,j) + (gam11(i,j)* 
     & (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* 
     & ds1(i,j)))/ dh1(i,j) - (2.0d0* (gam21(ihm1,jhm1)/ 
     & (dg2(ihm1,jhm1)* ds1(ihm1,j)) + gam21(i,jhm1)/ (dg2(i,jhm1)* 
     & ds1(i,j))))/ (dk1(i,jhm1)* ds2(i,jhm1)) - (2.0d0* 
     & sig22(i,j-1))/ ds2(i,jhm1)**2

       X(+1,-1)=
     & - sig21(i,jhm1)/ (2.0d0* ds1(i,jm1)* ds2(i,jhm1)) - (2.0d0* 
     &   gam12(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* dk1(i,jhm1)* 
     & ds2(ip1,jhm1)) - (2.0d0* gam21(i,jhm1))/ (dg2(i,jhm1)* 
     & dk1(i,jhm1)* ds1(i,jm1)* ds2(i,jhm1))

       X(+1,0)=
     & (gam11(i,j)/ (dh1(i,j)* ds1(i,j)**2) + gam11(i,j)/ 
     & (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + sig11(i,j)/ 
     & (dh1(i,j)* ds1(i,j)) - sig21(i,jhm1)/ (2.0d0* ds1(i,j)* 
     & ds2(i,jhm1)) + (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* dk1(i,jhm1)* ds2(ip1,jhm1)) + (2.0d0* 
     & gam21(i,jhm1))/ (dg2(i,jhm1)* dk1(i,jhm1)* ds1(i,j)* 
     & ds2(i,jhm1))

!! ================== BNDY 2,2
       elseif( 
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
     &       )  then
 
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
     & (gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)**2) + gam11(i,j)/ 
     & (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + 
     & sig11(ihm1,j)/ (dh1(i,j)* ds1(ihm1,j)) + (gam12(ihm1,jhm1)/ 
     & (dg1(ihm1,jhm1)* ds2(im1,jhm1)) + gam12(ihm1,j)/ 
     & (dg1(ihm1,j)* ds2(im1,j)))/ (dh1(i,j)* dk2(ihm1,j))

       X(-1,+1)=
     & - sig12(ihm1,j)/ (4.0d0* dh1(i,j)* ds2(im1,j)) - 
     &   sig21(ihm1,j)/ (4.0d0* dh2(i,j)* ds1(ihm1,jp1)) - 
     & gam12(ihm1,j)/ (dg1(ihm1,j)* dh1(i,j)* dk2(ihm1,j)* 
     & ds2(im1,j)) - gam21(ihm1,j)/ (dg2(ihm1,j)* dh2(i,j)* 
     & dk1(i,j)* ds1(ihm1,jp1))

       X(0,-1)=
     & (gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* dk2(ihm1,j)* 
     & ds2(i,jhm1)) + gam12(i,jhm1)/ (dg1(i,jhm1)* dk2(i,j)* 
     & ds2(i,jhm1)))/ dh1(i,j) + (sig12(ihm1,jhm1)/ (4.0d0* 
     & ds2(i,jhm1)) - sig12(i,jhm1)/ (4.0d0* ds2(i,jhm1)))/ 
     & dh1(i,j) - (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,jm1)) - 
     & sig21(i,jhm1)/ (4.0d0* ds1(i,jm1)))/ dh2(i,j) + (gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)**2) + gam22(i,j)/ (dh2(i,j)* 
     & ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* 
     & ds2(i,jhm1)) + (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* 
     & ds1(ihm1,jm1)) + gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,jm1)))/ 
     & (dh2(i,j)* dk1(i,jhm1))

       X(0,0)=
     & - (sig11(ihm1,j)/ ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ 
     &   dh1(i,j) - (sig22(i,j-1)/ ds2(i,jhm1) + sig22(i,j)/ 
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
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) - ((gam11(ip1,j)/ 
     & (dh1(ip1,j)* ds1(i,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 
     & 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(i,j) + (gam11(i,j)* 
     & (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* 
     & ds1(ihm1,j)))/ dh1(i,j) - ((gam22(i,jp1)/ (dh2(i,jp1)* 
     & ds2(i,j)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 1.0d0/ 
     & ds2(i,j)))/ dh2(i,j))/ ds2(i,j) + (gam22(i,j)* (1.0d0/ 
     & ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* ds2(i,jhm1)))/ 
     & dh2(i,j)

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

!! ================== BNDY Nx-1,2
       elseif( 
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
     &       )  then

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

       X(0,-1)=
     & (gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* dk2(ihm1,j)* 
     & ds2(i,jhm1)) + gam12(i,jhm1)/ (dg1(i,jhm1)* dk2(i,j)* 
     & ds2(i,jhm1)))/ dh1(i,j) + (sig12(ihm1,jhm1)/ (4.0d0* 
     & ds2(i,jhm1)) - sig12(i,jhm1)/ (4.0d0* ds2(i,jhm1)))/ 
     & dh1(i,j) - (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,jm1)) - 
     & sig21(i,jhm1)/ (4.0d0* ds1(i,jm1)))/ dh2(i,j) + (gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)**2) + gam22(i,j)/ (dh2(i,j)* 
     & ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* 
     & ds2(i,jhm1)) + (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* 
     & ds1(ihm1,jm1)) + gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,jm1)))/ 
     & (dh2(i,j)* dk1(i,jhm1))

       X(0,0)=
     & - (sig11(ihm1,j)/ ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ 
     &   dh1(i,j) - (sig22(i,j-1)/ ds2(i,jhm1) + sig22(i,j)/ 
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
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) - ((gam11(im1,j)/ 
     & (dh1(im1,j)* ds1(ihm1,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) 
     & + 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(ihm1,j) + (gam11(i,j)* 
     & (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* 
     & ds1(i,j)))/ dh1(i,j) - ((gam22(i,jp1)/ (dh2(i,jp1)* 
     & ds2(i,j)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 1.0d0/ 
     & ds2(i,j)))/ dh2(i,j))/ ds2(i,j) + (gam22(i,j)* (1.0d0/ 
     & ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* ds2(i,jhm1)))/ 
     & dh2(i,j)

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
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) + (gam11(i,j)/ 
     & (dh1(i,j)* ds1(i,j)**2) + gam11(i,j)/ (dh1(i,j)* 
     & ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + sig11(i,j)/ (dh1(i,j)* 
     & ds1(i,j)) + (gam12(i,jhm1)/ (dg1(i,jhm1)* ds2(ip1,jhm1)) + 
     & gam12(i,j)/ (dg1(i,j)* ds2(ip1,j)))/ (dh1(i,j)* dk2(i,j))

       X(+1,+1)=
     & sig12(i,j)/ (4.0d0* dh1(i,j)* ds2(ip1,j)) + sig21(i,j)/ 
     & (4.0d0* dh2(i,j)* ds1(i,jp1)) - gam12(i,j)/ (dg1(i,j)* 
     & dh1(i,j)* dk2(i,j)* ds2(ip1,j)) - gam21(i,j)/ (dg2(i,j)* 
     & dh2(i,j)* dk1(i,j)* ds1(i,jp1))

!! ================== BNDY 2,Ny-1
       elseif( 
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
     &       )  then
 
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
     & (gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)**2) + gam11(i,j)/ 
     & (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + 
     & sig11(ihm1,j)/ (dh1(i,j)* ds1(ihm1,j)) + (gam12(ihm1,jhm1)/ 
     & (dg1(ihm1,jhm1)* ds2(im1,jhm1)) + gam12(ihm1,j)/ 
     & (dg1(ihm1,j)* ds2(im1,j)))/ (dh1(i,j)* dk2(ihm1,j))

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
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & + (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,jm1)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,jm1)))/ (dh2(i,j)* 
     & dk1(i,jhm1))

       X(0,0)=
     & - (sig11(ihm1,j)/ ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ 
     &   dh1(i,j) - (sig22(i,j-1)/ ds2(i,jhm1) + sig22(i,j)/ 
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
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) - ((gam11(ip1,j)/ 
     & (dh1(ip1,j)* ds1(i,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 
     & 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(i,j) + (gam11(i,j)* 
     & (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* 
     & ds1(ihm1,j)))/ dh1(i,j) - ((gam22(i,jm1)/ (dh2(i,jm1)* 
     & ds2(i,jhm1)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 1.0d0/ 
     & ds2(i,j)))/ dh2(i,j))/ ds2(i,jhm1) + (gam22(i,j)* (1.0d0/ 
     & ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* ds2(i,j)))/ 
     & dh2(i,j)

       X(0,+1)=
     & (gam12(ihm1,j)/ (dg1(ihm1,j)* dk2(ihm1,j)* ds2(i,j)) + 
     & gam12(i,j)/ (dg1(i,j)* dk2(i,j)* ds2(i,j)))/ dh1(i,j) - 
     & (sig12(ihm1,j)/ (4.0d0* ds2(i,j)) - sig12(i,j)/ (4.0d0* 
     & ds2(i,j)))/ dh1(i,j) + (sig21(ihm1,j)/ (4.0d0* 
     & ds1(ihm1,jp1)) - sig21(i,j)/ (4.0d0* ds1(i,jp1)))/ dh2(i,j) 
     & + (gam22(i,j)/ (dh2(i,j)* ds2(i,j)**2) + gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j)/ 
     & (dh2(i,j)* ds2(i,j)) + (gam21(ihm1,j)/ (dg2(ihm1,j)* 
     & ds1(ihm1,jp1)) + gam21(i,j)/ (dg2(i,j)* ds1(i,jp1)))/ 
     & (dh2(i,j)* dk1(i,j))

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

!! ================== BNDY Nx-1,Ny-1
       elseif( 
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
     &       )  then
 
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
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & + (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,jm1)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,jm1)))/ (dh2(i,j)* 
     & dk1(i,jhm1))

       X(0,0)=
     & - (sig11(ihm1,j)/ ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ 
     &   dh1(i,j) - (sig22(i,j-1)/ ds2(i,jhm1) + sig22(i,j)/ 
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
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) - ((gam11(im1,j)/ 
     & (dh1(im1,j)* ds1(ihm1,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) 
     & + 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(ihm1,j) + (gam11(i,j)* 
     & (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* 
     & ds1(i,j)))/ dh1(i,j) - ((gam22(i,jm1)/ (dh2(i,jm1)* 
     & ds2(i,jhm1)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 1.0d0/ 
     & ds2(i,j)))/ dh2(i,j))/ ds2(i,jhm1) + (gam22(i,j)* (1.0d0/ 
     & ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* ds2(i,j)))/ 
     & dh2(i,j)

       X(0,+1)=
     & (gam12(ihm1,j)/ (dg1(ihm1,j)* dk2(ihm1,j)* ds2(i,j)) + 
     & gam12(i,j)/ (dg1(i,j)* dk2(i,j)* ds2(i,j)))/ dh1(i,j) - 
     & (sig12(ihm1,j)/ (4.0d0* ds2(i,j)) - sig12(i,j)/ (4.0d0* 
     & ds2(i,j)))/ dh1(i,j) + (sig21(ihm1,j)/ (4.0d0* 
     & ds1(ihm1,jp1)) - sig21(i,j)/ (4.0d0* ds1(i,jp1)))/ dh2(i,j) 
     & + (gam22(i,j)/ (dh2(i,j)* ds2(i,j)**2) + gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j)/ 
     & (dh2(i,j)* ds2(i,j)) + (gam21(ihm1,j)/ (dg2(ihm1,j)* 
     & ds1(ihm1,jp1)) + gam21(i,j)/ (dg2(i,j)* ds1(i,jp1)))/ 
     & (dh2(i,j)* dk1(i,j))

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
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) + (gam11(i,j)/ 
     & (dh1(i,j)* ds1(i,j)**2) + gam11(i,j)/ (dh1(i,j)* 
     & ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + sig11(i,j)/ (dh1(i,j)* 
     & ds1(i,j)) + (gam12(i,jhm1)/ (dg1(i,jhm1)* ds2(ip1,jhm1)) + 
     & gam12(i,j)/ (dg1(i,j)* ds2(ip1,j)))/ (dh1(i,j)* dk2(i,j))

       X(+1,+1)=
     & sig12(i,j)/ (4.0d0* dh1(i,j)* ds2(ip1,j)) + sig21(i,j)/ 
     & (4.0d0* dh2(i,j)* ds1(i,jp1)) - gam12(i,j)/ (dg1(i,j)* 
     & dh1(i,j)* dk2(i,j)* ds2(ip1,j)) - gam21(i,j)/ (dg2(i,j)* 
     & dh2(i,j)* dk1(i,j)* ds1(i,jp1))

!! ================== BNDY 1,j  2<=j<=Ny-1
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .ge. jst) .and. (j.le. jen)) then
 
       X(0,-2)=
     & -gam22(i,jm1)/ (dh2(i,jm1)* dh2(i,j)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-1)=
     & ((gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)) + (gam22(i,jm1)* 
     & (1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1)))/ dh2(i,jm1))/ 
     & ds2(i,jhm1) + gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* 
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & - sig12(i,jhm1)/ (2.0d0* ds1(i,j)* ds2(i,jhm1)) + (2.0d0* 
     & gam21(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* dk2(i,j)* 

     & ds1(i,jm1)) + (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dk2(i,j)* ds1(i,j)* ds2(i,jhm1))

       X(0,0)=
     & (sig12(i,jhm1)/ (2.0d0* ds2(i,jhm1)) - sig12(i,j)/ (2.0d0* 
     & ds2(i,j)))/ ds1(i,j) - (2.0d0* sig11(i,j))/ ds1(i,j)**2 - 
     & ((2.0d0* gam21(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* 
     & ds1(i,j)) + (2.0d0* gam21(i,j))/ (dg1(i,j)* dg2(i,j)* 
     & ds1(i,j)))/ dk2(i,j) - (sig22(i,j-1)/ ds2(i,jhm1) + 
     & sig22(i,j)/ ds2(i,j))/ dh2(i,j) - ((gam22(i,jm1)/ 
     & (dh2(i,jm1)* ds2(i,jhm1)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) 
     & + 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,jhm1) + (gam22(i,jp1)/ 
     & (dh2(i,jp1)* ds2(i,j)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 
     & 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,j))/ dh2(i,j) - (2.0d0* 
     & (gam12(i,jhm1)/ (dg1(i,jhm1)* ds2(i,jhm1)) + gam12(i,j)/ 
     & (dg1(i,j)* ds2(i,j))))/ (dk2(i,j)* ds1(i,j)) + gam11(ip1,j)/ 
     & (dh1(ip1,j)* ds1(i,j)**2.0d0* ds1(ihp1,j))

       X(0,+1)=
     & ((gam22(i,j)/ (dh2(i,j)* ds2(i,j)) + (gam22(i,jp1)* (1.0d0/ 
     & ds2(i,j) + 1.0d0/ ds2(i,jhp1)))/ dh2(i,jp1))/ ds2(i,j) + 
     & gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + 
     & sig22(i,j)/ (dh2(i,j)* ds2(i,j)) + sig12(i,j)/ (2.0d0* 
     & ds1(i,j)* ds2(i,j)) + (2.0d0* gam21(i,j))/ (dg1(i,j)* 
     & dg2(i,j)* dk2(i,j)* ds1(i,jp1)) + (2.0d0* gam12(i,j))/ 
     & (dg1(i,j)* dk2(i,j)* ds1(i,j)* ds2(i,j))

       X(0,+2)=
     & -gam22(i,jp1)/ (dh2(i,j)* dh2(i,jp1)* ds2(i,j)* ds2(i,jhp1))

       X(+1,-1)=
     & - sig12(i,jhm1)/ (2.0d0* ds1(i,j)* ds2(ip1,jhm1)) - (2.0d0* 
     &   gam21(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* dk2(i,j)* 
     & ds1(i,jm1)) - (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dk2(i,j)* ds1(i,j)* ds2(ip1,jhm1))

       X(+1,0)=
     & (2.0d0* sig11(i,j))/ ds1(i,j)**2 + ((2.0d0* gam21(i,jhm1))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,j)) + (2.0d0* gam21(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds1(i,j)))/ dk2(i,j) + (sig12(i,jhm1)/ 
     & (2.0d0* ds2(ip1,jhm1)) - sig12(i,j)/ (2.0d0* ds2(ip1,j)))/ 
     & ds1(i,j) + (2.0d0* (gam12(i,jhm1)/ (dg1(i,jhm1)* 
     & ds2(ip1,jhm1)) + gam12(i,j)/ (dg1(i,j)* ds2(ip1,j))))/ 
     & (dk2(i,j)* ds1(i,j)) - (gam11(ip1,j)* ((1.0d0/ ds1(i,j) + 
     & 1.0d0/ ds1(ihp1,j))/ dh1(ip1,j) + 1.0d0/ (dh1(ip2,j)* 
     & ds1(ihp1,j))))/ (ds1(i,j)* ds1(ihp1,j))

       X(+1,+1)=
     & sig12(i,j)/ (2.0d0* ds1(i,j)* ds2(ip1,j)) - (2.0d0* 
     & gam21(i,j))/ (dg1(i,j)* dg2(i,j)* dk2(i,j)* ds1(i,jp1)) - 
     & (2.0d0* gam12(i,j))/ (dg1(i,j)* dk2(i,j)* ds1(i,j)* 
     & ds2(ip1,j))

       X(+2,0)=
     & (gam11(ip1,j)* ((1.0d0/ ds1(ihp1,j) + 1.0d0/ ds1(ihp2,j))/ 
     & dh1(ip2,j) + 1.0d0/ (dh1(ip1,j)* ds1(ihp1,j))))/ (ds1(i,j)* 
     & ds1(ihp1,j))

       X(+3,0)=
     & -gam11(ip1,j)/ (dh1(ip2,j)* ds1(i,j)* ds1(ihp1,j)* 
     & ds1(ihp2,j))


!! ================== BNDY N,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin)   
     &        .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .ge. jst) .and. (j.le. jen)) then
 
       X(-3,0)=
     & -gam11(im1,j)/ (dh1(im2,j)* ds1(ihm3,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))


       X(-2,0)=
     & (gam11(im1,j)* ((1.0d0/ ds1(ihm3,j) + 1.0d0/ ds1(ihm2,j))/ 
     & dh1(im2,j) + 1.0d0/ (dh1(im1,j)* ds1(ihm2,j))))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j))

       X(-1,-1)=
     & sig12(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)* ds2(im1,jhm1)) - 
     & (2.0d0* gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & dk2(ihm1,j)* ds1(ihm1,jm1)) - (2.0d0* gam12(ihm1,jhm1))/ 
     & (dg1(ihm1,jhm1)* dk2(ihm1,j)* ds1(ihm1,j)* ds2(im1,jhm1))

       X(-1,0)=
     & (2.0d0* sig11(ihm1,j))/ ds1(ihm1,j)**2 + ((2.0d0* 
     & gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,j)) + (2.0d0* gam21(ihm1,j))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* ds1(ihm1,j)))/ dk2(ihm1,j) - (sig12(ihm1,jhm1)/ 
     & (2.0d0* ds2(im1,jhm1)) - sig12(ihm1,j)/ (2.0d0* 
     & ds2(im1,j)))/ ds1(ihm1,j) + (2.0d0* (gam12(ihm1,jhm1)/ 
     & (dg1(ihm1,jhm1)* ds2(im1,jhm1)) + gam12(ihm1,j)/ 
     & (dg1(ihm1,j)* ds2(im1,j))))/ (dk2(ihm1,j)* ds1(ihm1,j)) - 
     & (gam11(im1,j)* ((1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j))/ 
     & dh1(im1,j) + 1.0d0/ (dh1(im2,j)* ds1(ihm2,j))))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j))

       X(-1,+1)=
     & - sig12(ihm1,j)/ (2.0d0* ds1(ihm1,j)* ds2(im1,j)) - (2.0d0* 
     &   gam21(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* dk2(ihm1,j)* 
     & ds1(ihm1,jp1)) - (2.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* 
     & dk2(ihm1,j)* ds1(ihm1,j)* ds2(im1,j))

       X(0,-2)=
     & -gam22(i,jm1)/ (dh2(i,jm1)* dh2(i,j)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-1)=
     & ((gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)) + (gam22(i,jm1)* 
     & (1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1)))/ dh2(i,jm1))/ 
     & ds2(i,jhm1) + gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* 
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & + sig12(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)* ds2(i,jhm1)) + 
     & (2.0d0* gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & dk2(ihm1,j)* ds1(ihm1,jm1)) + (2.0d0* gam12(ihm1,jhm1))/ 
     & (dg1(ihm1,jhm1)* dk2(ihm1,j)* ds1(ihm1,j)* ds2(i,jhm1))

       X(0,0)=
     & gam11(im1,j)/ (dh1(im1,j)* ds1(ihm2,j)* ds1(ihm1,j)**2) - 
     & (2.0d0* sig11(ihm1,j))/ ds1(ihm1,j)**2 - ((2.0d0* 
     & gam21(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds1(ihm1,j)) + (2.0d0* gam21(ihm1,j))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* ds1(ihm1,j)))/ dk2(ihm1,j) - (sig22(i,j-1)/ 
     & ds2(i,jhm1) + sig22(i,j)/ ds2(i,j))/ dh2(i,j) - 
     & (sig12(ihm1,jhm1)/ (2.0d0* ds2(i,jhm1)) - sig12(ihm1,j)/ 
     & (2.0d0* ds2(i,j)))/ ds1(ihm1,j) - (2.0d0* (gam12(ihm1,jhm1)/ 
     & (dg1(ihm1,jhm1)* ds2(i,jhm1)) + gam12(ihm1,j)/ (dg1(ihm1,j)* 
     & ds2(i,j))))/ (dk2(ihm1,j)* ds1(ihm1,j)) - ((gam22(i,jm1)/ 
     & (dh2(i,jm1)* ds2(i,jhm1)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) 
     & + 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,jhm1) + (gam22(i,jp1)/ 
     & (dh2(i,jp1)* ds2(i,j)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 
     & 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,j))/ dh2(i,j)

       X(0,+1)=
     & ((gam22(i,j)/ (dh2(i,j)* ds2(i,j)) + (gam22(i,jp1)* (1.0d0/ 
     & ds2(i,j) + 1.0d0/ ds2(i,jhp1)))/ dh2(i,jp1))/ ds2(i,j) + 
     & gam22(i,j)/ (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + 
     & sig22(i,j)/ (dh2(i,j)* ds2(i,j)) - sig12(ihm1,j)/ (2.0d0* 
     & ds1(ihm1,j)* ds2(i,j)) + (2.0d0* gam21(ihm1,j))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* dk2(ihm1,j)* ds1(ihm1,jp1)) + 
     & (2.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* dk2(ihm1,j)* 
     & ds1(ihm1,j)* ds2(i,j))

       X(0,+2)=
     & -gam22(i,jp1)/ (dh2(i,j)* dh2(i,jp1)* ds2(i,j)* ds2(i,jhp1))


!! ================== BNDY i,1 2<=i<=Nx-1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .ge. ist) .and. (i.le. ien)) then

       X(-2,0)=
     & -gam11(im1,j)/ (dh1(im1,j)* dh1(i,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))

       X(-1,0)=
     & ((gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)) + (gam11(im1,j)* 
     & (1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j)))/ dh1(im1,j))/ 
     & ds1(ihm1,j) + gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)* 
     & ds1(i,j)))/ dh1(i,j) + sig11(ihm1,j)/ (dh1(i,j)* 
     & ds1(ihm1,j)) - sig21(ihm1,j)/ (2.0d0* ds1(ihm1,j)* ds2(i,j)) 
     & + (2.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* 
     & dk1(i,j)* ds2(im1,j)) + (2.0d0* gam21(ihm1,j))/ 
     & (dg2(ihm1,j)* dk1(i,j)* ds1(ihm1,j)* ds2(i,j))

       X(-1,+1)=
     & - sig21(ihm1,j)/ (2.0d0* ds1(ihm1,jp1)* ds2(i,j)) - (2.0d0* 
     &   gam12(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* dk1(i,j)* 
     & ds2(im1,j)) - (2.0d0* gam21(ihm1,j))/ (dg2(ihm1,j)* 
     & dk1(i,j)* ds1(ihm1,jp1)* ds2(i,j))

       X(0,0)=
     & (sig21(ihm1,j)/ (2.0d0* ds1(ihm1,j)) - sig21(i,j)/ (2.0d0* 
     & ds1(i,j)))/ ds2(i,j) - (2.0d0* sig22(i,j))/ ds2(i,j)**2 - 
     & ((2.0d0* gam12(ihm1,j))/ (dg1(ihm1,j)* dg2(ihm1,j)* 
     & ds2(i,j)) + (2.0d0* gam12(i,j))/ (dg1(i,j)* dg2(i,j)* 

     & ds2(i,j)))/ dk1(i,j) - (sig11(ihm1,j)/ ds1(ihm1,j) + 
     & sig11(i,j)/ ds1(i,j))/ dh1(i,j) - ((gam11(im1,j)/ 
     & (dh1(im1,j)* ds1(ihm1,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) 
     & + 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(ihm1,j) + (gam11(ip1,j)/ 
     & (dh1(ip1,j)* ds1(i,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 
     & 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(i,j))/ dh1(i,j) - (2.0d0* 
     & (gam21(ihm1,j)/ (dg2(ihm1,j)* ds1(ihm1,j)) + gam21(i,j)/ 
     & (dg2(i,j)* ds1(i,j))))/ (dk1(i,j)* ds2(i,j)) + gam22(i,jp1)/ 
     & (dh2(i,jp1)* ds2(i,j)**2.0d0* ds2(i,jhp1))

       X(0,+1)=
     & (2.0d0* sig22(i,j))/ ds2(i,j)**2 + ((2.0d0* gam12(ihm1,j))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* ds2(i,j)) + (2.0d0* gam12(i,j))/ 
     & (dg1(i,j)* dg2(i,j)* ds2(i,j)))/ dk1(i,j) + (sig21(ihm1,j)/ 
     & (2.0d0* ds1(ihm1,jp1)) - sig21(i,j)/ (2.0d0* ds1(i,jp1)))/ 
     & ds2(i,j) + (2.0d0* (gam21(ihm1,j)/ (dg2(ihm1,j)* 
     & ds1(ihm1,jp1)) + gam21(i,j)/ (dg2(i,j)* ds1(i,jp1))))/ 
     & (dk1(i,j)* ds2(i,j)) - (gam22(i,jp1)* ((1.0d0/ ds2(i,j) + 
     & 1.0d0/ ds2(i,jhp1))/ dh2(i,jp1) + 1.0d0/ (dh2(i,jp2)* 
     & ds2(i,jhp1))))/ (ds2(i,j)* ds2(i,jhp1))

       X(0,+2)=
     & (gam22(i,jp1)* ((1.0d0/ ds2(i,jhp1) + 1.0d0/ ds2(i,jhp2))/ 
     & dh2(i,jp2) + 1.0d0/ (dh2(i,jp1)* ds2(i,jhp1))))/ (ds2(i,j)* 
     & ds2(i,jhp1))

       X(0,+3)=
     & -gam22(i,jp1)/ (dh2(i,jp2)* ds2(i,j)* ds2(i,jhp1)* 
     & ds2(i,jhp2))

       X(+1,0)=
     & ((gam11(i,j)/ (dh1(i,j)* ds1(i,j)) + (gam11(ip1,j)* (1.0d0/ 
     & ds1(i,j) + 1.0d0/ ds1(ihp1,j)))/ dh1(ip1,j))/ ds1(i,j) + 
     & gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + 
     & sig11(i,j)/ (dh1(i,j)* ds1(i,j)) + sig21(i,j)/ (2.0d0* 
     & ds1(i,j)* ds2(i,j)) + (2.0d0* gam12(i,j))/ (dg1(i,j)* 
     & dg2(i,j)* dk1(i,j)* ds2(ip1,j)) + (2.0d0* gam21(i,j))/ 
     & (dg2(i,j)* dk1(i,j)* ds1(i,j)* ds2(i,j))

       X(+1,+1)=
     & sig21(i,j)/ (2.0d0* ds1(i,jp1)* ds2(i,j)) - (2.0d0* 
     & gam12(i,j))/ (dg1(i,j)* dg2(i,j)* dk1(i,j)* ds2(ip1,j)) - 
     & (2.0d0* gam21(i,j))/ (dg2(i,j)* dk1(i,j)* ds1(i,jp1)* 
     & ds2(i,j))

       X(+2,0)=
     & -gam11(ip1,j)/ (dh1(i,j)* dh1(ip1,j)* ds1(i,j)* ds1(ihp1,j))


!       write(5997,*) '***********i,j =',i,j
!       write(5997,*) jp1,jp2,jp3
!       write(5997,*) X(0,0), X(0,1), X(0,2), X(0,3)
!       write(5997,*) gam22(i,jp1)/ 
!     & (dh2(i,jp1)* ds2(i,j)**2.0d0* ds2(i,jhp1)),- 
!     & (gam22(i,jp1)* ((1.0d0/ ds2(i,j) + 1.0d0/ ds2(i,jhp1))/ 
!     & dh2(i,jp1) + 1.0d0/ (dh2(i,jp2)* ds2(i,jhp1))))/ (ds2(i,j)* 
!     & ds2(i,jhp1)),(gam22(i,jp1)* ((1.0d0/ ds2(i,jhp1) + 1.0d0/ 
!     & ds2(i,jhp2))/ dh2(i,jp2) + 1.0d0/ (dh2(i,jp1)* 
!     & ds2(i,jhp1))))/ (ds2(i,j)* ds2(i,jhp1)),
!     & -gam22(i,jp1)/ (dh2(i,jp2)* ds2(i,j)* ds2(i,jhp1)* 
!     & ds2(i,jhp2))


!! ================== BNDY i,N 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin)   
     &        .and. (BoundaryFibShellin(1,2) .eq. 1)  
     &    .and. (i .ge. ist) .and. (i.le. ien)) then

       X(-2,0)=
     & -gam11(im1,j)/ (dh1(im1,j)* dh1(i,j)* ds1(ihm2,j)* 
     & ds1(ihm1,j))

       X(-1,-1)=
     & sig21(ihm1,jhm1)/ (2.0d0* ds1(ihm1,jm1)* ds2(i,jhm1)) - 
     & (2.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 

     & dk1(i,jhm1)* ds2(im1,jhm1)) - (2.0d0* gam21(ihm1,jhm1))/ 
     & (dg2(ihm1,jhm1)* dk1(i,jhm1)* ds1(ihm1,jm1)* ds2(i,jhm1))

       X(-1,0)=
     & ((gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)) + (gam11(im1,j)* 
     & (1.0d0/ ds1(ihm2,j) + 1.0d0/ ds1(ihm1,j)))/ dh1(im1,j))/ 
     & ds1(ihm1,j) + gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)* 
     & ds1(i,j)))/ dh1(i,j) + sig11(ihm1,j)/ (dh1(i,j)* 
     & ds1(ihm1,j)) + sig21(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)* 
     & ds2(i,jhm1)) + (2.0d0* gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* 
     & dg2(ihm1,jhm1)* dk1(i,jhm1)* ds2(im1,jhm1)) + (2.0d0* 
     & gam21(ihm1,jhm1))/ (dg2(ihm1,jhm1)* dk1(i,jhm1)* 
     & ds1(ihm1,j)* ds2(i,jhm1))

       X(0,-3)=
     & -gam22(i,jm1)/ (dh2(i,jm2)* ds2(i,jhm3)* ds2(i,jhm2)* 
     & ds2(i,jhm1))

       X(0,-2)=
     & (gam22(i,jm1)* ((1.0d0/ ds2(i,jhm3) + 1.0d0/ ds2(i,jhm2))/ 
     & dh2(i,jm2) + 1.0d0/ (dh2(i,jm1)* ds2(i,jhm2))))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1))

       X(0,-1)=
     & (2.0d0* sig22(i,j-1))/ ds2(i,jhm1)**2 + ((2.0d0* 
     & gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds2(i,jhm1)) + (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds2(i,jhm1)))/ dk1(i,jhm1) - (sig21(ihm1,jhm1)/ 
     & (2.0d0* ds1(ihm1,jm1)) - sig21(i,jhm1)/ (2.0d0* 
     & ds1(i,jm1)))/ ds2(i,jhm1) + (2.0d0* (gam21(ihm1,jhm1)/ 
     & (dg2(ihm1,jhm1)* ds1(ihm1,jm1)) + gam21(i,jhm1)/ 
     & (dg2(i,jhm1)* ds1(i,jm1))))/ (dk1(i,jhm1)* ds2(i,jhm1)) - 
     & (gam22(i,jm1)* ((1.0d0/ ds2(i,jhm2) + 1.0d0/ ds2(i,jhm1))/ 
     & dh2(i,jm1) + 1.0d0/ (dh2(i,jm2)* ds2(i,jhm2))))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1))

       X(0,0)=
     & gam22(i,jm1)/ (dh2(i,jm1)* ds2(i,jhm2)* ds2(i,jhm1)**2) - 
     & (2.0d0* sig22(i,j-1))/ ds2(i,jhm1)**2 - ((2.0d0* 
     & gam12(ihm1,jhm1))/ (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* 
     & ds2(i,jhm1)) + (2.0d0* gam12(i,jhm1))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds2(i,jhm1)))/ dk1(i,jhm1) - (sig11(ihm1,j)/ 
     & ds1(ihm1,j) + sig11(i,j)/ ds1(i,j))/ dh1(i,j) - 
     & (sig21(ihm1,jhm1)/ (2.0d0* ds1(ihm1,j)) - sig21(i,jhm1)/ 
     & (2.0d0* ds1(i,j)))/ ds2(i,jhm1) - (2.0d0* (gam21(ihm1,jhm1)/ 
     & (dg2(ihm1,jhm1)* ds1(ihm1,j)) + gam21(i,jhm1)/ (dg2(i,jhm1)* 
     & ds1(i,j))))/ (dk1(i,jhm1)* ds2(i,jhm1)) - ((gam11(im1,j)/ 
     & (dh1(im1,j)* ds1(ihm1,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) 
     & + 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(ihm1,j) + (gam11(ip1,j)/ 
     & (dh1(ip1,j)* ds1(i,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 
     & 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(i,j))/ dh1(i,j)

       X(+1,-1)=
     & - sig21(i,jhm1)/ (2.0d0* ds1(i,jm1)* ds2(i,jhm1)) - (2.0d0* 
     &   gam12(i,jhm1))/ (dg1(i,jhm1)* dg2(i,jhm1)* dk1(i,jhm1)* 
     & ds2(ip1,jhm1)) - (2.0d0* gam21(i,jhm1))/ (dg2(i,jhm1)* 
     & dk1(i,jhm1)* ds1(i,jm1)* ds2(i,jhm1))

       X(+1,0)=
     & ((gam11(i,j)/ (dh1(i,j)* ds1(i,j)) + (gam11(ip1,j)* (1.0d0/ 
     & ds1(i,j) + 1.0d0/ ds1(ihp1,j)))/ dh1(ip1,j))/ ds1(i,j) + 
     & gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + 
     & sig11(i,j)/ (dh1(i,j)* ds1(i,j)) - sig21(i,jhm1)/ (2.0d0* 
     & ds1(i,j)* ds2(i,jhm1)) + (2.0d0* gam12(i,jhm1))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* dk1(i,jhm1)* ds2(ip1,jhm1)) + 
     & (2.0d0* gam21(i,jhm1))/ (dg2(i,jhm1)* dk1(i,jhm1)* ds1(i,j)* 
     & ds2(i,jhm1))

       X(+2,0)=
     & -gam11(ip1,j)/ (dh1(i,j)* dh1(ip1,j)* ds1(i,j)* ds1(ihp1,j))

!! ========================  row 2   2,j 2<=j<=Ny-1
       elseif((i.eq. 2) .and. (
     &                          (BoundaryFibShellin(2,1) .eq. 1)
     &                           .or. 
     &                          (BoundaryFibShellin(2,1) .eq. 2)
     &                         )
     &    .and. (j .ge. jst) .and. (j.le. jen)) then
 

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
     & (gam11(i,j)/ (dh1(i,j)* ds1(ihm1,j)**2) + gam11(i,j)/ 
     & (dh1(i,j)* ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + 
     & sig11(ihm1,j)/ (dh1(i,j)* ds1(ihm1,j)) + (gam12(ihm1,jhm1)/ 
     & (dg1(ihm1,jhm1)* ds2(im1,jhm1)) + gam12(ihm1,j)/ 
     & (dg1(ihm1,j)* ds2(im1,j)))/ (dh1(i,j)* dk2(ihm1,j))

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
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & + (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,jm1)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,jm1)))/ (dh2(i,j)* 
     & dk1(i,jhm1))

       X(0,0)=
     & - ((gam22(i,jm1)/ (dh2(i,jm1)* ds2(i,jhm1)) + (gam22(i,j)* 
     &   (1.0d0/ ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ dh2(i,j))/ 
     & ds2(i,jhm1) + (gam22(i,jp1)/ (dh2(i,jp1)* ds2(i,j)) + 
     & (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ 
     & dh2(i,j))/ ds2(i,j))/ dh2(i,j) - (sig11(ihm1,j)/ ds1(ihm1,j) 
     & + sig11(i,j)/ ds1(i,j))/ dh1(i,j) - (sig22(i,j-1)/ 
     & ds2(i,jhm1) + sig22(i,j)/ ds2(i,j))/ dh2(i,j) - 
     & ((gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* ds2(i,jhm1)) + 
     & gam12(ihm1,j)/ (dg1(ihm1,j)* ds2(i,j)))/ dk2(ihm1,j) + 
     & (gam12(i,jhm1)/ (dg1(i,jhm1)* ds2(i,jhm1)) + gam12(i,j)/ 
     & (dg1(i,j)* ds2(i,j)))/ dk2(i,j))/ dh1(i,j) - 
     & ((gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,j)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,j)))/ dk1(i,jhm1) + 
     & (gam21(ihm1,j)/ (dg2(ihm1,j)* ds1(ihm1,j)) + gam21(i,j)/ 
     & (dg2(i,j)* ds1(i,j)))/ dk1(i,j))/ dh2(i,j) - 
     & (sig12(ihm1,jhm1)/ (4.0d0* ds2(i,jhm1)) - sig12(ihm1,j)/ 
     & (4.0d0* ds2(i,j)) - sig12(i,jhm1)/ (4.0d0* ds2(i,jhm1)) + 
     & sig12(i,j)/ (4.0d0* ds2(i,j)))/ dh1(i,j) - 
     & (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,j)) - sig21(ihm1,j)/ 
     & (4.0d0* ds1(ihm1,j)) - sig21(i,jhm1)/ (4.0d0* ds1(i,j)) + 
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) - ((gam11(ip1,j)/ 
     & (dh1(ip1,j)* ds1(i,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 
     & 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(i,j) + (gam11(i,j)* 
     & (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* 
     & ds1(ihm1,j)))/ dh1(i,j)

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


!! ========================  row 2   Nx-1,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin-1) .and. (
     &                          (BoundaryFibShellin(2,2) .eq. 1)
     &                           .or. 
     &                          (BoundaryFibShellin(2,2) .eq. 2)
     &                         )
     &    .and. (j .ge. jst) .and. (j.le. jen)) then


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
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & + (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,jm1)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,jm1)))/ (dh2(i,j)* 
     & dk1(i,jhm1))

       X(0,0)=
     & - ((gam22(i,jm1)/ (dh2(i,jm1)* ds2(i,jhm1)) + (gam22(i,j)* 
     &   (1.0d0/ ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ dh2(i,j))/ 
     & ds2(i,jhm1) + (gam22(i,jp1)/ (dh2(i,jp1)* ds2(i,j)) + 
     & (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ 
     & dh2(i,j))/ ds2(i,j))/ dh2(i,j) - (sig11(ihm1,j)/ ds1(ihm1,j) 
     & + sig11(i,j)/ ds1(i,j))/ dh1(i,j) - (sig22(i,j-1)/ 
     & ds2(i,jhm1) + sig22(i,j)/ ds2(i,j))/ dh2(i,j) - 
     & ((gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* ds2(i,jhm1)) + 
     & gam12(ihm1,j)/ (dg1(ihm1,j)* ds2(i,j)))/ dk2(ihm1,j) + 
     & (gam12(i,jhm1)/ (dg1(i,jhm1)* ds2(i,jhm1)) + gam12(i,j)/ 
     & (dg1(i,j)* ds2(i,j)))/ dk2(i,j))/ dh1(i,j) - 

     & ((gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,j)) + 

     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,j)))/ dk1(i,jhm1) + 
     & (gam21(ihm1,j)/ (dg2(ihm1,j)* ds1(ihm1,j)) + gam21(i,j)/ 
     & (dg2(i,j)* ds1(i,j)))/ dk1(i,j))/ dh2(i,j) - 
     & (sig12(ihm1,jhm1)/ (4.0d0* ds2(i,jhm1)) - sig12(ihm1,j)/ 
     & (4.0d0* ds2(i,j)) - sig12(i,jhm1)/ (4.0d0* ds2(i,jhm1)) + 
     & sig12(i,j)/ (4.0d0* ds2(i,j)))/ dh1(i,j) - 
     & (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,j)) - sig21(ihm1,j)/ 
     & (4.0d0* ds1(ihm1,j)) - sig21(i,jhm1)/ (4.0d0* ds1(i,j)) + 
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) - ((gam11(im1,j)/ 
     & (dh1(im1,j)* ds1(ihm1,j)) + (gam11(i,j)* (1.0d0/ ds1(ihm1,j) 
     & + 1.0d0/ ds1(i,j)))/ dh1(i,j))/ ds1(ihm1,j) + (gam11(i,j)* 
     & (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ (dh1(i,j)* 

     & ds1(i,j)))/ dh1(i,j)

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
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) + (gam11(i,j)/ 
     & (dh1(i,j)* ds1(i,j)**2) + gam11(i,j)/ (dh1(i,j)* 
     & ds1(ihm1,j)* ds1(i,j)))/ dh1(i,j) + sig11(i,j)/ (dh1(i,j)* 
     & ds1(i,j)) + (gam12(i,jhm1)/ (dg1(i,jhm1)* ds2(ip1,jhm1)) + 
     & gam12(i,j)/ (dg1(i,j)* ds2(ip1,j)))/ (dh1(i,j)* dk2(i,j))

       X(+1,+1)=
     & sig12(i,j)/ (4.0d0* dh1(i,j)* ds2(ip1,j)) + sig21(i,j)/ 
     & (4.0d0* dh2(i,j)* ds1(i,jp1)) - gam12(i,j)/ (dg1(i,j)* 
     & dh1(i,j)* dk2(i,j)* ds2(ip1,j)) - gam21(i,j)/ (dg2(i,j)* 
     & dh2(i,j)* dk1(i,j)* ds1(i,jp1))

!! ========================  row 2   i,2 2<=i<=Nx-1
       elseif((j.eq. 2) .and. (
     &                          (BoundaryFibShellin(1,1) .eq. 1)
     &                           .or. 
     &                          (BoundaryFibShellin(1,1) .eq. 2)
     &                         )
     &    .and. (i .ge. ist) .and. (i.le. ien)) then

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

       X(0,-1)=
     & (gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* dk2(ihm1,j)* 
     & ds2(i,jhm1)) + gam12(i,jhm1)/ (dg1(i,jhm1)* dk2(i,j)* 
     & ds2(i,jhm1)))/ dh1(i,j) + (sig12(ihm1,jhm1)/ (4.0d0* 
     & ds2(i,jhm1)) - sig12(i,jhm1)/ (4.0d0* ds2(i,jhm1)))/ 
     & dh1(i,j) - (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,jm1)) - 
     & sig21(i,jhm1)/ (4.0d0* ds1(i,jm1)))/ dh2(i,j) + (gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)**2) + gam22(i,j)/ (dh2(i,j)* 
     & ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* 
     & ds2(i,jhm1)) + (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* 
     & ds1(ihm1,jm1)) + gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,jm1)))/ 
     & (dh2(i,j)* dk1(i,jhm1))

       X(0,0)=
     & - ((gam11(im1,j)/ (dh1(im1,j)* ds1(ihm1,j)) + (gam11(i,j)* 
     &   (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ dh1(i,j))/ 
     & ds1(ihm1,j) + (gam11(ip1,j)/ (dh1(ip1,j)* ds1(i,j)) + 
     & (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ 
     & dh1(i,j))/ ds1(i,j))/ dh1(i,j) - (sig11(ihm1,j)/ ds1(ihm1,j) 
     & + sig11(i,j)/ ds1(i,j))/ dh1(i,j) - (sig22(i,j-1)/ 
     & ds2(i,jhm1) + sig22(i,j)/ ds2(i,j))/ dh2(i,j) - 
     & ((gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* ds2(i,jhm1)) + 
     & gam12(ihm1,j)/ (dg1(ihm1,j)* ds2(i,j)))/ dk2(ihm1,j) + 
     & (gam12(i,jhm1)/ (dg1(i,jhm1)* ds2(i,jhm1)) + gam12(i,j)/ 
     & (dg1(i,j)* ds2(i,j)))/ dk2(i,j))/ dh1(i,j) - 
     & ((gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,j)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,j)))/ dk1(i,jhm1) + 
     & (gam21(ihm1,j)/ (dg2(ihm1,j)* ds1(ihm1,j)) + gam21(i,j)/ 
     & (dg2(i,j)* ds1(i,j)))/ dk1(i,j))/ dh2(i,j) - 
     & (sig12(ihm1,jhm1)/ (4.0d0* ds2(i,jhm1)) - sig12(ihm1,j)/ 
     & (4.0d0* ds2(i,j)) - sig12(i,jhm1)/ (4.0d0* ds2(i,jhm1)) + 
     & sig12(i,j)/ (4.0d0* ds2(i,j)))/ dh1(i,j) - 
     & (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,j)) - sig21(ihm1,j)/ 
     & (4.0d0* ds1(ihm1,j)) - sig21(i,jhm1)/ (4.0d0* ds1(i,j)) + 
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) - ((gam22(i,jp1)/ 
     & (dh2(i,jp1)* ds2(i,j)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) + 
     & 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,j) + (gam22(i,j)* 
     & (1.0d0/ ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* 
     & ds2(i,jhm1)))/ dh2(i,j)

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

!! ========================  row 2   i,Ny 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin-1) .and. (
     &                          (BoundaryFibShellin(1,2) .eq. 1)
     &                           .or. 
     &                          (BoundaryFibShellin(1,2) .eq. 2)
     &                         )
     &    .and. (i .ge. ist) .and. (i.le. ien)) then

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
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
     & + (gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,jm1)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,jm1)))/ (dh2(i,j)* 
     & dk1(i,jhm1))

       X(0,0)=
     & - ((gam11(im1,j)/ (dh1(im1,j)* ds1(ihm1,j)) + (gam11(i,j)* 
     &   (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ dh1(i,j))/ 
     & ds1(ihm1,j) + (gam11(ip1,j)/ (dh1(ip1,j)* ds1(i,j)) + 
     & (gam11(i,j)* (1.0d0/ ds1(ihm1,j) + 1.0d0/ ds1(i,j)))/ 
     & dh1(i,j))/ ds1(i,j))/ dh1(i,j) - (sig11(ihm1,j)/ ds1(ihm1,j) 
     & + sig11(i,j)/ ds1(i,j))/ dh1(i,j) - (sig22(i,j-1)/ 
     & ds2(i,jhm1) + sig22(i,j)/ ds2(i,j))/ dh2(i,j) - 
     & ((gam12(ihm1,jhm1)/ (dg1(ihm1,jhm1)* ds2(i,jhm1)) + 
     & gam12(ihm1,j)/ (dg1(ihm1,j)* ds2(i,j)))/ dk2(ihm1,j) + 
     & (gam12(i,jhm1)/ (dg1(i,jhm1)* ds2(i,jhm1)) + gam12(i,j)/ 
     & (dg1(i,j)* ds2(i,j)))/ dk2(i,j))/ dh1(i,j) - 
     & ((gam21(ihm1,jhm1)/ (dg2(ihm1,jhm1)* ds1(ihm1,j)) + 
     & gam21(i,jhm1)/ (dg2(i,jhm1)* ds1(i,j)))/ dk1(i,jhm1) + 
     & (gam21(ihm1,j)/ (dg2(ihm1,j)* ds1(ihm1,j)) + gam21(i,j)/ 
     & (dg2(i,j)* ds1(i,j)))/ dk1(i,j))/ dh2(i,j) - 
     & (sig12(ihm1,jhm1)/ (4.0d0* ds2(i,jhm1)) - sig12(ihm1,j)/ 
     & (4.0d0* ds2(i,j)) - sig12(i,jhm1)/ (4.0d0* ds2(i,jhm1)) + 
     & sig12(i,j)/ (4.0d0* ds2(i,j)))/ dh1(i,j) - 
     & (sig21(ihm1,jhm1)/ (4.0d0* ds1(ihm1,j)) - sig21(ihm1,j)/ 
     & (4.0d0* ds1(ihm1,j)) - sig21(i,jhm1)/ (4.0d0* ds1(i,j)) + 
     & sig21(i,j)/ (4.0d0* ds1(i,j)))/ dh2(i,j) - ((gam22(i,jm1)/ 
     & (dh2(i,jm1)* ds2(i,jhm1)) + (gam22(i,j)* (1.0d0/ ds2(i,jhm1) 
     & + 1.0d0/ ds2(i,j)))/ dh2(i,j))/ ds2(i,jhm1) + (gam22(i,j)* 
     & (1.0d0/ ds2(i,jhm1) + 1.0d0/ ds2(i,j)))/ (dh2(i,j)* 
     & ds2(i,j)))/ dh2(i,j)

       X(0,+1)=
     & (gam12(ihm1,j)/ (dg1(ihm1,j)* dk2(ihm1,j)* ds2(i,j)) + 
     & gam12(i,j)/ (dg1(i,j)* dk2(i,j)* ds2(i,j)))/ dh1(i,j) - 
     & (sig12(ihm1,j)/ (4.0d0* ds2(i,j)) - sig12(i,j)/ (4.0d0* 
     & ds2(i,j)))/ dh1(i,j) + (sig21(ihm1,j)/ (4.0d0* 
     & ds1(ihm1,jp1)) - sig21(i,j)/ (4.0d0* ds1(i,jp1)))/ dh2(i,j) 
     & + (gam22(i,j)/ (dh2(i,j)* ds2(i,j)**2) + gam22(i,j)/ 
     & (dh2(i,j)* ds2(i,jhm1)* ds2(i,j)))/ dh2(i,j) + sig22(i,j)/ 
     & (dh2(i,j)* ds2(i,j)) + (gam21(ihm1,j)/ (dg2(ihm1,j)* 
     & ds1(ihm1,jp1)) + gam21(i,j)/ (dg2(i,j)* ds1(i,jp1)))/ 
     & (dh2(i,j)* dk1(i,j))

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

!other points
       else

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
     & ds2(i,j)))/ dh2(i,j) + sig22(i,j-1)/ (dh2(i,j)* ds2(i,jhm1)) 
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
     & dh1(i,j) - (sig22(i,j-1)/ ds2(i,jhm1) + sig22(i,j)/ 
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

       endif
        do i2=-3,3
          do j2=-3,3
            if(abs(X(i2,j2)).gt. 0.d0) then
!              if((i2.eq.0).and.(j2.eq.0)) then  !TOBE MODIFIED OTHER PLACES 
                do k=1,3
                  nojtmp=(k-1)*nsnq+(i-1)*ns_ibmin+j
                  ForceOut(nojtmp)=ForceOut(nojtmp)-X(i2,j2)*
     &            (GX(k,i+i2,j+j2)+dGXdt(k,i+i2,j+j2))
!              else
!                  ForceOut(nojtmp)=ForceOut(nojtmp)+X(i2,j2)*
!     &            GX(k,i+i2,j+j2)
!              endif
               enddo
            endif
            X(i2,j2)=-X(i2,j2)*alpha_coefK(1)
          enddo
        enddo
       X(0,0)=X(0,0)+X00*scal_K

       Y=0.0       !KSHTEST 20170524-8:17S

!===nonlinear parts

       if((j .le. 2).and. (BoundaryFibShellin(1,1) .eq. 0))then
          Y(1:3,1:3,-3:3,-3:3)=0.0
       elseif((j .ge. ns_ibmin-1)
     &       .and. (BoundaryFibShellin(1,2) .eq. 0))then
          Y(1:3,1:3,-3:3,-3:3)=0.0
       elseif((i .le. 2).and. (BoundaryFibShellin(2,1) .eq. 0))then
          Y(1:3,1:3,-3:3,-3:3)=0.0
       elseif((i .ge. nq_ibmin-1)
     &       .and. (BoundaryFibShellin(2,2) .eq. 0))then
          Y(1:3,1:3,-3:3,-3:3)=0.0

       elseif((j .le. 1).and. (BoundaryFibShellin(1,1) .eq. 2))then
          Y(1:3,1:3,-3:3,-3:3)=0.0
       elseif((j .ge. ns_ibmin)
     &       .and. (BoundaryFibShellin(1,2) .eq. 2))then
          Y(1:3,1:3,-3:3,-3:3)=0.0
       elseif((i .le. 1).and. (BoundaryFibShellin(2,1) .eq. 2))then
          Y(1:3,1:3,-3:3,-3:3)=0.0
       elseif((i .ge. nq_ibmin)
     &       .and. (BoundaryFibShellin(2,2) .eq. 2))then
          Y(1:3,1:3,-3:3,-3:3)=0.0
!! ================== BNDY 1,1
       elseif(((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   .and.
     &    ((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   )then
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
              Y(i2,j2,0:1,0)=Y(i2,j2,0:1,0)+Siga_d(i,j,1)*DX_ip_j(j2)
     &                      *(/-1.0/ds1(i,j),1.0/ds1(i,j)/)
     &                  *-(2.0*(GX(i2,i,j) - GX(i2,ihp1,j)))/ds1(i,j)**2  
              enddo;enddo       
              DY_i_jp(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)=Y(i2,j2,0,0:1)+Siga_d(i,j,4)*DY_i_jp(j2)
     &                      *(/-1.0/ds2(i,j),1.0/ds2(i,j)/)
     &                  *-(2.0*(GX(i2,i,j) - GX(i2,i,jhp1)))/ds2(i,j)**2  
              enddo;enddo  
!! ================== BNDY Nx,1
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   )then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
            Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &                *(2.0*(GX(i2,ihm1,j)-GX(i2,i,j)))/ds1(ihm1,j)**2  
              enddo;enddo       
              DY_i_jp(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)
     &                =Y(i2,j2,0,0:1)
     &                +Siga_d(i,j,4)*DY_i_jp(j2)
     &                 *(/-1.0/ds2(i,j),1.0/ds2(i,j)/)
     &                 *-(2.0*(GX(i2,i,j) - GX(i2,i,jhp1)))/ds2(i,j)**2  
              enddo;enddo  
!! ================== BNDY 1,Ny
       elseif(((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
              Y(i2,j2,0:1,0)=
     &                Y(i2,j2,0:1,0)
     &                +Siga_d(i,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(i,j),1.0/ds1(i,j)/)
     &                *-(2.0*(GX(i2,i,j) - GX(i2,ihp1,j)))/ds1(i,j)**2  
              enddo;enddo       
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
     &                *(2.0*(GX(i2,i,jhm1)-GX(i2,i,j)))/ds2(i,jhm1)**2  
              enddo;enddo  
!! ================== BNDY Nx,Ny
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
            Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &                *(2.0*(GX(i2,ihm1,j)-GX(i2,i,j)))/ds1(ihm1,j)**2  
              enddo;enddo       
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
     &                *(2.0*(GX(i2,i,jhm1)-GX(i2,i,j)))/ds2(i,jhm1)**2  
              enddo;enddo  
!! ================== BNDY 1,2
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .eq. 2)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
              Y(i2,j2,0:1,0)=
     &                Y(i2,j2,0:1,0)
     &                +Siga_d(i,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(i,j),1.0/ds1(i,j)/)
     &                *-(2.0*(GX(i2,i,j) - GX(i2,ihp1,j)))/ds1(i,j)**2  
              enddo;enddo       

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
              temp_ibm3=Siga_d(i,j,5)*
     &                -( (GX(i2,i,j) - GX(i2,i,jhp1))/(2.0*ds2(i,j)) 
     &                + (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))
     &                  /(2.0*ds2(ip1,j)) )/ds1(i,j)
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
              temp_ibm3=Siga_d(i,jhm1,5)*
     &                -((GX(i2,i,jhm1) - GX(i2,i,j))/(2.0*ds2(i,jhm1)) 
     &                + (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))
     &                /(2.0*ds2(ip1,jhm1)))/ds1(i,j)

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
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)=
     &                Y(i2,j2,0,0:1)
     &                +Siga_d(i,j,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,j ),1.0/ds2(i,j )/)
     &             *-(GX(i2,i,j) - GX(i2,i,jhp1))/(dh2(i,j)*ds2(i,j))  
              enddo;enddo 
!! ================== BNDY 1,Ny-1
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .eq. ns_ibmin-1)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
              Y(i2,j2,0:1,0)=
     &                Y(i2,j2,0:1,0)
     &                +Siga_d(i,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(i,j),1.0/ds1(i,j)/)
     &                *-(2.0*(GX(i2,i,j) - GX(i2,ihp1,j)))/ds1(i,j)**2  
              enddo;enddo       

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
              temp_ibm3=Siga_d(i,j,5)*
     &                -( (GX(i2,i,j) - GX(i2,i,jhp1))/(2.0*ds2(i,j)) 
     &                + (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))
     &                  /(2.0*ds2(ip1,j)) )/ds1(i,j)
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
              temp_ibm3=Siga_d(i,jhm1,5)*
     &                -((GX(i2,i,jhm1) - GX(i2,i,j))/(2.0*ds2(i,jhm1)) 
     &                + (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))
     &                /(2.0*ds2(ip1,jhm1)))/ds1(i,j)

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
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)=
     &                Y(i2,j2,0,0:1)
     &                +Siga_d(i,j,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,j ),1.0/ds2(i,j )/)
     &             *-(GX(i2,i,j) - GX(i2,i,jhp1))/(dh2(i,j)*ds2(i,j))  
              enddo;enddo 
!! ================== BNDY Nx,2
       elseif((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .eq. 2)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
              Y(i2,j2,-1:0,0)=
     &                Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &                *(2.0*(GX(i2,ihm1,j) - GX(i2,i,j)))/ds1(ihm1,j)**2  
              enddo;enddo       

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
              temp_ibm3=Siga_d(ihm1,j,5)*
     &                -( (GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))
     &                   /(2.0*ds2(ihm1,jhm1)) 
     &                + (GX(i2,i,jhm1) - GX(i2,i,j))
     &                  /(2.0*ds2(i,jhm1)) )/ds1(ihm1,jhm1)
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
!sig12 -1,-1

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
              temp_ibm3=Siga_d(ihm1,jhm1,5)*
     &                ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))
     &                /(2.0*ds2(im1,jhm1)) 
     &                + (GX(i2,i,jhm1) - GX(i2,i,j))
     &                /(2.0*ds2(i,jhm1)))
     &                /ds1(ihm1,j)

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
              temp_ibm3=Siga_d(ihm1,j,5)*
     &                ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))
     &                /(2.0*ds2(im1,j)) 
     &                + (GX(i2,i,j) - GX(i2,i,jhp1))
     &                /(2.0*ds2(i,j)))
     &                /ds1(ihm1,j)
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
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)=
     &                Y(i2,j2,0,0:1)
     &                +Siga_d(i,j,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,j ),1.0/ds2(i,j )/)
     &             *-(GX(i2,i,j) - GX(i2,i,jhp1))/(dh2(i,j)*ds2(i,j))  
              enddo;enddo 
       elseif((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .eq. ns_ibmin-1)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
              Y(i2,j2,-1:0,0)=
     &                Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &                *(2.0*(GX(i2,ihm1,j) - GX(i2,i,j)))/ds1(ihm1,j)**2  
              enddo;enddo       

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
              temp_ibm3=Siga_d(ihm1,j,5)*
     &                -( (GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))
     &                  /(2.0*ds2(ihm1,jhm1)) 
     &                + (GX(i2,i,jhm1) - GX(i2,i,j))
     &                  /(2.0*ds2(i,jhm1)) )/ds1(ihm1,jhm1)
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
!sig12 -1,-1

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
              temp_ibm3=Siga_d(ihm1,jhm1,5)*
     &                ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))
     &                /(2.0*ds2(im1,jhm1)) 
     &                + (GX(i2,i,jhm1) - GX(i2,i,j))
     &                /(2.0*ds2(i,jhm1)))
     &                /ds1(ihm1,j)

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
              temp_ibm3=Siga_d(ihm1,j,5)*
     &                ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))
     &                /(2.0*ds2(im1,j)) 
     &                + (GX(i2,i,j) - GX(i2,i,jhp1))
     &                /(2.0*ds2(i,j)))
     &                /ds1(ihm1,j)
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
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)=
     &                Y(i2,j2,0,0:1)
     &                +Siga_d(i,j,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,j ),1.0/ds2(i,j )/)
     &             *-(GX(i2,i,j) - GX(i2,i,jhp1))/(dh2(i,j)*ds2(i,j))  
              enddo;enddo 
!! ================== BNDY 2,1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .eq. 2)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
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
              temp_ibm3=Siga_d(ihm1,j,5)*
     &                -((GX(i2,ihm1,j) - GX(i2,i,j))
     &                /(2.0*ds1(ihm1,j)) 
     &                + (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))
     &                /(2.0*ds1(ihm1,jp1)))/ds2(i,j)
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
              temp_ibm3=Siga_d(i,j,5)*
     &                -((GX(i2,i,j) - GX(i2,ihp1,j))
     &                /(2.0*ds1(i,j)) 
     &                + (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))
     &                /(2.0*ds1(i,jp1)))/ds2(i,j)
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
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)=Y(i2,j2,0,0:1)+Siga_d(i,j,4)*DY_i_jp(j2)
     &                      *(/-1.0/ds2(i,j),1.0/ds2(i,j)/)
     &                  *-(2.0*(GX(i2,i,j) - GX(i2,i,jhp1)))/ds2(i,j)**2  
              enddo;enddo 
!! ================== BNDY Nx-1,1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .eq. nq_ibmin-1)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
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
              temp_ibm3=Siga_d(ihm1,j,5)*
     &                -((GX(i2,ihm1,j) - GX(i2,i,j))
     &                /(2.0*ds1(ihm1,j)) 
     &                + (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))
     &                /(2.0*ds1(ihm1,jp1)))/ds2(i,j)
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
              temp_ibm3=Siga_d(i,j,5)*
     &                -((GX(i2,i,j) - GX(i2,ihp1,j))
     &                /(2.0*ds1(i,j)) 
     &                + (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))
     &                /(2.0*ds1(i,jp1)))/ds2(i,j)
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
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)=Y(i2,j2,0,0:1)+Siga_d(i,j,4)*DY_i_jp(j2)
     &                      *(/-1.0/ds2(i,j),1.0/ds2(i,j)/)
     &                  *-(2.0*(GX(i2,i,j) - GX(i2,i,jhp1)))/ds2(i,j)**2  
              enddo;enddo 
!! ================== BNDY 2,Ny
       elseif((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (i .eq. 2)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
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
              temp_ibm3=Siga_d(ihm1,jhm1,5)*
     &                ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))
     &                /(2.0*ds1(ihm1,jm1)) 
     &                + (GX(i2,ihm1,j) - GX(i2,i,j))
     &                /(2.0*ds1(ihm1,j)))/ds2(i,jhm1)
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
              temp_ibm3=Siga_d(i,jhm1,5)*
     &                ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))
     &                /(2.0*ds1(i,jm1)) 
     &                + (GX(i2,i,j) - GX(i2,ihp1,j))
     &                /(2.0*ds1(i,j)))/ds2(i,jhm1)
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
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                      Y(i2,j2,0,-1:0)
     &                      +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                      *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
     &               *(2.0*(GX(i2,i,jhm1) - GX(i2,i,j)))/ds2(i,jhm1)**2  
              enddo;enddo 
!! ================== BNDY Nx-1,Ny
       elseif((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (i .eq. nq_ibmin-1)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
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
              temp_ibm3=Siga_d(ihm1,jhm1,5)*
     &                ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))
     &                /(2.0*ds1(ihm1,jm1)) 
     &                + (GX(i2,ihm1,j) - GX(i2,i,j))
     &                /(2.0*ds1(ihm1,j)))/ds2(i,jhm1)
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
              temp_ibm3=Siga_d(i,jhm1,5)*
     &                ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))
     &                /(2.0*ds1(i,jm1)) 
     &                + (GX(i2,i,j) - GX(i2,ihp1,j))
     &                /(2.0*ds1(i,j)))/ds2(i,jhm1)
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
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                      Y(i2,j2,0,-1:0)
     &                      +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                      *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
     &               *(2.0*(GX(i2,i,jhm1) - GX(i2,i,j)))/ds2(i,jhm1)**2  
              enddo;enddo 
!! ================== BNDY 1,j  2<=j<=Ny-1
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .ge. jst) .and. (j.le. jen)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,0:1,0)
     &                =Y(i2,j2,0:1,0)
     &                +Siga_d(i,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(i,j),1.0/ds1(i,j)/)
     &            *-(2.0*(GX(i2,i,j) - GX(i2,ihp1,j)))/ds1(i,j)**2   
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
     &                /(2.0*ds2(i,jhm1)) 
     &                + (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))
     &                /(2.0*ds2(ip1,jhm1)))/ds1(i,j)
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
     &               /(2.0*ds2(i,j)) 
     &               + (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))
     &               /(2.0*ds2(ip1,j)))/ds1(i,j)
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
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
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
!! ================== BNDY N,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin)   
     &        .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .ge. jst) .and. (j.le. jen)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &            *(GX(i2,ihm1,j) - GX(i2,i,j))/(dh1(i,j)*ds1(ihm1,j))   
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
     &                /(2.0*ds2(i,jhm1)) 
     &                + (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))
     &                /(2.0*ds2(ip1,jhm1)))/ds1(i,j)
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
     &               /(2.0*ds2(i,j)) 
     &               + (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))
     &               /(2.0*ds2(ip1,j)))/ds1(i,j)
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
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
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
!! ================== BNDY i,1 2<=i<=Nx-1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .ge. ist) .and. (i.le. ien)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
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
     &                -((GX(i2,ihm1,j) - GX(i2,i,j))
     &               /(2.0*ds1(ihm1,j)) 
     &               + (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))
     &               /(2.0*ds1(ihm1,jp1)))/ds2(i,j)
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
     &               -((GX(i2,i,j) - GX(i2,ihp1,j))
     &               /(2.0*ds1(i,j)) 
     &               + (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))
     &               /(2.0*ds1(i,jp1)))/ds2(i,j)
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
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,0:1)=
     &                Y(i2,j2,0,0:1)
     &                +Siga_d(i,j,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,j ),1.0/ds2(i,j )/)
     &             *-(GX(i2,i,j) - GX(i2,i,jhp1))/(dh2(i,j)*ds2(i,j))  
              enddo;enddo 
!! ================== BNDY i,N 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin)   
     &        .and. (BoundaryFibShellin(1,2) .eq. 1)  
     &    .and. (i .ge. ist) .and. (i.le. ien)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &            *(GX(i2,ihm1,j) - GX(i2,i,j))/(dh1(i,j)*ds1(ihm1,j))   
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
     &                /(2.0*ds2(i,jhm1)) 
     &                + (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))
     &                /(2.0*ds2(ip1,jhm1)))/ds1(i,j)
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
     &               /(2.0*ds2(i,j)) 
     &               + (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))
     &               /(2.0*ds2(ip1,j)))/ds1(i,j)
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
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
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
       else

!! ================== BNDY 2,2
!! ================== BNDY Nx-1,2
!! ================== BNDY 2,Ny-1
!! ================== BNDY Nx-1,Ny-1
!! ========================  row 2   1,j
!! ========================  row 2   N,j
!! ========================  row 2   i,1
!!  ========================  row 2   i,N
!! ===============  Middle point

              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
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
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
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
       endif

      if(MaterialTypeIBM .eq. 1) then
       if((j .le. 2).and. (BoundaryFibShellin(1,1) .eq. 0))then

       elseif((j .ge. ns_ibmin-1)
     &       .and. (BoundaryFibShellin(1,2) .eq. 0))then

       elseif((i .le. 2).and. (BoundaryFibShellin(2,1) .eq. 0))then

       elseif((i .ge. nq_ibmin-1)
     &       .and. (BoundaryFibShellin(2,2) .eq. 0))then


       elseif((j .le. 1).and. (BoundaryFibShellin(1,1) .eq. 2))then

       elseif((j .ge. ns_ibmin)
     &       .and. (BoundaryFibShellin(1,2) .eq. 2))then

       elseif((i .le. 1).and. (BoundaryFibShellin(2,1) .eq. 2))then

       elseif((i .ge. nq_ibmin)
     &       .and. (BoundaryFibShellin(2,2) .eq. 2))then

!! ================== BNDY 1,1
       elseif(((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   .and.
     &    ((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   )then
!! +++++ GAM11
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ dh1(ihp1,j) - ((GX(i2,ihp1,j) 
     & - GX(i2,ihp2,j))/ ds1(ihp1,j) - (GX(i2,ihp2,j) - 
     & GX(i2,ihp3,j))/ ds1(ihp2,j))/ dh1(ihp2,j))/ (ds1(i,j)* 
     & ds1(ihp1,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(4.0d0* ((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - 
     & (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))/ ds2(ihp1,j)))/ 
     & (dg1(i,j)* dg2(i,j)* ds1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(4.0d0* ((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - 
     & (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))/ ds1(i,jhp1)))/ 
     & (dg1(i,j)* dg2(i,j)* ds2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ dh2(i,jhp1) - ((GX(i2,i,jhp1) 
     & - GX(i2,i,jhp2))/ ds2(i,jhp1) - (GX(i2,i,jhp2) - 
     & GX(i2,i,jhp3))/ ds2(i,jhp2))/ dh2(i,jhp2))/ (ds2(i,j)* 
     & ds2(i,jhp1)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 

!! ================== BNDY Nx,1
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1)) 
     &   )then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,ihm3,j) - GX(i2,ihm2,j))/ ds1(ihm3,j) - 
     & (GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j))/ dh1(ihm2,j) - 
     & ((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ dh1(ihm1,j))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 


!! +++++ GAM12
!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (4.0d0* ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j)))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* ds1(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (4.0d0* ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1)))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* ds2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM22
!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ dh2(i,jhp1) - ((GX(i2,i,jhp1) 
     & - GX(i2,i,jhp2))/ ds2(i,jhp1) - (GX(i2,i,jhp2) - 
     & GX(i2,i,jhp3))/ ds2(i,jhp2))/ dh2(i,jhp2))/ (ds2(i,j)* 
     & ds2(i,jhp1)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 


!! ================== BNDY 1,Ny
       elseif(((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
!! +++++ GAM11
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ dh1(ihp1,j) - ((GX(i2,ihp1,j) 
     & - GX(i2,ihp2,j))/ ds1(ihp1,j) - (GX(i2,ihp2,j) - 
     & GX(i2,ihp3,j))/ ds1(ihp2,j))/ dh1(ihp2,j))/ (ds1(i,j)* 
     & ds1(ihp1,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (4.0d0* ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1)))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* ds1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (4.0d0* ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j)))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* ds2(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,i,jhm3) - GX(i2,i,jhm2))/ ds2(i,jhm3) - 
     & (GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2))/ dh2(i,jhm2) - 
     & ((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ dh2(i,jhm1))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!! ================== BNDY Nx,Ny
       elseif(((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1)) 
     &   .and.
     &    ((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1)) 
     &   )then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,ihm3,j) - GX(i2,ihm2,j))/ ds1(ihm3,j) - 
     & (GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j))/ dh1(ihm2,j) - 
     & ((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ dh1(ihm1,j))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 


!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(4.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) 
     & - (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1)))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* ds1(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(4.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) 
     & - (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j)))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* ds2(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,i,jhm3) - GX(i2,i,jhm2))/ ds2(i,jhm3) - 
     & (GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2))/ dh2(i,jhm2) - 
     & ((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ dh2(i,jhm1))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!! ================== BNDY 1,2
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .eq. 2)) then
!! +++++ GAM11
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ dh1(ihp1,j) - ((GX(i2,ihp1,j) 
     & - GX(i2,ihp2,j))/ ds1(ihp1,j) - (GX(i2,ihp2,j) - 
     & GX(i2,ihp3,j))/ ds1(ihp2,j))/ dh1(ihp2,j))/ (ds1(i,j)* 
     & ds1(ihp1,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1)))/ 
     & (dg1(i,jhm1)* dk2(i,j)* ds1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - 
     & (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))/ ds2(ihp1,j)))/ 
     & (dg1(i,j)* dk2(i,j)* ds1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j)))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* dk2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - 
     & (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))/ ds1(i,jhp1)))/ 
     & (dg1(i,j)* dg2(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 

!! ================== BNDY 1,Ny-1
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .eq. ns_ibmin-1)) then
!! +++++ GAM11
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ dh1(ihp1,j) - ((GX(i2,ihp1,j) 
     & - GX(i2,ihp2,j))/ ds1(ihp1,j) - (GX(i2,ihp2,j) - 
     & GX(i2,ihp3,j))/ ds1(ihp2,j))/ dh1(ihp2,j))/ (ds1(i,j)* 
     & ds1(ihp1,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1)))/ 
     & (dg1(i,jhm1)* dk2(i,j)* ds1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - 
     & (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))/ ds2(ihp1,j)))/ 
     & (dg1(i,j)* dk2(i,j)* ds1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j)))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* dk2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - 
     & (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))/ ds1(i,jhp1)))/ 
     & (dg1(i,j)* dg2(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY Nx,2
       elseif((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .eq. 2)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,ihm3,j) - GX(i2,ihm2,j))/ ds1(ihm3,j) - 
     & (GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j))/ dh1(ihm2,j) - 
     & ((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ dh1(ihm1,j))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 


!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) 
     & - (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1)))/ 
     & (dg1(ihm1,jhm1)* dk2(ihm1,j)* ds1(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j)))/ (dg1(ihm1,j)* 
     & dk2(ihm1,j)* ds1(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) 
     & - (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j)))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1)))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM22
!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY Nx,Ny-1
       elseif((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .eq. ns_ibmin-1)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,ihm3,j) - GX(i2,ihm2,j))/ ds1(ihm3,j) - 
     & (GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j))/ dh1(ihm2,j) - 
     & ((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ dh1(ihm1,j))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 


!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) 
     & - (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1)))/ 
     & (dg1(ihm1,jhm1)* dk2(ihm1,j)* ds1(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j)))/ (dg1(ihm1,j)* 
     & dk2(ihm1,j)* ds1(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) 
     & - (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j)))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1)))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY 2,1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .eq. 2)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j)))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - 
     & (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))/ ds2(ihp1,j)))/ 
     & (dg1(i,j)* dg2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1)))/ 
     & (dg2(ihm1,j)* dk1(i,j)* ds2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - 
     & (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))/ ds1(i,jhp1)))/ 
     & (dg2(i,j)* dk1(i,j)* ds2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ dh2(i,jhp1) - ((GX(i2,i,jhp1) 
     & - GX(i2,i,jhp2))/ ds2(i,jhp1) - (GX(i2,i,jhp2) - 
     & GX(i2,i,jhp3))/ ds2(i,jhp2))/ dh2(i,jhp2))/ (ds2(i,j)* 
     & ds2(i,jhp1)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY Nx-1,1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .eq. nq_ibmin-1)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo


!! +++++ GAM12
!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j)))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - 
     & (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))/ ds2(ihp1,j)))/ 
     & (dg1(i,j)* dg2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1)))/ 
     & (dg2(ihm1,j)* dk1(i,j)* ds2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - 
     & (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))/ ds1(i,jhp1)))/ 
     & (dg2(i,j)* dk1(i,j)* ds2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ dh2(i,jhp1) - ((GX(i2,i,jhp1) 
     & - GX(i2,i,jhp2))/ ds2(i,jhp1) - (GX(i2,i,jhp2) - 
     & GX(i2,i,jhp3))/ ds2(i,jhp2))/ dh2(i,jhp2))/ (ds2(i,j)* 
     & ds2(i,jhp1)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY 2,Ny
       elseif((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (i .eq. 2)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) 
     & - (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1)))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1)))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) 
     & - (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j)))/ 
     & (dg2(ihm1,jhm1)* dk1(i,jhm1)* ds2(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j)))/ (dg2(i,jhm1)* 
     & dk1(i,jhm1)* ds2(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,i,jhm3) - GX(i2,i,jhm2))/ ds2(i,jhm3) - 
     & (GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2))/ dh2(i,jhm2) - 
     & ((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ dh2(i,jhm1))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY Nx-1,Ny
       elseif((j.eq. ns_ibmin) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (i .eq. nq_ibmin-1)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo


!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) 
     & - (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1)))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1)))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;


              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) 
     & - (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j)))/ 
     & (dg2(ihm1,jhm1)* dk1(i,jhm1)* ds2(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j)))/ (dg2(i,jhm1)* 
     & dk1(i,jhm1)* ds2(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,i,jhm3) - GX(i2,i,jhm2))/ ds2(i,jhm3) - 
     & (GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2))/ dh2(i,jhm2) - 
     & ((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ dh2(i,jhm1))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY 2,2
       elseif( 
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
     &       )  then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dg1(ihm1,jhm1)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j))/ (dg1(ihm1,j)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1))/ 
     & (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp1,jhp1))/ ds2(ihp1,j))/ (dg1(i,j)* dh1(i,j)* 
     & dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1))/ 
     & (dg2(ihm1,j)* dh2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j))/ (dg2(i,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,ihp1,jhp1))/ ds1(i,jhp1))/ (dg2(i,j)* dh2(i,j)* 
     & dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY Nx-1,2
       elseif( 
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
     &       )  then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo


!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dg1(ihm1,jhm1)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j))/ (dg1(ihm1,j)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1))/ 
     & (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp1,jhp1))/ ds2(ihp1,j))/ (dg1(i,j)* dh1(i,j)* 
     & dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1))/ 
     & (dg2(ihm1,j)* dh2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j))/ (dg2(i,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,ihp1,jhp1))/ ds1(i,jhp1))/ (dg2(i,j)* dh2(i,j)* 
     & dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY 2,Ny-1
       elseif( 
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
     &       )  then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dg1(ihm1,jhm1)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j))/ (dg1(ihm1,j)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1))/ 
     & (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp1,jhp1))/ ds2(ihp1,j))/ (dg1(i,j)* dh1(i,j)* 
     & dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1))/ 
     & (dg2(ihm1,j)* dh2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j))/ (dg2(i,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,ihp1,jhp1))/ ds1(i,jhp1))/ (dg2(i,j)* dh2(i,j)* 
     & dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY Nx-1,Ny-1
       elseif( 
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
     &       )  then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo


!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dg1(ihm1,jhm1)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j))/ (dg1(ihm1,j)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1))/ 
     & (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp1,jhp1))/ ds2(ihp1,j))/ (dg1(i,j)* dh1(i,j)* 
     & dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1))/ 
     & (dg2(ihm1,j)* dh2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j))/ (dg2(i,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,ihp1,jhp1))/ ds1(i,jhp1))/ (dg2(i,j)* dh2(i,j)* 
     & dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY 1,j  2<=j<=Ny-1
       elseif((i.eq. 1) .and. (BoundaryFibShellin(2,1) .eq. 1) 
     &    .and. (j .ge. jst) .and. (j.le. jen)) then
!! +++++ GAM11
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ dh1(ihp1,j) - ((GX(i2,ihp1,j) 
     & - GX(i2,ihp2,j))/ ds1(ihp1,j) - (GX(i2,ihp2,j) - 
     & GX(i2,ihp3,j))/ ds1(ihp2,j))/ dh1(ihp2,j))/ (ds1(i,j)* 
     & ds1(ihp1,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1)))/ 
     & (dg1(i,jhm1)* dk2(i,j)* ds1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - 
     & (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))/ ds2(ihp1,j)))/ 
     & (dg1(i,j)* dk2(i,j)* ds1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j)))/ (dg1(i,jhm1)* 
     & dg2(i,jhm1)* dk2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - 
     & (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))/ ds1(i,jhp1)))/ 
     & (dg1(i,j)* dg2(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 

!! ================== BNDY N,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin)   
     &        .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .ge. jst) .and. (j.le. jen)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,ihm3,j) - GX(i2,ihm2,j))/ ds1(ihm3,j) - 
     & (GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j))/ dh1(ihm2,j) - 
     & ((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ dh1(ihm1,j))/ 
     & (ds1(ihm2,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 


!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) 
     & - (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1)))/ 
     & (dg1(ihm1,jhm1)* dk2(ihm1,j)* ds1(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j)))/ (dg1(ihm1,j)* 
     & dk2(ihm1,j)* ds1(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) 
     & - (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j)))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1)))/ 
     & (dg1(ihm1,j)* dg2(ihm1,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY i,1 2<=i<=Nx-1
       elseif((j.eq. 1) .and. (BoundaryFibShellin(1,1) .eq. 1) 
     &    .and. (i .ge. ist) .and. (i.le. ien)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j)))/ (dg1(ihm1,j)* 
     & dg2(ihm1,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - 
     & (GX(i2,ihp1,j) - GX(i2,ihp1,jhp1))/ ds2(ihp1,j)))/ 
     & (dg1(i,j)* dg2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1)))/ 
     & (dg2(ihm1,j)* dk1(i,j)* ds2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - 
     & (GX(i2,i,jhp1) - GX(i2,ihp1,jhp1))/ ds1(i,jhp1)))/ 
     & (dg2(i,j)* dk1(i,j)* ds2(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ dh2(i,jhp1) - ((GX(i2,i,jhp1) 
     & - GX(i2,i,jhp2))/ ds2(i,jhp1) - (GX(i2,i,jhp2) - 
     & GX(i2,i,jhp3))/ ds2(i,jhp2))/ dh2(i,jhp2))/ (ds2(i,j)* 
     & ds2(i,jhp1)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ================== BNDY i,N 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin)   
     &        .and. (BoundaryFibShellin(1,2) .eq. 1)  
     &    .and. (i .ge. ist) .and. (i.le. ien)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) 
     & - (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1)))/ 
     & (dg1(ihm1,jhm1)* dg2(ihm1,jhm1)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1)))/ 
     & (dg1(i,jhm1)* dg2(i,jhm1)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(2.0d0* ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) 
     & - (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j)))/ 
     & (dg2(ihm1,jhm1)* dk1(i,jhm1)* ds2(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (2.0d0* ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j)))/ (dg2(i,jhm1)* 
     & dk1(i,jhm1)* ds2(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -(((GX(i2,i,jhm3) - GX(i2,i,jhm2))/ ds2(i,jhm3) - 
     & (GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2))/ dh2(i,jhm2) - 
     & ((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ dh2(i,jhm1))/ 
     & (ds2(i,jhm2)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 
!! ========================  row 2   2,j 2<=j<=Ny-1
       elseif((i.eq. 2) .and. (
     &                          (BoundaryFibShellin(2,1) .eq. 1)
     &                           .or. 
     &                          (BoundaryFibShellin(2,1) .eq. 2)
     &                         )
     &    .and. (j .ge. jst) .and. (j.le. jen)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dg1(ihm1,jhm1)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j))/ (dg1(ihm1,j)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1))/ 
     & (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp1,jhp1))/ ds2(ihp1,j))/ (dg1(i,j)* dh1(i,j)* 
     & dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1))/ 
     & (dg2(ihm1,j)* dh2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j))/ (dg2(i,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,ihp1,jhp1))/ ds1(i,jhp1))/ (dg2(i,j)* dh2(i,j)* 
     & dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ========================  row 2   Nx-1,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin-1) .and. (
     &                          (BoundaryFibShellin(2,2) .eq. 1)
     &                           .or. 
     &                          (BoundaryFibShellin(2,2) .eq. 2)
     &                         )
     &    .and. (j .ge. jst) .and. (j.le. jen)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo


!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dg1(ihm1,jhm1)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j))/ (dg1(ihm1,j)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1))/ 
     & (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp1,jhp1))/ ds2(ihp1,j))/ (dg1(i,j)* dh1(i,j)* 
     & dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1))/ 
     & (dg2(ihm1,j)* dh2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j))/ (dg2(i,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,ihp1,jhp1))/ ds1(i,jhp1))/ (dg2(i,j)* dh2(i,j)* 
     & dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ========================  row 2   i,2 2<=i<=Nx-1
       elseif((j.eq. 2) .and. (
     &                          (BoundaryFibShellin(1,1) .eq. 1)
     &                           .or. 
     &                          (BoundaryFibShellin(1,1) .eq. 2)
     &                         )
     &    .and. (i .ge. ist) .and. (i.le. ien)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dg1(ihm1,jhm1)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j))/ (dg1(ihm1,j)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1))/ 
     & (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp1,jhp1))/ ds2(ihp1,j))/ (dg1(i,j)* dh1(i,j)* 
     & dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1))/ 
     & (dg2(ihm1,j)* dh2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j))/ (dg2(i,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,ihp1,jhp1))/ ds1(i,jhp1))/ (dg2(i,j)* dh2(i,j)* 
     & dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j)) 
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3  
              enddo;enddo 
!! ========================  row 2   i,Ny 2<=i<=Nx-1
       elseif((j.eq. ns_ibmin-1) .and. (
     &                          (BoundaryFibShellin(1,2) .eq. 1)
     &                           .or. 
     &                          (BoundaryFibShellin(1,2) .eq. 2)
     &                         )
     &    .and. (i .ge. ist) .and. (i.le. ien)) then
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dg1(ihm1,jhm1)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j))/ (dg1(ihm1,j)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1))/ 
     & (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp1,jhp1))/ ds2(ihp1,j))/ (dg1(i,j)* dh1(i,j)* 
     & dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1))/ 
     & (dg2(ihm1,j)* dh2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j))/ (dg2(i,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,ihp1,jhp1))/ ds1(i,jhp1))/ (dg2(i,j)* dh2(i,j)* 
     & dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 
!other points
       else
!! +++++ GAM11
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihm1,j)-GX(1:3,ihm2,j))/ds1(ihm2,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihm1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm2,j) - GX(i2,ihm1,j))/ ds1(ihm2,j) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dh1(ihm1,j)* 
     & dh1(i,j)* ds1(ihm1,j)) 
                 Y(i2,j2,-2:0,0)
     &                =Y(i2,j2,-2:0,0)
     &                +gama_d(ihm1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm2,j)
     &                 ,-1.0/ds1(ihm2,j)-1.0/ds1(ihm1,j)
     &                  ,1.0/ds1(ihm1,j)        
     &                 /)/dh1(ihm1,j)
     &                * temp_ibm3  
              enddo;enddo 
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(ihm1,j)) + 
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - (GX(i2,i,j) - 
     & GX(i2,ihp1,j))/ ds1(i,j))/ (dh1(i,j)* ds1(i,j)))/ dh1(i,j) 
                 Y(i2,j2,-1:1,0)
     &                =Y(i2,j2,-1:1,0)
     &                +gama_d(i,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(ihm1,j)
     &                 ,-1.0/ds1(ihm1,j)-1.0/ds1(i,j)
     &                  ,1.0/ds1(i,j)        
     &                 /)/dh1(i,j)
     &                * temp_ibm3  
              enddo;enddo
!gama11(1,0)
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp2,j)-GX(1:3,ihp1,j))/ds1(ihp1,j)   !x+1/2,y
              DX_im_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dh1(ihp1,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp2,j))/ ds1(ihp1,j))/ (dh1(i,j)* dh1(ihp1,j)* 
     & ds1(i,j)) 
                 Y(i2,j2,0:2,0)
     &                =Y(i2,j2,0:2,0)
     &                +gama_d(ihp1,j,1)*DXX_i_j(j2)
     &                *(/1.0/ds1(i,j)
     &                 ,-1.0/ds1(i,j)-1.0/ds1(ihp1,j)
     &                  ,1.0/ds1(ihp1,j)        
     &                 /)/dh1(ihp1,j)
     &                * temp_ibm3  
              enddo;enddo



!! +++++ GAM12
!gama12(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))/ ds2(ihm1,jhm1) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dg1(ihm1,jhm1)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)= 0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))/ ds2(ihm1,j) - 
     & (GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j))/ (dg1(ihm1,j)* 
     & dh1(i,j)* dk2(ihm1,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - 
     & (GX(i2,ihp1,jhm1) - GX(i2,ihp1,j))/ ds2(ihp1,jhm1))/ 
     & (dg1(i,jhm1)* dh1(i,j)* dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama12( 0,0)
             i1s(1)= 0;i1s(2)=1
             i2s(1)= 0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;
              DY_ip_j(1:3)=
     &        (GX(1:3,iu,ju)-GX(1:3,iu,jv))
     &                     /ds2(iu,jv) 
              DY_im_j(1:3)=
     &         (GX(1:3,iv,ju)-GX(1:3,iv,jv))
     &                     /ds2(iv,jv) 

              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dg1(iv,jv)

              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,ihp1,j) - 
     & GX(i2,ihp1,jhp1))/ ds2(ihp1,j))/ (dg1(i,j)* dh1(i,j)* 
     & dk2(i,j)) 
                 Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(1))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iv,jv)
     &                   ,1.0/ds2(iv,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                =Y(i2,j2,i1s(1):i1s(2),i2s(2))
     &                +gama_d(iv,jv,3)*DXY_i_j(j2)
     &                *(/-1.0/ds2(iu,jv)
     &                   ,1.0/ds2(iu,jv)/)
     &                 /dg1(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM21
!gama21(-1,-1)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=-1;i2s(2)=0
             iu=i;iv=ihm1;iw=ihm2;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))/ ds1(ihm1,jhm1) - 
     & (GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j))/ (dg2(ihm1,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(-1,0)
             i1s(1)=-1;i1s(2)=0
             i2s(1)=0;i2s(2)=1
             iu=i;iv=ihm1;iw=ihm2;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,ihm1,j) - GX(i2,i,j))/ ds1(ihm1,j) - 
     & (GX(i2,ihm1,jhp1) - GX(i2,i,jhp1))/ ds1(ihm1,jhp1))/ 
     & (dg2(ihm1,j)* dh2(i,j)* dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo

!gama21(0,-1)
             i1s(1)=0;i1s(2)=1
             i2s(1)=-1;i2s(2)=0
             iu=ihp1;iv=i;iw=ihm1;
             ju=j;jv=jhm1;jw=jhm2;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))/ ds1(i,jhm1) - 
     & (GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j))/ (dg2(i,jhm1)* 
     & dh2(i,j)* dk1(i,jhm1)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 

!gama21(0,0)
             i1s(1)=0;i1s(2)=1
             i2s(1)=0;i2s(2)=1
             iu=ihp1;iv=i;iw=ihm1;
             ju=jhp1;jv=j;jw=jhm1;

              DX_i_jp(1:3)=
     &         (GX(1:3,iu,ju)-GX(1:3,iv,ju))
     &                     /ds1(i,ju) 
              DX_i_jm(1:3)=
     &         (GX(1:3,iu,jv)-GX(1:3,iv,jv))
     &                     /ds1(iv,jv) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dg2(iv,jv)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,ihp1,j))/ ds1(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,ihp1,jhp1))/ ds1(i,jhp1))/ (dg2(i,j)* dh2(i,j)* 
     & dk1(i,j)) 
                 Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(1),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,jv)
     &                   ,1.0/ds1(iv,jv)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3
                 Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                =Y(i2,j2,i1s(2),i2s(1):i2s(2))
     &                +gama_d(iv,jv,3)*DYX_i_j(j2)
     &                *(/-1.0/ds1(iv,ju)
     &                   ,1.0/ds1(iv,ju)/)
     &                 /dg2(iv,jv)
     &                * temp_ibm3  
              enddo;enddo 



!! +++++ GAM22
!gama22(0,-1)
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              DY_i_jm(1:3)=
     &       (GX(1:3,i,jhm1)-GX(1:3,i,jhm2))/ds2(i,jhm2) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhm1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,jhm2) - GX(i2,i,jhm1))/ ds2(i,jhm2) - 
     & (GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1))/ (dh2(i,jhm1)* 
     & dh2(i,j)* ds2(i,jhm1)) 
                 Y(i2,j2,0,-2:0)
     &                =Y(i2,j2,0,-2:0)
     &                +gama_d(i,jhm1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm2)
     &                 ,-1.0/ds2(i,jhm2)-1.0/ds2(i,jhm1)
     &                  ,1.0/ds2(i,jhm1)        
     &                 /)/dh2(i,jhm1)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,0)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,j)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & (((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,jhm1)) + 
     & ((GX(i2,i,jhm1) - GX(i2,i,j))/ ds2(i,jhm1) - (GX(i2,i,j) - 
     & GX(i2,i,jhp1))/ ds2(i,j))/ (dh2(i,j)* ds2(i,j)))/ dh2(i,j) 
                 Y(i2,j2,0,-1:1)
     &                =Y(i2,j2,0,-1:1)
     &                +gama_d(i,j,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,jhm1)
     &                 ,-1.0/ds2(i,jhm1)-1.0/ds2(i,j)
     &                  ,1.0/ds2(i,j)        
     &                 /)/dh2(i,j)
     &                * temp_ibm3  
              enddo;enddo 

!gama22(0,1)
              DY_i_jp(1:3)=
     &      (GX(1:3,i,jhp2)-GX(1:3,i,jhp1))/ds2(i,jhp1) !x,y+1/2
              DY_i_jm(1:3)=
     &      (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /dh2(i,jhp1)
              do i2=1,3;do j2=1,3
                 temp_ibm3=
     & -((GX(i2,i,j) - GX(i2,i,jhp1))/ ds2(i,j) - (GX(i2,i,jhp1) - 
     & GX(i2,i,jhp2))/ ds2(i,jhp1))/ (dh2(i,j)* dh2(i,jhp1)* 
     & ds2(i,j))
                 Y(i2,j2,0,0:2)
     &                =Y(i2,j2,0,0:2)
     &                +gama_d(i,jhp1,2)*DYY_i_j(j2)
     &                *(/1.0/ds2(i,j)
     &                 ,-1.0/ds2(i,j)-1.0/ds2(i,jhp1)
     &                  ,1.0/ds2(i,jhp1)        
     &                 /)/dh2(i,jhp1)
     &                * temp_ibm3
              enddo;enddo 
       endif
      endif

      
      do i3=1,3
            Y(i3,i3,-3:3,-3:3)=-alpha_coefK(1)*
     &           Y(i3,i3,-3:3,-3:3)+X(-3:3,-3:3)
      do j3=1,3
            Y(i3,j3,0,0)=Y(i3,j3,0,0)+(1.0-X00)*
     &           aamat_intmp(i3,j3)
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
       do i2=-3,3
        do j2=-3,3
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
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      Subroutine TAGMAPCross(TAG0ij,MAPTAG0ij)
      implicit none
       integer TAG0ij,MAPTAG0ij(3)
       logical, save :: FirstCallTAGMAP = .TRUE.
       integer,save :: TAGMAP(100,3)
      
      if(FirstCallTAGMAP)then
      TAGMAP(1,1:3) = (/-2,-2,1/)
      TAGMAP(2,1:3) = (/-1,-2,1/)
      TAGMAP(3,1:3) = (/0,-2,1/)
      TAGMAP(4,1:3) = (/1,-2,1/)
      TAGMAP(5,1:3) = (/2,-2,1/)
      TAGMAP(6,1:3) = (/-2,-1,1/)
      TAGMAP(7,1:3) = (/-1,-1,1/)
      TAGMAP(8,1:3) = (/0,-1,1/)
      TAGMAP(9,1:3) = (/1,-1,1/)
      TAGMAP(10,1:3) = (/2,-1,1/)
      TAGMAP(11,1:3) = (/-2,0,1/)
      TAGMAP(12,1:3) = (/-1,0,1/)
      TAGMAP(13,1:3) = (/0,0,1/)
      TAGMAP(14,1:3) = (/1,0,1/)
      TAGMAP(15,1:3) = (/2,0,1/)
      TAGMAP(16,1:3) = (/-2,1,1/)
      TAGMAP(17,1:3) = (/-1,1,1/)
      TAGMAP(18,1:3) = (/0,1,1/)
      TAGMAP(19,1:3) = (/1,1,1/)
      TAGMAP(20,1:3) = (/2,1,1/)
      TAGMAP(21,1:3) = (/-2,2,1/)
      TAGMAP(22,1:3) = (/-1,2,1/)
      TAGMAP(23,1:3) = (/0,2,1/)
      TAGMAP(24,1:3) = (/1,2,1/)
      TAGMAP(25,1:3) = (/2,2,1/)
      TAGMAP(26,1:3) = (/-2,-2,2/)
      TAGMAP(27,1:3) = (/-1,-2,2/)
      TAGMAP(28,1:3) = (/0,-2,2/)
      TAGMAP(29,1:3) = (/1,-2,2/)
      TAGMAP(30,1:3) = (/2,-2,2/)
      TAGMAP(31,1:3) = (/-2,-1,2/)
      TAGMAP(32,1:3) = (/-1,-1,2/)
      TAGMAP(33,1:3) = (/0,-1,2/)
      TAGMAP(34,1:3) = (/1,-1,2/)
      TAGMAP(35,1:3) = (/2,-1,2/)
      TAGMAP(36,1:3) = (/-2,0,2/)
      TAGMAP(37,1:3) = (/-1,0,2/)
      TAGMAP(38,1:3) = (/0,0,2/)
      TAGMAP(39,1:3) = (/1,0,2/)
      TAGMAP(40,1:3) = (/2,0,2/)
      TAGMAP(41,1:3) = (/-2,1,2/)
      TAGMAP(42,1:3) = (/-1,1,2/)
      TAGMAP(43,1:3) = (/0,1,2/)
      TAGMAP(44,1:3) = (/1,1,2/)
      TAGMAP(45,1:3) = (/2,1,2/)
      TAGMAP(46,1:3) = (/-2,2,2/)
      TAGMAP(47,1:3) = (/-1,2,2/)
      TAGMAP(48,1:3) = (/0,2,2/)
      TAGMAP(49,1:3) = (/1,2,2/)
      TAGMAP(50,1:3) = (/2,2,2/)
      TAGMAP(51,1:3) = (/-2,-2,3/)
      TAGMAP(52,1:3) = (/-1,-2,3/)
      TAGMAP(53,1:3) = (/0,-2,3/)
      TAGMAP(54,1:3) = (/1,-2,3/) 
      TAGMAP(55,1:3) = (/2,-2,3/)
      TAGMAP(56,1:3) = (/-2,-1,3/)
      TAGMAP(57,1:3) = (/-1,-1,3/)
      TAGMAP(58,1:3) = (/0,-1,3/)
      TAGMAP(59,1:3) = (/1,-1,3/)
      TAGMAP(60,1:3) = (/2,-1,3/)
      TAGMAP(61,1:3) = (/-2,0,3/)
      TAGMAP(62,1:3) = (/-1,0,3/)
      TAGMAP(63,1:3) = (/0,0,3/)
      TAGMAP(64,1:3) = (/1,0,3/)
      TAGMAP(65,1:3) = (/2,0,3/)
      TAGMAP(66,1:3) = (/-2,1,3/)
      TAGMAP(67,1:3) = (/-1,1,3/)
      TAGMAP(68,1:3) = (/0,1,3/)
      TAGMAP(69,1:3) = (/1,1,3/)
      TAGMAP(70,1:3) = (/2,1,3/)
      TAGMAP(71,1:3) = (/-2,2,3/)
      TAGMAP(72,1:3) = (/-1,2,3/)
      TAGMAP(73,1:3) = (/0,2,3/)
      TAGMAP(74,1:3) = (/1,2,3/)
      TAGMAP(75,1:3) = (/2,2,3/)
      TAGMAP(76,1:3) = (/-2,-2,4/)
      TAGMAP(77,1:3) = (/-1,-2,4/)
      TAGMAP(78,1:3) = (/0,-2,4/)
      TAGMAP(79,1:3) = (/1,-2,4/)
      TAGMAP(80,1:3) = (/2,-2,4/)
      TAGMAP(81,1:3) = (/-2,-1,4/)
      TAGMAP(82,1:3) = (/-1,-1,4/)
      TAGMAP(83,1:3) = (/0,-1,4/)
      TAGMAP(84,1:3) = (/1,-1,4/)
      TAGMAP(85,1:3) = (/2,-1,4/)
      TAGMAP(86,1:3) = (/-2,0,4/)
      TAGMAP(87,1:3) = (/-1,0,4/)
      TAGMAP(88,1:3) = (/0,0,4/)
      TAGMAP(89,1:3) = (/1,0,4/)
      TAGMAP(90,1:3) = (/2,0,4/)
      TAGMAP(91,1:3) = (/-2,1,4/)
      TAGMAP(92,1:3) = (/-1,1,4/)
      TAGMAP(93,1:3) = (/0,1,4/)
      TAGMAP(94,1:3) = (/1,1,4/)
      TAGMAP(95,1:3) = (/2,1,4/)
      TAGMAP(96,1:3) = (/-2,2,4/)
      TAGMAP(97,1:3) = (/-1,2,4/)
      TAGMAP(98,1:3) = (/0,2,4/)
      TAGMAP(99,1:3) = (/1,2,4/)
      TAGMAP(100,1:3) = (/2,2,4/)
      FirstCallTAGMAP = .FALSE.
      endif
      MAPTAG0ij(1:3)=TAGMAP(TAG0ij,1:3)
      return
      end subroutine
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      Subroutine TAGMAPddKapa(TAG0ij,MAPTAG0ij)
      implicit none
       integer TAG0ij,MAPTAG0ij(3)
       logical, save :: FirstCallTAGMAPddKapa = .TRUE.
       integer,save :: TAGMAPddKap(81,4)
       integer j1,ixua,j2
      
      if(FirstCallTAGMAPddKapa)then

      do j1=1,81
         ixua=int((j1-1)/(9*3))
         TAGMAPddKap(j1,1)=ixua+1

         j2=j1-ixua*(9*3)
         ixua=int((j2-1)/(9))
         TAGMAPddKap(j1,2)=ixua+1

         j2=j2-ixua*(9)
         ixua=int((j2-1)/(3))
         TAGMAPddKap(j1,3)=ixua-1
      
         ixua=j2-ixua*3
         TAGMAPddKap(j1,4)=ixua-1
      enddo
      FirstCallTAGMAPddKapa = .FALSE.
      endif
      MAPTAG0ij(1:4)=TAGMAPddKap(TAG0ij,1:4)
      return
      end subroutine
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|      

