c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine matfibrous_Huang_GALPHA(
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


      real*8 DX_i_j(3),DX_ip_j(3),DX_im_j(3)
     &      ,DY_i_j(3),DY_ip_j(3),DY_im_j(3)
     &      ,DX_i_jp(3),DX_ip_jp(3),DX_im_jp(3)
     &      ,DY_i_jp(3),DY_ip_jp(3),DY_im_jp(3)
     &      ,DX_i_jm(3),DX_ip_jm(3),DX_im_jm(3)
     &      ,DY_i_jm(3),DY_ip_jm(3),DY_im_jm(3)
     &      ,DXX_i_j(3),DYY_i_j(3),DXY_i_j(3),DYX_i_j(3)
     
      integer i,j,im1,im2,im3,jm1,jm2,jm3,ihm1,ihm2,ihm3,jhm1,jhm2,jhm3
      integer ip1,ip2,ip3,jp1,jp2,jp3,ihp1,ihp2,ihp3,jhp1,jhp2,jhp3,k
      real*8 X(-3:3,-3:3),coef,X00,Y(3,3,-3:3,-3:3),temp_ibm3
      integer iAdd(-3:3,-3:3,1:3)
     &  ,iiAdd(-3:3),jiAdd(-3:3),i2,j2,ntemi,noj
     &  ,ist,ien,jst,jen,iflagsimulation,nsnq,i3,j3,nojtmp

      integer i1s(2),i2s(2),iu,iv,iw,ju,jv,jw
!! ================== Inside i,j  2<i<nx-1; 2<j<ny-1
      ForceOut=0.0
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

         do i2=-3,3
           do j2=-3,3
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
        write(*,*) '1,1' 
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
        write(*,*) 'N,1' 
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
        write(*,*) '1,N' 
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
        write(*,*) 'N,N'
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
            X(i2,j2)=-X(i2,j2)*alpha_coefK
          enddo
        enddo
       X(0,0)=X(0,0)+X00*scal_K

       Y=0.0       !KSHTEST 20170524-8:17S
       !goto 1314  !KSHTEST 20170524-8:17S

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
       print*,'Nx==1'
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
       print*,'1==Ny'
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
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &            *(2.0*(GX(i2,ihm1,j) - GX(i2,i,j)))/(ds1(ihm1,j)**2)   
              enddo;enddo 
!sig12 -1,0

              DY_ip_j(1:3)=     !i,j+0.5
     &        (GX(1:3,i,jhp1)-GX(1:3,i,j))
     &                     /ds2(i,j) 
              DY_im_j(1:3)=
     &         (GX(1:3,ihm1,jhp1)-GX(1:3,ihm1,j))/ds2(ihm1,j)  !i-1,j+0.5
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))    !i-0.5,j+0.5

              DX_i_jp(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,ihm1,jhp1))
     &                     /ds1(ihm1,jhp1)                   !i-0.5,j+1
              DX_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j)) 
     &                     /ds1(ihm1,j)                      !i-0.5,j
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))    !i-0.5,j+0.5
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(ihm1,j,5)*(
     &                ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))
     &                /(2.0*ds2(ihm1,j)) 
     &                + (GX(i2,i,j) - GX(i2,i,jhp1))
     &                /(2.0*ds2(i,j)))/ds1(ihm1,j)
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
!sig12 -1,-1

              DY_ip_j(1:3)=
     &        (GX(1:3,ihm1,j)-GX(1:3,ihm1,jhm1))
     &                     /ds2(ihm1,jhm1)    !i-1,j-0.5
              DY_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1)   !i,j-0.5
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))  !i-0.5,j-0.5

              DX_i_jp(1:3)=
     &         (GX(1:3,i,jhm1)-GX(1:3,ihm1,jhm1))
     &                     /ds1(ihm1,jhm1)   !i-0.5,j-1
              DX_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))
     &                     /ds1(ihm1,j)      !i-0.5,j
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))  !i-0.5,j-0.5
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(ihm1,jhm1,5)*(
     &                ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))
     &                /(2.0*ds2(ihm1,jhm1)) 
     &                + (GX(i2,i,jhm1) - GX(i2,i,j))
     &                /(2.0*ds2(i,jhm1)))/ds1(ihm1,j)
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
!! ================== BNDY Nx,Ny-1
       elseif((i.eq. nq_ibmin) .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .eq. ns_ibmin-1)) then
              DX_ip_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))/ds1(ihm1,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,-1:0,0)
     &                =Y(i2,j2,-1:0,0)
     &                +Siga_d(ihm1,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(ihm1,j),1.0/ds1(ihm1,j)/)
     &            *(2.0*(GX(i2,ihm1,j) - GX(i2,i,j)))/(ds1(ihm1,j)**2)   
              enddo;enddo 
!sig12 -1,0

              DY_ip_j(1:3)=     !i,j+0.5
     &        (GX(1:3,i,jhp1)-GX(1:3,i,j))
     &                     /ds2(i,j) 
              DY_im_j(1:3)=
     &         (GX(1:3,ihm1,jhp1)-GX(1:3,ihm1,j))/ds2(ihm1,j)  !i-1,j+0.5
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))    !i-0.5,j+0.5

              DX_i_jp(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,ihm1,jhp1))
     &                     /ds1(ihm1,jhp1)                   !i-0.5,j+1
              DX_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j)) 
     &                     /ds1(ihm1,j)                      !i-0.5,j
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))    !i-0.5,j+0.5
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(ihm1,j,5)*(
     &                ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))
     &                /(2.0*ds2(ihm1,j)) 
     &                + (GX(i2,i,j) - GX(i2,i,jhp1))
     &                /(2.0*ds2(i,j)))/ds1(ihm1,j)
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
!sig12 -1,-1

              DY_ip_j(1:3)=
     &        (GX(1:3,ihm1,j)-GX(1:3,ihm1,jhm1))
     &                     /ds2(ihm1,jhm1)    !i-1,j-0.5
              DY_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1)   !i,j-0.5
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))  !i-0.5,j-0.5

              DX_i_jp(1:3)=
     &         (GX(1:3,i,jhm1)-GX(1:3,ihm1,jhm1))
     &                     /ds1(ihm1,jhm1)   !i-0.5,j-1
              DX_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))
     &                     /ds1(ihm1,j)      !i-0.5,j
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))  !i-0.5,j-0.5
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(ihm1,jhm1,5)*(
     &                ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))
     &                /(2.0*ds2(ihm1,jhm1)) 
     &                + (GX(i2,i,jhm1) - GX(i2,i,j))
     &                /(2.0*ds2(i,jhm1)))/ds1(ihm1,j)
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
     &            *-(2.0*(GX(i2,i,j) - GX(i2,ihp1,j)))/(ds1(i,j)**2)   
              enddo;enddo
!sig12 0,-1

              DY_ip_j(1:3)=
     &        (GX(1:3,ihp1,j)-GX(1:3,ihp1,jhm1))
     &                     /ds2(ihp1,jhm1)        !i+1,j-0.5
              DY_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1)    !i,j-0.5
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))  !i+0.5,j-0.5

              DX_i_jp(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))   !i+0.5,j
     &                     /ds1(i,j) 
              DX_i_jm(1:3)=
     &         (GX(1:3,ihp1,jhm1)-GX(1:3,i,jhm1))  
     &                     /ds1(i,jhm1)      !i+0.5,j-1
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))    !i+0.5,j-0.5
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
     &            *(2.0*(GX(i2,ihm1,j) - GX(i2,i,j)))/(ds1(ihm1,j)**2)   
              enddo;enddo 
!sig12 -1,0

              DY_ip_j(1:3)=     !i,j+0.5
     &        (GX(1:3,i,jhp1)-GX(1:3,i,j))
     &                     /ds2(i,j) 
              DY_im_j(1:3)=
     &         (GX(1:3,ihm1,jhp1)-GX(1:3,ihm1,j))/ds2(ihm1,j)  !i-1,j+0.5
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))    !i-0.5,j+0.5

              DX_i_jp(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,ihm1,jhp1))
     &                     /ds1(ihm1,jhp1)                   !i-0.5,j+1
              DX_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j)) 
     &                     /ds1(ihm1,j)                      !i-0.5,j
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))    !i-0.5,j+0.5
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(ihm1,j,5)*(
     &                ((GX(i2,ihm1,j) - GX(i2,ihm1,jhp1))
     &                /(2.0*ds2(ihm1,j)) 
     &                + (GX(i2,i,j) - GX(i2,i,jhp1))
     &                /(2.0*ds2(i,j)))/ds1(ihm1,j)
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
!sig12 -1,-1

              DY_ip_j(1:3)=
     &        (GX(1:3,ihm1,j)-GX(1:3,ihm1,jhm1))
     &                     /ds2(ihm1,jhm1)    !i-1,j-0.5
              DY_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1)   !i,j-0.5
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))  !i-0.5,j-0.5

              DX_i_jp(1:3)=
     &         (GX(1:3,i,jhm1)-GX(1:3,ihm1,jhm1))
     &                     /ds1(ihm1,jhm1)   !i-0.5,j-1
              DX_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))
     &                     /ds1(ihm1,j)      !i-0.5,j
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))  !i-0.5,j-0.5
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(ihm1,jhm1,5)*(
     &                ((GX(i2,ihm1,jhm1) - GX(i2,ihm1,j))
     &                /(2.0*ds2(ihm1,jhm1)) 
     &                + (GX(i2,i,jhm1) - GX(i2,i,j))
     &                /(2.0*ds2(i,jhm1)))/ds1(ihm1,j)
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

              DY_ip_j(1:3)=     !i,j+0.5
     &        (GX(1:3,i,jhp1)-GX(1:3,i,j))
     &                     /ds2(i,j) 
              DY_im_j(1:3)=
     &         (GX(1:3,ihm1,jhp1)-GX(1:3,ihm1,j))/ds2(ihm1,j)  !i-1,j+0.5
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))    !i-0.5,j+0.5

              DX_i_jp(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,ihm1,jhp1))
     &                     /ds1(ihm1,jhp1)                   !i-0.5,j+1
              DX_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j)) 
     &                     /ds1(ihm1,j)                      !i-0.5,j
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))    !i-0.5,j+0.5
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
     &                     /ds2(ihp1,j)   !i+1,j+0.5
              DY_im_j(1:3)=
     &         (GX(1:3,i,jhp1)-GX(1:3,i,j))/ds2(i,j)   !i,j+0.5 
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))

              DX_i_jp(1:3)=
     &         (GX(1:3,ihp1,jhp1)-GX(1:3,i,jhp1))
     &                     /ds1(i,jhp1)    !i+0.5,j+1
              DX_i_jm(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))
     &                     /ds1(i,j)    !i+0.5,j
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
     &             *-(2.0*(GX(i2,i,j)-GX(i2,i,jhp1)))/(ds2(i,j)**2)  
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
              DX_ip_j(1:3)=
     &         (GX(1:3,ihp1,j)-GX(1:3,i,j))/ds1(i,j)   !x+1/2,y
              do i2=1,3;do j2=1,3
                 Y(i2,j2,0:1,0)
     &                =Y(i2,j2,0:1,0)
     &                +Siga_d(i,j,1)*DX_ip_j(j2)
     &                *(/-1.0/ds1(i,j),1.0/ds1(i,j)/)
     &            *-(GX(i2,i,j) - GX(i2,ihp1,j))/(dh1(i,j)*ds1(i,j))   
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
     &                ((GX(i2,i,jhm1) - GX(i2,ihp1,jhm1))
     &                /(2.0*ds1(i,jhm1)) 
     &                + (GX(i2,i,j) - GX(i2,ihp1,j))
     &                /(2.0*ds1(i,j)))/ds2(i,jhm1)
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
!sig12 -1,-1

              DY_ip_j(1:3)=
     &        (GX(1:3,ihm1,j)-GX(1:3,ihm1,jhm1))
     &                     /ds2(ihm1,jhm1)    !i-1,j-0.5
              DY_im_j(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1)   !i,j-0.5
              DY_i_j(1:3)=0.5*(DY_ip_j(1:3)+DY_im_j(1:3))  !i-0.5,j-0.5

              DX_i_jp(1:3)=
     &         (GX(1:3,i,jhm1)-GX(1:3,ihm1,jhm1))
     &                     /ds1(ihm1,jhm1)   !i-0.5,j-1
              DX_i_jm(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,ihm1,j))
     &                     /ds1(ihm1,j)      !i-0.5,j
              DX_i_j(1:3)=0.5*(DX_i_jp(1:3)+DX_i_jm(1:3))  !i-0.5,j-0.5
              do i2=1,3;do j2=1,3
              temp_ibm3=Siga_d(ihm1,jhm1,5)*(
     &               ((GX(i2,ihm1,jhm1) - GX(i2,i,jhm1))
     &               /(2.0*ds1(ihm1,jhm1)) 
     &               + (GX(i2,ihm1,j) - GX(i2,i,j))
     &               /(2.0*ds1(ihm1,j)))/ds2(i,jhm1)
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
              DY_i_jp(1:3)=
     &         (GX(1:3,i,j)-GX(1:3,i,jhm1))/ds2(i,jhm1) !x,y+1/2
              do i2=1,3;do j2=1,3
              Y(i2,j2,0,-1:0)=
     &                Y(i2,j2,0,-1:0)
     &                +Siga_d(i,jhm1,4)*DY_i_jp(j2)
     &                *(/-1.0/ds2(i,jhm1),1.0/ds2(i,jhm1)/)
     &             *(2.0*(GX(i2,i,jhm1)-GX(i2,i,j)))/(ds2(i,jhm1)**2)  
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

 !1314 continue !KSHTEST 20170524-8:17S      
      do i3=1,3
         do j3=1,3
            Y(i3,j3,-3:3,-3:3)=-alpha_coefK*
     &           Y(i3,j3,-3:3,-3:3)
         enddo
            Y(i3,i3,-3:3,-3:3)=
     &           Y(i3,i3,-3:3,-3:3)+X(-3:3,-3:3)
         do j3=1,3
            Y(i3,j3,0,0)=Y(i3,j3,0,0)+(1.0-X00)*
     &           aamat_intmp(i3,j3)
         enddo
      enddo

      do i3=1,3
      noj=noj+1
      iacsrfibershell(noj)=ntemi+1
      do j3=1,3
       do i2=-3,3
        do j2=-3,3
            if(abs(Y(i3,j3,i2,j2)).gt. 1.d-17) then
              ntemi=ntemi+1
              jacsrfibershell(ntemi)=iAdd(i2,j2,j3)
              aacsrfibershell(ntemi)=Y(i3,j3,i2,j2)
            endif
        enddo
       enddo
       enddo !j3
      enddo !i3

! if piezoelectric
!      call PiezoforceF(
!     &       ds1(ihm2:ihp2,jhm2:jhp2)
!     &      ,ds2(ihm2:ihp2,jhm2:jhp2)
!     &      ,dg1(ihm2:ihp2,jhm2:jhp2)
!     &      ,dg2(ihm2:ihp2,jhm2:jhp2)
!     &      ,dh1(ihm2:ihp2,jhm2:jhp2)
!     &      ,dh2(ihm2:ihp2,jhm2:jhp2)
!     &      ,GX(1,ihm2:ihp2,jhm2:jhp2)
!     &      ,GX(2,ihm2:ihp2,jhm2:jhp2)
!     &      ,GX(3,ihm2:ihp2,jhm2:jhp2)
!     &      ,al(ihm2:ihp2,jhm2:jhp2)
!     &      ,Vpiezo(ihm2:ihp2,jhm2:jhp2)
!     &      ,VecPiezo(1:3,1:4,ihm1:ihp1,jhm1:jhp1) !NX
!     &      ,A0,TAG0)      


!              VECpiezo_tmp(1:3,26,jq,j)=normalpiezo(1:3)                !26: nrmalpiezo
!              VECpiezo_tmp(1:3,27,jq,j)=Rpiezo*(/Epiezo,Gpiezo,Fpiezo/) !27:[E R,G R,F R]
!              VECpiezo_tmp(1:3,28,jq,j)=Ppiezo*(/Lpiezo,Mpiezo,Npiezo/) !28:[E R,G R,F R]
!              VECpiezo_tmp(1:3,29,jq,j)=Ppiezo*(/Epiezo,Gpiezo,Fpiezo/)  !6:[E P,G P,F P]

!     &      ,VecPiezo(2,ihm1:ihp2,jhm2:jhp2) !NY
!     &      ,VecPiezo(3,ihm1:ihp2,jhm2:jhp2) !NZ
!     &      ,VecPiezo(4,ihm1:ihp2,jhm2:jhp2) !EP
!     &      ,VecPiezo(5,ihm1:ihp2,jhm2:jhp2) !ER
!     &      ,VecPiezo(6,ihm1:ihp2,jhm2:jhp2) !FP
!     &      ,VecPiezo(7,ihm1:ihp2,jhm2:jhp2) !FR
!     &      ,VecPiezo(8,ihm1:ihp2,jhm2:jhp2) !GP
!     &      ,VecPiezo(9,ihm1:ihp2,jhm2:jhp2) !GR
!     &      ,VecPiezo(10,ihm1:ihp2,jhm2:jhp2) !LP
!     &      ,VecPiezo(11,ihm1:ihp2,jhm2:jhp2) !NP
!     &      ,VecPiezo(12,ihm1:ihp2,jhm2:jhp2) !MP



      enddo !j
      enddo !i
      noj=noj+1
      iacsrfibershell(noj)=ntemi+1 

      if(ntemi .ge. ntermiMAX)  stop " ntemi is more than ntermiMAX"
      return 
      end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
