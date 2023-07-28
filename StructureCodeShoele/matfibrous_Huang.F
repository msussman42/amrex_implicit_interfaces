c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine matfibrous_Huang(
     &         ntermiMAXrow,ntermiMAX
     &        ,iacsrfibershell,jacsrfibershell,aacsrfibershell
     &        ,nq_ibmin,ns_ibmin
     &        ,dt2
     &        ,density_coefin,cs_ibmin,density_coefpin
     &        ,BoundaryFibShellin
     &        ,Stiff_link_fsh
     &        ,Sig11  
     &        ,Sig12 
     &        ,Sig21 
     &        ,Sig22           
     &        ,gam11  
     &        ,gam12
     &        ,gam21
     &        ,gam22
     &        ,dh1  
     &        ,ds1  
     &        ,dg1 
     &        ,dk1    
     &        ,dh2  
     &        ,ds2  
     &        ,dg2 
     &        ,dk2
     &        ,iflagsimulation    
     &        )
      implicit none
      integer nq_ibmin,ns_ibmin,ntermiMAXrow,ntermiMAX
      integer iacsrfibershell(ntermiMAXrow+1)
     &       ,jacsrfibershell(ntermiMAX)
      real*8  aacsrfibershell(ntermiMAX)
      real*8  dt2,density_coefin,cs_ibmin,density_coefpin,dt
      integer BoundaryFibShellin(2,2)
      real*8  Stiff_link_fsh(nq_ibmin,ns_ibmin)
     &        ,Sig11(nq_ibmin,ns_ibmin)
     &        ,Sig12(nq_ibmin,ns_ibmin) 
     &        ,Sig21(nq_ibmin,ns_ibmin) 
     &        ,Sig22(nq_ibmin,ns_ibmin)           
     &        ,gam11(nq_ibmin,ns_ibmin)  
     &        ,gam12(nq_ibmin,ns_ibmin)
     &        ,gam21(nq_ibmin,ns_ibmin)
     &        ,gam22(nq_ibmin,ns_ibmin)
     &        ,dh1(nq_ibmin,ns_ibmin)
     &        ,ds1(nq_ibmin,ns_ibmin)  
     &        ,dg1(nq_ibmin,ns_ibmin) 
     &        ,dk1(nq_ibmin,ns_ibmin)    
     &        ,dh2(nq_ibmin,ns_ibmin)  
     &        ,ds2(nq_ibmin,ns_ibmin)  
     &        ,dg2(nq_ibmin,ns_ibmin) 
     &        ,dk2(nq_ibmin,ns_ibmin)
     
      integer i,j,im1,im2,im3,jm1,jm2,jm3,ihm1,ihm2,ihm3,jhm1,jhm2,jhm3
      integer ip1,ip2,ip3,jp1,jp2,jp3,ihp1,ihp2,ihp3,jhp1,jhp2,jhp3
      real*8 X(-3:3,-3:3),coef,X00
      integer iAdd(-3:3,-3:3), iiAdd(-3:3),jiAdd(-3:3),i2,j2,ntemi,noj
     &  ,ist,ien,jst,jen,iflagsimulation


!! ================== Inside i,j  2<i<nx-1; 2<j<ny-1
      dt=sqrt(dt2)
      iacsrfibershell=0
      jacsrfibershell=0
      aacsrfibershell=0.0

      ntemi=0;noj=0
      do i=1,nq_ibmin
      do j=1,ns_ibmin
         noj=noj+1
         iacsrfibershell(noj)=ntemi+1  
         ist=1
         ien=nq_ibmin

         jst=1
         jen=ns_ibmin
         
         if(BoundaryFibShellin(1,1) .eq. 3) ist=1
         if(BoundaryFibShellin(1,2) .eq. 3) ien=nq_ibmin
            


         if(BoundaryFibShellin(2,1) .eq. 3) jst=1
         if(BoundaryFibShellin(2,2) .eq. 3) jen=ns_ibmin


         
         iAdd(0,0)=(i-1)*ns_ibmin+j

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
            iAdd(i2,j2)=(iiAdd(i2)-1)*ns_ibmin+jiAdd(j2)
         enddo
         enddo         

!         pause
         X=0.d0
         if(iflagsimulation==1)then
         coef=density_coefin+density_coefpin+cs_ibmin/2.0*dt
     &       +Stiff_link_fsh(i,j)
         else
         coef=0.0
         endif


       X00=coef
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
            X(i2,j2)=-dt2*X(i2,j2)
        enddo
       enddo
       X(0,0)=X(0,0)+X00
       do i2=-3,3
        do j2=-3,3
            if(abs(X(i2,j2)).gt. 0.d0) then
              ntemi=ntemi+1
              jacsrfibershell(ntemi)=iAdd(i2,j2)
              aacsrfibershell(ntemi)=X(i2,j2)
            endif
        enddo
       enddo

      enddo !j
      enddo !i
      noj=noj+1
      iacsrfibershell(noj)=ntemi+1 
      return 
      end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
