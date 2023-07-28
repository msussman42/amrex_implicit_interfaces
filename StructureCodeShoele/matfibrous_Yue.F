c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine matfibrous_Yue(
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
     &        ,SigB
     &        ,e4coefint 
     &        ,iflagsimulation
     &        ,FlagForceCalculation  
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
     &       ,jm5,jp5,jhm5,ihm5 
      integer im4,ip4,ihm4,ihp4,jm4,jp4,jhm4,jhp4,jhp5,ihp5,ip5,im5
      integer ip1,ip2,ip3,jp1,jp2,jp3,ihp1,ihp2,ihp3,jhp1,jhp2,jhp3
      real*8 X(-4:4,-4:4),coef,X00,Xp(-4:4,-4:4)
      real*8 SigB(4,nq_ibmin+ns_ibmin,3)
     &      ,e4coefint,coefX00,e5c
      integer FlagForceCalculation(nq_ibmin,ns_ibmin)
      integer iAdd(-4:4,-4:4), iiAdd(-4:4),jiAdd(-4:4),i2,j2,ntemi,noj
     &  ,ist,ien,jst,jen,iflagsimulation,iAddp(-4:4,-4:4),ksecond,idata
     &  ,IDatap(100)
      real*8 Datap(100),e6e
      integer formulationFlaga,formulationFlagb
      parameter(formulationFlaga=1)
      parameter(formulationFlagb=1)

!! ================== Inside i,j  2<i<nx-1; 2<j<ny-1
      dt=sqrt(dt2)
      iacsrfibershell=0
      jacsrfibershell=0
      aacsrfibershell=0.0

      if((e4coefint .gt. 50.0) .and. (e4coefint .lt. 100.0)) 
     &   e4coefint=e4coefint-50.0

      e5c=2.0-e4coefint 
      e6e=e4coefint

      ntemi=0;noj=0
      do i=1,nq_ibmin
      do j=1,ns_ibmin
         coefX00=1.0
         noj=noj+1
         iacsrfibershell(noj)=ntemi+1  
         ist=1
         ien=nq_ibmin

         jst=1
         jen=ns_ibmin
         
         if(BoundaryFibShellin(1,1) .eq. 3) ist=1
         if(BoundaryFibShellin(1,2) .eq. 3) ien=nq_ibmin
            
         FlagForceCalculation(i,j)=0

         if(BoundaryFibShellin(2,1) .eq. 3) jst=1
         if(BoundaryFibShellin(2,2) .eq. 3) jen=ns_ibmin

         ksecond=0
         
         iAdd(0,0)=(i-1)*ns_ibmin+j

         im1=mod(i-1+nq_ibmin,nq_ibmin)
         if(im1 .eq. 0) im1=nq_ibmin
         im2=mod(i-2+nq_ibmin,nq_ibmin)
         if(im2 .eq. 0) im2=nq_ibmin
         im3=mod(i-3+nq_ibmin,nq_ibmin)
         if(im3 .eq. 0) im3=nq_ibmin
         im4=mod(i-4+nq_ibmin,nq_ibmin)
         if(im4 .eq. 0) im4=nq_ibmin
         im5=mod(i-5+nq_ibmin,nq_ibmin)
         if(im5 .eq. 0) im5=nq_ibmin

         ip1=mod(i+1+nq_ibmin,nq_ibmin)
         if(ip1 .eq. 0) ip1=nq_ibmin
         ip2=mod(i+2+nq_ibmin,nq_ibmin)
         if(ip2 .eq. 0) ip2=nq_ibmin
         ip3=mod(i+3+nq_ibmin,nq_ibmin)
         if(ip3 .eq. 0) ip3=nq_ibmin
         ip4=mod(i+4+nq_ibmin,nq_ibmin)
         if(ip4 .eq. 0) ip4=nq_ibmin
         ip5=mod(i+5+nq_ibmin,nq_ibmin)
         if(ip5 .eq. 0) ip5=nq_ibmin


         ihm1=im1
         ihm2=im2
         ihm3=im3
         ihm4=im4
         ihm5=im5

         ihp1=ip1
         ihp2=ip2
         ihp3=ip3
         ihp4=ip4
         ihp5=ip5

         jm1=mod(j-1+ns_ibmin,ns_ibmin)
         if(jm1 .eq. 0) jm1=ns_ibmin
         jm2=mod(j-2+ns_ibmin,ns_ibmin)
         if(jm2 .eq. 0) jm2=ns_ibmin
         jm3=mod(j-3+ns_ibmin,ns_ibmin)
         if(jm3 .eq. 0) jm3=ns_ibmin
         jm4=mod(j-4+ns_ibmin,ns_ibmin)
         if(jm4 .eq. 0) jm4=ns_ibmin
         jm5=mod(j-5+ns_ibmin,ns_ibmin)
         if(jm5 .eq. 0) jm5=ns_ibmin

         jp1=mod(j+1+ns_ibmin,ns_ibmin)
         if(jp1 .eq. 0) jp1=ns_ibmin
         jp2=mod(j+2+ns_ibmin,ns_ibmin)
         if(jp2 .eq. 0) jp2=ns_ibmin
         jp3=mod(j+3+ns_ibmin,ns_ibmin)
         if(jp3 .eq. 0) jp3=ns_ibmin
         jp4=mod(j+4+ns_ibmin,ns_ibmin)
         if(jp4 .eq. 0) jp4=ns_ibmin
         jp5=mod(j+5+ns_ibmin,ns_ibmin)
         if(jp5 .eq. 0) jp5=ns_ibmin

         jhm1=jm1
         jhm2=jm2
         jhm3=jm3
         jhm4=jm4
         jhm5=jm5
      
         jhp1=jp1
         jhp2=jp2
         jhp3=jp3
         jhp4=jp4
         jhp5=jp5
         
         iiAdd(-4:4)=(/im4,im3,im2,im1,i,ip1,ip2,ip3,ip4/)
         jiAdd(-4:4)=(/jm4,jm3,jm2,jm1,j,jp1,jp2,jp3,jp4/)

         do i2=-4,4
         do j2=-4,4
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
       coefX00=0.0d0   
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
       coefX00=0.0d0   
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
       coefX00=0.0d0 

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
       coefX00=0.0d0 
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

 
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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
  
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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

    
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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

      
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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

    
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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

 
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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

  
      FlagForceCalculation(i,j)=1
      coefX00=1.0d0
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

 
      FlagForceCalculation(i,j)=1
      coefX00=1.0d0

!! ================== BNDY 2,N  
       elseif((i.eq. 2) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (j .eq. ns_ibmin)) then
!= ===============  2,N VER2
!CCCCCCCC+X-xCCCCCCCC

      X(-1,-2)=(0.25*sigB(4,im1,3))/(dh1(im1,j)*(0.5*ds2(im1,jm2) - 1.
     &5*ds2(im1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(im1,jm2)*d
     &s2(i,jm2)*ds2(i,jm1))

      X(-1,-1)=(0.5*sig12(im1,jm1))/(ds1(im1,jm1)*ds2(i,jm1)) - (1.0*s
     &igB(4,im1,3))/(dh1(im1,j)*(0.5*ds2(im1,jm2) - 1.5*ds2(im1,jm1)))
     & + (2.0*gam22(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) - (2
     &.0*e6e*gam21(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) + (2.
     &0*gam21(i,jm1))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm2)*ds2(i,jm1))

      X(-1,0)=gam11(i,j)/(dh1(i,j)**2*ds1(im1,j)**2) + sigB(4,im1,1)/(
     &dh1(i,j)*ds1(im1,j)) + (0.5*sig12(im1,jm1))/(ds1(im1,j)*ds2(i,jm
     &1)) + (0.75*sigB(4,im1,3))/(dh1(im1,j)*(0.5*ds2(im1,jm2) - 1.5*d
     &s2(im1,jm1))) + (gam22(i,j)*(e4coefint - 2.0))/(dh1(i,j)*ds1(im1
     &,j)*ds2(i,jm1)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*d
     &s2(i,jm1)**2)

      X(0,-3)=-(1.0*gam22(i,jm1))/(dh2(i,jm2)*ds2(i,jm3)*ds2(i,jm2)*ds
     &2(i,jm1))

      X(0,-2)=((0.25*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) 
     &- (0.25*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(im1,
     &j) + gam22(i,j)/(dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) + (1.0*gam
     &22(i,jm1)*((1.0/ds2(i,jm3) + 1.0/ds2(i,jm2))/dh2(i,jm2) + 1.0/(d
     &h2(i,jm1)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) + (2.0*gam21(i,j
     &m1)*(1.0/ds1(im1,jm2) + 1.0/ds1(i,jm2)))/(dh1(i,jm2)*ds2(i,jm2)*
     &ds2(i,jm1))

      X(0,-1)=(2.0*sig22(i,jm1))/ds2(i,jm1)**2 - ((1.0*sigB(4,im1,3))/
     &(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) - (1.0*sigB(4,i,3))/(0.5*ds2(i
     &,jm2) - 1.5*ds2(i,jm1)))/dh1(im1,j) - (2.0*((0.25*sig12(im1,jm1)
     &)/ds1(im1,jm1) - (0.25*sig12(i,jm1))/ds1(i,jm1)))/ds2(i,jm1) - (
     &gam22(i,j)*((2.0*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/dh1(i,jm1)
     & + (1.0/ds2(i,jm2) + 1.0/ds2(i,jm1))/dh2(i,jm1)))/ds2(i,jm1)**2 
     &- (1.0*gam22(i,jm1)*((1.0*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2
     &(i,jm1) + 1.0/(dh2(i,jm2)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) 
     &- (2.0*gam21(i,jm1)*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/(dh1(i,
     &jm1)*ds2(i,jm2)*ds2(i,jm1)) + (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,j
     &m1) + 1.0/ds1(i,jm1)))/(dh1(i,jm1)*ds2(i,jm1)**2)

      X(0,0)=((0.75*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) -
     & (0.75*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(im1,j
     &) - (2.0*((0.25*sig12(im1,jm1))/ds1(im1,j) - (0.25*sig12(i,jm1))
     &/ds1(i,j)))/ds2(i,jm1) - ((1.0*gam11(ip1,j))/(dh1(ip1,j)*ds1(i,j
     &)**2) + (gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*d
     &s1(im1,j)))/dh1(i,j) - (2.0*sig22(i,jm1))/ds2(i,jm1)**2 - ((1.0*
     &sigB(4,im1,1))/ds1(im1,j) + (1.0*sigB(4,i,1))/ds1(i,j))/dh1(i,j)
     & + (1.0*gam22(i,j))/(dh2(i,jm1)*ds2(i,jm1)**3) + (1.0*gam22(i,jm
     &1))/(dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) - (gam22(i,j)*(e4coefi
     &nt - 2.0)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds2(i,jm1)*
     &*2) - (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(
     &i,j)*ds2(i,jm1)**2)

      X(+1,-2)=- (0.25*sigB(4,i,3))/(dh1(im1,j)*(0.5*ds2(ip1,jm2) - 1.
     &5*ds2(ip1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(i,jm2)*ds2
     &(i,jm2)*ds2(i,jm1))

      X(+1,-1)=(1.0*sigB(4,i,3))/(dh1(im1,j)*(0.5*ds2(ip1,jm2) - 1.5*d
     &s2(ip1,jm1))) - (0.5*sig12(i,jm1))/(ds1(i,jm1)*ds2(i,jm1)) + (2.
     &0*gam22(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) - (2.0*e6e*g
     &am21(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) + (2.0*gam21(i,
     &jm1))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm2)*ds2(i,jm1))

      X(+1,0)=((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)) + (1.0*
     &gam11(ip1,j)*(1.0/ds1(i,j) + 1.0/ds1(ip1,j)))/(dh1(ip1,j)*ds1(i,
     &j)))/dh1(i,j) + (1.0*sigB(4,i,1))/(dh1(i,j)*ds1(i,j)) - (0.5*sig
     &12(i,jm1))/(ds1(i,j)*ds2(i,jm1)) - (0.75*sigB(4,i,3))/(dh1(im1,j
     &)*(0.5*ds2(ip1,jm2) - 1.5*ds2(ip1,jm1))) + (1.0*gam22(i,j)*(e4co
     &efint - 2.0))/(dh1(i,j)*ds1(i,j)*ds2(i,jm1)**2) + (2.0*e6e*gam21
     &(i,j))/(dh1(i,j)*ds1(i,j)*ds2(i,jm1)**2)

      X(+2,0)=-(1.0*gam11(ip1,j))/(dh1(i,j)*dh1(ip1,j)*ds1(i,j)*ds1(ip
     &1,j))

      FlagForceCalculation(i,j)=1
      coefX00=1.0d0
!! ================== BNDY N-1,N  
       elseif((i.eq. nq_ibmin-1) .and. (BoundaryFibShellin(1,2) .eq. 1) 
     &    .and. (j .eq. ns_ibmin)) then
!= ===============  N-1,N VER2
!CCCCCCCC+X-xCCCCCCCC

      X(-2,0)=-(1.0*gam11(im1,j))/(dh1(im1,j)*dh1(i,j)*ds1(im2,j)*ds1(
     &im1,j))

      X(-1,-2)=(0.25*sigB(4,im1,3))/(dh1(im1,j)*(0.5*ds2(im1,jm2) - 1.
     &5*ds2(im1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(im1,jm2)*d
     &s2(i,jm2)*ds2(i,jm1))

      X(-1,-1)=(0.5*sig12(im1,jm1))/(ds1(im1,jm1)*ds2(i,jm1)) - (1.0*s
     &igB(4,im1,3))/(dh1(im1,j)*(0.5*ds2(im1,jm2) - 1.5*ds2(im1,jm1)))
     & + (2.0*gam22(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) - (2
     &.0*e6e*gam21(i,j))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm1)**2) + (2.
     &0*gam21(i,jm1))/(dh1(i,jm1)*ds1(im1,jm1)*ds2(i,jm2)*ds2(i,jm1))

      X(-1,0)=(1.0*((1.0*gam11(i,j))/(dh1(i,j)*ds1(im1,j)*ds1(i,j)) + 
     &(gam11(im1,j)*(1.0/ds1(im2,j) + 1.0/ds1(im1,j)))/(dh1(im1,j)*ds1
     &(im1,j))))/dh1(i,j) + sigB(4,im1,1)/(dh1(i,j)*ds1(im1,j)) + (0.5
     &*sig12(im1,jm1))/(ds1(im1,j)*ds2(i,jm1)) + (0.75*sigB(4,im1,3))/
     &(dh1(im1,j)*(0.5*ds2(im1,jm2) - 1.5*ds2(im1,jm1))) + (gam22(i,j)
     &*(e4coefint - 2.0))/(dh1(i,j)*ds1(im1,j)*ds2(i,jm1)**2) + (2.0*e
     &6e*gam21(i,j))/(dh1(i,j)*ds1(im1,j)*ds2(i,jm1)**2)

      X(0,-3)=-(1.0*gam22(i,jm1))/(dh2(i,jm2)*ds2(i,jm3)*ds2(i,jm2)*ds
     &2(i,jm1))

      X(0,-2)=((0.25*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) 
     &- (0.25*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(im1,
     &j) + gam22(i,j)/(dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) + (1.0*gam
     &22(i,jm1)*((1.0/ds2(i,jm3) + 1.0/ds2(i,jm2))/dh2(i,jm2) + 1.0/(d
     &h2(i,jm1)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) + (2.0*gam21(i,j
     &m1)*(1.0/ds1(im1,jm2) + 1.0/ds1(i,jm2)))/(dh1(i,jm2)*ds2(i,jm2)*
     &ds2(i,jm1))

      X(0,-1)=(2.0*sig22(i,jm1))/ds2(i,jm1)**2 - ((1.0*sigB(4,im1,3))/
     &(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) - (1.0*sigB(4,i,3))/(0.5*ds2(i
     &,jm2) - 1.5*ds2(i,jm1)))/dh1(im1,j) - (2.0*((0.25*sig12(im1,jm1)
     &)/ds1(im1,jm1) - (0.25*sig12(i,jm1))/ds1(i,jm1)))/ds2(i,jm1) - (
     &gam22(i,j)*((2.0*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/dh1(i,jm1)
     & + (1.0/ds2(i,jm2) + 1.0/ds2(i,jm1))/dh2(i,jm1)))/ds2(i,jm1)**2 
     &- (1.0*gam22(i,jm1)*((1.0*(1.0/ds2(i,jm2) + 1.0/ds2(i,jm1)))/dh2
     &(i,jm1) + 1.0/(dh2(i,jm2)*ds2(i,jm2))))/(ds2(i,jm2)*ds2(i,jm1)) 
     &- (2.0*gam21(i,jm1)*(1.0/ds1(im1,jm1) + 1.0/ds1(i,jm1)))/(dh1(i,
     &jm1)*ds2(i,jm2)*ds2(i,jm1)) + (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,j
     &m1) + 1.0/ds1(i,jm1)))/(dh1(i,jm1)*ds2(i,jm1)**2)

      X(0,0)=((0.75*sigB(4,im1,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)) -
     & (0.75*sigB(4,i,3))/(0.5*ds2(i,jm2) - 1.5*ds2(i,jm1)))/dh1(im1,j
     &) - (2.0*((0.25*sig12(im1,jm1))/ds1(im1,j) - (0.25*sig12(i,jm1))
     &/ds1(i,j)))/ds2(i,jm1) - (2.0*sig22(i,jm1))/ds2(i,jm1)**2 - ((1.
     &0*sigB(4,im1,1))/ds1(im1,j) + (1.0*sigB(4,i,1))/ds1(i,j))/dh1(i,
     &j) - (1.0*((1.0*gam11(im1,j))/(dh1(im1,j)*ds1(im1,j)**2) + (1.0*
     &gam11(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*ds1(i,j)))
     &)/dh1(i,j) + (1.0*gam22(i,j))/(dh2(i,jm1)*ds2(i,jm1)**3) + (1.0*
     &gam22(i,jm1))/(dh2(i,jm1)*ds2(i,jm2)*ds2(i,jm1)**2) - (gam22(i,j
     &)*(e4coefint - 2.0)*(1.0/ds1(im1,j) + 1.0/ds1(i,j)))/(dh1(i,j)*d
     &s2(i,jm1)**2) - (2.0*e6e*gam21(i,j)*(1.0/ds1(im1,j) + 1.0/ds1(i,
     &j)))/(dh1(i,j)*ds2(i,jm1)**2)

      X(+1,-2)=- (0.25*sigB(4,i,3))/(dh1(im1,j)*(0.5*ds2(ip1,jm2) - 1.
     &5*ds2(ip1,jm1))) - (2.0*gam21(i,jm1))/(dh1(i,jm2)*ds1(i,jm2)*ds2
     &(i,jm2)*ds2(i,jm1))

      X(+1,-1)=(1.0*sigB(4,i,3))/(dh1(im1,j)*(0.5*ds2(ip1,jm2) - 1.5*d
     &s2(ip1,jm1))) - (0.5*sig12(i,jm1))/(ds1(i,jm1)*ds2(i,jm1)) + (2.
     &0*gam22(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) - (2.0*e6e*g
     &am21(i,j))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm1)**2) + (2.0*gam21(i,
     &jm1))/(dh1(i,jm1)*ds1(i,jm1)*ds2(i,jm2)*ds2(i,jm1))

      X(+1,0)=(1.0*gam11(i,j))/(dh1(i,j)**2*ds1(i,j)**2) + (1.0*sigB(4
     &,i,1))/(dh1(i,j)*ds1(i,j)) - (0.5*sig12(i,jm1))/(ds1(i,j)*ds2(i,
     &jm1)) - (0.75*sigB(4,i,3))/(dh1(im1,j)*(0.5*ds2(ip1,jm2) - 1.5*d
     &s2(ip1,jm1))) + (1.0*gam22(i,j)*(e4coefint - 2.0))/(dh1(i,j)*ds1
     &(i,j)*ds2(i,jm1)**2) + (2.0*e6e*gam21(i,j))/(dh1(i,j)*ds1(i,j)*d
     &s2(i,jm1)**2)

    
      FlagForceCalculation(i,j)=1
      coefX00=1.0d0
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

 
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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

     
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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

     
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0

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

    
       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
!! ================== BNDY N,j 2<=j<=Ny-1
       elseif((i.eq. nq_ibmin)   
     &        .and. (BoundaryFibShellin(2,2) .eq. 1) 
     &    .and. (j .gt. jst) .and. (j.lt. jen)) then
!= ===============  N,j VER2
!CCCCCCCC+X-xCCCCCCCC 
      print*, 'here'
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

 
      FlagForceCalculation(i,j)=1
      coefX00=1.0d0
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

      FlagForceCalculation(i,j)=1
      coefX00=1.0d0
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

       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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

       FlagForceCalculation(i,j)=1
       coefX00=1.0d0
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

       FlagForceCalculation(i,j)=1
       coefX00=1.0d0

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

       FlagForceCalculation(i,j)=1
       coefX00=1.0d0

!other points
       else
       FlagForceCalculation(i,j)=1
       coefX00=1.0
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

       endif
       if(       FlagForceCalculation(i,j) .eq. 1)then
         do i2=-4,4
          do j2=-4,4
              X(i2,j2)=-dt2*X(i2,j2)
          enddo
         enddo
         if(ksecond .eq. 1) then
           Datap(1:idata)=-dt2*Datap(1:idata)
         endif
       endif
       X(0,0)=X(0,0)+X00*coefX00
       if(ksecond .eq. 1) then 
          do i2=1,idata         
              ntemi=ntemi+1
              jacsrfibershell(ntemi)=IDatap(i2)
              aacsrfibershell(ntemi)=Datap(i2)
           enddo           
       else
         do i2=-4,4
          do j2=-4,4
            if(abs(X(i2,j2)).gt. 1.d-17) then
              ntemi=ntemi+1
              jacsrfibershell(ntemi)=iAdd(i2,j2)
              aacsrfibershell(ntemi)=X(i2,j2)
            endif
          enddo
         enddo
       endif

      enddo !j
      enddo !i
      noj=noj+1
      iacsrfibershell(noj)=ntemi+1 
  
!      do i=1,noj
!      write(151,*) iacsrfibershell(i)
!      enddo
!      do i=1,ntemi
!        write(152,*) jacsrfibershell(i),aacsrfibershell(i)
!      enddo
!      pause 
200   format(5(i24,' '),28(f20.9,' '))
      return 
      end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
