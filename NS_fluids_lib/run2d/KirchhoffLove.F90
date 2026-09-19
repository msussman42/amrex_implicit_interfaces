      program main
      IMPLICIT NONE
      real*8, dimension(:,:), allocatable :: w
      real*8, dimension(:,:), allocatable :: load
      real*8, dimension(:,:), allocatable :: Aw
      real*8, dimension(:,:), allocatable :: resid
      real*8, dimension(:,:), allocatable :: diagonal
      real*8 my_pi
      real*8 H,E,gamma,force
      integer polar
      integer i,j
      integer istride,jstride
      integer ibase,jbase
      real*8 prob_hi(2)
      real*8 dx(2)
      integer n_cell(2)
      integer n_cell_load(2)
      integer ngrow
      integer niter
      integer error_count
      real*8 load_area
      real*8 error
      real*8 error_tol

 
      my_pi=4.0d0*atan(1.0d0)
      H=1.0D0  !thickness
      E=1.0D+6  !Young's modulus
      gamma=0.25D0  !Poisson ratio
      force=1.0d0 !Force
      polar=1 !xy? or r?
      probhi(1)=1.0D0
      probhi(2)=1.0D0
      n_cell(1)=32 
      n_cell(2)=32 
      n_cell_load(1)=2
      n_cell_load(2)=2
      do dir=1,2
       dx(dir)=probhi(dir)/n_cell(dir)
      enddo
      ngrow=4
     
      allocate(w(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow))
      allocate(Aw(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow))
      allocate(resid(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow))
      allocate(load(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow))
      allocate(diagonal(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow))
      load=0.0d0
      w=0.0d0
      Aw=0.0d0
      resid=0.0d0
      diagonal=0.0d0
      
      if (polar.eq.1) then
       load_area=my_pi*(n_cell_load(1)*dx(1))**2 
      else if (polar.eq.0) then
       load_area=n_cell_load(1)*n_cell_load(2)*dx(1)*dx(2)
      else
       print *,"polar invalid"
       stop
      endif
      do i=0,n_cell_load(1)-1
      do j=0,n_cell_load(2)-1
       load(i,j)=force/load_area
      enddo
      enddo

      if (polar.eq.1) then
       do istride=0,3
        w=0.0d0
        do i=0,n_cell(1)-1
         ibase=i/4
         if (i-ibase*4.eq.istride) then
          w(i,0)=1.0d0
         endif
        enddo !i
        call apply(w,Aw,n_cell,dx,polar)
        do i=0,n_cell(1)-1
         ibase=i/4
         if (i-ibase*4.eq.istride) then
          diagonal(i,0)=Aw(i,0)
         endif
        enddo !i
       enddo !istride
      else if (polar.eq.0) then
       do istride=0,3
       do jstride=0,3
        w=0.0d0
        do i=0,n_cell(1)-1
        do j=0,n_cell(2)-1
         ibase=i/4
         jbase=j/4
         if ((i-ibase*4.eq.istride).and. &
             (j-jbase*4.eq.jstride)) then
          w(i,j)=1.0d0
         endif
        enddo !j
        enddo !i
        call apply(w,Aw,n_cell,dx,polar)
        do i=0,n_cell(1)-1
        do j=0,n_cell(2)-1
         ibase=i/4
         jbase=j/4
         if ((i-ibase*4.eq.istride).and. &
             (j-jbase*4.eq.jstride)) then
          diagonal(i,j)=Aw(i,j)
         endif
        enddo !i
        enddo !j
       enddo !jstride
       enddo !istride
      else
       print *,"polar invalid"
       stop
      endif
      error=1.0D+10
      error_tol=1.0D-8
      w=0.0d0
      niter=0
      do while (error.gt.error_tol)
       call residual(w,resid,load,n_cell,dx,polar) !load-Aw
       error=0.0d0
       error_count=0
       do i=0,n_cell(1)-1
       do j=0,n_cell(2)-1
        if (diagonal(i,j).gt.0.0d0) then
         w(i,j)=w(i,j)+resid(i,j)/diagonal(i,j)
        else
         print *,"diagonal invalid"
         stop
        endif
        error_count=error_count+1
        error=error+resid(i,j)**2
       enddo !j
       enddo !i
       error=sqrt(error/error_count)
       print *,"niter,error_count,error ",niter,error_count,error
       niter=niter+1
      enddo !while error>error_tol

      if (polar.eq.0) then
       print *,'VARIABLES="x","y","w"'
       print *,"zone i= ",n_cell(1)+2," j= ",n_cell(2)+2," f=point"
       print *,"solutiontime= 0.0  strandid=1"
       do i=-1,n_cell(1)
       do j=-1,n_cell(2)
        x=(i+0.5d0)*dx(1)  
        y=(j+0.5d0)*dx(2) 
        print *,x,y,w(i,j)
       enddo 
       enddo 
      else if (polar.eq.1) then
       do i=-1,n_cell(1)
        x=(i+0.5d0)*dx(1)  
        print *,x,w(i,0)
       enddo 
      else
       print *,"polar invalid"
       stop
      endif

      deallocate(w)
      deallocate(load)
      deallocate(Aw)
      deallocate(resid)
      deallocate(diagonal)

      return
      end

