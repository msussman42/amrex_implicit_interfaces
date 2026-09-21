
      subroutine set_ghost(w,n_cell,ngrow,polar)
      IMPLICIT NONE
      integer, intent(in) :: n_cell(2)
      integer, intent(in) :: ngrow
      integer, intent(in) :: polar
      real*8, intent(inout) :: &
             w(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1)
      integer ighost,jghost

      if (polar.eq.1) then
       do ighost=1,ngrow
        w(-ighost,0)=w(ighost-1,0)
        w(n_cell(1)+ighost-1,0)=-w(n_cell(1)-ighost,0)
       enddo
      else if (polar.eq.0) then
       do ighost=0,n_cell(1)-1
       do jghost=1,ngrow
        w(ighost,n_cell(2)+jghost-1)=w(ighost,n_cell(2)-jghost)
        w(ighost,-jghost)=w(ighost,jghost-1)
       enddo !jghost
       enddo !ighost
       do jghost=-ngrow,n_cell(2)+ngrow-1
       do ighost=1,ngrow
        w(n_cell(1)+ighost-1,jghost)=-w(n_cell(1)-ighost,jghost)
        w(-ighost,jghost)=w(ighost-1,jghost)
       enddo !ighost
       enddo !jghost
      else
       print *,"polar invalid"
       stop
      endif

      return
      end subroutine set_ghost


      subroutine apply_lap(w,lap_w,n_cell,ngrow,dx,polar)
      IMPLICIT NONE
      real*8, intent(in) :: dx(2)
      integer, intent(in) :: n_cell(2)
      integer, intent(in) :: ngrow
      integer, intent(in) :: polar
      real*8, intent(inout) :: &
             w(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1)
      real*8, intent(out) :: &
             lap_w(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1)
      integer i,j
      real*8 r,r_left,r_right,grad_left,grad_right

      call set_ghost(w,n_cell,ngrow,polar)
      if (polar.eq.1) then
       do i=0,n_cell(1)-1
        r=(i+0.5d0)*dx(1)
        r_left=r-0.5d0*dx(1)
        r_right=r+0.5d0*dx(1)
        if (r.gt.0.0d0) then
         !do nothing
        else
         print *,"r invalid"
         stop
        endif
        grad_left=(w(i,0)-w(i-1,0))/dx(1)
        grad_right=(w(i+1,0)-w(i,0))/dx(1)
        lap_w(i,0)=(r_right*grad_right-r_left*grad_left)/(r*dx(1))
       enddo
      else if (polar.eq.0) then
       do i=0,n_cell(1)-1
       do j=0,n_cell(2)-1
        lap_w(i,j)=(w(i+1,j)+w(i-1,j)-2.0d0*w(i,j))/(dx(1)**2)+ & 
                   (w(i,j+1)+w(i,j-1)-2.0d0*w(i,j))/(dx(2)**2) 
       enddo !j
       enddo !i
      else
       print *,"polar invalid"
       stop
      endif

      return
      end subroutine apply_lap


      subroutine apply(w,Aw,lap_w,n_cell,ngrow,dx,polar, &
                      acceleration_coefficient)
      IMPLICIT NONE
      real*8, intent(in) :: acceleration_coefficient
      real*8, intent(in) :: dx(2)
      integer, intent(in) :: n_cell(2)
      integer, intent(in) :: ngrow
      integer, intent(in) :: polar
      real*8, intent(inout) :: &
             w(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1)
      real*8, intent(out) :: &
             Aw(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1)
      real*8, intent(out) :: &
             lap_w(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1)
      integer i,j

      call set_ghost(w,n_cell,ngrow,polar)
      call apply_lap(w,lap_w,n_cell,ngrow,dx,polar)
      call apply_lap(lap_w,Aw,n_cell,ngrow,dx,polar)
      if (polar.eq.1) then
       do i=0,n_cell(1)-1
        Aw(i,0)=Aw(i,0)+acceleration_coefficient*w(i,0)
       enddo
      else if (polar.eq.0) then
       do i=0,n_cell(1)-1
       do j=0,n_cell(2)-1
        Aw(i,j)=Aw(i,j)+acceleration_coefficient*w(i,j)
       enddo !j
       enddo !i
      else
       print *,"polar invalid"
       stop
      endif

      return
      end subroutine apply

      program main
      IMPLICIT NONE
      real*8, dimension(:,:), allocatable :: w
      real*8, dimension(:,:), allocatable :: w_m1
      real*8, dimension(:,:), allocatable :: w_m2
      real*8, dimension(:,:), allocatable :: load
      real*8, dimension(:,:), allocatable :: Aw
      real*8, dimension(:,:), allocatable :: lap_w
      real*8, dimension(:,:), allocatable :: resid
      real*8, dimension(:,:), allocatable :: diagonal
      real*8 x,y
      real*8 my_pi
      real*8 H,E,D,gamma,force
      integer polar
      integer i,j
      integer dir
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
      real*8 time
      real*8 density
      real*8 drop_density
      real*8 drop_stop_time
      real*8 drop_initial_velocity
      real*8 drop_radius
      real*8 drop_mass
      real*8 acceleration
      real*8 acceleration_coefficient
      real*8 stop_time
      real*8 dt
      integer number_steps
      integer plot_int
      integer file_unit
      integer step_number
      character(len=64) :: filename

 
      my_pi=4.0d0*atan(1.0d0)
      density=1.03D0 ! density of PDMS
      H=0.1D0  !thickness=1mm=0.1cm
       !Shear modulus G=E/(2(1+gamma))
      E=1.0D+7  !Young's modulus 2.5E+6 Pascal=2.5E+7 dyne/cm^2
      gamma=0.25D0  !Poisson ratio
       !forward momentum of the drop is defeated in about 4ms
       !drop initial velocity is 72 cm/s
       !drop mass=density * volume
      drop_stop_time=4.0D-3
      stop_time=0.03 !30ms
      number_steps=100
      plot_int=10
      drop_initial_velocity=72.0d0
      drop_density=1.0d0
      drop_radius=0.28*0.5d0 !2.8 mm diameter
      drop_mass=drop_density*(4.0d0/3.0d0)*my_pi*(drop_radius**3)
      acceleration=drop_initial_velocity/drop_stop_time
      force=drop_mass*acceleration !Force
      print *,"drop density ",drop_density
      print *,"drop mass ",drop_mass
      print *,"drop_radius ",drop_radius
      print *,"drop_stop_time ",drop_stop_time
      print *,"drop_initial_velocity ",drop_initial_velocity
      print *,"acceleration ",acceleration

      D=2.0d0*(H**3)*E/(3.0d0*(1.0d0-gamma**2))
      if (D.gt.0.0D0) then
       !do nothing
      else
       print *,"D invalid"
       stop
      endif

      force=force/D
      acceleration_coefficient=density*H/D

      print *,"force ",force

      polar=0 !xy? or r?
      print *,"polar=",polar

      prob_hi(1)=0.5d0*4.0D0 !40 mm/2 (symmetry)
      prob_hi(2)=0.5d0*3.2D0 !32 mm/2 (symmetry)
      n_cell(1)=40 
      n_cell(2)=32 
      n_cell_load(1)=4
      n_cell_load(2)=4
      do dir=1,2
       dx(dir)=prob_hi(dir)/n_cell(dir)
      enddo
      ngrow=4
    
      if (polar.eq.1) then
       load_area=my_pi*(n_cell_load(1)*dx(1))**2 
       prob_hi(2)=1.0d0
       n_cell(2)=1
       n_cell_load(2)=1
       dx(2)=1.0d0
      else if (polar.eq.0) then
       ! solve quarter domain problem
       load_area=n_cell_load(1)*n_cell_load(2)*dx(1)*dx(2)*4.0d0
      else
       print *,"polar invalid"
       stop
      endif

      allocate(w(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1))
      allocate(w_m1(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1))
      allocate(w_m2(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1))
      allocate(Aw(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1))
      allocate(resid(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1))
      allocate(load(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1))
      allocate(diagonal(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1))
      allocate(lap_w(-ngrow:n_cell(1)+ngrow-1,-ngrow:n_cell(2)+ngrow-1))

      load=0.0d0
      w=0.0d0
      w_m1=0.0d0
      w_m2=0.0d0
      Aw=0.0d0
      resid=0.0d0
      diagonal=0.0d0
      lap_w=0.0d0
      
      time=0.0d0
      step_number=0
      dt=stop_time/number_steps
      if (dt.gt.0.0d0) then
       !do nothing
      else
       print *,"dt invalid"
       stop
      endif
      acceleration_coefficient=acceleration_coefficient/(dt*dt)

      do while (step_number.lt.number_steps)

       if (polar.eq.1) then
        do istride=0,3
         w=0.0d0
         do i=0,n_cell(1)-1
          ibase=i/4
          if (i-ibase*4.eq.istride) then
           w(i,0)=1.0d0
          endif
         enddo !i
         call apply(w,Aw,lap_w,n_cell,ngrow,dx,polar,acceleration_coefficient)
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
         call apply(w,Aw,lap_w,n_cell,ngrow,dx,polar,acceleration_coefficient)
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
       if (polar.eq.1) then
        if (n_cell(2).eq.1) then
         !do nothing
        else
         print *,"n_cell(2) invalid"
         stop
        endif
       else if (polar.eq.0) then
        !do nothing
       else
        print *,"polar invalid"
        stop
       endif
       load=0.0d0
       if (time.le.drop_stop_time) then
        do i=0,n_cell_load(1)-1
        do j=0,n_cell_load(2)-1
         load(i,j)=force/load_area
        enddo
        enddo
       endif

       do i=0,n_cell(1)-1
       do j=0,n_cell(2)-1
        load(i,j)=load(i,j)+acceleration_coefficient* &
             (2.0d0*w_m1(i,j)-w_m2(i,j))
       enddo !j
       enddo !i

       error=1.0D+10
       error_tol=1.0D-8
       w=0.0d0
       niter=0
       do while ((error.gt.error_tol).and.(niter.lt.200000))
        call apply(w,Aw,lap_w,n_cell,ngrow,dx,polar,acceleration_coefficient)
        do i=0,n_cell(1)-1
        do j=0,n_cell(2)-1
         resid(i,j)=load(i,j)-Aw(i,j)
        enddo !j
        enddo !i

        error=0.0d0
        error_count=0
        do i=0,n_cell(1)-1
        do j=0,n_cell(2)-1
         if (diagonal(i,j).gt.0.0d0) then
          w(i,j)=w(i,j)+0.5d0*resid(i,j)/diagonal(i,j)
         else
          print *,"diagonal invalid"
          stop
         endif
         error_count=error_count+1
         error=error+resid(i,j)**2
        enddo !j
        enddo !i
        error=sqrt(error/error_count)
        if (10000*(niter/10000).eq.niter) then
         print *,"niter,error_count,error ",niter,error_count,error
        endif
        niter=niter+1
       enddo !while error>error_tol

       if ((step_number/plot_int)*plot_int.eq.step_number) then

        write(filename,'(A,I3.3,A)') 'step_',step_number,'.tec'

        file_unit=10
        !file="output.txt" is an alternative
        open(unit=file_unit,file=trim(filename),status="replace",action="write")
        print *,"filename ",filename
        print *,"file_unit=",file_unit
        print *,"time=",time

        if (polar.eq.0) then
         write(file_unit,*) 'VARIABLES="x","y","w"'
         write(file_unit,*) "zone i= ",n_cell(1)+1," j= ",n_cell(2)+1," f=point"
         write(file_unit,*) "solutiontime=",time,"strandid=1"
         do j=0,n_cell(2)
         do i=0,n_cell(1)
          x=i*dx(1)  
          y=j*dx(2) 
          write(file_unit,*) x,y,0.25d0*(w(i,j)+w(i-1,j-1)+w(i-1,j)+w(i,j-1))
         enddo 
         enddo 
        else if (polar.eq.1) then
         do i=-1,n_cell(1)
          x=(i+0.5d0)*dx(1)  
          write(file_unit,*) x,w(i,0)
         enddo 
        else
         print *,"polar invalid"
         stop
        endif
   
        close(file_unit)
       endif

       step_number=step_number+1
       time=time+dt

       do i=0,n_cell(1)-1
       do j=0,n_cell(2)-1
        w_m2(i,j)=w_m1(i,j)
        w_m1(i,j)=w(i,j)
       enddo !j
       enddo !i
       w=0.0d0

      enddo  ! do while (step_number.lt.number_steps)

      deallocate(w)
      deallocate(w_m1)
      deallocate(w_m2)
      deallocate(load)
      deallocate(Aw)
      deallocate(resid)
      deallocate(diagonal)
      deallocate(lap_w)

      return
      end

