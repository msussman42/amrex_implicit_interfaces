 
      module variable_temperature_drop
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      REAL(kind=8) :: axisymmetric_disk_L
      REAL(kind=8) :: axisymmetric_disk_TSAT
      REAL(kind=8) :: axisymmetric_disk_TDIFF
      REAL(kind=8) :: axisymmetric_disk_thermal_diffusivity
      REAL(kind=8), dimension(2) :: axisymmetric_disk_radblob
      INTEGER :: axisymmetric_disk_stefan_flag
      INTEGER :: axisymmetric_disk_ireverse
      INTEGER :: axisymmetric_disk_isink
      REAL(kind=8) :: axisymmetric_disk_time

      REAL(kind=8),dimension(:,:),allocatable :: axisymmetric_disk_state

      INTEGER :: axisymmetric_disk_NCELL
      INTEGER :: axisymmetric_disk_NCYCLE
      INTEGER :: axisymmetric_polar_flag  ! 1=polar coord.  0=x coord.
      INTEGER :: axisymmetric_probtype

      contains

      subroutine get_dr(dr,time_index)
      IMPLICIT NONE
      
      REAL(kind=8) :: dr
      INTEGER  :: time_index

      if (axisymmetric_disk_NCELL.gt.0) then
       dr=axisymmetric_disk_radblob(time_index)/axisymmetric_disk_NCELL
      else
       print *,"axisymmetric_disk_NCELL invalid"
       stop
      endif

      end subroutine get_dr

       ! TDIFF=T_DISK_CENTER - TSAT 
      subroutine disk_eval_initial_temp(rval,temperature)
      IMPLICIT NONE

      REAL(kind=8) :: rval,temperature,RR,local_pi

      local_pi=4.0d0*atan(1.0d0)

      RR=axisymmetric_disk_radblob(1)
      if (RR.gt.0.0d0) then

       if (axisymmetric_probtype.eq.4) then

        if (rval.gt.RR) then
         temperature=axisymmetric_disk_TSAT
        else if ((rval.ge.0.0d0).and.(rval.le.RR)) then
           ! temperature(0)=TSAT+TDIFF
           ! temperature(RR)=TSAT
         temperature=axisymmetric_disk_TSAT+ &
          axisymmetric_disk_TDIFF*cos(0.5d0*rval*local_pi/RR)
        else
         print *,"rval invalid"
         stop
        endif
       else if (axisymmetric_probtype.eq.16) then
        if (rval.ge.0.2d0) then
         temperature=0.0d0
        else if ((rval.ge.0.0).and.(rval.le.0.2d0)) then
         temperature=5.0d0*cos(rval*local_pi/0.2d0)+5.0d0
        else if (rval.le.0.0d0) then
         temperature=10.0d0
        else
         print *,"rval invalid"
         stop
        endif
       else
        print *,"axisymmetric_probtype invalid"
        stop
       endif
      else
       print *,"RR invalid"
       stop
      endif
      if (temperature.ge.0.0d0) then
       ! do nothing
      else
       print *,"temperature invalid2"
       stop
      endif

      end subroutine disk_eval_initial_temp

      subroutine disk_interp_T(rval,time_index,temperature)
      IMPLICIT NONE

      REAL(kind=8) :: rval
      REAL(kind=8) :: temperature
      INTEGER :: time_index
      REAL(kind=8) :: RR,dr,ri,TLEFT,TRIGHT,RLEFT,RRIGHT,dr_interp
      REAL(kind=8) :: local_TWALL
      INTEGER :: i,ie
 
      RR=axisymmetric_disk_radblob(time_index)
      ie=axisymmetric_disk_NCELL-1

      if (rval.ge.RR) then
       temperature=axisymmetric_disk_TSAT
      else if (rval.lt.0.0d0) then
       if (axisymmetric_polar_flag.eq.0) then
        temperature=axisymmetric_disk_TSAT+axisymmetric_disk_TDIFF
       else
        print *,"axisymmetric_polar_flag or rval invalid"
        stop
       endif
      else if ((rval.ge.0.0d0).and.(rval.le.RR)) then
       call get_dr(dr,time_index)
       if (dr.gt.0.0d0) then
        ! (i+1/2)dr=rval
        i=NINT(rval/dr-0.5d0)
        if (i.eq.-1) then
         i=0
        else if (i.eq.ie+1) then
         i=ie
        else if ((i.ge.0).and.(i.le.ie)) then
         ! do nothing
        else
         print *,"i invalid"
         stop
        endif
        ri=(i+0.5d0)*dr

        if ((ri.gt.0.0d0).and.(ri.lt.RR)) then

         if (rval.le.ri) then
          if (i.eq.0) then
           if (axisymmetric_polar_flag.eq.1) then
            temperature=axisymmetric_disk_state(i,time_index)
           else if (axisymmetric_polar_flag.eq.0) then
            local_TWALL=axisymmetric_disk_TSAT+axisymmetric_disk_TDIFF
            temperature=(ri-rval)*local_TWALL/ri+ &
             rval*axisymmetric_disk_state(i,time_index)/ri
           else
            print *,"axisymmetric_polar_flag invalid"
            stop
           endif
          else if ((i.ge.1).and.(i.le.ie)) then
           TLEFT=axisymmetric_disk_state(i-1,time_index)
           TRIGHT=axisymmetric_disk_state(i,time_index)
           RLEFT=ri-dr
           RRIGHT=ri
           dr_interp=RRIGHT-RLEFT
           temperature=(RRIGHT-rval)*TLEFT/dr_interp + &
                       (rval-RLEFT)*TRIGHT/dr_interp
          else
           print *,"i invalid"
           stop
          endif
         else if (rval.ge.ri) then
          TLEFT=axisymmetric_disk_state(i,time_index)
          RLEFT=ri
          if (i.eq.ie) then
           TRIGHT=axisymmetric_disk_TSAT
           RRIGHT=RR 
          else if ((i.ge.0).and.(i.lt.ie)) then
           TRIGHT=axisymmetric_disk_state(i+1,time_index)
           RRIGHT=ri+dr
          else 
           print *,"i invalid"
           stop
          endif
          dr_interp=RRIGHT-RLEFT
          temperature=(RRIGHT-rval)*TLEFT/dr_interp + &
                      (rval-RLEFT)*TRIGHT/dr_interp
         else
          print *,"rval invalid"
          stop
         endif
        else
         print *,"ri invalid"
         stop
        endif
       else
        print *,"dr invalid"
        stop
       endif
      else
       print *,"rval invalid"
       stop
      endif
      if (temperature.ge.0.0d0) then
       ! do nothing
      else
       print *,"temperature invalid1"
       stop
      endif
 
      end subroutine disk_interp_T
 
      subroutine interp_T_to_grid(x,time_index,temperature)
      IMPLICIT NONE

      REAL(kind=8) :: x(2) 
      REAL(kind=8) :: temperature
      INTEGER :: time_index
      REAL(kind=8) :: dist

      if (axisymmetric_polar_flag.eq.1) then
       dist=sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)
       call disk_interp_T(dist,time_index,temperature) 
      else if (axisymmetric_polar_flag.eq.0) then
       dist=x(2)
       call disk_interp_T(dist,time_index,temperature) 
      else
       print *,"axisymmetric_polar_flag invalid"
       stop
      endif

      end subroutine interp_T_to_grid

      subroutine disk_get_speed(time_index,speed)
      IMPLICIT NONE

      REAL(kind=8) :: speed
      INTEGER :: time_index
      INTEGER :: ie
      REAL(kind=8) :: dr,TR

      ie=axisymmetric_disk_NCELL-1
      call get_dr(dr,time_index)
      if (dr.gt.0.0d0) then
       TR=(axisymmetric_disk_TSAT- &
           axisymmetric_disk_state(ie,time_index))/ &
          (0.5d0*dr)
       speed=abs(axisymmetric_disk_thermal_diffusivity*TR/ &
         axisymmetric_disk_L)
      else
       print *,"dr invalid"
       stop
      endif

      end subroutine disk_get_speed

       ! level set for material in the drop
      subroutine interp_LS_vel_to_grid(x,time_index,LS,vel)
      IMPLICIT NONE

      REAL(kind=8) :: x(2) 
      REAL(kind=8) :: vel(2) 
      REAL(kind=8) :: LS
      REAL(kind=8) :: RR,dist,speed
      INTEGER :: time_index

      RR=axisymmetric_disk_radblob(time_index)

      if (RR.gt.0.0d0) then
       vel(1)=0.0d0
       vel(2)=0.0d0
       if (axisymmetric_polar_flag.eq.1) then
        dist=sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)
        LS=RR-dist
        call disk_get_speed(time_index,speed)
        if (dist.gt.0.0d0) then
         vel(1)=speed*(x(1)-xblob)/dist
         vel(2)=speed*(x(2)-yblob)/dist
        else if (dist.eq.0.0d0) then
         ! do nothing
        else
         print *,"dist invalid"
         stop
        endif
       else if (axisymmetric_polar_flag.eq.0) then
        dist=x(2)
        LS=RR-dist
        call disk_get_speed(time_index,speed)
        vel(2)=speed
       else
        print *,"axisymmetric_polar_flag invalid"
        stop
       endif
      else
       print *,"RR invalid"
       stop
      endif

      return
      end subroutine interp_LS_vel_to_grid

      subroutine axisymmetric_disk_close()
      IMPLICIT NONE
 
      deallocate(axisymmetric_disk_state)

      return
      end subroutine axisymmetric_disk_close

       ! material 1 is inside, material 2 is outside
       ! if ireverse==0 then material 1 shrinks.
       ! if ireverse==1 then material 2 shrinks.
       ! TDIFF=T_DISK_CENTER - TSAT 
      subroutine axisymmetric_disk_init(L,TSAT,TDIFF, &
        thermal_diffusivity,RR, &
        stefan_flag,ireverse,isink,probtype,polar_flag)
      IMPLICIT NONE

      REAL(kind=8), intent(in) :: L,TSAT,TDIFF,thermal_diffusivity,RR
      INTEGER, intent(in) :: stefan_flag
      INTEGER, intent(in) :: ireverse
      INTEGER, intent(in) :: isink
      INTEGER, intent(in) :: probtype
      INTEGER, intent(in) :: polar_flag
      INTEGER      :: i,ie
      REAL(kind=8) :: dr,rval

      if ((L.ne.0.0d0).and. &
          (TSAT.ge.0.0d0).and. &
          (thermal_diffusivity.gt.0.0d0).and. &
          (RR.gt.0.0d0).and. &
          ((stefan_flag.eq.0).or.(stefan_flag.eq.1)).and. &
          ((ireverse.eq.0).or.(ireverse.eq.1)).and. &
          ((isink.eq.0).or.(isink.eq.1))) then
          ! do nothing
      else
       print *,"parameter bust axisymmetric_disk_init"
       stop
      endif
      axisymmetric_disk_isink=isink
      axisymmetric_disk_ireverse=ireverse
      axisymmetric_disk_L=L
      axisymmetric_disk_TSAT=TSAT
      axisymmetric_disk_TDIFF=TDIFF
      axisymmetric_disk_thermal_diffusivity=thermal_diffusivity
      axisymmetric_disk_radblob(1)=RR
      axisymmetric_disk_radblob(2)=RR
      axisymmetric_disk_stefan_flag=stefan_flag

      axisymmetric_polar_flag=polar_flag
      axisymmetric_probtype=probtype

      axisymmetric_disk_NCELL=1024   ! 1024
      axisymmetric_disk_NCYCLE=128   ! 128
      ie=axisymmetric_disk_NCELL-1
      allocate(axisymmetric_disk_state(0:ie,2))
      axisymmetric_disk_time=0.0d0

      call get_dr(dr,1)
      if (dr.gt.0.0d0) then
       do i=0,ie
        rval=(i+0.5d0)*dr
        call disk_eval_initial_temp(rval,axisymmetric_disk_state(i,1))
        call disk_eval_initial_temp(rval,axisymmetric_disk_state(i,2))
        if (1.eq.0) then
         print *,"r,tempinit ",rval,axisymmetric_disk_state(i,2)
        endif
       enddo
      else
       print *,"dr invalid"
       stop
      endif

      return
      end subroutine axisymmetric_disk_init
      
      subroutine axisymmetric_disk_advance(dt)
      IMPLICIT NONE

      REAL(kind=8) :: dt
      REAL(kind=8) :: dtsub
      REAL(kind=8) :: dt_i
      REAL(kind=8) :: interp_temp_i
      REAL(kind=8), dimension(:), allocatable :: L_array
      REAL(kind=8), dimension(:), allocatable :: U_array
      REAL(kind=8), dimension(:), allocatable :: D_array
      REAL(kind=8), dimension(:), allocatable :: RHS
      REAL(kind=8), dimension(:), allocatable :: SOLN
      INTEGER :: i,ie,nsub
      REAL(kind=8) :: speed,RR_old,RR_new,dr_new
      REAL(kind=8) :: ri_new,rplus,rminus,coef

      if ((dt.gt.0.0d0).and. &
          (axisymmetric_disk_NCYCLE.ge.1)) then
       ! do nothing
      else
       print *,"dt or axisymmetric_disk_NCYCLE invalid"
       stop
      endif
      ie=axisymmetric_disk_NCELL-1

      allocate(L_array(0:ie))
      allocate(U_array(0:ie))
      allocate(D_array(0:ie))
      allocate(RHS(0:ie))
      allocate(SOLN(0:ie))

      dtsub=dt/axisymmetric_disk_NCYCLE
      do nsub=1,axisymmetric_disk_NCYCLE
       call disk_get_speed(1,speed)
       RR_old=axisymmetric_disk_radblob(1)
       RR_new=RR_old
       if (axisymmetric_disk_stefan_flag.eq.1) then
        if (axisymmetric_disk_ireverse.eq.0) then !inside material shrinks
         RR_new=RR_new-dtsub*abs(speed)
        else if (axisymmetric_disk_ireverse.eq.1) then !inside material grows
         RR_new=RR_new+dtsub*abs(speed)
        else
         print *,"axisymmetric_disk_ireverse invalid"
         stop
        endif
       else if (axisymmetric_disk_stefan_flag.eq.0) then
        ! do nothing
       else
        print *,"axisymmetric_disk_stefan_flag invalid"
        stop
       endif
       axisymmetric_disk_radblob(2)=RR_new
       call get_dr(dr_new,2)

       if (dr_new.gt.0.0d0) then
        do i=0,ie
         ri_new=(i+0.5d0)*dr_new

         if (1.eq.0) then

          if ((ri_new.gt.0.0d0).and.(ri_new.le.RR_old)) then
           dt_i=dtsub
           call disk_interp_T(ri_new,1,interp_temp_i)
          else if ((ri_new.ge.RR_old).and.(ri_new.lt.RR_new)) then
           dt_i=dtsub*(RR_new-ri_new)/(RR_new-RR_old) 
           interp_temp_i=axisymmetric_disk_TSAT
          else
           print *,"ri_new invalid"
           stop
          endif

         else if (1.eq.1) then
          dt_i=dtsub
          if ((i.eq.ie-2).or.(i.eq.ie-1).or.(i.eq.ie)) then
           interp_temp_i=axisymmetric_disk_state(i,1)
          else
           call disk_interp_T(ri_new,1,interp_temp_i)
          endif
         endif
         rplus=ri_new+0.5d0*dr_new
         rminus=ri_new-0.5d0*dr_new

         if (axisymmetric_polar_flag.eq.1) then

          coef=axisymmetric_disk_thermal_diffusivity/ &
           (ri_new*(dr_new**2))

          if (i.eq.ie) then
           rplus=2.0d0*rplus
          else if ((i.ge.0).and.(i.lt.ie)) then 
           ! do nothing
          else
           print *,"i invalid"
           stop
          endif

          if (i.eq.0) then
           D_array(i)=1.0d0/dt_i+coef*rplus
           L_array(i)=0.0d0
          else if ((i.ge.1).and.(i.le.ie)) then 
           D_array(i)=1.0d0/dt_i+coef*(rplus+rminus)
           L_array(i)=-coef*rminus
          else
           print *,"i invalid"
           stop
          endif
          if (i.eq.ie) then
           U_array(i)=0.0d0
           RHS(i)=interp_temp_i/dt_i+coef*rplus*axisymmetric_disk_TSAT
          else if ((i.ge.0).and.(i.lt.ie)) then 
           U_array(i)=-coef*rplus
           RHS(i)=interp_temp_i/dt_i
          else
           print *,"i invalid"
           stop
          endif

         else if (axisymmetric_polar_flag.eq.0) then

          coef=axisymmetric_disk_thermal_diffusivity/(dr_new**2)

            ! u_{i-1} - 2 u_{i} + u_{i+1} 
            ! u_{i}-u_{w} = u_{w}-u_{i-1}
            ! u_{i-1}=2 u_{w} - u_{i}
          if (i.eq.0) then
           D_array(i)=1.0d0/dt_i+3.0d0*coef
           L_array(i)=0.0d0
           U_array(i)=-coef
           RHS(i)=interp_temp_i/dt_i+ &
             2.0d0*coef*(axisymmetric_disk_TSAT+axisymmetric_disk_TDIFF)
          else if ((i.ge.1).and.(i.lt.ie)) then 
           D_array(i)=1.0d0/dt_i+2.0d0*coef
           L_array(i)=-coef
           U_array(i)=-coef
           RHS(i)=interp_temp_i/dt_i
          else if (i.eq.ie) then
           D_array(i)=1.0d0/dt_i+3.0d0*coef
           U_array(i)=0.0d0
           L_array(i)=-coef
           RHS(i)=interp_temp_i/dt_i+2.0d0*coef*axisymmetric_disk_TSAT
          else
           print *,"i invalid"
           stop
          endif
         else 
          print *,"axisymmetric_polar_flag invalid"
          stop
         endif

         if (1.eq.0) then
          print *,"i,interp_temp,temp ",i,interp_temp_i, &
                   axisymmetric_disk_state(i,1)
         endif
        enddo ! i=0..ie
       else
        print *,"dr_new invalid"
        stop
       endif

       call tridiag_solve(L_array,U_array,D_array,ie+1,RHS,SOLN)
       if (axisymmetric_disk_isink.eq.0) then
        ! do nothing
       else if (axisymmetric_disk_isink.eq.1) then
        SOLN(0)=axisymmetric_disk_TSAT+axisymmetric_disk_TDIFF
       else
        print *,"axisymmetric_disk_isink invalid"
        stop
       endif
       axisymmetric_disk_radblob(1)=RR_new
       axisymmetric_disk_time=axisymmetric_disk_time+dtsub
       do i=0,ie
        axisymmetric_disk_state(i,1)=SOLN(i)
        axisymmetric_disk_state(i,2)=SOLN(i)
        ri_new=(i+0.5d0)*dr_new
        if (1.eq.0) then
         print *,"time,r,tempnew ",axisymmetric_disk_time, &
           ri_new,axisymmetric_disk_state(i,2)
        endif
       enddo
      enddo ! nsub=1..axisymmetric_disk_NCYCLE

      deallocate(L_array)
      deallocate(U_array)
      deallocate(D_array)
      deallocate(RHS)
      deallocate(SOLN)

      end subroutine axisymmetric_disk_advance

      end module variable_temperature_drop
