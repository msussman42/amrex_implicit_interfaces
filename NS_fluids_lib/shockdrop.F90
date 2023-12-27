#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

MODULE shockdrop
use amrex_fort_module, only : amrex_real

       real(amrex_real) shockdrop_P
       real(amrex_real) shockdrop_T
       real(amrex_real) shockdrop_DEN
       real(amrex_real) shockdrop_M
       real(amrex_real) shockdrop_VEL
       real(amrex_real) shockdrop_C

       real(amrex_real) shockdrop_P1
       real(amrex_real) shockdrop_T1
       real(amrex_real) shockdrop_DEN1
       real(amrex_real) shockdrop_M1
       real(amrex_real) shockdrop_VEL1
       real(amrex_real) shockdrop_C1

       real(amrex_real) shockdrop_gamma

CONTAINS

          ! units are cgs
       subroutine shockdrop_init()
       IMPLICIT NONE

!      shockdrop_M=1.4
       shockdrop_M=3.0

         ! if this value is not 1.39861, then material_type=5 parameters
         ! should be changed too:
       shockdrop_gamma=1.39861111111111d0
        ! if compressible liquid, this value should be the same
        ! as the equilibrium pressure of the liquid drop.
       shockdrop_P=1.0d+6 ! 1 g cm/s^2 /cm^2=g/(s^2 cm)=(1/10) kg/(s^2 m)
       shockdrop_T=278.0
       shockdrop_DEN=0.00123
       shockdrop_C=(shockdrop_gamma*shockdrop_P/shockdrop_DEN)**half
       shockdrop_VEL=shockdrop_M*shockdrop_C

       shockdrop_P1=shockdrop_P* &
        (two*shockdrop_gamma*(shockdrop_M**2)-shockdrop_gamma+one)/ &
        (shockdrop_gamma+one)

       shockdrop_T1=shockdrop_T* &
        ( two*shockdrop_gamma*(shockdrop_M**2)-shockdrop_gamma+one )* &
        ( (shockdrop_gamma-one)*(shockdrop_M**2)+two )/ &
        ( ((shockdrop_gamma+one)*shockdrop_M)**2 )

       shockdrop_DEN1=shockdrop_DEN* &
        ((shockdrop_gamma+one)*(shockdrop_M**2))/ &
        ((shockdrop_gamma-one)*(shockdrop_M**2)+two)

       shockdrop_C1=sqrt(shockdrop_gamma*shockdrop_P1/shockdrop_DEN1)
       shockdrop_M1= &
        ((shockdrop_gamma-one)*(shockdrop_M**2)+two)/ &
        (two*shockdrop_gamma*(shockdrop_M**2)-shockdrop_gamma+one)
       shockdrop_M1=sqrt(shockdrop_M1)
       shockdrop_VEL1=shockdrop_M1*shockdrop_C1

       print *,"shockdrop: upstream den,vel,T,M,C ", &
        shockdrop_DEN,shockdrop_VEL,shockdrop_T,shockdrop_M,shockdrop_C
       print *,"shockdrop: downstream den,vel,T,M,C ", &
        shockdrop_DEN1,shockdrop_VEL1,shockdrop_T1,shockdrop_M1,shockdrop_C

       return
       end subroutine shockdrop_init

       subroutine shockdrop_velocity(x,y,z,vel, &
         xblob,yblob,zblob,radblob,zblob2,axis_dir)
       IMPLICIT NONE
       integer axis_dir,dir
       real(amrex_real) x,y,z,xblob,yblob,zblob,radblob,zblob2,LS
       real(amrex_real) vel(SDIM)

       if (SDIM.eq.2) then
        if (abs(z-y).gt.1.0E-6) then
         print *,"expecting z=y"
         stop
        endif
       endif
       do dir=1,SDIM
        vel(dir)=zero
       enddo

       call shockdrop_dropLS(x,y,z,LS,xblob,yblob,zblob,radblob,axis_dir)
       if (LS.ge.zero) then
        ! do nothing
       else
        call shockdrop_shockLS(x,y,z,LS,zblob2,axis_dir)
        if (LS.ge.zero) then  ! upstream
         vel(SDIM)=zero
        else
         vel(SDIM)=shockdrop_VEL-shockdrop_VEL1
        endif
       endif 

       return
       end subroutine shockdrop_velocity


       subroutine shockdrop_maxvelocity(vel,axis_dir)
       IMPLICIT NONE
       integer axis_dir
       real(amrex_real) vel

       vel=max(abs(shockdrop_VEL),abs(shockdrop_VEL1))

       return
       end subroutine shockdrop_maxvelocity



       subroutine shockdrop_pressure(x,y,z,pres, &
         xblob,yblob,zblob,radblob,zblob2,axis_dir)
       IMPLICIT NONE
       integer axis_dir
       real(amrex_real) x,y,z,xblob,yblob,zblob,radblob,zblob2
       real(amrex_real) pres,LS

       if (SDIM.eq.2) then
        if (abs(z-y).gt.1.0E-6) then
         print *,"expecting z=y"
         stop
        endif
       endif

       call shockdrop_dropLS(x,y,z,LS,xblob,yblob,zblob,radblob,axis_dir)
       if (LS.ge.zero) then ! liquid
        pres=shockdrop_P
       else
        call shockdrop_shockLS(x,y,z,LS,zblob2,axis_dir)
        if (LS.ge.zero) then  ! upstream
         pres=shockdrop_P
        else
         pres=shockdrop_P1
        endif
       endif 

       return
       end subroutine shockdrop_pressure


       subroutine shockdrop_gas_density(x,y,z,den, &
         xblob,yblob,zblob,radblob,zblob2,axis_dir)
       IMPLICIT NONE
       integer axis_dir
       real(amrex_real) x,y,z,xblob,yblob,zblob,radblob,zblob2
       real(amrex_real) den,LS

       if (SDIM.eq.2) then
        if (abs(z-y).gt.1.0E-6) then
         print *,"expecting z=y"
         stop
        endif
       endif

       call shockdrop_dropLS(x,y,z,LS,xblob,yblob,zblob,radblob,axis_dir)
       if (LS.ge.zero) then ! liquid
        den=shockdrop_DEN
       else
        call shockdrop_shockLS(x,y,z,LS,zblob2,axis_dir)
        if (LS.ge.zero) then  ! upstream
         den=shockdrop_DEN
        else
         den=shockdrop_DEN1
        endif
       endif 

       return
       end subroutine shockdrop_gas_density


       subroutine shockdrop_gas_temperature(x,y,z,temp, &
         xblob,yblob,zblob,radblob,zblob2,axis_dir)
       IMPLICIT NONE
       integer axis_dir
       real(amrex_real) x,y,z,xblob,yblob,zblob,radblob,zblob2
       real(amrex_real) temp,LS

       if (SDIM.eq.2) then
        if (abs(z-y).gt.1.0E-6) then
         print *,"expecting z=y"
         stop
        endif
       endif

       call shockdrop_dropLS(x,y,z,LS,xblob,yblob,zblob,radblob,axis_dir)
       if (LS.ge.zero) then ! liquid
        temp=shockdrop_T
       else
        call shockdrop_shockLS(x,y,z,LS,zblob2,axis_dir)
        if (LS.ge.zero) then  ! upstream
         temp=shockdrop_T
        else
         temp=shockdrop_T1
        endif
       endif 

       return
       end subroutine shockdrop_gas_temperature


       ! probtype=1 in the inputs file
       ! axis_dir=150 shock drop
       ! axis_dir=151 shock column
       ! LS>0 upstream of the shock  z>zblob2
       ! LS<0 downstream of the shock z<zblob2
       SUBROUTINE shockdrop_shockLS(x,y,z,LS,zblob2,axis_dir)
       IMPLICIT NONE
       integer,INTENT(in) :: axis_dir
       real(amrex_real),INTENT(in) :: x,y,z,zblob2
       real(amrex_real),INTENT(out) :: LS

       if (SDIM.eq.2) then
        if (abs(z-y).gt.1.0E-6) then
         print *,"z<>y error"
         stop
        endif
       endif

       if ((axis_dir.eq.150).or.(axis_dir.eq.151)) then
        LS=z-zblob2
       else
        print *,"axis_dir invalid"
        stop
       endif

       return
       END SUBROUTINE shockdrop_shockLS

       ! probtype=1 in the inputs file
       ! axis_dir=150 shock drop
       ! axis_dir=151 shock column
       ! LS>0 in the drop
       subroutine shockdrop_dropLS(x,y,z,LS, &
         xblob,yblob,zblob,radblob,axis_dir)
       IMPLICIT NONE
       integer axis_dir
       real(amrex_real) x,y,z,LS,xblob,yblob,zblob,radblob,mag

       if (SDIM.eq.2) then
        if (abs(z-y).gt.1.0E-6) then
         print *,"z<>y error"
         stop
        endif
       endif

       if (axis_dir.eq.150) then ! shock drop
        mag=(x-xblob)**2+(y-yblob)**2
        if (SDIM.eq.3) then
         mag=mag+(z-zblob)**2
        endif
        LS=radblob-sqrt(mag) 
       else if (axis_dir.eq.151) then ! shock column
        if (SDIM.eq.2) then
         mag=(y-yblob)**2
         LS=radblob-sqrt(mag) 
        else if (SDIM.eq.3) then
         mag=(z-zblob)**2
         LS=radblob-sqrt(mag) 
        else 
         print *,"dimension bust"
         stop
        endif
       else
        print *,"axis_dir invalid in shockdrop_dropLS"
        stop
       endif
 
       return
       END SUBROUTINE shockdrop_dropLS



END MODULE shockdrop
