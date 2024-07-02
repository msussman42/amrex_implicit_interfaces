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

        !see run2d/NormalShockWaveNASA.F90
        !see run2d/inputs.shockdrop

       real(amrex_real) shockdrop_R
       real(amrex_real) shockdrop_cp
       real(amrex_real) shockdrop_cv

        !upstream supersonic relative to shock
       real(amrex_real) shockdrop_P0
       real(amrex_real) shockdrop_T0
       real(amrex_real) shockdrop_DEN0
       real(amrex_real) shockdrop_M0
       real(amrex_real) shockdrop_VEL0
       real(amrex_real) shockdrop_C0
       real(amrex_real) shockdrop_EE0

        !downstream subsonic relative to shock
       real(amrex_real) shockdrop_P1
       real(amrex_real) shockdrop_T1
       real(amrex_real) shockdrop_DEN1
       real(amrex_real) shockdrop_M1
       real(amrex_real) shockdrop_VEL1
       real(amrex_real) shockdrop_C1
       real(amrex_real) shockdrop_EE1

       real(amrex_real) shockdrop_gamma

CONTAINS

          ! units are cgs
       subroutine shockdrop_init()
       use probcommon_module
       IMPLICIT NONE

        !material_type=5 EOS_air
        !see PROBCOMMON.F90
       shockdrop_R=R_AIR_PARMS !ergs/(Kelvin g)
       shockdrop_cv=CV_AIR_PARMS !ergs/(Kelvin g)
       shockdrop_cp=shockdrop_cv+shockdrop_R !ergs/(Kelvin g)

!      shockdrop_M0=1.4
!      shockdrop_M0=3.0
       shockdrop_M0=1.17
       shockdrop_T0=278.0d0 !tempconst(2)
       shockdrop_DEN0=0.00125335272d0 !denconst(2)

       if (abs(shockdrop_DEN0-fort_denconst(2)).le.EPS6) then
        !do nothing
       else
        print *,"shockdrop_DEN0 ",shockdrop_DEN0
        print *,"fort_denconst(2) ",fort_denconst(2)
        print *,"mismatch"
        stop
       endif

       if (abs(shockdrop_T0-fort_tempconst(2)).le.EPS6) then
        !do nothing
       else
        print *,"shockdrop_T0 ",shockdrop_T0
        print *,"fort_tempconst(2) ",fort_tempconst(2)
        print *,"mismatch"
        stop
       endif
       if (abs(shockdrop_T0-fort_tempconst(1)).le.EPS6) then
        !do nothing
       else
        print *,"shockdrop_T0 ",shockdrop_T0
        print *,"fort_tempconst(1) ",fort_tempconst(1)
        print *,"mismatch"
        stop
       endif
    

        ! this value must be consistent with material_type=5 parameters
       shockdrop_gamma=shockdrop_cp/shockdrop_cv

        ! if compressible liquid, this value should be the same
        ! as the equilibrium pressure of the liquid drop.
       shockdrop_P0=(shockdrop_gamma-1.0d0)* &
        shockdrop_DEN0*shockdrop_cv*shockdrop_T0

       shockdrop_C0=(shockdrop_gamma*shockdrop_P0/shockdrop_DEN0)**half
       shockdrop_VEL0=shockdrop_M0*shockdrop_C0

       shockdrop_P1=shockdrop_P0* &
        (two*shockdrop_gamma*(shockdrop_M0**2)-shockdrop_gamma+one)/ &
        (shockdrop_gamma+one)

       shockdrop_T1=shockdrop_T0* &
        ( two*shockdrop_gamma*(shockdrop_M0**2)-shockdrop_gamma+one )* &
        ( (shockdrop_gamma-one)*(shockdrop_M0**2)+two )/ &
        ( ((shockdrop_gamma+one)*shockdrop_M0)**2 )

       shockdrop_DEN1=shockdrop_DEN0* &
        ((shockdrop_gamma+one)*(shockdrop_M0**2))/ &
        ((shockdrop_gamma-one)*(shockdrop_M0**2)+two)

       shockdrop_C1=sqrt(shockdrop_gamma*shockdrop_P1/shockdrop_DEN1)
       shockdrop_M1= &
        ((shockdrop_gamma-one)*(shockdrop_M0**2)+two)/ &
        (two*shockdrop_gamma*(shockdrop_M0**2)-shockdrop_gamma+one)
       shockdrop_M1=sqrt(shockdrop_M1)
       shockdrop_VEL1=shockdrop_M1*shockdrop_C1

       print *,"shockdrop: upstream den,approaching SPEED,T,M,C ", &
        shockdrop_DEN0,shockdrop_VEL0,shockdrop_T0,shockdrop_M0,shockdrop_C0
       print *,"shockdrop: downstream den,SPEED,T,M,C ", &
        shockdrop_DEN1,shockdrop_VEL1,shockdrop_T1,shockdrop_M1,shockdrop_C1

       return
       end subroutine shockdrop_init

       subroutine shockdrop_velocity(x,y,z,vel, &
         xblob,yblob,zblob,radblob,zblob2,axis_dir)
       IMPLICIT NONE
       integer axis_dir,dir
       real(amrex_real) x,y,z,xblob,yblob,zblob,radblob,zblob2,LS
       real(amrex_real) vel(SDIM)

       if (SDIM.eq.2) then
        if (abs(z-y).le.1.0E-6) then
         !do nothing
        else
         print *,"expecting z=y"
         stop
        endif
       endif

       do dir=1,SDIM
        vel(dir)=zero
       enddo

       call shockdrop_dropLS(x,y,z,LS,xblob,yblob,zblob,radblob,axis_dir)
       if (LS.ge.zero) then
        ! do nothing (drop is upstream from shock and
        ! stationary in the "upstream frame of reference")
        ! shock velocity > 0
        ! shock is approaching with speed: shockdrop_VEL0
       else
        call shockdrop_shockLS(x,y,z,LS,zblob2,axis_dir)
         ! in shock frame of reference:
         ! upstream: v=-shockdrop_VEL0
         ! downstream: v=-shockdrop_VEL1
        if (LS.ge.zero) then  ! upstream (above the shock)
         vel(SDIM)=zero
        else
         vel(SDIM)=shockdrop_VEL0-shockdrop_VEL1
        endif
       endif 

       return
       end subroutine shockdrop_velocity


       subroutine shockdrop_maxvelocity(vel,axis_dir)
       IMPLICIT NONE
       integer axis_dir
       real(amrex_real) vel

       vel=max(abs(shockdrop_VEL0),abs(shockdrop_VEL1))

       return
       end subroutine shockdrop_maxvelocity



       subroutine shockdrop_pressure(x,y,z,pres, &
         xblob,yblob,zblob,radblob,zblob2,axis_dir)
       IMPLICIT NONE
       integer axis_dir
       real(amrex_real) x,y,z,xblob,yblob,zblob,radblob,zblob2
       real(amrex_real) pres,LS

       if (SDIM.eq.2) then
        if (abs(z-y).le.1.0E-6) then
         !do nothing
        else
         print *,"expecting z=y"
         stop
        endif
       endif

       call shockdrop_dropLS(x,y,z,LS,xblob,yblob,zblob,radblob,axis_dir)
       if (LS.ge.zero) then ! liquid
        pres=shockdrop_P0
       else
        call shockdrop_shockLS(x,y,z,LS,zblob2,axis_dir)
        if (LS.ge.zero) then  ! upstream (above the approaching shock)
         pres=shockdrop_P0
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
        if (abs(z-y).le.1.0E-6) then
         !do nothing
        else
         print *,"expecting z=y"
         stop
        endif
       endif

       call shockdrop_dropLS(x,y,z,LS,xblob,yblob,zblob,radblob,axis_dir)
       if (LS.ge.zero) then ! liquid
        den=shockdrop_DEN0
       else
        call shockdrop_shockLS(x,y,z,LS,zblob2,axis_dir)
        if (LS.ge.zero) then  ! upstream (above the shock)
         den=shockdrop_DEN0
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
        if (abs(z-y).le.1.0E-6) then
         !do nothing
        else
         print *,"expecting z=y"
         stop
        endif
       endif

       call shockdrop_dropLS(x,y,z,LS,xblob,yblob,zblob,radblob,axis_dir)
       if (LS.ge.zero) then ! liquid
        temp=shockdrop_T0
       else
        call shockdrop_shockLS(x,y,z,LS,zblob2,axis_dir)
        if (LS.ge.zero) then  ! upstream (above the shock)
         temp=shockdrop_T0
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
        if (abs(z-y).le.1.0E-6) then
         !do nothing
        else
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
        if (abs(z-y).le.1.0E-6) then
         !do nothing
        else
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
