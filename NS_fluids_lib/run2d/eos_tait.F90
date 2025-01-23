      subroutine rhotait(p,rho)
      IMPLICIT NONE

      real*8 rho,p
      real*8 A_TAIT,B_TAIT,RHOBAR_TAIT,GAMMA_TAIT

      A_TAIT=1.0D+6
      B_TAIT=3.31D+9
      RHOBAR_TAIT=1.0D0
      GAMMA_TAIT=7.15D0

      rho=(p-A_TAIT)/B_TAIT+1.0D0
      rho=(rho**(1.0D0/GAMMA_TAIT))*RHOBAR_TAIT
      return
      end

      subroutine ptait(p,rho)
      IMPLICIT NONE

      real*8 rho,p
      real*8 A_TAIT,B_TAIT,RHOBAR_TAIT,GAMMA_TAIT

      A_TAIT=1.0D+6
      B_TAIT=3.31D+9
      RHOBAR_TAIT=1.0D0
      GAMMA_TAIT=7.15D0

      p=B_TAIT*( (rho/RHOBAR_TAIT)**GAMMA_TAIT - 1.0D0 ) + A_TAIT

      return
      end

      subroutine soundtait(p,rho,sound)
      IMPLICIT NONE

      real*8 rho,p,sound
      real*8 A_TAIT,B_TAIT,RHOBAR_TAIT,GAMMA_TAIT

      A_TAIT=1.0D+6
      B_TAIT=3.31D+9
      RHOBAR_TAIT=1.0D0
      GAMMA_TAIT=7.15D0

      call ptait(p,rho)
      sound=B_TAIT*GAMMA_TAIT*(  (rho/RHOBAR_TAIT)**(GAMMA_TAIT-1.0D0) )/ &
        RHOBAR_TAIT
      sound=sqrt(sound) 

      return
      end
 
      program main
      IMPLICIT NONE
      real*8 minimum_pressure,minimum_density
      real*8 minimum_sound

      minimum_pressure=220.2726d0
      call rhotait(minimum_pressure,minimum_density)
      print *,"minimum_pressure= ",minimum_pressure
      print *,"minimum_density= ",minimum_density
      call ptait(minimum_pressure,minimum_density)
      print *,"SANITY CHECK: "
      print *,"minimum_pressure= ",minimum_pressure
      print *,"minimum_density= ",minimum_density
 
      call soundtait(minimum_pressure,minimum_density,minimum_sound)
      print *,"minimum_sound= ",minimum_sound

      return
      end

