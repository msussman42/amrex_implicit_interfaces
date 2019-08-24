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

      real*8 diameter,pv,p,sigma,U,rho,rhowater
      real*8 RE,gravity,depth,rho_at_depth,visc
      real*8 dimensionless_slope1
      real*8 dimensionless_slope2
      real*8 dimensionless_slope3
      real*8 USTAR,UPART1,UPART2,TB,RS,T_DELAY,TPART1
      real*8 minimum_density,minimum_pressure,minimum_sound

      minimum_pressure=1.0D0
      call rhotait(minimum_pressure,minimum_density)
      print *,"minimum_pressure= ",minimum_pressure
      print *,"minimum_density= ",minimum_density
      call ptait(minimum_pressure,minimum_density)
      print *,"SANITY CHECK: "
      print *,"minimum_pressure= ",minimum_pressure
      print *,"minimum_density= ",minimum_density
 
      call soundtait(minimum_pressure,minimum_density,minimum_sound)
      print *,"minimum_sound= ",minimum_sound


      pv=1.0e+6
      rhowater=1.0
      if (1.eq.0) then
       diameter=1.5
       rho=1.41 ! Ertacetal
       sigma=11.7
       U=4.4D+2
      else if (1.eq.0) then
       diameter=1.5
       rho=7.8 ! stainless steel
       sigma=20.3
       U=3.4D+2
      else if (1.eq.0) then
       diameter=4.5
       rho=1.41
       sigma=35.37
       U=2.4D+2
      else if (1.eq.1) then ! figure 4, left (stainless steel)
       diameter=4.5
       RS=0.5*diameter
       rho=7.8
       sigma=37.9
       U=2.3D+2
           ! 0.1=113-52  blue=113-74
       dimensionless_slope1=(0.1*(113.0-74.0)/(113.0-52.0))/0.5
           ! blue=(0.1*(113-93)/(113-52))/(0.5*(242-183)/(220-130)) 
       dimensionless_slope2=(0.1*(113.0-93.0)/(113.0-52.0))/ &
                            (0.5*(242.0-183.0)/(220.0-130.0))
           ! blue=(0.1*(93-72)/(113-52))/(0.5*(569-242)/(220-130)) 
       dimensionless_slope3=(0.1*(93.0-72.0)/(113.0-52.0))/ &
                            (0.5*(569.0-242.0)/(220.0-130.0))
       USTAR=U/dimensionless_slope1
       UPART1=dimensionless_slope2*USTAR
       UPART2=dimensionless_slope3*USTAR
       TB=RS/USTAR
       T_DELAY=TB*(184.0-131.0)/(307.0-131.0)
       TPART1=TB*(0.5*(242.0-183.0)/(220.0-130.0))
      else if (1.eq.0) then ! figure 2 (stainless steel)
       diameter=4.5
       rho=7.8
       sigma=63.74
       U=1.8D+2
      endif
      visc=0.01
      RE=0.5*rhowater*diameter*U/visc
      p=0.5*rhowater*U*U*sigma+pv
      gravity=981.0
      depth=(p-pv)/(rhowater*gravity)
      call rhotait(p,rho_at_depth)

      print *,"rhowater=",rhowater
      print *,"diameter,rhosolid,sigma,U ",diameter,rho,sigma,U
      print *,"p,pv,gravity,depth ",p,pv,gravity,depth
      print *,"RE=",RE
      print *,"p-pv ",p-pv
      print *,"rho_at_depth ",rho_at_depth
      print *,"USTAR ",USTAR
      print *,"T_DELAY,UPART1,UPART2 ",T_DELAY,UPART1,UPART2
      print *,"TPART1 ",TPART1
      print *,"TB ",TB
      
      return
      end

