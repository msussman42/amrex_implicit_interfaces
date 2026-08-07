      program main
      IMPLICIT NONE

      real*8 delta_x,speed,artificial
      print *,"Kurganov Tadmor artificial viscosity is "
      print *,"the difference between central and upwind differencing"
      print *,"which is speed * \Delta x/2"
      delta_x=0.01
      speed=1.0E+5
      print *,"delta_x= ",delta_x
      print *,"speed=",speed
      artificial=speed*delta_x/2.0
      print *,"speed * \Delta x/2=artificial= ",artificial
      end
