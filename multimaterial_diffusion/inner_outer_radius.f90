      subroutine eff_rad(V1,V2,R1,R2,R1C,R2C)
      IMPLICIT NONE

      real*8 V1,V2,R1,R2,R1C,R2C
      real*8 my_pi
      
      my_pi=4.0d0*atan(1.0d0)
       ! pi r^2 =v
       ! r=sqrt(v/pi)
      R1=sqrt(V1/my_pi)
      R2=sqrt((V1+V2)/my_pi)

      print *,"V1,V2,R1,R2 ",V1,V2,R1,R2
      print *,"R1-R1C ",R1-R1C
      print *,"R2-R2C ",R2-R2C
      R1C=R1
      R2C=R2
      return
      end subroutine

      PROGRAM inner_outer

      IMPLICIT NONE

      real*8 V1,V2
      real*8 R1,R2
      real*8 R1C,R2C

       ! 32
      R1C=0.0d0
      R2C=0.0d0
      V1=0.143975d0
      V2=0.183520d0
      call eff_rad(V1,V2,R1,R2,R1C,R2C)
       ! 64
      V1=0.144011d0
      V2=0.196533d0
      call eff_rad(V1,V2,R1,R2,R1C,R2C)
       ! 128
      V1=0.144036d0
      V2=0.205831d0
      call eff_rad(V1,V2,R1,R2,R1C,R2C)
       ! 256
      V1=0.144030d0
      V2=0.211203d0
      call eff_rad(V1,V2,R1,R2,R1C,R2C)
 
      END PROGRAM
