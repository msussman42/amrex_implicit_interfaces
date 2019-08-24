      program main
      IMPLICIT NONE
 
      real*8 F1,F3,M1,M3
      real*8 F1end,F3end,M1end,M3end

      F1=0.719170084025
      F3=0.151385621791
      M1=0.719170084025
      M3=0.141394170753

      F1end=0.0731710581883
      F3end=0.852891760181
      M1end=0.0731710581851
      M3end=0.796600904009

      print *,"mass init= ",M1+M3
      print *,"volume init= ",F1+F3
      print *,"mass init end= ",M1end+M3end
      print *,"volume init end= ",F1end+F3end
      print *,"difference in volume: ",F1end+F3end-(F1+F3)
      print *,"difference in mass: ",M1end+M3end-(M1+M3)
      end

