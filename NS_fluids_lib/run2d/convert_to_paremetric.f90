      PROGRAM tecwhale


      IMPLICIT NONE

      Integer :: N,N2, i
      Real :: t1,t2,cenx,ceny,scale

      print *,"This program reads cenx and ceny"
      print *,"and outputs a gnuplot file: cenxy"
      open(unit=2, file= 'cenx')
      open(unit=4, file= 'ceny')
      open(unit=14, file= 'cenxy')

      print *,"reading number of time entries"
      read(2,*) N
      read(4,*) N2
      if (N.eq.N2) then
       print *,"number of entries: ",N
      else
       print *,"N <> N2"
       stop
      endif

      scale=0.585
      do i=1,N
       read(2,*) t1,cenx
       read(4,*) t2,ceny
       write(14,*) cenx/scale,ceny/scale
      enddo
      close(2)
      close(4)
      close(14)

      END PROGRAM
