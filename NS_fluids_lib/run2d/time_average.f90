      PROGRAM time_average

      IMPLICIT NONE

      integer, PARAMETER :: TIMESTEPSMAX=60000
      integer, PARAMETER :: MODEMAX=200
      Real*8, PARAMETER :: tstart=0.012
      Real*8, PARAMETER :: tend=0.022
      Real*8, PARAMETER :: slice_position=1.5
      Real*8, PARAMETER :: PI=3.141592654
      Real*8 :: xstart,xend
      Real*8 :: xdata(0:TIMESTEPSMAX)  ! array for x variable
      Real*8 :: tt_average(0:TIMESTEPSMAX)  ! array for x variable
      Real*8 :: tt_max(0:TIMESTEPSMAX)  ! array for x variable
      Real*8 :: tdata(0:TIMESTEPSMAX)  ! array for t variable
      Real*8 :: xx_average(0:TIMESTEPSMAX) ! array for t variable
      Real*8 :: slice(0:TIMESTEPSMAX,2) ! array for t variable
      character*256 garbage
      integer :: ispace,itime,itime_slice
      Real*8 :: prev_time,current_time
      Real*8 :: old_space_average
      Real*8 :: prev_space,current_space,interface_height
      integer :: nspace,ntime
      Real*8 :: time_weight
      Real*8 :: space_weight
      Real*8 :: probhix
      integer :: line_marker
      integer :: ncomp,ncomp_experiment
      Real*8 xraster,yraster,xactual,yactual
      Real*8 L,T1,T2,TH,H1,H2,HH,DT,det
      integer imode,imodemax,j,k
      Real*8 :: A(2,2)
      Real*8 :: B(2)
      Real*8 :: X(2)
      Real*8 :: AN(MODEMAX)
      Real*8 :: BN(MODEMAX)
      Real*8 :: CN(MODEMAX)

      xstart=0.5
      xend=1.5
      probhix=2.0
      ncomp_experiment=12

      print *,"This program finds time average profile vs x"
      print *,"This program also finds time max profile vs x"
      print *,"This program outputs experimental data"
      print *,"probhix= ",probhix
      print *,"xstart=",xstart
      print *,"xend=",xend
      print *,"This program also finds space average height vs t"
      print *,"xstart=",xstart
      print *,"xend=",xend
      print *,"tstart=",tstart
      print *,"tend=",tend
      print *,"This program also finds slice height vs t"
      print *,"slice_position=",slice_position

      print *,"input files: height_integral.txt, jet_height.txt"
      print *,"input files: experiment.txt"
      print *,"output: space_average.txt, time_average.txt"
      print *,"time_max.txt"
      print *,"slice.txt"
      print *,"mode.txt"
      print *,"experiment_out.txt"
      print *,"the lowest possible order of quadrature is used"

      open(unit=17, file= 'experiment_paper.txt')
      open(unit=18, file= 'experiment.txt')
      open(unit=19, file= 'height_integral.txt')
      open(unit=20, file= 'jet_height.txt')
      open(unit=21, file= 'space_average.txt')
      open(unit=22, file= 'time_average.txt')
      open(unit=23, file= 'time_max.txt')
      open(unit=24, file= 'experiment_out.txt')
      open(unit=25, file= 'experiment_paper_out.txt')
      open(unit=26, file= 'slice.txt')
      open(unit=27, file= 'mode.txt')

      do ncomp=1,ncomp_experiment
       read(18,*) xraster,yraster
       xactual=0.048*(xraster-6)/19.0+0.576
       yactual=0.048*(270.0-yraster)/19.0
       write(24,*) xactual,yactual
      enddo
       
      ncomp_experiment=8
      do ncomp=1,ncomp_experiment
       read(17,*) xraster,yraster
       xactual=(2.0/sqrt(3.0))*0.048*(xraster-68)/66.0+0.576
       yactual=0.048*(430.0-yraster)/66.0
       write(25,*) xactual,yactual
      enddo





      line_marker=0

      do ispace=0,TIMESTEPSMAX
       xdata(ispace)=0.0
       tt_average(ispace)=0.0
       tt_max(ispace)=0.0
      enddo
      do itime=0,TIMESTEPSMAX
       tdata(itime)=0.0
       xx_average(itime)=0.0
       slice(itime,1)=0.0
       slice(itime,2)=0.0
      enddo

      itime_slice=0
      prev_time=0.0
      current_time=0.0
      do while (current_time.lt.tstart)
       prev_time=current_time
       read(19,*) current_time,old_space_average

       current_space=0.0
       prev_space=0.0
       do while (current_space.lt.probhix-1.0E-12)
        prev_space=current_space
        read(20,*) current_space,interface_height
        if ((prev_space.lt.slice_position).and. &
            (current_space.ge.slice_position)) then
         slice(itime_slice,2)=interface_height
         slice(itime_slice,1)=current_time
        endif
        line_marker=line_marker+1
       enddo
       read(20,'(A)') garbage
       line_marker=line_marker+1
       print *,"reading garbage at line ",line_marker
       print *,"probhix= ",probhix
     
       itime_slice=itime_slice+1
 
      enddo  ! current_time<tstart

      itime=0

      time_weight=0.0

      do while (current_time.lt.tend)

       prev_time=current_time
       read(19,*) current_time,old_space_average

       print *,"in time loop, prev_time,current_time: ",prev_time, &
         current_time
       print *,"itime= ",itime

       tdata(itime)=prev_time      

       prev_space=0.0
       current_space=0.0

       print *,"CUR,XST ",current_space,xstart

       do while (current_space.lt.xstart)
        prev_space=current_space
        read(20,*) current_space,interface_height
        line_marker=line_marker+1
       enddo
  
       ispace=0

       space_weight=0.0

       print *,"in time loop, first current_space= ",current_space
       print *,"xstart=",xstart
       print *,"line_marker= ",line_marker

       do while (current_space.lt.probhix-1.0E-8)
        prev_space=current_space
        read(20,*) current_space,interface_height
        if ((prev_space.lt.slice_position).and. &
            (current_space.ge.slice_position)) then
         slice(itime_slice,2)=interface_height
         slice(itime_slice,1)=current_time
        endif
        line_marker=line_marker+1

        if (current_space.lt.xend) then
         xdata(ispace)=prev_space
         tt_average(ispace)=tt_average(ispace)+ &
          (current_time-prev_time)*interface_height
         if (interface_height.gt.tt_max(ispace)) then
          tt_max(ispace)=interface_height
         endif
         xx_average(itime)=xx_average(itime)+ &
          (current_space-prev_space)*interface_height
   
         space_weight=space_weight+(current_space-prev_space)

         ispace=ispace+1
        endif  ! current_space<xend

       enddo  ! current_space<probhix

       line_marker=line_marker+1
       print *,"reading garbage at line ",line_marker
       print *,"probhix= ",probhix
       print *,"xend= ",xend
       print *,"current_space=",current_space
       read(20,'(A)') garbage

       nspace=ispace

       xx_average(itime)=xx_average(itime)/space_weight

       time_weight=time_weight+(current_time-prev_time)

       itime=itime+1
       itime_slice=itime_slice+1
      enddo  ! while current_time<tend

      ntime=itime

      do ispace=0,nspace-1
       tt_average(ispace)=tt_average(ispace)/time_weight
       write(22,*) xdata(ispace),tt_average(ispace)
       write(23,*) xdata(ispace),tt_max(ispace)
      enddo

      do itime=0,ntime-1
       write(21,*) tdata(itime),xx_average(itime)
      enddo 
      do itime=0,itime_slice-1
       write(26,*) slice(itime,1),slice(itime,2)
      enddo 

      do j=1,2
      do k=1,2
       A(j,k)=0.0
      enddo
      enddo
      do j=1,2
       B(j)=0.0
       X(j)=0.0
      enddo

      do itime=0,itime_slice-1
       if ((slice(itime,1).gt.tstart).and. &
           (slice(itime,1).lt.tend).and. &
           (itime.gt.0)) then
        A(1,1)=A(1,1)+1.0
        A(1,2)=A(1,2)+slice(itime,1)
        A(2,1)=A(2,1)+slice(itime,1)
        A(2,2)=A(2,2)+slice(itime,1)**2
        B(1)=B(1)+slice(itime,2)
        B(2)=B(2)+slice(itime,1)*slice(itime,2)
       endif
      enddo
      det=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      X(1)=(A(2,2)*B(1)-A(1,2)*B(2))/det
      X(2)=(-A(2,1)*B(1)+A(1,1)*B(2))/det

      L=0.5*(tend-tstart)
      do imode=1,MODEMAX
       AN(imode)=0.0
       BN(imode)=0.0
       CN(imode)=0.0
      enddo
      do itime=0,itime_slice-1
       if ((slice(itime,1).gt.tstart).and. &
           (slice(itime,1).lt.tend).and. &
           (itime.gt.0)) then
        T1=slice(itime-1,1)
        H1=slice(itime-1,2)
        T2=slice(itime,1)
        H2=slice(itime,2)
        TH=0.5*(T1+T2)
        HH=0.5*(H1+H2)
        HH=HH-X(1)-X(2)*TH
        TH=TH-tstart-L
        DT=T2-T1
        do imode=1,MODEMAX
         AN(imode)=AN(imode)+HH*cos(imode*PI*TH/L)*DT 
         BN(imode)=BN(imode)+HH*sin(imode*PI*TH/L)*DT 
        enddo
       endif
      enddo

      do imode=1,MODEMAX
       CN(imode)=sqrt(AN(imode)**2+BN(imode)**2)
       print *,"imode,CN ",imode,CN(imode)
       write(27,*) imode,CN(imode)
      enddo
      imodemax=1
      do imode=2,MODEMAX
       if (CN(imode).gt.CN(imodemax)) then
        imodemax=imode
       endif
      enddo
      print *,"MODEMAX,imodemax ",MODEMAX,imodemax
      print *,"least squares tstart,tend,X(1),X(2) ", &
        tstart,tend,X(1),X(2)

      close(17)
      close(18)
      close(19)
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)

      END PROGRAM
