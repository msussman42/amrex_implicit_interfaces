c******  TEST Subroutine *****
      program Main
      INCLUDE   'grid_def'
      integer TIME_STEPin,nsout(Nr_IBM)
      integer N_TIME_STEPS,ngrab
      logical FIRST_TIMEin
      real*8 TIMEin,DELTA_Tin
      real*8 forcein(2,NS_IBM,Nr_IBM),velocityout(2,NS_IBM,Nr_IBM)
     &      ,positionout(2,NS_IBM,Nr_IBM)

      N_TIME_STEPS=10000
      DELTA_Tin =0.0001

      open (6,FILE='OUTPUT.txt',STATUS="unknown")

      forcein=0.0
      velocityout=0.0
      positionout=0.0
      TIME_STEPin=0
      TIMEin=0.0
      FIRST_TIMEin=.TRUE.

      call Main_FSI(TIME_STEPin,TIMEin,DELTA_Tin,FIRST_TIMEin
     &                      ,forcein,velocityout,positionout,nsout)


      ngrab=0

      CALL SAVE_STATS_2D(ngrab,TIME_STEPin,TIMEin
     &                  ,forcein,velocityout,positionout
     &                  ,nsout,Nr_ibm,Ns_IBM)
         
      FIRST_TIMEin=.FALSE.
      DO TIME_STEPin = TIME_STEPin+1, TIME_STEPin+N_TIME_STEPS
           WRITE(6,*) 'beginning TIME_STEP = ',TIME_STEP
           TIMEin=TIMEin+DELTA_Tin
           Do i=1,Nr_IBM
                 forcein(2,1:nsout(i),i)=-500.0*(1-exp(-TIMEin))
           enddo 
           call Main_FSI(TIME_STEPin,TIMEin,DELTA_Tin,FIRST_TIMEin
     &                  ,forcein,velocityout,positionout,nsout)

           IF (MOD(TIME_STEPin,50).EQ. 0) THEN
              ngrab=ngrab+1
              CALL SAVE_STATS_2D(ngrab,TIME_STEPin,TIMEin
     &                  ,forcein,velocityout,positionout
     &                  ,nsout,Nr_ibm,Ns_IBM)
           END IF
      ENDDO
      WRITE(6,*)
      WRITE(6,*) '        ****** END ******'
      Close(6)

      END
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_2D(ngrab,TIME_STEP,TIME
     &                  ,forcein,velocityout,positionout
     &                  ,nsout,Nr_ibm,Ns_IBM)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      integer Nr_ibm,Ns_IBM
      CHARACTER FNAME*40,supp*4,suppi*3
      integer ngrab
      integer i,j,k,n
      real*8 forcein(2,NS_IBM,Nr_IBM),velocityout(2,NS_IBM,Nr_IBM)
     &      ,positionout(2,NS_IBM,Nr_IBM)

      integer TIME_STEP,nsout(Nr_IBM)
      real*8 TIME

c Write out the mean statistics at each time
             write(supp,'(i4.4)') ngrab
            do i=1,Nr_IBM
             write(suppi,'(i3.3)') i
             Fname='STRUCT_'//suppi//'_'//supp//'.dat'
             !write(*,*) Fname
             IF (ngrab .eq. 1) then 
               open(4000+i,file=Fname,form='formatted',status='unknown')
             ELSE
               close(4000+i)
               open(4000+i,file=Fname,form='formatted',status='unknown')
             END IF
               write(4000+i,498) 'variables="X","Y","U","V","Fx","Fy"'
498         format(1(A37))
	         write(4000+i,599) TIME_STEP,TIME,nsout(i),nsout(i)-1
599            format(1x,'zone T="IBM',i7,F20.9,'" N=',i7,' E=',
     &             i7,' ,ZONETYPE=FELINESEG, DATAPACKING=POINT')
               do j=1,nsout(i)
                  write(4000+i,402) 
     +             positionout(1,j,i)
     +            ,positionout(2,j,i)
     +            ,velocityout(1,j,i)
     +            ,velocityout(2,j,i)
     +            ,forcein(1,j,i)
     +            ,forcein(2,j,i)
               enddo
               do j=1,nsout(i)-1
                    K=J+1
                    write(4000+i,*) k-1,k
               enddo
            enddo 
402      format(24(F20.9,' '))
      return
      end
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE Main_FSI(TIME_STEPin,TIMEin,DELTA_Tin,FIRST_TIMEin
     &                      ,forcein,velocityout,positionout,nsout)
      INCLUDE 'header'


      real*8 forcein(2,NS_IBM,Nr_IBM),velocityout(2,NS_IBM,Nr_IBM)
     &      ,positionout(2,NS_IBM,Nr_IBM)
      real*8 normStruct,tstruct
      INTEGER J,I,iiter,iiter2,nsout(Nr_IBM)
      integer TIME_STEPin
      logical FIRST_TIMEin
      real*8 TIMEin,DELTA_Tin
      real*8 denom



       TIME_STEP=TIME_STEPin
       TIME=TIMEin
       DELTA_T=DELTA_Tin
       FIRST_TIME=FIRST_TIMEin

	 k_flag=.true.
c
c MAIN STRUCTURE
c 
      PI=4.0*datan(1.d0)
      if(FIRST_TIME) then
           write(*,*) 'Opening Input file'
	 open (9997,FILE='inputstruct.dat',STATUS="unknown")
c
c Initialize the IBM part
c 
        CALL INITIALIZE_IBM
        nsout(1:Nr_IBM)= NS_IBM_r(1:Nr_IBM)
           write(*,*) 'Opening Input file'

      endif


      CALL BOUNDARY_IBM

c Formally 2nd-order method with solving structure and fluid seperately and
c match them by feedback control of position and velocity 
! if converged
      DO I=1,Nr_IBM
       DO J=1,NS_IBM_r(i)
          GX_IBMpre(I,J)=GX_IBM(I,J)
          GY_IBMpre(I,J)=GY_IBM(I,J)
          GX_IBM(I,J)=GX_IBM_MASSIVE(I,J)
          GY_IBM(I,J)=GY_IBM_MASSIVE(I,J)
          FK_MASS1O(I,J)=FK_MASS1(I,J)
          FK_MASS2O(I,J)=FK_MASS2(I,J)

       END DO
      END DO
      DO I=1,Nr_IBM
       DO J=1,NS_IBM_r(i)
          GX_IBMo1(I,J)=GX_IBM_MASSIVE(I,J)
          GY_IBMo1(I,J)=GY_IBM_MASSIVE(I,J)

       END DO
      END DO

c 
c Computing predicted X^*=2 * X^{n} - X^{n-1}
      DO I=1,Nr_IBM
       DO J=1,NS_IBM_r(i)
          GX_IBM_MASSIVEO(I,J)=2.0d0*GX_IBM(I,J)-GX_IBMpre(I,J)
          GY_IBM_MASSIVEO(I,J)=2.0d0*GY_IBM(I,J)-GY_IBMpre(I,J)
          GX_IBM_MASSIVE(I,J)=GX_IBM_MASSIVEO(I,J)
          GY_IBM_MASSIVE(I,J)=GY_IBM_MASSIVEO(I,J)
          FS_1_IBMo(I,J)=FS_1_IBM(I,J)
       END DO
      END DO 

      DO I=1,Nr_IBM
         DO J=1,NS_IBM_r(i)
          GX_IBMo1(I,J)=GX_IBM_MASSIVE(I,J)
          GY_IBMo1(I,J)=GY_IBM_MASSIVE(I,J)
         END DO
      END DO

      do iiter=1,maxiter  

      do iiter2=1,maxiterinner-1
          CALL IBM2_PER_PES_1ST(forcein)

          DO I=1,Nr_IBM
             DO J=1,NS_IBM_r(i)
               GX_IBMo1(I,J)=GX_IBM_MASSIVE(I,J)
               GY_IBMo1(I,J)=GY_IBM_MASSIVE(I,J)
             END DO
          END DO
      END DO

      CALL IBM2_PER_PES_1ST(forcein)

      tstruct=0.0
      DO I=1,Nr_IBM
       DO J=1,NS_IBM_r(i)
          normStruct=0.0
          normStruct=normStruct+(GX_IBMo1(I,J)-GX_IBM_MASSIVE(I,J))**2
          normStruct=normStruct+(GY_IBMo1(I,J)-GY_IBM_MASSIVE(I,J))**2
          tstruct=tstruct+normStruct !max(tstruct,sqrt(normStruct))
       END DO
      END DO
      denom=Nr_IBM*Ns_IBM
      normStruct=sqrt(tstruct)/denom !tstruct 

      DO I=1,Nr_IBM
       DO J=1,NS_IBM_r(i)
          GX_IBMo1(I,J)=GX_IBM_MASSIVE(I,J)
          GY_IBMo1(I,J)=GY_IBM_MASSIVE(I,J)
       END DO
      END DO
                     !error check and norm calculation
!      print*,'Conv.', iiter2, normstruct

       print*,'iS', iiter,tstruct,normstruct
      if(normstruct.lt. maxnstruct) then
         if(iiter .ge. 1)  exit
      endif

	end do

      do i=1,nr_ibm
       do j=1,ns_ibm_r(i)
          velocityout(1,j,i)=vibm1(i,j)
          velocityout(2,j,i)=vibm2(i,j)
          positionout(1,j,i)=gx_ibm_massive(i,j)
          positionout(2,j,i)=gy_ibm_massive(i,j)
       end do
      end do
c FINISH MAIN STRUCTURE
      return
      END


c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine create_grid(ds_ibm0,dr_ibm0)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c this subroutine creates the lagrangian grids and their mass properties
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      integer i,j,j2,i_fin,type_tmp1,i2,jstartkb,jendkb,jkb
      integer ixfin(2000),jstartdens,jenddens
      real*8 a_r,r,r0,s,s0,alpha,beta,gama,ds_ibm0,dr_ibm0,temp
     &     ,temp2,temp3,temp4,c_prime,deltatmp,ksh_tmp
     &     ,delta_tmp1,delta_tmp2,dx_tmp1
     &     ,delta_tmpx(2),delta_tmpy(2)
     &     ,delta_tmp3,delta_tmp4,lmax_f,a_r2,dabs_a_r2
      real*8 denom,numerator,parm1,parm2

      pi = 4. * atan(1.0)

      open (9987,file='rigidgrid.dat',status="unknown")
 
      write(*,*) 'start reading'
      read(9997,*) 
      read(9997,*) betaGlobal,gammaGlobal,deltaGlobal,alphaGlobal
      read(9997,*) K_MASSIVE_IBM,  KS_IBM,  KB_IBM, K_LINK
      read(9997,*) CLOSE_BD_flage, para_coor_flag, CREATE_IBM_flag
      read(9997,*) STRU_SIMPLIFY_FLAG
      read(9997,*) Delta_typeX, Delta_typeY, Delta_typeZ, Density_Coef
      read(9997,*) CS_IBM,CB_IBM,CS_IBM_target
      read(9997,*) min_grid_x,min_grid_y,time_max_df
      read(9997,*) 
      MASSIVE_IBM=.true.
      temp4=0.01
      alpha=alphaglobal
      gama= gammaglobal !gamma here length of body 
      deltatmp=deltaglobal
      beta=betaglobal 
      dr_ibm0=0.0
      denom=ns_ibm-1
      dx_tmp1=gammaglobal/denom

      if (nr_ibm .gt. 1) dr_ibm0=0.0
      ds_ibm0=1.d0/denom
      r0= dr_ibm0/2. 
      s0=0.0*ds_ibm0/2.
	  read(9997,*) nr_ibm_f,nr_ibm_fb
	  read(9997,*) maxnstruct,maxiter,maxiterinner,
     &              density_coefp

         lmax_f=gammaglobal

      do i=1,nr_ibm
         ns_ibm_r(i)=0
      enddo
      do i=0,nr_ibm_fb-1
        i_fix=1
    
         read(9997,*)
         type_tmp1=1.0
	  read(9997,*)  flexible_i(i+1)
            read(9997,*) free_bc(i),fix_bc(i)
      read(9997,*) ninteractionlist(i+1)
      if (ninteractionlist(i+1) .gt. 0) then
       read(9997,*) Interactionlist(i+1,1:ninteractionlist(i+1))
      endif
	   ! write(*,*)  flexible_i(i+1)
       ! write(*,*) ninteractionlist(i+1)
       ! write(*,*) Interactionlist(i+1,1:ninteractionlist(i+1))
         if (real(type_tmp1) .gt. 0.5) then
	  read(9997,*)  delta_tmp1,delta_tmp2,delta_tmp3,delta_tmp4
	   ! write(*,*)  delta_tmp1,delta_tmp2,delta_tmp3,delta_tmp4
         a_r  =sqrt( (delta_tmp1-delta_tmp3)**2+
     &               (delta_tmp2-delta_tmp4)**2)

        numerator=ns_ibm
        ns_ibm_r(i+1)=int(numerator/lmax_f*a_r+1d-8)+0

         ! write (*,*) 'Ns_ibm_r for ',i+1,'  is ',ns_ibm_r(i+1)
        denom=ns_ibm_r(i+1)-1
        ds_ibm0=1.d0/denom
        s0=0.0*ds_ibm0/2.
        r=r0+(i)*dr_ibm0
        read(9997,*) a_r2, cs_ibm_r(i+1), kb_ibm_r(i+1), cb_ibm_r(i+1)
       
        if (a_r2.ge.0.0) then
         dabs_a_r2=a_r2
        else 
         dabs_a_r2=-a_r2
        endif
        if(dabs_a_r2 .lt. 1.0d-14) then
           read(9997,*) jstartkb,jendkb
            ! write(*,*) jstartkb,jendkb
           do j=jstartkb,jendkb
              read(9997,*) jkb, temp2
               ! write(*,*) jkb, temp2
              kb_ibm_rs(i+1,jkb)=temp2*kb_ibm_r(i+1)
              cb_ibm_rs(i+1,jkb)=temp2*cb_ibm_r(i+1)
           enddo
           do j=1,jstartkb-1
               kb_ibm_rs(i+1,j)=kb_ibm_rs(i+1,jstartkb)
               cb_ibm_rs(i+1,j)=cb_ibm_rs(i+1,jstartkb)
           enddo
           do j=jendkb+1,ns_IBM_r(i+1)+1
               kb_ibm_rs(i+1,j)=kb_ibm_rs(i+1,jendkb)
               cb_ibm_rs(i+1,j)=cb_ibm_rs(i+1,jendkb)
           enddo
         elseif(a_r2 .gt. 1.0d-14) then
           a_r=a_r2
c variable stiffness
           cs_ibm_r(i+1)=cs_ibm*a_r !------------------------------added damping
           kb_ibm_r(i+1)=kb_ibm*a_r
           cb_ibm_r(i+1)=cb_ibm*a_r
           do j=1,ns_IBM_r(i+1)+1
               kb_ibm_rs(i+1,j)=kb_ibm_r(i+1)
               cb_ibm_rs(i+1,j)=cb_ibm_r(i+1)
           enddo
         else
           do j=1,ns_IBM_r(i+1)+1
               kb_ibm_rs(i+1,j)=kb_ibm_r(i+1)
               cb_ibm_rs(i+1,j)=cb_ibm_r(i+1)
           enddo
          endif
        read(9997,*) a_r2, density_coef_r(i+1),density_coefP_r(i+1)
         ! write(*,*) a_r2, density_coef_r(i+1),density_coefP_r(i+1)
        if (a_r2.ge.0.0) then
         dabs_a_r2=a_r2
        else 
         dabs_a_r2=-a_r2
        endif
        if(dabs_a_r2 .lt. 1.0d-14) then
           read(9997,*) jstartdens,jenddens
           do j=jstartdens,jenddens
              read(9997,*) jkb, temp2
               ! write(*,*) jkb, temp2
               density_coef_rs(i+1,jkb)=temp2* density_coef_r(i+1)
           enddo
           do j=1,jstartdens-1
               density_coef_rs(i+1,j)=density_coef_rs(i+1,jstartdens)
           enddo
           do j=jenddens+1,ns_IBM_r(i+1)+1
               density_coef_rs(i+1,j)=density_coef_rs(i+1,jenddens)
           enddo
         elseif(a_r2 .gt. 1.0d-14) then
           a_r=a_r2
c variable stiffness
           density_coef_r(i+1)=density_coef*a_r !------------------------------added damping
           density_coefP_r(i+1)=density_coefP*a_r !------------------------------added damping
           do j=1,ns_IBM_r(i+1)+1
                density_coef_rs(i+1,j)= density_coef_r(i+1)
           enddo
        else
           density_coef_r(i+1)=density_coef !------------------------------added damping
           density_coefP_r(i+1)=density_coefP !------------------------------added damping
           do j=1,ns_IBM_r(i+1)+1
                density_coef_rs(i+1,j)= density_coef_r(i+1)
           enddo
        endif

        ks_ibm_r(i+1)=ks_ibm*a_r
        k_massive_ibm_r(i+1)=k_massive_ibm*a_r
        temp=0.5*(alpha+beta)


        alpha=alphaglobal+delta_tmp1
        gama= gammaglobal !gamma here length of body 
        deltatmp=deltaglobal
        beta=betaglobal+delta_tmp2
          
          read(9997,*) gt0(i+1,1), gt0(i+1,2)
          read(9997,*) tramp(i+1)
          read(9997,*) ampx(i+1),freqx(i+1),phix(i+1)
          read(9997,*) ampy(i+1),freqy(i+1),phiy(i+1)
          read(9997,*) ampt(i+1),freqt(i+1),phit(i+1)           

 
          phix(i+1)=pi*  phix(i+1)
          phiy(i+1)=pi*  phiy(i+1)
          phit(i+1)=pi*  phit(i+1)


	  read(9997,*) frequency_fin(i+1)
	  read(9997,*) i_fin
	  read(9997,*) ixfin(1:i_fin)    
 
        do j=1,i_fin
              if(ixfin(j) .le. 0) then
              target_point_num(i+1,j)=ns_ibm_r(i+1)+ixfin(j)+0
              if(target_point_num(i+1,j) .gt. ns_ibm_r(i+1)-2) then
                write(*,*) 'Error in inputfin and targetpoints define'
                write(*,*) 'point is wrong: target',
     &          target_point_num(i+1,j),'> ns_ibm of bodyi-2',
     &          ns_ibm_r(i+1)-2,' for body',i+1,' ,point',
     &          j
                exit
              endif
              else
              target_point_num(i+1,j)=ixfin(j)+0
              endif
              i_fix=i_fix+1
        end do





        gx_ibm(i+1,1)=delta_tmp1
        gy_ibm(i+1,1)=delta_tmp2

        do j=0,ns_ibm_r(i+1)-1

           numerator=j
           s=s0+numerator*ds_ibm0
           gx_ibm(i+1,j+1)=delta_tmp1+(delta_tmp3-delta_tmp1)*s

           gy_ibm(i+1,j+1)=delta_tmp2+(delta_tmp4-delta_tmp2)*s  !dsin(2.0d0*pi*s)

           gx_bp(i+1,j+1)=gx_ibm(i+1,j+1)
           gy_bp(i+1,j+1)=gy_ibm(i+1,j+1)
		 if (j .eq. 0) then 
		    temp2=gx_bp(i+1,1)-gx_ibm(i+1,1)
		    temp3=gy_bp(i+1,1)-gy_ibm(i+1,1)
		 end if
		 gx_ibm(i+1,j+1)=temp2+gx_ibm(i+1,j+1)
		 gy_ibm(i+1,j+1)=temp3+gy_ibm(i+1,j+1)
           gx_bp(i+1,j+1)=gx_ibm(i+1,j+1)
           gy_bp(i+1,j+1)=gy_ibm(i+1,j+1)
           gx_ibmpre(i+1,j+1)=gx_ibm(i+1,j+1)
           gy_ibmpre(i+1,j+1)=gy_ibm(i+1,j+1)
           if (massive_ibm) mass_ibm(i+1,j+1)=1.0*ds_ibm0
        end do

!         print*, target_point_num
         target_num(i+1)=i_fix-1
          read(9997,*) target_k_link(i+1,1:target_num(i+1)) 
          read(9997,*) target_t_link(i+1,1:target_num(i+1))   
          read(9997,*) a_fin(i+1,1:target_num(i+1))    
          read(9997,*) phi_fin(i+1,1:target_num(i+1))   
           ! write(*,*) target_k_link(i+1,1:target_num(i+1)) 
           ! write(*,*) target_t_link(i+1,1:target_num(i+1))   
           ! write(*,*) a_fin(i+1,1:target_num(i+1))    
           ! write(*,*) phi_fin(i+1,1:target_num(i+1))   
          do j=1, target_num(i+1)
               phi_fin(i+1,j)=pi* phi_fin(i+1,j)
          enddo
          read(9997,*) flagfixed(i+1,1:target_num(i+1))   
           ! write(*,*) flagfixed(i+1,1:target_num(i+1))   

      else
	  read(9997,*)  delta_tmpx(1),delta_tmpy(1)
     &               ,delta_tmpx(2),delta_tmpy(2)
         call segment(delta_tmpx,delta_tmpy,dx_tmp1,i)
      endif
      enddo
        read(9997,*)  
        read(9997,*) i_fix
         ! write(*,*) i_fix


      do i=1,i_fix
        read(9987,*) i2,j, temp2, temp3
         ! write(*,*) i2,j, temp2, temp3
           gx_ibm(i2,j)=temp2
           gy_ibm(i2,j)=temp3
           parm1=ns_ibm_r(i2)
           parm2=j
           ns_ibm_r(i2)=int(max(parm1,parm2)+1.0d-5)
           gx_bp(i2,j)=gx_ibm(i2,j)
           gy_bp(i2,j)=gy_ibm(i2,j)
           gx_ibmpre(i2,j)=gx_ibm(i2,j)
           gy_ibmpre(i2,j)=gy_ibm(i2,j)
      enddo

      do i=nr_ibm_fb,nr_ibm-1
        a_r=1.0   ! +dsin(2*pi1*r-pi1/2.)
c variable stiffness
        ks_ibm_r(i+1)=ks_ibm*a_r
        cs_ibm_r(i+1)=cs_ibm*a_r !------------------------------added damping
        kb_ibm_r(i+1)=kb_ibm*a_r
        cb_ibm_r(i+1)=cb_ibm*a_r
        do j=1,ns_IBM_r(i+1)+1
               kb_ibm_rs(i+1,j)=kb_ibm_r(i+1)
               cb_ibm_rs(i+1,j)=cb_ibm_r(i+1)
        enddo
        k_massive_ibm_r(i+1)=k_massive_ibm*a_r
        Gt0(i+1,1)=gx_bp(i+1,1)
        Gt0(i+1,2)=gy_bp(i+1,1)
        Tramp(i+1)=0.0
        frequency_fin(i+1)=0.0
        ampx(i+1)=0.0
        freqx(i+1)=0.0
        phix(i+1)=0.0
        ampy(i+1)=0.0
        freqy(i+1)=0.0
        phiy(i+1)=0.0
        ampt(i+1)=0.0
        freqt(i+1)=0.0
        phit(i+1)=0.0          
          read(9997,*) i2
          read(9997,*) gt0(i2,1), gt0(i2,2)
          read(9997,*) tramp(i2)
          read(9997,*) ampx(i2),freqx(i2),phix(i2)
          read(9997,*) ampy(i2),freqy(i2),phiy(i2)
          read(9997,*) ampt(i2),freqt(i2),phit(i2)           

          phix(i2)=pi*  phix(i2)
          phiy(i2)=pi*  phiy(i2)
          phit(i2)=pi*  phit(i2)
      enddo

      do i=1,nr_ibm
         print*,'body',i,'with',ns_ibm_r(i),'points'!',kb_ibm_r(i)
       enddo

       close(9998)
      return
      end
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine segment(delta_tmpx,delta_tmpy,Dx_tmp1,i)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c This subroutine creates the lagrangian grids and their mass properties
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      real*8 delta_tmpx(2),delta_tmpy(2),Dx_tmp1,L_tmp,dx_tmp2
      integer i , ni,j
      real*8 denom,dj
           L_tmp=sqrt((delta_tmpx(2)-delta_tmpx(1))**2+
     &                (delta_tmpy(2)-delta_tmpy(1))**2)
           ni=min(int((L_tmp-DX_tmp1)/Dx_tmp1)+1,Ns_IBM)
           denom=ni
           dx_tmp2=(L_tmp)/denom

          do j=0,ni-1
              dj=j
		 GX_IBM(i+1,j+1)=delta_tmpx(1)
     &        +(dj+0.5)/denom
     &        *(delta_tmpx(2)-delta_tmpx(1))
		 GY_IBM(i+1,j+1)=delta_tmpy(1)
     &        +(dj+0.5)/denom
     &        *(delta_tmpy(2)-delta_tmpy(1))
           GX_BP(i+1,j+1)=GX_IBM(i+1,j+1)
           GY_BP(i+1,j+1)=GY_IBM(i+1,j+1)
           GX_IBMpre(i+1,j+1)=GX_IBM(i+1,j+1)
           GY_IBMpre(i+1,j+1)=GY_IBM(i+1,j+1)
           enddo
      Ns_IBM_r(i+1)=ni
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine initialize_ibm
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c this subroutine sets the initial position of immersed boundary and their mass
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      integer nr_ibm_t,ns_ibm_t
      integer i,j,k
      real*8  coor1(2),coor2(2),dist,ds_ibm0,dr_ibm0
      write(6,*) 'ibm grid size: nr_ibm =',nr_ibm
      write (6,*) 'coordinate of ibm points in r'
      write(6,*) 'ibm grid size: nr_ibm =',ns_ibm
      write (6,*) 'coordinate of ibm points in r'
c initialized ibm grid or reading saved grid file
       call create_grid(ds_ibm0,dr_ibm0)
      if (create_ibm_flag) then

        open (14,file='gridibm.txt',form='formatted',status='unknown')
        write(14,*) nr_ibm,ns_ibm
        do i=1,nr_ibm
           do j=1,ns_ibm
             fk_mass1(i,j)=0.0d0
             fk_mass2(i,j)=0.0d0
             write(14,*) gx_ibm(i,j),gy_ibm(i,j)
           end do
        end do
        close(14)
        if (massive_ibm) then 
         open (15,file='massibm.txt',form='formatted',status='unknown')
         write(15,*) nr_ibm,ns_ibm
         do i=1,nr_ibm
           do j=1,ns_ibm
             write(15,*) mass_ibm(i,j)
           end do
         end do
         close(15)
        end if
        open (12,file='gridibm.txt',form='formatted',status='old')
        read (12,*) nr_ibm_t,ns_ibm_t
      end if        

c check to make sure that ibm input file is the correct dimensions
         if (nr_ibm_t .ne. nr_ibm) then
             write(6,*) 'nr_ibm, nr_ibm_t',nr_ibm,nr_ibm_t
             stop 'error: r input in ibm wrong dimensions'
         else if (ns_ibm_t .ne. ns_ibm) then
             write(6,*) 'ns_ibm, ns_ibm_t',nr_ibm,nr_ibm_t
             stop 'error: s input in ibm wrong dimensions'
         end if
c read lagrangian coord. massless part
         do i=1,nr_ibm
           do j=1,ns_ibm
             read(12,*) gx_ibm(i,j),gy_ibm(i,j)
                 write(6,*) 'gx_ibm(',i,j,') = ',gx_ibm(i,j)
                 write(6,*) 'gy_ibm(',i,j,') = ',gy_ibm(i,j)
           end do
         end do
 
c initialize lagrangian coord. massive part
      if (massive_ibm) then
            do i=1,nr_ibm
              do j=1,ns_ibm
                gx_ibm_massive(i,j)=gx_ibm(i,j)
                gy_ibm_massive(i,j)=gy_ibm(i,j)
                gx_ibmpre(i,j)=gx_ibm(i,j)
                gy_ibmpre(i,j)=gy_ibm(i,j)
              end do
            end do
         open (13,file='massibm.txt',form='formatted',status='old')
         read(13,*) nr_ibm_t , ns_ibm_t
c check to make sure that ibm input file is the correct dimensions
         if (nr_ibm_t .ne. nr_ibm) then
             write(6,*) 'nr_ibm, nr_ibm_t',nr_ibm,nr_ibm_t
             stop 'error: r input in ibm wrong dimensions'
         else if (ns_ibm_t .ne. ns_ibm) then
             write(6,*) 'ns_ibm, ns_ibm_t',ns_ibm,ns_ibm_t
             stop 'error: s input in ibm wrong dimensions'
         end if
c read mass of each node in ibm
         do i=1,nr_ibm
           do j=1, ns_ibm
              read(13,*) mass_ibm(i,j)
                write(6,*) 'mass_ibm(',i,j,') = ',mass_ibm(i,j)
           end do
         end do
         close(13) 
      end if 
c satistfying periodic boundary for closed 2d surface
      if (close_bd_flage) then 
          do i=1,nr_ibm
            gx_ibm(i,0)       =  gx_ibm(i,ns_ibm)      
            gx_ibm(i,ns_ibm+1)=  gx_ibm(i,1) 
            gy_ibm(i,0)       =  gy_ibm(i,ns_ibm)      
            gy_ibm(i,ns_ibm+1)=  gy_ibm(i,1) 
            gx_bp(i,0)       =  gx_bp(i,ns_ibm)      
            gx_bp(i,ns_ibm+1)=  gx_bp(i,1) 
            gy_bp(i,0)       =  gy_bp(i,ns_ibm)      
            gy_bp(i,ns_ibm+1)=  gy_bp(i,1) 
            if (massive_ibm) then    
              gx_ibm_massive(i,0)       = gx_ibm_massive(i,ns_ibm)      
              gx_ibm_massive(i,ns_ibm+1)= gx_ibm_massive(i,1) 
              gy_ibm_massive(i,0)       = gy_ibm_massive(i,ns_ibm)      
              gy_ibm_massive(i,ns_ibm+1)= gy_ibm_massive(i,1) 
            end if
          end do
            s_start=0 
            s_end=ns_ibm+1
       else
          s_start=1 
          s_end=ns_ibm
       end if
c find the distance between ibm points and choose lagrangian coordinate

        if (para_coor_flag .eqv. .false.) then
          do i=1,nr_ibm
            do j=1,ns_ibm_r(i)-1
             coor1(1)=gx_bp(i,j)
             coor1(2)=gy_bp(i,j)
             coor2(1)=gx_bp(i,j+1)
             coor2(2)=gy_bp(i,j+1)
             dsf_ibm(i,j)=dist(coor1,coor2)
           enddo
           do j=ns_ibm_r(i),ns_ibm+2
              dsf_ibm(i,j)=dsf_ibm(i,j-1)  
           end do
           dsf_ibm(i,0)=dsf_ibm(i,1)
           dsf_ibm(i,-1)=dsf_ibm(i,1)
           dsf_ibm(i,-2)=dsf_ibm(i,1)
          end do
c find the staggered distance between ibm points  
          do i=1,nr_ibm
           do j=-1,ns_ibm
            ds_ibm(i,j)=(dsf_ibm(i,j-1) + dsf_ibm(i,j)) / 2.0d0
           end do
          end do
          s_start=1 
          s_end=ns_ibm
        else
         do i=1,nr_ibm
           do j=1,ns_ibm
             dsf_ibm(i,j)=ds_ibm0
             ds_ibm(i,j) =ds_ibm0
           end do
           dsf_ibm(i,0)=dsf_ibm(i,1)
           dsf_ibm(i,s_end)=dsf_ibm(i,s_end-1)
         end do
        end if
c here force for feedback algorithm is calculated 
        do i=1,nr_ibm
           do j=1,ns_ibm_r(i)-1
c computing (dx.dx)^{n}
              ds_ibm0=(gx_ibm(i,j+1)-gx_ibm(i,j))/dsf_ibm(i,j)
              dr_ibm0=(gy_ibm(i,j+1)-gy_ibm(i,j))/dsf_ibm(i,j) 
c saving computed term for using later time
              fs_1_ibm(i,j)=ds_ibm0**2.+dr_ibm0**2.
           end do
        end do
c initializing velocities for feedback algorithm
	print*, 'initialization of vibm1'
        vibm1=0.0d0
        vibm1_pre=0.0d0
        vibm2=0.0d0
        vibm2_pre=0.0d0

c determine the width of delta function
      call delta_width(del_x,del_y,del_z
     &                ,delta_typex,delta_typey,delta_typez)
      close(12) 
      return
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine ibm2_per_pes_1st(forcein)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      integer nrin,nsin
      real*8 forcein(2,NS_IBM,Nr_IBM)

      integer i,j,k,para,t_temp,n,j2, ni2Penalty,i2p
      integer ierrksh1,iksh1,iksh2,iksh3,i2
      real*8 aamat(-3:ns_ibm+3,-3:ns_ibm+3),bbmat(ns_ibm,2)
      real*8 temp1, temp2, temp3, temp4,ksh_tmp
      real*8 temp_ibm1,fn1,fn2,temp_ibm2
     &      ,temp_check1,temp_check2,temp_check3, temp2_ibm, temp3_ibm
      real*8 temp4_1_ibm(nr_ibm,0:ns_ibm),temp4_2_ibm(nr_ibm,0:ns_ibm)
     &      ,a(ns_ibm,-3:3) !     ,b(ns_ibm),c(ns_ibm)
     &      ,mat(4,4)
     &      ,tdx0_1(nr_ibm),tdx0_2(nr_ibm)
     &      ,tdxnm1_1(nr_ibm),tdxnm1_2(nr_ibm)
     &      ,coor1(2),coor2(2),
     &       dist,S_Alpha
      real*8       suminertia(2),sumviscous(2),sumflink(2),
     &           sumtension(2),sumfbend(2),sumfluid(2),
     &           ffluidsum1(nr_ibm,0:ns_ibm),
     &           ffluidsum2(nr_ibm,0:ns_ibm), timetemp2
      logical implicitflag,implicitflagvi
      integer itmax,precond,success,itconverge
      real*8 tol, kb_ibmt,kb_ibmtvi,dn,dkj2

      parameter (tol = 1.0d-11, precond=1, itmax=300)
c free:           free_bc =.true.    fix_bc =.false.
c Hinge:          free_bc =.false.   fix_bc =.false.
c clamp:          free_bc =.false.   fix_bc =.true.
      parameter( implicitflagvi=.true.,implicitflag=.true.) 
c time step
      temp4=delta_t
      timetemp2=(1.d0-exp(-time2/0.1))           
      eps=1.0d-8
      temp1=temp4**2.0

c calculating feedback force for matching bcs
      do i=1,nr_ibm
        do j=1,ns_ibm
           if(j .le. ns_ibm_r(i)) then
            fibm1(i,j)= forcein(1,j,i)
            fibm2(i,j)= forcein(2,j,i)
           endif
 
           vibm1_pre(i,j)=vibm1(i,j)
           vibm2_pre(i,j)=vibm2(i,j)

!add jan06_13 for force calculations
           ffluidsum1(i,j)=fibm1(i,j)
           ffluidsum2(i,j)=fibm2(i,j)

        end do
      end do

c-------------------------------------------------------------
c	temp_ibm1=0.0;temp_check2=0.0
c      do i=1,nr_ibm
c        do j=1,ns_ibm
c		temp_ibm1=temp_ibm1+fibm1(i,j)*ds_ibm(i,j)
c		temp_check2=temp_check2+fibm2(i,j)*ds_ibm(i,j)
c        end do
c      end do
c	 print*,temporary_sum,temp_ibm1;print*,temp_check1,
c       temp_check2;pause
c-------------------------------------------------------------
      n=1
      do t_temp=1,n   !t_temp
      dn=n
      temp4=delta_t/dn
      temp1=temp4**2.0
c imposing bcs
          call boundary_ibm
!      goto 440
c 
c initialize the bending force in ibm
      do i=1,nr_ibm
        do j=1,ns_ibm
         fb_1_ibm(i,j)=0.0d0
         fb_2_ibm(i,j)=0.0d0
         fbv_1_ibm(i,j)=0.0d0
         fbv_2_ibm(i,j)=0.0d0
         fbvi_1_ibm(i,j)=0.0d0
         fbvi_2_ibm(i,j)=0.0d0
         f_impuls1(i,j)=0.d0
         f_impuls2(i,j)=0.d0
        end do
      end do
c 
c compute bending force in ibm
c here we assume ( small deformation + neglecting 2nd order )
        f_link1=0.0
        f_link2=0.0

      do i=1,nr_ibm_f
        temp_ibm1=(target_points(i,2,1)-target_points(i,1,1))
        temp_ibm2=(target_points(i,2,2)-target_points(i,1,2))
        temp2_ibm=( temp_ibm1 ** 2. + temp_ibm2 ** 2.+eps ) ** 0.5
        temp3_ibm=temp_ibm1/temp2_ibm
        temp2_ibm=temp_ibm2/temp2_ibm


        do j=1,ns_ibm_r(i)
           kb_ibm=kb_ibm_rs(i,j) 
          if (j .eq. 1) then
c depending on fixed boundary condition we have different formulation
            if(.not.(free_bc(i))) then  
               if ( fix_bc(i) ) then
c (d_ssx)_1=[ dx_{1+1/2} - (1,0) } / 0.5 dsf_1
c also  dx_{1+1/2}=[x_2-x_1]/dsf_1
                  temp4_1_ibm(i,j)= kb_ibm *
     &          ( (gx_ibm_massive(i,j+1)-gx_ibm_massive(i,j))
     &            /dsf_ibm(i,j)-temp3_ibm )
     &            / ( 0.50d0 * dsf_ibm(i,j) )
                  temp4_2_ibm(i,j)= kb_ibm *
     &          ( (gy_ibm_massive(i,j+1)-gy_ibm_massive(i,j))
     &            /dsf_ibm(i,j)-temp2_ibm )
     &            / ( 0.50d0 * dsf_ibm(i,j) )
               else                 
                  temp4_1_ibm(i,j)=0.0d0
                  temp4_2_ibm(i,j)=0.0d0
               end if
            else
               temp4_1_ibm(i,j)=0.d0
               temp4_2_ibm(i,j)=0.d0
            endif
          elseif (j .eq. ns_ibm_r(i)) then 
            temp4_1_ibm(i,j)=0.0d0
            temp4_2_ibm(i,j)=0.0d0
          else            
            temp4_1_ibm(i,j)= kb_ibm *
     &    ( (gx_ibm_massive(i,j+1)-gx_ibm_massive(i,j))/dsf_ibm(i,j)
     &     -(gx_ibm_massive(i,j)-gx_ibm_massive(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
            temp4_2_ibm(i,j)=  kb_ibm *         
     &    ( (gy_ibm_massive(i,j+1)-gy_ibm_massive(i,j))/dsf_ibm(i,j)
     &     -(gy_ibm_massive(i,j)-gy_ibm_massive(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
          end if

        end do
           kb_ibm=kb_ibm_rs(i,1)
        tdx0_1(i)= kb_ibm*temp4_1_ibm(i,2)/dsf_ibm(i,1)
        tdx0_2(i)= kb_ibm*temp4_2_ibm(i,2)/dsf_ibm(i,1)

           kb_ibm=kb_ibm_rs(i,ns_ibm_r(i))
        tdxnm1_1(i)= -kb_ibm*temp4_1_ibm(i,ns_ibm_r(i)-1)           
     &    /dsf_ibm(i,ns_ibm_r(i)-1)
        tdxnm1_2(i)= -kb_ibm*temp4_2_ibm(i,ns_ibm_r(i)-1)           
     &    /dsf_ibm(i,ns_ibm_r(i)-1)


        do j=2,ns_ibm_r(i)-1
           fb_1_ibm(i,j)= - 1.0*
     &    ( (temp4_1_ibm(i,j+1)-temp4_1_ibm(i,j))/dsf_ibm(i,j)
     &     -(temp4_1_ibm(i,j)-temp4_1_ibm(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
           fb_2_ibm(i,j)= - 1.0 *
     &    ( (temp4_2_ibm(i,j+1)-temp4_2_ibm(i,j))/dsf_ibm(i,j)
     &     -(temp4_2_ibm(i,j)-temp4_2_ibm(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
         end do
         j = ns_ibm_r(i)
         fb_1_ibm(i,j)= 1.0*
     &     (2.0*temp4_1_ibm(i,j-1)-temp4_1_ibm(i,j-2))/dsf_ibm(i,j-1)
     &      / dsf_ibm(i,j-2)
         fb_2_ibm(i,j)= 1.0 *
     &     (2.0*temp4_2_ibm(i,j-1)-temp4_2_ibm(i,j-2))/dsf_ibm(i,j-1)
     &      / dsf_ibm(i,j-2)
         if(free_bc(i)) then
             j = 1
             fb_1_ibm(i,j)= 1.0 *
     &       (2.0*temp4_1_ibm(i,j+1)-temp4_1_ibm(i,j+2))/dsf_ibm(i,j+1)
     &          / dsf_ibm(i,j+2)
             fb_2_ibm(i,j)= 1.0 *
     &       (2.0*temp4_2_ibm(i,j+1)-temp4_2_ibm(i,j+2))/dsf_ibm(i,j+1)
     &         / dsf_ibm(i,j+2)
         else
             fb_1_ibm(i,1)= 0.0d0
             fb_2_ibm(i,1)= 0.0d0
        endif

!bending damping ++++++++++++++++++++++++++++++++++++++++++++++++++
         if (implicitflagvi) then
         do j=1,ns_ibm_r(i)
           cb_ibm=cb_ibm_rs(i,j) 
          if (j .eq. 1) then
c depending on fixed boundary condition we have different formulation
            if(.not.(free_bc(i))) then  
               if ( fix_bc(i) ) then
c (d_ssx)_1=[ dx_{1+1/2} - (1,0) } / 0.5 dsf_1
c also  dx_{1+1/2}=[x_2-x_1]/dsf_1
                  temp4_1_ibm(i,j)= cb_ibm *
     &          ( (gx_ibm(i,j+1)-gx_ibm(i,j))
     &            /dsf_ibm(i,j)-temp3_ibm*0.0 )
     &            / ( 0.50d0 * dsf_ibm(i,j) )/temp4
                  temp4_2_ibm(i,j)= cb_ibm *
     &          ( (gy_ibm(i,j+1)-gy_ibm(i,j))
     &            /dsf_ibm(i,j)-temp2_ibm*0.0 )
     &            / ( 0.50d0 * dsf_ibm(i,j) )/temp4
               else                 
                  temp4_1_ibm(i,j)=0.0d0
                  temp4_2_ibm(i,j)=0.0d0
               end if
            else
               temp4_1_ibm(i,j)=0.d0
               temp4_2_ibm(i,j)=0.d0
            endif
          elseif (j .eq. ns_ibm_r(i)) then 
            temp4_1_ibm(i,j)=0.0d0
            temp4_2_ibm(i,j)=0.0d0
          else            
            temp4_1_ibm(i,j)= cb_ibm *
     &    ( (gx_ibm(i,j+1)-gx_ibm(i,j))/dsf_ibm(i,j)
     &     -(gx_ibm(i,j)-gx_ibm(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)/temp4
            temp4_2_ibm(i,j)=  cb_ibm *         
     &    ( (gy_ibm(i,j+1)-gy_ibm(i,j))/dsf_ibm(i,j)
     &     -(gy_ibm(i,j)-gy_ibm(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)/temp4
          end if

        end do
           cb_ibm=cb_ibm_rs(i,1)

        do j=2,ns_ibm_r(i)-1
           fbvi_1_ibm(i,j)= - 1.0*
     &    ( (temp4_1_ibm(i,j+1)-temp4_1_ibm(i,j))/dsf_ibm(i,j)
     &     -(temp4_1_ibm(i,j)-temp4_1_ibm(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
           fbvi_2_ibm(i,j)= - 1.0 *
     &    ( (temp4_2_ibm(i,j+1)-temp4_2_ibm(i,j))/dsf_ibm(i,j)
     &     -(temp4_2_ibm(i,j)-temp4_2_ibm(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
         end do
         j = ns_ibm_r(i)
         fbvi_1_ibm(i,j)= 1.0*
     &     (2.0*temp4_1_ibm(i,j-1)-temp4_1_ibm(i,j-2))/dsf_ibm(i,j-1)
     &      / dsf_ibm(i,j-2)
         fbvi_2_ibm(i,j)= 1.0 *
     &     (2.0*temp4_2_ibm(i,j-1)-temp4_2_ibm(i,j-2))/dsf_ibm(i,j-1)
     &      / dsf_ibm(i,j-2)
         if(free_bc(i)) then
             j = 1
             fbvi_1_ibm(i,j)= 1.0 *
     &       (2.0*temp4_1_ibm(i,j+1)-temp4_1_ibm(i,j+2))/dsf_ibm(i,j+1)
     &          / dsf_ibm(i,j+2)
             fbvi_2_ibm(i,j)= 1.0 *
     &       (2.0*temp4_2_ibm(i,j+1)-temp4_2_ibm(i,j+2))/dsf_ibm(i,j+1)
     &         / dsf_ibm(i,j+2)
         else
             fbvi_1_ibm(i,1)= 0.0d0
             fbvi_2_ibm(i,1)= 0.0d0
        endif


        else
        do j=1,ns_ibm
         fbvi_1_ibm(i,j)=0.0d0
         fbvi_2_ibm(i,j)=0.0d0

        end do

        endif



        do j=1,ns_ibm_r(i)
           cb_ibm=cb_ibm_rs(i,j) 
          if (j .eq. 1) then
c depending on fixed boundary condition we have different formulation
            if(.not.(free_bc(i))) then  
               if ( fix_bc(i) ) then
c (d_ssx)_1=[ dx_{1+1/2} - (1,0) } / 0.5 dsf_1
c also  dx_{1+1/2}=[x_2-x_1]/dsf_1
                  temp4_1_ibm(i,j)= cb_ibm *
     &          ( (vibm1(i,j+1)-vibm1(i,j))
     &            /dsf_ibm(i,j)-temp3_ibm*0.0)
     &            / ( 0.50d0 * dsf_ibm(i,j) )
                  temp4_2_ibm(i,j)= cb_ibm *
     &          ( (vibm2(i,j+1)-vibm2(i,j))
     &            /dsf_ibm(i,j)-temp2_ibm*0.0)
     &            / ( 0.50d0 * dsf_ibm(i,j) )
               else                 
                  temp4_1_ibm(i,j)=0.0d0
                  temp4_2_ibm(i,j)=0.0d0
               end if
            else
               temp4_1_ibm(i,j)=0.d0
               temp4_2_ibm(i,j)=0.d0
            endif
          elseif (j .eq. ns_ibm_r(i)) then 
            temp4_1_ibm(i,j)=0.0d0
            temp4_2_ibm(i,j)=0.0d0
          else            
            temp4_1_ibm(i,j)= cb_ibm *
     &    ( (vibm1(i,j+1)-vibm1(i,j))/dsf_ibm(i,j)
     &     -(vibm1(i,j)-vibm1(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
            temp4_2_ibm(i,j)=  cb_ibm *         
     &    ( (vibm2(i,j+1)-vibm2(i,j))/dsf_ibm(i,j)
     &     -(vibm2(i,j)-vibm2(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
          end if

        end do
           cb_ibm=cb_ibm_rs(i,1)
        tdx0_1(i)= tdx0_1(i)+cb_ibm*temp4_1_ibm(i,2)/dsf_ibm(i,1)
        tdx0_2(i)= tdx0_2(i)+cb_ibm*temp4_2_ibm(i,2)/dsf_ibm(i,1)

           cb_ibm=cb_ibm_rs(i,ns_ibm_r(i))
        tdxnm1_1(i)=tdxnm1_1(i)-cb_ibm*temp4_1_ibm(i,ns_ibm_r(i)-1)           
     &    /dsf_ibm(i,ns_ibm_r(i)-1)
        tdxnm1_2(i)=tdxnm1_2(i)-cb_ibm*temp4_2_ibm(i,ns_ibm_r(i)-1)           
     &    /dsf_ibm(i,ns_ibm_r(i)-1)


        do j=2,ns_ibm_r(i)-1
           fbv_1_ibm(i,j)= - 1.0*
     &    ( (temp4_1_ibm(i,j+1)-temp4_1_ibm(i,j))/dsf_ibm(i,j)
     &     -(temp4_1_ibm(i,j)-temp4_1_ibm(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
           fbv_2_ibm(i,j)= - 1.0 *
     &    ( (temp4_2_ibm(i,j+1)-temp4_2_ibm(i,j))/dsf_ibm(i,j)
     &     -(temp4_2_ibm(i,j)-temp4_2_ibm(i,j-1))/dsf_ibm(i,j-1) )
     &      / ds_ibm(i,j)
         end do
         j = ns_ibm_r(i)
         fbv_1_ibm(i,j)= 1.0*
     &     (2.0*temp4_1_ibm(i,j-1)-temp4_1_ibm(i,j-2))/dsf_ibm(i,j-1)
     &      / dsf_ibm(i,j-2)
         fbv_2_ibm(i,j)= 1.0 *
     &     (2.0*temp4_2_ibm(i,j-1)-temp4_2_ibm(i,j-2))/dsf_ibm(i,j-1)
     &      / dsf_ibm(i,j-2)
         if(free_bc(i)) then
             j = 1
             fbv_1_ibm(i,j)= 1.0 *
     &       (2.0*temp4_1_ibm(i,j+1)-temp4_1_ibm(i,j+2))/dsf_ibm(i,j+1)
     &          / dsf_ibm(i,j+2)
             fbv_2_ibm(i,j)= 1.0 *
     &       (2.0*temp4_2_ibm(i,j+1)-temp4_2_ibm(i,j+2))/dsf_ibm(i,j+1)
     &         / dsf_ibm(i,j+2)
         else
             fbv_1_ibm(i,1)= 0.0d0
             fbv_2_ibm(i,1)= 0.0d0
        endif
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

             ni2Penalty=nInteractionList(i)
             
             if (ni2Penalty.gt. 0) then
             do j2=1,ns_ibm_r(i)
               do i2p=1,ni2Penalty
                    i2=InteractionList(i,i2p)
                     do j=1,ns_ibm_r(i2)
                        if((
     &    abs((gy_ibm_massive(i,j2)-gy_ibm_massive(i2,j))/min_grid_y)
     &    .le.4.0) .and.
     &    (abs((gx_ibm_massive(i,j2)-gx_ibm_massive(i2,j))/min_grid_x)
     &    .le.4.0) ) then
                         temp_ibm1= 
     &    abs((gy_ibm_massive(i,j2)-gy_ibm_massive(i2,j))/min_grid_y)
                           call delta_fun(delta_typey,temp_ibm1,fn2)
                        temp_ibm2= 
     &    abs((gx_ibm_massive(i,j2)-gx_ibm_massive(i2,j))/min_grid_x)
                           call delta_fun(delta_typex,temp_ibm2,fn1)
                        temp3_ibm=
     &                  sqrt((temp_ibm1*min_grid_y)**2
     &                      +(temp_ibm2*min_grid_x)**2)+1.0d-9
                 f_impuls1(i,j2)= f_impuls1(i,j2)+
     &                               fn1*fn2
     &                        *ds_ibm(i,j2)
     &                        /(min_grid_x*min_grid_y)
     &    *(gx_ibm_massive(i,j2)-gx_ibm_massive(i2,j))/temp3_ibm
                 f_impuls2(i,j2)= f_impuls2(i,j2)+
     &                               fn1*fn2
     &                        *ds_ibm(i,j2)
     &                        /(min_grid_x*min_grid_y)
     &    *(gy_ibm_massive(i,j2)-gy_ibm_massive(i2,j))/temp3_ibm
                         endif
                    
                   end do
               enddo
             end do 
          endif
       enddo

!             print*,maxval(f_impuls1),maxval(f_impuls2),temp3_ibm

      do i=1,nr_ibm_f
        density_coefp= density_coefp_r(i) 
        do j=1,ns_ibm_r(i)
           finer_1_ibm(i,j)= density_coefp *
     &    ( gx_ibmo1(i,j)-gx_ibm_massiveo(i,j) )/temp1
           finer_2_ibm(i,j)= density_coefp *
     &    ( gy_ibmo1(i,j)-gy_ibm_massiveo(i,j) )/temp1
         end do

         force_points=0.0
         do j=1,target_num(i)
            k=target_point_num(i,j)
             temp_ibm1=(gx_ibm_massive(i,k)-target_points(i,j,1))
             temp_ibm2=(gy_ibm_massive(i,k)-target_points(i,j,2))
             
             if(time2 .le. target_t_link(i,j)) then
             ksh_tmp=target_k_link(i,j)
             else
             ksh_tmp=0.0
             endif
             force_points(i,j,1)=k_link*temp_ibm1*ksh_tmp
             force_points(i,j,2)=k_link*temp_ibm2*ksh_tmp

             temp_ibm1=target_points_v(i,j,1)
             temp_ibm2=target_points_v(i,j,2)

             temp_ibm1=vibm1_pre(i,k)-temp_ibm1
             temp_ibm2=vibm2_pre(i,k)-temp_ibm2

             force_points(i,j,1)=force_points(i,j,1)+
     &              cs_ibm_target*k_link*temp_ibm1*ksh_tmp
             force_points(i,j,2)=force_points(i,j,2)+
     &              cs_ibm_target*k_link*temp_ibm2*ksh_tmp

             do j2=1,ns_ibm_r(i)
                 if ((k .le. 1) .or. (k .gt. ns_ibm_r(i))) then
                          fn1=0.d0
			              if (k .eq. j2) then
				                fn1=1.d0/ds_ibm(i,j2)
                    endif
                else
                    dkj2=k-j2
                    temp3_ibm= abs(dkj2/(1.d0))
                    call delta_fun(delta_typex,temp3_ibm,fn1)
                    fn1=fn1/ds_ibm(i,j2)
                 endif
                 f_link1(i,j2)= f_link1(i,j2)+
     &                               fn1*force_points(i,j,1)
                 f_link2(i,j2)= f_link2(i,j2)+
     &                               fn1*force_points(i,j,2)
             end do 
         enddo
         tdx0_1(i)= tdx0_1(i)+f_link1(i,1)*ds_ibm(i,1)
         tdx0_2(i)= tdx0_2(i)+f_link2(i,1)*ds_ibm(i,1)
         tdxnm1_1(i)= tdxnm1_1(i)-f_link1(i,ns_ibm_r(i))           
     &    *ds_ibm(i,ns_ibm_r(i))
         tdxnm1_2(i)= tdxnm1_2(i)-f_link2(i,ns_ibm_r(i))           
     &    *ds_ibm(i,ns_ibm_r(i))
      end do


c calculation of tension forces to impose inextensibility
      call cal_tension(temp1,tdx0_1,tdx0_2,tdxnm1_1,tdxnm1_2)
      if(implicitflag) then
          kb_ibmt=0.0
      else
          kb_ibmt=1.0
      endif

      if(implicitflagvi) then
          kb_ibmtvi=0.0
      else
          kb_ibmtvi=1.0
      endif

      do i=1,nr_ibm_f
        density_coefp= density_coefp_r(i) 
       do j=1,ns_ibm_r(i) 
        cs_ibm= cs_ibm_r(i) 
        density_coef= density_coef_rs(i,j) 
           if( (j.le.1) .or. (j .ge. ns_ibm_r(i)) )then
           fibm1(i,j)=(density_coef)*gx_ibm_massiveo(i,j)
     &        +density_coefp*gx_ibm_massiveo(i,j)
     &        +cs_ibm*gx_ibm(i,j)*temp4
     &        + temp1*(-fibm1(i,j)
     &                 + finer_1_ibm(i,j)  + f_impuls1(i,j) 
     &                 - 1.0*f_link1(i,j)*ds_ibm(i,1)+  
     &                  kb_ibmt  *fb_1_ibm(i,j)+
     &                  kb_ibmtvi*fbv_1_ibm(i,j)-
     &                  (1.0-kb_ibmtvi)*fbvi_1_ibm(i,j)+
     &                  density_coef*fr*cos(the_grav))
           fibm2(i,j)=(density_coef)*gy_ibm_massiveo(i,j)
     &        +density_coefp*gy_ibm_massiveo(i,j)
     &        +cs_ibm*gy_ibm(i,j)*temp4
     &        + temp1*(-fibm2(i,j) 
     &                 + finer_2_ibm(i,j) + f_impuls2(i,j)
     &                 - 1.0*f_link2(i,j)*ds_ibm(i,1)+ 
     &                  kb_ibmt*  fb_2_ibm(i,j)+
     &                  kb_ibmtvi*fbv_2_ibm(i,j)-
     &                  (1.0-kb_ibmtvi)*fbvi_2_ibm(i,j)+
     &                  density_coef*fr*sin(the_grav))
           else
           fibm1(i,j)=(density_coef)*gx_ibm_massiveo(i,j)
     &        +density_coefp*gx_ibm_massiveo(i,j)
     &        +cs_ibm*gx_ibm(i,j)*temp4
     &        + temp1*(-fibm1(i,j) - f_link1(i,j)
     &                 + finer_1_ibm(i,j) + f_impuls1(i,j)+
     &                  kb_ibmt*  fb_1_ibm(i,j)+
     &                  kb_ibmtvi*fbv_1_ibm(i,j)-
     &                  (1.0-kb_ibmtvi)*fbvi_1_ibm(i,j)+
     &                  density_coef*fr*cos(the_grav))
           fibm2(i,j)=(density_coef)*gy_ibm_massiveo(i,j)
     &        +density_coefp*gy_ibm_massiveo(i,j)
     &        +cs_ibm*gy_ibm(i,j)*temp4
     &        + temp1*(-fibm2(i,j) - f_link2(i,j)
     &                 + finer_2_ibm(i,j) + f_impuls2(i,j) +
     &                  kb_ibmt*  fb_2_ibm(i,j)+
     &                  kb_ibmtvi*fbv_2_ibm(i,j)-
     &                  (1.0-kb_ibmtvi)*fbvi_2_ibm(i,j)+
     &                  density_coef*fr*sin(the_grav))
            endif          

       end do
      end do

c prepare tridiagnal matrix for solving for x*
      do i=1,nr_ibm_f
!       write(*,*)  'kb', kb_ibm
        density_coefp= density_coefp_r(i) 
       a=0.0
       do j=1,ns_ibm_r(i)
         density_coef= density_coef_rs(i,j) 
         cs_ibm= cs_ibm_r(i) 
         kb_ibm=kb_ibm_rs(i,j)
         if(implicitflag) then
          kb_ibmt=1.0
         else
          kb_ibmt=0.0
         endif

         cb_ibm=cb_ibm_rs(i,j)
         if(implicitflagvi) then
          kb_ibmtvi=1.0
         else
          kb_ibmtvi=0.0
         endif


         if (j .eq. 1) then
          if(free_bc(i)) then
           j2=j+2
            mat(1,1)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2-1))
            mat(1,2)=-1.0/ds_ibm(i,j2)*(1.0/dsf_ibm(i,j2-1)
     &        +1.0/dsf_ibm(i,j2))
            mat(1,3)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2)) 
           j2=j+1
            mat(2,1)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2-1))
            mat(2,2)=-1.0/ds_ibm(i,j2)*(1.0/dsf_ibm(i,j2-1)
     &        +1.0/dsf_ibm(i,j2))
            mat(2,3)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2)) 
            temp3=fs_2_ibm(i,j)/dsf_ibm(i,j)/(dsf_ibm(i,j)*0.5d0)
            a(j,-1)=0.d0
            a(j,0)=(density_coef+density_coefp+cs_ibm*temp1**0.5)+temp3  
     &             +temp1*kb_ibm*(-mat(2,1))/dsf_ibm(i,j+1)
     &        *kb_ibmt
     &       /dsf_ibm(i,j)
     &             +temp1*cb_ibm/temp4*(-mat(2,1))/dsf_ibm(i,j+1)
     &        *kb_ibmtvi
     &       /dsf_ibm(i,j)

            a(j,1)=-temp3  +temp1*kb_ibm*(mat(1,1)-mat(2,2))
     &        *kb_ibmt
     &        /dsf_ibm(i,j+1)/dsf_ibm(i,j)
     &        +temp1*cb_ibm/temp4*(mat(1,1)-mat(2,2))
     &        *kb_ibmtvi
     &        /dsf_ibm(i,j+1)/dsf_ibm(i,j)

            a(j,2)= temp1*kb_ibm*(mat(1,2)-mat(2,3))
     &        *kb_ibmt
     &        /dsf_ibm(i,j+1)/dsf_ibm(i,j)
     &        +temp1*cb_ibm/temp4*(mat(1,2)-mat(2,3))
     &        *kb_ibmtvi
     &        /dsf_ibm(i,j+1)/dsf_ibm(i,j)

            a(j,3)= temp1*kb_ibm*(mat(1,3))
     &        *kb_ibmt
     &        /dsf_ibm(i,j+1)/dsf_ibm(i,j)
     &        + temp1*cb_ibm/temp4*(mat(1,3))
     &        *kb_ibmtvi
     &        /dsf_ibm(i,j+1)/dsf_ibm(i,j)

            fibm1(i,j)=fibm1(i,j)-temp1*2.0/dsf_ibm(i,j)*tdx0_1(i)
            fibm2(i,j)=fibm2(i,j)-temp1*2.0/dsf_ibm(i,j)*tdx0_2(i)
          else
            a(j,-1)=0.0d0
            a(j,0)=1.0d0/(dsf_ibm(i,3)*dsf_ibm(i,3))
            a(j,1)=0.0d0
          endif
        elseif(j .eq. 2) then
          if(fix_bc(i)) then
            a(j,-1)=0.0d0
            a(j,0)=1.0d0/(dsf_ibm(i,3)*dsf_ibm(i,3))
            a(j,1)=0.0d0
          else
          temp2=fs_2_ibm(i,j)/dsf_ibm(i,j)/ds_ibm(i,j)
          temp3=fs_2_ibm(i,j-1)/dsf_ibm(i,j-1)/ds_ibm(i,j)
         a(j,-1)=-2.0/dsf_ibm(i,j)**4*temp1*kb_ibm 
     &        *kb_ibmt+
     &           (-2.0)/dsf_ibm(i,j)**4*temp1*cb_ibm/temp4 
     &        *kb_ibmtvi
     &        -temp3
         a(j,0 )= 5.0/dsf_ibm(i,j)**4*temp1*kb_ibm 
     &        *kb_ibmt+
     &         5.0/dsf_ibm(i,j)**4*temp1*cb_ibm/temp4 
     &        *kb_ibmtvi
     &          +(temp2+temp3)
     &         +(density_coef+density_coefp+cs_ibm*temp1**0.5)
         a( j,1)=-4.0/dsf_ibm(i,j)**4*temp1*kb_ibm
     &        *kb_ibmt+
     &         (-4.0)/dsf_ibm(i,j)**4*temp1*cb_ibm/temp4
     &        *kb_ibmtvi
     &          -temp2
         a( j,2)= 1.0/dsf_ibm(i,j)**4*temp1*kb_ibm 
     &        *kb_ibmt+
     &         1.0/dsf_ibm(i,j)**4*temp1*cb_ibm/temp4 
     &        *kb_ibmtvi
         endif
        elseif(j .eq. ns_ibm_r(i)-1) then
          temp2=fs_2_ibm(i,j)/dsf_ibm(i,j)/ds_ibm(i,j)
          temp3=fs_2_ibm(i,j-1)/dsf_ibm(i,j-1)/ds_ibm(i,j)
         a(j,-2)= 1.0/dsf_ibm(i,j)**4*temp1*kb_ibm
     &         *kb_ibmt+
     &         1.0/dsf_ibm(i,j)**4*temp1*cb_ibm/temp4
     &         *kb_ibmtvi
         a(j,-1)=-4.0/dsf_ibm(i,j)**4*temp1*kb_ibm
     &         *kb_ibmt+
     &         (-4.0)/dsf_ibm(i,j)**4*temp1*cb_ibm/temp4
     &         *kb_ibmtvi
     &    -temp3
         a(j, 0)= 5.0/dsf_ibm(i,j)**4*temp1*kb_ibm 
     &         *kb_ibmt+
     &         (5.0)/dsf_ibm(i,j)**4*temp1*cb_ibm/temp4 
     &         *kb_ibmtvi
     &          +(temp2+temp3)
     &         +(density_coef+density_coefp+cs_ibm*temp1**0.5)
         a(j, 1)=-2.0/dsf_ibm(i,j)**4*temp1*kb_ibm
     &         *kb_ibmt+
     &         (-2.0)/dsf_ibm(i,j)**4*temp1*cb_ibm/temp4
     &         *kb_ibmtvi
     &         -temp2


         elseif (j .eq. ns_ibm_r(i)) then
           j2=j-2
            mat(1,1)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2-1))
            mat(1,2)=-1.0/ds_ibm(i,j2)*(1.0/dsf_ibm(i,j2-1)
     &        +1.0/dsf_ibm(i,j2))
            mat(1,3)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2)) 
           j2=j-1
            mat(2,1)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2-1))
            mat(2,2)=-1.0/ds_ibm(i,j2)*(1.0/dsf_ibm(i,j2-1)
     &        +1.0/dsf_ibm(i,j2))
            mat(2,3)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2))

            temp3=fs_2_ibm(i,j-1)/dsf_ibm(i,j-1)/(dsf_ibm(i,j-1)*0.5d0)
            a(j,-3)= temp1*kb_ibm*(mat(1,3))
     &         *kb_ibmt
     &        /dsf_ibm(i,j-1)/dsf_ibm(i,j)+
     &         temp1*cb_ibm/temp4*(mat(1,3))
     &         *kb_ibmtvi
     &        /dsf_ibm(i,j-1)/dsf_ibm(i,j)

            a(j,-2)= temp1*kb_ibm*(mat(1,2)-mat(2,3))
     &         *kb_ibmt
     &        /dsf_ibm(i,j-1)/dsf_ibm(i,j)+
     &         temp1*cb_ibm/temp4*(mat(1,2)-mat(2,3))
     &         *kb_ibmtvi
     &        /dsf_ibm(i,j-1)/dsf_ibm(i,j)

            a(j,-1)=-temp3 +
     &         temp1*kb_ibm*(mat(1,1)-mat(2,2))
     &         *kb_ibmt
     &        /dsf_ibm(i,j-1)/dsf_ibm(i,j)+
     &         temp1*cb_ibm/temp4*(mat(1,1)-mat(2,2))
     &         *kb_ibmtvi
     &        /dsf_ibm(i,j-1)/dsf_ibm(i,j)
            a(j,0)=(density_coef+density_coefp+cs_ibm*temp1**0.5)
     &            +temp3
     &            -temp1*kb_ibm
     &            *(mat(2,1))/dsf_ibm(i,j-1)/dsf_ibm(i,j)
     &         *kb_ibmt
     &            -temp1*cb_ibm/temp4
     &            *(mat(2,1))/dsf_ibm(i,j-1)/dsf_ibm(i,j)
     &         *kb_ibmtvi

            a(j,1)=0.0d0
            fibm1(i,j)=fibm1(i,j)+temp1*2.0/dsf_ibm(i,j-1)*tdxnm1_1(i)
            fibm2(i,j)=fibm2(i,j)+temp1*2.0/dsf_ibm(i,j-1)*tdxnm1_2(i)
         else 
           j2=j+1
            mat(1,1)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2-1))
            mat(1,2)=-1.0/ds_ibm(i,j2)*(1.0/dsf_ibm(i,j2-1)
     &        +1.0/dsf_ibm(i,j2))
            mat(1,3)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2)) 

           j2=j
            mat(2,1)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2-1))
            mat(2,2)=-1.0/ds_ibm(i,j2)*(1.0/dsf_ibm(i,j2-1)
     &        +1.0/dsf_ibm(i,j2))
            mat(2,3)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2))  
           j2=j-1
            mat(3,1)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2-1))
            mat(3,2)=-1.0/ds_ibm(i,j2)*(1.0/dsf_ibm(i,j2-1)
     &        +1.0/dsf_ibm(i,j2))
            mat(3,3)=1.0/(ds_ibm(i,j2)*dsf_ibm(i,j2)) 
          temp2=fs_2_ibm(i,j)/dsf_ibm(i,j)/ds_ibm(i,j)
          temp3=fs_2_ibm(i,j-1)/dsf_ibm(i,j-1)/ds_ibm(i,j)
          a(j,-2)= temp1*kb_ibm_rs(i,j-1)
     &         *kb_ibmt
     &        *(mat(3,1)/dsf_ibm(i,j-1))/dsf_ibm(i,j)+
     &         temp1*cb_ibm_rs(i,j-1)/temp4
     &         *kb_ibmtvi
     &        *(mat(3,1)/dsf_ibm(i,j-1))/dsf_ibm(i,j)

          a(j,-1)=-temp3 + temp1
     &        *(kb_ibm_rs(i,j-1)*mat(3,2)/dsf_ibm(i,j-1)
     &        -kb_ibm_rs(i,j)*mat(2,1)/dsf_ibm(i,j-1)
     &        -kb_ibm_rs(i,j)*mat(2,1)/dsf_ibm(i,j))
     &        /dsf_ibm(i,j)
     &         *kb_ibmt
     &         + temp1
     &        *(cb_ibm_rs(i,j-1)*mat(3,2)/dsf_ibm(i,j-1)
     &        -cb_ibm_rs(i,j)*mat(2,1)/dsf_ibm(i,j-1)
     &        -cb_ibm_rs(i,j)*mat(2,1)/dsf_ibm(i,j))/temp4
     &        /dsf_ibm(i,j)
     &         *kb_ibmtvi

          a(j,0)=(density_coef+density_coefp+cs_ibm*temp1**0.5)
     &         +(temp2+temp3)
     &         +temp1*(kb_ibm_rs(i,j+1)
     &        *mat(1,1)/dsf_ibm(i,j+1)+kb_ibm_rs(i,j-1)*mat(3,3)
     &        /dsf_ibm(i,j-1)-kb_ibm_rs(i,j)*mat(2,2)/dsf_ibm(i,j-1)
     &        -kb_ibm_rs(i,j)*mat(2,2)/dsf_ibm(i,j))/dsf_ibm(i,j)
     &         *kb_ibmt
     &         +temp1*(cb_ibm_rs(i,j+1)
     &        *mat(1,1)/dsf_ibm(i,j+1)+cb_ibm_rs(i,j-1)*mat(3,3)
     &        /dsf_ibm(i,j-1)-cb_ibm_rs(i,j)*mat(2,2)/dsf_ibm(i,j-1)
     &        -cb_ibm_rs(i,j)*mat(2,2)/dsf_ibm(i,j))/dsf_ibm(i,j)/temp4
     &         *kb_ibmtvi


          a(j,1)=-temp2 + temp1*
     &        (kb_ibm_rs(i,j+1)*mat(1,2)/dsf_ibm(i,j+1)
     &        -kb_ibm_rs(i,j)*mat(2,3)/dsf_ibm(i,j-1)
     &        -kb_ibm_rs(i,j)*mat(2,3)/dsf_ibm(i,j))
     &        /dsf_ibm(i,j)
     &         *kb_ibmt
     &        + temp1*
     &        (cb_ibm_rs(i,j+1)*mat(1,2)/dsf_ibm(i,j+1)
     &        -cb_ibm_rs(i,j)*mat(2,3)/dsf_ibm(i,j-1)
     &        -cb_ibm_rs(i,j)*mat(2,3)/dsf_ibm(i,j))/temp4
     &        /dsf_ibm(i,j)
     &         *kb_ibmtvi

          a(j,2)= temp1*kb_ibm_rs(i,j+1)
     &        *(mat(1,3)/dsf_ibm(i,j+1))/dsf_ibm(i,j)
     &         *kb_ibmt
     &          + temp1*cb_ibm_rs(i,j+1)/temp4
     &        *(mat(1,3)/dsf_ibm(i,j+1))/dsf_ibm(i,j)
     &         *kb_ibmtvi
!       write(*,*),'x',temp1,kb_ibm*mat(3,1)
          end if

        end do
        do j=1,ns_ibm_r(i) 
          do iksh1=-3,3
           a(j,iksh1)=dsf_ibm(i,3)*dsf_ibm(i,3)*a(j,iksh1)
          enddo
          fibm1(i,j)=dsf_ibm(i,3)*dsf_ibm(i,3)*fibm1(i,j)
          fibm2(i,j)=dsf_ibm(i,3)*dsf_ibm(i,3)*fibm2(i,j)
        enddo

        if(.not.(free_bc(i))) then
          if(fix_bc(i)) then
             fibm1(i,1)=target_points(i,1,1)        
             fibm2(i,1)=target_points(i,1,2)        

             fibm1(i,2)=target_points(i,2,1)        
             fibm2(i,2)=target_points(i,2,2)        
           else
             fibm1(i,1)=target_points(i,1,1)        
             fibm2(i,1)=target_points(i,1,2)        
           endif
          endif

       call bicgstab            
     &    (3,3,ns_ibm_r(i),a(1:ns_ibm_r(i),-3:3),fibm1(i,1:ns_ibm_r(i))
     &    ,gx_ibm_massive(i,1:ns_ibm_r(i))           
     &    ,gx_ibm_massive(i,1:ns_ibm_r(i)),
     &       itmax,tol,precond,success,itconverge)
         
!        print*,'x', success,itconverge
       iksh1=abs(success)
       call bicgstab           
     &     (3,3,ns_ibm_r(i),a(1:ns_ibm_r(i),-3:3),fibm2(i,1:ns_ibm_r(i))
     &      ,gy_ibm_massive(i,1:ns_ibm_r(i))          
     &      ,gy_ibm_massive(i,1:ns_ibm_r(i)),
     &       itmax,tol,precond,success,itconverge)


!        print*,'x', success,itconverge
!        success=-1



        iksh1=0     !iksh1+abs(success)
       if(iksh1 .ne. 0) then
        print*,'failed iterative'
       aamat=0.0d0
       bbmat=0.0d0

       do iksh1=1,ns_ibm_r(i)
             do iksh2= -3,3
                iksh3=iksh2+iksh1
                aamat(iksh1,iksh3)=a(iksh1,iksh2)


!                 write(1210,*) iksh1,iksh3,aamat(iksh1,iksh3)
             enddo
             bbmat(iksh1,1)=fibm1(i,iksh1)
             bbmat(iksh1,2)=fibm2(i,iksh1)
       enddo
!		pause

       call gaussj(aamat,ns_ibm_r(i),ns_ibm_r(i),bbmat,2,2,ierrksh1)
       do iksh1=1,ns_ibm_r(i)
             gx_ibm_massive(i,iksh1)=bbmat(iksh1,1)
             gy_ibm_massive(i,iksh1)=bbmat(iksh1,2)
       enddo
       endif


      end do
440   continue
      do i=1,nr_ibm
       do j=1,ns_ibm_r(i)
          vibm1(i,j)=(gx_ibm_massive(i,j)-gx_ibm(i,j))/temp4
          vibm2(i,j)=(gy_ibm_massive(i,j)-gy_ibm(i,j))/temp4
       end do
      end do

!      fibm1=0.0;fibm2=0.0
      do i=1,nr_ibm
        do j=1,ns_ibm
           if(j .le. ns_ibm_r(i)) then
           fibm1(i,j)=ffluidsum1(i,j)
           fibm2(i,j)=ffluidsum2(i,j)
	   endif
        end do
      end do
       
      end do
200   format(24(f20.9,'  '))
240   format(2(i5,'  '),24(f20.9,'  '))
 ! t_temp
      return 
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine cal_tension(dt2,tdx0_1,tdx0_2,tdxnm1_1,tdxnm1_2)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      integer i,j
      real*8 temp1,temp2,temp3,dt2,rhs(nr_ibm,ns_ibm)
      real*8 temp4_1_ibm(nr_ibm,0:ns_ibm),temp4_2_ibm(nr_ibm,0:ns_ibm)
     &      ,a(ns_ibm)
     &      ,b(ns_ibm),c(ns_ibm)
     &      ,tdx0_1(nr_ibm),tdx0_2(nr_ibm)
     &      ,tdxnm1_1(nr_ibm),tdxnm1_2(nr_ibm)
      temp3 = dt2  !here dt acually is delta_t^2 in input to subroutine
c here {dt^2} bring to t^{n+1/2} to have a good matrix shape
      do i=1,nr_ibm_f
        density_coefp= density_coefp_r(i) 
         do j=1,ns_ibm_r(i)-1
          density_coef= density_coef_rs(i,j) 
          cs_ibm= cs_ibm_r(i) 
c adding (dx.dx)^{n-1} to rhs  rhs <-- \frac{(1+(dx.dx)^{n-1})}
            rhs(i,j)=(density_coef+0.0*density_coefp)
     &      *(1.0d0+fs_1_ibmo(i,j))/2.0d0+density_coefp/2.0d0
c computing (dx.dx)^{n}

            temp1=(gx_ibm(i,j+1)-gx_ibm(i,j))/dsf_ibm(i,j)
            temp2=(gy_ibm(i,j+1)-gy_ibm(i,j))/dsf_ibm(i,j) 
c saving computed term for using later time
            fs_1_ibm(i,j)=temp1**2.+temp2**2.
c            fs_2_ibm(i,j)=ks_ibm*(fs_1_ibm(i,j)**0.5-1.0)*temp3
c rhs <-- rhs - \frac{ (dx.dx)^{n} }{dt^2}
            rhs(i,j)=rhs(i,j)-(density_coef+0.0*density_coefp)
     &      *fs_1_ibm(i,j)
c rhs <-- rhs + cs_ibm/2 \frac{ 1 - (dx.dx)^{n} }{dt} 
           rhs(i,j)=rhs(i,j)+cs_ibm/2.d0*(1.d0-fs_1_ibmo(i,j))
     &      *temp3**0.5
c rhs <---rhs - \(dv.dv)^{n}
            temp1=(vibm1(i,j+1)-vibm1(i,j))/dsf_ibm(i,j)
             temp2=(vibm2(i,j+1)-vibm2(i,j))/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)-(density_coef+0.0*density_coefp)
     &      *(temp1**2.+temp2**2.)*temp3
c computing (dx)^{*}
            temp4_1_ibm(i,j)=
     &          (gx_ibm_massive(i,j+1)-gx_ibm_massive(i,j))/dsf_ibm(i,j)
            temp4_2_ibm(i,j)=
     &          (gy_ibm_massive(i,j+1)-gy_ibm_massive(i,j))/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)-density_coefp
     &              *(temp4_1_ibm(i,j)**2+temp4_2_ibm(i,j)**2)/2.0d0

         end do
c      return
         do j=1,ns_ibm_r(i)-1
            if ((j .eq. 1).and. .not.(free_bc(i))) then
            temp2=0.0d0
            temp1=fbv_1_ibm(i,j+1) +fb_1_ibm(i,j+1) 
     &           +  f_impuls1(i,j+1)
     &           + 0.0*finer_1_ibm(i,j+1)-f_link1(i,j+1)*ds_ibm(i,2)
     &              -fibm1(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)-temp4_1_ibm(i,j)*temp2*temp3
            temp2=0.0d0
            temp1=fbv_2_ibm(i,j+1) +fb_2_ibm(i,j+1) 
     &           +   f_impuls2(i,j+1)
     &           +0.0*finer_2_ibm(i,j+1)-f_link2(i,j+1)*ds_ibm(i,2)
     &              -fibm2(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)-temp4_2_ibm(i,j)*temp2*temp3
            else if ((j .eq. 1).and. free_bc(i)) then
            temp2=fbv_1_ibm(i,j) +fb_1_ibm(i,j) 
     &           +   f_impuls1(i,j)
     &        +0.0*finer_1_ibm(i,j)-f_link1(i,j)*ds_ibm(i,1)-fibm1(i,j)
            temp1=fbv_1_ibm(i,j+1)  +fb_1_ibm(i,j+1)  
     &           +  f_impuls1(i,j+1)
     &           +0.0*finer_1_ibm(i,j+1)-f_link1(i,j+1)*ds_ibm(i,2)
     &              -fibm1(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)-temp4_1_ibm(i,j)*temp2*temp3
     &              -temp4_1_ibm(i,j)*tdx0_1(i)*temp3*2./dsf_ibm(i,j)**2
            temp2=fbv_2_ibm(i,j) +fb_2_ibm(i,j)  
     &           +  f_impuls2(i,j)
     &        +0.0*finer_2_ibm(i,j)-f_link2(i,j)*ds_ibm(i,1)-fibm2(i,j)
            temp1=fbv_2_ibm(i,j+1)  +fb_2_ibm(i,j+1)  
     &           +   f_impuls2(i,j+1)
     &           +0.0*finer_2_ibm(i,j+1)-f_link2(i,j+1)*ds_ibm(i,2)
     &              -fibm2(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)-temp4_2_ibm(i,j)*temp2*temp3
     &              -temp4_2_ibm(i,j)*tdx0_2(i)*temp3*2./dsf_ibm(i,j)**2
            elseif ((j .eq. ns_ibm_r(i)-1)) then  !check
            temp2=fbv_1_ibm(i,j) + fb_1_ibm(i,j)  
     &           +  f_impuls1(i,j)
     &      +0.0*finer_1_ibm(i,j)-f_link1(i,j)*ds_ibm(i,j)-fibm1(i,j)
            temp1=fbv_1_ibm(i,j+1)  +fb_1_ibm(i,j+1)  
     &           +   f_impuls1(i,j+1)
     &        +0.0*finer_1_ibm(i,j+1)-f_link1(i,j+1)*ds_ibm(i,j+1)
     &              -fibm1(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)
     &              -temp4_1_ibm(i,j)*temp2*temp3
     &            -temp4_1_ibm(i,j)*tdxnm1_1(i)*temp3*2./dsf_ibm(i,j)**2
            temp2=fbv_2_ibm(i,j) +fb_2_ibm(i,j)  
     &           +   f_impuls2(i,j)
     &       +0.0*finer_2_ibm(i,j)-f_link2(i,j)*ds_ibm(i,j)-fibm2(i,j)
            temp1=fbv_2_ibm(i,j+1) +fb_2_ibm(i,j+1)  
     &           +   f_impuls2(i,j+1)
     &       +0.0*finer_2_ibm(i,j+1)-f_link2(i,j+1)*ds_ibm(i,j+1)
     &              -fibm2(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)
     &              -temp4_2_ibm(i,j)*temp2*temp3
     &            -temp4_2_ibm(i,j)*tdxnm1_2(i)*temp3*2./dsf_ibm(i,j)**2
            else
            temp2=fbv_1_ibm(i,j) +fb_1_ibm(i,j)  
     &           +  f_impuls1(i,j)
     &           +0.0*finer_1_ibm(i,j)-f_link1(i,j)-fibm1(i,j)
            temp1=fbv_1_ibm(i,j+1) +fb_1_ibm(i,j+1)  
     &           +   f_impuls1(i,j+1)
     &           +0.0*finer_1_ibm(i,j+1)-f_link1(i,j+1)-fibm1(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)-temp4_1_ibm(i,j)*temp2*temp3
            temp2=fbv_2_ibm(i,j) +fb_2_ibm(i,j)  
     &           +   f_impuls2(i,j)
     &           +0.0*finer_2_ibm(i,j)-f_link2(i,j)-fibm2(i,j)
            temp1=fbv_2_ibm(i,j+1) +fb_2_ibm(i,j+1)  
     &           +   f_impuls2(i,j+1)
     &           +0.0*finer_2_ibm(i,j+1)-f_link2(i,j+1)-fibm2(i,j+1)
            temp2=(temp1-temp2)/dsf_ibm(i,j)
            rhs(i,j)=rhs(i,j)-temp4_2_ibm(i,j)*temp2*temp3
            end if
         end do
      end do
      do i=1,nr_ibm_f
       do j=1,ns_ibm_r(i)-1
         temp2=1.0d0/dsf_ibm(i,j)/ds_ibm(i,j)
         temp3=1.0d0/dsf_ibm(i,j)/ds_ibm(i,j+1)
         if ((j .eq. 1) .and. .not.(free_bc(i))) then
          c(j)=0.0d0
          a(j)=-temp3*( temp4_1_ibm(i,j)**2.0
     &          +                temp4_2_ibm(i,j)**2.0 )
          b(j)=temp3*( temp4_1_ibm(i,j+1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j+1)*temp4_2_ibm(i,j) )
         elseif ((j .eq. 1) .and. free_bc(i)) then
          temp1=2.0d0/dsf_ibm(i,j)/dsf_ibm(i,j)
          c(j)=0.0d0
          a(j)=-(temp3+1.0d0*temp1)*( temp4_1_ibm(i,j)**2.0
     &          +                temp4_2_ibm(i,j)**2.0 )
          b(j)=temp3*( temp4_1_ibm(i,j+1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j+1)*temp4_2_ibm(i,j) )
         else if (j .eq. ns_ibm_r(i)-1) then
          temp1=2.0d0/dsf_ibm(i,j)/dsf_ibm(i,j)
          c(j)=temp2*( temp4_1_ibm(i,j-1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j-1)*temp4_2_ibm(i,j) )
          a(j)=-(temp2+1.0d0*temp1)*( temp4_1_ibm(i,j)**2.0
     &          +                temp4_2_ibm(i,j)**2.0 )
          b(j)=0.0d0
         else 
          c(j)=temp2*( temp4_1_ibm(i,j-1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j-1)*temp4_2_ibm(i,j) )
          a(j)=-(temp2+temp3)*( temp4_1_ibm(i,j)**2.0
     &          +                temp4_2_ibm(i,j)**2.0 )
          b(j)=temp3*( temp4_1_ibm(i,j+1)*temp4_1_ibm(i,j)
     &          +         temp4_2_ibm(i,j+1)*temp4_2_ibm(i,j) )
         end if
       end do
       call thomas(ns_ibm_r(i)-1,a,b,c,rhs(i,1:ns_ibm_r(i)-1)
     &            ,fs_2_ibm(i,1:ns_ibm_r(i)-1),ns_ibm_r(i)-1)
      end do
      return
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine pentdag(a,b,c,d,e,f,u,n) 
      implicit none
      save

c..solves for a vector u of length n in the pentadiagonal linear system 
c.. a_i u_(i-2) + b_i u_(i-1) + c_i u_i + d_i u_(i+1) + e_i u_(i+2) = f_i 
c..input are the a, b, c, d, e, and f and they are not modified 

c..in its clearest incarnation, this algorithm uses three storage arrays 
c..called p, q and r. here, the solution vector u is used for r, cutting 
c..the extra storage down to two arrays. 

c..declare the pass
      integer          n
      double precision a(n),b(n),c(n),d(n),e(n),f(n),u(n)

c..local variables
      integer          nmax,i 

      integer        i2

      parameter        (nmax=500) 
      double precision p(nmax),q(nmax),bet,den 


!      do i2=1,n
!         write(2234,244) a(i2),b(i2),c(i2),d(i2),e(i2),f(i2),u(i2)
!      enddo
!      pause
244    format(24(f20.9,'  '))


c..initialize elimination and backsubstitution arrays 
      if (c(1) .eq. 0.0)  stop 'eliminate u2 trivially' 
      bet  = 1.0d0/c(1) 
      p(1) = -d(1) * bet 
      q(1) = -e(1) * bet 
      u(1) = f(1)  * bet 

      bet = c(2) + b(2)*p(1) 
      if (bet .eq. 0.0) stop 'singular 1 in pentdag' 
      bet = -1.0d0/bet 
      p(2) = (d(2) + b(2)*q(1)) * bet 
      q(2) = e(2) * bet 
      u(2) = (b(2)*u(1) - f(2)) * bet 


c..reduce to upper triangular 
      do i=3,n 
       bet = b(i) + a(i) * p(i-2) 
       den = c(i) + a(i)*q(i-2) + bet*p(i-1) 
       if (den .eq. 0.0) stop 'singular 2 in pentdag' 
       den = -1.0d0/den 
       p(i) = (d(i) + bet*q(i-1)) * den 
       q(i) = e(i) * den 
       u(i) = (a(i)*u(i-2) + bet*u(i-1) - f(i)) * den 
      enddo

c..backsubstitution 
      u(n-1) = u(n-1) + p(n-1) * u(n) 
      do i=n-2,1,-1 
       u(i) = u(i) + p(i) * u(i+1) + q(i) * u(i+2) 
      enddo
      return 
      end 




      subroutine gaussj(a,n,np,b,m,mp,ierr)

c  purpose: solution of the system of linear equations ax = b by
c     gauss-jordan elimination, where a is a matrix of order n and b is
c     an n x m matrix.  on output a is replaced by its matrix inverse
c     and b is preplaced by the corresponding set of solution vectors.

c  source: w.h. press et al, "numerical recipes," 1989, p. 28.

c  modifications: 
c     1. double  precision.
c     2. error parameter ierr included.  0 = no error. 1 = singular 
c        matrix encountered; no inverse is returned.

c  prepared by j. applequist, 8/17/91.

      implicit real*8(a-h,o-z)

c        set largest anticipated value of n.

      parameter (nmax=500)
      dimension a(-3:np+3,-3:np+3),b(np,mp),
     &          ipiv(nmax),indxr(nmax),indxc(nmax)
      real*8 dabs_a

      ierr=0
      do 11 j=1,n
      ipiv(j)=0
 11   continue
      do 22 i=1,n
      big=0.d0
      do 13 j=1,n
      if (ipiv(j).ne.1) then
      do 12 k=1,n
      if (ipiv(k).eq.0) then

      dabs_a=a(j,k)
      if (dabs_a.lt.0.0) then
       dabs_a=-dabs_a
      endif

      if (dabs_a.ge.big) then
      big=dabs_a
      irow=j
      icol=k
      endif
      else if (ipiv(k).gt.1) then
      ierr=1
      return
      endif
 12   continue
      endif
 13   continue
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
      do 14 l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
 14   continue
      do 15 l=1,m
      dum=b(irow,l)
      b(irow,l)=b(icol,l)
      b(icol,l)=dum
 15   continue
      endif
      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.0.d0) then
      ierr=1
      return
      endif
      pivinv=1.d0/a(icol,icol)
      a(icol,icol)=1.d0
      do 16 l=1,n
      a(icol,l)=a(icol,l)*pivinv
 16   continue
      do 17 l=1,m
      b(icol,l)=b(icol,l)*pivinv
 17   continue
      do 21 ll=1,n
      if (ll.ne.icol) then
      dum=a(ll,icol)
      a(ll,icol)=0.d0
      do 18 l=1,n
      a(ll,l)=a(ll,l)-a(icol,l)*dum
 18   continue
      do 19 l=1,m
      b(ll,l)=b(ll,l)-b(icol,l)*dum
 19   continue
      endif
 21   continue
 22   continue
      do 24 l=n,1,-1
      if (indxr(l).ne.indxc(l)) then
      do 23 k=1,n
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l))
      a(k,indxc(l))=dum
 23   continue
      endif
 24   continue
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      FUNCTION S_Alpha(i,ns,fix_bc)
      integer istartS,iendS,iraiseS,ns,i
      parameter(istartS=3,iraiseS=6)
      real*8 S_function
      logical fix_bc
      real*8 denom,term1,term2


      iendS=ns+1-istartS-iraiseS
      denom=iraiseS
      term1=i-iendS
      term2=i-istartS

      if(fix_bc) then
      S_Alpha=1.0
     &       -S_function(term1/denom)
      else
      S_Alpha=S_function(term2/denom)
     &       -S_function(term1/denom)
      endif
!       write(433,*) i,ns,S_function(dble(i-istartS)/dble(iraiseS))
!     &                    ,S_function(dble(i-iendS)/dble(iraiseS))

      return
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE DELTA_FUN(Delta_type,R,F)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      IMPLICIT NONE
      integer Delta_type 
      real*8 R,F,PI,FI_IB_6 
      PI = 4. * ATAN(1.0)
      F=0.0
      Select case (Delta_type)
        case(1)
           IF (R .LT. 2.D0) F=0.25D0*(1.0D0+DCOS(PI*R/2.D0))
        case(2)
           IF (R .LT. 1.D0) THEN
              F=0.125*(3.-2.*R+SQRT(1.+4.*R-4.*R**2.))
           ELSE IF (R .LT. 2.D0) THEN
              F=0.125*(5.-2.*R-SQRT(-7.+12.*R-4.*R**2.))
           END IF
        case(3)
           IF (R .LE. 1.D0) THEN
              F=FI_IB_6(R)
           ELSE IF (R .LE. 2.D0) THEN
              F=21./16.+7./12.*R-7./8.*R**2+1./6.*R**3-1.5*FI_IB_6(R-1.)
           ELSE IF (R .LE. 3.D0) THEN
              F=9./8.-23./12.*R+3./4.*R**2-1/12.*R**3+0.5*FI_IB_6(R-2.)   
           END IF
        case(4)
           IF (R .LT. 1.D0) THEN
              F=1.-0.5*R-R**2+0.5*R**3
           ELSE IF (R .LT. 2.D0) THEN
              F=1.-11./6.*R+R**2-1./6.*R**3
           END IF
      End Select
      return
      end 
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      real*8 Function FI_IB_6(R)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      IMPLICIT NONE
      real*8 R
      FI_IB_6=61./112.-11./42.*R-11./56.*R**2.+1./12.*R**3
     &         +SQRT(3.)/336.*SQRT(243.+1584.*R-748.*R**2-1560.*R**3
     &         +500.*R**4+336.*R**5-112.*R**6)
      return
      end
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      FUNCTION S_function(X)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Function for fringe region S(x)
C INPUTS: x value 
c OUTPUTS: value of function
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      REAL X
      IF (X .LE. 0.0) THEN 
         S_function=0.0
      ELSE IF (X .LT. 1.0) then
         S_function=1.0/(1.0+EXP(1.0/(X-1.0)+1.0/(X)))
      ELSE 
         S_function=1.0
      END IF
      RETURN
      END
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE Delta_width(Del_X,Del_Y,Del_Z
     &                ,Delta_typeX,Delta_typeY,Delta_typeZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c width of Delta fuction ( in grid) in each direction
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      IMPLICIT NONE
      integer Del_X,Del_Y,Del_Z,Delta_typeX,Delta_typeY,Delta_typeZ
      Select case (Delta_typeX)
        case(1)
           Del_X=3
        case(2)
           Del_X=3
        case(3)
           Del_X=4
        case(4)
           Del_X=3
      End Select
      Select case (Delta_typeY)
        case(1)
           Del_Y=3
        case(2)
           Del_Y=3
        case(3)
           Del_Y=4
        case(4)
           Del_Y=3
      End Select
      Select case (Delta_typeZ)
        case(1)
           Del_Z=3
        case(2)
           Del_Z=3
        case(3)
           Del_Z=4
        case(4)
           Del_Z=3
      End Select
      return
      end 
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	real*8 FUNCTION Dist(xa,xe)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c    Computes the distance between two points
c    with coordinates (xa,ya) and (xe,ye)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	IMPLICIT NONE
c Cartesian dimension
	INTEGER Cdim 
	PARAMETER(Cdim=2)      
c coords of point 1
	real*8 xa(Cdim)
c coords of point 2      
	real*8 xe(Cdim)      
	INTEGER   :: N
	real*8 :: SUMS
	SUMS= 0.0d0
	DO 10,  N=1,Cdim
         SUMS= SUMS + (xa(n)-xe(n))**2.d0
10	CONTINUE
	Dist= SQRT(SUMS)
	RETURN
	END 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine boundary_ibm
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c satisfying boundary condition for ibm
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      real*8 temp,period_temp,gama_temp
      real*8 a_transient
      integer i,j,i2,i_fixr,i_fixcoef
      real*8 kin_t,kin_amp

     &    ,kin_teta,kin_t_total,kin_time
     &    ,kin_ax,kin_ang,kin_omega,kin_flag
     &    ,kin_ax_v
      real*8 ttramp,tmpksh1,tmpksh2,kin_rx,kin_rx_v
     &    ,kin_ry,kin_ry_v
     &    ,kin_rt,kin_rt_v,kin_time2
      real*8 temp1,temp2,temp3,arot_temp,arot0_temp
      real*8 gridxtmpO,gridytmpO
      real*8 dflag



      pi=4.0d0*atan(1.d0)

      do i2=1,nr_ibm
      ttramp=tramp(i2)
      period_temp=1.0/frequency_fin(i2)
      gama_temp=1.0d0

      kin_t=period_temp
      kin_t_total=kin_t
      kin_omega=2.*(pi)/kin_t

      if (time .le. max(time_max_df,0.0005d0)) then
             temp1=0.d0
             temp2=0.d0
      else
         kin_time=mod(time-max(1.5*time_max_df,0.0005d0),kin_t_total)
         kin_time2=time-max(1.5*time_max_df,0.0005d0)
!---------------  change feb26-----------------------------
         temp1=1.d0-
     &   dexp(-(time-max(1.5*time_max_df,0.0005d0))/(ttramp*kin_t))
         temp2=1.d0/(ttramp*kin_t)
     &   *dexp(-(time-max(1.5*time_max_df,0.0005d0))/(ttramp*kin_t))
!-----------------------------------------------------------
      endif


      if (flexible_i(i2) == 1) then
      do i=1,target_num(i2)
!          print*, i, a_fin(i), phi_fin(i), flagfixed(i)
          kin_amp=a_fin(i2,i)
          kin_teta=phi_fin(i2,i)
          dflag=flagfixed(i2,i)
          kin_flag=dflag
         kin_ax=
     &      kin_amp*sin(kin_omega*kin_time+kin_teta)

         kin_ax_v=
     &      (kin_amp*kin_omega*
     &       cos(kin_omega*kin_time+kin_teta))
     &     *temp1
     &     +kin_ax*temp2
         kin_ax=kin_ax*temp1

           i_fix=target_point_num(i2,i)
           target_points(i2,i,1)=   0.0
 
           target_points_v(i2,i,1)= 0.0 


           target_points(i2,i,2)= kin_ax  

           target_points_v(i2,i,2)= kin_ax_v  
        enddo
      endif
!rigid body motions
         i=i2
         kin_rt=ampt(i)
     &     *sin(2.0*pi*freqt(i)*kin_time2+phit(i))
         kin_rt_v=ampt(i)*2.0*pi*freqt(i)
     &     *cos(2.0*pi*freqt(i)*kin_time2+phit(i))
     &     *temp1
     &     +kin_rt*temp2
         kin_rt=kin_rt*temp1


         kin_rx=ampx(i)
     &     *sin(2.0*pi*freqx(i)*kin_time2+phix(i))
         kin_rx_v=ampx(i)*2.0*pi*freqx(i)
     &     *cos(2.0*pi*freqx(i)*kin_time2+phix(i))
     &     *temp1
     &     +kin_rx*temp2
         kin_rx=kin_rx*temp1


         kin_ry=ampy(i)
     &     *sin(2.0*pi*freqy(i)*kin_time2+phiy(i))
         kin_ry_v=ampy(i)*2.0*pi*freqy(i)
     &     *cos(2.0*pi*freqy(i)*kin_time2+phiy(i))
     &     *temp1
     &     +kin_ry*temp2
         kin_ry=kin_ry*temp1




      if (flexible_i(i2) == 1) then
         do i=1,target_num(i2)

           tmpksh1=target_points(i2,i,2)
           tmpksh2=target_points_v(i2,i,2)


           i_fix=target_point_num(i2,i)
           target_points(i2,i,1)= (gx_bp(i2,i_fix)-gt0(i2,1))
     &     *cos(kin_rt)
     &     + gt0(i2,1) + kin_rx +
     &       (gy_bp(i2,i_fix)+tmpksh1-gt0(i2,2))*(-sin(kin_rt))

           target_points_v(i2,i,1)= (gx_bp(i2,i_fix)-gt0(i2,1))
     &     *(-sin(kin_rt))*kin_rt_v
     &      +kin_rx_v
     &      +(gy_bp(i2,i_fix)+tmpksh1-gt0(i2,2))*(-cos(kin_rt))
     &      *kin_rt_v
     &      -sin(kin_rt)*tmpksh2

           target_points(i2,i,2)= (gx_bp(i2,i_fix)-gt0(i2,1))
     &     *sin(kin_rt)
     &     + gt0(i2,2) +kin_ry+  
     &       (gy_bp(i2,i_fix)+tmpksh1-gt0(i2,2))*cos(kin_rt)

           target_points_v(i2,i,2)= (gx_bp(i2,i_fix)-gt0(i2,1))
     &     *cos(kin_rt)*kin_rt_v
     &      +kin_ry_v
     &      +(gy_bp(i2,i_fix)+tmpksh1-gt0(i2,2))*(-sin(kin_rt))
     &      *kin_rt_v
     &      +cos(kin_rt)*tmpksh2

          enddo
        else ! rigid bodies
         do i=1,NS_IBM_r(i2)

           tmpksh1=0.0   
           tmpksh2=0.0  



           i_fix=i
           gridxtmpO= (gx_bp(i2,i_fix)-gt0(i2,1))
     &     *cos(kin_rt)
     &     + gt0(i2,1) + kin_rx +
     &       (gy_bp(i2,i_fix)+tmpksh1-gt0(i2,2))*(-sin(kin_rt))


            !print*, target_points(i2,i,1),target_points_v(i2,i,1)
           gridytmpO= (gx_bp(i2,i_fix)-gt0(i2,1))
     &     *sin(kin_rt)
     &     + gt0(i2,2) +kin_ry+  
     &       (gy_bp(i2,i_fix)+tmpksh1-gt0(i2,2))*cos(kin_rt)


          gx_ibm_massive(i2,i_fix)=gridxtmpO              
          gy_ibm_massive(i2,i_fix)=gridytmpO             
            !print*, target_points(i2,i,2),target_points_v(i2,i,2)
 
      enddo

      endif
      enddo
!	close(1610)
620    format(1x,2(i5,2x),24(f8.4,2x))
      return
      end
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine thomas 
     +
     +  (N   ! matrix size
     +  ,a   ! diagonal 
     +  ,b   ! super-diagonal row
     +  ,c   ! sub-diagonal row      
     +  ,s   ! rhs
     +  ,x   ! solution 
     +  ,Nmax
     +  )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c Thomas algorithm for tridiagonal systems
c Coefficient matrix:
c
c  | a1 b1  0   0  ...  0   0    0    |
c  | c2 a2  b2  0  ...  0   0    0    |
c  | 0  c3  a3  b3 ...  0   0    0    |
c  | ..............................   |
c  | 0  0   0   0  ... cn-1 an-1 bn-1 |
c  | 0  0   0   0  ...  0   cn   an   |
c
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      Implicit real*8 (a-h,o-z)
      Integer Nmax
      Dimension a(Nmax),b(Nmax),c(Nmax),s(Nmax),x(Nmax)
      Dimension d(Nmax),y(Nmax)

      Parameter (tol=0.00000001D0)

c--------
c prepare
c--------
c      print*,'n', N;pause
c      print*,'a', a;pause
c      print*,'b', b;pause
c      print*,'c', c;pause
c      print*,'s', s;pause
c      print*,'x', x;pause

      Na = N-1

c------------------------------
c reduction to upper bidiagonal
c------------------------------

      d(1) = b(1)/a(1)
      y(1) = s(1)/a(1)

      Do i=1,Na
       i1 = i+1
       Den   = a(i1)-c(i1)*d(i)
       d(i1) = b(i1)/Den
       y(i1) = (s(i1)-c(i1)*y(i))/Den
      End Do

c------------------
c Back substitution
c------------------

      x(N) = y(N)

      Do i=Na,1,-1
        x(i)= y(i)-d(i)*x(i+1)
      End Do

c-----------------------
c Verification and alarm
c-----------------------

      Res = s(1)-a(1)*x(1)-b(1)*x(2)

      If(abs(Res).gt.tol) write (6,*) " thomas: alarm"

      Do i=2,Na
        Res = s(i)-c(i)*x(i-1)-a(i)*x(i)-b(i)*x(i+1)
        If(abs(Res).gt.tol) write (6,*) " thomas: alarm"
      End Do

      Res = s(N)-c(N)*x(N-1)-a(N)*x(N)

      If(abs(Res).gt.tol) write (6,*) " thomas: alarm"

c-----
c Done
c-----

 100  Format (1x,f15.10)



      Return
      End
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine thomas_pr (N,A,B,C,S,X,Nmax)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c Thomas algorithm for tridiagonal systems
c
c  T . x = s
c with a periodic condition
c
c Coefficient matrix:
c  | a1 b1  0   0  ...  0      0       c1   |
c  | c2 a2  b2  0  ...  0      0       0    |
c  | 0  c3  a3  b3 ...  0      0       0    |
c  | ....................................   |
c  | 0  0   0   0  ... c(n-1) a(n-1) b(n-1) |
c  | bn 0   0   0  ...  0      c(n)   a(n)  |
c
c------------------------------------------
	IMPLICIT NONE
	INTEGER N,Nmax,Na,Nb
	real*8  A(Nmax),B(Nmax),C(Nmax),S(Nmax),X(Nmax),tol
	real*8  D(N),Y(N),Den,Res,R1,R0,X0(Nmax)
	real*8  save1,saveNa
	INTEGER i,i1
      Parameter (tol=0.00000001D0)
c--------
c parameter
c--------
      Na = N-1
      Nb = N-2
      save1  = s(1)
      saveNa = s(Na)
c------------------------------
c First assume that x(N) = 0
c and solve the first N-1 equations
c neglecting the last column
c and the last row
c------------------------------
      x(N) = 0.0D0
c-- REGULAR THOMAS

c-----
c reduction to upper bidiagonal
c-----
      d(1) = b(1)/a(1)
      y(1) = s(1)/a(1)
      Do 11,i=1,Nb
       i1 = i+1
       Den   = a(i1)-c(i1)*d(i)
       d(i1) = b(i1)/Den
       y(i1) = (s(i1)-c(i1)*y(i))/Den
11	CONTINUE
c----
c Back substitution
c----
      x(Na) = y(Na)
      Do 12, i=Nb,1,-1
        x(i)= y(i)-d(i)*x(i+1)
12	CONTINUE
c-----
c compute the first residual:
c-----
      R0 = a(N)*x(N) + b(N)*x(1) + c(N)*x(Na) - s(N)
c-----
c save the solution
c-----
      Do 13, i=1,N
        x0(i) = x(i)
13	CONTINUE

c------------------------------
c Second, assume that x(N) = 1
c and solve the first N-1 equations
c with a modified RHS
c------------------------------

      x(N) = 1.0D0

      s(1)  = s(1)  - c(1)  * x(N)
      s(Na) = s(Na) - b(Na) * x(N)


c-- REGULAR THOMAS

c-----
c reduction to upper bidiagonal
c-----

      d(1) = b(1)/a(1)
      y(1) = s(1)/a(1)

      Do i=1,Nb
       i1 = i+1
       Den   = a(i1)-c(i1)*d(i)
       d(i1) = b(i1)/Den
       y(i1) = (s(i1)-c(i1)*y(i))/Den
      End Do

c-----
c back substitution
c-----

      x(Na) = y(Na)

      Do i=Nb,1,-1
        x(i)= y(i)-d(i)*x(i+1)
      End Do

c------
c compute the second residual:
c-----

      R1 = a(N)*x(N) + b(N)*x(1) + c(N)*x(Na) - s(N)

c----------------------------
c rectify the right-hand side
c----------------------------

      s(1)  = save1
      s(Na) = saveNa

c---------------------------------
c compute the correct value of x(N)
c---------------------------------

      x(N) = -R0/(R1-R0)

c----------------------------
c compose the solution vector
c----------------------------

      Do i=1,Na
       x(i) = (x(i)-x0(i)) * x(N) + x0(i)
      End Do

c-----------------------
c Verification and alarm
c-----------------------

      Res = s(1) - a(1)*x(1) -b(1)*x(2)-c(1)*x(N)

      If(abs(Res).gt.tol) write (6,*) " thomas_pr: alarm, 1",Res

      Do i=2,Na
        Res = s(i)-c(i)*x(i-1)-a(i)*x(i)-b(i)*x(i+1)
        If(abs(Res).gt.tol) write (6,*) " thomas_pr: alarm ",i,Res
      End Do

      Res = s(N)-c(N)*x(Na)-a(N)*x(N)-b(N)*x(1)

      If(abs(Res).gt.tol) write (6,*) " thomas_pr: alarm ",N,Res

c-----
c Done
c-----

 100  Format (1x,f15.10)

      Return
      End




	

