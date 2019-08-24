c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine compute_force(timewrite)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c find drag and lift forces in 2d and 3d
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      implicit none
      integer nfile,nfile2
      real*8 temp1,temp2,temp3,force_1(7),force_2(7)
      real*8 temp1_local(nr_ibm)
     &    ,temp2_local(nr_ibm)
     &    ,temp3_local(nr_ibm)

      real*8 temp1_localv(nr_ibm)
     &    ,temp2_localv(nr_ibm)
     &    ,temp3_localv(nr_ibm)

      real*8 temp1b_local(nr_ibm)
      real*8 temp2a_local(nr_ibm)
      real*8 temp2b_local(nr_ibm)
      real*8 temp2c_local(nr_ibm)
      real*8 temp2d_local(nr_ibm)
      real*8 tempPiezo_localp(nr_ibm)

      real*8 temp2a_localforce(nr_ibm)
      real*8 temp2b_localforce(nr_ibm)
      real*8 temp2c_localforce(nr_ibm)

      integer i2,isec,i2l,i2g,k2
      real*8 das

      real*8 coor1(3),coor2(3),dist
!     ,tempa(nr_ibm,0:ns_ibm)
      real*8 tempb,temp,tempc, temp1Piezo
      integer i,j,k,ii,jq,ibdyg
      logical timewrite
      character  num*3, filename*52
      character  num2*3, filename2*52
      real*8 temprecord(3)

      real*8 temp_ibm1,temp_ibm2,temp_ibm3,ksh_tmp
     &    ,target_timefn,time_alphaf   
      real*8  real_distibm(50),rtmpusrfn(50)
      integer int_distibm(50),itmpusrfn(50)


      temp1=0.0d0; temp2=0.0d0; temp3=0.0d0
      temp1_local=0.0d0; temp2_local=0.0d0

      if(timewrite) then
      do i=1,nr_ibm
         nfile=i*100+21
         write(num,'(i3.3)') i
        filename='forces_each-'//num(1:3)//'.txt'
        open(nfile,file = filename,form='formatted',
     *  action='write',position='append')
      enddo

c find the drag and lift forces by method f*dt = int(m*(u0j-u_j), omega)
      do i=1,nr_ibm
       temp=0.d0;tempb=0.d0; tempc=0.d0
       temp1=0.0d0;temp2=0.0d0; temp3=0.0d0
       temp1Piezo=0.0d0
        if(bodytype(i) ==fiber) then
         i2l=Iglbloc_fib(i)
          do j=1,ns_ibm_r_fib(i2l)
             das=ds_ibm_fib(i2l,j)
          do isec=1,nsec_IBM
            temp1=temp1
     &      - ffluidsum1_fib(i2l,isec,j)*das*(dsecy_IBM+dsecz_IBM)
            tempb=tempb
     &      - ffluidsum1_fib(i2l,isec,j)*das*(dsecy_IBM+dsecz_IBM)
     &        *vibm1_fib(i2l,j)
            temp2=temp2
     &      - ffluidsum2_fib(i2l,isec,j)*das*(dsecy_IBM+dsecz_IBM)
            temp=temp
     &      - ffluidsum2_fib(i2l,isec,j)*das*(dsecy_IBM+dsecz_IBM)
     &        *vibm2_fib(i2l,j)
            temp3=temp3
     &      - ffluidsum3_fib(i2l,isec,j)*das*(dsecy_IBM+dsecz_IBM)
            tempc=tempc
     &      - ffluidsum3_fib(i2l,isec,j)*das*(dsecy_IBM+dsecz_IBM)
     &        *vibm3_fib(i2l,j)
          enddo
         end do
        elseif(bodytype(i) ==shell) then
         i2l=Iglbloc_esh(i)
          do j=1,ns_ibm_r_esh(i2l)
          das=0.d0
          do i2=1,nposele_esh(i2l,j,1)
             das=das+area(i2l,nposele_esh(i2l,j,1+i2))
          end do
            temp1=temp1
     &      - ffluidsum1_esh(i2l,j)*das
            tempb=tempb
     &      - ffluidsum1_esh(i2l,j)*das
     &        *vibm1_esh(i2l,j)
            temp2=temp2
     &      - ffluidsum2_esh(i2l,j)*das
            temp=temp
     &      - ffluidsum2_esh(i2l,j)*das
     &        *vibm2_esh(i2l,j)
            temp3=temp3
     &      - ffluidsum3_esh(i2l,j)*das
            tempc=tempc
     &      - ffluidsum3_esh(i2l,j)*das
     &        *vibm3_esh(i2l,j)
          enddo
        elseif(bodytype(i) ==fabric) then
         i2l=Iglbloc_fbc(i)
          do j=1,ns_ibm_r_fbc(i2l)
            das=areapnt_fbc(i2l,j)
            temp1=temp1
     &      - ffluidsum1_fbc(i2l,j)*das
            tempb=tempb
     &      - ffluidsum1_fbc(i2l,j)*das
     &        *vibm1_fbc(i2l,j)
            temp2=temp2
     &      - ffluidsum2_fbc(i2l,j)*das
            temp=temp
     &      - ffluidsum2_fbc(i2l,j)*das
     &        *vibm2_fbc(i2l,j)
            temp3=temp3
     &      - ffluidsum3_fbc(i2l,j)*das
            tempc=tempc
     &      - ffluidsum3_fbc(i2l,j)*das
     &        *vibm3_fbc(i2l,j)
          enddo
        elseif(bodytype(i) ==fibrousshell) then
         i2l=Iglbloc_fsh(i)
          do jq=1,nq_ibm_r_fsh(i2l)
          do j=1,ns_ibm_r_fsh(i2l)
            das=ds_ibm_fsh(i2l,jq,j)*ds2_ibm_fsh(i2l,jq,j)

            temp1=temp1
     &      - ffluidsum1_fsh(i2l,jq,j)*das
            tempb=tempb
     &      - ffluidsum1_fsh(i2l,jq,j)*das
     &        *vibm1_fsh(i2l,jq,j)
            temp2=temp2
     &      - ffluidsum2_fsh(i2l,jq,j)*das
            temp=temp
     &      - ffluidsum2_fsh(i2l,jq,j)*das
     &        *vibm2_fsh(i2l,jq,j)
            temp3=temp3
     &      - ffluidsum3_fsh(i2l,jq,j)*das
            tempc=tempc
     &      - ffluidsum3_fsh(i2l,jq,j)*das
     &        *vibm3_fsh(i2l,jq,j)


            temp1Piezo=temp1Piezo
     &       +1.0/(piezo_Coef_fsh(i2l,jq,j,1)+1.d-9)
     &       *vPiezo_fsh(i2l,jq,j)**2*das

          enddo
          enddo
        endif

          temp1_local(i)=temp1
          temp2_local(i)=temp2
          temp3_local(i)=temp3
          temp1_localv(i)=tempb
          temp2_localv(i)=temp
          temp3_localv(i)=tempc

          tempPiezo_localp(i)=temp1Piezo


          temp2a_local(i)=0.d0
          temp2b_local(i)=0.d0
          temp2c_local(i)=0.d0
          temp2d_local(i)=0.d0
          temp2a_localforce(i)=0.d0
          temp2b_localforce(i)=0.d0
          temp2c_localforce(i)=0.d0

        if(bodytype(i) ==fiber) then
         i2l=Iglbloc_fib(i)
          do j=1,target_num(i)
              temp2a_local(i)=temp2a_local(i)
     &         +force_points_fib(i2l,j,1)*target_points_v_fib(i2l,j,1) 
              temp2b_local(i)=temp2b_local(i)
     &         +force_points_fib(i2l,j,2)*target_points_v_fib(i2l,j,2) 
              temp2c_local(i)=temp2c_local(i)
     &         +force_points_fib(i2l,j,3)*target_points_v_fib(i2l,j,3) 

              temp2a_localforce(i)=temp2a_localforce(i)
     &         +force_points_fib(i2l,j,1)
              temp2b_localforce(i)=temp2b_localforce(i)
     &         +force_points_fib(i2l,j,2)
              temp2c_localforce(i)=temp2c_localforce(i)
     &         +force_points_fib(i2l,j,3)
           enddo
        elseif(bodytype(i) ==shell) then
         i2l=Iglbloc_esh(i)
          do j=1,target_num(i)
              temp2a_local(i)=temp2a_local(i)
     &         +force_points_esh(i2l,j,1)*target_points_v_esh(i2l,j,1) 
              temp2b_local(i)=temp2b_local(i)
     &         +force_points_esh(i2l,j,2)*target_points_v_esh(i2l,j,2) 
              temp2c_local(i)=temp2c_local(i)
     &         +force_points_esh(i2l,j,3)*target_points_v_esh(i2l,j,3) 

              temp2a_localforce(i)=temp2a_localforce(i)
     &         +force_points_esh(i2l,j,1)
              temp2b_localforce(i)=temp2b_localforce(i)
     &         +force_points_esh(i2l,j,2)
              temp2c_localforce(i)=temp2c_localforce(i)
     &         +force_points_esh(i2l,j,3)
           enddo

        elseif(bodytype(i) ==fabric) then
         i2l=Iglbloc_fbc(i)
          do j=1,target_num(i)
              temp2a_local(i)=temp2a_local(i)
     &         +force_points_fbc(i2l,j,1)*target_points_v_fbc(i2l,j,1) 
              temp2b_local(i)=temp2b_local(i)
     &         +force_points_fbc(i2l,j,2)*target_points_v_fbc(i2l,j,2) 
              temp2c_local(i)=temp2c_local(i)
     &         +force_points_fbc(i2l,j,3)*target_points_v_fbc(i2l,j,3) 

              temp2a_localforce(i)=temp2a_localforce(i)
     &         +force_points_fbc(i2l,j,1)
              temp2b_localforce(i)=temp2b_localforce(i)
     &         +force_points_fbc(i2l,j,2)
              temp2c_localforce(i)=temp2c_localforce(i)
     &         +force_points_fbc(i2l,j,3)
           enddo

        elseif(bodytype(i) ==fibrousshell) then
         i2l=Iglbloc_fsh(i)
         ibdyg=i
         i2g=i
         do j=1,target_num(i)
             k=target_point_num_fsh(i2l,1,j)
             jq=target_point_num_fsh(i2l,2,j)
             k2=jq
             temp_ibm1=  
     &         (gx_ibm_massive_fsh(i2l,jq,k)-target_points_fsh(i2l,j,1))
             temp_ibm2=  
     &         (gy_ibm_massive_fsh(i2l,jq,k)-target_points_fsh(i2l,j,2))
             temp_ibm3=  
     &         (gz_ibm_massive_fsh(i2l,jq,k)-target_points_fsh(i2l,j,3))
             time_alphaf=time
             include "distIBM_target_k_rampfn.Fi"               
             if(time .le. target_t_link_fsh(i2l,j)) then
                 ksh_tmp=target_k_link_fsh(i2l,j)*target_timefn
             else
                 ksh_tmp=0.0
             endif
             force_points_fsh(i2l,j,1)=k_link*temp_ibm1*ksh_tmp
             force_points_fsh(i2l,j,2)=k_link*temp_ibm2*ksh_tmp
             if(ndim .gt. 2)  
     &          force_points_fsh(i2l,j,3)=k_link*temp_ibm3*ksh_tmp

             temp_ibm1=target_points_v_fsh(i2l,j,1)
             temp_ibm2=target_points_v_fsh(i2l,j,2)
             temp_ibm3=target_points_v_fsh(i2l,j,3)

             temp_ibm1=vibm1_pre_fsh(i2l,jq,k)-temp_ibm1
             temp_ibm2=vibm2_pre_fsh(i2l,jq,k)-temp_ibm2
             temp_ibm3=vibm3_pre_fsh(i2l,jq,k)-temp_ibm3

             force_points_fsh(i2l,j,1)=force_points_fsh(i2l,j,1)+
     &              cs_ibm_target(i2g)*k_link*temp_ibm1*ksh_tmp
             force_points_fsh(i2l,j,2)=force_points_fsh(i2l,j,2)+
     &              cs_ibm_target(i2g)*k_link*temp_ibm2*ksh_tmp
             if(ndim .gt. 2)  
     &          force_points_fsh(i2l,j,3)=force_points_fsh(i2l,j,3)+
     &              cs_ibm_target(i2g)*k_link*temp_ibm3*ksh_tmp
          enddo  !j



          do j=1,target_num(i)
              temp2a_local(i)=temp2a_local(i)
     &         +force_points_fsh(i2l,j,1)*target_points_v_fsh(i2l,j,1) 
              temp2b_local(i)=temp2b_local(i)
     &         +force_points_fsh(i2l,j,2)*target_points_v_fsh(i2l,j,2) 
              temp2c_local(i)=temp2c_local(i)
     &         +force_points_fsh(i2l,j,3)*target_points_v_fsh(i2l,j,3) 


              temp2a_localforce(i)=temp2a_localforce(i)
     &         +force_points_fsh(i2l,j,1)
              temp2b_localforce(i)=temp2b_localforce(i)
     &         +force_points_fsh(i2l,j,2)
              temp2c_localforce(i)=temp2c_localforce(i)
     &         +force_points_fsh(i2l,j,3)

           enddo
        endif
        temp2d_local(i)=
     &         temp2a_local(i)+temp2b_local(i)+temp2c_local(i)
      end do
200   format(2(i24,' '),28(f20.9,' '))

      do i=1,nr_ibm
	nfile=i*100+21
        write(nfile,200)
!     column 
     &       time_step,
!     1
     &       i,
!     2
     &       time,
!     3
     &       temp1_local(i),
!     4
     &       temp2_local(i),
!     5
     &       temp3_local(i),
!     6
     &       temp1_localv(i),
!     7
     &       temp2_localv(i),
!     8
     &       temp3_localv(i),
!     9
     &       temp1_localv(i)+temp2_localv(i)+temp3_localv(i),
!     10
     &       temp2a_local(i),
!11   
     &       temp2b_local(i),
!     12
     &       temp2c_local(i),
!     13
     &       temp2d_local(i),
!     14
     &       temp2a_localforce(i),
!     15
     &       temp2b_localforce(i),
!     16
     &       temp2c_localforce(i),
!     17
     &       tempPiezo_localp(i)
!     18
      enddo
      
      do i=1,nr_ibm
       if(nibmptrecord(i) .gt. 0) then
        do j=1,nibmptrecord(i)
         nfile2=3000+i*100+30+j
         write(num,'(i3.3)') i
         write(num2,'(i3.3)') j
         filename2='IBM_pos-'//num(1:3)//'_'//num2(1:3)//'.txt'
         open(nfile2,file = filename2,form='formatted',
     *        action='write',position='append')
         if(bodytype(i) ==fiber) then
          i2l=Iglbloc_fib(i)
          k=ibmptrecord(i,j,1)
          temprecord(1)=GX_IBM_MASSIVE_fib(i2l,k)
          temprecord(2)=GY_IBM_MASSIVE_fib(i2l,k)
          temprecord(3)=GZ_IBM_MASSIVE_fib(i2l,k)
         elseif(bodytype(i) ==fibrousshell) then
          i2l=Iglbloc_fsh(i)
          k=ibmptrecord(i,j,1)
          k2=ibmptrecord(i,j,2)
          temprecord(1)=GX_IBM_MASSIVE_fsh(i2l,k,k2)
          temprecord(2)=GY_IBM_MASSIVE_fsh(i2l,k,k2)
          temprecord(3)=GZ_IBM_MASSIVE_fsh(i2l,k,k2)
         elseif(bodytype(i) ==shell) then
          i2l=Iglbloc_esh(i)
          k=ibmptrecord(i,j,1)
          temprecord(1)=GX_IBM_MASSIVE_esh(i2l,k)
          temprecord(2)=GY_IBM_MASSIVE_esh(i2l,k)
          temprecord(3)=GZ_IBM_MASSIVE_esh(i2l,k)
         elseif(bodytype(i) ==fabric) then
          i2l=Iglbloc_fbc(i)
          k=ibmptrecord(i,j,1)
          temprecord(1)=GX_IBM_MASSIVE_fbc(i2l,k)
          temprecord(2)=GY_IBM_MASSIVE_fbc(i2l,k)
          temprecord(3)=GZ_IBM_MASSIVE_fbc(i2l,k)
         else
          stop
         endif
         write(nfile2,200)
!     column 
     &       time_step,
!     1    
     &       i,
!     2          
     &       time,
!     3                
     &    temprecord(1:3)
        enddo
!     Closeing file handles using nfile2 
        do j=1,nibmptrecord(i)
         nfile2=3000+i*100+30+j
         write(num,'(i3.3)') i
         write(num2,'(i3.3)') j
         filename2='IBM_pos-'//num(1:3)//'_'//num2(1:3)//'.txt'
         close(nfile2)
        enddo
       endif
      enddo

!     Closeing file handles using nfile
      do i=1,nr_ibm
       nfile=i*100+21
       write(num,'(i3.3)') i
       filename='forces_each-'//num(1:3)//'.txt'
       close(nfile)
      enddo
      endif
      return
      end   
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
