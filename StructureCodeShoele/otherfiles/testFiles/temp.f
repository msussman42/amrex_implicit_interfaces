
!CAUTION:need to be rewrite for fibrous shell
      if(iksh1 .ne. 0) then
        imaster=masterbdy(ibdyg)
        if(BodyType(imaster) .ne. shell) then 
           write(*,*)'The master should be shell'
          stop
        endif
        if(BodyType(ibdyg) .eq. fiber) then 
        do noi=1,ns_ibm_r_fib(ibdy)
           ntemi=npos_fib(ibdy,noi)
           do noj=1,ns_ibm_r_fib(ibdy)
              ntemj=npos_fib(ibdy,noj) 
              aamat_total(imaster,ntemi,ntemj)=
     &        aamat_total(imaster,ntemi,ntemj)
     &                   +aamat_fib(noi,noj)
           end do
           bbmat_total(imaster,ntemi,1)=bbmat_total(imaster,ntemi,1)
     &                                +fibm1_fib(i,noi)
           bbmat_total(imaster,ntemi,2)=bbmat_total(imaster,ntemi,2)
     &                                +fibm2_fib(i,noi)
           if(ndim .gt. 2) 
     &     bbmat_total(imaster,ntemi,3)=bbmat_total(imaster,ntemi,3)
     &                                +fibm3_fib(i,noi)
        end do
        endif

        if(BodyType(ibdyg) .eq. shell) then 
        do noi=1,ns_ibm_r_esh(ibdy)
           ntemi=npos_esh(ibdy,noi)
           do noj=1,ns_ibm_r_esh(ibdy)
              ntemj=npos_esh(ibdy,noj) 
              aamat_total(imaster,ntemi,ntemj)=
     &        aamat_total(imaster,ntemi,ntemj)
     &                   +aamat_esh(noi,noj)
           end do
           bbmat_total(imaster,ntemi,1)=bbmat_total(imaster,ntemi,1)
     &                                +fibm1_esh(i,noi)
           bbmat_total(imaster,ntemi,2)=bbmat_total(imaster,ntemi,2)
     &                                +fibm2_esh(i,noi)
           if(ndim .gt. 2) 
     &     bbmat_total(imaster,ntemi,3)=bbmat_total(imaster,ntemi,3)
     &                                +fibm3_esh(i,noi)
        end do
        endif
      end if
      

      end do  !ibdy


!CAUTION:need to be rewrite for fibrous shell

      if(iksh1 .ne. 0) then
         jq=1
         do imaster=1,Nr_IBM_mas
         i=imaster
        if(BodyType(i) .ne. shell) then 
           write(*,*)'The master should be shell'
          stop
        endif
         if(ntermiMAXm(imaster) .le. 1)  then      
              ntemi=0
              do noj=1,ns_ibm_rall(i)
                 iacsr(noj)=ntemi+1
                 do noi=1,ns_ibm_rall(i)
                    ksh_tmp=aamat_total(imaster,noj,noi)
                    if(abs(ksh_tmp) .gt. 1d-12) then
                       ntemi=ntemi+1
                    endif
                 enddo
              enddo
              ntermiMAXm(imaster)=ntemi+1 
!         k_flagMAT=.false.
         allocate(aacsrm(ntermiMAXm(imaster))
     &            ,jacsrm(ntermiMAXm(imaster)))
         aacsrm=0.0
         jacsrm=0
         iacsrm=0
         ntemi=0
         do noj=1,ns_ibm_rall(i)
           iacsrm(noj)=ntemi+1
           do noi=1,ns_ibm_rall(i)
             ksh_tmp=aamat_total(imaster,noj,noi)
             if(abs(ksh_tmp) .gt. 1d-12) then
               ntemi=ntemi+1
               aacsrm(ntemi)=ksh_tmp
               jacsrm(ntemi)=noi
             endif
           enddo
          enddo
          iacsrm(noj)=ntemi+1 
          ntermiMAXmaster=ntermiMAXm(imaster)
         endif

         i2l=Iglbloc_esh(imaster)
          call gmres_csr(
     &      aacsrm(1: ntermiMAXmaster)
     &     ,jacsrm(1: ntermiMAXmaster)
     &     ,iacsrm(1: ns_ibm_rall(i)+1)
     &     ,bbmat_total(imaster,1: ns_ibm_rall(i),1)
     &     ,gx_ibm_massive_esh(i2l,1: ns_ibm_rall(i))
     &     ,tol,ns_ibm_rall(i),ntermiMAXmaster)

          call gmres_csr(
     &      aacsrm(1: ntermiMAXmaster)
     &     ,jacsrm(1: ntermiMAXmaster)
     &     ,iacsrm(1: ns_ibm_rall(i)+1)
     &     ,bbmat_total(imaster,1: ns_ibm_rall(i),2)
     &     ,gy_ibm_massive_esh(i2l,1: ns_ibm_rall(i))
     &     ,tol,ns_ibm_rall(i),ntermiMAXmaster)
         if(ndim .gt. 2) 
     &      call gmres_csr(
     &      aacsrm(1: ntermiMAXmaster)
     &     ,jacsrm(1: ntermiMAXmaster)
     &     ,iacsrm(1: ns_ibm_rall(i)+1)
     &     ,bbmat_total(imaster,1: ns_ibm_rall(i),3)
     &     ,gz_ibm_massive_esh(i2l,1:ns_ibm_rall(i))
     &     ,tol,ns_ibm_rall(i),ntermiMAXmaster)
c initialized master matrices for later use
         aamat_total(imaster,1: ns_ibm_rall(i),1: ns_ibm_rall(i))=0.0
         bbmat_total(imaster,1: ns_ibm_rall(i),1:3)=0.0
      enddo
      do ibdy=1,nr_ibm
        i=ibdy
        imaster=masterbdy(ibdy)
        if(BodyType(i) .eq. fiber) then 
        i2l=Iglbloc_fib(i)
        do noi=1,ns_ibm_r_fib(i2l)
           ntemi=npos_fib(i2l,noi)
           gx_ibm_massive_fib(i2l,noi)=gx_ibm_massive_esh(imaster,ntemi)
           gy_ibm_massive_fib(i2l,noi)=gy_ibm_massive_esh(imaster,ntemi)
           if(ndim .gt. 2) 
     &    gz_ibm_massive_fib(i2l,noi)=gz_ibm_massive_esh(imaster,ntemi)
        end do
        elseif(BodyType(i) .eq. shell) then 
        i2l=Iglbloc_esh(i)
        do noi=1,ns_ibm_r_esh(i2l)
           ntemi=npos_esh(i2l,noi)
           gx_ibm_massive_esh(i2l,noi)=gx_ibm_massive_esh(imaster,ntemi)
           gy_ibm_massive_esh(i2l,noi)=gy_ibm_massive_esh(imaster,ntemi)
           if(ndim .gt. 2) 
     &     gz_ibm_massive_esh(i2l,noi)=gz_ibm_massive_esh(imaster,ntemi)
        end do
        endif
      end do ! ibdy
      end if !iksh1

      do ibdy=1,nr_ibm
        i=ibdy
        if(BodyType(i) .eq. fiber) then 
        i2l=Iglbloc_fib(i)
        do j=1,ns_ibm_r_fib(i2l)
          vibm1_fib(i2l,j)= 
     &        (gx_ibm_massive_fib(i2l,j)-gx_ibm_fib(i2l,j))/dt
          vibm2_fib(i2l,j)= 
     &        (gy_ibm_massive_fib(i2l,j)-gy_ibm_fib(i2l,j))/dt
          if(ndim .gt. 2) 
     &     vibm3_fib(i2l,j)= 
     &        (gz_ibm_massive_fib(i2l,j)-gz_ibm_fib(i2l,j))/dt
        end do

        elseif((BodyType(i) .eq. shell).and. 
     &         (Genalpha_timesolver(i) .le. 1)) then 
        i2l=Iglbloc_esh(i)
        do j=1,ns_ibm_r_esh(i2l)
          vibm1_esh(i2l,j)= 
     &        (gx_ibm_massive_esh(i2l,j)-gx_ibm_esh(i2l,j))/dt
          vibm2_esh(i2l,j)= 
     &        (gy_ibm_massive_esh(i2l,j)-gy_ibm_esh(i2l,j))/dt
          if(ndim .gt. 2) 
     &     vibm3_esh(i2l,j)= 
     &        (gz_ibm_massive_esh(i2l,j)-gz_ibm_esh(i2l,j))/dt
        end do
        elseif(BodyType(i) .eq. fibrousshell) then 
        i2l=Iglbloc_fsh(i)
        do jq=1,nq_ibm_r_fsh(i2l)
        do j=1,ns_ibm_r_fsh(i2l)
          vibm1_fsh(i2l,jq,j)= 
     &     (gx_ibm_massive_fsh(i2l,jq,j)-gx_ibm_fsh(i2l,jq,j))/dt
          vibm2_fsh(i2l,jq,j)= 
     &     (gy_ibm_massive_fsh(i2l,jq,j)-gy_ibm_fsh(i2l,jq,j))/dt
          if(ndim .gt. 2) 
     &     vibm3_fsh(i2l,jq,j)= 
     &     (gz_ibm_massive_fsh(i2l,jq,j)-gz_ibm_fsh(i2l,jq,j))/dt

        end do
        end do 
        endif 

      end do ! ibdy

      do i=1,nr_ibm_esh
       i2g=Ilocglb_esh(i)
      if((Genalpha_timesolver(i2g) .gt. 1)) then
       fibm1_esh(i,1:Ns_IBM_esh)=0.0
       fibm2_esh(i,1:Ns_IBM_esh)=0.0
       fibm3_esh(i,1:Ns_IBM_esh)=0.0

       if(FluidForceFlag(i2g)) then
        do j=1,ns_ibm_r_esh(i)
           fk_mass1_esh(i,j)=fk_mass1o_esh(i,j)
     &          +alpha_ibm*delta_t*(uibm1_esh(i,j)-vibm1_esh(i,j))
     &                   *timetemp2

            fk_mass2_esh(i,j)=fk_mass2o_esh(i,j)
     &          +alpha_ibm*delta_t*(uibm2_esh(i,j)-vibm2_esh(i,j))
     &                   *timetemp2

            ffluidsum1_esh(i,j)= fk_mass1_esh(i,j)
     &                   +beta_ibm*(uibm1_esh(i,j)-vibm1_esh(i,j))
     &                   *timetemp2
            ffluidsum2_esh(i,j)= fk_mass2_esh(i,j)
     &                   +beta_ibm*(uibm2_esh(i,j)-vibm2_esh(i,j))
     &                   *timetemp2
            if(ndim .gt. 2) then
            fk_mass3_esh(i,j)=fk_mass3o_esh(i,j)
     &          +alpha_ibm*delta_t*(uibm3_esh(i,j)-vibm3_esh(i,j))
     &                   *timetemp2

            ffluidsum3_esh(i,j)= fk_mass3_esh(i,j)
     &                   +beta_ibm*(uibm3_esh(i,j)-vibm3_esh(i,j))
     &                   *timetemp2
            endif
         end do

         do j=1,ns_ibm_r_esh(i)
              fibm1_esh(i,j)=fibm1_esh(i,j)+ffluidsum1_esh(i,j)
              fibm2_esh(i,j)=fibm2_esh(i,j)+ffluidsum2_esh(i,j)
              if(ndim .gt. 2) 
     &           fibm3_esh(i,j)=fibm3_esh(i,j)+ffluidsum3_esh(i,j)
          end do
        endif  
      endif  !Genalpha_timesolver(i2g)
      end do  !nr_ibm_esh

3000    format(1x,24(f8.4,2x))
3001    format(1x,1(i5,2x),24(f8.4,2x))
3002    format(1x,2(i5,2x),24(f8.4,2x))
3003    format(1x,3(i5,2x),24(f8.4,2x))
3004    format(1x,4(i5,2x),24(f8.4,2x))
3005    format(1x,5(i5,2x),24(f8.4,2x))
3006    format(1x,6(i5,2x),24(f8.4,2x))
3007    format(1x,7(i5,2x),24(f8.4,2x))
3008    format(1x,8(i5,2x),24(f8.4,2x))
3009    format(1x,9(i5,2x),24(f8.4,2x))
3010    format(1x,10(i5,2x),24(f8.4,2x))
3011    format(1x,11(i5,2x))

      return 
      end

