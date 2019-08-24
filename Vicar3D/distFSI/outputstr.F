c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!@c
      subroutine outputstr(nt,nnt)
!@T routine to output data of the structures in Tecplot format.
!Tecplot files are named as \\
![[ibm-X{2,3}d{Fib,Fsh,Esh}-B{Body num}_{time step}.dat]]
!@q
!
!   routine to output data of the structure
!
      USE HeaderFSI
      implicit none
      integer nt,nnt,md,i,j,k,jtype,itype,jq,jq2,i2,inode(2)
      integer isec
      real*8  ysecIBM,zsecIBM
!
      integer  flag,kflag,mk1,mk2,mk3,iadd1,iadd2	
      character  num*7, filename*52,numbdy*4,numk*4

!
      if(mod(nt,nnt).eq.0)then
        md = nt
	write(num,'(i7.7)') md/1

         do i=1,nr_ibm_fib
            write(numbdy,'(i4.4)') i
      if(ndim .gt.2) then
         filename='ibm-3dFib-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,497)
497       format('title = "3d immersed boundary data"')
      else
         filename='ibm-2dFib-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,2497)
2497      format('title = "2d immersed boundary data"')
      endif
332   format('variables= x,y,z,u,v,w,p,lev,rho,s')
      write(840,332)
               if(nsec_IBM .le. 1) then
                 write(840,499) md,ns_ibm_r_fib(i),ns_ibm_r_fib(i)-1
499              format(1x,'zone T="ibm_',i7,'", n=',i7,', e=',
     &             i7,', zonetype=felineseg, datapacking=point')
                do j=1,ns_ibm_r_fib(i)
                   if(ndim .gt.2) then
                    write(840,401) 
     +              gx_ibm_massive_fib(i,j)
     +             ,gy_ibm_massive_fib(i,j)
     +             ,gz_ibm_massive_fib(i,j)
     +             ,vibm1_fib(i,j)
     +             ,vibm2_fib(i,j)
     +             ,vibm3_fib(i,j)
     +             ,0.d0,0.d0,0.d0
     +             ,tsibm_fib(i,j)
                  else
                   write(840,401) 
     +              gx_ibm_massive_fib(i,j)
     +             ,gy_ibm_massive_fib(i,j)
     +             ,0.0
     +             ,vibm1_fib(i,j)
     +             ,vibm2_fib(i,j)
     +             ,0.0
     +             ,0.d0,0.d0,0.d0
     +             ,tsibm_fib(i,j)
                  endif
                enddo
                do j=1,ns_ibm_r_fib(i)-1
                    k=j+1
                    write(840,*) k-1,k
                enddo
              else
                 write(840,498) md,ns_ibm_r_fib(i)*nsec_IBM
     &                        ,(ns_ibm_r_fib(i)-1)*2*(nsec_IBM-1)


               do isec=1,nsec_IBM
                 ysecIBM=real(isec*2-nsec_IBM-1)/2.0*dsecy_IBM
                 zsecIBM=real(isec*2-nsec_IBM-1)/2.0*dsecz_IBM
               do j=1,ns_ibm_r_fib(i)
                  if(ndim .gt.2) then
                   write(840,401) 
     +             gx_ibm_massive_fib(i,j)
     +            ,gy_ibm_massive_fib(i,j)+ysecIBM
     +            ,gz_ibm_massive_fib(i,j)+zsecIBM
     +            ,vibm1_fib(i,j)
     +            ,vibm2_fib(i,j)
     +            ,vibm3_fib(i,j)
     +            ,0.d0,0.d0,0.d0
     +            ,tsibm_fib(i,j)
                 else
                   write(*,*) 'ndim=2 and nsecIBM>1 is not allowed' 
                   STOP
                 endif
               enddo
               enddo

               do isec=1,nsec_IBM-1
               do j=1,ns_ibm_r_fib(i)-1
                    k=(isec-1)*ns_ibm_r_fib(i)+j
                    write(840,*) k,k+1,k+ns_ibm_r_fib(i)
                    write(840,*) k+1,k+ns_ibm_r_fib(i)+1
     +                               ,k+ns_ibm_r_fib(i)
               enddo
               enddo
             endif
           close(840)
        enddo !Nr_ibm_fib

        do i=1,nr_ibm_fsh

	   write(numbdy,'(i4.4)') i
      if(ndim .gt.2) then
          filename='ibm-3dFsh-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,4971)
4971      format('title = "3d immersed boundary data"')
      else
          filename='ibm-2dFsh-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,2498)
2498      format('title = "2d immersed boundary data"')
      endif
3321   format('variables= x,y,z,u,v,w,p,lev,rho,s')

      write(840,3321)

                 iadd1=0;iadd2=0
                 if (BoundaryFibShell_fsh(i,1,1) .eq. 3) iadd2=1
                 if (BoundaryFibShell_fsh(i,2,1) .eq. 3) iadd1=1

                 write(840,479) md,ns_ibm_r_fsh(i)+iadd2
     &                             ,nq_ibm_r_fsh(i)+iadd1
479              format(1x,'zone t="ibm_',i7,'", i=',i7,' j=',
     &             i7,' ,f=point')

               do jq=1,nq_ibm_r_fsh(i)
               do j=1,ns_ibm_r_fsh(i)
                  if(ndim .gt.2) then
                   write(840,401) 
     +             gx_ibm_massive_fsh(i,jq,j)
     +            ,gy_ibm_massive_fsh(i,jq,j)
     +            ,gz_ibm_massive_fsh(i,jq,j)
     +            ,vibm1_fsh(i,jq,j)
     +            ,vibm2_fsh(i,jq,j)
     +            ,vibm3_fsh(i,jq,j)
     +            ,0.d0,0.d0,0.d0
     +            ,tsibm_fsh(i,jq,j)
                 else
                   write(840,401) 
     +             gx_ibm_massive_fsh(i,jq,j)
     +            ,gy_ibm_massive_fsh(i,jq,j)
     +            ,0.0
     +            ,vibm1_fsh(i,jq,j)
     +            ,vibm2_fsh(i,jq,j)
     +            ,0.0
     +            ,0.d0,0.d0,0.d0
     +            ,tsibm_fsh(i,jq,j)
                 endif
               enddo
               if (BoundaryFibShell_fsh(i,1,1) .eq. 3) then
                  j=1
                  if(ndim .gt.2) then
                   write(840,401) 
     +             gx_ibm_massive_fsh(i,jq,j)
     +            ,gy_ibm_massive_fsh(i,jq,j)
     +            ,gz_ibm_massive_fsh(i,jq,j)
     +            ,vibm1_fsh(i,jq,j)
     +            ,vibm2_fsh(i,jq,j)
     +            ,vibm3_fsh(i,jq,j)
     +            ,0.d0,0.d0,0.d0
     +            ,tsibm_fsh(i,jq,j)
                 else
                   write(840,401) 
     +             gx_ibm_massive_fsh(i,jq,j)
     +            ,gy_ibm_massive_fsh(i,jq,j)
     +            ,0.0
     +            ,vibm1_fsh(i,jq,j)
     +            ,vibm2_fsh(i,jq,j)
     +            ,0.0
     +            ,0.d0,0.d0,0.d0
     +            ,tsibm_fsh(i,jq,j)
                 endif
               endif
               enddo

               if (BoundaryFibShell_fsh(i,2,1) .eq. 3) then
               jq=1
               do j=1,ns_ibm_r_fsh(i)
                  if(ndim .gt.2) then
                   write(840,401) 
     +             gx_ibm_massive_fsh(i,jq,j)
     +            ,gy_ibm_massive_fsh(i,jq,j)
     +            ,gz_ibm_massive_fsh(i,jq,j)
     +            ,vibm1_fsh(i,jq,j)
     +            ,vibm2_fsh(i,jq,j)
     +            ,vibm3_fsh(i,jq,j)
     +            ,0.d0,0.d0,0.d0
     +            ,tsibm_fsh(i,jq,j)
                 else
                   write(840,401) 
     +             gx_ibm_massive_fsh(i,jq,j)
     +            ,gy_ibm_massive_fsh(i,jq,j)
     +            ,0.0
     +            ,vibm1_fsh(i,jq,j)
     +            ,vibm2_fsh(i,jq,j)
     +            ,0.0
     +            ,0.d0,0.d0,0.d0
     +            ,tsibm_fsh(i,jq,j)
                 endif
               enddo
               if (BoundaryFibShell_fsh(i,1,1) .eq. 3) then
                  j=1
                  if(ndim .gt.2) then
                   write(840,401) 
     +             gx_ibm_massive_fsh(i,jq,j)
     +            ,gy_ibm_massive_fsh(i,jq,j)
     +            ,gz_ibm_massive_fsh(i,jq,j)
     +            ,vibm1_fsh(i,jq,j)
     +            ,vibm2_fsh(i,jq,j)
     +            ,vibm3_fsh(i,jq,j)
     +            ,aibm1_fsh(i,jq,j)
     +            ,aibm2_fsh(i,jq,j)
     +            ,aibm3_fsh(i,jq,j)
     +            ,tsibm_fsh(i,jq,j)
                 else
                   write(840,401) 
     +             gx_ibm_massive_fsh(i,jq,j)
     +            ,gy_ibm_massive_fsh(i,jq,j)
     +            ,0.0
     +            ,vibm1_fsh(i,jq,j)
     +            ,vibm2_fsh(i,jq,j)
     +            ,0.0
     +            ,aibm1_fsh(i,jq,j)
     +            ,aibm2_fsh(i,jq,j)
     +            ,0.0
     +            ,tsibm_fsh(i,jq,j)
                 endif
               endif
             endif
             close(840)
        enddo !Nr_ibm_fsh

        do i=1,nr_ibm_esh

	   write(numbdy,'(i4.4)') i
      if(ndim .gt.2) then
          filename='ibm-3dEsh-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,4972)
4972      format('title = "3d immersed boundary data"')
      else
          filename='ibm-2dEsh-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,2499)
2499      format('title = "2d immersed boundary data"')
      endif
      if(flag_fractal_esh(i)) then
3522   format('variables= x,y,z,u,v,w,s11,s22,s12,s')
      write(840,3522)
      else
3322   format('variables= x,y,z,u,v,w,au,av,aw,s')
      write(840,3322)
      endif



              write(840,498) time_step,ns_ibm_r_esh(i),numelr(i)


498      format(1x,'zone t="ibm',i7,'" ,n=',i7,' e=',
     &             i7,' ,zonetype=fetriangle, datapacking=point')
      if(flag_fractal_esh(i)) then
          do j=1,ns_ibm_r_esh(i)
                  write(840,401) 
     +             gx_ibm_massive_esh(i,j)
     +            ,gy_ibm_massive_esh(i,j)
     +            ,gz_ibm_massive_esh(i,j)
     +            ,vibm1_esh(i,j)
     +            ,vibm2_esh(i,j)
     +            ,vibm3_esh(i,j)
     +            ,fibmsig_esh(i,1,j)
     +            ,fibmsig_esh(i,2,j)
     +            ,fibmsig_esh(i,3,j)
     +            ,tsibm_esh(i,j)
          enddo
      else
          do j=1,ns_ibm_r_esh(i)
                  write(840,401) 
     +             gx_ibm_massive_esh(i,j)
     +            ,gy_ibm_massive_esh(i,j)
     +            ,gz_ibm_massive_esh(i,j)
     +            ,vibm1_esh(i,j)
     +            ,vibm2_esh(i,j)
     +            ,vibm3_esh(i,j)
     +            ,aibm1_esh(i,j)
     +            ,aibm2_esh(i,j)
     +            ,aibm3_esh(i,j)
     +            ,tsibm_esh(i,j)
          enddo
      endif
          do j=1,numelr(i)
             jtype=neatype_esh(i,j)
             itype=jtypeall(jtype)
             if(itype .eq. 0) then
                write(840,*) nea_esh(i,j,4)
     +                       ,nea_esh(i,j,7)
     +                       ,nea_esh(i,j,8)
             else
                write(840,*)  nea_esh(i,j,1)
     +                       ,nea_esh(i,j,2)
     +                       ,nea_esh(i,j,itype+1)
             endif
          enddo
       close(840)
        enddo !nr_ibm_esh
        do i=1,nr_ibm_fbc

	   write(numbdy,'(i4.4)') i
      if(ndim .gt.2) then
          filename='ibm-3dFbc-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,4976)
4976      format('title = "3d immersed boundary data for Fabric"')
      else
          filename='ibm-2dFbc-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,2501)
2501      format('title = "2d immersed boundary data for Fabric"')
      endif
3323   format('variables= x,y,z,u,v,w,au,av,aw,s')

        write(840,3323)
        write(840,4991) time_step,ns_ibm_r_fbc(i),numelr(i)

4991      format(1x,'zone t="ibm',i7,'" ,n=',i7,' e=',
     &             i7,' ,zonetype=fetriangle, datapacking=point')
          do j=1,ns_ibm_r_fbc(i)
                  write(840,401) 
     +             gx_ibm_massive_fbc(i,j)
     +            ,gy_ibm_massive_fbc(i,j)
     +            ,gz_ibm_massive_fbc(i,j)
     +            ,vibm1_fbc(i,j)
     +            ,vibm2_fbc(i,j)
     +            ,vibm3_fbc(i,j)
     +            ,aibm1_fbc(i,j)
     +            ,aibm2_fbc(i,j)
     +            ,aibm3_fbc(i,j)
     +            ,tsibm_fbc(i,j)
          enddo
          do j=1,numelr(i)
                write(840,*) nea_fbc(i,j,1)
     +                       ,nea_fbc(i,j,2)
     +                       ,nea_fbc(i,j,3)
          enddo
       close(840)
       if(FabricMeshFlag_fbc(i)) then
         do k=1,FabricMesh_nMesh(i)
	         write(numk,'(i3.3)') k
           filename='ibm-3dFbc-B'//numbdy(1:4)//'_M'//numk(1:4)
     +                          //'_'//num(1:7)//'.dat'
           open(840,file = filename)
           write(840,3776)
3776       format('title = "3d immersed boundary data for Mesh Fabric"')
3777       format('variables= x,y,z,u,v,w,au,av,aw,s')
           write(840,3777)
           write(840,3778) i,k,time_step
     &              ,FabricMesh_npoint(i,k),FabricMesh_nLine(i,k)
3778       format(1x,'zone t="ibmMeshofB',i4,' M',i4,'" ,n=',i7,' e=',
     &             i7,' ,f=point')
            do j=1,FabricMesh_npoint(i,k)
                  jq=FabricMesh_address(i,k,j)
                  write(840,401) 
     +             FabricMesh_coord(i,jq,1:3)
     +           ,FabricMesh_v(i,jq,1:3)
     +           ,FabricMesh_a(i,jq,1:3)
     +           ,FabricMesh_t(i,jq)
            enddo
            do j=1,FabricMesh_nLine(i,k)
                  jq=FabricMeshL_address(i,k,j)
                  do jq2=1,2
                  i2=FabricMeshL_con(i,jq,jq2)
                  inode(jq2)=FabricMeshL_addressRev(i,i2,2)
                  enddo
                  write(840,*) inode(1:2)
            enddo
             close(840)
          enddo
        endif
        if(FabricLineFlag_fbc(i)) then
         do k=1,FabricLine_nLine(i)
	         write(numk,'(i3.3)') k
           filename='ibm-3dFbc-B'//numbdy(1:4)//'_L'//numk(1:4)
     +                          //'_'//num(1:7)//'.dat'
           open(840,file = filename)
           write(840,3676)
3676       format('title = "3d immersed boundary data for Line Fabric"')
3677       format('variables= x,y,z,u,v,w,au,av,aw,s')
           write(840,3677)
           write(840,3678) i,k,time_step
     &            ,FabricLine_npoint(i,k),FabricLine_npoint(i,k)-1
3678       format(1x,'zone t="ibmLineofB',i4,' M',i4,'" ,n=',i7,' e=',
     &             i7,' ,f=point')
            do j=1,FabricLine_npoint(i,k)
                  jq=FabricLine_address(i,k,j)
                  write(840,401) 
     +             FabricLine_coord(i,jq,1:3)
     +           ,FabricLine_v(i,jq,1:3)
     +           ,FabricLine_a(i,jq,1:3)
     +           ,FabricLine_t(i,jq)
            enddo
            do j=1,FabricLine_npoint(i,k)-1
                  write(840,*) j, j+1
            enddo
            close(840)
          enddo
        endif
        enddo !nr_ibm_fbc

      endif
401   format(20(f20.9,' '))
      return 
      end 


c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
