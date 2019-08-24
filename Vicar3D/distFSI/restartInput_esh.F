c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!@c
      subroutine restartInput_esh(nt0,nt,nnt)
!@T 
!routine to output gemtery of the FE shells in input format required 
! to run Vicar3D.
!@q
      USE HeaderFSI
      implicit none
      CHARACTER(LEN=60) :: Formatin,infile

      character  num*7, filename*60,numbdy*4
      integer nt,nt0,nnt,md,i,j,j2,k,ibdy,ibdyg,ne,jtype
     &  ,itype,ntype,nos,ntem,noj,ig,noi
     &  ,ItypecontELE(3),Rot_matPara(3)

      real*8 areapnt_esh(nr_ibm_esh,ns_ibm_esh)
     &  ,FiberWrite_esh(nr_ibm_esh,ns_ibm_esh,3)
     &  ,x(maxordern,3),x0(maxordern,3),x_FibVec(maxordern,3)
     &  ,s(1,maxordern), dsdv(2,maxordern),ddsdvdv(3,maxordern)
     &  ,xFibvecLoc0(3),xFibvecLoc(3),dxpos0(2,3),dxpos(2,3)
     &  ,tmpv1(3),tmpv2(3),hmat0(3,3),hmat(3,3),tmp
     &  ,jmat(2,2), amat1(2,maxordern),sx(2,maxordern)
      Rot_matPara(1)=1
      Rot_matPara(2)=2
      Rot_matPara(3)=3
      FiberWrite_esh=0.0

      if(mod(nt,nnt).eq.0)then
        md = nt0
	      write(num,'(i7.7)') md/1


        areapnt_esh=0.0
        do i=1,nr_ibm_esh
        ibdy=i
        ibdyg=Ilocglb_esh(i)
        do ne=1,numelr(ibdy)
         jtype=neatype_esh(ibdy,ne)
         itype=jtypeall(jtype)
         if(itype .eq. 0) then
           ntype=12  
         else
           ntype=6+itype
         end if
         if(itype .eq. 0) then
                  ItypecontELE(1)=4
                  ItypecontELE(2)=7
                  ItypecontELE(3)=8
         else              
                  ItypecontELE(1)=1
                  ItypecontELE(2)=2
                  ItypecontELE(3)=itype+1
         endif

         do nos=1,ntype
            ntem=nea_esh(ibdy,ne,nos)
            noj=1
               x(nos,noj)=gx_ibm_massive_esh(ibdy,ntem)
            noj=2
               x(nos,noj)=gy_ibm_massive_esh(ibdy,ntem)
            noj=3
               x(nos,noj)=gz_ibm_massive_esh(ibdy,ntem)

            noj=1
               x0(nos,noj)=gx_bp_esh0(ibdy,ntem)
            noj=2
               x0(nos,noj)=gy_bp_esh0(ibdy,ntem)
            noj=3
               x0(nos,noj)=gz_bp_esh0(ibdy,ntem)

         do noj=1,2
            sx(noj,nos)=coorcurv_esh(ibdy,ne,nos,noj)
         enddo
         enddo

         if(Fibrous_Esh(ibdy) .eq. 1) then
         do ig=1,ngau0
             s(1,1:ntype)=nn0(jtype,ig,1:ntype)
             dsdv(1,1:ntype)=dnndv0(jtype,ig,1:ntype)
             dsdv(2,1:ntype)=dnndw0(jtype,ig,1:ntype)

             ntem=nea_esh(ibdy,ne,ItypecontELE(ig))
             xFibvecLoc0(1:3)=Materialpara_esh(ibdy,ntem,1:3)

             do noi=1,2
               do noj=1,2
                  jmat(noj,noi)=
     &              dot_product(dsdv(noj,1:ntype),sx(noi,1:ntype))
               end do
             end do 

             call invert(jmat,2)
  
             amat1(1:2,1:ntype)=matmul(jmat,dsdv(1:2,1:ntype))

             do noj=1,3
                do noi=1,2
                   dxpos(noi,noj)=
     &               dot_product(amat1(noi,1:ntype),x(1:ntype,noj))
                   dxpos0(noi,noj)=
     &               dot_product(amat1(noi,1:ntype),x0(1:ntype,noj))
                end do
             end do

             call crossvec(tmpv1(1:3),dxpos0(1,1:3),dxpos0(2,1:3))  
             call crossvec(tmpv2(1:3),dxpos(1,1:3),dxpos(2,1:3))    
 

             hmat0(1:3,1)=dxpos0(1,1:3)
             hmat0(1:3,2)=dxpos0(2,1:3)
             hmat0(1:3,3)=tmpv1(1:3)

             call invert(hmat0,3)

             hmat(1:3,1)=dxpos(1,1:3)
             hmat(1:3,2)=dxpos(2,1:3)
             hmat(1:3,3)=tmpv2(1:3)

            hmat=matmul(hmat,hmat0)
            xFibvecLoc(1:3)=matmul(hmat,xFibvecLoc0)
            tmp=sqrt(
     &                xFibvecLoc(1)**2
     &               +xFibvecLoc(2)**2
     &               +xFibvecLoc(3)**2
     &               )

             xFibvecLoc(1:3)=
     &               xFibvecLoc(1:3)
     &               /tmp

             nos=ig

             ntem=nea_esh(ibdy,ne,ItypecontELE(nos))

             areapnt_esh(ibdy,ntem)=areapnt_esh(ibdy,ntem)
     &                                   +1.0

             noj=1
             FiberWrite_esh(ibdy,ntem,noj)=
     &       FiberWrite_esh(ibdy,ntem,noj)+ xFibvecLoc(noj)

             noj=2
             FiberWrite_esh(ibdy,ntem,noj)=
     &       FiberWrite_esh(ibdy,ntem,noj)+ xFibvecLoc(noj)

             noj=3
             FiberWrite_esh(ibdy,ntem,noj)=
     &       FiberWrite_esh(ibdy,ntem,noj)+ xFibvecLoc(noj)
        enddo
        endif
        enddo

        if(Fibrous_Esh(ibdy) .eq. 1) then             
        do j=1,ns_ibm_r_esh(ibdy)
         FiberWrite_esh(ibdy,j,1:3)=FiberWrite_esh(ibdy,j,1:3)
     &                            /areapnt_esh(ibdy,j)


         tmp=sqrt(
     &                FiberWrite_esh(ibdy,j,1)**2
     &               +FiberWrite_esh(ibdy,j,2)**2
     &               +FiberWrite_esh(ibdy,j,3)**2
     &               )

         FiberWrite_esh(ibdy,j,1:3)=
     &                FiberWrite_esh(ibdy,j,1:3)
     &                /tmp

        enddo
        endif
	      write(numbdy,'(i4.4)') ibdy
        infile='ShellRestart_'//numbdy(1:4)//'_'//num(1:7)//'.dat'
        open(9998,file=infile,status='unknown',form='formatted' )
        write(9998,'(a)' ) numbdy(1:4)
        write(9998,402)   
     &           ns_ibm_r(ibdyg) 
     &          ,ns_ibm_i(ibdyg) 
     &          ,ns_ibm_rall(ibdyg) 
     &          ,numelr(ibdyg) 
     &          ,numeli(ibdyg) 
     &          ,masterbdy(ibdyg)

        write(9998,402) Rot_matPara(1:3), n_matpara_esh(ibdy)
        write(9998,402) jtypeall(1:jtypemax)



       Formatin = "(2X,(I7,' '),3(F16.8,'  '),(I7,' '),20(F16.8,'  '))"

       j2=0
       do j=1,ns_ibm_r_esh(ibdy)
         j2=j2+1
        if(Fibrous_Esh(ibdy) .eq. 1) then
         if(n_matpara_esh(ibdy) .gt. 3) then
            write(9998,Formatin) j2   
     &                    ,gx_ibm_massive_esh(ibdy,j)   
     &                    ,gy_ibm_massive_esh(ibdy,j)   
     &                    ,gz_ibm_massive_esh(ibdy,j)   
     &                    ,target_ktype_esh(ibdy,j)
     &                    ,target_kvalue_esh(ibdy,j)
     &                    ,FiberWrite_esh(ibdy,j,1:3)
     &      ,Materialpara_esh(ibdy,j,4:n_matpara_esh(ibdy))
         else
            write(9998,Formatin) j2   
     &                    ,gx_ibm_massive_esh(ibdy,j)   
     &                    ,gy_ibm_massive_esh(ibdy,j)   
     &                    ,gz_ibm_massive_esh(ibdy,j)   
     &                    ,target_ktype_esh(ibdy,j)
     &                    ,target_kvalue_esh(ibdy,j)
     &                    ,FiberWrite_esh(ibdy,j,1:3)
         endif
        else
         if(n_matpara_esh(ibdy) .gt. 3) then
            write(9998,Formatin) j2   
     &                    ,gx_ibm_massive_esh(ibdy,j)   
     &                    ,gy_ibm_massive_esh(ibdy,j)   
     &                    ,gz_ibm_massive_esh(ibdy,j)   
     &                    ,target_ktype_esh(ibdy,j)
     &                    ,target_kvalue_esh(ibdy,j)
     &      ,Materialpara_esh(ibdy,j,1:n_matpara_esh(ibdy))
         else
            write(9998,Formatin) j2   
     &                    ,gx_ibm_massive_esh(ibdy,j)   
     &                    ,gy_ibm_massive_esh(ibdy,j)   
     &                    ,gz_ibm_massive_esh(ibdy,j)   
     &                    ,target_ktype_esh(ibdy,j)
     &                    ,target_kvalue_esh(ibdy,j)
         endif
       endif
       enddo

       Formatin = "(2X,(I7,' '),3(F16.8,'  '),5(I7,' '))"
       do j=1,ns_ibm_i_esh(ibdy)
          j2=j2+1
          write(9998,Formatin) j2,   
     &         gx_ibm_massive_esh(ibdy,j2),   
     &         gy_ibm_massive_esh(ibdy,j2),   
     &         gz_ibm_massive_esh(ibdy,j2), 
     &         bndycnd_esh(ibdy,j2),
     &         boundrylist_esh(ibdy,j,1),   
     &         boundrylist_esh(ibdy,j,2),   
     &         boundrylist_esh(ibdy,j,3)

       enddo
       Formatin = "(2X,(I7,' '),3(F16.8,'  '))"
       do j=ns_ibm_r_esh(ibdy)+ns_ibm_i_esh(ibdy)+1   
     &             ,ns_ibm_rall_esh(ibdy) 
            write(9998,Formatin) j
     &                    ,gx_ibm_massive_esh(ibdy,j)   
     &                    ,gy_ibm_massive_esh(ibdy,j)   
     &                    ,gz_ibm_massive_esh(ibdy,j)  
       enddo   

                    
       if (masterbdy(ibdyg) .gt. 0) then
           do j=1,ns_ibm_r_esh(ibdy)+ns_ibm_i_esh(ibdy)
               write(9998,402) j, npos_esh(ibdy,j) 
           end do  
       endif
       do j=1,ns_ibm_r_esh(ibdy)
            write(9998,402) nposele_esh(ibdy,j,1:maxordern-6+2)
       end do

       do j=1,numelr(ibdyg)
            write(9998,402) neatype_esh(ibdy,j)
     &                 , nea_esh(ibdy,j,1:maxordern)
       end do        
       close (9998)

      enddo
      endif



        do i=1,nr_ibm_esh

	   write(numbdy,'(i4.4)') i
      if(ndim .gt.2) then
          filename='ibmRes-3dEsh-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,4972)
4972      format('title = "3d immersed boundary data"')
      else
          filename='ibmRes-2dEsh-B'//numbdy(1:4)//'_'//num(1:7)//'.dat'
          open(840,file = filename)
          write(840,2499)
2499      format('title = "2d immersed boundary data"')
      endif
      if(flag_fractal_esh(i)) then
3522   format('variables= x,y,z,u,v,w,ax,ay,az,s11,s22,s12,s')
      write(840,3522)
      else
3322   format('variables= x,y,z,u,v,w,ax,ay,az,s')
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
     +            ,FiberWrite_esh(i,j,1)
     +            ,FiberWrite_esh(i,j,2)
     +            ,FiberWrite_esh(i,j,3)
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
     +            ,FiberWrite_esh(i,j,1)
     +            ,FiberWrite_esh(i,j,2)
     +            ,FiberWrite_esh(i,j,3)
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



401   format(20(f20.9,' '))
402   format(25(i7,' '))
403   format((i7,' '),20(f20.9,' '))
      return 
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
