c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine boundary_ibm_generalalpha(
     &       i2
     &      ,Xtar_np1_esh
     &      ,dXtar_np1_esh
     &      ,X0_np1_esh)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c satisfying boundary condition for ibm with gen-alpha time marching technqiue
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      implicit none

      integer i,j,i2,jq,i2l

      real*8 ttramp,stramp,kin_time,temp1,temp2
     &       ,temp1x(3),temp2x(3)
     &       ,temp1r(3),temp2r(3)
     &       ,ctmp1,ctmp2,ctmp3,stmp1,stmp2,stmp3
     &       ,tmpcoord(3),tmpcoordV1(3),tmpcoordV2(3)

      real*8  kin_t(3),kin_omega(3),kin_amp(3),kin_theta(3),kin_flag(3)
     &       ,kin_X0(3)
     &       ,kin_ax(3),kin_ax_v(3)
     &       ,kin_omegax(3),kin_omegar(3)
     &       ,kin_ampx(3),kin_ampr(3)
     &       ,kin_tx(3),kin_tr(3)
     &       ,kin_phix(3),kin_phir(3),time_next


      real*8 tmpangle(3),dtmpangle(3),tmpmotion(3),dtmpmotion(3)
     &      ,s_function,ds_function


      real*8 tmpRmat(3,3),dtmpRmat(3,3)
     &       ,dtmpRmat_da1(3,3),dtmpRmat_da2(3,3),dtmpRmat_da3(3,3)

!=== user external kinemtic model 
      real*8 usrmovefn(6), rtmpusrfn(50)
      integer itmpusrfn(50)
      integer kin_seq,iseqtmp,i0tmp,i1tmp,i2tmp,i3tmp
      real*8 Rmattmp(1:3,1:3,1:3),dRmattmp(1:3,1:3,1:3)
      real*8 Xtar_np1_esh(Ns_IBM_esh,3),dXtar_np1_esh(Ns_IBM_esh,3)
     &      ,X0_np1_esh(Ns_IBM_esh,3)


      pi1=4.0d0*atan(1.d0)

      ttramp=tramps(i2,2)
      stramp=tramps(i2,1)

      kin_t(1:3)=1.0/(frequency_fin(i2,1:3)+1.0d-9)
 
      kin_omega(1:3)=2.*(pi1)/kin_t(1:3)

      do j=1,3
      kin_omegax(j)=2.*(pi1)*freqx(i2,j)
      kin_omegar(j)=2.*(pi1)*freqt(i2,j)
      kin_ampr(j)=ampt(i2,j)
      kin_ampx(j)=ampx(i2,j)
      kin_tx(j)=1.0/(freqx(i2,j)+1.0d-9)
      kin_tr(j)=1.0/(freqt(i2,j)+1.0d-9)
      kin_phix(j)=phix(i2,j)
      kin_phir(j)=phit(i2,j)
      kin_X0(j)=Gt0(i2,j)
      kin_seq=rot_seq(i2)
      enddo

      stramp=max(stramp,0.0005d0)

      time_next=time
      kin_time=time_next

      if (time_next .le. stramp) then
             temp1=0.d0
             temp2=0.d0
             temp1x(1:3)=0.d0
             temp2x(1:3)=0.d0
             temp1r(1:3)=0.d0
             temp2r(1:3)=0.d0
      else
         kin_time=time_next-stramp
!---------------  change feb26-----------------------------
         temp1=1.d0-
     &   dexp(-(time_next-stramp)/(ttramp))
         temp2=1.d0/(ttramp)
     &   *dexp(-(time_next-stramp)/(ttramp))

      do j=1,3
         temp1x(j)=1.d0-
     &   dexp(-(time_next-stramp)/(trampx(i2,j)))
         temp2x(j)=1.d0/trampx(i2,j)
     &   *dexp(-(time_next-stramp)/(trampx(i2,j)))

         temp1r(j)=1.d0-
     &   dexp(-(time_next-stramp)/(trampt(i2,j)))
         temp2r(j)=1.d0/trampt(i2,j)
     &   *dexp(-(time_next-stramp)/(trampt(i2,j)))
      enddo
!-----------------------------------------------------------
      endif

      do j=1,3
!rotational RBM
          if(kin_omegar(j) .ge. 0.0)then
          tmpangle(j)=kin_ampr(j)
     &                *sin(kin_omegar(j)*kin_time+kin_phir(j))

          dtmpangle(j)=tmpangle(j)*temp2r(j)
     &                +kin_ampr(j)*kin_omegar(j)*temp1r(j)
     &                *cos(kin_omegar(j)*kin_time+kin_phir(j))
          tmpangle(j)=tmpangle(j)*temp1r(j)
          else
          ctmp1= abs(kin_omegar(j))/2.0/pi1
          stmp1= kin_phir(j)*180/pi1
          ctmp2= (kin_time-stmp1)/ctmp1
          stmp2= 1.0/ctmp1
          
          ctmp3= s_function(ctmp2)
          stmp3= stmp2*ds_function(ctmp2)           
          tmpangle(j)=kin_ampr(j)
     &                *ctmp3
          dtmpangle(j)=kin_ampr(j)
     &                *stmp3        
          endif
!translational RBM
          if(kin_omegax(j) .ge. 0.0)then
          tmpmotion(j)=kin_ampx(j)
     &                *sin(kin_omegax(j)*kin_time+kin_phix(j))

          dtmpmotion(j)=tmpmotion(j)*temp2x(j)
     &                +kin_ampx(j)*kin_omegax(j)*temp1x(j)
     &                *cos(kin_omegax(j)*kin_time+kin_phix(j))
          tmpmotion(j)=tmpmotion(j)*temp1x(j)
          else
          ctmp1= abs(kin_omegax(j))/2.0/pi1
          stmp1= kin_phix(j)*180/pi1
          ctmp2= (kin_time-stmp1)/ctmp1
          stmp2= 1.0/ctmp1
          
          ctmp3= s_function(ctmp2)
          stmp3= stmp2*ds_function(ctmp2)           

          tmpmotion(j)=kin_ampx(j)
     &                *ctmp3
          dtmpmotion(j)=kin_ampx(j)
     &                *stmp3        
          endif
      enddo


       ctmp1=cos(tmpangle(1))
       ctmp2=cos(tmpangle(2))
       ctmp3=cos(tmpangle(3))
       stmp1=sin(tmpangle(1))
       stmp2=sin(tmpangle(2))
       stmp3=sin(tmpangle(3))


       Rmattmp(1,1,1:3)=(/1.0d0,   0.0d0,   0.0d0/)
       Rmattmp(1,2,1:3)=(/0.0d0, ctmp1,-stmp1/)              
       Rmattmp(1,3,1:3)=(/0.0d0, stmp1, ctmp1/)


       Rmattmp(2,1,1:3)=(/ ctmp2,0.0d0, stmp2/)
       Rmattmp(2,2,1:3)=(/0.0d0   ,1.0d0,   0.0d0/)              
       Rmattmp(2,3,1:3)=(/-stmp2,0.0d0, ctmp2/)


       Rmattmp(3,1,1:3)=(/ ctmp3,-stmp3,0.0d0/)
       Rmattmp(3,2,1:3)=(/ stmp3, ctmp3,0.0d0/)             
       Rmattmp(3,3,1:3)=(/0.0d0   ,0.0d0   ,1.0d0/)

       iseqtmp=mod(kin_seq,10)

       i0tmp=int((kin_seq+1.0e-5)/10.0)
       i1tmp=int(i0tmp/100)
       i3tmp=mod(i0tmp,10)
       i2tmp=mod(int(i0tmp/10),10)


       if(iseqtmp.eq. 0)  then !current coordinate
       tmpRmat(1:3,1:3)=matmul(Rmattmp(i2tmp,1:3,1:3)
     &               ,Rmattmp(i3tmp,1:3,1:3))
       tmpRmat(1:3,1:3)=matmul(Rmattmp(i1tmp,1:3,1:3)
     &               ,tmpRmat(1:3,1:3))
       else !fixed coordinate
       tmpRmat(1:3,1:3)=matmul(Rmattmp(i2tmp,1:3,1:3)
     &               ,Rmattmp(i1tmp,1:3,1:3))
       tmpRmat(1:3,1:3)=matmul(Rmattmp(i3tmp,1:3,1:3)
     &               ,tmpRmat(1:3,1:3))
       endif


       dRmattmp(1,1,1:3)=(/0.0d0,   0.0d0,   0.0d0/)
       dRmattmp(1,2,1:3)=(/0.0d0,-stmp1,-ctmp1/)              
       dRmattmp(1,3,1:3)=(/0.0d0, ctmp1,-stmp1/)


       dRmattmp(2,1,1:3)=(/-stmp2,0.0d0, ctmp2/)
       dRmattmp(2,2,1:3)=(/0.0d0   ,0.0d0,   0.0d0/)              
       dRmattmp(2,3,1:3)=(/-ctmp2,0.0d0,-stmp2/)


       dRmattmp(3,1,1:3)=(/-stmp3,-ctmp3,0.0d0/)
       dRmattmp(3,2,1:3)=(/ ctmp3,-stmp3,0.0d0/)             
       dRmattmp(3,3,1:3)=(/0.0d0   ,0.0d0   ,0.0d0/)


       if(iseqtmp.eq. 0)  then !current coordinate
       dtmpRmat_da1(1:3,1:3)=matmul(Rmattmp(i2tmp,1:3,1:3)
     &               ,dRmattmp(i3tmp,1:3,1:3))
       dtmpRmat_da1(1:3,1:3)=matmul(Rmattmp(i1tmp,1:3,1:3)
     &               ,dtmpRmat_da1(1:3,1:3))
       dtmpRmat_da2(1:3,1:3)=matmul(dRmattmp(i2tmp,1:3,1:3)
     &               ,Rmattmp(i3tmp,1:3,1:3))
       dtmpRmat_da2(1:3,1:3)=matmul(Rmattmp(i1tmp,1:3,1:3)
     &               ,dtmpRmat_da2(1:3,1:3))
       dtmpRmat_da3(1:3,1:3)=matmul(Rmattmp(i2tmp,1:3,1:3)
     &               ,Rmattmp(i3tmp,1:3,1:3))
       dtmpRmat_da3(1:3,1:3)=matmul(dRmattmp(i1tmp,1:3,1:3)
     &               ,dtmpRmat_da2(1:3,1:3))
       else !fixed coordinate
       dtmpRmat_da1(1:3,1:3)=matmul(Rmattmp(i2tmp,1:3,1:3)
     &               ,Rmattmp(i1tmp,1:3,1:3))
       dtmpRmat_da1(1:3,1:3)=matmul(dRmattmp(i3tmp,1:3,1:3)
     &               ,dtmpRmat_da1(1:3,1:3))
       dtmpRmat_da2(1:3,1:3)=matmul(dRmattmp(i2tmp,1:3,1:3)
     &               ,Rmattmp(i1tmp,1:3,1:3))
       dtmpRmat_da2(1:3,1:3)=matmul(Rmattmp(i3tmp,1:3,1:3)
     &               ,dtmpRmat_da2(1:3,1:3))
       dtmpRmat_da3(1:3,1:3)=matmul(Rmattmp(i2tmp,1:3,1:3)
     &               ,dRmattmp(i1tmp,1:3,1:3))
       dtmpRmat_da3(1:3,1:3)=matmul(Rmattmp(i3tmp,1:3,1:3)
     &               ,dtmpRmat_da3(1:3,1:3))
       endif


      do j=1,3
         dtmpRmat(j,1:3)= 
     &        dtmpRmat_da1(j,1:3)*dtmpangle(1)
     &       +dtmpRmat_da2(j,1:3)*dtmpangle(2)
     &       +dtmpRmat_da3(j,1:3)*dtmpangle(3)
      enddo





      if(BodyType(i2)==shell) then
        i2l=Iglbloc_esh(i2)
      do i=1,target_num(i2)
          
          kin_amp(1:3)=a_fin_esh(i2l,i,1:3)
          kin_theta(1:3)=phi_fin_esh(i2l,i,1:3)/180.0d0*pi1
          kin_flag=dble(flagfixed_esh(i2l,i))

         do j=1,3
          if(kin_omega(j) .ge. 0.0)then
          kin_ax(j)=
     &      kin_amp(j)*sin(kin_omega(j)*kin_time+kin_theta(j))
          kin_ax_v(j)=
     &      kin_ax(j) * temp2
     &     +kin_amp(j)*kin_omega(j)* temp1
     &                *cos(kin_omega(j)*kin_time+kin_theta(j))
          kin_ax(j)=kin_ax(j)*temp1   
          else
          ctmp1= abs(kin_omega(j))/2.0/pi1
          stmp1= kin_theta(j)*180/pi1
          ctmp2= (kin_time-stmp1)/ctmp1
          stmp2= 1.0/ctmp1
          
          ctmp3= s_function(ctmp2)
          stmp3= stmp2*ds_function(ctmp2)           

          kin_ax(j)=kin_amp(j)
     &                *ctmp3
          kin_ax_v(j)=kin_amp(j)
     &                *stmp3        
          endif
 
         enddo

          i_fix=target_point_num_esh(i2l,i)

          tmpcoord(1:3)=
     &           (/gx_bp_esh0(i2l,i_fix)-kin_X0(1)+kin_ax(1)
     &           , gy_bp_esh0(i2l,i_fix)-kin_X0(2)+kin_ax(2)
     &           , gz_bp_esh0(i2l,i_fix)-kin_X0(3)+kin_ax(3) /)

          tmpcoord(1:3)=matmul(tmpRmat,tmpcoord(1:3))

          tmpcoordV1(1:3)=matmul(dtmpRmat,tmpcoord(1:3))

          tmpcoordV2(1:3)=
     &           (/kin_ax_v(1)
     &           , kin_ax_v(2)
     &           , kin_ax_v(3) /)
          tmpcoordV2(1:3)=matmul(tmpRmat,tmpcoordV2(1:3))

          tmpcoordV1(1:3)=tmpcoordV1(1:3)+tmpcoordV2(1:3)

!          include "distIBM_target_movefn.Fi"   
          do j=1,3
          Xtar_np1_esh(i,j)=kin_X0(j)+tmpmotion(j)+tmpcoord(j)
          dXtar_np1_esh(i,j)= dtmpmotion(j)+tmpcoordV1(j)
          enddo
       enddo
       do j=1,ns_ibm_r_esh(i2l)

          tmpcoord(1:3)=
     &           (/gx_bp_esh0(i2l,j)-kin_X0(1)
     &           , gy_bp_esh0(i2l,j)-kin_X0(2)
     &           , gz_bp_esh0(i2l,j)-kin_X0(3) /)

          tmpcoord(1:3)=matmul(tmpRmat,tmpcoord(1:3))

          X0_np1_esh(j,1:3)=kin_X0(1:3)+tmpmotion(1:3)+tmpcoord(1:3)
        enddo
      endif
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
