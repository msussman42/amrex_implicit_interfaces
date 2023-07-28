c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine Inplane_LineDelingette(Ftmp
     &                           ,x0,y0,z0
     &                           ,x1,y1,z1
     &                           ,x2,y2,z2
     &                           ,cf_coef
     &                           ,D2L,DL,L0,L
     &                           ,Lin_NonlinMode
     &                           ,KK
     &                           ,CC
     &                           ,Ktmp,KtmpFlag)

      real x0,y0,z0,x1,y1,z1,x2,y2,z2,cf_coef
      real Ftmp(3,3),Ktmp(3,3,3,3),KK(3),CC(3),D2L(3),DL(3),L0(3),L(3)
      logical KtmpFlag
      integer Lin_NonlinMode
      integer LINmode,NonLinmode
      integer IJ(3,2),IK(3,2),iQ(3,3)
      real Q(3,3),dQ(3,3),sQ(3,3),vec1(3),veca(6,3)
     &    ,ctmpi,ctmpj,ctmpk,ktmpi,ktmpj,ktmpk
      integer i,j,k,im,km,jm,j1,m
      PARAMETER(LINmode=0,NonLinmode=1)

       Ftmp(1:2,1:3)=0.0
       IJ(1,1)=2;IJ(1,2)=3
       IK(1,1)=3;IK(1,2)=2 

       IJ(2,1)=1;IJ(2,2)=3
       IK(2,1)=3;IK(2,2)=1 

       IJ(3,1)=1;IJ(3,2)=2
       IK(3,1)=2;IK(3,2)=1

       Q(1,1:3)=(/x0,x1,x2/) 
       Q(2,1:3)=(/y0,y1,y2/) 
       Q(3,1:3)=(/z0,z1,z2/) 

       dQ(1,1:3)=Q(3,1:3)-Q(2,1:3)
       dQ(2,1:3)=Q(3,1:3)-Q(1,1:3)
       dQ(3,1:3)=Q(2,1:3)-Q(1,1:3)

       iQ(2,1)=3;   sQ(2,1)= 1.0
       iQ(1,2)=3;   sQ(1,2)=-1.0

       iQ(3,1)=2;   sQ(3,1)= 1.0
       iQ(1,3)=2;   sQ(1,3)=-1.0

       iQ(3,2)=1;   sQ(3,2)= 1.0
       iQ(2,3)=1;   sQ(2,3)=-1.0


       if(Lin_NonlinMode == NonLinmode) then
         do i=1,3
            do m=1,2
               j=IJ(i,m)
               k=IK(i,m)

               j1=iQ(j,i)
               vec1(1:3)=dQ(j1,1:3)*sQ(j,i)
               Ftmp(i,1:3)=vec1(1:3)*D2L(k)*KK(k)
     &                    +vec1(1:3)*(CC(j)*D2L(i)+CC(i)*D2L(j))
            enddo
         enddo
       elseif(Lin_NonlinMode == Linmode) then  
         do i=1,3
            do m=1,2
               j=IJ(i,m)
               k=IK(i,m)

               k1=iQ(k,i)
               ktmpj=KK(j)*2.0*L0(i)**2
               ctmpk=CC(k)*2.0*L0(i)*L0(j)  
               ctmpi=CC(i)*2.0*L0(k)*L0(j)         
               vec1(1:3)=dQ(k1,1:3)*sQ(k,i)
               Ftmp(i,1:3)=vec1(1:3)/L(j)*ktmpj*DL(j)
     &                    +vec1(1:3)/L(j)*(ctmpk*DL(i)+ctmpi*DL(k))
            enddo
         enddo 
       endif

       if(KtmpFlag) then
       Ktmp(1:3,1:3,1:3,1:3)=0.0
       if(Lin_NonlinMode == NonLinmode) then
         do im=1,3
         do jm=1,3
               km=IJ(jm,1)
               if(km.eq. im) km=IJ(jm,2)
               veca(1,1:3)=dQ(iQ(im,jm),1:3)*sQ(im,jm)  !Qi-Qj
               veca(2,1:3)=dQ(iQ(im,km),1:3)*sQ(im,km)  !Qi-Qk
               veca(3,1:3)=dQ(iQ(jm,im),1:3)*sQ(jm,im)  !Qj-Qi
               veca(4,1:3)=dQ(iQ(jm,km),1:3)*sQ(jm,km)  !Qj-Qk
               veca(5,1:3)=dQ(iQ(km,im),1:3)*sQ(km,im)  !Qk-Qi
               veca(6,1:3)=dQ(iQ(km,jm),1:3)*sQ(km,jm)  !Qk-Qj

               vec1(1)=KK(km)*D2L(km)+CC(jm)*D2L(im)+CC(im)*D2L(jm)
               do i=1,3
                  j=i
                  Ktmp(im,i,jm,j)=vec1(1)
                  do j=1,3
                     Ktmp(im,i,jm,j)=Ktmp(im,i,jm,j)
     &               +2.0*CC(km)*veca(2,i)*veca(6,j)    
     &               +2.0*CC(im)*veca(2,i)*veca(1,j) 
     &               +2.0*CC(jm)*veca(3,i)*veca(4,j) 
     &               +2.0*KK(km)*veca(1,i)*veca(1,j) 
                 enddo
               enddo
          enddo
          enddo
       elseif(Lin_NonlinMode == Linmode) then  
         do im=1,3
         do jm=1,3
               km=IJ(jm,1)
               if(km.eq. im) km=IJ(jm,2)
               ktmpk=KK(km)*2.0*L0(km)**2
               ctmpk=CC(km)*2.0*L0(im)*L0(jm)  
               ctmpi=CC(im)*2.0*L0(km)*L0(jm)         
               ctmpj=CC(jm)*2.0*L0(km)*L0(im) 

               veca(1,1:3)=dQ(iQ(im,jm),1:3)*sQ(im,jm)  !Qi-Qj
               veca(2,1:3)=dQ(iQ(im,km),1:3)*sQ(im,km)  !Qi-Qk
               veca(3,1:3)=dQ(iQ(jm,im),1:3)*sQ(jm,im)  !Qj-Qi
               veca(4,1:3)=dQ(iQ(jm,km),1:3)*sQ(jm,km)  !Qj-Qk
               veca(5,1:3)=dQ(iQ(km,im),1:3)*sQ(km,im)  !Qk-Qi
               veca(6,1:3)=dQ(iQ(km,jm),1:3)*sQ(km,jm)  !Qk-Qj

               vec1(1)=(ktmpk*DL(km)+ctmpj*DL(im)+ctmpi*DL(jm))/L(km)
               do i=1,3
                  j=i
                  Ktmp(im,i,jm,j)=vec1(1)
                  do j=1,3
                     Ktmp(im,i,jm,j)=Ktmp(im,i,jm,j)
     &               +ctmpk/(L(im)*L(jm))*veca(5,i)*veca(4,j)    
     &               +ctmpi/(L(km)*L(jm))*veca(2,i)*veca(1,j)  
     &               +ctmpj/(L(im)*L(km))*veca(3,i)*veca(4,j)  
     &               +(ktmpk*L0(km)-ctmpj*DL(im)-ctmpi*DL(jm))/L(km)**3
     &                                   *veca(1,i)*veca(1,j)
  
                 enddo
               enddo
          enddo
          enddo
        endif
         do im=1,3
         do jm=1,3
               do i=1,3
                  do j=1,3
                     Ktmp(im,i,jm,j)=cf_coef*Ktmp(im,i,jm,j)
                 enddo
               enddo
          enddo
          enddo
         do im=1,3
               do i=1,3
                     Ftmp(im,i)=cf_coef*Ftmp(im,i)
               enddo
         enddo
      endif
      return
      end subroutine
!======================================================================
!======================================================================
      subroutine Inplane_LineForce(Ftmp
     &                            ,x0,y0,z0
     &                            ,x1,y1,z1
     &                            ,L_over_Lmax,z3
     &                            ,InplaneLinemode,mcoef
     &                            ,cf_coef
     &                           ,Ktmp,KtmpFlag)

      real x0,y0,z0,x1,y1,z1,z3,cf_coef,cf_coef2,t0
      real Ftmp(2,3),Ktmp(2,3,2,3)
      real TmpVar_fbc1,TmpVar_fbc2,TmpVar_fbc3
      real TmpVar_fbc4,TmpVar_fbc5
      real TmpVar_fbc6,TmpVar_fbc7,TmpVar_fbc8
      real TmpVar_fbc9,TmpVar_fbc10
      real TmpVar_fbc11,TmpVar_fbc12,TmpVar_fbc13
      real TmpVar_fbc14,TmpVar_fbc15
      real TmpVar_fbc16,TmpVar_fbc17,TmpVar_fbc18
      real TmpVar_fbc19,TmpVar_fbc20
      logical KtmpFlag
      integer InplaneLinemode,mcoef
      integer WLCmode,FENEmode,POWmode
      PARAMETER(WLCmode=0,FENEmode=1,POWmode=2)

       Ftmp(1:2,1:3)=0.0
       if(InplaneLinemode == WLCmode) then
       maxtmpx=L_over_Lmax-0.25+0.25/((1.-L_over_Lmax)**2)   !WLC -A.1

       !!maxtmpx=z3*L_over_Lmax/(1.-L_over_Lmax**2)          !FENE

!!maxtmpx=1.0/(z3*L_over_Lmax)**m
!     POW (m is req to be defined)



       cf_coef2=cf_coef*maxtmpx/(z3*L_over_Lmax)
       Ftmp(1,1)=(x1-x0)*cf_coef2                                     
       Ftmp(1,2)=(y1-y0)*cf_coef2                                     
       Ftmp(1,3)=(z1-z0)*cf_coef2                                     
       Ftmp(2,1)=-Ftmp(1,1)                                   
       Ftmp(2,2)=-Ftmp(1,2)                                   
       Ftmp(2,3)=-Ftmp(1,3)

       if(KtmpFlag) then
       Ktmp(1:2,1:3,1:2,1:3)=0.0
      t2 = x0-x1
      t3 = y0-y1
      t4 = z0-z1
      t5 = t2**2
      t6 = t3**2
      t7 = t4**2
      t8 = t5+t6+t7
      t9 = sqrt(t8)
      t10 = 1.0D0/z3
      t11 = 1.0D0/sqrt(t8)
      t12 = z3**2
      t13 = t9-z3
      t14 = x0*2.0D0
      t27 = x1*2.0D0
      t15 = t14-t27
      t16 = t9*t10
      t17 = 1.0D0/t13**3
      t18 = y0*2.0D0
      t31 = y1*2.0D0
      t19 = t18-t31
      t20 = 1.0D0/(t9-z3)**2
      t21 = t12*t20*(1.0D0/4.0D0)
      t22 = t16+t21-1.0D0/4.0D0
      t23 = 1.0D0/t8**(3.0D0/2.0D0)
      t24 = z0*2.0D0
      t35 = z1*2.0D0
      t25 = t24-t35
      t26 = t2*t10*t11
      t38 = t11*t12*t15*t17*(1.0D0/4.0D0)
      t28 = t26-t38
      t29 = cf_coef*t2*t15*t22*t23*(1.0D0/2.0D0)
      t30 = t3*t10*t11
      t40 = t11*t12*t17*t19*(1.0D0/4.0D0)
      t32 = t30-t40
      t33 = cf_coef*t2*t19*t22*t23*(1.0D0/2.0D0)
      t34 = t4*t10*t11
      t41 = t11*t12*t17*t25*(1.0D0/4.0D0)
      t36 = t34-t41
      t37 = cf_coef*t2*t22*t23*t25*(1.0D0/2.0D0)
      t39 = cf_coef*t11*t22
      t42 = cf_coef*t3*t15*t22*t23*(1.0D0/2.0D0)
      t43 = cf_coef*t3*t19*t22*t23*(1.0D0/2.0D0)
      t44 = cf_coef*t3*t22*t23*t25*(1.0D0/2.0D0)
      t45 = cf_coef*t4*t15*t22*t23*(1.0D0/2.0D0)
      t46 = cf_coef*t4*t19*t22*t23*(1.0D0/2.0D0)
      t47 = cf_coef*t4*t22*t23*t25*(1.0D0/2.0D0)
      t48 = cf_coef*t2*t11*t28
      t49 = -t29+t39+t48
      t50 = cf_coef*t2*t11*t32
      t51 = -t33+t50
      t52 = cf_coef*t2*t11*t36
      t53 = -t37+t52
      t54 = cf_coef*t3*t11*t28
      t55 = -t42+t54
      t56 = cf_coef*t3*t11*t32
      t57 = t39-t43+t56
      t58 = cf_coef*t3*t11*t36
      t59 = -t44+t58
      t60 = cf_coef*t4*t11*t28
      t61 = -t45+t60
      t62 = cf_coef*t4*t11*t32
      t63 = -t46+t62
      t64 = cf_coef*t4*t11*t36
      t65 = t39-t47+t64
      Ktmp(1,1,1,1) = t29-cf_coef*t11
     &*(t16+t12*1.0D0/t13**2*(1.0D0/4.0D0)-
     &1.0D0/4.0D0)-cf_coef*t2*t11*t28
      Ktmp(1,1,1,2) = t33-cf_coef*t2*t11*t32
      Ktmp(1,1,1,3) = t37-cf_coef*t2*t11*t36
      Ktmp(1,1,2,1) = t49
      Ktmp(1,1,2,2) = t51
      Ktmp(1,1,2,3) = t53
      Ktmp(1,2,1,1) = t42-cf_coef*t3*t11*t28
      Ktmp(1,2,1,2) = -t39+t43-cf_coef*t3*t11*t32
      Ktmp(1,2,1,3) = t44-cf_coef*t3*t11*t36
      Ktmp(1,2,2,1) = t55
      Ktmp(1,2,2,2) = t57
      Ktmp(1,2,2,3) = t59
      Ktmp(1,3,1,1) = t45-cf_coef*t4*t11*t28
      Ktmp(1,3,1,2) = t46-cf_coef*t4*t11*t32
      Ktmp(1,3,1,3) = -t39+t47-cf_coef*t4*t11*t36
      Ktmp(1,3,2,1) = t61
      Ktmp(1,3,2,2) = t63
      Ktmp(1,3,2,3) = t65
      Ktmp(2,1,1,1) = t49
      Ktmp(2,1,1,2) = t51
      Ktmp(2,1,1,3) = t53
      Ktmp(2,1,2,1) = t29-t39-t48
      Ktmp(2,1,2,2) = t33-t50
      Ktmp(2,1,2,3) = t37-t52
      Ktmp(2,2,1,1) = t55
      Ktmp(2,2,1,2) = t57
      Ktmp(2,2,1,3) = t59
      Ktmp(2,2,2,1) = t42-t54
      Ktmp(2,2,2,2) = -t39+t43-t56
      Ktmp(2,2,2,3) = t44-t58
      Ktmp(2,3,1,1) = t61
      Ktmp(2,3,1,2) = t63
      Ktmp(2,3,1,3) = t65
      Ktmp(2,3,2,1) = t45-t60
      Ktmp(2,3,2,2) = t46-t62
      Ktmp(2,3,2,3) = -t39+t47-t64
       endif
      elseif(InplaneLinemode == FENEmode) then
       maxtmpx=z3*L_over_Lmax/(1.-L_over_Lmax**2)          !FENE

!!maxtmpx=1.0/(z3*L_over_Lmax)**m
!     POW (m is req to be defined)



       cf_coef2=cf_coef*maxtmpx/(z3*L_over_Lmax)
       Ftmp(1,1)=(x1-x0)*cf_coef2                                     
       Ftmp(1,2)=(y1-y0)*cf_coef2                                     
       Ftmp(1,3)=(z1-z0)*cf_coef2                                     
       Ftmp(2,1)=-Ftmp(1,1)                                   
       Ftmp(2,2)=-Ftmp(1,2)                                   
       Ftmp(2,3)=-Ftmp(1,3)

       if(KtmpFlag) then
       Ktmp(1:2,1:3,1:2,1:3)=0.0
      t2 = x0-x1
      t3 = y0-y1

      t4 = z0-z1
      t5 = 1.0D0/z3**2
      t6 = t2**2
      t7 = t3**2
      t8 = t4**2
      t9 = t6+t7+t8
      t10 = t5*t9
      t11 = t10-1.0D0
      t12 = 1.0D0/t11**2
      t13 = 1.0D0/t11
      t14 = cf_coef*t13
      t15 = x0*2.0D0
      t21 = x1*2.0D0
      t16 = t15-t21
      t17 = y0*2.0D0
      t22 = y1*2.0D0
      t18 = t17-t22
      t19 = z0*2.0D0
      t23 = z1*2.0D0
      t20 = t19-t23
      t24 = cf_coef*t2*t5*t12*t16
      t25 = -t14+t24
      t26 = cf_coef*t2*t5*t12*t18
      t27 = cf_coef*t2*t5*t12*t20
      t28 = cf_coef*t3*t5*t12*t16
      t29 = cf_coef*t3*t5*t12*t18
      t30 = -t14+t29
      t31 = cf_coef*t3*t5*t12*t20
      t32 = cf_coef*t4*t5*t12*t16
      t33 = cf_coef*t4*t5*t12*t18
      t34 = cf_coef*t4*t5*t12*t20
      t35 = -t14+t34
      Ktmp(1,1,1,1) = t14-cf_coef*t2*t5*t12*t16
      Ktmp(1,1,1,2) = -cf_coef*t2*t5*t12*t18
      Ktmp(1,1,1,3) = -cf_coef*t2*t5*t12*t20
      Ktmp(1,1,2,1) = t25
      Ktmp(1,1,2,2) = t26
      Ktmp(1,1,2,3) = t27
      Ktmp(1,2,1,1) = -cf_coef*t3*t5*t12*t16
      Ktmp(1,2,1,2) = t14-cf_coef*t3*t5*t12*t18
      Ktmp(1,2,1,3) = -cf_coef*t3*t5*t12*t20
      Ktmp(1,2,2,1) = t28
      Ktmp(1,2,2,2) = t30
      Ktmp(1,2,2,3) = t31
      Ktmp(1,3,1,1) = -cf_coef*t4*t5*t12*t16
      Ktmp(1,3,1,2) = -cf_coef*t4*t5*t12*t18
      Ktmp(1,3,1,3) = t14-cf_coef*t4*t5*t12*t20
      Ktmp(1,3,2,1) = t32
      Ktmp(1,3,2,2) = t33
      Ktmp(1,3,2,3) = t35
      Ktmp(2,1,1,1) = t25
      Ktmp(2,1,1,2) = t26
      Ktmp(2,1,1,3) = t27
      Ktmp(2,1,2,1) = t14-t24
      Ktmp(2,1,2,2) = -t26
      Ktmp(2,1,2,3) = -t27
      Ktmp(2,2,1,1) = t28
      Ktmp(2,2,1,2) = t30
      Ktmp(2,2,1,3) = t31
      Ktmp(2,2,2,1) = -t28
      Ktmp(2,2,2,2) = t14-t29
      Ktmp(2,2,2,3) = -t31
      Ktmp(2,3,1,1) = t32
      Ktmp(2,3,1,2) = t33
      Ktmp(2,3,1,3) = t35
      Ktmp(2,3,2,1) = -t32
      Ktmp(2,3,2,2) = -t33
      Ktmp(2,3,2,3) = t14-t34
       endif
      elseif(InplaneLinemode == POWmode) then
         maxtmpx=1.0/(z3*L_over_Lmax)**m
!     POW (m is req to be defined)

       cf_coef2=cf_coef*maxtmpx/(z3*L_over_Lmax)
       Ftmp(1,1)=(x1-x0)*cf_coef2                                     
       Ftmp(1,2)=(y1-y0)*cf_coef2                                     
       Ftmp(1,3)=(z1-z0)*cf_coef2                                     
       Ftmp(2,1)=-Ftmp(1,1)                                   
       Ftmp(2,2)=-Ftmp(1,2)                                   
       Ftmp(2,3)=-Ftmp(1,3)

       if(KtmpFlag) then
       Ktmp(1:2,1:3,1:2,1:3)=0.0
      t2 = x0-x1
      t3 = y0-y1
      t4 = z0-z1
      t5 = t2**2
      t6 = t3**2
      t7 = t4**2
      t8 = t5+t6+t7
      t11 = m*(1.0D0/2.0D0)
      t9 = -t11-3.0D0/2.0D0
      t10 = t8**t9
      t12 = y0*2.0D0
      t20 = y1*2.0D0
      t13 = t12-t20
      t14 = z0*2.0D0
      t23 = z1*2.0D0
      t15 = t14-t23
      t16 = x0*2.0D0
      t18 = x1*2.0D0
      t17 = t16-t18
      t19 = cf_coef*m*t2*t10*t17*(1.0D0/2.0D0)
      t21 = cf_coef*t2*t10*t13*(1.0D0/2.0D0)
      t22 = cf_coef*m*t2*t10*t13*(1.0D0/2.0D0)
      t24 = cf_coef*t2*t10*t15*(1.0D0/2.0D0)
      t25 = cf_coef*m*t2*t10*t15*(1.0D0/2.0D0)
      t26 = -t11-1.0D0/2.0D0
      t27 = t8**t26
      t28 = cf_coef*t27
      t29 = cf_coef*t3*t10*t17*(1.0D0/2.0D0)
      t30 = cf_coef*m*t3*t10*t17*(1.0D0/2.0D0)
      t31 = cf_coef*m*t3*t10*t13*(1.0D0/2.0D0)
      t32 = cf_coef*t3*t10*t15*(1.0D0/2.0D0)
      t33 = cf_coef*m*t3*t10*t15*(1.0D0/2.0D0)
      t34 = cf_coef*t4*t10*t17*(1.0D0/2.0D0)
      t35 = cf_coef*m*t4*t10*t17*(1.0D0/2.0D0)
      t36 = cf_coef*t4*t10*t13*(1.0D0/2.0D0)
      t37 = cf_coef*m*t4*t10*t13*(1.0D0/2.0D0)
      t38 = cf_coef*m*t4*t10*t15*(1.0D0/2.0D0)
      t39 = -t19+t28-cf_coef*t2*t10*t17*(1.0D0/2.0D0)
      t40 = -t21-t22
      t41 = -t24-t25
      t42 = cf_coef*t5*t10
      t43 = t21+t22
      t44 = t24+t25
      t45 = -t29-t30
      t46 = t28-t31-cf_coef*t3*t10*t13*(1.0D0/2.0D0)
      t47 = -t32-t33
      t48 = t29+t30
      t49 = cf_coef*t6*t10
      t50 = -t28+t31+t49
      t51 = t32+t33
      t52 = -t34-t35
      t53 = -t36-t37
      t54 = t28-t38-cf_coef*t4*t10*t15*(1.0D0/2.0D0)
      t55 = t34+t35
      t56 = t36+t37
      t57 = cf_coef*t7*t10
      t58 = -t28+t38+t57
      Ktmp(1,1,1,1) = t19+t42-cf_coef
     &*t8**(m*(-1.0D0/2.0D0)-1.0D0/2.0D0)
      Ktmp(1,1,1,2) = t43
      Ktmp(1,1,1,3) = t44
      Ktmp(1,1,2,1) = t39
      Ktmp(1,1,2,2) = t40
      Ktmp(1,1,2,3) = t41
      Ktmp(1,2,1,1) = t48
      Ktmp(1,2,1,2) = t50
      Ktmp(1,2,1,3) = t51
      Ktmp(1,2,2,1) = t45
      Ktmp(1,2,2,2) = t46
      Ktmp(1,2,2,3) = t47
      Ktmp(1,3,1,1) = t55
      Ktmp(1,3,1,2) = t56
      Ktmp(1,3,1,3) = t58
      Ktmp(1,3,2,1) = t52
      Ktmp(1,3,2,2) = t53
      Ktmp(1,3,2,3) = t54
      Ktmp(2,1,1,1) = t39
      Ktmp(2,1,1,2) = t40
      Ktmp(2,1,1,3) = t41
      Ktmp(2,1,2,1) = t19-t28+t42
      Ktmp(2,1,2,2) = t43
      Ktmp(2,1,2,3) = t44
      Ktmp(2,2,1,1) = t45
      Ktmp(2,2,1,2) = t46
      Ktmp(2,2,1,3) = t47
      Ktmp(2,2,2,1) = t48
      Ktmp(2,2,2,2) = t50
      Ktmp(2,2,2,3) = t51
      Ktmp(2,3,1,1) = t52
      Ktmp(2,3,1,2) = t53
      Ktmp(2,3,1,3) = t54
      Ktmp(2,3,2,1) = t55
      Ktmp(2,3,2,2) = t56
      Ktmp(2,3,2,3) = t58
       endif
      endif
      return
      end subroutine
!======================================================================
      subroutine inplane_Cq_fbcForce(Ftmp
     &                            ,x0,y0,z0
     &                            ,x1,y1,z1
     &                            ,x2,y2,z2,cf2
     &                           ,Ktmp,KtmpFlag)

      real x0,y0,z0,x1,y1,z1,x2,y2,z2,cf2,t0
      real Ftmp(3,3),Ktmp(3,3,3,3)
      real TmpVar_fbc1,TmpVar_fbc2,TmpVar_fbc3
      real TmpVar_fbc4,TmpVar_fbc5
      real TmpVar_fbc6,TmpVar_fbc7,TmpVar_fbc8
      real TmpVar_fbc9,TmpVar_fbc10
      real TmpVar_fbc11,TmpVar_fbc12,TmpVar_fbc13
      real TmpVar_fbc14,TmpVar_fbc15
      real TmpVar_fbc16,TmpVar_fbc17,TmpVar_fbc18
      real TmpVar_fbc19,TmpVar_fbc20
      logical KtmpFlag

       Ftmp(1:3,1:3)=0.0

      !diff(Atot,x0)*(Atot*4)
       Ftmp(0+1,1+0)=
     & cf2*(z1*x2*z2 + z1*x2*z0 + z1*x1*z2 - z1*x1*z0 - 
     & z0*x2*z2 + x1*z0*z2 + x1*y2*y1 + x1*y0*y2 - 
     & x1*y0*y1 + y2*y1*x2 - y2*y0*x2 + y0*y1*x2 +
     & z1*z1*x0 + z2*z2*x0 + y2*y2*x0 + y1*y1*x0 - 
     & 2*y2*y1*x0 - z1*z1*x2 - 2*z1*x0*z2 - 
     & x1*z2*z2 - x1*y2*y2 - y1*y1*x2)
   
    
      !diff(Atot,y0)*(Atot*4)
       Ftmp(0+1,1+1)=
     & cf2*(y1*z2*z1 + y1*z0*z2 - y1*z0*z1 + z2*z1*y2 - 
     & z2*z0*y2 + z1*z0*y2 + x1*y2*x2 + x1*x0*y2 + 
     & x1*y1*x2 - x1*y1*x0 - x0*y2*x2 + y1*x0*x2 - 
     & y1*z2*z2 - 2*y0*z2*z1 - z1*z1*y2 + 
     & y0*z2*z2 + z1*z1*y0 + x1*x1*y0 + y0*x2*x2 - 
     & x1*x1*y2 - 2*x1*y0*x2 - y1*x2*x2)
       

      !diff(Atot,z0)*(Atot*4)
       Ftmp(0+1,1+2)=
     & cf2*(y1*z2*y2 + y1*y0*z2 + y1*z1*y2 - y1*z1*y0 - 
     & y0*z2*y2 + z1*y0*y2 + z1*x2*x1 + z1*x0*x2 - 
     & z1*x0*x1 + x2*x1*z2 - x2*x0*z2 + x1*x0*z2 - 
     & y1*y1*z2 - 2*y1*z0*y2 - z1*y2*y2 + 
     & y1*y1*z0 + y2*y2*z0 + x2*x2*z0 + x1*x1*z0 - 
     & z1*x2*x2 - 2*z0*x2*x1 - x1*x1*z2)

   
      !diff(Atot,x1)*(Atot*4)
       Ftmp(1+1,1+0)=
     & cf2*(-z1*x2*z2 + z1*x2*z0 + z1*x0*z2 - z1*x0*z0 + 
     & z0*x2*z2 + z0*x0*z2 - y2*y1*x2 + y2*y1*x0 + 
     & y2*y0*x2 + y0*x0*y2 + y0*y1*x2 - y0*y1*x0 + 
     & x1*z2*z2 + x1*z0*z0 + x1*y2*y2 + x1*y0*y0 - 
     & z0*z0*x2 - 2*x1*z0*z2 - z2*z2*x0 
     & - 2*x1*y0*y2 - y2*y2*x0 - y0*y0*x2)
    
      !diff(Atot,y1)*(Atot*4)
       Ftmp(1+1,1+1)=
     & cf2*(-z2*z1*y2 + z2*z1*y0 + z2*z0*y2 + z0*y0*z2 + 
     & z0*z1*y2 - z0*z1*y0 - x1*y2*x2 + x1*y2*x0 + 
     & x1*y0*x2 - x1*y0*x0 + x0*y2*x2 + x0*y0*x2 - 
     & 2*y1*z2*z0 - z2*z2*y0 - z0*z0*y2 + y1*z2*z2 + 
     & y1*z0*z0 + y1*x2*x2 + y1*x0*x0 - x0*x0*y2 - 
     & 2*y1*x2*x0 - x2*x2*y0)


    
      !diff(Atot,z1)*(Atot*4)
       Ftmp(1+1,1+2)=
     & cf2*(-y1*z2*y2 + y1*y0*z2 + y1*z0*y2 - y1*z0*y0 + 
     & y0*z2*y2 + y0*z0*y2 - x2*x1*z2 + z0*x2*x1 + 
     & x2*x0*z2 + x0*z0*x2 + x1*x0*z2 - x0*x1*z0 - 
     & y0*y0*z2 - 2*z1*y0*y2 - z0*y2*y2 + z1*y2*y2 + 
     & z1*y0*y0 + z1*x2*x2 + z1*x0*x0 - 2*z1*x0*x2 - 
     & z0*x2*x2 - x0*x0*z2)

    
      !diff(Atot,x2)*(Atot*4)
       Ftmp(2+1,1+0)=
     & cf2*(-z1*x1*z2 + z1*x1*z0 + z1*x0*z2 + z1*x0*z0 + 
     & x1*z0*z2 - z0*x0*z2 - x1*y2*y1 + x1*y0*y2 + 
     & x1*y0*y1 + y2*y1*x0 - y0*x0*y2 + y0*y1*x0 + 
     & z1*z1*x2 + z0*z0*x2 + y1*y1*x2 + y0*y0*x2 - 
     & z1*z1*x0 - 2*z1*x2*z0 - x1*z0*z0 - 
     & x1*y0*y0 - y1*y1*x0 - 2*y0*y1*x2)
      

      !diff(Atot,y2)*(Atot*4)
       Ftmp(2+1,1+1)=
     & cf2*(-y1*z2*z1 + y1*z2*z0 + y1*z0*z1 + z2*z1*y0 - 
     & z0*y0*z2 + z0*z1*y0 - x1*y1*x2 + x1*y1*x0 + 
     & x1*y0*x2 + x1*y0*x0 + x0*y1*x2 - x0*y0*x2 - 
     & y1*z0*z0 - z1*z1*y0 - 2*z0*z1*y2 + z1*z1*y2 + 
     & z0*z0*y2 + x1*x1*y2 + x0*x0*y2 - x0*x0*y1 - 
     & x1*x1*y0 - 2*x1*y2*x0)


      !diff(Atot,z2)*(Atot*4)
       Ftmp(2+1,1+2)=
     & cf2*(-y1*z1*y2 + y1*z1*y0 + y1*z0*y2 + y1*z0*y0 + 
     & z1*y0*y2 - y0*z0*y2 - z1*x2*x1 + z1*x0*x2 + 
     & z1*x0*x1 + z0*x2*x1 - x0*z0*x2 + x0*x1*z0 - 
     & y1*y1*z0 - 2*y1*y0*z2 - z1*y0*y0 + y1*y1*z2 + 
     & y0*y0*z2 + x1*x1*z2 + x0*x0*z2 - z1*x0*x0 - 
     & x1*x1*z0 - 2*x1*x0*z2)

       if(KtmpFlag) then
       Ktmp(1:3,1:3,1:3,1:3)=0.0
      t2 = x1-x2
      t3 = y2**2
      t4 = z2**2
      t5 = y0*y1
      t6 = z0*z1
      t7 = y1**2
      t8 = z1**2
      t9 = x0*y1*2.0D0
      t10 = x2*y0
      t11 = x0*z1*2.0D0
      t12 = x2*z0
      t13 = y1-y2
      t14 = z1-z2
      t15 = x2*y1
      t16 = x2*y2
      t17 = z1*z2
      t18 = x2**2
      t19 = x1*y0*2.0D0
      t20 = x0*y2
      t21 = x1*y1
      t22 = x1*y2
      t23 = x0*x1
      t24 = x1**2
      t25 = y0*z1*2.0D0
      t26 = y2*z0
      t27 = x2*z1
      t28 = x2*z2
      t29 = y2*z1
      t30 = y2*z2
      t31 = x1*x2
      t32 = y1*y2
      t33 = x1*z0*2.0D0
      t34 = x0*z2
      t35 = x1*z1
      t36 = x1*z2
      t37 = y1*z0*2.0D0
      t38 = y0*z2
      t39 = y1*z1
      t40 = y1*z2
      t41 = x0*y1
      t42 = x2*y0*2.0D0
      t43 = x0*z1
      t44 = x2*z0*2.0D0
      t45 = x0-x2
      t46 = y0**2
      t47 = z0**2
      t48 = x0*y2*2.0D0
      t49 = x0*z2*2.0D0
      t50 = x1*y0
      t51 = z0*z2
      t52 = y0*z1
      t53 = y2*z0*2.0D0
      t54 = y0-y2
      t55 = z0-z2
      t56 = x0*y0
      t57 = x0**2
      t58 = y0*z2*2.0D0
      t59 = x1*z0
      t60 = y1*z0
      t61 = x0*x2
      t62 = y0*y2
      t63 = x0*z0
      t64 = y0*z0
      t65 = x2*y1*2.0D0
      t66 = x2*z1*2.0D0
      t67 = t5+t6-t17-t32-t46-t47+t51+t62
      t68 = cf2*t67
      t69 = -t9+t20-t22-t42+t50+t56+t65
      t70 = cf2*t69
      t71 = -t11+t34-t36-t44+t59+t63+t66
      t72 = cf2*t71
      t73 = x0-x1
      t74 = x1*y2*2.0D0
      t75 = y2*z1*2.0D0
      t76 = t10-t15-t19+t41-t48+t56+t74
      t77 = cf2*t76
      t78 = t6-t17+t23-t31-t47+t51-t57+t61
      t79 = cf2*t78
      t80 = -t25+t38-t40-t53+t60+t64+t75
      t81 = cf2*t80
      t82 = y0-y1
      t83 = z0-z1
      t84 = x1*z2*2.0D0
      t85 = y1*z2*2.0D0
      t86 = t12-t27-t33+t43-t49+t63+t84
      t87 = cf2*t86
      t88 = t26-t29-t37+t52-t58+t64+t85
      t89 = cf2*t88
      t90 = t5+t23-t31-t32-t46-t57+t61+t62
      t91 = cf2*t90
      Ktmp(1,1,1,1) = cf2*(t3+t4+t7+t8-y1*y2*2.0D0-z1*z2*2.0D0)
      Ktmp(1,1,1,2) = -cf2*t2*t13
      Ktmp(1,1,1,3) = -cf2*t2*t14
      Ktmp(1,1,2,1) = -cf2*(t3+t4+t5+t6-y0*y2-y1*y2-z0*z2-z1*z2)
      Ktmp(1,1,2,2) = cf2
     &*(t9+t10+t16+t22-x1*y0-x0*y2*2.0D0-x2*y1*2.0D0)
      Ktmp(1,1,2,3) = cf2
     &*(t11+t12+t28+t36-x1*z0-x0*z2*2.0D0-x2*z1*2.0D0)
      Ktmp(1,1,3,1) = cf2*(t5+t6-t7-t8+t17+t32-y0*y2-z0*z2)
      Ktmp(1,1,3,2) = cf2*(-t9-t10+t15+t21+t48+t50-x1*y2*2.0D0)
      Ktmp(1,1,3,3) = cf2*(-t11-t12+t27+t35+t49+t59-x1*z2*2.0D0)
      Ktmp(1,2,1,1) = -cf2*t2*t13
      Ktmp(1,2,1,2) = cf2*(t4+t8+t18+t24-x1*x2*2.0D0-z1*z2*2.0D0)
      Ktmp(1,2,1,3) = -cf2*t13*t14
      Ktmp(1,2,2,1) = cf2
     &*(t15+t16+t19+t20-x0*y1-x2*y0*2.0D0-x1*y2*2.0D0)
      Ktmp(1,2,2,2) = -cf2
     &*(t4+t6-t17+t18+t23-x0*x2-x1*x2-z0*z2)
      Ktmp(1,2,2,3) = cf2
     &*(t25+t26+t30+t40-y1*z0-y0*z2*2.0D0-y2*z1*2.0D0)
      Ktmp(1,2,3,1) = cf2*(-t19-t20+t21+t22+t41+t42-x2*y1*2.0D0)
      Ktmp(1,2,3,2) = cf2*(t6-t8+t17+t23-t24+t31-x0*x2-z0*z2)
      Ktmp(1,2,3,3) = cf2*(-t25-t26+t29+t39+t58+t60-y1*z2*2.0D0)
      Ktmp(1,3,1,1) = -cf2*t2*t14
      Ktmp(1,3,1,2) = -cf2*t13*t14
      Ktmp(1,3,1,3) = cf2*(t3+t7+t18+t24-x1*x2*2.0D0-y1*y2*2.0D0)
      Ktmp(1,3,2,1) = cf2
     &*(t27+t28+t33+t34-x0*z1-x2*z0*2.0D0-x1*z2*2.0D0)
      Ktmp(1,3,2,2) = cf2
     &*(t29+t30+t37+t38-y0*z1-y2*z0*2.0D0-y1*z2*2.0D0)
      Ktmp(1,3,2,3) = -cf2*(t3+t5+t18+t23-t31-t32-x0*x2-y0*y2)
      Ktmp(1,3,3,1) = cf2*(-t33-t34+t35+t36+t43+t44-x2*z1*2.0D0)
      Ktmp(1,3,3,2) = cf2*(-t37-t38+t39+t40+t52+t53-y2*z1*2.0D0)
      Ktmp(1,3,3,3) = cf2*(t5-t7+t23-t24+t31+t32-x0*x2-y0*y2)
      Ktmp(2,1,1,1) = -cf2*(t3+t4+t5+t6-t17-t32-y0*y2-z0*z2)
      Ktmp(2,1,1,2) = cf2*(t15+t16+t19+t20-t41-t42-x1*y2*2.0D0)
      Ktmp(2,1,1,3) = cf2*(t27+t28+t33+t34-t43-t44-x1*z2*2.0D0)
      Ktmp(2,1,2,1) = cf2*(t3+t4+t46+t47-y0*y2*2.0D0-z0*z2*2.0D0)
      Ktmp(2,1,2,2) = -cf2*t45*t54
      Ktmp(2,1,2,3) = -cf2*t45*t55
      Ktmp(2,1,3,1) = t68
      Ktmp(2,1,3,2) = t77
      Ktmp(2,1,3,3) = t87
      Ktmp(2,2,1,1) = cf2*(t9+t10+t16+t22-t48-t50-x2*y1*2.0D0)
      Ktmp(2,2,1,2) = -cf2*(t4+t6-t17+t18+t23-t31-t51-x0*x2)
      Ktmp(2,2,1,3) = cf2*(t29+t30+t37+t38-t52-t53-y1*z2*2.0D0)
      Ktmp(2,2,2,1) = -cf2*t45*t54
      Ktmp(2,2,2,2) = cf2*(t4+t18+t47+t57-x0*x2*2.0D0-z0*z2*2.0D0)
      Ktmp(2,2,2,3) = -cf2*t54*t55
      Ktmp(2,2,3,1) = t70
      Ktmp(2,2,3,2) = t79
      Ktmp(2,2,3,3) = t89
      Ktmp(2,3,1,1) = cf2*(t11+t12+t28+t36-t49-t59-x2*z1*2.0D0)
      Ktmp(2,3,1,2) = cf2*(t25+t26+t30+t40-t58-t60-y2*z1*2.0D0)
      Ktmp(2,3,1,3) = -cf2*(t3+t5+t18+t23-t31-t32-t61-t62)
      Ktmp(2,3,2,1) = -cf2*t45*t55
      Ktmp(2,3,2,2) = -cf2*t54*t55
      Ktmp(2,3,2,3) = cf2*(t3+t18+t46+t57-x0*x2*2.0D0-y0*y2*2.0D0)
      Ktmp(2,3,3,1) = t72
      Ktmp(2,3,3,2) = t81
      Ktmp(2,3,3,3) = t91
      Ktmp(3,1,1,1) = cf2*(t5+t6-t7-t8+t17+t32-t51-t62)
      Ktmp(3,1,1,2) = cf2*(-t19-t20+t21+t22+t41+t42-t65)
      Ktmp(3,1,1,3) = cf2*(-t33-t34+t35+t36+t43+t44-t66)
      Ktmp(3,1,2,1) = t68
      Ktmp(3,1,2,2) = t70
      Ktmp(3,1,2,3) = t72
      Ktmp(3,1,3,1) = cf2*(t7+t8+t46+t47-y0*y1*2.0D0-z0*z1*2.0D0)
      Ktmp(3,1,3,2) = -cf2*t73*t82
      Ktmp(3,1,3,3) = -cf2*t73*t83
      Ktmp(3,2,1,1) = cf2*(-t9-t10+t15+t21+t48+t50-t74)
      Ktmp(3,2,1,2) = cf2*(t6-t8+t17+t23-t24+t31-t51-t61)
      Ktmp(3,2,1,3) = cf2*(-t37-t38+t39+t40+t52+t53-t75)
      Ktmp(3,2,2,1) = t77
      Ktmp(3,2,2,2) = t79
      Ktmp(3,2,2,3) = t81
      Ktmp(3,2,3,1) = -cf2*t73*t82
      Ktmp(3,2,3,2) = cf2*(t8+t24+t47+t57-x0*x1*2.0D0-z0*z1*2.0D0)
      Ktmp(3,2,3,3) = -cf2*t82*t83
      Ktmp(3,3,1,1) = cf2*(-t11-t12+t27+t35+t49+t59-t84)
      Ktmp(3,3,1,2) = cf2*(-t25-t26+t29+t39+t58+t60-t85)
      Ktmp(3,3,1,3) = cf2*(t5-t7+t23-t24+t31+t32-t61-t62)
      Ktmp(3,3,2,1) = t87
      Ktmp(3,3,2,2) = t89
      Ktmp(3,3,2,3) = t91
      Ktmp(3,3,3,1) = -cf2*t73*t83
      Ktmp(3,3,3,2) = -cf2*t82*t83
      Ktmp(3,3,3,3) = cf2*(t7+t24+t46+t57-x0*x1*2.0D0-y0*y1*2.0D0)
       endif
      return
      end subroutine
!======================================================================
      subroutine Local_AreaForce(Ftmp
     &                            ,x0,y0,z0
     &                            ,x1,y1,z1
     &                            ,x2,y2,z2,cf2
     &                           ,Ktmp2,KtmpFlag)

      real x0,y0,z0,x1,y1,z1,x2,y2,z2,cf2,t0
      real Ftmp(3,3),Ktmp(3,3,3,3),Ktmp2(3,3,3,3)
      real TmpVar_fbc1,TmpVar_fbc2,TmpVar_fbc3
      real TmpVar_fbc4,TmpVar_fbc5
      real TmpVar_fbc6,TmpVar_fbc7,TmpVar_fbc8
      real TmpVar_fbc9,TmpVar_fbc10
      real TmpVar_fbc11,TmpVar_fbc12,TmpVar_fbc13
      real TmpVar_fbc14,TmpVar_fbc15
      real TmpVar_fbc16,TmpVar_fbc17,TmpVar_fbc18
      real TmpVar_fbc19,TmpVar_fbc20
      logical KtmpFlag

       Ftmp(1:3,1:3)=0.0
       !diff(Atot,x0)*(Atot*4)
       Ftmp(1,1)=Ftmp(1,1)+
     &   cf2*(z1*x2*z2 + z1*x2*z0 + z1*x1*z2 - z1*x1*z0 -
     &     z0*x2*z2 + x1*z0*z2 + x1*y2*y1 + x1*y0*y2 -
     &     x1*y0*y1 + y2*y1*x2 - y2*y0*x2 + y0*y1*x2 +
     &     z1*z1*x0 + z2*z2*x0 + y2*y2*x0 + y1*y1*x0 -
     &     2*y2*y1*x0 - z1*z1*x2 - 2*z1*x0*z2 -
     &     x1*z2*z2 - x1*y2*y2 - y1*y1*x2)
   
      !diff(Atot,x0)*(Atot*4)
       Ftmp(0+1,1+0)=
     & cf2*(z1*x2*z2 + z1*x2*z0 + z1*x1*z2 - z1*x1*z0 - 
     & z0*x2*z2 + x1*z0*z2 + x1*y2*y1 + x1*y0*y2 - 
     & x1*y0*y1 + y2*y1*x2 - y2*y0*x2 + y0*y1*x2 +
     & z1*z1*x0 + z2*z2*x0 + y2*y2*x0 + y1*y1*x0 - 
     & 2*y2*y1*x0 - z1*z1*x2 - 2*z1*x0*z2 - 
     & x1*z2*z2 - x1*y2*y2 - y1*y1*x2)
    
       !diff(Atot,y0)*(Atot*4)
       Ftmp(1,2)=Ftmp(1,2)+
     &   cf2*(y1*z2*z1 + y1*z0*z2 - y1*z0*z1 + z2*z1*y2 -
     &     z2*z0*y2 + z1*z0*y2 + x1*y2*x2 + x1*x0*y2 +
     &     x1*y1*x2 - x1*y1*x0 - x0*y2*x2 + y1*x0*x2 -
     &     y1*z2*z2 - 2*y0*z2*z1 - z1*z1*y2 +
     &     y0*z2*z2 + z1*z1*y0 + x1*x1*y0 + y0*x2*x2 -
     &     x1*x1*y2 - 2*x1*y0*x2 - y1*x2*x2)
       

       !diff(Atot,z0)*(Atot*4)
       Ftmp(1,3)=Ftmp(1,3)+
     &   cf2*(y1*z2*y2 + y1*y0*z2 + y1*z1*y2 - y1*z1*y0 -
     &     y0*z2*y2 + z1*y0*y2 + z1*x2*x1 + z1*x0*x2 -
     &     z1*x0*x1 + x2*x1*z2 - x2*x0*z2 + x1*x0*z2 -
     &     y1*y1*z2 - 2*y1*z0*y2 - z1*y2*y2 +
     &     y1*y1*z0 + y2*y2*z0 + x2*x2*z0 + x1*x1*z0 -
     &     z1*x2*x2 - 2*z0*x2*x1 - x1*x1*z2)

   
       !diff(Atot,x1)*(Atot*4)
       Ftmp(2,1)=Ftmp(2,1)+
     &   cf2*(-z1*x2*z2 + z1*x2*z0 + z1*x0*z2 - z1*x0*z0 +
     &     z0*x2*z2 + z0*x0*z2 - y2*y1*x2 + y2*y1*x0 +
     &     y2*y0*x2 + y0*x0*y2 + y0*y1*x2 - y0*y1*x0 +
     &     x1*z2*z2 + x1*z0*z0 + x1*y2*y2 + x1*y0*y0 -
     &     z0*z0*x2 - 2*x1*z0*z2 - z2*z2*x0
     &     - 2*x1*y0*y2 - y2*y2*x0 - y0*y0*x2)
    
      !diff(Atot,y1)*(Atot*4)
       Ftmp(2,2)=Ftmp(2,2)+
     &   cf2*(-z2*z1*y2 + z2*z1*y0 + z2*z0*y2 + z0*y0*z2 +
     &     z0*z1*y2 - z0*z1*y0 - x1*y2*x2 + x1*y2*x0 +
     &     x1*y0*x2 - x1*y0*x0 + x0*y2*x2 + x0*y0*x2 -
     &     2*y1*z2*z0 - z2*z2*y0 - z0*z0*y2 + y1*z2*z2 +
     &     y1*z0*z0 + y1*x2*x2 + y1*x0*x0 - x0*x0*y2 -
     &     2*y1*x2*x0 - x2*x2*y0)


    
       !diff(Atot,z1)*(Atot*4)
       Ftmp(2,3)=Ftmp(2,3)+
     &   cf2*(-y1*z2*y2 + y1*y0*z2 + y1*z0*y2 - y1*z0*y0 +
     &     y0*z2*y2 + y0*z0*y2 - x2*x1*z2 + z0*x2*x1 +
     &     x2*x0*z2 + x0*z0*x2 + x1*x0*z2 - x0*x1*z0 -
     &     y0*y0*z2 - 2*z1*y0*y2 - z0*y2*y2 + z1*y2*y2 +
     &     z1*y0*y0 + z1*x2*x2 + z1*x0*x0 - 2*z1*x0*x2 -
     &     z0*x2*x2 - x0*x0*z2)

    
       !diff(Atot,x2)*(Atot*4)
       Ftmp(3,1)=Ftmp(3,1)+
     &   cf2*(-z1*x1*z2 + z1*x1*z0 + z1*x0*z2 + z1*x0*z0 +
     &     x1*z0*z2 - z0*x0*z2 - x1*y2*y1 + x1*y0*y2 +
     &     x1*y0*y1 + y2*y1*x0 - y0*x0*y2 + y0*y1*x0 +
     &     z1*z1*x2 + z0*z0*x2 + y1*y1*x2 + y0*y0*x2 -
     &     z1*z1*x0 - 2*z1*x2*z0 - x1*z0*z0 -
     &     x1*y0*y0 - y1*y1*x0 - 2*y0*y1*x2)
      

       !diff(Atot,y2)*(Atot*4)
       Ftmp(3,2)=Ftmp(3,2)+
     &   cf2*(-y1*z2*z1 + y1*z2*z0 + y1*z0*z1 + z2*z1*y0 -
     &     z0*y0*z2 + z0*z1*y0 - x1*y1*x2 + x1*y1*x0 +
     &     x1*y0*x2 + x1*y0*x0 + x0*y1*x2 - x0*y0*x2 -
     &     y1*z0*z0 - z1*z1*y0 - 2*z0*z1*y2 + z1*z1*y2 +
     &     z0*z0*y2 + x1*x1*y2 + x0*x0*y2 - x0*x0*y1 -
     &     x1*x1*y0 - 2*x1*y2*x0)


       !diff(Atot,z2)*(Atot*4)
       Ftmp(3,3)=Ftmp(3,3)+
     &   cf2*(-y1*z1*y2 + y1*z1*y0 + y1*z0*y2 + y1*z0*y0 +
     &     z1*y0*y2 - y0*z0*y2 - z1*x2*x1 + z1*x0*x2 +
     &     z1*x0*x1 + z0*x2*x1 - x0*z0*x2 + x0*x1*z0 -
     &     y1*y1*z0 - 2*y1*y0*z2 - z1*y0*y0 + y1*y1*z2 +
     &     y0*y0*z2 + x1*x1*z2 + x0*x0*z2 - z1*x0*x0 -
     &     x1*x1*z0 - 2*x1*x0*z2)
       if(KtmpFlag) then
       Ktmp(1:3,1:3,1:3,1:3)=0.0
      t2 = x1-x2
      t3 = y2**2
      t4 = z2**2
      t5 = y0*y1
      t6 = z0*z1
      t7 = y1**2
      t8 = z1**2
      t9 = x0*y1*2.0D0
      t10 = x2*y0
      t11 = x0*z1*2.0D0
      t12 = x2*z0
      t13 = y1-y2
      t14 = z1-z2
      t15 = x2*y1
      t16 = x2*y2
      t17 = z1*z2
      t18 = x2**2
      t19 = x1*y0*2.0D0
      t20 = x0*y2
      t21 = x1*y1
      t22 = x1*y2
      t23 = x0*x1
      t24 = x1**2
      t25 = y0*z1*2.0D0
      t26 = y2*z0
      t27 = x2*z1
      t28 = x2*z2
      t29 = y2*z1
      t30 = y2*z2
      t31 = x1*x2
      t32 = y1*y2
      t33 = x1*z0*2.0D0
      t34 = x0*z2
      t35 = x1*z1
      t36 = x1*z2
      t37 = y1*z0*2.0D0
      t38 = y0*z2
      t39 = y1*z1
      t40 = y1*z2
      t41 = x0*y1
      t42 = x2*y0*2.0D0
      t43 = x0*z1
      t44 = x2*z0*2.0D0
      t45 = x0-x2
      t46 = y0**2
      t47 = z0**2
      t48 = x0*y2*2.0D0
      t49 = x0*z2*2.0D0
      t50 = x1*y0
      t51 = z0*z2
      t52 = y0*z1
      t53 = y2*z0*2.0D0
      t54 = y0-y2
      t55 = z0-z2
      t56 = x0*y0
      t57 = x0**2
      t58 = y0*z2*2.0D0
      t59 = x1*z0
      t60 = y1*z0
      t61 = x0*x2
      t62 = y0*y2
      t63 = x0*z0
      t64 = y0*z0
      t65 = x2*y1*2.0D0
      t66 = x2*z1*2.0D0
      t67 = t5+t6-t17-t32-t46-t47+t51+t62
      t68 = cf2*t67
      t69 = -t9+t20-t22-t42+t50+t56+t65
      t70 = cf2*t69
      t71 = -t11+t34-t36-t44+t59+t63+t66
      t72 = cf2*t71
      t73 = x0-x1
      t74 = x1*y2*2.0D0
      t75 = y2*z1*2.0D0
      t76 = t10-t15-t19+t41-t48+t56+t74
      t77 = cf2*t76
      t78 = t6-t17+t23-t31-t47+t51-t57+t61
      t79 = cf2*t78
      t80 = -t25+t38-t40-t53+t60+t64+t75
      t81 = cf2*t80
      t82 = y0-y1
      t83 = z0-z1
      t84 = x1*z2*2.0D0
      t85 = y1*z2*2.0D0
      t86 = t12-t27-t33+t43-t49+t63+t84
      t87 = cf2*t86
      t88 = t26-t29-t37+t52-t58+t64+t85
      t89 = cf2*t88
      t90 = t5+t23-t31-t32-t46-t57+t61+t62
      t91 = cf2*t90
      Ktmp(1,1,1,1) = cf2*(t3+t4+t7+t8-y1*y2*2.0D0-z1*z2*2.0D0)
      Ktmp(1,1,1,2) = -cf2*t2*t13
      Ktmp(1,1,1,3) = -cf2*t2*t14
      Ktmp(1,1,2,1) = -cf2*(t3+t4+t5+t6-y0*y2-y1*y2-z0*z2-z1*z2)
      Ktmp(1,1,2,2) = cf2
     &   *(t9+t10+t16+t22-x1*y0-x0*y2*2.0D0-x2*y1*2.0D0)
      Ktmp(1,1,2,3) = cf2
     &   *(t11+t12+t28+t36-x1*z0-x0*z2*2.0D0-x2*z1*2.0D0)
      Ktmp(1,1,3,1) = cf2*(t5+t6-t7-t8+t17+t32-y0*y2-z0*z2)
      Ktmp(1,1,3,2) = cf2*(-t9-t10+t15+t21+t48+t50-x1*y2*2.0D0)
      Ktmp(1,1,3,3) = cf2*(-t11-t12+t27+t35+t49+t59-x1*z2*2.0D0)
      Ktmp(1,2,1,1) = -cf2*t2*t13
      Ktmp(1,2,1,2) = cf2*(t4+t8+t18+t24-x1*x2*2.0D0-z1*z2*2.0D0)
      Ktmp(1,2,1,3) = -cf2*t13*t14
      Ktmp(1,2,2,1) = cf2
     &   *(t15+t16+t19+t20-x0*y1-x2*y0*2.0D0-x1*y2*2.0D0)
      Ktmp(1,2,2,2) = -cf2*(t4+t6-t17+t18+t23-x0*x2-x1*x2-z0*z2)
      Ktmp(1,2,2,3) = cf2
     &   *(t25+t26+t30+t40-y1*z0-y0*z2*2.0D0-y2*z1*2.0D0)
      Ktmp(1,2,3,1) = cf2*(-t19-t20+t21+t22+t41+t42-x2*y1*2.0D0)
      Ktmp(1,2,3,2) = cf2*(t6-t8+t17+t23-t24+t31-x0*x2-z0*z2)
      Ktmp(1,2,3,3) = cf2*(-t25-t26+t29+t39+t58+t60-y1*z2*2.0D0)
      Ktmp(1,3,1,1) = -cf2*t2*t14
      Ktmp(1,3,1,2) = -cf2*t13*t14
      Ktmp(1,3,1,3) = cf2*(t3+t7+t18+t24-x1*x2*2.0D0-y1*y2*2.0D0)
      Ktmp(1,3,2,1) = cf2
     &   *(t27+t28+t33+t34-x0*z1-x2*z0*2.0D0-x1*z2*2.0D0)
      Ktmp(1,3,2,2) = cf2
     &   *(t29+t30+t37+t38-y0*z1-y2*z0*2.0D0-y1*z2*2.0D0)
      Ktmp(1,3,2,3) = -cf2*(t3+t5+t18+t23-t31-t32-x0*x2-y0*y2)
      Ktmp(1,3,3,1) = cf2*(-t33-t34+t35+t36+t43+t44-x2*z1*2.0D0)
      Ktmp(1,3,3,2) = cf2*(-t37-t38+t39+t40+t52+t53-y2*z1*2.0D0)
      Ktmp(1,3,3,3) = cf2*(t5-t7+t23-t24+t31+t32-x0*x2-y0*y2)
      Ktmp(2,1,1,1) = -cf2*(t3+t4+t5+t6-t17-t32-y0*y2-z0*z2)
      Ktmp(2,1,1,2) = cf2*(t15+t16+t19+t20-t41-t42-x1*y2*2.0D0)
      Ktmp(2,1,1,3) = cf2*(t27+t28+t33+t34-t43-t44-x1*z2*2.0D0)
      Ktmp(2,1,2,1) = cf2*(t3+t4+t46+t47-y0*y2*2.0D0-z0*z2*2.0D0)
      Ktmp(2,1,2,2) = -cf2*t45*t54
      Ktmp(2,1,2,3) = -cf2*t45*t55
      Ktmp(2,1,3,1) = t68
      Ktmp(2,1,3,2) = t77
      Ktmp(2,1,3,3) = t87
      Ktmp(2,2,1,1) = cf2*(t9+t10+t16+t22-t48-t50-x2*y1*2.0D0)
      Ktmp(2,2,1,2) = -cf2*(t4+t6-t17+t18+t23-t31-t51-x0*x2)
      Ktmp(2,2,1,3) = cf2*(t29+t30+t37+t38-t52-t53-y1*z2*2.0D0)
      Ktmp(2,2,2,1) = -cf2*t45*t54
      Ktmp(2,2,2,2) = cf2*(t4+t18+t47+t57-x0*x2*2.0D0-z0*z2*2.0D0)
      Ktmp(2,2,2,3) = -cf2*t54*t55
      Ktmp(2,2,3,1) = t70
      Ktmp(2,2,3,2) = t79
      Ktmp(2,2,3,3) = t89
      Ktmp(2,3,1,1) = cf2*(t11+t12+t28+t36-t49-t59-x2*z1*2.0D0)
      Ktmp(2,3,1,2) = cf2*(t25+t26+t30+t40-t58-t60-y2*z1*2.0D0)
      Ktmp(2,3,1,3) = -cf2*(t3+t5+t18+t23-t31-t32-t61-t62)
      Ktmp(2,3,2,1) = -cf2*t45*t55
      Ktmp(2,3,2,2) = -cf2*t54*t55
      Ktmp(2,3,2,3) = cf2*(t3+t18+t46+t57-x0*x2*2.0D0-y0*y2*2.0D0)
      Ktmp(2,3,3,1) = t72
      Ktmp(2,3,3,2) = t81
      Ktmp(2,3,3,3) = t91
      Ktmp(3,1,1,1) = cf2*(t5+t6-t7-t8+t17+t32-t51-t62)
      Ktmp(3,1,1,2) = cf2*(-t19-t20+t21+t22+t41+t42-t65)
      Ktmp(3,1,1,3) = cf2*(-t33-t34+t35+t36+t43+t44-t66)
      Ktmp(3,1,2,1) = t68
      Ktmp(3,1,2,2) = t70
      Ktmp(3,1,2,3) = t72
      Ktmp(3,1,3,1) = cf2*(t7+t8+t46+t47-y0*y1*2.0D0-z0*z1*2.0D0)
      Ktmp(3,1,3,2) = -cf2*t73*t82
      Ktmp(3,1,3,3) = -cf2*t73*t83
      Ktmp(3,2,1,1) = cf2*(-t9-t10+t15+t21+t48+t50-t74)
      Ktmp(3,2,1,2) = cf2*(t6-t8+t17+t23-t24+t31-t51-t61)
      Ktmp(3,2,1,3) = cf2*(-t37-t38+t39+t40+t52+t53-t75)
      Ktmp(3,2,2,1) = t77
      Ktmp(3,2,2,2) = t79
      Ktmp(3,2,2,3) = t81
      Ktmp(3,2,3,1) = -cf2*t73*t82
      Ktmp(3,2,3,2) = cf2*(t8+t24+t47+t57-x0*x1*2.0D0-z0*z1*2.0D0)
      Ktmp(3,2,3,3) = -cf2*t82*t83
      Ktmp(3,3,1,1) = cf2*(-t11-t12+t27+t35+t49+t59-t84)
      Ktmp(3,3,1,2) = cf2*(-t25-t26+t29+t39+t58+t60-t85)
      Ktmp(3,3,1,3) = cf2*(t5-t7+t23-t24+t31+t32-t61-t62)
      Ktmp(3,3,2,1) = t87
      Ktmp(3,3,2,2) = t89
      Ktmp(3,3,2,3) = t91
      Ktmp(3,3,3,1) = -cf2*t73*t83
      Ktmp(3,3,3,2) = -cf2*t82*t83
      Ktmp(3,3,3,3) = cf2*(t7+t24+t46+t57-x0*x1*2.0D0-y0*y1*2.0D0)

       Ktmp2(1:3,1:3,1:3,1:3)=Ktmp2(1:3,1:3,1:3,1:3)
     &                       +Ktmp(1:3,1:3,1:3,1:3)
       endif

      return
      end subroutine
!======================================================================
      subroutine Total_AreaForce(Ftmp
     &                            ,x0,y0,z0
     &                            ,x1,y1,z1
     &                            ,x2,y2,z2,cf2
     &                           ,Ktmp,KtmpFlag)

      real x0,y0,z0,x1,y1,z1,x2,y2,z2,cf2,t0
      real Ftmp(3,3),ktmp(3,3,3,3)
      real TmpVar_fbc1,TmpVar_fbc2,TmpVar_fbc3
      real TmpVar_fbc4,TmpVar_fbc5
      real TmpVar_fbc6,TmpVar_fbc7,TmpVar_fbc8
      real TmpVar_fbc9,TmpVar_fbc10
      real TmpVar_fbc11,TmpVar_fbc12,TmpVar_fbc13
      real TmpVar_fbc14,TmpVar_fbc15
      real TmpVar_fbc16,TmpVar_fbc17,TmpVar_fbc18
      real TmpVar_fbc19,TmpVar_fbc20
      logical Ktmpflag

       Ftmp(1:3,1:3)=0.0
      !diff(Atot,x0)*(Atot*4)
       Ftmp(1,1)=
     &   cf2*(z1*x2*z2 + z1*x2*z0 + z1*x1*z2 - z1*x1*z0 -
     &     z0*x2*z2 + x1*z0*z2 + x1*y2*y1 + x1*y0*y2 -
     &     x1*y0*y1 + y2*y1*x2 - y2*y0*x2 + y0*y1*x2 +
     &     z1*z1*x0 + z2*z2*x0 + y2*y2*x0 + y1*y1*x0 -
     &     2*y2*y1*x0 - z1*z1*x2 - 2*z1*x0*z2 -
     &     x1*z2*z2 - x1*y2*y2 - y1*y1*x2)
   
    
      !diff(Atot,y0)*(Atot*4)
       Ftmp(1,2)=
     &   cf2*(y1*z2*z1 + y1*z0*z2 - y1*z0*z1 + z2*z1*y2 -
     &     z2*z0*y2 + z1*z0*y2 + x1*y2*x2 + x1*x0*y2 +
     &     x1*y1*x2 - x1*y1*x0 - x0*y2*x2 + y1*x0*x2 -
     &     y1*z2*z2 - 2*y0*z2*z1 - z1*z1*y2 +
     &     y0*z2*z2 + z1*z1*y0 + x1*x1*y0 + y0*x2*x2 -
     &     x1*x1*y2 - 2*x1*y0*x2 - y1*x2*x2)
       

      !diff(Atot,z0)*(Atot*4)
       Ftmp(1,3)=
     &   cf2*(y1*z2*y2 + y1*y0*z2 + y1*z1*y2 - y1*z1*y0 -
     &     y0*z2*y2 + z1*y0*y2 + z1*x2*x1 + z1*x0*x2 -
     &     z1*x0*x1 + x2*x1*z2 - x2*x0*z2 + x1*x0*z2 -
     &     y1*y1*z2 - 2*y1*z0*y2 - z1*y2*y2 +
     &     y1*y1*z0 + y2*y2*z0 + x2*x2*z0 + x1*x1*z0 -
     &     z1*x2*x2 - 2*z0*x2*x1 - x1*x1*z2)

   
      !diff(Atot,x1)*(Atot*4)
       Ftmp(2,1)=
     &   cf2*(-z1*x2*z2 + z1*x2*z0 + z1*x0*z2 - z1*x0*z0 +
     &     z0*x2*z2 + z0*x0*z2 - y2*y1*x2 + y2*y1*x0 +
     &     y2*y0*x2 + y0*x0*y2 + y0*y1*x2 - y0*y1*x0 +
     &     x1*z2*z2 + x1*z0*z0 + x1*y2*y2 + x1*y0*y0 -
     &     z0*z0*x2 - 2*x1*z0*z2 - z2*z2*x0
     &     - 2*x1*y0*y2 - y2*y2*x0 - y0*y0*x2)
    
        !diff(Atot,y1)*(Atot*4)
       Ftmp(2,2)=
     &   cf2*(-z2*z1*y2 + z2*z1*y0 + z2*z0*y2 + z0*y0*z2 +
     &     z0*z1*y2 - z0*z1*y0 - x1*y2*x2 + x1*y2*x0 +
     &     x1*y0*x2 - x1*y0*x0 + x0*y2*x2 + x0*y0*x2 -
     &     2*y1*z2*z0 - z2*z2*y0 - z0*z0*y2 + y1*z2*z2 +
     &     y1*z0*z0 + y1*x2*x2 + y1*x0*x0 - x0*x0*y2 -
     &     2*y1*x2*x0 - x2*x2*y0)


    
      !diff(Atot,z1)*(Atot*4)
       Ftmp(2,3)=
     &   cf2*(-y1*z2*y2 + y1*y0*z2 + y1*z0*y2 - y1*z0*y0 +
     &     y0*z2*y2 + y0*z0*y2 - x2*x1*z2 + z0*x2*x1 +
     &     x2*x0*z2 + x0*z0*x2 + x1*x0*z2 - x0*x1*z0 -
     &     y0*y0*z2 - 2*z1*y0*y2 - z0*y2*y2 + z1*y2*y2 +
     &     z1*y0*y0 + z1*x2*x2 + z1*x0*x0 - 2*z1*x0*x2 -
     &     z0*x2*x2 - x0*x0*z2)

    
      !diff(Atot,x2)*(Atot*4)
       Ftmp(3,1)=
     &   cf2*(-z1*x1*z2 + z1*x1*z0 + z1*x0*z2 + z1*x0*z0 +
     &     x1*z0*z2 - z0*x0*z2 - x1*y2*y1 + x1*y0*y2 +
     &     x1*y0*y1 + y2*y1*x0 - y0*x0*y2 + y0*y1*x0 +
     &     z1*z1*x2 + z0*z0*x2 + y1*y1*x2 + y0*y0*x2 -
     &     z1*z1*x0 - 2*z1*x2*z0 - x1*z0*z0 -
     &     x1*y0*y0 - y1*y1*x0 - 2*y0*y1*x2)
      

      !diff(Atot,y2)*(Atot*4)
       Ftmp(3,2)=
     &   cf2*(-y1*z2*z1 + y1*z2*z0 + y1*z0*z1 + z2*z1*y0 -
     &     z0*y0*z2 + z0*z1*y0 - x1*y1*x2 + x1*y1*x0 +
     &     x1*y0*x2 + x1*y0*x0 + x0*y1*x2 - x0*y0*x2 -
     &     y1*z0*z0 - z1*z1*y0 - 2*z0*z1*y2 + z1*z1*y2 +
     &     z0*z0*y2 + x1*x1*y2 + x0*x0*y2 - x0*x0*y1 -
     &     x1*x1*y0 - 2*x1*y2*x0)


      !diff(Atot,z2)*(Atot*4)
       Ftmp(3,3)=
     &   cf2*(-y1*z1*y2 + y1*z1*y0 + y1*z0*y2 + y1*z0*y0 +
     &     z1*y0*y2 - y0*z0*y2 - z1*x2*x1 + z1*x0*x2 +
     &     z1*x0*x1 + z0*x2*x1 - x0*z0*x2 + x0*x1*z0 -
     &     y1*y1*z0 - 2*y1*y0*z2 - z1*y0*y0 + y1*y1*z2 +
     &     y0*y0*z2 + x1*x1*z2 + x0*x0*z2 - z1*x0*x0 -
     &     x1*x1*z0 - 2*x1*x0*z2)

       if(Ktmpflag) then
       Ktmp(1:3,1:3,1:3,1:3)=0.0
      t2 = x1-x2
      t3 = y2**2
      t4 = z2**2
      t5 = y0*y1
      t6 = z0*z1
      t7 = y1**2
      t8 = z1**2
      t9 = x0*y1*2.0D0
      t10 = x2*y0
      t11 = x0*z1*2.0D0
      t12 = x2*z0
      t13 = y1-y2
      t14 = z1-z2
      t15 = x2*y1
      t16 = x2*y2
      t17 = z1*z2
      t18 = x2**2
      t19 = x1*y0*2.0D0
      t20 = x0*y2
      t21 = x1*y1
      t22 = x1*y2
      t23 = x0*x1
      t24 = x1**2
      t25 = y0*z1*2.0D0
      t26 = y2*z0
      t27 = x2*z1
      t28 = x2*z2
      t29 = y2*z1
      t30 = y2*z2
      t31 = x1*x2
      t32 = y1*y2
      t33 = x1*z0*2.0D0
      t34 = x0*z2
      t35 = x1*z1
      t36 = x1*z2
      t37 = y1*z0*2.0D0
      t38 = y0*z2
      t39 = y1*z1
      t40 = y1*z2
      t41 = x0*y1
      t42 = x2*y0*2.0D0
      t43 = x0*z1
      t44 = x2*z0*2.0D0
      t45 = x0-x2
      t46 = y0**2
      t47 = z0**2
      t48 = x0*y2*2.0D0
      t49 = x0*z2*2.0D0
      t50 = x1*y0
      t51 = z0*z2
      t52 = y0*z1
      t53 = y2*z0*2.0D0
      t54 = y0-y2
      t55 = z0-z2
      t56 = x0*y0
      t57 = x0**2
      t58 = y0*z2*2.0D0
      t59 = x1*z0
      t60 = y1*z0
      t61 = x0*x2
      t62 = y0*y2
      t63 = x0*z0
      t64 = y0*z0
      t65 = x2*y1*2.0D0
      t66 = x2*z1*2.0D0
      t67 = t5+t6-t17-t32-t46-t47+t51+t62
      t68 = cf2*t67
      t69 = -t9+t20-t22-t42+t50+t56+t65
      t70 = cf2*t69
      t71 = -t11+t34-t36-t44+t59+t63+t66
      t72 = cf2*t71
      t73 = x0-x1
      t74 = x1*y2*2.0D0
      t75 = y2*z1*2.0D0
      t76 = t10-t15-t19+t41-t48+t56+t74
      t77 = cf2*t76
      t78 = t6-t17+t23-t31-t47+t51-t57+t61
      t79 = cf2*t78
      t80 = -t25+t38-t40-t53+t60+t64+t75
      t81 = cf2*t80
      t82 = y0-y1
      t83 = z0-z1
      t84 = x1*z2*2.0D0
      t85 = y1*z2*2.0D0
      t86 = t12-t27-t33+t43-t49+t63+t84
      t87 = cf2*t86
      t88 = t26-t29-t37+t52-t58+t64+t85
      t89 = cf2*t88
      t90 = t5+t23-t31-t32-t46-t57+t61+t62
      t91 = cf2*t90
      Ktmp(1,1,1,1) = cf2*(t3+t4+t7+t8-y1*y2*2.0D0-z1*z2*2.0D0)
      Ktmp(1,1,1,2) = -cf2*t2*t13
      Ktmp(1,1,1,3) = -cf2*t2*t14
      Ktmp(1,1,2,1) = -cf2*(t3+t4+t5+t6-y0*y2-y1*y2-z0*z2-z1*z2)
      Ktmp(1,1,2,2) = cf2
     &   *(t9+t10+t16+t22-x1*y0-x0*y2*2.0D0-x2*y1*2.0D0)
      Ktmp(1,1,2,3) = cf2
     &   *(t11+t12+t28+t36-x1*z0-x0*z2*2.0D0-x2*z1*2.0D0)
      Ktmp(1,1,3,1) = cf2*(t5+t6-t7-t8+t17+t32-y0*y2-z0*z2)
      Ktmp(1,1,3,2) = cf2*(-t9-t10+t15+t21+t48+t50-x1*y2*2.0D0)
      Ktmp(1,1,3,3) = cf2*(-t11-t12+t27+t35+t49+t59-x1*z2*2.0D0)
      Ktmp(1,2,1,1) = -cf2*t2*t13
      Ktmp(1,2,1,2) = cf2*(t4+t8+t18+t24-x1*x2*2.0D0-z1*z2*2.0D0)
      Ktmp(1,2,1,3) = -cf2*t13*t14
      Ktmp(1,2,2,1) = cf2
     &   *(t15+t16+t19+t20-x0*y1-x2*y0*2.0D0-x1*y2*2.0D0)
      Ktmp(1,2,2,2) = -cf2*(t4+t6-t17+t18+t23-x0*x2-x1*x2-z0*z2)
      Ktmp(1,2,2,3) = cf2
     &   *(t25+t26+t30+t40-y1*z0-y0*z2*2.0D0-y2*z1*2.0D0)
      Ktmp(1,2,3,1) = cf2*(-t19-t20+t21+t22+t41+t42-x2*y1*2.0D0)
      Ktmp(1,2,3,2) = cf2*(t6-t8+t17+t23-t24+t31-x0*x2-z0*z2)
      Ktmp(1,2,3,3) = cf2*(-t25-t26+t29+t39+t58+t60-y1*z2*2.0D0)
      Ktmp(1,3,1,1) = -cf2*t2*t14
      Ktmp(1,3,1,2) = -cf2*t13*t14
      Ktmp(1,3,1,3) = cf2*(t3+t7+t18+t24-x1*x2*2.0D0-y1*y2*2.0D0)
      Ktmp(1,3,2,1) = cf2
     &   *(t27+t28+t33+t34-x0*z1-x2*z0*2.0D0-x1*z2*2.0D0)
      Ktmp(1,3,2,2) = cf2
     &   *(t29+t30+t37+t38-y0*z1-y2*z0*2.0D0-y1*z2*2.0D0)
      Ktmp(1,3,2,3) = -cf2*(t3+t5+t18+t23-t31-t32-x0*x2-y0*y2)
      Ktmp(1,3,3,1) = cf2*(-t33-t34+t35+t36+t43+t44-x2*z1*2.0D0)
      Ktmp(1,3,3,2) = cf2*(-t37-t38+t39+t40+t52+t53-y2*z1*2.0D0)
      Ktmp(1,3,3,3) = cf2*(t5-t7+t23-t24+t31+t32-x0*x2-y0*y2)
      Ktmp(2,1,1,1) = -cf2*(t3+t4+t5+t6-t17-t32-y0*y2-z0*z2)
      Ktmp(2,1,1,2) = cf2*(t15+t16+t19+t20-t41-t42-x1*y2*2.0D0)
      Ktmp(2,1,1,3) = cf2*(t27+t28+t33+t34-t43-t44-x1*z2*2.0D0)
      Ktmp(2,1,2,1) = cf2*(t3+t4+t46+t47-y0*y2*2.0D0-z0*z2*2.0D0)
      Ktmp(2,1,2,2) = -cf2*t45*t54
      Ktmp(2,1,2,3) = -cf2*t45*t55
      Ktmp(2,1,3,1) = t68
      Ktmp(2,1,3,2) = t77
      Ktmp(2,1,3,3) = t87
      Ktmp(2,2,1,1) = cf2*(t9+t10+t16+t22-t48-t50-x2*y1*2.0D0)
      Ktmp(2,2,1,2) = -cf2*(t4+t6-t17+t18+t23-t31-t51-x0*x2)
      Ktmp(2,2,1,3) = cf2*(t29+t30+t37+t38-t52-t53-y1*z2*2.0D0)
      Ktmp(2,2,2,1) = -cf2*t45*t54
      Ktmp(2,2,2,2) = cf2*(t4+t18+t47+t57-x0*x2*2.0D0-z0*z2*2.0D0)
      Ktmp(2,2,2,3) = -cf2*t54*t55
      Ktmp(2,2,3,1) = t70
      Ktmp(2,2,3,2) = t79
      Ktmp(2,2,3,3) = t89
      Ktmp(2,3,1,1) = cf2*(t11+t12+t28+t36-t49-t59-x2*z1*2.0D0)
      Ktmp(2,3,1,2) = cf2*(t25+t26+t30+t40-t58-t60-y2*z1*2.0D0)
      Ktmp(2,3,1,3) = -cf2*(t3+t5+t18+t23-t31-t32-t61-t62)
      Ktmp(2,3,2,1) = -cf2*t45*t55
      Ktmp(2,3,2,2) = -cf2*t54*t55
      Ktmp(2,3,2,3) = cf2*(t3+t18+t46+t57-x0*x2*2.0D0-y0*y2*2.0D0)
      Ktmp(2,3,3,1) = t72
      Ktmp(2,3,3,2) = t81
      Ktmp(2,3,3,3) = t91
      Ktmp(3,1,1,1) = cf2*(t5+t6-t7-t8+t17+t32-t51-t62)
      Ktmp(3,1,1,2) = cf2*(-t19-t20+t21+t22+t41+t42-t65)
      Ktmp(3,1,1,3) = cf2*(-t33-t34+t35+t36+t43+t44-t66)
      Ktmp(3,1,2,1) = t68
      Ktmp(3,1,2,2) = t70
      Ktmp(3,1,2,3) = t72
      Ktmp(3,1,3,1) = cf2*(t7+t8+t46+t47-y0*y1*2.0D0-z0*z1*2.0D0)
      Ktmp(3,1,3,2) = -cf2*t73*t82
      Ktmp(3,1,3,3) = -cf2*t73*t83
      Ktmp(3,2,1,1) = cf2*(-t9-t10+t15+t21+t48+t50-t74)
      Ktmp(3,2,1,2) = cf2*(t6-t8+t17+t23-t24+t31-t51-t61)
      Ktmp(3,2,1,3) = cf2*(-t37-t38+t39+t40+t52+t53-t75)
      Ktmp(3,2,2,1) = t77
      Ktmp(3,2,2,2) = t79
      Ktmp(3,2,2,3) = t81
      Ktmp(3,2,3,1) = -cf2*t73*t82
      Ktmp(3,2,3,2) = cf2*(t8+t24+t47+t57-x0*x1*2.0D0-z0*z1*2.0D0)
      Ktmp(3,2,3,3) = -cf2*t82*t83
      Ktmp(3,3,1,1) = cf2*(-t11-t12+t27+t35+t49+t59-t84)
      Ktmp(3,3,1,2) = cf2*(-t25-t26+t29+t39+t58+t60-t85)
      Ktmp(3,3,1,3) = cf2*(t5-t7+t23-t24+t31+t32-t61-t62)
      Ktmp(3,3,2,1) = t87
      Ktmp(3,3,2,2) = t89
      Ktmp(3,3,2,3) = t91
      Ktmp(3,3,3,1) = -cf2*t73*t83
      Ktmp(3,3,3,2) = -cf2*t82*t83
      Ktmp(3,3,3,3) = cf2*(t7+t24+t46+t57-x0*x1*2.0D0-y0*y1*2.0D0)

       endif
      return
      end subroutine
!======================================================================
      subroutine VolumeForce(Ftmp
     &                            ,x0,y0,z0
     &                            ,x1,y1,z1
     &                            ,x2,y2,z2,cf
     &                           ,Ktmp,KtmpFlag)

      real x0,y0,z0,x1,y1,z1,x2,y2,z2,cf,t0
      real Ftmp(3,3),ktmp(3,3,3,3)
      real TmpVar_fbc1,TmpVar_fbc2,TmpVar_fbc3
      real TmpVar_fbc4,TmpVar_fbc5
      real TmpVar_fbc6,TmpVar_fbc7,TmpVar_fbc8
      real TmpVar_fbc9,TmpVar_fbc10
      real TmpVar_fbc11,TmpVar_fbc12,TmpVar_fbc13
      real TmpVar_fbc14,TmpVar_fbc15
      real TmpVar_fbc16,TmpVar_fbc17,TmpVar_fbc18
      real TmpVar_fbc19,TmpVar_fbc20
      logical KtmpFlag

      Ktmp(1:3,1:3,1:3,1:3)=0.0
      Ftmp(1:3,1:3)=0.0
      Ftmp(1,1) = cf*(y1*z2-y2*z1)*(1.0D0/6.0D0)
      Ftmp(1,2) = cf*(x1*z2-x2*z1)*(-1.0D0/6.0D0)
      Ftmp(1,3) = cf*(x1*y2-x2*y1)*(1.0D0/6.0D0)
      Ftmp(2,1) = cf*(y0*z2-y2*z0)*(-1.0D0/6.0D0)
      Ftmp(2,2) = cf*(x0*z2-x2*z0)*(1.0D0/6.0D0)
      Ftmp(2,3) = cf*(x0*y2-x2*y0)*(-1.0D0/6.0D0)
      Ftmp(3,1) = cf*(y0*z1-y1*z0)*(1.0D0/6.0D0)
      Ftmp(3,2) = cf*(x0*z1-x1*z0)*(-1.0D0/6.0D0)
      Ftmp(3,3) = cf*(x0*y1-x1*y0)*(1.0D0/6.0D0)
      if(Ktmpflag) then
      t2 = cf*z2*(1.0D0/6.0D0)
      t3 = cf*x2*(1.0D0/6.0D0)
      t4 = cf*y1*(1.0D0/6.0D0)
      t5 = cf*y2*(1.0D0/6.0D0)
      t6 = cf*z0*(1.0D0/6.0D0)
      t7 = cf*x0*(1.0D0/6.0D0)
      t8 = cf*z1*(1.0D0/6.0D0)
      t9 = cf*y0*(1.0D0/6.0D0)
      t10 = cf*x1*(1.0D0/6.0D0)
      Ktmp(1,1,2,2) = t2
      Ktmp(1,1,2,3) = cf*y2*(-1.0D0/6.0D0)
      Ktmp(1,1,3,2) = cf*z1*(-1.0D0/6.0D0)
      Ktmp(1,1,3,3) = t4
      Ktmp(1,2,2,1) = -t2
      Ktmp(1,2,2,3) = t3
      Ktmp(1,2,3,1) = t8
      Ktmp(1,2,3,3) = cf*x1*(-1.0D0/6.0D0)
      Ktmp(1,3,2,1) = t5
      Ktmp(1,3,2,2) = -t3
      Ktmp(1,3,3,1) = -t4
      Ktmp(1,3,3,2) = t10
      Ktmp(2,1,1,2) = -t2
      Ktmp(2,1,1,3) = t5
      Ktmp(2,1,3,2) = t6
      Ktmp(2,1,3,3) = cf*y0*(-1.0D0/6.0D0)
      Ktmp(2,2,1,1) = t2
      Ktmp(2,2,1,3) = -t3
      Ktmp(2,2,3,1) = -t6
      Ktmp(2,2,3,3) = t7
      Ktmp(2,3,1,1) = -t5
      Ktmp(2,3,1,2) = t3
      Ktmp(2,3,3,1) = t9
      Ktmp(2,3,3,2) = -t7
      Ktmp(3,1,1,2) = t8
      Ktmp(3,1,1,3) = -t4
      Ktmp(3,1,2,2) = -t6
      Ktmp(3,1,2,3) = t9
      Ktmp(3,2,1,1) = -t8
      Ktmp(3,2,1,3) = t10
      Ktmp(3,2,2,1) = t6
      Ktmp(3,2,2,3) = -t7
      Ktmp(3,3,1,1) = t4
      Ktmp(3,3,1,2) = -t10
      Ktmp(3,3,2,1) = -t9
      Ktmp(3,3,2,2) = t7
      endif

      return
      end subroutine
!======================================================================
      subroutine Bending_ForceNEW(Ftmp
     &                            ,costet0,sintet0,sinsign
     &                            ,x0,y0,z0
     &                            ,x1,y1,z1
     &                            ,x2,y2,z2
     &                            ,x3,y3,z3
     &                            ,cf
     &                           ,Ktmp,KtmpFlag)

      real sinsign,costet0,sintet0
      real x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3
      real cf
      real Ftmp(4,3),Ktmp(4,3,4,3)
      real betab,A0(4,3),B0(4,3,4,3)
      logical KtmpFlag

      t2 = x0-x1
      t3 = y0-y2
      t4 = t2*t3
      t5 = x0-x2
      t6 = y0-y1
      t26 = t5*t6
      t7 = t4-t26
      t9 = z0-z2
      t10 = z0-z1
      t22 = t2*t9
      t23 = t5*t10
      t8 = t22-t23
      t29 = t6*t9
      t30 = t3*t10
      t11 = -t29+t30
      t12 = x1-x3
      t13 = y2-y3
      t14 = t12*t13
      t15 = x2-x3
      t16 = y1-y3
      t32 = t15*t16
      t17 = t14-t32
      t19 = z2-z3
      t20 = z1-z3
      t24 = t12*t19
      t25 = t15*t20
      t18 = t24-t25
      t35 = t16*t19
      t36 = t13*t20
      t21 = -t35+t36
      t27 = t7**2
      t28 = t8**2
      t31 = t29-t30
      t33 = t17**2
      t34 = t18**2
      t37 = t35-t36
      t38 = t31**2
      t39 = t27+t28+t38
      t40 = 1.0D0/sqrt(t39)
      t41 = t37**2
      t42 = t33+t34+t41
      t43 = 1.0D0/sqrt(t42)
      t45 = t8*t18*t40*t43
      t46 = t31*t37*t40*t43
      t44 = cos(signsin*acos(-t45-t46-t7*t17*1.0D0/sqrt(t27+t28+t11**2)*
     &1.0D0/sqrt(t33+t34+t21**2)))
      t47 = -t45-t46-t7*t17*t40*t43
      t48 = acos(t47)
      t49 = signsin*t48
      Betab 
     &= cf/sqrt(-t44**2+1.0D0)*(costet0*sin(t49)-sintet0*cos(t49))

      t2 = x0-x1
      t3 = y0-y2
      t4 = t2*t3
      t5 = x0-x2
      t6 = y0-y1
      t26 = t5*t6
      t7 = t4-t26
      t9 = z0-z2
      t10 = z0-z1
      t22 = t2*t9
      t23 = t5*t10
      t8 = t22-t23
      t29 = t6*t9
      t30 = t3*t10
      t11 = -t29+t30
      t12 = x1-x3
      t13 = y2-y3
      t14 = t12*t13
      t15 = x2-x3
      t16 = y1-y3
      t32 = t15*t16
      t17 = t14-t32
      t19 = z2-z3
      t20 = z1-z3
      t24 = t12*t19
      t25 = t15*t20
      t18 = t24-t25
      t35 = t16*t19
      t36 = t13*t20
      t21 = -t35+t36
      t27 = t7**2
      t28 = t8**2
      t31 = t29-t30
      t33 = t17**2
      t34 = t18**2
      t37 = t35-t36
      t38 = t31**2
      t39 = t27+t28+t38
      t40 = 1.0D0/sqrt(t39)
      t41 = t37**2
      t42 = t33+t34+t41
      t43 = 1.0D0/sqrt(t42)
      t44 = y1-y2
      t45 = z1-z2
      t46 = 1.0D0/t39
      t47 = x1-x2
      t52 = t7*t17*t40*t43
      t53 = t8*t18*t40*t43
      t54 = t31*t37*t40*t43
      t48 = -t52-t53-t54
      t49 = acos(t48)
      t50 = signsin*t49
      t51 = cos(t50)
      t55 = 1.0D0/t42
      t56 = t17*t44
      t57 = t18*t45
      t58 = t56+t57
      t59 = t7*t44
      t60 = t8*t45
      t61 = t59+t60
      t62 = t17*t47
      t63 = t62-t37*t45
      t64 = t7*t47
      t65 = t64-t31*t45
      t66 = t18*t47
      t67 = t37*t44
      t68 = t66+t67
      t69 = t8*t47
      t70 = t31*t44
      t71 = t69+t70
      Ftmp(1,1) 
     &= t40*t43*t58+t46*t61*cos(signsin*acos(-t8*t18*t40*t43-t31
     &*t37*t40*t43-t7*t17*1.0D0/sqrt(t27+t28+t11**2)*1.0D0/sqrt(t33+t34+
     &t21**2)))
      Ftmp(1,2) 
     &= -t40*t43*t63-t46*t51*t65
      Ftmp(1,3) 
     &= -t40*t43*t68-t46*t51*t71
      Ftmp(2,1) 
     &= -t46*t51*(t3*t7+t8*t9)+t51*t55*(t13*t17+t18*t19)-t40*t43
     &*(t3*t17-t7*t13-t8*t19+t9*t18)
      Ftmp(2,2) 
     &= t46*t51*(t5*t7-t9*t31)-t51*t55*(t15*t17-t19*t37)+t40*t43
     &*(t5*t17-t7*t15-t9*t37+t19*t31)
      Ftmp(2,3) 
     &= t46*t51*(t5*t8+t3*t31)-t51*t55*(t15*t18+t13*t37)+t40*t43
     &*(t5*t18-t8*t15+t3*t37-t13*t31)
      Ftmp(3,1) 
     &= t46*t51*(t6*t7+t8*t10)-t51*t55*(t16*t17+t18*t20)+t40*t43
     &*(t6*t17-t7*t16-t8*t20+t10*t18)
      Ftmp(3,2) 
     &= -t46*t51*(t2*t7-t10*t31)+t51*t55*(t12*t17-t20*t37)-t40*t
     &43*(t2*t17-t7*t12-t10*t37+t20*t31)
      Ftmp(3,3) 
     &= -t46*t51*(t2*t8+t6*t31)+t51*t55*(t12*t18+t16*t37)-t40*t4
     &3*(t2*t18-t8*t12+t6*t37-t16*t31)
      Ftmp(4,1) 
     &= t40*t43*t61+t51*t55*t58
      Ftmp(4,2) 
     &= -t40*t43*t65-t51*t55*t63
      Ftmp(4,3) 
     &= -t40*t43*t71-t51*t55*t68

      if(Ktmpflag) then

      t2 = x0-x1
      t3 = y0-y2
      t4 = t2*t3
      t5 = x0-x2
      t6 = y0-y1
      t26 = t5*t6
      t7 = t4-t26
      t9 = z0-z2
      t10 = z0-z1
      t22 = t2*t9
      t23 = t5*t10
      t8 = t22-t23
      t29 = t6*t9
      t30 = t3*t10
      t11 = -t29+t30
      t12 = x1-x3
      t13 = y2-y3
      t14 = t12*t13
      t15 = x2-x3
      t16 = y1-y3
      t32 = t15*t16
      t17 = t14-t32
      t19 = z2-z3
      t20 = z1-z3
      t24 = t12*t19
      t25 = t15*t20
      t18 = t24-t25
      t35 = t16*t19
      t36 = t13*t20
      t21 = -t35+t36
      t27 = t7**2
      t28 = t8**2
      t31 = t29-t30
      t33 = t17**2
      t34 = t18**2
      t37 = t35-t36
      t38 = t31**2
      t39 = t27+t28+t38
      t40 = 1.0D0/sqrt(t39)
      t41 = t37**2
      t42 = t33+t34+t41
      t43 = 1.0D0/sqrt(t42)
      t45 = t8*t18*t40*t43
      t46 = t31*t37*t40*t43
      t44 = cos(signsin*acos(-t45-t46-t7*t17*1.0D0/sqrt(t27+t28+t11**2)*
     &1.0D0/sqrt(t33+t34+t21**2)))
      t54 = t7*t17*t40*t43
      t47 = t45+t46+t54
      t48 = y1-y2
      t49 = z1-z2
      t50 = t7*t48*2.0D0
      t51 = t8*t49*2.0D0
      t52 = t50+t51
      t53 = 1.0D0/t39**(3.0D0/2.0D0)
      t55 = t47**2
      t56 = -t55+1.0D0
      t57 = 1.0D0/sqrt(t56)
      t58 = t7*t17*t43*t52*t53*(1.0D0/2.0D0)
      t59 = t8*t18*t43*t52*t53*(1.0D0/2.0D0)
      t60 = t31*t37*t43*t52*t53*(1.0D0/2.0D0)
      t67 = t17*t40*t43*t48
      t68 = t18*t40*t43*t49
      t61 = t58+t59+t60-t67-t68
      t62 = -t45-t46-t54
      t63 = acos(t62)
      t64 = signsin*t63
      t65 = cos(t64)
      t66 = sin(t64)
      t69 = t65**2
      t70 = -t69+1.0D0
      t71 = x1-x2
      t72 = t7*t71*2.0D0
      t74 = t31*t49*2.0D0
      t73 = t72-t74
      t75 = t37*t40*t43*t49
      t76 = t7*t17*t43*t53*t73*(1.0D0/2.0D0)
      t77 = t8*t18*t43*t53*t73*(1.0D0/2.0D0)
      t78 = t31*t37*t43*t53*t73*(1.0D0/2.0D0)
      t83 = t17*t40*t43*t71
      t79 = t75+t76+t77+t78-t83
      t80 = 1.0D0/t70**(3.0D0/2.0D0)
      t81 = costet0*t66
      t92 = sintet0*t65
      t82 = t81-t92
      t84 = 1.0D0/sqrt(t70)
      t85 = t8*t71*2.0D0
      t86 = t31*t48*2.0D0
      t87 = t85+t86
      t88 = t7*t17*t43*t53*t87*(1.0D0/2.0D0)
      t89 = t8*t18*t43*t53*t87*(1.0D0/2.0D0)
      t90 = t31*t37*t43*t53*t87*(1.0D0/2.0D0)
      t93 = t18*t40*t43*t71
      t94 = t37*t40*t43*t48
      t91 = t88+t89+t90-t93-t94
      t95 = t3*t7*2.0D0
      t96 = t8*t9*2.0D0
      t97 = t95+t96
      t98 = t13*t17*2.0D0
      t99 = t18*t19*2.0D0
      t100 = t98+t99
      t101 = 1.0D0/t42**(3.0D0/2.0D0)
      t102 = t3*t17*t40*t43
      t103 = t9*t18*t40*t43
      t104 = t7*t17*t40*t100*t101*(1.0D0/2.0D0)
      t105 = t8*t18*t40*t100*t101*(1.0D0/2.0D0)
      t106 = t31*t37*t40*t100*t101*(1.0D0/2.0D0)
      t108 = t7*t13*t40*t43
      t109 = t8*t19*t40*t43
      t110 = t7*t17*t43*t53*t97*(1.0D0/2.0D0)
      t111 = t8*t18*t43*t53*t97*(1.0D0/2.0D0)
      t112 = t31*t37*t43*t53*t97*(1.0D0/2.0D0)
      t107 = t102+t103+t104+t105+t106-t108-t109-t110-t111-t112
      t113 = t5*t7*2.0D0
      t117 = t9*t31*2.0D0
      t114 = t113-t117
      t115 = t15*t17*2.0D0
      t118 = t19*t37*2.0D0
      t116 = t115-t118
      t119 = t5*t17*t40*t43
      t120 = t19*t31*t40*t43
      t121 = t7*t17*t40*t101*t116*(1.0D0/2.0D0)
      t122 = t8*t18*t40*t101*t116*(1.0D0/2.0D0)
      t123 = t31*t37*t40*t101*t116*(1.0D0/2.0D0)
      t125 = t7*t15*t40*t43
      t126 = t9*t37*t40*t43
      t127 = t7*t17*t43*t53*t114*(1.0D0/2.0D0)
      t128 = t8*t18*t43*t53*t114*(1.0D0/2.0D0)
      t129 = t31*t37*t43*t53*t114*(1.0D0/2.0D0)
      t124 = t119+t120+t121+t122+t123-t125-t126-t127-t128-t129
      t130 = t5*t8*2.0D0
      t131 = t3*t31*2.0D0
      t132 = t130+t131
      t133 = t15*t18*2.0D0
      t134 = t13*t37*2.0D0
      t135 = t133+t134
      t136 = t5*t18*t40*t43
      t137 = t3*t37*t40*t43
      t138 = t7*t17*t40*t101*t135*(1.0D0/2.0D0)
      t139 = t8*t18*t40*t101*t135*(1.0D0/2.0D0)
      t140 = t31*t37*t40*t101*t135*(1.0D0/2.0D0)
      t142 = t8*t15*t40*t43
      t143 = t13*t31*t40*t43
      t144 = t7*t17*t43*t53*t132*(1.0D0/2.0D0)
      t145 = t8*t18*t43*t53*t132*(1.0D0/2.0D0)
      t146 = t31*t37*t43*t53*t132*(1.0D0/2.0D0)
      t141 = t136+t137+t138+t139+t140-t142-t143-t144-t145-t146
      t147 = t6*t7*2.0D0
      t148 = t8*t10*2.0D0
      t149 = t147+t148
      t150 = t16*t17*2.0D0
      t151 = t18*t20*2.0D0
      t152 = t150+t151
      t153 = t6*t17*t40*t43
      t154 = t10*t18*t40*t43
      t155 = t7*t17*t40*t101*t152*(1.0D0/2.0D0)
      t156 = t8*t18*t40*t101*t152*(1.0D0/2.0D0)
      t157 = t31*t37*t40*t101*t152*(1.0D0/2.0D0)
      t159 = t7*t16*t40*t43
      t160 = t8*t20*t40*t43
      t161 = t7*t17*t43*t53*t149*(1.0D0/2.0D0)
      t162 = t8*t18*t43*t53*t149*(1.0D0/2.0D0)
      t163 = t31*t37*t43*t53*t149*(1.0D0/2.0D0)
      t158 = t153+t154+t155+t156+t157-t159-t160-t161-t162-t163
      t164 = t2*t7*2.0D0
      t168 = t10*t31*2.0D0
      t165 = t164-t168
      t166 = t12*t17*2.0D0
      t169 = t20*t37*2.0D0
      t167 = t166-t169
      t170 = t2*t17*t40*t43
      t171 = t20*t31*t40*t43
      t172 = t7*t17*t40*t101*t167*(1.0D0/2.0D0)
      t173 = t8*t18*t40*t101*t167*(1.0D0/2.0D0)
      t174 = t31*t37*t40*t101*t167*(1.0D0/2.0D0)
      t176 = t7*t12*t40*t43
      t177 = t10*t37*t40*t43
      t178 = t7*t17*t43*t53*t165*(1.0D0/2.0D0)
      t179 = t8*t18*t43*t53*t165*(1.0D0/2.0D0)
      t180 = t31*t37*t43*t53*t165*(1.0D0/2.0D0)
      t175 = t170+t171+t172+t173+t174-t176-t177-t178-t179-t180
      t181 = t2*t8*2.0D0
      t182 = t6*t31*2.0D0
      t183 = t181+t182
      t184 = t12*t18*2.0D0
      t185 = t16*t37*2.0D0
      t186 = t184+t185
      t187 = t2*t18*t40*t43
      t188 = t6*t37*t40*t43
      t189 = t7*t17*t40*t101*t186*(1.0D0/2.0D0)
      t190 = t8*t18*t40*t101*t186*(1.0D0/2.0D0)
      t191 = t31*t37*t40*t101*t186*(1.0D0/2.0D0)
      t193 = t8*t12*t40*t43

      t194 = t16*t31*t40*t43
      t195 = t7*t17*t43*t53*t183*(1.0D0/2.0D0)
      t196 = t8*t18*t43*t53*t183*(1.0D0/2.0D0)
      t197 = t31*t37*t43*t53*t183*(1.0D0/2.0D0)
      t192 = t187+t188+t189+t190+t191-t193-t194-t195-t196-t197
      t198 = t17*t48*2.0D0
      t199 = t18*t49*2.0D0
      t200 = t198+t199
      t201 = t7*t17*t40*t101*t200*(1.0D0/2.0D0)
      t202 = t8*t18*t40*t101*t200*(1.0D0/2.0D0)
      t203 = t31*t37*t40*t101*t200*(1.0D0/2.0D0)
      t205 = t7*t40*t43*t48
      t206 = t8*t40*t43*t49
      t204 = t201+t202+t203-t205-t206
      t207 = t37*t49*2.0D0
      t209 = t17*t71*2.0D0
      t208 = t207-t209
      t210 = t7*t40*t43*t71
      t211 = t7*t17*t40*t101*t208*(1.0D0/2.0D0)
      t212 = t8*t18*t40*t101*t208*(1.0D0/2.0D0)
      t213 = t31*t37*t40*t101*t208*(1.0D0/2.0D0)
      t215 = t31*t40*t43*t49
      t214 = t210+t211+t212+t213-t215
      t216 = t18*t71*2.0D0
      t217 = t37*t48*2.0D0
      t218 = t216+t217
      t219 = t7*t17*t40*t101*t218*(1.0D0/2.0D0)
      t220 = t8*t18*t40*t101*t218*(1.0D0/2.0D0)
      t221 = t31*t37*t40*t101*t218*(1.0D0/2.0D0)
      t223 = t8*t40*t43*t71
      t224 = t31*t40*t43*t48
      t222 = t219+t220+t221-t223-t224
      A0(1,1) = -1.0D0/sqrt(-t44**2+1.0D0)*(costet0*signsin*t57*t61*t65+
     &signsin*sintet0*t57*t61*t66)+signsin*t57*t61*t65*t66*t80*t82
      A0(1,2) = t84*(costet0*signsin*t57*t65*t79+signsin*sintet0*t57*t66
     &*t79)-signsin*t57*t65*t66*t79*t80*t82
      A0(1,3) = t84*(costet0*signsin*t57*t65*t91+signsin*sintet0*t57*t66
     &*t91)-signsin*t57*t65*t66*t80*t82*t91
      A0(2,1) = -t84*(costet0*signsin*t57*t65*t107+signsin*sintet0*t57*t
     &66*t107)+signsin*t57*t65*t66*t80*t82*(t102+t103+t104+t105+t106-t10
     &8-t109-t110-t111-t112)
      A0(2,2) = t84*(costet0*signsin*t57*t65*t124+signsin*sintet0*t57*t6
     &6*t124)-signsin*t57*t65*t66*t80*t82*t124
      A0(2,3) = t84*(costet0*signsin*t57*t65*t141+signsin*sintet0*t57*t6
     &6*t141)-signsin*t57*t65*t66*t80*t82*t141
      A0(3,1) = t84*(costet0*signsin*t57*t65*t158+signsin*sintet0*t57*t6
     &6*t158)-signsin*t57*t65*t66*t80*t82*t158
      A0(3,2) = -t84*(costet0*signsin*t57*t65*t175+signsin*sintet0*t57*t
     &66*t175)+signsin*t57*t65*t66*t80*t82*(t170+t171+t172+t173+t174-t17
     &6-t177-t178-t179-t180)
      A0(3,3) = -t84*(costet0*signsin*t57*t65*t192+signsin*sintet0*t57*t
     &66*t192)+signsin*t57*t65*t66*t80*t82*(t187+t188+t189+t190+t191-t19
     &3-t194-t195-t196-t197)
      A0(4,1) = -t84*(costet0*signsin*t57*t65*t204+signsin*sintet0*t57*t
     &66*t204)+signsin*t57*t65*t66*t80*t82*t204
      A0(4,2) = -t84*(costet0*signsin*t57*t65*t214+signsin*sintet0*t57*t
     &66*t214)+signsin*t57*t65*t66*t80*t82*(t210+t211+t212+t213-t215)
      A0(4,3) = t84*(costet0*signsin*t57*t65*t222+signsin*sintet0*t57*t6
     &6*t222)-signsin*t57*t65*t66*t80*t82*t222

      do j1=1,4
      do k1=1,3
         do j2=1,4
         do k2=1,4
            Ktmp(j1,k1,j2,k2)=Ftmp(j1,k1)*A0(j2,k2)*cf
         enddo
         enddo
      enddo
      enddo

      t2 = x0-x1
      t3 = y0-y2
      t4 = t2*t3
      t5 = x0-x2
      t6 = y0-y1
      t26 = t5*t6
      t7 = t4-t26
      t9 = z0-z2
      t10 = z0-z1
      t22 = t2*t9
      t23 = t5*t10
      t8 = t22-t23
      t29 = t6*t9
      t30 = t3*t10
      t11 = -t29+t30
      t12 = x1-x3
      t13 = y2-y3
      t14 = t12*t13
      t15 = x2-x3
      t16 = y1-y3
      t32 = t15*t16
      t17 = t14-t32
      t19 = z2-z3
      t20 = z1-z3
      t24 = t12*t19
      t25 = t15*t20
      t18 = t24-t25
      t35 = t16*t19
      t36 = t13*t20
      t21 = -t35+t36
      t27 = t7**2
      t28 = t8**2
      t31 = t29-t30
      t33 = t17**2
      t34 = t18**2
      t37 = t35-t36
      t38 = t31**2
      t39 = t27+t28+t38
      t40 = 1.0D0/sqrt(t39)
      t41 = t37**2
      t42 = t33+t34+t41
      t43 = 1.0D0/sqrt(t42)
      t44 = y1-y2
      t45 = z1-z2
      t46 = t7*t44*2.0D0
      t47 = t8*t45*2.0D0
      t48 = t46+t47
      t55 = t7*t17*t40*t43
      t56 = t8*t18*t40*t43
      t57 = t31*t37*t40*t43
      t49 = -t55-t56-t57
      t50 = acos(t49)
      t51 = signsin*t50
      t52 = t7*t44
      t53 = t8*t45
      t54 = t52+t53
      t58 = t55+t56+t57
      t59 = 1.0D0/t39
      t60 = 1.0D0/t39**(3.0D0/2.0D0)
      t61 = cos(t51)
      t62 = x1-x2
      t63 = 1.0D0/t39**2
      t64 = t17*t44
      t65 = t18*t45
      t66 = t64+t65
      t67 = t7*t62*2.0D0
      t73 = t31*t45*2.0D0
      t68 = t67-t73
      t69 = sin(t51)
      t70 = t58**2
      t71 = -t70+1.0D0
      t72 = 1.0D0/sqrt(t71)
      t74 = t8*t62*2.0D0
      t75 = t31*t44*2.0D0
      t76 = t74+t75
      t77 = t3*t7*2.0D0
      t78 = t8*t9*2.0D0
      t79 = t77+t78
      t80 = t13*t17*2.0D0
      t81 = t18*t19*2.0D0
      t82 = t80+t81
      t83 = 1.0D0/t42**(3.0D0/2.0D0)
      t84 = t5*t7*2.0D0
      t86 = t9*t31*2.0D0
      t85 = t84-t86
      t87 = t15*t17*2.0D0
      t89 = t19*t37*2.0D0
      t88 = t87-t89
      t90 = t5*t8*2.0D0
      t91 = t3*t31*2.0D0
      t92 = t90+t91
      t93 = t15*t18*2.0D0
      t94 = t13*t37*2.0D0
      t95 = t93+t94
      t96 = t6*t7*2.0D0
      t97 = t8*t10*2.0D0
      t98 = t96+t97
      t99 = t16*t17*2.0D0
      t100 = t18*t20*2.0D0
      t101 = t99+t100
      t102 = t2*t7*2.0D0
      t104 = t10*t31*2.0D0
      t103 = t102-t104
      t105 = t12*t17*2.0D0
      t107 = t20*t37*2.0D0
      t106 = t105-t107
      t108 = t2*t8*2.0D0
      t109 = t6*t31*2.0D0
      t110 = t108+t109
      t111 = t12*t18*2.0D0
      t112 = t16*t37*2.0D0
      t113 = t111+t112
      t114 = t44**2
      t115 = t45**2
      t116 = t114+t115
      t117 = t17*t44*2.0D0
      t118 = t18*t45*2.0D0
      t119 = t117+t118
      t120 = t17*t62*2.0D0
      t122 = t37*t45*2.0D0
      t121 = t120-t122
      t123 = t18*t62*2.0D0
      t124 = t37*t44*2.0D0
      t125 = t123+t124
      t126 = t7*t62
      t132 = t31*t45
      t127 = t126-t132
      t128 = t7*t17*t43*t48*t60*(1.0D0/2.0D0)
      t129 = t8*t18*t43*t48*t60*(1.0D0/2.0D0)
      t130 = t31*t37*t43*t48*t60*(1.0D0/2.0D0)
      t200 = t17*t40*t43*t44
      t201 = t18*t40*t43*t45
      t131 = t128+t129+t130-t200-t201
      t133 = t17*t62
      t140 = t37*t45
      t134 = t133-t140
      t135 = t37*t40*t43*t45
      t136 = t7*t17*t43*t60*t68*(1.0D0/2.0D0)
      t137 = t8*t18*t43*t60*t68*(1.0D0/2.0D0)
      t138 = t31*t37*t43*t60*t68*(1.0D0/2.0D0)
      t205 = t17*t40*t43*t62
      t139 = t135+t136+t137+t138-t205
      t141 = t7*t17*t43*t60*t76*(1.0D0/2.0D0)
      t142 = t8*t18*t43*t60*t76*(1.0D0/2.0D0)
      t143 = t31*t37*t43*t60*t76*(1.0D0/2.0D0)
      t206 = t18*t40*t43*t62
      t207 = t37*t40*t43*t44
      t144 = t141+t142+t143-t206-t207
      t145 = t3*t17*t40*t43
      t146 = t9*t18*t40*t43
      t147 = t7*t17*t40*t82*t83*(1.0D0/2.0D0)
      t148 = t8*t18*t40*t82*t83*(1.0D0/2.0D0)
      t149 = t31*t37*t40*t82*t83*(1.0D0/2.0D0)
      t208 = t7*t13*t40*t43
      t209 = t8*t19*t40*t43
      t210 = t7*t17*t43*t60*t79*(1.0D0/2.0D0)
      t211 = t8*t18*t43*t60*t79*(1.0D0/2.0D0)
      t212 = t31*t37*t43*t60*t79*(1.0D0/2.0D0)
      t150 = t145+t146+t147+t148+t149-t208-t209-t210-t211-t212
      t151 = t9*t45
      t152 = t19*t45
      t153 = t5*t17*t40*t43
      t154 = t19*t31*t40*t43
      t155 = t7*t17*t40*t83*t88*(1.0D0/2.0D0)
      t156 = t8*t18*t40*t83*t88*(1.0D0/2.0D0)
      t157 = t31*t37*t40*t83*t88*(1.0D0/2.0D0)
      t158 = t5*t18*t40*t43
      t159 = t3*t37*t40*t43
      t160 = t7*t17*t40*t83*t95*(1.0D0/2.0D0)
      t161 = t8*t18*t40*t83*t95*(1.0D0/2.0D0)
      t162 = t31*t37*t40*t83*t95*(1.0D0/2.0D0)
      t163 = t6*t17*t40*t43
      t164 = t10*t18*t40*t43
      t165 = t7*t17*t40*t83*t101*(1.0D0/2.0D0)
      t166 = t8*t18*t40*t83*t101*(1.0D0/2.0D0)
      t167 = t31*t37*t40*t83*t101*(1.0D0/2.0D0)
      t168 = t10*t45
      t169 = t20*t45
      t170 = t2*t17*t40*t43
      t171 = t20*t31*t40*t43
      t172 = t7*t17*t40*t83*t106*(1.0D0/2.0D0)
      t173 = t8*t18*t40*t83*t106*(1.0D0/2.0D0)
      t174 = t31*t37*t40*t83*t106*(1.0D0/2.0D0)

      t217 = t7*t12*t40*t43
      t218 = t10*t37*t40*t43
      t219 = t7*t17*t43*t60*t103*(1.0D0/2.0D0)
      t220 = t8*t18*t43*t60*t103*(1.0D0/2.0D0)
      t221 = t31*t37*t43*t60*t103*(1.0D0/2.0D0)
      t175 = t170+t171+t172+t173+t174-t217-t218-t219-t220-t221
      t176 = t2*t18*t40*t43
      t177 = t6*t37*t40*t43
      t178 = t7*t17*t40*t83*t113*(1.0D0/2.0D0)
      t179 = t8*t18*t40*t83*t113*(1.0D0/2.0D0)
      t180 = t31*t37*t40*t83*t113*(1.0D0/2.0D0)
      t226 = t8*t12*t40*t43
      t227 = t16*t31*t40*t43
      t228 = t7*t17*t43*t60*t110*(1.0D0/2.0D0)
      t229 = t8*t18*t43*t60*t110*(1.0D0/2.0D0)
      t230 = t31*t37*t43*t60*t110*(1.0D0/2.0D0)
      t181 = t176+t177+t178+t179+t180-t226-t227-t228-t229-t230
      t182 = t7*t17*t40*t83*t119*(1.0D0/2.0D0)
      t183 = t8*t18*t40*t83*t119*(1.0D0/2.0D0)
      t184 = t31*t37*t40*t83*t119*(1.0D0/2.0D0)
      t231 = t7*t40*t43*t44
      t232 = t8*t40*t43*t45
      t185 = t182+t183+t184-t231-t232
      t186 = t62**2
      t187 = t115+t186
      t188 = t31*t40*t43*t45
      t189 = t7*t17*t40*t83*t121*(1.0D0/2.0D0)
      t190 = t8*t18*t40*t83*t121*(1.0D0/2.0D0)
      t191 = t31*t37*t40*t83*t121*(1.0D0/2.0D0)
      t233 = t7*t40*t43*t62
      t192 = t188+t189+t190+t191-t233
      t193 = t7*t17*t40*t83*t125*(1.0D0/2.0D0)
      t194 = t8*t18*t40*t83*t125*(1.0D0/2.0D0)
      t195 = t31*t37*t40*t83*t125*(1.0D0/2.0D0)
      t235 = t8*t40*t43*t62
      t236 = t31*t40*t43*t44
      t196 = t193+t194+t195-t235-t236
      t197 = t8*t62
      t198 = t31*t44
      t199 = t197+t198
      t202 = t18*t62
      t203 = t37*t44
      t204 = t202+t203
      t213 = t5*t62
      t214 = t3*t44
      t215 = t15*t62
      t216 = t13*t44
      t222 = t2*t62
      t223 = t6*t44
      t224 = t12*t62
      t225 = t16*t44
      t234 = t114+t186
      t237 = t151+t214
      t238 = t152+t216
      t239 = t40*t43*t238
      t240 = t3*t7
      t241 = t8*t9
      t242 = t240+t241
      t243 = t13*t62
      t244 = t14-t32+t243
      t245 = t3*t62
      t246 = -t4+t26+t245
      t247 = t59*t61*t246
      t248 = t3*t17
      t249 = t9*t18
      t260 = t7*t13
      t261 = t8*t19
      t250 = t248+t249-t260-t261
      t251 = t13*t17
      t252 = t18*t19
      t253 = t251+t252
      t254 = 1.0D0/t42
      t255 = t19*t62
      t256 = t24-t25+t255
      t257 = t9*t62
      t258 = -t22+t23+t257
      t259 = t59*t61*t258
      t262 = 1.0D0/t42**2
      t297 = t7*t15*t40*t43
      t298 = t9*t37*t40*t43
      t299 = t7*t17*t43*t60*t85*(1.0D0/2.0D0)
      t300 = t8*t18*t43*t60*t85*(1.0D0/2.0D0)
      t301 = t31*t37*t43*t60*t85*(1.0D0/2.0D0)
      t263 = t153+t154+t155+t156+t157-t297-t298-t299-t300-t301
      t302 = t8*t15*t40*t43
      t303 = t13*t31*t40*t43
      t304 = t7*t17*t43*t60*t92*(1.0D0/2.0D0)
      t305 = t8*t18*t43*t60*t92*(1.0D0/2.0D0)
      t306 = t31*t37*t43*t60*t92*(1.0D0/2.0D0)
      t264 = t158+t159+t160+t161+t162-t302-t303-t304-t305-t306
      t307 = t7*t16*t40*t43
      t308 = t8*t20*t40*t43
      t309 = t7*t17*t43*t60*t98*(1.0D0/2.0D0)
      t310 = t8*t18*t43*t60*t98*(1.0D0/2.0D0)
      t311 = t31*t37*t43*t60*t98*(1.0D0/2.0D0)
      t265 = t163+t164+t165+t166+t167-t307-t308-t309-t310-t311
      t266 = t15*t44
      t267 = -t14+t32+t266
      t268 = t5*t44
      t269 = t4-t26+t268
      t270 = t59*t61*t269
      t271 = t5*t7
      t278 = t9*t31
      t272 = t271-t278
      t273 = t151+t213
      t274 = t152+t215
      t275 = t40*t43*t274
      t276 = t5*t17
      t277 = t19*t31
      t279 = t15*t17
      t286 = t19*t37
      t280 = t279-t286
      t281 = t19*t44
      t282 = t35-t36+t281
      t283 = t9*t44
      t284 = -t29+t30+t283
      t285 = t59*t61*t284
      t287 = t5*t13
      t288 = t3*t15
      t289 = t287+t288
      t290 = t40*t43*t289
      t295 = t7*t15
      t296 = t9*t37
      t291 = t276+t277-t295-t296
      t292 = t9*t19*2.0D0
      t293 = t9**2
      t294 = t19**2
      t312 = t9*t10
      t313 = t19*t20
      t314 = t10*t19
      t315 = t9*t20
      t316 = t15*t45
      t317 = -t24+t25+t316
      t318 = t5*t45
      t319 = t22-t23+t318
      t320 = t59*t61*t319
      t321 = t5*t8
      t322 = t3*t31
      t323 = t321+t322
      t324 = t13*t45
      t325 = -t35+t36+t324
      t326 = t3*t45
      t327 = t29-t30+t326
      t328 = t59*t61*t327
      t329 = t5*t18
      t330 = t3*t37
      t331 = t15*t18
      t332 = t13*t37
      t333 = t331+t332
      t334 = t213+t214
      t335 = t215+t216
      t336 = t40*t43*t335
      t337 = t5*t19
      t338 = t9*t15
      t339 = t337+t338
      t340 = t40*t43*t339
      t346 = t8*t15
      t347 = t13*t31
      t341 = t329+t330-t346-t347
      t342 = t3*t19
      t343 = t9*t13
      t344 = t342+t343
      t345 = t40*t43*t344
      t348 = t5*t15*2.0D0
      t349 = t3*t13*2.0D0
      t350 = t5**2
      t351 = t3**2
      t352 = t15**2
      t353 = t13**2
      t354 = t2*t5
      t355 = t3*t6
      t356 = t12*t15
      t357 = t13*t16
      t358 = t2*t15
      t359 = t5*t12
      t360 = t6*t13
      t361 = t3*t16
      t362 = t168+t223
      t363 = t59*t61*t362
      t364 = t169+t225
      t365 = t6*t7
      t366 = t8*t10
      t367 = t365+t366
      t368 = t16*t62
      t369 = t14-t32+t368
      t370 = t40*t43*t369
      t371 = t6*t62
      t372 = -t4+t26+t371
      t373 = t6*t17
      t374 = t10*t18
      t375 = t16*t17
      t376 = t18*t20
      t377 = t375+t376
      t378 = t20*t62
      t379 = t24-t25+t378
      t380 = t40*t43*t379
      t381 = t10*t62
      t382 = -t22+t23+t381
      t383 = t312+t355
      t384 = t313+t357
      t385 = t314+t315+t360+t361
      t386 = t40*t43*t385
      t393 = t7*t16
      t394 = t8*t20
      t387 = t373+t374-t393-t394
      t388 = t4-t5*t6*2.0D0
      t389 = t14-t15*t16*2.0D0
      t390 = t5*t16
      t391 = t6*t15
      t392 = t4+t14-t26-t32+t390+t391
      t395 = t22-t5*t10*2.0D0
      t396 = t24-t15*t20*2.0D0
      t397 = t5*t20
      t398 = t10*t15
      t399 = t22-t23+t24-t25+t397+t398
      t400 = t12*t44
      t401 = -t14+t32+t400
      t402 = t40*t43*t401
      t403 = t2*t44
      t404 = t4-t26+t403
      t405 = t2*t7
      t413 = t10*t31
      t406 = t405-t413
      t407 = t168+t222
      t408 = t59*t61*t407
      t409 = t169+t224
      t410 = t2*t17
      t411 = t20*t31
      t421 = t7*t12
      t422 = t10*t37
      t412 = t410+t411-t421-t422
      t414 = t12*t17
      t423 = t20*t37
      t415 = t414-t423
      t416 = t20*t44
      t417 = t35-t36+t416
      t418 = t40*t43*t417
      t419 = t10*t44
      t420 = -t29+t30+t419
      t424 = t26-t2*t3*2.0D0
      t425 = t32-t12*t13*2.0D0
      t426 = t2*t13
      t427 = t3*t12
      t428 = -t4-t14+t26+t32+t426+t427
      t429 = t312+t354
      t430 = t313+t356
      t431 = t314+t315+t358+t359
      t432 = t40*t43*t431
      t433 = t29-t3*t10*2.0D0
      t434 = t35-t13*t20*2.0D0
      t435 = t3*t20
      t436 = t10*t13
      t437 = t29-t30+t35-t36+t435+t436
      t438 = t2*t16
      t439 = t6*t12
      t440 = t438+t439
      t441 = t40*t43*t440
      t442 = t10*t20*2.0D0
      t443 = t10**2
      t444 = t20**2
      t445 = t12*t45
      t446 = -t24+t25+t445
      t447 = t40*t43*t446
      t448 = t2*t45
      t449 = t22-t23+t448
      t450 = t2*t8
      t451 = t6*t31
      t452 = t450+t451
      t453 = t16*t45
      t454 = -t35+t36+t453
      t455 = t40*t43*t454
      t456 = t6*t45
      t457 = t29-t30+t456
      t458 = t2*t18
      t459 = t6*t37
      t467 = t8*t12
      t468 = t16*t31
      t460 = t458+t459-t467-t468
      t461 = t12*t18
      t462 = t16*t37
      t463 = t461+t462
      t464 = t222+t223
      t465 = t59*t61*t464
      t466 = t224+t225
      t469 = t23-t2*t9*2.0D0
      t470 = t25-t12*t19*2.0D0
      t471 = t2*t19
      t472 = t9*t12
      t473 = -t22+t23-t24+t25+t471+t472
      t474 = t30-t6*t9*2.0D0
      t475 = t36-t16*t19*2.0D0
      t476 = t6*t19
      t477 = t9*t16
      t478 = -t29+t30-t35+t36+t476+t477
      t479 = t354+t355
      t480 = t356+t357
      t481 = t358+t359+t360+t361
      t482 = t40*t43*t481
      t483 = t2*t20
      t484 = t10*t12
      t485 = t483+t484
      t486 = t40*t43*t485
      t487 = t6*t20
      t488 = t10*t16
      t489 = t487+t488
      t490 = t40*t43*t489
      t491 = t2*t12*2.0D0
      t492 = t6*t16*2.0D0
      t493 = t2**2
      t494 = t6**2
      t495 = t12**2
      t496 = t16**2
      t497 = t40*t43*t116
      t498 = t61*t238*t254
      t499 = t40*t43*t269
      t500 = t40*t43*t319
      t501 = t40*t43*t362
      t502 = t61*t254*t401
      t503 = t61*t254*t446
      t504 = t40*t43*t187
      t505 = t40*t43*t246
      t506 = t61*t254*t274
      t507 = t40*t43*t327
      t508 = t61*t254*t369
      t509 = t40*t43*t407
      t510 = t61*t254*t454
      t511 = t40*t43*t234
      t512 = t40*t43*t258
      t513 = t40*t43*t284
      t514 = t61*t254*t335
      t515 = t61*t254*t379
      t516 = t61*t254*t417
      t517 = t40*t43*t464
      B0(1,1,1,1) = t59*t116*cos(signsin*acos(-t8*t18*t40*t43-t31*t37*t4
     &0*t43-t7*t17*1.0D0/sqrt(t27+t28+t11**2)*1.0D0/sqrt(t33+t34+t21**2)
     &))-t43*t48*t60*t66*(1.0D0/2.0D0)-t48*t54*t61*t63+signsin*t54*t59*t
     &69*t72*t131
      B0(1,1,1,2) = -t44*t59*t61*t62+t43*t60*t66*t68*(1.0D0/2.0D0)+t54*t
     &61*t63*t68-signsin*t54*t59*t69*t72*t139
      B0(1,1,1,3) = -t45*t59*t61*t62+t43*t60*t66*t76*(1.0D0/2.0D0)+t54*t
     &61*t63*t76-signsin*t54*t59*t69*t72*t144
      B0(1,1,2,1) = t239-t59*t61*t237+t43*t60*t66*t79*(1.0D0/2.0D0)+t54*
     &t61*t63*t79-t40*t66*t82*t83*(1.0D0/2.0D0)+signsin*t54*t59*t69*t72*
     &t150
      B0(1,1,2,2) = t270-t40*t43*t267-t43*t60*t66*t85*(1.0D0/2.0D0)-t54*
     &t61*t63*t85+t40*t66*t83*t88*(1.0D0/2.0D0)-signsin*t54*t59*t69*t72*
     &t263
      B0(1,1,2,3) = t320-t40*t43*t317-t43*t60*t66*t92*(1.0D0/2.0D0)-t54*
     &t61*t63*t92+t40*t66*t83*t95*(1.0D0/2.0D0)-signsin*t54*t59*t69*t72*
     &t264
      B0(1,1,3,1) = t363-t40*t43*t364-t43*t60*t66*t98*(1.0D0/2.0D0)-t54*
     &t61*t63*t98+t40*t66*t83*t101*(1.0D0/2.0D0)-signsin*t54*t59*t69*t72
     &*t265
      B0(1,1,3,2) = t402-t59*t61*t404+t43*t60*t66*t103*(1.0D0/2.0D0)+t54
     &*t61*t63*t103-t40*t66*t83*t106*(1.0D0/2.0D0)+signsin*t54*t59*t69*t
     &72*t175
      B0(1,1,3,3) = t447-t59*t61*t449+t43*t60*t66*t110*(1.0D0/2.0D0)+t54
     &*t61*t63*t110-t40*t66*t83*t113*(1.0D0/2.0D0)+signsin*t54*t59*t69*t
     &72*t181
      B0(1,1,4,1) = t497-t40*t66*t83*t119*(1.0D0/2.0D0)+signsin*t54*t59*
     &t69*t72*t185
      B0(1,1,4,2) = -t40*t43*t44*t62+t40*t66*t83*t121*(1.0D0/2.0D0)-sign
     &sin*t54*t59*t69*t72*t192
      B0(1,1,4,3) = -t40*t43*t45*t62+t40*t66*t83*t125*(1.0D0/2.0D0)-sign
     &sin*t54*t59*t69*t72*t196
      B0(1,2,1,1) = -t44*t59*t61*t62+t43*t48*t60*t134*(1.0D0/2.0D0)+t48*
     &t61*t63*t127-signsin*t59*t69*t72*t127*t131
      B0(1,2,1,2) = t59*t61*t187-t43*t60*t68*t134*(1.0D0/2.0D0)-t61*t63*
     &t68*t127+signsin*t59*t69*t72*t127*t139
      B0(1,2,1,3) = -t44*t45*t59*t61-t43*t60*t76*t134*(1.0D0/2.0D0)-t61*
     &t63*t76*t127+signsin*t59*t69*t72*t127*t144
      B0(1,2,2,1) = t247-t40*t43*t244-t43*t60*t79*t134*(1.0D0/2.0D0)-t61
     &*t63*t79*t127+t40*t82*t83*t134*(1.0D0/2.0D0)-signsin*t59*t69*t72*t
     &127*t150
      B0(1,2,2,2) = t275-t59*t61*t273+t43*t60*t85*t134*(1.0D0/2.0D0)+t61
     &*t63*t85*t127-t40*t83*t88*t134*(1.0D0/2.0D0)+signsin*t59*t69*t72*t
     &127*t263
      B0(1,2,2,3) = t328-t40*t43*t325+t43*t60*t92*t134*(1.0D0/2.0D0)+t61
     &*t63*t92*t127-t40*t83*t95*t134*(1.0D0/2.0D0)+signsin*t59*t69*t72*t
     &127*t264
      B0(1,2,3,1) = t370-t59*t61*t372+t43*t60*t98*t134*(1.0D0/2.0D0)+t61
     &*t63*t98*t127-t40*t83*t101*t134*(1.0D0/2.0D0)+signsin*t59*t69*t72*
     &t127*t265
      B0(1,2,3,2) = t408-t40*t43*t409-t43*t60*t103*t134*(1.0D0/2.0D0)-t6
     &1*t63*t103*t127+t40*t83*t106*t134*(1.0D0/2.0D0)-signsin*t59*t69*t7
     &2*t127*t175
      B0(1,2,3,3) = t455-t59*t61*t457-t43*t60*t110*t134*(1.0D0/2.0D0)-t6
     &1*t63*t110*t127+t40*t83*t113*t134*(1.0D0/2.0D0)-signsin*t59*t69*t7
     &2*t127*t181
      B0(1,2,4,1) = -t40*t43*t44*t62+t40*t83*t119*t134*(1.0D0/2.0D0)-sig
     &nsin*t59*t69*t72*t127*t185
      B0(1,2,4,2) = t504-t40*t83*t121*t134*(1.0D0/2.0D0)+signsin*t59*t69
     &*t72*t127*t192
      B0(1,2,4,3) = -t40*t43*t44*t45-t40*t83*t125*t134*(1.0D0/2.0D0)+sig
     &nsin*t59*t69*t72*t127*t196
      B0(1,3,1,1) = -t45*t59*t61*t62+t43*t48*t60*t204*(1.0D0/2.0D0)+t48*
     &t61*t63*t199-signsin*t59*t69*t72*t131*t199
      B0(1,3,1,2) = -t44*t45*t59*t61-t43*t60*t68*t204*(1.0D0/2.0D0)-t61*
     &t63*t68*t199+signsin*t59*t69*t72*t139*t199
      B0(1,3,1,3) = t59*t61*t234-t43*t60*t76*t204*(1.0D0/2.0D0)-t61*t63*
     &t76*t199+signsin*t59*t69*t72*t144*t199
      B0(1,3,2,1) = t259-t40*t43*t256-t43*t60*t79*t204*(1.0D0/2.0D0)-t61
     &*t63*t79*t199+t40*t82*t83*t204*(1.0D0/2.0D0)-signsin*t59*t69*t72*t
     &150*t199
      B0(1,3,2,2) = t285-t40*t43*t282+t43*t60*t85*t204*(1.0D0/2.0D0)+t61
     &*t63*t85*t199-t40*t83*t88*t204*(1.0D0/2.0D0)+signsin*t59*t69*t72*t
     &199*t263
      B0(1,3,2,3) = t336-t59*t61*t334+t43*t60*t92*t204*(1.0D0/2.0D0)+t61
     &*t63*t92*t199-t40*t83*t95*t204*(1.0D0/2.0D0)+signsin*t59*t69*t72*t
     &199*t264
      B0(1,3,3,1) = t380-t59*t61*t382+t43*t60*t98*t204*(1.0D0/2.0D0)+t61
     &*t63*t98*t199-t40*t83*t101*t204*(1.0D0/2.0D0)+signsin*t59*t69*t72*
     &t199*t265
      B0(1,3,3,2) = t418-t59*t61*t420-t43*t60*t103*t204*(1.0D0/2.0D0)-t6
     &1*t63*t103*t199+t40*t83*t106*t204*(1.0D0/2.0D0)-signsin*t59*t69*t7
     &2*t175*t199
      B0(1,3,3,3) = t465-t40*t43*t466-t43*t60*t110*t204*(1.0D0/2.0D0)-t6
     &1*t63*t110*t199+t40*t83*t113*t204*(1.0D0/2.0D0)-signsin*t59*t69*t7
     &2*t181*t199
      B0(1,3,4,1) = -t40*t43*t45*t62+t40*t83*t119*t204*(1.0D0/2.0D0)-sig
     &nsin*t59*t69*t72*t185*t199
      B0(1,3,4,2) = -t40*t43*t44*t45-t40*t83*t121*t204*(1.0D0/2.0D0)+sig
     &nsin*t59*t69*t72*t192*t199
      B0(1,3,4,3) = t511-t40*t83*t125*t204*(1.0D0/2.0D0)+signsin*t59*t69
     &*t72*t196*t199
      B0(2,1,1,1) = t239-t59*t61*t237+t43*t48*t60*t250*(1.0D0/2.0D0)+t48
     &*t61*t63*t242-signsin*t59*t69*t72*t131*t242+signsin*t69*t72*t131*t
     &253*t254
      B0(2,1,1,2) = t247-t40*t43*t244-t43*t60*t68*t250*(1.0D0/2.0D0)-t61
     &*t63*t68*t242+signsin*t59*t69*t72*t139*t242-signsin*t69*t72*t139*t
     &253*t254
      B0(2,1,1,3) = t259-t40*t43*t256-t43*t60*t76*t250*(1.0D0/2.0D0)-t61
     &*t63*t76*t242+signsin*t59*t69*t72*t144*t242-signsin*t69*t72*t144*t
     &253*t254
      B0(2,1,2,1) = -t40*t43*(t292+t349)+t59*t61*(t293+t351)+t61*t254*(t
     &294+t353)-t43*t60*t79*t250*(1.0D0/2.0D0)-t61*t63*t79*t242-t61*t82*
     &t253*t262+t40*t82*t83*(t248+t249-t260-t261)*(1.0D0/2.0D0)-signsin*
     &t59*t69*t72*t150*t242+signsin*t69*t72*t253*t254*(t145+t146+t147+t1
     &48+t149-t208-t209-t210-t211-t212)
      B0(2,1,2,2) = t290-t3*t5*t59*t61-t13*t15*t61*t254+t61*t63*t85*t242
     &-t40*t83*t88*t250*(1.0D0/2.0D0)+t61*t88*t253*t262+t43*t60*t85*(t24
     &8+t249-t260-t261)*(1.0D0/2.0D0)+signsin*t59*t69*t72*t242*t263-sign
     &sin*t69*t72*t253*t254*t263
      B0(2,1,2,3) = t340-t5*t9*t59*t61-t15*t19*t61*t254+t61*t63*t92*t242
     &-t40*t83*t95*t250*(1.0D0/2.0D0)+t61*t95*t253*t262+t43*t60*t92*(t24
     &8+t249-t260-t261)*(1.0D0/2.0D0)+signsin*t59*t69*t72*t242*t264-sign
     &sin*t69*t72*t253*t254*t264
      B0(2,1,3,1) = t386-t59*t61*t383-t61*t254*t384+t61*t63*t98*t242-t40
     &*t83*t101*t250*(1.0D0/2.0D0)+t61*t101*t253*t262+t43*t60*t98*(t248+
     &t249-t260-t261)*(1.0D0/2.0D0)+signsin*t59*t69*t72*t242*t265-signsi
     &n*t69*t72*t253*t254*t265
      B0(2,1,3,2) = -t40*t43*t428-t59*t61*t424-t61*t254*t425-t43*t60*t10
     &3*t250*(1.0D0/2.0D0)-t61*t63*t103*t242-t61*t106*t253*t262+t40*t83*
     &t106*(t248+t249-t260-t261)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t175*
     &t242+signsin*t69*t72*t253*t254*(t170+t171+t172+t173+t174-t217-t218
     &-t219-t220-t221)
      B0(2,1,3,3) = -t40*t43*t473-t59*t61*t469-t61*t254*t470-t43*t60*t11
     &0*t250*(1.0D0/2.0D0)-t61*t63*t110*t242-t61*t113*t253*t262+t40*t83*
     &t113*(t248+t249-t260-t261)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t181*
     &t242+signsin*t69*t72*t253*t254*(t176+t177+t178+t179+t180-t226-t227
     &-t228-t229-t230)
      B0(2,1,4,1) = t498-t40*t43*t237-t61*t119*t253*t262+t40*t83*t119*(t
     &248+t249-t260-t261)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t185*t242+si
     &gnsin*t69*t72*t185*t253*t254
      B0(2,1,4,2) = t505-t61*t244*t254-t40*t83*t121*t250*(1.0D0/2.0D0)+t
     &61*t121*t253*t262+signsin*t59*t69*t72*t192*t242-signsin*t69*t72*t1
     &92*t253*t254
      B0(2,1,4,3) = t512-t61*t254*t256-t40*t83*t125*t250*(1.0D0/2.0D0)+t
     &61*t125*t253*t262+signsin*t59*t69*t72*t196*t242-signsin*t69*t72*t1
     &96*t253*t254
      B0(2,2,1,1) = t270-t40*t43*t267-t43*t48*t60*t291*(1.0D0/2.0D0)-t48
     &*t61*t63*t272+signsin*t59*t69*t72*t131*t272-signsin*t69*t72*t131*t
     &254*t280
      B0(2,2,1,2) = t275-t59*t61*t273+t43*t60*t68*t291*(1.0D0/2.0D0)+t61
     &*t63*t68*t272-signsin*t59*t69*t72*t139*t272+signsin*t69*t72*t139*t
     &254*t280
      B0(2,2,1,3) = t285-t40*t43*t282+t43*t60*t76*t291*(1.0D0/2.0D0)+t61
     &*t63*t76*t272-signsin*t59*t69*t72*t144*t272+signsin*t69*t72*t144*t
     &254*t280
      B0(2,2,2,1) = t290-t3*t5*t59*t61-t13*t15*t61*t254+t43*t60*t79*t291
     &*(1.0D0/2.0D0)+t61*t63*t79*t272-t40*t82*t83*t291*(1.0D0/2.0D0)+t61
     &*t82*t262*t280-signsin*t69*t72*t150*t254*t280+signsin*t59*t69*t72*
     &t272*(t145+t146+t147+t148+t149-t208-t209-t210-t211-t212)
      B0(2,2,2,2) = -t40*t43*(t292+t348)+t59*t61*(t293+t350)+t61*t254*(t
     &294+t352)-t43*t60*t85*t291*(1.0D0/2.0D0)-t61*t63*t85*t272-t61*t88*
     &t262*t280+t40*t83*t88*(t276+t277-t295-t296)*(1.0D0/2.0D0)-signsin*
     &t59*t69*t72*t263*t272+signsin*t69*t72*t254*t280*(t153+t154+t155+t1
     &56+t157-t297-t298-t299-t300-t301)
      B0(2,2,2,3) = t345-t3*t9*t59*t61-t13*t19*t61*t254-t43*t60*t92*t291
     &*(1.0D0/2.0D0)-t61*t63*t92*t272-t61*t95*t262*t280+t40*t83*t95*(t27
     &6+t277-t295-t296)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t264*t272+sign
     &sin*t69*t72*t254*t280*(t158+t159+t160+t161+t162-t302-t303-t304-t30
     &5-t306)
      B0(2,2,3,1) = -t40*t43*t392-t59*t61*t388-t61*t254*t389-t43*t60*t98
     &*t291*(1.0D0/2.0D0)-t61*t63*t98*t272-t61*t101*t262*t280+t40*t83*t1
     &01*(t276+t277-t295-t296)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t265*t2
     &72+signsin*t69*t72*t254*t280*(t163+t164+t165+t166+t167-t307-t308-t
     &309-t310-t311)
      B0(2,2,3,2) = t432-t59*t61*t429-t61*t254*t430+t61*t63*t103*t272-t4
     &0*t83*t106*t291*(1.0D0/2.0D0)+t61*t106*t262*t280+t43*t60*t103*(t27
     &6+t277-t295-t296)*(1.0D0/2.0D0)-signsin*t69*t72*t175*t254*t280+sig
     &nsin*t59*t69*t72*t272*(t170+t171+t172+t173+t174-t217-t218-t219-t22
     &0-t221)
      B0(2,2,3,3) = -t40*t43*t478-t59*t61*t474-t61*t254*t475+t61*t63*t11
     &0*t272-t40*t83*t113*t291*(1.0D0/2.0D0)+t61*t113*t262*t280+t43*t60*
     &t110*(t276+t277-t295-t296)*(1.0D0/2.0D0)-signsin*t69*t72*t181*t254
     &*t280+signsin*t59*t69*t72*t272*(t176+t177+t178+t179+t180-t226-t227
     &-t228-t229-t230)
      B0(2,2,4,1) = t499-t61*t254*t267-t40*t83*t119*t291*(1.0D0/2.0D0)+t
     &61*t119*t262*t280+signsin*t59*t69*t72*t185*t272-signsin*t69*t72*t1
     &85*t254*t280
      B0(2,2,4,2) = t506-t40*t43*t273-t61*t121*t262*t280+t40*t83*t121*(t
     &276+t277-t295-t296)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t192*t272+si
     &gnsin*t69*t72*t192*t254*t280
      B0(2,2,4,3) = t513-t61*t254*t282-t61*t125*t262*t280+t40*t83*t125*(
     &t276+t277-t295-t296)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t196*t272+s
     &ignsin*t69*t72*t196*t254*t280
      B0(2,3,1,1) = t320-t40*t43*t317-t43*t48*t60*t341*(1.0D0/2.0D0)-t48
     &*t61*t63*t323+signsin*t59*t69*t72*t131*t323-signsin*t69*t72*t131*t
     &254*t333
      B0(2,3,1,2) = t328-t40*t43*t325+t43*t60*t68*t341*(1.0D0/2.0D0)+t61
     &*t63*t68*t323-signsin*t59*t69*t72*t139*t323+signsin*t69*t72*t139*t
     &254*t333
      B0(2,3,1,3) = t336-t59*t61*t334+t43*t60*t76*t341*(1.0D0/2.0D0)+t61
     &*t63*t76*t323-signsin*t59*t69*t72*t144*t323+signsin*t69*t72*t144*t
     &254*t333
      B0(2,3,2,1) = t340-t5*t9*t59*t61-t15*t19*t61*t254+t43*t60*t79*t341
     &*(1.0D0/2.0D0)+t61*t63*t79*t323-t40*t82*t83*t341*(1.0D0/2.0D0)+t61
     &*t82*t262*t333-signsin*t69*t72*t150*t254*t333+signsin*t59*t69*t72*
     &t323*(t145+t146+t147+t148+t149-t208-t209-t210-t211-t212)
      B0(2,3,2,2) = t345-t3*t9*t59*t61-t13*t19*t61*t254-t43*t60*t85*t341
     &*(1.0D0/2.0D0)-t61*t63*t85*t323-t61*t88*t262*t333+t40*t83*t88*(t32
     &9+t330-t346-t347)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t263*t323+sign
     &sin*t69*t72*t254*t333*(t153+t154+t155+t156+t157-t297-t298-t299-t30
     &0-t301)
      B0(2,3,2,3) = -t40*t43*(t348+t349)+t59*t61*(t350+t351)+t61*t254*(t
     &352+t353)-t43*t60*t92*t341*(1.0D0/2.0D0)-t61*t63*t92*t323-t61*t95*
     &t262*t333+t40*t83*t95*(t329+t330-t346-t347)*(1.0D0/2.0D0)-signsin*
     &t59*t69*t72*t264*t323+signsin*t69*t72*t254*t333*(t158+t159+t160+t1
     &61+t162-t302-t303-t304-t305-t306)
      B0(2,3,3,1) = -t40*t43*t399-t59*t61*t395-t61*t254*t396-t43*t60*t98
     &*t341*(1.0D0/2.0D0)-t61*t63*t98*t323-t61*t101*t262*t333+t40*t83*t1
     &01*(t329+t330-t346-t347)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t265*t3
     &23+signsin*t69*t72*t254*t333*(t163+t164+t165+t166+t167-t307-t308-t
     &309-t310-t311)
      B0(2,3,3,2) = -t40*t43*t437-t59*t61*t433-t61*t254*t434+t61*t63*t10
     &3*t323-t40*t83*t106*t341*(1.0D0/2.0D0)+t61*t106*t262*t333+t43*t60*
     &t103*(t329+t330-t346-t347)*(1.0D0/2.0D0)-signsin*t69*t72*t175*t254
     &*t333+signsin*t59*t69*t72*t323*(t170+t171+t172+t173+t174-t217-t218
     &-t219-t220-t221)
      B0(2,3,3,3) = t482-t59*t61*t479-t61*t254*t480+t61*t63*t110*t323-t4
     &0*t83*t113*t341*(1.0D0/2.0D0)+t61*t113*t262*t333+t43*t60*t110*(t32
     &9+t330-t346-t347)*(1.0D0/2.0D0)-signsin*t69*t72*t181*t254*t333+sig
     &nsin*t59*t69*t72*t323*(t176+t177+t178+t179+t180-t226-t227-t228-t22
     &9-t230)
      B0(2,3,4,1) = t500-t61*t254*t317-t40*t83*t119*t341*(1.0D0/2.0D0)+t
     &61*t119*t262*t333+signsin*t59*t69*t72*t185*t323-signsin*t69*t72*t1
     &85*t254*t333
      B0(2,3,4,2) = t507-t61*t254*t325-t61*t121*t262*t333+t40*t83*t121*(
     &t329+t330-t346-t347)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t192*t323+s
     &ignsin*t69*t72*t192*t254*t333
      B0(2,3,4,3) = t514-t40*t43*t334-t61*t125*t262*t333+t40*t83*t125*(t
     &329+t330-t346-t347)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t196*t323+si
     &gnsin*t69*t72*t196*t254*t333
      B0(3,1,1,1) = t363-t40*t43*t364-t43*t48*t60*t387*(1.0D0/2.0D0)-t48
     &*t61*t63*t367+signsin*t59*t69*t72*t131*t367-signsin*t69*t72*t131*t
     &254*t377
      B0(3,1,1,2) = t370-t59*t61*t372+t43*t60*t68*t387*(1.0D0/2.0D0)+t61
     &*t63*t68*t367-signsin*t59*t69*t72*t139*t367+signsin*t69*t72*t139*t
     &254*t377
      B0(3,1,1,3) = t380-t59*t61*t382+t43*t60*t76*t387*(1.0D0/2.0D0)+t61
     &*t63*t76*t367-signsin*t59*t69*t72*t144*t367+signsin*t69*t72*t144*t
     &254*t377
      B0(3,1,2,1) = t386-t59*t61*t383-t61*t254*t384+t43*t60*t79*t387*(1.
     &0D0/2.0D0)+t61*t63*t79*t367-t40*t82*t83*t387*(1.0D0/2.0D0)+t61*t82
     &*t262*t377-signsin*t69*t72*t150*t254*t377+signsin*t59*t69*t72*t367
     &*(t145+t146+t147+t148+t149-t208-t209-t210-t211-t212)
      B0(3,1,2,2) = -t40*t43*t392-t59*t61*t388-t61*t254*t389-t43*t60*t85
     &*t387*(1.0D0/2.0D0)-t61*t63*t85*t367-t61*t88*t262*t377+t40*t83*t88
     &*(t373+t374-t393-t394)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t263*t367
     &+signsin*t69*t72*t254*t377*(t153+t154+t155+t156+t157-t297-t298-t29
     &9-t300-t301)
      B0(3,1,2,3) = -t40*t43*t399-t59*t61*t395-t61*t254*t396-t43*t60*t92
     &*t387*(1.0D0/2.0D0)-t61*t63*t92*t367-t61*t95*t262*t377+t40*t83*t95
     &*(t373+t374-t393-t394)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t264*t367
     &+signsin*t69*t72*t254*t377*(t158+t159+t160+t161+t162-t302-t303-t30
     &4-t305-t306)
      B0(3,1,3,1) = -t40*t43*(t442+t492)+t59*t61*(t443+t494)+t61*t254*(t
     &444+t496)-t43*t60*t98*t387*(1.0D0/2.0D0)-t61*t63*t98*t367-t61*t101
     &*t262*t377+t40*t83*t101*(t373+t374-t393-t394)*(1.0D0/2.0D0)-signsi
     &n*t59*t69*t72*t265*t367+signsin*t69*t72*t254*t377*(t163+t164+t165+
     &t166+t167-t307-t308-t309-t310-t311)
      B0(3,1,3,2) = t441-t2*t6*t59*t61-t12*t16*t61*t254+t61*t63*t103*t36
     &7-t40*t83*t106*t387*(1.0D0/2.0D0)+t61*t106*t262*t377+t43*t60*t103*
     &(t373+t374-t393-t394)*(1.0D0/2.0D0)-signsin*t69*t72*t175*t254*t377
     &+signsin*t59*t69*t72*t367*(t170+t171+t172+t173+t174-t217-t218-t219
     &-t220-t221)
      B0(3,1,3,3) = t486-t2*t10*t59*t61-t12*t20*t61*t254+t61*t63*t110*t3
     &67-t40*t83*t113*t387*(1.0D0/2.0D0)+t61*t113*t262*t377+t43*t60*t110
     &*(t373+t374-t393-t394)*(1.0D0/2.0D0)-signsin*t69*t72*t181*t254*t37
     &7+signsin*t59*t69*t72*t367*(t176+t177+t178+t179+t180-t226-t227-t22
     &8-t229-t230)
      B0(3,1,4,1) = t501-t61*t254*t364-t40*t83*t119*t387*(1.0D0/2.0D0)+t
     &61*t119*t262*t377+signsin*t59*t69*t72*t185*t367-signsin*t69*t72*t1
     &85*t254*t377
      B0(3,1,4,2) = t508-t40*t43*t372-t61*t121*t262*t377+t40*t83*t121*(t
     &373+t374-t393-t394)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t192*t367+si
     &gnsin*t69*t72*t192*t254*t377
      B0(3,1,4,3) = t515-t40*t43*t382-t61*t125*t262*t377+t40*t83*t125*(t
     &373+t374-t393-t394)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t196*t367+si
     &gnsin*t69*t72*t196*t254*t377
      B0(3,2,1,1) = t402-t59*t61*t404+t43*t48*t60*t412*(1.0D0/2.0D0)+t48
     &*t61*t63*t406-signsin*t59*t69*t72*t131*t406+signsin*t69*t72*t131*t
     &254*t415
      B0(3,2,1,2) = t408-t40*t43*t409-t43*t60*t68*t412*(1.0D0/2.0D0)-t61
     &*t63*t68*t406+signsin*t59*t69*t72*t139*t406-signsin*t69*t72*t139*t
     &254*t415
      B0(3,2,1,3) = t418-t59*t61*t420-t43*t60*t76*t412*(1.0D0/2.0D0)-t61
     &*t63*t76*t406+signsin*t59*t69*t72*t144*t406-signsin*t69*t72*t144*t
     &254*t415
      B0(3,2,2,1) = -t40*t43*t428-t59*t61*t424-t61*t254*t425-t43*t60*t79
     &*t412*(1.0D0/2.0D0)-t61*t63*t79*t406-t61*t82*t262*t415+t40*t82*t83
     &*(t410+t411-t421-t422)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t150*t406
     &+signsin*t69*t72*t254*t415*(t145+t146+t147+t148+t149-t208-t209-t21
     &0-t211-t212)
      B0(3,2,2,2) = t432-t59*t61*t429-t61*t254*t430+t61*t63*t85*t406-t40
     &*t83*t88*t412*(1.0D0/2.0D0)+t61*t88*t262*t415+t43*t60*t85*(t410+t4
     &11-t421-t422)*(1.0D0/2.0D0)-signsin*t69*t72*t254*t263*t415+signsin
     &*t59*t69*t72*t406*(t153+t154+t155+t156+t157-t297-t298-t299-t300-t3
     &01)
      B0(3,2,2,3) = -t40*t43*t437-t59*t61*t433-t61*t254*t434+t61*t63*t92
     &*t406-t40*t83*t95*t412*(1.0D0/2.0D0)+t61*t95*t262*t415+t43*t60*t92
     &*(t410+t411-t421-t422)*(1.0D0/2.0D0)-signsin*t69*t72*t254*t264*t41
     &5+signsin*t59*t69*t72*t406*(t158+t159+t160+t161+t162-t302-t303-t30
     &4-t305-t306)
      B0(3,2,3,1) = t441-t2*t6*t59*t61-t12*t16*t61*t254+t61*t63*t98*t406
     &-t40*t83*t101*t412*(1.0D0/2.0D0)+t61*t101*t262*t415+t43*t60*t98*(t
     &410+t411-t421-t422)*(1.0D0/2.0D0)-signsin*t69*t72*t254*t265*t415+s
     &ignsin*t59*t69*t72*t406*(t163+t164+t165+t166+t167-t307-t308-t309-t
     &310-t311)
      B0(3,2,3,2) = -t40*t43*(t442+t491)+t59*t61*(t443+t493)+t61*t254*(t
     &444+t495)-t43*t60*t103*t412*(1.0D0/2.0D0)-t61*t63*t103*t406-t61*t1
     &06*t262*t415+t40*t83*t106*(t410+t411-t421-t422)*(1.0D0/2.0D0)-sign
     &sin*t59*t69*t72*t175*t406+signsin*t69*t72*t254*t415*(t170+t171+t17
     &2+t173+t174-t217-t218-t219-t220-t221)
      B0(3,2,3,3) = t490-t6*t10*t59*t61-t16*t20*t61*t254-t43*t60*t110*t4
     &12*(1.0D0/2.0D0)-t61*t63*t110*t406-t61*t113*t262*t415+t40*t83*t113
     &*(t410+t411-t421-t422)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t181*t406
     &+signsin*t69*t72*t254*t415*(t176+t177+t178+t179+t180-t226-t227-t22
     &8-t229-t230)
      B0(3,2,4,1) = t502-t40*t43*t404-t61*t119*t262*t415+t40*t83*t119*(t
     &410+t411-t421-t422)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t185*t406+si
     &gnsin*t69*t72*t185*t254*t415
      B0(3,2,4,2) = t509-t61*t254*t409-t40*t83*t121*t412*(1.0D0/2.0D0)+t
     &61*t121*t262*t415+signsin*t59*t69*t72*t192*t406-signsin*t69*t72*t1
     &92*t254*t415
      B0(3,2,4,3) = t516-t40*t43*t420-t40*t83*t125*t412*(1.0D0/2.0D0)+t6
     &1*t125*t262*t415+signsin*t59*t69*t72*t196*t406-signsin*t69*t72*t19
     &6*t254*t415
      B0(3,3,1,1) = t447-t59*t61*t449+t43*t48*t60*t460*(1.0D0/2.0D0)+t48
     &*t61*t63*t452-signsin*t59*t69*t72*t131*t452+signsin*t69*t72*t131*t
     &254*t463
      B0(3,3,1,2) = t455-t59*t61*t457-t43*t60*t68*t460*(1.0D0/2.0D0)-t61
     &*t63*t68*t452+signsin*t59*t69*t72*t139*t452-signsin*t69*t72*t139*t
     &254*t463
      B0(3,3,1,3) = t465-t40*t43*t466-t43*t60*t76*t460*(1.0D0/2.0D0)-t61
     &*t63*t76*t452+signsin*t59*t69*t72*t144*t452-signsin*t69*t72*t144*t
     &254*t463
      B0(3,3,2,1) = -t40*t43*t473-t59*t61*t469-t61*t254*t470-t43*t60*t79
     &*t460*(1.0D0/2.0D0)-t61*t63*t79*t452-t61*t82*t262*t463+t40*t82*t83
     &*(t458+t459-t467-t468)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t150*t452
     &+signsin*t69*t72*t254*t463*(t145+t146+t147+t148+t149-t208-t209-t21
     &0-t211-t212)
      B0(3,3,2,2) = -t40*t43*t478-t59*t61*t474-t61*t254*t475+t61*t63*t85
     &*t452-t40*t83*t88*t460*(1.0D0/2.0D0)+t61*t88*t262*t463+t43*t60*t85
     &*(t458+t459-t467-t468)*(1.0D0/2.0D0)-signsin*t69*t72*t254*t263*t46
     &3+signsin*t59*t69*t72*t452*(t153+t154+t155+t156+t157-t297-t298-t29
     &9-t300-t301)
      B0(3,3,2,3) = t482-t59*t61*t479-t61*t254*t480+t61*t63*t92*t452-t40
     &*t83*t95*t460*(1.0D0/2.0D0)+t61*t95*t262*t463+t43*t60*t92*(t458+t4
     &59-t467-t468)*(1.0D0/2.0D0)-signsin*t69*t72*t254*t264*t463+signsin
     &*t59*t69*t72*t452*(t158+t159+t160+t161+t162-t302-t303-t304-t305-t3
     &06)
      B0(3,3,3,1) = t486-t2*t10*t59*t61-t12*t20*t61*t254+t61*t63*t98*t45
     &2-t40*t83*t101*t460*(1.0D0/2.0D0)+t61*t101*t262*t463+t43*t60*t98*(
     &t458+t459-t467-t468)*(1.0D0/2.0D0)-signsin*t69*t72*t254*t265*t463+
     &signsin*t59*t69*t72*t452*(t163+t164+t165+t166+t167-t307-t308-t309-
     &t310-t311)
      B0(3,3,3,2) = t490-t6*t10*t59*t61-t16*t20*t61*t254-t43*t60*t103*t4
     &60*(1.0D0/2.0D0)-t61*t63*t103*t452-t61*t106*t262*t463+t40*t83*t106
     &*(t458+t459-t467-t468)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t175*t452
     &+signsin*t69*t72*t254*t463*(t170+t171+t172+t173+t174-t217-t218-t21
     &9-t220-t221)
      B0(3,3,3,3) = -t40*t43*(t491+t492)+t59*t61*(t493+t494)+t61*t254*(t
     &495+t496)-t43*t60*t110*t460*(1.0D0/2.0D0)-t61*t63*t110*t452-t61*t1
     &13*t262*t463+t40*t83*t113*(t458+t459-t467-t468)*(1.0D0/2.0D0)-sign
     &sin*t59*t69*t72*t181*t452+signsin*t69*t72*t254*t463*(t176+t177+t17
     &8+t179+t180-t226-t227-t228-t229-t230)
      B0(3,3,4,1) = t503-t40*t43*t449-t61*t119*t262*t463+t40*t83*t119*(t
     &458+t459-t467-t468)*(1.0D0/2.0D0)-signsin*t59*t69*t72*t185*t452+si
     &gnsin*t69*t72*t185*t254*t463
      B0(3,3,4,2) = t510-t40*t43*t457-t40*t83*t121*t460*(1.0D0/2.0D0)+t6
     &1*t121*t262*t463+signsin*t59*t69*t72*t192*t452-signsin*t69*t72*t19
     &2*t254*t463
      B0(3,3,4,3) = t517-t61*t254*t466-t40*t83*t125*t460*(1.0D0/2.0D0)+t
     &61*t125*t262*t463+signsin*t59*t69*t72*t196*t452-signsin*t69*t72*t1
     &96*t254*t463
      B0(4,1,1,1) = t497-t43*t48*t54*t60*(1.0D0/2.0D0)+signsin*t66*t69*t
     &72*t131*t254
      B0(4,1,1,2) = -t40*t43*t44*t62+t43*t54*t60*t68*(1.0D0/2.0D0)-signs
     &in*t66*t69*t72*t139*t254
      B0(4,1,1,3) = -t40*t43*t45*t62+t43*t54*t60*t76*(1.0D0/2.0D0)-signs
     &in*t66*t69*t72*t144*t254
      B0(4,1,2,1) = t498-t40*t43*t237+t43*t54*t60*t79*(1.0D0/2.0D0)-t40*
     &t54*t82*t83*(1.0D0/2.0D0)-t61*t66*t82*t262+signsin*t66*t69*t72*t25
     &4*(t145+t146+t147+t148+t149-t208-t209-t210-t211-t212)
      B0(4,1,2,2) = t499-t61*t254*t267-t43*t54*t60*t85*(1.0D0/2.0D0)+t40
     &*t54*t83*t88*(1.0D0/2.0D0)+t61*t66*t88*t262-signsin*t66*t69*t72*t2
     &54*t263
      B0(4,1,2,3) = t500-t61*t254*t317-t43*t54*t60*t92*(1.0D0/2.0D0)+t40
     &*t54*t83*t95*(1.0D0/2.0D0)+t61*t66*t95*t262-signsin*t66*t69*t72*t2
     &54*t264
      B0(4,1,3,1) = t501-t61*t254*t364-t43*t54*t60*t98*(1.0D0/2.0D0)+t40
     &*t54*t83*t101*(1.0D0/2.0D0)+t61*t66*t101*t262-signsin*t66*t69*t72*
     &t254*t265
      B0(4,1,3,2) = t502-t40*t43*t404+t43*t54*t60*t103*(1.0D0/2.0D0)-t40
     &*t54*t83*t106*(1.0D0/2.0D0)-t61*t66*t106*t262+signsin*t66*t69*t72*
     &t254*(t170+t171+t172+t173+t174-t217-t218-t219-t220-t221)
      B0(4,1,3,3) = t503-t40*t43*t449+t43*t54*t60*t110*(1.0D0/2.0D0)-t40
     &*t54*t83*t113*(1.0D0/2.0D0)-t61*t66*t113*t262+signsin*t66*t69*t72*
     &t254*(t176+t177+t178+t179+t180-t226-t227-t228-t229-t230)
      B0(4,1,4,1) = t61*t116*t254-t40*t54*t83*t119*(1.0D0/2.0D0)-t61*t66
     &*t119*t262+signsin*t66*t69*t72*t185*t254
      B0(4,1,4,2) = t40*t54*t83*t121*(1.0D0/2.0D0)-t44*t61*t62*t254+t61*
     &t66*t121*t262-signsin*t66*t69*t72*t192*t254
      B0(4,1,4,3) = t40*t54*t83*t125*(1.0D0/2.0D0)-t45*t61*t62*t254+t61*
     &t66*t125*t262-signsin*t66*t69*t72*t196*t254
      B0(4,2,1,1) = -t40*t43*t44*t62+t43*t48*t60*t127*(1.0D0/2.0D0)-sign
     &sin*t69*t72*t131*t134*t254
      B0(4,2,1,2) = t504-t43*t60*t68*t127*(1.0D0/2.0D0)+signsin*t69*t72*
     &t134*t139*t254
      B0(4,2,1,3) = -t40*t43*t44*t45-t43*t60*t76*t127*(1.0D0/2.0D0)+sign
     &sin*t69*t72*t134*t144*t254
      B0(4,2,2,1) = t505-t61*t244*t254-t43*t60*t79*t127*(1.0D0/2.0D0)+t4
     &0*t82*t83*t127*(1.0D0/2.0D0)+t61*t82*t134*t262-signsin*t69*t72*t13
     &4*t150*t254
      B0(4,2,2,2) = t506-t40*t43*t273+t43*t60*t85*t127*(1.0D0/2.0D0)-t40
     &*t83*t88*t127*(1.0D0/2.0D0)-t61*t88*t134*t262+signsin*t69*t72*t134
     &*t254*(t153+t154+t155+t156+t157-t297-t298-t299-t300-t301)
      B0(4,2,2,3) = t507-t61*t254*t325+t43*t60*t92*t127*(1.0D0/2.0D0)-t4
     &0*t83*t95*t127*(1.0D0/2.0D0)-t61*t95*t134*t262+signsin*t69*t72*t13
     &4*t254*(t158+t159+t160+t161+t162-t302-t303-t304-t305-t306)
      B0(4,2,3,1) = t508-t40*t43*t372+t43*t60*t98*t127*(1.0D0/2.0D0)-t40
     &*t83*t101*t127*(1.0D0/2.0D0)-t61*t101*t134*t262+signsin*t69*t72*t1
     &34*t254*(t163+t164+t165+t166+t167-t307-t308-t309-t310-t311)
      B0(4,2,3,2) = t509-t61*t254*t409-t43*t60*t103*t127*(1.0D0/2.0D0)+t
     &40*t83*t106*t127*(1.0D0/2.0D0)+t61*t106*t134*t262-signsin*t69*t72*
     &t134*t175*t254
      B0(4,2,3,3) = t510-t40*t43*t457-t43*t60*t110*t127*(1.0D0/2.0D0)+t4
     &0*t83*t113*t127*(1.0D0/2.0D0)+t61*t113*t134*t262-signsin*t69*t72*t
     &134*t181*t254
      B0(4,2,4,1) = t40*t83*t119*t127*(1.0D0/2.0D0)-t44*t61*t62*t254+t61
     &*t119*t134*t262-signsin*t69*t72*t134*t185*t254
      B0(4,2,4,2) = t61*t187*t254-t40*t83*t121*t127*(1.0D0/2.0D0)-t61*t1
     &21*t134*t262+signsin*t69*t72*t134*t192*t254
      B0(4,2,4,3) = t40*t83*t125*t127*(-1.0D0/2.0D0)-t44*t45*t61*t254-t6
     &1*t125*t134*t262+signsin*t69*t72*t134*t196*t254
      B0(4,3,1,1) = -t40*t43*t45*t62+t43*t48*t60*t199*(1.0D0/2.0D0)-sign
     &sin*t69*t72*t131*t204*t254
      B0(4,3,1,2) = -t40*t43*t44*t45-t43*t60*t68*t199*(1.0D0/2.0D0)+sign
     &sin*t69*t72*t139*t204*t254
      B0(4,3,1,3) = t511-t43*t60*t76*t199*(1.0D0/2.0D0)+signsin*t69*t72*
     &t144*t204*t254
      B0(4,3,2,1) = t512-t61*t254*t256-t43*t60*t79*t199*(1.0D0/2.0D0)+t4
     &0*t82*t83*t199*(1.0D0/2.0D0)+t61*t82*t204*t262-signsin*t69*t72*t15
     &0*t204*t254
      B0(4,3,2,2) = t513-t61*t254*t282+t43*t60*t85*t199*(1.0D0/2.0D0)-t4
     &0*t83*t88*t199*(1.0D0/2.0D0)-t61*t88*t204*t262+signsin*t69*t72*t20
     &4*t254*(t153+t154+t155+t156+t157-t297-t298-t299-t300-t301)
      B0(4,3,2,3) = t514-t40*t43*t334+t43*t60*t92*t199*(1.0D0/2.0D0)-t40
     &*t83*t95*t199*(1.0D0/2.0D0)-t61*t95*t204*t262+signsin*t69*t72*t204
     &*t254*(t158+t159+t160+t161+t162-t302-t303-t304-t305-t306)
      B0(4,3,3,1) = t515-t40*t43*t382+t43*t60*t98*t199*(1.0D0/2.0D0)-t40
     &*t83*t101*t199*(1.0D0/2.0D0)-t61*t101*t204*t262+signsin*t69*t72*t2
     &04*t254*(t163+t164+t165+t166+t167-t307-t308-t309-t310-t311)
      B0(4,3,3,2) = t516-t40*t43*t420-t43*t60*t103*t199*(1.0D0/2.0D0)+t4
     &0*t83*t106*t199*(1.0D0/2.0D0)+t61*t106*t204*t262-signsin*t69*t72*t
     &175*t204*t254
      B0(4,3,3,3) = t517-t61*t254*t466-t43*t60*t110*t199*(1.0D0/2.0D0)+t
     &40*t83*t113*t199*(1.0D0/2.0D0)+t61*t113*t204*t262-signsin*t69*t72*
     &t181*t204*t254
      B0(4,3,4,1) = -t45*t61*t62*t254+t40*t83*t119*t199*(1.0D0/2.0D0)+t6
     &1*t119*t204*t262-signsin*t69*t72*t185*t204*t254
      B0(4,3,4,2) = -t44*t45*t61*t254-t40*t83*t121*t199*(1.0D0/2.0D0)-t6
     &1*t121*t204*t262+signsin*t69*t72*t192*t204*t254
      B0(4,3,4,3) = t61*t234*t254-t40*t83*t125*t199*(1.0D0/2.0D0)-t61*t1
     &25*t204*t262+signsin*t69*t72*t196*t204*t254
                  
      do j1=1,4
      do k1=1,3
         do j2=1,4
         do k2=1,4
            Ktmp(j1,k1,j2,k2)=Ktmp(j1,k1,j2,k2)+B0(j1,k1,j2,k2)*betab
         enddo
         enddo
      enddo
      enddo

      endif  !(Ktmpflag)

      Ftmp(1:4,1:3)=betab*Ftmp(1:4,1:3)
      return
      end subroutine
!======================================================================
      subroutine  point2trimesh(
     &                 DD
     &                ,PP
     &                ,BB
     &                ,FF
     &                ,nfaces,nedges,nvertices,nQPoints
     &                ,faces,edges
     &                ,vertices_x,vertices_y,vertices_z
     &                ,normals,qPoints) !, max_distance,useSubSurface

       real vertices_x(nvertices)
     &     ,vertices_y(nvertices)
     &     ,vertices_z(nvertices)
       integer nfaces,nvertices,nedges,nQPoints
       real vertices(nvertices,3),normals(nfaces,3)
     &      ,qpoints(nQPoints,1:3)
     &      ,BB(nQPoints,3),PP(nQPoints,3),DD(nQPoints)
       integer faces(nfaces,3),FF(nQPoints),edges(nedges,2)
       
       real P(3),B(3),d,point(1:3)
       integer i,f
      vertices(1:nvertices,1)=vertices_x(1:nvertices)
      vertices(1:nvertices,2)=vertices_y(1:nvertices)
      vertices(1:nvertices,3)=vertices_z(1:nvertices)
      do i = 1,nQPoints
            ! Determine the surface points
            point(1:3) = qPoints(i,1:3) ! (1 x 3) query point
            call processPoint(d,p,f,b
     &                ,nfaces,nedges,nvertices
     &                ,faces,edges,vertices,point,normals)
            DD(i) = d
            PP(i,1:3) = p(1:3)
            BB(i,1:3)=b(1:3)
            FF(i)=f
       enddo
       return
       end subroutine
!======================================================================
      subroutine processPoint(Dd,Pp,FF,BB
     &                ,nfaces,nedges,nvertices
     &                ,edges,faces,vertices,point,normals)
       integer nfaces,nvertices,nedges
       real vertices(nvertices,3),normals(nfaces,3),point(1:3)
     &      ,BB(3),PP(3),DD
       integer faces(nfaces,3),V(3),FF,edges(nedges,2)

       real dmin, r1(3),r2(3),r3(3),vq(3),D(3),rD(3),P(3,3)
     &      ,r31r31,r21r21,r21r31,r31vq,r21vq,ds,B(3,3),vec(3),t,x0
       integer connectedFaces,i,iv0,iv1,iv2,ie0,ie1,ivplace(3),j,f(3)

      d(3) = 0.0     ! (distanceTypes x 1)
      p(1:3,1:3)=0.0 ! (distanceTypes x xyz)
      f(1:3)=0       ! (distanceTypes x 1)

      ! find nearest vertice
      call  distance_to_vertices(d(1),p(1,1:3),f(1),v(1),B(1,1:3)
     &                ,nfaces,nvertices
     &                ,faces,vertices,point,normals)
      ! d:  (1 x 1) signed distance to surface
      ! p:  (1 x 3) corresponding point on surface
      ! v:  (1 x 1) nearest vertex
      ! f:  (1 x 1) corresponding face
      ! b:  (1 x 3) Bary coordinate of the corresponding element

      ! find nearest point on all edges
      call  distance_to_edges(d(2),p(2,1:3),f(2),v(2),B(2,1:3)
     &                ,nfaces,nedges,nvertices
     &                ,faces,edges,vertices,point,normals)

      ! find nearest point on all surfaces
      call  distance_to_surfaces(d(3),p(3,1:3),f(3),v(3),B(3,1:3)
     &                ,nfaces,nvertices
     &                ,faces,vertices,point,normals)
      ! find minimum distance type
      I = MINLOC(d(1:3),1) 

      DD = d(I) 
      PP(1:3) = p(I,1:3)
      BB(1:3) = B(I,1:3)
      FF = f(I)
      return 
      end subroutine
!======================================================================
       subroutine distance_to_vertices(DD,PP,FF,V,BB
     &                ,nfaces,nvertices,faces,vertices,qPoint,normals)
       integer nfaces,nvertices
       real vertices(nvertices,3),normals(nfaces,3),qPoint(1:3)
     &      ,BB(3),PP(3),DD
       integer faces(nfaces,3),V,FF

       real dmin, r1(3),r2(3),r3(3),vq(3),D,rD(3),P(3),coefficients
     &      ,r31r31,r21r21,r21r31,r31vq,r21vq,ds,bary(3),vec(3),t,x0
       integer connectedFaces,i,iv0,iv1,iv2,ie0,ie1,ivplace(3),j

       ! find nearest vertex
       dmin=1.0e5
       connectedFaces=-9999
       do i=1,nvertices
          vec(1:3)=vertices(i,1:3) 
          vec(1:3)=vec(1:3)-qPoint(1:3)
          x0=dot_product(vec(1:3),vec(1:3))
          if(x0<dmin) then
            V=i
            DD=sqrt(x0)
            dmin=x0
          endif
       enddo
       PP=vertices(V,1:3)
       connectedFaces=0

       do i=1,nfaces
          if(any(faces(i,1:3).eq.V)) then
            connectedFaces=i
            r1(1:3) = vertices(V,1:3)
!     (#edges x 3) first point of every edge
              do j=1,3
                 r3(1:3)=vertices(faces(i,j),1:3)
                 vec(j) = dot_product(r3(1:3)-r1(1:3),r3(1:3)-r1(1:3))
!     (#faces x 1)
              enddo
              iv0 = MINLOC(vec(1:3),1)
              BB(1:3)=0.0
              BB(iv0)=1.0
            exit
          endif
       enddo
       if(connectedFaces .lt. 0) then
          stop
       endif
       coefficients
     &   =dot_product(normals(connectedFaces,1:3),qPoint(1:3)-P(1:3))
       FF=connectedFaces
       if(coefficients .lt. 0) DD=DD*coefficients
       return
       end subroutine
!======================================================================
       subroutine distance_to_edges(DD,PP,FF,V,BB
     &                ,nfaces,nedges,nvertices
     &                ,faces,edges,vertices,qPoint,normals)

       integer nfaces,nedges,nvertices
       real vertices(nvertices,3),normals(nfaces,3),qpoint(1:3)
     &      ,BB(3),PP(3),DD
       integer faces(nfaces,3),edges(nedges,2),V,FF

       real dmin, r1(3),r2(3),r3(3),vq(3),D,rD(3),P(3),coefficients
     &      ,r31r31,r21r21,r21r31,r31vq,r21vq,ds,bary(3),vec(3),t,x0
       integer connectedFaces,i,iv0,iv1,iv2,ie0,ie1,ivplace(3),j

       
       ! Intersection between tangent of edge lines and query point
       dmin=1.0e5
       connectedFaces=-9999
       do i=1,nedges
          r1(1:3) = vertices(edges(i,1),1:3)
!     (#edges x 3) first point of every edge
          r2(1:3) = vertices(edges(i,2),1:3)
!     (#edges x 3) second point of every edge
          t = dot_product( qPoint(1:3)-r1(1:3), r2(1:3)-r1(1:3))
          t =t/ dot_product(r2(1:3)-r1(1:3),r2(1:3)-r1(1:3))
!     (#edges x 1) location of intersection relative to r1 and r2
          if((t .le. 0) .or. (t .ge. 1) ) cycle 
            ! Distance between intersection and query point
          P(1:3) = r1(1:3) + t*(r2(1:3)-r1(1:3))
!     (#edges x 3) intersection points
          vec(1:3) = qPoint(1:3)-P(1:3)
!     (#edges x 3)
            x0=dot_product(vec(1:3),vec(1:3))
          if(x0<dmin) then
            V=i
            DD=sqrt(x0)
            dmin=x0
            PP(1:3)=P(1:3)
          endif
       enddo
       ie0=edges(V,1)
       ie1=edges(V,2)
!       ieplace=(/ie0,ie1/)
!       if(ie0>ie1) then
!          i=ie0
!          ie0=ie1
!          ie1=i
!          ieplace=(/ieplace(2),ieplace(1)/)
!       endif
       do i=1,nfaces
          iv0=faces(i,1)
          iv1=faces(i,2)   
          iv2=faces(i,3)
          ivplace=(/1,2,3/)
          if(any(ivplace .eq. ie0) .and. any(ivplace .eq. ie1)) then
              connectedFaces=i
              r1(1:3) = vertices(edges(V,1),1:3)
!     (#edges x 3) first point of every edge
              r2(1:3) = vertices(edges(V,2),1:3)
!     (#edges x 3) second point of every edge
              
              do j=1,3
                 r3(1:3)=vertices(faces(i,j),1:3)
                 vec(j) = dot_product(r3(1:3)-r1(1:3),r3(1:3)-r1(1:3))
!     (#faces x 1)
                 P(j)   = dot_product(r3(1:3)-r2(1:3),r3(1:3)-r2(1:3))
!     (#faces x 1)
              enddo
              iv0 = MINLOC(vec(1:3),1)
              iv1 = MINLOC(  P(1:3),1)
              iv2 = 6-iv0-iv1
              BB(iv0)=1.0-t
              BB(iv1)= t
              BB(iv2)=0.0
              exit
          endif
 
!          if(iv0>iv1) then
!           j=iv0
!           iv0=iv1
!           iv1=j
!           ivplace=(/ivplace(2),ivplace(1),ivplace(3)/)
!         endif     
!         if(iv0>iv2) then
!           j=iv0
!           iv0=iv2
!           iv2=j
!           ivplace=(/ivplace(3),ivplace(2),ivplace(1)/)
!         elseif(iv1>iv2) then
!           j=iv1
!           iv1=iv2
!           iv2=j
!           ivplace=(/ivplace(1),ivplace(3),ivplace(2)/)
!         endif  
!
!         if(ie0 .eq. iv0) then
!            if((ie1 .eq. iv1).or. (ie1 .eq. iv2)) then
!             connectedFaces=i
!             exit
!            endif
!         endif
        enddo
       if(connectedFaces .lt. 0) then
          stop
       endif
       FF=connectedFaces
    
       coefficients
     &   =dot_product(normals(connectedFaces,1:3),qPoint(1:3)-PP(1:3))
       if(coefficients .lt. 0) DD=DD*coefficients
       return
       end subroutine
!======================================================================
       subroutine distance_to_surfaces(DD,PP,FF,V,BB
     &                ,nfaces,nvertices
     &                ,faces,vertices,point,normals)
       integer nfaces,nvertices
       real vertices(nvertices,3),normals(nfaces,3),point(1:3)
     &      ,BB(3),PP(3),DD
       integer faces(nfaces,3),V,FF

       real dmin, r1(3),r2(3),r3(3),vq(3),D,rD(3),P(3)
     &      ,r31r31,r21r21,r21r31,r31vq,r21vq,ds,bary(3)
       integer connectedFaces,i
       dmin=1.0e8
       connectedFaces=-9999
       do i=1,nfaces
          r1(1:3) = vertices(faces(i,1),1:3)
!     (#faces x 3) ! 1st vertex of every face
          r2(1:3) = vertices(faces(i,2),1:3)
!     (#faces x 3) ! 2nd vertex of every face
          r3(1:3) = vertices(faces(i,3),1:3)
!     (#faces x 3) ! 3rd vertex of every face

         vq(1:3) = point(1:3)-r1(1:3)  ! (#faces x 3)
         D = dot_product(vq(1:3),normals(i,1:3))
!     (#faces x 1) distance to surface
         rD(1:3) = normals(i,1:3)*D
!     (#faces x 3) vector from surface to query point
         P(1:3) = point(1:3)-rD(1:3)
!     (#faces x 3) nearest point on surface; can be outside triangle

! find barycentric coordinates
         !(query point as linear combination of two edges)
         r31r31 = dot_product(r3(1:3)-r1(1:3),r3(1:3)-r1(1:3))
!     (#faces x 1)
         r21r21 = dot_product(r2(1:3)-r1(1:3),r2(1:3)-r1(1:3))
!     (#faces x 1)
         r21r31 = dot_product(r2(1:3)-r1(1:3),r3(1:3)-r1(1:3))
!     (#faces x 1)
         r31vq  = dot_product(r3(1:3)-r1(1:3),vq(1:3))
!     (#faces x 1) 
         r21vq  = dot_product(r2(1:3)-r1(1:3),vq(1:3))
!     (#faces x 1)

         ds = r31r31*r21r21 - r21r31**2               ! (#faces x 1)
         if(abs(ds) .le. 1.0d-6) cycle
         bary(1) = (r21r21*r31vq-r21r31*r21vq)/ds ! (#faces x 3)
         bary(2) = (r31r31*r21vq-r21r31*r31vq)/ds ! (#faces x 3)
         bary(3) = 1.0 - bary(1) - bary(2)        ! (#faces x 3)

         ! exclude intersections that are outside the triangle
         if((any(bary .le. 0.0)) .or. (any(bary .ge. 1.0)) )   cycle
!     (#faces x 1)
          if(ABS(ds)<dmin) then
            V=i
            dmin=ds
            BB(1:3)=bary(1:3)
            PP=P(1:3)
          endif
       enddo    
       DD=Dmin  
       FF=V
      end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
