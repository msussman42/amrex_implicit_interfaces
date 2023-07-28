           if(e4coef(i2g).gt.0.0) then
           do j=1,ns_ibm_r_fsh(i)
! boundary 1
           jq=1
              DX_im_j(1:3)=
     &         (   -GXtmp(1:3,jq+2,j)
     &         +4.0*GXtmp(1:3,jq+1,j)
     &         -3.0*GXtmp(1:3,jq,j)   )
     &         /(3.0*dsf_IBM_fsh(i,jq,j)-dsf_IBM_fsh(i,jq+1,j) )        !x,y

           jq=1
              DX_ip_j(1:3)=
     &         (   -GXtmp(1:3,jq+2,j+1)
     &         +4.0*GXtmp(1:3,jq+1,j+1)
     &         -3.0*GXtmp(1:3,jq,j+1)   )
     &         /(3.0*dsf_IBM_fsh(i,jq,j+1)-dsf_IBM_fsh(i,jq+1,j+1) )        !x,y+1

              DX_ip_j(1:3)=0.5*(DX_ip_j(1:3)+DX_im_j(1:3))

              DY_ip_j(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))/dsf2_IBM_fsh(i,jq,j) !x,y+1/2

              sigB(1,j,1)=0.0
              sigB(1,j,3)=Mem_Coef_fsh(i,jq,j,3)
     &              *(dot_product(DX_ip_j(1:3),DY_ip_j(1:3))
     &                  -TzeroB_fsh(i,1,j,3)  )
              SigB_d(1,j,3)=Mem_Coef_fsh(i,jq,j,3)

              if(e4coef(i2g) .lt. 50.0) then   
              sigB(1,j,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &            (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &        +e4coef(i2g)
     &           *(1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))        
     &                  -TzeroB_fsh(i,1,j,2)  )

              SigB_d(1,j,1)=Mem_Coef_fsh(i,jq,j,2)
     &              /sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))**3
              SigB_d(1,j,2)=Mem_Coef_fsh(i,jq,j,2)*e4coef(i2g)
     &              /sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))**3


              elseif(e4coef(i2g) .le. 100.0) then   
                     temp_ibm3= e4coef(i2g)-50.0          
                     sigB(1,j,2)=Mem_Coef_fsh(i,jq,j,1)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))-1.0)        
     &            /(sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))) )
     &                  -TzeroB_fsh(i,1,j,2)  )

                 SigB_d(1,j,1)=Mem_Coef_fsh(i,jq,j,2)
     &              /sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))**3
     &              +Mem_Coef_fsh(i,jq,j,2)*temp_ibm3
     &             *(sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))-1.0)        
     &             /(-sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))**3
                 SigB_d(1,j,2)=Mem_Coef_fsh(i,jq,j,2)*temp_ibm3
     &              /sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))
     &            /(sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))

              endif

! boundary 2
           jq=nq_ibm_r_fsh(i)
              DX_im_j(1:3)=
     &         (    GXtmp(1:3,jq-2,j)
     &         -4.0*GXtmp(1:3,jq-1,j)
     &         +3.0*GXtmp(1:3,jq,j)   )
     &         /(3.0*dsf_IBM_fsh(i,jq-1,j)- dsf_IBM_fsh(i,jq-2,j) )

              DY_i_jp(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))/dsf2_IBM_fsh(i,jq,j) !x,y+1/2

              DY_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &                     /dsf2_IBM_fsh(i,jq,j-1)                      !x,y-1/2

              DY_ip_j(1:3)=
     &        0.5*(DY_i_jp(1:3)+DY_i_jm(1:3))

              sigB(2,j,2)=0.0
              sigB(2,j,3)=Mem_Coef_fsh(i,jq,j,3)
     &              *(dot_product(DX_im_j(1:3),DY_ip_j(1:3))
     &                  -TzeroB_fsh(i,2,j,3)  )
              SigB_d(2,j,3)=Mem_Coef_fsh(i,jq,j,3)
              if(e4coef(i2g) .lt. 50.0) then

              sigB(2,j,1)=Mem_Coef_fsh(i,jq,j,1)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DX_im_j(1:3),DX_im_j(1:3))))
     &        +e4coef(i2g)
     &          *(1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))        
     &                  -TzeroB_fsh(i,2,j,1)  )

              SigB_d(2,j,1)=Mem_Coef_fsh(i,jq,j,1)
     &              /sqrt(dot_product(DX_im_j(1:3),DX_im_j(1:3)))**3
              SigB_d(2,j,2)=Mem_Coef_fsh(i,jq,j,1)*e4coef(i2g)
     &              /sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))**3

              elseif(e4coef(i2g) .le. 100.0) then   
                     temp_ibm3= e4coef(i2g)-50.0          
                     sigB(2,j,1)=Mem_Coef_fsh(i,jq,j,1)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DX_im_j(1:3),DX_im_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))-1.0)        
     &            /(sqrt(dot_product(DX_im_j(1:3),DX_im_j(1:3)))) )
     &                  -TzeroB_fsh(i,2,j,1)  )

                 SigB_d(2,j,1)=Mem_Coef_fsh(i,jq,j,1)
     &              /sqrt(dot_product(DX_im_j(1:3),DX_im_j(1:3)))**3
     &              +Mem_Coef_fsh(i,jq,j,1)*temp_ibm3
     &            *(sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))-1.0)        
     &             /(-sqrt(dot_product(DX_im_j(1:3),DX_im_j(1:3))))**3
                 SigB_d(2,j,2)=Mem_Coef_fsh(i,jq,j,1)*temp_ibm3
     &              /sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))
     &            /(sqrt(dot_product(DX_im_j(1:3),DX_im_j(1:3))))
              endif
           enddo


           do jq=1,nq_ibm_r_fsh(i)
! boundary 3
           j=1
              DY_i_jp(1:3)=
     &         (   -GXtmp(1:3,jq,j+2)
     &         +4.0*GXtmp(1:3,jq,j+1)
     &         -3.0*GXtmp(1:3,jq,j)   )
     &         /(3.0*dsf2_IBM_fsh(i,jq,j)-dsf2_IBM_fsh(i,jq,j+1))       !x+1/2,y

              DX_ip_j(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))/dsf_IBM_fsh(i,jq,j)  !x,y+1/2

              DX_im_j(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))
     &                     /dsf_IBM_fsh(i,jq-1,j)                      !x,y-1/2

              DX_i_jp(1:3)=
     &        0.5*(DX_ip_j(1:3)+DX_im_j(1:3))


              sigB(3,jq,1)=0.0
              sigB(3,jq,3)=Mem_Coef_fsh(i,jq,j,3)
     &              *(dot_product(DX_i_jp(1:3),DY_i_jp(1:3))
     &                  -TzeroB_fsh(i,3,jq,3)  )
              SigB_d(3,jq,3)=Mem_Coef_fsh(i,jq,j,3)
              if(e4coef(i2g) .lt. 50.0) then 
              sigB(3,jq,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))
     &        +e4coef(i2g)
     &          *(1.0-1.0/sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3))))        
     &                  -TzeroB_fsh(i,3,jq,2)  )

              SigB_d(3,jq,1)=Mem_Coef_fsh(i,jq,j,2)*e4coef(i2g)
     &              /sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3)))**3
              SigB_d(3,jq,2)=Mem_Coef_fsh(i,jq,j,1)
     &              /sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3)))**3

              elseif(e4coef(i2g) .le. 100.0) then   
                     temp_ibm3= e4coef(i2g)-50.0          
                     sigB(3,jq,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3)))-1.0)        
     &            /(sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3)))) )
     &                  -TzeroB_fsh(i,3,jq,2)  )

                 SigB_d(3,jq,1)=Mem_Coef_fsh(i,jq,j,2)
     &              /sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3)))**3
     &              +Mem_Coef_fsh(i,jq,j,2)*temp_ibm3
     &             *(sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3)))-1.0)        
     &             /(-sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))**3
                 SigB_d(3,jq,2)=Mem_Coef_fsh(i,jq,j,2)*temp_ibm3
     &              /sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3)))
     &            /(sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))
              endif
! boundary 4
           j=ns_ibm_r_fsh(i)
              DY_i_jm(1:3)=
     &         (    GXtmp(1:3,jq,j-2)
     &         -4.0*GXtmp(1:3,jq,j-1)
     &         +3.0*GXtmp(1:3,jq,j)   )
     &         /(3.0*dsf2_IBM_fsh(i,jq,j-1)-dsf2_IBM_fsh(i,jq,j-2))     !x+1/2,y

              DX_ip_j(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))/dsf_IBM_fsh(i,jq,j) !x,y+1/2

              DX_im_j(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))
     &                     /dsf_IBM_fsh(i,jq-1,j)                      !x,y-1/2

              DX_i_jp(1:3)=
     &        0.5*(DX_ip_j(1:3)+DX_im_j(1:3))


              sigB(4,jq,1)=0.0
              sigB(4,jq,3)=Mem_Coef_fsh(i,jq,j,3)
     &            *(dot_product(DX_i_jp(1:3),DY_i_jm(1:3))
     &                  -TzeroB_fsh(i,4,jq,3)  )
              SigB_d(4,jq,3)=Mem_Coef_fsh(i,jq,j,3)
              if(e4coef(i2g) .lt. 50.0) then 
              sigB(4,jq,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DY_i_jm(1:3),DY_i_jm(1:3))))
     &        +e4coef(i2g)
     &          *(1.0-1.0/sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3))))        
     &                  -TzeroB_fsh(i,4,jq,2)  )

              SigB_d(4,jq,1)=Mem_Coef_fsh(i,jq,j,2)*e4coef(i2g)
     &              /sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3)))**3
              SigB_d(4,jq,2)=Mem_Coef_fsh(i,jq,j,1)
     &              /sqrt(dot_product(DY_i_jm(1:3),DY_i_jm(1:3)))**3

              elseif(e4coef(i2g) .le. 100.0) then   
                     temp_ibm3= e4coef(i2g)-50.0          
                     sigB(4,jq,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DY_i_jm(1:3),DY_i_jm(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3)))-1.0)        
     &            /(sqrt(dot_product(DY_i_jm(1:3),DY_i_jm(1:3)))) )
     &                  -TzeroB_fsh(i,4,jq,2)  )

                 SigB_d(4,jq,1)=Mem_Coef_fsh(i,jq,j,2)
     &              /sqrt(dot_product(DY_i_jm(1:3),DY_i_jm(1:3)))**3
     &              +Mem_Coef_fsh(i,jq,j,2)*temp_ibm3
     &             *(sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3)))-1.0)        
     &             /(-sqrt(dot_product(DY_i_jm(1:3),DY_i_jm(1:3))))**3
                 SigB_d(4,jq,2)=Mem_Coef_fsh(i,jq,j,2)*temp_ibm3
     &              /sqrt(dot_product(DX_i_jp(1:3),DX_i_jp(1:3)))
     &            /(sqrt(dot_product(DY_i_jm(1:3),DY_i_jm(1:3))))
              endif

           enddo

           endif
