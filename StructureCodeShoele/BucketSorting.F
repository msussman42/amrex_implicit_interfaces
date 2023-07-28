c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine BucketSorting
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c Bucket calculation for contact model
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       USE HeaderFSI
       implicit none
       integer i2,jq,j,IBucketx,IBuckety,IBucketz,ibkt,k,i2l
     & ,ne2,jtype2,itype2,ntype2,ItypecontELE(3),ntem2,nos2,noj
       real*8 Xcont(3)
          if(nr_ibm_fib.ge. 1) then 
            IBucket_fib=0
            SBucket_fib=0
          endif
          if(nr_ibm_fsh.ge. 1) then 
            IBucket_fsh=0
            SBucket_fsh=0
          endif
          if(nr_ibm_esh.ge. 1) then 
            IBucket_esh=0
            SBucket_esh=0
          endif
          if(nr_ibm_fbc.ge. 1) then 
            IBucket_fbc=0
            SBucket_fbc=0
          endif

          do i2=1,nr_ibm
            NBucket(i2,1:NBucketnum)=0
           if (bodytype(i2) .eq. shell)then
                i2l=Iglbloc_esh(i2)
                if(contacttype .gt. 1) then
                do ne2=1,numelr(i2)
                     jtype2=neatype_esh(i2l,ne2)
                     itype2=jtypeall(jtype2)
                     if(itype2 .eq. 0) then
                        ntype2=12  
                     else
                        ntype2=6+itype2
                     end if   
                     if(itype2 .eq. 0) then
                       ItypecontELE(1)=4
                       ItypecontELE(2)=7
                       ItypecontELE(3)=8
                     else              
                       ItypecontELE(1)=1
                       ItypecontELE(2)=2
                       ItypecontELE(3)=itype2+1
                     endif
                     Xcont(1:3)=0.d0
                     do nos2=1,3
                        ntem2=nea_esh(i2l,ne2,ItypecontELE(nos2))
                        noj=1
                        Xcont(noj)=Xcont(noj)+gx_ibm_esh(i2l,ntem2)/3.d0
                        noj=2
                        xcont(noj)=Xcont(noj)+gy_ibm_esh(i2l,ntem2)/3.d0
                        noj=3
                        xcont(noj)=Xcont(noj)+gz_ibm_esh(i2l,ntem2)/3.d0
                     enddo

                    IBucketx=
     &        int((Xcont(1)-BucketXStart(1))/BucketdX(1))+1
             
                    IBuckety=
     &        int((Xcont(2)-BucketXStart(2))/BucketdX(2))+1

                    IBucketz=
     &        int((Xcont(3)-BucketXStart(3))/BucketdX(3))+1

                    IBucket_esh(i2l,ne2)=
     &                       IBucketx
     &                     +(IBuckety-1)*IBucketnum(1)
     &                     +(IBucketz-1)*IBucketnum(1)*IBucketnum(2)

                    NBucket(i2,IBucket_esh(i2l,ne2))=
     &              NBucket(i2,IBucket_esh(i2l,ne2))+1

                enddo
                else
                do ne2=1,ns_ibm_r_esh(i2l)
                        noj=1
                        Xcont(noj)=gx_ibm_esh(i2l,ne2)
                        noj=2
                        Xcont(noj)=gy_ibm_esh(i2l,ne2)
                        noj=3
                        Xcont(noj)=gz_ibm_esh(i2l,ne2)


                    IBucketx=
     &        int((Xcont(1)-BucketXStart(1))/BucketdX(1))+1
             
                    IBuckety=
     &        int((Xcont(2)-BucketXStart(2))/BucketdX(2))+1

                    IBucketz=
     &        int((Xcont(3)-BucketXStart(3))/BucketdX(3))+1

                    IBucket_esh(i2l,ne2)=
     &                       IBucketx
     &                     +(IBuckety-1)*IBucketnum(1)
     &                     +(IBucketz-1)*IBucketnum(1)*IBucketnum(2)

                    NBucket(i2,IBucket_esh(i2l,ne2))=
     &              NBucket(i2,IBucket_esh(i2l,ne2))+1

                enddo
              endif
           elseif (bodytype(i2) .eq. fabric)then
                i2l=Iglbloc_fbc(i2)
                if(contacttype .gt. 1) then
                do ne2=1,numelr(i2)
                     Xcont(1:3)=0.d0
                     do nos2=1,3
                        ntem2=nea_fbc(i2l,ne2,ItypecontELE(nos2))
                        noj=1
                        Xcont(noj)=Xcont(noj)+gx_ibm_fbc(i2l,ntem2)/3.d0
                        noj=2
                        xcont(noj)=Xcont(noj)+gy_ibm_fbc(i2l,ntem2)/3.d0
                        noj=3
                        xcont(noj)=Xcont(noj)+gz_ibm_fbc(i2l,ntem2)/3.d0
                     enddo

                    IBucketx=
     &        int((Xcont(1)-BucketXStart(1))/BucketdX(1))+1
             
                    IBuckety=
     &        int((Xcont(2)-BucketXStart(2))/BucketdX(2))+1

                    IBucketz=
     &        int((Xcont(3)-BucketXStart(3))/BucketdX(3))+1

                    IBucket_fbc(i2l,ne2)=
     &                       IBucketx
     &                     +(IBuckety-1)*IBucketnum(1)
     &                     +(IBucketz-1)*IBucketnum(1)*IBucketnum(2)

                    NBucket(i2,IBucket_fbc(i2l,ne2))=
     &              NBucket(i2,IBucket_fbc(i2l,ne2))+1

                enddo
                else
                do ne2=1,ns_ibm_r_fbc(i2l)
                        noj=1
                        Xcont(noj)=gx_ibm_fbc(i2l,ne2)
                        noj=2
                        Xcont(noj)=gy_ibm_fbc(i2l,ne2)
                        noj=3
                        Xcont(noj)=gz_ibm_fbc(i2l,ne2)


                    IBucketx=
     &        int((Xcont(1)-BucketXStart(1))/BucketdX(1))+1
             
                    IBuckety=
     &        int((Xcont(2)-BucketXStart(2))/BucketdX(2))+1

                    IBucketz=
     &        int((Xcont(3)-BucketXStart(3))/BucketdX(3))+1

                    IBucket_fbc(i2l,ne2)=
     &                       IBucketx
     &                     +(IBuckety-1)*IBucketnum(1)
     &                     +(IBucketz-1)*IBucketnum(1)*IBucketnum(2)

                    NBucket(i2,IBucket_fbc(i2l,ne2))=
     &              NBucket(i2,IBucket_fbc(i2l,ne2))+1

                enddo
              endif
            elseif (bodytype(i2) .eq. fiber)then
                i2l=Iglbloc_fib(i2)              
                 do j=1,ns_ibm_r_fib(i2l)
                    IBucketx=
     &        int((gx_ibm_fib(i2l,j)-BucketXStart(1))/BucketdX(1))+1
             
                    IBuckety=
     &        int((gy_ibm_fib(i2l,j)-BucketXStart(2))/BucketdX(2))+1

                    IBucketz=
     &        int((gz_ibm_fib(i2l,j)-BucketXStart(3))/BucketdX(3))+1
               
                    IBucket_fib(i2l,j)=
     &                       IBucketx
     &                     +(IBuckety-1)*IBucketnum(1)
     &                     +(IBucketz-1)*IBucketnum(1)*IBucketnum(2)
                    NBucket(i2,IBucket_fib(i2l,j))=
     &              NBucket(i2,IBucket_fib(i2,j))+1
                 enddo
            elseif (bodytype(i2) .eq. fibrousshell)then
                i2l=Iglbloc_fsh(i2)              
                 do jq=1,nq_ibm_r_fsh(i2l)               
                 do j=1,ns_ibm_r_fsh(i2l)
                    IBucketx=
     &       int((gx_ibm_fsh(i2l,jq,j)-BucketXStart(1))/BucketdX(1))+1
             
                    IBuckety=
     &       int((gy_ibm_fsh(i2l,jq,j)-BucketXStart(2))/BucketdX(2))+1

                    IBucketz=
     &       int((gz_ibm_fsh(i2l,jq,j)-BucketXStart(3))/BucketdX(3))+1
               
                    IBucket_fsh(i2l,jq,j)=
     &                       IBucketx
     &                     +(IBuckety-1)*IBucketnum(1)
     &                     +(IBucketz-1)*IBucketnum(1)*IBucketnum(2)
                    NBucket(i2,IBucket_fsh(i2l,jq,j))=
     &              NBucket(i2,IBucket_fsh(i2l,jq,j))+1
                 enddo
                 enddo
            endif

            PBucket(i2,1)=1
            do ibkt=2,NBucketnum+1
               PBucket(i2,ibkt)=PBucket(i2,ibkt-1)+NBucket(i2,ibkt-1)
            enddo
            NBucket(i2,1:NBucketnum)=0
            if (bodytype(i2) .eq. shell)then
                i2l=Iglbloc_esh(i2)
                if(contacttype .gt. 1) then
                do ne2=1,numelr(i2)
                    k=IBucket_esh(i2l,ne2)
                    NBucket(i2,k)= NBucket(i2,k)+1
                    sBucket_esh(i2l,k,NBucket(i2,k)+PBucket(i2,k)-1)=ne2
                enddo
                else
                do ne2=1,ns_ibm_r_esh(i2l)
                    k=IBucket_esh(i2l,ne2)
                    NBucket(i2,k)= NBucket(i2,k)+1
                    sBucket_esh(i2l,k,NBucket(i2,k)+PBucket(i2,k)-1)=ne2
!                    write(1818,200) ne2,ns_ibm_r_esh(i2l),k
!     &             ,NBucket(i2,k),PBucket(i2,k)
!     &             ,sBucket_esh(i2l,k,NBucket(i2,k)+PBucket(i2,k)-1)

                enddo
                endif
             elseif (bodytype(i2) .eq. fiber)then
                i2l=Iglbloc_fib(i2)
                 do j=1,ns_ibm_r_fib(i2l)
                    k=IBucket_fib(i2l,j)
                    NBucket(i2,k)= NBucket(i2,k)+1
                    sBucket_fib(i2l,k,NBucket(i2,k)+PBucket(i2,k)-1)=j

                enddo             
             elseif (bodytype(i2) .eq. fibrousshell)then
                i2l=Iglbloc_fsh(i2)
                 do jq=1,nq_ibm_r_fsh(i2l)             
                 do j=1,ns_ibm_r_fsh(i2l)
                    k=IBucket_fsh(i2l,jq,j)
                  NBucket(i2,k)= NBucket(i2,k)+1
                  sBucket_fsh(i2l,k,NBucket(i2,k)+PBucket(i2,k)-1,1)=j
                  sBucket_fsh(i2l,k,NBucket(i2,k)+PBucket(i2,k)-1,2)=jq

                 enddo
                 enddo
              endif
         enddo
200   format(20(i5,' '))
         return
         end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
