c  August 31, 2023 grid_def.F changed
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine shapeinitialize(jtypeall
     *   ,nn,dnndv,dnndw,ddnndvdv,ddnndwdw,ddnndvdw
     *   ,nnCont,dnndvCont,dnndwCont,
     *    ddnndvdvCont,ddnndwdwCont,ddnndvdwCont
     *   ,nn0,dnndv0,dnndw0,ddnndvdv0,ddnndwdw0,ddnndvdw0
     *   ,gau,ngau,wei
     *   ,gauCont,ngauCont,weiCont
     *   ,gau0   ,ngau0   ,wei0
     *   ,gau_h  ,ngau_h  ,wei_h
     *   ,n_jtypeall)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c initialize shape functions for FE shell calculations
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       include 'grid_def'
       integer i, j, jtype, itype, n_jtypeall, jtypeall(jtypemax)
       real v, w, l(jtypemax,maxordern)
     *       , iv(jtypemax,maxordern,maxordern) 
     *       , phi(jtypemax,0:2,maxordern,12)
       integer  ngau,ngauCont,ngau0,ngau_h
       real  wei(ngaumax), wei0(ngaumax), weiCont(ngaumax)
       real  gau(2,ngaumax), gau0(2,ngaumax), gauCont(2,ngaumax)
       real  gau_h(ngaumaxh), wei_h(ngaumaxh)
       real nn(jtypemax,ngaumax,maxordern), s(maxordern)
       real dnndv(jtypemax,ngaumax,maxordern), dsdv(maxordern) 
       real dnndw(jtypemax,ngaumax,maxordern), dsdw (maxordern)
       real ddnndvdv(jtypemax,ngaumax,maxordern), ddsdvdv(maxordern) 
       real ddnndwdw(jtypemax,ngaumax,maxordern), ddsdwdw(maxordern)
       real ddnndvdw(jtypemax,ngaumax,maxordern), ddsdvdw(maxordern)

       real nnCont(jtypemax,ngaumax,maxordern)
       real dnndvCont(jtypemax,ngaumax,maxordern)
       real dnndwCont(jtypemax,ngaumax,maxordern)
       real ddnndvdvCont(jtypemax,ngaumax,maxordern) 
       real ddnndwdwCont(jtypemax,ngaumax,maxordern)
       real ddnndvdwCont(jtypemax,ngaumax,maxordern)

       real nn0(jtypemax,1,maxordern)
       real dnndv0(jtypemax,1,maxordern) 
       real dnndw0(jtypemax,1,maxordern)
       real ddnndvdv0(jtypemax,1,maxordern)
       real ddnndwdw0(jtypemax,1,maxordern)
       real ddnndvdw0(jtypemax,1,maxordern)


       call quad2d3n(ngau,gau,wei,0)
       call quad2d3n(ngauCont,gauCont,weiCont,1)
       call quad2d3n(ngau0,gau0,wei0,-1)
       call quad1d(gau_h,wei_h,ngau_h)

        do i=1,n_jtypeall
           jtype=i
           itype=jtypeall(i)
           open(13,file = 'shapedata.in',  form = 'formatted')
           call matread(jtype,itype,l,iv,phi)
           close(13)
           do j=1,ngau
              v=gau(1,j)
              w=gau(2,j)
              call eval_n_dn_ddn(v,w, 
     *                             s,dsdv,dsdw,ddsdvdv,ddsdwdw,ddsdvdw, 
     *                             itype, phi, l, iv, jtype,
     *                             jtypemax,maxordern)

              nn(i,j,:)=s 
              dnndv(i,j,:)=dsdv 
              dnndw(i,j,:)=dsdw 
              ddnndvdv(i,j,:)=ddsdvdv
              ddnndwdw(i,j,:)=ddsdwdw
              ddnndvdw(i,j,:)=ddsdvdw
           end do

           do j=1,ngauCont
              v=gauCont(1,j)
              w=gauCont(2,j)
              call eval_n_dn_ddn(v,w, 
     *                             s,dsdv,dsdw,ddsdvdv,ddsdwdw,ddsdvdw, 
     *                             itype, phi, l, iv, jtype,
     *                             jtypemax,maxordern)

              nnCont(i,j,:)=s 
              dnndvCont(i,j,:)=dsdv 
              dnndwCont(i,j,:)=dsdw 
              ddnndvdvCont(i,j,:)=ddsdvdv
              ddnndwdwCont(i,j,:)=ddsdwdw
              ddnndvdwCont(i,j,:)=ddsdvdw
           end do

           do j=1,ngau0

              v=gau0(1,j)
              w=gau0(2,j)


              call eval_n_dn_ddn(v,w, 
     *                             s,dsdv,dsdw,ddsdvdv,ddsdwdw,ddsdvdw, 
     *                             itype, phi, l, iv, jtype,
     *                             jtypemax,maxordern)

              nn0(i,j,:)=s 
              dnndv0(i,j,:)=dsdv 
              dnndw0(i,j,:)=dsdw 
              ddnndvdv0(i,j,:)=ddsdvdv
              ddnndwdw0(i,j,:)=ddsdwdw
              ddnndvdw0(i,j,:)=ddsdvdw
           end do
        end do   

        return
        end   
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
