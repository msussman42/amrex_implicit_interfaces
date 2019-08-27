module SOLIDCouplerIO

use CouplerIO, only : CouplerIO_WriteHeader,CouplerIO_WriteNodes, &
                      CouplerIO_ReadForces

implicit none

logical(4) :: DebugCouple
integer(4) :: FlagDim
integer(4) :: IntElemDim,IntElemDimPaddle,IntElemDimPool 
integer(4) :: UnitCouple 
integer(4) :: UnitOut
integer(4) :: UnitDebug
integer(4) :: StdOut
integer(4) :: NumNodes,NumNodesPaddle,NumNodesPool,use_temp
integer(4) :: NumIntElems,NumIntElemsPaddle,NumIntElemsPool
integer(4), dimension(4) :: NumFlags
integer(4), dimension(:), pointer :: HdrFlag
integer(4), dimension(:), pointer :: Eul2IntNode
integer(4), dimension(:,:), pointer :: ElemData
integer(4), dimension(:,:), pointer :: IntElem
integer(4) :: istepB,sci_sdim,sci_istop,sci_istep,sci_probtype
real(8) :: sci_curtime,sci_dt
real(8), dimension(3) :: sci_problo,sci_probhi
real(8) :: timeB,tstart,tfinish
real(8) :: dtB
integer(4), dimension(:), pointer :: NptFlag
integer(4), dimension(:), pointer :: ForFlag
real(8), dimension(:,:), allocatable :: Node_old,Node_new,Node_current
real(8), dimension(:,:), allocatable :: NodeVel_old,NodeVel_new
real(8), dimension(:), allocatable :: NodeTemp_old,NodeTemp_new
real(8), dimension(:,:), allocatable :: NodeForce
logical(4), dimension(:), allocatable :: ActiveIntElem
real(8) :: radblob,denpaddle,dampingpaddle,TorquePos,TorqueVel,radblobwall
real(8) :: raddust,dendust,adheredust,floordust,tempdust
real(8), dimension(3) :: torquePosDust,torqueVelDust, &
 centerDust,centerVelDust,centerStartDust
real(8), dimension(3) :: xblob,newxblob,xblobwall,newxblobwall

contains



subroutine scihandoffset(ofs,time)
IMPLICIT NONE

real(8) :: ofs,time

 ofs=0.0
 if (time.le.0.4) then
  ofs=-time/8.0
 else if ((time.ge.0.4).and.(time.lt.0.8)) then
  ofs=-0.05
 else if ((time.ge.0.8).and.(time.lt.1.2)) then
  ofs=(time-0.9)/2.0
 else
  ofs=0.15
 endif

return
end subroutine scihandoffset

subroutine geominit(curtime,dt,ifirst,sdim,istop,istep,probtype)
integer(4) :: i,j,k,i1,j1,k1,itimecount,ifirst,sdim
real(8) :: curtime,dt
integer(4), dimension(0:1000) :: itimearr
real(8), dimension(0:1000) :: timearr
real(8), dimension(3) :: maxnode,minnode,xval,xtemp
real(8), dimension(3) :: maxnodebefore,minnodebefore
real(8), dimension(3) :: xvalbefore
real(8), dimension(3) :: invertfactor
real(8) :: tper,tcrit,theta,radblobpool
integer(4) :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
integer(4) :: gwave
integer(4) :: probtype
integer(4) :: shift_from_zero_node
integer(4) :: shift_from_zero_face
integer(4) :: override_IntElemDim
character(20) :: dwave,poolname
real(8) :: ofs
real(8) :: plungerfreq,plungeramp

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=dt

  if (ifirst.eq.1) then
   allocate(HdrFlag(NumFlags(1)))
   allocate(NptFlag(1))
   allocate(ForFlag(NumFlags(4))) 
  endif

  if (probtype.eq.58) then
   iread=1
  else if (probtype.eq.55) then
   if (ifirst.eq.1) then
    iread=1
   else
    iread=0
   endif
  else
 
   OPEN(unit=25,file="hcorresp.txt",access='sequential', &
     form="formatted",status='old')
   READ(25,*) itimecount
   print *,"itimecount= ",itimecount

   if (itimecount.gt.1000) then
    print *,"itimecoount too big"
    stop
   endif

   do i1=1,itimecount
    READ(25,*) itimearr(i1),timearr(i1)
    print *,"index,time ",itimearr(i1),timearr(i1)
   enddo
   close(25)

   timearr(0)=0.0
   itimearr(0)=itimearr(itimecount)

   tper=curtime/timearr(itimecount)
   iper=INT(tper)
   tcrit=curtime-iper*timearr(itimecount)
   if (tcrit.gt.timearr(itimecount)+1.0E-8) then
    print *,"tcrit invalid"
    stop
   endif

   icrit=0
   do i1=1,itimecount
    if ((tcrit.le.timearr(i1)+1.0E-8).and.(icrit.eq.0)) then
     icrit=i1
    endif
   enddo
   if (icrit.eq.0) then
    print *,"icrit invalid"
    stop
   endif

   if (ifirst.eq.1) then
    tstart=-1.0
    tfinish=-1.0
   endif

   if ((tcrit.ge.tstart).and.(tcrit.le.tfinish)) then
    iread=0
   else
    iread=1
    tstart=timearr(icrit-1)
    tfinish=timearr(icrit)
   endif
  endif

  if (iread.eq.1) then

   NumNodesPool=0
   NumIntElemsPool=0
   IntElemDimPool=0

   do i1=1,2 
    if ((probtype.eq.55).or.(probtype.eq.58)) then
     j=0
    else
     if (i1.eq.1) then
      j=itimearr(icrit-1)
     else
      j=itimearr(icrit)
     endif
    endif

    invertfactor(1)=1.0
    invertfactor(2)=1.0
    invertfactor(3)=1.0

    if (probtype.eq.55) then
     xblob(1)=0.0
     xblob(2)=0.0
     xblob(3)=0.0
     newxblob(1)=0.0
     newxblob(2)=0.0
     newxblob(3)=0.055
     radblob=600.0
     dwave="Sphere.txt"
    else if (probtype.eq.58) then
     plungerfreq=14.0
     plungeramp=0.015  ! 0.015
     xblob(1)=0.0
     xblob(2)=0.0
     xblob(3)=0.0
     newxblob(1)=0.0
     newxblob(2)=0.0
     newxblob(3)=0.0
     radblob=1.0
     newxblob(1)=plungeramp*sin(2.0*3.14159*plungerfreq*curtime)
     dwave="square.txt"
    else if (probtype.eq.52) then
     xblob(1)=260.0
     xblob(2)=19.0
     xblob(3)=-191.0
     newxblob(3)=0.5
     newxblob(1)=0.25
     newxblob(2)=0.25
     radblob=440.0

     if (j.lt.10) then
      dwave = "h/hand000"//char(48+j)//".txt"
     else if (j.lt.100) then
       ewave = j/10
       dwave = "h/hand00"//char(48+ewave)//char(48+j-ewave*10)//".txt"
     else
       ewave = j/100
       fwave = (j-ewave*100)/10
       dwave = "h/hand0"//char(48+ewave)//char(48+fwave)// &
              char(48+j-fwave*10-ewave*100)//".txt"
     endif
    else if (probtype.eq.57) then
!       maya geberated heart
     xblob(1)=-375.0
     xblob(2)=890.0
     xblob(3)=-70.0
!       vue generated heart
!     xblob(1)=-15.0
!     xblob(2)=-50.0
!     xblob(3)=0.0

     newxblob(1)=0.5
     newxblob(2)=0.5
     newxblob(3)=0.5
     radblob=550.0

     if (j.lt.10) then
      dwave = "h/heartleft000"//char(48+j)//".txt"
     else if (j.lt.100) then
       ewave = j/10
       dwave = "h/heartleft00"//char(48+ewave)//char(48+j-ewave*10)//".txt"
     else
       ewave = j/100
       fwave = (j-ewave*100)/10
       dwave = "h/heartleft0"//char(48+ewave)//char(48+fwave)// &
              char(48+j-fwave*10-ewave*100)//".txt"
     endif
    else if (probtype.eq.562) then
     xblob(1)=0.0
     xblob(2)=0.0
     xblob(3)=0.0
     newxblob(1)=0.5
     newxblob(2)=0.5
     newxblob(3)=1.0
! was 10 for whale file dated December, 2007
     radblob=12.5  
! whaletailup or whaletaildown
     xblob(1)=0.0 ! will become y
     xblob(2)=-1.0 ! will become -x
     invertfactor(2)=-1.0
     xblob(3)=0.0  ! will become z
     newxblob(1)=0.5
     newxblob(2)=1.0
     newxblob(3)=0.5
     radblob=12.5
! whalepregnant
     if (1.eq.0) then
     xblob(1)=0.0 ! will become y
     xblob(3)=-1.0 ! will become -x
     invertfactor(3)=-1.0
     xblob(2)=0.0  ! will become z
     newxblob(1)=0.5
     newxblob(3)=1.0
     newxblob(2)=0.5
     radblob=12.5
     endif
 
     if (j.lt.10) then
      dwave = "h/whale000"//char(48+j)//".txt"
     else if (j.lt.100) then
       ewave = j/10
       dwave = "h/whale00"//char(48+ewave)//char(48+j-ewave*10)//".txt"
     else
       ewave = j/100
       fwave = (j-ewave*100)/10
       dwave = "h/whale0"//char(48+ewave)//char(48+fwave)// &
              char(48+j-fwave*10-ewave*100)//".txt"
     endif
    else if ((probtype.eq.56).or.(probtype.eq.561)) then
     xblob(1)=0.2
     xblob(2)=0.25
     xblob(3)=0.5
     if (1.eq.0) then
      newxblob(3)=0.65
      newxblob(1)=0.375
      newxblob(2)=0.375
      radblob=40.0
      radblobpool=40.0
     else
      newxblob(3)=26
      newxblob(1)=15   
      newxblob(2)=15
      radblob=1.0
      radblobpool=1.0
     endif


! Viorel's test problem with a box.
     if (1.eq.0) then   
      xblob(1)=0.0
      xblob(2)=0.0
      xblob(3)=0.0
      radblob=3.0
     endif

! Viorel's ball and scythe problems 
     if (1.eq.0) then   
      xblob(1)=2.5
      xblob(2)=86.0
      xblob(3)=-1.0
      newxblob(3)=30.0
      newxblob(1)=30.0
      newxblob(2)=90.0
      radblob=1.0
     endif
!Viorel's diving problem     
     if (1.eq.1) then   
      xblob(1)=0.0
      xblob(2)=36.0
      xblob(3)=-0.5
      newxblob(3)=30.0
      newxblob(1)=30.0
      newxblob(2)=18.0
      radblob=1.0
     endif
!Viorel's paddle problem     
     if (1.eq.1) then   
      xblob(1)=0.0
      xblob(2)=0.0
      xblob(3)=0.0
      newxblob(3)=30.0
      newxblob(1)=60.0
      newxblob(2)=40.0
      radblob=1.0
     endif

     if (probtype.eq.561) then
      poolname="h/Pool.txt"
      print *,"opening ",poolname
      OPEN(unit=141,file=poolname,access='sequential', &
        form="formatted",status='old')
      READ(141,*) NumNodesPool
      print *,"NumNodesPool ",NumNodesPool
      READ(141,*) NumIntElemsPool
      print *,"NumIntElemsPool ",NumIntElemsPool
      READ(141,*) IntElemDimPool
      print *,"IntElemDimPool ",IntElemDimPool
     else
      NumNodesPool=0
      NumIntElemsPool=0
      IntElemDimPool=0
     endif
 
     if (j.lt.10) then
      dwave = "h/swimmer000"//char(48+j)//".txt"
     else if (j.lt.100) then
       ewave = j/10
       dwave = "h/swimmer00"//char(48+ewave)//char(48+j-ewave*10)//".txt"
     else if (j.lt.1000) then
       ewave = j/100
       fwave = (j-ewave*100)/10
       dwave = "h/swimmer0"//char(48+ewave)//char(48+fwave)// &
              char(48+j-fwave*10-ewave*100)//".txt"
     else
       ewave = j/1000
       fwave = (j-ewave*1000)/100
       gwave = (j-fwave*100-ewave*1000)/10
       dwave = "h/swimmer"//char(48+ewave)//char(48+fwave)// &
        char(48+gwave)//char(48+j-gwave*10-fwave*100-ewave*1000)//".txt"
     endif
    else
     print *,"probtype invalid"
     stop
    endif

    print *,"opening ",dwave
    OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

    READ(14,*) NumNodes
    print *,"NumNodes ",NumNodes
    READ(14,*) NumIntElems
    print *,"NumIntElems ",NumIntElems
    READ(14,*) IntElemDim
    print *,"IntElemDim ",IntElemDim
    NumNodes=NumNodes+NumNodesPool
    NumIntElems=NumIntElems+NumIntElemsPool
    if (IntElemDimPool.ne.0) then
     if (IntElemDim.ne.IntElemDimPool) then
      print *,"IntElemDim inconsistent"
      stop
     endif
    endif

    if ((ifirst.eq.1).and.(i1.eq.1)) then
     allocate(Eul2IntNode(NumNodes))
     allocate(ElemData(NumFlags(2),NumIntElems))
     allocate(IntElem(IntElemDim,NumIntElems))
     allocate(ActiveIntElem(NumIntElems))
     allocate(Node_old(3,NumNodes))
     allocate(Node_new(3,NumNodes))
     allocate(Node_current(3,NumNodes))
     allocate(NodeVel_old(3,NumNodes))
     allocate(NodeVel_new(3,NumNodes))
     allocate(NodeTemp_old(NumNodes))
     allocate(NodeTemp_new(NumNodes))
     allocate(NodeForce(3,NumNodes))

     do i=1,NumNodes
      Eul2IntNode(i)=i
      NodeTemp_old(i)=0.0
      NodeTemp_new(i)=0.0
      do dir=1,3
       Node_old(dir,i)=0.0
       Node_new(dir,i)=0.0
       Node_current(dir,i)=0.0
       NodeVel_old(dir,i)=0.0
       NodeVel_new(dir,i)=0.0
       NodeForce(dir,i)=0.0
      enddo
     enddo
     do i=1,NumIntElems
      ActiveIntElem(i)=.true.
     enddo
    endif

    HdrFlag(1)=sdim
    HdrFlag(2)=0     ! stationary body

    do dir=1,3
     maxnode(dir)=-1.0e+10
     minnode(dir)=1.0e+10
     maxnodebefore(dir)=-1.0e+10
     minnodebefore(dir)=1.0e+10
    enddo

    if (probtype.eq.561) then

     do i=1,NumNodesPool
      READ(141,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)

      do dir=1,3
       xval(dir)=invertfactor(dir)* &
        (xvalbefore(dir)-xblob(dir))/radblobpool + newxblob(dir)
      enddo

      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
    
      if (i1.eq.1) then
       do dir=1,3
        Node_old(dir,i)=xval(dir)
       enddo
      else if (i1.eq.2) then
       do dir=1,3
        Node_new(dir,i)=xval(dir)
       enddo
      else
       print *,"i1 invalid"
       stop
      endif
      if (i.ne.j) then
       print *,"vertex mismatch reading pool file"
       stop
      endif
      
     enddo
     READ(141,*) j
     if (j.ne.NumIntElemsPool) then
      print *,"face mismatch reading pool file"
      stop
     endif

     do i=1,NumIntElemsPool
      READ(141,*) j
      if (i.ne.j) then
       print *,"face mismatch reading pool file"
       stop
      endif
      READ(141,*) k
      if (k.gt.IntElemDim) then
       print *,"too many vertices k,IntElemDim ",k,IntElemDim
       stop
      endif
      ElemData(1,i)=k   ! number of nodes in element
      ElemData(2,i)=1   ! part number
      ElemData(3,i)=2   ! singly wetted, but do not call "fill" for these...
      do j1=1,k
       READ(141,*) IntElem(j1,i) 
      enddo
     enddo  ! i, looping faces
     close(141)
    endif ! probtype.eq.561 


    shift_from_zero_node=0
    shift_from_zero_face=0
    override_IntElemDim=0
    do i=NumNodesPool+1,NumNodes
     READ(14,*) j,xvalbefore(1),xvalbefore(2),xvalbefore(3)
     if ((j.eq.0).and.(shift_from_zero_node.eq.0)) then
      shift_from_zero_node=1
      print *,"nodes in file start at 0; shifting to start at 1"
     endif
     if (shift_from_zero_node.eq.1) then
      j=j+1
     endif

     do dir=1,3
      xval(dir)=invertfactor(dir)* &
       (xvalbefore(dir)-xblob(dir))/radblob + newxblob(dir)
     enddo

     if (probtype.eq.55) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(1)
      xval(2)=xtemp(2)
      xval(3)=xtemp(3)
     else if (probtype.eq.58) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(1)
      xval(2)=xtemp(2)
      xval(3)=xtemp(3)
     else if (probtype.eq.52) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
      if (i1.eq.1) then
       call scihandoffset(ofs,tstart)
      else if (i1.eq.2) then
       call scihandoffset(ofs,tfinish)
      else
       print *,"i1 invalid"
       stop
      endif
      xval(sdim)=xval(sdim)+ofs
     else if (probtype.eq.562) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
! whale file used at beginning of summer 2007
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
! whaletailup or whaletaildown
      xval(1)=xtemp(2)
      xval(2)=xtemp(1)
      xval(3)=xtemp(3)
! whale pregnant
      if (1.eq.0) then
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
      endif
     else if ((probtype.eq.56).or.(probtype.eq.561)) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
     else if (probtype.eq.57) then
      do dir=1,3
       xtemp(dir)=xval(dir)
      enddo
      xval(1)=xtemp(3)
      xval(2)=xtemp(1)
      xval(3)=xtemp(2)
     else
      print *,"probtype invalid"
      stop
     endif
     do dir=1,3
      if ((minnode(dir).gt.xval(dir)).or.(i.eq.1)) then
       minnode(dir)=xval(dir)
      endif
      if ((maxnode(dir).lt.xval(dir)).or.(i.eq.1)) then
       maxnode(dir)=xval(dir)
      endif
      if ((minnodebefore(dir).gt.xvalbefore(dir)).or.(i.eq.1)) then
       minnodebefore(dir)=xvalbefore(dir)
      endif
      if ((maxnodebefore(dir).lt.xvalbefore(dir)).or.(i.eq.1)) then
       maxnodebefore(dir)=xvalbefore(dir)
      endif
     enddo
    
     if (i1.eq.1) then
      do dir=1,3
       Node_old(dir,i)=xval(dir)
      enddo
     else if (i1.eq.2) then
      do dir=1,3
       Node_new(dir,i)=xval(dir)
      enddo
     else
      print *,"i1 invalid"
      stop
     endif
     if (i-NumNodesPool.ne.j) then
      print *,"vertex mismatch"
      print *,"NumNodesPool=",NumNodesPool
      print *,"expected node index ",i-NumNodesPool
      print *,"node index read from file ",j
      stop
     endif
      
    enddo
    READ(14,*) j
    if (j.ne.NumIntElems-NumIntElemsPool) then
     print *,"face mismatch"
     stop
    endif

    do i=1+NumIntElemsPool,NumIntElems
     READ(14,*) j
     if ((j.eq.0).and.(shift_from_zero_face.eq.0)) then
      shift_from_zero_face=1
      print *,"faces in file start at 0; shifting to start at 1"
     endif
     if (shift_from_zero_face.eq.1) then
      j=j+1
     endif

     if (i-NumIntElemsPool.ne.j) then
      print *,"face mismatch"
      print *,"NumIntElemsPool=",NumIntElemsPool
      print *,"expected face index ",i-NumIntElemsPool
      print *,"face index read from file ",j
      stop
     endif
     READ(14,*) k
     if (k.gt.IntElemDim) then
      if (override_IntElemDim.eq.0) then
       override_intElemDim=1 
       print *,"too many vertices k,IntElemDim ",k,IntElemDim
       print *,"will override further errors to be ",IntElemDim
      endif
      k=IntElemDim
     endif
     ElemData(1,i)=k   ! number of nodes in element
     ElemData(2,i)=1   ! part number
     ElemData(3,i)=0   ! singly wetted
     if (probtype.eq.57) then  ! 57 heart    56 swimmer
      ElemData(3,i)=1   ! doubly wetted
     endif
     if ((probtype.eq.56).and.(1.eq.1)) then  ! BOXSWIMMER
      ElemData(3,i)=1   ! doubly wetted
     endif
     if (probtype.eq.562) then  ! 562 whale
      ElemData(3,i)=0   ! 0=singly wetted 1=doubly wetted
     endif
     do j1=1,k
      READ(14,*) IntElem(j1,i) 
      IntElem(j1,i)=IntElem(j1,i)+NumNodesPool
      if (shift_from_zero_node.eq.1) then
       IntElem(j1,i)=IntElem(j1,i)+1
      endif
     enddo
    enddo  ! i, looping faces
    close(14)
 
    print *,"i1,minnodebefore,maxnodebefore ",i1, &
      minnodebefore(1),minnodebefore(2),minnodebefore(3), &
      maxnodebefore(1),maxnodebefore(2),maxnodebefore(3)
    print *,"i1,minnode,maxnode ",i1,minnode(1),minnode(2),minnode(3), &
      maxnode(1),maxnode(2),maxnode(3)
   enddo ! i1

   if ((probtype.eq.55).or.(probtype.eq.58)) then
    theta=1.0
   else if (tcrit.le.tstart+1.0e-8) then
    theta=0.0
   else if (tcrit.ge.tfinish-1.0e-8) then
    theta=1.0
   else
    theta=(tcrit-tstart)/(tfinish-tstart)
   endif

   do i=1,NumNodes
    if (probtype.eq.58) then 
     do dir=1,3
      NodeVel_old(dir,i)=0.0
     enddo
     NodeVel_old(1,i)=2.0*3.14159*plungeramp*plungerfreq*  &
       cos(2.0*3.14159*plungerfreq*curtime)
    else if (tfinish-tstart.gt.1.0e-8) then
     do dir=1,3
      NodeVel_old(dir,i)=(Node_new(dir,i)-Node_old(dir,i))/(tfinish-tstart)
     enddo
    else
     do dir=1,3
      NodeVel_old(dir,i)=0.0
     enddo
    endif
    do dir=1,3
     Node_current(dir,i)=theta*Node_new(dir,i)+(1.0-theta)*Node_old(dir,i)
     NodeVel_new(dir,i)=NodeVel_old(dir,i)
    enddo
   enddo  ! i=1,NumNodes
  endif ! iread=1



  use_temp=0
  if (probtype.eq.55) then
   use_temp=1
   do i=1,NumNodes
    NodeTemp_old(i)=530.0
    NodeTemp_new(i)=530.0
   enddo
  endif
    
  if (ifirst.eq.1) then
   print *,"before CouplerIO_WriteHeader"
   call CouplerIO_WriteHeader(DebugCouple,FlagDim,IntElemDim,UnitCouple, &
    UnitOut,UnitDebug,StdOut,NumNodes,NumIntElems,NumFlags,HdrFlag, &
    Eul2IntNode,ElemData,IntElem,use_temp)
   print *,"after writeheader NumNodes,NumIntElems ",NumNodes,NumIntElems
  endif ! ifirst.eq.1

  print *,"before CouplerIO_WriteNodes NumNodes=",NumNodes
  call CouplerIO_WriteNodes(DebugCouple,istep,curtime,dt, &
     NptFlag,Node_current,NodeVel_new,ActiveIntElem,NodeTemp_new)
  print *,"after writenodes"

    ! ReadForces updates curtime
  print *,"before CouplerIO_ReadForces"
  call CouplerIO_ReadForces(DebugCouple,istepcfd,curtime,dt, &
     ForFlag,NodeForce)  
  timeB=curtime
  dtB=dt
  print *,"after readforces"

  print *,"after geominit  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt,Theta_Dot,Theta ",timeB,dtB,TorqueVel,TorquePos
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
      xblob(1),xblob(3),newxblob(1),newxblob(3),radblob
  print *," xnew=cos(theta)*(x-x_0)+sin(theta)*(z-z_0)+x_0"
  print *," znew=-sin(theta)*(x-x_0)+cos(theta)*(z-z_0)+z_0"

91     FORMAT(I7)
92     FORMAT(f11.3)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)
125    FORMAT(I4,f20.10)

return
end subroutine geominit


subroutine initpaddle(curtime,dt,sdim,istop,istep,probtype, &
  paddle_pos,paddle_vel)
integer(4) :: i,j,k,i1,j1,k1,itimecount,sdim
integer(4) :: inode,iface
real(8) :: curtime,dt
real(8), dimension(3) :: maxnode,minnode,xval,xtemp
real(8) :: tper,tcrit,theta
integer(4) :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
integer(4) :: probtype
character(20) :: dwave,dwave2
real(8) :: ofs,xx,zz
real(8), intent(in) :: paddle_pos,paddle_vel

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=0.0
  TorquePos=paddle_pos
  TorqueVel=paddle_vel

  allocate(HdrFlag(NumFlags(1)))
  allocate(NptFlag(1))
  allocate(ForFlag(NumFlags(4))) 

  xblob(1)=-1.8
  xblob(3)=0.8
  xblob(2)=4.3  ! vertical
  xblobwall(1)=-1.8
  xblobwall(3)=0.8
  xblobwall(2)=4.3  ! vertical
  newxblob(1)=1.45  
  newxblob(3)=1.0
  newxblob(2)=1.0  ! vertical
  newxblobwall(1)=1.05  
  newxblobwall(3)=1.0
  newxblobwall(2)=1.0  ! vertical
  radblob=60.0
  radblobwall=60.0
  denpaddle=0.0035
  dampingpaddle=0.01

  dwave="h/paddle.txt"
  dwave2="h/wall.txt"

  print *,"opening ",dwave
  OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')
  print *,"opening ",dwave2
  OPEN(unit=15,file=dwave2,access='sequential',form="formatted",status='old')

  READ(14,*) NumNodesPaddle
  print *,"NumNodesPaddle ",NumNodesPaddle
  READ(14,*) NumIntElemsPaddle
  print *,"NumIntElemsPaddle ",NumIntElemsPaddle
  READ(14,*) IntElemDimPaddle
  print *,"IntElemDimPaddle ",IntElemDimPaddle

  READ(15,*) NumNodes
  print *,"NumNodes ",NumNodes
  READ(15,*) NumIntElems
  print *,"NumIntElems ",NumIntElems
  READ(15,*) IntElemDim
  print *,"IntElemDim ",IntElemDim

  NumNodes=NumNodes+NumNodesPaddle
  NumIntElems=NumIntElems+NumIntElemsPaddle
  if (IntElemDim.lt.IntElemDimPaddle) then
   IntElemDim=IntElemDimPaddle
  endif

  allocate(Eul2IntNode(NumNodes))
  allocate(ElemData(NumFlags(2),NumIntElems))
  allocate(IntElem(IntElemDim,NumIntElems))
  allocate(ActiveIntElem(NumIntElems))
  allocate(Node_old(3,NumNodes))
  allocate(Node_new(3,NumNodes))
  allocate(Node_current(3,NumNodes))
  allocate(NodeVel_old(3,NumNodes))
  allocate(NodeVel_new(3,NumNodes))
  allocate(NodeTemp_old(NumNodes))
  allocate(NodeTemp_new(NumNodes))
  allocate(NodeForce(3,NumNodes))

  do inode=1,NumNodes
   Eul2IntNode(inode)=inode
   NodeTemp_old(inode)=0.0
   NodeTemp_new(inode)=0.0
   do dir=1,3
    Node_old(dir,inode)=0.0
    Node_new(dir,inode)=0.0
    Node_current(dir,inode)=0.0
    NodeVel_old(dir,inode)=0.0
    NodeVel_new(dir,inode)=0.0
    NodeForce(dir,inode)=0.0
   enddo
  enddo
  do iface=1,NumIntElems
   ActiveIntElem(iface)=.true.
  enddo

  HdrFlag(1)=sdim
  HdrFlag(2)=0     ! stationary body

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
  enddo

  do inode=1,NumNodes
   if (inode.le.NumNodesPaddle) then
    READ(14,*) j,xval(1),xval(2),xval(3)
    if (inode.ne.j) then
     print *,"inode,j mismatch"
     stop
    endif
    do dir=1,3
     xval(dir)=(xval(dir)-xblob(dir))/radblob + newxblob(dir)
    enddo
   else
    READ(15,*) j,xval(1),xval(2),xval(3)
    if (inode.ne.j+NumNodesPaddle) then
     print *,"inode,j mismatch"
     stop
    endif
    do dir=1,3
     xval(dir)=(xval(dir)-xblobwall(dir))/radblobwall + newxblobwall(dir)
    enddo
   endif


   do dir=1,3
    xtemp(dir)=xval(dir)
   enddo
   xval(1)=xtemp(1)
   xval(2)=xtemp(3)
   xval(3)=xtemp(2)
   do dir=1,3
    if ((minnode(dir).gt.xval(dir)).or.(inode.eq.1)) then
     minnode(dir)=xval(dir)
    endif
    if ((maxnode(dir).lt.xval(dir)).or.(inode.eq.1)) then
     maxnode(dir)=xval(dir)
    endif
   enddo
    
   do dir=1,3
    Node_old(dir,inode)=xval(dir)
    Node_new(dir,inode)=xval(dir)
   enddo
      
  enddo  ! inode=1,NumNodes
 
  READ(14,*) j
  if (j.ne.NumIntElemsPaddle) then
   print *,"face mismatch"
   stop
  endif

  READ(15,*) j
  if (j.ne.NumIntElems-NumIntElemsPaddle) then
   print *,"face mismatch"
   stop
  endif

  do iface=1,NumIntElems
   if (iface.le.NumIntElemsPaddle) then
    READ(14,*) j
    if (iface.ne.j) then
     print *,"face mismatch"
     stop
    endif
    READ(14,*) k
    if (k.gt.IntElemDim) then
     print *,"14 too many vertices k,IntElemDim ",k,IntElemDim
     stop
    endif
    do j1=1,k
     READ(14,*) IntElem(j1,iface) 
    enddo
   else
    READ(15,*) j
    if (iface.ne.j+NumIntElemsPaddle) then
     print *,"face mismatch"
     stop
    endif
    READ(15,*) k
    if (k.gt.IntElemDim) then
     print *,"15 too many vertices k,IntElemDim ",k,IntElemDim
     stop
    endif
    do j1=1,k
     READ(15,*) IntElem(j1,iface) 
     IntElem(j1,iface)=IntElem(j1,iface)+NumNodesPaddle
    enddo
   endif

   ElemData(1,iface)=k   ! number of nodes in element
   ElemData(2,iface)=1   ! part number
   ElemData(3,iface)=0   ! singly wetted
   ElemData(3,iface)=1  ! doubly wetted
   if (iface.gt.NumIntElemsPaddle) then
    ElemData(3,iface)=1  ! doubly wetted
   endif
  enddo  ! iface, looping faces

  close(14)
  close(15)
 
  print *,"minnode,maxnode ",minnode(1),minnode(2),minnode(3), &
      maxnode(1),maxnode(2),maxnode(3)


  do inode=1,NumNodes
   do dir=1,3
    NodeVel_old(dir,inode)=0.0
   enddo
   do dir=1,3
    Node_current(dir,inode)=Node_new(dir,inode)
    NodeVel_new(dir,inode)=NodeVel_old(dir,inode)
   enddo
  enddo  ! inode=1,NumNodes

  do inode=1,NumNodesPaddle
   xx=Node_current(1,inode)-newxblob(1)
   zz=Node_current(3,inode)-newxblob(3)
   Node_new(1,inode)=xx*cos(TorquePos)+zz*sin(TorquePos)+newxblob(1)
   Node_new(3,inode)=-xx*sin(TorquePos)+zz*cos(TorquePos)+newxblob(3)
   xx=Node_new(1,inode)-newxblob(1)
   zz=Node_new(3,inode)-newxblob(3)
   NodeVel_new(1,inode)=TorqueVel*zz
   NodeVel_new(3,inode)=-TorqueVel*xx
  enddo

  use_temp=0
  print *,"before CouplerIO_WriteHeader"
  call CouplerIO_WriteHeader(DebugCouple,FlagDim,IntElemDim,UnitCouple, &
   UnitOut,UnitDebug,StdOut,NumNodes,NumIntElems,NumFlags,HdrFlag, &
    Eul2IntNode,ElemData,IntElem,use_temp)
   print *,"after writeheader NumNodes,NumIntElems ",NumNodes,NumIntElems

  print *,"before CouplerIO_WriteNodes"
  call CouplerIO_WriteNodes(DebugCouple,istep,curtime,dt, &
     NptFlag,Node_new,NodeVel_new,ActiveIntElem,NodeTemp_new)
  print *,"after writenodes"
  print *,"before CouplerIO_ReadForces"
  call CouplerIO_ReadForces(DebugCouple,istepcfd,curtime,dt, &
     ForFlag,NodeForce)  
  print *,"after readforces"

91     FORMAT(I7)
92     FORMAT(f11.3)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)
125    FORMAT(I4,f20.10)

return
end subroutine initpaddle

subroutine rotatexyz(xx,angles,newxx)
real(8), dimension(3), intent(in) :: xx,angles
real(8), dimension(3), intent(out) :: newxx

real(8), dimension(3) :: xxz,xxy

! counter-clockwise rotation about z axis
 xxz(1)=xx(1)*cos(angles(1))-xx(2)*sin(angles(1))
 xxz(2)=xx(1)*sin(angles(1))+xx(2)*cos(angles(1))
 xxz(3)=xx(3)

! counter-clockwise rotation about y axis
 xxy(1)=xxz(1)*cos(angles(2))-xxz(3)*sin(angles(2))
 xxy(3)=xxz(1)*sin(angles(2))+xxz(3)*cos(angles(2))
 xxy(2)=xxz(2)

! counter-clockwise rotation about x axis
 newxx(2)=xxy(2)*cos(angles(3))-xxy(3)*sin(angles(3))
 newxx(3)=xxy(2)*sin(angles(3))+xxy(3)*cos(angles(3))
 newxx(1)=xxy(1)

return
end subroutine rotatexyz


subroutine velrotatexyz(xx,velangles,velxx)
real(8), dimension(3), intent(in) :: xx,velangles
real(8), dimension(3), intent(out) :: velxx

real(8), dimension(3) :: xxz,xxy
integer(4) :: dir

 do dir=1,3
  velxx(dir)=0.0
 enddo

! counter-clockwise rotation about z axis
 velxx(1)=velxx(1)-xx(2)*velangles(1)
 velxx(2)=velxx(2)+xx(1)*velangles(1)

! counter-clockwise rotation about y axis
 velxx(1)=velxx(1)-xx(3)*velangles(2)
 velxx(3)=velxx(3)+xx(1)*velangles(2)

! counter-clockwise rotation about x axis
 velxx(2)=velxx(2)-xx(3)*velangles(3)
 velxx(3)=velxx(3)+xx(2)*velangles(3)

return
end subroutine velrotatexyz


subroutine initSteamCleaning(curtime,dt,sdim,istop,istep,probtype, &
  restartTorquePos,restartTorqueVel,restartCenterPos,restartCenterVel)
integer(4) :: i,j,k,sdim
integer(4) :: inode,iface
real(8) :: curtime,dt
real(8), dimension(3) :: maxnode,minnode,xval
real(8), dimension(3) :: xx,rotxx,velxx
real(8) :: tper,tcrit,theta
integer(4) :: dir,istep,istop,istepcfd
integer(4) :: probtype
integer(4) :: ntheta,npsi
real(8), dimension(3), intent(in) :: restartTorquePos,restartTorqueVel, &
  restartCenterPos,restartCenterVel

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=0.0
  do dir=1,3
   torquePosDust(dir)=restartTorquePos(dir)
   torqueVelDust(dir)=restartTorqueVel(dir)
   centerDust(dir)=restartCenterPos(dir)
   centerVelDust(dir)=restartCenterVel(dir)
   centerStartDust(dir)=0.0
  enddo
  centerStartDust(sdim)=raddust
  centerStartDust(1)=0.0
  centerStartDust(sdim-1)=0.0

  tempdust=0.0  ! temperature of the dust particle
  dendust=1.0
  floordust=0.0
  adheredust=1.0e+20

  allocate(HdrFlag(NumFlags(1)))
  allocate(NptFlag(1))
  allocate(ForFlag(NumFlags(4))) 

  if (sdim.eq.2) then
   IntElemDim=2
   ntheta=32
   npsi=1
   NumIntElems=ntheta
   NumNodes=ntheta+1
  else if (sdim.eq.3) then
   IntElemDim=4
   ntheta=32
   npsi=32
   NumIntElems=ntheta*npsi
   NumNodes=(ntheta+1)*(npsi+1)
  endif
    
  print *,"NumNodes ",NumNodes
  print *,"NumIntElems ",NumIntElems
  print *,"IntElemDim ",IntElemDim

  allocate(Eul2IntNode(NumNodes))
  allocate(ElemData(NumFlags(2),NumIntElems))
  allocate(IntElem(IntElemDim,NumIntElems))
  allocate(ActiveIntElem(NumIntElems))
  allocate(Node_old(3,NumNodes))
  allocate(Node_new(3,NumNodes))
  allocate(Node_current(3,NumNodes))
  allocate(NodeVel_old(3,NumNodes))
  allocate(NodeVel_new(3,NumNodes))
  allocate(NodeTemp_old(NumNodes))
  allocate(NodeTemp_new(NumNodes))
  allocate(NodeForce(3,NumNodes))

  do inode=1,NumNodes
   Eul2IntNode(inode)=inode
   NodeTemp_old(inode)=tempdust
   NodeTemp_new(inode)=tempdust
   do dir=1,3
    Node_old(dir,inode)=0.0
    Node_new(dir,inode)=0.0
    Node_current(dir,inode)=0.0
    NodeVel_old(dir,inode)=0.0
    NodeVel_new(dir,inode)=0.0
    NodeForce(dir,inode)=0.0
   enddo
  enddo
  do iface=1,NumIntElems
   ActiveIntElem(iface)=.true.
  enddo

  HdrFlag(1)=sdim
  HdrFlag(2)=0     ! stationary body

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
  enddo

  print *,"initializing dust particle: radus,x,y,z "
  print *,raddust,centerStartDust(1),centerStartDust(2), &
   centerStartDust(3)
  print *,"initial torque ",torquePosDust(1),torquePosDust(2), &
   torquePosDust(3)
  print *,"initial offset ",centerDust(1),centerDust(2), &
   centerDust(3)
! init Node_new and element map
  call circlegeom(centerStartDust(1),centerStartDust(2), &
   centerStartDust(3),raddust,sdim,ntheta,npsi)

  do inode=1,NumNodes
   do dir=1,3
    xval(dir)=Node_new(dir,inode)
   enddo
   do dir=1,3
    if ((minnode(dir).gt.xval(dir)).or.(inode.eq.1)) then
     minnode(dir)=xval(dir)
    endif
    if ((maxnode(dir).lt.xval(dir)).or.(inode.eq.1)) then
     maxnode(dir)=xval(dir)
    endif
   enddo
  enddo  ! inode=1,NumNodes
 
  print *,"minnode,maxnode ",minnode(1),minnode(2),minnode(3), &
      maxnode(1),maxnode(2),maxnode(3)


  do inode=1,NumNodes
   do dir=1,3
    NodeVel_old(dir,inode)=0.0
   enddo

! Node_current stores nodes prior to rotation or translation.

   do dir=1,3
    Node_current(dir,inode)=Node_new(dir,inode)
    NodeVel_new(dir,inode)=NodeVel_old(dir,inode)
   enddo
  enddo  ! inode=1,NumNodes

  do inode=1,NumNodes
   do dir=1,3
    xx(dir)=Node_current(dir,inode)-centerStartDust(dir)
   enddo
   call rotatexyz(xx,torquePosDust,rotxx)
   do dir=1,3
    Node_new(dir,inode)=rotxx(dir)+centerStartDust(dir)+ &
      centerDust(dir)
   enddo 

   call velrotatexyz(rotxx,torqueVelDust,velxx)
   do dir=1,3
    NodeVel_new(dir,inode)=velxx(dir)+centerVelDust(dir)
   enddo
  enddo

  use_temp=1
  print *,"before CouplerIO_WriteHeader"
  call CouplerIO_WriteHeader(DebugCouple,FlagDim,IntElemDim,UnitCouple, &
   UnitOut,UnitDebug,StdOut,NumNodes,NumIntElems,NumFlags,HdrFlag, &
    Eul2IntNode,ElemData,IntElem,use_temp)
   print *,"after writeheader NumNodes,NumIntElems ",NumNodes,NumIntElems

  print *,"before CouplerIO_WriteNodes"
  call CouplerIO_WriteNodes(DebugCouple,istep,curtime,dt, &
     NptFlag,Node_new,NodeVel_new,ActiveIntElem,NodeTemp_new)
  print *,"after writenodes"
  print *,"before CouplerIO_ReadForces"
  call CouplerIO_ReadForces(DebugCouple,istepcfd,curtime,dt, &
     ForFlag,NodeForce)  
  print *,"after readforces"

return
end subroutine initSteamCleaning



subroutine initship(curtime,dt,sdim,istop,istep,probtype, &
  paddle_pos,paddle_vel)
integer(4) :: i,j,k,i1,j1,k1,itimecount,sdim
integer(4) :: inode,iface
real(8) :: curtime,dt
real(8), dimension(3) :: maxnode,minnode,xval,xtemp
real(8) :: tper,tcrit,theta
integer(4) :: iper,icrit,iread,ewave,fwave,dir,istep,istop,istepcfd
integer(4) :: probtype,filler
character(40) :: dwave,dwave2
real(8) :: ofs,xx,zz
real(8), intent(in) :: paddle_pos,paddle_vel

  DebugCouple=.true.
  FlagDim=4
  UnitCouple=1
  UnitOut=2
  UnitDebug=3
  StdOut=6

  NumFlags(1)=2
  NumFlags(2)=3
  NumFlags(3)=0
  NumFlags(4)=4  ! 3 in the manual ???

  timeB=curtime
  dtB=0.0
  TorquePos=paddle_pos
  TorqueVel=paddle_vel

  allocate(HdrFlag(NumFlags(1)))
  allocate(NptFlag(1))
  allocate(ForFlag(NumFlags(4))) 

  xblob(1)=0.0
  xblob(2)=0.0
  xblob(3)=0.0
  newxblob(1)=0.0
  newxblob(2)=0.0
  newxblob(3)=0.0  
  radblob=1.0
  denpaddle=1.0
  dampingpaddle=0.01

  dwave="5415_froude4136.cas"

  print *,"opening ",dwave
  OPEN(unit=14,file=dwave,access='sequential',form="formatted",status='old')

  READ(14,*) NumNodes,NumIntElems
  IntElemDim=3

  allocate(Eul2IntNode(NumNodes))
  allocate(ElemData(NumFlags(2),NumIntElems))
  allocate(IntElem(IntElemDim,NumIntElems))
  allocate(ActiveIntElem(NumIntElems))
  allocate(Node_old(3,NumNodes))
  allocate(Node_new(3,NumNodes))
  allocate(Node_current(3,NumNodes))
  allocate(NodeVel_old(3,NumNodes))
  allocate(NodeVel_new(3,NumNodes))
  allocate(NodeTemp_old(NumNodes))
  allocate(NodeTemp_new(NumNodes))
  allocate(NodeForce(3,NumNodes))

  do inode=1,NumNodes
   Eul2IntNode(inode)=inode
   NodeTemp_old(inode)=0.0
   NodeTemp_new(inode)=0.0
   do dir=1,3
    Node_old(dir,inode)=0.0
    Node_new(dir,inode)=0.0
    Node_current(dir,inode)=0.0
    NodeVel_old(dir,inode)=0.0
    NodeVel_new(dir,inode)=0.0
    NodeForce(dir,inode)=0.0
   enddo
  enddo
  do iface=1,NumIntElems
   ActiveIntElem(iface)=.true.
  enddo

  HdrFlag(1)=sdim
  HdrFlag(2)=0     ! stationary body

  do dir=1,3
   maxnode(dir)=0.0
   minnode(dir)=0.0
  enddo

  do inode=1,NumNodes
   READ(14,*) xval(1),xval(2),xval(3)
   do dir=1,3
     xval(dir)=(xval(dir)-xblob(dir))/radblob + newxblob(dir)
   enddo

   do dir=1,3
    xtemp(dir)=xval(dir)
   enddo
   xval(1)=-xtemp(1)
   xval(2)=xtemp(2)
   xval(3)=xtemp(3)
   do dir=1,3
    if ((minnode(dir).gt.xval(dir)).or.(inode.eq.1)) then
     minnode(dir)=xval(dir)
    endif
    if ((maxnode(dir).lt.xval(dir)).or.(inode.eq.1)) then
     maxnode(dir)=xval(dir)
    endif
   enddo
    
   do dir=1,3
    Node_old(dir,inode)=xval(dir)
    Node_new(dir,inode)=xval(dir)
   enddo
      
  enddo  ! inode=1,NumNodes
 
  do iface=1,NumIntElems
   READ(14,*) IntElem(3,iface),IntElem(2,iface),IntElem(1,iface), &
    filler

   ElemData(1,iface)=3   ! number of nodes in element
   ElemData(2,iface)=1   ! part number
   ElemData(3,iface)=0   ! singly wetted
  enddo  ! iface, looping faces

  close(14)
 
  print *,"minnode,maxnode ",minnode(1),minnode(2),minnode(3), &
      maxnode(1),maxnode(2),maxnode(3)


  do inode=1,NumNodes
   do dir=1,3
    NodeVel_old(dir,inode)=0.0
   enddo
   do dir=1,3
    Node_current(dir,inode)=Node_new(dir,inode)
    NodeVel_new(dir,inode)=NodeVel_old(dir,inode)
   enddo
  enddo  ! inode=1,NumNodes

  use_temp=0
  print *,"before CouplerIO_WriteHeader"
  call CouplerIO_WriteHeader(DebugCouple,FlagDim,IntElemDim,UnitCouple, &
   UnitOut,UnitDebug,StdOut,NumNodes,NumIntElems,NumFlags,HdrFlag, &
    Eul2IntNode,ElemData,IntElem,use_temp)
   print *,"after writeheader NumNodes,NumIntElems ",NumNodes,NumIntElems

  print *,"before CouplerIO_WriteNodes"
  call CouplerIO_WriteNodes(DebugCouple,istep,curtime,dt, &
     NptFlag,Node_new,NodeVel_new,ActiveIntElem,NodeTemp_new)
  print *,"after writenodes"
  print *,"before CouplerIO_ReadForces"
  call CouplerIO_ReadForces(DebugCouple,istepcfd,curtime,dt, &
     ForFlag,NodeForce)  
  print *,"after readforces"

191    FORMAT(I12,I12)
91     FORMAT(I7)
92     FORMAT(f11.3)
193    FORMAT(E15.11,E15.11,E15.11)
194    FORMAT(I12,I12,I12,I12)
93     FORMAT(i7,f11.3,f11.3,f11.3)
100    FORMAT(f12.7)
115    FORMAT(I4)
125    FORMAT(I4,f20.10)

return
end subroutine initship

subroutine overall_solid_advance()

integer(4) :: ifirst

 ifirst=0
 sci_istep=sci_istep+1

 if ((sci_probtype.eq.52).or.(sci_probtype.eq.57).or. &
     (sci_probtype.eq.56).or.(sci_probtype.eq.561).or. &
     ((sci_probtype.eq.55).and.(raddust.eq.0.0)).or. &
     (sci_probtype.eq.58).or. &
     (sci_probtype.eq.562)) then
    ! sci_curtime updated via ReadForces
  call geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop,sci_istep, &
   sci_probtype)
 else
  call advance_solid(sci_sdim,sci_problo,sci_probhi,sci_curtime,sci_dt, &
    sci_istop,sci_istep,sci_probtype)
 endif

return
end subroutine overall_solid_advance

subroutine overall_solid_init()
integer(4) :: ifirst,dir
real(8) :: paddle_pos,paddle_vel
real(8), dimension(3) :: restartTorquePos,restartTorqueVel, &
  restartCenterPos,restartCenterVel


 sci_sdim=2       ! CHANGE ME

 raddust=0.0

 sci_probtype=36  ! hydrobulge
 sci_probtype=58  ! plunger
 sci_probtype=57  ! heart
 sci_probtype=55  ! bubble cloud with sphere 
 sci_probtype=55  ! boiling, with dust particle, raddust>0
 sci_probtype=56  ! swimmer
 sci_probtype=561 ! swimmer with pool
 sci_probtype=562 ! whale
 sci_probtype=52  ! hand
 sci_probtype=50  ! paddle
 sci_probtype=9   ! ship

 sci_probtype=56
 raddust=0.01

! sends header message to fluids code
! sends node message to fluids code 
! reads force message

 sci_curtime=0.0
 sci_dt=0.0

! --------------- MODIFY THESE VARIABLES IF RESTARTING! ---------------
 do dir=1,3
  restartTorquePos(dir)=0.0
  restartTorqueVel(dir)=0.0
  restartCenterPos(dir)=0.0
  restartCenterVel(dir)=0.0
 enddo    
 paddle_pos=0.0
 paddle_vel=0.0
 sci_curtime=0.0
 sci_dt=0.0
! --------------- END SECTION TO MODIFY IF RESTARTING! ---------------

 sci_istop=0 
 sci_istep=0
 ifirst=1


 if ((sci_probtype.eq.52).or.(sci_probtype.eq.57).or. &
     (sci_probtype.eq.56).or.(sci_probtype.eq.561).or. &
     ((sci_probtype.eq.55).and.(raddust.eq.0.0)).or. &
     (sci_probtype.eq.58).or. &
     (sci_probtype.eq.562)) then
  call geominit(sci_curtime,sci_dt,ifirst,sci_sdim,sci_istop, &
    sci_istep,sci_probtype)
 else if (sci_probtype.eq.50) then
  call initpaddle(sci_curtime,sci_dt,sci_sdim,sci_istop,sci_istep, &
    sci_probtype,paddle_pos,paddle_vel)
 else if (sci_probtype.eq.9) then
  call initship(sci_curtime,sci_dt,sci_sdim,sci_istop,sci_istep,sci_probtype, &
    paddle_pos,paddle_vel)
 else if ((sci_probtype.eq.55).and.(raddust.gt.0.0)) then
  call initSteamCleaning(sci_curtime,sci_dt,sci_sdim,sci_istop,sci_istep, &
    sci_probtype,restartTorquePos,restartTorqueVel,restartCenterPos, &
    restartCenterVel)
 else
  call initialize_solid(sci_probtype,sci_sdim,sci_problo,sci_probhi,sci_dt) 
 endif
 print *,"after initialize solid dt=",sci_dt

return
end subroutine overall_solid_init

subroutine circlegeom(xcen,ycen,zcen,radcircle,ndim,ntheta,npsi)
integer(4) :: i,j,k,ndim,ntheta,npsi
real(8)    :: xcen,ycen,zcen,radcircle
real(8)    :: theta,psi,Pi

 Pi=2.0*1.570796327

 if (ndim.eq.2) then
  if (ntheta.ne.NumIntElems) then
   print *,"ntheta <> NumIntElems"
   stop
  endif
  if (ntheta+1.ne.NumNodes) then
   print *,"ntheta+1 <> NumNodes"
   stop
  endif
   
  do i=1,ntheta+1
   theta=(i-1.0)*2.0*Pi/ntheta
   Node_new(1,i)=xcen+radcircle*cos(theta)
   Node_new(2,i)=ycen+radcircle*sin(theta)
   Node_new(3,i)=0.0
  enddo

  do i=1,NumIntElems
   ElemData(1,i)=2
   ElemData(2,i)=i
   ElemData(3,i)=0  ! singly wetted
   IntElem(1,i)=i
   IntElem(2,i)=i+1
  enddo
 else if (ndim.eq.3) then
  if (ntheta*npsi.ne.NumIntElems) then
   print *,"ntheta*npsi<>NumIntElems"
   stop
  endif
  if ((ntheta+1)*(npsi+1).ne.NumNodes) then
   print *,"(ntheta+1)*(npsi+1)<>NumNodes"
   stop
  endif
  k=1
  do i=1,ntheta+1
  do j=1,npsi+1
   theta=(i-1.0)*2.0*Pi/ntheta
   psi=(j-1.0)*Pi/npsi-0.5*Pi
   Node_new(1,k)=xcen+radcircle*cos(psi)*cos(theta)
   Node_new(2,k)=ycen+radcircle*cos(psi)*sin(theta)
   Node_new(3,k)=zcen+radcircle*sin(psi)
   k=k+1
  enddo
  enddo
  i=1
  j=1
  do k=1,NumIntElems
   ElemData(1,k)=4
   ElemData(2,k)=k
   ElemData(3,k)=0  ! singly wetted
   IntElem(4,k)=(i-1)*(npsi+1)+j
   IntElem(3,k)=i*(npsi+1)+j
   IntElem(2,k)=i*(npsi+1)+j+1
   IntElem(1,k)=(i-1)*(npsi+1)+j+1
   j=j+1
   if (j.gt.npsi) then
    j=1
    i=i+1
   endif
  enddo
 else
  print *,"ndim invalid"
  stop
 endif

return
end subroutine circlegeom

subroutine initialize_solid(probtype,sdim,problo,probhi,startdt)
integer(4) :: probtype,sdim
real(8), dimension(3) :: problo,probhi
real(8) :: xblob2,yblob2,zblob2,radblob2
real(8) :: adv_vel
integer(4) :: i,j,k,dir,istep
real(8) :: starttime,startdt,theta
integer(4) :: iflip

 timeB=0.0
 dtB=0.0

 DebugCouple=.true.
 FlagDim=4
 UnitCouple=1
 UnitOut=2
 UnitDebug=3
 StdOut=6
 NumFlags(1)=2
 NumFlags(2)=3
 NumFlags(3)=0
 NumFlags(4)=4  ! 3 in the manual ???

 allocate(HdrFlag(NumFlags(1)))
 allocate(NptFlag(1))
 allocate(ForFlag(NumFlags(4))) 

 do i=1,3
  problo(i)=0.0
  probhi(i)=0.0
 enddo

 if (probtype.eq.36) then
  problo(1)=0.0
  problo(2)=0.0
  problo(3)=0.0
  probhi(1)=5.715
  probhi(2)=5.715
  probhi(3)=11.43
  adv_vel=0.0
  xblob2=4.445
  yblob2=11.43
  radblob2=0.635 
  if (sdim.eq.2) then
   NumIntElems=1
   NumNodes=2
   IntElemDim=2
  else if (sdim.eq.3) then
   NumIntElems=64
   NumNodes=2+NumIntElems*2
   IntElemDim=4
  endif
  HdrFlag(2)=0  ! stationary body 
 else if (probtype.eq.42) then
  iflip=1

  if (sdim.eq.3) then

  problo(1)=0.0
  problo(2)=0.0
  problo(3)=0.0
  probhi(1)=50.0
  probhi(2)=50.0 
  probhi(3)=50.0 
  adv_vel=-150.0
  xblob2=15.0
  yblob2=12.0    ! in cfd code, it is 13.0
  radblob2=5.0  ! in cfd code, it is 5.08
  xblob(1)=25.0
  xblob(2)=25.0    ! in cfd code, position of bubble
  xblob(3)=25.0
  if (radblob2.eq.0.0) then
   NumIntElems=1
   NumNodes=4
  else if (radblob2.gt.0.0) then
   NumIntElems=6
   NumNodes=8
  else
   print *,"radblob2 invalid"
   stop
  endif
  IntElemDim=4

  else

  problo(1)=0.0
  problo(2)=0.0
  probhi(1)=50.0
  probhi(2)=50.0 
  adv_vel=-40.0
  xblob2=40.0
  yblob2=7.0    ! in cfd code, it is 13.0
  radblob2=0.0   ! in cfd code, it is 5.08
  xblob(1)=0.0
  xblob(2)=25.0    ! in cfd code, position of bubble
  if (radblob2.eq.0.0) then
   NumIntElems=1
   NumNodes=2
   if (iflip.eq.1) then  ! setup for radblob2=0 only
    NumIntElems=2
    NumNodes=3
   endif
  else if (radblob2.gt.0.0) then
   NumIntElems=4
   NumNodes=4
  else
   print *,"radblob2 invalid"
   stop
  endif
  IntElemDim=sdim

  endif

  HdrFlag(2)=0  ! stationary body
 else 
  print *,"probtype invalid"
  stop
 endif

 allocate(Eul2IntNode(NumNodes))
 allocate(ElemData(NumFlags(2),NumIntElems))
 allocate(IntElem(IntElemDim,NumIntElems))
 allocate(ActiveIntElem(NumIntElems))
 allocate(Node_old(3,NumNodes))
 allocate(Node_new(3,NumNodes))
 allocate(Node_current(3,NumNodes))
 allocate(NodeVel_old(3,NumNodes))
 allocate(NodeVel_new(3,NumNodes))
 allocate(NodeTemp_old(NumNodes))
 allocate(NodeTemp_new(NumNodes))
 allocate(NodeForce(3,NumNodes))
 HdrFlag(1)=sdim

 do i=1,NumIntElems
  ActiveIntElem(i)=.true.
 enddo

 do i=1,NumNodes
  Eul2IntNode(i)=i
  NodeTemp_old(i)=0.0
  NodeTemp_new(i)=0.0
  do dir=1,3
   Node_old(dir,i)=0.0
   Node_new(dir,i)=0.0
   Node_current(dir,i)=0.0
   NodeVel_old(dir,i)=0.0
   NodeVel_new(dir,i)=0.0
   NodeForce(dir,i)=0.0
  enddo
 enddo

 if (probtype.eq.36) then
  if (sdim.eq.2) then
   ElemData(1,1)=2   ! number of nodes in element
   ElemData(2,1)=1   ! part number  
   ElemData(3,1)=0   ! singly wetted
   Node_new(1,1)=xblob2
   Node_new(2,1)=yblob2
   Node_new(1,2)=xblob2
   Node_new(2,2)=0.0

   IntElem(1,1)=1
   IntElem(2,1)=2
  else
   do i=1,NumIntElems+1
    theta=(i-1.0)*(1.570796327)/NumIntElems
    Node_new(1,2*i-1)=xblob2*cos(theta)
    Node_new(2,2*i-1)=xblob2*sin(theta)
    Node_new(3,2*i-1)=0.0
    Node_new(1,2*i)=xblob2*cos(theta)
    Node_new(2,2*i)=xblob2*sin(theta)
    Node_new(3,2*i)=yblob2
   enddo
   do i=1,NumIntElems
    ElemData(1,i)=4
    ElemData(2,i)=i
    ElemData(3,i)=0  ! singly wetted
    IntElem(1,i)=2*i-1
    IntElem(2,i)=2*i
    IntElem(3,i)=2*i+2
    IntElem(4,i)=2*i+1
   enddo
  endif
 else if (probtype.eq.42) then

  if (sdim.eq.3) then

  if (radblob2.eq.0.0) then
   ElemData(1,1)=4   ! number of nodes in element
   ElemData(2,1)=1   ! part number  
   ElemData(3,1)=1   ! doubly wetted
   Node_new(1,1)=xblob(1)-xblob2
   Node_new(2,1)=xblob(2)-xblob2
   Node_new(3,1)=xblob(3)+yblob2

   Node_new(1,2)=xblob(1)+xblob2
   Node_new(2,2)=xblob(2)-xblob2
   Node_new(3,2)=xblob(3)+yblob2

   Node_new(1,3)=xblob(1)+xblob2
   Node_new(2,3)=xblob(2)+xblob2
   Node_new(3,3)=xblob(3)+yblob2

   Node_new(1,4)=xblob(1)-xblob2
   Node_new(2,4)=xblob(2)+xblob2
   Node_new(3,4)=xblob(3)+yblob2

   IntElem(1,1)=1
   IntElem(2,1)=2
   IntElem(3,1)=3
   IntElem(4,1)=4
  else if (radblob2.gt.0.0) then
   Node_new(1,1)=xblob(1)-xblob2
   Node_new(2,1)=xblob(2)-xblob2
   Node_new(3,1)=xblob(3)+yblob2

   Node_new(1,2)=xblob(1)+xblob2
   Node_new(2,2)=xblob(2)-xblob2
   Node_new(3,2)=xblob(3)+yblob2

   Node_new(1,3)=xblob(1)+xblob2
   Node_new(2,3)=xblob(2)+xblob2
   Node_new(3,3)=xblob(3)+yblob2

   Node_new(1,4)=xblob(1)-xblob2
   Node_new(2,4)=xblob(2)+xblob2
   Node_new(3,4)=xblob(3)+yblob2

   Node_new(1,5)=xblob(1)-xblob2
   Node_new(2,5)=xblob(2)-xblob2
   Node_new(3,5)=xblob(3)+yblob2+radblob2

   Node_new(1,6)=xblob(1)+xblob2
   Node_new(2,6)=xblob(2)-xblob2
   Node_new(3,6)=xblob(3)+yblob2+radblob2

   Node_new(1,7)=xblob(1)+xblob2
   Node_new(2,7)=xblob(2)+xblob2
   Node_new(3,7)=xblob(3)+yblob2+radblob2

   Node_new(1,8)=xblob(1)-xblob2
   Node_new(2,8)=xblob(2)+xblob2
   Node_new(3,8)=xblob(3)+yblob2+radblob2

! top
   ElemData(1,1)=4
   ElemData(2,1)=1
   ElemData(3,1)=0   ! singly wetted
   IntElem(1,1)=5
   IntElem(2,1)=6
   IntElem(3,1)=7
   IntElem(4,1)=8

! bottom
   ElemData(1,2)=4
   ElemData(2,2)=2
   ElemData(3,2)=0   ! singly wetted
   IntElem(1,2)=4
   IntElem(2,2)=3
   IntElem(3,2)=2
   IntElem(4,2)=1

! front
   ElemData(1,3)=4
   ElemData(2,3)=3
   ElemData(3,3)=0   ! singly wetted
   IntElem(1,3)=1
   IntElem(2,3)=2
   IntElem(3,3)=6
   IntElem(4,3)=5

! back
   ElemData(1,4)=4
   ElemData(2,4)=4
   ElemData(3,4)=0   ! singly wetted
   IntElem(1,4)=8
   IntElem(2,4)=7
   IntElem(3,4)=3
   IntElem(4,4)=4

! right
   ElemData(1,5)=4
   ElemData(2,5)=5
   ElemData(3,5)=0   ! singly wetted
   IntElem(1,5)=2
   IntElem(2,5)=3
   IntElem(3,5)=7
   IntElem(4,5)=6

! left
   ElemData(1,6)=4
   ElemData(2,6)=6
   ElemData(3,6)=0   ! singly wetted
   IntElem(1,6)=1
   IntElem(2,6)=5
   IntElem(3,6)=8
   IntElem(4,6)=4
  else 
   print *,"radblob2 invalid"
   stop
  endif

  if (adv_vel.ne.0.0) then
   do i=1,NumNodes
    NodeVel_new(3,i)=adv_vel
   enddo
  endif

  else

  if (radblob2.eq.0.0) then
   ElemData(1,1)=2   ! number of nodes in element
   ElemData(2,1)=1   ! part number  
   ElemData(3,1)=1   ! doubly wetted
   Node_new(1,1)=0.0
   Node_new(2,1)=xblob(2)+yblob2
   Node_new(1,2)=xblob2
   Node_new(2,2)=xblob(2)+yblob2
   if (iflip.eq.1) then
    Node_new(1,1)=xblob2
    Node_new(2,1)=xblob(2)-0.5*xblob2
    Node_new(1,2)=yblob2
    Node_new(2,2)=xblob(2)
    Node_new(1,3)=yblob2
    Node_new(2,3)=xblob(2)+0.5*xblob2
    ElemData(1,2)=2   ! number of nodes in element
    ElemData(2,2)=2   ! part number  
    ElemData(3,2)=1   ! doubly wetted
   endif

   IntElem(1,1)=1
   IntElem(2,1)=2
   if (iflip.eq.1) then
    IntElem(1,2)=2
    IntElem(2,2)=3
   endif
  else if (radblob2.gt.0.0) then
   Node_new(1,1)=0.0
   Node_new(2,1)=xblob(2)+yblob2
   Node_new(1,2)=xblob2
   Node_new(2,2)=xblob(2)+yblob2
   Node_new(1,3)=xblob2
   Node_new(2,3)=xblob(2)+yblob2+radblob2
   Node_new(1,4)=0.0
   Node_new(2,4)=xblob(2)+yblob2+radblob2

   ElemData(1,1)=2
   ElemData(2,1)=1
   ElemData(3,1)=0   ! singly wetted
   IntElem(1,1)=1
   IntElem(2,1)=2

   ElemData(1,2)=2
   ElemData(2,2)=2
   ElemData(3,2)=0   ! singly wetted
   IntElem(1,2)=2
   IntElem(2,2)=3

   ElemData(1,3)=2
   ElemData(2,3)=3
   ElemData(3,3)=0   ! singly wetted
   IntElem(1,3)=3
   IntElem(2,3)=4

   ElemData(1,4)=2
   ElemData(2,4)=4
   ElemData(3,4)=0   ! singly wetted
   IntElem(1,4)=4
   IntElem(2,4)=1
  else 
   print *,"radblob2 invalid"
   stop
  endif

  if (adv_vel.ne.0.0) then
   do i=1,NumNodes
    if (iflip.eq.1) then
     NodeVel_new(1,i)=adv_vel
    else
     NodeVel_new(2,i)=adv_vel
    endif
   enddo
  endif

  endif

 else 
  print *,"probtype invalid"
  stop
 endif

 do i=1,NumNodes
  do dir=1,3
   Node_old(dir,i)=Node_new(dir,i)
   Node_current(dir,i)=Node_new(dir,i)
   NodeVel_old(dir,i)=NodeVel_new(dir,i)
  enddo
 enddo

! print *,"j=1,NumIntElems   NumIntElems=",NumIntElems
 do j=1,NumIntElems
 do i=1,3
!  print *,"i,j,ElemData(i,j) ",i,j,ElemData(i,j)
 enddo
 enddo
! print *,"j=1,NumNodes   NumNodes=",NumNodes
 do j=1,NumNodes
 do i=1,3
!  print *,"i,j,Node_new(i,j) ",i,j,Node_new(i,j)
 enddo
 enddo
! print *,"i=1,NumIntElems  j=1,ElemData NumIntElems=",NumIntElems
 do i=1,NumIntElems
 do j=1,ElemData(1,i)
!  print *,"i,j,IntElem(j,i) ",i,j,IntElem(j,i)
 enddo
 enddo
 use_temp=0
 print *,"before CouplerIO_WriteHeader"
 call CouplerIO_WriteHeader(DebugCouple,FlagDim,IntElemDim,UnitCouple, &
   UnitOut,UnitDebug,StdOut,NumNodes,NumIntElems,NumFlags,HdrFlag, &
   Eul2IntNode,ElemData,IntElem,use_temp)
 print *,"after writeheader NumNodes,NumIntElems ",NumNodes,NumIntElems
 istep=0
 starttime=0.0
 startdt=0.0
 print *,"before CouplerIO_WriteNodes"
 call CouplerIO_WriteNodes(DebugCouple,istep,starttime,startdt, &
   NptFlag,Node_new,NodeVel_new,ActiveIntElem,NodeTemp_new)
 print *,"after writenodes"
 print *,"before CouplerIO_ReadForces"
 call CouplerIO_ReadForces(DebugCouple,istep,starttime,startdt, &
  ForFlag,NodeForce)  
 print *,"after readforces"
 
return
end subroutine initialize_solid

subroutine advance_solid(sdim,problo,probhi,curtime,dt,istop,istep,probtype)
integer(4) :: sdim
real(8), dimension(3) :: problo,probhi,xx
real(8) :: curtime,dt
integer(4) :: istop,probtype
integer(4) :: i,j,k,dir,istep,istepcfd
real(8), dimension(3) :: totaltorque,totalforce,rotxx,velxx
real(8) :: xxmag,Pi,masspart

 DebugCouple=.true.
 
 do i=1,NumNodes
  do dir=1,3
   Node_old(dir,i)=Node_new(dir,i)
   NodeVel_old(dir,i)=NodeVel_new(dir,i)
  enddo
 enddo

 print *,"before CouplerIO_WriteNodes"
 call CouplerIO_WriteNodes(DebugCouple,istep,curtime,dt, &
   NptFlag,Node_new,NodeVel_new,ActiveIntElem,NodeTemp_new)
 print *,"advance: after write nodes "
 print *,"before CouplerIO_ReadForces"
 call CouplerIO_ReadForces(DebugCouple,istepcfd,curtime,dt, &
  ForFlag,NodeForce)  
 print *,"advance: after readforces "

 dtB=0.0
 if (curtime.gt.timeB) then
  dtB=curtime-timeB
  timeB=curtime
 endif

 if ((probtype.eq.55).and.(raddust.gt.0.0)) then
  do dir=1,3
   totaltorque(dir)=0.0
   totalforce(dir)=0.0
  enddo
  do i=1,NumNodes
   xxmag=0.0
   do dir=1,sdim
    xx(dir)=Node_new(dir,i)-centerStartDust(dir)-centerDust(dir)
    xxmag=xxmag+xx(dir)**2
   enddo
   xxmag=sqrt(xxmag)

! x-y torque
   totaltorque(1)=totaltorque(1)+xx(1)*NodeForce(2,i)- &
     xx(2)*NodeForce(1,i)
! x-z torque
   totaltorque(2)=totaltorque(2)+xx(3)*NodeForce(1,i)- &
     xx(1)*NodeForce(3,i)
! y-z torque
   totaltorque(3)=totaltorque(3)+xx(2)*NodeForce(3,i)- &
     xx(3)*NodeForce(2,i)
! torque should really sum to zero since there is no viscosity!
   do dir=1,3
    totaltorque(dir)=0.0
   enddo
   
   do dir=1,3
    totalforce(dir)=totalforce(dir)+NodeForce(dir,i)
   enddo

  enddo  ! looping nodes

! torque should really sum to zero since there is no viscosity, so all
! force is normal to the (circular) body.

  Pi=2.0*1.570796327
  if (sdim.eq.2) then
   masspart=dendust*Pi*raddust**2
  else if (sdim.eq.3) then
   masspart=dendust*4.0*Pi*(raddust**3)/3.0
  else
   print *,"sdim invalid"
   stop
  endif

  if ((totalforce(sdim).lt.adheredust).and. &
      (centerDust(sdim).eq.0.0)) then
   do dir=1,3
    totalforce(dir)=0.0
   enddo
  endif

  do dir=1,3
   centerVelDust(dir)=centerVelDust(dir)+dtB*totalforce(dir)/masspart 
   if (dir.eq.sdim) then
    if (centerDust(dir)+dtB*centerVelDust(dir).lt.0.0) then
     centerDust(dir)=0.0
     centerVelDust(dir)=0.0
    endif
   endif
   centerDust(dir)=centerDust(dir)+dtB*centerVelDust(dir)
   torqueVelDust(dir)=torqueVelDust(dir)+ &
     dtB*totaltorque(dir)/(masspart*(raddust**2))
   torquePosDust(dir)=torquePosDust(dir)+ &
     dtB*torqueVelDust(dir)
  enddo

  do i=1,NumNodes
   do dir=1,3
    NodeForce(dir,i)=0.0
    NodeVel_new(dir,i)=0.0
   enddo
   do dir=1,3
    xx(dir)=Node_current(dir,i)-centerStartDust(dir)
   enddo
   call rotatexyz(xx,torquePosDust,rotxx)
   do dir=1,3
    Node_new(dir,i)=rotxx(dir)+centerStartDust(dir)+ &
      centerDust(dir)
   enddo 
   call velrotatexyz(rotxx,torqueVelDust,velxx)
   do dir=1,3
    NodeVel_new(dir,i)=velxx(dir)+centerVelDust(dir)
   enddo

  enddo  ! traversing nodes

  print *,"after advance_solid_dust_particle  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt ",timeB,dtB
  print *,"old x,y,z, offset x,y,z ", &
     centerStartDust(1),centerStartDust(2),centerStartDust(3), &
     centerDust(1),centerDust(2),centerDust(3)
  print *,"totaltorque ",totaltorque(1),totaltorque(2), &
   totaltorque(3)
  print *,"totalforce ",totalforce(1),totalforce(2), &
   totalforce(3)


 else if (probtype.eq.50) then
! x-z torque
  totaltorque(2)=0.0
  do i=1,NumNodesPaddle
   totaltorque(2)=totaltorque(2)+ &
        NodeForce(1,i)*(Node_new(3,i)-newxblob(3))- &
        NodeForce(3,i)*(Node_new(1,i)-newxblob(1))
  enddo
  totaltorque(2)=totaltorque(2)/(denpaddle*NumNodesPaddle) 
  do i=1,NumNodes
   do dir=1,3
    NodeForce(dir,i)=0.0
    NodeVel_new(dir,i)=0.0
   enddo
  enddo
  TorqueVel=(TorqueVel+dtB*totaltorque(2))/(1.0+dtB*dampingpaddle)
  TorquePos=TorquePos+dtB*TorqueVel
  do i=1,NumNodesPaddle
   xx(1)=Node_current(1,i)-newxblob(1)
   xx(3)=Node_current(3,i)-newxblob(3)
   Node_new(1,i)=xx(1)*cos(TorquePos)+xx(3)*sin(TorquePos)+newxblob(1)
   Node_new(3,i)=-xx(1)*sin(TorquePos)+xx(3)*cos(TorquePos)+newxblob(3)
   xx(1)=Node_new(1,i)-newxblob(1)
   xx(3)=Node_new(3,i)-newxblob(3)
   NodeVel_new(1,i)=TorqueVel*xx(3)
   NodeVel_new(3,i)=-TorqueVel*xx(1)
  enddo

  print *,"after advance_solid_paddle  curtime,dt,istep ",curtime, &
    dt,istep
  print *,"time,dt,Theta_Dot,Theta ",timeB,dtB,TorqueVel,TorquePos
  print *,"old xcenter,zcenter, new xcenter,zcenter,scale ", &
     xblob(1),xblob(3),newxblob(1),newxblob(3),radblob
  print *," xnew=cos(theta)*(x-x_0)+sin(theta)*(z-z_0)+x_0"
  print *," znew=-sin(theta)*(x-x_0)+cos(theta)*(z-z_0)+z_0"

 else if (probtype.eq.9) then
  do i=1,NumNodes
   do dir=1,3
    NodeVel_new(dir,i)=0.0
   enddo
  enddo
 else
  do i=1,NumNodes
   do dir=1,3
    Node_new(dir,i)=Node_new(dir,i)+dtB*NodeVel_new(dir,i)
   enddo
  enddo
 endif

return
end subroutine advance_solid

end module SOLIDCouplerIO

PROGRAM solid

use SOLIDCouplerIO, only : overall_solid_init, overall_solid_advance


  call overall_solid_init()

  do while (1.eq.1)
   call overall_solid_advance()
  enddo

end PROGRAM solid
