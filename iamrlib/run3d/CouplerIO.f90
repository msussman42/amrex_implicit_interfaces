module CouplerIO
! USE DFLIB, only : sleepqq  !unix
!USE IFPORT, only : sleepqq  !unix
! for pgf95, "sleep"

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!This file contains coupled subroutines that echange information between the Euler and
!structure code.
!
! subroutines:
!   CouplerIO_WriteHeader
!      Structure: write element conectivity
!   CouplerIO_ReadHeader
!      Fluid: read element conectivity
!   CouplerIO_WriteNodes
!      Structure: write node locations and velocities
!   CouplerIO_ReadNodes
!      Fluid: read node locations and velocities
!   CouplerIO_WriteForces
!      Fluid: write Node Loads
!   CouplerIO_ReadHeader
!      Structure: read Node Loads
!   SafeOpen
!      recovers from an opening error
!   WaitOnFile
!      waits for file to appear
!   WaitRemoveFile
!      waits for file to be removed
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

integer(4), private :: NumNodes,NumIntElems,UnitCouple,UnitOut,UnitDebug,StdOut
integer(4), parameter, private  :: WAITTIME = 1, MAXLOOP=1000000000 
integer(4), parameter, private  :: FLAG_DIM=4
integer(4), dimension(FLAG_DIM), private :: NumFlags
integer(4), private :: GlobalID,GlobalComm,FluidMaster,LagMaster

integer(4), private :: use_temp

contains

! +++++++++++++++++++++++++++++ CouplerIO_WriteHeader +++++++++++++++++++++++++++++++++++++
subroutine CouplerIO_WriteHeader(Debug,FlagDim,IntElemDim,UnitCoupleIn,UnitOutIn,UnitDebugIn, &
                                 StdOutIn,NumNodesIn,NumIntElemsIn,NumFlagsIn,HdrFlag,Eul2LagNode, &
                                 ElemData,IntElem,use_tempIn)
!  Called by Structure Code at start of run: 
!      reads surface element connectivity written by Finite Element code at start of run 

!  ----- arguments -----
   logical(4), intent(in) :: Debug
   integer(4), intent(in) :: FlagDim
   integer(4), intent(in) :: IntElemDim   
   integer(4), intent(in) :: UnitCoupleIn
   integer(4), intent(in) :: UnitOutIn
   integer(4), intent(in) :: UnitDebugIn
   integer(4), intent(in) :: StdOutIn
   integer(4), intent(in) :: NumNodesIn
   integer(4), intent(in) :: NumIntElemsIn
   integer(4), intent(in) :: use_tempIn
   integer(4), dimension(FlagDim), intent(in) :: NumFlagsIn
   integer(4), dimension(NumFlagsIn(1)), intent(in) :: HdrFlag
   integer(4), dimension(NumNodesIn), intent(in) :: Eul2LagNode
   integer(4), dimension(NumFlagsIn(2),NumIntElemsIn), intent(in) :: ElemData   
   integer(4), dimension(IntElemDim,NumIntElemsIn), intent(in) :: IntElem

!  ----- local variables -----
   integer(4) :: j,k
   character(40) :: filename

   use_temp=use_tempIn
   NumNodes=NumNodesIn
   NumIntElems=NumIntElemsIn
   Unitcouple=UnitCoupleIn
   UnitOut=UnitOutIn
   UnitDebug=UnitDebugIn
   StdOut=StdOutIn

   if(FLAG_DIM/=FlagDim)then
      write(UnitOut,*)'FLAG_DIM = ',FLAG_DIM
      write(StdOut,*)'FLAG_DIM = ',FLAG_DIM
      call CouplerIO_Abort('CouplerIO_WriteHeader:: FlagDim /= FLAG_DIM')
   endif
   if(NumFlagsIn(2)<1)then
      call CouplerIO_Abort('CouplerIO_WriteHeader:: Number of ElemFlags must be greater than 0')
   endif

   NumFlags(:)=NumFlagsIn(:)

   open(UnitCouple,file='elements.bin',form='unformatted')

   if (use_tempIn.eq.0) then
    write(UnitCouple)NumNodes,NumIntElems,IntElemDim,NumFlags
   else if (use_tempIn.eq.1) then
    write(UnitCouple)NumNodes,NumIntElems,IntElemDim,NumFlags,use_tempIn
   else
    call CouplerIO_Abort('CouplerIO_WriteHeader:: use_tempIn invalid')
   endif

   !NumNodes : number of interface nodes
   !NumIntElems : number of interface elements
   !idynamic : coupling mode (0=fixed body,1= presribed body motion,2 deformin body)
   !NumDim : number of dimension (i.e. 1-D, 2-D or 3-D)
   !IntElemDim : Maximum number of nodes in surface element

   if(NumFlags(1)>0)then
      write(UnitCouple)(HdrFlag(k),k=1,NumFlags(1))
   endif

   write(UnitCouple)(Eul2LagNode(k),k=1,NumNodes)
   !Eul2LagNode(k) : Lagrange element number for surface element number k

   do k=1,NumIntElems
      write(UnitCouple)(ElemData(j,k),j=1,NumFlags(2)),(IntElem(j,k),j=1,ElemData(1,k))
      ! NodesinElem : number of nodes in element k
      ! IntElem(k,n) : node number k of element n   
   enddo   

   close(UnitCouple)

   open(UnitCouple,file='nodeok.mpf')
   write(UnitCouple,*)'echo nodes written'
   close(UnitCouple)

   if(Debug)write(UnitDebug,*)' BODY: HEADER MESSAGE written'

   return
end subroutine CouplerIO_WriteHeader

! +++++++++++++++++++++++++++++ CouplerIO_ReadHeader +++++++++++++++++++++++++++++++++++++
subroutine CouplerIO_ReadHeader(Debug,FlagDim,IntElemDim,UnitcoupleIn,UnitOutIn,UnitDebugIn, &
                                StdOutIn,NumNodesOut,NumIntElemsOut,NumFlagsOut,HdrFlag, &
                                Eul2LagNode,ElemData,IntElem,use_tempIn)
!  Called by Fluid Code: 
!     reads surface element connectivity written by Finite Element code at start of run 

!  ----- arguments -----
   logical(4), intent(in) :: Debug
   integer(4), intent(in) :: FlagDim
   integer(4), intent(out) :: IntElemDim 
   integer(4), intent(in) :: UnitcoupleIn,UnitOutIn
   integer(4), intent(in) :: UnitDebugIn,StdOutIn,use_tempIn
   integer(4), intent(out) :: NumNodesOut,NumIntElemsOut
   integer(4), dimension(FlagDim), intent(inout) :: NumFlagsOut                  ! 0 static BC, 1 dynamic BC
   integer(4), dimension(:), pointer :: HdrFlag
   integer(4), dimension(:), pointer :: Eul2LagNode
   integer(4), dimension(:,:), pointer :: ElemData
   integer(4), dimension(:,:), pointer :: IntElem

!  ----- Local variables -----
   logical(4) :: exists
   integer(4) :: j,k,idum,LagDim,error,ElemDim
   integer(4) :: NumTriElems     ! FE code value of number of elements
   character(40) :: filename
   integer(4) :: use_tempOut
 
   UnitCouple=UnitcoupleIn
   UnitOut=UnitOutIn
   UnitDebug=UnitDebugIn
   StdOut=StdOutIn
   use_temp=use_tempIn

   if(FLAG_DIM/=FlagDim)then
      write(UnitOut,*)'FLAG_DIM = ',FLAG_DIM
      write(StdOut,*)'FLAG_DIM = ',FLAG_DIM
      call CouplerIO_Abort('CouplerIO_ReadHeader:: FlagDim /= FLAG_DIM')
   endif
            
   !     +++++ Read Header +++++
   write(STDOUT,*)
   write(STDOUT,*)'+++ Coupled Run; Waiting for Lagrange Data +++'
   if(Debug)write(UnitDebug,*)' FLUID: waiting for HEADER MESSAGE'   

   call WaitOnFile('nodeok.mpf','lag_stop.mpf')
   open(Unitcouple,file='elements.bin',form='unformatted')

   !    +++++ Header Specs +++++
   if (use_tempIn.eq.0) then
    read(UnitCouple)NumNodesOut,NumIntElemsOut,IntElemDim,NumFlagsOut
    use_tempOut=0
   else if (use_tempIn.eq.1) then
    read(UnitCouple)NumNodesOut,NumIntElemsOut,IntElemDim,NumFlagsOut, &
     use_tempOut
   else
    call CouplerIO_Abort('use_tempIn invalid')
   endif
   if (use_tempIn.ne.use_tempOut) then
    call CouplerIO_Abort('use_tempOut<>use_tempIn')
   endif 

   !NumNodesOut : number of interface nodes
   !NumIntElemsOut : number of interface elements
   !IntElemDim : Maximum number of nodes in surface element
   !NumFlagsOut : Number of message flags

   NumNodes=NumNodesOut   
   NumIntElems=NumIntElemsOut
   NumFlags(:)=NumFlagsOut(:) 
   
   if(NumFlags(1)>0)then
      allocate(HdrFlag(NumFlags(1)),STAT=error)
      if(error/=0)then
         call CouplerIO_Abort('CouplerIO_ElementHeader:: HdrFlag allocation failed')
      endif
      read(UnitCouple)(HdrFlag(k),k=1,NumFlags(1))
   endif  

   !    +++++ Allocate Storage +++++
   allocate(ElemData(NumFlags(2),NumIntElems),IntElem(IntElemDim,NumIntElems), &
            Eul2LagNode(NumNodes),STAT=error)
   if(error/=0)then
      call CouplerIO_Abort('CouplerIO_ElementHeader:: allocation failed')
   endif

   read(UnitCouple)(Eul2LagNode(k),k=1,NumNodes)
   !Eul2LagNode(k) : Lagrange element number for surface element number k

   do k=1,NumIntElems
   !  for each element read
      read(UnitCouple)(ElemData(j,k),j=1,NumFlags(2)),(IntElem(j,k),j=1,ElemData(1,k))
      ! NodesinElem : number of nodes in element k
      ! DWIIntElem : True if element is doubly wetted; false otherwise
      ! IntElemPartNo : depends on Euler code needs. Currently assigns elements to distinct parts
      ! IntElem(k,n) : node number k of element n
   enddo

   close(UnitCouple)
   write(STDOUT,*)'Lagrange Header and Elements Read'

   call SafeOpen(UnitCouple,'nodeok.mpf','lag_stop.mpf')
   close(UnitCouple,status='delete')

   if(Debug)write(UnitDebug,*)' FLUID: HEADER MESSAGE read'   
         
   return
  
end subroutine CouplerIO_ReadHeader

! +++++++++++++++++++++++++++++ CouplerIO_WriteNodes +++++++++++++++++++++++++++++++++++++
subroutine CouplerIO_WriteNodes(Debug,istep,time,dt,NptFlag,Node,NodeVel,ActiveIntElem,NodeTemp)
!  Called by Structure Code: 
!     write node locations and velocities

!  ----- arguments -----
   logical(4), intent(in) :: Debug
   integer(4), intent(in) :: istep
   real(8), intent(in) :: time
   real(8), intent(in) :: dt

   integer(4), dimension(:), pointer :: NptFlag
   real(8), dimension(3,NumNodes):: Node,NodeVel
   real(8), dimension(NumNodes):: NodeTemp
   logical(4), dimension(NumIntElems) :: ActiveIntElem

!  ----- local variables -----
   integer(4) :: k,j

   if(Debug)write(UnitDebug,*)'BODY: Waiting for NODE MESSAGE to be read'
   call WaitRemoveFile('lag_ok.mpf','eul_stop.mpf')

!  read file
   open(UnitCouple,file='lag_out.bin',form='unformatted')

   write(UnitCouple)istep,NumNodes,NumIntElems,time,dt
   ! istep : step  number
   ! NumNodes : number of nodes
   ! NumIntElems : Number of Lagrange elements
   ! time : time
   ! dt : time step

   if(NumFlags(3)>0)then
      write(UnitCouple)(NptFlag(k),k=1,NumFlags(3))
   endif
    
   !  For each Lagrange element write whether it is active or not (i.e. whether it has failed)
   write(UnitCouple)(ActiveIntElem(k),k=1,NumIntElems)

   do k=1,NumNodes
   !  For each node write
      if (use_temp.eq.0) then
       write(UnitCouple)Node(1,k),Node(2,k),Node(3,k),NodeVel(1,k),NodeVel(2,k),NodeVel(3,k)
      else if (use_temp.eq.1) then
       write(UnitCouple)Node(1,k),Node(2,k),Node(3,k),NodeVel(1,k),NodeVel(2,k),NodeVel(3,k),NodeTemp(k)
      else 
       call CouplerIO_Abort('CouplerIO_WriteNodes:: use_temp invalid')
      endif
 

      ! Node(:,k) - node locations
      ! NodeVel(:,k) ; Node velocities
   enddo

   close(UnitCouple)

   call SafeOpen(UnitCouple,'lag_ok.mpf','eul_stop.mpf')
   write(UnitCouple,*)'lag_out.mpf done'
   close(UnitCouple)
   if(Debug)write(UnitDebug,'(a,i6,e16.8)')'BODY: NODE MESSAGE read, Body step,time =',istep,time


   return
end subroutine CouplerIO_WriteNodes

! +++++++++++++++++++++++++++++ CouplerIO_ReadNodes +++++++++++++++++++++++++++++++++++++
subroutine CouplerIO_ReadNodes(Debug,step_lag,time_lag,dt_lag,NptFlag,Node,NodeVel,ActiveIntElem,NodeTemp)
!  Called by Fluid Code: 
!     write node locations and velocities

!  ----- arguments -----
   logical(4), intent(in) :: Debug
   integer(4), intent(out) :: step_lag
   real(8), intent(out) :: dt_lag,time_lag
   integer(4), dimension(:), pointer :: NptFlag
   logical(4), dimension(NumIntElems), intent(out) :: ActiveIntElem
   real(8), dimension(3,NumNodes), intent(out) :: Node,NodeVel
   real(8), dimension(NumNodes), intent(out) :: NodeTemp

!  ----- local variables -----
   logical(4) :: exists
   real(8), parameter :: EPS=1.e-12
   real(8) :: delta_time
   integer(4) :: k,n,nstep,nodes_in,isopen,elem_in

   if(Debug)write(UnitDebug,'(a)')'FLUID:  Waiting for NODE MESSAGE'
   call WaitOnFile('lag_ok.mpf','lag_stop.mpf')      

!  Structure code still running
   open(UnitCouple,FILE='lag_out.bin',form='unformatted',IOSTAT=isopen)

   read(UnitCouple)step_lag,nodes_in,elem_in,time_lag,dt_lag
   ! istep_lag : step  number
   ! nodes_in : number of nodes
   ! elem_in : Number of Lagrange elements
   ! time_lag : Lagrange code time
   ! dt_lag : Lagrange code time step 

   print *,"in ReadNode step_lag,nodesin,elemin,timelag,dtlag "
   print *,step_lag,nodes_in,elem_in,time_lag,dt_lag
   print *,"Debug ",Debug
   print *,"numflag1-4 ",NumFlags(1),NumFlags(2),NumFlags(3),NumFlags(4)
 
   if(NumFlags(3)>0)then
      read(UnitCouple)(NptFlag(k),k=1,NumFlags(3))
   endif

   if(nodes_in /= NumNodes)then
      print *,"nodes_in,NumNodes ",nodes_in,NumNodes
      call CouplerIO_Abort('CouplerIO_ReadNodes :: nodes_in/=NumNodes')
   endif

   if(Elem_in /= NumIntElems)then
      call CouplerIO_Abort('CouplerIO_ReadNodes :: elem_in/=NumIntElems')
   endif

   !  For Lagrange element read whether it is active or not (i.e. whether it has failed)
   read(UnitCouple)(ActiveIntElem(n),n=1,elem_in)

   do n=1,nodes_in
   !  For each node read
      if (use_temp.eq.0) then
       read(UnitCouple)Node(1,n),Node(2,n),Node(3,n),NodeVel(1,n),NodeVel(2,n),NodeVel(3,n)
       NodeTemp(n)=0.0
      else if (use_temp.eq.1) then
       read(UnitCouple)Node(1,n),Node(2,n),Node(3,n),NodeVel(1,n),NodeVel(2,n),NodeVel(3,n),NodeTemp(n)
      else 
       call CouplerIO_Abort('CouplerIO_ReadNodes:: use_temp invalid')
      endif

      ! Node(:,k) - node locations
      ! NodeVel(:,k) ; Node velocities
   enddo
   close(UnitCouple)

   call SafeOpen(UnitCouple,'lag_ok.mpf','lagstop.mpf')
   close(UnitCouple,status='delete')

   if(Debug)write(UnitDebug,'(a,i6,e16.8)')'FLUID: NODE MESSAGE read, Body step,time =',step_lag,time_lag

   return
end subroutine CouplerIO_ReadNodes

! +++++++++++++++++++++++++++++ CouplerIO_WriteLoads +++++++++++++++++++++++++++++++++++++
subroutine CouplerIO_WriteForces(Debug,step,nlast,time,dt,ForFlag,NodeForce)
!  Called by Fluid Code: 
!     write node forces
   logical(4), intent(in) :: Debug
   integer(4), intent(in) :: step
   integer(4), intent(in) :: nlast
   real(8), intent(in) :: time,dt
   integer(4),dimension(:), pointer:: ForFlag
   real(8), dimension(3,NumNodes), intent(in) :: NodeForce

   logical(4) :: exists
   real(8) :: deltatime,timeinit
   integer(4) :: k,n,isopen,i

   !  target time reached, write node loads
   if(Debug)write(UnitDebug,*)'FLUID: Waiting for old FORCE MESSAGE to be read'

   !  Waiting for FE code to read eul_out.bin and remove eul_ok.mpf file
   Call WaitRemoveFile('eul_ok.mpf','lag_stop.mpf')

   ! write node loads file, eul_out.bin
   open(UnitCouple,FILE='eul_out.bin',form='unformatted',IOSTAT=isopen)

   write(UnitCouple)step,NumNodes,time,dt
   ! step - step number
   ! NumNodes - number of nodes
   ! time - Fluid code time
   ! dt - Fluid code step size

   if(NumFlags(4)>0)then
      write(UnitCouple)(ForFlag(k),k=1,NumFlags(4))
   endif
    
   do n=1,NumNodes
   !  for every node write
      write(UnitCouple)(NodeForce(i,n),i=1,3)
      ! NodeForce(:,n) - node n forces 
   enddo
   close(UnitCouple)

   ! write Euler OK signal file eul_out.mpf
   call SafeOpen(UnitCouple,'eul_ok.mpf','lag_stop.mpf')
   write(UnitCouple,*)'EUL: cycle done'
   close(UnitCouple)
   
   if(Debug)write(UnitDebug,'(a,i6,e16.8)')'FLUID: FORCE MESSAGE written, FLUID step time ',step,time

   return
end subroutine CouplerIO_WriteForces

! +++++++++++++++++++++++++++++ CouplerIO_ReadForces +++++++++++++++++++++++++++++++++++++
subroutine CouplerIO_ReadForces(Debug,Step_Eul,time_Eul,dt_Eul,ForFlag,NodeForce)
!  Called by Structure Code: 
!     read node forces

!  ----- arguments -----
   logical(4), intent(in) :: Debug
   integer(4), intent(out) :: Step_Eul   
   real(8), intent(out) :: time_Eul,dt_Eul
   integer(4),dimension(:), pointer:: ForFlag
   real(8), dimension(3,NumNodes) :: NodeForce 

!  ----- local variables   
   LOGICAL :: exists
   integer(4) :: NumEulNodes,k,j

   if(Debug)write(UnitDebug,'(a)')'BODY:  Waiting for FORCE MESSAGE'
   call WaitOnFile('eul_ok.mpf','eul_stop.mpf')

   !        read and delete euler file
   open(UnitCouple,file='eul_out.bin',form='unformatted')

   read(UnitCouple)Step_Eul,NumEulNodes,time_Eul,dt_Eul
   ! Step_Eul - step number
   ! NumEulNodes - number of nodes
   ! time_Eul - Fluid code time
   ! dt_Eul - Fluid code step size

   if(NumFlags(4)>0)then
      read(UnitCouple)(ForFlag(k),k=1,NumFlags(4))
   endif

   do k=1,NumNodes
   !  for every node write
      read(UnitCouple)(NodeForce(j,k),j=1,3)
      ! NodeForce(:,k) - node k forces 

   enddo

   close(UnitCouple,status='delete')

   !        delete flag file
   call SafeOpen(UnitCouple,'eul_ok.mpf','eul_stop.mpf')
   close(UnitCouple,status='delete')

   if(Debug)write(UnitDebug,'(a,i6,e16.8)')'BODY: FORCE MESSAGE read, FLUID step time ',Step_Eul,time_Eul

   return    
end subroutine CouplerIO_ReadForces

! +++++++++++++++++++++++++++++ SafeOpen +++++++++++++++++++++++++++++++++++++
subroutine SafeOpen(filenum,filename,stopfile)
!recover for file open error and try again

!   ----- arguments -----
   integer(4), intent(in) :: filenum        ! unit number of file to be opened
   character(*), intent(in) :: filename    ! name of file to be opened
   character(*), intent(in) :: stopfile    ! FE code stopped file


!  ----- local variables -----
   logical(4) :: opened,done
   integer(4) :: ict

   ict=0
   done=.false.
   opened=.false.
   do while ((ict.le.MAXLOOP).AND.(.NOT.done).AND.(.NOT.opened))
      ict=ict+1
      open(filenum,file=filename,err=100)
      opened=.true.
      cycle

100   continue
      if(ict.eq.1)write(UnitOut,*)'open waiting for file ',filename
!         call sleepqq(waittime)  (commented previously)
      inquire(file=stopfile,exist=done)
   enddo

   if(.not.opened)then
      call CouplerIO_Abort('Coupler SafeOpen:: timed out waiting for file '//filename)
   else if(done)then
      call CouplerIO_Abort('Coupler SafeOpen:: Fluid Stop File Found; file name is  '//stopfile)
   endif

   return
end subroutine SafeOpen

!++++++++++++++++++++++++++++ WaitOnFile +++++++++++++++++++++++++++++++++ 
subroutine WaitOnFile(filename,stopfile)
!  Wait on File to appear

!   ----- arguments -----
   character(*), intent(in) :: filename    ! name of file to be opened
   character(*), intent(in) :: stopfile    ! FE code stopped file

!   ----- local variable -----
   logical(4) :: here,done
   integer(4) :: ict,iend
   character(len=80) :: line

   inquire(file=filename,exist=here)

   ict=0
   done=.false.
   do while((.NOT.here).AND.(ict.le.MAXLOOP).AND.(.NOT.done))
      ict=ict+1
      inquire(file=filename,exist=here)
      call sleep(waittime)   ! unix
      inquire(file=stopfile,exist=done)
   enddo

   if(.not.here)then
      call CouplerIO_Abort('CouplerIO WaitOnFile:: timed out waiting for file '//filename)
   else if(done)then
      call CouplerIO_Abort('CouplerIO WaitOnFile::  Stop File Found; file name is '//stopfile)
   endif

   return
end subroutine WaitOnFile

! +++++++++++++++++++++++++++++ WaitRemoveFile +++++++++++++++++++++++++++++++++++++
subroutine WaitRemoveFile(filename,stopfile)
!   Wait for file to be removed

!   ----- arguments -----
   character(*), intent(in) :: filename    ! name of file to be opened
   character(*), intent(in) :: stopfile    ! FE code stopped file

!   ----- local variable -----
   logical(4) :: here,done
   character(len=80) :: line
   integer(4) :: ict

   inquire(file=filename,exist=here)

   ict=0
   done=.false.
   do while((here).AND.(ict.le.MAXLOOP).AND.(.NOT.done))

      ict=ict+1
      inquire(file=filename,exist=here)
      call sleep(waittime)   ! unix
      inquire(file=stopfile,exist=done)

   enddo

   if(here)then
      call CouplerIO_Abort('CouplerIO WaitRemoveFile:: timed out waiting for file '//filename)
   else if(done)then
      call CouplerIO_Abort('CouplerIO WaitRemoveFile::  Stop File Found; file name is '//stopfile)
   endif

   return

end subroutine WaitRemoveFile

subroutine CouplerIO_Abort(string)
   character(*), intent(in) :: string

   write(StdOut,'(a)')string
   write(UnitOut,'(a)')string
   stop 'CouplerIO_Abort'
      
end subroutine CouplerIO_Abort


end module CouplerIO
