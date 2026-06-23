
!-----------------------------------------------
!   Does parallel universe Grover searches for !
!   different noise scenarios                  !
!                                              !
!   Frank Tabakin & Bruno Julia Diaz           !
!                                              !
!       Aug 30th, 2008                         !
!----------------------------------------------------------------------------
!   Universe is split into Multiuniverses(parallel universes)               |
!   Universe consists of all processors:  nprocU                            |
!   Universe processors are labelled by myidU =0 ...nprocU-1                |
!   The nprocU Processors are separated into NGROUPS                        |
!   The overall master processor is masterU=0                               |
!   The overall communicator is commU                                       |
!                                                                           | 
!   Each Multiuniverse consists of nprocM=nprocU/NGROUP processors          |
!   Within each Multiuniverse the processors are labelled by myidM          |
!   myidM=0, nprocM-1                                                       |
!   The Multiuniverses are labelled by rankM                                |
!   The master processor is thin each Multiuniverse is masterM=0            |
!   The communicator for each Multiuniverse is commM                        |
!----------------------------------------------------------------------------

Program main
Use mpi
Use mpi_vars
Use qcmpisubs
Implicit None

DOUBLE PRECISION RAND(2)

INTEGER masterU,commU
INTEGER masterM,commM
INTEGER rankM
INTEGER NGROUP
INTEGER nq,NBITS,NP,NG,IR,I,IC

masterU=0
commU=MPI_comm_world
!
!initialise MPI universe
!
  call MPI_init(ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_init'

!get size (total number of processors) within MPI Universe

  call MPI_comm_size(commU,nprocU,ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_comm_size'

!get my myidU on each different processor
!within MPI universe (myidU=0 to nprocU-1)

  call MPI_comm_rank(commU,myidU,ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_comm_rank' 

!-------------------------------------------------------
!  Use INPUT.txt for nq, NPROCU, NGROUP test if set OK |
!  Need to set NP=NPROCU=2^p, with p=even              |
!  Need to set NG=NGROUP=2^g, with g=0,1,2...p         |
!  Then prow=pcol= Sqrt(NP) on 2D processor grid       |
!  and prow X pcol = NP=NPROCU                         |
!-------------------------------------------------------

     open(unit=1,file="Input.txt",status="old",form="formatted")
     read(1,*)nq
     read(1,*)NP
     read(1,*)NG

     close(1)

     if (1.eq.1) then
       print *,"nq=",nq
       print *,"NP=",NP
       print *,"NP=",NG
     endif

     if (NP.ne.nprocU.or.NG.gt.NP ) then
        print *,"Values for NP and/or NG not selected correctly"
        print *,"Need to fix mpirun -np NP  setting"
        stop 
     endif
!------------------------------------------
! Number of Qubits nq from input yields: 
!------------------------------------------  	
  NBITS = 2**nq 

!-----------------------------------------|
!  Get random location for marked item    |
!-----------------------------------------|

if (myidU.eq.masterU) then

  CALL RANDX(RAND,1)
  IR=INT(RAND(1)*NBITS)

endif

!
!Broadcast "IR"  the marked item to all processes 
!

call MPI_Bcast(IR, 1, MPI_INTEGER, 0, commU,ierr)
if (ierr/=0) stop 'Error in MPI_Bcast of IR'

!--------------------------------------|
! Stipulate Number groups              |
!--------------------------------------|

   NGROUP=NG

   nprocM=nprocU/NGROUP

   if (nprocM.EQ.0) THEN 
    WRITE(*,*) 'ERROR, nprocU SMALLER THAN NGROUP'
   endif

!
!  Defines which processors belong to which group
!  according to "rankM" rankM labels the Multiuniverse.
!  rankM=0 is the first group of processors: 
!  Multiuniverse 0
!
!  rankM=NGROUP-1 is the last group of processors: 
!  Multiuniverse NGROUP-1
!
rankM=Int(myidU/NPROCM)

!
! creates group(MULTIUNIVERSE) communicators commM
!

call MPI_COMM_SPLIT(commU,rankM,myidU,commM,ierr)
if (ierr/=MPI_success) stop 'Error in MPI_COMM_split'


call Grover(commM,nq,rankM,IR)
  
! finalize MPI

  call MPI_finalize(ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_finalize'
contains
!  subgrover.f90                                                  |
!                                                                 |
!  MPI multi-universe version for nq qubits                       |
!                                                                 |
! This version distributes Psi over nprocM                        |
! and distributes Rho over NPRHO processors                       |
! to allow parallel treatment of density matrix                   |
!                                                                 |
!       Aug 11, 2008                                              |
!                                                                 |
!   Bruno Julia Diaz & Frank Tabakin                              |
!   Version usin HALL                                             |
!------------------------------------------------------------------
!                                                                 |
! Output:                                                         |
!      fort.10    Master of All groups                            |
!                                                                 |
!      fort.10+rankM*nprocM+myidM   each processor output         |
!                                                                 |
!------------------------------------------------------------------

Subroutine Grover(commM,nq,rankM,IR)

Use MPI
Use MPI_VARS

Implicit None

DOUBLE PRECISION RAND(1:2,1:1) 
COMPLEX, PARAMETER   :: One  =( 1.0, 0.0 )
COMPLEX, PARAMETER   :: Zero  =( 0.0, 0.0 )

COMPLEX, ALLOCATABLE :: Psi1(:)
COMPLEX, ALLOCATABLE :: Psi2(:) 
COMPLEX, ALLOCATABLE :: Psi3(:) 
REAL,    ALLOCATABLE :: R(:)

REAL const, NORM,FNORM,isend(2), irecv(2),probstop 
REAL S

DOUBLE PRECISION starttime, endtime,midtime

INTEGER nq,NBITS,IR
INTEGER i, j,x, js, ntry,section,seat,seatn
INTEGER eloc
INTEGER n, np, npl, NPART,NPRHO
INTEGER RESULT
INTEGER OUTFILE
!
! Multiuniverse variables
!
INTEGER masterU,commU
INTEGER masterM,commM
INTEGER rankM
INTEGER TGROUPS

COMPLEX temp,getsum,getsum2,getsum3
!
! MPI declaration and timing
!
INTEGER dest

INTEGER n1tag,MPISTAT(mpi_status_size)
COMPLEX tracerho,tracerho2,tracerho3
NBITS=2**nq

!
! MPI initialization
!

commU=MPI_comm_world

call MPI_COMM_SIZE(commU,nprocU,ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_comm_size'
call MPI_comm_rank(commU,myidU,ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_comm_rank'
!
! Now for the multiUniverse within the Universe
!
call MPI_COMM_SIZE( commM, nprocM, ierr) 
call MPI_COMM_RANK( commM, myidM, ierr)
!
!  masterU=master of all universe
!  masterM is master of one Multiuniverse within Universe.
!
masterU=0
masterM=0
!
!start timing
!
if (myidM.eq.masterM) starttime = MPI_WTIME()

!
! Outfile
!
OUTFILE=10+rankM*nprocM+myidM

TGROUPS=nprocU/nprocM
!
! Memory allocation for distributed vectors and arrays.
!
NPART=NBITS/nprocM

Allocate(psi1(NPART))
Allocate(psi2(NPART))
Allocate(R(NPART))

! Needs to also allocate these
! for Rho construction and distribution.


!
! Open output files
!

OPEN(unit=OUTFILE)

write(OUTFILE,*) '-----------------------------------------------' 
write(OUTFILE,*) '|    This program is for Grover               |' 
write(OUTFILE,*) '|    Searching Algorithm-Single Marked Item   |'
write(OUTFILE,*) '|    State Version  with Amplitude mapping.   |' 
write(OUTFILE,*) '-----------------------------------------------' 
if (myidU.eq.masterU) then
write(OUTFILE,*)  
write(OUTFILE,*)  '          MASTER of ALL GROUPS'
write(OUTFILE,*)  
write(OUTFILE,8999)  nprocU,nprocU/nprocM
8999 Format('     Total N. processors :',I4, ' Groups :',I4)
write(OUTFILE,*)  
write(OUTFILE,*) '-----------------------------------------------' 
endif


if (myidM.eq.masterM) then
write(OUTFILE,*)  
write(OUTFILE,9000)  rankM,TGROUPS
9000 Format('           MASTER  of Group Number=',I4,' out of ',I4,' Groups')
write(OUTFILE,*)  
write(OUTFILE,9001)  nprocM
9001 Format('     N. processors     :',I4)
write(OUTFILE,*)  
write(OUTFILE,*) '-----------------------------------------------' 
else



write(OUTFILE,*) ' '
write(OUTFILE,9002)  myidM,rankM,TGROUPS
9002 Format('     Output: processor=',I4,' Group=',I4,' out of ',I4, ' Groups')
endif

write(OUTFILE,*)
write(OUTFILE,9003) nq
9003 Format('     Qubits           :   ', I4)
write(OUTFILE,9004) NPART
9004 Format('     W.f. array size  : ', I8)


!------------------------------------------------------------------------
! Following is done for all processors                                  |
!                                                                       |               
! Grover Searching Algorithm-Single Marked Item.                        |
! Grover's search algorithm making use of states in subspace |x>.       |
! The space |y> is part of the search, but keeps the same state during  |
! the process and can therefore be suppressed.                          |
!                                                                       |
! Form  Psi1= |++++ ...> = const* Sum_n |n>                             |
! for the initialization of the nq qubits.                              |
!                                                                       |
! Note  the wavefunction is distributed over the processors             |
! to save space.                                                        |
! No need to construct Hadamard matrix HALL at this point.              |
! This version does first Hall directly.                                |
!------------------------------------------------------------------------


!
! Initial wave function
!
const= 1./2**(nq/2)
Do 20 i=1,NPART
Psi1(i)=const*One
20 Continue

Norm=0.
Do 18 i=1,NPART 
Norm=abs(Psi1(i)*conjg(Psi1(i))) + Norm
18 continue

call MPI_Reduce(NORM,FNORM,1,MPI_REAL,MPI_SUM,0, commM ,ierr)
NORM=FNORM

If(myidM.eq.masterM)Write(OUTFILE,9005) Norm
9005 format('     Normalization    :',F8.3)

!------------------------------------------------------------------------
!  IR= integer to mark entry in OracleG                                 |
!  obtain for myidU=0 and then BCAST it to all.                          |
!  Determine                                                            |
!  section/seat for where it is located in Psi1 distribution            |
!  Section==which processor; seat=;location within that processor       |
!-----------------------------------------------------------------------|
!  We find the section(processor)  and seat(location within processor)
!  for the marked item.  
!
  section= Int(IR/NPART) 
  seat= Mod(IR,NPART)
!
! Introduce noise at specific locations eloc
! We use three possible locations. First universe 
! is noiseless.
!
! Noise location depends on the universe at stake
!
If(rankM.gt.0.and.myidM.eq.masterM) then
CALL RANDX(RAND,1)
eloc =   mod(int(3.*RAND(2,1)+rankM)+1,3)+1
!write(*,*) rand(2,1),eloc

write(OUTFILE,*)  
write(OUTFILE,*) '-----------------------------------------------' 
write(OUTFILE,*)' Group subject to noise=',rankM
write(OUTFILE,*)' At location (1,2,3) ',eloc
endif
!
!Place noise at location eloc=1
!
If(rankM.gt.0.and.eloc.eq.1.and.myidM.eq.masterM) then	
write(OUTFILE,*)'             ' 
write(OUTFILE,*)' --->   Noise before building the Grover operator' 
call Noise(nq,rankM,psi1,NPART,commM)
endif
!  Continue Process
!------------------------------------------------------------------------
!Step 2: Building the Grover operator                                   !
!                                                                       !
!   Prepare the Grover operator G as the "inversion about the mean"     !
!   times the Oracle                                                    !
!                                                                       !
!   The Oracle defines an arbitrary single marked state                 !
!   When acting in just the |x> space, the Grover Oracle G changes the  !
!   sign of the marked state.                                           !
!                                                                       !
!   The Grover process, G, must be applied a number of times which      !
!   is estimated here to be ntry.                                       !
!-----------------------------------------------------------------------!
!
! Place noise at location eloc=2
!
If(rankM.gt.0.and.eloc.eq.2.and.myidM.eq.masterM) then
write(OUTFILE,*)'               '          
write(OUTFILE,*)'  --->  Noise after building the Grover operator'          
call Noise(nq,rankM,psi1,NPART,commM)
endif
ntry= 2**(nq/2) -1
ntry= (3.14/4)* (2**(nq/2))
! LIMIT MAXIMUM NTRY TO COMPLETE FOR LARGE NQ  Oct,2008
If(ntry.gt.200) then
Write(*,*)'NTRY WAS SET TO MAXIMUM OF 200 versus',ntry
ntry=200
endif
!  Cut down ntry to speed up code  7/30 ft
!  ntry=ntry/5
!  Do few` cycles  to test 8/2/08   ft
!ntry=3    

Write(OUTFILE,9006) ntry,section,seat
9006 Format('     Ntry, Section, Seat       : ',I4,2x,I6,2x,I6) 

Do 45 j=1,ntry

!
! Do Oracle action on psi
! Need to reverse Psi1 for marked item in right section /seat
!
If(myidM.eq.section)Psi1(seat+1)= - Psi1(seat+1) 

!------------------------------------------------------
! Do Hall on psi                                      !
! This version is for MPI distribution of Psi1 case   !
!------------------------------------------------------
call HALL( nq,psi1,NPART,commM )

!-------------------------
! END First Hall on psi  |
!-------------------------
!
!--------------------------------------------------------
! Step 3: Inversion                                     |
!                                                       |
!  Inversion step ( 2 |00000...><00000... | - 1 )       |
!--------------------------------------------------------

Do 50 i=1, NPART
If(myidM.ne.masterM) Psi1(i)= - Psi1(i)
If(i.ne.1.and.myidM.eq.masterM) Psi1(i)= - Psi1(i)
50 continue

!--------------------------------------------|
!  Do Hall on psi again                      |
!--------------------------------------------|

call HALL( nq,psi1,NPART,commM )
!--------------------------------------------|
!  End Hall on psi again                     |
!--------------------------------------------|


!--------------------------------------------|
! Result after Ntry                          |
!--------------------------------------------|

Do 53 i=1, NPART
R(i)= abs(psi1(i)*conjg(psi1(i)))
53 Continue

!
! Each group master writes its timing record
! For each ntry cycle
!

if (myidM.eq.masterM)then 
midtime = MPI_WTIME()
Write(OUTFILE,9411)j, midtime-starttime
9411 Format('ntry number',1I4,' Time .........',F9.4, ' seconds')
endif

45 continue  

!------------------------------------
! END of Grover action, ntry times  |
!------------------------------------

!
! Place noise at location eloc=3
!
If(rankM.gt.0.and.eloc.eq.3.and.myidM.eq.masterM) then
write(OUTFILE,*)'              ' 
write(OUTFILE,*)' --->  Noise after Grover action' 
call Noise(nq,rankM,psi1,NPART,commM)
endif


!
!----------------------
! Look at final state |
!----------------------

Do 54 i=1,NPART
If(R(i).eq.maxval(R)) js=i-1
54 Continue

isend(1) = maxval(R)
isend(2) = js+  myidM*NPART

write(OUTFILE,*) 
write(OUTFILE,*) '    Max Prob located :', js
write(OUTFILE,9007)  maxval(R)
9007 format('     Max Probability  :',F9.4)


!---------------------------------------------------------|
! Finds maximum probability of all registers              |
!---------------------------------------------------------|

CALL MPI_REDUCE(ISEND,IRECV,1,MPI_2REAL,MPI_MAXLOC,masterM,commM,IERR)
if (ierr.ne.0) write(OUTFILE,*) 'Error in reducing the maximum prob, myidM=',myidM

!---------------------------|
! Final output from master  |
!---------------------------|
IF (myidM == masterM) THEN

RESULT=irecv(2)

write(OUTFILE,*) ''
write(OUTFILE,*) '-----------------------------------------------' 
write(OUTFILE,*) '     Result of search in Group                 '
write(OUTFILE,*) '-----------------------------------------------' 
write(OUTFILE,*) ''
9008 Format('            Maxloc=', I8,'  /  Max Prob= ',F8.6 )
write(OUTFILE,9008) RESULT,irecv(1)
write(OUTFILE,*)
write(OUTFILE,9009)
9009 Format('               --  COMPARE --           ')
write(OUTFILE,*)
write(OUTFILE,9010) IR
9010 Format('     Actual Value:IR=', I8 )
Write(OUTFILE,*)
!
! Each group master writes its timing record
!
endtime = MPI_WTIME()

write(OUTFILE,9011) endtime-starttime
9011 Format('      Time .........',F9.4, ' seconds')
write(OUTFILE,*)
write(OUTFILE,*) '-----------------------------------------------' 

ENDIF

close (unit=OUTFILE)


deAllocate(psi1)
deAllocate(psi2)
deAllocate(R)
!
!
End subroutine

end program main
