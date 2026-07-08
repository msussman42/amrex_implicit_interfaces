
!-----------------------------------------------
!   Does parallel universe Shor factorization  !
!   for different noise scenarios              !
!                                              !
!   Frank Tabakin & Bruno Julia Diaz           !
!                                              !
!       May 13, 2008                         !
!-----------------------------------------------
!   Universe is split into Multiuniverses(parallel universes)
!   Universe consists of all processors:  nprocU
!   Universe processors are labelled by myidU =0 ...nprocU-1
!   The nprocU Processors are separated into NGROUPS
!   The overall master processor is masterU=0
!   The overall communicator is commU
!
!   Each Multiuniverse consists of nprocM=nprocU/NGROUP processors
!   Within each Multiuniverse the processors are labelled by myidM
!   myidM=0, nprocM-1
!   The Multiuniverses are labelled by rankM
!   The master processor is thin each Multiuniverse is masterM=0
!   The communicator for each Multiuniverse is commM
!
Program main
#ifdef USE_MPI
Use mpi
#endif
Use mpi_vars
Use qcmpisubs
Implicit None

DOUBLE PRECISION RAND(2)

INTEGER masterU,commU
INTEGER masterM,commM
INTEGER rankM
INTEGER  NGROUP
INTEGER nq,NBITS
INTEGER IR

INTEGER Q,N1,N2,OUTFILE

masterU=0
#ifdef USE_MPI
commU=MPI_comm_world
#endif

!
!initialise MPI universe
!

  nprocU=1
#ifdef USE_MPI
  call MPI_init(ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_init'

!get size (total number of processors) within MPI Universe

  call MPI_comm_size(commU,nprocU,ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_comm_size'
#endif
  if (1.eq.1) then
   print *,"nprocU=",nprocU
  endif

!get my myidU on each different processor
!within MPI universe (myidU=0 to nprocU-1)

  myidU=0

#ifdef USE_MPI
  call MPI_comm_rank(commU,myidU,ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_comm_rank' 
#endif


OUTFILE=10

!
! Number of Qubits
!   MAX ALLOWED seems to be 12 at PITT	
!  nq    =  12
!  NBITS = 2**nq 


!-----------------------------------------|
!  Number to be factorized by all groups
!-----------------------------------------|

!  IR=7*11
!  IR=7*13 !OK
  IR=7*17 !OK
!  IR=7*19 !133

!---------------------------------------------------------
!  We need the minimum number of qubits for 
!  both registers. Here we should be careful, 
!  a high n1 number will lengthen the computation.  
!  By default we give estimates for n1 and n2, 
!  based on the fact that n1 should be able to 
!  represent numbers up to 2 IR^2, while n2 should 
!  allow up to M. However, in many cases it is 
!  helpful to lower the value of n1 
!  [n2 is not relevant for the code's performance]
!  and check the resulting probabilities.
!
!  Note divide by Log 2 to convert from base 10 
!  to base 2 Log
!---------------------------------------------------------

 Q = Ceiling(Log(IR**2.)/Log(2.))
 Q=2**Q


!
!  n1 register 1-- determined to accomodate up to number 2 IR^2 
!  n2 register 2--determined to accomodate up to number IR 
!  n1=Ceiling(Log(Real(IR*IR))/Log(2.))
!
 n1=Ceiling(Log(Real(Q))/Log(2.))
 n2=Ceiling(Log(Real(IR))/Log(2.))
 nq=n1+n2
 NBITS = 2**nq 

 if (1.eq.1) then
  print *,"nq=",nq
 endif

if (myidU.eq.masterU) then 
write(OUTFILE,*) '-----------------------------------------------' 
 Write(OUTFILE,*)
 Write(OUTFILE,*) '       Size of registers'
 Write(OUTFILE,*)
 Write(OUTFILE,*) '    Q  = ',Q
Write(OUTFILE,8003) nq
8003 Format('       nq = ', I2, '   qubits in the system.')
Write(OUTFILE,8004) n1 
8004 Format('       n1 = ', I2, '   qubits in reg 1.')
Write(OUTFILE,8005) n2 
8005 Format('       n2 = ', I2, '   qubits in reg 2.')
endif

!--------------------------------------|
! Stipulate Number groups              |
!--------------------------------------|

   NGROUP=1
#ifdef USE_MPI
   NGROUP=4
#endif

   nprocM=nprocU/NGROUP

   if (nprocM.EQ.0) THEN 
    WRITE(*,*) 'ERROR, nprocM SMALLER THAN NGROUP'
   endif

!
!  Defines which processors belong to which group
!  according to "rankM"  rankM labels the Multiuniverse.
!
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
!call MPI_COMM_SPLIT(commU,rankM,nprocU,commM,ierr)
#ifdef USE_MPI
call MPI_COMM_SPLIT(commU,rankM,myidU,commM,ierr)
if (ierr/=MPI_success) stop 'Error in MPI_COMM_split'
#endif

call Shor(commM,n1,n2,rankM,IR)
  
! finalize MPI
#ifdef USE_MPI
  call MPI_finalize(ierr)
  if (ierr/=MPI_success) stop 'Error in MPI_finalize'
#endif

end program main
!contains

!------------------------------------------------------------------
!
!  May 13th, 2008
!
!  State MPI version of Shor's factoring algorithm
!  This version distributes psi over Nprocs
!  Code assumes Nproc= 2^p, where p=0,1,2,3....
!
!  Bruno Julia-Diaz, Frank Tabakin 
!  Note the subroutine CHPASEK  is the controlphase(k) 
!
!------------------------------------------------------------------


SUBROUTINE Shor(commM,n1,n2,rankM,MP)
#ifdef USE_MPI
Use mpi
#endif
Use mpi_vars
Use qcmpisubs
Implicit None

DOUBLE PRECISION RAND(1:2,1:1)

COMPLEX, PARAMETER :: One  =( 1.0, 0.0 )
COMPLEX, PARAMETER :: Zero =( 0.0, 0.0 )


Complex, ALLOCATABLE :: Psi(:)
Complex, ALLOCATABLE :: PsiO(:)
INTEGER, ALLOCATABLE :: PArray(:,:)
Integer, ALLOCATABLE :: B(:)
Integer, ALLOCATABLE :: B1(:)
Integer, ALLOCATABLE :: B2(:)
Integer, ALLOCATABLE :: result(:)

INTEGER i,iv,ivsection,ivseat
INTEGER j,jv,j1,j2,ic, M, MP, Q, iok
INTEGER iokf12,iokf1,iokf2,icum,if12pos,if1pos,if2pos
INTEGER n1,n2,nq,k,na,nb
INTEGER NPART1,NPART2,NPART,D1,D2,NPARTA
INTEGER ir,xguess,phi,jr,reg2M,nproj,x,r,top,bot

INTEGER F1,F2,FACTORS(1000,2),rr(1000)

REAL const,tmp,prob(1000),xprob(1000)
REAL  NORM, NORMF, tprob

LOGICAL Even,True,False
!
! MPI declaration
!
INTEGER master, dest,rankM
INTEGER OUTFILE,COMM

INTEGER eloc
INTEGER masterU,commU
INTEGER masterM,commM

INTEGER TGROUPS
!
! MPI initialization
!

#ifdef USE_MPI
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
#endif

!
!  masterU=master of all universe
!  masterM is master of one Multiuniverse within Universe.
!
masterU=0
masterM=0


OUTFILE=10+rankM*nprocM+myidM

TGROUPS=nprocU/nprocM
!
! total number of qubits
!
 nq=n1+n2


write(OUTFILE,*) '-----------------------------------------------' 
write(OUTFILE,*) '|    This program is for Shor Factoring       |' 
write(OUTFILE,*) '|    Algorithm.                               |'
write(OUTFILE,*) '|    State Version  with Amplitude mapping.   |' 
write(OUTFILE,*) '-----------------------------------------------' 

if (myidU.eq.masterU) then
write(OUTFILE,*) '-----------------------------------------------' 
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
9001 Format('     N. processors :',I4)
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


!----------------------------------------------------------------
! Introduce noise at specific locations eloc
! We use three possible locations. First universe 
! is noiseless.
!
! Noise location depends on the universe at stake
!----------------------------------------------------------------

If(rankM.gt.0.and.myidM.eq.masterM) then
CALL RANDX(RAND,1)
eloc =   mod(int(2.*RAND(2,1)+rankM)+1,2)+1
write(OUTFILE,*)  
write(OUTFILE,*) '-----------------------------------------------' 
write(OUTFILE,*)' Group subject to noise=',rankM
write(OUTFILE,*)' At location (1/2) ',eloc
endif


!----------------------------------------------------------------
!
!  Choose the number to be factorized: 15, 21, 33, 35, 39, 55, 77.... 
!  and determine size of registers
!
!  M is the number we are trying to factor 
!  MP is the number we will work with and possibly modify:
!
!  Used MP=13*17+1  in the even aspect test.
!  Used MP=2**30  in the power of 2 aspect test.
!----------------------------------------------------------------

   M = MP

   if (1.eq.1) then
    print *,"will try to factor ",MP
   endif


!
! Following is done for master processor
!
If(myidM.eq.masterM) then

If(2**n1 .lt.MP**2/4) then 
Write(OUTFILE,*)'CHECK YOUR n1, YOU MAY HAVE PROBLEMS GETTING THE PERIOD see Gerjuoy eq 54'
endif
! 
! Select an integer (xguess) coprime to MP 
! Number of integers that are coprime to MP 
!
call EulerPhi(MP,phi)

Write(OUTFILE,*) 
Write(OUTFILE,*) '    Random xguess info:   '
Write(OUTFILE,*) 
Write(OUTFILE,*) '    - MP,phi',MP,phi

!
!Generate Random Integer between 2 and phi (EulerPhi(MP))
!
CALL RANDX(RAND,1)
ir = (phi - 2)*RAND(2,1) + 2


99 continue
Write(OUTFILE,*) 
Write(OUTFILE,*)'     Random Integer in (2, EulerPhi(MP))',ir
!
! Test if random number ir is coprime with MP
!  If so, accept as xguess

Write(OUTFILE,*)
If(gcd(ir,MP).eq.1.and.ir.gt.1) then
  Write(OUTFILE,*)'     coprime to MP: OK'
  xguess=ir
 else
  Write(OUTFILE,*)'  .......coprime to MP: failed ... will try again'
  ir=ir+1
  Go to 99
endif

 Write(OUTFILE,*) 
 Write(OUTFILE,*) '    xguess:', xguess

 Write(OUTFILE,*) 
write(OUTFILE,*) '-----------------------------------------------' 


endif  

!
! BCAST  xguess to all processors
!
#ifdef USE_MPI
call MPI_BCAST(xguess,1,MPI_INTEGER,masterM,commM,ierr)
#endif

! evaluate NPART1,NPART2,NPART

NPART1=2**n1/NPROCM
NPART2=2**n2/NPROCM
NPART =2**nq/NPROCM

Write(OUTFILE,*) 
Write(OUTFILE,8006) NPART
8006 Format('     Wavefunction arrays are of size ', I8, ' on this processor.') 
Write(OUTFILE,*) 
Write(OUTFILE,8007) NPART1,NPART2
8007 Format('     NPART1: ', I8, ' , NPART2:',I8) 
!
! Construct f_j=Mod(xguess^j,M) as a function
!
if (myidM.eq.masterM) then 
Write(OUTFILE,*) ' f_j=Mod(xguess^j,M)'
Write(OUTFILE,*)
endif
Do 20 j=0,NPART1-1 
  jv=j+myidM*NPART1
  If(j.lt.20.and.myidM.eq.masterM) Write(OUTFILE,8008) jv, f(xguess,jv,MP)
8008 format('    j= ',I3,'  f_j=',I8)
20 Continue

! ---------------------------------------------
!  Load the first register with all the 
!  integers less than or equal to n1
!  That is, Hadamard all n1 qubits in 
!  register 1 to get a uniform state distribution
!  In second register put |f(n1)>
!  Question:  How is this done in reality?
!  Answer: Need to do modular exponentiation on a QC.
! ---------------------------------------------

!
! allocate appropriate sizes
!
ALLOCATE( Psi(NPART))
ALLOCATE( B(nq) )
ALLOCATE( B1(n1) )
ALLOCATE( B2(n2) )

! ----------------------------------------------
!  Build Full 2**nq wave function
!  In this MPI version, we get the full state so 
!  we can study errors that could occur in register 2 etc.
!  Distribute over NPROCM
!  Could record Just Locations of 2^n1 entries--as a sparse column.
!  Construct Part of Psi on processor myid
! ---------------------------------------------

const=1./sqrt(2.**n1)
Do 10 i=0,NPART-1
Psi(i+1)= zero
!
!  iv is the absolute value of the  component label
!
iv=i + myidM* NPART
!
!  Section locates processor and seat location 
!  within processor of the iv component
!
ivsection=Int(iv/NPART)
ivseat=Mod(iv,NPART)
!
!  B is the binary version if iv
!
call dectobin(nq,iv,B)
!
!  Now split iv into n1 and n2 parts.
!
Do 15 j1=1,n1
B1(j1)=B(j1)
15 Continue
Do 16 j2=1,n2
B2(j2)=B(j2+n1)
16 Continue

call bintodec(n1,B1,D1)
call bintodec(n2,B2,D2)

!
!  D1 and D2 are decimal version of B1 and B2
!
If(myidM.eq.ivsection) then

!  Following does function placement on register 2
!  f(D1+1)=PowerMod(xguess,D1,MP)
!  If( D2.eq.PowerMod(xguess,D1,MP)) then

  If( D2.eq.f(xguess,D1,MP)) then 
  Psi(ivseat+1)= const*one
  endif
endif

10 Continue
! ----------------------------------------------------
!  The full wave function after we build the 
!  second register with f_j reads:
!
!   Psi = 1/sqrt(2^n1)  Sum_(j=0)^(2^n1-1) |j>  |f_j>
!
!   with  f_j = Mod(xguess^j, M)  having a period r
! ----------------------------------------------------

!-----------------------------------------------------
!  Apply QFT to First register n1
!-----------------------------------------------------

!
!Place noise at location eloc=1
!
If(rankM.gt.0.and.eloc.eq.1) then	
write(OUTFILE,*)'                                    ' 
write(OUTFILE,*)' --->  Noise before performing QFT' 
call Noise(nq,rankM,psi,NPART,commM)
endif
!  Continue Process


call QFT(nq,n1,psi,NPART,commM)

!
!Place noise at location eloc=2
!
If(rankM.gt.0.and.eloc.eq.2) then	
write(OUTFILE,*)' --->  Noise after performing QFT' 
call Noise(nq,rankM,psi,NPART,commM)
endif
!  Continue Process


! ---------------------------------------------------
!  Measure register 2 and determine wavefunction 
!  of register 1.
!
!  Generate a random choice for the measurement jr and
!  hence of the state |f(jr)>_n2 
! ---------------------------------------------------

!
! Master processor obtaines the random jr
!
If(myidM.eq.masterM) then
!
! Construct f_j=Mod(xguess^j,M) for random j=jr.
!
CALL RANDX(RAND,1)
jr =  RAND(2,1)*( 2**n1 -1)



write(OUTFILE,*) '-----------------------------------------------' 
write(OUTFILE,*)
write(OUTFILE,*) '    Measure register 2 '
write(OUTFILE,*)
write(OUTFILE,*) '         Choose a random possible value between:'
write(OUTFILE,*) '               1   -  2**n1-1 '


write(OUTFILE,*) '   rnd=',RAND(2,1), ',   ( 2**n1 -1)=',( 2**n1 -1)
write(OUTFILE,*)
write(OUTFILE,*) '       jr=',jr

reg2M = f(xguess,jr,MP)

write(OUTFILE,*) ' Measured value:'
write(OUTFILE,*)
write(OUTFILE,*) '   reg2M=PowerMod(xguess,jr,MP)=',reg2M
write(OUTFILE,*)
endif

nproj=n2

!
! Size of register after measuring
!
NPARTA=2**(nq-nproj)/nprocM

If(nprocM.gt.2**(nq-nproj))then
write(*,*)
write(*,*) '   The number of processors exceeds  projected wf size'
write(*,*) '   Could possibly use fewer processors '
write(*,*) 
endif

!
! Broadcasts the value measured in Reg 2
!
#ifdef USE_MPI
call MPI_Bcast(reg2M, 1, MPI_INTEGER, masterM, commM, ierr)
#endif

ALLOCATE( PsiO(NPARTA) )
ALLOCATE( PArray(nproj,2) )

Do i= 1,n2
 PArray(i,1)=i+n1
 call dectobin(n2,reg2M,B2)
 PArray(i,2)=B2(i)
enddo

!
!=================================
!  On to Projection   
!=================================
!
call ProjA(nq,nproj,PArray,Psi,PsiO,NPART,NPART1,commM)


ALLOCATE( result(NPART))

ic=1
Do i=0, NPARTA-1
  tmp= psiO(i+1)* conjg(psiO(i+1))
!
! only adds the probabilities which are larger than 
! a certain threshold, 0.001 in this case.
!
If(tmp.ge.0.001) then
  iv=i + myidM* NPARTA
  result(ic)=iv
  prob(ic)=tmp
!  write(OUTFILE,*) 'iv,ic, abs(psi) after projection', iv,ic,prob(ic)
  ic=ic+1
endif
enddo

tmp=0
iok=0
icum=0
Do 87 i=1,ic-1 

If(result(i).ne.0) then

top= result(i)
bot=2**n1
x=0
call CF(top, bot,x)
!  Is period r even and nonzero?

r=x

!
! Only small periods are considered, very larger ones, e.g. 512 
! would result in integers difficult to handle.
!
if (r.gt.20) then 
goto 87
endif
!
!need  period to get results
!

 f1=GCD(xguess**(r/2) +1,MP) 
 f2=GCD(xguess**(r/2) -1,MP) 
!
! If at least one of them is not zero
! the factorization process was 
! succesful (the other can be obtained dividing, easy)
!
If(abs(f1*f2).ne.1.and.abs(f1).lt.MP.and.abs(f2).lt.MP) then
  icum=icum+1
  iok=10
  tmp=tmp+prob(i)
  factors(icum,1)=f1
  factors(icum,2)=f2
  xprob(icum)=prob(i)
  rr(icum)=r
  If (abs(f1).ne.1.and.abs(f2).ne.1) then 
       iokf12=10
       if12pos=icum
  else if (abs(f1).ne.1) then 
       iokf1=10
       if1pos=icum
  else if (abs(f2).ne.1) then 
       iokf2=10
       if2pos=icum
  endif
   
endif
endif

87 Continue

!----------------------------------------------------------------
!
! Analysis and final output
!
!----------------------------------------------------------------

tmp=0.

If(iok.ne.10) then 
write(OUTFILE,*)'Unsuccesful on processor:',myidM
write(*,*)'Unsuccesful on processor:',myidM
write(*,*)
endif 


If(iok.eq.10) then 
 write(*,*)'      You were lucky on proc:',myidM
 write(*,8009) OUTFILE
8009 format('      check fort.',I3,' for details.')
 write(*,*)
 write(OUTFILE,*)
 write(OUTFILE,*)'      one of the factors is probably:'
 write(OUTFILE,*)
  tmp=0
  do i=1,icum
   write(OUTFILE,8010) factors(i,1),factors(i,2),rr(i)
8010 format('f1=',I4,'  f2=',I4,2x,'Period=',I10)
   tmp=tmp+xprob(i)
  enddo
 
  write(*,*)' Total (local) probability:',tmp, ' in proc:',myidM,' in group:',rankM



If (iokf12.eq.10) then
write(OUTFILE,*) '-----------------------------------------------------'
write(OUTFILE,*) '           Both factors were found simultaneously    '
write(OUTFILE,*) ' The factors are:' 
write(OUTFILE,8011)  factors(if12pos,1)
8011 format('                f1:',I4)
write(OUTFILE,8012)  factors(if12pos,2)
8012 format('                f2:',I4)
Write(OUTFILE,8013)tmp
8013 format(' With (local) probability:' ,F9.3)
Write(OUTFILE,*)

else If (iokf1.eq.10) then
write(OUTFILE,*) '-------------------------------------------'
write(OUTFILE,*) '           One factor was found'
write(OUTFILE,8011)  factors(if1pos,1)
Write(OUTFILE,8013)tmp
Write(OUTFILE,*)

else if (iokf2.eq.10) then
write(OUTFILE,*) '-------------------------------------------'
write(OUTFILE,*) '           One factor was found'
write(OUTFILE,8012)  factors(if2pos,2)
Write(OUTFILE,8013)tmp
Write(OUTFILE,*)

endif
endif 


#ifdef USE_MPI
CALL MPI_REDUCE(tmp,tprob,1,MPI_REAL,MPI_SUM,masterM,commM,ierr)
#endif
#ifndef USE_MPI
tprob=tmp
#endif

  if (myidM.eq.masterM) write(OUTFILE,*)' Total (group) probability:',tprob, ' in group:',rankM
  if (myidM.eq.masterM) write(*,*)' Total (group) probability:',tprob, ' in group:',rankM

DEALLOCATE( PsiO )
DEALLOCATE( Psi )

!==================================
End SUBROUTINE SHOR 

!end program main
