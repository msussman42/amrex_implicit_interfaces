!------------------------------------------------
!
!
! qcmpisubs.f90, QCMPI ver 0.0
!
! Frank Tabakin and Bruno Julia-Diaz
!
!   version Oct 20, 2008
! contains, 
!
!     subroutines:
!
!       bintodec
!       dectobin
!       OneOpA
!       TwoOpA
!       EulerPhi
!       splitn
!       ProjA
!       Randx
!       QFT
!       CF
!       SWAP
!       CPHASEK
!       HALL
!       HALL2
!       Entropy
!       EntropyP
!
!     function
!
!       GCD
!       SH
!------------------------------------------------
Module qcmpisubs
Contains
Subroutine dectobin(nq,D,B) 
!a decimal number, D, comes in and its binary equivalent, B, is returned
   Implicit None
   Integer :: i, nq
   Integer :: B(nq),D  
   Logical :: bool
   Do i=0,nq-1
      bool = Btest(D, i)
      B(nq-i)=0
      If(bool) B(nq-i)=1
   End Do
End Subroutine dectobin
!------------------------------------------------
Subroutine bintodec(nq,B,D) 
   Implicit None
   Integer :: i, nq
   Integer :: B(nq),D  
   Integer :: temp
   temp=0
   Do i=1,nq
      temp=B(i)*(2**(nq-i))+temp
   End Do
   D=temp
End Subroutine bintodec

!Subroutine OneOpA(nq,is,Op,psi,NPART)
SUBROUTINE OneOpA(nq,is,Op,psi,NPART,comm)

Use MPI
Use MPI_VARS

Implicit None
Integer :: nq, is,n,B(nq),qmark,stride,NPART,ntag1,ntag2 

Integer :: nv,nvp,nvsection,nvseat,nvpsection,nvpseat
Complex :: Op(2,2),psi(NPART),chi(NPART)
Complex  :: tmp01,tmp10
Complex, Parameter :: Zero =( 0.0, 0.0 )
Integer :: comm
Integer :: mpistat(MPI_STATUS_SIZE)
Integer :: req(1)

!
!  Modify wavefunction amplitudes due to a 2x2 operator OP 
!  acting on qubit "is" out of nq qubits.
!
!  Put Vector in Partner Form for each "is" case
!  Partner Form means finding and working with the wf components
!  for the q_1..... a=0.....q_nq   and q_1..... a=1.....q_nq       partners.
!

!
!  stride gives the increase in number when q=0 <--->  q=1 for qubit "is"
!
stride=2**(nq-is)
!write(9+myid,*)'nq,is,stride', nq,is,stride

! For is > x,  where nprocs=2^x,  the needed psi components
! are on the same processor
If(2**is.Gt.nprocs) Then
Do n=0,NPART-1
   Call dectobin(nq,n+myid*NPART,B)
   If(B(is).Eq.0) Then
      !  Following does 2 X 2 Multiplication
      chi(n+1)        =Op(1,1)*psi(n+1) + Op(1,2)*psi(n+stride+1) 
      chi(n+stride+1) =Op(2,1)*psi(n+1) + Op(2,2)*psi(n+stride+1)
   Endif
End Do
psi=chi
Return
Endif
! For is < x,  where nprocs=2^x,  the needed psi components
! are on  different  processors and labelled sends and receives are invoked
! Note nv=absolute value              -->  q_1 q_2.....qmark.... ....q_nq
!  and nvp=absolute value- of partner--->  q_1 q_2..not(qmark).......q_nq
!  Where pairs are qmark=0, not(qmark)=1 
!  and             qmark=1  not(qmark)=0 
!  nvsection and nvpsection are associated processors
!  nvseat and nvpseat are associated locations within processors

Do 20 n=0,NPART-1
nv = n+myid*NPART
nvsection=Int(nv/NPART)
nvseat=Mod(nv,NPART)
Call dectobin(nq,nv,B)
qmark=1
If(B(is).Eq.0) qmark=0
! nvp  are the partners to nv
!  nv is on current myid processor,  nvp are on processors above or below myid
nvp = nv + ((-1)**qmark)*stride
nvpsection=Int(nvp/NPART)
nvpseat=Mod(nvp,NPART)
!  Following Sends/Receives  the non-myid psi components
If(qmark.Eq.0) Then
ntag1=19
ntag2=19
tmp10=zero
tmp01=psi(nvseat+1)
Call MPI_IRECV(tmp10,1,MPI_COMPLEX,nvpsection,ntag2,COMM,req(1),ierr)
Call MPI_SEND(tmp01,1,MPI_COMPLEX,nvpsection,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)

!Following does 2 X 2 Multiplication
chi(n+1) =Op(1,1)*psi(nvseat+1) + Op(1,2)*tmp10
Else If (qmark.Eq.1) Then
ntag1=19
ntag2=19
tmp01=zero
tmp10=psi(n+1)
Call MPI_IRECV(tmp01,1,MPI_COMPLEX,nvpsection,ntag1,COMM,req(1),ierr)
Call MPI_SEND(tmp10,1,MPI_COMPLEX,nvpsection,ntag2,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)

!Following does 2 X 2 Multiplication
chi(n+1) =Op(2,1)*tmp01 + Op(2,2)*psi(n+1)
Endif
20 Continue
psi=chi 
End Subroutine OneOpA
!This is  version of TwoOpA where OP12 is not a product operator Op1 X Op2

!---------------------------------------------
!   TWOOPA
!---------------------------------------------

SUBROUTINE TwoOpA(nq,is1,is2,Op12,psi,NPART,COMM)

Use MPI
Use MPI_VARS

Implicit None
INTEGER nq, is1,is2,n,B(nq),qmark1,qmark2,stride1,stride2,NPART
INTEGER iv1,iv2,iv3,iv4
INTEGER ntag1,ntag2, ntag3,ntag4 
INTEGER nv,nvsection,nvseat 
INTEGER nvp1,nvp1section,nvp1seat
INTEGER nvp2,nvp2section,nvp2seat
INTEGER nvp3,nvp3section,nvp3seat
COMPLEX Op12(4,4),psi(NPART),chi(NPART)
COMPLEX  tmp11,tmp12,tmp13,tmp14
COMPLEX, PARAMETER :: Zero =( 0.0, 0.0 )
INTEGER mpistat(mpi_status_size),COMM
integer status(MPI_STATUS_SIZE,1), req(1)

call MPI_COMM_SIZE(COMM, nprocs, ierr)
call MPI_COMM_RANK(COMM, myid, ierr)
!
!  Modify wavefunction amplitudes due to a 4x4 operator OP12 
!  acting on qubits "is1" and "is2" out of nq qubits.
!
!  Put Vector in Partner Form for each "is1" "is2" case
!
!  stride1 gives the increase in number when q=0 <--->  q=1 for qubit "is1"
!  stride2 gives the increase in number when q=0 <--->  q=1 for qubit "is2"
!
stride1=2**(nq-is1)
stride2=2**(nq-is2)

! For stride1+stride2 < NPART,   the needed psi components
! are on the same processor
!
If(stride1+stride2.lt.NPART) then
Do 10 n=0,NPART-1
call dectobin(nq,n+myid*NPART,B)
If(B(is1).eq.0.and.B(is2).eq.0) then
!  Following does 4 X 4 Multiplication
iv1=n+1
iv2=n+stride1+1
iv3=n+stride2+1
iv4=n+stride1+stride2+1

chi(iv1)=Op12(1,1)*psi(iv1) + Op12(1,2)*psi(iv2)+ Op12(1,3)*psi(iv3) + Op12(1,4)*psi(iv4) 
chi(iv2)=Op12(2,1)*psi(iv1) + Op12(2,2)*psi(iv2)+ Op12(2,3)*psi(iv3) + Op12(2,4)*psi(iv4) 
chi(iv3)=Op12(3,1)*psi(iv1) + Op12(3,2)*psi(iv2)+ Op12(3,3)*psi(iv3) + Op12(3,4)*psi(iv4) 
chi(iv4)=Op12(4,1)*psi(iv1) + Op12(4,2)*psi(iv2)+ Op12(4,3)*psi(iv3) + Op12(4,4)*psi(iv4) 
endif
10 Continue
psi=chi
return
endif
! For stride1+stride2 > NPART,   the needed psi components
! are on  different  processors and labelled sends and receives are invoked

!  Note nv=absolute value-->  q_1 q_2.....qmark1.. qmark2  ........q_nq
!  Partners are nvp1-->  q_1 q_2... ......qmark1.. Not(qmark2) ....q_nq
!  Partners are nvp2-->  q_1 q_2.....Not(qmark1).. qmark2 ... .....q_nq
!  Partners are nvp3-->  q_1 q_2.....Not(qmark1).. Not(qmark2) ....q_nq
!  For example we have first set as   00,01,10,11
!  For example we have second set as  01,00,11,10
!  For example we have third  set as  10,11,00,01
!  For example we have 4th    set as  11,10,01,00
!   Each denoted as                nv,nvp1,nvp2,nvp3
Do 20 n=0,NPART-1
nv = n+myid*NPART
nvsection=Int(nv/NPART)
nvseat=Mod(nv,NPART)
call dectobin(nq,nv,B)
qmark1=1
If(B(is1).eq.0) qmark1=0
qmark2=1
If(B(is2).eq.0) qmark2=0
! nvp1,nvp2,nvp3  are the three partners to nv
! nv is on current myid processor,  nvp1  nvp2  nvp3 are on processors above or below  myid
nvp1 = nv + ((-1)**qmark2)*stride2 
nvp2 = nv + ((-1)**qmark1)*stride1
nvp3 = nv + ((-1)**qmark1)*stride1+ ((-1)**qmark2)*stride2
!
!  determine the associated processors
!
nvp1section=Int(nvp1/NPART)
nvp2section=Int(nvp2/NPART)
nvp3section=Int(nvp3/NPART)
!
!  determine the location within associated processors
!
nvp1seat=Mod(nvp1,NPART)
nvp2seat=Mod(nvp2,NPART)
nvp3seat=Mod(nvp3,NPART)
!  Following Sends/Receives  the non-myid psi components
!========================
!   case of 00,01,10,11
!
If(qmark1.eq.0.and.qmark2.eq.0) then
ntag1=nv+nvp1+nvp2+nvp3
tmp11=psi(nvseat+1)
!=============RECEIVES/SENDS==========
If(nvp1section.ne.nvsection) then
call MPI_IRECV(tmp12,1,MPI_COMPLEX,nvp1section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp11,1,MPI_COMPLEX,nvp1section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif
If(nvp2section.ne.nvsection) then
call MPI_IRECV(tmp13,1,MPI_COMPLEX,nvp2section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp11,1,MPI_COMPLEX,nvp2section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif
If(nvp3section.ne.nvsection) then
call MPI_IRECV(tmp14,1,MPI_COMPLEX,nvp3section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp11,1,MPI_COMPLEX,nvp3section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif
!
!  Check if nvp1 is on current processor nvsection==myid.
!
If(nvp1section.eq.nvsection) then
tmp12= psi(nvp1seat+1)
endif

!  Check if nvp2 is on current processor nvsection==myid.
!
If(nvp2section.eq.nvsection) then
tmp13= psi(nvp2seat+1)
endif
!
!  Check if nvp3 is on current processor nvsection==myid.
!
If(nvp3section.eq.nvsection) then
tmp14= psi(nvp3seat+1)
endif
!Following does 4 X 4 Multiplication
chi(n+1) =Op12(1,1)*tmp11 + Op12(1,2)*tmp12+ Op12(1,3)*tmp13 + Op12(1,4)*tmp14 

!========================
!   case of 01,00,11,10
!
else if (qmark1.eq.0.and.qmark2.eq.1) then 
ntag1=nv+nvp1+nvp2+nvp3
tmp12=psi(nvseat+1)
!=============RECEIVES/SENDS==========
If(nvp1section.ne.nvsection) then
call MPI_IRECV(tmp11,1,MPI_COMPLEX,nvp1section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp12,1,MPI_COMPLEX,nvp1section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif
If(nvp2section.ne.nvsection) then
call MPI_IRECV(tmp14,1,MPI_COMPLEX,nvp2section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp12,1,MPI_COMPLEX,nvp2section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif
If(nvp3section.ne.nvsection) then
call MPI_IRECV(tmp13,1,MPI_COMPLEX,nvp3section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp12,1,MPI_COMPLEX,nvp3section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif

If(nvp1section.eq.nvsection) then
tmp11= psi(nvp1seat+1)
endif
If(nvp2section.eq.nvsection) then
tmp14= psi(nvp2seat+1)
endif
If(nvp3section.eq.nvsection) then
tmp13= psi(nvp3seat+1)
endif
!Following does 4 X 4 Multiplication
chi(n+1) =Op12(2,1)*tmp11 + Op12(2,2)*tmp12+ Op12(2,3)*tmp13 + Op12(2,4)*tmp14 
!
!========================
!   case of 10,11,00,01
!
else if (qmark1.eq.1.and.qmark2.eq.0) then
ntag1=nv+nvp1+nvp2+nvp3
tmp13=psi(nvseat+1)
!=============RECEIVES/SENDS==========
If(nvp1section.ne.nvsection) then
call MPI_IRECV(tmp14,1,MPI_COMPLEX,nvp1section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp13,1,MPI_COMPLEX,nvp1section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif
If(nvp2section.ne.nvsection) then
call MPI_IRECV(tmp11,1,MPI_COMPLEX,nvp2section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp13,1,MPI_COMPLEX,nvp2section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif
If(nvp3section.ne.nvsection) then
call MPI_IRECV(tmp12,1,MPI_COMPLEX,nvp3section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp13,1,MPI_COMPLEX,nvp3section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif


If(nvp1section.eq.nvsection) then
tmp14= psi(nvp1seat+1)
endif
If(nvp2section.eq.nvsection) then
tmp11= psi(nvp2seat+1)
endif
If(nvp3section.eq.nvsection) then
tmp12= psi(nvp3seat+1)
endif
!Following does 4 X 4 Multiplication
chi(n+1) =Op12(3,1)*tmp11 + Op12(3,2)*tmp12+ Op12(3,3)*tmp13 + Op12(3,4)*tmp14 

!========================
!   case of 11,10,01,00 
!
else if (qmark1.eq.1.and.qmark2.eq.1) then
ntag1=nv+nvp1+nvp2+nvp3
tmp14=psi(nvseat+1)
!=============RECEIVES/SENDS==========
If(nvp1section.ne.nvsection) then
call MPI_IRECV(tmp13,1,MPI_COMPLEX,nvp1section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp14,1,MPI_COMPLEX,nvp1section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif
If(nvp2section.ne.nvsection) then
call MPI_IRECV(tmp12,1,MPI_COMPLEX,nvp2section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp14,1,MPI_COMPLEX,nvp2section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif
If(nvp3section.ne.nvsection) then
call MPI_IRECV(tmp11,1,MPI_COMPLEX,nvp3section,ntag1,COMM,req(1),ierr)
call MPI_SEND(tmp14,1,MPI_COMPLEX,nvp3section,ntag1,COMM,ierr)
call MPI_Wait(req(1),mpistat,ierr)
endif


If(nvp1section.eq.nvsection) then
tmp13= psi(nvp1seat+1)
endif
If(nvp2section.eq.nvsection) then
tmp12= psi(nvp2seat+1)
endif
If(nvp3section.eq.nvsection) then
tmp11= psi(nvp3seat+1)
endif

!Following does 4 X 4 Multiplication
chi(n+1) =Op12(4,1)*tmp11 + Op12(4,2)*tmp12+ Op12(4,3)*tmp13 + Op12(4,4)*tmp14 
endif
20 continue
psi=chi 
END SUBROUTINE 
!bottom of TwoOpA =============================
!
!------------------------------------------------
Subroutine CNOTA(nq,is1,is2,psi,NPART)

Use MPI
Use MPI_VARS
Implicit None

Integer nq, is1,is2,n,B(nq),qmark1,qmark2,stride1,stride2,NPART
Integer iv1,iv2,iv3,iv4
Integer ntag
Integer nv,nvsection,nvseat 
Integer nvp1,nvp1section,nvp1seat
Integer nvp2,nvp2section,nvp2seat
Integer nvp3,nvp3section,nvp3seat
Complex psi(NPART),chi(NPART)
Complex  tmp1,tmp2,tmp3,tmp4
Complex, Parameter :: One =( 1.0, 0.0 )
Complex, Parameter :: Zero =( 0.0, 0.0 )
Integer ::mpistat(MPI_STATUS_SIZE)

!
!  Modify wavefunction amplitudes due to a 4x4 operator OP12 
!  acting on qubits "is1" and "is2" out of nq qubits.
!
!  Put Vector in Partner Form for each "is1" "is2" case
!  Partner Form means that the operator act on wfs in each processor
!
!  See TwoOP for meaning of variables
!
stride1=2**(nq-is1)
stride2=2**(nq-is2)

! For stride1+stride2 < NPART,   the needed psi components
! are on the same processor
If(stride1+stride2.Lt.NPART) Then
Do n=0,NPART-1
   Call dectobin(nq,n+myid*NPART,B)
   If(B(is1).Eq.0.And.B(is2).Eq.0) Then

      iv1=n+1  !  00 PART
      iv2=n+stride2+1  !  01 PART
      iv3=n+stride1+1  !  10 PART
      iv4=n+stride1+stride2+1  !  11 PART

      chi(iv1)= ONE*psi(iv1)  
      chi(iv2)= One*psi(iv2) 
      chi(iv3)= One*psi(iv4) 
      chi(iv4)= One*psi(iv3)
   Endif
End Do
psi=chi
Return
Endif
! For stride1+stride2 > NPART,   the needed psi components
! are on  different  processors and labelled sends and receives are invoked
!
Do 20 n=0,NPART-1
nv = n+myid*NPART
nvsection=Int(nv/NPART)
nvseat=Mod(nv,NPART)
Call dectobin(nq,nv,B)
qmark1=1
If(B(is1).Eq.0) qmark1=0
qmark2=1
If(B(is2).Eq.0) qmark2=0
! nvp1,nvp2,nvp3  are the three partners to nv
! nv is on current myid processor,
! nvp1  nvp2  nvp3 are on processors above or below myid
!  
nvp1 = nv + ((-1)**qmark2)*stride2 
nvp2 = nv + ((-1)**qmark1)*stride1
nvp3 = nv + ((-1)**qmark1)*stride1+ ((-1)**qmark2)*stride2
nvp1section=Int(nvp1/NPART)
nvp2section=Int(nvp2/NPART)
nvp3section=Int(nvp3/NPART)
nvp1seat=Mod(nvp1,NPART)
nvp2seat=Mod(nvp2,NPART)
nvp3seat=Mod(nvp3,NPART)
!  Following Sends/Receives  the non-myid psi components

If(qmark1.Eq.0.And.qmark2.Eq.0) Then
chi(n+1) = psi(n+1)  

Else If (qmark1.Eq.0.And.qmark2.Eq.1) Then 
chi(n+1) = psi(n+1)  

Else If (qmark1.Eq.1.And.qmark2.Eq.0) Then
ntag=nv+nvp1+nvp2+nvp3

! Components 11 and 10 are on same processor
If(nvsection.Eq.nvp1section) chi(n+1)=psi(nvp1seat+1)


! Components 11 and 10 are not on same processor
!  Need to Transfer 10 and 11 components
If(nvsection.Ne.nvp1section) Then
tmp3=psi(n+1)
Call MPI_SENDRECV(tmp3,1,MPI_COMPLEX,nvp1section,ntag,tmp4,1,MPI_COMPLEX,nvp1section,ntag,MPI_COMM_WORLD,mpistat,ierr)
chi(n+1) = tmp4 
End If

Else If (qmark1.Eq.1.And.qmark2.Eq.1) Then
ntag=nv+nvp1+nvp2+nvp3

! Components 11 and 10 are on same processor
If(nvsection.Eq.nvp1section) chi(n+1)=psi(nvp1seat+1)

! Components 11 and 10 are not on same processor
!  Need to Transfer 10 and 11 components
If(nvsection.Ne.nvp1section) Then
tmp4=psi(n+1)
Call MPI_SENDRECV(tmp4,1,MPI_COMPLEX,nvp1section,ntag,tmp3,1,MPI_COMPLEX,nvp1section,ntag,MPI_COMM_WORLD,mpistat,ierr)
chi(n+1) = tmp3 
End If

Endif
20 Continue
psi=chi 
End Subroutine cnota


!==================================   
   
SUBROUTINE EulerPhi(n,phi)
INTEGER n,phi
INTEGER  i,igcd
phi=0  
Do i=1,n
igcd=GCD(n,i)
If(igcd.eq.1) phi=phi+1
end do 
end subroutine EulerPhi



SUBROUTINE SPLITN(nq,nq1,n,na,nb)
IMPLICIT NONE
INTEGER nq,nq1,nq2,n,na,nb,Ba(nq1),Bb(nq-nq1),iq
Logical :: bool
!
!  Splits space |n>_nq  ---->  | na >_nq1  | nb >_nq2
!  where  nq =nq1 + nq2
!
nq2=nq-nq1
!  i=0 to nq2-1  space b
!  i=nq2 to nq-1  space a
do iq=0,nq-1
bool = btest(n, iq)
If(iq.le.nq2-1) then
Bb(nq2-iq)=0
If(bool) Bb(nq2-iq)=1
else if (iq.ge.nq2) then
Ba(nq-iq)=0
If(bool) Ba(nq-iq)=1
endif
end do

call bintodec(nq1,Ba,na)
call bintodec(nq2,Bb,nb)
END SUBROUTINE  SPLITN
!bottom of SPLITN =============================



SUBROUTINE ProjA(nq,nproj,PArray,Psi,PsiO,NPART,NPARTA,COMM)

Use MPI
Use MPI_VARS

Implicit None

INTEGER    nq,nproj,i,j,n,nv2,NPART,NPARTA
INTEGER    nv,nvsection,nvseat
INTEGER    nv1,nv1section,nv1seat
COMPLEX, PARAMETER :: Zero =( 0.0, 0.0 )
COMPLEX  Psi(NPART),PsiO(NPARTA),chi(NPARTA)
INTEGER  overlap,rank,ic,icc
INTEGER  B(nq),B1(nq-nproj),B2(nproj),PArray(nproj,2),BA(nq)
REAL  NORM, NORMF
INTEGER mpistat(MPI_STATUS_SIZE),ntag,COMM


call MPI_COMM_SIZE(COMM, nprocs, ierr)
call MPI_COMM_RANK(COMM, myid, ierr)

   Do n=0,NPARTA-1
   PsiO(n+1)= Zero
   chi(n+1)= Zero
   enddo

   Do 10 n=0,NPART-1
   nv=n+myid*NPART
   nvsection=Int(nv/NPART)
   nvseat= Mod(nv,NPART)
   call dectobin(nq,nv,B)

!  B are the binary bits corresponding to nv  
!  Identify projected qubits and their values 

   Do 20 i=1, nq
   BA(i)=B(i)
   Do 30 j=1, nproj
   If(i.eq.PArray(j,1))BA(i)=99 
   30 continue    
   20 continue
   ic=1
   icc=1
   Do 25 i=1, nq
   If(BA(i).ne.99) then
   B1(ic)=B(i) 
   ic=ic+1
   endif
   If(BA(i).eq.99) then
   B2(icc)=B(i) 
   icc=icc+1
   endif
   25 continue
   call bintodec(nq-nproj,B1,nv1)
   call bintodec(nproj,B2,nv2)
!  Locate PsiO
   nv1section=Int(nv1/NPARTA)
   nv1seat= Mod(nv1,NPARTA)
   overlap=1
   Do 45 j=1, nproj
   If(B2(j).eq.PArray(j,2))then
   overlap=1*overlap 
   else
   overlap=0*overlap 
   endif
   45 continue    
   If(overlap.eq.1) then
!PsiO(nv1seat+1)  = Psi(n+1)+PsiO(nv1seat+1)
   If(nv1section.eq.myid)PsiO(nv1seat+1)= Psi(n+1)+PsiO(nv1seat+1)
   If(nv1section.ne.myid) write(*,*) 'got nv1section .ne. myid case'
!If(nv1section.ne.myid)chi(nv1seat+1)= Psi(n+1)+chi(nv1seat+1)
endif
10 continue  !End of n loop

Norm=0
Do 50 i=1, NPARTA
Norm=PsiO(i)*conjg(PsiO(i)) + Norm
50 continue

call MPI_Reduce(NORM,NORMF,1,MPI_REAL,MPI_SUM,0,COMM,ierr)
call MPI_BCAST(NORMF,1,MPI_REAL,0,COMM,ierr)

Do 60 i=1, NPARTA
PsiO(i)= PsiO(i)/sqrt(NormF)
60 continue

END SUBROUTINE ProjA
!bottom of ProjA ============================

!---------------------------------------------
!   RANDOM NUMBERS
!---------------------------------------------

SUBROUTINE Randx(pos,np)
implicit none
! Adapted from Mark Probert see http://www-users.york.ac.uk/~mijp1/
! variables required for random number generator
integer,parameter :: dp = kind(1.0d0)
integer np
real(kind=dp), dimension(1:2,1:np) :: pos 
integer ierr
integer :: rand_size
integer, allocatable, dimension(:) :: rand_seed
logical :: fixed_seed = .false. ! for nonrepeatability
!logical :: fixed_seed = .true. ! for repeatability
!f90 intrinsic time call
 character(len=10) :: system_time           !length is crucial ...
 real (kind=dp)    :: rtime
! fill with random numbers - random number code pinched from simple_md.f90
  call random_seed(size=rand_size)
  allocate(rand_seed(1:rand_size),stat=ierr)
  if (ierr/=0) stop 'Error in allocating rand_seed'
  if (fixed_seed) then
     rand_seed=223257167
  else
     call date_and_time(time=system_time)  !character string hhmmss.xxx
     read (system_time,*) rtime            !convert to real
     rand_seed = int(rtime * 1000.0_dp)    !0<rtime<235959999.0 which fits within huge(1)
  end if
  !
  call random_seed(put=rand_seed)
  deallocate(rand_seed,stat=ierr)
  if (ierr/=0) stop 'Error in deallocating rand_seed'
  !now random number generator is ready to use ..
  CALL RANDOM_NUMBER(pos)
end SUBROUTINE RANDX


SUBROUTINE QFT(nq,n1,psi,NPART,COMM)

Use MPI
Use MPI_VARS

Implicit None

INTEGER nq,n1,NPART,i,ic,k
COMPLEX had(2,2),psi(NPART)
INTEGER mpistat(mpi_status_size),COMM
COMPLEX, PARAMETER :: One  =( 1.0, 0.0 )

call MPI_COMM_SIZE(COMM, nprocs, ierr)
call MPI_COMM_RANK(COMM, myid, ierr)
!=================================
! Apply QFT to First register n1
! Define Hadamard
had(1,1)= One/sqrt(2.)
had(1,2)= One/sqrt(2.)
had(2,1)= One/sqrt(2.)
had(2,2)= -One/sqrt(2.)
!
Do ic =1,n1-1
call OneOpA(nq,ic,had,psi,NPART,COMM)
Do k=ic+1,n1
call CPHASEK(nq,k,ic,k,psi,NPART,COMM)
enddo
enddo
! Final Hadamard
call OneOpA(nq,n1,had,psi,NPART,COMM)
! Reverse order using pair swaps
Do i=1,n1/2
call SWAP(nq,i,n1+1-i,Psi,NPART,COMM)
enddo
RETURN
END SUBROUTINE QFT

SUBROUTINE CF(up,down,rv)

integer top,bot,i,rv,up,down
integer xcf(100),xcfp(100),p(100),q(100)


top=up
bot=down
!
! start
! 

i=0
 
10 i=i+1

    xcf(i)=top/bot
    top=mod(top,bot)
    if (top.eq.0) goto 90

    i=i+1

    xcf(i)=bot/top
    bot=mod(bot,top)

    if(bot.eq.0) goto 90

    goto 10

90 continue

!
! builds the fractions
!

rv=0
do j=1,i
!write(*,*) j,xcf(j)
        if ( j .eq. 1 ) then
          p(j) = xcf(j) 
          q(j) =  1
        else if ( j .eq. 2 ) then
          p(j) = xcf(j) * p(j-1) + 1
          q(j) = xcf(j) * q(j-1) 
        else
          p(j) = xcf(j) * p(j-1) + p(j-2)
          q(j) = xcf(j) * q(j-1) + q(j-2)
        end if
      r=q(j)
      If(2*Int(r/2.).eq.r.and.r.gt.0) then 
         rv=r
         goto 100
      endif
enddo
100 continue
!write(*,*) "rv:",rv
end subroutine CF


!---------------------------------------------
!   SWAP
!---------------------------------------------


SUBROUTINE SWAP(nq,is1,is2,psi,NPART,COMM)

Use MPI
Use MPI_VARS

Implicit None
INTEGER nq, is1,is2,n,B(nq),qmark1,qmark2,stride1,stride2,NPART
INTEGER iv1,iv2,iv3,iv4
INTEGER ntag
INTEGER nv,nvsection,nvseat 
INTEGER nvp1,nvp1section,nvp1seat
INTEGER nvp2,nvp2section,nvp2seat
INTEGER nvp3,nvp3section,nvp3seat
COMPLEX psi(NPART),chi(NPART)
COMPLEX  tmp1,tmp2,tmp3,tmp4
COMPLEX, PARAMETER :: Zero =( 0.0, 0.0 )
INTEGER mpistat(mpi_status_size),COMM


call MPI_COMM_SIZE(COMM, nprocs, ierr)
call MPI_COMM_RANK(COMM, myid, ierr)
!
!  Modify wavefunction amplitudes due to a 4x4 operator OP12 
!  acting on qubits "is1" and "is2" out of nq qubits.
!
!  Put Vector in Partner Form for each "is1" "is2" case
!  Partner Form means that the operator act on wfs in each processor
!
!  See TwoOP for meaning of variables
!
stride1=2**(nq-is1)
stride2=2**(nq-is2)

! For stride1+stride2 < NPART,   the needed psi components
! are on the same processor
If(stride1+stride2.lt.NPART) then
Do 10 n=0,NPART-1
call dectobin(nq,n+myid*NPART,B)
If(B(is1).eq.0.and.B(is2).eq.0) then

iv1=n+1  !  00 PART
iv2=n+stride2+1  !  01 PART
iv3=n+stride1+1  !  10 PART
iv4=n+stride1+stride2+1  !  11 PART

chi(iv1)=psi(iv1)  
chi(iv2)=psi(iv3) 
chi(iv3)=psi(iv2) 
chi(iv4)=psi(iv4)
endif
10 Continue
psi=chi
return
endif
! For stride1+stride2 > NPART,   the needed psi components
! are on  different  processors and labelled sends and receives are invoked
!
Do 20 n=0,NPART-1
nv = n+myid*NPART
nvsection=Int(nv/NPART)
nvseat=Mod(nv,NPART)
call dectobin(nq,nv,B)
qmark1=1
If(B(is1).eq.0) qmark1=0
qmark2=1
If(B(is2).eq.0) qmark2=0
! nvp1,nvp2,nvp3  are the three partners to nv
! nv is on current myid processor,
! nvp1  nvp2  nvp3 are on processors above or below myid
!  
nvp1 = nv + ((-1)**qmark2)*stride2 
nvp2 = nv + ((-1)**qmark1)*stride1
nvp3 = nv + ((-1)**qmark1)*stride1+ ((-1)**qmark2)*stride2
nvp1section=Int(nvp1/NPART)
nvp2section=Int(nvp2/NPART)
nvp3section=Int(nvp3/NPART)
nvp1seat=Mod(nvp1,NPART)
nvp2seat=Mod(nvp2,NPART)
nvp3seat=Mod(nvp3,NPART)
!  Following Sends/Receives  the non-myid psi components

If(qmark1.eq.0.and.qmark2.eq.0) then !nv==00 case
chi(n+1) = psi(n+1)  

else if (qmark1.eq.0.and.qmark2.eq.1) then !nv==01 case
!  Now nvp1==00,  nvp2=11;nvp3=10
ntag=nv+nvp1+nvp2+nvp3

! If Components 01 and 10 are on same processor
If(nvsection.eq.nvp3section) chi(n+1)=psi(nvp3seat+1)
! If Components 01 and 10 are not on same processor
!  Need to Transfer 01 and 10 components
If(nvsection.ne.nvp3section) then
tmp3=psi(n+1)
call MPI_SENDRECV(tmp3,1,MPI_COMPLEX,nvp3section,ntag,tmp4,1,MPI_COMPLEX,nvp3section,ntag,COMM,mpistat,ierr)
chi(n+1) = tmp4 
end if

else if (qmark1.eq.1.and.qmark2.eq.0) then !nv=10 case
!  Now nvp1==11,  nvp2=00;nvp3=01
ntag=nv+nvp1+nvp2+nvp3

! If Components 10 and 01 are on same processor
If(nvsection.eq.nvp3section) chi(n+1)=psi(nvp3seat+1)
! If Components 10 and 01 are not on same processor
!  Need to Transfer 10 and 01 components
If(nvsection.ne.nvp3section) then
tmp4=psi(n+1)
call MPI_SENDRECV(tmp4,1,MPI_COMPLEX,nvp3section,ntag,tmp3,1,MPI_COMPLEX, nvp3section,ntag,COMM,mpistat,ierr)
chi(n+1) = tmp3 

else if (qmark1.eq.1.and.qmark2.eq.1) then 
chi(n+1) = psi(n+1)  
end if

endif
20 continue
psi=chi 
END SUBROUTINE SWAP
!bottom of SWAP =============================


!---------------------------------------------
!  CPHASEK
!---------------------------------------------

SUBROUTINE CPHASEK(nq,is1,is2,k,psi,NPART,COMM)

Use MPI
Use MPI_VARS

Implicit None
INTEGER nq, is1,is2,k,n,B(nq),qmark1,qmark2,stride1,stride2,NPART
INTEGER iv1,iv2,iv3,iv4
INTEGER ntag
INTEGER nv,nvsection,nvseat 
INTEGER nvp1,nvp1section,nvp1seat
INTEGER nvp2,nvp2section,nvp2seat
INTEGER nvp3,nvp3section,nvp3seat
COMPLEX psi(NPART),chi(NPART)
COMPLEX  tmp1,tmp2,tmp3,tmp4,phi
COMPLEX, PARAMETER :: One =( 1.0, 0.0 )
COMPLEX, PARAMETER :: Zero =( 0.0, 0.0 )
COMPLEX, PARAMETER :: Eye =( 0.0, 1.0 )
REAL*8 pi, th
INTEGER mpistat(mpi_status_size),COMM


call MPI_COMM_SIZE(COMM, nprocs, ierr)
call MPI_COMM_RANK(COMM, myid, ierr)
!
!  Modify wavefunction amplitudes due to a 4x4 operator OP12 
!  acting on qubits "is1" and "is2" out of nq qubits.
!
!  Put Vector in Partner Form for each "is1" "is2" case
!  Partner Form means that the operator act on wfs in each processor
!
!  See TwoOP for meaning of variables
!

pi=4.D0*datan(1.D0)
th=2*pi/2**k 
!write(9+myid,*) 'psi',psi
stride1=2**(nq-is1)
stride2=2**(nq-is2)
phi =  cos(th)*One + sin(th)*Eye
!write(9+myid,*) 'k,pi,phi', k,pi,phi
! For stride1+stride2 < NPART,   the needed psi components
! are on the same processor
If(stride1+stride2.lt.NPART) then
Do 10 n=0,NPART-1
call dectobin(nq,n+myid*NPART,B)
If(B(is1).eq.0.and.B(is2).eq.0) then

iv1=n+1  !  00 PART
iv2=n+stride2+1  !  01 PART
iv3=n+stride1+1  !  10 PART
iv4=n+stride1+stride2+1  !  11 PART

chi(iv1)= psi(iv1)  
chi(iv2)= psi(iv2) 
chi(iv3)= psi(iv3) 
chi(iv4)=psi(iv4)*phi
endif
10 Continue
psi=chi
return
endif
! For stride1+stride2 > NPART,   the needed psi components
! are on  different  processors and labelled sends and receives are invoked
!
Do 20 n=0,NPART-1
nv = n+myid*NPART
nvsection=Int(nv/NPART)
nvseat=Mod(nv,NPART)
call dectobin(nq,nv,B)
qmark1=1
If(B(is1).eq.0) qmark1=0
qmark2=1
If(B(is2).eq.0) qmark2=0
! nvp1,nvp2,nvp3  are the three partners to nv
! nv is on current myid processor,
! nvp1  nvp2  nvp3 are on processors above or below myid
!  
nvp1 = nv + ((-1)**qmark2)*stride2 
nvp2 = nv + ((-1)**qmark1)*stride1
nvp3 = nv + ((-1)**qmark1)*stride1+ ((-1)**qmark2)*stride2
nvp1section=Int(nvp1/NPART)
nvp2section=Int(nvp2/NPART)
nvp3section=Int(nvp3/NPART)
nvp1seat=Mod(nvp1,NPART)
nvp2seat=Mod(nvp2,NPART)
nvp3seat=Mod(nvp3,NPART)
!  Following puts phase on only is1,is2= 11 cases.
!  No need for send/receives
 

If(qmark1.eq.0.and.qmark2.eq.0) then
chi(n+1) = psi(n+1)  

else if (qmark1.eq.0.and.qmark2.eq.1) then 
chi(n+1) = psi(n+1)  

else if (qmark1.eq.1.and.qmark2.eq.0) then
chi(n+1) = psi(n+1)  

else if (qmark1.eq.1.and.qmark2.eq.1) then
chi(n+1) = phi*psi(n+1)  
endif
20 continue
psi=chi 
END SUBROUTINE CPHASEK 
!bottom of CPHASEK =============================


!------------------------------------------------------------------
!                                                                 |
!  Noise.f90                                                      |
!                                                                 |
!       Introduces random noise on a given qubit                  |                                                                                                               |
!                                                                 |
!       May 3, 2008                                               |
!                                                                 |
!  Frank Tabakin & Bruno Julia Diaz                               |
!                                                                 |
!------------------------------------------------------------------


SUBROUTINE Noise(nq,rankM,psi1,NPART,commM)

Use MPI
Use MPI_VARS

Implicit None

INTEGER nq, rankM,NPART,qhit

REAl ALPHA,BETA,GAMMA
COMPLEX D2(2,2),psi1(NPART) 
INTEGER commM,OUTFILE
REAL*8 pi
Double Precision RAND(1:2,1:2)  !random number

call MPI_COMM_SIZE(commM, nprocM, ierr)
call MPI_COMM_RANK(commM, myidM, ierr)
OUTFILE=10+rankM*nprocM+myidM
write(OUTFILE,*)
write(OUTFILE,*) '--------------------------------'
write(OUTFILE,*)'         Noise Details:      '
write(OUTFILE,*)'                                '

!
!Introduce noise at location eloc=1
!
pi=4.D0*datan(1.D0)

If(rankM.gt.0) then
!
! Pick random qubit that gets hit by 
! a general one qubit unitary operator D2
!
call Randx(RAND,1)
qhit = mod( int(nq*RAND(2,1)),nq)+1

!
write(OUTFILE,*)'    Qubit hit: qhit=',qhit
!
! Pick random unitary one qubit operator
! This will keep the trace of rho=1.
! Users can insert their own operators.
!
call Randx(RAND,2)
alpha = pi*RAND(2,1) 
beta = 2.0*pi*RAND(2,2) 
!
! Small effect Limit
! alpha = 0.0001 
! beta = 0.000 
!
gamma = 0.0  

call U2(alpha,beta,gamma,D2)

write(OUTFILE,1000)alpha,beta,gamma
1000 format('   Alpha, Beta, Gamma=',3(f9.6,2x))
write(OUTFILE,*)
write(OUTFILE,*)'   Random 2X2 operator hitting that qubit is: D2='
write(OUTFILE,*)
write(OUTFILE,1001) D2(1,1),D2(1,2)
write(OUTFILE,1001) D2(2,1),D2(2,2)
1001 format(' (',2(F9.6,2x),'  ',2(F9.6,2x), ' )')
write(OUTFILE,*)

!
!  Act on Psi with one qubit operator acting on qbit nos "qhit."
!
call OneOpA(nq,qhit,D2,psi1,NPART,commM)
!
!  psi1 input and output--check this is OK
!
endif
END SUBROUTINE NOISE 


subroutine U2(Alpha,Beta, Gamma,d2)
IMPLICIT none
COMPLEX d2(2,2)
REAL Alpha,Beta, Gamma
COMPLEX, PARAMETER :: Zero =( 0.0, 0.0 )
COMPLEX, PARAMETER :: One =( 1.0, 0.0 )
COMPLEX, PARAMETER :: Eye =( 0.0, 1.0 )
REAL cs,sn
complex phase(4)
cs=cos(Beta/2.)
sn=sin(Beta/2.)
phase(1)=Cos(-Alpha/2. - Gamma/2.)*One+sin(-Alpha/2. - Gamma/2.)*Eye
phase(2)=Cos(-Alpha/2. + Gamma/2.)*One+sin(-Alpha/2. + Gamma/2.)*Eye
phase(3)=Cos(+Alpha/2. - Gamma/2.)*One+sin(+Alpha/2. - Gamma/2.)*Eye
phase(4)=Cos(+Alpha/2. + Gamma/2.)*One+sin(+Alpha/2. + Gamma/2.)*Eye

 d2(1,1)= phase(1)* cs 
 d2(1,2)=-phase(2)* sn 
 d2(2,1)=+phase(3)* sn
 d2(2,2)=+phase(4)* cs
end subroutine U2 

subroutine U4(Alpha,Beta, Gamma,d4)
IMPLICIT none
COMPLEX d4(4,4)
REAL Alpha,Beta, Gamma
REAL cs,sn
complex phase(16)
COMPLEX, PARAMETER :: Zero =( 0.0, 0.0 )
COMPLEX, PARAMETER :: One =( 1.0, 0.0 )
COMPLEX, PARAMETER :: Eye =( 0.0, 1.0 )
cs=cos(Beta/2.)
sn=sin(Beta/2.)


phase(1)=Cos(-3*Alpha/2. - 3*Gamma/2.)*One+sin(-3*Alpha/2. - 3*Gamma/2.)*Eye
phase(2)=Cos(-3*Alpha/2. - 1*Gamma/2.)*One+sin(-3*Alpha/2. - 1*Gamma/2.)*Eye
phase(3)=Cos(-3*Alpha/2. + 1*Gamma/2.)*One+sin(-3*Alpha/2. + 1*Gamma/2.)*Eye
phase(4)=Cos(-3*Alpha/2. + 3*Gamma/2.)*One+sin(-3*Alpha/2. + 3*Gamma/2.)*Eye

phase(5)=Cos(-1*Alpha/2. - 3*Gamma/2.)*One+sin(-1*Alpha/2. - 3*Gamma/2.)*Eye
phase(6)=Cos(-1*Alpha/2. - 1*Gamma/2.)*One+sin(-1*Alpha/2. - 1*Gamma/2.)*Eye
phase(7)=Cos(-1*Alpha/2. + 1*Gamma/2.)*One+sin(-1*Alpha/2. + 1*Gamma/2.)*Eye
phase(8)=Cos(-1*Alpha/2. + 3*Gamma/2.)*One+sin(-1*Alpha/2. + 3*Gamma/2.)*Eye

phase(9)= Cos(+1*Alpha/2. - 3*Gamma/2.)*One+sin(+1*Alpha/2. - 3*Gamma/2.)*Eye
phase(10)=Cos(+1*Alpha/2. - 1*Gamma/2.)*One+sin(+1*Alpha/2. - 1*Gamma/2.)*Eye
phase(11)=Cos(+1*Alpha/2. + 1*Gamma/2.)*One+sin(+1*Alpha/2. + 1*Gamma/2.)*Eye
phase(12)=Cos(+1*Alpha/2. + 3*Gamma/2.)*One+sin(+1*Alpha/2. + 3*Gamma/2.)*Eye

phase(13)=Cos(+3*Alpha/2. - 3*Gamma/2.)*One+sin(+3*Alpha/2. - 3*Gamma/2.)*Eye
phase(14)=Cos(+3*Alpha/2. - 1*Gamma/2.)*One+sin(+3*Alpha/2. - 1*Gamma/2.)*Eye
phase(15)=Cos(+3*Alpha/2. + 1*Gamma/2.)*One+sin(+3*Alpha/2. + 1*Gamma/2.)*Eye
phase(16)=Cos(+3*Alpha/2. + 3*Gamma/2.)*One+sin(+3*Alpha/2. + 3*Gamma/2.)*Eye



d4(1,1)=(cs**3)*phase(1)
d4(1,2)=-Sqrt(3.)*(cs**2)*sn*phase(2)
d4(1,3)=+Sqrt(3.)*(sn**2)*cs*phase(3)
d4(1,4)=-(sn**3)*phase(4)
!
d4(2,1)=+Sqrt(3.)*(cs**2)*sn*phase(5)
d4(2,2)=cs*(3*cs**2-2)*phase(6)
d4(2,3)=sn*(3*sn**2-2)*phase(7)
d4(2,4)=+Sqrt(3.)*(sn**2)*cs*phase(8)
!
d4(3,1)=+Sqrt(3.)*(cs**2)*sn*phase(9)
d4(3,2)=cs*(3*cs**2-2)*phase(10)
d4(3,3)=sn*(3*sn**2-2)*phase(11)
d4(3,4)=+Sqrt(3.)*(sn**2)*cs*phase(12)
!
d4(4,1)=(sn**3)*phase(13)
d4(4,2)=-Sqrt(3.)*(sn**2)*cs*phase(14)
d4(4,3)=+Sqrt(3.)*(cs**2)*sn*phase(15)
d4(4,4)=(cs**3)*phase(16)

end subroutine U4

SUBROUTINE PICKN(nq,nq2,PICK,n,na,nb)
IMPLICIT NONE
INTEGER nq,nq1,nq2,n,na,nb,Ba(nq-nq2),Bb(nq2),B(nq)
INTEGER ica,icb,iq
INTEGER PICK(nq)
!
!  Divides space |n>_nq  ---->  | na >_n1  | nb >_n2
!  where  nq =nq1 + nq2,  where space 2 is specified by PICK
!  IF PICK(i)=0  qubit i goes to system A
!  IF PICK(i)=1  qubit i goes to system B
!
nq1=nq-nq2
call dectobin(nq,n,B)
ica=1
icb=1
do iq=1,nq
if(PICK(iq).ne.0) then
Bb(icb)=B(iq)
icb=icb+1
endif
if(PICK(iq).eq.0) then
Ba(ica)=B(iq)
ica=ica+1
endif
enddo
call bintodec(nq1,Ba,na)
call bintodec(nq2,Bb,nb)

write(*,*)'B=',B
write(*,*)'Ba=',Ba
write(*,*)'Bb=',Bb
END SUBROUTINE PICKN
!bottom of PICKN =============================



SUBROUTINE PROJN(nq,nq2,PROJ,n,na,val)
IMPLICIT NONE
INTEGER nq,nq1,nq2,n,na,Ba(nq-nq2),B(nq)
INTEGER ica,iq,val
INTEGER PROJ(nq)
!
!  Projects out space B,| nb >_n2 to get | na >_n1 
!  that is <nb | n>_nq ---> | na >_n1  
!  where  nq =nq1 + nq2, 
!  space B is specified by PROJ
!  IF PROJ(i)=0  qubit i is not projected out fed to system A
!  IF PROJ(i)=1  qubit i is projected out with its value=0
!  IF PROJ(i)=2  qubit i is projected out with its value=1
!
nq1=nq-nq2
call dectobin(nq,n,B)
ica=1
val=1
do iq=1,nq
if(PROJ(iq).ne.0)then
if(PROJ(iq)-1.ne.B(iq))val=0*val
endif
if(PROJ(iq).eq.0) then
Ba(ica)=B(iq)
ica=ica+1
endif
enddo
call bintodec(nq1,Ba,na)

write(*,*)'B=',B
write(*,*)'Ba=',Ba
END SUBROUTINE PROJN
!bottom of PROJN =============================

SUBROUTINE HALL(nq,psi,NPART,COMM)

Use MPI
Use MPI_VARS

Implicit None
INTEGER nq,NPART,ic
COMPLEX had(2,2),psi(NPART)
INTEGER mpistat(mpi_status_size),COMM
COMPLEX, PARAMETER :: One  =( 1.0, 0.0 )

call MPI_COMM_SIZE(COMM, nprocs, ierr)
call MPI_COMM_RANK(COMM, myid, ierr)
!=================================
! Define Hadamard
had(1,1)= One/sqrt(2.)
had(1,2)= One/sqrt(2.)
had(2,1)= One/sqrt(2.)
had(2,2)= -One/sqrt(2.)
!
Do ic =1,nq
call OneOpA(nq,ic,had,psi,NPART,COMM)
enddo
RETURN
END SUBROUTINE HALL

SUBROUTINE HALL2(nq,psi,NPART,COMM)

Use MPI
Use MPI_VARS

Implicit None
INTEGER nq,NPART,n,np,npl,dest,seatn,OUTFILE, SH
COMPLEX psi(NPART), psi2(NPART),temp,getsum
INTEGER mpistat(mpi_status_size),COMM
COMPLEX, PARAMETER :: ZERO  =( 0.0, 0.0 )
call MPI_COMM_SIZE(COMM, nprocs, ierr)
call MPI_COMM_RANK(COMM, myid, ierr)
!------------------------------------------------------
! Do Hall on psi                                      !
! This version is for MPI distribution of Psi1 case   !
! Get signs for Hall using function SH                !
!------------------------------------------------------
OUTFILE=10
Do 21 n=0,NPART*nprocs-1

temp=zero

Do 31 np=0, NPART-1
npl=np + myid*NPART
temp= temp+SH(nq,n,npl)*psi(np+1)/sqrt(2.**nq) 
31 continue

!
!  Now get processor (dest) and seat(seatn) for n value
!

dest= Int(n/NPART)
seatn= Mod(n,NPART)


CALL MPI_REDUCE(TEMP,GETSUM,1,MPI_COMPLEX,MPI_SUM,0,comm,IERR)
    if (ierr.ne.0) write(OUTFILE,*) 'Error in mpi_reduce, myid=',myid
CALL MPI_BCAST(GETSUM, 1, MPI_COMPLEX, 0, comm, IERR)
     if (ierr.ne.0) write(OUTFILE,*) 'Error in mpi_bcast, myid=',myid
!
!  Put the getsum into the right place
!

If(myid.eq.dest)psi2(seatn+1)=getsum
21 continue 

psi=psi2

!-------------------
! END Hall on psi  |
!-------------------
!                        
RETURN
END SUBROUTINE HALL2


!  Entropy.f90                                                    |
!                                                                 |
! This version collects  Rho into masterU=0 processor             |
!                                                                 |
!                                                                 |
!   October 29, 2008                                              |
!                                                                 |
!   Bruno Julia Diaz & Frank Tabakin                              |
!                                                                 |
!------------------------------------------------------------------
!                                                                 |
! Output:                                                         |
!      fort.10    Master of All groups                            |
!                                                                 |
!      fort.10+rankM*nprocM+myidM   each processor output         |
!                                                                 |
!------------------------------------------------------------------

Subroutine Entropy(nq,commM,rankM,NPART,PSI1)
Use MPI
Use MPI_VARS

Implicit None
INTEGER NPART,rankM
COMPLEX PSI1(NPART)
COMPLEX, PARAMETER   :: Zero = ( 0.0, 0.0 )
COMPLEX, PARAMETER   :: One  = ( 1.0, 0.0 )
COMPLEX, ALLOCATABLE :: Psi3(:) 
REAL, ALLOCATABLE :: W(:) 

REAL Groupprob
INTEGER nq,NBITS
INTEGER i,j,n,NPRHO,k
INTEGER M
INTEGER OUTFILE
!
!
!  Multiuniverse variables
!
INTEGER masterU,commU
INTEGER masterM,commM
INTEGER TGROUPS
!
INTEGER ntag,MPISTAT(mpi_status_size)
!COMPLEX temp, temp2,TR1,TR2
REAL temp, temp2,TR1,TR2
REAL S,epsilon, tracerho,tracerho2
!
! MPI declaration 
!
commU=MPI_comm_world
!
call MPI_COMM_SIZE(commU,nprocU,ierr)
call MPI_COMM_RANK(commU,myidU,ierr)
!
! Now for the multiUniverse within the Universe
!
call MPI_COMM_SIZE(commM,nprocM,ierr) 
call MPI_COMM_RANK(commM,myidM,ierr)
!
!  masterU=master of all universe
!  masterM is master of one Multiuniverse within Universe.
!
masterU=0
masterM=0
!
! Outfile
!
OUTFILE=10+rankM*nprocM+myidM
!
TGROUPS=nprocU/nprocM
!
! Memory allocation for distributed vectors and arrays.
!
!  NPRHO = Number of Processors over which RHO is distributed.
!  NPRHO = 1 puts full RHO on myidU=0 processor.
!  NBITS X NBITS= Array size for  RHO.
!
!     Set up process grid that is as close to square as possible
!
     NBITS=2**nq
!
 ALLOCATE(W(NBITS))
 ALLOCATE( Psi3(NBITS) ) 
!
!

!=============== send Psi's to last node===============

! universal tag
ntag=17
!
! Gather Psi1 from the commM processors into a full array Psi3 on myidM=0
!
call MPI_GATHER(PSI1,NPART,MPI_COMPLEX,PSI3,NPART,MPI_COMPLEX,0,commM,ierr)
 if (ierr/=MPI_SUCCESS) stop 'Error in Entropy mpi_gather '

If (myidM.eq.masterM) then
! Probability of this group
!
!  Make rankM=0,  the noiseless multiuniverse with probability 1-x
!
!  For example:
!
!        for x=.2;  the rankM=0 has weight 80%
!
!        The 20%  can be distributed in various ways
!        among the remaining NGROUPS-1 multiuniverses.
!
!        One way is to have them all with single one-qubit noise
!         with equal weights of x/(NGROUPS-1)
!
!  Taking x=0.2
!
If(rankM.eq.0)Groupprob =1.-0.2
If(rankM.gt.0)Groupprob= 0.2/(TGROUPS-1)
If(TGROUPS.eq.1)Groupprob =1.
!
!   Another possible choice is uniform distribution:
!   Groupprob=REAL(nprocM)/REAL(nprocU)
!
!   One could also distribute the 20% over a set of single one-qubit errors,
!   plus some double one-qubit errors, plus some single two qubit errors.
!   Such a study could probe algorithm noise sensitivity.


!
! Calculate Trace of the density matrix of the group
!
 tracerho =Groupprob
 tracerho2=Groupprob*Groupprob
  write(OUTFILE,*) '-----------------------------------------------' 
  Write(OUTFILE,*) 
  Write(OUTFILE,*) '   Properties of the Density matrix of my group'
  Write(OUTFILE,*) '          (from groupprob)                     '
  Write(OUTFILE,9019) Groupprob
  Write(OUTFILE,*) 
  9019 Format('    Probability =',F9.4)
  Write(OUTFILE,9020)  tracerho
  9020 Format('             TraceRho = (',F9.4,')')

  Write(OUTFILE,9021)  tracerho2
  9021 Format('             TraceRho2= (',F9.4,')')
  Write(OUTFILE,*) 
  write(OUTFILE,*) '-----------------------------------------------' 
  endif


!---------------------------------------------------------|
! Builds final density matrix, collects from all MASTERS  |
!---------------------------------------------------------|

!   Construct Psi3 for all rankM and place on myidM=0
!   Continue for all rankM values.  
!
!   When done the full RHO is constructed 
!   in VONN and then CHEEV is called
!   to get the full ensemble average density matrix.
!   Full RHO has global dimension NBITS X NBITS

If(myidM.eq.masterM) then
Psi3=  Psi3*sqrt(Groupprob)
endif

!---------------------------------------------------------|
! Get Eigenvalues of RHO
! CHEEV computes just eigenvalues 'N'
!
!---------------------------------------------------------|
CALL VONN(commU,commM,rankM,nbits,psi3,TGROUPS,W,M)

! Get VonNeuman entropy and trace[rho^2]=tr2 from eigenvalues of RHO
! Print results
! For Eigenvalues cut off

epsilon=5.682566E-08

write(OUTFILE,*) 'NBITS,M=',NBITS,M

!
! Prints out Eigenvalues
!
if (myidU.eq.0) then 
  write(OUTFILE,*) '-----------------------------------------------' 
  write(OUTFILE,*) '                                               ' 
write (outfile,*)  '        Non-zero Eigenvalues of FULL RHO'
  write(OUTFILE,*) '                                               ' 
Do i=1,M
If(abs(W(i)).ge.epsilon)write(OUTFILE,*) 'i,W(i)=',i,W(i)
enddo 
  write(OUTFILE,*) '                                               ' 
  write(OUTFILE,*) '-----------------------------------------------' 
endif


!-----------------------------------
! Builds Entropy and writes out
!-----------------------------------
temp2=0.0
Do 70 i=1,M
If(abs(W(i)).LT.epsilon)W(i)=0.0
If(W(i).GT.0.0) temp2=temp2- W(i)*log(W(i))
70 continue
S = temp2

!-----------------------------------
! Computes the trace of Rho
!-----------------------------------
temp2=0.0
Do 75 i=1,M
If(Abs(W(i)).LT.epsilon)W(i)=0.0
If(W(i).GT.0.0) temp2=temp2+ W(i)
75 continue
Tr1 = temp2

!-----------------------------------
! Computes the trace of Rho**2
!-----------------------------------
temp2=0.0
Do 80 i=1,M
If(Abs(W(i)).LT.epsilon)W(i)=0.0
If(W(i).GT.0.0) temp2=temp2+ W(i)*W(i)
80 continue
Tr2 = temp2


!=============Calculate  Trace =================
if (myidU.eq.masterU) then 
Write(OUTFILE,*) 
Write(OUTFILE,*) '   Properties of the FULL Density matrix '
Write(OUTFILE,*) 
Write(OUTFILE,9020) Tr1
Write(OUTFILE,9021) Tr2
Write(OUTFILE,9022) S
  9022 Format('             Entropy  = (',F9.4,')')
Write(OUTFILE,*) 
write(OUTFILE,*) '-----------------------------------------------' 
endif

DEALLOCATE( PSI3 ) 
DEALLOCATE( W ) 
Return
End subroutine ENTROPY


subroutine VONN(commU,commM,rankM,n,psi,TGROUPS,W,M)
! Eigenvalues of large matrix A Constructed from vector Psi
! n =NBITS=2**nq      number of rows/columns in matrix A

Use MPI
Use MPI_VARS

Implicit None
! Multiuniverse variables
!
      INTEGER masterU,commU
      INTEGER masterM,commM
      INTEGER ntag
      integer TGROUPS,k,ng,myunit
      integer  rankM,mpistat(mpi_status_size,TGROUPS-1),req(TGROUPS-1)
      integer :: n    ! problem size 
      integer :: i,j,M
      complex, dimension(:,:), allocatable :: A
      complex, dimension(:,:), allocatable :: A0
      complex Psi(n),PHI(n),PHI2(n)
 
      
!.. Parameters ..for CHEEV
INTEGER   INFO, LDA, LWORK
COMPLEX  WORK(n,n)
CHARACTER JOBZ, UPLO
REAL W(n), RWORK(max(1, 3*n-2)),epsilon
COMPLEX, PARAMETER   :: One  =( 1.0, 0.0 )
COMPLEX, PARAMETER   :: EYE =( 0.0, 1.0 )
COMPLEX, PARAMETER   :: Zero =( 0.0, 0.0 )
!
epsilon=5.682566E-08
JOBZ = 'N'
UPLO = 'U'
LDA = n
LWORK = max(1,2*n-1)
INFO = 0
!
! 
masterU=0
masterM=0
!
call MPI_COMM_SIZE(commU,nprocU,ierr)
call MPI_comm_rank(commU,myidU,ierr)
!
! Now for the multiUniverse within the Universe
!
call MPI_COMM_SIZE( commM, nprocM, ierr) 
call MPI_COMM_RANK( commM, myidM, ierr)
!
! Construct arrays
! Global structure:  matrix A of n rows and n columns
!                    matrix A0 of n rows and n column
!    
 
! Initialize arrays     
      allocate(A(n,n)) 
      allocate(A0(n,n)) 


!  Distribute PSI=Psi3 
!
! Prepare Full State vector (PHI) for rankM = 0
!

If(rankM.eq.0.and.myidM.eq.0) then
 PHI=psi
 do i=1,n
 do j=1,n
A0(i,j) =   PHI(i)*conjg(PHI(j))
 enddo
  enddo
endif
!
!  Note A0 is the density matrix for the noiseless channel
!  It could be used, along with A, to evaluate the Fidelity
!
A=A0
!
ntag=88
!
!  SEND PSI  for each noise group (rankM>0) to masterU receive as PHI2
!
 Do ng=1,TGROUPS-1   !   Loop over Groups rankM>0
 k=nprocM*ng
!Use IRECV method:
If(myidU.eq.0)call MPI_IRECV(PHI2,N,MPI_COMPLEX,k,ntag,commU,req(ng),ierr)
If(myidU.eq.k)call MPI_SEND(psi,N,MPI_COMPLEX,masterU,ntag,commU,ierr)
If(myidU.eq.0)call MPI_Wait(req(ng),mpistat,ierr)
!
!
!   Add rankM>0 density matrix to full density matrix
!   stored on myidU=0. 
!

 do i=1,n
    do j=1,n
A(i,j) =  A(i,j) + PHI2(i)*conjg(PHI2(j))
    enddo
enddo

enddo

!
! Call LAPACK library routine
! Ask CHEEV to compute:  just eigenvalues 'N'
!
!  Matrix A has global dimension N X N. 
!
!write(89,*)' myidU  ',myidU
!If(myidU.eq.0) write(89,*)' A  ',A

CALL CHEEV(JOBZ, UPLO,n, A, LDA, W, WORK, LWORK, RWORK,INFO)
! Count eigenvalues >epsilon
M=n
! Deallocate the arrays
 
      deallocate(A, A0)
 
return 
END SUBROUTINE VONN

!bottom of Entropy ====================================

!  EntropyP.f90                                                   |
!                                                                 |
! This version distributes Psi over nprocM                        |
! and distributes Rho over NPRHO processors                       |
! to allow parallel treatment of density matrix                   |
!                                                                 |
!      July 29, 2008                                              |
!                                                                 |
!   Bruno Julia Diaz & Frank Tabakin                              |
!                                                                 |
!------------------------------------------------------------------
!                                                                 |
! Output:                                                         |
!      fort.10    Master of All groups                            |
!                                                                 |
!      fort.10+rankM*nprocM+myidM   each processor output         |
!                                                                 |
!------------------------------------------------------------------

Subroutine EntropyP(nq,commM,rankM,NPART,PSI1)
Use MPI
Use MPI_VARS

Implicit None
INTEGER NPART,rankM
COMPLEX PSI1(NPART)
COMPLEX, PARAMETER   :: Zero = ( 0.0, 0.0 )
COMPLEX, PARAMETER   :: One  = ( 1.0, 0.0 )
COMPLEX, ALLOCATABLE :: Psi3(:) 
REAL, ALLOCATABLE :: W(:) 

REAL Groupprob
INTEGER nq,NBITS
INTEGER i,j,n,NPRHO,k
INTEGER M,NB,prow,pcol
INTEGER OUTFILE
!
!  BLACS INFO
!  pcol,prow are 2D processor numbers
!
!  Total number of processors in 2D processor grid is nb ==nprocU
!
!  Multiuniverse variables
!
INTEGER masterU,commU
INTEGER masterM,commM
INTEGER TGROUPS
!
INTEGER ntag,MPISTAT(mpi_status_size)
REAL temp, temp2,TR1,TR2
REAL S,epsilon, tracerho,tracerho2
!
! MPI declaration 
!
commU=MPI_comm_world
!
call MPI_COMM_SIZE(commU,nprocU,ierr)
call MPI_COMM_RANK(commU,myidU,ierr)
!
! Now for the multiUniverse within the Universe
!
call MPI_COMM_SIZE(commM,nprocM,ierr) 
call MPI_COMM_RANK(commM,myidM,ierr)
!
!  masterU=master of all universe
!  masterM is master of one Multiuniverse within Universe.
!
masterU=0
masterM=0
!
! Outfile
!
OUTFILE=10+rankM*nprocM+myidM
!
TGROUPS=nprocU/nprocM
!
! Memory allocation for distributed vectors and arrays.
!
!  NPRHO = Number of Processors over which RHO is distributed.
!  NPRHO = 1 puts full RHO on myidU=0 processor.
!  NBITS X NBITS= Array size for  RHO.
!
!     Set up process grid that is as close to square as possible
!
     NBITS=2**nq
     prow = INT( SQRT( REAL(nprocU) ) )
     pcol = prow
     nb=prow*pcol
     NPRHO= nb
!
 ALLOCATE(W(NBITS))
 ALLOCATE( Psi3(NBITS) ) 
!
!

!write(88,*)'npart,myidM,rankM,psi1=',npart,myidM,rankM,psi1
!=============== send Psi's to last node===============

! universal tag
ntag=17
!
! Gather Psi1 from the commM processors into a full array Psi3 on myidM=0
!
call MPI_GATHER(PSI1,NPART,MPI_COMPLEX,PSI3,NPART,MPI_COMPLEX,0,commM,ierr)
 if (ierr/=MPI_SUCCESS) stop 'Error in Entropy mpi_gather '

!If(myidM.eq.0)write(88,*)'myidM,rankM,psi3=',myidM,rankM,psi3
If (myidM.eq.masterM) then
! Probability of this group
!
!  Make rankM=0,  the noiseless multiuniverse with probability 1-x
!
!  For example:
!
!        for x=.2;  the rankM=0 has weight 80%
!
!        The 20%  can be distributed in various ways
!        among the remaining NGROUPS-1 multiuniverses.
!
!        One way is to have them all with single one-qubit noise
!         with equal weights of x/(NGROUPS-1)
!
!  Taking x=0.2
!
If(rankM.eq.0)Groupprob =1.-0.2
If(rankM.gt.0)Groupprob= 0.2/(TGROUPS-1)
If(TGROUPS.eq.1)Groupprob =1.
!
!   Another possible choice is uniform distribution:
!   Groupprob=REAL(nprocM)/REAL(nprocU)
!
!   One could also distribute the 20% over a set of single one-qubit errors,
!   plus some double one-qubit errors, plus some single two qubit errors.
!   Such a study could probe algorithm noise sensitivity.


!
! Calculate Trace of the density matrix of the group
!
 tracerho =Groupprob
 tracerho2=Groupprob*Groupprob
  write(OUTFILE,*) '-----------------------------------------------' 
  Write(OUTFILE,*) 
  Write(OUTFILE,*) '   Properties of the Density matrix of my group'
  Write(OUTFILE,*) '          (from groupprob)                     '
  Write(OUTFILE,9019) Groupprob
  Write(OUTFILE,*) 
  9019 Format('    Probability =',F9.4)
  Write(OUTFILE,9020)  tracerho
  9020 Format('             TraceRho = (',F9.4,')')

  Write(OUTFILE,9021)  tracerho2
  9021 Format('             TraceRho2= (',F9.4,')')
  Write(OUTFILE,*) 
  write(OUTFILE,*) '-----------------------------------------------' 
  endif


!---------------------------------------------------------|
! Builds final density matrix, collects from all MASTERS  |
!---------------------------------------------------------|

!   Construct Psi3 for all rankM and place on myidM=0
!   Continue for all rankM values.  
!
!   When done the full RHO is constructed 
!   in VONNP and then PCHEEV is called
!   to get the full ensemble average density matrix.
!   This procedure avoids storing large Rho matrix;
!   the local matrix size for Rho is prow X pcol
!   where prow = sqrt(nprocU)= pcol or less.   
!   Full RHO has global dimension NBITS X NBITS

If(myidM.eq.masterM) then
Psi3=  Psi3*sqrt(Groupprob)
endif

!---------------------------------------------------------|
! Get Eigenvalues of RHO
! PCHEEVX computes just eigenvalues 'N', for range 'V'
! Range of eigenvalues VL to VU =  0 to 1.
!
!---------------------------------------------------------|
CALL VONNP(commU,commM,rankM,nbits,prow,pcol,psi3,TGROUPS,W,M)

! Get VonNeuman entropy and trace[rho^2]=tr2 from eigenvalues of RHO
! Print results
! For Eigenvalues cut off

epsilon=5.682566E-08

write(OUTFILE,*) 'NBITS,M=',NBITS,M

!
! Prints out Eigenvalues
!
if (myidU.eq.0) then 
  write(OUTFILE,*) '-----------------------------------------------' 
  write(OUTFILE,*) '                                               ' 
write (outfile,*)  '        Non-zero Eigenvalues of FULL RHO'
  write(OUTFILE,*) '                                               ' 
Do i=1,M
If(abs(W(i)).ge.epsilon)write(OUTFILE,*) 'i,W(i)=',i,W(i)
enddo 
  write(OUTFILE,*) '                                               ' 
  write(OUTFILE,*) '-----------------------------------------------' 
endif


!-----------------------------------
! Builds Entropy and writes out
!-----------------------------------
temp2=0.0
Do 70 i=1,M
If(abs(W(i)).LT.epsilon)W(i)=0.0
If(W(i).GT.0.0) temp2=temp2- W(i)*log(W(i))
70 continue
S = temp2

!-----------------------------------
! Computes the trace of Rho
!-----------------------------------
temp2=0.0
Do 75 i=1,M
If(Abs(W(i)).LT.epsilon)W(i)=0.0
If(W(i).GT.0.0) temp2=temp2+ W(i)
75 continue
Tr1 = temp2

!-----------------------------------
! Computes the trace of Rho**2
!-----------------------------------
temp2=0.0
Do 80 i=1,M
If(Abs(W(i)).LT.epsilon)W(i)=0.0
If(W(i).GT.0.0) temp2=temp2+ W(i)*W(i)
80 continue
Tr2 = temp2


!=============Calculate  Trace  Use Reduce to sum=================
if (myidU.eq.masterU) then 
Write(OUTFILE,*) 
Write(OUTFILE,*) '   Properties of the FULL Density matrix '
Write(OUTFILE,*) 
Write(OUTFILE,9020) Tr1
Write(OUTFILE,9021) Tr2
Write(OUTFILE,9022) S
  9022 Format('             Entropy  = (',F9.4,')')
Write(OUTFILE,*) 
write(OUTFILE,*) '-----------------------------------------------' 
endif

DEALLOCATE( PSI3 ) 
DEALLOCATE( W ) 
Return
End subroutine ENTROPYP


subroutine VONNP(commU,commM,rankM,n,prow,pcol,psi,TGROUPS,W,M)
! Eigenvalues of large matrix A Constructed from vector Psi
 
!              prow    number of rows in proc grid
!              pcol    number of columns in proc grid
!              n       number of rows/columns in matrix A
!              nb      matrix distribution block size   

Use MPI
Use MPI_VARS

Implicit None
! Multiuniverse variables
!
      INTEGER masterU,commU
      INTEGER masterM,commM
      integer TGROUPS,k,ng,myunit
      INTEGER ntag
      integer rankM,mpistat(mpi_status_size,TGROUPS-1),req(TGROUPS-1)
!      integer rankM,mpistat(mpi_status_size)
      integer :: n, nb    ! problem size and block size
      integer :: myArows, myAcols   ! size of local subset of global array
      integer :: i,j, igrid,jgrid, iproc,jproc, myi,myj, p
      complex, dimension(:,:), allocatable :: A
      complex, dimension(:,:), allocatable :: A0
      complex, dimension(:,:), allocatable :: Z
      complex Psi(n),PHI(n),PHI2(n)
      real W(n)
 
      integer :: numroc   ! blacs routine
      integer :: me, procs, icontxt, prow, pcol, myrow, mycol  ! blacs data
      integer :: info    ! scalapack return value
      integer, dimension(9) :: IDESA, IDESZ ! scalapack array desc
      
!.. Parameters ..for PCHEEVX
      INTEGER  IZ,JZ,IA,JA,LIWORK,LRWORK,LWORK,M,NZ
      INTEGER ABSTOL, ORFAC	 ! Orthogonalization tolerances
      REAL    VL,VU		 ! Eigenvalue range
      REAL    Epsilon		 ! Eigenvalue cutoff
      INTEGER, dimension(:), allocatable :: ICLUSTR,IFAIL,IWORK
      REAL, dimension(:), allocatable :: GAP,RWORK
      COMPLEX, dimension(:), allocatable :: WORK

COMPLEX, PARAMETER   :: One  =( 1.0, 0.0 )
COMPLEX, PARAMETER   :: EYE =( 0.0, 1.0 )
COMPLEX, PARAMETER   :: Zero =( 0.0, 0.0 )

  epsilon=5.682566E-08
! problem description
 
      write(90,*)prow
      write(90,*)pcol
      write(90,*)n
      nb=pcol*prow
      write(90,*)nb

!
! performs checks on the sizes
! 
      if (((n/nb) < prow) .or. ((n/nb) < pcol)) then
      print *,"Problem size too small for processor set!"
!  Reduce size of processor grid
      pcol=pcol/2
      prow=pcol
      nb=pcol*prow
      write(90,*)'Reduced processor set size!'
      write(90,*)prow
      write(90,*)pcol
      write(90,*)n
      endif
      if (((n/nb) < prow) .or. ((n/nb) < pcol)) then
      print *,"Problem size still too small for processor set!"
      stop 
      endif
!
!  Following arrays are not addressed for eigenvalue only runs
      ALLOCATE(ICLUSTR(1)) 
      ALLOCATE(GAP(1)) 
!
!     ALLOCATE(ICLUSTR(2*N*N)) 
!     ALLOCATE(GAP(N*N)) 
! 
      ALLOCATE(IFAIL(N)) 
! 
masterU=0
masterM=0
!
call MPI_COMM_SIZE(commU,nprocU,ierr)
call MPI_comm_rank(commU,myidU,ierr)
!
! Now for the multiUniverse within the Universe
!
call MPI_COMM_SIZE( commM, nprocM, ierr) 
call MPI_COMM_RANK( commM, myidM, ierr)
!
! Initialize blacs processor grid
 
      call blacs_pinfo   (me,procs)
      call blacs_get     (0, 0, icontxt)
      call blacs_gridinit(icontxt, 'R', prow, pcol)
      call blacs_gridinfo(icontxt, prow, pcol, myrow, mycol)
      myunit = 10+me
      write(myunit,*)"--------"
      write(myunit,*)"Output for processor ",me," to unit ",myunit
      write(myunit,*)"Proc ",me,": myrow, mycol in p-array is ", &
         myrow, mycol
 
! Construct local arrays
! Global structure:  matrix A of n rows and n columns
!                    matrix Z of n rows and n column
!    
      myArows = numroc(n, nb, myrow, 0, prow)
      myAcols = numroc(n, nb, mycol, 0, pcol)
      write(myunit,*)"Size of global array is ",n," x ",n
      write(myunit,*)"Size of block is        ",nb," x ",nb
      write(myunit,*)"Size of local array is  ",myArows," x ",myAcols
! Initialize local arrays     
      allocate(A(myArows,myAcols)) 
     allocate(A0(myArows,myAcols)) 
     allocate(Z(1,1)) 
!      allocate(Z(myArows,myAcols)) 
!  Z  are the ON eigenvectors; not referenced for ev only 'N" case

!  Distribute PSI=Psi3 
!
! Prepare Full State vector (PHI) for rankM = 0
!

If(rankM.eq.0.and.myidM.eq.0) PHI=psi
!
!
CALL MPI_BCAST(PHI,N,MPI_COMPLEX,masterU,commU,IERR)

 do i=1,n
   call g2L(i,n,prow,nb,iproc,myi)
   if (myrow==iproc) then
 do j=1,n
    call g2L(j,n,pcol,nb,jproc,myj)
    if (mycol==jproc) then
A0(myi,myj) =   PHI(i)*conjg(PHI(j))
    endif
 enddo
    endif
  enddo
!
!  Note A0 is the density matrix for the noiseless channel
!  It could be used, along with A, to evaluate the Fidelity
!
!write(77+rankM,*)'myidM,rankM,myidU=',myidM,rankM,myidU
!write(77+rankM,*)'A0=',A0
!deallocate(A,A0,Z)
 
A=A0
!
ntag=88
!
!  SEND PSI  for each noise group (rankM>0) to masterU receive as PHI2
!
 Do ng=1,TGROUPS-1   !   Loop over Groups rankM>0
 k=nprocM*ng
!If(myidU.eq.0)call MPI_RECV(PHI2,N,MPI_COMPLEX,k,ntag,commU,mpistat,ierr)
!If(myidU.eq.k)call MPI_SEND(psi,N,MPI_COMPLEX,masterU,ntag,commU,ierr)
! Use IRECV method:
If(myidU.eq.0)call MPI_IRECV(PHI2,N,MPI_COMPLEX,k,ntag,commU,req(ng),ierr)
If(myidU.eq.k)call MPI_SEND(psi,N,MPI_COMPLEX,masterU,ntag,commU,ierr)
If(myidU.eq.0) call MPI_Wait(req(ng),mpistat,ierr)
!
!  Distribute PHI2  for each noise group (rankM>0) to all commU
!
CALL MPI_BCAST(PHI2,N,MPI_COMPLEX,masterU,commU,IERR)
!
!   Add rankM>0 density matrix to full density matrix
!   stored on each local array.  Distributed density matrix
!
!write(87,*)'PHI2=',PHI2
 do i=1,n
   call g2L(i,n,prow,nb,iproc,myi)
   if (myrow==iproc) then
    do j=1,n
     call g2L(j,n,pcol,nb,jproc,myj)
     if (mycol==jproc) then
A(myi,myj) =  A(myi,myj) + PHI2(i)*conjg(PHI2(j))
     endif
    enddo
   endif
enddo
!write(87,*)'ng A=', ng,A

enddo  !  end of group loop
!
! Prepare array descriptors for ScaLAPACK 
 
      IDESA(1) = 1         ! descriptor type
      IDESA(2) = icontxt   ! blacs context
      IDESA(3) = n         ! global number of rows
      IDESA(4) = n         ! global number of columns
      IDESA(5) = nb        ! row block size
      IDESA(6) = nb        ! column block size
      IDESA(7) = 1         ! initial process row
      IDESA(8) = 1         ! initial process column
      IDESA(9) = myArows   ! leading dimension of local array

      do i=1,9
        IDESZ(i) = IDESA(i)
      enddo
 
! Call ScaLAPACK library routine
! Ask PCHEEVX to compute:  just eigenvalues 'N', for range 'V'
! Range of eigenvalues VL to VU (0.1 to 1.01)
!
!  Matrix A has global dimension N X N; Local=praw X pcol
!  ABSTOL and ORFAC set accuracy of orthonalization for vectors.
!
VL=0.000000001;VU=1.1
!VL=0.1;VU=1.01
ABSTOL= -1.0E+0
ORFAC=ABSTOL
!
!write(87,*)'before PCHEEV A=',A
! first carry out a workspace query
! use  lwork = -1
! use  lrwork = -1
! use  liwork = -1
      ALLOCATE (WORK(1))
      ALLOCATE (RWORK(1))
      ALLOCATE(IWORK(1)) 
CALL PCHEEVX('N','V','U',N,A,1,1,IDESA,VL,VU,13,-13, ABSTOL, &
                          M, NZ, W, ORFAC, Z, 1, 1, IDESZ,  &
                          WORK, -1, RWORK, -1, IWORK, &
                           -1, IFAIL, ICLUSTR,GAP, INFO )
! NOW SET DIMENSIONS
       LWORK = INT(ABS(WORK(1)))
       DEALLOCATE(WORK)
       ALLOCATE(WORK(LWORK))
      
       LRWORK = INT(ABS(RWORK(1)))
       DEALLOCATE(RWORK)
       ALLOCATE(RWORK(LRWORK))

       LIWORK =INT (ABS(IWORK(1)))
       DEALLOCATE(IWORK)
       ALLOCATE(IWORK(LIWORK))
!  NOW DO EIGENVALUE CALL
CALL PCHEEVX('N','V','U',N,A,1,1,IDESA,VL,VU,13,-13, ABSTOL, &
                          M, NZ, W, ORFAC, Z, 1, 1, IDESZ,  &
                          WORK, LWORK, RWORK, LRWORK, IWORK, &
                           LIWORK, IFAIL, ICLUSTR,GAP, INFO )

! Output: M=total nos ev found ->N
! Output: NZ=total nos ev computed ->M;  not referenced 
! Output: W(M) array of  eigenvalues
! write(87,*)'M=', M 
! write(87,*)'W=', W 
! Deallocate the local arrays
 
      deallocate(A, A0, Z)
 
! End blacs for processors that are used
      DEALLOCATE(ICLUSTR,IFAIL,IWORK) 
      DEALLOCATE(GAP)
      DEALLOCATE (RWORK,WORK)
 
      call blacs_gridexit(icontxt)
      call blacs_exit(1)
return 
END SUBROUTINE VONNP


 
! convert global index to local index in block-cyclic distribution
 
   subroutine g2L(i,n,np,nb,p,il)
 
   implicit none
   integer :: i    ! global array index, input
   integer :: n    ! global array dimension, input
   integer :: np   ! processor array dimension, input
   integer :: nb   ! block size, input
   integer :: p    ! processor array index, output
   integer :: il   ! local array index, output
   integer :: im1   
 
   im1 = i-1
   p   = mod((im1/nb),np)
   il  = (im1/(np*nb))*nb + mod(im1,nb) + 1
   
   return
   end subroutine g2L 
 
! convert local index to global index in block-cyclic distribution
 
   subroutine L2g(il,p,n,np,nb,i)
 
   implicit none
   integer :: il   ! local array index, input
   integer :: p    ! processor array index, input
   integer :: n    ! global array dimension, input
   integer :: np   ! processor array dimension, input
   integer :: nb   ! block size, input
   integer :: i    ! global array index, output
   integer :: ilm1   
 
   ilm1 = il-1
   i    = (((ilm1/nb) * np) + p)*nb + mod(ilm1,nb) + 1
   
   return
   end subroutine L2g 

FUNCTION GCD(NA, NB)
        IA = NA
        IB = NB
1   IF (IB.NE.0) THEN
    ITEMP = IA
     IA = IB
     IB = MOD(ITEMP, IB)
   GOTO 1
   END IF
   GCD = IA
   RETURN
   END FUNCTION GCD

!SH  is the sign in HALL ============================
!

FUNCTION SH(nq,n,np)

Use MPI
Use MPI_VARS

Implicit None
INTEGER nq,n,np,i, S,SH
! In n X np Hadamard sign is determined
! By number of times bits are equal and equal 1
Logical :: bool1
Logical :: bool2 
S=1
do i=0,nq-1
bool1 = btest(n, i)
bool2 = btest(np, i)
If(bool1) then
If(bool1.eqv.bool2) S=-S
endif
end do
SH=S
END FUNCTION SH
!
!bottom of SH ============================

!=============================
FUNCTION f(xguess,j,MP)
!f(j+1)=PowerMod(xguess,j,MP)
IMPLICIT none
INTEGER f, xguess,j,MP,Powermod
f=PowerMod(xguess,j,MP)
END FUNCTION F

!=============================
Function PowerMod(Ia,Ik, N)
!  returns  Mod(a^k,N)
IMPLICIT none
INTEGER Powermod, Ia,Ik,N,IB,IC,I
        I = ik
        IB = 1
        IC=ia
    Do While(I.GT.0)
    Do While( Mod(I,2).eq.0)
    I=I/2
    IC= Mod(IC*IC,N)
    end  do
    I=I-1
    IB=IB*IC
    IB=MOD(IB,n)
    end do

   If(Ik.eq.0.and.Ia.gt.0) IB=1
   PowerMod = IB
   RETURN
 END FUNCTION POWERMOD

End Module qcmpisubs

