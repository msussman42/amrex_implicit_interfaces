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
INTEGER nq,NPART,n,np,npl,dest,seatn,OUTFILE
!integer SH
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

!=============================
FUNCTION f(xguess,j,MP)
!f(j+1)=PowerMod(xguess,j,MP)
IMPLICIT none
INTEGER f, xguess,j,MP
!integer Powermod
f=PowerMod(xguess,j,MP)
END FUNCTION F

End Module qcmpisubs

