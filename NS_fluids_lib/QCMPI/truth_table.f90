Program main
Use mpi_vars
Use qcmpisubs
Implicit None

INTEGER nq,NBITS

     ierr=0
     nq=5
     nprocU=1
     myidU=0
     ierr=0
     nprocs=1
     myid=0
     nprocM=1
     myidM=0
     print *,"nq=",nq

  NBITS = 2**nq 

call truthtable(nq)
  
contains

Subroutine truthtable(nq)
Use MPI_VARS
Implicit None

DOUBLE PRECISION RAND(1:2,1:1) 
COMPLEX, PARAMETER   :: One  =( 1.0, 0.0 )
COMPLEX, PARAMETER   :: Zero  =( 0.0, 0.0 )

COMPLEX, ALLOCATABLE :: Psi1(:)

INTEGER nq,NBITS
INTEGER i,j,k
INTEGER NPART
INTEGER OUTFILE
!
! Multiuniverse variables
!
INTEGER rankM
COMPLEX :: Op(2,2)
integer :: is,is1,is2,B1(nq),B2(nq),commU

commU=0

NBITS=2**nq

!
! Outfile
!
OUTFILE=10

NPART=NBITS

Allocate(psi1(NPART))

OPEN(unit=OUTFILE)

write(OUTFILE,*)
write(OUTFILE,9003) nq
9003 Format('     Qubits           :   ', I4)
write(OUTFILE,9004) NPART
9004 Format('     W.f. array size  : ', I8)


do i=1,NPART

 do j=1,NPART
  Psi1(j)=Zero
 enddo
 Psi1(i)=One

 Op(1,1)=One
 Op(1,2)=Zero
 Op(2,1)=Zero
 Op(2,2)=One
 is=1
 call OneOpA(nq,is,Op,psi1,NPART,commU)
 is1=2
 is2=3
 call CNOTA(nq,is1,is2,psi1,NPART)
 is1=4
 is2=5
 k=1
 call cphasek(nq,is1,is2,k,psi1,NPART,commU)
 call dectobin(nq,i-1,B1)
 do k=1,NPART
  if (psi1(k).ne.Zero) then
   call dectobin(nq,k-1,B2)
   write(OUTFILE,9030) B1(1),B1(2),B1(3),B1(4),B1(5), &
    B2(1),B2(2),B2(3),B2(4),B2(5),psi1(k)
   write(OUTFILE,9040)
  endif
 enddo
   

9030 Format(I4,I4,I4,I4,I4,' -> ',I4,I4,I4,I4,I4,(F6.2,2X,F6.2,"i")) 
9040 Format('----------------------------------------------------------')
enddo !i=1,NPART

close (unit=OUTFILE)


deAllocate(psi1)

End subroutine

end program main
