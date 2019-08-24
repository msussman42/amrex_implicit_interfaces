MODULE GeneralClass
USE probcommon_module
use MOF_routines_module
use geometry_intersect_module
IMPLICIT NONE

! r1=radcen-radeps
! r2=radcen+radeps
! set radcen,radeps in main.F90
real(kind=8) :: radcen
real(kind=8) :: radeps ! thick:0.1d0  thin:0.005d0

real(kind=8),parameter :: afrac_eps = 1.0d-8
real(kind=8),parameter :: pi  = 4.0d0 * atan(1.0d0)

! set these vars in main.F90
integer :: radial_variation
integer :: dirichlet_pentafoil
integer :: dirichlet_annulus
real(kind=8) :: pentaeps ! thick: 0.2  thin: 0.05 penta_foil
! for probtype=16   thickness of the thermal layer
! 0.001d0 is thin.  0.03d0 is thicker.
real(kind=8) :: thermal_delta
! filament_test_type==0 for irregular material 2
! filament_test_type==1 for circular material 2
integer     :: filament_test_type
! --------------------------------------


integer,parameter      :: pcurve_num = 2000
real(kind=8)           :: pcurve_ls(2,pcurve_num+1)
real(kind=8)           :: pcurve_ls2(2,pcurve_num+1)
real(kind=8)           :: pcurve_rad(pcurve_num+1)

real(kind=8)           :: crit_dist(2)

real(kind=8)           :: space_partition

! polar solver Np=128 Mp=256  or,
! polar solver Np=256 Mp=512 
integer,parameter  :: Np=256
integer,parameter  :: Mp=512
real(kind=8),dimension(0:Np,0:Mp) :: upolar
real(kind=8)         :: r_polar(0:Np)
real(kind=8)         :: z_polar(0:Mp)
real(kind=8)         :: dr_polar,dz_polar
real(kind=8)         :: pcenter(2)
real(kind=8) :: rlo
real(kind=8) :: rhi
real(kind=8),parameter :: BC_T1=2.0d0
real(kind=8),parameter :: BC_T2=3.0d0


!//////////////////////////////////////
!    TYPE DEFINE
!/////////////////////////////////////
TYPE::Node
  REAL(KIND=8)       :: VAL(2)
  TYPE(Node),POINTER :: Next
END TYPE Node

TYPE::POINTS
    REAL(KIND=8)       :: VAL(2)
END TYPE

TYPE:: SIDETYPE
    TYPE(POINTS)       :: PT(2)
END TYPE 

type:: IFSEG
    integer             :: flag
    TYPE(POINTS)        :: PT(2)
END TYPE

TYPE::polygon
    INTEGER                  :: NodesNum
!    INTEGER                  :: TAG
    TYPE(POINTS)             :: CENTER
    TYPE(POINTS),ALLOCATABLE :: Nodes(:)
!    INTEGER                  :: nsig(4)         ! (1 for pos   -1 for neg)
    TYPE(SIDETYPE)           :: SIDES(4)
    TYPE(POINTS)             :: side_cen(4)
!    INTEGER                  :: csig(4)   
    TYPE(points)             :: centroid      
END TYPE

TYPE:: LINE
    INTEGER                :: TAG
    REAL(KIND=8)           :: VAL(3)
END TYPE

TYPE:: PLG_LIST
    INTEGER                       :: TAG
    INTEGER                       :: NUM
    TYPE(POLYGON),ALLOCATABLE     :: PLG(:) 
END TYPE



!TYPE::Polygon
!     INTEGER                :: NodesNum
!     TYPE(POINTS),ALLOCATABLE :: Nodes(:)
!     INTEGER                :: IDS(2) 
!END TYPE

contains
!////////////////////////////////////////////////
!----------------------------------------------
!    <SUBROUTINE> 4th-Order Runge-Kutta
!----------------------------------------------
SUBROUTINE f_RK(xk,yk,tk,tau,xkout,ykout)
IMPLICIT NONE

INTEGER,PARAMETER       :: N=1
INTEGER                 :: i
REAL(KIND=8),INTENT(IN) :: xk,yk,tk,tau
REAL(KIND=8),INTENT(OUT):: xKout,yKout
REAL(KIND=8)            :: hr
REAL(KIND=8),external   :: Uprescribe,Vprescribe
REAL(KIND=8)            :: k0,k1,k2,k3,l0,l1,l2,l3
REAL(KIND=8),ALLOCATABLE:: t(:),X(:),y(:)

ALLOCATE (x(N+1),y(N+1),t(N+1))

hr=tau/N

t(1)=tk
x(1)=xk
y(1)=yk


do i=1,N

K0=hr*Uprescribe(-1,-1,t(i),x(i),y(i))
l0=hr*Vprescribe(-1,-1,t(i),x(i),y(i))


k1=hr*Uprescribe(-1,-1,t(i)+0.5d0*hr,x(i)+0.5d0*k0,y(i)+0.5d0*l0)
l1=hr*Vprescribe(-1,-1,t(i)+0.5d0*hr,x(i)+0.5d0*k0,y(i)+0.5d0*l0)
k2=hr*Uprescribe(-1,-1,t(i)+0.5d0*hr,x(i)+0.5d0*k1,y(i)+0.5d0*l1)
l2=hr*Vprescribe(-1,-1,t(i)+0.5d0*hr,x(i)+0.5d0*k1,y(i)+0.5d0*l1)
k3=hr*Uprescribe(-1,-1,t(i)+hr,x(i)+k2,y(i)+l2)
l3=hr*Vprescribe(-1,-1,t(i)+hr,x(i)+k2,y(i)+l2)

x(i+1)=x(i)+1.0d0/6.0d0*(k0+2.0d0*k1+2.0d0*k2+k3)
y(i+1)=y(i)+1.0d0/6.0d0*(l0+2.0d0*l1+2.0d0*l2+l3)

t(i+1)=t(i)+hr

enddo

xKout=x(N+1)
yKout=y(N+1)


DEALLOCATE (X,Y,T) 

END SUBROUTINE f_RK

!-----------------------------------------------
!   <SUBROUTINE> polygoncopy
!----------------------------------------------
subroutine polygoncopy(p_id,plga,plgb)             
implicit none
integer, intent(in) :: p_id
TYPE(POLYGON),INTENT(IN)   ::PLGA
TYPE(POLYGON),INTENT(OUT)  ::PLGB
INTEGER                    :: i

if(plga%nodesnum .le. 2) then
  print *,"invalid polygon in polygon copy nodesum=",plga%nodesnum
  print *,"p_id= ",p_id
  stop
endif 

PLGB%NODESNUM=PLGA%NODESNUM

if(plgb%nodesnum .ne. 0) then
allocate(PLGB%nodes(PLGB%nodesnum))

 do i=1,plgb%nodesnum
   plgb%nodes(i)%VAL=plga%nodes(i)%VAL
 enddo

endif
!PLGA%NODESNUM=0
!deallocate(plga%nodes)

end subroutine polygoncopy

!------------------------------------------------
!     <subroutine> POLYGON DELETE
!------------------------------------------------
SUBROUTINE PLGDEL(plga)
IMPLICIT NONE

TYPE(POLYGON)    :: PLGA

if (plga%nodesnum .ne. 0) then
   PLGA%NODESNUM=0
   DEALLOCATE (PLGA%NODES)
endif

plga%nodesnum = 0

END SUBROUTINE PLGDEL
!---------------------------------------------------
subroutine outputplg(plg)
implicit none

type(polygon)        :: plg
integer              :: i

if (plg%nodesnum .eq. 0) then
   print *, "0 polygon"
else
    do i=1,plg%nodesnum
       print *, i, plg%nodes(i)%val
    enddo
endif

end subroutine
!------------------------------------------------
!    <SUBROUTINE> FINDINTERSECTION_VERTICAL
!------------------------------------------------
Subroutine FindIntersection(x,y,clc,zv)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: x,y
TYPE(POINTS),INTENT(OUT) :: zv
REAL(KIND=8),INTENT(IN):: clc

 zv%VAL(1)=clc
 zv%VAL(2)=(clc-y%VAL(1))/(x%VAL(1)-y%VAL(1))*x%VAL(2)+&
          &(clc-x%VAL(1))/(y%VAL(1)-x%VAL(1))*y%VAL(2)

return

END Subroutine FindIntersection
!------------------------------------------------
!     <SUBROUTINE> FINDINTERSECTION_HORIZONTAL
!------------------------------------------------
Subroutine FindIntersection_H(x,y,clc,zh)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: x,y
TYPE(POINTS),INTENT(OUT) :: zh
REAL(KIND=8),INTENT(IN):: clc

 zh%VAL(2)=clc
 zh%VAL(1)=(clc-y%VAL(2))/(x%VAL(2)-y%VAL(2))*x%VAL(1)+&
          &(clc-x%VAL(2))/(y%VAL(2)-x%VAL(2))*y%VAL(1)

return

END Subroutine FindIntersection_H
!----------------------------------------------------
!     <SUBROUTINE> DETERMINESIDE_VERTICAL
!----------------------------------------------------
!   CALL DETERMINESIDE_V(PLG%Nodes(i),PLG%Nodes(i+1),CL,PL,PR)
!-----------------------------------------------------------
SUBROUTINE DETERMINESIDE_V(a,b,clc,PL,PR)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: a,b
TYPE(POINTS)             :: Z
Real(kind=8),INTENT(IN)  :: clc
TYPE(Node),POINTER       :: PL,PR
!TYPE(Node),POINTER       :: PLHEAD,PRHEAD

IF(((a%VAL(1).GE. clc) .and. (b%VAL(1).GT. clc)).OR.&
     &((a%VAL(1).GT.clc).and. (b%VAL(1) .GE. clc))) THEN
   IF( maxval ( abs ( a%VAL - PR%VAL) ) < afrac_eps ) Then
!     ALLOCATE(PR%NEXT)
!     PR=> PR%NEXT
!     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ElSE
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ENDIF
ELSEIF (((a%VAL(1).LE. clc) .and. (b%VAL(1).LT. clc)).OR.&
     &((a%VAL(1).LT.clc).and. (b%VAL(1) .LE. clc))) THEN
    IF( maxval(abs( a%VAL - PL%VAL)) < afrac_eps ) Then
!    ALLOCATE(PL%NEXT)
!    PL=> PL%NEXT
!    PL%VAL=a%VAL
    ALLOCATE(PL%NEXT)
    PL=> PL%NEXT
    PL%VAL=b%VAL
 ElSE
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=a%VAL
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=b%VAL
   ENDIF 

ELSEIF((a%VAL(1).GT. clc) .and. (b%VAL(1).LT. clc)) THEN
  CALL FINDINTERSECTION(a,b,clc,z)
    IF( maxval(abs( a%VAL - PR%VAL)) < afrac_eps ) Then
!    ALLOCATE(PR%NEXT)
!    PR=>PR%NEXT
!    PR%VAL=a%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=z%VAL
    ELSE
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=a%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=z%VAL
    ENDIF
     
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=z%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=b%VAL
 

ELSEIF((a%VAL(1).LT.clc ).and. (b%VAL(1).GT.clc)) THEN
  CALL FINDINTERSECTION(a,b,clc,z)
    IF( maxval(abs( a%VAL - PL%VAL)) < afrac_eps ) Then
!    ALLOCATE(PL%NEXT)
!    PL=>PL%NEXT
!    PL%VAL=A%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=Z%VAL
    ELSE
    ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=A%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=Z%VAL
    ENDIF
 
   ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=Z%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=B%VAL

ELSEIF( abs(a%VAL(1) - clc) < afrac_eps  .and. abs( b%VAL(1) - clc) < afrac_eps )  THEN
    IF(a%VAL(2).GT.b%VAL(2))THEN
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=a%VAL
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=b%VAL    

    ELSEIF(a%VAL(2).LT. b%VAL(2))THEN
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL    
    ENDIF

ELSE
   !print *, 'ERROR5'

ENDIF
END SUBROUTINE DetermineSide_V

!###################################################
!#########################################SUBROUTINE HORIZONTAL_DETERMINESIDE
SUBROUTINE DetermineSide_H(a,b,CLC,PL,PR)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)   :: a,b
TYPE(POINTS)              :: Z
Real(kind=8),INTENT(IN)   :: clc
TYPE(Node),POINTER        :: PL,PR
!TYPE(Node),POINTER       :: PLHEAD,PRHEAD

IF(((a%VAL(2).GE. clc) .and. (b%VAL(2).GT. clc)).OR.&
     &((a%VAL(2).GT.clc).and. (b%VAL(2) .GE. clc))) THEN
   !IF(a%VAL(1).eq.PR%VAL(1) .and. a%VAL(2).eq.PR%VAL(2)) Then
   IF( maxval(abs(a%VAL - PR%VAL)) < afrac_eps  ) Then
!     ALLOCATE(PR%NEXT)
!     PR=> PR%NEXT
!     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ElSE
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ENDIF
ELSEIF (((a%VAL(2).LE. clc) .and. (b%VAL(2).LT. clc)).OR.&
     &((a%VAL(2).LT.clc).and. (b%VAL(2) .LE. clc))) THEN
   ! IF(a%VAL(1).eq.PL%VAL(1) .and. a%VAL(2).eq.PL%VAL(2)) Then
   IF( maxval(abs(a%VAL - PL%VAL)) < afrac_eps  ) Then
!    ALLOCATE(PL%NEXT)
!    PL=> PL%NEXT
!    PL%VAL=a%VAL
    ALLOCATE(PL%NEXT)
    PL=> PL%NEXT
    PL%VAL=b%VAL
 ElSE
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=a%VAL
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=b%VAL
   ENDIF 

ELSEIF((a%VAL(2).GT. clc) .and. (b%VAL(2).LT. clc)) THEN
  CALL FINDINTERSECTION_H(a,b,clc,z)
   ! IF(a%VAL(1).eq.PR%VAL(1) .and. a%VAL(2).eq.PR%VAL(2)) Then
   IF( maxval(abs(a%VAL - PR%VAL)) < afrac_eps  ) Then
!    ALLOCATE(PR%NEXT)
!    PR=>PR%NEXT
!    PR%VAL=a%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=z%VAL
    ELSE
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=a%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=z%VAL
    ENDIF
     
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=z%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=b%VAL
 

ELSEIF((a%VAL(2).LT.clc .and. b%VAL(2).GT.clc)) THEN
  CALL FINDINTERSECTION_H(a,b,clc,z)
   ! IF(a%VAL(1).eq.PL%VAL(1) .and. a%VAL(2).eq.PL%VAL(2)) Then
   IF( maxval(abs(a%VAL - PL%VAL)) < afrac_eps  ) Then
!    ALLOCATE(PL%NEXT)
!    PL=>PL%NEXT
!    PL%VAL=A%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=Z%VAL
    ELSE
    ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=A%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=Z%VAL
    ENDIF
 
   ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=Z%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=B%VAL

ELSEIF( abs(a%VAL(2) - clc) < afrac_eps .and. abs(b%VAL(2) - clc) < afrac_eps ) THEN
    IF(a%VAL(1).GT.b%VAL(1))THEN
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL    

    ELSEIF(a%VAL(1).LT. b%VAL(1))THEN
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=a%VAL
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=b%VAL    
    ENDIF

ELSE
   !print *, 'ERROR2'
ENDIF

END SUBROUTINE DetermineSide_H


!-----------------------------------------------------
!         <SUBROUTINE> CUTPOLYGON_VERTICAL
!-------------------------------------------------------
SUBROUTINE CUTPOLYGON_V(PLG,CL,PLGL,PLGR)        !SUBROUTINE CUTPOLYGON_V
IMPLICIT NONE

TYPE(POLYGON)             :: PLG
TYPE(POLYGON),INTENT(OUT) :: PLGL,PLGR
REAL(KIND=8),INTENT(IN)   :: CL
TYPE(Node),POINTER        :: PLHEAD,PL,PR,PRHEAD
INTEGER                   :: i,N
REAL(KIND=8)              :: P(2)

ALLOCATE(PRHEAD)
ALLOCATE(PLHEAD)

PLhead%val(1) = 1.0e+8
PLhead%val(2) = 1.0e+8
PRhead%val(1) = 1.0e+8
PRhead%val(2) = 1.0e+8

PL => PLhead
PR => PRhead

DO i=1,PLG%NodesNum-1                  !not allocate memory PLG!
   CALL DETERMINESIDE_V(PLG%Nodes(i),PLG%Nodes(i+1),CL,PL,PR)
ENDDO
   CALL DETERMINESIDE_V(PLG%Nodes(PLG%NODESNUM),PLG%Nodes(1),CL,PL,PR)

nullify(pl%next,pr%next)


!!!!!!!!!!!!!!!COUNT!!!!!!!!!!!!!!!!!!!1

PR=>PRHEAD
PR=>PR%NEXT
P=PR%VAL
N=1
do while(ASSOCIATED(PR%NEXT))
   PR=>PR%NEXT
   N=N+1
ENDDO
!IF((P(1).EQ.PR%VAL(1)).AND.(P(2).EQ.PR%VAL(2))) THEN
IF( maxval(abs(P - PR%VAL)) < afrac_eps ) THEN
   PLGR%NODESNUM=N-1
ELSE 
   PLGR%NODESNUM=N
ENDIF



!!********************

PL=>PLHEAD
PL=>PL%NEXT
P=PL%VAL
N=1

do while(ASSOCIATED(PL%NEXT))

   PL=>PL%NEXT
   N=N+1
ENDDO
!1F((P(1).EQ.PL%VAL(1)).AND.(P(2).EQ.PL%VAL(2))) THEN
IF( maxval(abs(P-PL%VAL)) < afrac_eps  ) then 

  PLGL%NODESNUM=N-1
ELSE 
   PLGL%NODESNUM=N
ENDIF

ALLOCATE(PLGL%NODES(PLGL%NODESNUM),PLGR%NODES(PLGR%NODESNUM))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PL=>PLHEAD
do i=1,PLGL%NODESNUM
   PL=>PL%NEXT
   PLGL%NODES(i)%VAL=PL%VAL
ENDDO


PR=>PRHEAD
do i=1,PLGR%NODESNUM
   PR=>PR%NEXT
   PLGR%NODES(i)%VAL=PR%VAL
ENDDO

                                 !NOT DEALLOCATE PLGL PLGR MEMORY YET!
!NULLIFY(PL,PR)
NULLIFY(PL%NEXT,PR%NEXT)
NULLIFY(PLHEAD,PRHEAD)
!DEALLOCATE(PLG%NODES)
CALL PLGDEL(PLG)
END SUBROUTINE CUTPOLYGON_V

!---------------------------------------------------
!       <SUBROUTINE> CUTPOLYGON_HORIZONTAL
!---------------------------------------------------
SUBROUTINE CUTPOLYGON_H(PLG,CL,PLGL,PLGR) !SUBROUTINE CUTPOLYGON_H
IMPLICIT NONE

TYPE(POLYGON)             :: PLG
TYPE(POLYGON),INTENT(OUT) :: PLGL,PLGR
REAL(KIND=8),INTENT(IN)   :: CL
TYPE(Node),POINTER        :: PLHEAD1,PL1,PR1,PRHEAD1
INTEGER                   :: i,N
REAL(KIND=8)              :: P(2)

ALLOCATE(PRHEAD1)
ALLOCATE(PLHEAD1)

PLhead1%val(1) = 1.0e+8
PLhead1%val(2) = 1.0e+8
PRhead1%val(1) = 1.0e+8
PRhead1%val(2) = 1.0e+8

PL1 => PLhead1
PR1 => PRhead1

DO i=1,PLG%NodesNum-1                  !not allocate memory PLG!
   CALL DETERMINESIDE_H(PLG%Nodes(i),PLG%Nodes(i+1),CL,PL1,PR1)
ENDDO
   CALL DETERMINESIDE_H(PLG%Nodes(PLG%NODESNUM),PLG%Nodes(1),CL,PL1,PR1)

nullify(pl1%next,pr1%next)

!!!!!!!!!!!!!!!COUNT!!!!!!!!!!!!!!!!!!!1

PR1=>PRHEAD1
PR1=>PR1%NEXT
P=PR1%VAL
N=1
do while(ASSOCIATED(PR1%NEXT))
   PR1=>PR1%NEXT
   N=N+1
ENDDO
!IF( abs( P(1) - PR1%VAL(1)) < afrac_eps .AND. (P(2).EQ.PR1%VAL(2))) THEN

IF( maxval(abs(P - PR1%VAL)) < afrac_eps ) THEN
   PLGR%NODESNUM=N-1
ELSE 
   PLGR%NODESNUM=N
ENDIF


PL1=>PLHEAD1
PL1=>PL1%NEXT
P=PL1%VAL
N=1
do while(ASSOCIATED(PL1%NEXT))
   PL1=>PL1%NEXT
   N=N+1
ENDDO

!IF((P(1).EQ.PL1%VAL(1)).AND.(P(2).EQ.PL1%VAL(2))) THEN
IF( maxval(abs(P - PL1%VAL)) < afrac_eps ) THEN
   PLGL%NODESNUM=N-1
ELSE 
   PLGL%NODESNUM=N
ENDIF

ALLOCATE(PLGL%NODES(PLGL%NODESNUM),PLGR%NODES(PLGR%NODESNUM))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PL1=>PLHEAD1
do i=1,PLGL%NODESNUM
   PL1=>PL1%NEXT
   PLGL%NODES(i)%VAL=PL1%VAL
ENDDO


PR1=>PRHEAD1
do i=1,PLGR%NODESNUM
   PR1=>PR1%NEXT
   PLGR%NODES(i)%VAL=PR1%VAL
ENDDO
                                 !NOT DEALLOCATE PLGL PLGR MEMORY YET!
!NULLIFY(PL1,PR1)
NULLIFY(PL1%NEXT,PR1%NEXT)
NULLIFY(PLHEAD1,PRHEAD1)
!DEALLOCATE(PLG%NODES)
CALL PLGDEL(PLG)
END SUBROUTINE CUTPOLYGON_H
!--------------------------------------------------------


subroutine PinL(mx,my,alpha,X,VOUT)
implicit none

REAL(KIND=8),INTENT(IN)   :: mx,my,alpha
TYPE(POINTS)              :: x
REAL(KIND=8)              :: VOUT

  VOUT= mx*x%val(1)+my*x%val(2)-alpha

end subroutine PinL
!----------------------------------------------------------
subroutine LxL(a,b,mx,my,alpha,z)
implicit none

type(points),intent(in)  :: a,b
real(kind=8),intent(in)  :: mx,my,alpha
real(kind=8)             :: slope

type(points),intent(out) :: z

slope = 0.0d0

if(abs(a%val(1) - b%val(1)) .lt. afrac_eps) then
   z%val(1) = a%val(1)
   z%val(2) = -(mx*z%val(1) - alpha)/my
elseif(abs(a%val(2) - b%val(2)) .lt. afrac_eps) then
   z%val(2) = a%val(2)
   z%val(1) = -(my*z%val(2)-alpha)/mx
else
   slope = (b%val(2)-a%val(2))/(b%val(1)-a%val(1))
   z%val(1) = (my*(slope*a%val(1)-a%val(2))+alpha)/(mx+my*slope)
   z%val(2) = slope*(z%val(1)-a%val(1))+a%val(2)
endif

end subroutine LxL
!--------------------------------------------------------
SUBROUTINE DETERMINESIDE(a,b,mx,my,alpha,PL,PR)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: a,b
real(kind=8),intent(in)  :: mx,my,alpha
real(kind=8)             :: nx,ny,nalpha


TYPE(POINTS)             :: Z
TYPE(Node),POINTER       :: PL,PR
real(kind=8)             :: aval,bval


nx = mx
ny = my
nalpha=alpha


if (nx .lt. 0.0d0) then
   nx = -1.0d0*nx
   ny = -1.0d0*ny
   nalpha = -1.0d0*nalpha
endif

if(nx .eq. 0.0d0  .and. ny .lt. 0.0d0) then
   ny = -1.0d0*ny
   nalpha = -1.0d0*nalpha
endif



call PinL(nx,ny,nalpha,a,aval)
call PinL(nx,ny,nalpha,b,bval)


IF(((aval.GE. 0.0d0) .and. (bval.GT. 0.0d0)).OR.&
     &((aval .GT. 0.0d0).and. (bval .GE. 0.0d0))) THEN
   IF( maxval ( abs ( a%VAL - PR%VAL) ) < afrac_eps ) Then
!     ALLOCATE(PR%NEXT)
!     PR=> PR%NEXT
!     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ElSE
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ENDIF
ELSEIF (((aval .LE. 0.0d0) .and. (bval.LT. 0.0d0)).OR.&
     &((aval .LT. 0.0d0).and. (bval .LE. 0.0d0))) THEN
    IF( maxval(abs( a%VAL - PL%VAL)) < afrac_eps ) Then
!    ALLOCATE(PL%NEXT)
!    PL=> PL%NEXT
!    PL%VAL=a%VAL
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=b%VAL
    ElSE
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=a%VAL
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=b%VAL
    ENDIF 

ELSEIF((aval.GT. 0.0d0) .and. (bval .LT. 0.0d0)) THEN
  CALL LxL(a,b,mx,my,alpha,z)
    IF( maxval(abs( a%VAL - PR%VAL)) < afrac_eps ) Then
!    ALLOCATE(PR%NEXT)
!    PR=>PR%NEXT
!    PR%VAL=a%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=z%VAL
    ELSE
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=a%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=z%VAL
    ENDIF
     
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=z%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=b%VAL
 

ELSEIF((aval.LT. 0.0d0 ).and. (bval.GT. 0.0d0)) THEN
  CALL LxL(a,b,mx,my,alpha,z)
    IF( maxval(abs( a%VAL - PL%VAL)) < afrac_eps ) Then
!    ALLOCATE(PL%NEXT)
!    PL=>PL%NEXT
!    PL%VAL=A%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=Z%VAL
    ELSE
    ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=A%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=Z%VAL
    ENDIF
 
   ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=Z%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=B%VAL
!------------------------------------------------------------
ELSEIF( aval .eq. 0.0d0  .and. bval .eq. 0.0d0 )  THEN
  if (ny .eq. 0.0d0 .and. nx .ne. 0.0d0) then
    
   IF(a%val(2) .GT. b%val(2))THEN
       IF( maxval ( abs ( a%VAL - PR%VAL) ) < afrac_eps ) Then
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL
       else
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=a%VAL
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL    
       endif
   ELSEIF(a%val(2) .LT. b%val(2))THEN
     IF( maxval ( abs ( a%VAL - PL%VAL) ) < afrac_eps ) Then
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL 
     else 
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL
     endif  
   else
     write(*,*) 'cutpolygon:ERROR1'
   endif


 elseif(nx .eq. 0.0d0 .and. ny .ne. 0.0d0 ) then
     IF(a%val(1) .lT. b%val(1))THEN
       IF( maxval ( abs ( a%VAL - PR%VAL) ) < afrac_eps ) Then
                ALLOCATE(PR%NEXT)
                PR=>PR%NEXT
                PR%VAL=b%VAL
       else
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=a%VAL
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL    
       endif
     ELSEIF(a%val(1) .gT. b%val(1))THEN
      IF( maxval ( abs ( a%VAL - PL%VAL) ) < afrac_eps ) Then
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL 
      else 
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL
      endif  
     else
       write(*,*) 'cutpolygon: ERROR2'
     endif
 elseif(nx .ne. 0.0d0 .and. ny .ne. 0.0d0)then
   if(nx*ny .gt. 0.0d0) then
     IF(a%val(1) .lT. b%val(1))THEN
       IF( maxval ( abs ( a%VAL - PR%VAL) ) < afrac_eps ) Then
                ALLOCATE(PR%NEXT)
                PR=>PR%NEXT
                PR%VAL=b%VAL
       else
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=a%VAL
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL    
       endif
     ELSEIF(a%val(1) .gT. b%val(1))THEN
      IF( maxval ( abs ( a%VAL - PL%VAL) ) < afrac_eps ) Then
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL 
      else 
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL
      endif  
     else
       write(*,*) 'cutpolygon: ERROR2'
     endif
   ELSE
     IF(a%val(1) .GT. b%val(1))THEN
       IF( maxval ( abs ( a%VAL - PR%VAL) ) < afrac_eps ) Then
                ALLOCATE(PR%NEXT)
                PR=>PR%NEXT
                PR%VAL=b%VAL
       else
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=a%VAL
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL    
       endif
     ELSEIF(a%val(1) .LT. b%val(1))THEN
      IF( maxval ( abs ( a%VAL - PL%VAL) ) < afrac_eps ) Then
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL 
      else 
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL
      endif  
     else
       write(*,*) 'cutpolygon: ERROR2'
     endif
   endif
 else
   print *, "cutpolygon: nx = 0, ny = 0"
   
ENDIF
!------------------------------------------------------------------
ELSE
   print *, 'ERROR3'
ENDIF

END SUBROUTINE DetermineSide


!--------------------------------------------------------
subroutine cutpolygon(a,b,alpha,PLG,PLGL,PLGR)
!----------------------------------------------------------
!   ax+by = alpha
!--------------------------------------------------------
implicit none

REAL(KIND=8),intent(in)    :: a,b,alpha
type(polygon)              :: plg
type(polygon),intent(out)  :: plgl,plgr

TYPE(Node),POINTER         :: PLHEAD1,PL1,PR1,PRHEAD1
INTEGER                    :: i,N
REAL(KIND=8)               :: P(2)

!print *, "mx=", a, "my=",b,"alpha=",alpha

ALLOCATE(PRHEAD1)
ALLOCATE(PLHEAD1)

PLhead1%val(1) = 1.0e+8
PLhead1%val(2) = 1.0e+8
PRhead1%val(1) = 1.0e+8
PRhead1%val(2) = 1.0e+8

PL1 => PLhead1
PR1 => PRhead1


if (a .eq. 0.0d0 .and. b .eq. 0.0d0) then
WRITE(*,*) " a and b both 0"
endif

DO i=1,PLG%NodesNum-1                  !not allocate memory PLG!
   CALL DETERMINESIDE(PLG%Nodes(i),PLG%Nodes(i+1),a,b,alpha,PL1,PR1)
ENDDO
   CALL DETERMINESIDE(PLG%Nodes(PLG%NODESNUM),PLG%Nodes(1),a,b,alpha,PL1,PR1)


nullify(pl1%next,pr1%next)

!!!!!!!!!!!!!!!COUNT!!!!!!!!!!!!!!!!!!!1

PR1=>PRHEAD1
if (.not. associated(pr1%next)) then
  PLGR%NODESNUM = 0
else
PR1=>PR1%NEXT
P=PR1%VAL
N=1
do while(ASSOCIATED(PR1%NEXT))
   PR1=>PR1%NEXT
   N=N+1
ENDDO
!IF( abs( P(1) - PR1%VAL(1)) < afrac_eps .AND. (P(2).EQ.PR1%VAL(2))) THEN

IF( maxval(abs(P - PR1%VAL)) < afrac_eps ) THEN
   PLGR%NODESNUM=N-1
ELSE 
   PLGR%NODESNUM=N
ENDIF

endif


PL1=>PLHEAD1
if(.not. associated(PL1%next)) then
   PLGL%NODESNUM = 0 
ELSE
PL1=>PL1%NEXT
P=PL1%VAL
N=1
do while(ASSOCIATED(PL1%NEXT))
   PL1=>PL1%NEXT
   N=N+1
ENDDO

!IF((P(1).EQ.PL1%VAL(1)).AND.(P(2).EQ.PL1%VAL(2))) THEN
IF( maxval(abs(P - PL1%VAL)) < afrac_eps ) THEN
   PLGL%NODESNUM=N-1
ELSE 
   PLGL%NODESNUM=N
ENDIF

endif

!-------------------------------------------------
IF(PLGL%NODESNUM .ne. 0) then
ALLOCATE(PLGL%NODES(PLGL%NODESNUM))
PL1=>PLHEAD1
do i=1,PLGL%NODESNUM
   PL1=>PL1%NEXT
   PLGL%NODES(i)%VAL=PL1%VAL
ENDDO

endif

if(PLGR%NODESNUM .ne. 0) then
ALLOCATE(PLGR%NODES(PLGR%NODESNUM))
PR1=>PRHEAD1
do i=1,PLGR%NODESNUM
   PR1=>PR1%NEXT
   PLGR%NODES(i)%VAL=PR1%VAL
ENDDO

endif 

                             !NOT DEALLOCATE PLGL PLGR MEMORY YET!
!NULLIFY(PL1,PR1)
NULLIFY(PL1%NEXT,PR1%NEXT)
NULLIFY(PLHEAD1,PRHEAD1)
!DEALLOCATE(PLG%NODES)

end subroutine cutpolygon

!//////////////////////////////////////////////////////////////
SUBROUTINE l_DetermineSide(alpha,beta,zeta,a,b,PL,PR)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: a,b
TYPE(POINTS)             :: Z
TYPE(NODE),POINTER       :: PL,PR
real(kind=8)             :: vala,valb
real(kind=8),intent(in)  :: alpha, beta, zeta

call cut_DIFF(alpha,beta,zeta,a,vala)
call cut_DIFF(alpha,beta,zeta,b,valb)


IF(((vala .GE. 0) .and. (valb.GT. 0)).OR.&
     &((vala.GT. 0).and. (valb .GE. 0))) THEN
   IF( maxval ( abs ( a%VAL - PR%VAL) ) < afrac_eps ) Then
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ElSE
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ENDIF
ELSEIF (((vala.LE. 0) .and. (valb.LT. 0)).OR.&
     &((vala.LT. 0).and. (valb .LE. 0))) THEN
    IF( maxval(abs( a%VAL - PL%VAL)) < afrac_eps ) Then
      ALLOCATE(PL%NEXT)
      PL=> PL%NEXT
      PL%VAL=b%VAL
    ElSE
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=a%VAL
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=b%VAL
   ENDIF 

ELSEIF((vala.GT. 0) .and. (valb.LT. 0)) THEN
  CALL l_FINDINTERSECTION(alpha,beta,zeta,a,b,z)                           !!!!!!!!!!!!!!!!
    IF( maxval(abs( a%VAL - PR%VAL)) < afrac_eps ) Then
     ALLOCATE(PR%NEXT)
     PR=>PR%NEXT
     PR%VAL=z%VAL
    ELSE
     ALLOCATE(PR%NEXT)
     PR=>PR%NEXT
     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=>PR%NEXT
     PR%VAL=z%VAL
    ENDIF
     
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=z%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=b%VAL
 

ELSEIF((vala .LT. 0 ).and. (valb .GT. 0)) THEN
  CALL l_FINDINTERSECTION(alpha,beta,zeta,a,b,z)
   IF( maxval(abs( a%VAL - PL%VAL)) < afrac_eps ) Then
      ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=Z%VAL
   ELSE
      ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=A%VAL
      ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=Z%VAL
   ENDIF
 
   ALLOCATE(PR%NEXT)
   PR=>PR%NEXT
   PR%VAL=Z%VAL
   ALLOCATE(PR%NEXT)
   PR=>PR%NEXT
   PR%VAL=B%VAL

ELSEIF(vala .EQ. 0 .and. valb .EQ. 0) THEN
    IF(a%VAL(2).GT.b%VAL(2))THEN
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=a%VAL
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=b%VAL    
    ELSEIF(a%VAL(2).LT. b%VAL(2))THEN
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL    
    ENDIF

ELSE
   !print *, 'ERROR4'

ENDIF
END SUBROUTINE l_DetermineSide

!----------------------------------------------------------------
SUBROUTINE l_CUTPOLYGON(PLG,a,b,d,PLGL,PLGR)             
IMPLICIT NONE

TYPE(POLYGON)             :: PLG
TYPE(POLYGON),INTENT(OUT) :: PLGL,PLGR
TYPE(Node),POINTER        :: PLHEAD,PL,PR,PRHEAD
INTEGER                   :: i,N
REAL(KIND=8)              :: P(2)
REAL(KIND=8),INTENT(IN)   :: a,b,d

plg%nodesnum = 4

ALLOCATE(PRHEAD)
ALLOCATE(PLHEAD)

PLhead%val(1) = 1.0e+8
PLhead%val(2) = 1.0e+8
PRhead%val(1) = 1.0e+8
PRhead%val(2) = 1.0e+8

PL => PLhead
PR => PRhead



DO i=1,PLG%NodesNum-1                  !not allocate memory PLG!
   CALL L_DETERMINESIDE(a,b,d,PLG%Nodes(i),PLG%Nodes(i+1),PL,PR)   
ENDDO
   CALL L_DETERMINESIDE(a,b,d,PLG%Nodes(PLG%NODESNUM),PLG%Nodes(1),PL,PR)
NULLIFY(PL%NEXT,PR%NEXT)

!--------count--------------------

PR=>PRHEAD
if (.not. associated(pr%next)) then
  PLGR%NODESNUM = 0
else
PR=>PR%NEXT
P=PR%VAL
N=1
do while(ASSOCIATED(PR%NEXT))
   PR=>PR%NEXT
   N=N+1
ENDDO
!IF( abs( P(1) - PR1%VAL(1)) < afrac_eps .AND. (P(2).EQ.PR1%VAL(2))) THEN

IF( maxval(abs(P - PR%VAL)) < afrac_eps ) THEN
   PLGR%NODESNUM=N-1
ELSE 
   PLGR%NODESNUM=N
ENDIF

endif
!-------------------------------
PL=>PLHEAD
if(.not. associated(PL%next)) then
   PLGL%NODESNUM = 0 
ELSE
PL=>PL%NEXT
P=PL%VAL
N=1
do while(ASSOCIATED(PL%NEXT))
   PL=>PL%NEXT
   N=N+1
ENDDO

!IF((P(1).EQ.PL1%VAL(1)).AND.(P(2).EQ.PL1%VAL(2))) THEN
IF( maxval(abs(P - PL%VAL)) < afrac_eps ) THEN
   PLGL%NODESNUM=N-1
ELSE 
   PLGL%NODESNUM=N
ENDIF

endif
!--------------------------------
IF(PLGL%NODESNUM .NE. 0) THEN
  ALLOCATE(PLGL%NODES(PLGL%NODESNUM))
  PL=>PLHEAD
  do i=1,PLGL%NODESNUM
    PL=>PL%NEXT
    PLGL%NODES(i)%VAL=PL%VAL
  ENDDO
ENDIF

IF(PLGR%NODESNUM .NE. 0) THEN
  ALLOCATE(PLGR%NODES(PLGR%NODESNUM))
  PR=>PRHEAD
  do i=1,PLGR%NODESNUM
    PR=>PR%NEXT
    PLGR%NODES(i)%VAL=PR%VAL
  ENDDO
ENDIF

                                 !NOT DEALLOCATE PLGL PLGR MEMORY YET!
NULLIFY(PL%NEXT,PR%NEXT)
NULLIFY(PLHEAD,PRHEAD)

END SUBROUTINE l_CUTPOLYGON

Subroutine L_FindIntersection(a,b,d,x1,x2,zout)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: x1,x2
TYPE(POINTS),INTENT(OUT) :: zout
REAL(KIND=8)             :: m
REAL(KIND=8),INTENT(IN)  :: a,b,d

if (x2%val(1) .ne. x1%val(1)) THEN
   m= (x2%val(2)-x1%val(2))/(x2%val(1)-x1%val(1))

   IF((b .NE. 0.0d0) .and. (a .NE. 0.0d0)) then
     zout%val(1)= (m*x1%val(1)-x1%val(2)-d/b)/(m+a/b)
     zout%val(2)= -a/b*zout%val(1)-d/b

   ELSEIF((b .EQ. 0.0d0) .and. (a .NE. 0.0d0)) then
    zout%val(1)= -d/a
    zout%val(2)= x1%val(2)+m*(-d/a-x1%val(1)) 
   ELSE 
    zout%val(2)= -d/b
    zout%val(1)= (-d/b-x1%val(2))/m+x1%val(1) 
   ENDIF
ELSE
   zout%val(1)= x1%val(1)
   zout%val(2)= (-d-a*x1%val(1))/b
   
ENDIF

END Subroutine l_FindIntersection

!////////////////////////////////////////////////////////////////////////

!----------------------------------------------------------------
!            <SUBROUTINE> DETERMINE_RANGE
!----------------------------------------------------------------
SUBROUTINE DetermineRange(POLYGONIN,X,Y,N,M,il,ir,id,iu)
IMPLICIT NONE
TYPE(POLYGON),INTENT(IN)              :: POLYGONIN
REAL(KIND=8),INTENT(IN)               :: X(N+1),Y(M+1)
REAL(KIND=8)                          :: CLD,CLU,CLL,CLR
REAL(KIND=8),ALLOCATABLE              :: TRANS(:)
INTEGER,INTENT(OUT)                   :: il,ir,iu,id         
INTEGER                               :: i,N,M

print *,"i=0,N-1 and j=0,N-1"
stop

ALLOCATE(TRANS(POLYGONIN%NODESNUM))
do i=1,polygonin%nodesnum
   trans(i)=polygonin%nodes(i)%val(2)
ENDDO  
!do i=1,polygonin%nodesnum
!    print*, trans(i)
!enddo
!!print *, trans
 CLD=MINVAL(trans,POLYGONIN%NODESNUM)
 CLU=MAXVAL(trans,POLYGONIN%NODESNUM)
!print *, "down=",cld,"up=", clu
!DEALLOCATE(TRANS)

!ALLOCATE(TRANS(POLYGONIN%NODESNUM)) 
do i=1,polygonin%nodesnum
   trans(i)=polygonin%nodes(i)%val(1)
ENDDO  
 CLL=MINVAL(trans,POLYGONIN%NODESNUM)
 CLR=MAXVAL(trans,POLYGONIN%NODESNUM)
!print *, "left=",cll,"right=", clr
DEALLOCATE(TRANS)

print *,"i=0,M-1 now ?"
stop

do i=1,M
   if ((CLD .ge. (Y(i)-afrac_eps)) .and. (CLD .LT. Y(i+1))) then
      id=i 
      exit  
   ENDIF
enddo

print *,"i=0,M-1 now ?"
stop

do i=1,M
   if (CLU.GT. Y(i) .and. CLU.LE.(Y(i+1)+afrac_eps)) then
      iu=i+1
      exit
   ENDIF
!print*,iu
enddo

print *,"i=0,N-1 now ?"
stop

do i=1,N
   if (CLL.ge. (x(i)-afrac_eps) .and. CLL .LT. X(i+1) )  then
      il=i
      exit
   ENDIF
enddo

print *,"i=0,N-1 now ?"
stop

do i=1,N
   if (CLR.GT. X(i) .and.  CLR .LE.(x(i+1)+afrac_eps)) then
      ir=i+1
      exit
   ENDIF
enddo
!!print *, il,ir,id,iu

END SUBROUTINE DetermineRange

!-------------------------------------------------------------------
!                    MESH_CUT   
!-------------------------------------------------------------------
SUBROUTINE HorizontalMeshCut(polygonIn,X,Y,NX,MY,polygontemp,id,iu)
IMPLICIT NONE
!allocate x(:) y(:) before CALL

TYPE(POLYGON)                         :: POLYGONIN
TYPE(POLYGON),ALLOCATABLE,Intent(out) :: POLYGONTEMP(:)
TYPE(POLYGON)                         :: PLGL,PLGR,PLGRTEMP
INTEGER                  ,INTENT(IN)  :: NX,MY
REAL(KIND=8)             ,INTENT(IN)  :: X(NX+1),Y(MY+1)
!INTEGER,ALLOCATABLE                  :: H_ID(:)
INTEGER                               :: RY
INTEGER                               :: i
INTEGER                               :: il,ir,iu,id
!INTEGER                              :: il1,ir1,iu1,id1
!INTEGER ,ALLOCATABLE                 :: YCOUNT(:)

!horizontal cut first
CALL DetermineRange(POLYGONIN,X,Y,NX,MY,il,ir,id,iu) 
!!print *, "bound of the mesh_cut(horizontal)", id, iu

RY= iu-id
!ALLOCATE(H_ID(RY))
ALLOCATE(POLYGONTEMP(RY))

IF(iu.eq.id+1) THEN
!  H_ID(RY)= id
  CALL polygoncopy(10,POLYGONIN,POLYGONTEMP(RY))
  CALL PLGDEL(POLYGONIN)
ELSE
     CALL polygoncopy(11,PolygonIn,PLGRTEMP)
     CALL PLGDEL(POLYGONIN)
     DO i=id+1,iu-1
        CALL CUTPOLYGON_H(PLGRTEMP,Y(i),PLGL,PLGR)
        CALL polygoncopy(12,PLGL,POLYGONTEMP(i-id))
        CALL PLGDEL(PLGL)
        CALL polygoncopy(13,PLGR,PLGRTEMP)
        CALL PLGDEL(PLGR)
     ENDDO
     CALL polygoncopy(14,PLGRTEMP,POLYGONTEMP(RY))
     CALL PLGDEL(PLGRTEMP)
ENDIF
 
!################test polygontemp###################
!  !print *, "After horizontal cut, we have:"
!do k=1,RY
!  !print *, "row= ",k
!do j=1,polygontemp(K)%nodesnum
!  !print *, polygontemp(k)%nodes(j)%VAL(1),polygontemp(k)%nodes(j)%VAL(2)
!  enddo
!enddo
!####################################################
END SUBROUTINE HorizontalMeshCut


!----------------------------------------------------------
SUBROUTINE VerticalMeshCut(PolygonIn, X,Y,NX,MY,il,ir,PolygonOut)

TYPE(POLYGON)                         :: PolygonIn
! YANG LIU: POINTER?
TYPE(POLYGON),INTENT(OUT),ALLOCATABLE :: PolygonOut(:)
TYPE(POLYGON)                         :: PLGL,PLGR,PLGRTEMP
! YANG LIU: PASS THE DIMENSIONS OF X and Y as a parameter.
! get rid of "ALLOCATABLE" (repeat for other violating cases)
INTEGER                  ,INTENT(IN)  :: NX,MY
REAL(KIND=8)             ,INTENT(IN)  :: X(NX+1),Y(MY+1)
INTEGER                               :: il,ir,id,iu
INTEGER                               :: RX
INTEGER                               :: j


!(plgs(i,j)%plg(ii)ertical cut

  !print *, "After vertical cut, we have:"
  CALL DETERMINERANGE(PolygonIn,X,Y,NX,MY,il,ir,id,iu)
  RX= ir-il
  !print *, "left and right bound for row",i,"is :", il, ir
  ALLOCATE(PolygonOut(RX))
  IF(ir.eq.il+1) THEN
    CALL polygoncopy(15,PolygonIn,polygonout(1)) 
    CALL PLGDEL(POLYGONIN)
  ELSE      
    CALL polygoncopy(16,POLYGONIN,PLGRTEMP)
    CALL PLGDEL(POLYGONIN)
    do j=il+1,ir-1
       CALL CUTPOLYGON_V(PLGRTEMP,X(j),PLGL,PLGR)
!!print *, "* ********************************"
!!print *, " in SUBROUTINE VerticalMeshCut(PolygonIn, X,Y,NX,MY,il,ir,PolygonOut) " 
!!print *, "* ********************************"

       CALL polygoncopy(17,PLGL,PolygonOut(j-il))
       CALL PLGDEL(PLGL)
       CALL polygoncopy(18,PLGR,PLGRTEMP)
       CALL PLGDEL(PLGR)
    enddo
    CALL polygoncopy(19,PLGRTEMP,PolygonOut(RX))
    CALL PLGDEL(PLGRTEMP)
  ENDIF

!   DO i=1,RX
!      !print *, "col", i
!      !print *, "number of nodes",polygonout(i)%nodesnum
!      do j=1,PolygonOut(i)%NODESNUM
!         !print *,  PolygonOut(i)%NODES(j)%VAL
!      ENDDO
!   enddo

!do i=1,RX
!   DEALLOCATE(PolygonOut(i)%NODES)
!   PolygonOut(i)%NODESNUM=0
!ENDDO
!DEALLOCATE (PolygonOut)


END SUBROUTINE VerticalMeshCut

!////////////////////////////////////////////////////////////////




!--------------------------------------------------------
!       triangle area , outer product
!---------------------------------------------------------
SUBROUTINE TRI_AREA(V1,V2,V3,AREA)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)                  :: V1,V2,V3
REAL(KIND=8),INTENT(OUT)                 :: AREA

 AREA = 0.5* ABS(V1%VAL(1)*(V2%VAL(2)-V3%VAL(2)) &
           & - V2%VAL(1)*(V1%VAL(2)-V3%VAL(2)) &
           & + V3%VAL(1)*(V1%VAL(2)-V2%VAL(2)))


END SUBROUTINE



!--------------------------------------------------------
subroutine init_cell(N,h,nox,noy,i,j,PIN)
implicit none

integer                                   :: N
real(kind=8)                              :: h
real(kind=8),dimension(-1:N)             :: nox,noy
INTEGER                                  :: i,j
TYPE(polygon)                            :: PIN

pin%nodesnum = 4
allocate(PIN%nodes(4))
PIN%CENTER%VAL(1) = nox(i)
PIN%CENTER%VAL(2) = noy(j)

!PIN%tag = 0

PIN%NODES(1)%VAL(1)= nox(i)-h/2.0d0
PIN%NODES(1)%VAL(2)= noy(j)-h/2.0d0
PIN%NODES(2)%VAL(1)= nox(i)+h/2.0d0
PIN%NODES(2)%VAL(2)= noy(j)-h/2.0d0
PIN%NODES(3)%VAL(1)= nox(i)+h/2.0d0
PIN%NODES(3)%VAL(2)= noy(j)+h/2.0d0
PIN%NODES(4)%VAL(1)= nox(i)-h/2.0d0
PIN%NODES(4)%VAL(2)= noy(j)+h/2.0d0

PIN%SIDES(1)%PT(1)= PIN%NODES(1)
PIN%SIDES(1)%PT(2)= PIN%NODES(2)
PIN%SIDES(2)%PT(1)= PIN%NODES(2)
PIN%SIDES(2)%PT(2)= PIN%NODES(3)
PIN%SIDES(3)%PT(1)= PIN%NODES(3)
PIN%SIDES(3)%PT(2)= PIN%NODES(4)
PIN%SIDES(4)%PT(1)= PIN%NODES(4)
PIN%SIDES(4)%PT(2)= PIN%NODES(1)

PIN%side_cen(1)%VAL(1)= 0.5D0*(PIN%SIDES(1)%PT(1)%VAL(1)+ PIN%SIDES(1)%PT(2)%VAL(1))
PIN%side_cen(1)%VAL(2)= 0.5D0*(PIN%SIDES(1)%PT(1)%VAL(2)+ PIN%SIDES(1)%PT(2)%VAL(2))
PIN%side_cen(2)%VAL(1)= 0.5D0*(PIN%SIDES(2)%PT(1)%VAL(1)+ PIN%SIDES(2)%PT(2)%VAL(1))
PIN%side_cen(2)%VAL(2)= 0.5D0*(PIN%SIDES(2)%PT(1)%VAL(2)+ PIN%SIDES(2)%PT(2)%VAL(2))
PIN%side_cen(3)%VAL(1)= 0.5D0*(PIN%SIDES(3)%PT(1)%VAL(1)+ PIN%SIDES(3)%PT(2)%VAL(1))
PIN%side_cen(3)%VAL(2)= 0.5D0*(PIN%SIDES(3)%PT(1)%VAL(2)+ PIN%SIDES(3)%PT(2)%VAL(2))
PIN%side_cen(4)%VAL(1)= 0.5D0*(PIN%SIDES(4)%PT(1)%VAL(1)+ PIN%SIDES(4)%PT(2)%VAL(1))
PIN%side_cen(4)%VAL(2)= 0.5D0*(PIN%SIDES(4)%PT(1)%VAL(2)+ PIN%SIDES(4)%PT(2)%VAL(2))


end subroutine init_cell



!------------------------------------------------------
! Linear reconstruction normal and initial guess
!        ax+by+d=0
!  input:  G ,X, Y
!------------------------------------------------------
SUBROUTINE N_COMP(N,h,iIN,jIN,G,X,Y,aout,bout,dout)
IMPLICIT NONE

integer,INTENT(IN)                                                          :: N
real(kind=8),intent(in)                                                     :: h
INTEGER     , INTENT(IN)                                                    :: iIN, jIN
REAL(KIND=8), INTENT(IN),DIMENSION(N,N)                                 :: G
REAL(KIND=8), INTENT(IN),DIMENSION(N)                                     :: X,Y
REAL(KIND=8), DIMENSION(iIN-1:iIN+1,jIN-1:jIN+1)                            :: w
REAL(KIND=8)                                                                :: w1,w2
REAL(KIND=8)                                                 :: aa,bb,cc,dd,ee,ff,gg,hh
REAL(KIND=8)                                                                :: aout, bout, dout
INTEGER                                                                     :: i,j
REAL(KIND=8), EXTERNAL                                                      :: NORM_2d



!calculate the weight
do i= iIN-1, iIN+1
  do j= jIN-1, jIN+1
    CALL SDD(h/2, g(i,j)/h, w1)               !-----------------------check  kin
    CALL SDD(h/2, NORM_2d(X(i),y(j),x(iIN),y(jIN))/h, w2)
!    !print *, "w1",w1,"w2",w2
    w(i,j) = w1*w2 
  enddo
enddo

!do i = iin-1,iin+1 
!  !print *, w(i,:) 
!enddo

 aa=0.0d0
 bb=0.0d0
 cc=0.0d0
 dd=0.0d0
 ee=0.0d0
 ff=0.0d0
 gg=0.0d0
 hh=0.0d0


do i= iIN-1, iIN+1
  do j= jIN-1, jIN+1
     aa = aa + 2.0d0*w(i,j)*((x(i)-x(iIN))**2.0d0)
     bb = bb + 2.0d0*w(i,j)*(x(i)-x(iIN))*(y(j)-y(jIN))
     cc = cc + 2.0d0*w(i,j)*(x(i)-x(iIN))
     dd = dd + 2.0d0*w(i,j)*((y(j)-y(jIN))**2.0d0)
     ee = ee + 2.0d0*w(i,j)*(y(j)-y(jIN))
     ff = ff + 2.0d0*w(i,j)*(x(i)-x(iIN))*G(i,j)
     gg = gg + 2.0d0*w(i,j)*(y(j)-y(jIN))*G(i,j)
     hh = hh + 2.0d0*w(i,j)*G(i,j)
  enddo
enddo

!!print *, aa,bb,cc,dd,ee,ff,gg,hh

if(bb .eq. 0.0d0  .and. ee .eq. 0.0d0) THEN
   bout = gg/dd
   aout = (ff-cc*hh)/(aa-cc*cc)
   dout = hh - cc*aout
else
  dout = ((aa*hh-ff*cc)/(aa*ee-bb*cc)-(aa*gg-bb*ff)/(aa*dd-bb*bb))/ &
         &((aa-cc*cc)/(aa*ee-bb*cc)-(aa*ee-bb*cc)/(aa*dd-bb*bb))
  bout = (aa*gg-bb*ff)/(aa*dd-bb*bb)-((aa*ee-bb*cc)/(aa*dd-bb*bb))*dout
  aout = (ff-cc*dout-bb*bout)/aa
endif
!dout = ((aa*hh-ff*cc)*(aa*dd-bb*bb)-(aa*gg-bb*ff)/(aa*ee-bb*cc))/ &
!    &((aa-cc*cc)*(aa*dd-bb*bb)-(aa*ee-bb*cc)*(aa*ee-bb*cc))
!bout = (aa*gg-bb*ff)/(aa*dd-bb*bb)-((aa*ee-bb*cc)/(aa*dd-bb*bb))*dout
!aout = (ff-cc*dout-bb*bout)/aa



END SUBROUTINE N_COMP
!---------------------------------------------------------
Subroutine SDD(a,x,delta)
IMPLICIT NONE

REAL(KIND=8), INTENT(IN) :: a
REAL(KIND=8), INTENT(IN) :: X
REAL(KIND=8), Intent(out):: delta

 delta = (1.0d0/(a*sqrt(pi)))*exp(-x**2.0d0/(a**2.0d0)) 

end subroutine SDD



!-----------------------------------
subroutine deter_min(a1,a2,a3,a4,b1,b2,b3,b4,pout)
implicit none

real(kind=8)         :: a1,a2,a3,a4,b1,b2,b3,b4
real(kind=8)         :: p1,p2
real(kind=8)         :: pout

call pp_com(a1,a2,p1)
call pp_com(p1,a3,p2)
call pp_com(p2,a4,p1)
call pp_com(p1,b1,p2)
call pp_com(p2,b2,p1)
call pp_com(p1,b3,p2)
call pp_com(p2,b4,p1)

pout = p1

end subroutine deter_min
!---------------------------------
subroutine pp_com(v1,v2,vout)
implicit none

real(kind=8)       :: v1,v2
real(kind=8)       :: vout

IF(abs(v1) .lt. abs(v2))THEN
  vout = v1
ELSE
  vout = v2
endif


end subroutine pp_com
!--------------------------------


!---------------------------------------------------------
SUBROUTINE RECON_L(a,b,c,X,VOUT)
IMPLICIT NONE

REAL(KIND=8),INTENT(IN)   :: a,b,c
TYPE(POINTS)              :: x
REAL(KIND=8)              :: VOUT

  VOUT= a*x%val(1)+b*x%val(2)+c

END SUBROUTINE RECON_L



!----------------------------------------------------------------
subroutine cell_line_inter(cellin,a,b,c,pout)
implicit none

type(polygon),intent(in)    :: cellin
real(kind=8),intent(in)  :: a,b,c
integer                  :: i,flag1,flag2
type(points)             :: p
type(points)             :: pout(2)

!do i = 1,cellin%nodesnum
!   !!print *, "cellin",i, cellin%nodes(i)%val 
!enddo
            
flag2 = 0

do i=1,4
!   !!print *, i,"cellside",cellin%sides(i)%PT
   call cell_side_inter(a,b,c,i,cellin%sides(i),flag1,p)
   if (flag1 .eq. 1 .and. flag2 .eq. 0) then
     pout(1)%val(1) = p%val(1)
     pout(1)%val(2) = p%val(2)
     flag2 = 1
   elseif(flag1 .eq. 1 .and. flag2 .eq. 1) then
     pout(2)%val(1) = p%val(1)
     pout(2)%val(2) = p%val(2) 
     exit          !>>>>>>>>>>>>>>>>>>>>>>>>
   endif
enddo

!   !!print *, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

end subroutine CELL_LINE_INTER
!-----------------------------------------
!------------------------------------------------
!   line cell-side intersection
!--------------------------------------------------
subroutine cell_side_inter(a,b,c,i,sidein,flag,pout)
implicit none

real(kind=8)             :: a,b,c
integer                  :: i
type(SIDETYPE)           :: sidein 
type(points)             :: pout
real(kind=8)             :: x1,y1,x2,y2,x0,y0
integer,intent(out)      :: flag

flag = 0
x1 = sidein%pt(1)%val(1)
x2 = sidein%pt(2)%val(1)
y1 = sidein%pt(1)%val(2)
y2 = sidein%pt(2)%val(2)

!write(*,*) "a=",a,"b=",b,"i=",i

IF(a .ne. 0.0d0 .and. b.ne. 0.0d0)then
  IF((i .eq. 1) .or. (i .eq. 3)) THEN
     x0 = (-b*y1-c)/a

!      !!print *,  "a",a, "b",b,"c",c
!      !!print *, "y1",y1,"y2",y2
!      !!print *, "x0",x0,"x1",x1,"x2",x2
     
  
!     if((x0 .ge. x1 .and. x0 .le. x2) .or. &
!       & (x0 .le. x1 .and. x0 .ge. x2)) then
     if( x0 .ge. min(x1,x2)  .and. x0 .le. max(x1,x2) ) then
        flag = 1
        pout%val(1) = x0
        pout%val(2) = y1
     else
        pout%val(1) = 0.0d0
        pout%val(2) = 0.0d0
     endif
  elseif((i .eq. 2) .or. (i .eq. 4)) THEN
     y0 = (-a*x1-c)/b
 !    !!print *,  "a",a, "b",b,"c",c
 !    !!print *, "x1",x1,"x2",x2
 !    !!print *, "y0",y0,"y1",y1,"y2",y2  

!     if((y0.ge. y1 .and. y0 .le. y2) .or. &
!         & (y0 .le. y1 .and. y0 .ge. y2)) then
     if( y0 .ge. min(y1,y2)  .and. y0 .le. max(y1,y2) ) then
        flag = 1
        pout%val(1) = x1
        pout%val(2) = y0
     else
        pout%val(1) = 0.0d0
        pout%val(2) = 0.0d0
     endif
  else
    !print *, "cell_side_inter: wrong numbers of vertex"    
  endif
elseif(a .eq. 0.0d0 .and. b .ne. 0.0d0)then
  if(i.eq. 2 .or. i .eq. 4) then
    flag = 1
    pout%val(1)= x1
    pout%val(2)= -c/b
!  else
!    write(*,*) "error, parallel"
  endif
elseif(a .ne. 0.0d0 .and. b .eq. 0.0d0)then
   if(i.eq. 1 .or. i .eq. 3) then
    flag = 1
    pout%val(1)= -c/a
    pout%val(2)= y1
!   else
!    write(*,*) "error, parallel"
   endif
endif

end subroutine cell_side_inter
!---------------------------------
subroutine point_inter(a,b,c,rg,pin,dist_out)
implicit none

real(kind=8)              :: a,b,c
type(points)              :: rg(2), pin
type(points)              :: inter
real(kind=8)              :: dist1,dist2,dist_out


IF(a .ne. 0.0d0) then
   inter%val(1) = ((b**2.0d0/a)*pin%val(1) - b*pin%val(2) -c)/(a+(b**2.0d0/a))
   inter%val(2) = (-a*inter%val(1)-c)/b
else
   inter%val(1) = pin%val(1)
   inter%val(2) = -c/b   
endif

iF( (inter%val(1) .le. max(rg(1)%val(1),rg(2)%val(1))) .and. &
   & ( inter%val(1) .ge. min(rg(1)%val(1),rg(2)%val(1)))&
   & .and. ( inter%val(2) .le. max(rg(1)%val(2),rg(2)%val(2))) .and. &
   &  (inter%val(2) .ge. min(rg(1)%val(2),rg(2)%val(2)))) then
  call pp_dist(pin,inter,dist_out)
else
  call pp_dist(pin,rg(1),dist1)
  call pp_dist(pin,rg(2),dist2)
  dist_out = min(dist1,dist2)
endif

end subroutine POINT_INTER



!------------------------------------------------------------
SUBROUTINE PP_DIST(P1,P2,DIST)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)         :: P1,P2
REAL(KIND=8),INTENT(OUT)        :: DIST

DIST = sqrt((P1%VAL(1)-P2%VAL(1))**2.0D0+ (P1%VAL(2)-P2%VAL(2))**2.0D0)


END SUBROUTINE PP_DIST



!---------------------------------------------------------------------
subroutine zalesakdist(cas,dist,xx,yy)
IMPLICIT NONE

integer,intent(in) :: cas
REAL(KIND=8) :: dist,xx,yy,x,y
REAL(KIND=8) :: dist1,dist2
REAL(KIND=8) :: xtarget(2)

!      if (probtype.ne.28) then
!       print *,"probtype invalid"
!       stop
!      endif



      x=xx
      y=yy
      xtarget(1)=xx
      xtarget(2)=yy

      if (cas .eq. 0) then
       dist=sqrt((x-50.0)**2+(y-75.0)**2)-15.0
       if ((x.ge.47.5).and.(x.le.52.5)) then
        if (y.le.60.0) then
         if (x.lt.50.0) then
          dist=sqrt( (y-60.0)**2+(x-47.5)**2 )
         else
          dist=sqrt( (y-60.0)**2+(x-52.5)**2 )
         endif
        else if (y.le.85.0) then
         if (x.lt.50.0) then
          dist1=x-47.5
         else
          dist1=52.5-x
         endif
         dist2=85.0-y
         dist=min(dist1,dist2)
        else if ((y.le.90.0).and.(dist.le.0.0d0)) then
         dist=max(dist,85.0-y)
        endif
       else if ((dist.lt.0.0d0).and.(x.lt.47.5)) then
        if (y.le.85.0) then
         dist=max(dist,(x-47.5))
        else
         dist=max(dist,-sqrt( (x-47.5)**2+(y-85.0)**2 ) )
        endif
       else if ((dist.lt.0.0d0).and.(x.gt.52.5)) then
        if (y.le.85.0) then
         dist=max(dist,(52.5-x))
        else
         dist=max(dist,-sqrt( (x-52.5)**2+(y-85.0)**2 ) )

        endif
       endif

      else if (cas.eq.1) then
       dist=sqrt((x-50.0)**2+(y-75.0)**2)-15.0
      else
       print *,"axis_dir invalid zalesakdist"
       stop
      endif


      return
      end subroutine


!////////////////////////////////////////////////////////////////


!-----------------------------------------------------------------
SUBROUTINE VOL_FRAC_CAL(POLYGONIN,a,b,c,POS_PLG,NEG_PLG)
IMPLICIT NONE

TYPE(POLYGON),INTENT(IN)    :: POLYGONIN
REAL(KIND=8)                :: a,b,c
TYPE(POLYGON)               :: PLGL,PLGR
REAL(KIND=8)                :: VAL,VAL1
!REAL(KIND=8)                :: TT_AREA, POS_AREA
INTEGER                     :: LR_FLAG     ! 0 for L    1 for R
TYPE(POLYGON),INTENT(OUT)   :: POS_PLG, NEG_PLG             


CALL l_cutpolygon(POLYGONIN,a,b,c,PLGL,PLGR)

CALL CUT_DIFF(a,b,c,POLYGONIN%NODES(1),val)
  IF(VAL .GE. 0.0D0) THEN
    LR_FLAG = 1   !POINT(1) AT PLGR
  ELSE
    LR_FLAG = 0   !POINT(1) AT PLGL
  ENDIF
                                                              !!!!!!!!!!!!!!!!!!!!!!
CALL RECON_L(a,b,c,POLYGONIN%NODES(1),val1)                  !!!!!!!!pos to neg!!!!!!!
IF(LR_FLAG .eq. 1 .AND. val1 .GT. 0) THEN                          !!!!!!!!!!!!!!!!!!!!!!
   CALL polygoncopy(20,PLGL,POS_PLG)
   CALL polygoncopy(21,PLGR,NEG_PLG)
ELSEIF(LR_FLAG .eq. 1 .AND. val1 .LT. 0)THEN
   CALL polygoncopy(22,PLGR,POS_PLG)
   CALL polygoncopy(23,PLGL,NEG_PLG)
ELSEIF(LR_FLAG .eq. 0 .AND. val1 .LT. 0)THEN
   CALL polygoncopy(24,PLGL,POS_PLG)
   CALL polygoncopy(25,PLGR,NEG_PLG)
ELSEIF(LR_FLAG .eq. 0 .AND. val1 .GT. 0)THEN
   CALL polygoncopy(26,PLGR,POS_PLG)
   CALL polygoncopy(27,PLGL,NEG_PLG)
ELSE
   !print *, "THE RECONSTRUCTION LINE GOES THROUGHT ONE VERTEX"
ENDIF


CALL PLGDEL(PLGL)
CALL PLGDEL(PLGR)
END SUBROUTINE VOL_FRAC_CAL

!-------------------------------------------------------
!-----------------------------------------------------




SUBROUTINE INIT_F0(NMAT,N,h,nox,noy,phi,cell,vf,centroid)
IMPLICIT NONE

INTEGER,INTENT(IN)      :: NMAT
integer,intent(in)      :: N
real(kind=8)            :: h
real(kind=8),intent(in) :: nox(N),noy(N)
type(polygon)           :: cell(N,N)
real(kind=8),intent(in) :: phi(N,N)
INTEGER                 :: i,j
real(kind=8)            :: aout,bout,dout,C
TYPE(POLYGON)           :: POS_PLG,NEG_PLG
REAL(KIND=8)            :: POS_VOL
INTEGER                 :: FLAG

real(kind=8)            :: vf(N,N,NMAT)
type(ifseg)             :: seg(N,N)
TYPE(POINTS),DIMENSION(N,N,NMAT):: CENTROID

print *,"i=0,N-1 and j=0,N-1"
stop

do i=1,N
   vf(i,1,1) = 0.0d0
   vf(i,N,1) = 0.0d0
   vf(i,1,2) = 1.0d0
   vf(i,N,2) = 1.0d0
   centroid(i,1,1)%val(1) = 0.0d0
   centroid(i,1,1)%val(2) = 0.0d0  
   centroid(i,N,1)%val(1) = 0.0d0
   centroid(i,N,1)%val(2) = 0.0d0  
   centroid(i,1,2)%VAL= cell(i,1)%center%VAL 
   centroid(i,N,2)%VAL= cell(i,N)%center%VAL 
enddo

print *,"i=0,N-1 and j=0,N-1"
stop

do i=1,N
   vf(1,i,1) = 0.0d0
   vf(N,i,1) = 0.0d0
   vf(1,i,2) = 1.0d0
   vf(N,i,2) = 1.0d0
   centroid(1,i,1)%val(1) = 0.0d0
   centroid(1,i,1)%val(2) = 0.0d0  
   centroid(N,i,1)%val(1) = 0.0d0
   centroid(N,i,1)%val(2) = 0.0d0
   centroid(1,i,2)%VAL= cell(1,i)%center%VAL 
   centroid(N,i,2)%VAL= cell(N,i)%center%VAL   
enddo

POS_PLG%NODESNUM = 0
NEG_PLG%NODESNUM = 0

print *,"i=1,N-2 and j=1,N-2"
stop

do i=2,N-1
   do j=2,N-1
    print *, i,j
    CALL NB_CHECK(i,j,N,flag,phi)
    print *, "flag=",flag
    print *, ">>>>>>>>>>>1<<<<<<<<<<<<<<"
    IF (flag .eq. 1) then
       print *, ">>>>>>>>>>>>2<<<<<<<<<<<<<"
       CALL N_comp(N,h,i,j,phi,nox,noy,aout,bout,dout)          ! ???????????????
        c = -aout*nox(i)-bout*noy(j)+dout
       
       call cell_line_inter(cell(i,j),aout,bout,c,seg(i,j)%PT)
       print *, seg(i,j)%pt(1),seg(i,j)%pt(2)


       CALL VOL_FRAC_CAL(cell(i,j),aout,bout,c,POS_PLG,NEG_PLG)
       print *, "in cell", i,j
       call outputplg(pos_plg)
       call outputplg(neg_plg)
       
       CALL POLY_CENTROID(POS_PLG,centroid(i,j,1))
       CALL POLY_CENTROID(NEG_PLG,centroid(i,j,2))
       
!       !print *, "initF=============================="
       CALL poly_area(POS_PLG,POS_VOL)
      
       VF(i,j,1)= POS_VOL/(h**2.0d0)
       vf(i,j,2) = 1-vf(i,j,1)       

       CALL PLGDEL(POS_PLG)
       CALL PLGDEL(NEG_PLG)

    ELSEIF(FLAG .EQ. -1  .AND. PHI(I,J) .LT. 0.0D0) THEN
       VF(I,J,1) = 1.0D0
       vf(i,j,2) = 0.0d0
       centroid(i,j,1)%val = cell(i,j)%center%val 
       centroid(i,j,2)%val(1) = 0.0d0 
       centroid(i,j,2)%val(2) = 0.0d0
       
    ELSEIF(FLAG .EQ. -1 .AND. PHI(I,J) .Ge. 0.0D0) THEN
       VF(I,J,1) = 0.0D0
       vf(i,j,2) = 1.0d0

       centroid(i,j,2)%val = cell(i,j)%center%val 
       centroid(i,j,1)%val(1) = 0.0d0 
       centroid(i,j,1)%val(2) = 0.0d0

    ENDIF
  enddo
enddo 

!do i = 1,N
!   do j =1,N
!      if(seg(i,j)%flag .eq. 1) then
!         print *, i,j, seg(i,j)%pt(1),seg(i,j)%pt(2) 
!      endif
!   enddo
!enddo

!!print *, "============================================="
! close(11)

END SUBROUTINE INIT_F0

!--------------------------------------------
SUBROUTINE nb_check(xi,yj,N,flag,phi) 
IMPLICIT NONE

integer               :: N
INTEGER,INTENT(IN)    :: xi,yj
real(kind=8)          :: phi(N,N)
integer               :: i,j
integer               :: flag


FLAG = -1

do i=-1,1,2
   do j= -1,1,2
      IF (phi(xi,yj)* phi(xi+i,yj+j) .Lt. 0.0D0  .AND. &
         & phi(xi,yj)* (phi(xi+i,yj+j)+phi(xi,yj)) .Lt. 0.0D0 ) THEN
        FLAG =1
      ENDIF
   enddo
enddo


END SUBROUTINE

!-----------------------------------------------------------------------

SUBROUTINE CUT_DIFF(a,b,d,x,val)
IMPLICIT NONE

TYPE(POINTS)            :: x
REAL(KIND=8),intent(out):: val
REAL(KIND=8),intent(in) :: a,b,d

  IF (a .gt. 0.0d0) THEN
    val= a*x%val(1)+ b*x%val(2)+ d
  ELSEIF(a .lt. 0.0d0) THEN
    val= -a*x%val(1)- b*x%val(2)- d
  ELSE
    IF(b .gt. 0.0d0) THEN
      VAL = b*x%val(2)+ d
    ELSEIF(b .lt. 0.0d0) THEN
      VAL = -b*x%val(2)- d
    ELSE
      !print *, "a and b are both 0"
    ENDIF
  ENDIF

END SUBROUTINE CUT_DIFF
!----------------------------------------------------------------------
!------------------------------------------------------------------------
!-----------------------NEW---------------------------------------------
subroutine parker_young(iin,jin,n,X,Y,f,h,mx,my)
implicit none

integer, intent(in)           :: iin,jin
integer, intent(in)           :: n
real(kind=8)                  :: x(n),y(n)
real(kind=8)                  :: f(n,n)
real(kind=8),intent(in)       :: h
real(kind=8)                  :: m_x(4),m_y(4)
real(kind=8),intent(out)      :: mx,my


m_x(1) = 0.5d0*(f(iin+1,jin+1)+f(iin+1,jin)-f(iin,jin+1)-f(iin,jin))/h
m_y(1) = 0.5d0*(f(iin+1,jin+1)+f(iin,jin+1)-f(iin+1,jin)-f(iin,jin))/h

m_x(2) = 0.5d0*(f(iin,jin+1)+f(iin,jin)-f(iin-1,jin+1)-f(iin-1,jin))/h
m_y(2) = 0.5d0*(f(iin,jin+1)+f(iin-1,jin+1)-f(iin,jin)-f(iin-1,jin))/h

m_x(3) = 0.5d0*(f(iin,jin)+f(iin,jin-1)-f(iin-1,jin)-f(iin-1,jin-1))/h
m_y(3) = 0.5d0*(f(iin,jin)+f(iin-1,jin)-f(iin,jin-1)-f(iin-1,jin-1))/h

m_x(4) = 0.5d0*(f(iin+1,jin)+f(iin+1,jin-1)-f(iin,jin)-f(iin,jin-1))/h
m_y(4) = 0.5d0*(f(iin+1,jin)+f(iin,jin)-f(iin+1,jin-1)-f(iin,jin-1))/h

mx = 0.25d0*(m_x(1)+m_x(2)+m_x(3)+m_x(4))
my = 0.25d0*(m_y(1)+m_y(2)+m_y(3)+m_y(4))



end subroutine parker_young
!-------------------------------------------------------------

!---------------------------------------------------------------
subroutine poly_area(polygonin,area_out)
implicit none

type(polygon)                :: polygonin
real(kind=8)                 :: area,area_out
integer                      :: i

area_out = 0.0d0

if (polygonin%nodesnum .lt. 3 .and. polygonin%nodesnum .gt. 0) then
      !print *, "poly_area: wrong numbers of vertex"
elseif(polygonin%nodesnum .gt. 2) then
      do i = 2,polygonin%nodesnum-1
         call tri_area(polygonin%nodes(1),polygonin%nodes(i),polygonin%nodes(i+1),area)
         area_out = area_out + area
      enddo
else
   area_out = 0.0d0
endif

end subroutine poly_area

!-----------------------------------------------------------
Subroutine f_diff(v_f,mx,my,mes,alpha,CellIn,fout)
implicit none

real(kind=8), intent(in)  :: v_f , mx , my, alpha, mes
type(polygon)             :: CellIn
real(kind=8)              :: area
TYPE(POLYGON)             :: plgl,plgr
real(kind=8)              :: fout

!print *,  "f_diff==========begin=================="
!print *, "mx=", mx, "my=",my,"alpha=",alpha
!print *,  "f_diff============end================"
!print *,  "f_diff$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"

if(mx .eq. 0.0d0 .and. my .eq. 0.0d0)then
   write(*,*) "mx=", mx, "my=",my,"alpha=",alpha
endif

call cutpolygon(mx,my,alpha,cellIn,plgl,plgr)

!print *, "plgl:"
call outputplg(plgl)
!print *, "plgr:"
call outputplg(plgr)

if( mx.lt. 0.0d0 .and. my .ne. 0.0d0)then
     call poly_area(plgl,area)
elseif( mx .gt. 0.0d0 .and. my .ne. 0.0d0) then  
  call poly_area(plgr,area) 
elseif(mx .eq. 0.0d0 .and. my .gt. 0.0d0) then 
     call poly_area(plgr,area)  
elseif(mx .eq. 0.0d0 .and. my .lt. 0.0d0) then                                
     call poly_area(plgl,area)        
elseif(my .eq. 0.0d0 .and. mx .gt. 0.0d0) then
      call poly_area(plgr,area) 
elseif(my .eq. 0.0d0 .and. mx .lt. 0.0d0) then
     call poly_area(plgl,area)
else
    !print *, "mx,my,both 0"
    area = 0.0d0
endif

 !print *, "area=",area,"mes=",mes,"v_f=",v_f

 fout = v_f - area/mes

!print *, "f_out=",fout
!print *,  "f_diff$$$$$$$$$$$$$$end$$$$$$$$$$$$$$$$$"
!print *, "v_f=", v_f, "area=", area, "mes=", mes
!call outputplg(plgl)
!!print *, "plgr:"
!call outputplg(plgr)


!!print *, "***********************************"

call plgdel(plgl)
call plgdel(plgr)

end subroutine f_diff


!-----------------------------------------------------------------
subroutine Brents_Method(a,b,h,v_f,mx,my,mes,CellIn,root)
implicit none

real(kind=8),intent(in)        :: v_f,h,mx,my
real(kind=8)                   :: a,b
type(polygon)                  :: cellin
real(kind=8)                   :: mes
real(kind=8)                   :: a1,b1
real(kind=8)                   :: fa,fb,s,c,fc,fs,d,fa1,fb1

real(kind=8),intent(out)       :: root
integer                        :: mflag

!print *, "Brents Method1========begin=============="
!print *, "========================================="           

!Write(10,*) "a=" , a , "b=" ,b

!print *, "volume fraction = ", v_f 
!print *, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"  
!print *, "calc_fa>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
if(mx .eq. 0.0d0 .and. my .eq. 0.0d0)then
  write(*,*) "=====In Brents ==============="
endif
call f_diff(v_f,mx,my,mes,a,CellIn,fa)
!print *, "calc_fb>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
call f_diff(v_f,mx,my,mes,b,CellIn,fb)

!print *, "Brents Method1======end================"
!print *, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"  



!print *, "fa=", fa, "fb=", fb


if (fa*fb .ge. 0.0d0 ) then 
   print *, "out of cell"
endif

if (abs(fa) .lt. abs(fb)) then
  a1 = b
  b1 = a 
  a = a1
  b = b1
  fa1 = fb
  fb1 = fa
  fa = fa1
  fb = fb1
endif

 c = a 
!print *, "=====fc========begin=============="
!print *, "=========================================" 
 call f_diff(v_f,mx,my,mes,c,CellIn,fc)
!print *, "c=",c,"fc=",fc
!print *, "B2222222222222222222222========end=============="
!print *, "================end  =========================" 
 mflag = 1

do while( (abs(fb).gt. 0.0d0) .and. &
          (abs(fs).gt. 0.0d0) .and. &
          (abs(b-a).gt.afrac_eps))
   if ((fa .ne. fc) .and. (fb .ne. fc)) then
       s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) &
           & + c*fa*fb/((fc-fa)*(fc-fb))
   else
       s = b - fb*(b-a)/(fb-fa)
   endif

 !print *, "s=",s
   
   if ((( s  - 0.25d0*(3.0d0*a+b))*(s-b) .gt. 0.0d0 ) .or. &
       & (mflag .eq. 1 .and. abs(s-b) .ge. 0.5d0*abs(b-c)) .or. &
       & (mflag .eq. 0 .and. abs(s-b) .ge. 0.5d0*abs(c-d)) .or. &
       & (mflag .eq. 1 .and. abs(b-c) .lt. afrac_eps) .or. &
       & (mflag .eq. 0 .and. abs(c-d) .lt. afrac_eps) ) then
      s = 0.5d0*(a+b)
      mflag = 1
   else
      mflag = 0
   endif

!print *, "33333333333333========begin=============="
!print *, "=========================================" 
!print *, "@@@@@@@@@@@@@@@@@@fsfsfs@@@@@@@@@@@@@@@@@@@@@@@@@@@"
   call f_diff(v_f,mx,my,mes,s,CellIn,fs)
   
!print *, "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
!print *, "33333333333333========end=============="
!print *, "=========================================" 

   d = c
   c = b

   if (fa*fs .lt. 0.0d0 ) then
      b = s
   else
      a = s
   endif
   
   if(abs(fa) .lt. abs(fb)) then
      a1 = b
      b1 = a
      a = a1
      b = b1
     fa1 = fb
     fb1 = fa
     fa = fa1
     fb = fb1
   endif  

enddo

if ((abs(fb).lt. afrac_eps).or. &
    (abs(b-a).lt.afrac_eps)) then
   root = b
elseif(abs(fs) .lt. afrac_eps) then
   root = s
else
  !print *, "Neither b nor s is the root" 
endif


end subroutine Brents_Method

!---------------------------------------------------------------------
subroutine vol_matching(mx,my,alpha,h,N,iin,jin,cell,mes,vf,alpha_out,plgout)
implicit none

real(kind=8),intent(in)              :: vf(N,N)
real(kind=8),intent(in)              :: mx,my,h
real(kind=8),dimension(4)            :: alpha
real(kind=8)                         :: alpha_min,alpha_max
integer     ,intent(in)              :: N, iin,jin
type(polygon),intent(in)             :: cell(N,N)
REAL(KIND=8),INTENT(IN)              :: mes(N,N)
type(polygon)                        :: plgl,plgr


real(kind=8)                         :: alpha_out
type(polygon)                        :: plgout


 alpha_min = minval(alpha)
 alpha_max = maxval(alpha)

!print *, "min=" , alpha_min , "max=" ,alpha_max

if(mx .eq. 0.0d0 .and. my .eq. 0.0d0)then
  write(*,*) "=====In volm ==============="
endif

call Brents_Method(alpha_min,alpha_max,h,vf(iin,jin),mx,my,  &
       & mes(iin,jin),Cell(iin,jin),alpha_out) 

!print *, "new_alpha", alpha_out,"*******************************"




call cutpolygon(mx,my,alpha_out,cell(iin,jin),plgl,plgr)

!print *, ".>>>>>>>>>>>>>>>>>end<<<<<<<<<<<<<<<<<<<"

if( mx.lt. 0.0d0 .and. my .ne. 0.0d0)then
     call polygoncopy(28,plgl,plgout)
elseif( mx .gt. 0.0d0 .and. my .ne. 0.0d0) then  
     call polygoncopy(29,plgr,plgout)
elseif(mx .eq. 0.0d0 .and. my .gt. 0.0d0) then 
     call polygoncopy(30,plgr,plgout) 
elseif(mx .eq. 0.0d0 .and. my .lt. 0.0d0) then                                
    call polygoncopy(31,plgl,plgout)        
elseif(my .eq. 0.0d0 .and. mx .gt. 0.0d0) then
    call polygoncopy(32,plgr,plgout) 
elseif(my .eq. 0.0d0 .and. mx .lt. 0.0d0) then
     call polygoncopy(33,plgl,plgout)
else
    !print *, "mx,my,both 0"
    plgout%nodesnum = 0
endif



!if( mx.lt. 0.0d0)then 
!  call polygoncopy(10,plgl,plgout)
!elseif( mx .gt. 0.0d0) then  
!  call polygoncopy(10,plgr,plgout) 
!else 
!  if( my .gt. 0.0d0) then
!     call polygoncopy(10,plgr,plgout)
!  else
!     call polygoncopy(10,plgl,plgout)
!  endif
!endif

call plgdel(plgl)
call plgdel(plgr)

end subroutine vol_matching

SUBROUTINE INIT_V(N,xCC,yCC,local_probtype,iten, &
                scomp,local_sdim,vel)
IMPLICIT NONE

integer,intent(in)             :: N,scomp,local_sdim,local_probtype,iten
real(kind=8),intent(in)        :: xCC(-1:N),yCC(-1:N)
real(kind=8),intent(out)       :: vel(-1:N,-1:N,scomp+local_sdim-1)
INTEGER                  :: i,j
INTEGER                  :: scomp_offset
INTEGER                  :: nten
INTEGER                  :: nmat
real(kind=8), external  :: Uprescribe,Vprescribe 
real(kind=8) :: local_time

if (local_sdim.eq.2) then

 nmat=num_materials
 nten=((nmat-1)*(nmat-1)+nmat-1)/2

 scomp_offset=scomp-nmat-1
 if ((scomp_offset.ge.0).and.(scomp_offset.le.(nten-1)*local_sdim)) then
  if ((scomp_offset/local_sdim)*local_sdim.eq.scomp_offset) then
   local_time=0.0d0
   do i=-1,N
   do j=-1,N
     ! Uprescribe, Vprescribe in: multimat_FVM.F90
    vel(i,j,scomp)= Uprescribe(local_probtype,iten,local_time,xCC(i),yCC(j))
    vel(i,j,scomp+1)= Vprescribe(local_probtype,iten,local_time,xCC(i),yCC(j))
   enddo
   enddo
  else
   print *,"scomp_offset invalid"
   stop
  endif
 else
  print *,"scomp_offset invalid"
  stop
 endif

else
 print *,"local_sdim invalid"
 stop
endif

END SUBROUTINE INIT_V

SUBROUTINE INIT_LS(N,xCC,yCC,local_dx,local_probtype,local_nmat, &
                scomp,local_sdim,LS,ncomp_LS)
IMPLICIT NONE

integer,intent(in)             :: N,scomp,local_sdim,local_probtype
integer,intent(in)             :: local_nmat,ncomp_LS
real(kind=8),intent(in)        :: local_dx(local_sdim)
real(kind=8),intent(in)        :: xCC(-1:N),yCC(-1:N)
real(kind=8),intent(out)       :: LS(-1:N,-1:N,ncomp_LS)
INTEGER                  :: i,j,ils
real(kind=8) :: local_time
real(kind=8) :: local_LS(local_nmat*(local_sdim+1))

local_time=0.0d0

do i=-1,N
   do j=-1,N 
     ! in: multimat_FVM.F90
    call LSprescribe(local_probtype,local_nmat,local_sdim, &
            local_time,xCC(i),yCC(j),local_dx,local_LS)
    do ils=1,local_nmat*(local_sdim+1)
     LS(i,j,scomp+ils-1)=local_LS(ils)
    enddo
   enddo
enddo

END SUBROUTINE INIT_LS

!---------------------------------------------------------

!--------------------------------------------------------
subroutine poly_centroid(polygonin,c)
implicit none

type(polygon)              :: polygonin
integer                    :: i
real(kind=8)               :: A
type(points),intent(out)   :: c
integer                    :: N

 N= polygonin%nodesnum
 if (N .le. 2) then
    print *, "invalid polygon  in poly_centroid",N
    stop
 else
!print *, " N=", N
A = 0.0d0
do  i = 1,N-1
   A = A + 0.5d0*(polygonin%nodes(i)%val(1)*polygonin%nodes(i+1)%val(2)  &
           - polygonin%nodes(i+1)%val(1)*polygonin%nodes(i)%val(2))
enddo
   A = A + 0.5d0*(polygonin%nodes(N)%val(1)*polygonin%nodes(1)%val(2)  &
           - polygonin%nodes(1)%val(1)*polygonin%nodes(N)%val(2))

! print *, "A=", A 

 c%val(1) = 0
 c%val(2) = 0

do  i = 1,N-1
   c%val(1) = c%val(1) + (1.0d0/(6.0d0*A))*(polygonin%nodes(i)%val(1)+polygonin%nodes(i+1)%val(1))*  &
          & (polygonin%nodes(i)%val(1)*polygonin%nodes(i+1)%val(2) - &
          & polygonin%nodes(i+1)%val(1)*polygonin%nodes(i)%val(2))
   c%val(2) = c%val(2) + (1.0d0/(6.0d0*A))*(polygonin%nodes(i)%val(2)+polygonin%nodes(i+1)%val(2))* &
          & (polygonin%nodes(i)%val(1)*polygonin%nodes(i+1)%val(2) - &
          & polygonin%nodes(i+1)%val(1)*polygonin%nodes(i)%val(2))
enddo

   c%val(1) = c%val(1) + (1.0d0/(6.0d0*A))*(polygonin%nodes(N)%val(1)+polygonin%nodes(1)%val(1))*  &
          & (polygonin%nodes(N)%val(1)*polygonin%nodes(1)%val(2) - &
          & polygonin%nodes(1)%val(1)*polygonin%nodes(N)%val(2))
   c%val(2) = c%val(2) + (1.0d0/(6.0d0*A))*(polygonin%nodes(N)%val(2)+polygonin%nodes(1)%val(2))*   &
          & (polygonin%nodes(N)%val(1)*polygonin%nodes(1)%val(2) - &
          &  polygonin%nodes(1)%val(1)*polygonin%nodes(N)%val(2))

endif

end subroutine
!////////////////////////////////////////////////////////////
subroutine output_multi_polygon(mx,my,intercept,h,N,iin,jin,cell,plgin,mes,&
                                 plgout_d,plgout_l,alpha)
implicit none

!TYPE(POINTS),intent(in)              :: centroid
!real(kind=8),intent(in)              :: vf(N,N)
real(kind=8),intent(in)              :: mx,my,h
real(kind=8)                         :: nx, ny
real(kind=8)                         :: alpha
integer     ,intent(in)              :: N, iin,jin
type(polygon),intent(in)             :: cell(N,N)
type(polygon)                        :: plgin
REAL(KIND=8),INTENT(IN)              :: mes(N,N)
type(polygon)                        :: plgl,plgr
real(kind=8)                         :: intercept

type(polygon)                        :: plgout_d , plgout_l


print *,"i=0,N-1 and j=0,N-1 output_multi_polygon"
stop

nx = mx
ny = my


if(abs(mx) .lt. afrac_eps) then
  nx = 0.0d0
endif

if(abs(my) .lt. afrac_eps) then
  ny = 0.0d0
endif

alpha = nx*cell(iin,jin)%center%val(1) + ny*cell(iin,jin)%center%val(2) - intercept

call cutpolygon(nx,ny,alpha,plgin,plgl,plgr)

if( nx.lt. 0.0d0 .and. ny .ne. 0.0d0)then
     call polygoncopy(41,plgl,plgout_d)
     call polygoncopy(42,plgr,plgout_l)
elseif( nx .gt. 0.0d0 .and. ny .ne. 0.0d0) then  
     call polygoncopy(43,plgr,plgout_d)
     call polygoncopy(44,plgl,plgout_l)
elseif(nx .eq. 0.0d0 .and. ny .gt. 0.0d0) then 
     call polygoncopy(45,plgr,plgout_d)
     call polygoncopy(46,plgl,plgout_l) 
elseif(nx .eq. 0.0d0 .and. ny .lt. 0.0d0) then                                
    call polygoncopy(47,plgl,plgout_d)
     call polygoncopy(48,plgr,plgout_l)        
elseif(ny .eq. 0.0d0 .and. nx .gt. 0.0d0) then
    call polygoncopy(49,plgr,plgout_d)
     call polygoncopy(50,plgl,plgout_l) 
elseif(ny .eq. 0.0d0 .and. nx .lt. 0.0d0) then
     call polygoncopy(51,plgl,plgout_d)
     call polygoncopy(52,plgr,plgout_l)
else
    !print *, "mx,my,both 0"
!    plgout%nodesnum = 0
endif


call plgdel(plgl)
call plgdel(plgr)

end subroutine output_multi_polygon








!--------------------------------------------------------------------
subroutine plglistdel(plgs)
implicit none

type(plg_list)    :: plgs
integer           :: i

!if(plgs%num .ne. 0)then
   plgs%tag = 0
   do i = 1,plgs%num
     call plgdel(plgs%plg(i))
   enddo
   plgs%num = 0
   deallocate(plgs%plg)
!else
!   write(*,*) "the polygon list is already deleted"

!endif

end subroutine
!------------------------------------------------------------
subroutine plglistcopy(plgs1,plgs2)
implicit none

type(plg_list)      :: plgs1,plgs2
integer             :: i

plgs2%num = plgs1%num
plgs2%tag = plgs1%tag

allocate(plgs2%plg(plgs2%num))
do i = 1,plgs2%num
   call polygoncopy(60,plgs1%plg(i), plgs2%plg(i))
enddo

end subroutine
!---------------------------------------------------
subroutine output_plglist(plgs)
implicit none

type(plg_list)           :: plgs
integer                  :: i


!print *, "------------output polygon list-----------"
if (plgs%num .ne. 0) then
do i = 1, plgs%num
   !print *,  "plgs(1-plgs%num):", i
   call outputplg(plgs%plg(i))   
enddo
endif
!print *, "----------------list end------------------"


end subroutine
!----------------------------------------------------------
SUBROUTINE ult_cutpolygon(l_recon,plg,PLGs) 
! input l_recon,   plg
! output plg_list(plgs)

IMPLICIT NONE

TYPE(line),intent(in)                 :: L_recoN 
!INTEGER,intent(in)                    :: order
type(polygon)                         :: plg
integer                               :: i,j
real(kind=8),allocatable              :: vval(:)
type(plg_list)                       :: PLGs
type(polygon)                         :: plgl,plgr
integer                               :: flag,flag1
!  not deallocate  plg in yet

flag = 0
flag1 = 0

allocate(vval(plg%nodesnum))
do i = 1,plg%nodesnum
   call PinL(l_recon%val(1),l_recon%val(2),l_recon%val(3),plg%nodes(i),vval(i))
enddo

do i = 1,plg%nodesnum
   do j = 1,plg%nodesnum
     if(vval(i)*vval(j) .lt. 0.0d0) then
       flag = 1
       call cutpolygon(l_recon%val(1),l_recon%val(2),l_recon%val(3),plg,plgl,plgr)
       plgs%num = 2
       allocate(plgs%plg(2))
       call polygoncopy(71,plgl,plgs%plg(1))
       call polygoncopy(72,plgr,plgs%plg(2))
       call plgdel(plgl)
       call plgdel(plgr)
       flag1 = 1
       exit
     endif
   enddo
   if(flag1 .eq. 1) exit
enddo

if(flag .eq. 0) then
   plgs%num = 1
   allocate(plgs%plg(1))
   call polygoncopy(73,plg,plgs%plg(1))
   call plgdel(plg) 
endif


deallocate(vval)
END SUBROUTINE ult_cutpolygon




!-----------------------------------------------------------
subroutine plglist_to_trilist(list1,list2,N)
implicit none

type(plg_list),intent(in)     :: list1
type(plg_list),intent(out)    :: list2
integer                       :: i,j
integer,intent(out)           :: N
integer                       :: M,m1,m2


!initialize list 2
N = 0
do i = 1,list1%num
   if(list1%plg(i)%nodesnum .le. 2) then
      write(*,*) "no such polygon"
   endif
   N = N + (list1%plg(i)%nodesnum - 2)
enddo

list2%num = N
allocate(list2%plg(N))
do i = 1,list2%num
   list2%plg(i)%nodesnum = 3
   allocate(list2%plg(i)%nodes(3))
enddo
!--------------------------------

m = 0
do i = 1,list1%num
   m1 = list1%plg(i)%nodesnum - 2

   m2 = 2
   do j = m+1,m+m1
      list2%plg(j)%nodes(1)%val = list1%plg(i)%nodes(1)%val
      list2%plg(j)%nodes(2)%val = list1%plg(i)%nodes(m2)%val
      list2%plg(j)%nodes(3)%val = list1%plg(i)%nodes(m2+1)%val
      m2 = m2+1
   enddo
   m = m + m1

enddo



end subroutine

!-------------------------------------------------------------
!subroutine plg_to_trilist(plgin,listout)
!implicit none

!type(polygon),intent(in)     :: plgin
!type(plg_list)               :: listout
!integer                      :: i
!integer                      :: m

!if(plgin%nodesnum .le. 2) then
!  write(*,*) "no such polygon"

!else

!  listout%num = plgin%nodesnum - 2
!  allocate(listout%plg(listout%num))
!  do i = 1, listout%num  
!     listout%plg(i)%nodesnum = 3
!     allocate(listout%plg(i)%nodes(3))
!  enddo

!  m = 2
!  do i = 1,listout%num
!      listout%plg(i)%nodes(1)%val = plgin%nodes(1)%val
!      listout%plg(i)%nodes(2)%val = plgin%nodes(m)%val
!      listout%plg(i)%nodes(3)%val = plgin%nodes(m+1)%val
!      m = m+1
!  enddo

!endif

!end subroutine




!--------------------------------------------------------
subroutine plg_line_inter(plg,a,b,c,pout)
implicit none

type(polygon),intent(in)    :: PLG
real(kind=8),intent(in)   :: a,b,c
integer                   :: i,flag1,flag2,dupflag,flag3
type(points)              :: p
type(points)              :: pout(2)

flag2 = 0
flag3 = 0

do i = 1,plg%nodesnum-1 
   call line_line_inter(plg%nodes(i),plg%nodes(i+1),a,b,c,flag1,dupflag,p)
   if(flag1 .eq. 1 .and. flag2 .eq. 0 ) then
     pout(1)%val = p%val
     flag2 = 1
     flag3 = dupflag
   elseif(flag1 .eq. 1 .and. flag2 .eq. 1 .and. flag3 .eq. 1) then
     flag3 = 0
   elseif(flag1 .eq. 1 .and. flag2 .eq. 1 .and. flag3 .ne. 1)then
     pout(2)%val = p%val     
     exit
   endif
enddo

if(flag2 .eq. 0) then
  write(*,*) "warning: polygon_line intercept at only one point"
endif

end subroutine plg_line_inter


!-------------------------------------------------------
SUBROUTINE line_line_inter(p1,p2,a,b,c,flag,dupflag,pout)
implicit none

type(points),intent(in)    :: p1,p2
real(kind=8)               :: a, b, c
type(points),intent(out)   :: pout
real(kind=8)               :: x1,y1,x2,y2,x,y
real(kind=8)               :: m                        !slope
integer,intent(out)        :: flag,dupflag

x1 = p1%val(1)
y1 = p1%val(2)
x2 = p2%val(1)
y2 = p2%val(2)
flag = 0
dupflag = 0


if(abs(x2-x1) .gt. 0.0d0) then         ! not use afrac_eps as the difference
  m = (y2-y1)/(x2-x1)
  !-------------------------------------
  if(b .ne. 0.0d0) then
     if((m+a/b) .eq. 0.0d0) then
        write(*,*) "warning, denominator=0"
     endif
     x = (m*x1 - c/b - y1)/(m+a/b)
     y = (-a*x-c)/b
  elseif(a .ne. 0.0d0 .and. b .eq. 0.0d0) then
     x = -c/a
     y = y1 + m*(x - x1)       
  else 
     write(*,*) "a and b both 0 in line_line_inter"
  endif
else    ! where x1 = x2
  if(b .ne. 0.0d0) then
    x = x1
    y = (-c - a*x1)/b
  else
    write(*,*) "warning: the cutline is coincide with one edge"
  endif
endif

!--------------------------------------------
if( (x-x1)*(x-x2) .le. 0.0  .and.  (x-x2).ne. 0.0d0 ) then
    flag = 1
    pout%val(1) = x 
    pout%val(2) = y
elseif( (x-x1)*(x-x2) .le. 0.0  .and.  (x-x2).eq. 0.0d0 ) then
    dupflag = 1                                        ! dupflag = 1 
    flag = 1
    pout%val(1) = x 
    pout%val(2) = y    
else
    flag = 0
    pout%val = 0.0d0
endif
!--------------------------------------
end subroutine line_line_inter

!////////////////////////////////////////////////////////////
subroutine selec_m(mofdata,dmsn,nmat,order,a,b,c,m)
implicit none

integer     , intent(in) :: dmsn, nmat, order
real(kind=8), intent(in) :: mofdata((dmsn*2+3)*nmat)
integer                  :: i
real(kind=8)             :: a,b,c
integer,intent(out)      :: m

do i = 1, nmat
   IF(order .eq. mofdata((2+dmsn)+(i-1)*(dmsn*2+3)))then
     a = mofdata((2+dmsn+1)+(i-1)*(dmsn*2+3))
     b = mofdata((2+dmsn+2)+(i-1)*(dmsn*2+3))
     c = mofdata((2+dmsn+3)+(i-1)*(dmsn*2+3))

     m = i
     exit
   endif
enddo

end subroutine selec_m

!-------------------------------------------------------------------
SUBROUTINE INIT_F01(N,h,iin,jin,nox,noy,phi1,phi2,cell,vf,centroid)
IMPLICIT NONE

integer,parameter       :: nmat = 3
integer,parameter       :: dmsn = 2
integer,intent(in)      :: N
real(kind=8)            :: h
real(kind=8),intent(in) :: nox(N),noy(N)
integer                 :: iin,jin 
type(polygon)           :: cell(N,N)
real(kind=8),intent(in) :: phi1(N,N),phi2(N,N)
!INTEGER                 :: i,j
real(kind=8)            :: aout,bout,dout,C
TYPE(POLYGON)           :: POS_PLG,NEG_PLG,SUBPOS,SUBNEG
REAL(KIND=8)            :: VOL1,VOL2,vol3
INTEGER                 :: flag1,flag2
!type(POINTS)            :: centroid_1,centroid_2,centroid_3


!real(kind=8)            :: vfdata(M*(2+dmsn))
!   vfdata(...) = { tag, vol_frac, centroid    }
!type(ifseg)             :: seg1(N,N),seg2(N,N)
type(ifseg)              :: seg1,seg2
TYPE(POINTS),DIMENSION(N,N,nmat)     :: CENTROID
real(kind=8)                         :: vf(N,N,NMAT)

! initial vfdata
!vfdata = 0.0d0

 centroid(iin,jin,:)%val(1) = 0.0d0
 centroid(iin,jin,:)%val(2) = 0.0d0
  vf(iin,jin,:) = 0.0d0


POS_PLG%NODESNUM = 0
NEG_PLG%NODESNUM = 0
subpos%NODESNUM = 0
subneg%NODESNUM = 0
!centroid_1%VAL = 0.0d0
!centroid_2%VAL = 0.0d0
!centroid_3%VAL = 0.0d0

flag1 = 0
flag2 = 0

IF(iin .eq. 1  .or. iin .eq. N .or. jin .eq. 1 .or. jin .eq. N)THEN
   vf(iin,jin,3) = 1.0d0
   centroid(iin,jin,3)%val(1) = cell(iin,jin)%center%val(1)
   centroid(iin,jin,3)%val(2) = cell(iin,jin)%center%val(2)
ELSE

    CALL NB_CHECK(iin,jin,N,flag1,phi1)
    CALL NB_CHECK(iin,jin,N,flag2,phi2)

    IF (flag1 .eq. 1 .AND. FLAG2 .eq. 1) then
     ! write(*,*) "1" 
    ! first material
       CALL N_comp(N,h,iin,jin,phi1,nox,noy,aout,bout,dout)
        c = -aout*nox(iin)-bout*noy(jin)+dout

       call cell_line_inter(cell(iin,jin),aout,bout,c,seg1%PT)
!       print *, seg1%pt(1),seg1%pt(2)
                     
       CALL VOL_FRAC_CAL(cell(iin,jin),aout,bout,c,POS_PLG,NEG_PLG)
       CALL POLY_CENTROID(POS_PLG,centroid(iin,jin,1))       

     ! second material
     ! aout bout dout using here again       
       CALL N_comp(N,h,iin,jin,phi2,nox,noy,aout,bout,dout)
        c = -aout*nox(iin)-bout*noy(jin)+dout

       call cell_line_inter(CELL(IIN,JIN),aout,bout,c,seg2%PT)
!       print *, seg2%pt(1),seg2%pt(2)
                     
       CALL VOL_FRAC_CAL(NEG_PLG,aout,bout,c,SUBPOS,SUBNEG)
        write(*,*) "pos",subpos%nodesnum
        write(*,*) "neg",subneg%nodesnum
         
       CALL POLY_CENTROID(SUBPOS,centroid(iin,jin,2))         
       CALL POLY_CENTROID(SUBNEG,centroid(iin,jin,3))         
!       vfdata(3) = centroid_1%val(1)
!       vfdata(4) = centroid_1%val(2)
!       vfdata(7) = centroid_2%val(1)
!       vfdata(8) = centroid_2%val(2)
!       vfdata(11) = centroid_3%val(1)
!       vfdata(12) = centroid_3%val(2)
  
       CALL poly_area(POS_PLG,vol1)
       CALL poly_area(subpos,vol2)
       call poly_area(subneg,vol3)
       vf(iin,jin,1)= VOL1/(h**2.0d0)
       vf(iin,jin,2) = vol2/(h**2.0d0)
       vf(iin,jin,3) = vol3/(h**2.0d0)
       
!       vfdata(1)  = 1.0d0
!       vfdata(5)  = 1.0d0
!       vfdata(9)  = 1.0d0

       CALL PLGDEL(POS_PLG)
       CALL PLGDEL(NEG_PLG)
       CALL PLGDEL(SUBPOS)
       CALL PLGDEL(SUBNEG) 
    ELSEIF(FLAG1 .EQ. -1 .AND. PHI1(IIN,JIN) .LT. 0.0D0) THEN
!       vfdata(1) = 1.0D0
       vf(iin,jin,1) = 1.0d0
       centroid(iin,jin,1)%val(1) = cell(iin,jin)%center%val(1)
       centroid(iin,jin,1)%val(2) = cell(iin,jin)%center%val(2)
       
    ELSEIF(FLAG1 .EQ. -1 .AND. PHI1(IIN,JIN) .GT. 0.0D0 &
            &.and.  flag2 .eq. 1) THEN
     !  write(*,*)  "2" 
       CALL N_comp(N,h,iin,jin,phi2,nox,noy,aout,bout,dout)
        c = -aout*nox(iin)-bout*noy(jin)+dout

       call cell_line_inter(CELL(IIN,JIN),aout,bout,c,seg2%PT)
!       print *, seg2%pt(1),seg2%pt(2)                    

       CALL VOL_FRAC_CAL(CELL(IIN,JIN),aout,bout,c,POS_PLG,NEG_PLG)
       CALL POLY_CENTROID(POS_PLG,centroid(iin,jin,2))         
       CALL POLY_CENTROID(NEG_PLG,centroid(iin,jin,3))         
!       vfdata(7) = centroid_1%val(1)
!       vfdata(8) = centroid_1%val(2)
!       vfdata(11) = centroid_2%val(1)
!       vfdata(12) = centroid_2%val(2)
 
       CALL poly_area(POS_PLG,vol1)
       CALL poly_area(neg_plg,vol2)
       vf(iin,jin,2)= VOL1/(h**2.0d0)
       vf(iin,jin,3) = vol2/(h**2.0d0)
      
!       vfdata(5) = 1.0d0
!       vfdata(9) = 1.0d0
       CALL PLGDEL(POS_PLG)
       CALL PLGDEL(NEG_PLG)
     
    ELSEIF(FLAG1 .EQ. 1 .and.  flag2 .eq. -1) THEN
     !  write(*,*)  "3" 
       CALL N_comp(N,h,iin,jin,phi1,nox,noy,aout,bout,dout)
        c = -aout*nox(iin)-bout*noy(jin)+dout

       call cell_line_inter(CELL(IIN,JIN),aout,bout,c,seg1%PT)
       print *, seg1%pt(1),seg1%pt(2)                    

       CALL VOL_FRAC_CAL(CELL(IIN,JIN),aout,bout,c,POS_PLG,NEG_PLG)
       CALL POLY_CENTROID(POS_PLG,centroid(iin,jin,1))         
       CALL POLY_CENTROID(NEG_PLG,centroid(iin,jin,2))         
!       vfdata(3) = centroid_1%val(1)
!       vfdata(4) = centroid_1%val(2)
!       vfdata(7) = centroid_2%val(1)
!       vfdata(8) = centroid_2%val(2)
 
       CALL poly_area(POS_PLG,vol1)
       CALL poly_area(neg_plg,vol2)
       vf(iin,jin,1)= VOL1/(h**2.0d0)
       vf(iin,jin,2) = vol2/(h**2.0d0)
      
!       vfdata(1) = 1.0d0
!       vfdata(5) = 1.0d0
       CALL PLGDEL(POS_PLG)
       CALL PLGDEL(NEG_PLG)

    ELSEIF(FLAG1 .EQ. -1 .AND. PHI1(IIN,JIN) .GT. 0.0D0 &
            &.and.  flag2 .eq. -1 .AND. PHI2(IIN,JIN) .GT. 0.0D0) THEN
!       vfdata(9) = 1.0d0
       vf(iin,jin,3) = 1.0d0
       centroid(iin,jin,3)%val(1) = cell(iin,jin)%center%val(1)
       centroid(iin,jin,3)%val(2) = cell(iin,jin)%center%val(2)
    ENDIF
 ENDIF


END SUBROUTINE INIT_F01

!-----------------------------------------------------------
SUBROUTINE circledist(dist,x,y,r)
IMPLICIT NONE


REAL(KIND=8) :: dist,x,y,r


   dist = sqrt( (x-50.0d0)**2.0d0 + (y-50.0d0)**2.0d0) - r



END SUBROUTINE circledist



!////////////////////////////////////////////////////////////
!   <SUBROUTINE> SPACE_LOOP
SUBROUTINE SPACE_LOOP(dmsn,nmat,vf,centroid,CELL,X,Y,T_current,N,M,TAU,h)
IMPLICIT NONE


integer                   ,intent(in) :: nmat                  ! number of material
integer                   ,intent(in) :: dmsn                  ! dimension
integer                               :: datalen

INTEGER                  ,INTENT(IN) :: N,M                   ! num of cell(i direction), time step
REAL(KIND=8)             ,INTENT(IN) :: TAU,h                  
REAL(KIND=8)             ,INTENT(IN) :: X(N+1),Y(N+1)
REAL(KIND=8)             ,INTENT(IN) :: T_current
TYPE(POLYGON)            ,INTENT(IN) :: CELL(N,N)


REAL(KIND=8)                         :: mes(N,N)
type(polygon)                        :: fic
real(kind=8)                         :: vol(N,N,nmat)
real(kind=8)                         :: Cxtemp(n,n,nmat), Cytemp(n,n,nmat)

type(polygon)                        :: comp_area(N,N,nmat)
!real(kind=8)                         :: area_temp
TYPE(POLYGON),ALLOCATABLE            :: POLYGONTEMP(:)
TYPE(POLYGON),ALLOCATABLE            :: PolygonOut(:)
INTEGER                              :: RX,RY
INTEGER                              :: i,j,k,i1,i2,j1,im
real(kind=8)                         :: a,b,c
integer                              :: ii,jj
INTEGER                              :: lowbound,upbound
INTEGER                              :: leftbound,rightbound
real(kind=8)                         :: area_temp
real(kind=8)                         :: alpha
type(points)                         :: pout(2)
type(points)                         :: centroid_temp
!integer                              :: cen_flag(n,n)

TYPE(POINTS),DIMENSION(N,N,nmat)     :: CENTROID
real(kind=8)                         :: vf(N,N,NMAT)

! MOF value transfer
real(kind=8)                          :: xtetlist(dmsn+1,dmsn,2000)
real(kind=8)                          :: mofdata((dmsn*2+3)*nmat)
type(polygon)                         :: plgout_d,plgout_l
real(kind=8)                          :: dx(dmsn)
integer                               :: mof_verbose

TYPE(line)                            :: L_recon(N,N,nmat) 
TYPE(PLG_LIST)                        :: PLGS(N,N,nmat)
TYPE(PLG_LIST)                        :: PLGS_TEMP1,PLGS_TEMP2,TRI_TEMP
real(kind=8)                          :: area_test, area_test_tot
integer                               :: finecut_flag
integer                               :: n_temp


type(polygon)                         :: tbd_plg
real(kind=8)                          :: vf_temp
integer                               :: c_order
real(kind=8)                          :: vflist(nmat) 
integer                               :: orderlist(nmat)
integer                               :: order_out
real(kind=8)                          :: max_out
integer                               :: mat


datalen = dmsn*2+3
xtetlist(:,:,:) = 0.0d0 
dx(:) = h
mes(:,:) = h**(2.0d0)


!***********************************************
 Cxtemp(:,:,:) = 0.0d0
 Cytemp(:,:,:) = 0.0d0
 vol(:,:,:) = 0.0d0
 comp_area(:,:,:)%nodesnum = 0
 l_recon(:,:,:)%tag = 0
 

print *,"i=0,N-1 and j=0,N-1 space loop"
stop

!____________________________________________________________
!print *, "step1"
do i = 1,N
  do j = 1, N
 !     print *, centroid_d(i,j)
 !   mofdata 
 !   material 1  F
 !               centroid(2)
 !               order
 !               slope(2)
 !               intercept
 !     n*(x - centroid(2)) + intercept = 0
     print *, "i=",i,"j=",j

     ! initial mofdata
  
     do im = 1,nmat
       mofdata(1+(im-1)*datalen) = vf(i,j,im)
       mofdata(2+(im-1)*datalen) = centroid(i,j,im)%val(1)-cell(i,j)%center%val(1)
       mofdata(3+(im-1)*datalen) = centroid(i,j,im)%val(2)-cell(i,j)%center%val(2)
       mofdata(4+(im-1)*datalen) = 0.0d0
       mofdata(5+(im-1)*datalen) = 0.0d0
       mofdata(6+(im-1)*datalen) = 0.0d0
       mofdata(7+(im-1)*datalen) = 0.0d0 
     enddo

     print *, "%%%%%%%%%%%%%%%%%%%%%%%%"
     print *, "mofdata"
     print *, mofdata
     print *, "%%%%%%%%%%%%%%%%%%%%%%%%"

      mof_verbose = 0
     ! call multimaterial_MOF( &
     !                         mof_verbose, &
     !                         xtetlist, &
     !                         xtetlist, &
     !                         1000, &
     !                         cell(i,j)%center%val, & 
     !                         dx, &
     !                         mofdata,&
     !                         multi_centroidA, &  
     !                         1.0d0, &
     !                         0, &                      ! continuous MOF = 0
     !                         cell(i,j)%center%val,&
     !                         dx, &
     !                         0, &
     !                         nmat, &
     !                         2) 
      ! mofdata updated
      ! 
      ! 
      do k = 1,nmat
         comp_area(i,j,k)%nodesnum = 0              
      enddo
      !initialize
      call polygoncopy(74,cell(i,j),tbd_plg)                 ! deallocate     

      vf_temp = 1.0d0
      c_order = 1
      mat = 0

      do  i1 =1,nmat
         vflist(i1) = mofdata(1+(i1-1)*datalen)                   ! 2-D
      enddo
      print *, "vflist:"
      print *, vflist
 
      do i1 = 1,nmat
         orderlist(i1) = mofdata((2+dmsn)+(i1-1)*datalen)
      enddo
      print *, "orderlist:"
      print *, orderlist

      CALL findmax(nmat,vflist,order_out,max_out)     
      print *, "order_out"
      print *, order_out
      print *, "max_out"
      print *, max_out

     
      if(max_out .ge. 1.0d0-afrac_eps)then
         call polygoncopy(75,cell(i,j), comp_area(i,j,order_out))             
         call plgdel(tbd_plg)
         vf_temp = 0.0d0
      else

         do while(tbd_plg%nodesnum .ne. 0)
            call selec_m(mofdata,dmsn,nmat,c_order,a,b,c,mat)
            print *, "c_order",c_order
            print *, "a=",a,"b=",b,"c=",c
            print *, "mat=",mat

           if(a .ne. 0.0d0 .or. b .ne. 0.0d0) then
            print *, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

            call output_multi_polygon(a,b,c,h,N,i,j,cell,tbd_plg,mes,&
                                plgout_d,plgout_l,alpha)
            print *, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
            
            call plg_line_inter(tbd_plg,a,b,-alpha,pout)
            print *, pout(1)%val,pout(2)%val

            l_recon(i,j,mat)%tag = 1      
            l_recon(i,j,mat)%val(1) = a
            l_recon(i,j,mat)%val(2) = b
            l_recon(i,j,mat)%val(3) = alpha
            
            call polygoncopy(76,plgout_d,comp_area(i,j,mat))
            call plgdel(plgout_d)
            call plgdel(tbd_plg) 

              call poly_area(comp_area(i,j,mat),area_temp)
              vf_temp = vf_temp - area_temp/mes(i,j)
  
              if(plgout_l%nodesnum .gt. 2) then                                     
                call polygoncopy(77,plgout_l,tbd_plg)
                call plgdel(plgout_l)
              elseif(plgout_l%nodesnum .eq. 2 .or. plgout_l%nodesnum .eq. 1) then
                call plgdel(plgout_l)
                write(*,*) "suspicious number of nodes"
          !    else
          !      exit
              endif
            endif

            c_order = c_order + 1

           ! if(c_order .gt. nmat) then
           !   exit
           ! endif

          enddo
       endif


     ! print *, "i=",i,"j=",j,"volume fraction left:",vf_temp
     
    ! print *, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
         !   call polygoncopy(10,plgout_d, comp_area(i,j))



       !print *, "MOF reconstruction>>>>>>>>>>>>>>"     
       !print *, "for cell", i,j
       !print *, "volume fraction=", vf(i,j)
       !print *, "cell_center", cell(i,j)%center%val
       !print *, "centroid_dark", centroid_d(i,j)%val
       !print *, "centroid_light",centroid_l(i,j)%val

       !print *, "$$$$$$$$$$$$$$$$$$$"
      ! do k = 1,14
          !print *, "mofdata",mofdata(k)
      ! enddo
      !print *, "$$$$$$$$$$$$$$$$$4$"
      !print *, "mx=",mofdata(5) 
      !print *,  "my=", mofdata(6) 
      !print *,  "intercept", mofdata(7)
   ENDDO
 ENDDO

     !print *, "l_recon:", i,j
     !print *, l_recon(i,j)%val
!do i = 1,N
!   print *, cell_tag(i,:)
   
!enddo

area_test_tot = 0.0d0

print *,"i=0,N-1 and j=0,N-1 space loop"
stop

!print *, "sum of step 1:"
do i=1,N
  do j = 1,N
     do k = 1,nmat
  
    !print *, "i=",i,"j=",j   
!    call outputplg(comp_area(i,j))
         if(comp_area(i,j,k)%nodesnum .ne. 0.0d0)then
          call poly_area(comp_area(i,j,k),area_test)
         endif
    !print *, "comp_area,before cut:", area_test
        area_test_tot = area_test_tot + area_test
     enddo 
  enddo
enddo

print *, "total area:",  area_test_tot

print *,"i=0,N-1 and j=0,N-1 space loop"
stop

do i=1,N
   do j = 1,N
     do k = 1,nmat
      plgs(i,j,k)%num = 0
     enddo
   enddo
enddo

print *,"i=1,N-2 and j=1,N-2 space loop"
stop

do i=2,N-1
  do j =2,N-1
    do k = 1,nmat
     if(comp_area(i,j,k)%nodesnum .ne. 0.0) then
       !print *, "cell working on:",i,j
       PLGS(i,j,k)%num = 1                        ! start from plgs(i,j)%plg(1)
       allocate(plgs(i,j,k)%plg(1))
       call polygoncopy(78,comp_area(i,j,k),plgs(i,j,k)%plg(1))   
                                               ! not yet deallocate comp_area


       do ii=-1,1                             ! check the eight neighbour cells
          do jj=-1,1
             do im = 1,nmat
             PLGS_TEMP2%NUM = 0
             ALLOCATE(PLGS_TEMP2%PLG(60*nmat))
             finecut_flag = 0
  
            if(ii .ne. 0 .and. jj .ne. 0)then
                
               !----------------------------------------------
                if(l_recon(i+ii,j+jj,im)%tag .eq. 1) then
                   finecut_flag = 1
                   do i1 = 1,plgs(i,j,k)%num             !! for each polygon list 
                      call ult_cutpolygon(L_RECON(i+ii,j+jj,im), &
                                  plgs(i,j,k)%plg(i1),PLGS_TEMP1)
                      call plglist_to_trilist(plgs_temp1,tri_temp,n_temp)
                      do j1 = 1, tri_temp%num     
                         call polygoncopy(79,tri_temp%plg(j1),plgs_temp2%plg(plgs_temp2%num+1))
                         PLGS_temp2%num = plgs_temp2%num + 1
                      enddo
                     !----------------------------------------
!                      IF(PLGS_TEMP1%NUM .EQ. 1) THEN                    
!                         CALL polygoncopy(10,PLGS_TEMP1%PLG(1), &
!                                  PLGS_temp2%plg(plgs_temp2%num+1))  
!                         PLGS_temp2%num = plgs_temp2%num + 1
!                      ELSEIF(Plgs_temp1%num .eq. 2) then
!                         CALL polygoncopy(10,PLGS_TEMP1%PLG(1), &
!                                   PLGS_temp2%plg(plgs_temp2%num+1))
!                         CALL polygoncopy(10,PLGS_TEMP1%PLG(2), &
!                                   PLGS_temp2%plg(plgs_temp2%num+2))
!                         PLGS_temp2%num = plgs_temp2%num + 2
!                      ELSE
!                         WRITE(*,*) "Wrong num of polygons"                          
!                      ENDIF
                      
                       call plglistdel(plgs_temp1)
                       call plglistdel(tri_temp)
                     !-----------------------------------
                    enddo               
                  endif
                !---------------------------------------------

                if(finecut_flag .eq. 1) then
                 call plglistdel(plgs(i,j,k))
                 call plglistcopy(plgs_temp2,plgs(i,j,k))
                endif

            endif

             call plglistdel(plgs_temp2)
            enddo
          enddo
       enddo


     endif
   enddo
  enddo
enddo



!area_test_tot = 0.0d0
print *,"i=0,N-1 and j=0,N-1 space loop"
stop

!do i= 1,N
!   do j = 1,N
      !print *, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      !print *, "polygon_list", "i=",i,"j",j
!      do i1 = 1,plgs(i,j)%num
!         call poly_area(plgs(i,j)%plg(i1), area_test)
!         area_test_tot = area_test_tot + area_test
!      enddo
  !print *, "tot area after further cut for each comp area:", &
      !print *, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
!   enddo
!enddo
!print *, "end of step 2"





!-------------------------------------------------------------------
!      output triangles of plgs(i,j)
!-------------------------------------------------------------------
print *,"i=0,N-1 and j=0,N-1 space loop"
stop
do i = 1,N
  do j = 1,N
    do k = 1,nmat
      do i1 = 1,plgs(i,j,k)%num
         print *, plgs(i,j,k)%plg(i1)%nodes(1), &
               plgs(i,j,k)%plg(i1)%nodes(2), &
               plgs(i,j,k)%plg(i1)%nodes(3)
      enddo
    enddo
  enddo
enddo


!--------------------------------------------------------------------







!print *, "step 3 "

!__________________________________________________________
!>>>>advect each polygon in plgs(i,j) for each (i,j)<<<<<
!__________________________________________________________

DO im = 1,nmat

print *,"i=0,N-1 and j=0,N-1 space loop"
stop
do i= 1,N
  do j= 1,N
      upbound =0
      lowbound = 0
      leftbound = 0
      rightbound = 0     
    do ii = 1,plgs(i,j,im)%num 
      if (plgs(i,j,im)%plg(ii)%nodesnum .ne. 0.0) then
        fic%nodesnum = plgs(i,j,im)%plg(ii)%nodesnum                     !deallocate
        allocate(fic%nodes(fic%nodesnum))
        do k = 1, fic%nodesnum
          CALL f_RK(plgs(i,j,im)%plg(ii)%NODES(k)%VAL(1),plgs(i,j,im)%plg(ii)%NODES(k)%VAL(2)&
                    ,T_current,tau,FIC%NODES(k)%VAL(1)&
                    ,FIC%NODES(k)%VAL(2))
        enddo
  
        CALL HorizontalMeshCut(FIC,X,Y,N,N,POLYGONTEMP,lowbound,upbound)
        RY= upbound-lowbound
        do i1 = 1,RY
          CALL VerticalMeshCut(polygontemp(i1),X,Y,N,N,leftbound,rightbound,PolygonOut)
          RX= rightbound-leftbound
          do j1=1,RX
            call poly_area(polygonout(j1),area_temp)
            vol(leftbound+j1-1,lowbound+i1-1,im) = vol(leftbound+j1-1,lowbound+i1-1,im)+ area_temp
            call poly_centroid(polygonout(j1),centroid_temp)
             Cxtemp(leftbound+j1-1,lowbound+i1-1,im) = &
               Cxtemp(leftbound+j1-1,lowbound+i1-1,im)+ area_temp*centroid_temp%val(1)
             Cytemp(leftbound+j1-1,lowbound+i1-1,im) = &
              Cytemp(leftbound+j1-1,lowbound+i1-1,im)+ area_temp*centroid_temp%val(2)       
          enddo
          DO i2=1,RX
             CALL PLGDEL(POLYGONOUT(i2))
          ENDDO
          DEALLOCATE (PolygonOut)     
        enddo
        call plgdel(fic)
      endif
    enddo

  enddo
enddo

enddo
print *,"i=0,N-1 and j=0,N-1"
stop

do im = 1,nmat
print *,"i=0,N-1 and j=0,N-1 space loop"
stop
do i=1,N
   do j=1,N
     if(vol(i,j,im) .ge. afrac_eps) then
!      print *, i,j,"mes=",mes(i,j),"volume=",vol(i,j)
      vf(i,j,im) = vol(i,j,im)/mes(i,j)
      Centroid(i,j,im)%val(1) = Cxtemp(i,j,im)/vol(i,j,im)
      Centroid(i,j,im)%val(2) = Cytemp(i,j,im)/vol(i,j,im)

!      Centroid_l(i,j,im)%val(1) = (mes(i,j)*cell(i,j,im)%center%val(1) - &
!                 centroid_d(i,j)%val(1)*vol(i,j))/(mes(i,j)-vol(i,j))
!      Centroid_l(i,j)%val(2) = (mes(i,j)*cell(i,j)%center%val(2) - &
!                 centroid_d(i,j)%val(2)*vol(i,j))/(mes(i,j)-vol(i,j))      

!       print *, "Cxtemp=",Cxtemp(i,j),"Cytemp=", Cytemp(i,j)
!       print *, centroid_d(i,j)%val
       !print *, "centroid_dark", centroid_d(i,j)%val
       !print *, "centroid_light",centroid_l(i,j)%val


      else
        vf(i,j,im) = 0.0d0
        centroid(i,j,im)%val(1) = 0.0d0
        centroid(i,j,im)%val(2) = 0.0d0
      endif

   enddo
enddo
enddo




!print *, ">>>>>>>>>>>>>>>"
!do  i = 1,N
!   print *, centroid_l(i,:)
!enddo

print *,"i=0,N-1 and j=0,N-1 space loop"
stop


do i=1,N
   do j=1,N
      do k =1,nmat
        call plgdel(comp_area(i,j,k))
      enddo
   enddo
enddo

print *,"i=0,N-1 and j=0,N-1 space loop"
stop
!do i = 1, N
!  print *, vf(i,:)
!enddo

!print *, "end"

END SUBROUTINE SPACE_LOOP

!----------------------------------------------------------
subroutine findmax(N,x,order_out,max_out)
implicit none

integer             :: N
real(kind=8)        :: x(N)
integer             :: i
integer             :: order
real(kind=8)        :: val
integer             :: order_out
real(kind=8)        :: max_out

  
  order = 1
  val = x(1)
print *,"i=1,N-1 findmax?"
stop
  do  i = 2,N
    if(x(i) .gt. x(i-1))then
      order = i
      val = x(i)
    endif       
  enddo
  order_out = order
  max_out = val  

end subroutine findmax




END MODULE GeneralCLass








