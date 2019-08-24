      SUBROUTINE blstr2int(iarr, n, str)
      CHARACTER*(*) str
      INTEGER n, i, j
      INTEGER iarr(n)
      INTEGER EOS
      PARAMETER (EOS=-1)
 
      IF ( n .LE. len(str) ) THEN
          call bl_abort("blstr2int: str to large for iarr")
      END IF
      DO J = 1, N
          iarr(J) = ichar(' ')
      END DO
      j = 1
      DO i = 1, len(str)
          iarr(j) = ichar(str(i:i))
          j = j + 1
      END DO
      iarr(j) = EOS
 
      END
      SUBROUTINE blint2str(str, iarr, n)
      CHARACTER*(*) str
      INTEGER n
      INTEGER iarr(n)
      INTEGER EOS
      PARAMETER (EOS=-1)
      INTEGER i
 
      DO i = 1, LEN(str)
          str(i:i) = ' '
      END DO
      DO i = 1, n
          IF ( i .GT. LEN(str) ) then
             call bl_abort("blint2str: iarr to large for str")
          end if
          IF ( iarr(i) .EQ. EOS ) GO TO 100
          str(i:i) = char(iarr(i))
      END DO
 100  CONTINUE
 
      END
