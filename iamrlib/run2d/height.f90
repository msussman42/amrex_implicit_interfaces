      subroutine get_bottom_height(x,height)

      IMPLICIT NONE
      INTEGER,PARAMETER :: DP=SELECTED_REAL_KIND(15)
      REAL(KIND=DP),INTENT(IN) :: x
      REAL(KIND=DP) :: h,l
      REAL(KIND=DP),INTENT(OUT) :: height
      !x/h=[-52,44] z/h=[0,5]
      h=0.1143
      l=2.5*h
      height=h*(1-2*x**2/l**2+x**4/l**4)

      end subroutine get_bottom_height








