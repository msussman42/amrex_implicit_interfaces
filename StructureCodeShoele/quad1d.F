c  August 31, 2023 grid_def.F changed
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE quad1d(Cor,Wi,Intord)
c------------------------------------
c 	 Returns Gauss coordinates and Weights for up to 8 Gauss points
c------------------------------------
       IMPLICIT NONE
       include 'grid_def'
       REAL*8  Wi(8),Cor(8)      !       weigths
       INTEGER Intord	 !	 integration order
      print*, igauh
       Intord=igauh
       SELECT CASE (Intord)
       CASE (1)
          Cor(1)= 0.
          Wi(1) = 2.0
       CASE(2)
          Cor(1)= .577350269	; Cor(2)= -Cor(1)
          Wi(1) = 1.0 ;  Wi(2) = Wi(1)
       CASE(3)
          Cor(1)= .774596669	; Cor(2)= 0.0 ; Cor(3)= -Cor(1)
          Wi(1) = .555555555	; Wi(2) = .888888888 ; Wi(3) = Wi(1)
       CASE(4)
          Cor(1)= .861136311 ; Cor(2)= .339981043
          Cor(3)= -Cor(2) ; Cor(4)= -Cor(1)
          Wi(1) = .347854845 ; Wi(2) = .652145154
           Wi(3) = Wi(2) ; Wi(4) = Wi(1)
       CASE(5)
          Cor(1)= .9061798459 ; Cor(2)= .5384693101
           Cor(3)= .0 ; Cor(4)= -Cor(2)
          Cor(5)= -Cor(1)
          Wi(1) = .236926885 ; Wi(2) = .478628670
          Wi(3) = .568888888 ; Wi(4) = Wi(2)
          Wi(5) = Wi(1)
       CASE(6)
          Cor(1)= .932469514 ; Cor(2)= .661209386 ; Cor(3)= .238619186
          Cor(4)= -Cor(3) ;  Cor(5)= -Cor(2) ; Cor(6)= -Cor(1)
          Wi(1) = .171324492 ; Wi(2) = .360761573 ; Wi(3) = .467913934
          Wi(4) = Wi(3) ; Wi(5) = Wi(2) ; Wi(6) = Wi(1)
       CASE(7)
          Cor(1)= .949107912 ; Cor(2)= .741531185 ; Cor(3)= .405845151
          Cor(4)= 0.
          Cor(5)= -Cor(3) ;Cor(6)= -Cor(2) ;Cor(7)= -Cor(1)
          Wi(1) = .129484966 ; Wi(2) = .279705391 ; Wi(3) = .381830050
          Wi(4) = .417959183
          Wi(5) = Wi(3) ; Wi(6) = Wi(2) ; Wi(7) = Wi(1)
		   CASE(8)
          Cor(1)= .960289856 ; Cor(2)= .796666477
          Cor(3)= .525532409 ; Cor(4)= .183434642
          Cor(5)= -Cor(4) ; Cor(6)= -Cor(3)
          Cor(7)= -Cor(2) ; Cor(8)= -Cor(1)
          Wi(1) = .101228536 ; Wi(2) = .222381034
          Wi(3) = .313706645 ;Wi(4) = .362683783
          Wi(5) = Wi(4) ; Wi(6) = Wi(3) ; Wi(7) = Wi(2) ; Wi(8) = Wi(1)
	     CASE DEFAULT
       write(*,*) 'Gauss points not in range 1-8'
       write(*,*) 'Please chaneg igauh to <8'       
       STOP
       END SELECT
       RETURN
       END 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
