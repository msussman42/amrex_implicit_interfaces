      program main
      IMPLICIT NONE

      real*8,parameter      :: macrolayer = 2.8E-4
      real*8,parameter      :: microlayer = 1.0E-9
      real*8 factor

      factor=log(macrolayer/microlayer)/(macrolayer-microlayer)
      print *,"macrolayer=",macrolayer
      print *,"microlayer=",microlayer
      print *,"factor= ",factor
 
      return
      end

