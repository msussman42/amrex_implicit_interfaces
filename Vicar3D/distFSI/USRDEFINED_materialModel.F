c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine USRDEFINED_materialModel(
     &            ibdyg,Phi_mem_coef,gmetric_con0 
     &           ,CauchyGreen_C,CauchyGreen_Cbar
     &           ,detC,tracC,detC_inplane,tracC_inplane
     &           ,Fibrous_Esh,MatposFib_con
     &           ,DiffPhi_C,DiffPhi_C33
     &           ,Diff2Phi_C,Diff2Phi_CC33,Diff2Phi_C33 )    
      
      implicit none
      integer ibdyg
      real*8      Phi_mem_coef(15)
     &     ,gmetric_con0(2,2)
     &     ,CauchyGreen_C(3,3)
     &     ,CauchyGreen_Cbar(3,3)
     &     ,detC,tracC,detC_inplane,tracC_inplane
     &     ,Fibrous_Esh
     &     ,MatposFib_con(3,3)
     &     ,DiffPhi_C(3)
!     OUTPUT
     &     ,DiffPhi_C33
!     OUTPUT
     &     ,Diff2Phi_C(3,3)
!     OUTPUT
     &     ,Diff2Phi_CC33(3)
!     OUTPUT
     &     ,Diff2Phi_C33
!     OUTPUT
! Ref: J. Kiendl et al. / Comput. Methods Appl. Mech. Engrg. 291 (2015) 280–303
      return
      end  subroutine 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
