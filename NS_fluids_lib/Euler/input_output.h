      integer fid,fid_l,fid_r,fid_xr,fid_xm,fid_xiL,fid_xiR,fid_p
      integer fid_rm,fid_iLiR,fid_pu,fid_pr,fid_xv,fid_xp
      integer fid_xiS
      parameter( fid=11, fid_l=12, fid_r=13, fid_xr=14,
     &     fid_xm=15, fid_xiL=16, fid_xiR=17, fid_p=18,
     &     fid_rm=21, fid_iLiR=22, fid_pu=23, fid_pr=24,
     &     fid_xv=25, fid_xp=26, fid_xiS=29 )

c-----------------------------------------------------------------------
c     print flags  
c-----------------------------------------------------------------------
c
      logical    PRINT_NEWTON
      logical    PRINT_NEWTON_CONV
      logical    PRINT_RP
      logical    PRINT_WC
      logical    PRINT_INIT_WC
      common /print_flags/ PRINT_NEWTON,
     &                     PRINT_NEWTON_CONV,
     &                     PRINT_RP,
     &                     PRINT_WC, 
     &                     PRINT_INIT_WC
