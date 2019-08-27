      real*8 timeInit
      integer nStateInit
      real*8 cInit(0:N_INIT_STATE,1:N_EQUATIONS)
      real*8 uInit(0:N_INIT_STATE,1:N_EXT_VARS)
      real*8 xInit(0:N_INIT_DISC)
      real*8 wInit(0:N_INIT_DISC)
      real*8 dxInit
      common /exact_pars/ timeInit,
     &                    cInit,uInit,xInit,wInit,dxInit,
     &                    nStateInit
      integer lCSt, rCSt
      common /const_st_indx/ lCSt, rCSt
