      integer    N_EQUATIONS
      parameter( N_EQUATIONS=3 )

      integer    N_WAVES
      parameter( N_WAVES=3 )

      integer    N_EXT_VARS
      parameter( N_EXT_VARS=4 )

      integer    N_STATES
      parameter( N_STATES=N_WAVES+1 )

      integer    N_EDGES
      parameter( N_EDGES=2*N_WAVES )

c------------------------------------------------------------------------

c     Directions.
      integer    DIR_L,    DIR_R,   DIR_M
      parameter( DIR_L=-1, DIR_R=+1,  DIR_M=0)

c     Input Side
      integer    BEHIND, AHEAD
      parameter( BEHIND=-1, AHEAD=1 )

c     Significance
      integer    INSIGNIFICANT,   SIGNIFICANT
      parameter( INSIGNIFICANT=0, SIGNIFICANT=1 )

c------------------------------------------------------------------------

      integer WAVE_TO_DIR(N_WAVES)
      common /WAVE_TO_DIRECTION/ WAVE_TO_DIR

      integer WAVE_TO_EDGE_L(N_WAVES), WAVE_TO_EDGE_R(N_WAVES)
      common /wave_to_edge/ WAVE_TO_EDGE_L, WAVE_TO_EDGE_R

      integer WAVE_TO_EDGE_A(N_WAVES), WAVE_TO_EDGE_B(N_WAVES)
      common /wave_to_edge_ab/ WAVE_TO_EDGE_A, WAVE_TO_EDGE_B

      integer WAVE_TO_STATE_L(N_WAVES), WAVE_TO_STATE_R(N_WAVES)
      common /wave_to_state/ WAVE_TO_STATE_L, WAVE_TO_STATE_R

      integer WAVE_TO_STATE_A(N_WAVES), WAVE_TO_STATE_B(N_WAVES)
      common /wave_to_state_ab/ WAVE_TO_STATE_A, WAVE_TO_STATE_B

c-------------------------------------------------------------------------
c     Used Only in slope2state.f
c-------------------------------------------------------------------------

c     Number of edges plus one.
      integer    N_EDGES_P1
      parameter( N_EDGES_P1=N_EDGES+1 )

      logical IS_L_EDGE(N_EDGES_P1)
      common /left_edges/ IS_L_EDGE

      integer EDGE_TO_L_STATE(N_EDGES_P1)
      common /left_states/ EDGE_TO_L_STATE

      integer EDGE_TO_L_WAVE(N_EDGES_P1)
      common /left_waves_GRP/ EDGE_TO_L_WAVE 

c-------------------------------------------------------------------------
c     Wave Types:
c-------------------------------------------------------------------------

      integer    N_WAVE_TYPES
      parameter( N_WAVE_TYPES=3 )

c     Wave type names.  
      integer    RFN,    SFR,    CD
      parameter( RFN=1,  SFR=2,  CD=3 )

c-----------------------------------------------------------------------
c     Names of primative variables:
c 
      integer    IRHO,   IVEL,   IPRE,   ISND
      parameter( IRHO=1, IVEL=2, IPRE=3, ISND=4 )

c-----------------------------------------------------------------------
c     Names of conservative variables:
c
      integer    IMOM,   IENE
      parameter( IMOM=2, IENE=3 )

c--------------------------------------------------------------------------

      character*8 WAVE_TYPE_STR(N_WAVE_TYPES)
      common /wave_types/ WAVE_TYPE_STR

c     Xmgr colors used in plotting waves:
      character*8 WV_COLOR(N_WAVE_TYPES)
      common /WV_COLOR_CM/ WV_COLOR

c     Xmgr linewidths used in plotting waves:
      character*11 WV_LINEWIDTH(N_WAVE_TYPES)
      common /WV_LINEWIDTH_CM/ WV_LINEWIDTH
