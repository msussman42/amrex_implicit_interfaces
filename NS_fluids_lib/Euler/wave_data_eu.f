      block data wave_data

      include 'wave_constants.h'

c--------------------------------------------------------------
c--------------------------------------------------------------
c
c      WAVE_TO_DIR:
c      maps wave index to its Direction. 
c      Direction being DIR_L=-1 for left-facing wave,
c                      DIR_M= 0 for center wave,  
c                      DIR_R=+1 for right-facing wave
c
      data WAVE_TO_DIR/ DIR_L, DIR_M, DIR_R /
c--------------------------------------------------------------
c--------------------------------------------------------------
c
c      WAVE_TO_EDGE_L:
c      maps wave index to the index of its LEFT edge
c
c      WAVE_TO_EDGE_R:
c      maps wave index to the index of its RIGHT edge
c
      data WAVE_TO_EDGE_L /1,3,5/
      data WAVE_TO_EDGE_R /2,4,6/
c--------------------------------------------------------------
c--------------------------------------------------------------
c
c      WAVE_TO_EDGE_A:
c      maps wave index to the index of its AHEAD edge
c
c      WAVE_TO_EDGE_B:
c      maps wave index to the index of its BEHIND edge
c
      data WAVE_TO_EDGE_A /1,0,6/
      data WAVE_TO_EDGE_B /2,0,5/
c NOTE: WAVE_TO_EDGE_A(2) and WAVE_TO_EDGE_B(2)
c--------------------------------------------------------------
c--------------------------------------------------------------
c
c      WAVE_TO_STATE_L:
c      maps wave index to the index of its LEFT state
c
c      WAVE_TO_STATE_R:
c      maps wave index to the index of its RIGHT state
c
      data WAVE_TO_STATE_L /1,2,3/
      data WAVE_TO_STATE_R /2,3,4/
c--------------------------------------------------------------
c--------------------------------------------------------------
c
c     WAVE_TO_STATE_A:
c     maps wave index to the index of the state AHEAD of the wave
c
c     WAVE_TO_STATE_B:
c     maps wave index to the index of the state BEHIND of the wave
c
c     assumes 0 if the AHEAD/BEHIND direction is not defined
c
      data WAVE_TO_STATE_A /1,0,4/
      data WAVE_TO_STATE_B /2,0,3/
c--------------------------------------------------------------
c--------------------------------------------------------------
c     Used Only in slope2state.f
c--------------------------------------------------------------
c--------------------------------------------------------------
c  
c     IS_L_EDGE:
c     assume .TRUE. when the edge is a LEFT edge of the wave 
c             (note that fake (N_EDGES+1)th edge is a LEFT edge)
c
      data IS_L_EDGE /.true., .false., .true., .false., 
     &     .true., .false., .true./
c--------------------------------------------------------------
c--------------------------------------------------------------
c      
c      EDGE_TO_L_STATE:
c      assumes 0 if the edge is a RIGHT edge
c      assumes (1...N_STATES) the index of the state LEFT to the edge,
c                if the edge is a LEFT edge 
c      (note: N_EDGES+1 elements)
c
      data EDGE_TO_L_STATE /1, 0, 2, 0, 3, 0, 4/
c--------------------------------------------------------------
c--------------------------------------------------------------
c      
c      EDGE_TO_L_WAVE:
c      assumes 0 if the edge is a LEFT edge
c      assumes (1...N_WAVES) the index of the wave LEFT to the edge,
c                if the edge is a RIGHT edge               
c      (note: N_EDGES+1 elements)
c
      data EDGE_TO_L_WAVE /0, 1, 0, 2, 0, 3, 0/ 
c--------------------------------------------------------------
c--------------------------------------------------------------
c     Wave Types:
c-------------------------------------------------------------------------
c-------------------------------------------------------------------------
c
c     WAVE_TYPE_STR ( wave type string ):
c
c     1  rarefaction fan
c     2  shock front
c
      data WAVE_TYPE_STR /'  RFn   ', '  SFr   ','   CD   '/
c--------------------------------------------------------------
c--------------------------------------------------------------
c
c     WV_COLOR:  Xmgr colors used in plotting.
c     
c
      data WV_COLOR /'color 4 ','color 1 ','color 2 '/
c--------------------------------------------------------------
c--------------------------------------------------------------
c
c     WV_LINEWIDTH:  Xmgr linewidths used in plotting.
c     
c
      data WV_LINEWIDTH /'linewidth 3','linewidth 5','linewidth 5'/
c
      end
