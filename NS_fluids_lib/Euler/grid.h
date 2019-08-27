c---------------------------------------------------------------
c---------------------------------------------------------------
c                 ARRAY BOUNDARIES:
c---------------------------------------------------------------
c---------------------------------------------------------------
c                 Primary
c---------------------------------------------------------------
      integer    I_FIRST_CELL
      parameter( I_FIRST_CELL=1 )

      integer    I_LAST_CELL
      parameter( I_LAST_CELL=1000 )

      integer    N_GHOST_CELLS
      parameter( N_GHOST_CELLS=1 )

c--------------------------------------------------------------- 
c                 Derived
c---------------------------------------------------------------
      integer    I_FISRT_GHOST_CELL
      parameter( I_FISRT_GHOST_CELL=I_FIRST_CELL-N_GHOST_CELLS )

      integer    I_LAST_GHOST_CELL
      parameter( I_LAST_GHOST_CELL =I_LAST_CELL +N_GHOST_CELLS )

      integer    I_FISRT_GHOST_NODE
      parameter( I_FISRT_GHOST_NODE=I_FISRT_GHOST_CELL )

      integer    I_LAST_GHOST_NODE
      parameter( I_LAST_GHOST_NODE =I_LAST_GHOST_CELL+1 )
c---------------------------------------------------------------
c---------------------------------------------------------------
c                 CELL-NODE CONNECTIONS:
c---------------------------------------------------------------
c---------------------------------------------------------------
      integer leftCell(I_FISRT_GHOST_NODE:I_LAST_GHOST_NODE)
      integer rghtCell(I_FISRT_GHOST_NODE:I_LAST_GHOST_NODE)
      common /node_to_cell/ leftCell,rghtCell

      integer leftNode(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL)
      integer rghtNode(I_FISRT_GHOST_CELL:I_LAST_GHOST_CELL)
      common /cell_to_node/ leftNode,rghtNode

