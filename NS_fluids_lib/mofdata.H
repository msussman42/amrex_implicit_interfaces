      integer order_algorithm(1000)
      integer renormalize_order_algorithm(1000)
      integer rank_algorithm(1000)
      integer min_rank
      integer nonzero_ranks
      integer MOFITERMAX
      integer MOFITERMAX_AFTER_PREDICT

      common /mofcommon/ &
        order_algorithm, &
        renormalize_order_algorithm, &
        rank_algorithm, &
	min_rank, &
	nonzero_ranks, &
        MOFITERMAX, &
	MOFITERMAX_AFTER_PREDICT

