
#ifndef LOCAL_THREAD_CLASS_H
#define LOCAL_THREAD_CLASS_H

#include <iostream>
#include <functional>
#include <string>
#include <cstdio>
#include <vector>
#include <memory>


// SUSSMAN
class thread_class {

 public:
  static void Initialize();
  static void init_d_numPts(double BA_d_numPts);
  static void sync_tile_d_numPts();
  static int nthreads;
  static double number_mfiter_loops;
  static std::vector<double> tile_d_numPts;
  static double boxarray_d_numPts;
  static void reconcile_d_numPts(int caller_loop_id,
    const std::string& caller_string);
};

#endif /*LOCAL_THREAD_CLASS_H*/
