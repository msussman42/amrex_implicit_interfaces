#include <local_thread_class.H>
#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Box.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BLFort.H>
#include <AMReX_Utility.H>
#include <AMReX_Random.H>
#include <AMReX_Print.H>
#include <AMReX_Arena.H>
#include <AMReX_BLBackTrace.H>
#include <AMReX_MemPool.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_Machine.H>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <csignal>
#include <cfenv>
#include <iostream>
#include <iomanip>
#include <new>
#include <stack>
#include <limits>
#include <vector>
#include <algorithm>

int thread_class::nthreads;
double thread_class::number_mfiter_loops=0.0;
std::vector<double> thread_class::tile_d_numPts;
double thread_class::boxarray_d_numPts=0.0;

void 
thread_class::Initialize() {

 nthreads=1;

#ifdef _OPENMP
 int tid=0;
 int nthreads_local;
#endif

#ifdef _OPENMP
#pragma omp parallel private(nthreads_local, tid)
{

 /* Obtain thread number */
 tid = omp_get_thread_num();

 /* Only master thread does this */
 if (tid == 0) {
  nthreads_local = omp_get_num_threads();
  nthreads=nthreads_local;
 }
}  /* All threads join master thread and disband */
#endif

 if (nthreads<1)
  amrex::Error("nthreads invalid");

 tile_d_numPts.resize(nthreads);
 for (int tid_local=0;tid_local<nthreads;tid_local++)
  tile_d_numPts[tid_local]=0.0;
 number_mfiter_loops=0.0;

} // end subroutine thread_class::Initialize

void 
thread_class::reconcile_d_numPts(int caller_loop_id) {

 number_mfiter_loops=number_mfiter_loops+1.0;

 if (tile_d_numPts[0]==boxarray_d_numPts) {
  // do nothing
 } else {
  amrex::ParallelDescriptor::Barrier();
  std::fflush(NULL);
  amrex::ParallelDescriptor::Barrier();
  for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
   amrex::ParallelDescriptor::Barrier();
   if (amrex::ParallelDescriptor::MyProc()==pid) {
    std::fflush(NULL);
    std::cout << "on processor: " << pid << '\n';
    std::cout << "tile_d_numPts[0]= " << tile_d_numPts[0] << '\n';
    std::cout << "boxarray_d_numPts= " << boxarray_d_numPts << '\n';
    std::cout << "number_mfiter_loops= " << number_mfiter_loops << '\n';
    std::cout << "caller_loop_id= " << caller_loop_id << '\n';
    std::fflush(NULL);
   }
  }  // pid=0..NProcs-1
  amrex::ParallelDescriptor::Barrier();
  std::fflush(NULL);
  amrex::ParallelDescriptor::Barrier();
  amrex::Error("reconcile_d_numPts found tile sum <> boxarray sum\n");
 }

} // end subroutine thread_class::reconcile_d_numPts()

void
thread_class::init_d_numPts(double BA_d_numPts) {

 for (int tid_local=0;tid_local<nthreads;tid_local++) {
  tile_d_numPts[tid_local] = 0.0;
 }
 boxarray_d_numPts=BA_d_numPts;

} // end subroutine thread_class::init_d_numPts

void 
thread_class::sync_tile_d_numPts() {

 for (int tid_local=1;tid_local<nthreads;tid_local++) {
  tile_d_numPts[0]+=tile_d_numPts[tid_local];
 }

} // end subroutine thread_class::sync_tile_d_numPts()
