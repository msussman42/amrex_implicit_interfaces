#include <local_thread_class.H>
#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#include <AMReX_AmrCore.H>

using namespace amrex;

namespace amrex{

extern void fortran_parameters();
extern void fortran_deallocate_parameters();

}

void
fork_job(int fork_id) {


 std::fflush(NULL);
 amrex::ParallelDescriptor::Barrier();

 for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
  amrex::ParallelDescriptor::Barrier();
  if (amrex::ParallelDescriptor::MyProc()==pid) {
   std::fflush(NULL);
   std::cout << 
     "fork_id= " << fork_id << " proc = " << 
     amrex::ParallelDescriptor::MyProc() << "\n";
   std::cout << 
     "fork_id= " << fork_id << " NProcs()= " << 
     amrex::ParallelDescriptor::NProcs() << '\n';
   std::cout << 
     "fork_id= " << fork_id << " PROC= " << 
         amrex::ParallelDescriptor::MyProc() << 
         " thread_class::nthreads= " << 
         thread_class::nthreads << '\n';
   std::fflush(NULL);
  }
 }  // pid=0..NProcs-1

 amrex::ParallelDescriptor::Barrier();

 const double run_strt = amrex::ParallelDescriptor::second();

 int  wait_for_key=0;
 int  max_step;  // int is 4 bytes
 amrex::Real strt_time;
 amrex::Real stop_time;

 ParmParse pp;

 max_step  = -1;    
 strt_time =  0.0;  
 stop_time = -1.0;  

 double sleepsec=0.0;

 pp.queryAdd("wait_for_key",wait_for_key);
 pp.queryAdd("max_step",max_step);
 pp.queryAdd("strt_time",strt_time);
 pp.queryAdd("stop_time",stop_time);

 pp.queryAdd("sleepsec",sleepsec);

 if (strt_time < 0.0)
  amrex::Abort("MUST SPECIFY a non-negative strt_time");

 if (max_step < 0 && stop_time < 0.0)
 {
     amrex::Abort(
         "Exiting because neither max_step nor stop_time is non-negative.");
 }

  // NavierStokes.cpp (fortran_parameters) ->
  // PROB_CPP_PARMS.F90 (fort_override) ->
  // PROB_CPP_PARMS.F90 (SUB_INIT_MODULE) 
 fortran_parameters();

   // 1. call AmrMesh()
   //   a. AmrMesh() calls Geometry::Setup() 
   //     (rb==nullptr, coord=-1, is_per=nullptr)
   //   b. AmrMesh() calls InitAmrMesh
   //     i. geom[i].define(index_domain)  i=0..max_level
   //                      
   // 2. Initialize()
   // 3. InitAmr() 
   //   a.  levelbld = getLevelBld();
 AmrCore* amrptr = new AmrCore();

 amrex::ParallelDescriptor::Barrier();

   // AmrCore::init 
 amrptr->init(strt_time,stop_time);

 amrex::ParallelDescriptor::Barrier();

   // if not subcycling then levelSteps(level) is independent of "level"
   // initially, cumTime()==0.0
 while ( amrptr->okToContinue()           &&
        (amrptr->levelSteps(0) < max_step || max_step < 0) &&
        (amrptr->cumTime() < stop_time || stop_time < 0.0) )
 {
  amrex::ParallelDescriptor::Barrier();
  std::fflush(NULL);
  BL_PROFILE_INITIALIZE();
  std::fflush(NULL);

   // coarseTimeStep is in amrlib/AMReX_AmrCore.cpp
  amrptr->coarseTimeStep(stop_time); // synchronizes internally

  amrex::ParallelDescriptor::Barrier();
  std::fflush(NULL);
  BL_PROFILE_FINALIZE();
  std::fflush(NULL);
  std::cout << "TIME= " << amrptr->cumTime() << " PROC= " <<
    amrex::ParallelDescriptor::MyProc() << " sleepsec= " << sleepsec << '\n';
  std::fflush(NULL);
  amrex::Sleep(sleepsec);
  amrex::ParallelDescriptor::Barrier();
 }
 amrex::ParallelDescriptor::Barrier();

 delete amrptr;

 if (CArena* arena = dynamic_cast<CArena*>(amrex::The_Arena()))
 {
  //
  // We're using a CArena -- output some FAB memory stats.
  // This'll output total # of bytes of heap space in the Arena.
  // It's actually the high water mark of heap space required by FABs.
  //
  char buf[256];

  sprintf(buf,
          "CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %ld",
          amrex::ParallelDescriptor::MyProc(),
          arena->heap_space_used());

  std::cout << buf << '\n';
 }

 const int IOProc   = amrex::ParallelDescriptor::IOProcessorNumber();
 double    run_stop = amrex::ParallelDescriptor::second() - run_strt;

 amrex::ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

 if (amrex::ParallelDescriptor::IOProcessor())
     std::cout << "Run time = " << run_stop << '\n';

  // NavierStokes.cpp (fortran_deallocate_parameters) ->
  // PROB_CPP_PARMS.F90 (fort_deallocate_module) ->
  // PROB_CPP_PARMS.F90 (SUB_DEALLOCATE_MODULE) 
 fortran_deallocate_parameters();

}

int
main (int   argc,
      char* argv[])
{

    std::cout.imbue(std::locale("C"));

    amrex::Initialize(argc,argv);  

    thread_class::Initialize();

    std::fflush(NULL);
    amrex::ParallelDescriptor::Barrier();

    for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
     amrex::ParallelDescriptor::Barrier();
     if (amrex::ParallelDescriptor::MyProc()==pid) {
      std::fflush(NULL);
      std::cout << 
	"Multimaterial SUPERMESH/SPECTRAL, 06/12/23, 19:10 on proc " << 
        amrex::ParallelDescriptor::MyProc() << "\n";
      std::cout << "NProcs()= " << 
        amrex::ParallelDescriptor::NProcs() << '\n';
      std::cout << "PROC= " << amrex::ParallelDescriptor::MyProc() << 
	    " thread_class::nthreads= " << 
	    thread_class::nthreads << '\n';
      std::fflush(NULL);
     }
    }  // pid=0..NProcs-1

    amrex::ParallelDescriptor::Barrier();

    if (amrex::ParallelDescriptor::IOProcessor()) {
     std::cout << "after the barrier on IO processor " << 
	    amrex::ParallelDescriptor::MyProc() << "\n";
     int double_size=sizeof(double);
     if (sizeof(amrex::Real)!=double_size) 
      amrex::Error("expecting amrex::Real and double_size to be equal");
     if (double_size!=8)
      amrex::Error("expecting double_size == 8 ");
     int int_size=sizeof(int);
     if ((int_size==4)||(int_size==8)) {
      // do nothing
     } else {
      std::cout << "int_size= " << int_size << '\n';
      amrex::Error("expecting int_size == 4 or 8");
     }
     std::cout << "double_size= " << double_size << '\n';
     std::cout << "int_size= " << int_size << '\n';
    } //IOProcessor==TRUE

    int fork_id=0;
    fork_job(fork_id);

    amrex::Finalize();

    return 0;

}

