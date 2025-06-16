#include <local_thread_class.H>
#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#include <AMReX_AmrCore.H>

extern "C" void cpp_reduce_real_sum(int n,double sync_data[]);
extern "C" void main_cpp_keyboard();

using namespace amrex;

namespace amrex{

extern void fortran_parameters();
extern void fortran_deallocate_parameters();

}

void main_cpp_keyboard() {

  std::cout << 
   "press any number then enter (called from SOLIDFLUID_F90_KEYBOARD) \n";
  int n_input;
  std::cin >> n_input;

} //end subroutine main_cpp_keyboard

void cpp_reduce_real_sum(int n,double sync_data[]) {

 int debug_reduce=0;

 if (debug_reduce==1) {
  std::fflush(NULL);
  amrex::ParallelDescriptor::Barrier();

  for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
   amrex::ParallelDescriptor::Barrier();
   if (amrex::ParallelDescriptor::MyProc()==pid) {
    std::fflush(NULL);
    for (int i=0;i<n;i++) {
     std::cout << "BEFORE proc= " << pid << " i= " << i << " data= " <<
      sync_data[i] << '\n';
    }
    std::fflush(NULL);
   }
  }  // pid=0..NProcs-1
 }


 amrex::ParallelDescriptor::Barrier();
 amrex::ParallelDescriptor::ReduceRealSum(sync_data,n);
 amrex::ParallelDescriptor::Barrier();

 if (debug_reduce==1) {
  std::fflush(NULL);
  amrex::ParallelDescriptor::Barrier();

  for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
   amrex::ParallelDescriptor::Barrier();
   if (amrex::ParallelDescriptor::MyProc()==pid) {
    std::fflush(NULL);
    for (int i=0;i<n;i++) {
     std::cout << "AFTER proc= " << pid << " i= " << i << " data= " <<
      sync_data[i] << '\n';
    }
    std::fflush(NULL);
   }
  }  // pid=0..NProcs-1
 }

} // end subroutine cpp_reduce_real_sum

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

 int  max_step;  // int is 4 bytes
 amrex::Real strt_time;
 amrex::Real stop_time;

 ParmParse pp;
 ParmParse ppamr("amr");
 ParmParse ppns("ns");

  // LSA = Linear Stability Analysis
  // ABEL OKOJUNU
  // LSA_nsteps_power_method=number of power method iterations.
 int local_LSA_nsteps_power_method=0;
 ppamr.queryAdd("LSA_nsteps_power_method",local_LSA_nsteps_power_method);
 if (local_LSA_nsteps_power_method>=0) {
  //do nothing
 } else
  amrex::Error("expecting local_LSA_nsteps_power_method>=0");


 max_step  = -1;    
 strt_time =  0.0;  
 stop_time = -1.0;  

 double sleepsec=0.0;

 pp.queryAdd("max_step",max_step);

 Real local_fixed_dt=0.0;

 if (local_LSA_nsteps_power_method==0) {
  //do nothing
 } else if (local_LSA_nsteps_power_method>0) {

   //ABEL OKOJUNO
   //for LSA, max_step>=1, and 
   //LSA_steps=(max_step-initial_levelSteps)>=1
   //(initial_levelSteps=amrptr->levelSteps(0))
   //probably best to have 
   //LSA_steps=(max_step-initial_levelSteps)>1 so that
   //the first step is the "forcing step" and the next step(s) is/are the
   //resultant "perturbed" state.
  if (max_step>=1) {
   //do nothing
  } else
   amrex::Error("expecting 1<=max_step");

   //time step must be fixed if LSA.
  ppns.queryAdd("fixed_dt",local_fixed_dt);
  if (local_fixed_dt>0.0) {
   //do nothing
  } else
   amrex::Error("expecting ns.fixed_dt>0.0 if LSA");

 } else
  amrex::Error("expecting local_LSA_nsteps_power_method>=0");

 pp.queryAdd("strt_time",strt_time);
 pp.queryAdd("stop_time",stop_time);

 pp.queryAdd("sleepsec",sleepsec);

 if (strt_time < 0.0)
  amrex::Abort("MUST SPECIFY a non-negative strt_time");

 if (max_step < 0 && stop_time < 0.0) {
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
   //  a) if !restart_file.empty() and restart_file!="init" then:
   //     "restart(restsrt_file)"
   //  b) otherwise: "initialInit(strt_time,stop_time)"
   //
 amrptr->init(strt_time,stop_time);

 amrex::ParallelDescriptor::Barrier();

 Real initial_cumTime=amrptr->cumTime();
 int initial_levelSteps=amrptr->levelSteps(0);
 int LSA_steps=max_step-initial_levelSteps;
 Real compare_time_scale=stop_time-initial_cumTime;

 if (amrex::ParallelDescriptor::IOProcessor()) {
  std::cout << "stop_time= " << stop_time <<'\n';
  std::cout << "initial_cumTime= " << initial_cumTime <<'\n';
  std::cout << "compare_time_scale= " << compare_time_scale <<'\n';
  std::cout << "initial_levelSteps= " << initial_levelSteps <<'\n';
  std::cout << "max_step= " << max_step <<'\n';
  std::cout << "LSA_steps=max_step-initial_levelSteps= " << LSA_steps <<'\n';
 }

 if (local_LSA_nsteps_power_method==0) {
  //do nothing
 } else if (local_LSA_nsteps_power_method>0) {

  if (initial_cumTime>0.0) {
   //do nothing
  } else
   amrex::Error("LSA: expecting initial_cumTime>0.0");
 
//ABEL OKOJUNO
//LSA_steps=max_step-initial_levelSteps
  if ((LSA_steps>0)&&(LSA_steps<9999)) {
   //do nothing
  } else
   amrex::Error("LSA: expecting 0<LSA_steps<9999");

  Real time_scale=local_fixed_dt*LSA_steps;
//ABEL OKOJUNO: compare_time_scale=stop_time-initial_cumTime;
  if (std::abs(compare_time_scale-time_scale)<=1.0e-4*time_scale) {
   //do nothing
  } else
   amrex::Error("LSA: need |compare_time_scale-time_scale|<eps");

 } else
  amrex::Error("expecting local_LSA_nsteps_power_method>=0");

//ABEL OKOJUNO:
//1. read in the "steady data" at time t0: data(t0)
//2. with no perturbations run to time t1: data(t1)
//3. let x^{(0)} initialized
//4. k=0
//5. using data(t0)+eps * x^{(k)} run to time t1: data(t1)^{(k)}
//6. let x^{(k+1)}=normalized(data(t1)^{(k)}-data(t1))
//7. k=k+1
//8. go back to step 5.
 for (int LSA_current_step=0;
      LSA_current_step<=local_LSA_nsteps_power_method;
      LSA_current_step++) {

  if (LSA_current_step==0) {
   //do nothing
  } else if ((LSA_current_step>=1)&&
             (LSA_current_step<=local_LSA_nsteps_power_method)) {
   amrex::ParallelDescriptor::Barrier();
   amrptr->rewindTimeStep(stop_time,LSA_current_step,
    initial_cumTime,initial_levelSteps);
   amrex::ParallelDescriptor::Barrier();
  } else
   amrex::Error("LSA_current_step invalid");

   // if not subcycling then levelSteps(level) is independent of "level"
   // initially, cumTime()==0.0
   // if LSA_current_step>=1, then reset:
   // cumTime,levelSteps,initial state
  while ( amrptr->okToContinue()           &&
         (amrptr->levelSteps(0) < max_step || max_step < 0) &&
         (amrptr->cumTime() < stop_time || stop_time < 0.0) ) {
   amrex::ParallelDescriptor::Barrier();
   std::fflush(NULL);
   BL_PROFILE_INITIALIZE();
   std::fflush(NULL);

   // coarseTimeStep is in amrlib/AMReX_AmrCore.cpp
   // timeStep is in amrlib/AMReX_AmrCore.cpp
   amrptr->coarseTimeStep(stop_time,
     LSA_current_step,initial_levelSteps);//synchronizes internally

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

 }  //LSA_current_step=0 .... local_LSA_nsteps_power_method

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
	"Multimaterial SUPERMESH/SPECTRAL, June 16, 2025, 18:01pm on proc " << 
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

     if ((sizeof(amrex::Real)==sizeof(double))||
         (sizeof(amrex::Real)==sizeof(float))) {
      //do nothing
     } else {
      amrex::Error("amrex::Real invalid size");
     }
     int int_size=sizeof(int);
     if ((int_size==4)||(int_size==8)) {
      // do nothing
     } else {
      std::cout << "int_size= " << int_size << '\n';
      amrex::Error("expecting int_size == 4 or 8");
     }
     std::cout << "real size= " << sizeof(amrex::Real) << '\n';
     std::cout << "int size= " << int_size << '\n';
    } //IOProcessor==TRUE

    int fork_id=0;
    fork_job(fork_id);

     //AmrCore::regrid_ba.clear();
     //AmrCore::initial_ba.clear();
    amrex::Finalize();

    return 0;

}

