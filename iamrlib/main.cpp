#define OUTPUT_CONTOUR 1

#include <winstd.H>
#include <algorithm>
#include <vector>

#if defined(BL_OLD_STL)
#include <stdio.h>
#include <math.h>
#else
#include <cstdio>
#include <cmath>
#endif

#include <CArena.H>
#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>

extern void fortran_parameters();

int
main (int   argc,
      char* argv[])
{

    std::cout.imbue(std::locale("C"));
    BoxLib::Initialize(argc,argv);  // mpi initialization.
    ParallelDescriptor::Barrier();
    std::cout << "Multimaterial SPECTRAL, 8/24/19, 11:00 on proc " << 
	    ParallelDescriptor::MyProc() << "\n";
    ParallelDescriptor::Barrier();
    if (ParallelDescriptor::IOProcessor())
     std::cout << "after the barrier on IO processor " << 
	    ParallelDescriptor::MyProc() << "\n";

    const Real run_strt = ParallelDescriptor::second();

    int  wait_for_key=0;
    int  max_step;  // int is 4 bytes
    Real strt_time;
    Real stop_time;

    ParmParse pp;

    max_step  = -1;    
    strt_time =  0.0;  
    stop_time = -1.0;  

    pp.query("wait_for_key",wait_for_key);
    pp.query("max_step",max_step);
    pp.query("strt_time",strt_time);
    pp.query("stop_time",stop_time);

    if (strt_time < 0.0)
        BoxLib::Abort("MUST SPECIFY a non-negative strt_time");

    if (max_step < 0 && stop_time < 0.0)
    {
        BoxLib::Abort(
            "Exiting because neither max_step nor stop_time is non-negative.");
    }

     
    Amr* amrptr = new Amr();

    ParallelDescriptor::Barrier();

    fortran_parameters();

      // Amr::init 
    amrptr->init(strt_time,stop_time);

    ParallelDescriptor::Barrier();

      // if not subcycling then levelSteps(level) is independent of "level"
      // initially, cumTime()==0.0
    while ( amrptr->okToContinue()           &&
           (amrptr->levelSteps(0) < max_step || max_step < 0) &&
           (amrptr->cumTime() < stop_time || stop_time < 0.0) )
    {
        ParallelDescriptor::Barrier();

         // coarseTimeStep is in amrlib/Amr.cpp
        amrptr->coarseTimeStep(stop_time); // synchronizes internally

        ParallelDescriptor::Barrier();
    }
    ParallelDescriptor::Barrier();

    delete amrptr;

    if (CArena* arena = dynamic_cast<CArena*>(BoxLib::The_Arena()))
    {
        //
        // We're using a CArena -- output some FAB memory stats.
        // This'll output total # of bytes of heap space in the Arena.
        // It's actually the high water mark of heap space required by FABs.
        //
        char buf[256];

        sprintf(buf,
                "CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %ld",
                ParallelDescriptor::MyProc(),
                arena->heap_space_used());

        std::cout << buf << '\n';
    }

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Run time = " << run_stop << '\n';

    BoxLib::Finalize();

    return 0;
}
