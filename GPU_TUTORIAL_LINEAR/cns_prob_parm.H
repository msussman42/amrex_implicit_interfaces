#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    int num_state_variables=0;
    Vector< amrex::Real > eigenvalues;
    Vector< amrex::Real > Reigenvector1;
    Vector< amrex::Real > Reigenvector2;
    Vector< amrex::Real > Leigenvector1;
    Vector< amrex::Real > Leigenvector2;
};

#endif
