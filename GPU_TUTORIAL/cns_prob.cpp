
#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>
#include "cns_prob_parm.H"
#include "CNS.H"

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {

        amrex_real test_float_size=1.0/10.0;
        std::cout << "sizeof(test_float_size)=" << 
         sizeof(test_float_size) << '\n';

        amrex::ParmParse pp("prob");

        pp.query("p_l", CNS::h_prob_parm->p_l);
        pp.query("p_r", CNS::h_prob_parm->p_r);
        pp.query("rho_l", CNS::h_prob_parm->rho_l);
        pp.query("rho_r", CNS::h_prob_parm->rho_r);
        pp.query("u_l", CNS::h_prob_parm->u_l);
        pp.query("u_r", CNS::h_prob_parm->u_r);

#ifdef AMREX_USE_CUDA
        amrex::Gpu::htod_memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#else
	std::memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#endif
    }
}
