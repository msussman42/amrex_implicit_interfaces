
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

        pp.get("sound_speed", CNS::h_prob_parm->sound_speed);
	 //A P^{-1} = P^{-1} \Lambda
	 //A = P^{-1} \Lambda P
	 //right eigenvectors are columns of P^{-1}
	 //left eigenvectors are rows of P
	 //P A P^{-1} = P P^{-1} \Lambda = \Lambda
	 //P A = \Lambda P

#ifdef AMREX_USE_CUDA
        amrex::Gpu::htod_memcpy(CNS::d_prob_parm, CNS::h_prob_parm, 
		sizeof(ProbParm));
#else
        std::memcpy(CNS::d_prob_parm, CNS::h_prob_parm,sizeof(ProbParm));
#endif
    }
}
