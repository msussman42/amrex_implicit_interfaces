
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

        pp.get("num_state_variables", CNS::h_prob_parm->num_state_variables);
	int local_nstate=CNS::h_prob_parm->num_state_variables;
	BL_ASSERT(local_nstate==2);
        pp.getarr("eigenvalues", CNS::h_prob_parm->eigenvalues,0,local_nstate);
	 //A P^{-1} = P^{-1} \Lambda
	 //A = P^{-1} \Lambda P
	 //right eigenvectors are columns of P^{-1}
	 //left eigenvectors are rows of P
	 //P A P^{-1} = P P^{-1} \Lambda = \Lambda
	 //P A = \Lambda P
        pp.query("Reigenvector1", CNS::h_prob_parm->Reigenvector1,
	  0,local_nstate);
        pp.query("Reigenvector2", CNS::h_prob_parm->Reigenvector2,
	  0,local_nstate);
        pp.query("Leigenvector1", CNS::h_prob_parm->Leigenvector1,
	  0,local_nstate);
        pp.query("Leigenvector2", CNS::h_prob_parm->Leigenvector2,
	  0,local_nstate);

        amrex::Gpu::htod_memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
    }
}
