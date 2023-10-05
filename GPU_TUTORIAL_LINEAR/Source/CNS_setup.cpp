
#include <CNS.H>

using namespace amrex;

int CNS::num_state_data_types = 0;
Parm* CNS::h_parm = nullptr;
Parm* CNS::d_parm = nullptr;
ProbParm* CNS::h_prob_parm = nullptr;
ProbParm* CNS::d_prob_parm = nullptr;

//static Box the_same_box (const Box& b) { return b; }
//static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall
//
/*
static int scalar_bc[] =
{
  BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, 
  BCType::reflect_even, BCType::reflect_even
};
*/

static int pres_bc[] =
{
  BCType::int_dir, BCType::reflect_even, BCType::ext_dir, BCType::reflect_even, 
  BCType::reflect_even, BCType::reflect_even
};

//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall
static int norm_vel_bc[] =
{
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_odd,  BCType::reflect_odd,  BCType::reflect_odd
};

static int tang_vel_bc[] =
{
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_odd
};

/*
static
void
set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,scalar_bc[lo_bc[i]]);
        bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}
*/

static
void
set_pres_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,pres_bc[lo_bc[i]]);
        bc.setHi(i,pres_bc[hi_bc[i]]);
    }
}


static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
#endif
}

void
CNS::variableSetUp ()
{
    h_parm = new Parm{}; // This is deleted in CNS::variableCleanUp().
    h_prob_parm = new ProbParm{};
    d_parm = (Parm*)The_Arena()->alloc(sizeof(Parm));
    d_prob_parm = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm));

    read_params();

    AMREX_ALWAYS_ASSERT(NUM_STATE==4);
    AMREX_ALWAYS_ASSERT(NUM_GROW==2);
   
    bool state_data_extrap = false;
    bool store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,NUM_GROW,NUM_STATE,
                           &cell_cons_interp,
			   state_data_extrap,store_in_checkpoint);

    Vector<BCRec>       bcs(NUM_STATE);
    Vector<std::string> name(NUM_STATE);
    BCRec bc;
    int cnt = 0;
    set_pres_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "P";
    cnt++; set_x_vel_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "U";
    cnt++; set_y_vel_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "V";
    cnt++; set_z_vel_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "W";

    StateDescriptor::BndryFunc bndryfunc(cns_bcfill);
    bndryfunc.setRunOnGPU(true);  // I promise the bc function will launch gpu kernels.

    int start_comp=0;
    desc_lst.setComponent(State_Type,
                          start_comp,
                          name,
                          bcs,
                          bndryfunc);

    num_state_data_types = desc_lst.size();

}

void
CNS::variableCleanUp ()
{
    delete h_parm;
    delete h_prob_parm;
    The_Arena()->free(d_parm);
    The_Arena()->free(d_prob_parm);
    desc_lst.clear();
    derive_lst.clear();
}
