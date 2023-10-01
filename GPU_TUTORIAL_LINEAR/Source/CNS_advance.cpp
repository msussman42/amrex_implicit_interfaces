
#include "CNS.H"
#include "CNS_hydro_K.H"
#include "CNS_K.H"

using namespace amrex;

Real
CNS::advance (Real time, Real dt, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("CNS::advance()");

    for (int i = 0; i < num_state_data_types; ++i) {
        state[i].allocOldData();
        state[i].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab dSdt(grids,dmap,NUM_STATE,0,MFInfo(),Factory());
    MultiFab Sborder(grids,dmap,NUM_STATE,NUM_GROW,MFInfo(),Factory());

    FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);
    compute_dSdt(Sborder, dSdt, dt);
    Real alpha=Real(1.0);
    Real beta=dt;
    int ngrow=0;
    int scomp=0;

    // U^{n+1} = alpha * U^n + beta*dUdt^n
    MultiFab::LinComb(S_new,alpha,Sborder,scomp,beta,dSdt,scomp,scomp, 
       NUM_STATE,ngrow);

    return dt;
}

void
CNS::compute_dSdt (const MultiFab& S, MultiFab& dSdt, Real dt)
{
    BL_PROFILE("CNS::compute_dSdt()");

    const auto dx = geom.CellSizeArray();
    const int ncomp = NUM_STATE;
    const int nchar = NUM_STATE;

    Parm const* lparm = d_parm;

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& sfab = S.array(mfi);
        auto const& dsdtfab = dSdt.array(mfi);

        const Box& bxg1 = amrex::grow(bx,1);

	FArrayBox wtmp,wnew;
        wtmp.resize(bxg1, nchar);
        wnew.resize(bx, nchar);

	 // "Eli" = "Extend Life" 
	 // (otherwise, host might delete the device data before the
	 // device is finished)
        Elixir weli = wtmp.elixir();
        auto const& w = wtmp.array();

        Elixir wneweli = wnew.elixir();
        auto const& wnew_d = wnew.array();

        amrex::ParallelFor(bxg1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_ctochar(i, j, k, sfab, w, *lparm);
        });

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_lax_wendroff(i, j, k, wnew_d, w, *lparm, dx);
        });

        amrex::ParallelFor(bx, 
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_chartoc(i, j, k, sfab, dsdtfab, wnew_d, dxinv);
        });

        // don't have to do this, but we could
        weli.clear(); // don't need them anymore
        wneweli.clear();
    }


}


