
#include <CNS.H>
#include <CNS_K.H>
#include <CNS_tagging.H>
#include <CNS_parm.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <climits>

using namespace amrex;

constexpr int CNS::NUM_GROW;
constexpr int CNS::State_Type;

int CNS::NUM_STATE=0;

BCRec     CNS::phys_bc;

int       CNS::verbose = 0;
IntVect   CNS::hydro_tile_size {AMREX_D_DECL(1024,16,16)};
Real      CNS::cfl       = 0.3;

CNS::CNS ()
{}

CNS::CNS (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          const DistributionMapping& dm,
          Real            time)
    : AmrLevel(papa,lev,level_geom,bl,dm,time)
{

    buildMetrics();
}

CNS::~CNS ()
{}

void
CNS::init (AmrLevel& old)
{
    auto& oldlev = dynamic_cast<CNS&>(old);

    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev.state[State_Type].curTime();
    Real prev_time = oldlev.state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    FillPatch(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
}

void
CNS::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();
    Real dt_old = (cur_time - prev_time)/static_cast<Real>(parent->MaxRefRatio(level-1));
    setTimeLevel(cur_time,dt_old,dt);

    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);
}

void
CNS::initData ()
{
    BL_PROFILE("CNS::initData()");

    const auto geomdata = geom.data();
    const auto dx=geom.CellSizeArray();
    const auto prob_lo=geom.ProbLoArray();

    MultiFab& S_new = get_new_data(State_Type);

    Parm const* lparm = d_parm;
    ProbParm const* lprobparm = d_prob_parm;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        auto sfab = S_new.array(mfi);

        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_initdata(i, j, k, 
			 sfab, 
			 geomdata,
			 dx,prob_lo, 
			 *lparm, 
			 *lprobparm);
        });
    }
}

void
CNS::computeInitialDt (int                    finest_level,
                       int                    /*sub_cycle*/,
                       Vector<int>&           n_cycle,
                       const Vector<IntVect>& /*ref_ratio*/,
                       Vector<Real>&          dt_level,
                       Real                   stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0) {
        return;
    }

    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
CNS::computeNewDt (int                    finest_level,
                   int                    /*sub_cycle*/,
                   Vector<int>&           n_cycle,
                   const Vector<IntVect>& /*ref_ratio*/,
                   Vector<Real>&          dt_min,
                   Vector<Real>&          dt_level,
                   Real                   stop_time,
                   int                    post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0) {
        return;
    }

    for (int i = 0; i <= finest_level; i++)
    {
        dt_min[i] = getLevel(i).estTimeStep();
    }

    if (post_regrid_flag == 1)
    {
        //
        // Limit dt's by pre-regrid dt
        //
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],dt_level[i]);
        }
    }
    else
    {
        //
        // Limit dt's by change_max * old dt
        //
        static Real change_max = 1.1;
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
        }
    }

    //
    // Find the minimum over all levels
    //
    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
CNS::post_regrid (int /*lbase*/, int /*new_finest*/)
{
}

void
CNS::post_timestep (int /*iteration*/)
{
    BL_PROFILE("post_timestep");

    if (level < parent->finestLevel()) {
        avgDown();
    }
}

void
CNS::postCoarseTimeStep (Real time)
{
    BL_PROFILE("postCoarseTimeStep()");

    // This only computes sum on level 0
    if (verbose >= 2) {
        printTotal(time);
    }
}


void
CNS::compute_errors (MultiFab& Error_Analysis,const Real time) {
 BL_PROFILE("CNS::compute_errors()");

 const auto dx = geom.CellSizeArray();
 const auto prob_lo=geom.ProbLoArray();

   //AMReX_BLassert.H
 AMREX_ALWAYS_ASSERT(NUM_STATE==AMREX_SPACEDIM+1);
 AMREX_ALWAYS_ASSERT(Error_Analysis.nGrow()==0);

 //Parm const* lparm = d_parm;
 ProbParm const* lprob_parm = d_prob_parm;

 MultiFab& S_new = get_new_data(State_Type);

 for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
  const Box& bx = mfi.tilebox();

  AMREX_ALWAYS_ASSERT(bx.ixType()==IndexType::TheCellType());
  AMREX_ALWAYS_ASSERT(Error_Analysis[mfi].box().ixType()==IndexType::TheCellType());

  auto snewfab = S_new.array(mfi);
  auto errorfab = Error_Analysis.array(mfi);

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
   cns_compute_errors(i,j,k,
                      snewfab,
                      errorfab,
                      dx,prob_lo,
                      time);
  });

 }


}

void
CNS::printTotal (amrex::Real time) 
{
    const MultiFab& S_new = get_new_data(State_Type);
    int ncomp=S_new.nComp();
    AMREX_ALWAYS_ASSERT(ncomp==AMREX_SPACEDIM+1);

    MultiFab Error_Analysis(grids,dmap,2,0,MFInfo(),Factory());
    compute_errors(Error_Analysis,time);

    std::array<Real,AMREX_SPACEDIM+1> tot;
    std::array<Real,2> err_tot;

    for (int comp = 0; comp < ncomp; ++comp) {
        tot[comp] = S_new.sum(comp,true) * geom.ProbSize();
    }
    for (int comp = 0; comp < 2; ++comp) {
        err_tot[comp] = Error_Analysis.sum(comp,true);
    }

#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelDescriptor::ReduceRealSum(tot.data(), ncomp, ParallelDescriptor::IOProcessorNumber());
#if (AMREX_SPACEDIM==1)
            amrex::Print().SetPrecision(17) << "\n[CNS] Total p      is " << tot[0] << "\n"
                                            <<   "      Total u      is " << tot[1] << "\n";
#endif
#if (AMREX_SPACEDIM==2)
            amrex::Print().SetPrecision(17) << "\n[CNS] Total p      is " << tot[0] << "\n"
                                            <<   "      Total u      is " << tot[1] << "\n"
                                            <<   "      Total v      is " << tot[2] << "\n";
#endif
#if (AMREX_SPACEDIM==3)
            amrex::Print().SetPrecision(17) << "\n[CNS] Total p      is " << tot[0] << "\n"
                                            <<   "      Total u      is " << tot[1] << "\n"
                                            <<   "      Total v      is " << tot[2] << "\n"
                                            <<   "      Total w      is " << tot[3] << "\n";
#endif
#ifdef BL_LAZY
        });
#endif


#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelDescriptor::ReduceRealSum(err_tot.data(), 2, ParallelDescriptor::IOProcessorNumber());
            amrex::Print().SetPrecision(17) << "\n[CNS] Total vol    is " << err_tot[0] << "\n"
                                            <<   "      Total err    is " << err_tot[1]/err_tot[0] << "\n";
#ifdef BL_LAZY
        });
#endif



}

void
CNS::post_init (Real /*stop_time*/)
{
    if (level > 0) return;
    for (int k = parent->finestLevel()-1; k >= 0; --k) {
        getLevel(k).avgDown();
    }

    if (verbose >= 2) {
        amrex::Real time=0.0;
        printTotal(time);
    }
}

void
CNS::post_restart ()
{
}

void
CNS::errorEst (TagBoxArray& tags, int, int, Real /*time*/, int, int)
{
    BL_PROFILE("CNS::errorEst()");

    const MultiFab& S_new = get_new_data(State_Type);
    const Real cur_time = state[State_Type].curTime();
    MultiFab p(S_new.boxArray(), S_new.DistributionMap(), 1, 1);
    FillPatch(*this, p, p.nGrow(), cur_time, State_Type, 0, 1, 0);

    const char   tagval = TagBox::SET;
//  const char clearval = TagBox::CLEAR;
    const Real p_threshold = 1.0e+20;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(p,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
     const Box& bx = mfi.tilebox();

     const auto pfab = p.array(mfi);
     auto tag = tags.array(mfi);

     amrex::ParallelFor(bx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
        cns_tag_perror(i, j, k, tag, pfab, p_threshold, tagval);
       });
    }
}

void
CNS::read_params ()
{
    ParmParse pp("cns");

    pp.get("NUM_STATE",NUM_STATE);

    pp.query("v", verbose);

    Vector<int> tilesize(AMREX_SPACEDIM);
    if (pp.queryarr("hydro_tile_size", tilesize, 0, AMREX_SPACEDIM))
    {
        for (int i=0; i<AMREX_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
    }

    pp.query("cfl", cfl);

    Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
    pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    h_parm->Initialize();
#ifdef AMREX_USE_CUDA
    amrex::Gpu::htod_memcpy(d_parm, h_parm, sizeof(Parm));
#else
    std::memcpy(d_parm, h_parm, sizeof(Parm));
#endif

}

void
CNS::avgDown ()
{
    BL_PROFILE("CNS::avgDown()");

    if (level == parent->finestLevel()) return;

    auto& fine_lev = getLevel(level+1);

    MultiFab& S_crse =          get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    amrex::average_down(S_fine, S_crse, fine_lev.geom, geom,
                        0, S_fine.nComp(), parent->refRatio(level));

}

void
CNS::buildMetrics ()
{
    // make sure dx == dy == dz
    const Real* dx = geom.CellSize();
#if (AMREX_SPACEDIM==3)
    if (std::abs(dx[0]-dx[1]) > Real(1.e-12)*dx[0] || std::abs(dx[0]-dx[2]) > Real(1.e-12)*dx[0]) {
        amrex::Abort("CNS: must have dx == dy == dz\n");
    }
#endif
#if (AMREX_SPACEDIM==2)
    if (std::abs(dx[0]-dx[1]) > Real(1.e-12)*dx[0]) {
        amrex::Abort("CNS: must have dx == dy\n");
    }
#endif
}

Real
CNS::estTimeStep ()
{
    BL_PROFILE("CNS::estTimeStep()");

    const auto dx = geom.CellSizeArray();
    const MultiFab& S = get_new_data(State_Type);
    Parm const* lparm = d_parm;
    ProbParm const* lprob_parm = d_prob_parm;

    Real estdt = amrex::ReduceMin(S, 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& fab) -> Real
    {
        return cns_estdt(bx, fab, dx, *lparm,*lprob_parm);
    });

    estdt *= cfl;
    ParallelDescriptor::ReduceRealMin(estdt);

    return estdt;
}

Real
CNS::initialTimeStep ()
{
    return estTimeStep();
}

