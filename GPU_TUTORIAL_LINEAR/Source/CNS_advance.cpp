#include "CNS.H"
#include "CNS_K.H"
#include "DSDT_F.H"

using namespace amrex;

Real
CNS::advance (Real time, Real dt, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("CNS::advance()");

    for (int i = 0; i < num_state_data_types; ++i) {
        state[i].allocOldData();
        state[i].swapTimeLevels(dt);
    }

    AMREX_ALWAYS_ASSERT(NUM_GROW==2);
    AMREX_ALWAYS_ASSERT(NUM_STATE==AMREX_SPACEDIM+1);

    MultiFab& S_new = get_new_data(State_Type);
    //MultiFab& S_old = get_old_data(State_Type);
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
CNS::compute_dSdt (const MultiFab& S, MultiFab& dSdt, Real dt) {
 BL_PROFILE("CNS::compute_dSdt()");

 //const int ncomp = NUM_STATE;
 //const int nchar = NUM_STATE;

   //AMReX_BLassert.H
 AMREX_ALWAYS_ASSERT(NUM_STATE==AMREX_SPACEDIM+1);
 AMREX_ALWAYS_ASSERT(dSdt.nGrow()==0);
 AMREX_ALWAYS_ASSERT(S.nGrow()==2);
 AMREX_ALWAYS_ASSERT(dt>Real(0.0));

 //Parm const* lparm = d_parm;
 ProbParm const* lprob_parm = d_prob_parm;

 for (MFIter mfi(S); mfi.isValid(); ++mfi) {
  const Box& bx = mfi.tilebox();

  AMREX_ALWAYS_ASSERT(bx.ixType()==IndexType::TheCellType());
  AMREX_ALWAYS_ASSERT(dSdt[mfi].box().ixType()==IndexType::TheCellType());

  if (1==0) {
   const Real* dx=geom.CellSize();
   Real c=lprob_parm->sound_speed;

   const FArrayBox& sfab=S[mfi];
   const FArrayBox& dsdtfab=dSdt[mfi];
   const int* lo=bx.loVect();
   const int* hi=bx.hiVect();
   fort_dsdtfab(
     dsdtfab.dataPtr(),
     ARLIM(dsdtfab.loVect()),
     ARLIM(dsdtfab.hiVect()),
     sfab.dataPtr(),
     ARLIM(sfab.loVect()),
     ARLIM(sfab.hiVect()),
     dx,
     &dt,
     &c,
     lo,hi);

  } else {

  const auto dx = geom.CellSizeArray();
  auto const& sfab = S.array(mfi);
  auto const& dsdtfab = dSdt.array(mfi);

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
   Real c=lprob_parm->sound_speed;
   Real c2=c*c;

   AMREX_ASSERT(c2>0.0);
   AMREX_ASSERT(c>0.0);

   Real h=dx[0];
   Real h2=h*h;
   AMREX_ASSERT(h>0.0);
   AMREX_ASSERT(h2>0.0);

   int comp=0;
   Real local_dsdt=Real(0.0);

   local_dsdt=
    -c2*(sfab(i+1,j,k,1)-sfab(i-1,j,k,1))/(Real(2.0)*h)+
    Real(0.5)*dt*(sfab(i+1,j,k,comp)-
     Real(2.0)*sfab(i,j,k,comp)+
     sfab(i-1,j,k,comp))/h2;
#if ((AMREX_SPACEDIM==2)||(AMREX_SPACEDIM==3))
   local_dsdt=local_dsdt-
    c2*(sfab(i,j+1,k,2)-sfab(i,j-1,k,2))/(Real(2.0)*h)+
    Real(0.5)*dt*(sfab(i,j+1,k,comp)-
     Real(2.0)*sfab(i,j,k,comp)+
     sfab(i,j-1,k,comp))/h2;
#if (AMREX_SPACEDIM==3) 	
   local_dsdt=local_dsdt-
    c2*(sfab(i,j,k+1,3)-sfab(i,j,k-1,3))/(Real(2.0)*h)+
    Real(0.5)*dt*(sfab(i,j,k+1,comp)-
     Real(2.0)*sfab(i,j,k,comp)+
     sfab(i,j,k-1,comp))/h2;
#endif
#endif
   dsdtfab(i,j,k,comp)=local_dsdt;

   AMREX_ASSERT(dsdtfab(i,j,k,comp)==dsdtfab(i,j,k,comp));

   comp=1;
   dsdtfab(i,j,k,comp)=
    -(sfab(i+1,j,k,0)-sfab(i-1,j,k,0))/(Real(2.0)*h)+
    Real(0.5)*dt*(sfab(i+1,j,k,comp)-
     Real(2.0)*sfab(i,j,k,comp)+
     sfab(i-1,j,k,comp))/h2;

   AMREX_ASSERT(dsdtfab(i,j,k,comp)==dsdtfab(i,j,k,comp));

#if ((AMREX_SPACEDIM==2)||(AMREX_SPACEDIM==3))
   comp=2;
   dsdtfab(i,j,k,comp)=
    -(sfab(i,j+1,k,0)-sfab(i,j-1,k,0))/(Real(2.0)*h)+
    Real(0.5)*dt*(sfab(i,j+1,k,comp)-
     Real(2.0)*sfab(i,j,k,comp)+
     sfab(i,j-1,k,comp))/h2;
   AMREX_ASSERT(dsdtfab(i,j,k,comp)==dsdtfab(i,j,k,comp));

#if (AMREX_SPACEDIM==3) 	
   comp=3;
   dsdtfab(i,j,k,comp)=
    -(sfab(i,j,k+1,0)-sfab(i,j,k-1,0))/(Real(2.0)*h)+
    Real(0.5)*dt*(sfab(i,j,k+1,comp)-
     Real(2.0)*sfab(i,j,k,comp)+
     sfab(i,j,k-1,comp))/h2;

   AMREX_ASSERT(dsdtfab(i,j,k,comp)==dsdtfab(i,j,k,comp));
#endif
#endif


  });

  }

 }


}

