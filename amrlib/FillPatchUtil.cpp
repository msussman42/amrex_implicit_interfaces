
#include <AMReX_Utility.H>
#include <FillPatchUtil.H>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex
{

    void FillPatchSingleLevel (
      int level,
      MultiFab& mf, 
      Real time, 
      const MultiFab& smf, 
      int scomp, 
      int dcomp, 
      int ncomp,
      const Geometry& geom, 
      PhysBCFunctBaseSUSSMAN& physbcf,
      Vector<int> scompBC_map,
      int bfact)
    {
     BL_PROFILE("FillPatchSingleLevel");

     BL_ASSERT(scomp+ncomp <= smf.nComp());
     BL_ASSERT(dcomp+ncomp <= mf.nComp());
     if (scompBC_map.size()!=ncomp)
      amrex::Error("scompBC_map has invalid size");

      // src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period
     mf.ParallelCopy(smf, scomp, dcomp, ncomp, IntVect{0}, 
        mf.nGrowVect(), geom.periodicity());

     if (1==0) {
      std::cout << "scomp=" << scomp << " dcomp= " << dcomp << " ncomp= " <<
       ncomp << '\n';
      for (int n=0;n<ncomp;n++) {
       Real smf_norm=smf.norm1(scomp+n,0);
       Real mf_norm=mf.norm1(dcomp+n,0);
       std::cout << "n= " << n << " BCmap(n)= " << scompBC_map[n] <<
        " smf_norm= " << smf_norm << " mf_norm= " << mf_norm << '\n';
      }
     }

     physbcf.FillBoundary(level,mf,time,dcomp,scompBC_map,ncomp,bfact);
    }


    void FillPatchTwoLevels (
     MultiFab& mf,   // target
     Real time,
     MultiFab& cmf,  // coarse
     MultiFab& fmf,  // fine
     int scomp, 
     int dcomp, 
     int ncomp,
     const Geometry& cgeom, 
     const Geometry& fgeom, 
     PhysBCFunctBaseSUSSMAN& cbc, 
     PhysBCFunctBaseSUSSMAN& fbc,
     Interpolater* mapper, 
     const Vector<BCRec>& global_bcs,
     Vector<int> scompBC_map,
     int levelc,int levelf,
     int bfactc,int bfactf) {

     BL_PROFILE("FillPatchTwoLevels");

     if ((levelc<0)||(levelc!=levelf-1))
      amrex::Error("levelc or levelf invalid");
     if (scompBC_map.size()!=ncomp)
      amrex::Error("scompBC_map has invalid size");
     if (global_bcs.size()<ncomp) {
      std::cout << "time= " << time << '\n';
      std::cout << "scomp,dcomp,ncomp= " << scomp << ' ' << dcomp <<
       ' ' << ncomp << '\n';
      std::cout << "levelc, levelf, bfactc, bfactf = " << levelc << ' ' << 
       levelf << ' ' << bfactc << ' ' << bfactf << '\n';
      std::cout << "global_bcs.size() " << global_bcs.size() << '\n';
      amrex::Error("global_bcs has invalid size");
     }

#ifdef AMREX_USE_EB
     amrex::Error("AMREX_USE_EB should never be defined");
#else
     EB2::IndexSpace const* index_space=nullptr;   
#endif

     IntVect ratio_vec(D_DECL(2,2,2));

     int ngrow = mf.nGrow();
     const IntVect& ngrow_vec=mf.nGrowVect();	   
     
     int do_the_interp=0;
     if (ngrow>0) {
      do_the_interp=1;
     } else if (ngrow==0) {
      // do nothing
     } else {
      amrex::Error("ngrow invalid");
     }
     if (do_the_interp==0) {
      if (mf.boxArray()!=fmf.boxArray()) {
       do_the_interp=1;
      }
      if (do_the_interp==0) {
       if (mf.DistributionMap()!=fmf.DistributionMap()) {
        do_the_interp=1;
       }
      } else if (do_the_interp==1) {
       // do nothing
      } else
       amrex::Error("do_the_interp invalid");
     } else if (do_the_interp==1) {
      // do nothing
     } else
      amrex::Error("do_the_interp invalid");

     if (do_the_interp==1) {

      const InterpolaterBoxCoarsener& coarsener = 
        mapper->BoxCoarsener(bfactc,bfactf);
	    
      Box fdomain = fgeom.Domain();
      fdomain.convert(mf.boxArray().ixType());

      Box fdomain_g(fdomain);
      for (int i = 0; i < AMREX_SPACEDIM; ++i) {
       if (fgeom.isPeriodic(i)) {
        fdomain_g.grow(i,ngrow);
       }
      } // i

        // find coarsen( mf intersect complement(fmf) within fdomain_g ).
      const FabArrayBase::FPinfo& fpc = 
       FabArrayBase::TheFPinfo(fmf, mf, fdomain_g, ngrow_vec, coarsener,
        amrex::coarsen(fgeom.Domain(),ratio_vec),index_space);

      bool empty_flag=fpc.ba_crse_patch.empty();

      if (empty_flag) {
	      // do nothing, no data on the coarse level is needed to
	      // fill the fine level target.
      } else if (! empty_flag) {

       MultiFab mf_crse_patch(fpc.ba_crse_patch,fpc.dm_crse_patch,ncomp,0,
         MFInfo().SetTag("mf_crse_patch"),FArrayBoxFactory());

        // This data will be interpolated next.       
       FillPatchSingleLevel(
        levelc,
        mf_crse_patch, // target data (coarse)
        time, 
        cmf,  // coarse data source
        scomp, 
        0,   // dst_comp
        ncomp, 
        cgeom, 
        cbc,
        scompBC_map,
        bfactc);

       Vector< BCRec > local_bcs;
       local_bcs.resize(ncomp);
       for (int i=0;i<ncomp;i++)
        local_bcs[i]=global_bcs[scompBC_map[i]]; 
	
       bool cc = fpc.ba_crse_patch.ixType().cellCentered();
       if ((cc!=true)&&(cc!=false))
        amrex::Error("cc bust");

#ifdef _OPENMP
#pragma omp parallel if (cc)
#endif
       for (MFIter mfi(mf_crse_patch); mfi.isValid(); ++mfi) {
        int li = mfi.LocalIndex();
        int gi = fpc.dst_idxs[li];		
        const Box& dbx = fpc.dst_boxes[li];
		    
        Vector<BCRec> bcr(ncomp);
        int src_comp_bcs=0;
        int dest_comp_bcr=0;
        amrex::setBC(dbx,fdomain,src_comp_bcs,dest_comp_bcr,ncomp,
          local_bcs,bcr);
		    
        mapper->interp(time,
                       mf_crse_patch[mfi], // source
	               0,
		       mf[gi], //dest; does not overwrite existing fine.
		       dcomp,
		       ncomp,
		       dbx,
   		       cgeom,
		       fgeom,
		       bcr,
                       levelc,levelf,bfactc,bfactf);
       } // mfi
      } else {
       amrex::Error("empty_flag invalid");
      }
     } else if (do_the_interp==0) {
      // do nothing, mf and fmf have same boxarray's and dmaps and
      // ngrow==0
     } else
      amrex::Error("do_the_interp invalid");

     FillPatchSingleLevel(
      levelf,
      mf,  // mf is now init with coarse data.
      time, 
      fmf, 
      scomp, 
      dcomp, 
      ncomp, 
      fgeom, 
      fbc,
      scompBC_map,
      bfactf);

    }  //FillPatchTwoLevels

}  // namespace amrex
