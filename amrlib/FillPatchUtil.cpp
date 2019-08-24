
#include <FillPatchUtil.H>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace BoxLib
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
      PhysBCFunctBase& physbcf,
      Array<int> scompBC_map,
      int bfact)
    {
     BL_PROFILE("FillPatchSingleLevel");

     BL_ASSERT(scomp+ncomp <= smf.nComp());
     BL_ASSERT(dcomp+ncomp <= mf.nComp());
     if (scompBC_map.size()!=ncomp)
      BoxLib::Error("scompBC_map has invalid size");

      // src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period
     mf.copy(smf, scomp, dcomp, ncomp, 0, mf.nGrow(), geom.periodicity());

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
     MultiFab& mf, 
     Real time,
     MultiFab& cmf, 
     MultiFab& fmf, 
     int scomp, 
     int dcomp, 
     int ncomp,
     const Geometry& cgeom, 
     const Geometry& fgeom, 
     PhysBCFunctBase& cbc, 
     PhysBCFunctBase& fbc,
     Interpolater* mapper, 
     const Array<BCRec>& global_bcs,
     Array<int> scompBC_map,
     int levelc,int levelf,
     int bfactc,int bfactf) {

     BL_PROFILE("FillPatchTwoLevels");

     if ((levelc<0)||(levelc!=levelf-1))
      BoxLib::Error("levelc or levelf invalid");
     if (scompBC_map.size()!=ncomp)
      BoxLib::Error("scompBC_map has invalid size");
     if (global_bcs.size()<ncomp) {
      std::cout << "time= " << time << '\n';
      std::cout << "scomp,dcomp,ncomp= " << scomp << ' ' << dcomp <<
       ' ' << ncomp << '\n';
      std::cout << "levelc, levelf, bfactc, bfactf = " << levelc << ' ' << 
       levelf << ' ' << bfactc << ' ' << bfactf << '\n';
      std::cout << "global_bcs.size() " << global_bcs.size() << '\n';
      BoxLib::Error("global_bcs has invalid size");
     }

     int ngrow = mf.nGrow();
	    
     if ((ngrow > 0)||(mf.getBDKey() != fmf.getBDKey())) {

      const InterpolaterBoxCoarsener& coarsener = 
        mapper->BoxCoarsener(bfactc,bfactf);
	    
      Box fdomain = fgeom.Domain();
      fdomain.convert(mf.boxArray().ixType());

      Box fdomain_g(fdomain);
      for (int i = 0; i < BL_SPACEDIM; ++i) {
       if (fgeom.isPeriodic(i)) {
        fdomain_g.grow(i,ngrow);
       }
      } // i

        // find coarsen( mf intersect complement(fmf) within fdomain_g ).
      const FabArrayBase::FPinfo& fpc = 
       FabArrayBase::TheFPinfo(fmf, mf, fdomain_g, ngrow, coarsener);

      if ( ! fpc.ba_crse_patch.empty()) {
       MultiFab mf_crse_patch(fpc.ba_crse_patch, ncomp, 0, 
        fpc.dm_crse_patch);
		
       FillPatchSingleLevel(
        levelc,
        mf_crse_patch, 
        time, 
        cmf, 
        scomp, 
        0,   // dst_comp
        ncomp, 
        cgeom, 
        cbc,
        scompBC_map,
        bfactc);

       Array< BCRec > local_bcs;
       local_bcs.resize(ncomp);
       for (int i=0;i<ncomp;i++)
        local_bcs[i]=global_bcs[scompBC_map[i]]; 
	
       bool cc = fpc.ba_crse_patch.ixType().cellCentered();
       if ((cc!=true)&&(cc!=false))
        BoxLib::Error("cc bust");

#ifdef _OPENMP
#pragma omp parallel if (cc)
#endif
       for (MFIter mfi(mf_crse_patch); mfi.isValid(); ++mfi) {
        int li = mfi.LocalIndex();
        int gi = fpc.dst_idxs[li];		
        const Box& dbx = fpc.dst_boxes[li];
		    
        Array<BCRec> bcr(ncomp);
        int src_comp_bcs=0;
        int dest_comp_bcr=0;
        BoxLib::setBC(dbx,fdomain,src_comp_bcs,dest_comp_bcr,ncomp,
          local_bcs,bcr);
		    
        mapper->interp(time,
                       mf_crse_patch[mfi],
	               0,
		       mf[gi],
		       dcomp,
		       ncomp,
		       dbx,
   		       cgeom,
		       fgeom,
		       bcr,
                       levelc,levelf,bfactc,bfactf);
       } // mfi
      } // some boxes need to be interpolated
     }  // coarse data needed

     FillPatchSingleLevel(
      levelf,
      mf, 
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
}
