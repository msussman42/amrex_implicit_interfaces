
#include <climits>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <Interpolater.H>
#include <INTERP_F.H>

namespace amrex {

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
PCInterpNull              pc_interp_null;
PCInterp                  pc_interp;
LSHOInterp                ls_ho_interp_LOW_PARM;
LSHOInterp                ls_ho_interp_HIGH_PARM;
SEMInterp                 sem_interp_DEFAULT;
SEMInterp                 sem_interp_LOW_PARM;
SEMInterp                 sem_interp_HIGH_PARM;
multiMOFInterp            multi_mof_interp;
multiEXTMOFInterp         multi_extmof_interp;
BurnVelInterp             burnvel_interp;
BurnVelInterp             tsat_interp;
UMACInterp                umac_interp;
VMACInterp                vmac_interp;
WMACInterp                wmac_interp;
UMACInterp                xd_mac_interp;
VMACInterp                yd_mac_interp;
WMACInterp                zd_mac_interp;
maskSEMInterp             mask_sem_interp;


static
Vector<int>
GetBCArray (const Vector<BCRec>& bcr)
{
    Vector<int> bc(2*AMREX_SPACEDIM*bcr.size());

    for (int n = 0; n < bcr.size(); n++)
    {
        const int* b_rec = bcr[n].vect();

        for (int m = 0; m < 2*AMREX_SPACEDIM; m++)
        {
            bc[2*AMREX_SPACEDIM*n + m] = b_rec[m];
        }
    }

    return bc;
}

InterpolaterBoxCoarsener::InterpolaterBoxCoarsener(
 Interpolater* mapper_,int bfactc_,int bfactf_) {

 mapper=mapper_;
 bfactc=bfactc_;
 bfactf=bfactf_;

}

InterpolaterBoxCoarsener
Interpolater::BoxCoarsener (int bfactc,int bfactf)
{
    return InterpolaterBoxCoarsener(this, bfactc,bfactf);
}

Box
InterpolaterBoxCoarsener::doit (const Box& fine) const
{
    return mapper->CoarseBox(fine, bfactc,bfactf);
}

BoxConverter*
InterpolaterBoxCoarsener::clone () const
{
    return new InterpolaterBoxCoarsener(mapper, bfactc,bfactf);
}


Interpolater::~Interpolater () {}


multiMOFInterp::~multiMOFInterp () {}

Box
multiMOFInterp::CoarseBox (const Box& fine,int bfactc,int bfactf)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 Box crse = amrex::coarsen(fine,2);
 
 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;       
}

void
multiMOFInterp::interp (Real time,
  const FArrayBox& crse,
  int               crse_comp,
  FArrayBox&        fine,
  int               fine_comp,
  int               ncomp,
  const Box&        fine_region,
  const Geometry&   crse_geom,
  const Geometry&   fine_geom,
  Vector<BCRec>&     bcr,
  int levelc,int levelf,
  int bfactc,int bfactf)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo = crse.box().loVect();
    const int* chi = crse.box().hiVect();
    const int* flo = fine.loVect();
    const int* fhi = fine.hiVect();
    const int* lo  = fine_region.loVect();
    const int* hi  = fine_region.hiVect();

    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);

    const Real* prob_lo=fine_geom.ProbLo();
    const Real* dxf = fine_geom.CellSize();
    const Real* dxc = crse_geom.CellSize();

    int nmat=multiMOFInterp_nmat;
    int ngeom_raw=multiMOFInterp_ngeom_raw;
    int ngeom_recon=multiMOFInterp_ngeom_recon;

    if (ngeom_raw!=AMREX_SPACEDIM+1)
     amrex::Error("ngeom_raw invalid");
    if (ngeom_recon!=2*AMREX_SPACEDIM+3)
     amrex::Error("ngeom_recon invalid");

    if (nmat<1)
     amrex::Error("nmat invalid in multi mof interp");

    if (ncomp!=nmat*ngeom_raw) {
     std::cout << "ncomp " << ncomp << '\n';
     amrex::Error("must interpolate all multiMOF data at once");
    }

    Box reconbox(crse.box());
    FArrayBox* reconfab=new FArrayBox(reconbox,nmat*ngeom_recon);

     // 2 loops:
     // 1. MOF reconstruction on coarse level
     // 2. interpolate from coarse to fine (traverse fine grid)
    FORT_MULTIMOFINTERP(
     &time,
     cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
     clo,chi,
     fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
     lo,hi,
     reconfab->dataPtr(),
     AMREX_ARLIM(reconfab->loVect()),AMREX_ARLIM(reconfab->hiVect()),
     prob_lo,dxf,dxc,&nmat,
     &ngeom_recon,&ngeom_raw,
     &levelc,&levelf,
     &bfactc,&bfactf);

    delete reconfab;
}



multiEXTMOFInterp::~multiEXTMOFInterp () {}

Box
multiEXTMOFInterp::CoarseBox (const Box& fine,int bfactc,int bfactf)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 Box crse = amrex::coarsen(fine,2);
 
 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;       
}

void
multiEXTMOFInterp::interp (Real time,
  const FArrayBox& crse,
  int               crse_comp,
  FArrayBox&        fine,
  int               fine_comp,
  int               ncomp,
  const Box&        fine_region,
  const Geometry&   crse_geom,
  const Geometry&   fine_geom,
  Vector<BCRec>&     bcr,
  int levelc,int levelf,
  int bfactc,int bfactf)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo = crse.box().loVect();
    const int* chi = crse.box().hiVect();
    const int* flo = fine.loVect();
    const int* fhi = fine.hiVect();
    const int* lo  = fine_region.loVect();
    const int* hi  = fine_region.hiVect();

    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);

    const Real* prob_lo=fine_geom.ProbLo();
    const Real* dxf = fine_geom.CellSize();
    const Real* dxc = crse_geom.CellSize();

    int nmat=multiMOFInterp_nmat;
    int ngeom_raw=multiMOFInterp_ngeom_raw;
    int ngeom_recon=multiMOFInterp_ngeom_recon;

    if (ngeom_raw!=AMREX_SPACEDIM+1)
     amrex::Error("ngeom_raw invalid");
    if (ngeom_recon!=2*AMREX_SPACEDIM+3)
     amrex::Error("ngeom_recon invalid");

    if (nmat<1)
     amrex::Error("nmat invalid in multi ext mof interp");

    if (ncomp!=nmat*ngeom_recon) {
     std::cout << "ncomp " << ncomp << '\n';
     amrex::Error("must interpolate all multiEXTMOF data at once");
    }
    // in NavierStokes::VOF_Recon
    // 1. get MOF data with 1 ghost cell (so that CMOF can be chosen)
    // 2. reconstruct interior cells only.
    // 3. do extended filpatch; MOF used for coarse/fine and ext_dir cells.
    FORT_MULTIEXTMOFINTERP(
     &time,
     cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
     fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
     lo,hi,
     prob_lo,
     dxf,dxc,&nmat,
     &ngeom_recon,&ngeom_raw,
     &levelc,&levelf,
     &bfactc,&bfactf);

}



BurnVelInterp::~BurnVelInterp () {}

Box
BurnVelInterp::CoarseBox (const Box& fine,int bfactc,int bfactf)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 Box crse = amrex::coarsen(fine,2);
 
 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;       
}

void
BurnVelInterp::interp (Real time,
  const FArrayBox& crse,
  int               crse_comp,
  FArrayBox&        fine,
  int               fine_comp,
  int               ncomp,
  const Box&        fine_region,
  const Geometry&   crse_geom,
  const Geometry&   fine_geom,
  Vector<BCRec>&     bcr,
  int levelc,int levelf,
  int bfactc,int bfactf)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo = crse.box().loVect();
    const int* chi = crse.box().hiVect();
    const int* flo = fine.loVect();
    const int* fhi = fine.hiVect();
    const int* lo  = fine_region.loVect();
    const int* hi  = fine_region.hiVect();

    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);

    if ((crse_comp>=0)&&(fine_comp>=0)) {
     // do nothing
    } else
     amrex::Error("crse_comp or fine_comp invalid");

    const Real* prob_lo=fine_geom.ProbLo();
    const Real* dxf = fine_geom.CellSize();
    const Real* dxc = crse_geom.CellSize();

    int nmat=burnvel_nmat;
    int nten=burnvel_nten;
    int ncomp_check=nten+nten*burnvel_ncomp_per;

    int velflag=0;

    if (burnvel_ncomp_per==2) { // interface temperature, mass fraction
     velflag=0;
    } else if (burnvel_ncomp_per==AMREX_SPACEDIM) {
     velflag=1;
    } else
     amrex::Error("burnvel_ncomp_per invalid");

    if ((crse.nComp()>=ncomp_check+crse_comp)&&
        (fine.nComp()>=ncomp_check+fine_comp)) {
     // do nothing
    } else
     amrex::Error("crse.nComp() or fine.nComp() invalid");

    if (nten!=((nmat-1)*(nmat-1)+nmat-1)/2) 
     amrex::Error("nten invalid");

    if (nmat<1)
     amrex::Error("nmat invalid in burnvel interp");

    if (ncomp!=ncomp_check) {
     std::cout << "ncomp " << ncomp << '\n';
     amrex::Error("must interpolate all burnvel data at once");
    }
     // first nmat components are the status.
     // next sdim * nmat components are the burning velocities.
    FORT_EXT_BURNVEL_INTERP(
     &velflag,
     &time,
     cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
     fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
     lo,hi,
     prob_lo,
     dxf,dxc,
     &nmat,
     &nten,
     &ncomp,
     &levelc,&levelf,
     &bfactc,&bfactf);

} // end subroutine BurnVelInterp::interp




PCInterp::~PCInterp () {}

Box
PCInterp::CoarseBox (const Box& fine,int bfactc,int bfactf)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 Box crse = amrex::coarsen(fine,2);

 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;
}

void
PCInterp::interp (
 Real time,
 const FArrayBox& crse,
 int              crse_comp,
 FArrayBox&       fine,
 int              fine_comp,
 int              ncomp,
 const Box&       fine_region,
 const Geometry&  crse_geom,
 const Geometry&  fine_geom,
 Vector<BCRec>&     bcr,
 int levelc,int levelf,
 int bfactc,int bfactf)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo  = crse.box().loVect();
    const int* chi  = crse.box().hiVect();
    const int* flo  = fine.loVect();
    const int* fhi  = fine.hiVect();
    const int* fblo = fine_region.loVect();
    const int* fbhi = fine_region.hiVect();

    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);
    int zapflag=0;

    FORT_PCINTERP (
     &zapflag,
     cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
     fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
     fblo,fbhi,
     &ncomp,
     &levelc,&levelf,
     &bfactc,&bfactf);
}


LSHOInterp::~LSHOInterp () {}

Box
LSHOInterp::CoarseBox (const Box& fine,int bfactc,int bfactf)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 Box crse = amrex::coarsen(fine,2);

 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;
}

void
LSHOInterp::interp (
 Real time,
 const FArrayBox& crse,
 int              crse_comp,
 FArrayBox&       fine,
 int              fine_comp,
 int              ncomp,
 const Box&       fine_region,
 const Geometry&  crse_geom,
 const Geometry&  fine_geom,
 Vector<BCRec>&     bcr,
 int levelc,int levelf,
 int bfactc,int bfactf)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo  = crse.box().loVect();
    const int* chi  = crse.box().hiVect();
    const int* flo  = fine.loVect();
    const int* fhi  = fine.hiVect();
    const int* fblo = fine_region.loVect();
    const int* fbhi = fine_region.hiVect();

    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);

    const Real* prob_lo=fine_geom.ProbLo();
    const Real* dxf = fine_geom.CellSize();
    const Real* dxc = crse_geom.CellSize();

    int nmat=LSHOInterp_nmat;

    if ((LSHOInterp_LO!=0)&&(LSHOInterp_LO!=1))
     amrex::Error("LSHOInterp_LO invalid");

    if (nmat<1)
     amrex::Error("nmat invalid in ls ho interp");
    if (ncomp!=(AMREX_SPACEDIM+1)*nmat) {
     std::cout << "ncomp " << ncomp << '\n';
     amrex::Error("must interpolate all ls ho data at once");
    }

    FORT_LSHOINTERP (
     &LSHOInterp_LO,
     cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
     clo,chi,
     fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
     fblo,fbhi,
     prob_lo,
     dxf,dxc,
     &nmat,
     &ncomp,
     &levelc,&levelf,
     &bfactc,&bfactf);

} //subroutine LSHOInterp::interp



SEMInterp::~SEMInterp () {}

Box
SEMInterp::CoarseBox (const Box& fine,int bfactc,int bfactf)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 Box crse = amrex::coarsen(fine,2);

 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;
}

void
SEMInterp::interp (
 Real time,
 const FArrayBox& crse,
 int              crse_comp,
 FArrayBox&       fine,
 int              fine_comp,
 int              ncomp,
 const Box&       fine_region,
 const Geometry&  crse_geom,
 const Geometry&  fine_geom,
 Vector<BCRec>&     bcr,
 int levelc,int levelf,
 int bfactc,int bfactf)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo  = crse.box().loVect();
    const int* chi  = crse.box().hiVect();
    const int* flo  = fine.loVect();
    const int* fhi  = fine.hiVect();
    const int* fblo = fine_region.loVect();
    const int* fbhi = fine_region.hiVect();

    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);

    const Real* dxf = fine_geom.CellSize();
    const Real* dxc = crse_geom.CellSize();

      // enable_spectral:
      // 0 - low order
      // 1 - space/time spectral
      // 2 - space spectral only
      // 3 - time spectral only
    FORT_SEMINTERP (
     &interp_enable_spectral,
     dxc,dxf, 
     cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
     fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),fblo,fbhi,
     &ncomp,
     &levelc,&levelf,
     &bfactc,&bfactf);
}

maskSEMInterp::~maskSEMInterp () {}

Box
maskSEMInterp::CoarseBox (const Box& fine,int bfactc,int bfactf)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

  // the smallest coarse box whose refinement contains "fine"
 Box crse = amrex::coarsen(fine,2);

 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;
}

void
maskSEMInterp::interp (
 Real time,
 const FArrayBox& crse,
 int              crse_comp,
 FArrayBox&       fine,
 int              fine_comp,
 int              ncomp,
 const Box&       fine_region,
 const Geometry&  crse_geom,
 const Geometry&  fine_geom,
 Vector<BCRec>&     bcr,
 int levelc,int levelf,
 int bfactc,int bfactf)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo  = crse.box().loVect();
    const int* chi  = crse.box().hiVect();
    const int* flo  = fine.loVect();
    const int* fhi  = fine.hiVect();
    const int* fblo = fine_region.loVect();
    const int* fbhi = fine_region.hiVect();


    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);

    FORT_MASKINTERPPC (
     cdat,
     AMREX_ARLIM(clo),AMREX_ARLIM(chi),
     fdat,
     AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
     fblo,fbhi,
     &ncomp,
     &levelc,&levelf,
     &bfactc,&bfactf);
}




PCInterpNull::~PCInterpNull () {}

Box
PCInterpNull::CoarseBox (const Box& fine,int bfactc,int bfactf)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 Box crse=amrex::coarsen(fine,2);

 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;
}

void
PCInterpNull::interp (
 Real time,
 const FArrayBox& crse,
 int              crse_comp,
 FArrayBox&       fine,
 int              fine_comp,
 int              ncomp,
 const Box&       fine_region,
 const Geometry&  crse_geom,
 const Geometry&  fine_geom,
 Vector<BCRec>&     bcr,
 int levelc,int levelf,
 int bfactc,int bfactf)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo  = crse.box().loVect();
    const int* chi  = crse.box().hiVect();
    const int* flo  = fine.loVect();
    const int* fhi  = fine.hiVect();
    const int* fblo = fine_region.loVect();
    const int* fbhi = fine_region.hiVect();


    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);

    int zapflag=1;
    FORT_PCINTERP (
     &zapflag,
     cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
     fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),fblo,fbhi,
     &ncomp,
     &levelc,&levelf,
     &bfactc,&bfactf);
}

UMACInterp::~UMACInterp() {}


Box
UMACInterp::CoarseBox(const Box& fine,int bfactc,int bfactf)
{

 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (fine.ixType()!=IndexType::TheUMACType())
  amrex::Error("error in CoarseBox");

 Box crse = amrex::coarsen(fine,2);
 if (crse.ixType()!=IndexType::TheUMACType())
  amrex::Error("error in CoarseBox: crse");

 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;
}

void
UMACInterp::interp(
 Real time,
 const FArrayBox& crse, int crse_comp,
 FArrayBox& fine, int fine_comp,
 int ncomp,
 const Box& fine_region, 
 const Geometry& crse_geom,
 const Geometry& fine_geom,
 Vector<BCRec>&     bcr,
 int levelc,int levelf,
 int bfactc,int bfactf)
{
 if ((ncomp<1)||(ncomp>9999))
  amrex::Error("invalid ncomp umac interp");

 BL_ASSERT(bcr.size() >= ncomp);

 Vector<int> bcfine = GetBCArray(bcr);

 IndexType typ(fine_region.ixType());

 if (typ!=IndexType::TheUMACType())
  amrex::Error("fine_region has incorrect type");

 if (fine.box().ixType()!=typ)
  amrex::Error("fine box invalid");
 if (crse.box().ixType()!=typ)
  amrex::Error("crse box invalid");

 Box fine_bx = fine_region & fine.box();

 Box crse_bx(CoarseBox(fine_bx,bfactc,bfactf));
 if (crse_bx.ixType()!=typ)
  amrex::Error("crse_bx invalid");

 const Real* prob_lo=fine_geom.ProbLo();
 const Real* dxf = fine_geom.CellSize();
 const Real* dxc = crse_geom.CellSize();

 int dir=0;

  // enable_spectral:
  // 0 - low order
  // 1 - space/time spectral
  // 2 - space spectral only
  // 3 - time spectral only
 FORT_EDGEINTERP(
   &interp_enable_spectral,
   &dir,
   crse.dataPtr(crse_comp),
   AMREX_ARLIM(crse.loVect()),AMREX_ARLIM(crse.hiVect()),
   crse_bx.loVect(),crse_bx.hiVect(),
   fine.dataPtr(fine_comp),
   AMREX_ARLIM(fine.loVect()),AMREX_ARLIM(fine.hiVect()),
   fine_bx.loVect(),fine_bx.hiVect(),
   prob_lo,dxf,dxc,
   &ncomp,
   &levelc,&levelf,
   &bfactc,&bfactf);
}




VMACInterp::~VMACInterp() {}


Box
VMACInterp::CoarseBox(const Box& fine,int bfactc,int bfactf)
{


 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (fine.ixType()!=IndexType::TheVMACType())
  amrex::Error("error in VMAC CoarseBox");

 Box crse = amrex::coarsen(fine,2);
 if (crse.ixType()!=IndexType::TheVMACType())
  amrex::Error("error in VMAC CoarseBox: crse");

 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;

}

void
VMACInterp::interp(
 Real time,
 const FArrayBox& crse, int crse_comp,
 FArrayBox& fine, int fine_comp,
 int ncomp,
 const Box& fine_region, 
 const Geometry& crse_geom,
 const Geometry& fine_geom,
 Vector<BCRec>&     bcr,
 int levelc,int levelf,
 int bfactc,int bfactf)
{
 if ((ncomp<1)||(ncomp>9999))
  amrex::Error("invalid ncomp vmac interp");

 BL_ASSERT(bcr.size() >= ncomp);

 Vector<int> bcfine = GetBCArray(bcr);

 IndexType typ(fine_region.ixType());

 if (typ!=IndexType::TheVMACType())
  amrex::Error("fine_region has incorrect type");

 if (fine.box().ixType()!=typ)
  amrex::Error("fine box invalid");
 if (crse.box().ixType()!=typ)
  amrex::Error("crse box invalid");

 Box fine_bx = fine_region & fine.box();
 Box crse_bx(CoarseBox(fine_bx,bfactc,bfactf));
 if (crse_bx.ixType()!=typ)
  amrex::Error("crse_bx invalid");

 const Real* prob_lo=fine_geom.ProbLo();
 const Real* dxf = fine_geom.CellSize();
 const Real* dxc = crse_geom.CellSize();

 int dir=1;

  // enable_spectral:
  // 0 - low order
  // 1 - space/time spectral
  // 2 - space spectral only
  // 3 - time spectral only
 FORT_EDGEINTERP(
   &interp_enable_spectral,
   &dir,
   crse.dataPtr(crse_comp),
   AMREX_ARLIM(crse.loVect()),AMREX_ARLIM(crse.hiVect()),
   crse_bx.loVect(),crse_bx.hiVect(),
   fine.dataPtr(fine_comp),
   AMREX_ARLIM(fine.loVect()),AMREX_ARLIM(fine.hiVect()),
   fine_bx.loVect(),fine_bx.hiVect(),
   prob_lo,dxf,dxc,
   &ncomp,
   &levelc,&levelf,
   &bfactc,&bfactf);
}


WMACInterp::~WMACInterp() {}


Box
WMACInterp::CoarseBox(const Box& fine,int bfactc,int bfactf)
{

 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (fine.ixType()!=IndexType::TheWMACType())
  amrex::Error("error in WMAC CoarseBox");

 Box crse = amrex::coarsen(fine,2);
 if (crse.ixType()!=IndexType::TheWMACType())
  amrex::Error("error in WMAC CoarseBox: crse");

 if (bfactc>1) {
  Box e_crse = amrex::coarsen(crse,bfactc);
  e_crse.refine(bfactc);
  crse=e_crse;
 } else if (bfactc==1) {
  // do nothing
 } else
  amrex::Error("bfactc invalid");

 return crse;

}

void
WMACInterp::interp(
 Real time,
 const FArrayBox& crse, int crse_comp,
 FArrayBox& fine, int fine_comp,
 int ncomp,
 const Box& fine_region,
 const Geometry& crse_geom,
 const Geometry& fine_geom,
 Vector<BCRec>&     bcr,
 int levelc,int levelf,
 int bfactc,int bfactf)
{
 if ((ncomp<1)||(ncomp>9999))
  amrex::Error("invalid ncomp wmac interp");

 BL_ASSERT(bcr.size() >= ncomp);

 Vector<int> bcfine = GetBCArray(bcr);

 IndexType typ(fine_region.ixType());

 if (typ!=IndexType::TheWMACType())
  amrex::Error("fine_region has incorrect type");

 if (fine.box().ixType()!=typ)
  amrex::Error("fine box invalid");
 if (crse.box().ixType()!=typ)
  amrex::Error("crse box invalid");

 Box fine_bx = fine_region & fine.box();
 Box crse_bx(CoarseBox(fine_bx,bfactc,bfactf));
 if (crse_bx.ixType()!=typ)
  amrex::Error("crse_bx invalid");

 const Real* prob_lo=fine_geom.ProbLo();
 const Real* dxf = fine_geom.CellSize();
 const Real* dxc = crse_geom.CellSize();

 int dir=2;

  // enable_spectral:
  // 0 - low order
  // 1 - space/time spectral
  // 2 - space spectral only
  // 3 - time spectral only
 FORT_EDGEINTERP(
   &interp_enable_spectral,
   &dir,
   crse.dataPtr(crse_comp),
   AMREX_ARLIM(crse.loVect()),AMREX_ARLIM(crse.hiVect()),
   crse_bx.loVect(),crse_bx.hiVect(),
   fine.dataPtr(fine_comp),
   AMREX_ARLIM(fine.loVect()),AMREX_ARLIM(fine.hiVect()),
   fine_bx.loVect(),fine_bx.hiVect(),
   prob_lo,dxf,dxc,
   &ncomp,
   &levelc,&levelf,
   &bfactc,&bfactf);
}

} // namespace amrex
