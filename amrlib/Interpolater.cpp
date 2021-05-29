
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
UMACInterp                xd_mac_interp;
UMACInterp                xd_mac_lo_interp;
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

// CELL - CELL - CELL  (grid_type=-1)
//    x_i = (i+1/2)*dx y_j=(j+1/2)*dy z_k=(k+1/2)*dz
//   xflux: NODE - CELL - CELL ("umac grid")  (grid_type=0)
//   yflux: CELL - NODE - CELL ("vmac grid")  (grid_type=1)
//   zflux: CELL - CELL - NODE ("wmac grid")  (grid_type=2)
// NODE - CELL - CELL (grid_type 0)
//   xflux: CELL - CELL - CELL (grid_type -1)
//   yflux: NODE - NODE - CELL (grid_type= 3)
//   zflux: NODE - CELL - NODE (grid_type= 4)
// CELL - NODE - CELL (grid_type 1)
//   zflux: CELL - NODE - NODE (grid_type= 5)
// constructor
InterpolaterBoxCoarsener::InterpolaterBoxCoarsener(
 Interpolater* mapper_,int bfactc_,int bfactf_,int grid_type_) {

 if ((grid_type_==-1)||
     ((grid_type_>=0)&&(grid_type_<=5))) {
  // do nothing
 } else
  amrex::Error("grid_type_ invalid");

 mapper=mapper_;
 bfactc=bfactc_;
 bfactf=bfactf_;
 grid_type=grid_type_;

}

InterpolaterBoxCoarsener
Interpolater::BoxCoarsener (int bfactc,int bfactf,int grid_type)
{
 if ((grid_type==-1)||
     ((grid_type>=0)&&(grid_type<=5))) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

 return InterpolaterBoxCoarsener(this, bfactc,bfactf,grid_type);
}

Box
InterpolaterBoxCoarsener::doit (const Box& fine) const
{
    return mapper->CoarseBox(fine,bfactc,bfactf,grid_type);
}

BoxConverter*
InterpolaterBoxCoarsener::clone () const
{
    return new InterpolaterBoxCoarsener(mapper,bfactc,bfactf,grid_type);
}


Interpolater::~Interpolater () {}


multiMOFInterp::~multiMOFInterp () {}

Box
multiMOFInterp::CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
  int bfactc,int bfactf,
  int grid_type)
{

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
multiEXTMOFInterp::CoarseBox (const Box& fine,int bfactc,int bfactf,
		int grid_type)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
  int bfactc,int bfactf,
  int grid_type)
{

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
BurnVelInterp::CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
  int bfactc,int bfactf,
  int grid_type)
{

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
PCInterp::CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 Box crse = amrex::coarsen(fine,2);

  // CELL - CELL - CELL
 if (grid_type==-1) {
  if ((fine.ixType()==IndexType::TheCellType())&&
      (crse.ixType()==IndexType::TheCellType)) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");

  // NODE - CELL - CELL
 } else if (grid_type==0) {
  if ((fine.ixType()==IndexType::TheUMACType())&&
      (crse.ixType()==IndexType::TheUMACType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");

   // CELL - NODE - CELL
 } else if (grid_type==1) {
  if ((fine.ixType()==IndexType::TheVMACType())&&
      (crse.ixType()==IndexType::TheVMACType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");

   // CELL - CELL - NODE
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((fine.ixType()==IndexType::TheWMACType())&&
      (crse.ixType()==IndexType::TheWMACType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
   // NODE - NODE - CELL 
 } else if (grid_type==3) {
  if ((fine.ixType()==IndexType::TheYUMACType())&&
      (crse.ixType()==IndexType::TheYUMACType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
   // NODE - CELL - NODE 
 } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
  if ((fine.ixType()==IndexType::TheZUMACType())&&
      (crse.ixType()==IndexType::TheZUMACType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
   // CELL - NODE - NODE
 } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
  if ((fine.ixType()==IndexType::TheZVMACType())&&
      (crse.ixType()==IndexType::TheZVMACType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
 } else
  amrex::Error("grid_type invalid");

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
 int bfactc,int bfactf,
 int grid_type)
{

 if ((ncomp<1)||(ncomp>9999))
  amrex::Error("invalid ncomp PCInterp::interp");

 BL_ASSERT(bcr.size() >= ncomp);

 IndexType typ(fine_region.ixType());

 Box fine_bx = fine_region & fine.box();

 Box crse_bx(CoarseBox(fine_bx,bfactc,bfactf,grid_type));

  // CELL - CELL - CELL
 if (grid_type==-1) {
  if ((typ==IndexType::TheCellType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
  // NODE - CELL - CELL
 } else if (grid_type==0) {
  if ((typ==IndexType::TheUMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - NODE - CELL
 } else if (grid_type==1) {
  if ((typ==IndexType::TheVMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - CELL - NODE
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((typ==IndexType::TheWMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // NODE - NODE - CELL 
 } else if (grid_type==3) {
  if ((typ==IndexType::TheYUMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // NODE - CELL - NODE 
 } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
  if ((typ==IndexType::TheZUMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - NODE - NODE
 } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
  if ((typ==IndexType::TheZVMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
 } else
  amrex::Error("grid_type invalid");

 const Real* prob_lo=fine_geom.ProbLo();
 const Real* dxf = fine_geom.CellSize();
 const Real* dxc = crse_geom.CellSize();

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
  &grid_type,
  &zapflag,
  cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
  crse_bx.loVect(),crse_bx.hiVect(),
  fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
  fine_bx.loVect(),fine_bx.hiVect(),
  fblo,fbhi,
  prob_lo,dxf,dxc,
  &ncomp,
  &levelc,&levelf,
  &bfactc,&bfactf);
}


LSHOInterp::~LSHOInterp () {}

Box
LSHOInterp::CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
 int bfactc,int bfactf,
 int grid_type)
{

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
SEMInterp::CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
 int bfactc,int bfactf,
 int grid_type)
{

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
maskSEMInterp::CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
 int bfactc,int bfactf,
 int grid_type)
{

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
PCInterpNull::CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type)
{
 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 if (grid_type==-1) {
  // do nothing
 } else
  amrex::Error("grid_type invalid");

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
 int bfactc,int bfactf,
 int grid_type)
{

 if ((ncomp<1)||(ncomp>9999))
  amrex::Error("invalid ncomp PCInterpNull::interp");

 BL_ASSERT(bcr.size() >= ncomp);

 IndexType typ(fine_region.ixType());

 Box fine_bx = fine_region & fine.box();

 Box crse_bx(CoarseBox(fine_bx,bfactc,bfactf,grid_type));

  // CELL - CELL - CELL
 if (grid_type==-1) {
  if ((typ==IndexType::TheCellType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
  // NODE - CELL - CELL
 } else if (grid_type==0) {
  if ((typ==IndexType::TheUMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - NODE - CELL
 } else if (grid_type==1) {
  if ((typ==IndexType::TheVMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - CELL - NODE
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((typ==IndexType::TheWMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // NODE - NODE - CELL 
 } else if (grid_type==3) {
  if ((typ==IndexType::TheYUMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // NODE - CELL - NODE 
 } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
  if ((typ==IndexType::TheZUMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - NODE - NODE
 } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
  if ((typ==IndexType::TheZVMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
 } else
  amrex::Error("grid_type invalid");

 const Real* prob_lo=fine_geom.ProbLo();
 const Real* dxf = fine_geom.CellSize();
 const Real* dxc = crse_geom.CellSize();

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
  &grid_type,
  &zapflag,
  cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
  crse_bx.loVect(),crse_bx.hiVect(),
  fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
  fblo,fbhi,
  prob_lo,dxf,dxc,
  &ncomp,
  &levelc,&levelf,
  &bfactc,&bfactf);
}

UMACInterp::~UMACInterp() {}


Box
UMACInterp::CoarseBox(const Box& fine,int bfactc,int bfactf,int grid_type)
{

 if ((bfactc<1)||(bfactf<1))
  amrex::Error("bfactc or bfactf invalid");
 if (bfactf>bfactc)
  amrex::Error("cannot have bfactf>bfactc");

 Box crse = amrex::coarsen(fine,2);

 if (grid_type==0) {
  if ((fine.ixType()==IndexType::TheUMACType())&&
      (crse.ixType()==IndexType::TheUMACType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
 } else if (grid_type==1) {
  if ((fine.ixType()==IndexType::TheVMACType())&&
      (crse.ixType()==IndexType::TheVMACType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((fine.ixType()==IndexType::TheWMACType())&&
      (crse.ixType()==IndexType::TheWMACType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
 } else
  amrex::Error("grid_type invalid");

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
 const FArrayBox& crse, 
 int crse_comp,
 FArrayBox& fine, 
 int fine_comp,
 int ncomp,
 const Box& fine_region, 
 const Geometry& crse_geom,
 const Geometry& fine_geom,
 Vector<BCRec>&     bcr,
 int levelc,int levelf,
 int bfactc,int bfactf,
 int grid_type)
{

 if ((ncomp<1)||(ncomp>9999))
  amrex::Error("invalid ncomp umac interp");

 BL_ASSERT(bcr.size() >= ncomp);

 IndexType typ(fine_region.ixType());

 Box fine_bx = fine_region & fine.box();

 Box crse_bx(CoarseBox(fine_bx,bfactc,bfactf,grid_type));

 if (grid_type==0) {
  if ((typ==IndexType::TheUMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
 } else if (grid_type==1) {
  if ((typ==IndexType::TheVMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((typ==IndexType::TheWMACType())&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
 } else
  amrex::Error("grid_type invalid");

 const Real* prob_lo=fine_geom.ProbLo();
 const Real* dxf = fine_geom.CellSize();
 const Real* dxc = crse_geom.CellSize();

  // enable_spectral:
  // 0 - low order
  // 1 - space/time spectral
  // 2 - space spectral only
  // 3 - time spectral only
 FORT_EDGEINTERP(
   &interp_enable_spectral,
   &grid_type,
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
