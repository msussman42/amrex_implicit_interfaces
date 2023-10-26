
#include <climits>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <Interpolater.H>
#include <INTERP_F.H>
#include <INDEX_TYPE_MACROS.H>
#include <EXTRAP_COMP.H>

namespace amrex {

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
PCInterpNull              pc_interp_null;
PCInterp                  pc_interp;
LSInterp                  ls_interp;
SEMInterp                 sem_interp_DEFAULT;
SEMInterp                 sem_interp_LOW_PARM;
SEMInterp                 sem_interp_HIGH_PARM;
multiMOFInterp            multi_mof_interp;
multiEXTMOFInterp         multi_extmof_interp;
BurnVelInterp             burnvel_interp;
BurnVelInterp             tsat_interp;
BurnVelInterp             drag_interp;
UMACInterp                umac_interp;
UMACInterp                xd_mac_interp;
UMACInterp                xd_mac_lo_interp;
PCInterp                  tensor_pc_interp;
maskSEMInterp             mask_sem_interp;


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
  Vector<BCRec>& /* bcr */,
  int levelc,int levelf,
  int bfactc,int bfactf,
  int grid_type)
{

 if (time>=0.0) {
  //do nothing
 } else 
  amrex::Error("time invalid");

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

 int local_nmat=multiMOFInterp_nmat;
 int ngeom_raw=multiMOFInterp_ngeom_raw;
 int ngeom_recon=multiMOFInterp_ngeom_recon;

 if (ngeom_raw!=AMREX_SPACEDIM+1)
  amrex::Error("ngeom_raw invalid");
 if (ngeom_recon!=2*AMREX_SPACEDIM+3)
  amrex::Error("ngeom_recon invalid");

 if (local_nmat<1)
  amrex::Error("local_nmat invalid in multi mof interp");

 if (ncomp!=local_nmat*ngeom_raw) {
  std::cout << "ncomp " << ncomp << '\n';
  amrex::Error("must interpolate all multiMOF data at once");
 }

 Box reconbox(crse.box());
 FArrayBox* reconfab=new FArrayBox(reconbox,local_nmat*ngeom_recon);

  // 2 loops:
  // 1. MOF reconstruction on coarse level
  // 2. interpolate from coarse to fine (traverse fine grid)
 fort_multimofinterp(
  &time,
  cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
  clo,chi,
  fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
  lo,hi,
  reconfab->dataPtr(),
  AMREX_ARLIM(reconfab->loVect()),AMREX_ARLIM(reconfab->hiVect()),
  prob_lo,dxf,dxc,
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
  Vector<BCRec>& /* bcr */,
  int levelc,int levelf,
  int bfactc,int bfactf,
  int grid_type)
{

 if (time>=0.0) {
  //do nothing
 } else 
  amrex::Error("time invalid");

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

 int local_nmat=multiMOFInterp_nmat;
 int ngeom_raw=multiMOFInterp_ngeom_raw;
 int ngeom_recon=multiMOFInterp_ngeom_recon;

 if (ngeom_raw!=AMREX_SPACEDIM+1)
  amrex::Error("ngeom_raw invalid");
 if (ngeom_recon!=2*AMREX_SPACEDIM+3)
  amrex::Error("ngeom_recon invalid");

 if (local_nmat<1)
  amrex::Error("local_nmat invalid in multi ext mof interp");

 if (ncomp!=local_nmat*ngeom_recon) {
  std::cout << "ncomp " << ncomp << '\n';
  amrex::Error("must interpolate all multiEXTMOF data at once");
 }
 // in NavierStokes::VOF_Recon
 // 1. get MOF data with 1 ghost cell (so that CMOF can be chosen)
 // 2. reconstruct interior cells only.
 // 3. do extended filpatch; MOF used for coarse/fine and ext_dir cells.
 fort_multiextmofinterp(
  &time,
  cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
  fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
  lo,hi,
  prob_lo,
  dxf,dxc,
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
  Vector<BCRec>& /* bcr */,
  int levelc,int levelf,
  int bfactc,int bfactf,
  int grid_type)
{

 if (time>=0.0) {
  //do nothing
 } else 
  amrex::Error("time invalid");

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

 int local_nmat=burnvel_nmat;
 int local_nten=burnvel_nten;
 int ncomp_check=local_nten+local_nten*burnvel_ncomp_per;

 if (burnvel_ncomp_per==0) {
  int num_materials=local_nmat; //the macro N_DRAG depends on "num_materials"
  if (num_materials>=2) {
   ncomp_check=N_DRAG;
  } else
   amrex::Error("BurnVelInterp::interp => need num_materials>=2");
 } else if (burnvel_ncomp_per>0) {
  // do nothing
 } else
  amrex::Error("burnvel_ncomp_per invalid");

 if (ncomp_check==burnvel_ncomp) {
  // do nothing
 } else
  amrex::Error("ncomp_check invalid");

 int velflag=0;

  //interface temperature,mass fraction
 if (burnvel_ncomp_per==EXTRAP_PER_TSAT) {
  velflag=0;
 } else if ((burnvel_ncomp_per==EXTRAP_PER_BURNING)&&
            (EXTRAP_PER_TSAT!=EXTRAP_PER_BURNING)) {
  velflag=1;
 } else if (burnvel_ncomp_per==0) {
  velflag=2;
 } else
  amrex::Error("burnvel_ncomp_per invalid");

 if ((crse.nComp()>=ncomp_check+crse_comp)&&
     (fine.nComp()>=ncomp_check+fine_comp)) {
  // do nothing
 } else
  amrex::Error("crse.nComp() or fine.nComp() invalid");

 if (local_nten!=((local_nmat-1)*(local_nmat-1)+local_nmat-1)/2) 
  amrex::Error("local_nten invalid");

 if (local_nmat<1)
  amrex::Error("local_nmat invalid in burnvel interp");

 if (ncomp!=ncomp_check) {
  std::cout << "ncomp " << ncomp << '\n';
  amrex::Error("must interpolate all burnvel data at once");
 }
  // if (velflag=0 or 1):
  // first num_materials components are the status.
  // next sdim * num_materials components are the burning velocities.
  
 fort_ext_burnvel_interp(
  &velflag,
  &time,
  cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
  fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
  lo,hi,
  prob_lo,
  dxf,dxc,
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
      (crse.ixType()==IndexType::TheCellType())) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");

  // NODE - CELL - CELL
 } else if (grid_type==0) {
  if ((fine.ixType()==TheUMACType)&&
      (crse.ixType()==TheUMACType)) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");

   // CELL - NODE - CELL
 } else if (grid_type==1) {
  if ((fine.ixType()==TheVMACType)&&
      (crse.ixType()==TheVMACType)) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");

   // CELL - CELL - NODE
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((fine.ixType()==TheWMACType)&&
      (crse.ixType()==TheWMACType)) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
   // NODE - NODE - CELL 
 } else if (grid_type==3) {
  if ((fine.ixType()==TheYUMACType)&&
      (crse.ixType()==TheYUMACType)) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
   // NODE - CELL - NODE 
 } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
  if ((fine.ixType()==TheZUMACType)&&
      (crse.ixType()==TheZUMACType)) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
   // CELL - NODE - NODE
 } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
  if ((fine.ixType()==TheZVMACType)&&
      (crse.ixType()==TheZVMACType)) {
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

 if (time>=0.0) {
  //do nothing
 } else 
  amrex::Error("time invalid");

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
  if ((typ==TheUMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - NODE - CELL
 } else if (grid_type==1) {
  if ((typ==TheVMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - CELL - NODE
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((typ==TheWMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // NODE - NODE - CELL 
 } else if (grid_type==3) {
  if ((typ==TheYUMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // NODE - CELL - NODE 
 } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
  if ((typ==TheZUMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - NODE - NODE
 } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
  if ((typ==TheZVMACType)&&
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
// const int* fblo = fine_region.loVect();
// const int* fbhi = fine_region.hiVect();

 const Real* cdat  = crse.dataPtr(crse_comp);
 Real*       fdat  = fine.dataPtr(fine_comp);
 int zapflag=0;

   //fine_bx = fine_region & fine.box();
 fort_pcinterp (
  &grid_type,
  &zapflag,
  cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
  crse_bx.loVect(),crse_bx.hiVect(),
  fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
  fine_bx.loVect(),fine_bx.hiVect(),
  prob_lo,dxf,dxc,
  &ncomp,
  &levelc,&levelf,
  &bfactc,&bfactf);
}


LSInterp::~LSInterp () {}

Box
LSInterp::CoarseBox (const Box& fine,int bfactc,int bfactf,int grid_type)
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
LSInterp::interp (
 Real time,
 const FArrayBox& crse,
 int              crse_comp,
 FArrayBox&       fine,
 int              fine_comp,
 int              ncomp,
 const Box&       fine_region,
 const Geometry&  crse_geom,
 const Geometry&  fine_geom,
 Vector<BCRec>& /*bcr */,
 int levelc,int levelf,
 int bfactc,int bfactf,
 int grid_type)
{

 if (time>=0.0) {
  //do nothing
 } else 
  amrex::Error("time invalid");

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

 int local_nmat=LSInterp_nmat;

 if (local_nmat<1)
  amrex::Error("local_nmat invalid in ls interp");
 if (ncomp!=(AMREX_SPACEDIM+1)*local_nmat) {
  std::cout << "ncomp " << ncomp << '\n';
  amrex::Error("must interpolate all ls data at once");
 }

 AMREX_ALWAYS_ASSERT(crse.nComp()>=ncomp);
 AMREX_ALWAYS_ASSERT(fine.nComp()>=ncomp);

 fort_lsinterp (
  cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
  clo,chi,
  fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
  fblo,fbhi,
  prob_lo,
  dxf,dxc,
  &ncomp,
  &levelc,&levelf,
  &bfactc,&bfactf);

} //subroutine LSInterp::interp



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
 Vector<BCRec>& /* bcr */,
 int levelc,int levelf,
 int bfactc,int bfactf,
 int grid_type)
{

 if (time>=0.0) {
  //do nothing
 } else 
  amrex::Error("time invalid");

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
 fort_seminterp (
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
 const Geometry&  /*crse_geom*/,
 const Geometry&  /*fine_geom*/,
 Vector<BCRec>& /* bcr */,
 int levelc,int levelf,
 int bfactc,int bfactf,
 int grid_type)
{

 if (time>=0.0) {
  //do nothing
 } else 
  amrex::Error("time invalid");

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

 fort_maskinterppc (
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

 if (time>=0.0) {
  //do nothing
 } else 
  amrex::Error("time invalid");

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
  if ((typ==TheUMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - NODE - CELL
 } else if (grid_type==1) {
  if ((typ==TheVMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - CELL - NODE
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((typ==TheWMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // NODE - NODE - CELL 
 } else if (grid_type==3) {
  if ((typ==TheYUMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // NODE - CELL - NODE 
 } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
  if ((typ==TheZUMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
   // CELL - NODE - NODE
 } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
  if ((typ==TheZVMACType)&&
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
// const int* fblo = fine_region.loVect();
// const int* fbhi = fine_region.hiVect();

 const Real* cdat  = crse.dataPtr(crse_comp);
 Real*       fdat  = fine.dataPtr(fine_comp);

   //fine_bx = fine_region & fine.box();
 int zapflag=1;
 fort_pcinterp (
  &grid_type,
  &zapflag,
  cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
  crse_bx.loVect(),crse_bx.hiVect(),
  fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
  fine_bx.loVect(),fine_bx.hiVect(),
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
  if ((fine.ixType()==TheUMACType)&&
      (crse.ixType()==TheUMACType)) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
 } else if (grid_type==1) {
  if ((fine.ixType()==TheVMACType)&&
      (crse.ixType()==TheVMACType)) {
   // do nothing
  } else
   amrex::Error("fine or crse box has wrong grid_type");
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((fine.ixType()==TheWMACType)&&
      (crse.ixType()==TheWMACType)) {
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

 if (time>=0.0) {
  //do nothing
 } else 
  amrex::Error("time invalid");

 if ((ncomp<1)||(ncomp>9999))
  amrex::Error("invalid ncomp umac interp");

 BL_ASSERT(bcr.size() >= ncomp);

 IndexType typ(fine_region.ixType());

 Box fine_bx = fine_region & fine.box();

 Box crse_bx(CoarseBox(fine_bx,bfactc,bfactf,grid_type));

 if (grid_type==0) {
  if ((typ==TheUMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
 } else if (grid_type==1) {
  if ((typ==TheVMACType)&&
      (fine.box().ixType()==typ)&&
      (crse.box().ixType()==typ)&&
      (crse_bx.ixType()==typ)) {
   // do nothing
  } else
   amrex::Error("typ invalid");
 } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
  if ((typ==TheWMACType)&&
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
 fort_edgeinterp(
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
