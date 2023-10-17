
#include <iostream>
#include <algorithm>

#include <unistd.h>
#include <cmath>

#include <AMReX_RealBox.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <StateData.H>
#include <StateDescriptor.H>
#include <INDEX_TYPE_MACROS.H>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <INTERP_F.H>

#include <EXTRAP_COMP.H>

namespace amrex {

const Real INVALID_TIME = -1.0e200;

// for each AmrLevel, there is a class of type "StateData" and within StateData
// there is an array of MultiFabs called "new_data"
// 1. class AmrCore
// 2. Vector<std::unique_ptr<AmrLevel> > amr_level; 
// 3. Vector<StateData> state;  
//    cell centered data, MAC data, stress tensor data 
// 4. Vector< MultiFab* > new_data; 
//    for low order in time: new_data[0]  (previous data)
//                           new_data[1]  (new data)
//    spectral deferred correction (a SLAB) in time:
//                           new_data[0] ....
//                           new_data[order+1]
StateData::StateData () 
{
   StateData_level=0;

   StateData_MAX_NUM_SLAB=33;
   StateData_slab_dt_type=0;

   desc = 0;
   descGHOST = 0;

   time_array.resize(StateData_MAX_NUM_SLAB);
   new_data.resize(StateData_MAX_NUM_SLAB);

   for (int i=0;i<StateData_MAX_NUM_SLAB;i++) {
    new_data[i]=0;
    time_array[i] = INVALID_TIME;
   }
   bfact_time_order=0;
}

StateData::StateData (
  AmrCore& papa,
  int level,
  int max_level,
  const Box& p_domain,
  const BoxArray& grds,
  const DistributionMapping& dm,
  const StateDescriptor* d,
  const StateDescriptor* dGHOST,
  Real cur_time,
  Real dt,
  int time_order,
  int slab_dt_type, // 0=SEM 1=evenly spaced
  int MAX_NUM_SLAB)
{
    define(papa,level,max_level,
	   p_domain, grds, dm, *d, *dGHOST, cur_time, dt,
           time_order,slab_dt_type,MAX_NUM_SLAB);
}

void
StateData::define (
  AmrCore& papa,
  int level,
  int max_level,
  const Box& p_domain,
  const BoxArray& grds,
  const DistributionMapping& dm,
  const StateDescriptor& d,
  const StateDescriptor& dGHOST,
  Real time,
  Real dt,
  int time_order,
  int slab_dt_type, // 0=SEM 1=evenly spaced
  int MAX_NUM_SLAB)
{

    AmrCore* parent=&papa;

    if (time_order==parent->Time_blockingFactor()) {
     //do nothing
    } else
     amrex::Error("expecting time_order==parent->Time_blockingFactor()");

    if (max_level>=0) {
     // do nothing
    } else
     amrex::Error("max_level>=0 violated");

    if (dt<=0.0) {
     std::cout << "dt = " << dt << '\n';
     amrex::Error("dt invalid in StateData define");
    }
    if (time<0.0) {
     std::cout << "time = " << time << '\n';
     amrex::Error("time invalid in StateData define");
    }
   
    StateData_level=level;

    StateData_MAX_NUM_SLAB=MAX_NUM_SLAB;
    if (StateData_MAX_NUM_SLAB<33)
     amrex::Error("StateData_MAX_NUM_SLAB too small");

    StateData_slab_dt_type=slab_dt_type;
    if ((StateData_slab_dt_type!=0)&&
        (StateData_slab_dt_type!=1))
     amrex::Error("StateData_slab_dt_type invalid");
 
    bfact_time_order=time_order;
    if ((bfact_time_order>StateData_MAX_NUM_SLAB)||
        (bfact_time_order<1)) {
     std::cout << "bfact_time_order= " << bfact_time_order << '\n';
     amrex::Error("bfact_time_order invalid in define");
    }
    time_array.resize(StateData_MAX_NUM_SLAB);
    new_data.resize(StateData_MAX_NUM_SLAB);

    domain = p_domain;
    desc = &d;
    descGHOST = &dGHOST;
    grids = grds;
    dmap = dm;
    //
    // Convert to proper type.
    //
    IndexType typ(desc->getType());
    if (!typ.cellCentered())
    {
        domain.convert(typ);
        grids.convert(typ);
    }
      //time_array[0]=time-dt
      //time_array[bfact_time_order]=time
    Real slablow=time-dt;
    Real slabhigh=time;
    int do_scale_time=0;

    if (dt>=1.0) {
     do_scale_time=1;
     slablow=time/dt-1.0;
     slabhigh=time/dt;
    } else if ((dt>0.0)&&(dt<1.0)) {
     // do nothing
    } else {
     amrex::Error("dt invalid");
    }

    fort_gl_slab(time_array.dataPtr(),&StateData_slab_dt_type,
                 &bfact_time_order,&slablow,&slabhigh);

    if (do_scale_time==1) {
     for (int islab=0;islab<=bfact_time_order;islab++)
      time_array[islab]=time_array[islab]*dt;
    } else if (do_scale_time==0) {
     // do nothing
    } else {
     amrex::Error("do_scale_time invalid");
    }

    for (int i=0;i<=bfact_time_order-1;i++) {
     if (time_array[i]>time_array[i+1]) {
      std::cout << "i= " << i << " time_array " << time_array[i];
      std::cout << "i+1= " << i+1 << " time_array " << time_array[i+1];
      amrex::Error("time_array corruption");
     }
    }
 
    int ncomp = desc->nComp();
    int state_holds_data = desc->get_state_holds_data();

    for (int i=0;i<=bfact_time_order;i++) {

     if (state_holds_data==1) {
      new_data[i]=new MultiFab(grids,dmap,ncomp,desc->nExtra(),
       MFInfo().SetTag("new_data"),FArrayBoxFactory());
     } else if (state_holds_data==0) {
      new_data[i]=nullptr;
     } else
      amrex::Error("state_holds_data invalid");
    }// for (int i=0;i<=bfact_time_order;i++) 

    buildBC();
} // end subroutine StateData::define 

// this constructs a variable of type "StateData" from checkpoint data.
void
StateData::restart (
  AmrCore& papa,
  int time_order,
  int slab_dt_type, // 0=SEM 1=evenly spaced
  int MAX_NUM_SLAB,
  int level,
  int max_level,
  std::istream& is,
  const Box& p_domain,
  const BoxArray&        grds,
  const DistributionMapping& dm,
  const StateDescriptor& d,
  const StateDescriptor& dGHOST,
  const std::string&     chkfile)
{
    AmrCore* parent=&papa;

    if (time_order==parent->Time_blockingFactor()) {
     //do nothing
    } else
     amrex::Error("expecting time_order==parent->Time_blockingFactor()");

    if (max_level>=0) {
     // do nothing
    } else
     amrex::Error("max_level>=0 violated");

    StateData_MAX_NUM_SLAB=MAX_NUM_SLAB;
    if (StateData_MAX_NUM_SLAB<33)
     amrex::Error("StateData_MAX_NUM_SLAB too small");

    StateData_slab_dt_type=slab_dt_type;
    if ((StateData_slab_dt_type!=0)&&
        (StateData_slab_dt_type!=1))
     amrex::Error("StateData_slab_dt_type invalid");

    time_array.resize(StateData_MAX_NUM_SLAB);
    new_data.resize(StateData_MAX_NUM_SLAB);

    bfact_time_order=time_order;

    if ((bfact_time_order<1)||
        (bfact_time_order>StateData_MAX_NUM_SLAB)) {
     std::cout << "bfact_time_order= " << bfact_time_order << '\n';
     amrex::Error("bfact_time_order invalid in restart");
    }

    StateData_level=level;

    desc = &d;
    descGHOST = &dGHOST;

    domain = p_domain;
    grids = grds;
    dmap = dm;

    // Convert to proper type.
    IndexType typ(desc->getType());
    if (!typ.cellCentered()) {
        domain.convert(typ);
        grids.convert(typ);
    }

    Box domain_in;
    BoxArray grids_in;

    is >> domain_in;
    grids_in.readFrom(is);

    if (domain_in!=domain) {
     std::cout << "domain: " << domain;
     std::cout << "domain_in: " << domain_in;
     amrex::Error("domain_in invalid");
    }

    if (! amrex::match(grids_in,grids))
     amrex::Error("grids_in invalid");     

    for (int i=0;i<=bfact_time_order;i++) {
     is >> time_array[i];
    }

    int nsets;
    is >> nsets;
    if (nsets!=bfact_time_order)
     amrex::Error("all slab data should be checkpointed");

    for (int i=0;i<StateData_MAX_NUM_SLAB;i++)
     new_data[i]=nullptr;

    std::string mf_name;
    std::string FullPathName;

    int state_holds_data = desc->get_state_holds_data();

    for (int i=0;i<=bfact_time_order;i++) {

     if (state_holds_data==1) {
      new_data[i]=new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(),
        MFInfo().SetTag("new_data"),FArrayBoxFactory());
     } else if (state_holds_data==0) {
      new_data[i]=nullptr;
     } else
      amrex::Error("state_holds_data invalid");


       // read the file name from the header file.
     is >> mf_name;
      //
      // Note that mf_name is relative to the Header file:
      // mf_name=Level_<level num>/SD_<state index>_New_MF<slab index>
      //
      // We need to prepend the name of the chkfile directory.
      //
      // e.g. chkfile=./chk<nsteps>
      //
     if (state_holds_data==1) {

      FullPathName = chkfile;
      if (!chkfile.empty() && chkfile[chkfile.length()-1] != '/')
       FullPathName += '/';
      FullPathName += mf_name;
       // read from a file other than the header file.
      VisMF::Read(*new_data[i], FullPathName);

     } else if (state_holds_data==0) {
      new_data[i]=nullptr;
     } else
      amrex::Error("state_holds_data invalid");

    }  // i=0 ... bfact_time_order

    buildBC();
}

void
StateData::buildBC ()
{
    int ncomp = desc->nComp();
    bc.resize(ncomp);
    for (int i = 0; i < ncomp; i++)
    {
        bc[i].resize(grids.size());
        for (int j = 0; j < grids.size(); j++)
        {
            BCRec bcr;
            amrex::setBC(grids[j],domain,desc->getBC(i),bcr);
            bc[i][j]=bcr;
        }
    }

    int ncompGHOST = descGHOST->nComp();
    bcGHOST.resize(ncompGHOST);
    for (int i = 0; i < ncompGHOST; i++)
    {
        bcGHOST[i].resize(grids.size());
        for (int j = 0; j < grids.size(); j++)
        {
            BCRec bcr;
            amrex::setBC(grids[j],domain,descGHOST->getBC(i),bcr);
            bcGHOST[i][j]=bcr;
        }
    }
}

StateData::~StateData() {

 desc = 0;
 descGHOST = 0;
 for (int i=0;i<=bfact_time_order;i++) {

  int state_holds_data=1;
  if (new_data[i]==nullptr)
   state_holds_data=0;

  if (state_holds_data==1) {
   delete new_data[i];
  } else if (state_holds_data==0) {
   // do nothing
  } else
   amrex::Error("state_holds_data invalid");

 } // i=0..bfact_time_order
} // end subroutine StateData::~StateData() 

const StateDescriptor*
StateData::descriptor () const
{
    return desc;
}


const StateDescriptor*
StateData::descriptorGHOST () const
{
    return descGHOST;
}


int StateData::get_bfact_time_order() const
{
 return bfact_time_order;
}

const Box&
StateData::getDomain () const
{
    return domain;
}

Real
StateData::slabTime (int slab_index) const
{
 if ((slab_index<0)||(slab_index>bfact_time_order))
  amrex::Error("slab_index invalid");

 return time_array[slab_index];
}


MultiFab&
StateData::newData (int slab_index)
{
 int project_slab_index=slab_index;
 if (project_slab_index==-1)
  project_slab_index=0;
 if (project_slab_index==bfact_time_order+1)
  project_slab_index=bfact_time_order;
 if ((project_slab_index<0)||
     (project_slab_index>bfact_time_order)) {
  std::cout << "bfact_time_order= " << bfact_time_order << '\n';
  std::cout << "project_slab_index= " << project_slab_index << '\n';
  amrex::Error("project_slab_index invalid1");
 }
 if (new_data[project_slab_index] == nullptr)
  amrex::Error("new_data[project_slab_index]==nullptr");

 return *new_data[project_slab_index];
}

const MultiFab&
StateData::newData (int slab_index) const
{
 int project_slab_index=slab_index;
 if (project_slab_index==-1)
  project_slab_index=0;
 if (project_slab_index==bfact_time_order+1)
  project_slab_index=bfact_time_order;
 if ((project_slab_index<0)||
     (project_slab_index>bfact_time_order)) {
  std::cout << "bfact_time_order= " << bfact_time_order << '\n';
  std::cout << "project_slab_index= " << project_slab_index << '\n';
  amrex::Error("project_slab_index invalid2");
 }
 if (new_data[project_slab_index] == nullptr)
  amrex::Error("new_data[project_slab_index]==nullptr");

 return *new_data[project_slab_index];
}

Vector<BCRec>&
StateData::getBCs (int comp)
{
    return bc[comp];
}

const BCRec&
StateData::getBC (int comp, int i) const
{
    return bc[comp][i];
}


Vector<BCRec>&
StateData::getBCsGHOST (int comp)
{
    return bcGHOST[comp];
}

const BCRec&
StateData::getBCGHOST (int comp, int i) const
{
    return bcGHOST[comp][i];
}


// dt will be modified if t_new>0 and t_new-dt<0
void
StateData::setTimeLevel (Real time,Real& dt)
{
  //time_array[0]=time-dt
  //time_array[bfact_time_order]=time

 if (time==0.0) {
  // do nothing
 } else if ((time>0.0)&&(time-dt<0.0)) {
  dt=time;
 } else if ((time>0.0)&&(time-dt>=0.0)) {
  // do nothing
 } else
  amrex::Error("time or dt bust");

 if (dt<1.0e-12) {
  std::cout << "dt= " << dt << '\n';
  amrex::Error("dt<1e-12 in setTimeLevel StateData");
 }
 if (dt<1.0e-99) {
  std::cout << "dt= " << dt << '\n';
  amrex::Error("dt<1e-99 in setTimeLevel StateData");
 }
 if (time<0.0)
  amrex::Error("time<0 in setTimeLevel StateData");

 if ((StateData_slab_dt_type!=0)&&
     (StateData_slab_dt_type!=1)) 
  amrex::Error("StateData_slab_dt_type invalid");

 Real slablow=time-dt;
 Real slabhigh=time;
 int do_scale_time=0;

 if (dt>=1.0) {
  do_scale_time=1;
  slablow=time/dt-1.0;
  slabhigh=time/dt;
 } else if ((dt>0.0)&&(dt<1.0)) {
  // do nothing
 } else {
  amrex::Error("dt invalid");
 }

 fort_gl_slab(time_array.dataPtr(),&StateData_slab_dt_type,
              &bfact_time_order,&slablow,&slabhigh);

 if (do_scale_time==1) {
  for (int islab=0;islab<=bfact_time_order;islab++)
   time_array[islab]=time_array[islab]*dt;
 } else if (do_scale_time==0) {
  // do nothing
 } else {
  amrex::Error("do_scale_time invalid");
 }

 if (do_scale_time==1) {
  if (std::abs(time_array[0]/dt-slablow)>1.0e-10)
   amrex::Error("time_array[0] inv in setTimeLevel StateData");
  if (std::abs(time_array[bfact_time_order]/dt-slabhigh)>1.0e-10)
   amrex::Error("time_array[bfact_time_order] inv setTimeLevel StateData");
 } else if (do_scale_time==0) {
  if (std::abs(time_array[0]-slablow)>1.0e-10*dt)
   amrex::Error("time_array[0] inv in setTimeLevel StateData");
  if (std::abs(time_array[bfact_time_order]-slabhigh)>1.0e-10*dt)
   amrex::Error("time_array[bfact_time_order] inv setTimeLevel StateData");
 } else {
  amrex::Error("dt invalid");
 }
  
} // subroutine setTimeLevel

void
StateData::get_grid_type(IndexType local_typ,int& grid_type) {

 grid_type=-1;

  // CELL - CELL - CELL
 if (local_typ.cellCentered()) {
  grid_type=-1;
  // CELL - CELL - CELL
 } else if (local_typ==TheXUMACType) {
  grid_type=-1;
  amrex::Error("expecting TheXUMACType==cellCentered");
  // CELL - CELL - CELL
 } else if (local_typ==TheYVMACType) {
  grid_type=-1;
  amrex::Error("expecting TheYVMACType==cellCentered");
  // CELL - CELL - CELL
 } else if ((local_typ==TheZWMACType)&&
            (AMREX_SPACEDIM==3)) {
  grid_type=-1;
  amrex::Error("expecting TheZWMACType==cellCentered");
  // NODE - CELL - CELL
 } else if (local_typ==TheUMACType) {
  grid_type=0;
  // CELL - NODE - CELL
 } else if (local_typ==TheVMACType) {
  grid_type=1;
  // CELL - CELL - NODE 
 } else if ((local_typ==TheWMACType)&&
            (AMREX_SPACEDIM==3)) {
  grid_type=2;
  // NODE - NODE - CELL
 } else if (local_typ==TheYUMACType) {
  grid_type=3;
  // NODE - NODE - CELL
 } else if (local_typ==TheXVMACType) {
  grid_type=3;
  amrex::Error("expecting TheXVMACType==TheYUMACType");
  // NODE - CELL - NODE 
 } else if ((local_typ==TheZUMACType)&&
            (AMREX_SPACEDIM==3)) {
  grid_type=4;
  // NODE - CELL - NODE 
 } else if ((local_typ==TheXWMACType)&&
            (AMREX_SPACEDIM==3)) {
  grid_type=4;
  amrex::Error("expecting TheZUMACType==TheXWMACType");
 } else if ((local_typ==TheZVMACType)&&
            (AMREX_SPACEDIM==3)) {
  grid_type=5;
  // CELL - NODE - NODE 
 } else if ((local_typ==TheYWMACType)&&
            (AMREX_SPACEDIM==3)) {
  grid_type=5;
  amrex::Error("expecting TheZVMACType==TheYWMACType");
 } else
  amrex::Error("grid_type not supported");
} // end subroutine StateData::get_grid_type()

void
StateData::FillBoundary (
 int level,
 FArrayBox& dest,
 Real time,
 const Real* dx,
 const RealBox& prob_domain,
 int dcomp,
 Vector<int> scompBC_map,
 int ncomp,
 int bfact)
{

    IndexType local_typ(desc->getType());

    if (dest.box().ixType()==local_typ) {
     // do nothing
    } else
     amrex::Error("(desc) dest.box().ixType()!=local_typ");

    int grid_type=-1;
    get_grid_type(local_typ,grid_type);

    if (domain.contains(dest.box())) return;

    const Box& bx  = dest.box();
    const int* dlo = dest.loVect();
    const int* dhi = dest.hiVect();
    const int* plo = domain.loVect();
    const int* phi = domain.hiVect();

    Vector<int> bcrs;

    Real xlo[AMREX_SPACEDIM];
    BCRec bcr;
    const Real* problo = prob_domain.lo();

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(dlo[i]-plo[i]);
    }
    int dc_offset=0;
    for (int i = 0; i < ncomp; )
    {
        const int dc  = dcomp+dc_offset;
        const int sc  = scompBC_map[i];
        Real*     dat = dest.dataPtr(dc);

        if (desc->master(sc))
        {
            int groupsize = desc->groupsize(sc);

            BL_ASSERT(groupsize != 0);

            if (groupsize+i <= ncomp)
            {
                //
                // Can do the whole group at once.
                //
                bcrs.resize(2*AMREX_SPACEDIM*groupsize);
                int* bci  = bcrs.dataPtr();

		int kbase=0;

                for (int j = 0; j < groupsize; j++)
                {
                    dc_offset+=1;

                    amrex::setBC(bx,domain,desc->getBC(sc+j),bcr);

                    const int* local_bc = bcr.vect();

                    for (int k = 0; k < 2*AMREX_SPACEDIM; k++)
                        bci[kbase+k] = local_bc[k];

                    kbase += 2*AMREX_SPACEDIM;
                }
		if (kbase==2*AMREX_SPACEDIM*groupsize) {
		 // do nothing
		} else
		 amrex::Error("kbase invalid");

                //
                // Use the "group" boundary fill routine.
                //
                desc->bndryFill(sc)(
                  &grid_type,
                  &level,
                  dat,dlo,dhi,plo,phi,dx,xlo,
                  &time,bcrs.dataPtr(),&sc,&groupsize,&bfact,true);

                i += groupsize;
            }
            else
            {
                int single_ncomp=1;
                dc_offset+=1;

                amrex::setBC(bx,domain,desc->getBC(sc),bcr);
                desc->bndryFill(sc)( 
                  &grid_type,
                  &level,
                  dat,dlo,dhi,plo,phi,dx,xlo,
                  &time,bcr.vect(),&sc,&single_ncomp,&bfact);
                i++;
            }
        }
        else
        {
            int single_ncomp=1;
            dc_offset+=1;

            amrex::setBC(bx,domain,desc->getBC(sc),bcr);
            desc->bndryFill(sc)(
              &grid_type,
              &level,
              dat,dlo,dhi,plo,phi,dx,xlo,
              &time,bcr.vect(),&sc,&single_ncomp,&bfact);
            i++;
        }
    }
} //  StateData::FillBoundary 

void
StateData::FillBoundaryGHOST (
 int level,
 FArrayBox& dest,
 Real time,
 const Real* dx,
 const RealBox& prob_domain,
 int dcomp,
 Vector<int> scompBC_map,
 int ncomp,
 int bfact)
{
    IndexType local_typ(descGHOST->getType());

    if (dest.box().ixType()==local_typ) {
     // do nothing
    } else
     amrex::Error("(descGHOST) dest.box().ixType()!=local_typ");

    int grid_type=-1;
    get_grid_type(local_typ,grid_type);

    if (domain.contains(dest.box())) return;

    const Box& bx  = dest.box();
    const int* dlo = dest.loVect();
    const int* dhi = dest.hiVect();
    const int* plo = domain.loVect();
    const int* phi = domain.hiVect();

    Vector<int> bcrs;

    Real xlo[AMREX_SPACEDIM];
    BCRec bcr;
    const Real* problo = prob_domain.lo();

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(dlo[i]-plo[i]);
    }
    int dc_offset=0;
    for (int i = 0; i < ncomp; )
    {
        const int dc  = dcomp+dc_offset;
        const int sc  = scompBC_map[i];
        Real*     dat = dest.dataPtr(dc);

        if (descGHOST->master(sc))
        {
            int groupsize = descGHOST->groupsize(sc);

            BL_ASSERT(groupsize != 0);

            if (groupsize+i <= ncomp)
            {
                //
                // Can do the whole group at once.
                //
                bcrs.resize(2*AMREX_SPACEDIM*groupsize);
                int* bci  = bcrs.dataPtr();

		int kbase=0;

                for (int j = 0; j < groupsize; j++)
                {
                    dc_offset+=1;

                    amrex::setBC(bx,domain,descGHOST->getBC(sc+j),bcr);

                    const int* local_bc = bcr.vect();

                    for (int k = 0; k < 2*AMREX_SPACEDIM; k++)
                        bci[kbase+k] = local_bc[k];

                    kbase += 2*AMREX_SPACEDIM;
                }
		if (kbase==2*AMREX_SPACEDIM*groupsize) {
		 // do nothing
		} else
		 amrex::Error("kbase invalid");

                //
                // Use the "group" boundary fill routine.
                //
                descGHOST->bndryFill(sc)(
                  &grid_type,
                  &level,
                  dat,dlo,dhi,plo,phi,dx,xlo,
                  &time,bcrs.dataPtr(),&sc,&groupsize,&bfact,true);

                i += groupsize;
            }
            else
            {
                int single_ncomp=1;
                dc_offset+=1;

                amrex::setBC(bx,domain,descGHOST->getBC(sc),bcr);
                descGHOST->bndryFill(sc)(
                  &grid_type,
                  &level,
                  dat,dlo,dhi,plo,phi,dx,xlo,
                  &time,bcr.vect(),&sc,&single_ncomp,&bfact);
                i++;
            }
        }
        else
        {
            int single_ncomp=1;
            dc_offset+=1;

            amrex::setBC(bx,domain,descGHOST->getBC(sc),bcr);
            descGHOST->bndryFill(sc)(
              &grid_type,
              &level,
              dat,dlo,dhi,plo,phi,dx,xlo,
              &time,bcr.vect(),&sc,&single_ncomp,&bfact);
            i++;
        }
    }
} //  StateData::FillBoundaryGHOST 


void
StateData::get_time_index(Real time,Real &nudge_time,int& best_index) {

 nudge_time=time;
 if (bfact_time_order<1)
  amrex::Error("bfact_time_order invalid");

 if ((time_array[0]<=0.0)&&(time_array[bfact_time_order]==0.0)) {
  best_index=bfact_time_order;
  nudge_time=0.0;
 } else if ((time_array[0]>=0.0)&&(time_array[bfact_time_order]>0.0)) {

  Real time_scale_begin=time_array[0]/time_array[bfact_time_order];
  Real time_scale=time_scale_begin;

  Real eps=(1.0-time_scale_begin)*1.0e-10;
  if (eps<=0.0) {
   std::cout << "time_scale_begin= " << time_scale_begin << '\n';
   std::cout << "time_array(0)= " << time_array[0] << '\n';
   std::cout << "time_array(bfact_time_order)= " << 
    time_array[bfact_time_order] << '\n';
   std::cout << "bfact_time_order= " << bfact_time_order << '\n';
   std::cout << "time= " << time << '\n';
   std::cout << "desc->nComp() " << desc->nComp() << '\n';

   for (int i=0;i<desc->nComp();i++)
    std::cout << "i= " << i << " name= " << desc->name(i) << '\n';

   std::cout << "desc->nExtra() " << desc->nExtra() << '\n';

   amrex::Error("eps invalid (a)");
  }
  time_scale=time/time_array[bfact_time_order];
  if (time_scale<time_scale_begin-eps)
   amrex::Error("time_scale too small");
  if (time_scale>1.0+eps)
   amrex::Error("time_scale too big");

  best_index=bfact_time_order;
  nudge_time=time_array[bfact_time_order];

  for (int slab_step=bfact_time_order-1;slab_step>=0;slab_step--) {
   Real time_scale_current=time_array[slab_step]/time_array[bfact_time_order];
   Real time_scale_best=time_array[best_index]/time_array[bfact_time_order];
   if (std::abs(time_scale-time_scale_current)<
       std::abs(time_scale-time_scale_best)) {
    best_index=slab_step;
    nudge_time=time_array[slab_step];
   }
  } // slabstep

 } else {
  std::cout << "bfact_time_order= " << bfact_time_order << '\n';
  std::cout << "time= " << time << '\n';
  for (int itime=0;itime<=bfact_time_order;itime++) {
   std::cout << "itime= " << itime << "time_array[itime]= " <<
    time_array[itime] << '\n';
  }
  amrex::Error("time_array bust");
 }

} // get_time_index

void
StateData::get_time_bounding_box(Real time,Real &nudge_time,
  int &start_index) {

 nudge_time=time;
 if (bfact_time_order<1)
  amrex::Error("bfact_time_order invalid");

 if ((time_array[0]<=0.0)&&(time_array[bfact_time_order]==0.0)) {
  start_index=bfact_time_order-1;
  nudge_time=0.0;
 } else if ((time_array[0]>=0.0)&&(time_array[bfact_time_order]>0.0)) {

  Real time_scale_begin=time_array[0]/time_array[bfact_time_order];
  Real time_scale=time_scale_begin;

  Real eps=(1.0-time_scale_begin)*1.0e-10;
  if (eps<=0.0) {
   std::cout << "time_scale_begin= " << time_scale_begin << '\n';
   std::cout << "time_array(0)= " << time_array[0] << '\n';
   std::cout << "time_array(bfact_time_order)= " << 
    time_array[bfact_time_order] << '\n';
   std::cout << "bfact_time_order= " << bfact_time_order << '\n';
   std::cout << "time= " << time << '\n';
   std::cout << "desc->nComp() " << desc->nComp() << '\n';

   for (int i=0;i<desc->nComp();i++)
    std::cout << "i= " << i << " name= " << desc->name(i) << '\n';

   std::cout << "desc->nExtra() " << desc->nExtra() << '\n';

   amrex::Error("eps invalid (b)");
  }

  time_scale=time/time_array[bfact_time_order];
  if (time_scale<time_scale_begin-eps)
   amrex::Error("time_scale too small");
  if (time_scale>1.0+eps)
   amrex::Error("time_scale too big");

  Real time_scale_current=time_array[1]/time_array[bfact_time_order];
  Real time_scale_current_alt=
   time_array[bfact_time_order-1]/time_array[bfact_time_order];

  if (time_scale<=time_scale_current) {
   if (time_scale<time_scale_begin)
    nudge_time=time_array[0];
   start_index=0;
  } else if (time_scale>=time_scale_current_alt) {
   if (time_scale>1.0)
    nudge_time=time_array[bfact_time_order];
   start_index=bfact_time_order-1;
  } else {
   start_index=0;
   Real time_scale_loop= 
	  time_array[start_index+1]/time_array[bfact_time_order];
   while (time_scale>time_scale_loop) {
    start_index++;
    if (start_index>=bfact_time_order)
     amrex::Error("start_index invalid");
    time_scale_loop= 
      time_array[start_index+1]/time_array[bfact_time_order];
   }
  }
 }

} // get_time_bounding_box

void
StateData::CopyNewToOld(int level,int max_level) {

 if (max_level>=0) {
  // do nothing
 } else
  amrex::Error("max_level>=0 violated");

 if (level>=0) {
  // do nothing
 } else
  amrex::Error("level>=0 violated");

 int state_holds_data = desc->get_state_holds_data();

 if (state_holds_data==1) {

  MultiFab & newmulti = *new_data[bfact_time_order];
  int ncomp=newmulti.nComp();
  int ngrow=newmulti.nGrow();

  for (int i=0;i<bfact_time_order;i++) {
   MultiFab & oldmulti = *new_data[i];
   MultiFab::Copy(oldmulti,newmulti,0,0,ncomp,ngrow);
  }
 } else if (state_holds_data==0) {
  for (int i=0;i<bfact_time_order;i++) {
   if (new_data[i]!=nullptr)
    amrex::Error("new_data[i]!=nullptr");
  }
 } else
  amrex::Error("state_holds_data invalid");

} // end subroutine CopyNewToOld

void
StateData::CopyOldToNew(int level,int max_level) {

 if (max_level>=0) {
  // do nothing
 } else
  amrex::Error("max_level>=0 violated");

 if (level>=0) {
  // do nothing
 } else
  amrex::Error("level>=0 violated");

 int state_holds_data = desc->get_state_holds_data();

 if (state_holds_data==1) {

  MultiFab & oldmulti = *new_data[0];

  int ncomp=oldmulti.nComp();
  int ngrow=oldmulti.nGrow();

  for (int i=1;i<=bfact_time_order;i++) {
   MultiFab & newmulti = *new_data[i];
   MultiFab::Copy(newmulti,oldmulti,0,0,ncomp,ngrow);
  }

 } else if (state_holds_data==0) {
  for (int i=1;i<=bfact_time_order;i++) {
   if (new_data[i]!=nullptr)
    amrex::Error("new_data[i]!=nullptr");
  }
 } else
  amrex::Error("state_holds_data invalid");

} // end subroutine CopyOldToNew

void
StateData::checkPoint (const std::string& name,
                       const std::string& fullpathname,
                       std::ostream&  os,
		       int level,int max_level)
{

    if (level>=0) {
     // do nothing
    } else
     amrex::Error("level invalid");

    if (max_level>=0) {
     // do nothing
    } else
     amrex::Error("max_level invalid");

     // on input:  
     //   name=Level_<level num>/SD_<state index>
     //
     // mf_name=Level_<level num>/SD_<state index>_New_MF<slab index>
    Vector<std::string> NewSuffix;
    NewSuffix.resize(StateData_MAX_NUM_SLAB);
    Vector<std::string> mf_name;
    mf_name.resize(StateData_MAX_NUM_SLAB);
    for (int i=0;i<=bfact_time_order;i++) {
     std::stringstream slab_string_stream(std::stringstream::in |
      std::stringstream::out);
     slab_string_stream << i;
     std::string slab_string=slab_string_stream.str();
     NewSuffix[i]="_New_MF" + slab_string;
     mf_name[i]=name;
     mf_name[i]+=NewSuffix[i];
    }

    if (ParallelDescriptor::IOProcessor()) {

     os << domain << '\n';

     grids.writeOn(os);

     for (int i=0;i<=bfact_time_order;i++)  {
       os << time_array[i] << '\n';
     }

         
      // output to the header file:
     os << bfact_time_order << '\n';
     for (int i=0;i<=bfact_time_order;i++) {
       os << mf_name[i] << '\n';
     }  // i

    }  // IOProcessor ?

    int state_holds_data = desc->get_state_holds_data();

    for (int i=0;i<=bfact_time_order;i++) {
     if (state_holds_data==1) {

      if (new_data[i]==nullptr)
       amrex::Error("new_data not allocated");

       // mf_fullpath_new=./chk<nsteps>/Level_<level num>/SD_<state index>
       //   _New_MF<slab index> 
       // fullpathname=./chk<nsteps>/Level_<level num>/SD_<state index>
      std::string mf_fullpath_new = fullpathname; 
      mf_fullpath_new += NewSuffix[i];
      // this file is not the header file
      VisMF::Write(*new_data[i],mf_fullpath_new);

     } else if (state_holds_data==0) {
      if (new_data[i]!=nullptr)
       amrex::Error("new_data[i]!=nullptr");
     } else
      amrex::Error("state_holds_data invalid");

    }  // i=0..bfact_time_order

} // end subroutine StateData::checkPoint

void
StateData::printTimeInterval (std::ostream &os) const
{
    os << '['
       << time_array[0]
       << "] ["
       << time_array[bfact_time_order]
       << ']'
       << '\n';
}

StateDataPhysBCFunct::StateDataPhysBCFunct (
 StateData& sd, 
 const Geometry& geom_)
    : statedata(&sd),
      geom(geom_)
{ }

// dx,prob_domain come from variables that are local to this class.
void
StateDataPhysBCFunct::FillBoundary (
 int level,
 MultiFab& mf, 
 Real time,
 int dcomp, 
 Vector<int> scompBC_map,
 int ncomp, 
 int bfact)
{
 BL_PROFILE("StateDataPhysBCFunct::FillBoundary");

 const Box&     domain      = statedata->getDomain();
 const int*     domainlo    = domain.loVect();
 const int*     domainhi    = domain.hiVect();
 const Real*    dx          = geom.CellSize();
 const RealBox& prob_domain = geom.ProbDomain();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mf.boxArray().d_numPts());

#ifdef CRSEGRNDOMP
#ifdef _OPENMP
#pragma omp parallel
#endif
#endif
{
  for (MFIter mfi(mf,false); mfi.isValid(); ++mfi) {
   const Box& tilegrid=mfi.tilebox();

   int tid_current=0;
#ifdef _OPENMP
   tid_current = omp_get_thread_num();
#endif
   if ((tid_current>=0)&&(tid_current<thread_class::nthreads)) {
    // do nothing
   } else
    amrex::Error("tid_current invalid");

   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FArrayBox& dest = mf[mfi];
   const Box& bx = dest.box();
           
   bool has_phys_bc = false;
   bool is_periodic = false;
   for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    bool touch = bx.smallEnd(i) < domainlo[i] || bx.bigEnd(i) > domainhi[i];
    if (geom.isPeriodic(i)) {
     is_periodic = is_periodic || touch;
    } else {
     has_phys_bc = has_phys_bc || touch;
    }
   } // i
           
   if (has_phys_bc) {
    statedata->FillBoundary(
     level,
     dest, 
     time, 
     dx, 
     prob_domain, 
     dcomp, 
     scompBC_map, 
     ncomp,
     bfact);
       	
    if (is_periodic) { // fix corner
     Box GrownDomain = domain;
     for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      if (!geom.isPeriodic(dir)) {
       const int lo = domainlo[dir] - bx.smallEnd(dir);
       const int hi = bx.bigEnd(dir) - domainhi[dir];
       if (lo > 0) GrownDomain.growLo(dir,lo);
       if (hi > 0) GrownDomain.growHi(dir,hi);
      }
     }
        	    
     for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      if (!geom.isPeriodic(dir)) continue;
        		
      Box lo_slab = bx;
      Box hi_slab = bx;
      lo_slab.shift(dir, domain.length(dir));
      hi_slab.shift(dir,-domain.length(dir));
      lo_slab &= GrownDomain;
      hi_slab &= GrownDomain;
      if (lo_slab.ok()) {
       lo_slab.shift(dir,-domain.length(dir));
       FArrayBox tmp;
       tmp.resize(lo_slab,ncomp);

       Array4<Real> const& dest_array=dest.array();
       Array4<Real> const& tmp_array=tmp.array();

       Box dest_box=dest.box();
       Box tmp_box=tmp.box();
       dest_box &= tmp_box;

       const Dim3 lo3a=amrex::lbound(dest_box);
       const Dim3 hi3a=amrex::ubound(dest_box);
       for (int n=0;n<ncomp;++n) {
       for (int z=lo3a.z;z<=hi3a.z;++z) {
       for (int y=lo3a.y;y<=hi3a.y;++y) {
       for (int x=lo3a.x;x<=hi3a.x;++x) {
        tmp_array(x,y,z,n)=dest_array(x,y,z,n+dcomp);
       }
       }
       }
       }

       tmp.shift(dir,domain.length(dir));

       int dcomp_tmp=0;

       statedata->FillBoundary(
        level,
        tmp, 
        time, 
        dx, 
        prob_domain, 
        dcomp_tmp, 
        scompBC_map, 
        ncomp,
        bfact);
        		    
       tmp.shift(dir,-domain.length(dir));

       dest_box=dest.box();
       tmp_box=tmp.box();
       dest_box &= tmp_box;

       const Dim3 lo3b=amrex::lbound(dest_box);
       const Dim3 hi3b=amrex::ubound(dest_box);
       for (int n=0;n<ncomp;++n) {
       for (int z=lo3b.z;z<=hi3b.z;++z) {
       for (int y=lo3b.y;y<=hi3b.y;++y) {
       for (int x=lo3b.x;x<=hi3b.x;++x) {
        dest_array(x,y,z,n+dcomp)=tmp_array(x,y,z,n);
       }
       }
       }
       }

      } // lo_slab
      if (hi_slab.ok()) {
       hi_slab.shift(dir,domain.length(dir));
       FArrayBox tmp;
       tmp.resize(hi_slab,ncomp);

       Array4<Real> const& dest_array=dest.array();
       Array4<Real> const& tmp_array=tmp.array();

       Box dest_box=dest.box();
       Box tmp_box=tmp.box();
       dest_box &= tmp_box;

       const Dim3 lo3a=amrex::lbound(dest_box);
       const Dim3 hi3a=amrex::ubound(dest_box);
       for (int n=0;n<ncomp;++n) {
       for (int z=lo3a.z;z<=hi3a.z;++z) {
       for (int y=lo3a.y;y<=hi3a.y;++y) {
       for (int x=lo3a.x;x<=hi3a.x;++x) {
        tmp_array(x,y,z,n)=dest_array(x,y,z,n+dcomp);
       }
       }
       }
       }

       tmp.shift(dir,-domain.length(dir));

       int dcomp_tmp=0;

       statedata->FillBoundary(
        level,
        tmp, 
        time, 
        dx, 
        prob_domain, 
        dcomp_tmp, 
        scompBC_map, 
        ncomp,
        bfact);
        		    
       tmp.shift(dir,domain.length(dir));

       dest_box=dest.box();
       tmp_box=tmp.box();
       dest_box &= tmp_box;

       const Dim3 lo3b=amrex::lbound(dest_box);
       const Dim3 hi3b=amrex::ubound(dest_box);
       for (int n=0;n<ncomp;++n) {
       for (int z=lo3b.z;z<=hi3b.z;++z) {
       for (int y=lo3b.y;y<=hi3b.y;++y) {
       for (int x=lo3b.x;x<=hi3b.x;++x) {
        dest_array(x,y,z,n+dcomp)=tmp_array(x,y,z,n);
       }
       }
       }
       }

      } // hi_slab
     } // dir
    } // periodic?
   } // has_phys_bc?
  } // mfi
} // omp
  thread_class::sync_tile_d_numPts();
  ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
  thread_class::reconcile_d_numPts(LOOP_FILLBOUNDARY,"FillBoundary");

} // StateDataPhysBCFunct::FillBoundary


StateDataPhysBCFunctGHOST::StateDataPhysBCFunctGHOST (
 StateData& sd, 
 const Geometry& geom_)
    : statedata(&sd),
      geom(geom_)
{ }

// dx,prob_domain come from variables that are local to this class.
void
StateDataPhysBCFunctGHOST::FillBoundary (
 int level,
 MultiFab& mf, 
 Real time,
 int dcomp, 
 Vector<int> scompBC_map,
 int ncomp, 
 int bfact)
{
 BL_PROFILE("StateDataPhysBCFunctGHOST::FillBoundary");

 const Box&     domain      = statedata->getDomain();
 const int*     domainlo    = domain.loVect();
 const int*     domainhi    = domain.hiVect();
 const Real*    dx          = geom.CellSize();
 const RealBox& prob_domain = geom.ProbDomain();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mf.boxArray().d_numPts());

#ifdef CRSEGRNDOMP
#ifdef _OPENMP
#pragma omp parallel
#endif
#endif
{
  for (MFIter mfi(mf,false); mfi.isValid(); ++mfi) {
   const Box& tilegrid=mfi.tilebox();

   int tid_current=0;
#ifdef _OPENMP
   tid_current = omp_get_thread_num();
#endif
   if ((tid_current>=0)&&(tid_current<thread_class::nthreads)) {
    // do nothing
   } else
    amrex::Error("tid_current invalid");

   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FArrayBox& dest = mf[mfi];
   const Box& bx = dest.box();
           
   bool has_phys_bc = false;
   bool is_periodic = false;
   for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    bool touch=((bx.smallEnd(i)<domainlo[i])||(bx.bigEnd(i)>domainhi[i]));
    if (geom.isPeriodic(i)) {
     is_periodic = is_periodic || touch;
    } else {
     has_phys_bc = has_phys_bc || touch;
    }
   } // i
           
   if (has_phys_bc) {

    statedata->FillBoundaryGHOST(
     level,
     dest, 
     time, 
     dx, 
     prob_domain, 
     dcomp, 
     scompBC_map, 
     ncomp,
     bfact);
       	
    if (is_periodic) { // fix corner
     Box GrownDomain = domain;
     for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      if (!geom.isPeriodic(dir)) {
       const int lo = domainlo[dir] - bx.smallEnd(dir);
       const int hi = bx.bigEnd(dir) - domainhi[dir];
       if (lo > 0) GrownDomain.growLo(dir,lo);
       if (hi > 0) GrownDomain.growHi(dir,hi);
      }
     } // dir
        	    
     for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      if (!geom.isPeriodic(dir)) continue;
        		
      Box lo_slab = bx;
      Box hi_slab = bx;
      lo_slab.shift(dir, domain.length(dir));
      hi_slab.shift(dir,-domain.length(dir));
      lo_slab &= GrownDomain;
      hi_slab &= GrownDomain;
      if (lo_slab.ok()) {
         // lo_slab is in the GrownDomain.
       lo_slab.shift(dir,-domain.length(dir));
       FArrayBox tmp;
       tmp.resize(lo_slab,ncomp);

       Array4<Real> const& dest_array=dest.array();
       Array4<Real> const& tmp_array=tmp.array();

       Box dest_box=dest.box();
       Box tmp_box=tmp.box();
       dest_box &= tmp_box;

       const Dim3 lo3a=amrex::lbound(dest_box);
       const Dim3 hi3a=amrex::ubound(dest_box);
       for (int n=0;n<ncomp;++n) {
       for (int z=lo3a.z;z<=hi3a.z;++z) {
       for (int y=lo3a.y;y<=hi3a.y;++y) {
       for (int x=lo3a.x;x<=hi3a.x;++x) {
        tmp_array(x,y,z,n)=dest_array(x,y,z,n+dcomp);
       }
       }
       }
       }

        // tmp.Box() is outside the Growndomain.
       tmp.shift(dir,domain.length(dir));

       if (1==0) {
        std::cout << "dir= " << dir << '\n';
        std::cout << "ncomp= " << ncomp << '\n';
        std::cout << "dcomp= " << dcomp << '\n';
        std::cout << "GrownDomain= " << GrownDomain << '\n';
        std::cout << "tmp.box() " << tmp.box() << '\n';
       }

       int dcomp_tmp=0;

       statedata->FillBoundaryGHOST(
        level,
        tmp, 
        time, 
        dx, 
        prob_domain, 
        dcomp_tmp, 
        scompBC_map,
        ncomp,
        bfact);
        		   
        // tmp.Box() is inside the Growndomain 
       tmp.shift(dir,-domain.length(dir));

       dest_box=dest.box();
       tmp_box=tmp.box();
       dest_box &= tmp_box;

       const Dim3 lo3b=amrex::lbound(dest_box);
       const Dim3 hi3b=amrex::ubound(dest_box);
       for (int n=0;n<ncomp;++n) {
       for (int z=lo3b.z;z<=hi3b.z;++z) {
       for (int y=lo3b.y;y<=hi3b.y;++y) {
       for (int x=lo3b.x;x<=hi3b.x;++x) {
        dest_array(x,y,z,n+dcomp)=tmp_array(x,y,z,n);
       }
       }
       }
       }

      } // lo_slab

      if (hi_slab.ok()) {
       hi_slab.shift(dir,domain.length(dir));

       FArrayBox tmp;
       tmp.resize(hi_slab,ncomp);

       Array4<Real> const& dest_array=dest.array();
       Array4<Real> const& tmp_array=tmp.array();

       Box dest_box=dest.box();
       Box tmp_box=tmp.box();
       dest_box &= tmp_box;

       const Dim3 lo3a=amrex::lbound(dest_box);
       const Dim3 hi3a=amrex::ubound(dest_box);
       for (int n=0;n<ncomp;++n) {
       for (int z=lo3a.z;z<=hi3a.z;++z) {
       for (int y=lo3a.y;y<=hi3a.y;++y) {
       for (int x=lo3a.x;x<=hi3a.x;++x) {
        tmp_array(x,y,z,n)=dest_array(x,y,z,n+dcomp);
       }
       }
       }
       }

       tmp.shift(dir,-domain.length(dir));

       int dcomp_tmp=0;

       statedata->FillBoundaryGHOST(
        level,
        tmp, 
        time, 
        dx, 
        prob_domain, 
        dcomp_tmp, 
        scompBC_map, 
        ncomp,
        bfact);
        		    
       tmp.shift(dir,domain.length(dir));

       dest_box=dest.box();
       tmp_box=tmp.box();
       dest_box &= tmp_box;

       const Dim3 lo3b=amrex::lbound(dest_box);
       const Dim3 hi3b=amrex::ubound(dest_box);
       for (int n=0;n<ncomp;++n) {
       for (int z=lo3b.z;z<=hi3b.z;++z) {
       for (int y=lo3b.y;y<=hi3b.y;++y) {
       for (int x=lo3b.x;x<=hi3b.x;++x) {
        dest_array(x,y,z,n+dcomp)=tmp_array(x,y,z,n);
       }
       }
       }
       }

      } // hi_slab
     } // dir
    } // periodic?
   } // has_phys_bc?
  } // mfi
} // omp
  thread_class::sync_tile_d_numPts();
  ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
  thread_class::reconcile_d_numPts(LOOP_FILLBOUNDARY_GHOST,"FillBoundary");

} // StateDataPhysBCFunctGHOST::FillBoundary

} // namespace amrex
