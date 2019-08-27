
#include <algorithm>
#include <string>
#include <iostream>

#include <StateDescriptor.H>
#include <Interpolater.H>
#include <AMReX_BCRec.H>

namespace amrex {

StateDescriptor::BndryFunc::BndryFunc ()
    :
    m_func(0),
    m_gfunc(0)
{}

StateDescriptor::BndryFunc::BndryFunc (BndryFuncDefault inFunc)
    :
    m_func(inFunc),
    m_gfunc(0)
{}

StateDescriptor::BndryFunc::BndryFunc (BndryFuncDefault inFunc,
                                       BndryFuncDefault gFunc)
    :
    m_func(inFunc),
    m_gfunc(gFunc)
{}

StateDescriptor::BndryFunc*
StateDescriptor::BndryFunc::clone () const
{
    return new BndryFunc(*this);
}

StateDescriptor::BndryFunc::~BndryFunc () {}

// FILL
void
StateDescriptor::BndryFunc::operator () (
 int* level,
 Real* data,const int* lo,const int* hi,
 const int* dom_lo, const int* dom_hi,
 const Real* dx, const Real* grd_lo,
 const Real* time, const int* bc,
 const int* scomp,int* ncomp,int* bfact) const
{
    BL_ASSERT(m_func != 0);

    m_func(level,
      data,AMReX_ARLIM(lo),AMReX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,bc,
      scomp,ncomp,bfact);
}

// GROUP FILL
void
StateDescriptor::BndryFunc::operator () (
 int* level,
 Real* data,const int* lo,const int* hi,
 const int* dom_lo, const int* dom_hi,
 const Real* dx, const Real* grd_lo,
 const Real* time, const int* bc, 
 const int* scomp,int* ncomp,int* bfact,bool) const
{
    BL_ASSERT(m_gfunc != 0);

    m_gfunc(
      level,
      data,AMReX_ARLIM(lo),AMReX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,bc,
      scomp,ncomp,bfact);
}

// this function cannot throw an exception.
DescriptorList::DescriptorList () noexcept
{}

void
DescriptorList::clear ()
{
    desc.clear();
    desc.resize(0);
}

int
DescriptorList::size () const
{
    return desc.size();
}

void
DescriptorList::saveMapper(int indx,int comp,Interpolater*& interp) {

 desc[indx].saveMapper(comp,interp);
}

void 
DescriptorList::resetMapper(int indx,int comp,Interpolater* interp) {

 desc[indx].resetMapper(comp,interp);
}

void
DescriptorList::save_bcrecs_statedesc(int indx,int comp,BCRec& bcr) {
    
 desc[indx].save_bcrecs_statedesc(comp,bcr);
}   

void
DescriptorList::reset_bcrecs(int indx,int comp,BCRec bcr) {
    
 desc[indx].reset_bcrecs(comp,bcr);
}   




void
DescriptorList::resetComponentBCs (int                               indx,
                                   int                               comp,
                                   const BCRec&                      bc,
                                   const StateDescriptor::BndryFunc& func)
{
    desc[indx].resetComponentBCs(comp,bc,func);
}

void
DescriptorList::setComponent (int        indx,
       int                               comp,
       const std::string&                nm,
       const BCRec&                      bc,
       const StateDescriptor::BndryFunc& func,
       Interpolater*                     interp,
       int                               max_map_start_comp,
       int                               min_map_end_comp)
{
    desc[indx].setComponent(comp,nm,bc,func,interp,max_map_start_comp,
      min_map_end_comp);
}

void
DescriptorList::setComponent (int    indx,
   int                               comp,
   const Array<std::string>&         nm,
   const Array<BCRec>&               bc,
   const StateDescriptor::BndryFunc& func,
   Interpolater*                     interp)
{
    for (int i = 0; i < nm.size(); i++)
    {
        const bool master = (i == 0) ? true : false;

        desc[indx].setComponent(comp+i,nm[i],bc[i],func,interp,
         master,nm.size());
    }
}

const StateDescriptor&
DescriptorList::operator[] (int k) const
{
    return desc[k];
}

void
DescriptorList::addDescriptor (int indx,
      IndexType                   typ,
      int                         nextra,
      int                         num_comp, 
      Interpolater*               interp,
      bool                        store_in_checkpoint)
{
    if (indx >= desc.size())
        desc.resize(indx+1);
    desc.set(indx,new StateDescriptor(typ,indx,nextra,num_comp,
     interp,store_in_checkpoint));
}  

StateDescriptor::StateDescriptor () noexcept
    :
    id(-1),
    ncomp(0),
    ngrow(0),
    mapper(0),
    m_store_in_checkpoint(true)
{}

StateDescriptor::StateDescriptor (IndexType btyp,
        int                         ident,
        int                         nextra, 
        int                         num_comp,
        Interpolater*               interp,
        bool                        store_in_checkpoint)
    :
    type(btyp),
    id(ident),
    ncomp(num_comp),
    ngrow(nextra),
    mapper(interp),
    m_store_in_checkpoint(store_in_checkpoint)
{
    BL_ASSERT (num_comp > 0);
   
    names.resize(num_comp);
    bc.resize(num_comp);
    bc_func.resize(num_comp);
    mapper_comp.resize(num_comp);
    m_master.resize(num_comp);
    m_groupsize.resize(num_comp);
    max_map_start_comp.resize(num_comp);
    min_map_end_comp.resize(num_comp);
}

StateDescriptor::~StateDescriptor ()
{
    mapper = 0;
}

void
StateDescriptor::saveMapper(int comp,Interpolater*& interp) {

 interp=mapper_comp[comp];

}

void
StateDescriptor::resetMapper(int comp,Interpolater* interp) {

 mapper_comp[comp] = interp;

}

void
StateDescriptor::save_bcrecs_statedesc(int comp,BCRec& bcr) {

 bcr=bc[comp];
 
}   

void
StateDescriptor::reset_bcrecs(int comp,BCRec bcr) {
   
 bc[comp] = bcr;

}



void
StateDescriptor::resetComponentBCs (int              comp,
                                    const BCRec&     bcr,
                                    const BndryFunc& func)
{
    BL_ASSERT(comp >= 0 && comp < ncomp);

    bc_func.clear(comp);
    bc_func.set(comp,func.clone());
    bc[comp] = bcr;
}

IndexType
StateDescriptor::getType () const
{
    return type;
}


int
StateDescriptor::nComp () const
{
    return ncomp;
}


int
StateDescriptor::nExtra () const
{
    return ngrow;
}

Interpolater*
StateDescriptor::interp () const
{
    return mapper;
}

Interpolater*
StateDescriptor::interp (int i) const
{
    return mapper_comp[i] == 0 ? mapper : mapper_comp[i];
}

const std::string&
StateDescriptor::name (int i) const
{
    return names[i];
}

const BCRec&
StateDescriptor::getBC (int i) const
{
    return bc[i];
}

const Array<BCRec>&
StateDescriptor::getBCs () const
{
    return bc;
}


bool
StateDescriptor::store_in_checkpoint () const
{
    return m_store_in_checkpoint;
}


const StateDescriptor::BndryFunc&
StateDescriptor::bndryFill (int i) const
{
    return bc_func[i];
}

void
StateDescriptor::check_inRange (Array<int> scompBC_map, int nc) const
{
 if (scompBC_map.size()!=nc)
  amrex::Error("scompBC_map has invalid size");
 
 for (int i=0;i<nc;i++) {
  if ((scompBC_map[i]<0)||(scompBC_map[i]>=ncomp))
   amrex::Error("scompBC_map is corrupt");
 }
}

void
StateDescriptor::define (IndexType btyp,
      int                         ident,
      int                         nextra,
      int                         num_comp,
      Interpolater*               interp,
      bool                        store_in_checkpoint)
{
    type     = btyp;
    id       = ident;
    ngrow    = nextra;
    ncomp    = num_comp;
    mapper   = interp;
    m_store_in_checkpoint = store_in_checkpoint;

    BL_ASSERT (num_comp > 0);
   
    names.resize(num_comp);
    bc.resize(num_comp);
    bc_func.resize(num_comp);
    mapper_comp.resize(num_comp);
    m_master.resize(num_comp);
    m_groupsize.resize(num_comp);
    max_map_start_comp.resize(num_comp);
    min_map_end_comp.resize(num_comp);
}

void
StateDescriptor::setComponent (int    comp,
    const std::string&                nm,
    const BCRec&                      bcr,
    const StateDescriptor::BndryFunc& func,
    Interpolater*                     interp, 
    int                               max_map_start_comp_,
    int                               min_map_end_comp_)
{
    bc_func.clear(comp);
    bc_func.set(comp,func.clone());

    names[comp]       = nm;
    bc[comp]          = bcr;
    mapper_comp[comp] = interp;
    m_master[comp]    = false;
    m_groupsize[comp] = 0;

    if (max_map_start_comp_>=0 && min_map_end_comp_>=0)
    {
        BL_ASSERT(comp >= max_map_start_comp_ &&
                  comp <= min_map_end_comp_   &&
                  min_map_end_comp_ < ncomp);
        max_map_start_comp[comp] = max_map_start_comp_;
        min_map_end_comp[comp]   = min_map_end_comp_;
    }
    else
    {
        max_map_start_comp[comp] = comp;
        min_map_end_comp[comp]   = comp;
    }
}


void
StateDescriptor::setComponent (int                               comp,
                               const std::string&                nm,
                               const BCRec&                      bcr,
                               const StateDescriptor::BndryFunc& func,
                               Interpolater*                     interp,
                               bool                              master,
                               int                               groupsize)
{
    setComponent(comp,nm,bcr,func,interp,-1,-1);

    m_master[comp]    = master;
    m_groupsize[comp] = groupsize;
}

void
StateDescriptor::dumpNames (std::ostream& os,
                            int           start_comp,
                            int           num_comp) const
{
    BL_ASSERT(start_comp >= 0 && start_comp+num_comp <= ncomp);

    for (int k = 0; k < num_comp; k++)
    {
        os << names[start_comp+k] << ' ';
    }
}

void
StateDescriptor::setUpMaps (int&                use_default_map,
                            const Interpolater* default_map,
                            int                 start_comp,
                            int                 num_comp,
                            Interpolater**&     maps, 
                            int&                nmaps,
                            int*&               map_start_comp,
                            int*&               map_num_comp, 
                            int*&               max_start_comp,
                            int*&               min_end_comp) const
{
    BL_ASSERT(start_comp>=0 && start_comp+num_comp-1 < ncomp && num_comp>0);

    maps           = 0;
    map_start_comp = 0;
    map_num_comp   = 0;
    max_start_comp = 0;
    min_end_comp   = 0;
    //
    // First, count number of interpolaters needed and allocate.
    //
    Interpolater* map = mapper_comp[start_comp];
    if (!map) map = (Interpolater*) default_map;
    nmaps = 1; 
    int icomp = start_comp+1;

    use_default_map = 1;
    while (icomp < start_comp+num_comp)
    {
        Interpolater* mapper_icomp = mapper_comp[icomp];
        if (!mapper_icomp)
        {
            mapper_icomp = (Interpolater *) default_map;
        }
        else
        {
            use_default_map = 0;
        }
        if (map != mapper_icomp)
        {
            map = mapper_icomp;
            nmaps++;
        }
        icomp++;
    }

    if (use_default_map) return;

    maps           = new Interpolater*[nmaps];
    map_start_comp = new int[nmaps];
    map_num_comp   = new int[nmaps];
    min_end_comp   = new int[nmaps];
    max_start_comp = new int[nmaps];
    //
    // Now fill the slots.
    //
    int imap             = 0;

    if (mapper_comp[start_comp])
    {
        maps[imap] = mapper_comp[start_comp];
    }
    else
    {
        maps[imap] = (Interpolater *) default_map;
    }

    icomp                = start_comp+1;
    map_start_comp[imap] = start_comp;
    map_num_comp[imap]   = 1;
    min_end_comp[imap]   = min_map_end_comp[start_comp];
    max_start_comp[imap] = max_map_start_comp[start_comp];

    while (icomp < start_comp+num_comp)
    {
        Interpolater* mapper_icomp = mapper_comp[icomp];

        if (!mapper_icomp)
            mapper_icomp = (Interpolater *) default_map;

        if (maps[imap] != mapper_icomp)
        {
            imap++;

            BL_ASSERT (imap < nmaps);

            maps[imap]           = mapper_icomp;
            map_start_comp[imap] = icomp;
            map_num_comp[imap]   = 1;
            min_end_comp[imap]   = min_map_end_comp[icomp];
            max_start_comp[imap] = max_map_start_comp[icomp];

        }
        else
        {
            map_num_comp[imap]++;
            min_end_comp[imap]   = std::max(min_end_comp[imap],min_map_end_comp[icomp]);
            max_start_comp[imap] = std::min(max_start_comp[imap],max_map_start_comp[icomp]);
        }
        icomp++;
    }
}

void
StateDescriptor::cleanUpMaps (Interpolater**& maps, 
                              int*&           map_start_comp,
                              int*&           map_num_comp,
                              int*&           max_start_comp,
                              int*&           min_end_comp) const
{
    delete [] maps;
    delete [] map_start_comp;
    delete [] map_num_comp;
    delete [] max_start_comp;
    delete [] min_end_comp;
}

std::vector< std::pair<int,int> >
StateDescriptor::sameInterps (Array<int> scompBC_map,
                              int ncomp) const
{
    if (ncomp<1)
     amrex::Error("ncomp<1");

    if (scompBC_map.size()!=ncomp)
     amrex::Error("scompBC_map has invalid size");

    std::vector< std::pair<int,int> > range;

    Interpolater* map = interp(scompBC_map[0]);

    int SComp = 0;
    int NComp = 1;

    for (int i = 1; i < ncomp; i++)
    {
        if (map == interp(scompBC_map[i]))
        {
            NComp++;
        }
        else
        {
            range.push_back(std::pair<int,int>(SComp,NComp));

            map   = interp(scompBC_map[i]);
            SComp = i;
            NComp = 1;
        }
    }

    range.push_back(std::pair<int,int>(SComp,NComp));

#ifndef NDEBUG
    int local_sum = 0;
    for (int i = 0; i < range.size(); i++)
        local_sum += range[i].second;
    BL_ASSERT(local_sum == ncomp);
#endif

    return range;
}

} // namespace amrex
