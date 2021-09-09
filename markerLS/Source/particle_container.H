#ifndef PARTICLE_CONTAINER_H_
#define PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

struct IntData
{
    enum {
        i,
		j,
		k,
		ncomps
    };
};

struct IntData_interface
{
    enum {
        sorted,
        i,
		j,
		k,
		ncomps
    };
};

//ParticleContainer<NStructReal,NStructInt,NArrayReal,NarrayInt> (AoS real&int, SoA real&int)

//inner narrow band is computational region for most operations (need extra cell layer for redistancing)
class PC_MaskInner : public amrex::ParticleContainer<0,0,0, IntData::ncomps>
{
    public:

        using MyParIter = amrex::ParIter<0,0,0, IntData::ncomps>;
		
        PC_MaskInner (const amrex::Geometry            & a_geom,
                      const amrex::DistributionMapping & a_dmap,
                      const amrex::BoxArray            & a_ba);

        
        void writeParticles(int n);
        //void clearParticles();
};

//!!SoA real data is for length of coeff vector : tricubic interp =4x4x4, current least squares a_coeffs:4 for polyOrder 1, 10 for polyOrder 2, 20 for polyOrder 3
//outer narrow band includes inner band (where advection is done) + extra layer outside inner band (need extra layer for redist.)
class PC_MaskOuter : public amrex::ParticleContainer<0,0,0, IntData::ncomps>
{
    public:

        using MyParIter = amrex::ParIter<0,0,0, IntData::ncomps>;
		
        PC_MaskOuter (const amrex::Geometry            & a_geom,
                      const amrex::DistributionMapping & a_dmap,
                      const amrex::BoxArray            & a_ba);

        
        void writeParticles(int n);
};

//interface particles are seeded uniformly in cells where LS changes sign, they are then attracted onto the interface
class PC_Interface : public amrex::ParticleContainer<4,4,0,3> //AoS:phi, G_Ptc(3), i,j,k_cell, sorted   SoA:i,j,k_cell !!check if need cell index in SoA or AoS
{
    public:
        using MyParIter = amrex::ParIter<4,4,0,3>;
		
        PC_Interface (const amrex::Geometry            & a_geom,
                      const amrex::DistributionMapping & a_dmap,
                      const amrex::BoxArray            & a_ba);

        void writeParticles(int n);
        
        void ReBin();
        //void MoveParticles();
        int SumCellVectors();
        int numUnsorted();
        int numWrongCell();
        void visitAllParticles();
        
        void UpdateCellVectors();
        void UpdateFortranStructures();
    //protected:

        //void UpdateCellVectors();
        //void UpdateFortranStructures();

        // used to store vectors of particle indices on a cell-by-cell basis
        std::map<int, amrex::BaseFab<std::vector<int> > > m_cell_vectors;

        // primitive data for passing the above into Fortran
        std::map<int, amrex::BaseFab<int> > m_vector_size;
        std::map<int, amrex::BaseFab<int*> > m_vector_ptrs;
    protected:
        // used to keep track of whether we need to rebuild the above or not
        bool m_vectors_initialized = false;
        amrex::BoxArray::RefID m_BARef;
        amrex::DistributionMapping::RefID m_DMRef;

    private:
        void correctCellVectors(int old_index, int new_index, int grid, const ParticleType& p) override;
  
};

//init inner narrow band mask
void InitParticles_inner(const amrex::Geometry& geom, amrex::MultiFab& a_phi, amrex::iMultiFab& a_mask, 
                        PC_MaskInner& NarrowBand_inner, PC_MaskOuter& NarrowBand_outer);

//init outer narrow band and interface particles
void InitParticles_outer(const amrex::Geometry& geom, amrex::MultiFab& a_phi, amrex::iMultiFab& a_mask, 
                         PC_MaskInner& NarrowBand_inner, PC_MaskOuter& NarrowBand_outer, PC_Interface& InterfaceParticles, bool init_IntfPtc,
                         const int polyOrder, int numrows, int numcols, amrex::Real* w, amrex::Real* P, amrex::Real* A_inv, amrex::Real* a_coeffs);
//outside of class since used to init particles from two diff classes at once



#endif
