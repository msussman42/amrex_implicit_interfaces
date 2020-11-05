
#include "myfunc.H"
#include "myfunc_F.H"

#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

using namespace amrex;

void getGridNormals(MultiFab& phi,
              PC_MaskInner& NarrowBand_inner,
              const Geometry& geom,
              MultiFab& grid_normals){
    const Real* dx = geom.CellSize();
    
    for (MFIter mfi(phi); mfi.isValid(); ++mfi){
        const Box& bx = mfi.validbox();
        int lev = 0;
        auto& particles = NarrowBand_inner.GetParticles(lev)[std::make_pair(mfi.index(),mfi.LocalTileIndex())];
        PC_MaskInner::SoA& NBinnerParticles = particles.GetStructOfArrays();
        //mask index, narrow band
        PC_MaskInner::IntVector& mindex_i = NBinnerParticles.GetIntData(0);
        PC_MaskInner::IntVector& mindex_j = NBinnerParticles.GetIntData(1);
        PC_MaskInner::IntVector& mindex_k = NBinnerParticles.GetIntData(2);
        int* mi_i = mindex_i.dataPtr();
        int* mi_j = mindex_j.dataPtr();
        int* mi_k = mindex_k.dataPtr();
        int num_NB_inner = NBinnerParticles.size();
        
        GridNormals(BL_TO_FORTRAN_ANYD(grid_normals[mfi]),
                    BL_TO_FORTRAN_ANYD(phi[mfi]),
                    mi_i, mi_j, mi_k, 
                    &num_NB_inner, dx);
    }
    grid_normals.FillBoundary(geom.periodicity()); //not needed?
}//end getGridNormals


void advance_LS (MultiFab& phi_old,
              MultiFab& phi,
              iMultiFab& mask,
			  Array<MultiFab, AMREX_SPACEDIM>& flux,
			  Array<MultiFab, AMREX_SPACEDIM>& vel,
              Real dt,
			  Real time,
              const Geometry& geom,
              PC_MaskInner& NarrowBand_inner,
              PC_MaskOuter& NarrowBand_outer,
              MultiFab& G){
    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    phi_old.FillBoundary(geom.periodicity());

    const Real* dx = geom.CellSize();

    const Box& domain_bx = geom.Domain();
    
    //for RK steps --! any way to store as 1d array (length of mask size) on each grid?
    G.setVal(0.0);
    
    for (MFIter mfi(phi); mfi.isValid(); ++mfi){
        const Box& bx = mfi.validbox();
        init_vel(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_BOX(domain_bx),
                     BL_TO_FORTRAN_ANYD(vel[0][mfi]),
                     BL_TO_FORTRAN_ANYD(vel[1][mfi]),
                     BL_TO_FORTRAN_ANYD(vel[AMREX_SPACEDIM-1][mfi]),
                     dx, geom.ProbLo(), geom.ProbHi());
    }
	// Init vel //!!can move into main if vel field const
    
    for (int iter = 1; iter<=3; iter++){ //RK intervals (3 for RK3 Williamson)
        for ( MFIter pti(phi); pti.isValid(); ++pti ){
        //for ( PC_MaskInner::MyParIter pti(NarrowBand_inner,0); pti.isValid(); ++pti ){
            const Box& bx = pti.validbox();
            //std::cout << "box : " << bx << "\n" << "domainbox : " << domain_bx << "\n";
            
            /*
            // Compute fluxes one grid at a time
            compute_flux(BL_TO_FORTRAN_BOX(bx),
                         BL_TO_FORTRAN_BOX(domain_bx),
                         BL_TO_FORTRAN_ANYD(phi_old[pti]),
                         BL_TO_FORTRAN_ANYD(flux[0][pti]),
                         BL_TO_FORTRAN_ANYD(flux[1][pti]),
                        #if (AMREX_SPACEDIM == 3) 
                         BL_TO_FORTRAN_ANYD(flux[AMREX_SPACEDIM-1][pti]),
                        #endif 
                         BL_TO_FORTRAN_ANYD(vel[0][pti]),
                         BL_TO_FORTRAN_ANYD(vel[1][pti]),
                        #if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_ANYD(vel[AMREX_SPACEDIM-1][pti]),
                        #endif
                         dx);
        
            // Advance the solution one grid at a time
            
            //Finite diff
             update_phi(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(phi_old[pti]),
                       BL_TO_FORTRAN_ANYD(phi[pti]),
                       BL_TO_FORTRAN_ANYD(flux[0][pti]),
                       BL_TO_FORTRAN_ANYD(flux[1][pti]),
                       BL_TO_FORTRAN_ANYD(flux[2][pti]),
                       dx, &dt);
            */

            int lev = 0;
            auto& particles = NarrowBand_inner.GetParticles(lev)[std::make_pair(pti.index(),pti.LocalTileIndex())];
            PC_MaskInner::SoA& NBinnerParticles = particles.GetStructOfArrays();
            //mask index, narrow band
            PC_MaskInner::IntVector& mindex_i = NBinnerParticles.GetIntData(0);
            PC_MaskInner::IntVector& mindex_j = NBinnerParticles.GetIntData(1);
            PC_MaskInner::IntVector& mindex_k = NBinnerParticles.GetIntData(2);
            int* mi_i = mindex_i.dataPtr();
            int* mi_j = mindex_j.dataPtr();
            int* mi_k = mindex_k.dataPtr();
            //Vector<int>& NB_index;
            //for (int i=0; i<AMREX_SPACEDIM; i++){
            //    NB_index[i] = NBinnerParticles.GetIntData(i);
            //}
            //auto& NBinnerParticles = NarrowBand_inner.GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
            int num_NB_inner = NBinnerParticles.size();
            //std::cout << "num_NB_inner : " << num_NB_inner << "\n";
            
            //update phi, rk3 time step
            RK3_Williamson(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(phi_old[pti]),
                       BL_TO_FORTRAN_ANYD(phi[pti]),
                       BL_TO_FORTRAN_ANYD(mask[pti]),
                       BL_TO_FORTRAN_ANYD(vel[0][pti]),
                       BL_TO_FORTRAN_ANYD(vel[1][pti]),
                       BL_TO_FORTRAN_ANYD(vel[AMREX_SPACEDIM-1][pti]),
                       dx, &dt, &time, mi_i, mi_j, mi_k, &num_NB_inner, &iter, BL_TO_FORTRAN_ANYD(G[pti]));

        }//end MyParIter pti
        phi.FillBoundary(geom.periodicity());
    }//endfor iter
        
}//end advance_LS

//!!need to put in remove particles function wherever particles move out of cell
void attract_PtcIntf (PC_Interface& InterfaceParticles,
                  PC_MaskInner& NarrowBand_inner,
                  MultiFab& phi,
                  MultiFab& grid_normals,
                  Array<MultiFab, AMREX_SPACEDIM>& vel,
                  const Geometry& geom,
                  const int polyOrder, int numrows, int numcols, Real* w, Real* P, Real* A_inv){
    //InterfaceParticles.UpdateCellVectors();//moved to main
    
    const Real* dx = geom.CellSize();
    const Real* p_lo = geom.ProbLo();
    for ( MFIter pti(phi); pti.isValid(); ++pti ){
        const Box& bx = pti.validbox();
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
   
        int lev = 0;
        auto& particles = InterfaceParticles.GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        PC_Interface::AoS& particles_data = particles.GetArrayOfStructs();
        
        int nstride = particles_data.dataShape().first;
        int num_IntfPtc = particles_data.size();
        
        auto& maskparticles = NarrowBand_inner.GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        PC_MaskInner::SoA& NBinnerParticles = maskparticles.GetStructOfArrays();
        //mask index, narrow band
        PC_MaskInner::IntVector& mindex_i = NBinnerParticles.GetIntData(0);
        PC_MaskInner::IntVector& mindex_j = NBinnerParticles.GetIntData(1);
        PC_MaskInner::IntVector& mindex_k = NBinnerParticles.GetIntData(2);
        int* mi_i = mindex_i.dataPtr();
        int* mi_j = mindex_j.dataPtr();
        int* mi_k = mindex_k.dataPtr();
        int maskSize = NBinnerParticles.size();
        
        //send interface particles onto zero LS
        attractToIntf(particles_data.data(), &num_IntfPtc,
                   InterfaceParticles.m_vector_ptrs[grid_id].dataPtr(),
                   InterfaceParticles.m_vector_size[grid_id].dataPtr(),
                   InterfaceParticles.m_vector_ptrs[grid_id].loVect(),
                   InterfaceParticles.m_vector_ptrs[grid_id].hiVect(),
                   BL_TO_FORTRAN_ANYD(phi[pti]),
                   BL_TO_FORTRAN_ANYD(grid_normals[pti]),
                   BL_TO_FORTRAN_ANYD(vel[0][pti]),
                   BL_TO_FORTRAN_ANYD(vel[1][pti]),
                   BL_TO_FORTRAN_ANYD(vel[AMREX_SPACEDIM-1][pti]),
                   &maskSize, mi_i, mi_j, mi_k,
                   dx, p_lo,
                   &polyOrder, &numrows, &numcols, w, P, A_inv); 
                   
    }//end MyParIter pti
    InterfaceParticles.Redistribute(); //!?check that this should be here and not inside pti loop
    InterfaceParticles.UpdateCellVectors();//InterfaceParticles.ReBin();
}// end attract_PtcIntf



//!!need to put in remove particles function wherever particles move out of cell
void advance_PTC (PC_Interface& InterfaceParticles,
                  MultiFab& phi,
                  Array<MultiFab, AMREX_SPACEDIM>& vel,
                  const Geometry& geom,
                  Real dt, Real time){
                      
    const Real* dx = geom.CellSize();
    const Real* p_lo = geom.ProbLo();
    
    for (MFIter mfi(phi); mfi.isValid(); ++mfi){//!!!this is already in advance_LS, move
        const Box& bx = mfi.validbox();
        init_vel(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_BOX(geom.Domain()),
                     BL_TO_FORTRAN_ANYD(vel[0][mfi]),
                     BL_TO_FORTRAN_ANYD(vel[1][mfi]),
                     BL_TO_FORTRAN_ANYD(vel[AMREX_SPACEDIM-1][mfi]),
                     dx, geom.ProbLo(), geom.ProbHi());
    }
    
    //G_Ptc.setVal(0.0);
    for (int iter = 1; iter<=3; iter++){ //RK intervals
        for ( MFIter pti(phi); pti.isValid(); ++pti ){
        //for ( PC_Interface::MyParIter pti(InterfaceParticles,0); pti.isValid(); ++pti ){
            const Box& bx = pti.validbox();
       
            int lev = 0;
            auto& particles = InterfaceParticles.GetParticles(lev)[std::make_pair(pti.index(),pti.LocalTileIndex())];
            PC_Interface::AoS& particles_data = particles.GetArrayOfStructs();
            
            //int nstride = particles_data.dataShape().first;
            //std::cout << nstride << '\n';
            int num_IntfPtc = particles_data.size();//check if correct num of particles (yes), or if need to divide by nstride
            //std::cout << num_IntfPtc << '\n';
            
            PC_Interface::SoA& intfPtc_SOA = particles.GetStructOfArrays();
            PC_Interface::IntVector& pindex_i = intfPtc_SOA.GetIntData(0);
            PC_Interface::IntVector& pindex_j = intfPtc_SOA.GetIntData(1);
            PC_Interface::IntVector& pindex_k = intfPtc_SOA.GetIntData(2);
            int* pi_i = pindex_i.dataPtr();
            int* pi_j = pindex_j.dataPtr();
            int* pi_k = pindex_k.dataPtr();
            //std::cout << pindex_i.size() << '\n'; //(same as num_IntfPtc)
            
            if (num_IntfPtc > 0){
            //update particle position, rk3 time step
            RK3_Williamson_Ptc(particles_data.data(),
                       BL_TO_FORTRAN_ANYD(vel[0][pti]),
                       BL_TO_FORTRAN_ANYD(vel[1][pti]),
                       BL_TO_FORTRAN_ANYD(vel[AMREX_SPACEDIM-1][pti]),
                       &dt, &time, &num_IntfPtc,
                       pi_i, pi_j, pi_k,
                       dx, p_lo,
                       &iter);
            }          
            
            
            
        }//end MyParIter pti
        //std::cout << " end iter adv ptc \n";
        
        //update particle cell
        InterfaceParticles.Redistribute();
        
        InterfaceParticles.UpdateCellVectors();
        //if a particle is removed from a cell, ReBin() takes those flagged as removed and adds them to correct cell sorted map 
      //  InterfaceParticles.ReBin(); //not currently doing anything | need to add in function to flag particles as removed from cell
    }//endfor iter 
}//end advance_PTC

void redist (MultiFab& phi_old,
             MultiFab& dist,
             MultiFab& phi,
             iMultiFab& mask,
             MultiFab& grid_normals,
			 Array<MultiFab, AMREX_SPACEDIM>& flux,
			 Array<MultiFab, AMREX_SPACEDIM>& vel,
             Real dt,
			 Real time,
             const Geometry& geom,
             PC_MaskInner& NarrowBand_inner,
             PC_MaskOuter& NarrowBand_outer,
             MultiFab& G,
             PC_Interface& InterfaceParticles,
             const int polyOrder, int numrows, int numcols, Real* weights, Real* P, Real* A_inv){

    int redist_method = 1; //1: pde based (Sussman, Smereka, Osher), 2: crossing time (Sethian)
    const Real* dx = geom.CellSize();
    
    correctGridLS(InterfaceParticles, NarrowBand_outer, phi, geom, polyOrder, numrows, numcols, weights, P, A_inv);
    
  //  getGridNormals(phi, NarrowBand_inner, geom, grid_normals); //some normals return 0 ?? otherwise, correct values
    
    MultiFab::Copy(phi_old, phi, 0, 0, 1, 1); //needed for crossing time
    
    amrex::Real tau = 0.0; //pseudotime
    amrex::Real dtau = dx[0]/3; //redistancing time step restriction to be monotone: (dtau/dx)|s_ij|<=1/2
    const int max_iter = 4;//42; //!!placeholder
    int iter_redist = 1; //should only take one or two iterations within tube
    while (iter_redist <= max_iter){
        //amrex::Print() << "  redist iter: " << iter_redist << " \n";   
        for (MFIter mfi(phi); mfi.isValid(); ++mfi){
            const Box& bx = mfi.validbox();
            int lev = 0;
            auto& particles_outer = NarrowBand_outer.GetParticles(lev)[std::make_pair(mfi.index(),mfi.LocalTileIndex())];
            PC_MaskOuter::SoA& NBouterParticles = particles_outer.GetStructOfArrays();
            //mask index, narrow band
            PC_MaskOuter::IntVector& mindex_i = NBouterParticles.GetIntData(0);
            PC_MaskOuter::IntVector& mindex_j = NBouterParticles.GetIntData(1);
            PC_MaskOuter::IntVector& mindex_k = NBouterParticles.GetIntData(2);
            int* mi_i = mindex_i.dataPtr();
            int* mi_j = mindex_j.dataPtr();
            int* mi_k = mindex_k.dataPtr();
            int num_NB_outer = NBouterParticles.size();
            
            //!! add in if num_NB_out > 0
            //amrex::Print() << "  numNBouter: " << num_NB_outer << " \n"; 
            
            if (redist_method == 1){
            Redistance(BL_TO_FORTRAN_BOX(bx),
                    BL_TO_FORTRAN_ANYD(phi[mfi]),
                    BL_TO_FORTRAN_ANYD(mask[mfi]),
                    dx, &dtau, &tau, 
                    mi_i, mi_j, mi_k, 
                    &num_NB_outer); 
            }else if (redist_method == 2){ 
            //auto& innerparticles = NarrowBand_inner.GetParticles(0)[std::make_pair(mfi.index(),mfi.LocalTileIndex())];
            //PC_MaskInner::SoA& NBinnerParticles = innerparticles.GetStructOfArrays();
            //int num_NB_inner = NBinnerParticles.size();
            
            /*GridNormals(BL_TO_FORTRAN_ANYD(grid_normals[mfi]),
                    BL_TO_FORTRAN_ANYD(phi[mfi]),
                    mi_i, mi_j, mi_k, 
                    &num_NB_outer, dx); */
            crossingTime(BL_TO_FORTRAN_BOX(bx), dx,
                        BL_TO_FORTRAN_ANYD(phi[mfi]),
                        BL_TO_FORTRAN_ANYD(phi_old[mfi]), 
                        BL_TO_FORTRAN_ANYD(dist[mfi]), 
                        BL_TO_FORTRAN_ANYD(mask[mfi]), 
                        BL_TO_FORTRAN_ANYD(grid_normals[mfi]),
                        mi_i, mi_j, mi_k, &num_NB_outer,
                        &tau, &dtau,
                        &iter_redist, &max_iter);
            }
        }
        if (iter_redist == 21){
            MultiFab::Copy(phi_old, phi, 0, 0, 1, 1); 
            tau = 0.0;
        }
        phi.FillBoundary(geom.periodicity());
        phi_old.FillBoundary(geom.periodicity());
        
        
        correctGridLS(InterfaceParticles, NarrowBand_outer, phi, geom, polyOrder, numrows, numcols, weights, P, A_inv);
        iter_redist = iter_redist + 1;
        
        //if (time>=.292 && time<=.293){
        //    const std::string& pltfile = amrex::Concatenate("plt",iter_redist,5);
        //    WriteSingleLevelPlotfile(pltfile, phi, {"phi"}, geom, tau, iter_redist);
        //}
        ///break; //!!break for crossingTime, temporary
    }//end while
}//end redist




        
void correctGridLS (PC_Interface& InterfaceParticles,
                    PC_MaskOuter& NarrowBand_outer,
                    MultiFab& phi,
                    const Geometry& geom,
                    const int polyOrder, int numrows, int numcols, Real* w, Real* P, Real* A_inv){
                        
    std::cout << " correcting GridLS using Interface Particles \n";
    
    //amrex::Print() << " Total particles " << InterfaceParticles.TotalNumberOfParticles() << std::endl;
    int intf_ptc_count = 0;
    
    //InterfaceParticles.UpdateCellVectors();//moved to main
    const Real* dx = geom.CellSize();
    const Real* p_lo = geom.ProbLo();
    //for ( MFIter pti(phi); pti.isValid(); ++pti ){
    for ( PC_Interface::MyParIter pti(InterfaceParticles,0); pti.isValid(); ++pti ){
        const Box& bx = pti.validbox(); //tilebox()
   
        int lev = 0;
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        
        auto& particle_tile = InterfaceParticles.GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        PC_Interface::AoS& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();
        //std::cout << " np:  " << np << std::endl;
        
        intf_ptc_count = intf_ptc_count + np;
        /*
        for(int pindex = 0; pindex < np; ++pindex){
            auto& p = particles[pindex];
            
            //int nstride = particles.dataShape().first;
            //int num_IntfPtc = particles.size();
            const IntVect& iv = InterfaceParticles.Index(p, lev);
            auto test = InterfaceParticles.m_cell_vectors[pti.index()](iv);
            amrex::Print() << " Total particles " << InterfaceParticles.TotalNumberOfParticles() << std::endl;
            amrex::Print() << " Number of particles in cell vectors is " << test[1] << std::endl;
        }
        */
        
        PC_Interface::SoA& intfPtc_SOA = particle_tile.GetStructOfArrays();
        PC_Interface::IntVector& pindex_i = intfPtc_SOA.GetIntData(0);
        PC_Interface::IntVector& pindex_j = intfPtc_SOA.GetIntData(1);
        PC_Interface::IntVector& pindex_k = intfPtc_SOA.GetIntData(2);
        int* pi_i = pindex_i.dataPtr(); //!!not used anymore, keeping temporarily but will prob remove (already stored in AoS)
        int* pi_j = pindex_j.dataPtr();
        int* pi_k = pindex_k.dataPtr();
        
        
        //!using outerband since easier for a_coeffs (can change later)
        auto& outerNBmask = NarrowBand_outer.GetParticles(lev)[std::make_pair(pti.index(),pti.LocalTileIndex())];
        PC_MaskOuter::SoA& NBouterParticles = outerNBmask.GetStructOfArrays();
        //mask index, narrow band
        PC_MaskOuter::IntVector& mindex_i = NBouterParticles.GetIntData(0);
        PC_MaskOuter::IntVector& mindex_j = NBouterParticles.GetIntData(1);
        PC_MaskOuter::IntVector& mindex_k = NBouterParticles.GetIntData(2);
        int* mi_i = mindex_i.dataPtr();
        int* mi_j = mindex_j.dataPtr();
        int* mi_k = mindex_k.dataPtr();
        int maskSize = NBouterParticles.size();
        /* //finding a_coeffs inside func, don't need to store
        //!!!find how to put this into an array of vectors?
        PC_MaskOuter::RealVector& coeffs_mask0 = NBouterParticles.GetRealData(0);
        PC_MaskOuter::RealVector& coeffs_mask1 = NBouterParticles.GetRealData(1);
        PC_MaskOuter::RealVector& coeffs_mask2 = NBouterParticles.GetRealData(2);
        PC_MaskOuter::RealVector& coeffs_mask3 = NBouterParticles.GetRealData(3);
        Real* a_coeffs_mask0 = coeffs_mask0.dataPtr();
        Real* a_coeffs_mask1 = coeffs_mask1.dataPtr();
        Real* a_coeffs_mask2 = coeffs_mask2.dataPtr();
        Real* a_coeffs_mask3 = coeffs_mask3.dataPtr();
        */
        
        minimizeError(particles.data(), &np,
                      bx.loVect(), bx.hiVect(),
                      BL_TO_FORTRAN_ANYD(phi[pti]),
                      InterfaceParticles.m_vector_ptrs[grid_id].dataPtr(),
                      InterfaceParticles.m_vector_size[grid_id].dataPtr(),
                      InterfaceParticles.m_vector_ptrs[grid_id].loVect(),
                      InterfaceParticles.m_vector_ptrs[grid_id].hiVect(),
                      mi_i, mi_j, mi_k, &maskSize,
                      p_lo, dx,
                      &polyOrder, &numrows, &numcols, w, P, A_inv); 
                  
    }//end MyParIter pti
    
    //amrex::Print() << " particle count " << intf_ptc_count << std::endl;
    
    phi.FillBoundary(geom.periodicity());
}// end correctGridLS

