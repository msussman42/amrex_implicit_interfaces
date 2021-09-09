
#include "myfunc.H"
#include "myfunc_F.H"

#include "particle_container.H"

using namespace amrex;

/*
namespace {
    
    void get_position_unit_cell(Real* r, const IntVect& N, int i_part){
        //uniformly place particles in cell, N in each dir
        //N^3 particles per cell
        int Nppc = pow(N,3);
        
    }
    
}
*/

PC_MaskInner::PC_MaskInner(const Geometry  & a_geom,
             const DistributionMapping     & a_dmap,
             const BoxArray                & a_ba)
    : ParticleContainer<0,0,0, IntData::ncomps> (a_geom, a_dmap, a_ba){}
    
PC_MaskOuter::PC_MaskOuter(const Geometry  & a_geom,
             const DistributionMapping     & a_dmap,
             const BoxArray                & a_ba)
    : ParticleContainer<0,0,0, IntData::ncomps> (a_geom, a_dmap, a_ba){}
	
PC_Interface::PC_Interface(const Geometry  & a_geom,
             const DistributionMapping     & a_dmap,
             const BoxArray                & a_ba)
    : ParticleContainer<4,4,0,3> (a_geom, a_dmap, a_ba){} //!!can this move to particle_operations.cpp instead?
    

void InitParticles_inner(const Geometry& geom, MultiFab& phi_mf, iMultiFab& mask, 
                            PC_MaskInner& NarrowBand_inner, PC_MaskOuter& NarrowBand_outer){
    //Sweep through domain, mask cells within dist. gamma of interface as inner narrow band (computation region)
    //init particle mask locations 
    
	const int lev = 0;
	//const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* p_lo = geom.ProbLo();
	//const Real* p_hi = geom.ProbHi();
	
	double gamma = 6*dx[0]; //6dx for 5th order WENO, 4dx for 3rd order WENO
	
    //!! turning on tiling is causing issues
    for (MFIter mfi(phi_mf,false); mfi.isValid(); ++mfi){ 
		//'particles' starts off empty
		
		//tiling
		const Box& tile_box  = mfi.tilebox();
		Dim3 lo = lbound(tile_box);
		Dim3 hi = ubound(tile_box);
		//auto const& phi = phi_mf.array(mfi);
		const FArrayBox& phi = phi_mf[mfi]; //!does this include ghost cells? -yes -need FillBoundary()
        //const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
		
		const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
		auto& particles = NarrowBand_inner.GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles_outerband = NarrowBand_outer.GetParticles(0)[std::make_pair(grid_id,tile_id)];
		
		//iterate over tile (C++11 range-based for loop)
		for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)){
			int i = iv[0];
			int j = iv[1];
			int k = iv[2];
            
			if (fabs(phi(iv)) < gamma){
				PC_MaskInner::ParticleType p_NBinner;
				p_NBinner.id()  = PC_MaskInner::ParticleType::NextID();
				p_NBinner.cpu() = ParallelDescriptor::MyProc();
				
				p_NBinner.pos(0) = p_lo[0] + (iv[0]+0.5)*dx[0];
				p_NBinner.pos(1) = p_lo[1] + (iv[1]+0.5)*dx[1];
				p_NBinner.pos(2) = p_lo[2] + (iv[2]+0.5)*dx[2];
			
				// AoS int data
				//NarrowBand_inner.idata(0) = 2; //inner narrow band
				//SoA int data
				std::array<int, 3> int_attribs;
				int_attribs[IntData::i] = i;
				int_attribs[IntData::j] = j;
				int_attribs[IntData::k] = k;
				
				//AMREX_ASSERT(this->Index(NarrowBand_inner, lev) == iv);
				//append new particles
				particles.push_back(p_NBinner);
				particles.push_back_int(int_attribs);
                
                //full narrow band includes inner + 1 layer of cells outside
                particles_outerband.push_back(p_NBinner);
				particles_outerband.push_back_int(int_attribs);
                
                //mask -> 3: contains interface, 2: inner narrow band, 1: outer narrow band, 0: not in computational region
                mask[mfi](iv) = 2;
			}else{
                mask[mfi](iv) = 0;
            }
		}//endfor IntVect iv
    }//endfor MFIter mfi
	mask.FillBoundary(geom.periodicity());
    //NarrowBand_inner.Redistribute(); //not needed since '.GetParticles(lev)[std::make_pair(grid_id,tile_id)]' statement  ??
}//END InitParticles_inner()



void PC_MaskInner::writeParticles(int n)
{
    BL_PROFILE("PC_MaskInner::writeParticles");
    const std::string& pltfile = amrex::Concatenate("InnerNarrowBand", n, 5);
    WriteAsciiFile(pltfile);
}
/*
void ParticleContainer_Mask::clearParticles(){
	const int lev = 0;
	using MyParIter = ParIter<0,1,0,3>;
	for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
		*this.id()=-1;
	}
}
*/









void InitParticles_outer(const Geometry& geom, MultiFab& phi_mf, iMultiFab& mask_imf, 
                         PC_MaskInner& NarrowBand_inner, PC_MaskOuter& NarrowBand_outer, PC_Interface& InterfaceParticles, bool init_IntfPtc,
                         const int polyOrder, int numrows, int numcols, Real* w, Real* P, Real* A_inv, Real* a_coeffs){
    //sweep through innerband and flag cells just outside as outer narrow band (use one extra cell outside for redist.)
    
	const int lev = 0;
	//const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* p_lo = geom.ProbLo();
	//const Real* p_hi = geom.ProbHi();
    
    int i,j,k; //indices
    int p,q,r; //iterators around indices
    
    // Add to level 0, grid 0, tile 0   (lev)[grid,tile]
    // Redistribute() will move them to the proper place.
    auto& particles_outerband = NarrowBand_outer.GetParticles(0)[std::make_pair(0,0)];
    auto& particles_interface = InterfaceParticles.GetParticles(0)[std::make_pair(0,0)];
    
    bool interfaceCell_init = 0;
	
    //iterate through inner narrowband particles, flag cells just outside this as outer narrowband
	//using MyParIter = ParIter<0,0,0, IntData::ncomps>;
    for (PC_MaskInner::MyParIter pti(NarrowBand_inner, lev); pti.isValid(); ++pti){
        const Box& bx_mask = pti.validbox();
        //const Box& box = fabArray.box((*mask_imf)[pti]);
        Dim3 lo = lbound(bx_mask);
		Dim3 hi = ubound(bx_mask);
        //mask_imf.FillBoundary(geom.periodicity());
        IArrayBox& mask = mask_imf[pti];
        phi_mf.FillBoundary(geom.periodicity());//!check if still needed
        const FArrayBox& phi = phi_mf[pti];
        
        auto& particle_attributes = pti.GetStructOfArrays();
        PC_MaskInner::IntVector& i_index = particle_attributes.GetIntData(0);
        PC_MaskInner::IntVector& j_index = particle_attributes.GetIntData(1);
        PC_MaskInner::IntVector& k_index = particle_attributes.GetIntData(2);
        //auto& i_index = particle_attributes.GetIntData(0);
        //auto& j_index = particle_attributes.GetIntData(1);
        //auto& k_index = particle_attributes.GetIntData(2);
        
        //iterate through indices (SoA data)
        for (int nb_inner = 0; nb_inner < pti.numParticles(); ++nb_inner) {
            i = i_index[nb_inner];
            j = j_index[nb_inner];
            k = k_index[nb_inner];
            interfaceCell_init = 0;
            
            if (mask({i,j,k}) != 2){
                amrex::Error("Error -- Mask multifab & particle disagreement");
            }
            //iterate in region immediately around inner band cells
            for (r=k-1; r<=k+1; r++){
            for (q=j-1; q<=j+1; q++){
            for (p=i-1; p<=i+1; p++){
                if (p!=i || q!=j || r!=k){
                    //if just outside of inner narrowband, flag as outer narrowband
                    if (mask({p,q,r}) == 0){
                        mask({p,q,r}) = 1; 
                        //3: contains interface, 2:inner narrow band, 1: outer narrow band, 0: not in computational region
                        
                        PC_MaskOuter::ParticleType p_NBouter;
                        p_NBouter.id()  = PC_MaskOuter::ParticleType::NextID();
                        p_NBouter.cpu() = ParallelDescriptor::MyProc();
                        
                        p_NBouter.pos(0) = p_lo[0] + (p+0.5)*dx[0];
                        p_NBouter.pos(1) = p_lo[1] + (q+0.5)*dx[1];
                        p_NBouter.pos(2) = p_lo[2] + (r+0.5)*dx[2];
                    
                        // AoS int data
                        //p_NBouter.idata(0) = 1; //outer narrow band
                        //SoA int data
                        std::array<int, 3> int_attribs;
                        int_attribs[IntData::i] = p;
                        int_attribs[IntData::j] = q;
                        int_attribs[IntData::k] = r;
                        
                        //AMREX_ASSERT(this->Index(NarrowBand_outer, lev) == IntVect({i,j,k}));
                        //append new particles
                        particles_outerband.push_back(p_NBouter);
                        particles_outerband.push_back_int(int_attribs);
                    }//endif mask  
                    //if phi changes sign, flag as interface cell
                    if ((phi({i,j,k})*phi({p,q,r}) <= 0) &&
                        (p==i || q==j || r==k) &&//if not counting corners
                        interfaceCell_init != 1 &&
                        init_IntfPtc == true //check if interface particles need to be re-initialized
                        ){
                        mask({i,j,k}) = 3;
                        interfaceCell_init = 1; //mark cell as initialized as interface cell so as to not double count
                        
                        int N = 2; // N^3 particles per interface cell !hardcoded, change to input file
                        if (N<=0 || N>5){
                            amrex::Error("Error -- Fix num interface particles per cell");
                        }
                        //uniformly place particles in cell, N in each dir
                        //N^3 particles per cell
                        //int Nppc = pow(N,3);
                        
                        int currentIndex[3] = {i,j,k};
                        //get polynomial approximation to phi within cell // interpolation of phi
                        C_getPolyInterpCoeffs(&polyOrder, &numrows, &numcols, dx, w, P, A_inv, BL_TO_FORTRAN_ANYD(phi_mf[pti]), currentIndex, a_coeffs);
                        
                        std::array<Real, 4> a_coeffs_temp;
                        for (int a=0; a<4; a++){
                            a_coeffs_temp[a] = a_coeffs[a];
                        }
                        
                        for (int kk = 1; kk<=N; kk++){
                        for (int jj = 1; jj<=N; jj++){
                        for (int ii = 1; ii<=N; ii++){
                            
                            PC_Interface::ParticleType p_interface;
                            p_interface.id()  = PC_Interface::ParticleType::NextID();
                            p_interface.cpu() = ParallelDescriptor::MyProc();
                            
                            double cellCenter_x = p_lo[0] + (i+0.5)*dx[0];
                            double cellCenter_y = p_lo[1] + (j+0.5)*dx[1];
                            double cellCenter_z = p_lo[2] + (k+0.5)*dx[2];
                            double ptcPos_x = p_lo[0] + (i+ii/(N+1.0))*dx[0];
                            double ptcPos_y = p_lo[1] + (j+jj/(N+1.0))*dx[1];
                            double ptcPos_z = p_lo[2] + (k+kk/(N+1.0))*dx[2];
                            
                            p_interface.pos(0) = ptcPos_x;
                            p_interface.pos(1) = ptcPos_y;
                            p_interface.pos(2) = ptcPos_z;
                            
                            
                            //// AoS real data
                            //phi
                            double polyval=0;
                            int index_a=0;
                            for (int r_p = 0; r_p<=polyOrder; r_p++){
                            for (int t_p = 0; t_p<=polyOrder; t_p++){
                            for (int s_p = 0; s_p<=polyOrder; s_p++){
                                if (s_p + t_p + r_p > polyOrder){
                                    //   %do nothing, cycle
                                }else{
                                polyval = polyval + a_coeffs[index_a]*pow(ptcPos_x-cellCenter_x, s_p)
                                                        *pow(ptcPos_y-cellCenter_y, t_p)
                                                        *pow(ptcPos_z-cellCenter_z, r_p);
                                index_a = index_a+1;
                                } //endif
                            }
                            }
                            }
                            
                            
                            p_interface.rdata(0) = polyval;//phi({i,j,k}); //!!placeholder
                            //G_Ptc //!!prob move to SoA
                            p_interface.rdata(1) = 0.0; //for RK-iterations
                            p_interface.rdata(2) = 0.0;
                            p_interface.rdata(3) = 0.0;
                            
                            p_interface.idata(0) = 0; //sorted
                            p_interface.idata(1) = i; //cell indices
                            p_interface.idata(2) = j; 
                            p_interface.idata(3) = k; 
                            
                            
                            
                            //SoA int data
                            // cell index//!!currently also stored in AoS (prob need AoS & not SoA for this?)
                            std::array<int, 3> int_attribs;
                            int_attribs[IntData::i] = i;
                            int_attribs[IntData::j] = j;
                            int_attribs[IntData::k] = k;
                        
                            
                            //AMREX_ASSERT(this->Index(InterfaceParticles, lev) == IntVect({i,j,k}));
                            //append new particles
                            particles_interface.push_back(p_interface);
                            particles_interface.push_back_int(int_attribs);
                            
                        }//endfor ii
                        }//endfor jj
                        }//endfor kk
                    }//endif phi changes sign 
                }//endif p,q,r != i,j,k
            }//endfor p
            }//endfor q
            }//endfor r
        }//endfor nb_inner
        
    }//endfor MyParIter pti
	//put new particles in correct container
	NarrowBand_outer.Redistribute();
    InterfaceParticles.Redistribute();
    mask_imf.FillBoundary(geom.periodicity());
}//END InitParticles_outer()

void PC_MaskOuter::writeParticles(int n)
{
    BL_PROFILE("PC_MaskOuter::writeParticles");
    const std::string& pltfile = amrex::Concatenate("OuterNarrowBand", n, 5);
    WriteAsciiFile(pltfile);
}


void PC_Interface::writeParticles(int n)
{
    BL_PROFILE("PC_MaskOuter::writeParticles");
    const std::string& pltfile = amrex::Concatenate("InterfaceParticles", n, 5);
    WriteAsciiFile(pltfile);
}
