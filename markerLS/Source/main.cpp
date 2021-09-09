
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "myfunc.H"
#include "myfunc_F.H"
#include "INIT_MARCHING_F.H"

#include "particle_container.H"

using namespace amrex;




int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    main_main();
    
    amrex::Finalize();
    return 0;
}

void main_main ()
{
    // start clock for total runtime
    Real strt_time = amrex::second();

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, nsteps, plot_int;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default
    
    bool plotEulerian, plotParticles;

    // inputs parameters
    {
        // ParmParse : read from the inputs file
        ParmParse pp;

        // n_cell is number of cells on each side of square/cubic domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Default nsteps to 10, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        pp.queryarr("is_periodic", is_periodic);
        
        
        pp.query("plotEulerian", plotEulerian);
        pp.query("plotParticles", plotParticles);
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

       // This defines the physical box, [0,1] in each direction.
        Real local_xlo[3];
	local_xlo[0]=0.0;
	local_xlo[1]=0.0;
	local_xlo[2]=0.0;
        RealBox real_box({AMREX_D_DECL(local_xlo[0],local_xlo[1],local_xlo[2])},
                         {AMREX_D_DECL( 1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array 
    int Nghost = 3;
    // Ncomp = number of components for each array
    int Ncomp  = 1;
  
    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two phi multifabs; one will store the old state, the other the new.
    MultiFab phi(ba, dm, Ncomp, Nghost);
    MultiFab grid_normals(ba, dm, 3, 1); //normal directions for each Eulerian grid cell
    iMultiFab phi_NBmask(ba, dm, Ncomp, 1); //! make local?
    
    MultiFab phi_old(ba, dm, 1, 1);
    MultiFab dist(ba, dm, 1, 1);
    
    MultiFab G(ba, dm, 1, 0);
    //MultiFab G_Ptc(ba, dm, 3, 0);

    // Initialize phi by calling a Fortran routine.
    // MFIter = MultiFab Iterator
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi ){
        const Box& bx = mfi.validbox();

        init_phi(BL_TO_FORTRAN_BOX(bx),
                 BL_TO_FORTRAN_ANYD(phi[mfi]),
                 geom.CellSize(), geom.ProbLo(), geom.ProbHi());
        phi.FillBoundary(geom.periodicity());
    }

    const amrex::Real* xlo_fort=geom.ProbLo();
    int n_cell_fort[3];
    n_cell_fort[0]=n_cell;
    n_cell_fort[1]=n_cell;
    n_cell_fort[2]=n_cell;
    fort_init_marching(xlo_fort,geom.CellSize(),n_cell_fort);

    // compute the time step
    const amrex::Real* dx=geom.CellSize();
    Real dt = dx[0]/8;// 0.9*dx[0]*dx[0] / (2.0*AMREX_SPACEDIM);
    amrex::Print() << "dt: " << dt << "\n";
    amrex::Print() << "dx: " << dx[0] << "\n";

    // time = starting time in the simulation
    Real time = 0.0;

    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0 && plotEulerian){
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
        WriteSingleLevelPlotfile(pltfile, phi, {"phi"}, geom, time, 0);
    }

    // build the flux multifabs
    Array<MultiFab, AMREX_SPACEDIM> flux;
	Array<MultiFab, AMREX_SPACEDIM> vel;
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++){
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);
        flux[dir].define(edge_ba, dm, 1, 0);
		vel[dir].define(edge_ba, dm, 1, 0);
    }
	
	
	const int polyOrder = 3;
    int numrows, numcols;
    polyMatrixDim(&numrows, &numcols, &polyOrder); //find numrows/numcols for given polyOrder
    //std::cout << "numrows: " << numrows << ", numcols: " << numcols << std::endl;
    Real* weights = new Real[numrows];
    //P dim: numrows x numcols
    Real* P = new Real[numrows*numcols];
    //A_inv dim: numcols x numcols
    Real* A_inv = new Real[numcols*numcols];
    Real* a_coeffs = new Real[numcols];
    
    init_Pmatrix(&polyOrder, &numrows, &numcols, dx, weights, P);
    init_LeastSquaresMatrix(&numrows, &numcols, weights, P, A_inv);
    
	
    //create inner narrow band
	PC_MaskInner NarrowBand_inner(geom, dm, ba);
    
    //create outer narrow band
    PC_MaskOuter NarrowBand_outer(geom, dm, ba);
    PC_Interface InterfaceParticles(geom, dm, ba);
    amrex::Print() << "init inner \n";
    InitParticles_inner(geom, phi, phi_NBmask, NarrowBand_inner, NarrowBand_outer);
    bool init_IntfPtc = true;
    amrex::Print() << "init outer \n";
    InitParticles_outer(geom, phi, phi_NBmask, NarrowBand_inner, NarrowBand_outer, InterfaceParticles, init_IntfPtc,
                        polyOrder, numrows, numcols, weights, P, A_inv, a_coeffs);
    amrex::Print() << "particles initialized \n";
    InterfaceParticles.UpdateCellVectors();//initialize cell sorting vector
    
    //move particles to interface
    getGridNormals(phi, NarrowBand_inner, geom, grid_normals);
    attract_PtcIntf(InterfaceParticles, NarrowBand_inner, phi, grid_normals, vel, geom, polyOrder, numrows, numcols, weights, P, A_inv);
    //!!!^put this back in!
    
    //!!!TESTING
    //InterfaceParticles.UpdateCellVectors();
    //correctGridLS (InterfaceParticles, NarrowBand_outer, phi, geom, polyOrder, numrows, numcols, weights, P, A_inv);
    
    if (plotParticles){
        NarrowBand_inner.writeParticles(0);
        NarrowBand_outer.writeParticles(0);
        InterfaceParticles.writeParticles(0);
    }
	//const std::string& pltmask = "plt_m";
    //WriteSingleLevelPlotfile(pltmask, GetVecOfConstPtrs(narrow_band_pc), {"mask"}, geom,0,0);
	
    init_IntfPtc = false; //reinit interface particles

    for (int n = 1; n <= nsteps; ++n){
        if (time>2.1){break;}
        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "step " << n << " time: "<< time << "\n";
        
        MultiFab::Copy(phi_old, phi, 0, 0, 1, 1);
		//MultiFab::Copy(mfdst, mfsrc, sc, dc, nc, ng); // Copy from mfsrc to mfdst
		// MultiFab mfdst: destination
		// MultiFab mfsrc: source
		// int      sc   : starting component index in mfsrc for this operation
		// int      dc   : starting component index in mfdst for this operation
		// int      nc   : number of components for this operation
		// int      ng   : number of ghost cells involved in this operation
		//                 mfdst and mfsrc may have more ghost cells
        
        if (n%100 == 0){
            //init_IntfPtc = true;
        }else{
            init_IntfPtc = false;
        }
        

        // new_phi = old_phi + dt * (something)
        amrex::Print() << " advect particles \n";
        advance_PTC(InterfaceParticles, phi, vel, geom, dt, time);//, G_Ptc);
        amrex::Print() << " advect LS \n";
        advance_LS(phi_old, phi, phi_NBmask, flux, vel, dt, time, geom, NarrowBand_inner, NarrowBand_outer, G);    
        amrex::Print() << " redistancing \n";        
        redist(phi_old, dist, phi, phi_NBmask, grid_normals, flux, vel, dt, time, geom, NarrowBand_inner, NarrowBand_outer, G, InterfaceParticles, polyOrder, numrows, numcols, weights, P, A_inv);
        //phi.FillBoundary(geom.periodicity());//!check if needed
            //!! moved two lines below into redist...
        //InterfaceParticles.UpdateCellVectors();
        //correctGridLS(InterfaceParticles, NarrowBand_outer, phi, geom, polyOrder, numrows, numcols, weights, P, A_inv);
        time = time + dt;
        
        phi.FillBoundary(geom.periodicity()); //check if still needed
        //destroy all mask particles
		NarrowBand_inner.~PC_MaskInner();
		new(&NarrowBand_inner) PC_MaskInner(geom, dm, ba);
        NarrowBand_outer.~PC_MaskOuter();
		new(&NarrowBand_outer) PC_MaskOuter(geom, dm, ba);
        if (init_IntfPtc){
            InterfaceParticles.~PC_Interface();
            new(&InterfaceParticles) PC_Interface(geom, dm, ba);
        }
        amrex::Print() << " re-initializing mask \n";
        //re-initialize mask particles //!can avoid this every iteration w/ AdalsteinssonSethian1995 narrow band update? --NO reinit every iter instead
		InitParticles_inner(geom, phi, phi_NBmask, NarrowBand_inner, NarrowBand_outer);
        InitParticles_outer(geom, phi, phi_NBmask, NarrowBand_inner, NarrowBand_outer, InterfaceParticles, init_IntfPtc,
                            polyOrder, numrows, numcols, weights, P, A_inv, a_coeffs);
        
        /*
        //interpolate vel to cell centers? or adjust geom
        if (n==1){ //vel field const. for now
            amrex::Print() << " -> writing velocity field \n";   
            std::string pltfile = "u";
            WriteSingleLevelPlotfile(pltfile, vel[0], {"u"}, geom, time, n);
            pltfile = "v";
            WriteSingleLevelPlotfile(pltfile, vel[1], {"v"}, geom, time, n);
            pltfile = "w";
            WriteSingleLevelPlotfile(pltfile, vel[2], {"w"}, geom, time, n);
        } */
        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if ((plot_int > 0 && plotEulerian) && n%plot_int == 0){
            amrex::Print() << " -> writing Eulerian data \n";   
            const std::string& pltfile = amrex::Concatenate("plt",n,5);
            WriteSingleLevelPlotfile(pltfile, phi, {"phi"}, geom, time, n);
            
            //const std::string& pltfile_mask = amrex::Concatenate("mask",n,5);
            //WriteSingleLevelPlotfile(pltfile_mask, phi_innermask, {"mask"}, geom, time, n);
        }
        
        if (plotParticles && n%plot_int == 0){
            amrex::Print() << " -> writing Lagrangian data \n";  
            NarrowBand_inner.writeParticles(n);
            NarrowBand_outer.writeParticles(n);
            InterfaceParticles.writeParticles(n);
        }
        //amrex::Print() << "number of mask particles " << narrow_band_pc.GetParticles() << "\n";
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = amrex::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
    
    
    
    
    //delete arrays initialized w/ new[]
    //for(int i = 0; i < numrows; i++){
    //    delete[] P[i];
    //}
     
    //for(int i = 0; i < numcols; i++){
    //    delete[] A_inv[i];
    //}
    delete[] P, A_inv, weights, a_coeffs;
}
