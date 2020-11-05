#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

#include "particle_container.H"

void main_main ();

void getGridNormals(amrex::MultiFab& phi,
              PC_MaskInner& NarrowBand_inner,
              const amrex::Geometry& geom,
              amrex::MultiFab& grid_normals);
              
void advance_LS (amrex::MultiFab& phi_old,
              amrex::MultiFab& phi_new,
              amrex::iMultiFab& mask,
			  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& flux,
			  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& vel,
			  amrex::Real dt,
              amrex::Real time,
              const amrex::Geometry& geom,
              PC_MaskInner& NarrowBand_inner,
              PC_MaskOuter& NarrowBand_outer,
              amrex::MultiFab& G);
              
void attract_PtcIntf (PC_Interface& InterfaceParticles,
                  PC_MaskInner& NarrowBand_inner,
                  amrex::MultiFab& phi,
                  amrex::MultiFab& grid_normals,
                  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& vel,
                  const amrex::Geometry& geom,
                  const int polyOrder, int numrows, int numcols, amrex::Real* w, amrex::Real* P, amrex::Real* A_inv);
              
void advance_PTC (PC_Interface& InterfaceParticles,
                  amrex::MultiFab& phi,
                  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& vel,
                  const amrex::Geometry& geom,
                  amrex::Real dt,
                  amrex::Real time);
              
void redist (amrex::MultiFab& phi_old,
              amrex::MultiFab& phi_update,
              amrex::MultiFab& phi,
              amrex::iMultiFab& mask,
              amrex::MultiFab& grid_normals,
			  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& flux,
			  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& vel,
			  amrex::Real dt,
              amrex::Real time,
              const amrex::Geometry& geom,
              PC_MaskInner& NarrowBand_inner,
              PC_MaskOuter& NarrowBand_outer,
              amrex::MultiFab& G,
              PC_Interface& InterfaceParticles,
              const int polyOrder, int numrows, int numcols, amrex::Real* w, amrex::Real* P, amrex::Real* A_inv);
              
void correctGridLS (PC_Interface& InterfaceParticles,
                    PC_MaskOuter& NarrowBand_outer,
                    amrex::MultiFab& phi,
                    const amrex::Geometry& geom,
                    const int polyOrder, int numrows, int numcols, amrex::Real* w, amrex::Real* P, amrex::Real* A_inv);
              
#endif
