#ifndef _IFRAME_H_
#define _IFRAME_H_

#if (BL_SPACEDIM == 2 || BL_SPACEDIM == 3)
#include <CArena.H>
#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>

// describes linear transformation between integer coords and real coords
// in 2d.  useful for IGraphTool
class IFrame {
	private:
		int ilo,ihi,jlo,jhi;
		Real xlo,xhi,ylo,yhi;
		Real delx,dely;
	public:
		IFrame(){}
		IFrame(const Box&,Real xmin,Real xmax, 
				Real ymin, Real ymax);
		~IFrame(){}
		
		IntVect toIV( const Real* coord) const;
		Real toX(const IntVect) const;
		Real toY(const IntVect) const;
};
#endif
#endif
