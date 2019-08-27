#if (HAS_WINDOWS)

#include <winstd.H>

#include <algorithm>
#include <vector>

#if defined(BL_OLD_STL)
#include <stdio.h>
#include <math.h>
#include <string.h>
#else
#include <cstdio>
#include <cmath>
#include <cstring>
#endif

#include "PSfile.H"

#define MULTIFRAME 1

void PSfile::flush() {file << '\n'; }

// this is a width of a linewidth=1 line in points
static const Real WIDTH = 0.24; 

// this is the X dimension of units of WIDTH points
static const int    FACTOR = 3196;

// these margins place the entire image within the printable
// region for the laser writer.
static const int    XMARGIN = 52;
static const int    YMARGIN = 76;

#define DRAWPATH if (npts>0) {file << "S\n";}; npts = 0;

void PSfile::PSheader(int nstep)
{
char snum[8];

//   sprintf(snum,"%04d",framecount);
   sprintf(snum,"%06d",nstep);
   strcpy(filename,save_name);
   if (MULTIFRAME==1)
    strcat(filename,snum);

   strcat(filename,".ps");
   file.open(filename,std::ios::out);
   if (!file) {
	   std::cerr << "cant open postscript file: " << filename << '\n';
      abort();
   };

   // write out header
   file << "%!PS\n";
   file << "%  " << save_wid << "   " << save_high << '\n';
   file << WIDTH << ' ' << WIDTH << " scale\n";
   file << "-90 rotate\n";
   file << -FACTOR-XMARGIN << ' ' << YMARGIN << " translate\n";
   file << "1 setlinecap\n";
   file << "1 setlinejoin\n";
   file << "/M {moveto} def\n";
   file << "/L {lineto} def\n";
   file << "/S {stroke} def\n";
   file << '\n';
} // PSheader()

PSfile::PSfile(int wid, int high, const char *name,int nstep)
{

   // create postscript file
   framecount=1;
   save_wid=wid;
   save_high=high;

   save_name = new char[strlen(name)+ 5];
   strcpy(save_name,name);
   filename = new char[strlen(name) + 10];
   strcpy(filename,save_name);
   

   PSheader(nstep);
   // init scaling factors
   Real dmax = Real( std::max(wid,high) );
   Real xlen = Real(wid)/dmax;
   Real ylen = Real(high)/dmax;
   xfactor = int( xlen * FACTOR );
   yfactor = int( ylen * FACTOR );
   xcur = ycur = 0;
   npts = 0;
}

PSfile::~PSfile()
{
   DRAWPATH;
   file << "showpage\n";
   file << '\n';
   file.close();
   delete filename;
   delete save_name;
}

void PSfile::newPage(int nstep)
{
   DRAWPATH;
   file << "copypage\n";
   file << "erasepage\n";

   framecount++;
   if (MULTIFRAME==1) {
    file << "showpage\n";
    file << '\n';
    file.close();
    PSheader(nstep);
   }

}

void PSfile::movePen(Real x, Real y)
{
   DRAWPATH;
   xcur = int( x*xfactor );
   ycur = int( y*yfactor );
   npts++;
}

void PSfile::drawLine(Real x, Real y, int lev)
{
   if (npts == 1) {
      file << xcur << ' ' << ycur << " M\n";
   };
   xcur = int( x*xfactor );
   ycur = int( y*yfactor );
   if (npts > 0) {
      file << xcur << ' ' << ycur << " L\n";
      if(lev == 0) {
        file << "[]" << " 0 setdash " << "\n";
      } else {
//        file << "[" << 9*(4-lev) << "]" << " 0 setdash " << "\n";
        if(lev == 1) {
           file << "[" << 2*(3-lev) << " " << 5*(3-lev) << "]"
                << " 0 setdash " << "\n";
        } else {
          if(lev == 2) {
              file << "[" << 1*(3-lev) << " " << 5*(3-lev) << "]"
                << " 0 setdash " << "\n";
          } else {
		  std::cout << " drawLine handles 2 levels of refinement "
                 << " You lose... " << "\n";
            exit(0);
          }
	}
      }
   } else {
      file << xcur << ' ' << ycur << " M\n";
   };
   npts++;
}

void PSfile::setLineWidth(int lw)
{
   if (npts > 0) {
      DRAWPATH;
      npts = 1;
   };
   file << lw << " setlinewidth\n";
}

#undef DRAWPATH
#endif

