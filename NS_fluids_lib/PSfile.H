#ifndef _PSfile_H_
#define _PSfile_H_ 

#include <winstd.H>

#include <fstream>
#include <algorithm>
#include <vector>

#if defined(BL_OLD_STL)
#include <stdio.h>
#include <math.h>
#else
#include <cstdio>
#include <cmath>
#endif

#include <CArena.H>
#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>

 
class PSfile {
   public:
      PSfile(int wid, int high, const char *name,int nstep);
      ~PSfile();

      void PSheader(int nstep);
      void newPage(int nstep);
      void movePen(Real x, Real y);
      void drawLine(Real x, Real y, int lev=0);
      void setLineWidth(int lw);
      void flush();

   private:
      int      save_wid,save_high;
      int      framecount;
      int      npts;
      int      xfactor;
      int      yfactor;
      int      xcur;
      int      ycur;
      std::ofstream file;
      char     *filename;
      char     *save_name;

}; 

#endif

