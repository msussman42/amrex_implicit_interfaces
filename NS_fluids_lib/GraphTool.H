#ifndef _GraphTool_H_
#define _GraphTool_H_  

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

#include <CArena.H>
#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>

#if (HAS_WINDOWS)
#include "XWindow.H"
#endif
#include "PSfile.H"
#include "IFrame.H"

typedef unsigned short GTDevice;
const GTDevice noDevice = 0;
const GTDevice xWinDevice = 1;
const GTDevice psDevice = 2;

// This class is a C++ version of Louis Howell's GraphTools library.
// It has the ability to do simple line drawing graphics.
// The devices supported are X-windows and postscript files and 
// an object of this class can have either one or both active at the
// same time.
// Many of the operations, such as querring the mouse and setting color
// and fonts are not implemented for the postscript device and
// become no-ops.
// The user can change the currently active devices at any time
// through a variety of member functions and all subsequent
// graphics operations will only be directed to those devices.
// One sets the active devices by ORing together any combination
// of the above constants (noDevice, XWinDevice, psDevice).
class GraphTool {
   public:
      GraphTool(Real x1, Real y1, Real x2, Real y2, 
                const char *str, int wid, int high,
		GTDevice device = xWinDevice);
      // create a graphtool:
      //   (x1,y1) X (x2,y2) physical coordinates of window
      //   str               window/postscript file name
      //   (wid,high)        width and height (in pixels) of the window
      //   device            selected devices           

      GraphTool(const Real *lo_pt, const Real *hi_pt,
                const char *str, int wid, int high,
		GTDevice device = xWinDevice);
      // create a graphtool:
      //   lo_pt, hi_pt      arrays specifying physical coordinates of window
      //   str               window/postscript file name
      //   (wid,high)        width and height (in pixels) of the window
      //   device            selected devices           

      GraphTool(const Box &bx,
                const char *str, int wid, int high,
		GTDevice device = xWinDevice);
      // create a graphtool:
      //   bx	             physical coordinates specified via box
      //   str               window/postscript file name
      //   (wid,high)        width and height (in pixels) of the window
      //   device            selected devices           

      GraphTool(int nstep,Real x1, Real y1, Real x2, Real y2,
                const char *str, int maxwinsize = 600,
		GTDevice  device = xWinDevice);
      // greate a graphtool with correct aspec ratio and
      // maximum length having maxwinsize pixels.

      GraphTool(const Real *lo_pt, const Real *hi_pt,
                const char *str, int maxwinsize = 600,
		GTDevice  device = xWinDevice);
      // greate a graphtool with correct aspec ratio and
      // maximum length having maxwinsize pixels.

      GraphTool(const Box &bx,
                const char *str, int maxwinsize = 600,
                GTDevice  device = xWinDevice);
      // greate a graphtool with correct aspec ratio and
      // maximum length having maxwinsize pixels.

      ~GraphTool();
      // distructor

      GraphTool& setDevice(GTDevice device);
      // set these as the active devices

      GraphTool& setFrame(const IFrame &fr);
      // set the current int to Real mapping

      GTDevice getDevice() const;
      // returns currently active devices

      IFrame getFrame() const;
      // returns currently active int to Real mapping

      GraphTool& addDevice(GTDevice device);
      // add these devices to the currently active devices

      GraphTool& rmDevice(GTDevice device);
      // remove these devices from the currently active devices

      GraphTool& newPage(int nstep);
      // clear XWindow if active
      // end (not start new)  page of PostScript file if active

      GraphTool& movePen(Real x, Real y);
      // move to given location without drawing line


      GraphTool& movePen(const IntVect &v);
      // move to given location without drawing line

      GraphTool& drawLine(Real x, Real y, int lev = 0);
      // draw line from current location to this location


      GraphTool& drawLine(const IntVect &v);
      // draw line from current location to this location

      GraphTool& drawBox(Real xlo, Real ylo, Real xhi, Real yhi,
		int lev=0);
      // draw a rectangle with the given diagional locations

      GraphTool& drawBox(const Box &bx);
      // draw a rectangle with the given diagional locations

      GraphTool& setClipRegion(Real xlo, Real ylo,
                               Real xhi, Real yhi);
      // set the clipping region so that nothing is drawn outside
      // this this box.


      GraphTool& setClipRegion(const Box &bx);
      // set the clipping region so that nothing is drawn outside
      // this this box.

      GraphTool& setLineWidth(int lw);
      // change the width of the lines drawn.
      
      GraphTool& setFont(char *);
      // change the font to the given XWindow font.
      // NoOp for PostScript device.
      
      GraphTool& putString(Real x, Real y, const char *str);
      // output the given string at the given location
      // NoOp for PostScript device.

      GraphTool& putString(const IntVect &v, const char *str);
      // output the given string at the given location
      // NoOp for PostScript device.

      GraphTool& setfgColor(const char *colorname);
      // change the XWindow foreground color to the named color.
      // NoOp for PostScript device.

      GraphTool& setbgColor(const char *colorname);
      // change the XWindow background color to the named color.
      // NOTE: will also blank the screen!
      // NoOp for PostScript device.

      GraphTool& defineCmap( unsigned short *red, unsigned short *green,
                             unsigned short *blue, int num);
      // NoOp for PostScript device.

      GraphTool& setfgColor(int colorTableEntry);
      // change the XWindow foreground color to that specified
      // by the colorTableEntry.
      // NoOp for PostScript device.

      GraphTool& setbgColor(int colorTableEntry);
      // change the XWindow background color to that specified
      // by the colorTableEntry.
      // NOTE: will also blank the screen!
      // NoOp for PostScript device.

      int        getMouse(Real &x, Real &y) const;
      // query the mouse and return the number of the button
      // pushed and the (x,y) coordinates (in physical coordinates)
      // of where the button was pressed.
      // NoOp for PostScript device.

      int        getMouse(IntVect &v) const;
      // query the mouse and return the number of the button
      // pushed and the position (mapped Real to IntVect)
      // of where the button was pressed.
      // NoOp for PostScript device.

   protected:
      int       width, height;      // width and height of window
      Real    xlo,ylo,xhi,yhi;    // physical domain
      Real    x0, y0;             // current position
      Real    xcL,xcU,ycL,ycU;    // clipping coordinates
      int       was_out;            // was last point out of range?
      GTDevice  dev;                // encoding of active devices
      char      *name;		    // name of window/file
#if (HAS_WINDOWS)
      XWindow   *win;               // pointer to XWindow object
#endif
      PSfile    *ps;		    // pointer to postscript object

      IFrame  frame;		    // current int to Real mapping

      Real    clipRatio(Real,Real,Real,Real);
      void      initGT(int,int,Real,Real,Real,Real,
                       const char*,GTDevice,int nstep=0);

};

#endif
