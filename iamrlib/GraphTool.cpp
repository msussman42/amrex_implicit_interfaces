#if (HAS_WINDOWS)
#include "GraphTool.H"

GraphTool& GraphTool::setFrame(const IFrame &fr)
{
   frame = fr;
   return *this;
}

GTDevice GraphTool::getDevice() const 
{
   return dev;
}

IFrame GraphTool::getFrame() const
{
   return frame;
}

GraphTool& GraphTool::movePen(const IntVect &v)
{
   return movePen(frame.toX(v),frame.toY(v));
}

GraphTool& GraphTool::drawLine(const IntVect &v)
{
   return drawLine(frame.toX(v),frame.toY(v));
}


GraphTool& GraphTool::drawBox(const Box &b)
{
   Real x1 = frame.toX(b.smallEnd());
   Real x2 = frame.toX(b.bigEnd());
   Real y1 = frame.toY(b.smallEnd());
   Real y2 = frame.toY(b.bigEnd());
   return drawBox(x1,y1,x2,y2);
}

GraphTool& GraphTool::putString(const IntVect &v, const char *str)
{
   return putString(frame.toX(v),frame.toY(v),str);
}

#define XTOINT(x) int( width* ((x-xlo)/(xhi-xlo)) )
#define YTOINT(y) int( height*((y-ylo)/(yhi-ylo)) )
#define PS_XSCALE(x) ((x)-xlo)/(xhi-xlo)
#define PS_YSCALE(y) ((y)-ylo)/(yhi-ylo)

// private function to init a GraphTool
void GraphTool::initGT(int wid, int high, Real x1, Real y1,
                  Real x2, Real y2, const char *str, GTDevice device,
                  int nstep)
{
   name = new char[strlen(str) + 1];
   strcpy(name,str);
   dev = device;
   if ( (x1 >= x2) || (y1 >= y2) ) {
	   std::cerr << "GraphTool: invalid domain" << '\n';
      abort();
   };
   xcL = xlo = x1;     ycL = ylo = y1;
   xcU = xhi = x2;     ycU = yhi = y2;
   x0 = xlo;     y0 = ylo;
   was_out = 0;
   width = wid;
   height = high;

#if (BL_SPACEDIM == 2)
   IntVect v_lo(0,0);
   IntVect v_hi(wid-1,high-1);
#else
   IntVect v_lo(0,0,0);
   IntVect v_hi(wid-1,high-1,0);
#endif
   Box     bx(v_lo,v_hi);
   IFrame  fr(bx,xlo,xhi,ylo,yhi);
   frame = fr;


#if (HAS_WINDOWS)
   win = NULL;
   ps  = NULL;
   if ((dev&xWinDevice) == xWinDevice) {
      win = new XWindow(width,height,str);
   };
   if ((dev&psDevice) == psDevice) {
      ps = new PSfile(width,height,str,nstep);
   };
#endif
   
}

GraphTool::GraphTool(const Box &bx,
                     const char *str, int wid, int high, GTDevice device)
{
   const int* lo_d = bx.smallEnd().getVect();
   const int* hi_d = bx.bigEnd().getVect();
   Real x1 = (Real) lo_d[0];
   Real x2 = (Real) hi_d[0];
   Real y1 = (Real) lo_d[1];
   Real y2 = (Real) hi_d[1];
   initGT(wid,high,x1,y1,x2,y2,str,device);
}

GraphTool::GraphTool(int nstep,Real x1, Real y1, Real x2, Real y2,
                     const char *str, int maxwinsize,
		     GTDevice  device)
{
   Real maxlen = std::max(x2-x1,y2-y1);
   int    wid    = int( maxwinsize*(x2-x1)/maxlen );
   int    high   = int( maxwinsize*(y2-y1)/maxlen );
   initGT(wid,high,x1,y1,x2,y2,str,device,nstep);
}		     

GraphTool::GraphTool(Real x1, Real y1, Real x2, Real y2,
                     const char *str, int wid, int high, GTDevice device)
{
   initGT(wid,high,x1,y1,x2,y2,str,device);
}		     

GraphTool::GraphTool(const Real *lo_pt, const Real *hi_pt,
                     const char *str, int wid, int high, GTDevice device)
{
   initGT(wid,high,lo_pt[0],lo_pt[1],hi_pt[0],hi_pt[1],str,device);
}	  
		

GraphTool::GraphTool(const Real *lo_pt, const Real *hi_pt,
                     const char *str, int maxwinsize,
		     GTDevice  device)
{
   Real x1 = lo_pt[0];
   Real x2 = hi_pt[0];
   Real y1 = lo_pt[1];
   Real y2 = hi_pt[1];
   Real maxlen = std::max(x2-x1,y2-y1);
   int    wid    = int( maxwinsize*(x2-x1)/maxlen );
   int    high   = int( maxwinsize*(y2-y1)/maxlen );
   initGT(wid,high,x1,y1,x2,y2,str,device);
}


GraphTool::GraphTool(const Box &bx,
                     const char *str, int maxwinsize,
		     GTDevice  device)
{
   const int* lo_d = bx.smallEnd().getVect();
   const int* hi_d = bx.bigEnd().getVect();
   Real x1 = (Real) lo_d[0];
   Real x2 = (Real) hi_d[0];
   Real y1 = (Real) lo_d[1];
   Real y2 = (Real) hi_d[1];
   Real maxlen = std::max(x2-x1,y2-y1);
   int    wid    = int( maxwinsize*(x2-x1)/maxlen );
   int    high   = int( maxwinsize*(y2-y1)/maxlen );
   initGT(wid,high,x1,y1,x2,y2,str,device);
   IFrame dummy(bx,x1,x2,y1,y2);
   frame = dummy;
}

GraphTool::~GraphTool()
{
   delete name;
#if (HAS_WINDOWS)
   delete win;
   delete ps;
#endif
}

GraphTool& GraphTool::rmDevice(GTDevice device)
{
  if (dev & device) dev ^= device;
  return *this;
}

GraphTool& GraphTool::setDevice(GTDevice device)
{
   dev = device;
#if (HAS_WINDOWS)
   if ( ((dev&xWinDevice) == xWinDevice) && (win == NULL) ) {
      win = new XWindow(width,height,name);
   };
   if ( ((dev&psDevice) == psDevice) && (ps == NULL) ) {
      ps = new PSfile(width,height,name,0);
   };
#endif
   return *this;
}

GraphTool& GraphTool::addDevice(GTDevice device)
{
   dev |= device;
#if (HAS_WINDOWS)
   if ( ((dev&xWinDevice) == xWinDevice) && (win == NULL) ) {
      win = new XWindow(width,height,name);
   };
   if ( ((dev&psDevice) == psDevice) && (ps == NULL) ) {
      ps = new PSfile(width,height,name,0);
   };
#endif
   return *this;
}

GraphTool& GraphTool::newPage(int nstep)
{
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
      win->newPage();
   };
   if ((dev&psDevice) == psDevice) {
      ps->newPage(nstep);
   };
#endif
   return *this;
}

GraphTool& GraphTool::movePen(Real x, Real y)
{
   int out = ( (x<xcL) || (x>xcU) || (y<ycL) || (y>ycU) );
   x0 = x;
   y0 = y;
   was_out = out;
   if (!out) {
#if (HAS_WINDOWS)
      if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(x);
         int yi = YTOINT(y);
         win->movePen(xi,yi);
      };
      if ((dev&psDevice) == psDevice) {
         ps->movePen(PS_XSCALE(x),PS_YSCALE(y));
      };
#endif
   };
   return *this;
}


GraphTool& GraphTool::setClipRegion(const Box &b)
{
   Real x1 = frame.toX(b.smallEnd());
   Real x2 = frame.toX(b.bigEnd());
   Real y1 = frame.toY(b.smallEnd());
   Real y2 = frame.toY(b.bigEnd());
   return setClipRegion(x1,y1,x2,y2);
}
                                    
GraphTool& GraphTool::drawLine(Real x, Real y, int lev)
{
   int out = ( (x<xcL) || (x>xcU) || (y<ycL) || (y>ycU) );
   if ( (!out) && (!was_out) ) {
      // old point and new point are in clip region
#if (HAS_WINDOWS)
      if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(x);
         int yi = YTOINT(y);
         win->drawLine(xi,yi);
      };
      if ((dev&psDevice) == psDevice) {
         ps->drawLine(PS_XSCALE(x),PS_YSCALE(y), lev);
      };
#endif
   } else if ((!was_out) && out) {
      // old point was in range, new point out of range
      // find exit point and only draw to there
      Real c = clipRatio(x0,y0,x,y);
      Real xc = x0 + c*(x-x0);
      Real yc = y0 + c*(y-y0);
#if (HAS_WINDOWS)
      if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(xc);
         int yi = YTOINT(yc); 
         win->drawLine(xi,yi);
      };
      if ((dev&psDevice) == psDevice) {
         ps->drawLine(PS_XSCALE(xc),PS_YSCALE(yc), lev);
      };
#endif
      
   } else if (was_out && (!out)) {
      // old point was out of range, new point is in range
      // find entry point and only draw from there
      Real c = clipRatio(x,y,x0,y0);
      Real xc = x + c*(x0-x);
      Real yc = y + c*(y0-y);
#if (HAS_WINDOWS)
      if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(xc);
         int yi = YTOINT(yc); 
         win->movePen(xi,yi);
         xi = XTOINT(x);
         yi = YTOINT(y); 
	 win->drawLine(xi,yi);
      };
      if ((dev&psDevice) == psDevice) {
         ps->movePen(PS_XSCALE(xc),PS_YSCALE(yc));
         ps->drawLine(PS_XSCALE(x),PS_YSCALE(y), lev);
      };
#endif
   } else {
      // both points outside range, does line intersect
      // box at all?  if so, draw on intersection
   };
   x0 = x;
   y0 = y;
   was_out = out;
   return *this;
}

GraphTool& GraphTool::drawBox(Real x1, Real y1,
                              Real x2, Real y2, int lev) 
{
   movePen(x1,y1);
   drawLine(x2,y1,lev);
   drawLine(x2,y2,lev);
   drawLine(x1,y2,lev);
   drawLine(x1,y1,lev);
   return *this;
}

GraphTool& GraphTool::setClipRegion(Real x1, Real y1,
                                    Real x2, Real y2)
{
   if ( (x1 >= x2) || (y1 >= y2) ) {
	   std::cout << "Invalid Clip region, ignoring..." << '\n';
   } else {
      xcL = x1;
      xcU = x2;
      ycL = y1;
      ycU = y2;
   };
   return *this;
}

GraphTool& GraphTool::setLineWidth(int lw)
{
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
      win->setLineWidth(lw);
   };
   if ((dev&psDevice) == psDevice) {
      ps->setLineWidth(lw);
   };
#endif
   return *this;
}

GraphTool& GraphTool::setFont(char *font_name)
{
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
      win->setFont(font_name);
   };
   if ((dev&psDevice) == psDevice) {
   };
#endif
   return *this;
}

GraphTool& GraphTool::putString(Real x, Real y, const char *str)
{
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(x);
         int yi = YTOINT(y);
         win->putString(xi,yi,str);
   };
   if ((dev&psDevice) == psDevice) {
   };
#endif
   return *this;
}

GraphTool& GraphTool::setfgColor(const char* color)
{
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
      win->setfgColor(color);
   };
   if ((dev&psDevice) == psDevice) {
   };
#endif
   return *this;
}

GraphTool& GraphTool::setbgColor(const char* color)
{
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
      win->setbgColor(color);
   };
   if ((dev&psDevice) == psDevice) {
   };
#endif
   return *this;
}

GraphTool& GraphTool::setfgColor(int color)
{
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
      win->setfgColor(color);
   };
   if ((dev&psDevice) == psDevice) {
   };
#endif
   return *this;
}

GraphTool& GraphTool::setbgColor(int color)
{
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
      win->setbgColor(color);
   };
   if ((dev&psDevice) == psDevice) {
   };
#endif
   return *this;
}

GraphTool& GraphTool::defineCmap(unsigned short *red, unsigned short *green,
                                 unsigned short *blue, int num)
{
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
      win->defineCmap(red, green, blue, num);
   };
   if ((dev&psDevice) == psDevice) {
   };
#endif
   return *this;
}


int GraphTool::getMouse(IntVect &v) const
{
   Real x[2];
   int but = getMouse(x[0],x[1]);
   v = frame.toIV(x);
   return but;
}

int GraphTool::getMouse(Real &x, Real &y) const
{
   int button = -1;
#if (HAS_WINDOWS)
   if ((dev&xWinDevice) == xWinDevice) {
      int i,j;
      button = win->getMouse(i,j);
      x = xlo + (xhi-xlo)*( Real(i)/Real(width) );
      y = ylo + (yhi-ylo)*( Real(j)/Real(height) );
   };
#endif
   return button;
}

//  Private function: Determines fraction of distance from
//  (x1,y1) to (x2,y2) at which boundary is intersected
Real GraphTool::clipRatio(Real x1, Real y1, Real x2, Real y2)
{
   Real r = 1.0;
   if (x2<xcL) r = std::min(r,(x1-xcL)/(x1-x2));
   if (y2<ycL) r = std::min(r,(y1-ycL)/(y1-y2));
   if (x2>xcU) r = std::min(r,(xcU-x1)/(x2-x1));
   if (y2>ycU) r = std::min(r,(ycU-y1)/(y2-y1));
   return 0.9999*r;
}

#undef  XTOINT
#undef  YTOINT
#undef  PS_XSCALE
#undef  PS_YSCALE
#endif

