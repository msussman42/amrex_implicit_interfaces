#!MC 1400
# Created by Tecplot 360 build 14.0.1.26249
$!VarSet |MFBD| = '/home/mof/geometry6.0MOF/iamrlib/run2d/testbubble'
$!READDATASET  '"|MFBD|/nddata000000.tec" "|MFBD|/nddata000010.tec" "|MFBD|/nddata000020.tec" "|MFBD|/nddata000030.tec" "|MFBD|/nddata000040.tec" "|MFBD|/nddata000050.tec" "|MFBD|/nddata000060.tec" "|MFBD|/nddata000070.tec" "|MFBD|/nddata000080.tec" "|MFBD|/nddata000090.tec" "|MFBD|/nddata000100.tec" "|MFBD|/nddata000110.tec" "|MFBD|/nddata000120.tec" "|MFBD|/nddata000130.tec" "|MFBD|/nddata000140.tec" "|MFBD|/nddata000150.tec" "|MFBD|/nddata000160.tec" "|MFBD|/nddata000170.tec" "|MFBD|/nddata000180.tec" "|MFBD|/nddata000190.tec" "|MFBD|/nddata000200.tec" "|MFBD|/nddata000210.tec" "|MFBD|/nddata000220.tec" "|MFBD|/nddata000230.tec" "|MFBD|/nddata000240.tec" "|MFBD|/nddata000250.tec" "|MFBD|/nddata000260.tec" "|MFBD|/nddata000270.tec" "|MFBD|/nddata000280.tec" "|MFBD|/nddata000290.tec" "|MFBD|/nddata000300.tec" "|MFBD|/nddata000310.tec" "|MFBD|/nddata000320.tec" "|MFBD|/nddata000330.tec" "|MFBD|/nddata000340.tec" "|MFBD|/nddata000350.tec" "|MFBD|/nddata000360.tec" "|MFBD|/nddata000370.tec" "|MFBD|/nddata000380.tec" "|MFBD|/nddata000390.tec" "|MFBD|/nddata000400.tec" "|MFBD|/nddata000410.tec" "|MFBD|/nddata000420.tec" "|MFBD|/nddata000430.tec" "|MFBD|/nddata000440.tec" "|MFBD|/nddata000450.tec" "|MFBD|/nddata000460.tec" "|MFBD|/nddata000470.tec" "|MFBD|/nddata000480.tec" "|MFBD|/nddata000490.tec" "|MFBD|/nddata000500.tec" "|MFBD|/nddata000510.tec" "|MFBD|/nddata000520.tec" "|MFBD|/nddata000530.tec" "|MFBD|/nddata000540.tec" "|MFBD|/nddata000550.tec" "|MFBD|/nddata000560.tec" "|MFBD|/nddata000570.tec" "|MFBD|/nddata000580.tec" "|MFBD|/nddata000590.tec" "|MFBD|/nddata000600.tec" "|MFBD|/nddata000610.tec" "|MFBD|/nddata000620.tec" "|MFBD|/nddata000630.tec" "|MFBD|/nddata000640.tec" "|MFBD|/nddata000650.tec" "|MFBD|/nddata000660.tec" "|MFBD|/nddata000670.tec" "|MFBD|/nddata000680.tec" "|MFBD|/nddata000690.tec" "|MFBD|/nddata000700.tec" "|MFBD|/nddata000710.tec" "|MFBD|/nddata000716.tec" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"X" "Y" "U" "V" "U2" "V2" "vof1" "vof2" "vof3" "vof4" "lsdist1" "lsdist2" "lsdist3" "lsdist4" "pres" "den1" "den2" "den3" "den4"'
$!GLOBALCONTOUR 1  VAR = 11
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!FIELDLAYERS SHOWCONTOUR = YES
$!GLOBALCONTOUR 1  VAR = 7
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!CONTOURLEVELS NEW
  CONTOURGROUP = 1
  RAWDATA
1
0.5
$!REDRAWALL 
$!EXPORTSETUP EXPORTFORMAT = FLASH
$!EXPORTSETUP IMAGEWIDTH = 792
$!EXPORTSETUP EXPORTFNAME = '/home/mof/geometry6.0MOF/iamrlib/run2d/testbubble/export.swf'
$!ANIMATETIME 
  STARTTIME = 0
  ENDTIME = 10
  SKIP = 1
  CREATEMOVIEFILE = YES
  LIMITSCREENSPEED = NO
  MAXSCREENSPEED = 12
$!RemoveVar |MFBD|
