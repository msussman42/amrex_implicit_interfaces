#!MC 1410
$!LOOP 99
$!VarSet |MFBD| = '/home/mark/geometrySPECTRAL/iamrlib/run2d/nucleate'
$!VarSet |NDD| = '/nddata00'
$!VarSet |ZZ| = '00'
$!VarSet |PLT| = '.plt'
$!VarSet |FileBase| = '|MFBD||NDD||LOOP%02d||ZZ||PLT|'
$!READDATASET '|FileBase|'
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  VARNAMELIST = '"X" "Y" "U01" "V01" "PMG" "PEOS" "DIV" "DIVDT" "MCH" "F01" "F02" "F03" "L0101" "L0202" "L0303" "L0102" "L0103" "L0203" "NX0101" "NY0101" "NX0202" "NY0202" "NX0303" "NY0303" "NX0102" "NY0102" "NX0103" "NY0103" "NX0203" "NY0203" "D01" "T01" "D02" "T02" "D03" "T03" "MU01" "MU02" "MU03" "DT01" "TR01" "TRT01" "TRTF01" "VORT01" "DT02" "TR02" "TRT02" "TRTF02" "VORT02" "DT03" "TR03" "TRT03" "TRTF03" "VORT03"'
$!RemoveVar |MFBD|
$!RemoveVar |NDD|
$!RemoveVar |ZZ|
$!RemoveVar |PLT|
$!RemoveVar |FileBase|
$!EndLOOP
