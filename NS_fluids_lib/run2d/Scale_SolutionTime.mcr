#!MC 1410

$!EXTENDEDCOMMAND COMMANDPROCESSORID='extend time mcr'
	COMMAND='QUERY.NUMTIMESTEPS NUMTIMESTEPS'
$!AlterData
  Equation = '{TSTAR}=0.0'

$!LOOP |NUMTIMESTEPS|
	$!EXTENDEDCOMMAND COMMANDPROCESSORID='extend time mcr'
		COMMAND='SET.CURTIMESTEP |LOOP|'
	$!EXTENDEDCOMMAND COMMANDPROCESSORID='extend time mcr'
		COMMAND='QUERY.TIMEATSTEP |LOOP| CURTIME'
	$!VARSET|LOCALTIME|=(|CURTIME|/0.0005664)

	$!EXTENDEDCOMMAND COMMANDPROCESSORID='extendmcr'
		COMMAND='QUERY.ACTIVEZONES ZNUM'
	$!AlterData  [|ZNUM|]
	  Equation = '{TSTAR}=|LOCALTIME|'
$!ENDLOOP

$!GETVARNUMBYNAME |numTSTAR|
	NAME = "TSTAR"

$!AttachText 
  AnchorPos
    {
    X = 10
    Y = 90
    }
  TextShape
    {
    IsBold = No
    }
  Text = 'Time: &(MAXVAR[|numTSTAR|])'