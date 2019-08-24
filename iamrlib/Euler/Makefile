FC     = gfortran
#FC     = ifort
FFLAGS = -g -O -I. -I..
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib


#----------------------------------------------------------------------
# 	Default action is to build all needed components.
#----------------------------------------------------------------------

all: libne.a

HEADERS=grid.h initial.h                 \
exact.h mat_constants.h parameters.h wave_constants.h           \
rarFan.h referenceValues.h 
# HEADERS=../grid.h ../initial.h                 \
# exact.h mat_constants.h parameters.h wave_constants.h           \
# ../rarFan.h ../referenceValues.h 

#----------------------------------------------------------------------
# 	Newton code.
#----------------------------------------------------------------------

NEWTON = insideFan_eu.o inv2st_eu.o readEOSConst_eu.o maxTime_eu.o newton_eu.o \
onWaveCurve_eu.o riemann_eu.o rpSnapshot_eu.o \
slope2state_eu.o \
st2inv_eu.o \
stIndex_eu.o states_eu.o statesO_eu.o valueAndSlope_eu.o wave_eu.o wave_data_eu.o wvIndex_eu.o \
thermSt2snd_eu.o thermSt2invEntr_eu.o veljump2vel_eu.o slopeInFan2velNsnd_eu.o primExtSt2W_eu.o \
compPower_eu.o pre2veljumpExpan_eu.o pre2veljumpCompr_eu.o stPrimExt2stPrimExt_eu.o pre2densCompr_eu.o \
pre2densExpan_eu.o waveCD_eu.o setVacMidState_eu.o pre2absDpreDvelCompr_eu.o pre2absDpreDvelExpan_eu.o \
pre2relMassFluxCompr_eu.o setRegMidState_eu.o
#rho2DmomDrhoCompr_eu.o rho2DmomDrhoExpan_eu.o rho2momCompr_eu.o rho2momExpan_eu.o slope2rho_eu.o 
#st2W_eu.o st2char_eu.o st2st_eu.o

libne.a:${NEWTON}
	${ARCH} ${ARCHFLAGS} $@ ${NEWTON} 
	${RANLIB} $@ 

#----------------------------------------------------------------------
#----------------------------------------------------------------------

insideFan_eu.o:insideFan_eu.f  ${HEADERS}

inv2st_eu.o:inv2st_eu.f  ${HEADERS}

readEOSConst_eu.o:readEOSConst_eu.f  ${HEADERS}

maxTime_eu.o:maxTime_eu.f  ${HEADERS}

newton_eu.o:newton_eu.f  ${HEADERS}

onWaveCurve_eu.o:onWaveCurve_eu.f  ${HEADERS}

#rho2DmomDrhoCompr_eu.o:rho2DmomDrhoCompr_eu.f  ${HEADERS}

#rho2DmomDrhoExpan_eu.o:rho2DmomDrhoExpan_eu.f  ${HEADERS}

#rho2momCompr_eu.o:rho2momCompr_eu.f  ${HEADERS}

#rho2momExpan_eu.o:rho2momExpan_eu.f  ${HEADERS}

#slope2rho_eu.o:slope2rho_eu.f  ${HEADERS}

#st2W_eu.o:st2W_eu.f  ${HEADERS}

#st2char_eu.o:st2char_eu.f  ${HEADERS}

#st2st_eu.o:st2st_eu.f  ${HEADERS}

riemann_eu.o:riemann_eu.f  ${HEADERS}

rpSnapshot_eu.o:rpSnapshot_eu.f  ${HEADERS}

slope2state_eu.o:slope2state_eu.f  ${HEADERS}

st2inv_eu.o:st2inv_eu.f  ${HEADERS}

stIndex_eu.o:stIndex_eu.f ${HEADERS}

states_eu.o:states_eu.f  ${HEADERS}

statesO_eu.o:statesO_eu.f  ${HEADERS}

valueAndSlope_eu.o:valueAndSlope_eu.f  ${HEADERS}

wave_eu.o:wave_eu.f  ${HEADERS}

wave_data_eu.o:wave_data_eu.f  ${HEADERS}

wvIndex_eu.o:wvIndex_eu.f ${HEADERS}

compPower_eu.o: compPower_eu.f  ${HEADERS}

pre2veljumpExpan_eu.o: pre2veljumpExpan_eu.f ${HEADERS}

pre2veljumpCompr_eu.o: pre2veljumpCompr_eu.f ${HEADERS}

stPrimExt2stPrimExt_eu.o: stPrimExt2stPrimExt_eu.f ${HEADERS}

pre2densCompr_eu.o: pre2densCompr_eu.f ${HEADERS}

pre2densExpan_eu.o: pre2densExpan_eu.f ${HEADERS}

waveCD_eu.o: waveCD_eu.f ${HEADERS}

setVacMidState_eu.o: setVacMidState_eu.f ${HEADERS}

pre2absDpreDvelCompr_eu.o: pre2absDpreDvelCompr_eu.f ${HEADERS}

pre2absDpreDvelExpan_eu.o: pre2absDpreDvelExpan_eu.f ${HEADERS}

pre2relMassFluxCompr_eu.o: pre2relMassFluxCompr_eu.f ${HEADERS}

pre2relMassFluxCompr_eu.o: pre2relMassFluxCompr_eu.f ${HEADERS}

setRegMidState_eu.o: setRegMidState_eu.f ${HEADERS}

#-------------------------------------------------------

compAver_eu.o:compAver_eu.f  ${HEADERS}

compInitialLevel_eu.o:compInitialLevel_eu.f  ${HEADERS}

nonDimExactPars_eu.o:nonDimExactPars_eu.f  ${HEADERS}

readInitialLevel_eu.o:readInitialLevel_eu.f  ${HEADERS}

setExactPars_eu.o:setExactPars_eu.f  ${HEADERS} 

#-------------------------------------------------------

compAver_RF_eu.o:compAver_RF_eu.f  ${HEADERS}

compInitialLevel_RF_eu.o:compInitialLevel_RF_eu.f  ${HEADERS}

nonDimExactPars_RF_eu.o:nonDimExactPars_RF_eu.f  ${HEADERS}

readInitialLevel_RF_eu.o:readInitialLevel_RF_eu.f  ${HEADERS}

setExactPars_RF_eu.o:setExactPars_RF_eu.f  ${HEADERS} 

#-------------------------------------------------------


.PHONY: clean all
clean: 
	rm -f *_eu.o *.out *.dat core
