#!/bin/bash
for file in FABRIC_DROP.F90 FABRIC_DROP_MEHDI.F90 GENERAL_PHASE_CHANGE.F90 ICE_ON_SUBSTRATE.F90 MEHDI_SHOCK_DROP.F90 MITSUHIRO_MELTING.F90 SIMPLE_KASSEMI.F90 SIMPLE_PALMORE_DESJARDINS.F90 YAOHONG_INKJET.F90 ZEYU_droplet_impact.F90 flexible_plate_impact.F90 helix.F90 rigid_FSI.F90 sample_user_defined.F90 sinking_particle.F90 thermalspray.F90 wavy_channel.F90
do
echo "Processing $file"
sed 's/ibase+2/ibase+ENUM_TEMPERATUREVAR+1/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
