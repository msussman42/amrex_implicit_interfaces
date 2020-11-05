#!/bin/bash
echo "converting to .vtp"
python amrex_particles_to_vtp.py 00000 00200 InnerNarrowBand 50 &
python amrex_particles_to_vtp.py 00000 00200 OuterNarrowBand 50 &
python amrex_particles_to_vtp.py 00000 00200 InterfaceParticles 50 &
wait
echo "creating VisIt sequence"
ls -1 InnerNarrowBand*.vtp | tee movie_innerNB.visit &
ls -1 OuterNarrowBand*.vtp | tee movie_outerNB.visit &
ls -1 InterfaceParticles*.vtp | tee movie_interface.visit &
wait
echo "done"