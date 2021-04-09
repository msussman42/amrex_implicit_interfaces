plot "aratio" with lines
replot "lev3/aratio" with lines
replot "lev4/aratio" with lines
replot "numerical_spherical_1D" using 2:3 with lines

plot "KASSEMI/aratio" with lines
replot "KASSEMI/sealed/aratio" with lines
replot "PD/aratio" with lines
replot "numerical_spherical_1D_KASSEMI" using 2:3 with lines
replot "numerical_spherical_1D_PD" using 2:3 with lines

plot "aratio_Borodulin_lev2" with lines
replot "aratio_Borodulin_lev3" with lines
replot "aratio_Borodulin_lev4" with lines
replot "Borodulin_KASSEMI_1D" using 2:3 with lines
replot (1000.0-x)/1000.0 with lines

plot "aratio_Villegas_KASSEMI_lev2" with lines
replot "aratio_Villegas_KASSEMI_lev2_sealed" with lines
replot "aratio_Villegas_PD_lev2" with lines
replot "Villegas_PD_1D" using 2:3 with lines
replot "Villegas_PD_analytical" using 1:3 with lines


