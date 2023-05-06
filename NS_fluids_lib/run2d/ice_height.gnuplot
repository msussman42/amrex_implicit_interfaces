set xlabel "time"
set ylabel "centerline ice height"
plot "32x32LiquidIceSLICE" with lines
replot "64x64LiquidIceSLICE" with lines
replot "128x128LiquidIceSLICE" with lines
set terminal eps
set output "ICE_HEIGHT.eps"
replot

