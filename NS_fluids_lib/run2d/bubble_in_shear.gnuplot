plot "critical_Reynolds_Number"
replot "138.3*exp(-1.41*x)" with lines
set terminal png
set output "Bubble_Deformation_Chart.png"
replot

