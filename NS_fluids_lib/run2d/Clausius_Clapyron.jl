# https://julialang.org/
# https://julialang.org/downloads/
#
# FOR UBUNTU LINUX:
# 1. wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/
#           julia-1.5.0-linux-x86_64.tar.gz
# 2. tar zxvf julia-1.5.0-linux-x86_64.tar.gz
# 3. place following line in the ~/.bashrc file:
#    export PATH="$PATH:/path/to/<Julia directory>/bin"
# 4. julia
# 5. julia> import Pkg;         # first time only
# 6a. julia> Pkg.add("Plots")    # first time only
# 6b. julia> Pkg.add("PyPlot")    # first time only
# 7. julia> include("taylor_series_simple.jl")
# 8. julia> exit()
#
# FOR MacBook:  (this approach is easy, bypasses Jupyter.  
# not recommended to use Jupyter unless
# you know it well)
# 
# 1. https://julialang.org/downloads/
# 2. LTS release
# 3. macOS, 64 bit, click to download
# 4. click on the dmg file in the Downloads tab in the lower right hand corner
#    to install.
# 5. Double Click on the Julia icon to open; it will fail, then go to 
#    "security and privacy" and click "Open Anyway" where it says
#    "Julia-1.0" was blocked ...
# 6. Double Click on the Julia icon to open; the computer will give
#    you more warnings, but go ahead and open it anyway (my computer
#    still works)
# 7. pwd()
# 8. cd("tex/MAD3703FALL2020")   (this is where my julia code resides)
# 9. include("taylor_series_simple.jl")

function Tgamma_from_Y(Y)
 
 WA=4.0e-3
 WV=2.01e-3
 R=8.31446261815324e+7  # ergs/(Kelvin mol)
 R=8.31446261815324  # J/(Kelvin mol)
 X=WA*Y/(WV*(1.0-Y)+WA*Y) 
 Tgamma=-log(X)*R/(L*WV)
 Tgamma=Tgamma+1.0/TSAT
 Tgamma=1.0/Tgamma
    
 return Tgamma
 
end 

function do_the_plot(PD_name,color_name)

tz_array=readdlm(PD_name,Float64)
Nplot_nodes=size(tz_array,1)
Nplot_intervals=Nplot_nodes-1
println("Nplot_nodes=",Nplot_nodes)
println("Nplot_intervals=",Nplot_intervals)

t_array=Array{Float64,1}(undef,Nplot_nodes)
y_array=Array{Float64,1}(undef,Nplot_nodes)
z_array=Array{Float64,1}(undef,Nplot_nodes)

a=0.0 
b=0.04
h=(b-a)/Nplot_intervals
l=0.11564438522003612
lambda=0.1
tphys=0.46733780447863266
xblob2=0.025

for iplot in 0:Nplot_intervals

 t_array[iplot+1]=tz_array[iplot+1,1]
 y_array[iplot+1]=exact_interface_location(t_array[iplot+1],l,lambda,tphys,xblob2)
 z_array[iplot+1]=tz_array[iplot+1,2]*10.0

end
plot(t_array,y_array,color="blue",linewidth=2.0,linestyle="--")
plot(t_array,z_array,color=color_name,linewidth=2.0,linestyle="--")

end



using PyPlot
using DelimitedFiles

do_the_plot("PD_compute","red")
do_the_plot("PD_compute1","green")
do_the_plot("PD_compute2","orange")
do_the_plot("PD_compute3","purple")
