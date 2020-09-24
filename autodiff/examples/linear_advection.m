function linear_advection
  
Ncells=64
Nnodes=65
x_array=zeros(1,Nnodes);
uexact_array=zeros(1,Nnodes);
vn_array=zeros(1,Nnodes);
vnp1_array=zeros(1,Nnodes);

xlo=0.0
xhi=1.0
stop_time=0.5
h=(xhi-xlo)/Ncells
a=1.0
k=0.5*h/a

for i=0:Ncells
	xi=xlo+i*h;
	x_array(i+1)=xi;
        vnp1_array(i+1)=u0(xi);
        uexact_array(i+1)=uexact(xi-stop_time*a);
end
disp('plotting initial conditions')
plot(x_array,vnp1_array)

current_time=0.0
nsteps=0

while (current_time<stop_time)

 for i=0:Ncells
  vn_array(i+1)=vnp1_array(i+1);
 end

 for i=1:Ncells
  DM=vn_array(i+1)-vn_array(i);
  vnp1_array(i+1)=vn_array(i+1)-(k/h)*DM;
 end
 vnp1_array(1)=vnp1_array(Ncells+1);

 nsteps=nsteps+1
 if (current_time+k>=stop_time-1.0e-10)
  k=stop_time-current_time+1.0e-10
 end
 current_time=current_time+k
 nsteps=nsteps+1

end

plot(x_array,uexact_array,x_array,vnp1_array)


endfunction

function u0x=u0(x)

 ulocal=0.0;
 if (x<0.5)
	 ulocal=-1.0;
 elseif (x>=0.5)
	 ulocal=1.0;
 else
	 disp('x became corrupt')
	 x
	 quit
 end
 u0x=ulocal;
end


function uexactx=uexact(x)

 xlocal=x;
 while (xlocal<0.0)
  xlocal=xlocal+1.0;
 end

 while (1.0-xlocal<0.0)
  xlocal=xlocal-1.0;
 end

 uexactx=u0(xlocal);

end

