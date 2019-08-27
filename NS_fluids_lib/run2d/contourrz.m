finedata=load('visual00001');
D=size(finedata)
if (D(2)~=10) 
 sprintf('invalid number of columns\n')
end

numrows=D(1)
rimax=0;
rjmax=0;
for I=1:numrows,
 if (finedata(I,1)>rimax) 
  rimax=finedata(I,1);
 end
 if (finedata(I,2)>rjmax) 
  rjmax=finedata(I,2);
 end
end
imax=round(rimax)+1
jmax=round(rjmax)+1

imax2=2*imax
volume1=zeros(jmax,imax2);
U=zeros(jmax,imax2);
V=zeros(jmax,imax2);
P=zeros(jmax,imax2);
viscosity=zeros(jmax,imax2);
shear=zeros(jmax,imax2);
 
for I=1:numrows,
 II=round(finedata(I,1))+1;
 JJ=round(finedata(I,2))+1;
 U(JJ,II+imax)=finedata(I,3);
 V(JJ,II+imax)=finedata(I,4);
 P(JJ,II+imax)=finedata(I,8);
 viscosity(JJ,II+imax)=finedata(I,9);
 shear(JJ,II+imax)=finedata(I,10);
 volume1(JJ,II+imax)=finedata(I,5);
 U(JJ,imax-II+1)=-finedata(I,3);
 V(JJ,imax-II+1)=finedata(I,4);
 P(JJ,imax-II+1)=finedata(I,8);
 viscosity(JJ,imax-II+1)=finedata(I,9);
 shear(JJ,imax-II+1)=finedata(I,10);
 volume1(JJ,imax-II+1)=finedata(I,5);
end
hold off;
contour(volume1,[0]);
hold on;
contour(P);
if (1==0)
 hold off;
 contour(shear);
end
