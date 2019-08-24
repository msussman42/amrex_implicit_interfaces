%function [m1,m2]=contourxyz(k)
function [m1,m2]=conframe(k)

if k<10,
   finedata=load(['visual0000',int2str(k)]);
elseif k<100
     finedata=load(['visual000',int2str(k)]);
else
     finedata=load(['visual00',int2str(k)]);
end

D=size(finedata)
if (D(2)~=12) 
 sprintf('invalid number of columns\n')
end

numrows=D(1)
rimax=0
rjmax=0
rkmax=0
for I=1:numrows,
 if (finedata(I,1)>rimax) 
  rimax=finedata(I,1);
 end
 if (finedata(I,2)>rjmax) 
  rjmax=finedata(I,2);
 end
 if (finedata(I,3)>rkmax) 
  rkmax=finedata(I,3);
 end
end
imax=round(rimax)+1
jmax=round(rjmax)+1
kmax=round(rkmax)+1
x=1:1:imax;
y=1:1:jmax;
z=1:1:kmax;

volume1=zeros(imax,jmax,kmax);
 
for I=1:numrows,
 II=round(finedata(I,1))+1;
 JJ=round(finedata(I,2))+1;
 KK=round(finedata(I,3))+1;
 volume1(II,JJ,KK)=finedata(I,7);
end

figure,hold off
p = patch(isosurface(volume1, 0));
m1=get(p,'Vertices');
m2=get(p,'Faces');
isonormals(volume1, p)
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
daspect([1 1 1])
view(3)
camlight; lighting phong;
%axis([-6 imax+6 -6 jmax+6 -6 kmax+6])
axis off
view(142,22.0);
camlight;
alpha(0.6);
%M = getframe;hold off

