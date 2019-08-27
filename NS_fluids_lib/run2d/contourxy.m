finedata=load('visual00010');
D=size(finedata)
if (D(2)~=7) 
 sprintf('invalid number of columns\n')
end

numrows=D(1)
imax=24
jmax=64

x=zeros(imax,jmax);
y=zeros(imax,jmax);

volume1=zeros(imax,jmax);

II=1
JJ=1 
for I=1:numrows,
 x(II,JJ)=finedata(I,1);
 y(II,JJ)=finedata(I,2);
 volume1(II,JJ)=finedata(I,5);
 JJ=JJ+1;
 if (JJ>jmax),
  JJ=1;
  II=II+1;
 end
end
hold off;
%%[C,LEV]=contourc(x,y,volume1,[0,0]);
contour(x,y,volume1,[0,0]);
hold on;
