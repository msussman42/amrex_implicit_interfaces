finedata=load('visual00010');
D=size(finedata)
if (D(2)~=9) 
 sprintf('invalid number of columns\n')
end

numrows=D(1)
imax=24
jmax=32
kmax=64

x=zeros(imax,jmax,kmax);
y=zeros(imax,jmax,kmax);
z=zeros(imax,jmax,kmax);

volume1=zeros(imax,jmax,kmax);
% U=zeros(imax,jmax,kmax);
% V=zeros(imax,jmax,kmax);
% W=zeros(imax,jmax,kmax);

II=1
JJ=1
KK=1 
for I=1:numrows,
 x(II,JJ,KK)=finedata(I,1);
 y(II,JJ,KK)=finedata(I,2);
 z(II,JJ,KK)=finedata(I,3);
%  U(II,JJ,KK)=finedata(I,4);
%  V(II,JJ,KK)=finedata(I,5);
%  W(II,JJ,KK)=finedata(I,6);
 volume1(II,JJ,KK)=finedata(I,7);
 KK=KK+1;
 if (KK>kmax),
  KK=1;
  JJ=JJ+1;
  if (JJ>jmax),
   JJ=1;
   II=II+1;
  end
 end   
 end

[F V]=isosurface(z, x, y, volume1, 0.0);
p = patch('Faces',F,'Vertices',V);
%isonormals(x,y,z,volume1, p)
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
daspect([1 1 1])
view(3)
%%camlight; 
lighting phong;
%axis([x(1,1,1) x(1,1,imax) y(1,1,1) y(1,1,jmax) z(1,1,1) z(1,1,kmax)])
%axis off
%view(-37.5,10.0);
camlight
