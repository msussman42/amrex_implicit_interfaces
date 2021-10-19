z=zeros(33,33,2);
for i=1:33
for j=1:33
 x=(i-1)/32.0;
 y=(j-1)/32.0;
 z(i,j,1)=(x-0.5)^2+(y-0.5)^2;
 z(i,j,2)=x+y;
end
end
contour(z(:,:,1))
hold on
contour(z(:,:,2))
quiver(z(:,:,1),z(:,:,2))
%contour (uniform* files), quiver (uniform* files), isosurface (uniform* files), scatter (PARCON*)
