A=zeros(25,25);
rhs=zeros(25,1);
coeff=zeros(25,1);
lder=1
mder=1
idx1=1
for l=0:4,
for m=0:4,
 rhs(idx1)=0.0;
 if (l==lder) & (m==mder)
  rhs(idx1)=factorial(l)*factorial(m);
 end
 idx1=idx1+1;
end
end 
x=0.0;
y=0.0;
idx1=1;
for l=0:4,
for m=0:4,
 idx2=1;
 for i=-2:2,
 for j=-2:2,
  x=i; 
  y=j;
  h=power(x,l)*power(y,m);
  if (l<2)
   hxx=0.0;
  else
   hxx=l*(l-1)*power(x,l-2)*power(y,m);
  end
  if (m<2)
   hyy=0.0;
  else
   hyy=m*(m-1)*power(y,m-2)*power(x,l);
  end
  if (l<4)
   hxxxx=0.0;
  else
   hxxxx=l*(l-1)*(l-2)*(l-3)*power(x,l-4)*power(y,m);
  end
  if (m<4)
   hyyyy=0.0;
  else
   hyyyy=m*(m-1)*(m-2)*(m-3)*power(y,m-4)*power(x,l);
  end
  if ((l<2) | (m<2))
   hxxyy=0.0;
  else
   hxxyy=l*(l-1)*m*(m-1)*power(x,l-2)*power(y,m-2);
  end
  temp=h+hxx/24.0+hyy/24.0+hxxyy/(24.0*24.0);
  temp=temp+hxxxx/(16.0*120.0)+hyyyy/(16.0*120.0);
  A(idx1,idx2)=temp;
  idx2=idx2+1;
 end
 end
 idx1=idx1+1;
end
end
format long e

coeff=inv(A)*rhs
coeff1d=zeros(5,1);
idx2=1;
for i=-2:2,
 coeff1d(i+3)=0.0;
 for j=-2:2,
  coeff1d(i+3)=coeff1d(i+3)+coeff(idx2);
  idx2=idx2+1;
 end
end
coeff1d
