N=32
A=zeros(N,N);
D=zeros(N+2);
for i=1:N+1,
 D(i)=1.0;
 if ((i>=9)&(i<=17))
  D(i)=1.0;
 end
end
for i=1:N,
 A(i,i)=1.0/D(i)+1.0/D(i+1);
 if (i>1) 
  A(i,i-1)=1.0/D(i);
 end
 if (i<N)
  A(i,i+1)=1.0/D(i+1);
 end
end

acond=cond(A)
