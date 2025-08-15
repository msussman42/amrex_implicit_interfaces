%drive.matlab.com -> upload -> import data -> all files 
%y=1.1 e^(ct)
%log((y-1.0)/0.1)=c t
data=load("ampfine");
time=data(1:end,1);
y=data(1:end,2);
ymod=log((y-1.0)/0.1);
plot(time,ymod);
