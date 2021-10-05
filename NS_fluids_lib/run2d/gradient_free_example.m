clear;
clc;
% these two should actually be equivalent.  "0" is the seed.
rng('default')
%rng(0,'twister')

% the dimension of x is implicitly defined by how "costFunc" is defined.
J = @(x) costFunc(x);
%https://www.mathworks.com/help/gads/patternsearch.html#buxdit7-options
%options = optimoptions('patternsearch','PlotFcn','psplotbestf', 'MaxIter', 10000, 'MaxFunEvals',10000, 'pollmethod', 'MADSPositiveBasis2N');
options = optimoptions('patternsearch', 'MaxIter', 10000, 'MaxFunEvals',10000, 'pollmethod', 'MADSPositiveBasis2N');

global N_DATA X_DATA Y_DATA CRITERION

%CRITERION=0  % Chebyshev
CRITERION=2  % Least Squares

N_DATA=5
N_MODEL=3
X_DATA=zeros([1 N_DATA])  %row vector
Y_DATA=zeros([1 N_DATA]) 

X_DATA(1)=0.1
X_DATA(2)=0.2
X_DATA(3)=0.3
X_DATA(4)=0.4
X_DATA(5)=0.5
Y_DATA(1)=0.06
Y_DATA(2)=0.12
Y_DATA(3)=0.36
Y_DATA(4)=0.65
Y_DATA(5)=0.95

u_IC=zeros([1 N_MODEL]) % row vector
for i=1:N_MODEL
 u_IC(i)=0.0
end

[u_optimal,fval,exitflag,output] = patternsearch(J,u_IC,[],[],[],[],[],[],[],options) 

N_RANGE=100;
%hold off;
plot(X_DATA,Y_DATA,'.b','MarkerSize',15);
hold on;
X_RANGE=linspace(X_DATA(1),X_DATA(N_DATA),N_RANGE);
Y_RANGE=zeros(1,length(X_RANGE));
for i=1:length(X_RANGE)
 X_MODEL=X_RANGE(i);
 Y_MODEL=u_optimal(1)*X_MODEL^2+u_optimal(2)*X_MODEL+u_optimal(3);
 Y_RANGE(i)=Y_MODEL;
end
plot(X_RANGE,Y_RANGE,'--','LineWidth',2);

function Jreturn = costFunc(u)

 global N_DATA N_MODEL X_DATA Y_DATA CRITERION

 local_cost=0.0;

 for i=1:N_DATA
  X_MODEL=X_DATA(i);
  Y_MODEL=u(1)*X_MODEL^2+u(2)*X_MODEL+u(3);
  Y_OBS=Y_DATA(i);
  DEVIATION=abs(Y_MODEL-Y_OBS);
  if (CRITERION==0) % Chebyshev criterion
   if (DEVIATION>local_cost)
    local_cost=DEVIATION;
   end
  elseif (CRITERION==2) % least squares
   local_cost=local_cost+DEVIATION^2;
  else
   disp('CRITERION invalid')
   quit
  end
 end %for i=1:N_DATA

 Jreturn=local_cost;

end % function costFunc

