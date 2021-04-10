% Nesterov's algorithm appling to exercise 2

clear;
clc;

%% Settings of parameter matrices and vectors
n=10;
BB = zeros(n,n);
S = zeros(n,n);
e_1=zeros(n,1);
e_n=zeros(n,1);
%x=zeros(n,1);

%define matrix S
S(1,1)=1;
S(1,2)=-1;
S(n,n)=1;
S(n,n-1)=-1;
for i=2:(n-1)
	S(i,i-1)=-1;
	S(i,i)=2;
	S(i,i+1)=-1;
end

%define matrix BB
BB(1,1)=2;
BB(1,2)=-1;
BB(n,n)=1;
BB(n,n-1)=-1;
for i=2:(n-1)
	BB(i,i-1)=-1;
	BB(i,i)=2;
	BB(i,i+1)=-1;
end

%define vector e_n
e_n(n,1)=1;

%define vector e_1
e_1(1,1)=1;

%% Settings of initial values

xup=zeros(n,1);
x=xup; 
xstar=[1:n]';
maxIter=500;
L = 100;
M=L;
%h = 1/L;%h>1/L

  format long
 
%% Nesterov's method algorithm in original function form
% define f function 
  f = @(x)((L/80)*x'*BB*x-(L/40)*e_n'*x);
  gradf =@(x)((L/40)*BB*x-(L/40)*e_n); 
    
% define g function
  g = @(x)((M/80)*x'*S*x-(M/40)*(e_n-e_1)'*x);
  gradg =@(x)((M/40)*S*x-(M/40)*(e_n-e_1)); 
  
% define phi function
  phi =@(x) (f(x) + g(x));
  gradphi =@(x) (gradf(x) + gradg(x));
  
  vec_fx = zeros(maxIter,1);
  vec_ec = zeros(maxIter,1);
  
%iteration for original function
tic
  for iter =1:maxIter
      lambda = 2/(iter+1);
      gamma = lambda;
      beta = (2*L)/iter; %Works when (20*L)/iter
      xlw=(1-lambda)*xup + lambda*x;
      d = gradphi(xlw);
      x = x - (1/beta)*d;
      xup=(1-gamma)*xup + gamma*x;
      %save obj value and eclapse time of each iter
      fx = phi(xup);
      vec_fx(iter) = fx;
      vec_ec(iter)= toc;
      iter;
      max(xup); % increases at an increase rate 
      f(xstar)-f(xup); %diverges at iter=4   
  end  
  xup;
  max(xup)
  f(xup)-f(xstar);
  vec_fx = vec_fx(1:iter);
  plot(vec_fx);
  vec_ec =vec_ec(1:iter);
  plot(vec_ec,vec_fx)

%% Remark:
% Setting function Lip half doesn't work without beta *10, 
% solution is to shrink funtion Lip by 10 times to trade off beta change.
% Now it even works with n=1000, iter=50000
 