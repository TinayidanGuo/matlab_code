% Nesterov's algorithm appling to exercise 2

clear;
clc;

%% Settings of parameter matrices and vectors
n=2;
BB = zeros(n,n);
S = zeros(n,n);
e_1=zeros(n,1);
e_n=zeros(n,1);
x=zeros(n,1);
%define matrix S
%S = zeros(n,n);
for i=1
	S(i,i)=1;
	S(i,i+1)=-1;
end
for i=2:(n-1)
	S(i,i-1)=-1;
	S(i,i)=2;
	S(i,i+1)=-1;
end
for i = n
	S(i,i)=1;
	S(i,i-1)=-1;
end
%S
%define matrix Bsqr(=BB)
%Bsqr = zeros(n,n);
for i=1
	BB(i,i)=2;
	BB(i,i+1)=-1;
end
for i=2:(n-1)
	BB(i,i-1)=-1;
	BB(i,i)=2;
	BB(i,i+1)=-1;
end
for i = n
	BB(i,i)=1;
	BB(i,i-1)=-1;
end
%BB
%find the eigenvalue of Bsqr
%eig(BB)
%V=right eigenvector of Bsqr, D=eigenvalue matrix
% [V,D]=eig(BB)
%  V'*S*V

%define vector e_n
%e_n=zeros(n,1);
for i=n
	e_n(i,1)=1;
end
%define vector e_1
%e_1=zeros(n,1);
for i=1
	e_1(i,1)=1;
end

%% Settings of initial values
h = .09; 
xup=[.5 1.5]';
x=xup; 
xstar=[1 2]';
maxIter=500;
lambda = 2/(n+1);
gamma = lambda;
beta = (2*L)/n;
L = 1/h;
M = L;
  format long
 
%% Nesterov's method algorithm in original function form
% define f function 
  f = @(x)((L/4)*x'*BB*x-(L/2)*e_n'*x);
  gradf =@(x)((L/2)*BB*x-(L/2)*e_n); 
    
% define g function
  g = @(x)((M/4)*x'*S*x-(M/2)*(e_n-e_1)'*x);
  gradg =@(x)((M/2)*S*x-(M/2)*(e_n-e_1)); 
  
% define phi function
  phi =@(x) f(x) + g(x);
  gradphi =@(x) gradf(x) + gradg(x);
  
%iteration for original function
  for iter =1:maxIter
      x=xup;
      xlw=(1-lambda)*xup + lambda*x;
      gradphilw = gradphi(xlw);
      x = x - (1/beta)*gradphilw;
      xup=(1-gamma)*xup + (gamma)*x;
  end  
  xup
 %% Define functions with U transform 
 
  %%With U being eigenvector decomposition of S,
 [U,D]=eig(S);
  A=U'*BB*U;
  D;%=U'*S*U
  
  % define fu function 
  fu = @(x)((L/4)*x'*A*x-(L/2)*(Ux)'*e_n);
  gradfu =@(x)((L/2)*A*x-(L/2)*U'*e_n); 
    
% define gu function
  gu = @(x)((M/4)*x'*D*x-(M/2)*(Ux)'*(e_n-e_1));
  gradgu =@(x)((M/2)*D*x-(M/2)*U'*(e_n-e_1)); 
  
% define phiu function
  phiu =@(x) fu(x) + gu(x);
  gradphiu =@(x) gradfu(x) + gradgu(x);
  
%   lambda = 2/(n+1);
%   gamma = lambda;
%   beta = (2*L)/n;
% iteration for function with U transformation  
  for iter =1:maxIter
      xlw=(1-lambda)*xup + lambda*x; %update 'x_lower bar'
      gradphiuxlw = gradphiu(xlw); %grad phiu evaluated at 'x_lowerbar'
      x = x - (1/beta)*gradphiuxlw; %update 'x'
      xup=(1-gamma)*xup + gamma*x; %update 'x_upper bar'
  end
  xup
  U\xup %to testify that inv(U)*xup is the solution of phi(with no U transf)
  