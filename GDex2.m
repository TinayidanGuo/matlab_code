% gradient Descent method appling to exercise 2

clear;
clc;

%% Settings of initial values
x=[1.9;2.3;3.5;4.4;5.9];
h = .001; 
xstar=[1:5]';
maxIter=50000;
L = 1/(1.9*h);
M=L;

n=5;

%% Settings of parameter matrices and vectors

%define matrix S
S = zeros(n,n);
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

%define matrix BB
BB = zeros(n,n);
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
%find the eigenvalue of BB
%eig(BB)
%V=right eigenvector of BB, D=eigenvalue matrix
% [V,D]=eig(BB)
%  V'*S*V
% [V,D]=eig(S);
%  V'*BB*V;

%define vector e_n
e_n=zeros(n,1);
for i=n
	e_n(i,1)=1;
end

%define vector e_1
e_1=zeros(n,1);
for i=1
	e_1(i,1)=1;
end

%% GD with No orthogonal matrix transformation
format long
% define f function 
  f = @(x)((L/4)*x'*BB*x-(L/2)*e_n'*x);
  gradf =@(x)((L/2)*BB*x-(L/2)*e_n); 
   
% define g function
  e = e_n-e_1;
  g = @(x)((M/4)*x'*S*x-(M/2)*e'*x);
  gradg =@(x)((M/2)*S*x-(M/2)*e); 
  
% define phi function
  phi =@(x) (f(x) + g(x));
  gradphi =@(x) (gradf(x) + gradg(x));
  %gradphi =@(x)((L/2)*BB*x+(M/2)*S*x-(M/2)*e-(L/2)*e_n);

%gradient decent method for exercise 2
  for iter =1:maxIter   
      df = gradf(x);
      dg = gradg(x);
      d = df+dg;
      %gradphix = gradphi(x);
      %d = gradphix;
%       iter;
      x = x - h*d;
  end
x 
%% GD with U being eigenvector matrix of S,
  [U,D]=eig(S);
  A=U'*BB*U;
  D;%=U'*S*U
  
  % define fu function 
  fu = @(x)((L/4)*x'*A*x-(L/2)*(Ux)'*e_n); 
  %optimal value = U'*x, where x is optimal value of f
  gradfu =@(x)((L/2)*A*x-(L/2)*U'*e_n); 
    
  % define gu function
  gu = @(x)((M/4)*x'*D*x-(M/2)*(Ux)'*(e_n-e_1));
  gradgu =@(x)((M/2)*D*x-(M/2)*U'*(e_n-e_1)); 
  
  % define phiu function
  phiu =@(x) fu(x) + gu(x);
  gradphiu =@(x) gradfu(x) + gradgu(x);

  for iter =1:maxIter 
      dfu = gradfu(x);
      dgu = gradgu(x);
      du = dfu+dgu;
      x = x - h*du;      
  end
 x
 U\x
%% theoretical solution of phi function:
% v = L*BB+M*S;
% b = M*e+L*e_n;
% v\b

% v = L*A+M*D;
% b = M*U'*e+L*U'*e_n;
% v\b