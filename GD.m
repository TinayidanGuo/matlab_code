% gradient Descent method appling to exercise 2

clear;
clc;

%% Settings of initial values
n=100; %works when n is up to 5
x=zeros(n,1);
 
%xstar=[1:n]';
maxIter=50000;
L = 10;
M=L;
h = 1/(2*L);%h>1/L

%% Settings of parameter matrices and vectors
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
BB(1,1)= 2;
BB(1,2)= -1;
BB(n,n)= 1;
BB(n,n-1)= -1;
for i=2:(n-1)
	BB(i,i-1)=-1;
	BB(i,i)=2;
	BB(i,i+1)=-1;
end

%define vector e_n
e_n(n,1)=1;

%define vector e_1
e_1(1,1)=1;

%% GD with original function form
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
  vec_fx = zeros(maxIter,1);
  vec_ec = zeros(maxIter,1);
%gradient decent method for exercise 2
tic
  for iter =1:maxIter   
      df = gradf(x);
      dg = gradg(x);
      d = df+dg;
      %gradphix = gradphi(x);
      %d = gradphix;
      %iter;
      x = x - h*d;
      fx = phi(x);
      vec_fx(iter) = fx;
      vec_ec(iter)= toc;
      max(x);
  end
x';
  vec_fx = vec_fx(1:iter);
  vec_ec =vec_ec(1:iter);
  plot(vec_ec,vec_fx)
  plot(vec_fx)
