% AGS algorithm appling to exercise 2

clear;
clc;

%% Settings of parameter matrices and vectors
n=200; %works when n is < 20, solution loses accuracy when n=30 and up
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
%h=1e-5;% decrease h pushes up accuracy for n upto 30
xup=zeros(n,1); % set initial value of 'x_upper bar'
x=xup;
xstar=[1:n]';

maxIter=5000; %maxIter=5000 doesn't work 
%L = 1/(2*h);
%L = 1e10; in debug, use a large L value to expell other mistakes
L = 100;
M = L*L;

  format long
 
%% Define functions

% define f function 
  f = @(x)((L/8)*x'*BB*x-(L/4)*e_n'*x);
  gradf =@(x)((L/4)*BB*x-(L/4)*e_n); 
    
% define g function
  g = @(x)((M/8)*x'*S*x-(M/4)*(e_n-e_1)'*x);
  gradg =@(x)((M/4)*S*x-(M/4)*(e_n-e_1)); 
  
% define phi function
  phi =@(x) f(x) + g(x);
  gradphi =@(x) gradf(x) + gradg(x);

%% AGS method algorithm iterations
  %fx = f(xup);
  vec_fx = zeros(maxIter,1);
  vec_ec = zeros(maxIter,1);
  %T = ceil(sqrt(M/L)); %set T to be square root of the ratio between M n L
  T=1;
  tic;
  for iter =1:maxIter
      gamma = 2/(iter+1);
      lambda = gamma;
      beta = 2*L/iter;
      xlw=(1-gamma)*xup + gamma*x;
      gradfxlw = gradf(xlw); %gradf evaluated at 'x_low'
      d = gradfxlw;
      utld = xup;%set initial value of 'u_tilde' = 'x_upperbar'
      u = x; %set initial value of 'u' = 'x'      
      for t = 1:T
          if t==1
              alpha = 1; %set initial value of alpha
              Lambda = 1; %set initial value of Lambda
          else
              alpha =sqrt(4*Lambda+Lambda^2)/2 -Lambda/2; %update alpha
              Lambda = Lambda*(1-alpha); %update Lambda
          end
          p = 1/alpha -1;% find p at iteration t
          q = (2*M*alpha)/iter; %find q at iteration t
          ulw=(1-gamma)*xup + gamma*(1-alpha)*utld + gamma*alpha*u; % the expression of 'u_lowerbar'
          ud = beta*(1+p)+q;% define the denorminator of the expression of u to be 'ud' 
          un = beta*x + (beta*p+q)*u - gradg(ulw)-d; % define the numerator of the expression of u to be 'un' 
          u = un/ud; %find the expression of u
          utld=(1-alpha)*utld + alpha*u; %find 'u_tilde' 
          up=(1-gamma)*xup + gamma*utld; %update 'u_upperbar'
      end
      x = u; % update x as the value of u at iteration T 
      xup = up; % update xup as the value of up at iteration T
      fx = phi(xup);
      vec_fx(iter) = fx;
%       %vec_ec(iter)= cumsum(toc);
      vec_ec(iter)= toc;
      max(xup);
      phi(xup)- phi(xstar);% decreases until n=46, then increase
  end 
  xup';
  max(xup)
  phi(xup)- phi(xstar)
  toc;
  
  vec_fx = vec_fx(1:iter);
  vec_ec = vec_ec(1:iter);
  tiledlayout(1,2);
  nexttile
  plot(vec_ec,vec_fx)
  nexttile
  plot(vec_fx)
  
  %% print mutiple figures in one window
%   tiledlayout(1,3);
%   nexttile
%   plot(vec_ec,vec_fx)
%   nexttile
%   plot(vec_fx)
%   nexttile
%   plot(vec_ec)

 %% Remark
% When debugging, following strategy might helpful
% figure
% plot(vec_fx(1:100)) % plot the first 100 "vec_fx" values
% set algorithm Lipschitz constant value L = function Lip value (L/4*BB)
% find(vec_fx > 0,1) % find the iteration that vec_fx >0 at first time
% if algorithm converges firstly but then quickly diverges and blow up,
%   possbilly step size is too large.

