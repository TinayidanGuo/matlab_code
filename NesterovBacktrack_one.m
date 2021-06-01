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
x=zeros(n,1); 
xstar=[1:n]';
maxIter=500;
L_prev = 5;
L=100;
M=40; %L starting at 100, find L=M=40
%M=100; %L starting at 10, find L=0
%M=.1; %L starting at 10, find L=5
gamma=1;


  format long
 
%% Define functions
% define f function 
  f = @(x)((L/8)*x'*BB*x-(L/4)*e_n'*x);
  gradf =@(x)((L/4)*BB*x-(L/4)*e_n); 
    
% define g function
  g = @(x)((M/8)*x'*S*x-(M/4)*(e_n-e_1)'*x);
  gradg =@(x)((M/4)*S*x-(M/4)*(e_n-e_1)); 
  
% define phi function
  phi =@(x) (f(x) + g(x));
  gradphi =@(x) (gradf(x) + gradg(x));
  

  
%% Backtracking linesearch for L
% Define lhs = phi(xup)- phi(xlw)- (xup-xlw)'*gradphi(xlw);
% and rhs = norm(xup-xlw)^2;
% then L >= 2*lhs/rhs, true;
% vec_x = zeros(n,maxIter);
% vec_xup = zeros(n,maxIter);
% vec_gamma = zeros(maxIter);
% vec_fx = zeros(maxIter,1);
% vec_ec = zeros(maxIter,1);

vec_L = zeros(maxIter,1);
vec_ratio = zeros(maxIter,1);

tic
for   iter =1:maxIter    
      xup_prev= xup;
      x_prev =x;
      gamma_prev = gamma;
      L = .5*L_prev; 
      while true         
          beta = 2*(L+M)/iter;
          r = (L_prev+M)/(L+M);
          gamma = (-r*gamma_prev^2+sqrt(r^2*gamma_prev^4 + 4*r*gamma_prev^2))/2;        
          xlw=(1-gamma)*xup_prev + gamma*x_prev;
          d = gradphi(xlw);
          x = x_prev - (1/beta)*d;
          xup=(1-gamma)*xup_prev + gamma*x;
          lhs = phi(xup)- phi(xlw)- (xup-xlw)'*d;
          rhs = norm(xup-xlw)^2;
          ratio = 2*lhs/rhs;   
          vec_L(iter) = L;
          vec_ratio(iter)=ratio;       
          if lhs <= (L+M)*rhs/2
              break;
          else
              L = 2*L; %L_k <=2*L
          end         
      end  
      L_prev = L;
%      vec_ec(iter)= toc;
end  
%Now xup=xlw is the optimal solution
%L =10 is found out
  %% plot L value and objective value
  figure
  tiledlayout(2,1);
  nexttile
  vec_L =vec_L(1:iter);
  %legend({'y = sin(x)','y = cos(x)'},'Location','southwest')
  plot(vec_L) 
  title('L vs iteration')


  nexttile
  vec_ratio =vec_ratio(1:iter);
  %legend({'y = sin(x)','y = cos(x)'},'Location','southwest')
  plot(vec_ratio) 
  title('ratio vs iteration')

%%
%   nexttile
%   vec_fx = vec_fx(1:iter);
%   plot(vec_fx);
%   title('objective value vs iteration')
%   
%   nexttile
%   vec_ec =vec_ec(1:iter);
%   plot(vec_ec,vec_fx)
%   title('objective value vs eclapse time')
%% Debug
% plot(vec_L(1:200))
% find(vec_L > 100,1)
% find(vec_ratio < vec_L,1)  
 