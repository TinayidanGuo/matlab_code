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
maxIter=5000; 
L_prev = 1;
L=5; 
M_prev = 5;
M=100;
gamma=1;
T=1;
%with this setting, maxIter=500, found L=128, M=160
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
  

  
%% Backtracking linesearch for L & M
vec_L = zeros(maxIter,1);
vec_ratio = zeros(maxIter,1);

vec_M = zeros(T,1);
vec_ratio1 = zeros(T,1);

tic
for   iter =1:maxIter      
      xup_prev= xup;
      x_prev =x;
      gamma_prev = gamma;
      L = .5*L_prev; 
      while true 
          %T = ceil(sqrt(M/L));
          beta = 2*L/iter;
          
          r = L_prev/L;
          gamma = (-r*gamma_prev^2+sqrt(r^2*gamma_prev^4 + 4*r*gamma_prev^2))/2;   
                    
          xlw=(1-gamma)*xup_prev + gamma*x_prev;
          df = gradf(xlw);
          dphi = gradphi(xlw);
          for t = 1:T
              utld = xup_prev;
              u = x_prev;  
              M = .5*M_prev;  
              while true
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
                  
                  dg =gradg(ulw);
                  ud = beta*(1+p)+q;% define the denorminator of the expression of u to be 'ud' 
                  un = beta*x + (beta*p+q)*u - gradg(ulw)-df; % define the numerator of the expression of u to be 'un' 
                  u = un/ud; %find the expression of u
                  utld=(1-alpha)*utld + alpha*u; %find 'u_tilde' 
                  up=(1-gamma)*xup + gamma*utld; %update 'u_upperbar'   

                  ls = g(up)- g(ulw)- (up-ulw)'*dg;
                  rs = norm(up-ulw)^2;         
                  ratio1 = 2*ls/rs;  
                  vec_M(t) = M;
                  vec_ratio1(t)=ratio1;        
                  if ls <= M*rs/2
                      break;
                  else
                      M = 2*M; %M_k <=2*M
                  end   
              end
              M_prev = M;
              T = ceil(sqrt(M/L));
          end
          x = u; % update x as the value of u at iteration T 
          xup = up; % update xup as the value of up at iteration T         
          lhs = phi(xup)- phi(xlw)- (xup-xlw)'*dphi;
          rhs = norm(xup-xlw)^2;         
          ratio = 2*lhs/rhs;  
          vec_L(iter) = L;
          vec_ratio(iter)=ratio;        
          if lhs <= L*rhs/2
              break;
          else
              L = 2*L; %L_k <=2*L
          end         
      end  
      L_prev = L;
end  
phi(xup)- phi(xstar)
%Now xup=xlw is the optimal solution
%L =10 is found out
  %% plot L value and objective value
  figure
  tiledlayout(2,2);
  nexttile
  vec_L =vec_L(1:iter);
  plot(vec_L) 
  title('L vs iteration')

  nexttile
  vec_M =vec_M(1:T);
  plot(vec_M) 
  title('M vs iteration')
  
  nexttile
  vec_ratio =vec_ratio(1:iter);
  plot(vec_ratio) 
  title('ratio(L) vs iteration')
  
  nexttile
  vec_ratio1 =vec_ratio1(1:T);
  plot(vec_ratio1) 
  title('ratio1(M) vs iteration')
 

%% Debug
% plot(vec_L(1:200))
% find(vec_L > 100,1)
% find(vec_ratio < vec_L,1)  
 
