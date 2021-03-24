% AGS algorithm appling to exercise 2

clear;
clc;

%% Settings of parameter matrices and vectors
n=2;
BB = zeros(n,n);
S = zeros(n,n);
e_1=zeros(n,1);
e_n=zeros(n,1);
%x=zeros(n,1);

%define matrix S
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

%define vector e_n
for i=n
	e_n(i,1)=1;
end

%define vector e_1
for i=1
	e_1(i,1)=1;
end

%% Let's solve for alpha first
% t=50;
% alpha_0 = 1;
% Lambda = alpha_0;
% for j=1:t
%   j;
%   alpha =sqrt(5*Lambda)/2 -Lambda/2;
%   Lambda = Lambda*(1-alpha);
% end 

%% Settings of initial values
h = .09; 
xup=[.5 1.5]'; % set initial value of 'x_upper bar'
x=xup;
xstar=[1 2]';
maxIter=500;
L = 1/h;
M = L*L;
gamma = 2/(n+1);
lambda = gamma;
beta = 2*L/n;

  format long
 
%% AGS method algorithm in original function form
% define f function 
  f = @(x)((L/4)*x'*BB*x-(L/2)*e_n'*x);
  gradf =@(x)((L/2)*BB*x-(L/2)*e_n); 
    
% define g function
  g = @(x)((M/4)*x'*S*x-(M/2)*(e_n-e_1)'*x);
  gradg =@(x)((M/2)*S*x-(M/2)*(e_n-e_1)); 
  
% define phi function
  phi =@(x) f(x) + g(x);
  gradphi =@(x) gradf(x) + gradg(x);
  
% AGS's method for exercise 2

%%un-comment the following code to perform AGS algorithm with no
%%transformation

  utld = xup;%set initial value of 'u_tilde' = 'x_upperbar'
  u = x; %set initial value of 'u' = 'x'
  T = sqrt(M/L); %set T to be square root of the ratio between M n L
  
  for iter =1:maxIter
      xlw=(1-gamma)*xup + gamma*x;
      gradfxlw = gradf(xlw); %gradf evaluated at 'x_low'
      d = gradfxlw;
      for t = 1:T 
          if t==1
              alpha = 1; %set initial value of alpha
              Lambda = alpha; %set initial value of Lambda
          else
          for j=1:t-1 % find alpha at iteration t
              j;
              alpha =sqrt(5*Lambda)/2 -Lambda/2; %update alpha
              Lambda = Lambda*(1-alpha); %update Lambda
          end
          end
          p = 1/alpha -1;% find p at iteration t
          q = (2*M*alpha)/n; %find q at iteration t
          ulw=(1-gamma)*xup + gamma*(1-alpha)*utld + gamma*alpha*u; % the expression of 'u_lowerbar'
          ud = beta*(1+p)+q;% define the denorminator of the expression of u to be 'ud' 
          un = beta*x + (beta*p+q)*u - gradg(ulw)-d; % define the numerator of the expression of u to be 'un' 
          u = un/ud; %find the expression of u
          utld=(1-alpha)*utld + alpha*u; %find 'u_tilde' 
          up=(1-gamma)*xup + gamma*utld; %update 'u_upperbar'
      end
      x = u; % update x as the value of u at iteration T 
      xup = up; % update xup as the value of up at iteration T 
  end  
  xup
  
 %% Define functions with U transform 
 
  %%With U being eigenvector matrix of S,
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
    % phiu =@(x) fu(x) + gu(x);
    % gradphiu =@(x) gradfu(x) + gradgu(x);
  
  % AGS's method for functions with U transformation
  
    utld = xup;%set initial value of 'u_tilde' = 'x_upperbar'
    u = x; %set initial value of 'u' = 'x'
    
  for iter =1:maxIter
      xlw=(1-gamma)*xup + gamma*x;
      gradfuxlw = gradfu(xlw); %gradf evaluated at 'x_low'
      d = gradfuxlw;
      T = sqrt(M/L); %set T to be square root of the ratio between M n L
      for t = 1:T 
          if t==1
              alpha = 1; %set initial value of alpha
              Lambda = alpha; %set initial value of Lambda
          else
          for j=1:t-1 % find alpha at iteration t
              j;
              alpha =sqrt(5*Lambda)/2 -Lambda/2; %update alpha
              Lambda = Lambda*(1-alpha); %update Lambda
          end
          end
          p = 1/alpha -1;% find p at iteration t
          q = (2*M*alpha)/n; %find q at iteration t
          ulw=(1-gamma)*xup + gamma*(1-alpha)*utld + gamma*alpha*u; % the expression of 'u_lowerbar'
          ud = beta*(1+p)+q;% define the denorminator of the expression of u to be 'ud' 
          un = beta*x + (beta*p+q)*u - gradgu(ulw)-d; % define the numerator of the expression of u to be 'un' 
          u = un/ud; %find the expression of u
          utld=(1-alpha)*utld + alpha*u; %find 'u_tilde' 
          up=(1-gamma)*xup + gamma*utld; %update 'u_upperbar'
      end
      x = u; % update x as the value of u at iteration T 
      xup = up; % update xup as the value of up at iteration T 
  end  
  xup
  U\xup %to testify that inv(U)*xup is the solution of phi(with no U transf)
  