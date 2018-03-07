function [ y ] = New_Trial_Kummer_R_approch_mu(Mu)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Pf = 0.1
std_dev = 10.^(0.8./10);

mu = 10.^(Mu./10)
%Pf = 0.1;
u = 5;
%lambda = ones(1, length(Pf));

  %  syms lamb;
   % a = check(u, lamb) == Pf*gamma(u);
   % lambda = solve(a, lamb);
%lambda(1) = 15.98;

lamb = sym('lamb');

% == Pf*gamma(u);
%y = solve(a, lamb);
%r = a == P*gamma(U);
lambda = solve(base, lamb);
 
%mu = 10.^(-3.914./10);
c = 0.5;
A = gamma(1+2./c).^(c./2);

time  = [-1.2247 0 1.2247];
W = [0.2954 1.1816 0.2954];
di = [0.66 0.1667 0.1667];

%alpha = zeros(1, length(time));
%for i = 1:length(time)    
 %    alpha(i) = exp(-1*c./2.*((sqrt(2)*std_dev.*(time(i))) + mu)); % mu == mean
%end
     
%Y = zeros(1, length(H));

Y = zeros(1,length(mu));
     
for j = 1: length(mu)
X = [mu(j)  mu(j)+sqrt(3)*std_dev  mu(j)-sqrt(3)*std_dev];

g = zeros(1,1);
for h = 1:length(X)
%y = zeros(1,1);
% for alpha i
     

              out1 = zeros(1,1);
              for L = 0:u-1
              out1 = (((lambda(1))^L)*(exp(-1*(lambda(1)./2))))./((gamma(L+1))*(2^L)) + out1;
              end
   out1 = out1*di(h);
              out2  = zeros(1,1);
              for k = 0:101
              out2  = ((lambda(1).^u)*((-1).^k)*gamma(1+c*k/2)*((A.^k)*kummer(1+c*k./2 , u+1, lambda(1)./2)))./((2.^(u))*(exp(X(h)).^(c*k/2))*gamma(k+1)*(exp(lambda(1)./2))*gamma(u+1)) + out2 ;
              end
out2 = out2*di(h);
g = out1 + out2 +g
end
y = g;
Y(j) = y;

end
plot(Mu, Y)
%function [out] = sum1()
%out = 1;
%for L = 0:u-1
%   out = (((lambda)^L)*(exp(-1*(lambda./2))))./((gamma(L+1))*(2^L)) + out;
%end
%end

%function [out] = sum2(alpha)
    
    %   out  = 1;
    %for k = 0:inf
     %  out  = ((lambda^u)*((-1)^k)*2*gamma(1+c*k/2)*((A.*alpha)^k)*kummer(1+c*k/2, u+1, lambda/2))/((2^u)*(2^(c*k/2))*gamma(k+1)*(exp(lambda/2))*gamma(u+1)) + out ;
    %end
%end

function [f] = kummer(a,b,x)

% This function estimates the Kummer function with the specified tolerance
% the generalized hypergeometric series, noted below.  This solves Kummer's
% differential equation:
%
%       x*g''(x) + (b - x)*g'(x) - a*g(x) = 0

% Default tolerance is tol = 1e-10.  Feel free to change this as needed.
tol = 1e-10;

% Estimates the value by summing powers of the generalized hypergeometric
% series:
%
%       sum(n=0-->Inf)[(a)_n*x^n/{(b)_n*n!}
%
% until the specified tolerance is acheived.

term = x*a/b;
f = 1 + term;
n = 1;
an = a;
bn = b;
nmin = 10;
while(n < nmin)||max(abs(term) > tol)
  n = n + 1;
  an = an + 1;
  bn = bn + 1;
  term = x.*term*an/bn/n;
  f = f + term;
end

% VERSION INFORMATION
% v1 - Written to support only scalar inputs for x
% v2 - Changed to support column inputs for x by using the repmat
% command and using matrix multiplication to achieve the desired sum
%
% v3 - Credit goes to Ben Petschel for making this suggestion.
%    The previous method of creating vectors for multiplication to
%    produce the sum was replaced by a while loop that executes
%    until a certain tolerance is achieved.  My previous thinking
%    was avoiding a loop would produce a code that would execute
%    faster.  Ben pointed out this is not necessarily true, and not
%    true in this case.  Not only does the while loop used execute
%    faster for this calculation, but it is also more accurate.
end


end
