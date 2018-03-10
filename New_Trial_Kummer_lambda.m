function [ y ] = New_Trial_Kummer_lambda(Pf)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


std_dev = 10.^(0.806./10);
%Pf = 0.1;
u = 5;
lambda = ones(1, length(Pf));
for i = 1:length(Pf)
lambda(i) = 2*(gammaincinv(1-Pf(i), u./2));
end
lambda
mu = 10.^(-3.914./10);
c = 1.5;
A = gamma(1+2./c).^(c./2);

time  = [-1.2247 0 1.2247];
W = [0.2954 1.1816 0.2954];

alpha = zeros(1, length(time));
for i = 1:length(time)    
     alpha(i) = exp(-1*c./2.*((sqrt(2)*std_dev.*(time(i))) + mu)); % mu == mean
end
     
Y = zeros(1, length(mu));
for j = 1: length(Pf)


% for alpha i
     

     y = zeros(1,1);
     for i= 1:length(time)

              out1 = zeros(1,1);
              for L = 0:u-1
              out1 = (((lambda(j))^L)*(exp(-1*(lambda(j)./2))))./((gamma(L+1))*(2^L)) + out1;
              end
   
              out2  = zeros(1,1);
              for k = 0:10
              out2  = ((lambda(j).^u)*((-1).^k)*gamma(1+c*k/2)*((A*alpha(i)).^k)*kummer((1+c*k./2, u+1, lambda(j)/2))/((2.^(u/2))*(2.^(c*k/2))*gamma(k+1)*(exp(lambda(j)/2))*gamma(u+1)) + out2 ;
              end

      y = W(i)*(out1 + out2) + y;
     end

      y = (1./sqrt(pi)).*y;

   Y(j) = y;

end
plot(Pf, Y)
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
