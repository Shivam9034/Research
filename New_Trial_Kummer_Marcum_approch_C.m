function [ y ] = New_Trial_Kummer_Marcum_approch_C(c)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%Pf = 0.1;

%lambda = ones(1, length(Pf));

  %  syms lamb;
   % a = check(u, lamb) == Pf*gamma(u);
   % lambda = solve(a, lamb);


%lamb = sym('lamb');

% == Pf*gamma(u);
%y = solve(a, lamb);
%r = a == P*gamma(U);
%lambda = solve(base, lamb);
 
%mu = 10.^(-3.914./10);
%c = 0.5;

%time  = [-1.2247 0 1.2247];
%W = [0.2954 1.1816 0.2954];


%alpha = zeros(1, length(time));
%for i = 1:length(time)    
 %    alpha(i) = exp(-1*c./2.*((sqrt(2)*std_dev.*(time(i))) + mu)); % mu == mean
%end
     
%Y = zeros(1, length(H));

Y = zeros(1,length(c));
syms gama
for j = 1: length(c)
fun =  summe(gama, c(j));
y = int(fun,gama , 0, 50);
Y(j) = y;
end

plot(c, Y)
end
function [out] = summe(x, d)
std_dev = 10.^(5./10);  % ==1.20

mu = 10.^(8./10); % ==0.406
lambda = 15.98;
%Pf = 0.1;
u = 5;
X = [mu  mu+sqrt(3)*std_dev  mu-sqrt(3)*std_dev];
di = [0.66 0.1667 0.1667];
A = gamma(1+2./d).^(d./2);
A
g = zeros(1,1);
for h = 1:length(X)
%y = zeros(1,1);
% for alpha i
     
g =  di(h).*Marcum_Q(sqrt(2.*x), sqrt(lambda),u).*(d./2).*((gamma(1+2./d)./exp(X(h))).^(d./2)).*(x.^(d./2 -1)).*exp((-1.*x.*gamma(1+2./d)./exp(X(h))).^(d./2)) + g; 

end
out = g;      
  
end

