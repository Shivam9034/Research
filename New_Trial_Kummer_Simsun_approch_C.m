function [ y ] = New_Trial_Kummer_Simsun_approch_C(c)
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
std_dev = 10.^(5./10);  % ==1.20

mu = 10.^(8./10); % ==0.406
lambda = 15.98;
%Pf = 0.1;
u = 5;
X = [mu  mu+sqrt(3)*std_dev  mu-sqrt(3)*std_dev];
D = [0.66 0.1667 0.1667];

for j = 1: length(c)
%A = gamma(1+2./c(j)).^(c(j)./2);

out1 = zeros(1,1);
out2 = zeros(1,1);
out3 = zeros(1,1);
for h = 1:length(X)
%y = zeros(1,1);
% for alpha i
     
out1 =  D(h).*Marcum_Q(sqrt(2.*10), sqrt(lambda),u).*(c(j)./2).*((gamma(1+2./c(j))./exp(X(h))).^(c(j)./2)).*(10.^(c(j)./2 -1)).*exp((-1.*10.*gamma(1+2./c(j))./exp(X(h))).^(c(j)./2)) + out1; 
out1
out2 =  4*D(h).*Marcum_Q(sqrt(11), sqrt(lambda),u).*(c(j)./2).*((gamma(1+2./c(j))./exp(X(h))).^(c(j)./2)).*((11./2).^(c(j)./2 -1)).*exp((-1.*(11./2).*gamma(1+2./c(j))./exp(X(h))).^(c(j)./2)) + out2; 
out2
out3 =  D(h).*Marcum_Q(sqrt(2.*1), sqrt(lambda),u).*(c(j)./2).*((gamma(1+2./c(j))./exp(X(h))).^(c(j)./2)).*(1.^(c(j)./2 -1)).*exp((-1.*1.*gamma(1+2./c(j))./exp(X(h))).^(c(j)./2)) + out3; 
out3
end

y = 9./6*(out1+out2+out3);
Y(j) = y;
end

plot(c, Y)
end

