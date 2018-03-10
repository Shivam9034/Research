function [ y ] = New_trial_Pd_Equation( Pf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
std_dev = 10.^(1./10);
%Pf = 0.1;
u = 5;

lambda = zeros(1, length(Pf));
    for i = 1:length(Pf)
    lambda(i) = 2*(gammaincinv(1-Pf(i), u./2));
    end
lambda
mu =4;
c = 2;
A = gamma(1+2./c).^(c./2);
A
time  = [-1.2247 0 1.2247];
W = [0.2954 1.1816 0.2954];

alpha = zeros(1, length(time));
   for i = 1:3   %==3
    alpha(i) = exp(-1*c./2.*((sqrt(2)*std_dev.*(time(i))) + mu)); % mu == mean
   end
    alpha



Y = zeros(1, length(lambda));
for j = 1:length(lambda)
    
z = zeros(1,1);
syms gama;
for i = 1:length(alpha)
 
    z = alpha(i)*W(i)*((gama)^((c./2) -1))*exp(-1*((gama)^c./2)*A*alpha(i)) + z;
    
end

z = (c*A./(2*sqrt(pi)))*Marcum_Q(sqrt(2*gama), sqrt(lambda(j)), u)*z;

fun  =  z;

y = int(fun,gama,  0, 100) ;
Y(j) = y;

end
plot(Pf, Y)
end
function [out] = Marcum_Q(a, b, M) % M = u 
out = zeros(1,1);

for k = 1-M:10
    out  = ((a./b)^k)*(besseli(k, a*b)) + out;
end
out = exp(-(a^2 + b^2)./2)*out ;

end


