function [ out ] = New_trial_Pd_Exact_PDvsPF( Pf )
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
mu =10.^(0/10);
c = 1.2;
A = gamma(1+2./c).^(c./2);
A
time  = [-1.2247 0 1.2247];
W = [0.2954 1.1816 0.2954];

Y = zeros(1, length(lambda));
for j = 1:length(lambda)
out  = zeros(1,1);
syms w                                                          %w=SNR for log normal 
syms y

p= A*(1/w)^(c/2);
p1=A*((1/w).*y).^(c/2);
p2=exp(-p1);
LN=exp(-((10*log10(w)-mu)^2)/(2*(std_dev)^2))*(1/(sqrt(2*pi)*std_dev*w));      % Log Normal Distribution
f=4.3429*(c/2)*p*p2*LN.*(y.^((c/2)-1));                        %Composite weibull log normal
l=int(f,w,0,50);                         %PDF for weibull over log normal

out = l*Marcum_Q(2*sqrt(y), sqrt(lambda), u);

out = int(out, y, 0, 10);
Y(j) = out;
end
end
function [out] = Marcum_Q(a, b, M) % M = u 
out = zeros(1,1);

for k = 1-M:10
    out  = ((a./b)^k)*(besseli(k, a*b)) + out;
end
out = exp(-(a^2 + b^2)./2)*out ;

end


