function [ y ] = New_Trial_Kummer_Marcum_approch_C(c)

Y = zeros(1,length(c));

for j = 1: length(c)
fun = @(gama) summe(gama, c(j));
g = integral(fun , 0, 50);
y = g;
Y(j) = y;
end

plot(c, Y)
end
function [out] = summe(x, d)
std_dev = 10.^(5./10);  % ==3.1623

mu = 10.^(8./10); % ==6.3096
lambda = 15.98;
%Pf = 0.1;
u = 5;
X = [mu  mu+sqrt(3)*std_dev  mu-sqrt(3)*std_dev];
D = [0.66 0.1667 0.1667];
A = gamma((1+2./d).^(d./2)); %==approx 0.9

g = zeros(1,1);
for h = 1:length(X)
%y = zeros(1,1);
% for alpha i
     
x1 = D(h).*Marcum_Q(sqrt(2.*x), sqrt(lambda),u).*d.*A.*(x.^(d./2 -1))  ;
x2  = 2.*(exp(X(h)*d./2)).*exp(((x./exp(X(h))).^(d./2)).*(A));
g = x1./x2 +g;
end
g
out = g;      
  
end

