function [ y ] = New_trial_Pd_Equation( Pf )

std_dev = 1;
u = 5;
mu =0;
c = 2.5;
A = gamma(1+2./c).^(c./2);
lambda = [15.98 13.44 11.78 10.47 9.34 8.29 7.26 6.17 4.86 0];
time  = [-1.2247 0 1.2247];
%time  = [-5.38748089	-4.60368245	-3.94476404	-3.347854567	-2.788806058	-2.254974002	-1.738537712	-1.234076215	-0.737473729	-0.245340708	0.245340708	0.737473729	1.234076215	1.738537712	2.254974002	2.788806058	3.347854567	3.94476404	4.60368245	5.38748089];

W = [0.2954 1.1816 0.2954];
%W = [2.23E-13	4.40E-10	1.09E-07	7.80E-06	2.28E-04	0.003243773	0.024810521	0.109017206	0.286675505	0.46224367	0.46224367	0.286675505	0.109017206	0.024810521	0.003243773	2.28E-04	7.80E-06	1.09E-07	4.40E-10	2.23E-13];


alpha = zeros(1, length(time));
   for i = 1:length(time)   %==3
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

z = (c*A./(2*sqrt(pi)))*Marcum_Q(sqrt(2*gama), sqrt(lambda(j)), u).*z;

fun  =  z;

y = int(fun,gama,  0, 100) ;
Y(j) = y;

end
plot(Pf, Y)
end
function [out] = Marcum_Q(a, b, M) % M = u 
out = zeros(1,1);

for k = 1-M:100
    out  = ((a./b).^k)*(besseli(k, a.*b)) + out;
end
out = exp(-(a.^2 + b.^2)./2).*out ;

end


