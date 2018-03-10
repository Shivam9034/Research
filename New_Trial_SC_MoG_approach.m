function [ y ] = New_Trial_SC_MoG_approach(Gamath)

gama = 10.^(Gamath./10);
std_dev_e = 1;  % ==1.20
std_dev_m = 2;  % ==1.20

mu_e = 0; % ==0.406
mu_m = 0;
c = 2;
W = [2.23E-13	4.40E-10	1.09E-07	7.80E-06	2.28E-04	0.003243773	0.024810521	0.109017206	0.286675505	0.46224367	0.46224367	0.286675505	0.109017206	0.024810521	0.003243773	2.28E-04	7.80E-06	1.09E-07	4.40E-10	2.23E-13];
time  = [-5.38748089	-4.60368245	-3.94476404	-3.347854567	-2.788806058	-2.254974002	-1.738537712	-1.234076215	-0.737473729	-0.245340708	0.245340708	0.737473729	1.234076215	1.738537712	2.254974002	2.788806058	3.347854567	3.94476404	4.60368245	5.38748089];
alpha_e = zeros(1, length(time));
for i = 1:length(time)    
    alpha_e(i) = exp(-1*c./2.*((sqrt(2)*std_dev_e.*(time(i))) + mu_e)); % mu == mean
end
alpha_e

alpha_m = zeros(1, length(time));
for i = 1:length(time)    
    alpha_m(i) = exp(-1*c./2.*((sqrt(2)*std_dev_m.*(time(i))) + mu_m)); % mu == mean
end
alpha_m

Y = zeros(1,length(gama));

for j = 1: length(gama)

x2 = zeros(1,1);
for h = 1:length(W)
     for k =1:length(W)
out1 =  W(k).*W(h);
out1
out2 = alpha_e(h) + alpha_m(k).*(1+gama(j)).^(c./2); 
     
     x2 = (1./(2.*pi)).*(out1./out2) +x2; 
  
     end
     
end

x1 = 1./sqrt(pi).*sum(W);
x1
y = x1 - x2;
y
Y(j) = y;
end

plot(Gamath, Y)
end

