function [ y ] = New_Trial_SC_R_approach(Gamath)

gama = 10.^(Gamath./10);

std_dev_e = 0;  % ==1.20
std_dev_m = 0;  % ==1.20

c = 2;
mu_e = 1; % ==0.406
mu_m = 2;
Xi = [mu_m  mu_m+sqrt(3)*std_dev_m  mu_m-sqrt(3)*std_dev_m];
Yj = [mu_e  mu_e+sqrt(3)*std_dev_e  mu_e-sqrt(3)*std_dev_e];

Di = [0.66 0.1667 0.1667];
Dj = [0.66 0.1667 0.1667];

Y = zeros(1,length(gama));
for j = 1: length(gama)

y = zeros(1,1);
for h = 1:length(Xi)
     for k =1:length(Yj)
out1 =  Di(h).*Dj(k).*((gamma(1+2./c)./exp(Xi(h))).^(c./2)).*(1+gama(j)).^(c./2); 
out1
out2 = ((gamma(1+2./c)./exp(Yj(k))).^(c./2)) + ((gamma(1+2./c)./exp(Xi(h))).^(c./2)).*(1+gama(j)).^(c./2); 
out2     
    y = out1./out2 +y; 
     end
end
Y(j) = y;
end

plot(Gamath, Y)
end

