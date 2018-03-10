function [ y ] = New_Trial_MOG_Wl_exact(Gama)

gama = 10.^(Gama./10);

std_dev = 0.806;  % ==3.1623

mu = -3.914; % ==6.3096

u = 5;
c = 2.5;
A = gamma(1+2./c).^(c./2); %==approx 0.9

Y = zeros(1,length(gama));

syms w
for j = 1: length(gama)
x1 = c.*A.*(gama(j).^(c./2 -1))  ;
x2  = 2.*(w.^(c./2)).*exp(((gama(j)./w).^(c./2)).*A).*exp(((log(w) - mu).^(2))./(2.*std_dev.^(2))).*std_dev.*w.*sqrt(2.*pi);
out = x1./x2;      
y = int(out ,w, 0, 50);
Y(j) = y;
end

plot(Gama, Y)
end
