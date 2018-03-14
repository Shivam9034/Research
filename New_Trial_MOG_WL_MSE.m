function [ y ] = New_Trial_MOG_WL_MSE(Gama)
gama = 10.^(Gama./10);
std_dev = 0.806;  % ==3.1623
 mu = -3.914;
c = 2.5;
lambda = [15.98 13.44 11.78 10.47 9.34 8.29 7.26 6.17 4.86 0]; % for u == 5
%time  = [-2.020 -0.9585 0 0.9585 2.020];
time  = [-5.38748089	-4.60368245	-3.94476404	-3.347854567	-2.788806058	-2.254974002	-1.738537712	-1.234076215	-0.737473729	-0.245340708	0.245340708	0.737473729	1.234076215	1.738537712	2.254974002	2.788806058	3.347854567	3.94476404	4.60368245	5.38748089];
%W = [0.01995 0.3936 0.9453 0.3936 0.01995];
%W = [0 0.4038 1 0.4038 0];
W = [2.23E-13	4.40E-10	1.09E-07	7.80E-06	2.28E-04	0.003243773	0.024810521	0.109017206	0.286675505	0.46224367	0.46224367	0.286675505	0.109017206	0.024810521	0.003243773	2.28E-04	7.80E-06	1.09E-07	4.40E-10	2.23E-13];
A = gamma(1+2./c).^(c./2); %==approx 0.9

alpha = zeros(1, length(time));
for i = 1:length(time)    
    alpha(i) = exp(-1*c./2.*((sqrt(2)*std_dev.*(time(i))) + mu)); % mu == mean
end
alpha     

Y = zeros(1,length(gama));

for j = 1: length(gama)

    g = zeros(1,1);
    for i = 1:length(alpha)
    x1 = alpha(i).*W(i).*gama(j).^(c./2 -1);
    
    
x2 = exp((gama(j).^(c./2)).*A.*alpha(i));
    
    g = x1./x2 +g;
    
    end
    y = c.*A ./(2.*sqrt(pi)).*g;
Y(j) = y;

end
Y
X = zeros(1,length(gama));

syms w
for j = 1: length(gama)
x1 = c.*A.*(gama(j).^(c./2 -1))  ;
x2  = 2.*(w.^(c./2)).*exp(((gama(j)./w).^(c./2)).*A).*exp(((log(w) - mu).^(2))./(2.*std_dev.^(2))).*std_dev.*w.*sqrt(2.*pi);
out = x1./x2;      

y = int(out ,w, 0, 50);

X(j) = y;
end
X

MSE = ((X - Y).^(2))./length(gama);
MSE = sum(MSE);
MSE
end

