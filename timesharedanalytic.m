%PDF for MRC time shared multipath/shadowed and unshadowed fading
%Mean MM=1dB
%Variance VV=0.25dB
%Shape parameter c
%Rice factor kk(dB)
%SNR YY(dB)

c=0.5;
MM=1;
m=db2mag(MM);
VV=0.25;
v=db2mag(VV);
YY=0:0.1:20;
y=db2mag(YY);
kk=0;
k=db2mag(kk);

%Holtzman approximation of composite weibull lognormal fading
a=gamma(1+(2/c));
a1=exp(m);
a2=exp(m+(sqrt(3)*v));
a3=exp(m-(sqrt(3)*v));
P=(2/3)*(c/2)*((a/a1)^(c/2)).*((y.^((c/2)-1))).*exp(-((a/a1).*y).^(c/2));
P1=(1/6)*(c/2)*((a/a2)^(c/2)).*((y.^((c/2)-1))).*exp(-((a/a2).*y).^(c/2));
P2=(1/6)*(c/2)*((a/a3)^(c/2)).*((y.^((c/2)-1))).*exp(-((a/a3).*y).^(c/2));
R=P+P1+P2;
%Unshadowed Rician fading
pu=(1+k)/10;
ru=exp(-k);
tu=exp(-pu.*y);
z=2*sqrt(k*pu.*y);
I=besseli(0, z);
T=pu*ru*tu.*I;%Rician fading

Rup=(0.4*T)+(0.6*R);                                            %Time shared combined PDF


%PDF Analytic for MRC time shared multipath/shadowed and unshadowed fading


syms w                                                          %w=SNR for log normal 

p=(a/w)^(c/2);
p1=((a/w).*y).^(c/2);
p2=exp(-p1);
LN=exp(-((10*log10(w)-m)^2)/(2*v^2))*(1/(sqrt(2*pi)*v*w));      % Log Normal Distribution
f=4.3429*(c/2)*p*p2*LN.*(y.^((c/2)-1));                        %Composite weibull log normal
l=int(f,w,0,50);                                                %PDF for weibull over log normal


Rup1=(0.4*T)+(0.6*l);                                           %Time shared combined PDF(analytic)


semilogy(YY,Rup,'r',YY,Rup1,'b')