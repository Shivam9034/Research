
function [ y ] = New_trial_NL_H_exact( Gama )

gama=10.^(Gama./10);

m=0.5;

wo=10.^(5./10);

beta=0.5;

y = zeros(1,1);

syms w

for j = 1:length(gama)
    
    
    x1=m.^m*(gama(j).^(m-1)).*exp(-(m.*gama(j)./w)).*(w.^(beta-1)).*exp(-(w./wo));
    x2=gamma(m).*gamma(beta).*(wo.^beta).*(w.^m);
    
    g=x1./x2;
    g
    
    y=int(g, w, 0, 50);
   
    
Y(j) = y;
          
end    
   
plot(Gama, Y)
end



