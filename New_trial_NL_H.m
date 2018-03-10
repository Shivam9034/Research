function [ y ] = New_trial_NL_H( Gama )

gama=10.^(Gama./10);

m=2;

wo=10.^(10./10);

beta=2;

y = zeros(1,1);

for j = 1:length(gama)
    
    z=(m.*gama(j)./wo);
    x1=m.^m*(gama(j).^(m-1)).*H_Function_NewTrial([1-beta+m],[-1],[], [],[0], [1], [], [], z);
    x2=gamma(m).*gamma(beta).*wo.^m;
    
    y=x1./x2;
    y
    
Y(j) = y;
          
end    
   
plot(Gama, Y)
end
