function [ y ] = New_trial_kummer2_Pf( Pf )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

std_dev = 10.^(1./10);
%Pf = 0.1;
u = 5;

lambda = zeros(1, length(Pf));
    for i = 1:length(Pf)
    lambda(i) = 2*(gammaincinv(1-Pf(i), u./2));
    end
lambda
mu =4;
c = 2;
A = gamma(1+2./c).^(c./2);
A
time  = [-1.2247 0 1.2247];
W = [0.2954 1.1816 0.2954];

alpha = zeros(1, length(time));
   for i = 1:3   %==3
    alpha(i) = exp(-1*c./2.*((sqrt(2)*std_dev.*(time(i))) + mu)); % mu == mean
   end
    alpha
Y = zeros(1, length(lambda));
for j = 1:length(lambda)
    
y = zeros(1,1);

    for k = 0:20 % summision of k
          
          for i = 1:3
        
            for L = 0:20
                       
            y1 = (check(u+k, lambda(j))*alpha(i)*W(i)*((-1).^L)*((A*alpha(i)).^L)*gamma(k+(c./2)+(c*L./2)));
            y1
             y2 = (gamma(u+k)*gamma(k+1)*gamma(L+1));
            y2
             y = y1./y2 +y;
            end
            
          end
        
    end
y = y*(c*A./(2*sqrt(pi)))
Y(j) = y;
end

plot(Pf, Y)

end

function [out1]  = check(n, lambda)
          out1 = zeros(1,1);
          for t = 0:n -1
          out1 = (lambda./2).^t/(gamma(t+1)) + out1;
          end
          out1 = out1*gamma(n)*exp(-lambda./2);
end