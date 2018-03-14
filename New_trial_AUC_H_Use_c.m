function [ y ] = New_trial_AUC_H_Use_c( c )

std_dev = 1;

u = 5;
mu =0;
%c = 2;

%time  = [-1.2247 0 1.2247];
%W = [0.2954 1.1816 0.2954];
time  = [-5.38748089	-4.60368245	-3.94476404	-3.347854567	-2.788806058	-2.254974002	-1.738537712	-1.234076215	-0.737473729	-0.245340708	0.245340708	0.737473729	1.234076215	1.738537712	2.254974002	2.788806058	3.347854567	3.94476404	4.60368245	5.38748089];


%W = [0.1667 0.032 0.1259 ];   % normalize value..
W = [2.23E-13	4.40E-10	1.09E-07	7.80E-06	2.28E-04	0.003243773	0.024810521	0.109017206	0.286675505	0.46224367	0.46224367	0.286675505	0.109017206	0.024810521	0.003243773	2.28E-04	7.80E-06	1.09E-07	4.40E-10	2.23E-13];


   
Y = zeros(1, length(c));
for j = 1:length(c)
A = gamma(1+2./c(j)).^(c(j)./2);
    alpha = zeros(1, length(time));
   for i = 1:length(time)   %==3
    alpha(i) = exp(-1*c(j)./2.*((sqrt(2)*std_dev.*(time(i))) + mu)); % mu == mean
   end   
y = zeros(1,1);
  
          for i = 1:length(time)
        
            for L = 1:10
            H = H_Function_NewTrial(1-L-c(j)./2,c(j)./2,[], [],0, 1, [], [], A.*alpha(i));
            K = kummer_21(1,2*u + L,u+1,1/2);
            y1 = alpha(i)*W(i)*gamma(2*u +L)*((1/2).^(L+u))*H*K;
            
            y2 = gamma(u+L)*gamma(L+1)*(2.^(u-1))*gamma(u)*u;
            
            y = y1./y2 +y;
            end
            
          end
y = y*(c(j)*A./(2*sqrt(pi)));
y
Y(j) = y;
          
end    
   
plot(c, Y)
end

function [out1]  = check(n, lambda)
          out1 = zeros(1,1);
          for t = 0:n -1
          out1 = (lambda./2).^t/(gamma(t+1)) + out1;
          end
          out1 = out1*gamma(n)*exp(-lambda./2);
end
