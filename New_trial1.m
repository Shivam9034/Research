function [ y ] = New_trial1( mu)

%alpha == matrix 1xlength(time) == no need to add parameter
%A == constant ==NP
%c ==constant ==shape parameter ==NP
%std_dev == constant 
%time == matrix 1xlength(time)
%mu  = constant == mean ==prameter
%L ==matrix 1xinf   on deifned length == no need to ad parameter
%u ==constant ==order
%lambda == constant === therhold
%W == matrix 1xlength(time) == gauss harmite parameter

std_dev = 10.^(1./10);
Pf = 0.1;
u = 5;
lambda = 2*(gammaincinv(1-Pf, u./2));

c = 1.5;
A = gamma(1+2./c).^(c./2);

time  = [-1.2247 0 1.2247];
W = [0.2954 1.1816 0.2954];

Y = zeros(1, length(mu));   % see it is capital letter Y, New Y so that we can point our output toY for the plot..
for n = 1: length(mu)

alpha = zeros(1, length(time));   % initializing alpha
for i = 1:length(time)
    
alpha(i) = exp(-1*c./2.*((sqrt(2)*std_dev.*(time(i))) + mu(n))); % mu == mean

end

d = (c./2).*(A./sqrt(pi));
y = zeros(1,1);
for i = 1:length(time)
    exp1 = zeros(1,1);
    for L  = 0:15       %% we need to put some value instead of infnity

       exp1 = gamma_comb(L, u, lambda).*(foxy_h(L, c, A.*alpha(i))) +exp1;
       
    end
    y = alpha(i).*W(i).*exp1 +y;
end    

y = d.*(y);
        
Y(1,n) = y;

end
Y
plot(mu, Y);
xlabel('Mean value');
ylabel('probabiity');

function [ O ] = gamma_comb(p, x, h)
 O = (gammainc(p+(x./2) , h./2))./(gamma(p+1).*(gamma(p+x)));
end

% making the fox h function with H(m, n, p, q ==1,1,1,1)
% here p = L, x = c, z = A.*alpha
function [ out ] = foxy_h(p,x, z)
fun = @(s)(gamma(-s).*gamma(p + (x./2) + ((x./2).*s))).*(z.^s);
out = (1./2*pi*1i).*integral(fun, 3, 4);      %changes are require in limit
end
end
