function [out] = Marcum_Q(a, b, M) % M = u 
out = zeros(1,1);

for k = 1-M:10
    out  = ((a./b)^k)*(besseli(k, a*b)) + out;
end
out = exp(-(a^2 + b^2)./2)*out ;

end

