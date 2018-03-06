function [ Out_m ] = for_m(m,bq, Bq)
Out_m = 1;
syms s;
for k = 1:m
   
    Out_m = (gamma(bq(k) - Bq(k).*s)).*(Out_m);
    
end

end