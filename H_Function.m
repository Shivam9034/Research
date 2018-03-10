function [ y ] = H_Function( m, n, ap, Ap, bq, Bq, Z )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
syms s;
fun  = (((for_m(m,bq, Bq)).*(for_n(n,ap, Ap)).*(Z.^s))./((for_q(m,bq, Bq)).*(for_p(n,ap, Ap))));
y = int(fun,s, 3 , 4 );

end

function [ Out_m ] = for_m(m,bq, Bq)

if m == 0
    Out_m = 1;
else
    Out_m = 1;
    syms s;
    for k = 1:m
    Out_m = (gamma(bq(k) + Bq(k).*s)).*(Out_m); 
    end
end
end

function [ Out_m ] = for_n(n,ap, Ap)

if n == 0
    Out_m =1;
else
    Out_m = 1;
    syms s;
    for k = 1:n
    Out_m = (gamma(1 - ap(k) + Ap(k).*s)).*(Out_m); 
    end
end
end

function [ Out_m ] = for_q(m,bq, Bq)

if m+1 >= length(bq)
    Out_m =1;
else
Out_m = 1;
syms s;
for k = m +1 : length(bq)
    Out_m = (gamma(1 - bq(k) +  Bq(k).*s)).*(Out_m); 
end
end
end

function [ Out_m ] = for_p(n,ap, Ap)

if n+1 >= length(ap)
    Out_m = 1;
else
Out_m = 1;
syms s;
for k = n+1 : length(ap)
    Out_m = (gamma(ap(k) + Ap(k).*s)).*(Out_m); 
end
end
end
