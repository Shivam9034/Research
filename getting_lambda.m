function [y] = getting_lambda(Pf, u)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%lambda = ones(1, length(Pf));

[xy] = deal({lamb});
   
a = check(u, xy);% == Pf*gamma(u);
%y = solve(a, lamb);
r = a == Pf*gamma(u);
y = solve(r, xy);
 
end