function [ y ] = gammaa(z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
syms t;
Funn = (exp(-t))*(t.^(z-1));
y = int(Funn, t, 0, 50);

end

