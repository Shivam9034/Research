function [ out ] = foxy_h(p,x, z)


L = 2;
c = 1.5;
        %% Contour preparation:
epsilon = 10^1.2;
Sups = min((L+ c./2)./(c./2)); Infs = max(0);
if(isempty(Sups) && isempty(Infs))
WPx=1;
elseif(isempty(Sups) && ~isempty(Infs))
WPx = Infs +epsilon;
elseif(~isempty(Sups) && isempty(Infs))
WPx = Sups -epsilon;
else
WPx = (Sups + Infs)/2;% s between Sups and Infs
end


fun = @(s)(gamma(-s).*gamma(p + (x./2) + ((x./2).*s))).*(z.^s);
out = (1./2*pi*1i).*integral(fun, WPx -1i*10, WPx +1i*10);
end
