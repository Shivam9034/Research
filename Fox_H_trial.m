function out = Fox_H_trial(an, An, ap, Ap,bm, Bm ,bq, Bq,z)
% an, An, ap, Ap,bm, Bm ,bq, Bq,z
%% Integrand definition
F = @(s)(GammaProd(bm,Bm,s).* GammaProd(1-an,-An,s).* z.^-s )./ (GammaProd(1-bq,-Bq,s).* GammaProd(ap,Ap,s));
%% Contour preparation:
epsilon = 10^1.2;
Sups = min((1-an)./An); Infs = max(-bm./Bm);
if(isempty(Sups) && isempty(Infs))
WPx=1;
elseif(isempty(Sups) && ~isempty(Infs))
WPx = Infs +epsilon;
elseif(~isempty(Sups) && isempty(Infs))
WPx = Sups -epsilon;
else
WPx = (Sups + Infs)/2;% s between Sups and Infs
end
%% integration:
infity = 10;
out = (1/(2i*pi))*integral(F,WPx-1i*infity, WPx+1i*infity);
return
%% ***** GammaProd subfunction *****
function fnl = GammaProd(p,x,s)
   fnl  = 1;
    for i = 1:length(p)
       fnl  = gamma(p(i) - x(i).*s)*fnl;
    end
end 
end
    
