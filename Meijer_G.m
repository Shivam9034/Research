function out = Meijer_G(an, ap, bm, bq, z)
%***** Integrand definition *****
syms s;
F = @(s) (GammaProd(bm,s).*GammaProd(1-an,-s).*(z.^(-s)))./(GammaProd(ap,s).*GammaProd(1-bq,-s));

%***** Contour definition *****
Sups = min(bm);
Infs = -max(1-an); % cs
cs = (Sups + Infs)/2;% s between Sups and Infs
W = 50; % W
%***** Bivariate Meijer G *****
out = real((1/(2*pi*1i))*integral( F,cs-1i*W,cs+1i*W)); 

function output = GammaProd(p,z)
[pp, zz] = meshgrid(p,z);
if (isempty(p)) 
  output = ones(size(z));
else
  output = reshape(prod(gammac(pp+zz),2),size(z));
end
end
end