function [ y ] = H_Function_NewTrial(an,An,ap, Ap,bm, Bm, bq, Bq, Z )

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Contour preparation:
epsilon = 1;
Sups = max((1-an)./An); Infs = min(-bm./Bm);
if(isempty(Sups) && isempty(Infs))
WPx=1;
elseif(isempty(Sups) && ~isempty(Infs))
WPx = Infs -epsilon;
elseif(~isempty(Sups) && isempty(Infs))
WPx = Sups +epsilon;
else
WPx = (Sups + Infs)/2;% s between Sups and Infs
end

infity = 10;
for i = 1:length(Z)
fun  = @(s) ((for_m(s, bm, Bm)).*(for_n(s,an, An)).*(Z(i).^(-s))./((for_q(s,bq, Bq)).*(for_p(s,ap, Ap)))); 
y = (1./(2*pi*1i))*integral(fun, WPx-1i*100, WPx+1i*100);
end


function [ Out_m ] = for_m(x, bm, Bm)
m = length(bm);
if m == 0
    Out_m = 1;
else
    Out_m = 1;
    %syms s;
    for k = 1:m
    Out_m = (gammaZ(bm(k) + Bm(k).*x)).*(Out_m); 
    end
end
Out_m
end

function [ Out_n ] = for_n(x, an, An)
n = length(an);
if n == 0
    Out_n =1;
else
    Out_n = 1;
    %syms s;
    for k = 1:n
    Out_n = (gammaZ(1 - an(k) - An(k).*x)).*(Out_n); 
    end
end
Out_n
end

function [ Out_q ] = for_q(x,bq, Bq)
q = length(bq);
if q == 0
    Out_q =1;
else
Out_q = 1;
%syms s;
for k = 1 : q
    Out_q = (gammaZ(1 - bq(k) - Bq(k).*x)).*(Out_q); 
end
end
Out_q
end

function [ Out_p ] = for_p(x, ap, Ap)
p = length(ap);
if p == 0
    Out_p = 1;
else
Out_p = 1;
%syms s;
for k = 1:p
    Out_p = (gammaZ(ap(k) + Ap(k).*x)).*(Out_p); 
end
end
Out_p
end


end

function [f] = gammaZ(z)
% GAMMA  Gamma function valid in the entire complex plane.
%        Accuracy is 15 significant digits along the real axis
%        and 13 significant digits elsewhere.
%        This routine uses a superb Lanczos series
%        approximation for the complex Gamma function.
%
%        z may be complex and of any size.
%        Also  n! = prod(1:n) = gamma(n+1)
%
%usage: [f] = gamma(z)
%       
%tested on versions 6.0 and 5.3.1 under Sun Solaris 5.5.1
%
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931-944
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996
%            W. J. Cody "An Overview of Software Development for Special
%            Functions", 1975
%
%see also:   GAMMA GAMMALN GAMMAINC PSI
%see also:   mhelp GAMMA
%
%Paul Godfrey
%pgodfrey@intersil.com
%http://winnie.fit.edu/~gabdo/gamma.txt
%Sept 11, 2001

siz = size(z);
z=z(:);
zz=z;

f = 0.*z; % reserve space in advance

p=find(real(z)<0);
if ~isempty(p)
   z(p)=-z(p);
end

% 15 sig. digits for 0<=real(z)<=171
% coeffs should sum to about g*g/2+23/24

g=607/128; % best results when 4<=g<=5

c = [  0.99999999999999709182;
      57.156235665862923517;
     -59.597960355475491248;
      14.136097974741747174;
      -0.49191381609762019978;
        .33994649984811888699e-4;
        .46523628927048575665e-4;
       -.98374475304879564677e-4;
        .15808870322491248884e-3;
       -.21026444172410488319e-3;
        .21743961811521264320e-3;
       -.16431810653676389022e-3;
        .84418223983852743293e-4;
       -.26190838401581408670e-4;
        .36899182659531622704e-5];

%Num Recipes used g=5 with 7 terms
%for a less effective approximation

z=z-1;
zh =z+0.5;
zgh=zh+g;
%trick for avoiding FP overflow above z=141
zp=zgh.^(zh*0.5);

ss=0.0;
for pp=size(c,1)-1:-1:1
    ss=ss+c(pp+1)./(z+pp);
end

%sqrt(2Pi)
sq2pi=  2.5066282746310005024157652848110;
f=(sq2pi*(c(1)+ss)).*((zp.*exp(-zgh)).*zp);

f(z==0 | z==1) = 1.0;

%adjust for negative real parts
if ~isempty(p)
   f(p)=-pi./(zz(p).*f(p).*sin(pi*zz(p)));
end

%adjust for negative poles
p=find(round(zz)==zz & imag(zz)==0 & real(zz)<=0);
if ~isempty(p)
   f(p)=Inf;
end

f=reshape(f,siz);

end





