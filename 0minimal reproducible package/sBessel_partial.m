function Bsl_p = sBessel_partial(n,x,type)
%%
% Obtain the differentiate of Spherical Bessel Function 
% n--order
% x--independent variable (SHOULD be row vector)
% type--the first kind or the second kind
% Based on Eq.(6.69) at E.G.W@1999
%%

Bsl_p = sBessel(n-1,x,type) - (n+1)*sBessel(n,x,type)/(x+eps);

% if n ==0 && x == 0 && type == 1 %% NOTE: from figure 6.6 of E.G.W@1999
%     Bsl_p = 0;
% end
if type == 1 %% NOTE: from figure 6.6 of E.G.W@1999
    Bsl_p(find(n == 0 & x == 0)) = 0;
end

% checking
checking = 0;

if checking == 1
    delta_x = x/100000 + eps;
    Bs1 = sBessel(n,x,type);
    Bs2 = sBessel(n,x+delta_x,type);
    Bsl_p_check = (Bs2 - Bs1) / delta_x;
    Bsl_p_check
end
