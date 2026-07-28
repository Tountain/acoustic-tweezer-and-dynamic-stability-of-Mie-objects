function Bsl= sBessel(n,x,type)
%%
% Spherical Bessel Function
% n--order
% x--independent variable (SHOULD be row vector)
% type--the first kind or the second kind
%%

if type == 1    % Eq.(6.57)
    Bsl = (pi/2./(x+eps)).^0.5.*besselj(n+0.5, x);
elseif type == 2
    Bsl = (pi/2./(x+eps)).^0.5.*bessely(n+0.5, x);
end


% if type == 1 %% NOTE: from figure 6.6 of E.G.W@1999
%     if (n == 0 && x == 0) 
%         Bsl = 1;
%     end
% end
if type == 1 %% NOTE: from figure 6.6 of E.G.W@1999
    Bsl(find(n == 0 & x == 0)) = 1;
end

% Bsl(find(Bsl <= -10^5)) = -10^5;

% check table can be found in Eq.(6.60) and Eq.(6.61) in BOOK E.G.Williams@1999@

%% for n<0, the results of 'besselj' and 'bessely' can be validated by Eqs. (4.30) and (4.31) in BOOK E.G.Williams@1999@
% if (n+0.5) >= 0
%     if type == 1    % Eq.(6.57)
%         Bsl = (pi/2./(x+eps)).^0.5.*besselj(n+0.5, x);
%     elseif type == 2
%         Bsl = (pi/2./(x+eps)).^0.5.*bessely(n+0.5, x);
%     end
% else
%     order = -(n+0.5);
%     besselj_temp = cos(order*pi)*besselj(order, x) - sin(order*pi)*besselj(order, x);
%     bessely_temp = cos(order*pi)*bessely(order, x) - sin(order*pi)*bessely(order, x);
%     if type == 1   
%         Bsl = (pi/2./(x+eps)).^0.5.*besselj_temp;
%     elseif type == 2
%         Bsl = (pi/2./(x+eps)).^0.5.*bessely_temp;
%     end
% end