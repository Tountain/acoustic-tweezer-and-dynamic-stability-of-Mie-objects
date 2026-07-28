function [Y_nm_partial_theta, Y_nm_partial_phi] = sHarmonics_partial(n,m,theta,phi)
%%
% Obtain the deviation of spherical harmonics function 
% n--lower order
% m--upper order
% theta--polar angle, radian mesurement (SHOULD be row vector)   [0,pi]
% phi--azimuthal angle, radian mesurement (SHOULD be row vector) [0,2*pi]
% Y_nm-- row:length(phi), column:length(theta)
%%

P_nm = sLegendre(n,m,theta);        % 1 * length(theta)
exp_phi = (exp(1i * m * phi'));      % length(phi) * 1
P_nm_partial = sLegendre_partial(n,m,theta);
coeff = sqrt(((2*n+1)/(4*pi)) * (factorial(n-m)/factorial(n+m)));

Y_nm_partial_theta = coeff * exp_phi * (P_nm_partial .* sin(-theta));
Y_nm_partial_phi = coeff * exp_phi * P_nm * 1i*m;


% checking
checking = 0;

if checking == 1
    delta_theta = theta/1000 + eps; delta_phi = phi/1000 + eps;
    Ynm1 = sHarmonics(n,m, theta, phi);
    Ynm2 = sHarmonics(n,m, theta+delta_theta, phi);
    Ynm3 = sHarmonics(n,m, theta, phi+delta_phi);
    Y_nm_partial_theta_check = (Ynm2 - Ynm1) / delta_theta;
    Y_nm_partial_phi_check = (Ynm3 - Ynm1) / delta_phi;
    Y_nm_partial_theta_check
    Y_nm_partial_phi_check
end
% Note that checking should not work in theta == 0, while from theoretical
% deviation, Ynm_partial_theta == 0 when theta == 0, regardless with n and
% m.