function [Snmvu_partial_r, Snmvu_partial_theta, Snmvu_partial_phi] = Snmvu_coeff_partial(n, m, nu, mu, k, r, theta, phi)
%%
% this function is used to compute the partial derivative of Translation
% coefficient "Snm_vu" of the first kind ONLY.
% Based on definition Eq.(3.80) for Translation coefficient "Snm_vu" of the
% first kind at P.A.MARTIN@2006@Multiple_Scattering.
%
% Snmvu_coeff -> S_(n,v)^(m,u) or S_(n,nu)^(m,mu); (v -> nu; u -> mu)
%
% Gaunt coefficient "gaunt value" and Spherical Harmonics "Ynm" are
% obtained by sub-function "gaunt_coeff" and "sHarmonics".
%%

kr = k*r;
coeff = 4*pi*1i^(nu-n);
q_lower = abs(n-nu);
q_upper = n+nu;

sum_partial_r = 0;
sum_partial_theta = 0;
sum_partial_phi = 0;

for qq = q_lower:q_upper
    if abs(mu-m)<=qq
        temp = (1i)^qq * (-1)^m * gaunt_coeff(n,m,nu,-mu,qq);
        [Y_nm_partial_theta, Y_nm_partial_phi] = sHarmonics_partial(qq,mu-m,theta,phi);
        sum_partial_r = sum_partial_r + temp * ( conj(sHarmonics(qq,mu-m,theta,phi)) * sBessel_partial(qq,kr,1) * k );
        sum_partial_theta = sum_partial_theta + temp * ( conj(Y_nm_partial_theta) * sBessel(qq,kr,1) );
        sum_partial_phi = sum_partial_phi + temp * ( conj(Y_nm_partial_phi) * sBessel(qq,kr,1) );
    end 
end
Snmvu_partial_r = coeff * sum_partial_r;
Snmvu_partial_theta = coeff * sum_partial_theta;
Snmvu_partial_phi = coeff * sum_partial_phi;


% check the correction of Translation coefficient "Snmvu" can use
% Eq.(3.82), Eq.(3.89), Eq.(3.93), Eq.(3.94), Eq.(3.95), Eq.(3.98),
% Eq.(3.99), Eq.(3.100) and Lamma 3.29 at
% P.A.MARTIN@2006@Multiple_Scattering


% checking
checking = 0;

if checking == 1
    delta_r = r/100000 + eps; delta_theta = theta/1000 + eps; delta_phi = phi/1000 + eps;
    Snmvu1 = Snmvu_coeff(n, m, nu, mu, k*r, theta, phi);
    Snmvu2 = Snmvu_coeff(n, m, nu, mu, k*(r+delta_r), theta, phi);
    Snmvu3 = Snmvu_coeff(n, m, nu, mu, k*r, theta+delta_theta, phi);
    Snmvu4 = Snmvu_coeff(n, m, nu, mu, k*r, theta, phi+delta_phi);
    Snmvu_partial_r_check = (Snmvu2 - Snmvu1) / delta_r;
    Snmvu_partial_theta_check = (Snmvu3 - Snmvu1) / delta_theta;
    Snmvu_partial_phi_check = (Snmvu4 - Snmvu1) / delta_phi;
    Snmvu_partial_r_check
    Snmvu_partial_theta_check
    Snmvu_partial_phi_check
end
% Note that checking should not work in theta == 0, while from theoretical
% deviation, Ynm_partial_theta == 0 when theta == 0, regardless with n and
% m.
