function [db_bs_coeff_translation, db_bs_coeff_partial_r, db_bs_coeff_partial_theta, db_bs_coeff_partial_phi] = ...
            beam_shape_coeff_partial_using_database(nn, mm, db_bs_coeff, Snmvu_1, Snmvu_partial_r, Snmvu_partial_theta, Snmvu_partial_phi)
%%
% This function is used to determine the translated beam-shape coefficients
% of the standard beam-shape coefficients and the partial derivatives of
% the translated beam-shape coefficients using database
% "database_tranlation_coeff_partial.m".  
% The standard beam-shape coefficients 'db_bs_coeff' is calculated from
% function "beam_shape_coeff.m".
% Code structure refers to code "beam_shape_coeff_partial.m"
%%

[N, ~] = size(db_bs_coeff); N = N - 1;

%% ONLY for single layer media and single particle 

% 
db_bs_coeff_translation = 0;
for nu = 0:N
    for mu = -nu:nu
        %[Snmvu_1, ~] = Snmvu_coeff(nu, mu, nn, mm, k*r, theta, phi);
        db_bs_coeff_translation = db_bs_coeff_translation + ...
            db_bs_coeff(nu+1, nu+mu+1) * Snmvu_1(nu+1, nu+mu+1, nn+1, nn+mm+1);
    end
end

%
db_bs_coeff_partial_r = 0;
db_bs_coeff_partial_theta = 0;
db_bs_coeff_partial_phi = 0;

for nu = 0:N
    for mu = -nu:nu
        %[Snmvu_partial_r, Snmvu_partial_theta, Snmvu_partial_phi] = Snmvu_coeff_partial(nu, mu, nn, mm, k, r, theta, phi);
        db_bs_coeff_partial_r = db_bs_coeff_partial_r + db_bs_coeff(nu+1, mu+nu+1) * Snmvu_partial_r(nu+1, nu+mu+1, nn+1, nn+mm+1);
        db_bs_coeff_partial_theta = db_bs_coeff_partial_theta + db_bs_coeff(nu+1, mu+nu+1) * Snmvu_partial_theta(nu+1, nu+mu+1, nn+1, nn+mm+1);
        db_bs_coeff_partial_phi = db_bs_coeff_partial_phi + db_bs_coeff(nu+1, mu+nu+1) * Snmvu_partial_phi(nu+1, nu+mu+1, nn+1, nn+mm+1);
    end
end

%%