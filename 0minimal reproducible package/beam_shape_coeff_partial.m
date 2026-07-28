function [db_bs_coeff_translation, db_bs_coeff_partial_r, db_bs_coeff_partial_theta, db_bs_coeff_partial_phi] ...
                = beam_shape_coeff_partial(nn, mm, db_bs_coeff, k, r, theta, phi)
%%
% This function is used to determine the translated beam-shape coefficients
% of the standard beam-shape coefficients and the partial derivatives of
% the translated beam-shape coefficients. 
% The standard beam-shape coefficients 'db_bs_coeff' is calculated from
% function "beam_shape_coeff.m".
%%

[N, ~] = size(db_bs_coeff); N = N - 1;


%% ONLY for single layer media and single particle 

% 
db_bs_coeff_translation = 0;
for nu = 0:N
    for mu = -nu:nu
        [Snmvu_1, ~] = Snmvu_coeff(nu, mu, nn, mm, k*r, theta, phi);
        db_bs_coeff_translation = db_bs_coeff_translation + ...
            db_bs_coeff(nu+1, nu+mu+1) * Snmvu_1;
    end
end

%
db_bs_coeff_partial_r = 0;
db_bs_coeff_partial_theta = 0;
db_bs_coeff_partial_phi = 0;

for nu = 0:N
    for mu = -nu:nu
        [Snmvu_partial_r, Snmvu_partial_theta, Snmvu_partial_phi] = Snmvu_coeff_partial(nu, mu, nn, mm, k, r, theta, phi);
        db_bs_coeff_partial_r = db_bs_coeff_partial_r + db_bs_coeff(nu+1, mu+nu+1) * Snmvu_partial_r;
        db_bs_coeff_partial_theta = db_bs_coeff_partial_theta + db_bs_coeff(nu+1, mu+nu+1) * Snmvu_partial_theta;
        db_bs_coeff_partial_phi = db_bs_coeff_partial_phi + db_bs_coeff(nu+1, mu+nu+1) * Snmvu_partial_phi;
    end
end
 
 
% checking
checking = 0;

if checking == 1
    delta_r = r/100000 + eps; delta_theta = theta/1000 + eps; delta_phi = phi/1000 + eps;
    db_bs_coeff1 = db_bs_coeff_translation;
    
    db_bs_coeff2 = 0;
    for nu = 0:N
        for mu = -nu:nu
            [Snmvu_1, ~] = Snmvu_coeff(nu, mu, nn, mm, k*(r+delta_r), theta, phi);
            db_bs_coeff2 = db_bs_coeff2 + ...
                db_bs_coeff(nu+1, nu+mu+1) * Snmvu_1;
        end
    end
    
    db_bs_coeff3 = 0;
    for nu = 0:N
        for mu = -nu:nu
            [Snmvu_1, ~] = Snmvu_coeff(nu, mu, nn, mm, k*r, theta+delta_theta, phi);
            db_bs_coeff3 = db_bs_coeff3 + ...
                db_bs_coeff(nu+1, nu+mu+1) * Snmvu_1;
        end
    end
    
    db_bs_coeff4 = 0;
    for nu = 0:N
        for mu = -nu:nu
            [Snmvu_1, ~] = Snmvu_coeff(nu, mu, nn, mm, k*r, theta, phi+delta_phi);
            db_bs_coeff4 = db_bs_coeff4 + ...
                db_bs_coeff(nu+1, nu+mu+1) * Snmvu_1;
        end
    end
    
    db_bs_coeff_partial_r_check = (db_bs_coeff2 - db_bs_coeff1) / delta_r;
    db_bs_coeff_partial_theta_check = (db_bs_coeff3 - db_bs_coeff1) / delta_theta;
    db_bs_coeff_partial_phi_check = (db_bs_coeff4 - db_bs_coeff1) / delta_phi;
    db_bs_coeff_partial_r_check
    db_bs_coeff_partial_theta_check
    db_bs_coeff_partial_phi_check
end
% Note that checking should not work in theta == 0, while from theoretical
% deviation, Ynm_partial_theta == 0 when theta == 0, regardless with n and
% m.
