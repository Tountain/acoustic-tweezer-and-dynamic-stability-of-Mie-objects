function [db_bs_coeff_rot_translation, db_bs_coeff_rot_partial_alpha, db_bs_coeff_rot_partial_beta, db_bs_coeff_rot_partial_gamma] ...
                = beam_shape_coeff_rot_partial(nn, mm, db_bs_coeff, alpha_pd, beta_pd, gamma_pd)
%%
% This function is used to determine the translated beam-shape coefficients
% of the standard beam-shape coefficients and the partial derivatives of
% the translated beam-shape coefficients. 
% The standard beam-shape coefficients 'db_bs_coeff' is calculated from
% function "beam_shape_coeff.m".
%%

[N, ~] = size(db_bs_coeff); N = N - 1;


%% ONLY for single layer media and single particle 

db_bs_coeff_rot_translation = 0;
for ll = -nn:nn
    D_nml = rotation_matrix([-gamma_pd, -beta_pd, -alpha_pd], nn, ll, mm); %  theta_rotation
    db_bs_coeff_rot_translation = db_bs_coeff_rot_translation + db_bs_coeff(nn+1, ll+nn+1) * D_nml;
end


%
db_bs_coeff_rot_partial_alpha = 0;
db_bs_coeff_rot_partial_beta = 0;
db_bs_coeff_rot_partial_gamma = 0;

for ll = -nn:nn
    [D_nml_partial_alpha, D_nml_partial_beta, D_nml_partial_gamma] = rotation_matrix_partial([-gamma_pd, -beta_pd, -alpha_pd], nn, ll, mm); %  theta_rotation
    db_bs_coeff_rot_partial_alpha = db_bs_coeff_rot_partial_alpha + db_bs_coeff(nn+1, ll+nn+1) * D_nml_partial_alpha;
    db_bs_coeff_rot_partial_beta = db_bs_coeff_rot_partial_beta + db_bs_coeff(nn+1, ll+nn+1) * D_nml_partial_beta;
    db_bs_coeff_rot_partial_gamma = db_bs_coeff_rot_partial_gamma + db_bs_coeff(nn+1, ll+nn+1) * D_nml_partial_gamma;
end
 
 
% checking
checking = 0;

if checking == 1
    delta_alpha = alpha_pd/10000000000 + 1/10000000000; delta_beta = beta_pd/10000000000 + eps; delta_gamma = gamma_pd/10000000000 + 1/10000000000;
    db_bs_coeff1 = db_bs_coeff_rot_translation;
    
    db_bs_coeff2 = 0;
    for ll = -nn:nn
        D_nml = rotation_matrix([-gamma_pd+delta_alpha, -beta_pd, -alpha_pd], nn, ll, mm); %  theta_rotation
        db_bs_coeff2 = db_bs_coeff2 + db_bs_coeff(nn+1, ll+nn+1) * D_nml;
    end
    
    db_bs_coeff3 = 0;
    for ll = -nn:nn
        D_nml = rotation_matrix([-gamma_pd, -beta_pd+delta_beta, -alpha_pd], nn, ll, mm); %  theta_rotation
        db_bs_coeff3 = db_bs_coeff3 + db_bs_coeff(nn+1, ll+nn+1) * D_nml;
    end
    
    db_bs_coeff4 = 0;
    for ll = -nn:nn
        D_nml = rotation_matrix([-gamma_pd, -beta_pd, -alpha_pd+delta_gamma], nn, ll, mm); %  theta_rotation
        db_bs_coeff4 = db_bs_coeff4 + db_bs_coeff(nn+1, ll+nn+1) * D_nml;
    end
    
    db_bs_coeff_rot_partial_alpha_check = (db_bs_coeff2 - db_bs_coeff1) / delta_alpha;
    db_bs_coeff_rot_partial_beta_check = (db_bs_coeff3 - db_bs_coeff1) / delta_beta;
    db_bs_coeff_rot_partial_gamma_check = (db_bs_coeff4 - db_bs_coeff1) / delta_gamma;
    db_bs_coeff_rot_partial_alpha_check
    db_bs_coeff_rot_partial_beta_check
    db_bs_coeff_rot_partial_gamma_check
end
% Note that checking should not work in theta == 0, while from theoretical
% deviation, Ynm_partial_theta == 0 when theta == 0, regardless with n and
% m.
