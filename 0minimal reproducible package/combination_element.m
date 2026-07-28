function combi_coeff = combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn, mm, vv, uu)
%%
% This function is used to calculate the combination coefficients.
% The combination coefficients involve the contributions from the translate
% beam-shape coefficients and the scalar scattering coefficients.
%%

% Nt = size(trans_bs_coeff, 2);

%%

alpha_ii_nm = real(trans_bs_coeff{ii}(nn, mm));
beta_ii_nm = imag(trans_bs_coeff{ii}(nn, mm));

alpha_jj_vu = real(trans_bs_coeff{jj}(vv, uu));
beta_jj_vu = imag(trans_bs_coeff{jj}(vv, uu));

sigma_nm = real(db_s_coeff(nn, mm));
epsilon_nm = imag(db_s_coeff(nn, mm));

sigma_vu = real(db_s_coeff(vv, uu));
epsilon_vu = imag(db_s_coeff(vv, uu));


tau_ii_nm = alpha_ii_nm * sigma_nm - beta_ii_nm * epsilon_nm;
eta_ii_nm = beta_ii_nm * sigma_nm + alpha_ii_nm * epsilon_nm;
tau_jj_vu = alpha_jj_vu * sigma_vu - beta_jj_vu * epsilon_vu;
eta_jj_vu = beta_jj_vu * sigma_vu + alpha_jj_vu * epsilon_vu;


%%

Real = (alpha_ii_nm + tau_ii_nm) * tau_jj_vu + (beta_ii_nm + eta_ii_nm) * eta_jj_vu;
Imag = (beta_ii_nm + eta_ii_nm) * tau_jj_vu - (alpha_ii_nm + tau_ii_nm) * eta_jj_vu;

combi_coeff = Real + Imag * 1i;

%%