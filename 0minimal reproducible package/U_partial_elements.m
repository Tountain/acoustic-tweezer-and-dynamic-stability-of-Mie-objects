function U_partial = U_partial_elements(trans_bs_coeff1, trans_bs_coeff2, ii, jj, db_s_coeff1, db_s_coeff2, nn, mm, n_s, m_s)
%%
% this code is used to calculate the real and imaginary parts of ' anm^(i)*
% an'm'^(j) * snm * sn'm' ', which includes 7 different formats, and
% summarize all of them form the radiation effect of i and j transducers: '
% U_nm^n'm' (ij) '.   
%%

alpha_ii_nm = real(trans_bs_coeff1{ii}(nn, mm));
beta_ii_nm = imag(trans_bs_coeff1{ii}(nn, mm));

alpha_jj_nms = real(trans_bs_coeff2{jj}(n_s, m_s));
beta_jj_nms = imag(trans_bs_coeff2{jj}(n_s, m_s));

sigma_nm = real(db_s_coeff1(nn, mm));
epsilon_nm = imag(db_s_coeff1(nn, mm));

sigma_nms = real(db_s_coeff2(n_s, m_s));
epsilon_nms = imag(db_s_coeff2(n_s, m_s));

%%

A = (alpha_ii_nm * alpha_jj_nms + beta_ii_nm * beta_jj_nms) * sigma_nms + (beta_ii_nm * alpha_jj_nms - alpha_ii_nm * beta_jj_nms) * epsilon_nms;
B = (beta_ii_nm * alpha_jj_nms - alpha_ii_nm * beta_jj_nms) * sigma_nms - (alpha_ii_nm * alpha_jj_nms + beta_ii_nm * beta_jj_nms) * epsilon_nms;

%%

Real = A * sigma_nm - B * epsilon_nm;
Imag = A * epsilon_nm + B * sigma_nm;

U_partial = Real + Imag * 1i;

%%
