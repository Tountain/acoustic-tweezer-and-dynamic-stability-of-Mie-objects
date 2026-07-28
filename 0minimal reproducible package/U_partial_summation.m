function U_partial = U_partial_summation(trans_bs_coeff, trans_bs_coeff_partial, db_s_coeff, db_s_coeff_partial, ii,jj, nn, mm, n_s, m_s)
%%
% Sum up the 7 element from "U_partial_elements.m", and return the real and
% imaginary parts of the partial derivatives of radiation effects.       
%%

[rol, col] = size(db_s_coeff);

part1 = U_partial_elements(trans_bs_coeff_partial, trans_bs_coeff, ii, jj, ones(rol, col), db_s_coeff, nn, mm, n_s, m_s);
part2 = U_partial_elements(trans_bs_coeff, trans_bs_coeff_partial, ii, jj, ones(rol, col), db_s_coeff, nn, mm, n_s, m_s);
part3 = U_partial_elements(trans_bs_coeff, trans_bs_coeff, ii, jj, ones(rol, col), db_s_coeff_partial, nn, mm, n_s, m_s);
part4 = U_partial_elements(trans_bs_coeff_partial, trans_bs_coeff, ii, jj, db_s_coeff, db_s_coeff, nn, mm, n_s, m_s);
part5 = U_partial_elements(trans_bs_coeff, trans_bs_coeff, ii, jj, db_s_coeff_partial, db_s_coeff, nn, mm, n_s, m_s);
part6 = U_partial_elements(trans_bs_coeff, trans_bs_coeff_partial, ii, jj, db_s_coeff, db_s_coeff, nn, mm, n_s, m_s);
part7 = U_partial_elements(trans_bs_coeff, trans_bs_coeff, ii, jj, db_s_coeff, db_s_coeff_partial, nn, mm, n_s, m_s);

U_partial = part1 + part2 + part3 + part4 + part5 + part6 + part7;

% U_partial_real = real(U_partial);
% U_partial_imag = imag(U_partial);

%%
