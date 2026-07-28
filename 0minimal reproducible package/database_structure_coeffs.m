function [db_filename_str] = database_structure_coeffs(N)
%%
% build the database for 2 kinds of structure coefficients (Gam and Amd),
% and their radial partial derivative (Gam_partial_u and Amd_partial_u).
% Therefore, saving the time for repeating calculating these
% two types of structure coefficients.
% 
% NOTE: this database is created only for multiple irregular particle
% system 'irregular_body == 1' meanwhile 'multi_particle == 1'.
%%

parameters_names;

if exist([db_filename_str, '.mat']) ~= 0        % if already exist the database, 
    return;                                     % then do not create again for saving time.
end


%% the Gam. Gam_partial_u, Amd, and Amd_partial_u coefficients

particles_Cartesian_data;

% Gam_nmvu_ii -> Gam_(n,m)^(v,u) of the ii-th particle
Gam_nmvu_ii = zeros(N+1, 2*N+1, N+1, 2*N+1, particle_number);       %% 'ii' means the ii-th particle index of the structure coefficients
Gam_u_nmvu_ii = zeros(N+1, 2*N+1, N+1, 2*N+1, particle_number);
Amd_nmvu_ii = zeros(N+1, 2*N+1, N+1, 2*N+1, particle_number);
Amd_u_nmvu_ii = zeros(N+1, 2*N+1, N+1, 2*N+1, particle_number);

for ii = 1 : particle_number
    for n_s = 0 : N
        for m_s = -n_s : n_s

            [Gam_nm_ii, Gam_u_nm_ii, Amd_nm_ii, Amd_u_nm_ii] = ...
                    structual_functions_multi_bodies(n_s, m_s, Cn_multi{ii}, N);    % (n_s, m_s) or (vv, uu) [corresponding to (n', m') in text] are the upper indice of the structure coefficients

            Gam_nmvu_ii(:, :, n_s+1, n_s+m_s+1, ii) = Gam_nm_ii;  
            Gam_u_nmvu_ii(:, :, n_s+1, n_s+m_s+1, ii) = Gam_u_nm_ii;
            Amd_nmvu_ii(:, :, n_s+1, n_s+m_s+1, ii) = Amd_nm_ii;
            Amd_u_nmvu_ii(:, :, n_s+1, n_s+m_s+1, ii) = Amd_u_nm_ii;
 
        end
    end
    fprintf('Structure Coefficients Database Preparing %d%% (%d-th of %d particles). \n', ...
        round(100*ii/particle_number), ii, particle_number);
end
  
  
% % ============= Parallel version 5 for speed up =============
% for nn = 0:N
%     indices_1 = nn+1;
%     particle_number = particle_number;
%     parfor mm = -nn:nn
%         
%         temp_Snmvu_1 = zeros(N+1, 2*N+1, particle_number, particle_number);
%         temp_Snmvu_2 = zeros(N+1, 2*N+1, particle_number, particle_number);   
%         for nu = 0:N
%             for mu = -nu:nu
%                 indices_3 = nu+1;
%                 indices_4 = nu+mu+1;
%                 for ii = 1:particle_number
%                     for jj = 1:particle_number
%                         if ii == jj
%                             temp_Snmvu_1(indices_3, indices_4, ii, jj) = 0;
%                             temp_Snmvu_2(indices_3, indices_4, ii, jj) = 0;
%                         else
%                             [temp_Snmvu_1(indices_3, indices_4, ii, jj), temp_Snmvu_2(indices_3, indices_4, ii, jj)] = ...
%                                 Snmvu_coeff(nn, mm, nu, mu, kr_ij(ii, jj), theta_ij(ii,jj), phi_ij(ii,jj));
%                         end
%                     end
%                 end
%             end
%         end
%         Snmvu_1(indices_1, indices_1 + mm, :, :, :, :) = temp_Snmvu_1;
%         Snmvu_2(indices_1, indices_1 + mm, :, :, :, :) = temp_Snmvu_2;
%         
%     end
%     fprintf('Translation Coefficients Database Preparing %d%% \n', ...
%             round(100*nn/N));
% end

%% save the database by meaningful name

save([db_filename_str, '.mat'], 'Gam_nmvu_ii', 'Gam_u_nmvu_ii', 'Amd_nmvu_ii', 'Amd_u_nmvu_ii', 'particle_number');

%%