function s_vui = scattering_coeff_irregular_multi_bodies(BC, db_bs_coeff)
%%
% this function calculate the scattering coefficients of multi irregular
% bodies "s_nmi" under Hard wall (Neumann) Boundary Conditions.
% For irregular body, the scattering coefficients are related to the
% beam-shape coefficient 'db_bs_coeff'.
% The irregular bodies maps to spherical body using mapping coefficients
% 'Cn'.
%
% Note: this function is only valid for 'irregular_body == 1' meanwhile
% 'multi_particle == 1'.
%%

[row, ~] = size(db_bs_coeff);
N = row-1;

%% Building the database of the separate translation coefficients between i-th and j-th objects 'Snmvu_ij'
% using "database_tranlation_coeffs.m"
[db_filename_ts] = database_translation_coeffs(N);
load([db_filename_ts, '.mat'], 'Snmvu_1', 'Snmvu_2', 'r_ij');   % Snmvu -> S_(n,v)^(m,u) and i-th toward j-th

% delete([db_filename_ts, '.mat']);
  
%% Building the database of the structure coefficients for i-th object 'Gam_nmvu_i' and 'Amd_nmvu_i'
% referring to "structure_coefficients.m"
% Gam_nmvu_ii -> Gam_(n,m)^(v,u) of the ii-th particle
[db_filename_str] = database_structure_coeffs(N);
load([db_filename_str, '.mat'], 'Gam_nmvu_ii', 'Gam_u_nmvu_ii', 'Amd_nmvu_ii', 'Amd_u_nmvu_ii', 'particle_number');

% delete([db_filename_str, '.mat']);
   
%% Building the database of the beam-shape-structure coefficients for i-th object 'A_vui'
% referring to "multi_particle_beam_shape_coeff.m"
% Note: the probe particle is set to 'ii == 1', the other indice (ii ~= 1)
% as the source particles 
a_nmi_col = multi_particle_beam_shape_coeff_irregular(Snmvu_1, N, db_bs_coeff);       % the beam shape coefficients of (n,m) for i-th particle have been rearranged as a column vector
a_nmi = reshape(a_nmi_col, (N+1)^2, particle_number);       % rearrange the beam shape coefficients as a matrix, the i-th column represents the translated beam shape coefficients of the i-th particle.
if strcmp(BC, 'rigid') ~= 1
    A_vui = multi_particle_beam_structure_coeff(a_nmi, Gam_nmvu_ii);    
elseif strcmp(BC, 'rigid') == 1
    A_vui = multi_particle_beam_structure_coeff(a_nmi, Gam_u_nmvu_ii);    
end

% rearrange the matrix format of beam-shape-structure coefficients to a
% single column 'A'
A = zeros(particle_number * (N+1)^2, 1);
index = 0;
for ii = 1 : particle_number
    for n_s = 0 : N
        for m_s = -n_s : n_s
            index = index + 1;
            A(index) = A_vui(n_s+1, n_s+m_s+1, ii);     % a column vector, for matrix calculation convenience.
        end
    end
end

%% The coefficient matrix 'F' for the scalar scattering coefficients of multi irregular bodies 's_vuj'
% Note: the probe particle is set to 'ii == 1', the other indice (ii ~= 1)
% as the source particles

% Referred to my note, the indice (v,u), (n',m'), (n,m) here represent by
% (v,u)->(nu,mu), (n',m')->(n_s,m_s), (n,m)->(nn,mm).
% for coefficients B
% ===================== single thread version =============================
% B_nm_nsms_jl_numu = zeros(N+1, 2*N+1, N+1, 2*N+1, particle_number, particle_number, N+1, 2*N+1);
% for nn = 0 : N
%     for mm = -nn : nn
%         for n_s = 0 : N
%             for m_s = -n_s : n_s
%                 for ll = 1 : particle_number
%                     for jj = 1 : particle_number        
%                         
%                         if jj ~= ll             % 'll' means the probe particle currently is set to the ll-th particle
%                             for nu = 0 : N
%                                 for mu = -nu : nu
%                                     if strcmp(BC, 'rigid') ~= 1
%                                           B_nm_nsms_jl_numu(nn+1, nn+mm+1, n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1) = ...
%                                               Gam_nmvu_ii(nn+1, nn+mm+1, n_s+1, n_s+m_s+1, ll) * Snmvu_2(nu+1, nu+mu+1, nn+1, nn+mm+1, jj, ll);
%                                     elseif strcmp(BC, 'rigid') == 1
%                                           B_nm_nsms_jl_numu(nn+1, nn+mm+1, n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1) = ...
%                                               Gam_u_nmvu_ii(nn+1, nn+mm+1, n_s+1, n_s+m_s+1, ll) * Snmvu_2(nu+1, nu+mu+1, nn+1, nn+mm+1, jj, ll);
%                                     end
%                                 end
%                             end
%                         end
%                         
%                     end
%                 end
%             end
%         end
%     end
%     fprintf('Coefficients B Preparing %d%% \n',round(100*nn/N));
% end
% ===================== single thread version =============================

% ========================= parallel version ==============================
% B_nm_nsms_jl_numu -> B_(n,m)^{(n',m')(jl)}(v,u)
B_nm_nsms_jl_numu = zeros(N+1, 2*N+1, N+1, 2*N+1, particle_number, particle_number, N+1, 2*N+1);
for nn = 0 : N
    index_nn = nn+1;
    particle_number = particle_number;
    parfor mm = -nn : nn
        B_nsms_jl_numu = zeros(N+1, 2*N+1, particle_number, particle_number, N+1, 2*N+1);
        for n_s = 0 : N
            for m_s = -n_s : n_s
                for ll = 1 : particle_number
                    for jj = 1 : particle_number        
                        
                        if jj ~= ll             % 'll' means the probe particle currently is set to the ll-th particle
                            for nu = 0 : N
                                for mu = -nu : nu
                                    if strcmp(BC, 'rigid') ~= 1
                                        B_nsms_jl_numu(n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1) = ...
                                            Gam_nmvu_ii(nn+1, nn+mm+1, n_s+1, n_s+m_s+1, ll) * Snmvu_2(nu+1, nu+mu+1, nn+1, nn+mm+1, jj, ll);
                                    elseif strcmp(BC, 'rigid') == 1
                                        B_nsms_jl_numu(n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1) = ...
                                            Gam_u_nmvu_ii(nn+1, nn+mm+1, n_s+1, n_s+m_s+1, ll) * Snmvu_2(nu+1, nu+mu+1, nn+1, nn+mm+1, jj, ll);
                                    end
                                end
                            end
                        end
                        
                    end
                end
            end
        end
        B_nm_nsms_jl_numu(index_nn, index_nn+mm, :, :, :, :, :, :) = B_nsms_jl_numu;
    end
    fprintf('Coefficients B Preparing %d%% \n',round(100*nn/N));
end
% ========================= parallel version ==============================

% for coefficients D and E
% DE_nsms_jl_numu -> if j ~= l, D^{(n',m')(jl)}(v,u), if j == l, E^{(n',m')(ll)}(v,u). 
DE_nsms_jl_numu = zeros(N+1, 2*N+1, particle_number, particle_number, N+1, 2*N+1);
for n_s = 0 : N
    for m_s = -n_s : n_s
        for ll = 1 : particle_number
            for jj = 1 : particle_number        
                for nu = 0 : N
                    for mu = -nu : nu
                        
                        if jj ~= ll         % 'll' means the probe particle currently is set to the ll-th particle
                            % for coefficients D
                            for nn = 0 : N
                                for mm = -nn : nn
                                    DE_nsms_jl_numu(n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1) = ...
                                        DE_nsms_jl_numu(n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1) + B_nm_nsms_jl_numu(nn+1, nn+mm+1, n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1);
                                end
                            end
                        else
                            % for coefficients E
                            if strcmp(BC, 'rigid') ~= 1
                                DE_nsms_jl_numu(n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1) = ...
                                    Amd_nmvu_ii(nu+1, nu+mu+1, n_s+1, n_s+m_s+1, ll);
                            elseif strcmp(BC, 'rigid') == 1
                                DE_nsms_jl_numu(n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1) = ...
                                    Amd_u_nmvu_ii(nu+1, nu+mu+1, n_s+1, n_s+m_s+1, ll);
                            end
                        end
                        
                    end
                end
            end
        end
    end
    fprintf('Coefficients D and E Preparing %d%% \n',round(100*n_s/N));
end

% % for coefficients E
% E_nsms_l_numu = zeros(N+1, 2*N+1, particle_number, N+1, 2*N+1);
% for n_s = 0 : N
%     for m_s = -n_s : n_s
%         for ll = 1 : particle_number        % 'll' means the probe particle currently is set to the ll-th particle
%             
%             for nu = 0 : N
%                 for mu = -nu : nu
%                     E_nsms_l_numu(n_s+1, n_s+m_s+1, ll, nu+1, nu+mu+1) = ...
%                         Amd_nmvu_ii(nu+1, nu+mu+1, n_s+1, n_s+m_s+1, ll);
%                 end
%             end
%             
%         end
%     end
% end



%% Solving the governing equations for the scalar scattering coefficients 'G'
% referring to "interaction_beam_shape_coeff.m"

% rearrange the matrix format of coefficients and DE to a single row vector
% 'F_nsms_ii' for each indice (n',m',ll).
F_nsms_ii = zeros(N+1, 2*N+1, particle_number, particle_number * (N+1)^2);
for ll = 1 : particle_number
    for n_s = 0 : N
        for m_s = -n_s : n_s
            
            index = 0;
            for jj = 1 : particle_number
                for nu = 0 : N
                    for mu = -nu : nu
                        index = index + 1;
                        F_nsms_ii(n_s+1, n_s+m_s+1, ll, index) = DE_nsms_jl_numu(n_s+1, n_s+m_s+1, jj, ll, nu+1, nu+mu+1);
                    end
                end
            end
            
        end
    end
end

% for coefficients matrix 'F': a matrix with 'particle_number * (N+1)^2' column and 'particle_number * (N+1)^2' row 
F = zeros(particle_number * (N+1)^2, particle_number * (N+1)^2);
index = 0;
for ll = 1 : particle_number
    for n_s = 0 : N
        for m_s = -n_s : n_s
            index = index + 1;
            F(index, :) = F_nsms_ii(n_s+1, n_s+m_s+1, ll, :);
        end
    end
end

% figure(1); 
% pclr = pcolor(abs((F(:,:,1))));
% set(pclr, 'LineStyle','none');colorbar;
% currentCmap = colormap(hot);
% colormap(flipud(currentCmap));
% caxis([-0,10^5]);
% % axis([0,147,0,147]);
% % axis([0,192,0,192]);
% axis([0,243,0,243]);
% % axis([0,300,0,300]);
% % axis([0,363,0,363]);
% axis equal;

% figure(2);
% pclr = surf(abs(inv(F(:,:,1))));
% set(pclr, 'LineStyle','none');colorbar; 
% caxis([-10,10]);

% figure(3);
% pclr = plot(abs(A),'o');
% axis([0,1000]);

fprintf('The condition number of the coefficient matrix F: %e \n', cond(F));
if rcond(F) < 10^(-10)
    fprintf('Warning: the condition number of the coefficient matrix F is too high.\n');
end
   
G = F \ (-A);             % equivalent to 'G = inv(F) * A'
% G = lsqminnorm(F, -A);      % equivalent to 'G = inv(F) * A'
% G = pinv(F)*(-A);      % equivalent to 'G = inv(F) * A'
  
% rearrange the "G" based on the database projecting rule
s_vui = zeros(N+1, 2*N+1, particle_number);
index = 0;
% for ii = 1 : particle_number
%     for nn = 0 : N
%         for mm = -nn : nn
%             index = index + 1;
%             if a_nmi_col(index) == 0
%                 s_vui(nn+1, nn+mm+1, ii) = 0;
%             else
%                 s_vui(nn+1, nn+mm+1, ii) = G(index)/a_nmi_col(index);
%             end
%         end
%     end
% end
for ii = 1 : particle_number
    for nn = 0 : N
        for mm = -nn : nn
            index = index + 1;
            s_vui(nn+1, nn+mm+1, ii) = G(index);    % scattering coefficients, not scalar scattering coefficients.
        end
    end
end
% figure(2); 
% pclr = pcolor(abs(s_vui(:, :, 1)));
% set(pclr, 'LineStyle','none');
% caxis([-5,5]); 

%%
