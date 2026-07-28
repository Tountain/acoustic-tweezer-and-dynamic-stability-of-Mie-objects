function [db_filename] = database_beam_scattering_coeffs()
%%
% build the database for "beam_shape_coeffcient" and
% "scattering_coefficient" to avoid repeating call function
% "beam_shape_coeff" and "scattering_coefficient", which will cost tons of
% time. 
% 
% database will automatically be created if current folder has not a
% corresponding database. (check the database file exists or not by the
% filename)
% database maximum size: db_size_nn = 30 (seeing "parameters.m");
%%

parameters_names;

if exist([db_filename, '.mat']) ~= 0        % if already exist the database, 
    return;                                 % then do not create again for saving time.
end

% definiation the size of the coefficients
db_bs_coeff = zeros(db_size_nn+1, 2*db_size_nn+1);
if multi_particle == 0      % single particle system
    if irregular_body == 0      % single spherical object
        db_s_coeff = zeros(db_size_nn+1, 1);
    elseif irregular_body == 1  % single nonspherical object
        db_s_coeff = zeros(db_size_nn+1, 2*db_size_nn+1);
    end
else                        % multi-particle system
    particles_Cartesian_data;
    if irregular_body == 0          % multi spherical objects
        db_s_coeff = zeros(db_size_nn+1, particle_number);
    elseif irregular_body == 1      % multi axisymmetric objects
        db_s_coeff = zeros(db_size_nn+1, 2*db_size_nn+1, particle_number);
    end
end

A_delay = 0;
phi_delay = 0;

%% determine the Truncation number of expansion terms: N

N = 0;          % Truncation number (based on "sn") for single-particle system
if multi_particle == 0
%     error = 1;
    % "cut-off theory" based on G.T.S@2011@IEEE p.299
%     while (error > 0.01) && (N < db_size_nn)
%         if irregular_body == 0
%             error = abs(db_s_coeff(N+1+1) / db_s_coeff(N+0+1));      
%             N = N + 1;
%         elseif irregular_body == 1
%             error = abs(sum(db_s_coeff(N+1+1, :)) / sum(db_s_coeff(N+0+1, :)));
%             N = N + 1;
%         end
%     end
    % get the small "N" from "N=ka+6" and "cut-off theory"
%     if N > round(fluid_k*particle_radius + 8)  
%         N = round(fluid_k*particle_radius + 8);
%     end
    N = db_size_nn;            % a larger "N", a "wider" plain wave
else
    % maximum "database size for beam-shape coeff, espacially plain wave,
    % as small database size will lose the plain wave information in
    % farfield, but for multi-particle system, we need farfield
    % information.
    N = db_size_nn;            % based on J.H@2016@IEEE, NOTE: the particle radius "ka<4"     
end

% avoid too many expansion terms, "db_size" is the maximum size of
% "beam_shape_coefficient" and "scattering_coefficient" in the database.
% if N > db_size_nn       
%     N = db_size_nn;
% end


%% database for "beam_shape_coeffcient" single layer media: db_bs_coeff
% mapping relation of "db_bs_coeff" and "beam_shape_coeffcient"
%   1.first subscript for "scattering_coefficient" is "nn" 
%     but "nn+1" for "db_s_coeff"
%   2.second subscript for "scattering_coefficient" is "mm" 
%     but "mm+nn+1" for "db_s_coeff"
%
% particularly, beam_shape_coeffcient(NN, MM)
% then, db_bs_coeff(NN+1, NN+MM+1)
% where NN = nn + num_C, NN = mm + num_D
%
% NOTE: after finishing position translation(num_C or num_D), the new
%       position NN and MM should meet "abs(MM) <= abs(NN)", because for
%       "beam_shape_coeffcient", we require "abs(mm) <= abs(nn)".

if multi_layer == 0
    
    if rot_matrix ~= 1
        
        for nn = 0:N
            for mm = -nn:nn
                db_bs_coeff(nn+1, mm+nn+1) = beam_shape_coeff(nn, mm);
            end
            fprintf('Beam-Shape and Scattering Coefficients Database Preparing %d%% \n',round(100*nn/N));
        end
        
    elseif rot_matrix == 1
          
        for nn = 0:N
            for mm = -nn:nn
                db_bs_coeff(nn+1, mm+nn+1) = beam_shape_coeff(nn, mm);
            end
            fprintf('Beam-Shape and Scattering Coefficients Database Preparing %d%% \n',round(100*nn/N));
        end
          
        if strcmp(wave_type, 'plain') == 1  % for plain wave ONLY, for phase array the equivalent beam-shape is given below.
            db_bs_coeff_rot = zeros(db_size_nn+1, 2*db_size_nn+1);
            %rotation matrix is introduced to rotate the beam-shape 
            %coefficients (which is equivalent to a rotational
            %transformation (from OCS to CCS).) 
            for nn = 0:N
                for mm = -nn:nn
                    sum_bs_ll = 0;
                    for ll = -nn:nn
                        D_nml = rotation_matrix([-theta_rotation(3), -theta_rotation(2), -theta_rotation(1)], nn, ll, mm); %  theta_rotation
                        sum_bs_ll = sum_bs_ll + db_bs_coeff(nn+1, ll+nn+1) * D_nml;
                    end
                    db_bs_coeff_rot(nn+1, mm+nn+1) = sum_bs_ll;
                end
            end   
        end  
        
%         db_bs_coeff = db_bs_coeff_rot;  % test the addition rotation theorem
           
    end
    
end  

%% database for "beam_shape_coeffcient" multi-layer medium: db_bs_coeff
% mapping relation of "db_bs_coeff" and "beam_shape_coeffcient"
%   1.first subscript for "scattering_coefficient" is "nn" 
%     but "nn+1" for "db_s_coeff"
%   2.second subscript for "scattering_coefficient" is "mm" 
%     but "mm+nn+1" for "db_s_coeff"
%
% particularly, beam_shape_coeffcient(NN, MM)
% then, db_bs_coeff(NN+1, NN+MM+1)
% where NN = nn + num_C, NN = mm + num_D
%
% NOTE: after finishing position translation(num_C or num_D), the new
%       position NN and MM should meet "abs(MM) <= abs(NN)", because for
%       "beam_shape_coeffcient", we require "abs(mm) <= abs(nn)".

if multi_layer == 1
    
    % prepare the amplitude field on the integration sperical surface 'U'
%     [U0, T] = total_amplitude_distribution_all_plain_wave_component();
    U = ASA_amplitude_field(U0, T);
%     load('inter_dist_source_layer50_AIR-rho2o5c680-PLANE_v2_t1.mat')
    
    for nn = 0:N
        for mm = -nn:nn
            
            db_bs_coeff(nn+1, mm+nn+1) = beam_shape_coeff_multi_layer(U, nn, mm);
%             db_bs_coeff(nn+1, mm+nn+1) = beam_shape_coeff_multi_layer(U_8_plane_30bauto08, nn, mm);
            
        end
        fprintf('Beam-Shape and Scattering Coefficients Database Preparing %d%% \n',round(100*nn/N));
    end
    
end
% for ii = 1:N
%     for jj = 1:2*N+1
%         if abs(db_bs_coeff(ii, jj)) > 2*abs(db_bs_coeff(1, 1))
%             db_bs_coeff(ii, jj) = 0;
%         end
%     end
% end

%% database for "equivalent beam_shape_coeffcient": equivalent_bs_coeff -> db_bs_coeff

if multi_layer == 0
    
    if rot_matrix ~= 1
        
        % translation of beam-shape coefficients for the Bessel beams.
        if (strcmp(wave_type, 'zero-Bessel') == 1 || strcmp(wave_type, 'non-zero-Bessel') == 1) && translated_Bessel_beam == 1
            [r_t, theta_t, phi_t] = ...
            coords_system_relative_positions_general(beam_source, [0, 0, 0]);
            if strcmp(fluid, 'co2') == 1
                fluid_k_bessel = real(fluid_k) + 1i * imag(fluid_k) / cos(BETA);     % for BETA is 75 degree
            end
            kr_t = fluid_k_bessel * r_t;

            bs_translation = zeros(size(db_bs_coeff));
            for nn = 0:N
                indices_1 = nn + 1;
                parfor mm = -nn:nn          % 'mm' is sliced variable
                    temp_bs_t = 0;
                    for nu = 0:N
                        for mu = -nu:nu
                            [Snmvu_1, ~] = Snmvu_coeff(nu, mu, nn, mm, kr_t, theta_t, phi_t);
                            temp_bs_t = temp_bs_t + ...
                                db_bs_coeff(nu+1, nu+mu+1) * Snmvu_1;
                        end
                    end
                    bs_translation(indices_1, indices_1 + mm) = temp_bs_t;
                end
                fprintf('Translated Beam-Shape Coefficient of Bessel Beams %d%% \n', ...
                    round(100*nn/N));
            end

            db_bs_coeff = bs_translation;
        end
    
        % equivalent beam-shape coefficients for the phase array.
        if strcmp(wave_type, 'phase_array_transducer') == 1 
            [equivalent_bs_q_coeff, transducer, A_delay, phi_delay] = phase_array_beam_shape_coeff(wave_type, N, db_bs_coeff, fluid_k, transducer, transducer_number, theta_rotation);
            db_bs_coeff = equivalent_bs_q_coeff;      % equivalent_bs_coeff -> db_bs_coeff
        elseif strcmp(wave_type, 'phase_array_transducer2') == 1 
            [~, transducer, A_delay, phi_delay] = phase_array_beam_shape_coeff(wave_type, N, db_bs_coeff, fluid_k, transducer, transducer_number, theta_rotation);
        else
            transducer = [0,0,0];
        end
        
    elseif rot_matrix == 1 
         
        if strcmp(wave_type, 'phase_array_transducer') == 1 
            [equivalent_bs_q_coeff, transducer, A_delay, phi_delay] = phase_array_beam_shape_coeff(wave_type, N, db_bs_coeff, fluid_k, transducer, transducer_number, [0,0,0]);
            db_bs_coeff = equivalent_bs_q_coeff;      % equivalent_bs_coeff -> db_bs_coeff
            %rotation of equivalent beam-shape coefficients
            db_bs_coeff_rot = zeros(db_size_nn+1, 2*db_size_nn+1);
            for nn = 0:N 
                for mm = -nn:nn
                    sum_bs_ll = 0;
                    for ll = -nn:nn
                        D_nml = rotation_matrix([-theta_rotation(3), -theta_rotation(2), -theta_rotation(1)], nn, ll, mm); %  theta_rotation
                        sum_bs_ll = sum_bs_ll + db_bs_coeff(nn+1, ll+nn+1) * D_nml;
                    end
                    db_bs_coeff_rot(nn+1, mm+nn+1) = sum_bs_ll;
                end
            end
        end
        
    end 
    
    fprintf('Equivalent Beam-Shape Coefficient is prepared well! \n');
    
elseif multi_layer == 1
    
    % left for future 
     
end


%% position dependance of beam-shape coefficients and its partial derivatives for the any wave-kinds. 

if derivative_field_is_needed == 1
     
    [r_pd, theta_pd, phi_pd] = ...
            coords_system_relative_positions_general([derivativeX, derivativeY, derivativeZ], [0, 0, 0]); % '- \vec{r}'
    db_bs_coeff_translation = zeros(db_size_nn+1, 2*db_size_nn+1);
    db_bs_coeff_partial_r = zeros(db_size_nn+1, 2*db_size_nn+1);
    db_bs_coeff_partial_theta = zeros(db_size_nn+1, 2*db_size_nn+1);
    db_bs_coeff_partial_phi = zeros(db_size_nn+1, 2*db_size_nn+1);
    fluid_k = fluid_k;
     
    for nn = 0:N
        indices_1 = nn + 1;
        parfor mm = -nn:nn          % 'mm' is sliced variable
            [db_bs_coeff_translation(indices_1, indices_1+mm), db_bs_coeff_partial_r(indices_1, indices_1+mm), db_bs_coeff_partial_theta(indices_1, indices_1+mm), db_bs_coeff_partial_phi(indices_1, indices_1+mm)] ...
                        = beam_shape_coeff_partial(nn, mm, db_bs_coeff, fluid_k, r_pd, theta_pd, phi_pd);
        end 
        fprintf('Translated Beam-Shape Coefficient and its partial derivatives (translation) %d%% \n', round(100*nn/N));
    end  
    
    db_bs_coeff = db_bs_coeff_translation;
      
end

if derivative_rotationfield_is_needed == 1
    
    % there lack a code to map (derivative_thetaX, derivative_thetaY,
    % derivative_thetaZ) to (alpha, beta, gamma), while if
    % derivative_thetaX == derivative_thetaZ == 0, alpha == 0, beta ==
    % derivative_thetaZ, gamma == 0. so if (derivative_thetaX,
    % derivative_thetaY, derivative_thetaZ) = (0,t,0), 
    % (alpha, beta, gamma) = (0,t,0). 
    if (derivative_thetaX == 0 && derivative_thetaZ == 0)
        alpha_pd = derivative_thetaX; beta_pd = derivative_thetaY; gamma_pd = derivative_thetaZ;
    else
        error('(derivative_thetaX, derivative_thetaY, derivative_thetaZ) != (alpha, beta, gamma)');
    end
    db_bs_coeff_rot_translation = zeros(db_size_nn+1, 2*db_size_nn+1);
    db_bs_coeff_rot_partial_alpha = zeros(db_size_nn+1, 2*db_size_nn+1);
    db_bs_coeff_rot_partial_beta = zeros(db_size_nn+1, 2*db_size_nn+1);
    db_bs_coeff_rot_partial_gamma = zeros(db_size_nn+1, 2*db_size_nn+1);
        
    for nn = 0:N 
        indices_1 = nn + 1;
        parfor mm = -nn:nn          % 'mm' is sliced variable
            [db_bs_coeff_rot_translation(indices_1, indices_1+mm), db_bs_coeff_rot_partial_alpha(indices_1, indices_1+mm), db_bs_coeff_rot_partial_beta(indices_1, indices_1+mm), db_bs_coeff_rot_partial_gamma(indices_1, indices_1+mm)] ...
                        = beam_shape_coeff_rot_partial(nn, mm, db_bs_coeff_rot, alpha_pd, beta_pd, gamma_pd);
        end  
        fprintf('Translated Beam-Shape Coefficient and its partial derivatives (rotation) %d%% \n', round(100*nn/N));
    end 
    
    db_bs_coeff_rot2 = db_bs_coeff_rot_translation;
    
end

%% database for "scattering_coefficient": db_s_coeff
% mapping relation of "db_s_coeff" and "scattering_coefficient"
%   subscript for "scattering_coefficient" is "nn" 
%   but "nn+1" for "db_s_coeff" 
%
% particularly, scattering_coefficient(NN)
% then, db_s_coeff(NN+1)
% where NN = nn + num_C

if multi_particle == 0      % single particle system
    
    if irregular_body == 0
        
        for nn = 0:db_size_nn   % single spherical object
            db_s_coeff(nn+1) = scattering_coefficient(nn); 
        end
        
    elseif irregular_body == 1  % single nonspherical object
        
        if rot_matrix ~= 1 && derivative_field_is_needed ~= 1
            % means that we does not introduce addition rotation theorem to
            % return the CCS to OCS, and another rotation tranformation is
            % needed before visulization the pressure field or the
            % radiation force/torque.
            snm_ib = scattering_coeff_irregular_body(BC, db_bs_coeff);
%             snm_ib_simplified = scattering_coeff_irregular_body_simplified(BC, db_bs_coeff);
            for nn = 0:db_size_nn
                for mm = -nn:nn
                    db_s_coeff(nn+1, mm+nn+1) = snm_ib(nn+1, mm+nn+1);
%                     db_s_coeff_simplified(nn+1, mm+nn+1) = snm_ib_simplified(nn+1, mm+nn+1);
                end
            end
            
        elseif rot_matrix ~= 1 && derivative_field_is_needed == 1
            % means that we does not introduce addition rotation theorem to
            % return the CCS to OCS, and another rotation tranformation is
            % needed before visulization the pressure field or the
            % radiation force/torque.
            snm_ib = scattering_coeff_irregular_body(BC, db_bs_coeff);
%             snm_ib_simplified = scattering_coeff_irregular_body_simplified(BC, db_bs_coeff);
            for nn = 0:db_size_nn
                for mm = -nn:nn
                    db_s_coeff(nn+1, mm+nn+1) = snm_ib(nn+1, mm+nn+1);
%                     db_s_coeff_simplified(nn+1, mm+nn+1) = snm_ib_simplified(nn+1, mm+nn+1);
                end
            end
             
            % position dependance of scalar scattering coefficients and its partial derivatives for the any wave-kinds.
            [db_s_coeff_partial_r, db_s_coeff_partial_theta, db_s_coeff_partial_phi] ...
                = scattering_coeff_irregular_body_partial(BC, db_bs_coeff, db_bs_coeff_partial_r, db_bs_coeff_partial_theta, db_bs_coeff_partial_phi, db_s_coeff);
            
        elseif rot_matrix == 1 && derivative_rotationfield_is_needed ~= 1
            snm_ib_rot = scattering_coeff_irregular_body(BC, db_bs_coeff_rot);
            db_s_coeff_rot = snm_ib_rot;
            %rotation matrix is introduced to rotate the scalar scattering
            %coefficients (which is equivalent to a rotational
            %transformation (from CCS to OCS).) 
            for nn = 0:db_size_nn  
                for mm = -nn:nn
                    sum_bss_coeff = 0;
                    for ll = -nn:nn
                        D_nml = rotation_matrix([theta_rotation(1), theta_rotation(2), theta_rotation(3)], nn, ll, mm); %  (theta_z, theta_y, theta_x)
                        sum_bss_coeff = sum_bss_coeff + db_bs_coeff_rot(nn+1, ll+nn+1) * db_s_coeff_rot(nn+1, ll+nn+1) * D_nml;
                    end
                    db_s_coeff_original(nn+1, mm+nn+1) = sum_bss_coeff;         % scattering coefficients
                    if db_bs_coeff(nn+1, mm+nn+1) == 0          % avoid 'db_bs_coeff(nn+1, mm+nn+1) == 0' and 'db_s_coeff(nn+1, mm+nn+1) == inf'.
                        db_s_coeff(nn+1, mm+nn+1) = 0;          % The reason why not '(db_bs_coeff(nn+1, mm+nn+1) + eps)' is because it will enlarge 'db_s_coeff(nn+1, mm+nn+1)', and make the pressure field and radiation force/torque mismatch.
                    else
                        db_s_coeff(nn+1, mm+nn+1) = sum_bss_coeff / db_bs_coeff(nn+1, mm+nn+1);     % scalar scattering coefficients
                    end
                    % if isinf(db_s_coeff(nn+1, mm+nn+1))         % avoid 'db_bs_coeff(nn+1, mm+nn+1) == 0' and 'db_s_coeff(nn+1, mm+nn+1) == inf'.
                    %     db_s_coeff(nn+1, mm+nn+1) = 0;          % The reason why not '(db_bs_coeff(nn+1, mm+nn+1) + eps)' is because it will enlarge 'db_s_coeff(nn+1, mm+nn+1)', and make the pressure field and radiation force/torque mismatch.
                    % end
                end
            end
            
        elseif rot_matrix == 1 && derivative_rotationfield_is_needed == 1
            snm_ib_rot2 = scattering_coeff_irregular_body(BC, db_bs_coeff_rot2);
            db_s_coeff_rot2 = snm_ib_rot2;
            
            % rotation dependance of scalar scattering coefficients and its partial derivatives for the any wave-kinds. 
            [db_s_coeff_rot_partial_alpha, db_s_coeff_rot_partial_beta, db_s_coeff_rot_partial_gamma] ...
                = scattering_coeff_irregular_body_partial(BC, db_bs_coeff_rot2, db_bs_coeff_rot_partial_alpha, db_bs_coeff_rot_partial_beta, db_bs_coeff_rot_partial_gamma, db_s_coeff_rot2);
        
            % NOTE: considering that 'theta_x+derivative_thetaX == 0',
            % 'theta_y+derivative_thetaY == 0', 'theta_z+derivative_thetaZ
            % == 0', so not another rotation matrix is needed in
            % 'db_bs_coeff_rot2' and 'db_s_coeff_rot2'.
            if (theta_x + derivative_thetaX ~= 0) || (theta_y + derivative_thetaY ~= 0) || (theta_z + derivative_thetaZ ~= 0)
                error('Remain the object is in its standard orientation. \n');
            else
                db_bs_coeff = db_bs_coeff_rot2;
                db_s_coeff = db_s_coeff_rot2;
            end
            
        end
        
    end
    
elseif multi_particle == 1                        % multi-particle system
    
    if irregular_body == 0
        for nn = 0:db_size_nn   % multiple spherical objects
            db_s_coeff(nn+1, :) = scattering_coefficient(nn);   % scalar scattering coefficients
        end
    elseif irregular_body == 1  % multiple nonspherical object
        if rot_matrix ~= 1
            error('For irregular multiple object problem, addition rotation theorem is REQUIRED (rotational transformation cannot handle)!\n');
        elseif rot_matrix == 1
            if theta_rotation(1) == 0 && theta_rotation(2) == 0 && theta_rotation(3) == 0
                % in this moment, we firstly restrict the problem as NO
                % particle rotation phenamena (theta_rotation = [0, 0, 0]). 
                s_nmi = scattering_coeff_irregular_multi_bodies(BC, db_bs_coeff);   % NOTE: this function includes many database building, such as database for translation coefficients and database for structure coefficients.
                db_s_coeff = s_nmi;      % scattering coefficients, NOT scalar scattering coefficients!!! Note: the probe particle is set to 'ii == 1'
            else
                
                % left for future
                
            end
        end
    end
    
end

%% save these databases by meaningful name

% r = max(max(abs(db_s_coeff)))/max(max(abs(db_bs_coeff)))
% f_ka  

% % visualiza the coefficients matrices
% sum_bs = bs_visualize_colormap(db_bs_coeff);
% sum_bs;
% figure(1); 
% figure;
% heatmap(abs(db_bs_coeff(:,:,1)));%,'GridVisible','off');
% caxis([0,12]); 
% bs_visualize_colormap;
% sum_bs = sum(sum(abs(db_bs_coeff(:,:,1))))
% figure(2); 
% heatmap(abs(db_bs_coeff_rot(:,:,1)));%,'GridVisible','off');
% caxis([0,12]);
% bs_visualize_colormap;
% sum_bs_rot = sum(sum(abs(db_bs_coeff_rot(:,:,1))))
% figure(3);
% figure; 
% heatmap(abs(db_s_coeff(:,:,1)));
% caxis([0,0.2]);
% figure; 
% heatmap(abs(db_s_coeff_simplified(:,:,1)));
% caxis([0,0.2]);
% figure(4);
% pclr = pcolor(abs(db_s_coeff_original(:,:,1)));
% set(pclr, 'LineStyle','none');  
% caxis([0,5]);  
    
if derivative_field_is_needed == 1
    save([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_bs_coeff_partial_r', 'db_bs_coeff_partial_theta', 'db_bs_coeff_partial_phi', ...
        'db_s_coeff', 'db_s_coeff_partial_r', 'db_s_coeff_partial_theta', 'db_s_coeff_partial_phi', 'transducer', 'A_delay', 'phi_delay');
elseif derivative_rotationfield_is_needed == 1
    save([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_bs_coeff_rot_partial_alpha', 'db_bs_coeff_rot_partial_beta', 'db_bs_coeff_rot_partial_gamma', ...
        'db_s_coeff', 'db_s_coeff_rot_partial_alpha', 'db_s_coeff_rot_partial_beta', 'db_s_coeff_rot_partial_gamma', 'transducer', 'A_delay', 'phi_delay');
else
    save([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff', 'transducer', 'A_delay', 'phi_delay');
end

%%