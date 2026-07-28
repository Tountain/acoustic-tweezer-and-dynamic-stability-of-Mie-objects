function [Torque_x_partial,Torque_y_partial,Torque_z_partial,transducer] = radiation_torque_based_Analyses_partial()
%%
% this function is used to obtain the Radiation Torque from the pressure
% distribution surrounding the particle.
%
%
% NOTE: This function IS ABLE TO work for conditions of incident wave 
% along different direction. (Therefore, SUITS for multi-particle system)
%%

parameters_names;

% build the database, if already exist, then it will jump to next sentense.
if multi_particle == 0      % for single spherical and axisymmetric object
    if derivative_rotationfield_is_needed == 1
        [db_filename] = database_beam_scattering_coeffs();      
        load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_bs_coeff_rot_partial_alpha', 'db_bs_coeff_rot_partial_beta', 'db_bs_coeff_rot_partial_gamma', ...
                                         'db_s_coeff', 'db_s_coeff_rot_partial_alpha', 'db_s_coeff_rot_partial_beta', 'db_s_coeff_rot_partial_gamma', 'transducer');
    elseif derivative_rotationfield_is_needed ~= 1
        error('The code is desighed for partial derivatives (rotation) ONLY!\n');
    end
else
    error('The code is temporally vaild for single object system ONLY!\n');
end

if multi_layer == 1
    error('The code is temporally vaild for single layer system ONLY!\n');
end

 
%% obtain the torque constant "torque coeff", "BB_p" and "BB_n"


coeff_torque = - (1)^2 / (2 * fluid_rho * fluid_c^2);

BB_p = @(nn,mm) sqrt((nn-mm) * (nn+mm+1));
BB_n = @(nn,mm) sqrt((nn+mm) * (nn-mm+1));

% [r_pd, theta_pd, phi_pd] = ...
%             coords_system_relative_positions_general([derivativeX, derivativeY, derivativeZ], [0, 0, 0]);
% [dr_dx, dr_dy, dr_dz, dtheta_dx, dtheta_dy, dtheta_dz, dphi_dx, dphi_dy, dphi_dz] ...
%                     = SphericaltoCartesian_partial(r_pd, theta_pd, phi_pd);
dalpha_dthetax = -cos(derivative_thetaX) / tan(derivative_thetaY);
dalpha_dthetay = -sin(derivative_thetaX) / tan(derivative_thetaY);
dalpha_dthetaz = 1;
dbeta_dthetax = -sin(derivative_thetaX);
dbeta_dthetay = cos(derivative_thetaX);
dbeta_dthetaz = 0;
dgamma_dthetax = cos(derivative_thetaX) / sin(derivative_thetaY);
dgamma_dthetay = sin(derivative_thetaX) / sin(derivative_thetaY);
dgamma_dthetaz = 0;

%% obtain the radiation torque rotating with X-axis

sum_1_alpha = 0;
sum_1_beta = 0;
sum_1_gamma = 0;
sum_2_alpha = 0;
sum_2_beta = 0;
sum_2_gamma = 0;
for nn = 0:(N)    % (n), (m-1)
    for mm = -nn:(nn)
       
        m_minus_1 = mm - 1;
        if abs(m_minus_1) <= nn    % make sure the adjusting position within the database
            if irregular_body == 0
                % left blank
            elseif irregular_body == 1
                sum_1_alpha = sum_1_alpha + ...
                    BB_n(nn+1,mm+1) * radiation_effect_unit_partial(nn, mm, nn, m_minus_1, db_bs_coeff, db_bs_coeff_rot_partial_alpha, db_s_coeff, db_s_coeff_rot_partial_alpha);
                sum_1_beta = sum_1_beta + ...
                    BB_n(nn+1,mm+1) * radiation_effect_unit_partial(nn, mm, nn, m_minus_1, db_bs_coeff, db_bs_coeff_rot_partial_beta, db_s_coeff, db_s_coeff_rot_partial_beta);
                sum_1_gamma = sum_1_gamma + ...
                    BB_n(nn+1,mm+1) * radiation_effect_unit_partial(nn, mm, nn, m_minus_1, db_bs_coeff, db_bs_coeff_rot_partial_gamma, db_s_coeff, db_s_coeff_rot_partial_gamma);
            end
        end
        
    end
end

for nn = 0:(N)      % (n), (m+1)
    for mm = -nn:(nn)
        
        m_plus_1 = mm + 1;
        if abs(m_plus_1) <= nn    % make sure the adjusting position within the database
            if irregular_body == 0
                % left blank
            elseif irregular_body == 1
                sum_2_alpha = sum_2_alpha + ...
                    BB_p(nn+1,mm+1) * radiation_effect_unit_partial(nn, mm, nn, m_plus_1, db_bs_coeff, db_bs_coeff_rot_partial_alpha, db_s_coeff, db_s_coeff_rot_partial_alpha);
                sum_2_beta = sum_2_beta + ...
                    BB_p(nn+1,mm+1) * radiation_effect_unit_partial(nn, mm, nn, m_plus_1, db_bs_coeff, db_bs_coeff_rot_partial_beta, db_s_coeff, db_s_coeff_rot_partial_beta);
                sum_2_gamma = sum_2_gamma + ...
                    BB_p(nn+1,mm+1) * radiation_effect_unit_partial(nn, mm, nn, m_plus_1, db_bs_coeff, db_bs_coeff_rot_partial_gamma, db_s_coeff, db_s_coeff_rot_partial_gamma);
            end
        end
        
    end
end


Torque_x_partial_alpha = coeff_torque / (fluid_k^3) * 1/2 * real(sum_2_alpha + sum_1_alpha);
Torque_x_partial_beta = coeff_torque / (fluid_k^3) * 1/2 * real(sum_2_beta + sum_1_beta);
Torque_x_partial_gamma = coeff_torque / (fluid_k^3) * 1/2 * real(sum_2_gamma + sum_1_gamma);

Torque_x_partial = Torque_x_partial_alpha * dalpha_dthetax + Torque_x_partial_beta * dbeta_dthetax + Torque_x_partial_gamma * dgamma_dthetax;

    
%% obtain the radiation torque rotating with Y-axis


Torque_y_partial_alpha = coeff_torque / (fluid_k^3) * 1/2 * real(sum_2_alpha - sum_1_alpha);
Torque_y_partial_beta = coeff_torque / (fluid_k^3) * 1/2 * real(sum_2_beta - sum_1_beta);
Torque_y_partial_gamma = coeff_torque / (fluid_k^3) * 1/2 * real(sum_2_gamma - sum_1_gamma);

Torque_y_partial = Torque_y_partial_alpha * dalpha_dthetay + Torque_y_partial_beta * dbeta_dthetay + Torque_y_partial_gamma * dgamma_dthetay;


%% obtain the radiation torque rotating with Z-axis (parallel with incident wave)

sum_1_alpha = 0;
sum_1_beta = 0;
sum_1_gamma = 0;
for nn = 0:(N)     % (n), m
    for mm = -nn:(nn)

        if irregular_body == 0
            % left blank
        elseif irregular_body == 1
            sum_1_alpha = sum_1_alpha + ...
                mm * radiation_effect_unit_partial(nn, mm, nn, mm, db_bs_coeff, db_bs_coeff_rot_partial_alpha, db_s_coeff, db_s_coeff_rot_partial_alpha);
            sum_1_beta = sum_1_beta + ...
                mm * radiation_effect_unit_partial(nn, mm, nn, mm, db_bs_coeff, db_bs_coeff_rot_partial_beta, db_s_coeff, db_s_coeff_rot_partial_beta);
            sum_1_gamma = sum_1_gamma + ...
                mm * radiation_effect_unit_partial(nn, mm, nn, mm, db_bs_coeff, db_bs_coeff_rot_partial_gamma, db_s_coeff, db_s_coeff_rot_partial_gamma);
        end
        
    end
end


Torque_z_partial_alpha = coeff_torque / (fluid_k^3) * real(sum_1_alpha);
Torque_z_partial_beta = coeff_torque / (fluid_k^3) * real(sum_1_beta);
Torque_z_partial_gamma = coeff_torque / (fluid_k^3) * real(sum_1_gamma);

Torque_z_partial = Torque_z_partial_alpha * dalpha_dthetaz + Torque_z_partial_beta * dbeta_dthetaz + Torque_z_partial_gamma * dgamma_dthetaz;


%% anti-rotation of the particle coordinate system for irregular particles

if irregular_body == 1
    if rot_matrix ~= 1
        Rx=inv([1 0 0;
                0 cos(theta_rotation(1)) -sin(theta_rotation(1)); 
                0 sin(theta_rotation(1)) cos(theta_rotation(1))]);
        Ry=inv([cos(theta_rotation(2)) 0 sin(theta_rotation(2));
                0 1 0; 
                -sin(theta_rotation(2)) 0 cos(theta_rotation(2))]);
        Rz=inv([cos(theta_rotation(3)) -sin(theta_rotation(3)) 0; 
                sin(theta_rotation(3)) cos(theta_rotation(3)) 0;
                0 0 1]);
        Rzyx = Rz * Ry * Rx;
    %     Torque_anti = Rzyx_anti * [Torque_x; Torque_y; Torque_z];
        Torque_anti_partial = [Torque_x_partial, Torque_y_partial, Torque_z_partial] * Rzyx;
        Torque_x_partial = Torque_anti_partial(1);
        Torque_y_partial = Torque_anti_partial(2);
        Torque_z_partial = Torque_anti_partial(3);
    end
end


%% saving file  ('A.mat' for analyses results)

E0 = p_inlet^2 / (2 * fluid_rho * fluid_c^2);   % J/m3
pi_a3_E0 = (pi*particle_radius^3) * E0;     % J or N*m
% N_x = Torque_x / pi_a3_E0;        % dimensionaless radiation torque function x
% N_y = Torque_y / pi_a3_E0;        % dimensionaless radiation torque function y
% N_z = Torque_z / pi_a3_E0;        % dimensionaless radiation torque function z

% save([dir_file_torques, 'A.mat'], 'E0', 'pi_a3_E0', 'Torque_x', 'Torque_y', 'Torque_z');
if (strcmp(wave_type, 'single_transducer') || strcmp(wave_type, 'phase_array_transducer') || strcmp(wave_type, 'phase_array_transducer2')) ~= 1
    if multi_particle == 0
%         if multi_layer == 0
            save([dir_file_torques_partial, 'A.mat'], 'E0', 'pi_a3_E0', 'Torque_x_partial', 'Torque_y_partial', 'Torque_z_partial');
%         elseif multi_layer == 1
%             save([dir_file_torques, 'A.mat'], 'E0', 'pi_a3_E0', 'Torque_x', 'Torque_y', 'Torque_z', 'c_layer', 'density_layer' , 'depth_fun');
%         end
%     else
%         particles_Cartesian_data;
%         if irregular_body == 1
%             dir_file_torques_multi = [dir_file_torques_multi, '_Cn_multi'];
%             for ii =1:particle_number
%                 dir_file_torques_multi = [dir_file_torques_multi, num2str(length(Cn_multi{ii}))];
%             end
%         end
%         dir_file_torques_multi = [dir_file_torques_multi, '_Coords('];
%         particle = particle / particle_radius;
%         for ii =1:particle_number
%             particle_code = [num2str(particle(ii,1)), num2str(particle(ii,2)), num2str(particle(ii,3))];
%             dir_file_torques_multi = [dir_file_torques_multi, particle_code, ','];
%         end
%         dir_file_torques_multi(end) = ')';
%         if multi_layer == 0
%             save([dir_file_torques_multi, 'A.mat'], 'E0', 'pi_a3_E0', ...
%                 'Torque_x', 'Torque_y', 'Torque_z', 'particle', 'multi_particle_radius');
%         elseif multi_layer == 1
%             save([dir_file_torques_multi, 'A.mat'], 'E0', 'pi_a3_E0', ...
%                 'Torque_x', 'Torque_y', 'Torque_z', 'particle', 'multi_particle_radius', 'c_layer', 'density_layer' , 'depth_fun');
%         end
    end
else        % saving for 'single_transducer' and 'phase_array_transducer' cases
    if multi_particle == 0
%         if multi_layer == 0
            save([dir_file_torques_trans_partial, 'A.mat'], 'E0', 'pi_a3_E0', 'transducer', 'Torque_x_partial', 'Torque_y_partial', 'Torque_z_partial');
%         elseif multi_layer == 1
%             save([dir_file_torques_trans, 'A.mat'], 'E0', 'pi_a3_E0', 'transducer', 'Torque_x', 'Torque_y', 'Torque_z', 'c_layer', 'density_layer' , 'depth_fun');
%         end
%     else
%         particles_Cartesian_data;
%         dir_file_torques_trans = [dir_file_torques_trans, '_Coords('];
%         particle = particle / particle_radius;
%         for ii =1:particle_number
%             particle_code = [num2str(particle(ii,1)), num2str(particle(ii,2)), num2str(particle(ii,3))];
%             dir_file_torques_trans = [dir_file_torques_trans, particle_code, ','];
%         end
%         dir_file_torques_trans(end) = ')';
%         if multi_layer == 0
%             save([dir_file_torques_trans, 'A.mat'], 'E0', 'pi_a3_E0', 'transducer', ...
%                 'Torque_x', 'Torque_y', 'Torque_z', 'particle', 'multi_particle_radius');
%         elseif multi_layer == 1
%             save([dir_file_torques_trans, 'A.mat'], 'E0', 'pi_a3_E0', 'transducer', ...
%                 'Torque_x', 'Torque_y', 'Torque_z', 'particle', 'multi_particle_radius', 'c_layer', 'density_layer' , 'depth_fun');
%         end
    end
end

%% clear database for beam and scattering coefficient

if multi_particle == 0
    delete([db_filename, '.mat']);
% else
%     delete([db_filename, '.mat']);
%     delete([db_filename_ts, '.mat']);
%     delete([db_filename_ibs, '.mat']);
%     delete([db_filename_str, '.mat']);
end


%%

