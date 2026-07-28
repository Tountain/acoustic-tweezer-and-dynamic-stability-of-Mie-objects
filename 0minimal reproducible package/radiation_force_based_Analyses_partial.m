function [Frad_x_partial,Frad_y_partial,Frad_z_partial,transducer] = radiation_force_based_Analyses_partial()
%%
% this function is used to obtain the partial derivatives of Radiation
% Force from the pressure distribution surrounding the particle.
%
% NOTE: the radiation force is given in "radiation_force_based_Analyses.m",
% The partial derivative is based on the expression of radiation force,
% referred to my note POINT 12, PAGE 3, Eq.(VI).  
%%

parameters_names;

% build the database, if already exist, then it will jump to next sentense.
if multi_particle == 0      % for single spherical and axisymmetric object
    if derivative_field_is_needed == 1
        [db_filename] = database_beam_scattering_coeffs();      
        load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_bs_coeff_partial_r', 'db_bs_coeff_partial_theta', 'db_bs_coeff_partial_phi', ...
                                         'db_s_coeff', 'db_s_coeff_partial_r', 'db_s_coeff_partial_theta', 'db_s_coeff_partial_phi', 'transducer', 'A_delay', 'phi_delay');
    elseif derivative_field_is_needed ~= 1
        error('The code is desighed for partial derivatives ONLY!\n');
    end
else
    error('The code is temporally vaild for single object system ONLY!\n');
end

if multi_layer == 1
    error('The code is temporally vaild for single layer system ONLY!\n');
end

%% obtain the force constant "force coeff", "AA", "BB", "CC", "DD", "EE" and "FF"

  
coeff_force = - (1)^2 / (2 * fluid_rho * fluid_c^2);

AA = @(nn,mm) (-1) * sqrt(((nn+mm-1) * (nn+mm)) / ((2*nn-1) * (2*nn+1)));
BB = @(nn,mm) sqrt(((nn-mm+2) * (nn-mm+1)) / ((2*nn+1) * (2*nn+3)));
CC = @(nn,mm) sqrt(((nn-mm-1) * (nn-mm)) / ((2*nn-1) * (2*nn+1)));
DD = @(nn,mm) (-1) * sqrt(((nn+mm+2) * (nn+mm+1)) / ((2*nn+1) * (2*nn+3)));
EE = @(nn,mm) sqrt(((nn-mm) * (nn+mm)) / ((2*nn-1) * (2*nn+1)));
FF = @(nn,mm) sqrt(((nn-mm+1) * (nn+mm+1)) / ((2*nn+1) * (2*nn+3)));

[r_pd, theta_pd, phi_pd] = ...
            coords_system_relative_positions_general([derivativeX, derivativeY, derivativeZ], [0, 0, 0]);
[dr_dx, dr_dy, dr_dz, dtheta_dx, dtheta_dy, dtheta_dz, dphi_dx, dphi_dy, dphi_dz] ...
                    = SphericaltoCartesian_partial(r_pd, theta_pd, phi_pd);

%% obtain the radiation force along X-axis


sum_1_r = 0;
sum_1_theta = 0;
sum_1_phi = 0;
sum_2_r = 0;
sum_2_theta = 0;
sum_2_phi = 0;
sum_3_r = 0;
sum_3_theta = 0;
sum_3_phi = 0;
sum_4_r = 0;
sum_4_theta = 0;
sum_4_phi = 0;
for nn = 0:(N-1)    % (n+1), (m+1)
    for mm = -nn:(nn)
        
        n_plus_1 = nn + 1;      % adjust the projecting relation of database and beam-shape and scattering coefficents
        m_plus_1 = mm + 1;
        if abs(m_plus_1) <= n_plus_1    % make sure the adjusting position within the database
            if irregular_body == 0
                % left blank
            elseif irregular_body == 1
                sum_1_r = sum_1_r + ...
                    AA(nn+1,mm+1) * radiation_effect_unit_partial(nn, mm, n_plus_1, m_plus_1, db_bs_coeff, db_bs_coeff_partial_r, db_s_coeff, db_s_coeff_partial_r);
                sum_1_theta = sum_1_theta + ...
                    AA(nn+1,mm+1) * radiation_effect_unit_partial(nn, mm, n_plus_1, m_plus_1, db_bs_coeff, db_bs_coeff_partial_theta, db_s_coeff, db_s_coeff_partial_theta);
                sum_1_phi = sum_1_phi + ...
                    AA(nn+1,mm+1) * radiation_effect_unit_partial(nn, mm, n_plus_1, m_plus_1, db_bs_coeff, db_bs_coeff_partial_phi, db_s_coeff, db_s_coeff_partial_phi);
            end
        end
    
    end
end
for nn = 1:(N)      % (n-1), (m+1)
    for mm = -nn:(nn)
        
        n_minus_1 = nn - 1;
        m_plus_1 = mm + 1;
        if abs(m_plus_1) <= n_minus_1
            if irregular_body == 0
                % left blank
            elseif irregular_body == 1
                sum_2_r = sum_2_r + ...
                    BB(nn-1,mm+1) * radiation_effect_unit_partial(nn, mm, n_minus_1, m_plus_1, db_bs_coeff, db_bs_coeff_partial_r, db_s_coeff, db_s_coeff_partial_r);
                sum_2_theta = sum_2_theta + ...
                    BB(nn-1,mm+1) * radiation_effect_unit_partial(nn, mm, n_minus_1, m_plus_1, db_bs_coeff, db_bs_coeff_partial_theta, db_s_coeff, db_s_coeff_partial_theta);
                sum_2_phi = sum_2_phi + ...
                    BB(nn-1,mm+1) * radiation_effect_unit_partial(nn, mm, n_minus_1, m_plus_1, db_bs_coeff, db_bs_coeff_partial_phi, db_s_coeff, db_s_coeff_partial_phi);
            end
        end
    
    end
end
for nn = 0:(N-1)    % (n+1), (m-1)
    for mm = (-nn):(nn)
       
        n_plus_1 = nn + 1;
        m_minus_1 = mm - 1;
        if abs(m_minus_1) <= n_plus_1
            if irregular_body == 0
                % left blank
            elseif irregular_body == 1
                sum_3_r = sum_3_r + ...
                    CC(nn+1,mm-1) * radiation_effect_unit_partial(nn, mm, n_plus_1, m_minus_1, db_bs_coeff, db_bs_coeff_partial_r, db_s_coeff, db_s_coeff_partial_r);
                sum_3_theta = sum_3_theta + ...
                    CC(nn+1,mm-1) * radiation_effect_unit_partial(nn, mm, n_plus_1, m_minus_1, db_bs_coeff, db_bs_coeff_partial_theta, db_s_coeff, db_s_coeff_partial_theta);
                sum_3_phi = sum_3_phi + ...
                    CC(nn+1,mm-1) * radiation_effect_unit_partial(nn, mm, n_plus_1, m_minus_1, db_bs_coeff, db_bs_coeff_partial_phi, db_s_coeff, db_s_coeff_partial_phi);
            end
        end
    
    end
end
for nn = 1:(N)      % (n-1), (m-1)
    for mm = (-nn):(nn)
        
        n_minus_1 = nn - 1;
        m_minus_1 = mm - 1;
        if abs(m_minus_1) <= n_minus_1
            if irregular_body == 0
                % left blank
            elseif irregular_body == 1
                sum_4_r = sum_4_r + ...
                    DD(nn-1,mm-1) * radiation_effect_unit_partial(nn, mm, n_minus_1, m_minus_1, db_bs_coeff, db_bs_coeff_partial_r, db_s_coeff, db_s_coeff_partial_r);
                sum_4_theta = sum_4_theta + ...
                    DD(nn-1,mm-1) * radiation_effect_unit_partial(nn, mm, n_minus_1, m_minus_1, db_bs_coeff, db_bs_coeff_partial_theta, db_s_coeff, db_s_coeff_partial_theta);
                sum_4_phi = sum_4_phi + ...
                    DD(nn-1,mm-1) * radiation_effect_unit_partial(nn, mm, n_minus_1, m_minus_1, db_bs_coeff, db_bs_coeff_partial_phi, db_s_coeff, db_s_coeff_partial_phi);
            end
        end
    
    end
end

 
Frad_x_partial_r = coeff_force * 1/2 * real(1i * (sum_1_r - sum_2_r + sum_3_r - sum_4_r) / (fluid_k^2));     
Frad_x_partial_theta = coeff_force * 1/2 * real(1i * (sum_1_theta - sum_2_theta + sum_3_theta - sum_4_theta) / (fluid_k^2));  
Frad_x_partial_phi = coeff_force * 1/2 * real(1i * (sum_1_phi - sum_2_phi + sum_3_phi - sum_4_phi) / (fluid_k^2));  

Frad_x_partial = Frad_x_partial_r * dr_dx + Frad_x_partial_theta * dtheta_dx + Frad_x_partial_phi * dphi_dx;


%% obtain the radiation force along Y-axis


Frad_y_partial_r = coeff_force * 1/2 * real((sum_1_r - sum_2_r - sum_3_r + sum_4_r) / (fluid_k^2)); 
Frad_y_partial_theta = coeff_force * 1/2 * real((sum_1_theta - sum_2_theta - sum_3_theta + sum_4_theta) / (fluid_k^2)); 
Frad_y_partial_phi = coeff_force * 1/2 * real((sum_1_phi - sum_2_phi - sum_3_phi + sum_4_phi) / (fluid_k^2)); 

Frad_y_partial = Frad_y_partial_r * dr_dy + Frad_y_partial_theta * dtheta_dy + Frad_y_partial_phi * dphi_dy;


%% obtain the radiation force along Z-axis (parallel with incident wave)

sum_1_r = 0;
sum_1_theta = 0;
sum_1_phi = 0;
sum_2_r = 0;
sum_2_theta = 0;
sum_2_phi = 0;
for nn = 0:(N-1)     % (n+1), m
    for mm = -nn:(nn)
        
        n_plus_1 = nn + 1;
        if abs(mm) <= n_plus_1
            if irregular_body == 0
                % left blank
            elseif irregular_body == 1
                sum_1_r = sum_1_r + ...
                    EE(nn+1,mm) * radiation_effect_unit_partial(nn, mm, n_plus_1, mm, db_bs_coeff, db_bs_coeff_partial_r, db_s_coeff, db_s_coeff_partial_r);
                sum_1_theta = sum_1_theta + ...
                    EE(nn+1,mm) * radiation_effect_unit_partial(nn, mm, n_plus_1, mm, db_bs_coeff, db_bs_coeff_partial_theta, db_s_coeff, db_s_coeff_partial_theta);
                sum_1_phi = sum_1_phi + ...
                    EE(nn+1,mm) * radiation_effect_unit_partial(nn, mm, n_plus_1, mm, db_bs_coeff, db_bs_coeff_partial_phi, db_s_coeff, db_s_coeff_partial_phi);
            end
        end
        
    end
end
for nn = 1:(N)      % (n-1), m
    for mm = -(nn):(nn)
        
        n_minus_1 = nn - 1;
        if abs(mm) <= n_minus_1
            if irregular_body == 0
                % left blank
            elseif irregular_body == 1
                sum_2_r = sum_2_r + ...
                    FF(nn-1,mm) * radiation_effect_unit_partial(nn, mm, n_minus_1, mm, db_bs_coeff, db_bs_coeff_partial_r, db_s_coeff, db_s_coeff_partial_r);
                sum_2_theta = sum_2_theta + ...
                    FF(nn-1,mm) * radiation_effect_unit_partial(nn, mm, n_minus_1, mm, db_bs_coeff, db_bs_coeff_partial_theta, db_s_coeff, db_s_coeff_partial_theta);
                sum_2_phi = sum_2_phi + ...
                    FF(nn-1,mm) * radiation_effect_unit_partial(nn, mm, n_minus_1, m_plus_1, db_bs_coeff, db_bs_coeff_partial_phi, db_s_coeff, db_s_coeff_partial_phi);
            end
        end
        
    end
end


Frad_z_partial_r = coeff_force * real(1i * (sum_1_r - sum_2_r) / (fluid_k^2));      
Frad_z_partial_theta = coeff_force * real(1i * (sum_1_theta - sum_2_theta) / (fluid_k^2));    
Frad_z_partial_phi = coeff_force * real(1i * (sum_1_phi - sum_2_phi) / (fluid_k^2));   

Frad_z_partial = Frad_z_partial_r * dr_dz + Frad_z_partial_theta * dtheta_dz + Frad_z_partial_phi * dphi_dz;


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
    %     Frad_anti = Rzyx_anti * [Frad_x; Frad_y; Frad_z];
        Frad_anti_partial = [Frad_x_partial, Frad_y_partial, Frad_z_partial] * Rzyx;
        Frad_x_partial = Frad_anti_partial(1);
        Frad_y_partial = Frad_anti_partial(2); 
        Frad_z_partial = Frad_anti_partial(3);
    end
end


%% saving file  ('A.mat' for analyses results)

E0 = p_inlet^2 / (2 * fluid_rho * fluid_c^2);
pi_a2_E0 = (pi*particle_radius^2) * E0;
% Y_x = Frad_x / pi_a2_E0;        % dimensionaless radiation force function x
% Y_y = Frad_y / pi_a2_E0;        % dimensionaless radiation force function y
% Y_z = Frad_z / pi_a2_E0;        % dimensionaless radiation force function z
 
if (strcmp(wave_type, 'single_transducer') || strcmp(wave_type, 'phase_array_transducer') || strcmp(wave_type, 'phase_array_transducer2')) ~= 1
    if multi_particle == 0
%         if multi_layer == 0
            save([dir_file_forces_partial, 'A.mat'], 'E0', 'pi_a2_E0', 'Frad_x_partial', 'Frad_y_partial', 'Frad_z_partial');
%         elseif multi_layer == 1
%             save([dir_file_forces, 'A.mat'], 'E0', 'pi_a2_E0', 'Frad_x', 'Frad_y', 'Frad_z', 'c_layer', 'density_layer', 'depth_fun');
%         end
%     else
%         particles_Cartesian_data;
%         if irregular_body == 1
%             dir_file_forces_multi = [dir_file_forces_multi, '_Cn_multi'];
%             for ii =1:particle_number
%                 dir_file_forces_multi = [dir_file_forces_multi, num2str(length(Cn_multi{ii}))];
%             end
%         end
%         dir_file_forces_multi = [dir_file_forces_multi, '_Coords('];
%         particle = particle / particle_radius;
%         for ii =1:particle_number
%             particle_code = [num2str(particle(ii,1)), num2str(particle(ii,2)), num2str(particle(ii,3))];
%             dir_file_forces_multi = [dir_file_forces_multi, particle_code, ','];
%         end
%         dir_file_forces_multi(end) = ')';
%         if multi_layer == 0
%             save([dir_file_forces_multi, 'A.mat'], 'E0', 'pi_a2_E0', ...
%                 'Frad_x', 'Frad_y', 'Frad_z', 'particle', 'multi_particle_radius');
%         elseif multi_layer == 1
%             save([dir_file_forces_multi, 'A.mat'], 'E0', 'pi_a2_E0', ...
%                 'Frad_x', 'Frad_y', 'Frad_z', 'particle', 'multi_particle_radius', 'c_layer', 'density_layer' ,'depth_fun');
%         end
    end
else        % saving for 'single_transducer' and 'phase_array_transducer' cases
    if multi_particle == 0
%         if multi_layer == 0
            save([dir_file_forces_trans_partial, 'A.mat'], 'E0', 'pi_a2_E0', 'transducer', 'A_delay', 'phi_delay', 'Frad_x_partial', 'Frad_y_partial', 'Frad_z_partial');
%         elseif multi_layer == 1
%             save([dir_file_forces_trans, 'A.mat'], 'E0', 'pi_a2_E0', 'transducer', 'Frad_x', 'Frad_y', 'Frad_z', 'c_layer', 'density_layer' , 'depth_fun');
%         end
%     else
%         particles_Cartesian_data;
%         dir_file_forces_trans = [dir_file_forces_trans, '_Coords('];
%         particle = particle / particle_radius;
%         for ii =1:particle_number
%             particle_code = [num2str(particle(ii,1)), num2str(particle(ii,2)), num2str(particle(ii,3))];
%             dir_file_forces_trans = [dir_file_forces_trans, particle_code, ','];
%         end
%         dir_file_forces_trans(end) = ')';
%         if multi_layer == 0
%             save([dir_file_forces_trans, 'A.mat'], 'E0', 'pi_a2_E0', 'transducer', ...
%                 'Frad_x', 'Frad_y', 'Frad_z', 'particle', 'multi_particle_radius');
%         elseif multi_layer == 1
%             save([dir_file_forces_trans, 'A.mat'], 'E0', 'pi_a2_E0', 'transducer', ...
%                 'Frad_x', 'Frad_y', 'Frad_z', 'particle', 'multi_particle_radius', 'c_layer', 'density_layer' , 'depth_fun');
%         end
    end
end


%% clear database for beam and scattering coefficient

if multi_particle == 0
%     delete([db_filename, '.mat']);
% else
%     delete([db_filename, '.mat']);
%     delete([db_filename_ts, '.mat']);
%     delete([db_filename_ibs, '.mat']);
%     delete([db_filename_str, '.mat']);
end


%%
