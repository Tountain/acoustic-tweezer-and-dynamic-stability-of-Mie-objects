function [delay_amplitude, delay_phase, F] = design_transducer_parameters(A_delay, phi_delay, coeff_force, coeff_torque, Fs, Ts, dFs, mode_type, amp_fix)
%%
% This function is used to design the transducer parameter (phase and/or
% amplitude) of a transducer array. In this way, we need to specify the
% desired motions of a levitated particle, then we retrieval the phase
% and/or amplitude distributions.
%
% initial value of the transducer parameters: 'A_delay' and 'phi_delay'.
%%
 
parameters_names;
Nt = transducer_number;

if strcmp(wave_type, 'phase_array_transducer') ~= 1
    error('This function is designed for transducer array!\n');
end
if irregular_body == 0
    error('"irregular_body == 0" can be replaced with "irregular_body == 1" and "Cn = [a, 0]"!\n');
end

%% build the database for equivalent "beam_shape_coeffcient" 
db_bs_coeff = zeros(db_size_nn+1, 2*db_size_nn+1);
if multi_particle == 0      % single particle system
    % single nonspherical object
    db_s_coeff = zeros(db_size_nn+1, 2*db_size_nn+1);
end

N = db_size_nn; 
% database for "beam_shape_coeffcient": db_bs_coeff
for nn = 0:N
    for mm = -nn:nn
        
        db_bs_coeff(nn+1, mm+nn+1) = beam_shape_coeff(nn, mm);  % Eq. (6) of my manuscript for engineering the force and torque

    end
    fprintf('Beam-Shape and Scattering Coefficients Database Preparing %d%% \n',round(100*nn/N));
end


% the equivalent beam-shape coeffcients, that involve the contributions
% from the transducers ('A_delay' and 'phase_delay')
% equivalent_bs_q_coeff = db_bs_coeff;
[equivalent_bs_q_coeff] = ...
    transducer_transform_beam_shape_coeff(N, db_bs_coeff, fluid_k, transducer, transducer_number, theta_rotation, 2, A_delay, phi_delay);

% database for "scalar scattering_coefficient": db_s_coeff
% Note that for non-spherical object, the scalr scattering coefficients are
% dependent to the transducer parameters through 'equivalent_bs_q_coeff'.
if multi_particle == 0      % single particle system
    % single nonspherical object
    snm_ib = scattering_coeff_irregular_body(BC, equivalent_bs_q_coeff);
    for nn = 0:db_size_nn
        for mm = -nn:nn
            db_s_coeff(nn+1, mm+nn+1) = snm_ib(nn+1, mm+nn+1);
        end
    end
end

%% build the characteristic database; if already exist, then it will jump to load the data.
[db_filename_CharMat] = database_characteristic_matrix(N, db_bs_coeff, db_s_coeff);      % Pls make sure the 'wave_type' is 'phase_array_transducer' or 'phase_array_transducer2'
load([db_filename_CharMat, '.mat'], 'trans_bs_coeff', 'db_M_Fx', 'db_N_Fx', 'db_M_Fy', 'db_N_Fy', 'db_M_Fz', 'db_N_Fz', ...
                                'db_M_Tx', 'db_N_Tx', 'db_M_Ty', 'db_N_Ty', 'db_M_Tz', 'db_N_Tz');

% [db_filename_EQ_VIII] = database_EQ_VIII(N, db_bs_coeff, db_s_coeff);      % Pls make sure the 'wave_type' is 'phase_array_transducer' or 'phase_array_transducer2'
% load([db_filename_EQ_VIII, '.mat'], 'db_AA', 'db_BB', 'db_CC', 'db_DD', 'db_EE', 'db_FF', ...
%                                 'db_GGn', 'db_GGp', 'db_GG');
%                             
% [db_filename_EQ_X] = database_EQ_X(N, db_bs_coeff, db_s_coeff);      % Pls make sure the 'wave_type' is 'phase_array_transducer' or 'phase_array_transducer2'
% load([db_filename_EQ_X, '.mat'], 'db_AA_X', 'db_BB_X', 'db_CC_X', 'db_DD_X', 'db_EE_X', 'db_FF_X', ...
%                                 'db_GGn_X', 'db_GGp_X', 'db_GG_X');
%                             
% [db_filename_EQ_XII] = database_EQ_XII(N, db_bs_coeff, db_s_coeff, phi_delay);      % Pls make sure the 'wave_type' is 'phase_array_transducer' or 'phase_array_transducer2'
% load([db_filename_EQ_XII, '.mat'], 'db_AA_XII', 'db_BB_XII', 'db_CC_XII', 'db_DD_XII', 'db_EE_XII', 'db_FF_XII', ...
%                                 'db_GGn_XII', 'db_GGp_XII', 'db_GG_XII');
                   

%% scalar scattering_coefficient and partial derivatives of scalar scattering coefficients

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
                    = beam_shape_coeff_partial(nn, mm, equivalent_bs_q_coeff, fluid_k, r_pd, theta_pd, phi_pd);
    end 
    fprintf('Translated Equivalent Beam-Shape Coefficient and its partial derivatives (translation) %d%% \n', round(100*nn/N));
end  


% database for "scalar scattering_coefficient": db_s_coeff
% Note that for non-spherical object, the scalr scattering coefficients are
% dependent to the transducer parameters through 'equivalent_bs_q_coeff'.
if multi_particle == 0      % single particle system
    % single nonspherical object
    snm_ib = scattering_coeff_irregular_body(BC, db_bs_coeff_translation);
    for nn = 0:db_size_nn
        for mm = -nn:nn
            db_s_coeff(nn+1, mm+nn+1) = snm_ib(nn+1, mm+nn+1);
        end
    end
end

[db_s_coeff_partial_r, db_s_coeff_partial_theta, db_s_coeff_partial_phi] ...
       = scattering_coeff_irregular_body_partial(BC, db_bs_coeff_translation, db_bs_coeff_partial_r, db_bs_coeff_partial_theta, db_bs_coeff_partial_phi, db_s_coeff);


%% partial derivatives of transducer beam-shape coefficients

trans_bs_coeff_translation = cell(1,Nt);
trans_bs_coeff_partial_r = cell(1,Nt);
trans_bs_coeff_partial_theta = cell(1,Nt);
trans_bs_coeff_partial_phi = cell(1,Nt);

[r_pd, theta_pd, phi_pd] = ...
            coords_system_relative_positions_general([derivativeX, derivativeY, derivativeZ], [0, 0, 0]); % '- \vec{r}'
% establish database of partial derivatives of seperation matrix of first
% kind: avoiding repeatly calculation of 'Snmvu_partial'.
[db_filename_ts_partial] = database_tranlation_coeff_partial(N, fluid_k, r_pd, theta_pd, phi_pd);
load([db_filename_ts_partial, '.mat'], 'Snmvu_1', 'Snmvu_partial_r', 'Snmvu_partial_theta', 'Snmvu_partial_phi');

for ii = 1:Nt
    db_bs_coeff_translation = zeros(db_size_nn+1, 2*db_size_nn+1);
    db_bs_coeff_partial_r = zeros(db_size_nn+1, 2*db_size_nn+1);
    db_bs_coeff_partial_theta = zeros(db_size_nn+1, 2*db_size_nn+1);
    db_bs_coeff_partial_phi = zeros(db_size_nn+1, 2*db_size_nn+1);
    % initial for parfor
    fluid_k = fluid_k;
    trans_bs_coeff = trans_bs_coeff;
    Snmvu_1 = Snmvu_1;
    Snmvu_partial_r = Snmvu_partial_r;
    Snmvu_partial_theta = Snmvu_partial_theta;
    Snmvu_partial_phi = Snmvu_partial_phi;
    
    for nn = 0:N
        indices_1 = nn + 1;
        parfor mm = -nn:nn          % 'mm' is sliced variable
%             [db_bs_coeff_translation(indices_1, indices_1+mm), db_bs_coeff_partial_r(indices_1, indices_1+mm), db_bs_coeff_partial_theta(indices_1, indices_1+mm), db_bs_coeff_partial_phi(indices_1, indices_1+mm)] ...
%                         = beam_shape_coeff_partial(nn, mm, trans_bs_coeff{ii}, fluid_k, r_pd, theta_pd, phi_pd);
            [db_bs_coeff_translation(indices_1, indices_1+mm), db_bs_coeff_partial_r(indices_1, indices_1+mm), db_bs_coeff_partial_theta(indices_1, indices_1+mm), db_bs_coeff_partial_phi(indices_1, indices_1+mm)] ...
                        = beam_shape_coeff_partial_using_database(nn, mm, trans_bs_coeff{ii}, Snmvu_1, Snmvu_partial_r, Snmvu_partial_theta, Snmvu_partial_phi);
        end 
        fprintf('Translated Beam-Shape Coefficient of %d-th transducers and its partial derivatives (translation) %d%% \n', ii, round(100*nn/N));
    end  
    
    trans_bs_coeff_translation{ii} = db_bs_coeff_translation;
    trans_bs_coeff_partial_r{ii} = db_bs_coeff_partial_r;
    trans_bs_coeff_partial_theta{ii} = db_bs_coeff_partial_theta;
    trans_bs_coeff_partial_phi{ii} = db_bs_coeff_partial_phi;
    
    %trans_bs_coeff{ii} = trans_bs_coeff_translation{ii};
end


%% build the partial characteristic database; if already exist, then it will jump to load the data.

[db_filename_CharMat_partial] = database_characteristic_matrix_partial(N, ...
                trans_bs_coeff_translation, trans_bs_coeff_partial_r, trans_bs_coeff_partial_theta, trans_bs_coeff_partial_phi, ...
                db_s_coeff, db_s_coeff_partial_r, db_s_coeff_partial_theta, db_s_coeff_partial_phi);
load([db_filename_CharMat_partial, '.mat'], 'db_M_Fx_r', 'db_N_Fx_r', 'db_M_Fy_r', 'db_N_Fy_r', 'db_M_Fz_r', 'db_N_Fz_r', ...
                                    'db_M_Fx_theta', 'db_N_Fx_theta', 'db_M_Fy_theta', 'db_N_Fy_theta', 'db_M_Fz_theta', 'db_N_Fz_theta', ...
                                    'db_M_Fx_phi', 'db_N_Fx_phi', 'db_M_Fy_phi', 'db_N_Fy_phi', 'db_M_Fz_phi', 'db_N_Fz_phi');


%% solving the non-linear equations for the transducer parameters

mode_system = mode_type;        % for the governing system

% x0 = [ones(1,Nt) zeros(1,Nt)];
Amp_pha = [A_delay, phi_delay];
if mode_system == 1
    x0 = Amp_pha;
elseif mode_system == 2
    x0 = Amp_pha(1:Nt);
elseif mode_system == 3
%     x0 = [Amp_pha(1+Nt:2*Nt) Amp_pha(1)];
    x0 = [Amp_pha(1+Nt:2*Nt) amp_fix];      % 'amp_fix' is also the optimization parameter
elseif mode_system == 4
    x0 = [Amp_pha(1+Nt:2*Nt) amp_fix];      % 'amp_fix' is not the optimization parameter
end
% check_consistency();
F = evalEqs(x0, theta_rotation, Nt, db_M_Fx, db_N_Fx, db_M_Fy, db_N_Fy, db_M_Fz, db_N_Fz, db_M_Tx, db_N_Tx, db_M_Ty, db_N_Ty, db_M_Tz, db_N_Tz, coeff_force, coeff_torque, Fs, Ts, mode_system, x0(end));
F
% F_EQ_VIII = evalEqs_EQ_VIII(x0, theta_rotation, Nt, db_AA, db_BB, db_CC, db_DD, db_EE, db_FF, db_GGn, db_GGp, db_GG, coeff_force, coeff_torque, Fs, Ts, mode_system, x0(end));
% F_EQ_VIII
% F_EQ_X = evalEqs_EQ_X(x0, theta_rotation, Nt, db_AA_X, db_BB_X, db_CC_X, db_DD_X, db_EE_X, db_FF_X, db_GGn_X, db_GGp_X, db_GG_X, coeff_force, coeff_torque, Fs, Ts, mode_system, x0(end));
% F_EQ_X
% F_EQ_XII = evalEqs_EQ_XII(x0, theta_rotation, Nt, db_AA_XII, db_BB_XII, db_CC_XII, db_DD_XII, db_EE_XII, db_FF_XII, db_GGn_XII, db_GGp_XII, db_GG_XII, coeff_force, coeff_torque, Fs, Ts, mode_system, x0(end));
% F_EQ_XII
% F_EQ_XII2 = evalEqs_EQ_XII2(x0, theta_rotation, Nt, db_AA_XII, db_BB_XII, db_CC_XII, db_DD_XII, db_EE_XII, db_FF_XII, db_GGn_XII, db_GGp_XII, db_GG_XII, coeff_force, coeff_torque, Fs, Ts, mode_system, x0(end));
% F_EQ_XII2
F_partial = evalEqs_partial(x0, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                        db_M_Fx_r, db_N_Fx_r, db_M_Fy_r, db_N_Fy_r, db_M_Fz_r, db_N_Fz_r, ...
                        db_M_Fx_theta, db_N_Fx_theta, db_M_Fy_theta, db_N_Fy_theta, db_M_Fz_theta, db_N_Fz_theta, ...
                        db_M_Fx_phi, db_N_Fx_phi, db_M_Fy_phi, db_N_Fy_phi, db_M_Fz_phi, db_N_Fz_phi, ...
                        coeff_force, coeff_torque, Fs, Ts, mode_system, x0(end));
F_partial

% [delay_amplitude, delay_phase] = solving_nonlinear_system(x0, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
%                         db_M_Fx, db_N_Fx, db_M_Fy, db_N_Fy, db_M_Fz, db_N_Fz, db_M_Tx, db_N_Tx, db_M_Ty, db_N_Ty, db_M_Tz, db_N_Tz, ...
%                         db_M_Fx_r, db_N_Fx_r, db_M_Fy_r, db_N_Fy_r, db_M_Fz_r, db_N_Fz_r, ...
%                         db_M_Fx_theta, db_N_Fx_theta, db_M_Fy_theta, db_N_Fy_theta, db_M_Fz_theta, db_N_Fz_theta, ...
%                         db_M_Fx_phi, db_N_Fx_phi, db_M_Fy_phi, db_N_Fy_phi, db_M_Fz_phi, db_N_Fz_phi, ...
%                         coeff_force, coeff_torque, Fs, Ts, dFs, mode_system);
[delay_amplitude, delay_phase] = optimize_nonlinear_system(x0, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                        db_M_Fx, db_N_Fx, db_M_Fy, db_N_Fy, db_M_Fz, db_N_Fz, db_M_Tx, db_N_Tx, db_M_Ty, db_N_Ty, db_M_Tz, db_N_Tz, ...
                        db_M_Fx_r, db_N_Fx_r, db_M_Fy_r, db_N_Fy_r, db_M_Fz_r, db_N_Fz_r, ...
                        db_M_Fx_theta, db_N_Fx_theta, db_M_Fy_theta, db_N_Fy_theta, db_M_Fz_theta, db_N_Fz_theta, ...
                        db_M_Fx_phi, db_N_Fx_phi, db_M_Fy_phi, db_N_Fy_phi, db_M_Fz_phi, db_N_Fz_phi, ...
                        coeff_force, coeff_torque, Fs, Ts, dFs, mode_system);
if mode_system == 2
    delay_phase = Amp_pha(1+Nt:2*Nt);
end
 
%%

x = [delay_amplitude delay_phase];
if mode_system == 1
    x = x;
elseif mode_system == 2
    x = x(1:Nt);
elseif mode_system == 3
    x = [x(1+Nt:2*Nt) x(1)];
elseif mode_system == 4
    x = [x(1+Nt:2*Nt) x(1)];
end
F = evalEqs(x, theta_rotation, Nt, db_M_Fx, db_N_Fx, db_M_Fy, db_N_Fy, db_M_Fz, db_N_Fz, db_M_Tx, db_N_Tx, db_M_Ty, db_N_Ty, db_M_Tz, db_N_Tz, coeff_force, coeff_torque, [0 0 0], [0 0 0], mode_system, x0(end));
F
F(1:3) 
F(4:6)
F_partial = evalEqs_partial(x, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                        db_M_Fx_r, db_N_Fx_r, db_M_Fy_r, db_N_Fy_r, db_M_Fz_r, db_N_Fz_r, ...
                        db_M_Fx_theta, db_N_Fx_theta, db_M_Fy_theta, db_N_Fy_theta, db_M_Fz_theta, db_N_Fz_theta, ...
                        db_M_Fx_phi, db_N_Fx_phi, db_M_Fy_phi, db_N_Fy_phi, db_M_Fz_phi, db_N_Fz_phi, ...
                        coeff_force, coeff_torque, [0 0 0], [0 0 0], mode_system, x0(end));
F_partial

%% clear database for characteristic matrix

delete([db_filename_CharMat, '.mat']);
delete([db_filename_CharMat_partial, '.mat']);
delete([db_filename_ts_partial, '.mat']);
% delete([db_filename_EQ_VIII, '.mat']);
% delete([db_filename_EQ_X, '.mat']);
% delete([db_filename_EQ_XII, '.mat']);

%%