function [A_delay, phi_delay, error_A_record, error_P_record, error_record] = iterative_system()
%%
% This function is used to iteratively call
% "design_transducer_parameters.m", which is used to solve the system of
% nonlinear equation based on the inputting transducer parameters. And the
% object of this function is find a convergent solution of the transducer
% parameters.
%%

parameters;  % getting the transducer number


A_delay = 2*[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];      % radius = 10mm
phi_delay = [0     pi/4     pi/4      pi/4      pi/4     pi/4     pi/4      pi/4      pi/4     7*pi/4      7*pi/4      7*pi/4     7*pi/4     7*pi/4      7*pi/4     7*pi/4      7*pi/4   7*pi/4     7*pi/4      7*pi/4     7*pi/4     7*pi/4      7*pi/4      7*pi/4     7*pi/4];

% convergent parameter
relax_factor = 1;
max_step = 50;
max_error = 0.001; % 0.005 for nonsphere
step = 1;

% desired radiation force and torque 
gravity_a = 9.8;
particle_mass = 0.00001;
% motion_inertia = (2/5) * particle_mass * particle_radius^2;     % for sphere rotates around the central axis
G = particle_mass * gravity_a;
coeff_force = - 1 / (2 * fluid_rho * fluid_c^2) / fluid_k^2;
coeff_torque = - 1 / (2 * fluid_rho * fluid_c^2) / fluid_k^3;
  
dimensionaless_design_Fx = 0;
dimensionaless_design_Fy = 0;
dimensionaless_design_Fz = 1;
dimensionaless_design_Tx = 0;
dimensionaless_design_Ty = 0;
dimensionaless_design_Tz = 0; 
Fx = dimensionaless_design_Fx * G; 
Fy = dimensionaless_design_Fy * G;
Fz = dimensionaless_design_Fz * G;
Tx = dimensionaless_design_Tx * G / fluid_k;        % make sure the desirable 'Tx' is within reasonable range.
Ty = dimensionaless_design_Ty * G / fluid_k;        % make sure the desirable 'Ty' is within reasonable range.
Tz = dimensionaless_design_Tz * G / fluid_k;        % make sure the desirable 'Tz' is within reasonable range.
Fx = Fx + 0;
Fy = Fy + 0;
Fz = Fz + 0;
Tx = Tx + 0;
Ty = Ty + 0;
Tz = Tz + 0;
Fs = [Fx; Fy; Fz];
Ts = [Tx; Ty; Tz];

dFs = [-0.2, -0.2, -0.1];
  
%% the iterative process for phase_delay only
 
error = +inf;
error_A = error;
error_P = error;
error_record = [];
error_A_record = [];
error_P_record = [];
amp_fix = A_delay(1);
mode_type = 4; 
while (error > max_error) && (step <= max_step) 
     
    [delay_amplitude, delay_phase, F] = design_transducer_parameters(A_delay, phi_delay, coeff_force, coeff_torque, Fs, Ts, dFs, mode_type, amp_fix);
    for amp_ii = 1:length(delay_amplitude)
        if delay_amplitude(amp_ii) < 0
            delay_amplitude = abs(delay_amplitude);
            delay_phase = delay_phase + pi;
        end
    end
    delay_phase = mod(delay_phase, 2*pi);     % limit the range in [0,2*pi]
    
    error_A = (delay_amplitude - A_delay) ./ (relax_factor);
    error_P = (delay_phase - phi_delay) ./ (2*pi);
    error = max(max(abs([error_A; error_P])));
    error_A_record = [error_A_record error_A];
    error_P_record = [error_P_record error_P];
    error_record = [error_record error];

    % the 'error' is defined as the maximum deviation between the desired
    % radiation force and torque and the calculated radiation force and
    % torque based on 'delay_amplitude' and 'delay_phase' in current
    % iterative step.  
    % A_delay = abs(delay_amplitude);
    A_delay = delay_amplitude;
%     A_delay(find(A_delay < 0)) = 0;
    phi_delay = delay_phase;
    phi_delay = mod(phi_delay, 2*pi);     % limit the range in [0,2*pi]
    

      
     
    % A_delay_disp and phi_delay_disp: 
    indices = round(transducer(:,1:2)/(2*trans_radius)) - min(min(round(transducer(:,1:2)/(2*trans_radius)))) + 1;
    for ii = 1 : transducer_number
        A_delay_disp(indices(ii,1), indices(ii,2)) = A_delay(ii);
        phi_delay_disp(indices(ii,1), indices(ii,2)) = phi_delay(ii);
    end
    fprintf('\n');
    fprintf('Step %d for error: %f.\n', step, error);
    fprintf('Amplitude distribution:\n'); disp(roundn(A_delay_disp,-3));
    fprintf('Phase distribution:\n'); disp(roundn(phi_delay_disp*180/pi,-1));
    step = step + 1;
       
    % evaluate the radiation force and torque under current transducer parameters
    if mod(step,5) == 0
        replace_Amp_delay_phi_delay(A_delay, phi_delay);    % rewrite the 'A_delay' and 'phi_delay' into "phase_array_beam_shape_coeff.m".
        [Frad_x,Frad_y,Frad_z,~] = radiation_force_based_Analyses();
        [Torque_x,Torque_y,Torque_z,~] = radiation_torque_based_Analyses();
        Frad_z - F(3)
    end
    
end


%% visualization of the distributions

figure(1);
% visualize the amplitude delay of the current step
[X,Y]=meshgrid([1:size(A_delay_disp,1)+1],[1:size(A_delay_disp,2)+1]);
% adding a surrounding zeros for visualization
A_delay_vis = [A_delay_disp; zeros(1,size(A_delay_disp,2))];
A_delay_vis(1,size(A_delay_disp,1)+1) = 0;
pclr = pcolor(X, Y, A_delay_vis);
xlabel('y-axis [mm]');
ylabel('x-axis [mm]');
% set(pclr, 'LineStyle','none');
%set(pclr, 'box', 'on');
colormap('hot');
%colormap('hsv');
colorbar;
caxis([0,3]);

figure(2);
% visualize the phase delay of the current step
[X,Y]=meshgrid([1:size(phi_delay_disp,1)+1],[1:size(phi_delay_disp,2)+1]);
% adding a surrounding zeros for visualization
phi_delay_vis = [phi_delay_disp*180/pi; zeros(1,size(phi_delay_disp,2))];
phi_delay_vis(1,size(phi_delay_disp,1)+1) = 0;
pclr = pcolor(X, Y, phi_delay_vis);
xlabel('y-axis [mm]');
ylabel('x-axis [mm]');
% set(pclr, 'LineStyle','none');
%set(pclr, 'box', 'on');
%colormap('hot');
colormap('hsv');
colorbar;
caxis([0,360]);

%%