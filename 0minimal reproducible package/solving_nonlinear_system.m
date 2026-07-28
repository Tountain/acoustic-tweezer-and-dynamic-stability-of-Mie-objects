function [delay_amplitude, delay_phase] = solving_nonlinear_system(x0, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                db_M_Fx, db_N_Fx, db_M_Fy, db_N_Fy, db_M_Fz, db_N_Fz, db_M_Tx, db_N_Tx, db_M_Ty, db_N_Ty, db_M_Tz, db_N_Tz, ...
                db_M_Fx_r, db_N_Fx_r, db_M_Fy_r, db_N_Fy_r, db_M_Fz_r, db_N_Fz_r, ...
                db_M_Fx_theta, db_N_Fx_theta, db_M_Fy_theta, db_N_Fy_theta, db_M_Fz_theta, db_N_Fz_theta, ...
                db_M_Fx_phi, db_N_Fx_phi, db_M_Fy_phi, db_N_Fy_phi, db_M_Fz_phi, db_N_Fz_phi, ...
                coeff_force, coeff_torque, Fs, Ts, dFs, mode)
%%
% This is function is used to solve a possible root of a given non-linear
% system of equations.
%%

options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'MaxFunEvals', 10000, 'MaxIter', 5000, 'TolFun', 1e-15);
% options = optimset('Display','iter');
Eqs = @(x) systemEqs(x, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                db_M_Fx, db_N_Fx, db_M_Fy, db_N_Fy, db_M_Fz, db_N_Fz, db_M_Tx, db_N_Tx, db_M_Ty, db_N_Ty, db_M_Tz, db_N_Tz, ...
                db_M_Fx_r, db_N_Fx_r, db_M_Fy_r, db_N_Fy_r, db_M_Fz_r, db_N_Fz_r, ...
                db_M_Fx_theta, db_N_Fx_theta, db_M_Fy_theta, db_N_Fy_theta, db_M_Fz_theta, db_N_Fz_theta, ...
                db_M_Fx_phi, db_N_Fx_phi, db_M_Fy_phi, db_N_Fy_phi, db_M_Fz_phi, db_N_Fz_phi, ...
                coeff_force, coeff_torque, Fs, Ts, dFs, mode, x0(end));

%% the system of non-linear equations

if mode == 1
%     x0 = [ones(1, Nt) zeros(1, Nt)];
%     x0 = x0;
    x = fsolve(Eqs, x0, options);
    
    delay_amplitude = x(1 : Nt);
    delay_phase = x(Nt+1 : end);
end

%% all transducer in-phase

if mode == 2
%     x0 = ones(1, Nt);
%     x0 = x0(1:Nt);
    x = fsolve(Eqs, x0, options);
    
    delay_amplitude = x;
    delay_phase = zeros(1, Nt);
end

%% all transducer in-amplitude
 
if mode == 3
%     x0 = [zeros(1, Nt) 1];
%     x0 = [x0(1+Nt:2*Nt) x0(1)];
    x = fsolve(Eqs, x0, options);
    
    delay_amplitude = x(Nt+1) * ones(1, Nt);
    delay_phase = x(1:Nt);
end

%% all transducer in-amplitude and amplitudes are the same and unchanged

if mode == 4
%     x0 = [zeros(1, Nt) 1];
%     x0 = [x0(1+Nt:2*Nt) x0(1)];
    x = fsolve(Eqs, x0(1:Nt), options);
    
    delay_amplitude = x0(Nt+1) * ones(1, Nt);
    delay_phase = x(1:Nt);
end

%%