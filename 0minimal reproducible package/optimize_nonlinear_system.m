function [delay_amplitude, delay_phase] = optimize_nonlinear_system(x0, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                db_M_Fx, db_N_Fx, db_M_Fy, db_N_Fy, db_M_Fz, db_N_Fz, db_M_Tx, db_N_Tx, db_M_Ty, db_N_Ty, db_M_Tz, db_N_Tz, ...
                db_M_Fx_r, db_N_Fx_r, db_M_Fy_r, db_N_Fy_r, db_M_Fz_r, db_N_Fz_r, ...
                db_M_Fx_theta, db_N_Fx_theta, db_M_Fy_theta, db_N_Fy_theta, db_M_Fz_theta, db_N_Fz_theta, ...
                db_M_Fx_phi, db_N_Fx_phi, db_M_Fy_phi, db_N_Fy_phi, db_M_Fz_phi, db_N_Fz_phi, ...
                coeff_force, coeff_torque, Fs, Ts, dFs, mode)
%%
% This is function is used to solve a possible root of a given non-linear
% system of equations.
%%

options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations', 1e5);  % for sphere
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations', 1e5);   % for nonsphere 

F_Eqs = @(x) systemEqs(x, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                db_M_Fx, db_N_Fx, db_M_Fy, db_N_Fy, db_M_Fz, db_N_Fz, db_M_Tx, db_N_Tx, db_M_Ty, db_N_Ty, db_M_Tz, db_N_Tz, ...
                db_M_Fx_r, db_N_Fx_r, db_M_Fy_r, db_N_Fy_r, db_M_Fz_r, db_N_Fz_r, ...
                db_M_Fx_theta, db_N_Fx_theta, db_M_Fy_theta, db_N_Fy_theta, db_M_Fz_theta, db_N_Fz_theta, ...
                db_M_Fx_phi, db_N_Fx_phi, db_M_Fy_phi, db_N_Fy_phi, db_M_Fz_phi, db_N_Fz_phi, ...
                coeff_force, coeff_torque, Fs, Ts, dFs, mode, x0(end));

% contraints by Lyapunov
F_partial = @(x) evalEqs_partial(x, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                        db_M_Fx_r, db_N_Fx_r, db_M_Fy_r, db_N_Fy_r, db_M_Fz_r, db_N_Fz_r, ...
                        db_M_Fx_theta, db_N_Fx_theta, db_M_Fy_theta, db_N_Fy_theta, db_M_Fz_theta, db_N_Fz_theta, ...
                        db_M_Fx_phi, db_N_Fx_phi, db_M_Fy_phi, db_N_Fy_phi, db_M_Fz_phi, db_N_Fz_phi, ...
                        coeff_force, coeff_torque, Fs, Ts, mode, x0(end));

Eqs = @(x) sum(F_Eqs(x).^2) + sum((F_partial(x).*[1,1,1]).^2);
                    
myConstraints = @(x) LyapunovConstraints(x, F_partial);

%% the system of non-linear equations

if mode == 1
%     x0 = [ones(1, Nt) zeros(1, Nt)];
%     x0 = x0;
    lb = zeros(1, 2*Nt);
    ub = ones(1, Nt) * 10;
    ub = ones(Nt+1, 2*Nt) * 2*pi;
    [x, fval] = fmincon(Eqs, x0, [],[],[],[], lb, ub, myConstraints, options);
    
    delay_amplitude = x(1 : Nt);
    delay_phase = x(Nt+1 : end);
end

%% all transducer in-phase

if mode == 2
%     x0 = ones(1, Nt);
%     x0 = x0(1:Nt);
    lb = zeros(1, Nt);
    ub = ones(1, Nt) * 10;
    [x, fval] = fmincon(Eqs, x0, [],[],[],[], lb, ub, myConstraints, options);
    
    delay_amplitude = x;
    delay_phase = zeros(1, Nt);
end

%% all transducer in-amplitude
 
if mode == 3
%     x0 = [zeros(1, Nt) 1];
%     x0 = [x0(1+Nt:2*Nt) x0(1)];
    lb = zeros(1, Nt+1);
    ub = ones(1, Nt) * 2*pi;
    ub = [ub, 10];
    [x, fval] = fmincon(Eqs, x0, [],[],[],[], lb, ub, myConstraints, options);
    
    delay_amplitude = x(Nt+1) * ones(1, Nt);
    delay_phase = x(1:Nt);
end

%% all transducer in-amplitude and amplitudes are the same and unchanged

if mode == 4
%     x0 = [zeros(1, Nt) 1];
%     x0 = [x0(1+Nt:2*Nt) x0(1)];
    lb = zeros(1, Nt);
    ub = ones(1, Nt) * 2*pi;
    [x, fval] = fmincon(Eqs, x0, [],[],[],[], lb, ub, myConstraints, options);
    
    delay_amplitude = x0(Nt+1) * ones(1, Nt);
    delay_phase = x(1:Nt);
end

%%