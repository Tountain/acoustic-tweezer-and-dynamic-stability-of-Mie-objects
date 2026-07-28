function F = evalEqs_partial(x, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                            M_Fx_r, N_Fx_r, M_Fy_r, N_Fy_r, M_Fz_r, N_Fz_r, ...
                            M_Fx_theta, N_Fx_theta, M_Fy_theta, N_Fy_theta, M_Fz_theta, N_Fz_theta, ...
                            M_Fx_phi, N_Fx_phi, M_Fy_phi, N_Fy_phi, M_Fz_phi, N_Fz_phi, ...
                            coeff_force, coeff_torque, Fs, Ts, mode, amp)
%%
% Evaluate the system of partial non-linear equations.
%%

F = zeros(1,3);
F_r = zeros(1,3);
F_theta = zeros(1,3);
F_phi = zeros(1,3);

% G = (G / coeff_force);
% F(3) = F(3) - G;


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

[dr_dx, dr_dy, dr_dz, dtheta_dx, dtheta_dy, dtheta_dz, dphi_dx, dphi_dy, dphi_dz] ...
                    = SphericaltoCartesian_partial(r_pd, theta_pd, phi_pd);
                
%% the system of non-linear equations


if mode == 1
    % partial r
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_r(1) = F_r(1) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx_r(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx_r(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy_r(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy_r(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz_r(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz_r(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_r(2) = F_r(2) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx_r(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx_r(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy_r(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy_r(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz_r(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz_r(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_r(3) = F_r(3) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx_r(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx_r(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy_r(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy_r(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz_r(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz_r(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_r(1) = F_r(1) * coeff_force;
    F_r(2) = F_r(2) * coeff_force;
    F_r(3) = F_r(3) * coeff_force;
%     F(4) = F(4) * coeff_torque;
%     F(5) = F(5) * coeff_torque;
%     F(6) = F(6) * coeff_torque;

    % partial theta
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_theta(1) = F_theta(1) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx_theta(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx_theta(ii,jj)) * Rzyx(1,1) ...
                                    + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy_theta(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy_theta(ii,jj)) * Rzyx(2,1) ...
                                    + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz_theta(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz_theta(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_theta(2) = F_theta(2) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx_theta(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx_theta(ii,jj)) * Rzyx(1,2) ...
                                    + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy_theta(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy_theta(ii,jj)) * Rzyx(2,2) ...
                                    + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz_theta(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz_theta(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_theta(3) = F_theta(3) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx_theta(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx_theta(ii,jj)) * Rzyx(1,3) ...
                                    + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy_theta(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy_theta(ii,jj)) * Rzyx(2,3) ...
                                    + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz_theta(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz_theta(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_theta(1) = F_theta(1) * coeff_force;
    F_theta(2) = F_theta(2) * coeff_force;
    F_theta(3) = F_theta(3) * coeff_force;
    
    % partial phi
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_phi(1) = F_phi(1) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx_phi(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx_phi(ii,jj)) * Rzyx(1,1) ...
                                    + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy_phi(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy_phi(ii,jj)) * Rzyx(2,1) ...
                                    + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz_phi(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz_phi(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_phi(2) = F_phi(2) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx_phi(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx_phi(ii,jj)) * Rzyx(1,2) ...
                                    + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy_phi(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy_phi(ii,jj)) * Rzyx(2,2) ...
                                    + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz_phi(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz_phi(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_phi(3) = F_phi(3) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx_phi(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx_phi(ii,jj)) * Rzyx(1,3) ...
                                    + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy_phi(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy_phi(ii,jj)) * Rzyx(2,3) ...
                                    + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz_phi(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz_phi(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_phi(1) = F_phi(1) * coeff_force;
    F_phi(2) = F_phi(2) * coeff_force;
    F_phi(3) = F_phi(3) * coeff_force;
end


%% all transducer in-phase


if mode == 2
    % partial r
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_r(1) = F_r(1) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx_r(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx_r(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy_r(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy_r(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz_r(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz_r(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_r(2) = F_r(2) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx_r(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx_r(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy_r(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy_r(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz_r(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz_r(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_r(3) = F_r(3) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx_r(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx_r(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy_r(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy_r(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz_r(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz_r(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_r(1) = F_r(1) * coeff_force;
    F_r(2) = F_r(2) * coeff_force;
    F_r(3) = F_r(3) * coeff_force;
%     F(4) = F(4) * coeff_torque;
%     F(5) = F(5) * coeff_torque;
%     F(6) = F(6) * coeff_torque;
    
    % partial theta
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_theta(1) = F_theta(1) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx_theta(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx_theta(ii,jj)) * Rzyx(1,1) ...
                                    + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy_theta(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy_theta(ii,jj)) * Rzyx(2,1) ...
                                    + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz_theta(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz_theta(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_theta(2) = F_theta(2) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx_theta(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx_theta(ii,jj)) * Rzyx(1,2) ...
                                    + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy_theta(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy_theta(ii,jj)) * Rzyx(2,2) ...
                                    + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz_theta(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz_theta(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_theta(3) = F_theta(3) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx_theta(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx_theta(ii,jj)) * Rzyx(1,3) ...
                                    + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy_theta(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy_theta(ii,jj)) * Rzyx(2,3) ...
                                    + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz_theta(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz_theta(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_theta(1) = F_theta(1) * coeff_force;
    F_theta(2) = F_theta(2) * coeff_force;
    F_theta(3) = F_theta(3) * coeff_force;
    
    % partial phi
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_phi(1) = F_phi(1) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx_phi(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx_phi(ii,jj)) * Rzyx(1,1) ...
                                    + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy_phi(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy_phi(ii,jj)) * Rzyx(2,1) ...
                                    + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz_phi(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz_phi(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_phi(2) = F_phi(2) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx_phi(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx_phi(ii,jj)) * Rzyx(1,2) ...
                                    + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy_phi(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy_phi(ii,jj)) * Rzyx(2,2) ...
                                    + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz_phi(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz_phi(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_phi(3) = F_phi(3) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx_phi(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx_phi(ii,jj)) * Rzyx(1,3) ...
                                    + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy_phi(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy_phi(ii,jj)) * Rzyx(2,3) ...
                                    + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz_phi(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz_phi(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_phi(1) = F_phi(1) * coeff_force;
    F_phi(2) = F_phi(2) * coeff_force;
    F_phi(3) = F_phi(3) * coeff_force;
end


%% all transducer in-amplitude


if mode == 3
    % partial r
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_r(1) = F_r(1) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx_r(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx_r(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy_r(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy_r(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz_r(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz_r(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_r(2) = F_r(2) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx_r(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx_r(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy_r(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy_r(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz_r(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz_r(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_r(3) = F_r(3) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx_r(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx_r(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy_r(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy_r(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz_r(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz_r(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_r(1) = F_r(1) * coeff_force;
    F_r(2) = F_r(2) * coeff_force;
    F_r(3) = F_r(3) * coeff_force;
%     F(4) = F(4) * coeff_torque;
%     F(5) = F(5) * coeff_torque;
%     F(6) = F(6) * coeff_torque;

    % partial theta
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_theta(1) = F_theta(1) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx_theta(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx_theta(ii,jj)) * Rzyx(1,1) ...
                                    + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy_theta(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy_theta(ii,jj)) * Rzyx(2,1) ...
                                    + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz_theta(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz_theta(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_theta(2) = F_theta(2) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx_theta(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx_theta(ii,jj)) * Rzyx(1,2) ...
                                    + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy_theta(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy_theta(ii,jj)) * Rzyx(2,2) ...
                                    + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz_theta(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz_theta(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_theta(3) = F_theta(3) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx_theta(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx_theta(ii,jj)) * Rzyx(1,3) ...
                                    + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy_theta(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy_theta(ii,jj)) * Rzyx(2,3) ...
                                    + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz_theta(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz_theta(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_theta(1) = F_theta(1) * coeff_force;
    F_theta(2) = F_theta(2) * coeff_force;
    F_theta(3) = F_theta(3) * coeff_force;
    
    % partial phi
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_phi(1) = F_phi(1) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx_phi(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx_phi(ii,jj)) * Rzyx(1,1) ...
                                    + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy_phi(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy_phi(ii,jj)) * Rzyx(2,1) ...
                                    + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz_phi(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz_phi(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_phi(2) = F_phi(2) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx_phi(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx_phi(ii,jj)) * Rzyx(1,2) ...
                                    + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy_phi(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy_phi(ii,jj)) * Rzyx(2,2) ...
                                    + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz_phi(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz_phi(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_phi(3) = F_phi(3) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx_phi(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx_phi(ii,jj)) * Rzyx(1,3) ...
                                    + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy_phi(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy_phi(ii,jj)) * Rzyx(2,3) ...
                                    + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz_phi(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz_phi(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_phi(1) = F_phi(1) * coeff_force;
    F_phi(2) = F_phi(2) * coeff_force;
    F_phi(3) = F_phi(3) * coeff_force;
end


%% all transducer in-amplitude and amplitude is the same and unchanged


if mode == 4
    % partial r
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_r(1) = F_r(1) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx_r(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx_r(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy_r(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy_r(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz_r(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz_r(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_r(2) = F_r(2) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx_r(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx_r(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy_r(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy_r(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz_r(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz_r(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_r(3) = F_r(3) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx_r(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx_r(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy_r(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy_r(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz_r(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz_r(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_r(1) = F_r(1) * coeff_force;
    F_r(2) = F_r(2) * coeff_force;
    F_r(3) = F_r(3) * coeff_force;
%     F(4) = F(4) * coeff_torque;
%     F(5) = F(5) * coeff_torque;
%     F(6) = F(6) * coeff_torque;

    % partial theta
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_theta(1) = F_theta(1) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx_theta(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx_theta(ii,jj)) * Rzyx(1,1) ...
                                    + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy_theta(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy_theta(ii,jj)) * Rzyx(2,1) ...
                                    + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz_theta(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz_theta(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_theta(2) = F_theta(2) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx_theta(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx_theta(ii,jj)) * Rzyx(1,2) ...
                                    + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy_theta(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy_theta(ii,jj)) * Rzyx(2,2) ...
                                    + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz_theta(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz_theta(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_theta(3) = F_theta(3) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx_theta(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx_theta(ii,jj)) * Rzyx(1,3) ...
                                    + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy_theta(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy_theta(ii,jj)) * Rzyx(2,3) ...
                                    + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz_theta(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz_theta(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_theta(1) = F_theta(1) * coeff_force;
    F_theta(2) = F_theta(2) * coeff_force;
    F_theta(3) = F_theta(3) * coeff_force;
    
    % partial phi
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F_phi(1) = F_phi(1) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx_phi(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx_phi(ii,jj)) * Rzyx(1,1) ...
                                    + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy_phi(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy_phi(ii,jj)) * Rzyx(2,1) ...
                                    + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz_phi(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz_phi(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F_phi(2) = F_phi(2) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx_phi(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx_phi(ii,jj)) * Rzyx(1,2) ...
                                    + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy_phi(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy_phi(ii,jj)) * Rzyx(2,2) ...
                                    + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz_phi(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz_phi(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F_phi(3) = F_phi(3) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx_phi(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx_phi(ii,jj)) * Rzyx(1,3) ...
                                    + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy_phi(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy_phi(ii,jj)) * Rzyx(2,3) ...
                                    + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz_phi(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz_phi(ii,jj)) * Rzyx(3,3) );
            
%             % for radiation torque in x-axis
%             F(4) = F(4) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                            + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                            + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                            + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                            + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                            + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                            + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F_phi(1) = F_phi(1) * coeff_force;
    F_phi(2) = F_phi(2) * coeff_force;
    F_phi(3) = F_phi(3) * coeff_force;
end

%%

% partial Fx / partial x
F(1) = F_r(1) * dr_dx + F_theta(1) * dtheta_dx + F_phi(1) * dphi_dx;
% partial Fy / partial y
F(2) = F_r(2) * dr_dy + F_theta(2) * dtheta_dy + F_phi(2) * dphi_dy;
% partial Fz / partial z
F(3) = F_r(3) * dr_dz + F_theta(3) * dtheta_dz + F_phi(3) * dphi_dz;

%%