function F = systemEqs(x, theta_rotation, Nt, r_pd, theta_pd, phi_pd, ...
                M_Fx, N_Fx, M_Fy, N_Fy, M_Fz, N_Fz, M_Tx, N_Tx, M_Ty, N_Ty, M_Tz, N_Tz, ...
                M_Fx_r, N_Fx_r, M_Fy_r, N_Fy_r, M_Fz_r, N_Fz_r, ...
                M_Fx_theta, N_Fx_theta, M_Fy_theta, N_Fy_theta, M_Fz_theta, N_Fz_theta, ...
                M_Fx_phi, N_Fx_phi, M_Fy_phi, N_Fy_phi, M_Fz_phi, N_Fz_phi, ...
                coeff_force, coeff_torque, Fs, Ts, dFs, mode, amp)
%%
% The system of non-linear equations.
% NOTE: consider the characteristic matrix for radiation force is much
% larger than the characteristic matrix for radiation torque, we need to
% normalize the characteristic matrix for radiation force and make the both
% types of characteristic matrix into a similar magnitude.
%%

enhanec_factor = 10;

% F = zeros(1,5);
F = zeros(1,9);
F_r = zeros(1,3);
F_theta = zeros(1,3);
F_phi = zeros(1,3);

% G = (G / coeff_force);
% F(1) = F(1) - 0;
% F(2) = F(2) - 0;
% F(3) = F(3) - G;
% 
% design_Tz = (design_Tz / coeff_torque);
% F(4) = F(4) - 0;
% F(5) = F(5) - 0;
% F(6) = F(6) - design_Tz;
Fs = (Fs / coeff_force);
F(1) = Fs(1);
F(2) = Fs(2);
F(3) = Fs(3);

Ts = (Ts / coeff_torque);
F(4) = Ts(1);
F(5) = Ts(2);
F(6) = Ts(3);

dFs = (dFs / coeff_force);
F(7) = dFs(1);
F(8) = dFs(2);
F(9) = dFs(3);

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
                
% for optimization problem, partial F_x,y,z / partial x,y,z are not include
% here, so F(7), F(8), F(9) are all return zero.
dr_dx = 0; dr_dy = 0; dr_dz = 0; dtheta_dx = 0; dtheta_dy = 0; dtheta_dz = 0; dphi_dx = 0; dphi_dy = 0; dphi_dz = 0;
F(7) = 0;
F(8) = 0;
F(9) = 0;


%% the system of non-linear equations

% if mode == 1
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx(ii,jj));
%             
%             % for radiation force in y-axis
%             F(2) = F(2) - 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy(ii,jj));
%             
%             % for radiation force in z-axis
%             F(3) = F(3) + 1*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz(ii,jj));
%             
%             % for radiation torque in x-axis
%             F(4) = F(4) - 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj));
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) - 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj));
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) - 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj));
%         end
%     end
% end
if mode == 1
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F(1) = F(1) - (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F(2) = F(2) - (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F(3) = F(3) - (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz(ii,jj)) * Rzyx(3,3) );
            
            % for radiation torque in x-axis
            F(4) = F(4) - (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
                           + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation torque in y-axis
            F(5) = F(5) - (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
                           + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation torque in z-axis
            F(6) = F(6) - (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
                           + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    
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
        end
    end 
    % partial Fx / partial x
    F(7) = F(7) - (F_r(1) * dr_dx + F_theta(1) * dtheta_dx + F_phi(1) * dphi_dx);
    % partial Fy / partial y
    F(8) = F(8) - (F_r(2) * dr_dy + F_theta(2) * dtheta_dy + F_phi(2) * dphi_dy);
    % partial Fz / partial z
    F(9) = F(9) - (F_r(3) * dr_dz + F_theta(3) * dtheta_dz + F_phi(3) * dphi_dz);
end

            
            
%% all transducer in-phase

% if mode == 2
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) + 0.5*x(ii)*x(jj) * N_Fx(ii,jj);
%             
%             % for radiation force in y-axis
%             F(2) = F(2) - 0.5*x(ii)*x(jj) * M_Fy(ii,jj);
%             
%             % for radiation force in z-axis
%             F(3) = F(3) - 1*x(ii)*x(jj) * N_Fz(ii,jj);
%             
%             % for radiation torque in x-axis
%             F(4) = F(4) + 0.5*x(ii)*x(jj) * M_Tx(ii,jj);
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + 0.5*x(ii)*x(jj) * N_Ty(ii,jj);
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + 1*x(ii)*x(jj) * M_Tz(ii,jj);
%         end
%     end
% end
if mode == 2
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F(1) = F(1) - (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F(2) = F(2) - (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F(3) = F(3) - (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz(ii,jj)) * Rzyx(3,3) );
            
            % for radiation torque in x-axis
            F(4) = F(4) - (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,1) ...
                           + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation torque in y-axis
            F(5) = F(5) - (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,2) ...
                           + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation torque in z-axis
            F(6) = F(6) - (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,3) ...
                           + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    
    
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
        end
    end
    % partial Fx / partial x
    F(7) = F(7) - (F_r(1) * dr_dx + F_theta(1) * dtheta_dx + F_phi(1) * dphi_dx);
    % partial Fy / partial y
    F(8) = F(8) - (F_r(2) * dr_dy + F_theta(2) * dtheta_dy + F_phi(2) * dphi_dy);
    % partial Fz / partial z
    F(9) = F(9) - (F_r(3) * dr_dz + F_theta(3) * dtheta_dz + F_phi(3) * dphi_dz);
end
% if mode == 2
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) + x(ii)*x(jj) * N_Fx(ii,jj) * Rzyx(1,1) ...
%                         - x(ii)*x(jj) * M_Fy(ii,jj) * Rzyx(2,1) ...
%                         - x(ii)*x(jj) * N_Fz(ii,jj) * Rzyx(3,1);
%             
%             % for radiation force in y-axis
%             F(2) = F(2) + x(ii)*x(jj) * N_Fx(ii,jj) * Rzyx(1,2) ...
%                         - x(ii)*x(jj) * M_Fy(ii,jj) * Rzyx(2,2) ...
%                         - x(ii)*x(jj) * N_Fz(ii,jj) * Rzyx(3,2);
%             
%             % for radiation force in z-axis
%             F(3) = F(3) + x(ii)*x(jj) * N_Fx(ii,jj) * Rzyx(1,3) ...
%                         - x(ii)*x(jj) * M_Fy(ii,jj) * Rzyx(2,3) ...
%                         - x(ii)*x(jj) * N_Fz(ii,jj) * Rzyx(3,3);
%             
%             % for radiation torque in x-axis
%             F(4) = F(4) + x(ii)*x(jj) * M_Tx(ii,jj) * Rzyx(1,1) ...
%                         + x(ii)*x(jj) * N_Ty(ii,jj) * Rzyx(2,1) ...
%                         + x(ii)*x(jj) * M_Tz(ii,jj) * Rzyx(3,1);
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + x(ii)*x(jj) * M_Tx(ii,jj) * Rzyx(1,2) ...
%                         + x(ii)*x(jj) * N_Ty(ii,jj) * Rzyx(2,2) ...
%                         + x(ii)*x(jj) * M_Tz(ii,jj) * Rzyx(3,2);
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + x(ii)*x(jj) * M_Tx(ii,jj) * Rzyx(1,3) ...
%                         + x(ii)*x(jj) * N_Ty(ii,jj) * Rzyx(2,3) ...
%                         + x(ii)*x(jj) * M_Tz(ii,jj) * Rzyx(3,3);
%         end
%     end
% end


%% all transducer in-amplitude

% if mode == 3
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) + 0.5*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj));
%             
%             % for radiation force in y-axis
%             F(2) = F(2) - 0.5*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj));
%             
%             % for radiation force in z-axis
%             F(3) = F(3) + 1*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj));
%             
%             % for radiation torque in x-axis
%             F(4) = F(4) - 0.5*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj));
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) - 0.5*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj));
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) - 1*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj));
%         end
%     end
% end
if mode == 3
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F(1) = F(1) - (  (-0.5)*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F(2) = F(2) - (  (-0.5)*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F(3) = F(3) - (  (-0.5)*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,3) );
            
            % for radiation torque in x-axis
            F(4) = F(4) - (  0.5*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
                           + 1*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation torque in y-axis
            F(5) = F(5) - (  0.5*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
                           + 1*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation torque in z-axis
            F(6) = F(6) - (  0.5*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(Nt+1)^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
                           + 1*(x(Nt+1)^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    
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
        end
    end        
    % partial Fx / partial x
    F(7) = F(7) - (F_r(1) * dr_dx + F_theta(1) * dtheta_dx + F_phi(1) * dphi_dx);
    % partial Fy / partial y
    F(8) = F(8) - (F_r(2) * dr_dy + F_theta(2) * dtheta_dy + F_phi(2) * dphi_dy);
    % partial Fz / partial z
    F(9) = F(9) - (F_r(3) * dr_dz + F_theta(3) * dtheta_dz + F_phi(3) * dphi_dz);        
end
% if mode == 3
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) + (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Fx(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Fx(ii,jj)) * Rzyx(1,1) ...
%                         - (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Fy(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Fy(ii,jj)) * Rzyx(2,1) ...
%                         - (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Fz(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Fz(ii,jj)) * Rzyx(3,1);
%             
%             % for radiation force in y-axis
%             F(2) = F(2) + (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Fx(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Fx(ii,jj)) * Rzyx(1,2) ...
%                         - (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Fy(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Fy(ii,jj)) * Rzyx(2,2) ...
%                         - (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Fz(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Fz(ii,jj)) * Rzyx(3,2);
%             
%             % for radiation force in z-axis
%             F(3) = F(3) + (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Fx(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Fx(ii,jj)) * Rzyx(1,3) ...
%                         - (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Fy(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Fy(ii,jj)) * Rzyx(2,3) ...
%                         - (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Fz(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Fz(ii,jj)) * Rzyx(3,3);
%             
%             % for radiation torque in x-axis
%             F(4) = F(4) + (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Tx(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                         + (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Ty(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                         + (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Tz(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Tz(ii,jj)) * Rzyx(3,1);
%                     
%             % for radiation torque in y-axis
%             F(5) = F(5) + (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Tx(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                         + (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Ty(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                         + (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Tz(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Tz(ii,jj)) * Rzyx(3,2);
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Tx(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                         + (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Ty(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                         + (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Tz(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Tz(ii,jj)) * Rzyx(3,3);
%         end
%     end
% end

% F;
%% all transducer in-amplitude and amplitude is the same and unchanged

% if mode == 4
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) + 0.5*(amp^2 * sin(x(ii)-x(jj)) * M_Fx(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Fx(ii,jj));
%             
%             % for radiation force in y-axis
%             F(2) = F(2) - 0.5*(amp^2 * cos(x(ii)-x(jj)) * M_Fy(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Fy(ii,jj));
%             
%             % for radiation force in z-axis
%             F(3) = F(3) - 1*(amp^2 * sin(x(ii)-x(jj)) * M_Fz(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Fz(ii,jj));
%             
%             % for radiation torque in x-axis
%             F(4) = F(4) + (amp^2 * cos(x(ii)-x(jj)) * M_Tx(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Tx(ii,jj));
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + (amp^2 * sin(x(ii)-x(jj)) * M_Ty(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Ty(ii,jj));
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (amp^2 * cos(x(ii)-x(jj)) * M_Tz(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Tz(ii,jj));
%         end
%     end
% end
if mode == 4
    for ii = 1 : Nt
        for jj = 1 : Nt
            % for radiation force in x-axis
            F(1) = F(1) - (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F(2) = F(2) - (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F(3) = F(3) - (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,3) );
            
            % for radiation torque in x-axis
            F(4) = F(4) - (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
                           + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation torque in y-axis
            F(5) = F(5) - (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
                           + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation torque in z-axis
            F(6) = F(6) - (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
                           + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F(1) = enhanec_factor * F(1);
    F(2) = enhanec_factor * F(2);
    F(4) = enhanec_factor * F(4);
    F(5) = enhanec_factor * F(5);
    F(6) = enhanec_factor * F(6);
    
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
        end
    end    
    % partial Fx / partial x
    F(7) = F(7) - (F_r(1) * dr_dx + F_theta(1) * dtheta_dx + F_phi(1) * dphi_dx);
    % partial Fy / partial y
    F(8) = F(8) - (F_r(2) * dr_dy + F_theta(2) * dtheta_dy + F_phi(2) * dphi_dy);
    % partial Fz / partial z
    F(9) = F(9) - (F_r(3) * dr_dz + F_theta(3) * dtheta_dz + F_phi(3) * dphi_dz);              
end
% if mode == 4
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) + (amp^2 * sin(x(ii)-x(jj)) * M_Fx(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Fx(ii,jj)) * Rzyx(1,1) ...
%                         - (amp^2 * cos(x(ii)-x(jj)) * M_Fy(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Fy(ii,jj)) * Rzyx(2,1) ...
%                         - (amp^2 * sin(x(ii)-x(jj)) * M_Fz(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Fz(ii,jj)) * Rzyx(3,1);
%             
%             % for radiation force in y-axis
%             F(2) = F(2) + (amp^2 * sin(x(ii)-x(jj)) * M_Fx(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Fx(ii,jj)) * Rzyx(1,2) ...
%                         - (amp^2 * cos(x(ii)-x(jj)) * M_Fy(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Fy(ii,jj)) * Rzyx(2,2) ...
%                         - (amp^2 * sin(x(ii)-x(jj)) * M_Fz(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Fz(ii,jj)) * Rzyx(3,2);
%             
%             % for radiation force in z-axis
%             F(3) = F(3) + (amp^2 * sin(x(ii)-x(jj)) * M_Fx(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Fx(ii,jj)) * Rzyx(1,3) ...
%                         - (amp^2 * cos(x(ii)-x(jj)) * M_Fy(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Fy(ii,jj)) * Rzyx(2,3) ...
%                         - (amp^2 * sin(x(ii)-x(jj)) * M_Fz(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Fz(ii,jj)) * Rzyx(3,3);
%             
%             % for radiation torque in x-axis
%             F(4) = F(4) + (amp^2 * cos(x(ii)-x(jj)) * M_Tx(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
%                         + (amp^2 * sin(x(ii)-x(jj)) * M_Ty(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
%                         + (amp^2 * cos(x(ii)-x(jj)) * M_Tz(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Tz(ii,jj)) * Rzyx(3,1);
%                     
%             % for radiation torque in y-axis
%             F(5) = F(5) + (amp^2 * cos(x(ii)-x(jj)) * M_Tx(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
%                         + (amp^2 * sin(x(ii)-x(jj)) * M_Ty(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
%                         + (amp^2 * cos(x(ii)-x(jj)) * M_Tz(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Tz(ii,jj)) * Rzyx(3,2);
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (amp^2 * cos(x(ii)-x(jj)) * M_Tx(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
%                         + (amp^2 * sin(x(ii)-x(jj)) * M_Ty(ii,jj) + amp^2 * cos(x(ii)-x(jj)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
%                         + (amp^2 * cos(x(ii)-x(jj)) * M_Tz(ii,jj) - amp^2 * sin(x(ii)-x(jj)) * N_Tz(ii,jj)) * Rzyx(3,3);
%         end
%     end
% end

%%
