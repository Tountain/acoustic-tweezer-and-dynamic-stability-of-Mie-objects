function F = evalEqs(x, theta_rotation, Nt, M_Fx, N_Fx, M_Fy, N_Fy, M_Fz, N_Fz, M_Tx, N_Tx, M_Ty, N_Ty, M_Tz, N_Tz, coeff_force, coeff_torque, Fs, Ts, mode, amp)
%%
% Evaluate the system of non-linear equations.
%%

% F = zeros(1,5);
F = zeros(1,6);

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

%% the system of non-linear equations

Fx1 = 0;
Fx2 = 0;
Fy1 = 0;
Fy2 = 0;
Fz1 = 0;
Fz2 = 0;
Tx1 = 0;
Tx2 = 0;
Ty1 = 0;
Ty2 = 0;
Tz1 = 0;
Tz2 = 0;
% if mode == 1
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             Fx1 = Fx1 + (x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx(ii,jj));
%             
%             Fx2 = Fx2 + (x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx(ii,jj));
% 
%             % for radiation force in y-axis
%             Fy1 = Fy1 + (x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy(ii,jj));
%             
%             Fy2 = Fy2 + (x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy(ii,jj));
% 
%             % for radiation force in z-axis
%             Fz1 = Fz1 + (x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz(ii,jj));
%             
%             Fz2 = Fz2 + (x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz(ii,jj));
% 
%             % for radiation torque in x-axis
%             Tx1 = Tx1 + (x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj));
%             
%             Tx2 = Tx2 + (x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj));
% 
%             % for radiation torque in y-axis
%             Ty1 = Ty1 + (x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj));
%            
%             Ty2 = Ty2 + (x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj));
%         
%             % for radiation torque in z-axis
%             Tz1 = Tz1 + (x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj));
%             
%             Tz2 = Tz2 + (x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj));
%             
%         end
%     end
%     
%     % for radiation force in x-axis
%     F(1) = - 0.5 * coeff_force * (Fx1 + Fx2);
% 
%     % for radiation force in y-axis
%     F(2) = + 0.5 * coeff_force * (Fy1 - Fy2);
% 
%     % for radiation force in z-axis
%     F(3) = - 1 * coeff_force * (Fz1 + Fz2);
% 
%     % for radiation torque in x-axis
%     F(4) = + 0.5 * coeff_torque * (Tx1 - Tx2);
% 
%     % for radiation torque in y-axis
%     F(5) = + 0.5 * coeff_torque * (Ty1 + Ty2);
% 
%     % for radiation torque in z-axis
%     F(6) = + 1 * coeff_torque * (Tz1 - Tz2);  
%     
% end
if mode == 1
    for ii = 1 : Nt
        for jj = 1 : Nt
            F(1) = F(1) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F(2) = F(2) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F(3) = F(3) + (  (-0.5)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Fy(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Fz(ii,jj)) * Rzyx(3,3) );
            
            % for radiation torque in x-axis
            F(4) = F(4) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
                           + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation torque in y-axis
            F(5) = F(5) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
                           + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation torque in z-axis
            F(6) = F(6) + (  0.5*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
                           + 1*(x(ii)*x(jj)*cos(x(jj+Nt)-x(ii+Nt)) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(x(jj+Nt)-x(ii+Nt)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F(1) = F(1) * coeff_force;
    F(2) = F(2) * coeff_force;
    F(3) = F(3) * coeff_force;
    F(4) = F(4) * coeff_torque;
    F(5) = F(5) * coeff_torque;
    F(6) = F(6) * coeff_torque;
end


%% all transducer in-phase

% if mode == 2
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) - x(ii)*x(jj) * N_Fx(ii,jj);
% 
%             % for radiation force in y-axis
%             F(2) = F(2) + x(ii)*x(jj) * M_Fy(ii,jj);
% 
%             % for radiation force in z-axis
%             F(3) = F(3) - x(ii)*x(jj) * N_Fz(ii,jj);
% 
%             % for radiation torque in x-axis
%             F(4) = F(4) + x(ii)*x(jj) * M_Tx(ii,jj);
% 
%             % for radiation torque in y-axis
%             F(5) = F(5) + x(ii)*x(jj) * N_Ty(ii,jj);
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + x(ii)*x(jj) * M_Tz(ii,jj);
%         end
%     end
% end
% if mode == 2
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) + 0.5*x(ii)*x(jj) * N_Fx(ii,jj) * Rzyx(1,1) ...
%                         - 0.5*x(ii)*x(jj) * M_Fy(ii,jj) * Rzyx(2,1) ...
%                         - 0.5*x(ii)*x(jj) * N_Fz(ii,jj) * Rzyx(3,1);
%             
%             % for radiation force in y-axis
%             F(2) = F(2) + 0.5*x(ii)*x(jj) * N_Fx(ii,jj) * Rzyx(1,2) ...
%                         - 0.5*x(ii)*x(jj) * M_Fy(ii,jj) * Rzyx(2,2) ...
%                         - 0.5*x(ii)*x(jj) * N_Fz(ii,jj) * Rzyx(3,2);
%             
%             % for radiation force in z-axis
%             F(3) = F(3) + x(ii)*x(jj) * N_Fx(ii,jj) * Rzyx(1,3) ...
%                         - x(ii)*x(jj) * M_Fy(ii,jj) * Rzyx(2,3) ...
%                         - x(ii)*x(jj) * N_Fz(ii,jj) * Rzyx(3,3);
%             
%             % for radiation torque in x-axis
%             F(4) = F(4) + 0.5*x(ii)*x(jj) * M_Tx(ii,jj) * Rzyx(1,1) ...
%                         + 0.5*x(ii)*x(jj) * N_Ty(ii,jj) * Rzyx(2,1) ...
%                         + 0.5*x(ii)*x(jj) * M_Tz(ii,jj) * Rzyx(3,1);
%             
%             % for radiation torque in y-axis
%             F(5) = F(5) + 0.5*x(ii)*x(jj) * M_Tx(ii,jj) * Rzyx(1,2) ...
%                         + 0.5*x(ii)*x(jj) * N_Ty(ii,jj) * Rzyx(2,2) ...
%                         + 0.5*x(ii)*x(jj) * M_Tz(ii,jj) * Rzyx(3,2);
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + x(ii)*x(jj) * M_Tx(ii,jj) * Rzyx(1,3) ...
%                         + x(ii)*x(jj) * N_Ty(ii,jj) * Rzyx(2,3) ...
%                         + x(ii)*x(jj) * M_Tz(ii,jj) * Rzyx(3,3);
%         end
%     end
% end
if mode == 2
    for ii = 1 : Nt
        for jj = 1 : Nt
            F(1) = F(1) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F(2) = F(2) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F(3) = F(3) + (  (-0.5)*(x(ii)*x(jj)*sin(0) * M_Fx(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*cos(0) * M_Fy(ii,jj) - x(ii)*x(jj)*sin(0) * N_Fy(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(x(ii)*x(jj)*sin(0) * M_Fz(ii,jj) + x(ii)*x(jj)*cos(0) * N_Fz(ii,jj)) * Rzyx(3,3) );
            
            % for radiation torque in x-axis
            F(4) = F(4) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,1) ...
                           + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation torque in y-axis
            F(5) = F(5) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,2) ...
                           + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation torque in z-axis
            F(6) = F(6) + (  0.5*(x(ii)*x(jj)*cos(0) * M_Tx(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(ii)*x(jj)*sin(0) * M_Ty(ii,jj) + x(ii)*x(jj)*cos(0) * N_Ty(ii,jj)) * Rzyx(2,3) ...
                           + 1*(x(ii)*x(jj)*cos(0) * M_Tz(ii,jj) - x(ii)*x(jj)*sin(0) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F(1) = F(1) * coeff_force;
    F(2) = F(2) * coeff_force;
    F(3) = F(3) * coeff_force;
    F(4) = F(4) * coeff_torque;
    F(5) = F(5) * coeff_torque;
    F(6) = F(6) * coeff_torque;
end


%% all transducer in-amplitude

% if mode == 3
%     for ii = 1 : Nt
%         for jj = 1 : Nt
%             % for radiation force in x-axis
%             F(1) = F(1) - (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Fx(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Fx(ii,jj));
% 
%             % for radiation force in y-axis
%             F(2) = F(2) + (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Fy(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Fy(ii,jj));
% 
%             % for radiation force in z-axis
%             F(3) = F(3) - (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Fz(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Fz(ii,jj));
% 
%             % for radiation torque in x-axis
%             F(4) = F(4) + (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Tx(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Tx(ii,jj));
% 
%             % for radiation torque in y-axis
%             F(5) = F(5) + (x(Nt+1)^2 * sin(x(ii)-x(jj)) * M_Ty(ii,jj) + x(Nt+1)^2 * cos(x(ii)-x(jj)) * N_Ty(ii,jj));
%             
%             % for radiation torque in z-axis
%             F(6) = F(6) + (x(Nt+1)^2 * cos(x(ii)-x(jj)) * M_Tz(ii,jj) - x(Nt+1)^2 * sin(x(ii)-x(jj)) * N_Tz(ii,jj));
%         end
%     end
% end
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
if mode == 3
    for ii = 1 : Nt
        for jj = 1 : Nt
            F(1) = F(1) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F(2) = F(2) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F(3) = F(3) + (  (-0.5)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fx(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Fy(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Fz(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,3) );
            
            % for radiation torque in x-axis
            F(4) = F(4) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
                           + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation torque in y-axis
            F(5) = F(5) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
                           + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation torque in z-axis
            F(6) = F(6) + (  0.5*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tx(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(x(Nt+1)^2*sin(x(jj)-x(ii)) * M_Ty(ii,jj) + x(Nt+1)^2*cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
                           + 1*(x(Nt+1)^2*cos(x(jj)-x(ii)) * M_Tz(ii,jj) - x(Nt+1)^2*sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F(1) = F(1) * coeff_force;
    F(2) = F(2) * coeff_force;
    F(3) = F(3) * coeff_force;
    F(4) = F(4) * coeff_torque;
    F(5) = F(5) * coeff_torque;
    F(6) = F(6) * coeff_torque;
end


%% all transducer in-amplitude and amplitude is the same and unchanged

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
if mode == 4
    for ii = 1 : Nt
        for jj = 1 : Nt
            F(1) = F(1) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,1) ...
                           + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation force in y-axis
            F(2) = F(2) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,2) ...
                           + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation force in z-axis
            F(3) = F(3) + (  (-0.5)*(amp^2 * sin(x(jj)-x(ii)) * M_Fx(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Fy(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Fy(ii,jj)) * Rzyx(2,3) ...
                           + (-1)*(amp^2 * sin(x(jj)-x(ii)) * M_Fz(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Fz(ii,jj)) * Rzyx(3,3) );
            
            % for radiation torque in x-axis
            F(4) = F(4) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,1) ...
                           + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,1) ...
                           + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,1) );
            
            % for radiation torque in y-axis
            F(5) = F(5) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,2) ...
                           + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,2) ...
                           + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,2) );
            
            % for radiation torque in z-axis
            F(6) = F(6) + (  0.5*(amp^2 * cos(x(jj)-x(ii)) * M_Tx(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tx(ii,jj)) * Rzyx(1,3) ...
                           + 0.5*(amp^2 * sin(x(jj)-x(ii)) * M_Ty(ii,jj) + amp^2 * cos(x(jj)-x(ii)) * N_Ty(ii,jj)) * Rzyx(2,3) ...
                           + 1*(amp^2 * cos(x(jj)-x(ii)) * M_Tz(ii,jj) - amp^2 * sin(x(jj)-x(ii)) * N_Tz(ii,jj)) * Rzyx(3,3) );
        end
    end
    F(1) = F(1) * coeff_force;
    F(2) = F(2) * coeff_force;
    F(3) = F(3) * coeff_force;
    F(4) = F(4) * coeff_torque;
    F(5) = F(5) * coeff_torque;
    F(6) = F(6) * coeff_torque;
end

%%