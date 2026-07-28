function [db_filename_CharMat] = database_characteristic_matrix(N, db_bs_coeff, db_s_coeff)
%%
% build the database for "characteristic matrix" to avoid repeating call
% function "transducer_transform_beam_shape_coeff.m",
% "scattering_coefficient.m", "real_element.m", and "imag_element.m", which
% will cost tons of time. 
% 
% database will automatically be created if current folder has not a
% corresponding database. (check the database file exists or not by the
% filename)
% database size: [Nt, Nt] (Nt is the transducer number:
% 'transducer_number' in "parameters.m"); 
%%

parameters_names;

if strcmp(wave_type, 'phase_array_transducer') ~= 1
    error('Database characteristic matrix cannot be created for NOT phase array wave type!\n');
end

if exist([db_filename_CharMat, '.mat']) ~= 0        % if already exist the database, 
    return;                                 % then do not create again for saving time.
end

Nt = transducer_number;


%% obtain the force and torque constant "AA", "BB", "CC", "DD", "EE", "FF", and "GG" 

% coeff_force = - (1)^2 / (2 * fluid_rho * fluid_c^2);

AA = @(nn,mm) (-1) * sqrt(((nn+mm-1) * (nn+mm)) / ((2*nn-1) * (2*nn+1)));
BB = @(nn,mm) sqrt(((nn-mm+2) * (nn-mm+1)) / ((2*nn+1) * (2*nn+3)));
CC = @(nn,mm) sqrt(((nn-mm-1) * (nn-mm)) / ((2*nn-1) * (2*nn+1)));
DD = @(nn,mm) (-1) * sqrt(((nn+mm+2) * (nn+mm+1)) / ((2*nn+1) * (2*nn+3)));
EE = @(nn,mm) sqrt(((nn-mm) * (nn+mm)) / ((2*nn-1) * (2*nn+1)));
FF = @(nn,mm) sqrt(((nn-mm+1) * (nn+mm+1)) / ((2*nn+1) * (2*nn+3)));

% coeff_torque = - (1)^2 / (2 * fluid_rho * fluid_c^2);

GG_p = @(nn,mm) sqrt((nn-mm) * (nn+mm+1));
GG_n = @(nn,mm) sqrt((nn+mm) * (nn-mm+1));

%% prepare the translate beam-shape coefficients

[trans_bs_coeff] = ...
    transducer_transform_beam_shape_coeff(N, db_bs_coeff, fluid_k, transducer, transducer_number, theta_rotation, 1, inf, inf);

%% establish the characteristic matrix for radiation force

db_M_Fx = zeros(Nt, Nt);
db_N_Fx = zeros(Nt, Nt);
db_M_Fy = zeros(Nt, Nt);
db_N_Fy = zeros(Nt, Nt);
db_M_Fz = zeros(Nt, Nt);
db_N_Fz = zeros(Nt, Nt);
for ii = 1:Nt
    for jj = 1:Nt
        
        sum_1_r = 0; sum_1_i = 0;
        sum_2_r = 0; sum_2_i = 0;
        sum_3_r = 0; sum_3_i = 0;
        sum_4_r = 0; sum_4_i = 0;
        sum_5_r = 0; sum_5_i = 0;
        sum_6_r = 0; sum_6_i = 0;
        for nn = 0:(N-1)    % (n+1), (m+1)
            for mm = -nn:(nn)
                n_plus_1 = nn + 1;      % adjust the projecting relation of database and beam-shape and scattering coefficents
                m_plus_1 = mm + 1;
                if abs(m_plus_1) <= n_plus_1    % make sure the adjusting position within the database
                    R = real(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_plus_1+1, n_plus_1+m_plus_1+1));
                    I = imag(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_plus_1+1, n_plus_1+m_plus_1+1));
                    sum_1_r = sum_1_r + AA(nn+1,mm+1) * R;
                    sum_1_i = sum_1_i + AA(nn+1,mm+1) * I;
                end
            end
        end
        
        for nn = 1:(N)      % (n-1), (m+1)
            for mm = -nn:(nn)
                n_minus_1 = nn - 1;
                m_plus_1 = mm + 1;
                if abs(m_plus_1) <= n_minus_1
                    R = real(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_minus_1+1, n_minus_1+m_plus_1+1));
                    I = imag(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_minus_1+1, n_minus_1+m_plus_1+1));
                    sum_2_r = sum_2_r + BB(nn-1,mm+1) * R;
                    sum_2_i = sum_2_i + BB(nn-1,mm+1) * I;
                end
            end
        end
        
        for nn = 0:(N-1)    % (n+1), (m-1)
            for mm = (-nn):(nn)
                n_plus_1 = nn + 1;
                m_minus_1 = mm - 1;
                if abs(m_minus_1) <= n_plus_1
                    R = real(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_plus_1+1, n_plus_1+m_minus_1+1));
                    I = imag(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_plus_1+1, n_plus_1+m_minus_1+1));
                    sum_3_r = sum_3_r + CC(nn+1,mm-1) * R;
                    sum_3_i = sum_3_i + CC(nn+1,mm-1) * I;
                end
            end
        end
        
        for nn = 1:(N)      % (n-1), (m-1)
            for mm = (-nn):(nn)
                n_minus_1 = nn - 1;
                m_minus_1 = mm - 1;
                if abs(m_minus_1) <= n_minus_1
                    R = real(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_minus_1+1, n_minus_1+m_minus_1+1));
                    I = imag(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_minus_1+1, n_minus_1+m_minus_1+1));
                    sum_4_r = sum_4_r + DD(nn-1,mm-1) * R;
                    sum_4_i = sum_4_i + DD(nn-1,mm-1) * I;
                end
            end
        end

        for nn = 0:(N-1)     % (n+1), m
            for mm = -nn:(nn)
                n_plus_1 = nn + 1;
                if abs(mm) <= n_plus_1
                    R = real(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_plus_1+1, n_plus_1+mm+1));
                    I = imag(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_plus_1+1, n_plus_1+mm+1));
                    sum_5_r = sum_5_r + EE(nn+1,mm) * R;
                    sum_5_i = sum_5_i + EE(nn+1,mm) * I;
                end
            end
        end
        
        for nn = 1:(N)      % (n-1), m
            for mm = -(nn):(nn)
                n_minus_1 = nn - 1;
                if abs(mm) <= n_minus_1
                    R = real(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_minus_1+1, n_minus_1+mm+1));
                    I = imag(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, n_minus_1+1, n_minus_1+mm+1));
                    sum_6_r = sum_6_r + FF(nn-1,mm) * R;
                    sum_6_i = sum_6_i + FF(nn-1,mm) * I;
                end
            end
        end
       
        db_M_Fx(ii, jj) = sum_1_r - sum_2_r + sum_3_r - sum_4_r;
        db_N_Fx(ii, jj) = sum_1_i - sum_2_i + sum_3_i - sum_4_i;
        db_M_Fy(ii, jj) = sum_1_r - sum_2_r - sum_3_r + sum_4_r;
        db_N_Fy(ii, jj) = sum_1_i - sum_2_i - sum_3_i + sum_4_i;
        db_M_Fz(ii, jj) = sum_5_r - sum_6_r;
        db_N_Fz(ii, jj) = sum_5_i - sum_6_i;
        
    end
    
    fprintf('Characteristic Matrix Database Preparing %d%% \n', ...
        roundn(100*ii/Nt,0)/2);
    
end


%% establish the characteristic matrix for radiation torque

db_M_Tx = zeros(Nt, Nt);
db_N_Tx = zeros(Nt, Nt);
db_M_Ty = zeros(Nt, Nt);
db_N_Ty = zeros(Nt, Nt);
db_M_Tz = zeros(Nt, Nt);
db_N_Tz = zeros(Nt, Nt);
for ii = 1:Nt
    for jj = 1:Nt
        
        sum_1_r = 0; sum_1_i = 0;
        sum_2_r = 0; sum_2_i = 0;
        sum_3_r = 0; sum_3_i = 0;
        for nn = 0:(N)    % (n), (m-1)
            for mm = -nn:(nn)
                % adjust the projecting relation of database and beam-shape and scattering coefficents
                m_minus_1 = mm - 1;
                if abs(m_minus_1) <= nn    % make sure the adjusting position within the database
                    R = real(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, nn+1, nn+m_minus_1+1));
                    I = imag(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, nn+1, nn+m_minus_1+1));
                    sum_1_r = sum_1_r + GG_n(nn,mm) * R;
                    sum_1_i = sum_1_i + GG_n(nn,mm) * I;
                end
            end
        end

        for nn = 0:(N)      % (n), (m+1)
            for mm = -nn:(nn)
                % adjust the projecting relation of database and beam-shape and scattering coefficents
                m_plus_1 = mm + 1;
                if abs(m_plus_1) <= nn    % make sure the adjusting position within the database
                    R = real(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, nn+1, nn+m_plus_1+1));
                    I = imag(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, nn+1, nn+m_plus_1+1));
                    sum_2_r = sum_2_r + GG_p(nn,mm) * R;
                    sum_2_i = sum_2_i + GG_p(nn,mm) * I;
                end
            end
        end

        for nn = 0:(N)     % (n), m
            for mm = -nn:(nn)
                R = real(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, nn+1, nn+mm+1));
                I = imag(combination_element(trans_bs_coeff, db_s_coeff, ii, jj, nn+1, nn+mm+1, nn+1, nn+mm+1));
                sum_3_r = sum_3_r + mm * R;
                sum_3_i = sum_3_i + mm * I;
            end
        end

        db_M_Tx(ii, jj) = sum_1_r + sum_2_r;
        db_N_Tx(ii, jj) = sum_1_i + sum_2_i;
        db_M_Ty(ii, jj) = sum_2_r - sum_1_r;
        db_N_Ty(ii, jj) = sum_2_i - sum_1_i;
        db_M_Tz(ii, jj) = sum_3_r;
        db_N_Tz(ii, jj) = sum_3_i;
        
    end
    
    fprintf('Characteristic Matrix Database Preparing %d%% \n', ...
        50 + roundn(100*ii/Nt,0)/2);
    
end

%% save these databases by meaningful name

% r = max(max(abs(db_s_coeff)))/max(max(abs(db_bs_coeff)))
% f_ka 

save([db_filename_CharMat, '.mat'], 'trans_bs_coeff', 'db_M_Fx', 'db_N_Fx', 'db_M_Fy', 'db_N_Fy', 'db_M_Fz', 'db_N_Fz', ...
                                'db_M_Tx', 'db_N_Tx', 'db_M_Ty', 'db_N_Ty', 'db_M_Tz', 'db_N_Tz');


%%