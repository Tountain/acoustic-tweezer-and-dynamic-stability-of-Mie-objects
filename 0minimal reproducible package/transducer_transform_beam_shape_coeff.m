function [coeff] = ...
    transducer_transform_beam_shape_coeff(N, db_bs_coeff, fluid_k, transducer, transducer_number, theta_rotation, mode, A_delay, phi_delay)
%%
% This function is used to calculate the translate beam-shape coefficients
% of every transducer based on the relative position relationship.
%
% Writing structure is mainly referred to code
% "phase_array_beam_shape_coeff.m".
%%

% the phase array information
% transducer(1, :) = [0, 0, 0] * (2*trans_radius);         % probe transducer
% transducer(2, :) = [1, 0, 0] * (2*trans_radius);         % source transducer(s)
% transducer(3, :) = [-1, 0, 0] * (2*trans_radius);
% transducer(4, :) = [0, 1, 0] * (2*trans_radius);
% transducer(5, :) = [0, -1, 0] * (2*trans_radius);
% transducer(6, :) = [1, 1, 0] * (2*trans_radius); 
% transducer(7, :) = [1, -1, 0] * (2*trans_radius); 
% transducer(8, :) = [-1, 1, 0] * (2*trans_radius);
% transducer(9, :) = [-1, -1, 0] * (2*trans_radius);
% transducer(1, :) = [0, 0, 0] * (2*trans_radius);         % probe transducer
% transducer(2, :) = [1, 0, 0] * (2*trans_radius);         % source transducer(s)
% transducer(3, :) = [1, 1, 0] * (2*trans_radius);
% transducer(4, :) = [0, 1, 0] * (2*trans_radius);
% transducer(5, :) = [-1, 1, 0] * (2*trans_radius);
% transducer(6, :) = [-1, 0, 0] * (2*trans_radius); 
% transducer(7, :) = [-1, -1, 0] * (2*trans_radius); 
% transducer(8, :) = [0, -1, 0] * (2*trans_radius);
% transducer(9, :) = [1, -1, 0] * (2*trans_radius);
% transducer(1, :) = [0, 0, 0] * (2*trans_radius);         % probe transducer
% transducer(2, :) = [1, 0, 0] * (2*trans_radius);         % source transducer(s)
% transducer(3, :) = [0, 1, 0] * (2*trans_radius);
% transducer(4, :) = [-1, 0, 0] * (2*trans_radius); 
% transducer(5, :) = [0, -1, 0] * (2*trans_radius);

    
if transducer_number ~= size(transducer, 1)
    error('The transducer number is wrong.\n');
end

% rotation of relative positions among the transducers
% here 'theta_rotation' should take a opposite direction
Rx=[1 0 0;
    0 cos(theta_rotation(1)) -sin(theta_rotation(1)); 
    0 sin(theta_rotation(1)) cos(theta_rotation(1))];
Ry=[cos(theta_rotation(2)) 0 sin(theta_rotation(2));
    0 1 0; 
    -sin(theta_rotation(2)) 0 cos(theta_rotation(2))];
Rz=[cos(theta_rotation(3)) -sin(theta_rotation(3)) 0; 
    sin(theta_rotation(3)) cos(theta_rotation(3)) 0; 
    0 0 1];
Rxyz = Rx * Ry * Rz;
for ii = 1:transducer_number
    if ii == 1
        continue;
    end
    transducer(ii, :) = (transducer(ii, :) - transducer(1, :)) * Rxyz;
end


%% calculating the 'translate beam-shape coefficient' of each transducer

if mode == 1
   
    equ_A = cell(1, transducer_number);
    for ii = 1 : transducer_number
        equ_A{ii} = db_bs_coeff * 1;
    end

    % the relative position of source transducers r.w.t the probe transducer (ref. "database_translation_coeffs.m") 
    r_iq = zeros(transducer_number, 1);     % "iq" means particle "i" toward probe transducer "q"
    theta_iq = zeros(transducer_number, 1);
    phi_iq = zeros(transducer_number, 1);
    % number "ii == 1" represents the probe transducer 'q'
    for ii = 1:transducer_number
        if ii == 1
            continue;
        end
        [r_iq(ii,1), theta_iq(ii,1), phi_iq(ii,1)] = ...
            coords_system_relative_positions_general(transducer(ii, :), transducer(1, :));
    end
    kr_iq = fluid_k * r_iq;

    % the transform beam-shape coefficient 3D matrix: bs_nmi_q (ref. "database_translation_coeffs.m")
    trans_bs_coeff = cell(1, transducer_number);
    for ii = 1:transducer_number
        if ii == 1                  % except to probe transducer 'ii==1'
            trans_bs_coeff{ii} = equ_A{ii};
        else
            trans_bs_coeff{ii} = transform_beam_shape_matrix(equ_A{ii}, kr_iq(ii), theta_iq(ii), phi_iq(ii), N, ii);
        end
        fprintf('Transform Beam-Shape Coefficients Database Preparing %d%% \n', ...
            round(100*ii/transducer_number));
    end
    
    coeff = trans_bs_coeff;

%% calculating the 'equivalent beam-shape coefficient' of the phase array
elseif mode == 2
    
    phase_delay_vector = A_delay .* exp(1i * -phi_delay);              % loading the 'phase_delay_vector'
    
    if transducer_number ~= length(A_delay)
        error('The size of the transducer parameter is wrong.\n');
    end
    
    equ_A = cell(1, transducer_number);
    for ii = 1 : transducer_number
        equ_A{ii} = db_bs_coeff * phase_delay_vector(ii);
    end
    
    r_iq = zeros(transducer_number, 1);     % "iq" means particle "i" toward probe transducer "q"
    theta_iq = zeros(transducer_number, 1);
    phi_iq = zeros(transducer_number, 1);
    % number "ii == 1" represents the probe transducer 'q'
    for ii = 1:transducer_number
        if ii == 1
            continue;
        end
        [r_iq(ii,1), theta_iq(ii,1), phi_iq(ii,1)] = ...
            coords_system_relative_positions_general(transducer(ii, :), transducer(1, :));
    end
    kr_iq = fluid_k * r_iq;

    % the transform beam-shape coefficient 3D matrix: bs_nmi_q (ref. "database_translation_coeffs.m")
    bs_nmi_q = cell(1, transducer_number);
    for ii = 1:transducer_number
        if ii == 1                  % except to probe transducer 'ii==1'
            bs_nmi_q{ii} = equ_A{ii};
        else
            bs_nmi_q{ii} = transform_beam_shape_matrix(equ_A{ii}, kr_iq(ii), theta_iq(ii), phi_iq(ii), N, ii);
        end
        fprintf('Transform Beam-Shape Coefficients Database Preparing %d%% \n', ...
            round(100*ii/transducer_number));
    end

    % the equivalent beam-shape coefficient: equivalent_bs_coeff
    % summing up all 2D matrix layers in the 3D matrix 'bs_nmi_q' in element by
    % element.
    equivalent_bs_q_coeff = zeros(size(bs_nmi_q{1}));
    for ii = 1:transducer_number
        equivalent_bs_q_coeff = equivalent_bs_q_coeff + bs_nmi_q{ii};
    end
    
    coeff = equivalent_bs_q_coeff;
    
else
    
    error('The ''mode'' should be 1 or 2!\n');
    
end
%%