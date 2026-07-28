%%
% save the parameters, which uses to control the properties of the system.
%%

%%
% ================ fundimental parameters ==================
 
irregular_body = 1;                 % 1: the particle is axisymmetric arbitrary nonspherical object; 0: the particle is spherical object.

if irregular_body == 0
    rot_matrix = 0;                 % default value.
    particle_radius = 5*10^-3;   % particle radius (a); 30 um;
    Cn = [particle_radius, 0];
    theta_x = 0;
    theta_y = 0;
    theta_z = 0;
    theta_rotation = [theta_x, theta_y, theta_z];
elseif irregular_body == 1
    rot_matrix = 0;         % 1: activate the addition rotation theorem; 0: inactivate the addition rotation theorem. 
    particle_radius = 10*10^-3;                   % particle radial mean-radius (a); 
    %particle_major_length = 2*particle_radius * 2;    % the major length (along axisymmetric axis) of the irregular particle
    Cn = [particle_radius, 0];                    % spherical
    theta_x = 0;                % particle counter-clock wise rotation (position Tx) along x-axis for positive 'theta_x'
    theta_y = 0;                % particle counter-clock wise rotation (position Ty) along y-axis for positive 'theta_y'
    theta_z = 0;                % particle counter-clock wise rotation (position Tz) along z-axis for positive 'theta_z'
    theta_rotation = [theta_x, theta_y, theta_z];
    % 'theta_rotation' follow right-hand principle to determine the
    % particle rotating direction.
    
    % NOTE: rot_matrix ~= 1 and rot_matrix == 1, they follow different
    % rotational principles. 
    %
    % For rot_matrix ~= 1, we employ the rotational transformation, and the
    % rotational principle is the object rotates anti-clockwise by about
    % 'theta_x' along the x-axis to Oxy'z'; then, the object rotates
    % anti-clockwise by 'theta_y' along the y'-axis to Ox'y'z''; finally,
    % the object rotates anti-clockwise by 'theta_z' along the z''-axis to
    % Ox''y''z''.
    % For rot_matrix == 1, we employ the addition rotation theorem, and the
    % rotational principle is the object rotates anti-clockwise by about
    % 'theta_x' (here, 'theta_x' means rotates along z-axis) along the
    % z-axis to Ox'y'z; then, the object rotates anti-clockwise by
    % 'theta_y' along the y'-axis to Ox''y'z'; finally, the object rotates
    % anti-clockwise by 'theta_z' along the z'-axis to Ox'''y''z'.
    %
    % The different between these two rotational principles is the first
    % rotation process, one rotates along x-axis, while the other along
    % z-axis. In other words, if 'theta_x' remain 'theta_x=0', the above
    % two rotational process will be consistent.
end

% derivative_field_is_needed == 1 is designed for phase array
derivative_field_is_needed = 1;     % 1: if the partial derivative (translation) of beam-shape and scalar scattering coefficients, as well as derivative of radiation force/torque for judgement of Lyapunov stability, are needed.
if derivative_field_is_needed == 1
    if rot_matrix == 0 %&& theta_x == 0 && theta_y == 0 && theta_z == 0
        derivativeX = 0 - 1*particle_radius/1000000;      % need to set deviationX + 1*particle_radius/10000
        derivativeY = 0 - 1*particle_radius/1000000;      % need to set deviationY + 1*particle_radius/10000
        derivativeZ = 0 - 1*particle_radius/1000000;      % need to set deviationZ + 1*particle_radius/10000
    else
        error('Need to inactive rotation matrix for partial derivatives (translation) of coefficients and radiation effects.\n');
    end
end

% derivative_rotationfield_is_needed == 1 is designed for phase array
derivative_rotationfield_is_needed = 0;     % 1: if the partial derivative (rotation) of beam-shape and scalar scattering coefficients, as well as derivative of radiation force/torque for judgement of Lyapunov stability, are needed.
if derivative_rotationfield_is_needed == 1
    if rot_matrix == 1
        derivative_thetaX = 0;
        derivative_thetaY = 0 - pi/4;                   % need to set theta_y + pi/4
        derivative_thetaZ = 0;
    else
        error('Need to active rotation matrix for partial derivatives (rotation) of coefficients and radiation effects.\n');
    end
    
    if (theta_x + derivative_thetaX ~= 0) || (theta_y + derivative_thetaY ~= 0) || (theta_z + derivative_thetaZ ~= 0)
        error('Remain the object is in its standard orientation. \n');
    end
end

freq = 40000;               % 40 kHz
omega = freq*2*pi;

%sound_intensity = 33*10^3;         % w/m2

% transducers parameters used in calculation
v0 = 1;                                                                    % transducer vibration velocity amplitude, m/s, assuming v0 = 1 m/s
inter_dist = 0.05;                                                         % m, distance between the transducer and the particle center (i.e., the origin of the coordinate system)

% ==========================================================
%%

%%
% ============ Database&Grid control parameters ============

multi_particle = 0;          % 1: multi-particle system; 0: single-particle system "pressure_contour.m" and "radiation_force_based_Analyses.m"
multi_particle_force_torque = 0*multi_particle;     % function valid only if 'multi_particle == 1'
                             % 2: obtain force&torque of all particle of multi-particle system; 
                             %      Use for "multi_radiation_force_based_Analyses.m"
                             % 1: obtain force&torque of singe probe particle of multi-particle system; 
                             %      Use for "radiation_force_based_Analyses.m" and "radiation_torque_based_Analyses.m" 
                             % 0: obtain pressure field of multi-particle system (if multi_particle == 1).
                             %      Use for "pressure_contour.m" and "pressure_contour_vertical.m" 

vector_amination_visualization = 0;                 % function valid for both 'multi_particle == 0' and 'multi_particle == 1'
                             % 1: for visualizing the force vectors or making animation for particle's moving; 
                             %      Use for "visual_radiation_force_vector.m" and "animation.m" 
                             % 0: for visualizing the pressure contour or getting single position's forces and torques. 
                             %      Use for "pressure_contour.m", "pressure_contour_vertical.m", "radiation_force_based_Analyses.m" and "radiation_torque_based_Analyses.m" 

range_r_coeff = 5;          % the range enlarge coefficient of radius
grid_resolution = 150;       % the grid notes along r-axis and theta-axis

if multi_particle == 0       % database maximum size (larger than Truncation requirement) for outsidest root
    if length(Cn) >= 3
        db_size_nn = ceil(omega * (Cn(1) + abs(Cn(3))) / 340)+6;         % for single-particle system, maximum size of database larger than 30 will waste the computational resource.  
    else
        db_size_nn = 15; % 13
    end
    if db_size_nn < 12
        db_size_nn = 15; % 13
    end
    if db_size_nn > 15
        db_size_nn = 15; % 15
    end
%     db_size_nn = 12;    % phase array a=2mm
%     db_size_nn = 20;    % phase array a=5mm
    db_size_nn = 22;    % phase array a=10mm
else
    % db_size_nn = 25;         % for multi spherical particle system, maximum size of database larger than 20 will out of MATLAB computational memory.
    db_size_nn = 25;         % for multi-axisymmetric objects, if db_size_nn > 8, the condition number of the eigenmatric will be very small (smaller than 10^-20).
end

% ==========================================================
%%
 
%% 
% =================== particle properties ==================

BC = 'rigid';               % Boundary Conditions "rigid" or "compressible".
% BC = 'soft';
% BC = 'olive_oil';             % Note: for irregular_body == 1, if BC ~= 'rigid', then it corresponds to sound soft BC='soft'. Only for irregular_body == 0, the compressible BC can work. 
% BC = 'benzene';
% BC = 'polyurethane';

[particle_rho, particle_k] = material_property(BC, freq);
% particle_rho = 30; particle_k = 0;  % rigid EPS
% particle_mass = 4/3 * pi * particle_radius^3 * particle_rho;
% motion_inertia = (2/5) * particle_mass * particle_radius^3;

% ==========================================================
%%

%%
% =================== medium properties ====================

fluid = 'air';                  % medium types "air" or "water" and others
%fluid = 'mixture';              % "mixture" means considering attenuation  
% fluid = 'water';
% fluid = 'co2';                  % medium types "air" or "water" and others

%======= define the Mie scattering crystals' moving velocity ======
if multi_particle == 1
    velocity = [0, 0, 0] * particle_radius/1;    % particle radius per second
end
%==================================================================

% for "mixture (attenuation)" fluid only, the ratio of particles' volume to the whole
% calculating domain volume
particles_volume_fraction = 0.064;   
% (Details of parameters setting refers "absorption_manual.txt")
% =========== AIR, rigid, 100um, 1MHz (kr=1.85) ============
% particles_volume_fraction = 0.065  -> particles_occupy_per_meter = 50.0%, Dis = 4.0r cr N = 11 
% particles_volume_fraction = 0.049  -> particles_occupy_per_meter = 45.5%, Dis = 4.4r    N = 12  
% particles_volume_fraction = 0.043  -> particles_occupy_per_meter = 43.5%, Dis = 4.6r cr  
% particles_volume_fraction = 0.034  -> particles_occupy_per_meter = 40.0%, Dis = 5.0r    N = 13
% particles_volume_fraction = 0.025  -> particles_occupy_per_meter = 36.4%, Dis = 5.5r    N = 14  
% particles_volume_fraction = 0.023  -> particles_occupy_per_meter = 35.1%, Dis = 5.7r cr  
% particles_volume_fraction = 0.019  -> particles_occupy_per_meter = 33.3%, Dis = 6.0r    N = 15 
% particles_volume_fraction = 0.012  -> particles_occupy_per_meter = 28.6%, Dis = 7.0r    N = 18 
% particles_volume_fraction = 0.008  -> particles_occupy_per_meter = 25.0%, Dis = 8.0r cr N = 20 
% =========== AIR, rigid, 100um, 1MHz (kr=1.85) ============   

[fluid_c, fluid_k, fluid_rho, fluid_viscosity] = fluid_property(fluid, omega); 

f_ka = real(fluid_k)*particle_radius;

% ==========================================================
%%

%%
% ===== compressible particle BC for scattering coeffs======

% for multi layer media, the 'fluid_rho' and 'fluid_k' need to revise
gama = (fluid_rho * particle_k) / (particle_rho * fluid_k);


% ==========================================================
%%

%%
% ========= incident progress Directions&Deviation =========

wave_type = 'phase_array_transducer';                 % based on addition theorem


if strcmp(wave_type, 'standing_plain') == 1
    %direction = 'X';
    %direction = 'Y';
    direction = 'Z';
    %direction = 'XY'; (functionless)
    %direction = 'XZ';
    
%     symbol = 'positive';                % DON'T CHANGE! for two-direction (X&Z-axis) standing wave, we define symbol always "positive", then saving file looks order-well
    symbol = 'arbitrary';               % for arbitrary incidence of standing wave at potential well ONLY.
    if strcmp(symbol, 'arbitrary') == 1 % 'direction' must be 'Z'
        direction = 'Z';                % (theta_inc, phi_inc) = (0, 0) for along Z-axis
        theta_inc = 6*pi/12+0.005;                  % (theta_inc, phi_inc) = (pi/2, 0) for along X-axis
        phi_inc = 0;                    % (theta_inc, phi_inc) = (pi/2, pi/2) for along Y-axis
        if theta_inc < 0                % theta_inc in [0, pi]
            theta_inc = theta_inc+pi;
        elseif theta_inc > pi
            theta_inc = theta_inc-pi;
        end
    end
    
    %================= attenuation adjustion ===============
    % as determination of attenuation coefficient of mixture
    % medium is based on single direction plain 
    % (Wei_ronghao@1965@). therefore, for standing plain wave
    % or multi-plaine waves, we need to turn back to single
    % plain wave for attenuation coefficient "absp_factor". 
    % once "absp_factor ~=0" means the coefficient already be
    % determined, then below code will not run.
    global absp_factor          % "absorption_factor.m"
    if strcmp(fluid, 'mixture') == 1 && absp_factor == 0
        wave_type = 'plain';    % temporally change "standing_plain" as "plain"
        direction = 'Z';
    end
    %================= attenuation adjustion ===============
    
elseif strcmp(wave_type, 'plain') == 1
    %direction = 'X';
    %direction = 'Y';
    direction = 'Z';
    
    %symbol = 'positive';
    %symbol = 'negative';
    symbol = 'arbitrary';               % for arbitrary incidence of plane wave mode ONLY.
    if strcmp(symbol, 'arbitrary') == 1 % 'direction' must be 'Z'
        direction = 'Z';                % (theta_inc, phi_inc) = (0, 0) for along Z-axis
        theta_inc = 0;                  % (theta_inc, phi_inc) = (pi/2, 0) for along X-axis
        phi_inc = 0;                    % (theta_inc, phi_inc) = (pi/2, pi/2) for along Y-axis
        if theta_inc < 0                % theta_inc in [0, pi]
            theta_inc = theta_inc+pi;
        elseif theta_inc > pi
            theta_inc = theta_inc-pi;
        end
    end
elseif strcmp(wave_type, 'phase_array_transducer') == 1 || strcmp(wave_type, 'phase_array_transducer2') == 1
    direction = 'Z';
    symbol = 'positive';
else  
    %direction = 'X';
    %direction = 'Y';
    direction = 'Z';
    
    symbol = 'positive';
    %symbol = 'negative';
end

deviationX = 0.0;
deviationY = 0;
deviationZ = 0.0;
if vector_amination_visualization == 1 || multi_particle_force_torque == 2
    global X_vec Y_vec Z_vec        % for "visual_radiation_force_vector.m", "animation.m" and "multi_radiation_force_based_Analyses.m"
    deviationX = X_vec;             % particle deviation from wave input positions (if "XZ" case, specific for "X")
    deviationY = Y_vec;             % particle deviation from wave input positions (if "XY" case, specific for "Y")??
    deviationZ = Z_vec;             % particle deviation from wave input positions (if "XZ" case, specific for "Z")
end

[dir_sign, x_translation, y_translation, z_translation] = ...
    wave_direction_particle_deviation(wave_type, direction, symbol, deviationX, deviationY, deviationZ);

% ==========================================================
%%

%%
% ==== incident progress wave Amplitude&Phase of Micros ====

%p_inlet = round(sqrt(2 * sound_intensity * fluid_c * fluid_rho));      % initial amplitude given based on input intensity
p_inlet = 1;               % initial amplitude of proba particle [Pa] if no attenuation effects
phase_inlet = 0;            % initial phase of proba particle [rad]

% adjusting amplitude of proba particle [Pa] to keep microphones' amplitude
% consistent no matter exist attenuation or not.
[p_0] = proba_particle_input_amplitude(p_inlet, fluid, fluid_k, particle_radius, range_r_coeff);

% ==========================================================
%%

%%
% =========== incident progress wave Type input ============
% The 'p_i' will be used in "beam_shape_coeff.m" for the numerical
% beam-shape coefficient. However, for the transducer scenarios, as the
% beam-shape coefficient of the wave functions can not be numerical
% decomposite independently with the position vector vec{r}. So, for the
% wave functions other than 'plain', 'zero-Bessel' and 'non-zero-Bessel'
% types, we have to use Fourier expansion to simplify the wave functions as
% a series of plain waves, and each of them can be express by
% position-independent beam-shape coefficient. 
% Therefore, for complex wave function, the below codes are useless, while
% a new counterpart function "G_i.m" is used to obtain the transerve shape
% function for further Fourier expansion.

% Bessel beam parameters
BN = 0;
BETA = 1.309;
X0 = 0;
Y0 = 0; 
Z0 = 0;

translated_Bessel_beam = 1;
beam_source = [0, 0, -0.00] + eps;

% transducer array parameters
trans_radius = 0.005;           % NOTE: keep same in "total_effects_all_plain_wave_component.m". radius of transducer is 5mm
transducer_number = 25;          % NOTE: keep same in "total_effects_all_plain_wave_component.m". for phase array
% the phase array information
% transducer(1, :) = [0, 0, 0] * (2*trans_radius);         % probe transducer
% transducer(2, :) = [1, 0, 0] * (2*trans_radius);         % source transducer(s)
% transducer(3, :) = [-1, 0, 0] * (2*trans_radius);
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
transducer(1, :) = [0, 0, 0] * (2*trans_radius);         % probe transducer
transducer(2, :) = [1, 0, 0] * (2*trans_radius);         % source transducer(s)
transducer(3, :) = [1, 1, 0] * (2*trans_radius);
transducer(4, :) = [0, 1, 0] * (2*trans_radius);
transducer(5, :) = [-1, 1, 0] * (2*trans_radius);
transducer(6, :) = [-1, 0, 0] * (2*trans_radius); 
transducer(7, :) = [-1, -1, 0] * (2*trans_radius); 
transducer(8, :) = [0, -1, 0] * (2*trans_radius);
transducer(9, :) = [1, -1, 0] * (2*trans_radius);
transducer(10, :) = [2, -1, 0] * (2*trans_radius);
transducer(11, :) = [2, 0, 0] * (2*trans_radius);
transducer(12, :) = [2, 1, 0] * (2*trans_radius);
transducer(13, :) = [2, 2, 0] * (2*trans_radius);
transducer(14, :) = [1, 2, 0] * (2*trans_radius);
transducer(15, :) = [0, 2, 0] * (2*trans_radius);
transducer(16, :) = [-1, 2, 0] * (2*trans_radius);
transducer(17, :) = [-2, 2, 0] * (2*trans_radius);
transducer(18, :) = [-2, 1, 0] * (2*trans_radius);
transducer(19, :) = [-2, 0, 0] * (2*trans_radius);
transducer(20, :) = [-2, -1, 0] * (2*trans_radius);
transducer(21, :) = [-2, -2, 0] * (2*trans_radius);
transducer(22, :) = [-1, -2, 0] * (2*trans_radius);
transducer(23, :) = [0, -2, 0] * (2*trans_radius);
transducer(24, :) = [1, -2, 0] * (2*trans_radius);
transducer(25, :) = [2, -2, 0] * (2*trans_radius);

theta_Fourier = 1.2496;         % used for plane wave expansion feature 
phi_Fourier = 0.7854;           % used for plane wave expansion feature 
Pr_0 = - fluid_rho * fluid_c * fluid_k * trans_radius^2 * v0 * 1i / 2;      % the transducer output power for 'Dims_Fourier = 3'.  

% wavefront function
[p_i] = ...
    wave_function(p_0, wave_type, direction, dir_sign, ...
    fluid_k, x_translation, y_translation, z_translation, trans_radius, transducer_number, Pr_0, inter_dist, theta_rotation, rot_matrix);

if (strcmp(wave_type, 'plain') == 1 || strcmp(wave_type, 'standing_plain') == 1) ...
        && BETA ~= 0 && X0 ~= 0 && Y0 ~= 0 && Z0 ~= 0
    error('No wave deviations for plain wave!\n');
end
if strcmp(wave_type, 'zero-Bessel') == 1 && BN ~= 0
    error('Bn == 0 for zero Bessel beam!\n');
end
if strcmp(wave_type, 'non-zero-Bessel') == 1 && BN < 1
    error('Bn >= 1 for non-zero Bessel beam!\n');
end


%%
% other settings for integration for the beam-shape coefficients

if (strcmp(wave_type, 'phase_array_transducer') == 1 || strcmp(wave_type, 'phase_array_transducer2') == 1)
    % for the transducer array cases, b is set to '0.6*2*transducer_radius' and normally 'transducer_radius is 0.005'; ref to TANG@2022@PRE of Fig. 4.  
    %b = 2.4 * 2 * trans_radius;
    %b = 2.4 * 1 * trans_radius;     % have to perform an independent analysis of 'b', and please try to avoid numerically determine the beam-shape coefficient!!!
    b = 21/4 * pi / fluid_k;        % kb = 13/4 * pi; in this way, kb >> 1 and sin(kb-1/2*n) != 0, so that j_n(kb) -> sin(kb-1/2*n)/kb != 0.
else
    % for plane wave cases, b is better to set as a slightly larger sphere that can surround the particle.
    b = 3 * particle_radius; 
    b = 21/4 * pi / fluid_k;
end

% integral radial distance (b) for beam-shape coeffcient, EMPIRICAL: MUST follow  b > a and b < inter_dist BUT not b >> a
if (strcmp(wave_type, 'phase_array_transducer') == 1 || strcmp(wave_type, 'phase_array_transducer2') == 1) && b > inter_dist / 1.1
    error('Integral radial distance (b) for beam-shape coeffcient must not include the source terms!\n');
end
if b < particle_radius * 1.3
    error('Integral radial distance (b) for beam-shape coeffcient must include the particle!\n');
end

% R = 1000 * particle_radius;         % farfield integral distance (R) MUST follow  R >> a, used in "radiation_force_based_NumIntPres.m" 
R = 2000 * b;


%% 
% multi-layer medium with irregular interface

multi_layer = 0;

if multi_layer == 1
    
    NMP = 8;                                % Fourier expansion series truncated to a finite number of terms 'NMP'
    inter_dist_source_layer = 0.05;        % m, vertical distance between the transducer and the irregular layer (the origin of the coordinate system is given in center of the irregular layer)
%     interface_type = 'plane';             % define the interface, used in "define_interface.m"
    interface_type = 'circle';             % define the interface, used in "define_interface.m"
    c_layer = 680;                         % m/s, speed of sound of the irregular layer (self define)
    density_layer = 2.5;                   % kg/m3, density of the irregular layer (self define)
%     c_layer = 340;                         % m/s, speed of sound of the irregular layer (air)
%     density_layer = 1.224;                   % kg/m3, density of the irregular layer (air)
%     c_layer = 1000;                         % m/s, speed of sound of the irregular layer (expanded polystyrene foam)
%     density_layer = 11.5;                   % kg/m3, density of the irregular layer (expanded polystyrene foam)
%     c_layer = 1500;                         % m/s, speed of sound of the irregular layer (water)
%     density_layer = 1000;                   % kg/m3, density of the irregular layer (water)
%     n_layer = c_layer/fluid_c;
    n_layer = 1;            % the acoustic refractive index of the irregular layer (expanded polystyrene foam)
    depth_fun = define_interface(interface_type);
    
    b = (0.8*pi) / (2*pi*freq / c_layer);   % optimation of 'b', making j_n(layer_kb) is local maximum point. Otherwise, if j_n(layer_kb) -> 0, the beam-shape coefficients will tend to infinity and make the result error.
    % integral radial distance (b) for beam-shape coeffcient, EMPIRICAL: MUST follow  b > a and b < inter_dist BUT not b >> a
    if (strcmp(wave_type, 'phase_array_transducer') == 1 || strcmp(wave_type, 'phase_array_transducer2') == 1) && b > inter_dist / 1.1
        error('Integral radial distance (b) for beam-shape coeffcient must not include the source terms!\n');
    end
    if b < particle_radius * 1.3
        error('Integral radial distance (b) for beam-shape coeffcient must include the particle!\n');
    end
    
    R = 2000 * b;
    
end

% ==========================================================
%%
