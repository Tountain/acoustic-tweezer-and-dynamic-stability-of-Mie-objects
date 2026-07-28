function [Anm, anm] = beam_shape_coeff_multi_layer(U, n, m)
%%
% This function is used to determine the beam shape coefficient "Anm" based
% on complex incident pressure amplitude field "U", which is calculated by
% sub-function "ASA_amplitude_field.m", and represents a sperical surface
% amplitude field, whose center is consistent with the particle weigh
% center and the radius is 'b'.
%
% Similar to "beam_shape_coeff.m", the numerical integration based on
% "trapz" function. 
%%

parameters;

theta = linspace(0, pi, 100*2);   % Empirical Relation: the interval of 
phi = linspace(0, 2*pi, 100*4);   % "theta" and "phi" should larger than "kb".

conj_Ynm = conj(sHarmonics(n, m, theta, phi));
inner_int = zeros(length(phi), 1);

U = conj(U);            % !!?

%% for multi-layer medium

if multi_layer == 0
    
    error('Try to run a sub-function ''beam_shape_coeff_multi_layer.m'', which is design for multi-layer medium only!\n');
    
elseif multi_layer == 1
    
    for ii = 1:length(phi)
        obj_func = U(ii, :) .* conj_Ynm(ii, :) .* sin(theta); 
        inner_int(ii) = trapz(theta, obj_func);
    end
%     figure(1);plot(real(U(ii, :) .* conj_Ynm(ii, :) .* sin(theta)));
%     figure(2);plot(imag(U(ii, :) .* conj_Ynm(ii, :) .* sin(theta)));
%     figure(3);plot(real(inner_int));
%     figure(4);plot(imag(inner_int));
    integral = trapz(phi, inner_int);

    layer_kb = 2*pi*freq / c_layer * b;         % wavenumber should corresponse to the fluid that involved the objects
    anm = integral / ((2*n + 1) * 1i^n * (sBessel(n, layer_kb, 1)+eps));   % G.T.S@2011@Off-axis, Eq.(8)
    Anm = anm * (2*n + 1) * 1i^n;                                    % more general explaination
     
end

%%
