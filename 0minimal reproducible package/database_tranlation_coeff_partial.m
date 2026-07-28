function [db_filename_ts_partial] = database_tranlation_coeff_partial(N, k, r, theta, phi)
%%
% This function is used to save the partial derivative of separation matrix
% of first kind for vector position '-(r, theta, phi)'.
%%

parameters_names;

if exist([db_filename_ts_partial, '.mat']) ~= 0        % if already exist the database, 
    return;                                     % then do not create again for saving time.
end

%%

Snmvu_1 = zeros(N+1, 2*N+1, N+1, 2*N+1);
Snmvu_partial_r = zeros(N+1, 2*N+1, N+1, 2*N+1);
Snmvu_partial_theta = zeros(N+1, 2*N+1, N+1, 2*N+1);
Snmvu_partial_phi = zeros(N+1, 2*N+1, N+1, 2*N+1);
for nn = 0:N
    indices_1 = nn+1;
    parfor mm = -nn:nn
        temp_Snmvu_1 = zeros(N+1, 2*N+1);
        temp_Snmvu_partial_r = zeros(N+1, 2*N+1);
        temp_Snmvu_partial_theta = zeros(N+1, 2*N+1);
        temp_Snmvu_partial_phi = zeros(N+1, 2*N+1);
        for nu = 0:N
            for mu = -nu:nu
                indices_3 = nu+1;
                indices_4 = nu+mu+1;
                [temp_Snmvu_1(indices_3, indices_4), ~] = Snmvu_coeff(nn, mm, nu, mu, k*r, theta, phi);
                [temp_Snmvu_partial_r(indices_3, indices_4), temp_Snmvu_partial_theta(indices_3, indices_4), temp_Snmvu_partial_phi(indices_3, indices_4)] = ...
                            Snmvu_coeff_partial(nn, mm, nu, mu, k, r, theta, phi);       
            end
        end
        Snmvu_1(indices_1, indices_1 + mm, :, :) = temp_Snmvu_1;
        Snmvu_partial_r(indices_1, indices_1 + mm, :, :) = temp_Snmvu_partial_r;
        Snmvu_partial_theta(indices_1, indices_1 + mm, :, :) = temp_Snmvu_partial_theta;
        Snmvu_partial_phi(indices_1, indices_1 + mm, :, :) = temp_Snmvu_partial_phi;
    end
    fprintf('Partial Translation Coefficients Database Preparing %d%% \n', ...
            round(100*nn/N));
end

%% save the database by meaningful name

save([db_filename_ts_partial, '.mat'], 'Snmvu_1', 'Snmvu_partial_r', 'Snmvu_partial_theta', 'Snmvu_partial_phi');

%%