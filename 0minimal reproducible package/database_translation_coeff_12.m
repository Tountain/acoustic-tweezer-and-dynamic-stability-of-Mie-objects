function [db_filename_Snmvu] = database_translation_coeff_12()
%%
% build the database for "translation coeffcient type 1" and
% "translation coeffcient type 2" to avoid repeating call function
% "Snmvu_coeff", "gaunt_value" and "wigner3j", which will cost tons of
% time.
% 
% database will automatically be created if current folder has not a
% corresponding database. (check the database file exists or not by the
% filename)
% database maximum size: db_size_nn = 30 (seeing "parameters_names.m");
%%

%% do not need to create this database!!

parameters_names;

if exist([db_filename_Snmvu, '.mat'], 'file') ~= 0        % if already exist the database, 
    return;                                               % then do not create again for saving time.
end

particles_Cartesian_data;           % input particles' data

% "db_Snmvu_1_coeff" and "db_Snmvu_2_coeff" both need 5 indexes,
db_Snmvu_1_coeff = zeros();
db_Snmvu_2_coeff = zeros();

%% do not need to create this database!!




