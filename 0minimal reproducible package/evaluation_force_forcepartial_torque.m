function evaluation_force_forcepartial_torque(A_delay, phi_delay, deviations)
%%
% for jj = 1:2:length(A_delay_dataset)
%     evaluation_force_forcepartial_torque(A_delay_dataset{jj}, phi_delay_dataset{jj}, posture{jj});
% end
%%

figure_FontSize = 20;
line_width = 2;

% A_delay, phi_delay is from 'sphere_stable_trap-10mm.mat'
replace_Amp_delay_phi_delay(A_delay, phi_delay);

sample_number = 21;
delta_x = linspace(-1, 1, sample_number) * 0.01;
delta_y = linspace(-1, 1, sample_number) * 0.01;
delta_z = linspace(-1, 1, sample_number) * 0.01;

FX = [];
FPX = [];
for ii = 1:sample_number    % X-axis
    
    position = [delta_x(ii) + deviations(1), deviations(2), deviations(3)];
    rotation = [0, 0, 0];
    replace_position_theta(position, rotation);     % object POSTURE revise.
%     % renew the names for 'dir_file_forces_trans', 'dir_file_forces_trans_partial', and 'dir_file_torques_trans'.
%     parameters_names;
%     
%     if exist([dir_file_forces_trans, 'A.mat'], 'file') ~=0 && exist([dir_file_forces_trans_partial, 'A.mat'], 'file') ~=0 && exist([dir_file_torques_trans, 'A.mat'], 'file') ~=0
%         load([dir_file_forces_trans, 'A.mat'], 'Frad_x', 'Frad_y', 'Frad_z');
%         load([dir_file_forces_trans_partial, 'A.mat'], 'Frad_x_partial', 'Frad_y_partial', 'Frad_z_partial');
%         load([dir_file_torques_trans, 'A.mat'], 'Torque_x', 'Torque_y', 'Torque_z');
%     else
%         [Frad_x,Frad_y,Frad_z,~] = radiation_force_based_Analyses();
%         [Frad_x_partial,Frad_y_partial,Frad_z_partial,~] = radiation_force_based_Analyses_partial();
%         [Torque_x,Torque_y,Torque_z,~] = radiation_torque_based_Analyses();
%     end
    [Frad_x,Frad_y,Frad_z,~] = radiation_force_based_Analyses();
    [Frad_x_partial,Frad_y_partial,Frad_z_partial,~] = radiation_force_based_Analyses_partial();
    [Torque_x,Torque_y,Torque_z,~] = radiation_torque_based_Analyses();
    
    FX = [FX Frad_x];
    FPX = [FPX Frad_x_partial];

end

figure(1); box on; hold on;
plot(delta_x, FX, '-', 'linewidth', line_width);
xlabel('\rm{\fontname{Times new roman}\rm{{\it{x}}-axis displacement }\rm{[mm]}}');
ylabel('\rm{\fontname{Times new roman}{\it{F_x}}\rm{ [{\mu}N]}}');
set(gca, 'FontName', 'Times new roman');
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);


FY = [];
FPY = []; 
for ii = 1:sample_number    % Y-axis
    
    position = [deviations(1), delta_y(ii) + deviations(2), deviations(3)];
    rotation = [0, 0, 0];
    replace_position_theta(position, rotation);     % object POSTURE revise.
%     % renew the names for 'dir_file_forces_trans', 'dir_file_forces_trans_partial', and 'dir_file_torques_trans'.
%     parameters_names; 
%     
%     if exist([dir_file_forces_trans, 'A.mat'], 'file') ~=0 && exist([dir_file_forces_trans_partial, 'A.mat'], 'file') ~=0 && exist([dir_file_torques_trans, 'A.mat'], 'file') ~=0
%         load([dir_file_forces_trans, 'A.mat'], 'Frad_x', 'Frad_y', 'Frad_z');
%         load([dir_file_forces_trans_partial, 'A.mat'], 'Frad_x_partial', 'Frad_y_partial', 'Frad_z_partial');
%         load([dir_file_torques_trans, 'A.mat'], 'Torque_x', 'Torque_y', 'Torque_z');
%     else
%         [Frad_x,Frad_y,Frad_z,~] = radiation_force_based_Analyses();
%         [Frad_x_partial,Frad_y_partial,Frad_z_partial,~] = radiation_force_based_Analyses_partial();
%         [Torque_x,Torque_y,Torque_z,~] = radiation_torque_based_Analyses();
%     end
    [Frad_x,Frad_y,Frad_z,~] = radiation_force_based_Analyses();
    [Frad_x_partial,Frad_y_partial,Frad_z_partial,~] = radiation_force_based_Analyses_partial();
    [Torque_x,Torque_y,Torque_z,~] = radiation_torque_based_Analyses();
    
    FY = [FY Frad_y];
    FPY = [FPY Frad_y_partial];
    
end

figure(2); box on; hold on;
plot(delta_y, FY, '-', 'linewidth', line_width);
xlabel('\rm{\fontname{Times new roman}\rm{{\it{y}}-axis displacement }\rm{[mm]}}');
ylabel('\rm{\fontname{Times new roman}{\it{F_y}}\rm{ [{\mu}N]}}');
set(gca, 'FontName', 'Times new roman');
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);


FZ = [];
FPZ = [];
for ii = 1:sample_number    % Z-axis
    
    position = [deviations(1), deviations(2), delta_z(ii) + deviations(3)];
    rotation = [0, 0, 0];
    replace_position_theta(position, rotation);     % object POSTURE revise.
%     % renew the names for 'dir_file_forces_trans', 'dir_file_forces_trans_partial', and 'dir_file_torques_trans'.
%     parameters_names; 
%     
%     if exist([dir_file_forces_trans, 'A.mat'], 'file') ~=0 && exist([dir_file_forces_trans_partial, 'A.mat'], 'file') ~=0 && exist([dir_file_torques_trans, 'A.mat'], 'file') ~=0
%         load([dir_file_forces_trans, 'A.mat'], 'Frad_x', 'Frad_y', 'Frad_z');
%         load([dir_file_forces_trans_partial, 'A.mat'], 'Frad_x_partial', 'Frad_y_partial', 'Frad_z_partial');
%         load([dir_file_torques_trans, 'A.mat'], 'Torque_x', 'Torque_y', 'Torque_z');
%     else
%         [Frad_x,Frad_y,Frad_z,~] = radiation_force_based_Analyses();
%         [Frad_x_partial,Frad_y_partial,Frad_z_partial,~] = radiation_force_based_Analyses_partial();
%         [Torque_x,Torque_y,Torque_z,~] = radiation_torque_based_Analyses();
%     end
    [Frad_x,Frad_y,Frad_z,~] = radiation_force_based_Analyses();
    [Frad_x_partial,Frad_y_partial,Frad_z_partial,~] = radiation_force_based_Analyses_partial();
    [Torque_x,Torque_y,Torque_z,~] = radiation_torque_based_Analyses();

    FZ = [FZ Frad_z];
    FPZ = [FPZ Frad_z_partial];
    
end

figure(3); box on; hold on;
plot(delta_z, FZ, '-', 'linewidth', line_width);
xlabel('\rm{\fontname{Times new roman}\rm{{\it{z}}-axis displacement }\rm{[mm]}}');
ylabel('\rm{\fontname{Times new roman}{\it{F_z}}\rm{ [{\mu}N]}}');
set(gca, 'FontName', 'Times new roman');
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);

%%

save(['Forces&SpatialDeviation_X', num2str(deviations(1)), '_Y', num2str(deviations(2)), '_Z', num2str(deviations(3)), '.mat'], ...
        'delta_x', 'FX', 'FPX', 'delta_y', 'FY', 'FPY', 'delta_z', 'FZ', 'FPZ');
    
%%