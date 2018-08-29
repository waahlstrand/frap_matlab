%% Initialization.

clear
clc
close all hidden

%% Experimental and simulation parameters.

exp_sim_param = struct();

exp_sim_param.pixel_size = 7.5e-07; % m
exp_sim_param.number_of_pixels = 256;

exp_sim_param.number_of_prebleach_frames = 0;
exp_sim_param.number_of_bleach_frames = 1;
exp_sim_param.number_of_postbleach_frames = 100;
exp_sim_param.delta_t = 0.2; % s

exp_sim_param.number_of_pad_pixels = 128;

exp_sim_param.bleach_region.shape = "rectangle";%"circle";
exp_sim_param.bleach_region.x = 128; % pixels
exp_sim_param.bleach_region.y = 128; % pixels 
exp_sim_param.bleach_region.r = 15e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.lx = 20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.ly = 20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.upsampling_factor = 16;

%% System parameters.

D_SI = 5e-11; % m^2/s
D = D_SI / exp_sim_param.pixel_size^2; % pixels^2 / s
mobile_fraction = 1.0; % dimensionless
C0 = 1.0; % a.u. original concentration
alpha = 0.6; % a.u.  bleach factor
beta = 1.0; % a.u. imaging bleach factor

sys_param = [D, mobile_fraction, C0, alpha, beta];

%% Simulate.

tic
[C_prebleach, C_postbleach] = signal_d(sys_param, exp_sim_param);
toc

[rc_prebleach, rc_postbleach] = recovery_curve(C_prebleach, C_postbleach, exp_sim_param);

plot([rc_prebleach ; rc_postbleach])


%% Plot concentration.

% for current_frame = 1:exp_sim_param.number_of_prebleach_frames
%     figure, imagesc([C_prebleach(:, :, current_frame), C_prebleach_sim(:, :, current_frame), C_prebleach(:, :, current_frame) - C_prebleach_sim(:, :, current_frame)]), axis 'equal'
% end
% % for current_frame = 1:exp_sim_param.number_of_postbleach_frames
% %     figure, imagesc([C_postbleach(:, :, current_frame), C_postbleach_sim(:, :, current_frame), C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)]), axis 'equal'
% % end
% for current_frame = 1:exp_sim_param.number_of_prebleach_frames
%     figure, imagesc(C_prebleach(:, :, current_frame) - C_prebleach_sim(:, :, current_frame)), axis 'equal'
% end
% for current_frame = 1:exp_sim_param.number_of_postbleach_frames
%     figure, imagesc(C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)), axis 'equal'
% end
% max(abs(C_prebleach(:)-C_prebleach_sim(:)))
% max(abs(C_postbleach(:)-C_postbleach_sim(:)))

% figure, hold on, plot(C_postbleach(:, 128), 'k.-'), plot(C_postbleach_sim(:, 128), 'r.-')
% figure, plot(C_postbleach(128, :) - C_postbleach_sim(128, :))
