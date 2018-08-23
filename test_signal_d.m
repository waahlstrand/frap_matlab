%% Initialization.

clear
clc
close all hidden

%% Experimental and simulation parameters.

exp_sim_param = struct();

exp_sim_param.pixel_size = 7.5e-07; % m
exp_sim_param.number_of_pixels = 256;

exp_sim_param.number_of_prebleach_frames = 20;
exp_sim_param.number_of_bleach_frames = 4;
exp_sim_param.number_of_postbleach_frames = 100;
exp_sim_param.delta_t = 0.2; % s

exp_sim_param.number_of_pad_pixels = 128;

exp_sim_param.bleach_region.shape = 'rectangular';
exp_sim_param.bleach_region.x = 128; % pixels
exp_sim_param.bleach_region.y = 128; % pixels 
exp_sim_param.bleach_region.r = 15e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.lx = 20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.ly = 20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.upsampling_factor  = 3;

%% System parameters.

D_SI = 5e-11; % m^2/s
D = D_SI / exp_sim_param.pixel_size^2; % pixels^2 / s
mobile_fraction = 1.0; % dimensionless
C0 = 1.0; % a.u. original concentration
alpha = 0.9; % a.u.  bleach factor

sys_param = [D, mobile_fraction, C0, alpha];

%% Simulate.

tic
C = signal_d(sys_param, exp_sim_param);
toc

sigma_noise = 0;
C = C + sigma_noise * randn(size(C));

%% Plot solution.
% figure
% hold on
% imagesc(reshape(C, [number_of_pixels, number_of_pixels * number_of_images]))
% % axis 'equal'
% axis([0 number_of_images*number_of_pixels 0 number_of_pixels])
% axis off
% hold off
% 
% %% Plot recovery curve.
% [X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
% X = X - 0.5;
% Y = Y - 0.5;
% 
% if numel(param_bleach) == 3 % Circular.
%     ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
% else % Rectangular.
%     ind = find( X >= x_bleach - 0.5 * lx_bleach & X <= x_bleach + 0.5 * lx_bleach & Y >= y_bleach - 0.5 * ly_bleach & Y <= y_bleach + 0.5 * ly_bleach );
% end
% ind = ind(:);
% 
% recovery_curve = zeros(1, number_of_images);
% for current_image = 1:number_of_images
%     slice = C(:, :, current_image);
%     recovery_curve(current_image) = mean(slice(ind));
% end
% figure
% plot(delta_t:delta_t:number_of_images*delta_t, recovery_curve)
