%% Initialization.
clear
clc
close all hidden

random_seed = round(sum(1e6 * clock()));
random_stream = RandStream('mt19937ar', 'Seed', random_seed);
RandStream.setGlobalStream(random_stream);

%% Experimental and simulation parameters.

exp_sim_param = struct();

exp_sim_param.pixel_size = 7.5e-07; % m
exp_sim_param.number_of_pixels = 256;

exp_sim_param.number_of_prebleach_frames = 10;
exp_sim_param.number_of_bleach_frames = 2;
exp_sim_param.number_of_postbleach_frames = 50;
exp_sim_param.delta_t = 0.2; % s

exp_sim_param.number_of_pad_pixels = 128;

exp_sim_param.bleach_region.shape = "circle";
exp_sim_param.bleach_region.x = 128; % pixels
exp_sim_param.bleach_region.y = 128; % pixels 
exp_sim_param.bleach_region.r = 15e-6 / exp_sim_param.pixel_size; % pixels
% exp_sim_param.bleach_region.lx = 20e-6 / exp_sim_param.pixel_size; % pixels
% exp_sim_param.bleach_region.ly = 20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.upsampling_factor = 16;

%% System parameters.

D_SI = 5e-10; % m^2/s
D = D_SI / exp_sim_param.pixel_size^2; % pixels^2 / s
mobile_fraction = 1.0; % dimensionless
C0 = 1.0; % a.u. original concentration
alpha = 0.6; % a.u.  bleach factor
beta = 1.0; % a.u. imaging bleach factor

sys_param = [D, mobile_fraction, C0, alpha, beta];

%% Generate data.

[C_prebleach, C_postbleach] = signal_d(sys_param, exp_sim_param);
sigma = 0.1;
C_prebleach = C_prebleach + sigma * randn(size(C_prebleach));
C_postbleach = C_postbleach + sigma * randn(size(C_postbleach));

%% Fit parameters.

exp_sim_param.number_of_bleach_frames = 2;

fit_param = struct();

fit_param.mode = "pixel";%"recovery-curve"; % "pixel"
fit_param.use_parallel = false;
fit_param.number_of_fits = 1;
fit_param.guess = [D, mobile_fraction, C0, alpha, beta];
fit_param.lower_bound = [0.5 * D, 0, 0, 0, 0];
fit_param.upper_bound = [2 * D, 1, 2 * C0, 1, 1];

%% Fit.

[sys_param_hat, ss] = fit_d(C_prebleach, C_postbleach, exp_sim_param, fit_param);

%% Plot results.

[C_prebleach_model, C_postbleach_model] = signal_d(sys_param_hat, exp_sim_param);
C_model = cat(3, C_prebleach_model, C_postbleach_model);
C = cat(3, C_prebleach, C_postbleach);
R = C_model(:) - C(:);

[rc_prebleach_model, rc_postbleach_model] = recovery_curve(C_prebleach_model, C_postbleach_model, exp_sim_param);
rc_model = [rc_prebleach_model ; rc_postbleach_model];
[rc_prebleach, rc_postbleach] = recovery_curve(C_prebleach, C_postbleach, exp_sim_param);
rc = [rc_prebleach ; rc_postbleach];

figure
hold on
plot(rc, 'k.')
plot(rc_model, 'r')


figure, imagesc([reshape(C, [exp_sim_param.number_of_pixels, exp_sim_param.number_of_pixels * (exp_sim_param.number_of_prebleach_frames + exp_sim_param.number_of_postbleach_frames)]) ; ...
                 reshape(C_model, [exp_sim_param.number_of_pixels, exp_sim_param.number_of_pixels * (exp_sim_param.number_of_prebleach_frames + exp_sim_param.number_of_postbleach_frames)])]);
figure, imagesc(reshape(C_model - C, [exp_sim_param.number_of_pixels, exp_sim_param.number_of_pixels * (exp_sim_param.number_of_prebleach_frames + exp_sim_param.number_of_postbleach_frames)]));





% D = param_hat(1) * pixel_size^2
% mf = param_hat(2);
% Ib = param_hat(3);
% Iu = param_hat(4);
% 
% 
% %% Show results.
% 
% model = signal_d( ...
%     param_hat(1), ...
%     param_hat(2), ...
%     param_hat(3), ...
%     param_hat(4), ...
%     param_bleach, ...
%     delta_t, ...
%     number_of_pixels, ...
%     number_of_images, ...
%     number_of_pad_pixels);
%                 
% figure, imagesc([reshape(data, [number_of_pixels, number_of_pixels * number_of_images]) ; reshape(model, [number_of_pixels, number_of_pixels * number_of_images])])
% figure, imagesc(reshape(data - model, [number_of_pixels, number_of_pixels * number_of_images]))
% 
% rc_data = zeros(1, number_of_images_prebleach + number_of_images);
% for current_image = 1:number_of_images_prebleach
%     slice = data_prebleach(:, :, current_image);
%     rc_data(current_image) = mean(slice(ind));
% end
% for current_image = 1:number_of_images
%     slice = data(:, :, current_image);
%     rc_data(number_of_images_prebleach + current_image) = mean(slice(ind));
% end
% 
% rc_model = zeros(1, number_of_images_prebleach + number_of_images);
% for current_image = 1:number_of_images_prebleach
%     rc_model(current_image) = Iu;
% end
% for current_image = 1:number_of_images
%     slice = model(:, :, current_image);
%     rc_model(number_of_images_prebleach + current_image) = mean(slice(ind));
% end
% 
% figure, hold on
% plot([(-number_of_images_prebleach:-1)*delta_t (1:number_of_images)*delta_t], rc_data, 'ro');
% plot([(-number_of_images_prebleach:-1)*delta_t (1:number_of_images)*delta_t], rc_model, 'k-');
% 
% 
% 
