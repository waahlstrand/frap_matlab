%% Initialization.
clear
clc
close all hidden

%% Load data.
load('simulated_data_zero_noise.mat');
number_of_pixels = size(image_data_post_bleach, 1);

number_of_post_bleach_images = 2;
image_data_post_bleach = image_data_post_bleach(:, :, 1:number_of_post_bleach_images);

delta_t = 0.250; % s.
number_of_time_points_fine_per_coarse = 500; % dimensionless
number_of_pad_pixels = 128;
x_bleach = number_of_pixels / 2; % pixels
y_bleach = number_of_pixels / 2; % pixels

%% Add noise.
sigma_noise = 0.01;
image_data_post_bleach = image_data_post_bleach + sigma_noise * randn(size(image_data_post_bleach));

%% Parameter estimation pre-work.

% Set parameter bounds.
lb_SI = [1e-10, 0, 0, 0.8, 5, 0, 0]; % (D, kon, koff, mobile_fraction, r_bleach, intensity_inside_bleach_region, intensity_outside_bleach_region)
ub_SI = [1e-9, 10, 10, 1, 0.5*number_of_pixels, 1, 1]; % (D, kon, koff, mobile_fraction, r_bleach, intensity_inside_bleach_region, intensity_outside_bleach_region)

lb = lb_SI;
lb(1) = lb(1) / pixel_size^2;
ub = ub_SI;
ub(1) = ub(1) / pixel_size^2;

% Initial guess
param_guess_SI = [2.5e-10, 1.0, 1.0, 0.9, 32, 0.6, 0.9]; % (D, kon, koff, mf, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
param_guess = param_guess_SI;
param_guess(1) = param_guess(1) / pixel_size^2;

%% Least-squares optimization.

% fun = @(param)residual_diffusion_and_binding(   param(1), ...
%                                                 param(2), ...
%                                                 param(3), ...
%                                                 param(4), ...
%                                                 x_bleach, ...
%                                                 y_bleach, ...
%                                                 param(5), ...
%                                                 param(6), ...
%                                                 param(7), ...
%                                                 delta_t, ...
%                                                 number_of_time_points_fine_per_coarse, ...
%                                                 number_of_pad_pixels, ...
%                                                 image_data_post_bleach);
% options = optimoptions(@lsqnonlin);
% options.Algorithm = 'levenberg-marquardt';%'trust-region-reflective';
% options.Display = 'iter';
% options.TolFun = 1e-6;
% options.TolCon = 1e-6;
% options.TolX = 1e-6;
% options.UseParallel = true;
% 
% [param_hat, ss] = lsqnonlin(fun, param_guess, lb, ub, options);

%% Pattern-search optimization.
fun = @(param)sum ( residual_diffusion_and_binding( param(1), ...
                                                    param(2), ...
                                                    param(3), ...
                                                    param(4), ...
                                                    x_bleach, ...
                                                    y_bleach, ...
                                                    param(5), ...
                                                    param(6), ...
                                                    param(7), ...
                                                    delta_t, ...
                                                    number_of_time_points_fine_per_coarse, ...
                                                    number_of_pad_pixels, ...
                                                    image_data_post_bleach).^2 );
options = psoptimset(@patternsearch);
options.Display = 'iter';
options.UseParallel = true;

[param_hat, ss] = patternsearch(fun, param_guess, [], [], [], [], lb, ub, [], options);

