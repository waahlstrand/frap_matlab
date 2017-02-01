%% Initialization.
clear
clc
close all hidden

%% Load data.
load('simulated_data_zero_noise.mat');
number_of_pixels = size(image_data_post_bleach, 1);

number_of_images_post_bleach = 2;
image_data_post_bleach = image_data_post_bleach(:, :, 1:number_of_images_post_bleach);

%% Add noise.
sigma_noise = 0.01;
image_data_post_bleach = image_data_post_bleach + sigma_noise * randn(size(image_data_post_bleach));

%% Least squares optimization.
% % 
% % % Set parameter bounds.
% % % lb_SI = [1e-10 0 0 0.8]; % (D, kon, koff, mf)
% % % ub_SI = [1e-9 10 10 1]; % (D, kon, koff, mf)
% % % lb_SI = [1e-10, 0, 0, 0.8, 0.25*number_of_pixels, 0.25*number_of_pixels, 5, 0, 0]; % (D, kon, koff, mf, xc, yc, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
% % % ub_SI = [1e-9, 10, 10, 1, 0.75*number_of_pixels, 0.75*number_of_pixels, 0.5*number_of_pixels, 1, 1]; % (D, kon, koff, mf, xc, yc, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
% % lb_SI = [1e-10, 0, 0, 0.8, 5, 0, 0]; % (D, kon, koff, mf, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
% % ub_SI = [1e-9, 10, 10, 1, 0.5*number_of_pixels, 1, 1]; % (D, kon, koff, mf, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
% % 
% % lb = lb_SI;
% % lb(1) = lb(1) / pixel_size^2;
% % 
% % ub = ub_SI;
% % ub(1) = ub(1) / pixel_size^2;
% % 
% % % Initial guess
% % % param_guess_SI = [2e-10 0.5 0.5 1]; % (D, kon, koff, mf)
% % % param_guess_SI = [2.5e-10, 1.0, 1.0, 0.9, 0.5*number_of_pixels, 0.5*number_of_pixels, 30, 0.6, 0.9]; % (D, kon, koff, mf, xc, yc, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
% % param_guess_SI = [2.5e-10, 1.0, 1.0, 0.9, 32, 0.6, 0.9]; % (D, kon, koff, mf, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
% % param_guess = param_guess_SI;
% % param_guess(1) = param_guess(1) / pixel_size^2;
% % 
% % % Optimization routine options.
% % options = optimoptions(@lsqnonlin);
% % options.Algorithm = 'levenberg-marquardt';%'trust-region-reflective';
% % options.Display = 'iter';
% % options.TolFun = 1e-6;
% % options.TolCon = 1e-6;
% % options.TolX = 1e-6;
% % options.UseParallel = true;
% % 
% % % Function handle.
% % % fun = @(param)residual(param(1), param(2), param(3), param(4), xc, yc, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region, delta_t, image_data_post_bleach);
% % % fun = @(param)residual(param(1), param(2), param(3), param(4), param(5), param(6), param(7), param(8), param(9), delta_t, image_data_post_bleach);
% % fun = @(param)residual(param(1), param(2), param(3), param(4), xc, yc, param(5), param(6), param(7), delta_t, image_data_post_bleach);
% %     
% % % [param_hat, ss] = lsqnonlin(fun, param_guess, lb, ub, options);
% % [param_hat, ss] = lsqnonlin(fun, param_guess, [], [], options);

%% Pattern-search optimization.

lb_SI = [2.25e-10, 0.75, 0.75, 0.8, 24, 0.55, 0.8]; % (D, kon, koff, mf, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
ub_SI = [2.75e-10, 1.2, 1.25, 1, 40, 0.70, 1]; % (D, kon, koff, mf, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
lb = lb_SI;
lb(1) = lb(1) / pixel_size^2;
ub = ub_SI;
ub(1) = ub(1) / pixel_size^2;
param_guess_SI = [2.5e-10, 1.0, 1.0, 0.9, 32, 0.6, 0.9]; % (D, kon, koff, mf, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
param_guess = param_guess_SI;
param_guess(1) = param_guess(1) / pixel_size^2;

options = psoptimset(@patternsearch);
options.Algorithm = 'levenberg-marquardt';%'trust-region-reflective';
options.Display = 'iter';
options.CheckGradients = false;
options.SpecifyObjectiveGradient = false;
options.FunctionTolerance = 1e-6;
options.StepTolerance = 1e-6;
options.OptimalityTolerance = 1e-6;

fun = @(param)sum( residual(param(1), param(2), param(3), param(4), xc, yc, param(5), param(6), param(7), delta_t, image_data_post_bleach).^2 );

[param_hat, ss] = patternsearch(fun, param_guess, [], [], [], [], lb, ub, [], options);

