%% Initialization.
clear
clc
close all hidden

%% Load data.
file_path = '../../../data/data_binding_early_wood_6A/FRAP_002.mat';
raw_data = load(file_path);
image_data_pre_bleach = raw_data.imdata{1};
image_data_post_bleach = raw_data.imdata{3};
clear file_path raw_data

bit_depth = 16;
pixel_size = 7.5980e-07; % m.
delta_t = 0.2650; % s.

number_of_pixels = size(image_data_post_bleach, 1);
number_of_post_bleach_images = 10;

%% Extract desired numbers of images/frames to include in analysis.
image_data_post_bleach = image_data_post_bleach(:, :, 1:number_of_post_bleach_images);

%% Convert image data to double and rescale to [0, 1] range.
image_data_pre_bleach = double(image_data_pre_bleach);
image_data_pre_bleach = image_data_pre_bleach / (2^bit_depth - 1);

image_data_post_bleach = double(image_data_post_bleach);
image_data_post_bleach = image_data_post_bleach / (2^bit_depth - 1);

%% Background subtraction.
if do_background_subtraction
    image_data_post_bleach = subtract_background(image_data_pre_bleach, image_data_post_bleach);
end

%% Parameter estimation pre-work.

% Set parameter bounds.
lb_SI = [1e-10, 0, 0];%, 0.8, 5, 0, 0]; % (D, kon, koff, mobile_fraction, r_bleach, intensity_inside_bleach_region, intensity_outside_bleach_region)
ub_SI = [1e-9, 10, 10];%, 1, 0.5*number_of_pixels, 1, 1]; % (D, kon, koff, mobile_fraction, r_bleach, intensity_inside_bleach_region, intensity_outside_bleach_region)

lb = lb_SI;
lb(1) = lb(1) / pixel_size^2;
ub = ub_SI;
ub(1) = ub(1) / pixel_size^2;

% Initial guess
param_guess_SI = [2.5e-10, 0.5, 0.5];%, 0.9, 32, 0.6, 0.9]; % (D, kon, koff, mf, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region)
param_guess = param_guess_SI;
param_guess(1) = param_guess(1) / pixel_size^2;

%% Least-squares optimization.

fun = @(param)residual_diffusion_and_binding(   param(1), ...
                                                param(2), ...
                                                param(3), ...
                                                mobile_fraction, ...
                                                x_bleach, ...
                                                y_bleach, ...
                                                r_bleach, ...
                                                intensity_inside_bleach_region, ...
                                                intensity_outside_bleach_region, ...
                                                delta_t, ...
                                                number_of_time_points_fine_per_coarse, ...
                                                number_of_pad_pixels, ...
                                                image_data_post_bleach);
options = optimoptions(@lsqnonlin);
options.Algorithm = 'trust-region-reflective';
options.Display = 'iter';
options.FunctionTolerance = 1e-6;
options.OptimalityTolerance = 1e-6;
options.StepTolerance = 1e-6;
% options.UseParallel = true;

[param_hat, ss] = lsqnonlin(fun, param_guess, lb, ub, options);

%% Global optimization.
% fun = @(param)sum ( residual_diffusion_and_binding( param(1), ...
%                                                     param(2), ...
%                                                     param(3), ...
%                                                     param(4), ...
%                                                     x_bleach, ...
%                                                     y_bleach, ...
%                                                     param(5), ...
%                                                     param(6), ...
%                                                     param(7), ...
%                                                     delta_t, ...
%                                                     number_of_time_points_fine_per_coarse, ...
%                                                     number_of_pad_pixels, ...
%                                                     image_data_post_bleach).^2 );
% fun = @(param)sum ( residual_diffusion_and_binding( param(1), ...
%                                                     param(2), ...
%                                                     param(3), ...
%                                                     mobile_fraction, ...
%                                                     x_bleach, ...
%                                                     y_bleach, ...
%                                                     r_bleach, ...
%                                                     intensity_inside_bleach_region, ...
%                                                     intensity_outside_bleach_region, ...
%                                                     delta_t, ...
%                                                     number_of_time_points_fine_per_coarse, ...
%                                                     number_of_pad_pixels, ...
%                                                     image_data_post_bleach).^2 );
% delete(gcp('nocreate'))
% c = parcluster('local');
% c.NumWorkers = 8;
% parpool(c, c.NumWorkers);
% options = psoptimset(@patternsearch);
% options.Display = 'iter';
% options.UseParallel = true;
% [param_hat, ss] = patternsearch(fun, param_guess, [], [], [], [], lb, ub, [], options);

% options = saoptimset('simulannealbnd');
% options.Display = 'iter';
% options.ObjectiveLimit = 0;
% options.DisplayInterval = 1;
% options.TimeLimit = Inf;
% [param_hat, ss] = simulannealbnd(fun, param_guess, lb, ub, options);





