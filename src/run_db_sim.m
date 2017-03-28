%% Initialization.
clear
clc
close all hidden

addpath('signal_pde');

%% Init reandom stream.
random_seed = sum( 1e6 * clock() );
random_stream = RandStream('mt19937ar', 'Seed', random_seed);
RandStream.setGlobalStream(random_stream);

%% Simulate data.
D_SI = 2.5e-10; % m^2/s
pixel_size = 7.598e-07; % m
D = D_SI / pixel_size^2; % pixels^2 / s
k_on = 0.5;%0.05; % 1/s
k_off = 1.0;%0.01; % 1/s
mf = 0.9; % dimensionless
Ib = 0.6; % a.u.
Iu = 1.0; % a.u.
x_bleach = 128; % pixels
y_bleach = 128; % pixels
r_bleach = 32; % pixels

delta_t = 0.2650; % s
number_of_pixels = 256;
number_of_images = 20;
number_of_pad_pixels = 128;

data = signal_db(   D, ...
                    k_on, ...
                    k_off, ...
                    mf, ...
                    Ib, ...
                    Iu, ...
                    x_bleach, ...
                    y_bleach, ...
                    r_bleach, ...
                    delta_t, ...
                    number_of_pixels, ...
                    number_of_images, ...
                    number_of_pad_pixels);

sigma_noise = 0.025;
data = data + sigma_noise * randn(size(data));

%% Estimate parameters.

lb_D_SI = 1e-12;
ub_D_SI = 1e-8;
lb_D = lb_D_SI / pixel_size^2;
ub_D = ub_D_SI / pixel_size^2;

lb_k_on = 0;
ub_k_on = 10;

lb_k_off = 0;
ub_k_off = 10;

lb_

% Set parameter bounds.
lb = [0.0, 0.0, 1.0]; % mf, Ib, Iu
ub = [1.0, 1.0, 1.0];

% Initial guess for first estimation.
param_hat_1 = [mf, Ib, Iu];

% Set parameter bounds for second estimation.
lb_2_SI = [1e-12, 0, 0];
ub_2_SI = [1e-8, 10, 10];

lb_2 = lb_2_SI;
lb_2(1) = lb_2(1) / pixel_size^2;
ub_2 = ub_2_SI;
ub_2(1) = ub_2(1) / pixel_size^2;

% Initial guess for second estimation.
param_hat_2_SI = [D_SI, k_on, k_off];
param_hat_2 = param_hat_2_SI;
param_hat_2(1) = param_hat_2(1) / pixel_size^2;

% Optimization.
options_1 = optimoptions(@lsqnonlin);
options_1.Algorithm = 'trust-region-reflective';
options_1.Display = 'iter';
options_1.FunctionTolerance = 1e-6;
options_1.OptimalityTolerance = 1e-6;
options_1.StepTolerance = 1e-6;
options_1.CheckGradients = false;
options_1.SpecifyObjectiveGradient = true;

options_2 = optimoptions(@lsqnonlin);
options_2.Algorithm = 'trust-region-reflective';
options_2.Display = 'iter';
options_2.FunctionTolerance = 1e-6;
options_2.OptimalityTolerance = 1e-6;
options_2.StepTolerance = 1e-6;
options_2.UseParallel = true;

param_hat = [];
ss = [];

param_hat = [param_hat ; [param_hat_2 param_hat_1]];
ss = [ss ; Inf];

is_converged = false;
while ~is_converged
    [image_data_post_bleach_model_unscaled, initial_condition_model_unscaled] = signal_diffusion_and_binding(   param_hat_2(1), ...
                                                                                                                param_hat_2(2), ...
                                                                                                                param_hat_2(3), ...
                                                                                                                1.0, ...
                                                                                                                x_bleach, ...
                                                                                                                y_bleach, ...
                                                                                                                r_bleach, ...
                                                                                                                0.5, ...
                                                                                                                1.0, ...
                                                                                                                delta_t, ...
                                                                                                                number_of_time_points_fine_per_coarse, ...
                                                                                                                number_of_pixels, ...
                                                                                                                number_of_post_bleach_images, ...
                                                                                                                number_of_pad_pixels);
    fun_1 = @(param)residual_diffusion_and_binding_partial( param(1), ...
                                                            param(2), ...
                                                            param(3), ...
                                                            data, ...
                                                            image_data_post_bleach_model_unscaled, ...
                                                            initial_condition_model_unscaled);

    [param_hat_1, ss_1] = lsqnonlin(fun_1, param_hat_1, lb_1, ub_1, options_1);
    
    disp([param_hat_2 param_hat_1])
    
    fun_2 = @(param)residual_diffusion_and_binding_full(param(1), ...
                                                        param(2), ...
                                                        param(3), ...
                                                        param_hat_1(1), ...
                                                        x_bleach, ...
                                                        y_bleach, ...
                                                        r_bleach, ...
                                                        param_hat_1(2), ...
                                                        param_hat_1(3), ...
                                                        delta_t, ...
                                                        number_of_time_points_fine_per_coarse, ...
                                                        number_of_pad_pixels, ...
                                                        data);

    [param_hat_2, ss_2] = lsqnonlin(fun_2, param_hat_2, lb_2, ub_2, options_2);
    
    param_hat = [param_hat ; [param_hat_2 param_hat_1]];
    ss = [ss ; ss_2];
    
    disp([param_hat_2 param_hat_1])
    
    if size(param_hat, 1) >= 2
        max_jump = max(abs(param_hat(end, :) - param_hat(end - 1, :)));
        disp(max_jump)
        if max_jump < 1e-3
            is_converged = true;
        end
    end
            
end

param_hat_pixelbased = param_hat;
ss_pixelbased = ss;

%% Recovery curve-based fitting.

% Initial guess for first estimation.
param_hat_1 = [mf, Ib, Iu];

% Initial guess for second estimation.
param_hat_2_SI = [D_SI, k_on, k_off];
param_hat_2 = param_hat_2_SI;
param_hat_2(1) = param_hat_2(1) / pixel_size^2;

% Optimization.
options_1 = optimoptions(@lsqnonlin);
options_1.Algorithm = 'trust-region-reflective';
options_1.Display = 'iter';
options_1.FunctionTolerance = 1e-6;
options_1.OptimalityTolerance = 1e-6;
options_1.StepTolerance = 1e-6;
options_1.CheckGradients = false;
options_1.SpecifyObjectiveGradient = false; % SIC!!!

options_2 = optimoptions(@lsqnonlin);
options_2.Algorithm = 'trust-region-reflective';
options_2.Display = 'iter';
options_2.FunctionTolerance = 1e-6;
options_2.OptimalityTolerance = 1e-6;
options_2.StepTolerance = 1e-6;
options_2.UseParallel = true;

param_hat = [];
ss = [];

param_hat = [param_hat ; [param_hat_2 param_hat_1]];
ss = [ss ; Inf];

is_converged = false;
while ~is_converged
    [image_data_post_bleach_model_unscaled, initial_condition_model_unscaled] = signal_diffusion_and_binding(   param_hat_2(1), ...
                                                                                                                param_hat_2(2), ...
                                                                                                                param_hat_2(3), ...
                                                                                                                1.0, ...
                                                                                                                x_bleach, ...
                                                                                                                y_bleach, ...
                                                                                                                r_bleach, ...
                                                                                                                0.5, ...
                                                                                                                1.0, ...
                                                                                                                delta_t, ...
                                                                                                                number_of_time_points_fine_per_coarse, ...
                                                                                                                number_of_pixels, ...
                                                                                                                number_of_post_bleach_images, ...
                                                                                                                number_of_pad_pixels);
    fun_1 = @(param)residual_recoverycurve_diffusion_and_binding_partial(   param(1), ...
                                                                            param(2), ...
                                                                            param(3), ...
                                                                            x_bleach, ...
                                                                            y_bleach, ...
                                                                            r_bleach, ...
                                                                            data, ...
                                                                            image_data_post_bleach_model_unscaled, ...
                                                                            initial_condition_model_unscaled);

    [param_hat_1, ss_1] = lsqnonlin(fun_1, param_hat_1, lb_1, ub_1, options_1);
    
    disp([param_hat_2 param_hat_1])
    
    fun_2 = @(param)residual_recoverycurve_diffusion_and_binding_full(  param(1), ...
                                                                        param(2), ...
                                                                        param(3), ...
                                                                        param_hat_1(1), ...
                                                                        x_bleach, ...
                                                                        y_bleach, ...
                                                                        r_bleach, ...
                                                                        param_hat_1(2), ...
                                                                        param_hat_1(3), ...
                                                                        delta_t, ...
                                                                        number_of_time_points_fine_per_coarse, ...
                                                                        number_of_pad_pixels, ...
                                                                        data);

    [param_hat_2, ss_2] = lsqnonlin(fun_2, param_hat_2, lb_2, ub_2, options_2);
    
    param_hat = [param_hat ; [param_hat_2 param_hat_1]];
    ss = [ss ; ss_2];
    
    disp([param_hat_2 param_hat_1])
    
    if size(param_hat, 1) >= 2
        max_jump = max(abs(param_hat(end, :) - param_hat(end - 1, :)));
        disp(max_jump)
        if max_jump < 1e-3
            is_converged = true;
        end
    end
            
end

param_hat_recoverycurvebased = param_hat;
ss_recoverycurvebased = ss;

save(['est_' num2str(random_seed) '.mat'], 'param_hat_pixelbased', 'ss_pixelbased', 'param_hat_recoverycurvebased', 'ss_recoverycurvebased')
