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

sigma_noise = 0.0;%0.025;
data = data + sigma_noise * randn(size(data));

%% Estimate parameters.

% Set parameter bounds.
lb_D_SI = 1e-12;
ub_D_SI = 1e-8;
lb_D = lb_D_SI / pixel_size^2;
ub_D = ub_D_SI / pixel_size^2;

lb_k_on = 0;
ub_k_on = 10;

lb_k_off = 0;
ub_k_off = 10;

lb_mf = 0.0;
ub_mf = 1.0;

lb_Ib = 0.0;
ub_Ib = 1.0;

lb_Iu = 0.0;
ub_Iu = 1.0;

lb = [lb_D, lb_k_on, lb_k_off, lb_mf, lb_Ib, lb_Iu]; 
ub = [ub_D, ub_k_on, ub_k_off, ub_mf, ub_Ib, ub_Iu]; 

% Optimization options.
options = optimoptions(@lsqnonlin);
options.Algorithm = 'trust-region-reflective';
options.Display = 'iter';
options.FunctionTolerance = 1e-6;
options.OptimalityTolerance = 1e-6;
options.StepTolerance = 1e-6;
options.UseParallel = true;

% Optimization.
% param_hat = [];
% ss = [];
param_guess = lb + (ub - lb) .* rand(size(lb));

fun = @(param)residual_db(  param(1), ...
                            param(2), ...
                            param(3), ...
                            param(4), ...
                            param(5), ...
                            param(6), ...
                            x_bleach, ...
                            y_bleach, ...
                            r_bleach, ...
                            delta_t, ...
                            number_of_pixels, ...
                            number_of_images, ...
                            number_of_pad_pixels, ...
                            data);

[param_hat, ss] = lsqnonlin(fun, param_guess, lb, ub, options);
   