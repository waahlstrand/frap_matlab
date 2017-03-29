%% Initialization.
clear
clc
close all hidden

addpath('signal_pde');

%% Init parallel pool.
delete(gcp('nocreate'))
c = parcluster('local');
c.NumWorkers = 8;
parpool(c, c.NumWorkers);

%% Simulate data.
pixel_size = 7.598e-07; % m
delta_t = 0.2650; % s
number_of_pixels = 256;
number_of_images = 100
number_of_pad_pixels = 128;

D_SI = 2.5e-10; % m^2/s
D = D_SI / pixel_size^2; % pixels^2 / s
k_on = 0.5;%0.05; % 1/s
k_off = 1.0;%0.01; % 1/s
mf = 0.9; % dimensionless
Ib = 0.6; % a.u.
Iu = 0.9; % a.u.

x_bleach = 128; % pixels
y_bleach = 128; % pixels
r_bleach = 32; % pixels

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

param_true = [D, k_on, k_off, mf, Ib, Iu];
param_guess = param_true;
number_of_fits = 1;

% Loop forever in parallel.
parfor i = 1:100000
    random_seed = sum( 1e6 * clock() );
    random_stream = RandStream('mt19937ar', 'Seed', random_seed);
    RandStream.setGlobalStream(random_stream);

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

    sigma_noise = randsample([0.001 0.002 0.005 0.01 0.02 0.05], 1);
    data = data + sigma_noise * randn(size(data));

    [param_hat_px, ss_px] = estimate_db_px( data, ...
                                            x_bleach, ...
                                            y_bleach, ...
                                            r_bleach, ...
                                            delta_t, ...
                                            number_of_pixels, ...
                                            number_of_images, ...
                                            number_of_pad_pixels, ...
                                            lb, ...
                                            ub, ...
                                            param_guess, ...
                                            number_of_fits);

    [param_hat_rc, ss_rc] = estimate_db_rc( data, ...
                                            x_bleach, ...
                                            y_bleach, ...
                                            r_bleach, ...
                                            delta_t, ...
                                            number_of_pixels, ...
                                            number_of_images, ...
                                            number_of_pad_pixels, ...
                                            lb, ...
                                            ub, ...
                                            param_guess, ...
                                            number_of_fits);

    mat_file = matfile(['est_' num2str(random_seed) '.mat'], 'writable', true);
    mat_file.param_true = param_true;
    mat_file.sigma_noise = sigma_noise;
    mat_file.param_hat_px = param_hat_px;
    mat_file.ss_px = ss_px;
    mat_file.param_hat_rc = param_hat_rc;
    mat_file.ss_rc = ss_rc;
end