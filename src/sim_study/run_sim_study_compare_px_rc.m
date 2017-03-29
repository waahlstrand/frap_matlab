%% Initialization.
clear
clc
close all hidden

addpath('..');
addpath('../signal_pde');

%% Init parallel pool.
delete(gcp('nocreate'))
c = parcluster('local');
c.NumWorkers = 8;
parpool(c, c.NumWorkers);

%% Run simulation study.
pixel_size = 7.598e-07; % m
delta_t = 0.2650; % s
number_of_pixels = 256;
number_of_images = 200;
number_of_pad_pixels = 128;

mf = 0.9; % dimensionless
Ib = 0.6; % a.u.
Iu = 0.9; % a.u.

x_bleach = 128; % pixels
y_bleach = 128; % pixels
% r_bleach = 32; % pixels
r_bleach = 15e-6 / pixel_size; % pixels corresponding to 15 µm radius (30 µm diameter)

% Set parameter bounds.
lb_D_SI = 1e-12;
ub_D_SI = 1e-9;
lb_D = lb_D_SI / pixel_size^2;
ub_D = ub_D_SI / pixel_size^2;

lb_k_on = 0;
ub_k_on = 150;

lb_k_off = 0;
ub_k_off = 150;

lb_mf = 0.0;
ub_mf = 1.0;

lb_Ib = 0.0;
ub_Ib = 1.0;

lb_Iu = 0.0;
ub_Iu = 1.0;

lb = [lb_D, lb_k_on, lb_k_off, lb_mf, lb_Ib, lb_Iu]; 
ub = [ub_D, ub_k_on, ub_k_off, ub_mf, ub_Ib, ub_Iu]; 

number_of_fits = 1;

D_SI_VECTOR = [5e-12, 1e-11, 5e-11, 1e-10, 5e-10];
K_ON_OFF_MATRIX = [];
for k_on_exp = -2:2
    for k_off_exp = -2:2
        if k_on_exp <= k_off_exp + 1
            K_ON_OFF_MATRIX = [ K_ON_OFF_MATRIX ; 10^k_on_exp , 10^k_off_exp ];
        end
    end
end
    
% Loop forever in parallel.
parfor i = 1:100000
    random_seed = sum( 1e6 * clock() );
    random_stream = RandStream('mt19937ar', 'Seed', random_seed);
    RandStream.setGlobalStream(random_stream);
    
    % Randomize true parameters.
    D_SI = randsample(D_SI_VECTOR, 1); % m^2/s
    D = D_SI / pixel_size^2 % pixels^2 / s
    
    ind = randsample(1:size(K_ON_OFF_MATRIX, 1), 1);
    k_on = K_ON_OFF_MATRIX(ind, 1) % 1/s
    k_off = K_ON_OFF_MATRIX(ind, 2) % 1/s
    
    param_true = [D, k_on, k_off, mf, Ib, Iu];
    param_guess = param_true;

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