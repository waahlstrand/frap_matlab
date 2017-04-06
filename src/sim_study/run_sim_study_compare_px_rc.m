%% Initialization.
clear
clc
close all hidden

addpath('..');
addpath('../signal_pde');

random_seed = sum( 1e6 * clock() );
random_stream = RandStream('mt19937ar', 'Seed', random_seed);
RandStream.setGlobalStream(random_stream);

WALL_TIME = 3600 * 12;
time_start = tic();

%% Run simulation study.
pixel_size = 7.5e-07; % m
delta_t = 0.2; % s
number_of_pixels = 256;
number_of_images = 100;
number_of_pad_pixels = 128;

mf = 1.0; % dimensionless
Ib = 0.6; % a.u.
Iu = 1.0; % a.u.

x_bleach = 128; % pixels
y_bleach = 128; % pixels
r_bleach = 15e-6 / pixel_size; % pixels corresponding to 15 µm radius (30 µm diameter)

% Set parameter bounds.
lb_D_SI = 1e-12;
ub_D_SI = 1e-9;
lb_D = lb_D_SI / pixel_size^2;
ub_D = ub_D_SI / pixel_size^2;

lb_k_on = 0;
ub_k_on = 100;

lb_k_off = 0;
ub_k_off = 100;

lb_mf = 1.0;
ub_mf = 1.0;

lb_Ib = 0.0;
ub_Ib = 1.0;

lb_Iu = 0.0;
ub_Iu = 1.1;

lb = [lb_D, lb_k_on, lb_k_off, lb_mf, lb_Ib, lb_Iu]; 
ub = [ub_D, ub_k_on, ub_k_off, ub_mf, ub_Ib, ub_Iu]; 

number_of_fits = 1;

D_SI_VECTOR = [5e-12, 1e-11, 5e-11, 1e-10, 5e-10];
K_ON_VECTOR = [0.05, 0.1, 0.5, 1, 5];
K_OFF_VECTOR = [0.05, 0.1, 0.5, 1, 5];

PARAM_TRUE = [];
PARAM_HAT_PX = [];
PARAM_HAT_RC = [];
SS_PX = [];
SS_RC = [];
SIGMA_NOISE = [];

% Loop forever.
while toc(time_start) < WALL_TIME - 600
    % Randomize true parameters.
    D_SI = randsample(D_SI_VECTOR, 1); % m^2/s
    D = D_SI / pixel_size^2; % pixels^2 / s
    
    k_on = randsample(K_ON_VECTOR, 1); % 1/s
    k_off = randsample(K_OFF_VECTOR, 1); % 1/s
    
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

    sigma_noise = randsample([0.1 0.2 0.5], 1);
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
                                        
    PARAM_TRUE = [PARAM_TRUE ; param_true];
    PARAM_HAT_PX = [PARAM_HAT_PX ; param_hat_px];
    PARAM_HAT_RC = [PARAM_HAT_RC ; param_hat_rc];
    SS_PX = [SS_PX ; ss_px];
    SS_RC = [SS_RC ; ss_rc];
    SIGMA_NOISE = [SIGMA_NOISE ; sigma_noise];
end
save(['est_' num2str(random_seed) '.mat'], 'PARAM_TRUE', 'PARAM_HAT_PX', 'PARAM_HAT_RC', 'SS_PX', 'SS_RC', 'SIGMA_NOISE');
