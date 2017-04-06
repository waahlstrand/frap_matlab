clear
clc
close all hidden

load('results.mat');

pixel_size = 7.5e-07; % m
D_SI_VECTOR = [5e-12, 1e-11, 5e-11, 1e-10, 5e-10];
D_VECTOR = D_SI_VECTOR / pixel_size^2;
K_ON_VECTOR = [0.05, 0.1, 0.5, 1, 5];
K_OFF_VECTOR = [0.05, 0.1, 0.5, 1, 5];
SIGMA_VECTOR = [0.1, 0.2, 0.5];

NUMBER_OF_SIMULATIONS = [];
THETA_TRUE = [];
MSE_PX = [];
MSE_RC = [];
RE = [];


for ind_D = 1:5
    for ind_k_on = 1:5
        for ind_k_off = 1:5
            for ind_sigma = 1%1:3
                disp([ind_D ind_k_on ind_k_off ind_sigma])
                ind = abs(PARAM_TRUE(:, 1) - D_VECTOR(ind_D)) < 1e-3 & abs(PARAM_TRUE(:, 2) - K_ON_VECTOR(ind_k_on)) < 1e-3 & abs(PARAM_TRUE(:, 3) - K_OFF_VECTOR(ind_k_off)) < 1e-3 & abs(SIGMA_NOISE - SIGMA_VECTOR(ind_sigma)) < 1e-3;
                ind = find(ind);
                NUMBER_OF_SIMULATIONS = [NUMBER_OF_SIMULATIONS ; numel(ind)];
                THETA_TRUE = [THETA_TRUE ; D_VECTOR(ind_D) K_ON_VECTOR(ind_k_on) K_OFF_VECTOR(ind_k_off) SIGMA_VECTOR(ind_sigma)];
                MSE_PX = [MSE_PX ; mean( (PARAM_HAT_PX(ind, 1:3) - [D_VECTOR(ind_D) K_ON_VECTOR(ind_k_on) K_OFF_VECTOR(ind_k_off)]).^2, 1)];
                MSE_RC = [MSE_RC ; mean( (PARAM_HAT_RC(ind, 1:3) - [D_VECTOR(ind_D) K_ON_VECTOR(ind_k_on) K_OFF_VECTOR(ind_k_off)]).^2, 1)];
                RE = [RE ; MSE_PX(end, :) ./ MSE_RC(end, :)];
            end
        end
    end
end

[THETA_TRUE RE]