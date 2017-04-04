clear
clc
close all hidden

load('results.mat');

D_VECTOR = sort(unique(param_true(:, 1)));
K_ON_OFF_MATRIX = [];
for k_on_exp = -2:1%2
    for k_off_exp = -2:1%2
        if k_on_exp <= k_off_exp + 1 && k_on_exp >= k_off_exp - 1 % SIC SIC SIC
            K_ON_OFF_MATRIX = [ K_ON_OFF_MATRIX ; 10^k_on_exp , 10^k_off_exp ];
        end
    end
end
sigma_vector = [0.001, 0.01, 0.1];

NUMBER_OF_SIMULATIONS = [];
THETA_TRUE = [];
MSE_PX = [];
MSE_RC = [];
RE = [];


for ind_D = 1:5
    for ind_k = 1:size(K_ON_OFF_MATRIX, 1)
        for ind_sigma = 3%1:3
%             disp([ind_D ind_k ind_sigma]) 
            ind = abs(param_true(:, 1) - D_VECTOR(ind_D)) < 1e-3 & abs(param_true(:, 2) - K_ON_OFF_MATRIX(ind_k, 1)) < 1e-3 & abs(param_true(:, 3) - K_ON_OFF_MATRIX(ind_k, 2)) < 1e-3 & abs(sigma_noise - sigma_vector(ind_sigma)) < 1e-3;
            ind = find(ind);
            NUMBER_OF_SIMULATIONS = [NUMBER_OF_SIMULATIONS ; numel(ind)];
            THETA_TRUE = [THETA_TRUE ; D_VECTOR(ind_D) K_ON_OFF_MATRIX(ind_k, :) sigma_vector(ind_sigma)];
            MSE_PX = [MSE_PX ; mean( (param_hat_px(ind, 1:3) - [D_VECTOR(ind_D) K_ON_OFF_MATRIX(ind_k, :)]).^2, 1)];
            MSE_RC = [MSE_RC ; mean( (param_hat_rc(ind, 1:3) - [D_VECTOR(ind_D) K_ON_OFF_MATRIX(ind_k, :)]).^2, 1)];
            RE = [RE ; MSE_PX(end, :) ./ MSE_RC(end, :)];
        end
    end
end

[THETA_TRUE RE]


% ind_D = 1
% ind_k = 5
% ind_sigma = 2
% ind = abs(param_true(:, 1) - D_VECTOR(ind_D)) < 1e-3 & abs(param_true(:, 2) - K_ON_OFF_MATRIX(ind_k, 1)) < 1e-3 & abs(param_true(:, 3) - K_ON_OFF_MATRIX(ind_k, 2)) < 1e-3 & abs(sigma_noise - sigma_vector(ind_sigma)) < 1e-3;
% ind = find(ind);
% [param_hat_px(ind, 1:3) param_hat_rc(ind, 1:3)]

%             
% 
% std(param_hat_px(ind, :), [], 1) ./ param_true(ind(1), :)
% std(param_hat_rc(ind, :), [], 1) ./ param_true(ind(1), :)
% 
% a = 0;
% for ind_D = 1:5
%     for ind_k = 1:19
%         for ind_sigma = 1:6
%             ind = param_true(:, 1) == D_VECTOR(ind_D) & param_true(:, 2) == K_ON_OFF_MATRIX(ind_k, 1) & param_true(:, 3) == K_ON_OFF_MATRIX(ind_k, 2) & sigma_noise == sigma_vector(ind_sigma);
%             a = a + sum(ind);
%         end
%     end
% end