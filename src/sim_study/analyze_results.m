clear
clc
close all hidden

load('results.mat');

D_VECTOR = sort(unique(param_true(:, 1)));
K_ON_OFF_MATRIX = [];
for k_on_exp = -2:2
    for k_off_exp = -2:2
        if k_on_exp <= k_off_exp + 1
            K_ON_OFF_MATRIX = [ K_ON_OFF_MATRIX ; 10^k_on_exp , 10^k_off_exp ];
        end
    end
end
sigma_vector = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05];


ind_D = 1;
ind_k = 5;
ind_sigma = 6;
ind = param_true(:, 1) == D_VECTOR(ind_D) & param_true(:, 2) == K_ON_OFF_MATRIX(ind_k, 1) & param_true(:, 3) == K_ON_OFF_MATRIX(ind_k, 2) & sigma_noise == sigma_vector(ind_sigma);
sum(ind)
ind = find(ind);

std(param_hat_px(ind, :), [], 1) ./ param_true(ind(1), :)
std(param_hat_rc(ind, :), [], 1) ./ param_true(ind(1), :)
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