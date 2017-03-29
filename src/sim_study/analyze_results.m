clear
clc
close all hidden

load('results.mat');

sigma_vector = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05];

i = 6;
std(param_hat_px(sigma_noise == sigma_vector(i), :), [], 1) ./ param_true(1, :)
std(param_hat_rc(sigma_noise == sigma_vector(i), :), [], 1) ./ param_true(1, :)