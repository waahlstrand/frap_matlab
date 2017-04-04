clear
clc
close all hidden

folder = 'results_20170403';
files = dir([folder '\*.mat']);
number_of_files = numel(files);

param_true = nan(number_of_files, 6);
param_hat_px = nan(number_of_files, 6);
param_hat_rc = nan(number_of_files, 6);
sigma_noise = nan(number_of_files, 1);

ok = true(number_of_files, 1);

for current_file = 1:number_of_files
    %if mod(current_file, 100) == 0
        disp(current_file)
    %end
    
    file_path = [folder '\' files(current_file).name];
    
    try
        file_data = load(file_path);

        param_true(current_file, :) = file_data.param_true;
        param_hat_px(current_file, :) = file_data.param_hat_px;
        param_hat_rc(current_file, :) = file_data.param_hat_rc;
        sigma_noise(current_file, :) = file_data.sigma_noise;
    catch
        disp(['Deleting file # ' num2str(current_file)])
        ok(current_file) = false;
        delete(file_path)
    end
end
param_true = param_true(ok, :);
param_hat_px = param_hat_px(ok, :);
param_hat_rc = param_hat_rc(ok, :);
sigma_noise = sigma_noise(ok);

clear files number_of_files current_file file_path file_data ok

save('results.mat')
    