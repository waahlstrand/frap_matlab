%% Initialization.
clear
clc
close all hidden

%% Read file.
file_id = fopen('simulated_stochastic_data.bin');

D = fread(file_id, 1, 'float64');
k_on = fread(file_id, 1, 'float64');
k_off = fread(file_id, 1, 'float64');
mobile_fraction = fread(file_id, 1, 'float64');
alpha = fread(file_id, 1, 'float64');
r_bleach = fread(file_id, 1, 'float64');
number_of_pixels = fread(file_id, 1, 'int64');
number_of_prebleach_frames = fread(file_id, 1, 'int64');
number_of_bleach_frames = fread(file_id, 1, 'int64');
number_of_postbleach_frames = fread(file_id, 1, 'int64');
number_of_pad_pixels = fread(file_id, 1, 'int64');
delta_t = fread(file_id, 1, 'float64');
pixel_size = fread(file_id, 1, 'float64');
number_of_particles = fread(file_id, 1, 'int64');
t_exec = fread(file_id, 1, 'float64');
data = reshape(fread(file_id, 'int64'), [number_of_pixels, number_of_pixels, number_of_prebleach_frames + number_of_postbleach_frames]);

C_prebleach_sim = double(data(:, :, 1:number_of_prebleach_frames));
C_postbleach_sim = double(data(:, :, number_of_prebleach_frames+1:end));

fclose(file_id);

clear file_id 

%% Normalize count data.

C0 = 1;
norm_factor = number_of_particles / (number_of_pixels + 2 * number_of_pad_pixels)^2; % Average count in pixels.
C_prebleach_sim = C_prebleach_sim / norm_factor;
C_postbleach_sim = C_postbleach_sim / norm_factor;
%% Save.
% save('simulated_frap_data.mat')








% % % %% Plot data.
% % % k = 1 / (Ib*pi*r_bleach^2 + Iu * ((number_of_pixels+2*number_of_pad_pixels)^2 - pi*r_bleach^2));
% % % %k*Ib*pi*r_bleach^2+k*Iu*((number_of_pixels+2*number_of_pad_pixels)^2-pi*r_bleach^2)
% % % %expected = number_of_particles * Iu * k;
% % % data = data / (number_of_particles * k);
% % % return
% % % %M = load('../inference/simulated_data_diffusion_and_binding.mat');
% % % 
% % % data_deterministic = load('../../../../sim_study_stochastic_simulation/results_stochastic/simulated_stochastic_data_1.0e-10_0.5_1.0.mat');
% % % 
% % % figure
% % % hold on
% % % %imagesc([reshape(M.data, [number_of_pixels, number_of_pixels*number_of_images]) ; reshape(data-M.data, [number_of_pixels, number_of_pixels*number_of_images])]);
% % % %imagesc(reshape(data-M.data, [number_of_pixels, number_of_pixels*number_of_images]));
% % % imagesc(reshape(data-data_deterministic.data_stochastic, [number_of_pixels, number_of_pixels*number_of_images]), [-0.1,0.1]);
% % % %axis 'equal'
% % % axis([0 number_of_images*number_of_pixels 0 number_of_pixels])
% % % % axis off
% % % hold off
% % % 
% % % %figure
% % % %hold on
% % % %hist(data(:) - M.data(:), 100)
% % % %hold off
% % % 
% % % %sum(M.data(:))
% % % %sum(data(:))
% % % %std(data(:) - M.data(:))
% % % 
% % % %max(data(:))
% % % %min(data(:))
% % % % mean(data(:)) / mean(data(:,1,1))