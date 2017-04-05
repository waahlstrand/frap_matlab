%% Initialization.
%clear
clc
close all hidden

%% Read file.
file_id = fopen('simulated_stochastic_data_1.0e-10_5.0_5.0.dat');

D = fread(file_id, 1, 'float64');
k_on = fread(file_id, 1, 'float64');
k_off = fread(file_id, 1, 'float64');
mf = fread(file_id, 1, 'float64');
Ib = fread(file_id, 1, 'float64');
Iu = fread(file_id, 1, 'float64');
r_bleach = fread(file_id, 1, 'float64');
number_of_pixels = fread(file_id, 1, 'int64');
number_of_images = fread(file_id, 1, 'int64');
number_of_pad_pixels = fread(file_id, 1, 'int64');
delta_t = fread(file_id, 1, 'float64');
pixel_size = fread(file_id, 1, 'float64');
number_of_particles = fread(file_id, 1, 'int64');
data = reshape(fread(file_id, 'int64'), [number_of_pixels, number_of_pixels, number_of_images]);

fclose(file_id);

clear file_id 

%% Save.
save('simulated_frap_data.mat')

%% Plot data.
k = 1 / (Ib*pi*r_bleach^2 + Iu * ((number_of_pixels+2*number_of_pad_pixels)^2 - pi*r_bleach^2));
%k*Ib*pi*r_bleach^2+k*Iu*((number_of_pixels+2*number_of_pad_pixels)^2-pi*r_bleach^2)
%expected = number_of_particles * Iu * k;
data = data / (number_of_particles * k);

%M = load('../inference/simulated_data_diffusion_and_binding.mat');

figure
hold on
%imagesc([reshape(M.data, [number_of_pixels, number_of_pixels*number_of_images]) ; reshape(data-M.data, [number_of_pixels, number_of_pixels*number_of_images])]);
%imagesc(reshape(data-M.data, [number_of_pixels, number_of_pixels*number_of_images]));
imagesc(reshape(data, [number_of_pixels, number_of_pixels*number_of_images]));
axis 'equal'
axis([0 number_of_images*number_of_pixels 0 number_of_pixels])
axis off
hold off

%figure
%hold on
%hist(data(:) - M.data(:), 100)
%hold off

%sum(M.data(:))
%sum(data(:))
%std(data(:) - M.data(:))

%max(data(:))
%min(data(:))
% mean(data(:)) / mean(data(:,1,1))