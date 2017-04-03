%% Initialization.
%clear
clc
close all hidden

%% Read file.
file_id = fopen('simulated_frap_data_5.0e-10_1.0_10.0.dat');

D = fread(file_id, 1, 'float64');
k_on = fread(file_id, 1, 'float64');
k_off = fread(file_id, 1, 'float64');
mobile_fraction = fread(file_id, 1, 'float64');
delta_t = fread(file_id, 1, 'float64');
r_bleach = fread(file_id, 1, 'float64');
pixel_size = fread(file_id, 1, 'float64');
number_of_pixels = fread(file_id, 1, 'int64');
number_of_post_bleach_images = fread(file_id, 1, 'int64');
number_of_pad_pixels = fread(file_id, 1, 'int64');
intensity_inside_bleach_region = fread(file_id, 1, 'float64');
intensity_outside_bleach_region = fread(file_id, 1, 'float64');
number_of_particles = fread(file_id, 1, 'int64');

image_data_post_bleach = reshape(fread(file_id, 'int64'), [number_of_pixels, number_of_pixels, number_of_post_bleach_images]);

fclose(file_id);

clear file_id 

%% Save.
save('simulated_frap_data.mat')

%% Plot data.
k = 1 / (intensity_inside_bleach_region*pi*r_bleach^2 + intensity_outside_bleach_region * ((number_of_pixels+2*number_of_pad_pixels)^2 - pi*r_bleach^2));
%k*intensity_inside_bleach_region*pi*r_bleach^2+k*intensity_outside_bleach_region*((number_of_pixels+2*number_of_pad_pixels)^2-pi*r_bleach^2)
%expected = number_of_particles * intensity_outside_bleach_region * k;
image_data_post_bleach = image_data_post_bleach / (number_of_particles * k);

%M = load('../inference/simulated_data_diffusion_and_binding.mat');

figure
hold on
%imagesc([reshape(M.image_data_post_bleach, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]) ; reshape(image_data_post_bleach-M.image_data_post_bleach, [number_of_pixels, number_of_pixels*number_of_post_bleach_images])]);
%imagesc(reshape(image_data_post_bleach-M.image_data_post_bleach, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]));
imagesc(reshape(image_data_post_bleach, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]));
axis 'equal'
axis([0 number_of_post_bleach_images*number_of_pixels 0 number_of_pixels])
axis off
hold off

%figure
%hold on
%hist(image_data_post_bleach(:) - M.image_data_post_bleach(:), 100)
%hold off

%sum(M.image_data_post_bleach(:))
%sum(image_data_post_bleach(:))
%std(image_data_post_bleach(:) - M.image_data_post_bleach(:))

%max(image_data_post_bleach(:))
%min(image_data_post_bleach(:))
% mean(image_data_post_bleach(:)) / mean(image_data_post_bleach(:,1,1))