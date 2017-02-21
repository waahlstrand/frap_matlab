%% Initialization.
%clear
clc
close all hidden

%% Read file.
file_id = fopen('simulated_frap_data.dat');

D = fread(file_id, 1, 'float64');
k_on = fread(file_id, 1, 'float64');
k_off = fread(file_id, 1, 'float64');
mobile_fraction = fread(file_id, 1, 'float64');
delta_t = fread(file_id, 1, 'float64');
r_bleach = fread(file_id, 1, 'float64');
pixel_size = fread(file_id, 1, 'float64');
number_of_pixels = fread(file_id, 1, 'int64');
number_of_post_bleach_images = fread(file_id, 1, 'int64');

image_data_post_bleach = reshape(fread(file_id, 'int64'), [number_of_pixels, number_of_pixels, number_of_post_bleach_images]);

fclose(file_id);

clear file_id 

%% Save.
save('simulated_frap_data.mat')

%% Plot data.
figure
imagesc(reshape(image_data_post_bleach, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]));
axis 'equal'
axis([0 number_of_post_bleach_images*number_of_pixels 0 number_of_pixels])
axis off

%max(image_data_post_bleach(:))
%min(image_data_post_bleach(:))
mean(image_data_post_bleach(:)) / mean(image_data_post_bleach(:,1,1))