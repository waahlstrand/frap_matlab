%% Initialization.
clear
clc
close all hidden

%% Read file.
file_id = fopen('simulated_frap_data.dat');
file_data = fread(file_id, '*int64');
fclose(file_id);

%% Extract parameters and image data.
number_of_pixels = file_data(1);
number_of_post_bleach_images = file_data(2);
image_data_post_bleach = reshape(file_data(3:end), [number_of_pixels, number_of_pixels, number_of_post_bleach_images]);

%% Plot data.
figure
imagesc(reshape(image_data_post_bleach, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]));
axis 'equal'
axis([0 number_of_post_bleach_images*number_of_pixels 0 number_of_pixels])
axis off