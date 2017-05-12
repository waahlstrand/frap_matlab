clear
clc
close all hidden

folder = '\\sp.se\FB\FBs\SM2Open\Tommy\CLSM\20170425_Sticky_90Glu';
prefix_prebleach = '20170425_Sticky_90Glu_FRAP Pre Series71_t';
prefix_postbleach = '20170425_Sticky_90Glu_FRAP Pb1 Series71_t';
postfix = '_ch00';
number_of_images_prebleach = 50;
number_of_images_postbleach = 70;

number_of_pixels = 256;
bit_depth = 8;

delta_t = 0.500; % s
pixel_size = 0.36621e-6; % m
x_bleach = number_of_pixels / 2; % pixels
y_bleach = number_of_pixels / 2; % pixels
lx_bleach = 15e-6 / pixel_size; % m to pixels conversion
ly_bleach = lx_bleach;

%% Load prebleach data.
data_prebleach = zeros(number_of_pixels, number_of_pixels, number_of_images_prebleach);
for current_image = 1:number_of_images_prebleach
    index_str = num2str(current_image - 1);
    while numel(index_str) < 3
        index_str = ['0' index_str];
    end
    file_path = [folder '\' prefix_prebleach index_str postfix '.tif'];
    
    data_prebleach(:, :, current_image) = double(imread(file_path));
end
data_prebleach = data_prebleach / (2^bit_depth - 1);

%% Load postbleach data.
data = zeros(number_of_pixels, number_of_pixels, number_of_images_postbleach);

for current_image = 1:number_of_images_postbleach
    index_str = num2str(current_image - 1);
    while numel(index_str) < 3
        index_str = ['0' index_str];
    end
    file_path = [folder '\' prefix_postbleach index_str postfix '.tif'];
    
    data(:, :, current_image) = double(imread(file_path));
end
data = data / (2^bit_depth - 1);

number_of_images = number_of_images_postbleach;

%% Preprocessing.
data_prebleach_avg = mean(data_prebleach, 3);
% return
data = data - repmat(data_prebleach_avg, [1, 1, number_of_images_postbleach]) + mean(data_prebleach_avg(:));
% warning('no preprocessing')

clear file_path folder index_str postfix prefix current_image
clear data_prebleach data_prebleach_avg
clear number_of_images_prebleach number_of_images_postbleach

%% Save data.
save('data_sp2.mat');
