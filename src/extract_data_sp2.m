clear
clc
close all hidden

folder = '\\sp.se\FB\FBs\SM2Open\Hannes\18.4\fitcdextran50ppmnospheres';

number_of_pixels = 256;
bit_depth = 8;

%% Load prebleach data.
prefix = 'fitcdextran50ppmnospheres_FRAP Pre Series64_t';
postfix = '_ch00';
number_of_images = 30;

data_prebleach = zeros(number_of_pixels, number_of_pixels, number_of_images);
for current_image = 1:number_of_images
    index_str = num2str(current_image - 1);
    while numel(index_str) < 3
        index_str = ['0' index_str];
    end
    file_path = [folder '\' prefix index_str postfix '.tif'];
    
    data_prebleach(:, :, current_image) = double(imread(file_path));
end

data_prebleach = data_prebleach / (2^bit_depth - 1);

%% Load postbleach data.
prefix = 'fitcdextran50ppmnospheres_FRAP Pb1 Series64_t';
postfix = '_ch00';

number_of_images = 75;

delta_t = 0.500; % s
pixel_size = 0.1831e-6; % m
x_bleach = number_of_pixels / 2; % pixels
y_bleach = number_of_pixels / 2; % pixels
lx_bleach = 10e-6 / pixel_size; % µm to pixels conversion
ly_bleach = lx_bleach;
param_bleach = [lx_bleach ly_bleach];

data = zeros(number_of_pixels, number_of_pixels, number_of_images);

for current_image = 1:number_of_images
    index_str = num2str(current_image - 1);
    while numel(index_str) < 3
        index_str = ['0' index_str];
    end
    file_path = [folder '\' prefix index_str postfix '.tif'];
    
    data(:, :, current_image) = double(imread(file_path));
end

data = data / (2^bit_depth - 1);

%% Preprocessing.
data_prebleach_avg = mean(data_prebleach, 3);
data = data - repmat(data_prebleach_avg, [1, 1, number_of_images]) + mean(data_prebleach_avg(:));

clear file_path folder index_str postfix prefix current_image
clear data_prebleach data_prebleach_avg

%% Save data.
save('data_sp2.mat');
