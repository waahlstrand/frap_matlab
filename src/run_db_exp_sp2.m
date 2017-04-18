clear
clc
close all hidden

addpath('signal_pde');

%% Load data.
folder = '\\sp.se\FB\FBs\SM2Open\Hannes\13.4\pss002naf';
prefix = 'pss002naf_FRAP Pb1 Series41_t';
postfix = '_ch00';

number_of_images = 50;
number_of_pixels = 256;
bit_depth = 12;
delta_t = 0.500; % s
pixel_size = 0.732422e-6; % m
x_bleach = number_of_pixels / 2; % pixels
y_bleach = number_of_pixels / 2; % pixels
r_bleach = 25e-6 / pixel_size; % µm to pixels conversion

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

clear file_path folder index_str postfix prefix

%% Estimate parameters.
number_of_pad_pixels = 128;

lb_D_SI = 1e-11;
ub_D_SI = 1e-8;
lb_D = lb_D_SI / pixel_size^2;
ub_D = ub_D_SI / pixel_size^2;

lb_k_on = 0;
ub_k_on = 10;

lb_k_off = 0;
ub_k_off = 10;

lb_mf = 0.0;
ub_mf = 1.0;

lb_Ib = 0.0;
ub_Ib = 1.0;

lb_Iu = 0.0;
ub_Iu = 1.2;

lb = [lb_D, lb_k_on, lb_k_off, lb_mf, lb_Ib, lb_Iu]; 
ub = [ub_D, ub_k_on, ub_k_off, ub_mf, ub_Ib, ub_Iu]; 

param_guess = [];
number_of_fits = 1;

[param_hat_px, ss_px] = estimate_db_px( data, ...
                                        x_bleach, ...
                                        y_bleach, ...
                                        r_bleach, ...
                                        delta_t, ...
                                        number_of_pixels, ...
                                        number_of_images, ...
                                        number_of_pad_pixels, ...
                                        lb, ...
                                        ub, ...
                                        param_guess, ...
                                        number_of_fits);
                            
%% Show images.
model = signal_db(  param_hat_px(1), ...
                    param_hat_px(2), ...
                    param_hat_px(3), ...
                    param_hat_px(4), ...
                    param_hat_px(5), ...
                    param_hat_px(6), ...
                    x_bleach, ...
                    y_bleach, ...
                    r_bleach, ...
                    delta_t, ...
                    number_of_pixels, ...
                    number_of_images, ...
                    number_of_pad_pixels);
                
figure, imagesc([reshape(data, [number_of_pixels, number_of_pixels * number_of_images]) ; reshape(model, [number_of_pixels, number_of_pixels * number_of_images])])
figure, imagesc(reshape(data - model, [number_of_pixels, number_of_pixels * number_of_images]))

%% Show recovery curve.
[X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
X = X - 0.5;
Y = Y - 0.5;
ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
ind = ind(:);

rc_data = zeros(1, number_of_images);
for current_image = 1:number_of_images
    slice = data(:, :, current_image);
    rc_data(current_image) = mean(slice(ind));
end

rc_model = zeros(1, number_of_images);
for current_image = 1:number_of_images
    slice = model(:, :, current_image);
    rc_model(current_image) = mean(slice(ind));
end

figure, hold on
plot((1:number_of_images)*delta_t, rc_data, 'ro');
plot((1:number_of_images)*delta_t, rc_model, 'k-');
