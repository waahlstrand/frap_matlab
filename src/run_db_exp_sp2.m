clear
clc
close all hidden

addpath('signal');
addpath('estimation');

%% Load SP2 data.
load('data_sp2.mat');

number_of_images = 25;
data = data(:, :, 1:number_of_images);

if exist('r_bleach', 'var') % Circular.
    param_bleach = [x_bleach, y_bleach, r_bleach];
else % Rectangular.
    param_bleach = [x_bleach, y_bleach, lx_bleach, ly_bleach];
end

%% Estimate parameters.
number_of_pad_pixels = 128;

lb_D_SI = 1e-12;
ub_D_SI = 1e-9;
lb_D = lb_D_SI / pixel_size^2;
ub_D = ub_D_SI / pixel_size^2;

lb_k_on = 0;
ub_k_on = 10;

lb_k_off = 0;
ub_k_off = 10;

lb_mf = 0.6;
ub_mf = 1.0;

lb_Ib = 0.0;
ub_Ib = 1.0;

lb_Iu = 0.0;
ub_Iu = 1.2;

lb = [lb_D, lb_k_on, lb_k_off, lb_mf, lb_Ib, lb_Iu]; 
ub = [ub_D, ub_k_on, ub_k_off, ub_mf, ub_Ib, ub_Iu]; 

param_guess = [];
number_of_fits = 1;

[param_hat, ss] = estimate_db_px( ...
    data, ...
    param_bleach, ...
    delta_t, ...
    number_of_pixels, ...
    number_of_images, ...
    number_of_pad_pixels, ...
    lb, ...
    ub, ...
    param_guess, ...
    number_of_fits)
                            
%% Show images.
model = signal_db( ...
    param_hat(1), ...
    param_hat(2), ...
    param_hat(3), ...
    param_hat(4), ...
    param_hat(5), ...
    param_hat(6), ...
    param_bleach, ...
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

if numel(param_bleach) == 3 % Circular.
    ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
else % Rectangular.
    ind = find( X >= x_bleach - 0.5 * lx_bleach & X <= x_bleach + 0.5 * lx_bleach & Y >= y_bleach - 0.5 * ly_bleach & Y <= y_bleach + 0.5 * ly_bleach );
end
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

D = param_hat(1) * pixel_size^2
k_on = param_hat(2)
k_off = param_hat(3)
