clear
clc
close all hidden

addpath('signal');
addpath('estimation');

%% Load data.
% load('\\sp.se\FB\FBs\SM\samdata\Annika K\Lantmännen\FRAP 170418 Viscoferm inm test\buffert\frap_002.mat');
load('\\sp.se\FB\FBs\SM\samdata\Annika K\Lantmännen\FRAP 170418 Viscoferm inm test\Viscoferm\frap_002.mat');

data = experiment.postbleach.image_data;%(:, :, 1:number_of_images);
data = double(data);
data = data / (2^experiment.postbleach.bit_depth - 1);
number_of_images = size(data, 3);
% 
% data_prebleach = experiment.prebleach.image_data;
% data_prebleach = double(data_prebleach);
% data_prebleach = data_prebleach / (2^experiment.prebleach.bit_depth - 1);
% data_prebleach_avg = mean(data_prebleach, 3);
% data = data - repmat(data_prebleach_avg, [1, 1, number_of_images]) + mean(data_prebleach_avg(:));

number_of_pixels = size(experiment.postbleach.image_data, 1);
x_bleach = - experiment.bleach.bleach_position_x / experiment.postbleach.pixel_size_x + number_of_pixels / 2;
y_bleach = experiment.bleach.bleach_position_y / experiment.postbleach.pixel_size_y + number_of_pixels / 2;
if isequal(experiment.bleach.bleach_type, 'circle')
    r_bleach = experiment.bleach.bleach_size_x / experiment.postbleach.pixel_size_x;
    param_bleach = [x_bleach, y_bleach, r_bleach];
elseif isequal(experiment.bleach.bleach_type, 'rectangle')
    lx_bleach = experiment.bleach.bleach_size_x / experiment.postbleach.pixel_size_x;
    ly_bleach = experiment.bleach.bleach_size_y / experiment.postbleach.pixel_size_y;
    param_bleach = [x_bleach, y_bleach, lx_bleach, ly_bleach];
end

pixel_size = experiment.postbleach.pixel_size_x;
delta_t = experiment.postbleach.time_frame;



%% Estimate parameters.
number_of_pad_pixels = 128;

lb_D_SI = 1e-13;
ub_D_SI = 2e-9;
lb_D = lb_D_SI / pixel_size^2;
ub_D = ub_D_SI / pixel_size^2;

lb_mf = 0.6;
ub_mf = 1.0;

lb_Ib = 0.0;
ub_Ib = 1.0;

lb_Iu = 0.0;
ub_Iu = 1.2;

lb = [lb_D, lb_mf, lb_Ib, lb_Iu]; 
ub = [ub_D, ub_mf, ub_Ib, ub_Iu]; 

param_guess = [];
number_of_fits = 1;

[param_hat, ss] = estimate_d_px( ...
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
model = signal_d( ...
    param_hat(1), ...
    param_hat(2), ...
    param_hat(3), ...
    param_hat(4), ...
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
