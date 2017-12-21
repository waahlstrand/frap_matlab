%% Initialization.

clear
clc
close all hidden

addpath('signal');
addpath('estimation');

random_seed = round(sum(1e6 * clock()));
random_stream = RandStream('mt19937ar','Seed',random_seed);
RandStream.setGlobalStream(random_stream);

%% Load data.

load('\\sp.se\FB\FBs\SM\samdata\Annika K\Lantmännen\FRAP 171109\FITC samt inm viscoferm 5 punkter\frap_001.mat');

%% Settings.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit type, either pixel-based ('px') or recovery curve-based ('rc').
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimation_mode = 'px';
% estimation_mode = 'rc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bleaching (and laser fluctuation) correction. Specified by a list of 
% indices, the mean intensity in each image of which to divide each image 
% by. As described in e.g. Deschout (2010) and Xiong (2016). If variable is
% [] then no bleaching correction is performed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bleaching_correction_indices = [];

number_of_pixels = size(experiment.postbleach.image_data, 1);
[subX, subY] = ndgrid(1:number_of_pixels, 1:number_of_pixels);
ind = ones(number_of_pixels, number_of_pixels);
ind(2:end-1, 2:end-1) = 0;
subX = subX(ind(:) == 1);
subY = subY(ind(:) == 1);
bleaching_correction_indices = sub2ind([number_of_pixels number_of_pixels], subX, subY);
clear subX subY ind

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Background correction, either 'none', 'subtraction' of the average pre-
% bleach image (as in Jonasson (2010) or 'division' of the average pre-
% bleach image (as in Deschout (2010) and Xiong (2016)). Also, a kernel
% size for median prefiltering needs to be specified (with 1 being no
% prefiltering).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

background_correction = 'none';
% background_correction = 'subtraction';
% background_correction = 'division';

background_smoothing = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel executation in fitting.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_parallel = true;

%% Preprocess data.

[data, data_prebleach] = preprocess(experiment, bleaching_correction_indices, background_correction, background_smoothing);

%% Perform background subtraction.

number_of_pixels = size(data, 1);
number_of_images = size(data, 3);
number_of_images_prebleach = size(data_prebleach, 3);

x_bleach =  - experiment.bleach.bleach_position_x / experiment.postbleach.pixel_size_x + number_of_pixels / 2;
y_bleach = experiment.bleach.bleach_position_y / experiment.postbleach.pixel_size_y + number_of_pixels / 2;

if isequal(experiment.bleach.bleach_type, 'circle')
    r_bleach = 0.5 * experiment.bleach.bleach_size_x / experiment.postbleach.pixel_size_x;
%     r_bleach = experiment.bleach.bleach_size_x / experiment.postbleach.pixel_size_x
    param_bleach = [x_bleach, y_bleach, r_bleach];
elseif isequal(experiment.bleach.bleach_type, 'rectangle')
    lx_bleach = experiment.bleach.bleach_size_x / experiment.postbleach.pixel_size_x;
    ly_bleach = experiment.bleach.bleach_size_y / experiment.postbleach.pixel_size_y;
    param_bleach = [x_bleach, y_bleach, lx_bleach, ly_bleach];
end

pixel_size = experiment.postbleach.pixel_size_x;
delta_t = experiment.postbleach.time_frame;

%% Compute bleach region indices.

[X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
X = X - 0.5;
Y = Y - 0.5;

if isequal(experiment.bleach.bleach_type, 'circle')
    bleach_region = (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2;
elseif isequal(experiment.bleach.bleach_type, 'rectangle')
    bleach_region = X >= x_bleach - 0.5 * lx_bleach & X <= x_bleach + 0.5 * lx_bleach & Y >= y_bleach - 0.5 * ly_bleach & Y <= y_bleach + 0.5 * ly_bleach;
end
ind = find(bleach_region(:));

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

[param_hat, ss] = estimate_d( ...
    data_prebleach, ...
    data, ...
    param_bleach, ...
    delta_t, ...
    number_of_pad_pixels, ...
    lb, ...
    ub, ...
    param_guess, ...
    number_of_fits, ...
    estimation_mode, ...
    use_parallel);

D = param_hat(1) * pixel_size^2
mf = param_hat(2);
Ib = param_hat(3);
Iu = param_hat(4);


%% Show results.

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

rc_data = zeros(1, number_of_images_prebleach + number_of_images);
for current_image = 1:number_of_images_prebleach
    slice = data_prebleach(:, :, current_image);
    rc_data(current_image) = mean(slice(ind));
end
for current_image = 1:number_of_images
    slice = data(:, :, current_image);
    rc_data(number_of_images_prebleach + current_image) = mean(slice(ind));
end

rc_model = zeros(1, number_of_images_prebleach + number_of_images);
for current_image = 1:number_of_images_prebleach
    rc_model(current_image) = Iu;
end
for current_image = 1:number_of_images
    slice = model(:, :, current_image);
    rc_model(number_of_images_prebleach + current_image) = mean(slice(ind));
end

figure, hold on
plot([(-number_of_images_prebleach:-1)*delta_t (1:number_of_images)*delta_t], rc_data, 'ro');
plot([(-number_of_images_prebleach:-1)*delta_t (1:number_of_images)*delta_t], rc_model, 'k-');



