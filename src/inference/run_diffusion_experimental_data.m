%% Initialization.
clear
clc
close all hidden

%% Load data.

% Load.
% file_path = '../../../data/data_binding_early_wood_6A/FRAP_002.mat';
file_path = '../../../data/data_binding_early_wood_27A/FRAP_002.mat';
raw_data = load(file_path);

% Extract image data.
image_data_pre_bleach = raw_data.imdata{1};
image_data_post_bleach = raw_data.imdata{3};

% Extract pixel size.
pixel_size = raw_data.imfo.frapPb1Series22.hardwaresettinglist.scannersettingrecord.variant{37};

% Extract bit depth.
bit_depth = raw_data.imfo.frapPb1Series22.hardwaresettinglist.scannersettingrecord.variant{61};

% Time lag between frames.
delta_t = raw_data.imfo.frapPb1Series22.contextdescription.block_frap_time_info.time{3};

% Extract bleaching profile information.
mask = raw_data.imdata{2}(:,:,1) == 2^16 - 1;
area = sum(mask(:)) * raw_data.imfo.frapBleachSeries20.hardwaresettinglist.scannersettingrecord.variant{37}^2;
r_bleach_SI = sqrt(area/pi);
r_bleach = r_bleach_SI / pixel_size;

clear file_path raw_data mask area

% bit_depth = 16;
% pixel_size = 7.5980e-07; % m.
% delta_t = 0.2650; % s.

number_of_pixels = size(image_data_post_bleach, 1);
number_of_post_bleach_images = 2;

x_bleach = 128;
y_bleach = 128;
% r_bleach = 0.5 * 50e-6 / pixel_size % The 0.5 => diameter -> radius

number_of_time_points_fine_per_coarse = [];%500;
number_of_pad_pixels = 128;
number_of_iterations = 1000;

%% Extract desired numbers of images/frames to include in analysis.
image_data_post_bleach = image_data_post_bleach(:, :, 1:number_of_post_bleach_images);

%% Convert image data to double and rescale to [0, 1] range.
image_data_pre_bleach = double(image_data_pre_bleach);
image_data_pre_bleach = image_data_pre_bleach / (2^bit_depth - 1);

image_data_post_bleach = double(image_data_post_bleach);
image_data_post_bleach = image_data_post_bleach / (2^bit_depth - 1);

%% Background subtraction.
image_data_post_bleach = subtract_background(image_data_pre_bleach, image_data_post_bleach);

%% Compute recovery curve.
[X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
X = X - 0.5;
Y = Y - 0.5;
ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
ind = ind(:);
recovery_curve = zeros(1, number_of_post_bleach_images);
for current_image_post_bleach = 1:number_of_post_bleach_images
    slice = image_data_post_bleach(:, :, current_image_post_bleach);
    recovery_curve(current_image_post_bleach) = mean(slice(ind));
end

%% Parameter estimation pre-work.

% Set parameter bounds for first estimation.
lb_1 = [0, min(image_data_post_bleach(:)), min(image_data_post_bleach(:))]; % mobile_fraction, intensity_inside_bleach_region, intensity_outside_bleach_region
ub_1 = [1, max(image_data_post_bleach(:)), max(image_data_post_bleach(:))];

% Initial guess for first estimation.
mobile_fraction_hat = (max(recovery_curve)-min(recovery_curve))/(mean(image_data_pre_bleach(:))-min(recovery_curve));
intensity_inside_bleach_region_hat = min(recovery_curve);
intensity_outside_bleach_region_hat = mean(image_data_pre_bleach(:));
param_hat_1 = [mobile_fraction_hat, intensity_inside_bleach_region_hat, intensity_outside_bleach_region_hat];

% Set parameter bounds for second estimation.
lb_2_SI = 1e-12;
ub_2_SI = 1e-9;
lb_2 = lb_2_SI / pixel_size^2;
ub_2 = ub_2_SI / pixel_size^2;

% Initial guess for second estimation.
param_hat_2_SI = 1e-10;
param_hat_2 = param_hat_2_SI / pixel_size^2;

%% Least-squares optimization.

options_1 = optimoptions(@lsqnonlin);
options_1.Algorithm = 'trust-region-reflective';
options_1.Display = 'iter';
options_1.FunctionTolerance = 1e-6;
options_1.OptimalityTolerance = 1e-6;
options_1.StepTolerance = 1e-6;
options_1.CheckGradients = false;%true;
options_1.SpecifyObjectiveGradient = true;

options_2 = optimoptions(@lsqnonlin);
options_2.Algorithm = 'trust-region-reflective';
options_2.Display = 'iter';
options_2.FunctionTolerance = 1e-6;
options_2.OptimalityTolerance = 1e-6;
options_2.StepTolerance = 1e-6;
options_2.MaxIterations = 1;
options_2.UseParallel = true;

param_hat = zeros(number_of_iterations, 4);
ss = zeros(number_of_iterations, 1);

for current_iteration = 1:number_of_iterations
    [image_data_post_bleach_model_unscaled, initial_condition_model_unscaled] = signal_diffusion(   param_hat_2, ...
                                                                                                    1.0, ...
                                                                                                    x_bleach, ...
                                                                                                    y_bleach, ...
                                                                                                    r_bleach, ...
                                                                                                    0.5, ...
                                                                                                    1.0, ...
                                                                                                    delta_t, ...
                                                                                                    number_of_time_points_fine_per_coarse, ...
                                                                                                    number_of_pixels, ...
                                                                                                    number_of_post_bleach_images, ...
                                                                                                    number_of_pad_pixels);
    fun_1 = @(param)residual_diffusion_partial( param(1), ...
                                                param(2), ...
                                                param(3), ...
                                                image_data_post_bleach, ...
                                                image_data_post_bleach_model_unscaled, ...
                                                initial_condition_model_unscaled);

    [param_hat_1, ss_1] = lsqnonlin(fun_1, param_hat_1, lb_1, ub_1, options_1);
    
    disp([param_hat_2 param_hat_1])
    
    fun_2 = @(param)residual_diffusion_full(param, ...
                                            param_hat_1(1), ...
                                            x_bleach, ...
                                            y_bleach, ...
                                            r_bleach, ...
                                            param_hat_1(2), ...
                                            param_hat_1(3), ...
                                            delta_t, ...
                                            number_of_time_points_fine_per_coarse, ...
                                            number_of_pad_pixels, ...
                                            image_data_post_bleach);

    [param_hat_2, ss_2] = lsqnonlin(fun_2, param_hat_2, lb_2, ub_2, options_2);
    
    param_hat(current_iteration, :) = [param_hat_2 param_hat_1];
    ss(current_iteration) = ss_2;
    
    disp([param_hat_2 param_hat_1])
end