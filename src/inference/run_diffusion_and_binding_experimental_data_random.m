%% Initialization.
clear
clc
close all hidden

%% Random stream.
random_seed = round( sum( 1e6 * clock() ) );
random_stream = RandStream('mt19937ar', 'Seed', random_seed);
RandStream.setGlobalStream(random_stream);

%% Load data.

folder = randsample([6, 27], 1);
file = randsample(2:7, 1);

file_path = ['../../../data/data_binding_early_wood_' num2str(folder) 'A/FRAP_00' num2str(file) '.mat'];
raw_data = load(file_path);

% Extract image data.
image_data_pre_bleach = raw_data.imdata{1};
image_data_post_bleach = raw_data.imdata{3};

% Extract pixel size, bit depth, and time lag between frames. Start by
% finding the name of the structure field that's called something with
% 'frapPb1Series...'.
list_field_names = fieldnames(raw_data.imfo);
ind_field = find(~cellfun(@isempty, strfind(list_field_names, 'frapPb1Series')));

pixel_size = raw_data.imfo.(list_field_names{ind_field}).hardwaresettinglist.scannersettingrecord.variant{37};
bit_depth = raw_data.imfo.(list_field_names{ind_field}).hardwaresettinglist.scannersettingrecord.variant{61};
delta_t = raw_data.imfo.(list_field_names{ind_field}).contextdescription.block_frap_time_info.time{3};

% Extract bleaching profile information (radius). Start by finding the name 
% of the structure field that's called something with 'frapBleachSeries...'.
ind_field = find(~cellfun(@isempty, strfind(list_field_names, 'frapBleachSeries')));
pixel_size_bleach = raw_data.imfo.(list_field_names{ind_field}).hardwaresettinglist.scannersettingrecord.variant{37};
mask = raw_data.imdata{2}(:,:,1) == 2^16 - 1;
area = sum(mask(:)) * pixel_size_bleach^2;
r_bleach_SI = sqrt(area/pi);
r_bleach = r_bleach_SI / pixel_size;

clear file_path raw_data mask area r_bleach_SI list_field_names ind_field

number_of_pixels = size(image_data_post_bleach, 1);
number_of_post_bleach_images = 100;

x_bleach = 128;
y_bleach = 128;

number_of_time_points_fine_per_coarse = 500;
number_of_pad_pixels = 128;

%% Extract desired numbers of images/frames to include in analysis.
image_data_post_bleach = image_data_post_bleach(:, :, 1:number_of_post_bleach_images);

%% Convert image data to double and rescale to [0, 1] range.
image_data_pre_bleach = double(image_data_pre_bleach);
image_data_pre_bleach = image_data_pre_bleach / (2^bit_depth - 1);

image_data_post_bleach = double(image_data_post_bleach);
image_data_post_bleach = image_data_post_bleach / (2^bit_depth - 1);

%% Background subtraction.
image_data_post_bleach = subtract_background(image_data_pre_bleach, image_data_post_bleach);

%% Parameter estimation pre-work.

% Set parameter bounds for first estimation.
lb_1 = [0, min(image_data_post_bleach(:)), min(image_data_post_bleach(:))];
ub_1 = [1, max(image_data_post_bleach(:)), max(image_data_post_bleach(:))];

% Initial guess for first estimation.
param_hat_1 = lb_1 + (ub_1 - lb_1) .* rand(size(lb_1));
param_hat_1(2:3) = sort(param_hat_1(2:3), 'ascend'); % Make sure the two intensity levels are not switched.

% Set parameter bounds for second estimation.
lb_2_SI = [1e-12, 0, 0];
ub_2_SI = [1e-9, 50, 50];
lb_2 = lb_2_SI;
lb_2(1) = lb_2(1) / pixel_size^2;
ub_2 = ub_2_SI;
ub_2(1) = ub_2(1) / pixel_size^2;

% Initial guess for second estimation.
param_hat_2 = lb_2 + (ub_2 - lb_2) .* rand(size(lb_2));

%% Least-squares optimization.

options_1 = optimoptions(@lsqnonlin);
options_1.Algorithm = 'trust-region-reflective';
options_1.Display = 'iter';
options_1.FunctionTolerance = 1e-6;
options_1.OptimalityTolerance = 1e-6;
options_1.StepTolerance = 1e-6;
options_1.CheckGradients = false;
options_1.SpecifyObjectiveGradient = true;

options_2 = optimoptions(@lsqnonlin);
options_2.Algorithm = 'trust-region-reflective';
options_2.Display = 'iter';
options_2.FunctionTolerance = 1e-6;
options_2.OptimalityTolerance = 1e-6;
options_2.StepTolerance = 1e-6;
options_2.MaxIterations = 1;

param_hat = [];
ss = [];

param_hat = [param_hat ; [param_hat_2 param_hat_1]];
ss = [ss ; Inf];

is_converged = false;
while ~is_converged
    [image_data_post_bleach_model_unscaled, initial_condition_model_unscaled] = signal_diffusion_and_binding(param_hat_2(1), ...
                                                                                        param_hat_2(2), ...
                                                                                        param_hat_2(3), ...
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
    fun_1 = @(param)residual_diffusion_and_binding_partial( param(1), ...
                                                            param(2), ...
                                                            param(3), ...
                                                            image_data_post_bleach, ...
                                                            image_data_post_bleach_model_unscaled, ...
                                                            initial_condition_model_unscaled);

    [param_hat_1, ss_1] = lsqnonlin(fun_1, param_hat_1, lb_1, ub_1, options_1);
    
    disp([param_hat_2 param_hat_1])
    
    fun_2 = @(param)residual_diffusion_and_binding_full(param(1), ...
                                                        param(2), ...
                                                        param(3), ...
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
    
    param_hat = [param_hat ; [param_hat_2 param_hat_1]];
    ss = [ss ; ss_2];
    
    disp([param_hat_2 param_hat_1])
    
    if size(param_hat, 1) >= 2
        if max(abs(param_hat(end, :) - param_hat(end - 1, :))) < 1e-3
            is_converged = true;
        end
    end
            
end

save(['est_' num2str(folder) '_' num2str(file) '_' num2str(random_seed) '.mat'], 'param_hat', 'ss'); 