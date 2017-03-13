%% Initialization.
clear
clc
close all hidden

delete(gcp('nocreate'))
c = parcluster('local');
c.NumWorkers = 32;
parpool(c, c.NumWorkers);

%% Experiment parameters.
delta_t = 0.25; % s
number_of_post_bleach_images = 40;
number_of_pixels = 256;
number_of_pad_pixels = 128;
r_bleach = 32; % pixels
x_bleach = 128;
y_bleach = 128;
intensity_inside_bleach_region = 0.6;
intensity_outside_bleach_region = 0.9;

number_of_time_points_fine_per_coarse = [];

%% Particle parameters.
D = 200; % pixels^2 / s
k_on = 0.1; % 1/s
k_off = 1.0; % 1/s
mobile_fraction = 1.0;%0.8;

%% Simulate true data set.
image_data_post_bleach_true = signal_diffusion_and_binding(  D, ...
                                                        k_on, ...
                                                        k_off, ...
                                                        mobile_fraction, ...
                                                        x_bleach, ...
                                                        y_bleach, ...
                                                        r_bleach, ...
                                                        intensity_inside_bleach_region, ...
                                                        intensity_outside_bleach_region, ...
                                                        delta_t, ...
                                                        number_of_time_points_fine_per_coarse, ...
                                                        number_of_pixels, ...
                                                        number_of_post_bleach_images, ...
                                                        number_of_pad_pixels);
                                                    
                                                    
[X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
X = X - 0.5;
Y = Y - 0.5;
ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
ind = ind(:);
recovery_curve_true = zeros(1, number_of_post_bleach_images);
for current_image_post_bleach = 1:number_of_post_bleach_images
    slice = image_data_post_bleach_true(:, :, current_image_post_bleach);
    recovery_curve_true(current_image_post_bleach) = mean(slice(ind));
end                                    

%% Sweep parameters.
number_of_grid_points = 4+3*3;
[K_ON, K_OFF] = ndgrid(logspace(-2, 1, number_of_grid_points), logspace(-2, 1, number_of_grid_points));
SS_pixel_based = zeros(number_of_grid_points, number_of_grid_points);
SS_curve_based = zeros(number_of_grid_points, number_of_grid_points);

parfor i = 1:number_of_grid_points^2
    disp([K_ON(i) K_OFF(i)])
    image_data_post_bleach_model = signal_diffusion_and_binding(D, ...
                                                                K_ON(i), ...
                                                                K_OFF(i), ...
                                                                mobile_fraction, ...
                                                                x_bleach, ...
                                                                y_bleach, ...
                                                                r_bleach, ...
                                                                intensity_inside_bleach_region, ...
                                                                intensity_outside_bleach_region, ...
                                                                delta_t, ...
                                                                number_of_time_points_fine_per_coarse, ...
                                                                number_of_pixels, ...
                                                                number_of_post_bleach_images, ...
                                                                number_of_pad_pixels);
    
    SS_pixel_based(i) = sum( ( image_data_post_bleach_model(:) - image_data_post_bleach_true(:) ).^2 );
    
    recovery_curve_model = zeros(1, number_of_post_bleach_images);
    for current_image_post_bleach = 1:number_of_post_bleach_images
        slice = image_data_post_bleach_model(:, :, current_image_post_bleach);
        recovery_curve_model(current_image_post_bleach) = mean(slice(ind));
    end
    
    SS_curve_based(i) = sum( ( recovery_curve_model - recovery_curve_true ).^2 ); 
    
end
   
figure, hold on
surf(log10(K_ON), log10(K_OFF), log10(SS_pixel_based))

figure, hold on
surf(log10(K_ON), log10(K_OFF), log10(SS_curve_based))