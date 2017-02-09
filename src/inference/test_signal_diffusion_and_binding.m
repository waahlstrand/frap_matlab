%% Initialization.
clear
clc
close all hidden

%% Measurement parameters.
delta_t = 0.250; % s.
number_of_pixels = 256;
number_of_post_bleach_images = 1%200;
xc = number_of_pixels / 2;
yc = number_of_pixels / 2;
r_bleach_region = 32; % pixels
intensity_inside_bleach_region = 0.6;
intensity_outside_bleach_region = 0.85;

sigma_noise = 0.0;

%% Particle parameters.
D_SI = 2.5e-10; % m^2/s
pixel_size = 7.5e-07; % m.
D = D_SI / pixel_size^2; % pixels^2 / s
k_on = 1; % 1/s
k_off = 1; % 1/s
mobile_fraction = 0.90;%0.8;

% profile on
tic
image_data_post_bleach = signal(D, k_on, k_off, mobile_fraction, xc, yc, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region, delta_t, number_of_pixels, number_of_post_bleach_images);
toc
% profile off
% profile viewer
return

image_data_post_bleach = image_data_post_bleach + sigma_noise * randn(size(image_data_post_bleach));

%% Save data.
% save('simulated_data_zero_noise.mat', 'image_data_post_bleach', 'delta_t', 'pixel_size', 'D_SI', 'D', 'k_on', 'k_off', 'mobile_fraction', 'xc', 'yc', 'r_bleach_region', 'intensity_inside_bleach_region', 'intensity_outside_bleach_region');

%% Plot.
result_pde = [];
for current_image_post_bleach = 1:number_of_post_bleach_images
    result_pde = [result_pde, image_data_post_bleach(:, :, current_image_post_bleach)];
end
figure
imagesc(result_pde)
axis 'equal'
axis([0 number_of_post_bleach_images*number_of_pixels 0 number_of_pixels])
axis off

[X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
X = X - 0.5;
Y = Y - 0.5;
ind = find( (X - xc).^2 + (Y - yc).^2 <= r_bleach_region^2 );
ind = ind(:);
recovery_curve = zeros(1, number_of_post_bleach_images);
for current_image_post_bleach = 1:number_of_post_bleach_images
    slice = image_data_post_bleach(:, :, current_image_post_bleach);
    recovery_curve(current_image_post_bleach) = mean(slice(ind));
end
figure
plot(delta_t:delta_t:number_of_post_bleach_images*delta_t, recovery_curve)

