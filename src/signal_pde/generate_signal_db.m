%% Initialization.
clear
clc
close all hidden

%% Parameters.
D_SI = 2.5e-10; % m^2/s
pixel_size = 7.598e-07; % m
D = D_SI / pixel_size^2; % pixels^2 / s
k_on = 1.0; % 1/s
k_off = 1.0; % 1/s
mobile_fraction = 1.0; % dimensionless

delta_t = 0.2650; % s
number_of_pixels = 256;
number_of_images = 100;
number_of_pad_pixels = 128;
intensity_inside_bleach_region = 0.6; % a.u.
intensity_outside_bleach_region = 0.9; % a.u.
x_bleach = 128; % pixels
y_bleach = 128; % pixels
r_bleach = 32; % pixels

%% Simulate.
tic
signal = signal_db( D, ...
                    k_on, ...
                    k_off, ...
                    mobile_fraction, ...
                    intensity_inside_bleach_region, ...
                    intensity_outside_bleach_region, ...
                    x_bleach, ...
                    y_bleach, ...
                    r_bleach, ...
                    delta_t, ...
                    number_of_pixels, ...
                    number_of_images, ...
                    number_of_pad_pixels);
toc

%% Plot solution.
figure
hold on
imagesc(reshape(signal, [number_of_pixels, number_of_pixels * number_of_images]))
axis 'equal'
axis([0 number_of_images*number_of_pixels 0 number_of_pixels])
axis off
hold off

%% Plot recovery curve.
[X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
X = X - 0.5;
Y = Y - 0.5;
ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
ind = ind(:);
recovery_curve = zeros(1, number_of_images);
for current_image = 1:number_of_images
    slice = signal(:, :, current_image);
    recovery_curve(current_image) = mean(slice(ind));
end
figure
plot(delta_t:delta_t:number_of_images*delta_t, recovery_curve)

%% Save data.
% clear current_image X Y ind recovery_curve slice
% save('simulated_signal_db.mat')
