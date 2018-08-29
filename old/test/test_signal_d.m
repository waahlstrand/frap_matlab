%% Initialization.
clear
clc
close all hidden

addpath('../src/signal');

%% Parameters.
D_SI = 5e-11; % m^2/s
pixel_size = 7.5e-07; % m
D = D_SI / pixel_size^2; % pixels^2 / s
mf = 1.0; % dimensionless

delta_t = 0.2; % s
number_of_pixels = 256;
number_of_images = 100;
number_of_pad_pixels = 128;
Ib = 0.6; % a.u.
Iu = 1.0;%0.9; % a.u.
x_bleach = 128; % pixels
y_bleach = 128; % pixels
r_bleach = 20e-6 / pixel_size; % pixels corresponding to 15 �m radius (30 �m diameter)
lx_bleach = r_bleach;
ly_bleach = r_bleach;

% param_bleach = [x_bleach, y_bleach, r_bleach]; % Circular.
param_bleach = [x_bleach, y_bleach, lx_bleach, ly_bleach]; % Rectangular.

%% Simulate.
tic
data = signal_d( ...
    D, ...
    mf, ...
    Ib, ...
    Iu, ...
    param_bleach, ...
    delta_t, ...
    number_of_pixels, ...
    number_of_images, ...
    number_of_pad_pixels);
toc
return
sigma_noise = 0;
data = data + sigma_noise * randn(size(data));

%% Plot solution.
figure
hold on
imagesc(reshape(data, [number_of_pixels, number_of_pixels * number_of_images]))
% axis 'equal'
axis([0 number_of_images*number_of_pixels 0 number_of_pixels])
axis off
hold off

%% Plot recovery curve.
[X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
X = X - 0.5;
Y = Y - 0.5;

if numel(param_bleach) == 3 % Circular.
    ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
else % Rectangular.
    ind = find( X >= x_bleach - 0.5 * lx_bleach & X <= x_bleach + 0.5 * lx_bleach & Y >= y_bleach - 0.5 * ly_bleach & Y <= y_bleach + 0.5 * ly_bleach );
end
ind = ind(:);

recovery_curve = zeros(1, number_of_images);
for current_image = 1:number_of_images
    slice = data(:, :, current_image);
    recovery_curve(current_image) = mean(slice(ind));
end
figure
plot(delta_t:delta_t:number_of_images*delta_t, recovery_curve)
