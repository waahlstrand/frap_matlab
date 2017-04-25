%% Initialization.
clear
clc
close all hidden

addpath('../src/signal');
addpath('../src/signal_fourier');

%% Parameters.
D_SI = 5e-11; % m^2/s
pixel_size = 7.5e-07; % m
D = D_SI / pixel_size^2; % pixels^2 / s
k_on = 1e-1; % 1/s
k_off = 1e1; % 1/s
mf = 1.0; % dimensionless

delta_t = 0.2; % s
number_of_pixels = 256;
number_of_images = 1;
Ib = 0.6; % a.u.
Iu = 1.0;%0.9; % a.u.
x_bleach = 128; % pixels
y_bleach = 128; % pixels
r_bleach = 25e-6 / pixel_size; % pixels corresponding to 15 µm radius (30 µm diameter)
lx_bleach = r_bleach;
ly_bleach = r_bleach;

number_of_pad_pixels = 128;

% param_bleach = [x_bleach, y_bleach, r_bleach]; % Circular.
param_bleach = [x_bleach, y_bleach, lx_bleach, ly_bleach]; % Rectangular.

%% Simulate.
tic
data_fourier = signal_fourier_db( ...
    D, ...
    k_on, ...
    k_off, ...
    mf, ...
    Ib, ...
    Iu, ...
    param_bleach, ...
    delta_t, ...
    number_of_pixels, ...
    number_of_images);
toc

tic
data = signal_db( ...
    D, ...
    k_on, ...
    k_off, ...
    mf, ...
    Ib, ...
    Iu, ...
    param_bleach, ...
    delta_t, ...
    number_of_pixels, ...
    number_of_images, ...
    number_of_pad_pixels);
toc

for i = 1:number_of_images
    data_f(:, :, i) = fft2(data(:, :, i));
end


%% Plot solution.
figure
hold on
imagesc(reshape(real(data_fourier), [number_of_pixels, number_of_pixels * number_of_images]))
% axis 'equal'
axis([0 number_of_images*number_of_pixels 0 number_of_pixels])
axis off
hold off

figure
hold on
imagesc(reshape(abs(data_f), [number_of_pixels, number_of_pixels * number_of_images]))
% axis 'equal'
axis([0 number_of_images*number_of_pixels 0 number_of_pixels])
axis off
hold off