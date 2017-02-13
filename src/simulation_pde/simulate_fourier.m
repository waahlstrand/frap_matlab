%% Initialization.
clear
clc
close all hidden

%% Measurement parameters.
delta_t = 0.25; % s
number_of_post_bleach_images = 20;
number_of_pixels = 256;
number_of_pad_pixels = 256;
r_bleach_region = 32; % pixels

intensity_inside_bleach_region = 0.6;
intensity_outside_bleach_region = 0.9;

%% Particle parameters.
D = 400; % pixels^2 / s
k_on = 0%0.05; % 1/s
k_off = 0.05; % 1/s
mobile_fraction = 0.5;%0.8;

p_free_marginal = k_off / ( k_on + k_off );
p_bound_marginal = k_on / ( k_on + k_off );

%% Numerical Fourier transform of initial condition.
upsampling_factor = 3;
[X, Y] = meshgrid(1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels), 1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels));
X = X - 0.5;
Y = Y - 0.5;
xc = number_of_pad_pixels + number_of_pixels / 2;
yc = number_of_pad_pixels + number_of_pixels / 2;
U0 = zeros(size(X));
U0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 <= (upsampling_factor * r_bleach_region)^2 ) = intensity_inside_bleach_region;
U0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 > (upsampling_factor * r_bleach_region)^2 ) = intensity_outside_bleach_region;
clear X Y
% Distribute bound and free according to their marginal distribution,
% defined by p_free and p_bound.
B0 = p_bound_marginal * U0;
U0 = p_free_marginal * U0;

B0 = imresize(B0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);
U0 = imresize(U0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);

% figure, imagesc(U0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels)), axis 'equal'



% return
