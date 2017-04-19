%% Initialization.
clear
clc
close all hidden

%% Measurement parameters.
delta_t = 0.25; % s
number_of_post_bleach_images = 100;
number_of_pixels = 256;
number_of_pad_pixels = 0;
r_bleach_region = 32; % pixels

intensity_inside_bleach_region = 0.6;
intensity_outside_bleach_region = 0.9;

%% Particle parameters.
D = 400; % pixels^2 / s

%% Generate post bleach image data for pure diffusion.

% Initial condition. Create a high resolution initial condition which is
% then downsampled to avoid too sharp edges.

X = 1:(number_of_pixels + 2 * number_of_pad_pixels);
X = X - 0.5;
xc = number_of_pad_pixels + number_of_pixels / 2;

U0 = zeros(size(X));
U0( (X - xc).^2 <= ( r_bleach_region)^2 ) = intensity_inside_bleach_region;
U0( (X - xc).^2 > ( r_bleach_region)^2 ) = intensity_outside_bleach_region;

% sigma = 15;
% U0 = intensity_outside_bleach_region - (intensity_outside_bleach_region - intensity_inside_bleach_region) * exp( -(X - xc).^2 / (2*sigma^2));

% figure, plot(X, U0+B0)
% return

%% FFT of initial conditions.

F_U0 = fftshift(fft(U0));


%% FFT space time evolution of PDE system.

F_image_data_post_bleach_u = zeros(number_of_pixels, number_of_post_bleach_images);
image_data_post_bleach_u = zeros(number_of_pixels, number_of_post_bleach_images);

% XSI = (0:number_of_pixels-1)/number_of_pixels;
% XSISQ = XSI.^2;
XSI = -number_of_pixels/2:number_of_pixels/2-1;
XSI = XSI / number_of_pixels;
XSISQ = XSI.^2;

T = delta_t * (1:number_of_post_bleach_images);

for t = 1:number_of_post_bleach_images
    disp(t)
    A = - D * XSISQ;
    F_U = exp( A * T(t) ) .* F_U0;
    
    F_image_data_post_bleach_u(:, t) = F_U;
end

for i = 1:number_of_post_bleach_images
    image_data_post_bleach_u(:, i) = abs(ifft(ifftshift(F_image_data_post_bleach_u(:, i))));
end

figure, semilogy(X, abs(F_image_data_post_bleach_u))
figure, semilogy(X, abs(image_data_post_bleach_u))