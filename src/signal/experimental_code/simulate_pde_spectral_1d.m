%% Initialization.
clear
clc
close all hidden

%% Measurement parameters.
delta_t = 0.25; % s
number_of_post_bleach_images = 100%200;
number_of_pixels = 256;
number_of_pad_pixels = 128;
r_bleach_region = 32; % pixels

intensity_inside_bleach_region = 0.6;
intensity_outside_bleach_region = 0.9;

%% Particle parameters.
D = 700; % pixels^2 / s
k_on = 0.2; % 1/s
k_off = 3.0; % 1/s

p_free = k_off / ( k_on + k_off );
p_bound = k_on / ( k_on + k_off );

%% Generate post bleach image data for pure diffusion.

% Initial condition. Create a high resolution initial condition which is 
% then downsampled to avoid too sharp edges.


X = 1:(number_of_pixels + 2 * number_of_pad_pixels);
X = X - 0.5;
xc = number_of_pad_pixels + number_of_pixels / 2;

U0 = zeros(size(X));
U0( (X - xc).^2 <= ( r_bleach_region)^2 ) = intensity_inside_bleach_region;
U0( (X - xc).^2 > ( r_bleach_region)^2 ) = intensity_outside_bleach_region;

% sigma = 7.5;
% U0 = intensity_outside_bleach_region - (intensity_outside_bleach_region - intensity_inside_bleach_region) * exp( -(X - xc).^2 / (2*sigma^2));

% Distribute bound and free according to their marginal distribution,
% defined by p_free and p_bound.
B0 = p_bound * U0;
U0 = p_free * U0;

% figure, plot(X, U0+B0)
% return

% clear X


%% FFT of initial conditions.

F_U0 = fftshift(fft(U0));
F_B0 = fftshift(fft(B0));

%% FFT space time evolution of PDE system.

F_image_data_post_bleach_u = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
F_image_data_post_bleach_b = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
image_data_post_bleach_u = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
image_data_post_bleach_b = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);

XSI = -(number_of_pixels + 2 * number_of_pad_pixels)/2:(number_of_pixels + 2 * number_of_pad_pixels)/2-1;
XSI = XSI / (number_of_pixels + 2 * number_of_pad_pixels);
XSISQ = XSI.^2;

T = delta_t * (1:number_of_post_bleach_images);

for t = 1:number_of_post_bleach_images
    for i = 1:(number_of_pixels + 2 * number_of_pad_pixels)
        disp([t, i])
        A = [- D * XSISQ(i) - k_on, k_off ; k_on, - k_off];
        c_vector_hat = expm( A * T(t) ) * [F_U0(i) ; F_B0(i)];
        
        F_image_data_post_bleach_u(i, t) = c_vector_hat(1);
        F_image_data_post_bleach_b(i, t) = c_vector_hat(2);
    end
end

for i = 1:number_of_post_bleach_images
    image_data_post_bleach_u(:, i) = abs(ifft(ifftshift(F_image_data_post_bleach_u(:, i))));
    image_data_post_bleach_b(:, i) = abs(ifft(ifftshift(F_image_data_post_bleach_b(:, i))));
end

image_data_post_bleach_u = image_data_post_bleach_u(number_of_pad_pixels+1:end-number_of_pad_pixels, :);
image_data_post_bleach_b = image_data_post_bleach_b(number_of_pad_pixels+1:end-number_of_pad_pixels, :);

FRAP = image_data_post_bleach_u + image_data_post_bleach_b;
semilogy(FRAP)

% figure, imagesc(abs(image_data_post_bleach_u(:,:,i))+abs(image_data_post_bleach_b(:,:,i)))

% imagesc(reshape(FRAP, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]))
% axis 'equal'



% 
