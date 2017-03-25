%% Initialization.
clear
clc
close all hidden

%% Measurement parameters.
delta_t = 0.25; % s
number_of_post_bleach_images = 5;
number_of_pixels = 256;
number_of_pad_pixels = 0;
r_bleach_region = 32; % pixels

intensity_inside_bleach_region = 0.6;
intensity_outside_bleach_region = 0.9;

%% Particle parameters.
D = 400; % pixels^2 / s
k_on = 0.5; % 1/s
k_off = 1.0; % 1/s
mobile_fraction = 1.0;%0.8;

p_free = k_off / ( k_on + k_off );
p_bound = k_on / ( k_on + k_off );

%% Generate post bleach image data

% Initial condition. Create a high resolution initial condition which is 
% then downsampled to avoid too sharp edges.

upsampling_factor = 3;


[X, Y] = meshgrid(1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels), 1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels));
X = X - 0.5;
Y = Y - 0.5;
xc = number_of_pad_pixels + number_of_pixels / 2;
yc = number_of_pad_pixels + number_of_pixels / 2;

U0 = zeros(size(X));
U0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 <= (upsampling_factor * r_bleach_region)^2 ) = intensity_inside_bleach_region;
U0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 > (upsampling_factor * r_bleach_region)^2 ) = intensity_outside_bleach_region;

% Distribute bound and free according to their marginal distribution,
% defined by p_free and p_bound.
B0 = p_bound * U0;
U0 = p_free * U0;

B0 = imresize(B0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);
U0 = imresize(U0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);

% figure, imagesc(U0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels)), axis 'equal'

clear X Y

%% FFT of initial conditions.

F_U0 = fft2(U0);
% F_U0 = fftshift(F_U0);

F_B0 = fft2(B0);
% F_B0 = fftshift(F_B0);

%% FFT space time evolution of PDE system.

F_image_data_post_bleach_u = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);
F_image_data_post_bleach_b = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);
image_data_post_bleach_u = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);
image_data_post_bleach_b = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);

% [XSI1, XSI2] = meshgrid(-number_of_pixels/2:number_of_pixels/2-1, -number_of_pixels/2:number_of_pixels/2-1);
% XSISQ = XSI1.^2 + XSI2.^2;
[XSI1, XSI2] = meshgrid(0:number_of_pixels-1, 0:number_of_pixels-1);
XSISQ = XSI1.^2 + XSI2.^2;

% For single point.
% i = 1;
% j = 1;
% F_U = zeros(number_of_post_bleach_images, 1);
% F_B = zeros(number_of_post_bleach_images, 1);
T = delta_t * (1:number_of_post_bleach_images);
% 
% for i = 1:number_of_post_bleach_images
%     F_image_data_post_bleach_u(:, :, i) = exp( ( - D * XSISQ - k_on) * T(i) ) .* F_U0 + exp( k_off * T(i) ) .* F_B0;
%     F_image_data_post_bleach_b(:, :, i) = exp( k_on * T(i) ) .* F_U0 + exp( - k_off * T(i) ) .* F_B0;
% %     F_image_data_post_bleach_u(:, :, i) = exp( ( - D * XSISQ - k_on) .* T(i) .* F_U0 ) + exp( k_off * T(i) * F_B0 ); 
% %     F_image_data_post_bleach_b(:, :, i) = exp( k_on * T(i) * F_U0 ) + exp( - k_off * T(i) * F_B0 ); 
% 
%     image_data_post_bleach_u(:, :, i) = ifftshift(F_image_data_post_bleach_u(:, :, i));
%     image_data_post_bleach_u(:, :, i) = ifft2(image_data_post_bleach_u(:, :, i));
%     image_data_post_bleach_b(:, :, i) = ifftshift(F_image_data_post_bleach_b(:, :, i));
%     image_data_post_bleach_b(:, :, i) = ifft2(image_data_post_bleach_b(:, :, i));
% end
% 
% FRAP = abs(image_data_post_bleach_u) + abs(image_data_post_bleach_b);
% % figure, imagesc(abs(image_data_post_bleach_u(:,:,i))+abs(image_data_post_bleach_b(:,:,i)))
% 
% imagesc(reshape(FRAP, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]))
% axis 'equal'


for t = 1:number_of_post_bleach_images
    for i = 1:number_of_pixels
        disp([t, i])
        for j = 1:number_of_pixels
            A = [- D * XSISQ(i, j) - k_on, k_off ; k_on, - k_off];
            c_vector_hat = exp( A * T(t) ) * [F_U0(i, j) ; F_B0(i, j)];
            
            F_image_data_post_bleach_u(i, j, t) = c_vector_hat(1);
            F_image_data_post_bleach_b(i, j, t) = c_vector_hat(2);
        end
    end
end

% F_FRAP = abs(F_image_data_post_bleach_u) + abs(F_image_data_post_bleach_b);
% imagesc(reshape(log10(F_FRAP), [number_of_pixels, number_of_pixels*number_of_post_bleach_images]))
% axis 'equal'

for i = 1:number_of_post_bleach_images
%     image_data_post_bleach_u(:, :, i) = ifftshift(F_image_data_post_bleach_u(:, :, i));
    image_data_post_bleach_u(:, :, i) = abs(ifft2(F_image_data_post_bleach_u(:, :, i)));
%     image_data_post_bleach_b(:, :, i) = ifftshift(F_image_data_post_bleach_b(:, :, i));
    image_data_post_bleach_b(:, :, i) = abs(ifft2(F_image_data_post_bleach_b(:, :, i)));
end

FRAP = image_data_post_bleach_u + image_data_post_bleach_b;
% figure, imagesc(abs(image_data_post_bleach_u(:,:,i))+abs(image_data_post_bleach_b(:,:,i)))

imagesc(reshape(FRAP, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]))
axis 'equal'



% 
