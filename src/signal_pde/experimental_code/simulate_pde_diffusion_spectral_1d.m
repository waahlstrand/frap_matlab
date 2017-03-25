%% Initialization.
clear
clc
close all hidden

%% Measurement parameters.
delta_t = 0.25; % s
number_of_post_bleach_images = 200;
number_of_pixels = 1024%256;
number_of_pad_pixels = 0;
r_bleach_region = 128%32; % pixels

intensity_inside_bleach_region = 0.6;
intensity_outside_bleach_region = 0.9;

%% Particle parameters.
D = 400; % pixels^2 / s
k_on = 0.05; % 1/s
k_off = 0.05; % 1/s

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

F_U0 = fft(U0);
F_B0 = fft(B0);

%% FFT space time evolution of PDE system.

F_image_data_post_bleach_u = zeros(number_of_pixels, number_of_post_bleach_images);
F_image_data_post_bleach_b = zeros(number_of_pixels, number_of_post_bleach_images);
image_data_post_bleach_u = zeros(number_of_pixels, number_of_post_bleach_images);
image_data_post_bleach_b = zeros(number_of_pixels, number_of_post_bleach_images);

XSI = 0:number_of_pixels-1;
XSISQ = XSI.^2;

T = delta_t * (1:number_of_post_bleach_images);
% 
% for i = 1:number_of_post_bleach_images
%     F_image_data_post_bleach_u(:, :, i) = exp( ( - D * XSISQ - k_on) * T(i) ) .* F_U0 + exp( k_off * T(i) ) .* F_B0;
%     F_image_data_post_bleach_b(:, :, i) = exp( k_on * T(i) ) .* F_U0 + exp( - k_off * T(i) ) .* F_B0;
% %     F_image_data_post_bleach_u(:, :, i) = exp( ( - XSISQ - k_on) .* T(i) .* F_U0 ) + exp( k_off * T(i) * F_B0 ); 
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
        A = [- D * XSISQ(i) - k_on, k_off ; k_on, - k_off];
%         A = [- k_on, k_off ; k_on, - k_off];
        c_vector_hat = exp( A * T(t) ) * [F_U0(i) ; F_B0(i)];
        
        F_image_data_post_bleach_u(i, t) = c_vector_hat(1);
        F_image_data_post_bleach_b(i, t) = c_vector_hat(2);
    end
end

for i = 1:number_of_post_bleach_images
    image_data_post_bleach_u(:, i) = abs(ifft(F_image_data_post_bleach_u(:, i)));
    image_data_post_bleach_b(:, i) = abs(ifft(F_image_data_post_bleach_b(:, i)));
end

FRAP = image_data_post_bleach_u + image_data_post_bleach_b;
semilogy(X, FRAP)

% figure, imagesc(abs(image_data_post_bleach_u(:,:,i))+abs(image_data_post_bleach_b(:,:,i)))

% imagesc(reshape(FRAP, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]))
% axis 'equal'



% 
