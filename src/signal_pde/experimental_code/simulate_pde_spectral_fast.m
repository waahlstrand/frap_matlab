%% Initialization.
clear
clc
close all hidden

%% Measurement parameters.
delta_t = 0.25; % s
number_of_post_bleach_images = 10;
number_of_pixels = 256;
number_of_pad_pixels = 128;
r_bleach_region = 32; % pixels

intensity_inside_bleach_region = 0.6;
intensity_outside_bleach_region = 0.9;

%% Particle parameters.
D = 4000%1200; % pixels^2 / s
k_on = 0.2; % 1/s
k_off = 3.0; % 1/s

p_free = k_off / ( k_on + k_off );
p_bound = k_on / ( k_on + k_off );

%% Initial condition. 
% Create a high resolution initial condition which is then downsampled to 
% avoid too sharp edges.

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

clear X Y

%% FFT of initial conditions.

F_U0 = fftshift(fft2(U0));
F_B0 = fftshift(fft2(B0));

% figure, imagesc(log10(abs(F_U0)))
% figure, imagesc(log10(abs(F_B0)))

% return

%% FFT space time evolution of PDE system.

F_image_data_post_bleach_u = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
F_image_data_post_bleach_b = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
image_data_post_bleach_u = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
image_data_post_bleach_b = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);

[XSI1, XSI2] = meshgrid(-(number_of_pixels + 2 * number_of_pad_pixels)/2:(number_of_pixels + 2 * number_of_pad_pixels)/2-1, ...
                        -(number_of_pixels + 2 * number_of_pad_pixels)/2:(number_of_pixels + 2 * number_of_pad_pixels)/2-1);
XSI1 = XSI1 / (number_of_pixels + 2 * number_of_pad_pixels);
XSI2 = XSI2 / (number_of_pixels + 2 * number_of_pad_pixels);
XSISQ = XSI1.^2 + XSI2.^2;

% Diagonalization, excluding time t which DD will be multiplied with.
PP11 = -(k_on - k_off + (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*k_on);
PP12 = -(k_on - k_off - (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*k_on);
PP21 = 1;
PP22 = 1;

DD11 = -((k_on + k_off + (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ))./2;
DD22 = -((k_on + k_off - (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ))./2;

PPinv11 = -k_on./(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2);
PPinv12 = -(k_on - k_off - (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2));
PPinv21 = k_on./(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2);
PPinv22 = (k_on - k_off + (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2));
% 
% 
% tic
% for t = 1:number_of_post_bleach_images
% %     disp(t)
%     T = t * delta_t;
%     F_image_data_post_bleach_u(:, :, t) = (PP11 .* PPinv11 .* exp(DD11 * T) + PP12 .* PPinv21 .* exp(DD22 * T)) .* F_U0 + (PP11 .* PPinv12 .* exp(DD11 * T) + PP12 .* PPinv22 .* exp(DD22 * T)) .* F_B0;
%     F_image_data_post_bleach_b(:, :, t) = (PP21 .* PPinv11 .* exp(DD11 * T) + PP22 .* PPinv21 .* exp(DD22 * T)) .* F_U0 + (PP21 .* PPinv12 .* exp(DD11 * T) + PP22 .* PPinv22 .* exp(DD22 * T)) .* F_B0;
% end
% toc

tic
CONST11 = PP11 .* (PPinv11 .* F_U0 + PPinv12 .* F_B0);
CONST12 = PP12 .* (PPinv21 .* F_U0 + PPinv22 .* F_B0);
CONST21 = PP21 .* (PPinv11 .* F_U0 + PPinv12 .* F_B0);
CONST22 = PP22 .* (PPinv21 .* F_U0 + PPinv22 .* F_B0);
for t = 1:number_of_post_bleach_images
%     disp(t)
    T = t * delta_t;
    F_image_data_post_bleach_u(:, :, t) = CONST11 .* exp(DD11 * T) + CONST12 .* exp(DD22 * T);
    F_image_data_post_bleach_b(:, :, t) = CONST21 .* exp(DD11 * T) + CONST22 .* exp(DD22 * T);
end
toc
   
for t = 1:number_of_post_bleach_images
    disp(t)
    image_data_post_bleach_u(:, :, t) = abs(ifft2(ifftshift(F_image_data_post_bleach_u(:, :, t))));
    image_data_post_bleach_b(:, :, t) = abs(ifft2(ifftshift(F_image_data_post_bleach_b(:, :, t))));
end

image_data_post_bleach_u = image_data_post_bleach_u(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels, :);
image_data_post_bleach_b = image_data_post_bleach_b(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels, :);

image_data_post_bleach = image_data_post_bleach_u + image_data_post_bleach_b;

imagesc(reshape(image_data_post_bleach, [number_of_pixels, number_of_pixels*number_of_post_bleach_images]))
axis 'equal'