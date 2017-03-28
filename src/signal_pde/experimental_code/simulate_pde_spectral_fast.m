%% Initialization.
clear
clc
close all hidden

%% Measurement parameters.
delta_t = 0.25; % s
number_of_post_bleach_images = 100;
number_of_pixels = 256;
number_of_pad_pixels = 128;
r_bleach = 32; % pixels

intensity_inside_bleach_region = 0.6;
intensity_outside_bleach_region = 0.9;

%% Particle parameters.
D = 350;
k_on = 1;
k_off = 1;
mobile_fraction = 0.5;

p_free = k_off / ( k_on + k_off );
p_bound = k_on / ( k_on + k_off );

tic
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
U0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 <= (upsampling_factor * r_bleach)^2 ) = intensity_inside_bleach_region;
U0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 > (upsampling_factor * r_bleach)^2 ) = intensity_outside_bleach_region;

% Distribute bound and free according to their marginal distribution,
% defined by p_free and p_bound.
B0 = p_bound * U0;
U0 = p_free * U0;

B0 = imresize(B0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);
U0 = imresize(U0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);
clear X Y

%% FFT of initial conditions.
F_U0 = fft2(U0);
F_B0 = fft2(B0);

%% FFT space time evolution of PDE system.
F_image_data_post_bleach_u = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
F_image_data_post_bleach_b = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
image_data_post_bleach_u = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
image_data_post_bleach_b = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);

[XSI1, XSI2] = meshgrid(-(number_of_pixels + 2 * number_of_pad_pixels)/2:(number_of_pixels + 2 * number_of_pad_pixels)/2-1, ...
                        -(number_of_pixels + 2 * number_of_pad_pixels)/2:(number_of_pixels + 2 * number_of_pad_pixels)/2-1);

XSI1 = XSI1 * 2 * pi / (number_of_pixels + 2 * number_of_pad_pixels);
XSI2 = XSI2 * 2 * pi / (number_of_pixels + 2 * number_of_pad_pixels);
XSISQ = XSI1.^2 + XSI2.^2;
XSISQ = ifftshift(XSISQ);

% Diagonalization, excluding time t which DD will be multiplied with.
PP11 = -(k_on - k_off + (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*k_on);
PP12 = -(k_on - k_off - (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*k_on);

DD11 = -((k_on + k_off + (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ))./2;
DD22 = -((k_on + k_off - (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ))./2;

PPinv11 = -k_on./(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2);
PPinv12 = -(k_on - k_off - (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2));
PPinv21 = -PPinv11;%k_on./(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2);
PPinv22 = (k_on - k_off + (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2));

CONST11 = PP11 .* (PPinv11 .* F_U0 + PPinv12 .* F_B0);
CONST12 = PP12 .* (PPinv21 .* F_U0 + PPinv22 .* F_B0);
CONST21 = PPinv11 .* F_U0 + PPinv12 .* F_B0;
CONST22 = PPinv21 .* F_U0 + PPinv22 .* F_B0;

for t = 1:number_of_post_bleach_images
    T = t * delta_t;
    CONST1 = exp(DD11 * T);
    CONST2 = exp(DD22 * T);
    F_image_data_post_bleach_u(:, :, t) = CONST11 .* CONST1 + CONST12 .* CONST2;
    F_image_data_post_bleach_b(:, :, t) = CONST21 .* CONST1 + CONST22 .* CONST2;
end

image_data_post_bleach = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_post_bleach_images);
for t = 1:number_of_post_bleach_images
%     image_data_post_bleach_u(:, :, t) = abs(ifft2(F_image_data_post_bleach_u(:, :, t)));
%     image_data_post_bleach_b(:, :, t) = abs(ifft2(F_image_data_post_bleach_b(:, :, t)));
    image_data_post_bleach(:, :, t) = abs(ifft2(F_image_data_post_bleach_u(:, :, t) + F_image_data_post_bleach_b(:, :, t)));
end
image_data_post_bleach = image_data_post_bleach(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels, :);
% image_data_post_bleach_u = image_data_post_bleach_u(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels, :);
% image_data_post_bleach_b = image_data_post_bleach_b(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels, :);

% image_data_post_bleach = image_data_post_bleach_u + image_data_post_bleach_b;

% Take (im)mobile fraction into account.
image_data_post_bleach = mobile_fraction * (image_data_post_bleach) + ...
    (1 - mobile_fraction) * (U0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels) + ...
    B0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels));

toc
%% Plot.
figure, imagesc(reshape(image_data_post_bleach, [number_of_pixels, number_of_pixels * number_of_post_bleach_images]))
axis 'equal'
axis([0 number_of_post_bleach_images*number_of_pixels 0 number_of_pixels])
axis off

x_bleach = 128;
y_bleach = 128;

[X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
X = X - 0.5;
Y = Y - 0.5;
ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
ind = ind(:);
recovery_curve = zeros(1, number_of_post_bleach_images);
for current_image_post_bleach = 1:number_of_post_bleach_images
    slice = image_data_post_bleach(:, :, current_image_post_bleach);
    recovery_curve(current_image_post_bleach) = mean(slice(ind));
end
figure
plot(delta_t:delta_t:number_of_post_bleach_images*delta_t, recovery_curve)