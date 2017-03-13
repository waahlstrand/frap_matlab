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
k_on = 0.05; % 1/s
k_off = 0.05; % 1/s
mobile_fraction = 0.5;%0.8;

p_free_marginal = k_off / ( k_on + k_off );
p_bound_marginal = k_on / ( k_on + k_off );

%% Generate post bleach image data for pure diffusion.

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
B0 = p_bound_marginal * U0;
U0 = p_free_marginal * U0;

B0 = imresize(B0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);
U0 = imresize(U0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);

% figure, imagesc(U0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels)), axis 'equal'

clear X Y

% return

%% Time evolution of PDE system.

image_data_post_bleach_u = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);
image_data_post_bleach_b = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);

U = U0;
B = B0;

% number_of_time_points_fine_per_coarse = 350;
number_of_time_points_fine_per_coarse = 700;
dt = delta_t / number_of_time_points_fine_per_coarse;

laplacian_operator = [ 0 1 0 ; 1 -4 1 ; 0 1 0];
tic
for current_post_bleach_image = 1:number_of_post_bleach_images
    disp(current_post_bleach_image)
    
    for current_time_fine = 1:number_of_time_points_fine_per_coarse
        
        U(1, :) = p_free_marginal * intensity_outside_bleach_region;
        U(end, :) = p_free_marginal * intensity_outside_bleach_region;
        U(:, 1) = p_free_marginal * intensity_outside_bleach_region;
        U(:, end) = p_free_marginal * intensity_outside_bleach_region;
        
        diffusion_term = D * conv2(U, laplacian_operator, 'same');
        
        U_new = ( diffusion_term - k_on * U + k_off * B ) * dt + U;
        B = ( k_on * U - k_off * B ) * dt + B;
        U = U_new;
    end
    
    image_data_post_bleach_u(:, :, current_post_bleach_image) = U(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);
    image_data_post_bleach_b(:, :, current_post_bleach_image) = B(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);
end

% Take (im)mobile fraction into account.
image_data_post_bleach = mobile_fraction * (image_data_post_bleach_u + image_data_post_bleach_b) + ...
    (1 - mobile_fraction) * (U0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels) + ...
    B0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels));
toc

%% Plot.
result_pde = [];
for current_image_post_bleach = 1:number_of_post_bleach_images
    result_pde = [result_pde, image_data_post_bleach(:, :, current_image_post_bleach)];
end
figure
imagesc(result_pde)
axis 'equal'
axis([0 number_of_post_bleach_images*number_of_pixels 0 number_of_pixels])
axis off
% 
% [X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
% X = X - 0.5;
% Y = Y - 0.5;
% xc = number_of_pixels / 2;
% yc = number_of_pixels / 2;
% ind = find( (X - xc).^2 + (Y - yc).^2 <= r_bleach_region^2 );
% ind = ind(:);
% recovery_curve = zeros(1, number_of_post_bleach_images);
% for current_image_post_bleach = 1:number_of_post_bleach_images
%     slice = image_data_post_bleach(:, :, current_image_post_bleach);
%     recovery_curve(current_image_post_bleach) = mean(slice(ind));
% end
% figure
% plot(delta_t:delta_t:number_of_post_bleach_images*delta_t, recovery_curve)

