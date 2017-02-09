%% Initialization.
% clear
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

C0 = zeros(size(X));
C0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 <= (upsampling_factor * r_bleach_region)^2 ) = intensity_inside_bleach_region;
C0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 > (upsampling_factor * r_bleach_region)^2 ) = intensity_outside_bleach_region;

C0 = imresize(C0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);

clear X Y

%% Time evolution of PDE system.

C = C0;

% number_of_time_points_fine_per_coarse = 350;
number_of_time_points_fine_per_coarse = 100;
dt = delta_t / number_of_time_points_fine_per_coarse;

[XP, YP] = meshgrid(-3:3, -3:3);
free_propagator = normpdf(XP, 0, sqrt(2*D*dt)) .* normpdf(YP, 0, sqrt(2*D*dt));
free_propagator = free_propagator / sum(free_propagator (:));
bound_propagator = zeros(size(XP));
bound_propagator(ceil(size(XP, 1) / 2), ceil(size(XP, 2) / 2)) = 1;

propagator = p_free_marginal * free_propagator + p_bound_marginal * bound_propagator

image_data_post_bleach = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);

tic
for current_post_bleach_image = 1:number_of_post_bleach_images
    disp(current_post_bleach_image)
    
    for current_time_fine = 1:number_of_time_points_fine_per_coarse
        
        C(1, :) = p_free_marginal * intensity_outside_bleach_region;
        C(end, :) = p_free_marginal * intensity_outside_bleach_region;
        C(:, 1) = p_free_marginal * intensity_outside_bleach_region;
        C(:, end) = p_free_marginal * intensity_outside_bleach_region;
        
        C = conv2(C, propagator, 'same');
    end
    
    image_data_post_bleach(:, :, current_post_bleach_image) = C(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);
end

% Take (im)mobile fraction into account.
image_data_post_bleach = mobile_fraction * image_data_post_bleach + ...
    (1 - mobile_fraction) * C0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);
toc

%% Plot.
result_pde_special = [];
for current_image_post_bleach = 1:number_of_post_bleach_images
    result_pde_special = [result_pde_special, image_data_post_bleach(:, :, current_image_post_bleach)];
end
figure
imagesc(result_pde_special)
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

