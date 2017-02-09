function image_data_post_bleach = signal_diffusion_and_binding(D, k_on, k_off, mobile_fraction, xc, yc, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region, delta_t, number_of_pixels, number_of_images_post_bleach)

disp([D, k_on, k_off, mobile_fraction, xc, yc, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region])

number_of_pad_pixels = 128;

p_free_marginal = k_off / ( k_on + k_off );
p_bound_marginal = k_on / ( k_on + k_off );

% Initial condition. Create a high resolution initial condition which is 
% then downsampled to avoid too sharp edges. Distribute bound and free 
% according to their marginal distribution, defined by p_free and p_bound.

upsampling_factor = 3;

[X, Y] = meshgrid(1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels), 1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels));
X = X - 0.5;
Y = Y - 0.5;
xc = number_of_pad_pixels + xc;
yc = number_of_pad_pixels + yc;

U0 = zeros(size(X));
U0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 <= (upsampling_factor * r_bleach_region)^2 ) = intensity_inside_bleach_region;
U0( (X - upsampling_factor * xc).^2 + (Y - upsampling_factor * yc).^2 > (upsampling_factor * r_bleach_region)^2 ) = intensity_outside_bleach_region;

U0 = imresize(U0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);

B0 = p_bound_marginal * U0;
U0 = p_free_marginal * U0;

% Time evolution of PDE system.

image_data_post_bleach_u = zeros(number_of_pixels, number_of_pixels, number_of_images_post_bleach);
image_data_post_bleach_b = zeros(number_of_pixels, number_of_pixels, number_of_images_post_bleach);

U = U0;
B = B0;

number_of_time_points_fine_per_coarse = 600;
dt = delta_t / number_of_time_points_fine_per_coarse;

laplacian_operator = [ 0 1 0 ; 1 -4 1 ; 0 1 0];

n = number_of_pixels + 2 * number_of_pad_pixels; 
x_edges = [1:n n*ones(1, n-2) n:-1:1 n*ones(1, n-2)];
y_edges = [ones(1, n-1) 1:n  n*ones(1, n-2) n:-1:2];
ind_edges = sub2ind([n n], x_edges, y_edges);

for current_post_bleach_image = 1:number_of_images_post_bleach
    
    for current_time_fine = 1:number_of_time_points_fine_per_coarse
        
        U(ind_edges) = p_free_marginal * intensity_outside_bleach_region;
      
        del_U = conv2(U, laplacian_operator, 'same');
        
        U_new = ( D * del_U - k_on * U + k_off * B ) * dt + U;
        B = ( k_on * U - k_off * B ) * dt + B;
        U = U_new;
    end
    
    image_data_post_bleach_u(:, :, current_post_bleach_image) = U(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);
    image_data_post_bleach_b(:, :, current_post_bleach_image) = B(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);
end

% Take (im)mobile fraction into account and add the free and bound
% contribution to the fluorescence. Even though U0 and B0 are 2-dim
% matrices, they are actually added to each slice in dim 3 of the 3-D
% matrix by writing like this.
image_data_post_bleach = mobile_fraction * (image_data_post_bleach_u + image_data_post_bleach_b) + ...
    (1 - mobile_fraction) * (U0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels) + ...
    B0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels));

end

