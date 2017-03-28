function [image_data_post_bleach, initial_condition] = signal_diffusion(D, ...
                                                                        mobile_fraction, ...
                                                                        x_bleach, ...
                                                                        y_bleach, ...
                                                                        r_bleach, ...
                                                                        intensity_inside_bleach_region, ...
                                                                        intensity_outside_bleach_region, ...
                                                                        delta_t, ...
                                                                        number_of_time_points_fine_per_coarse, ...
                                                                        number_of_pixels, ...
                                                                        number_of_post_bleach_images, ...
                                                                        number_of_pad_pixels)

% If time step not given, compute the 'adaptive' time step.
if isempty(number_of_time_points_fine_per_coarse)
    number_of_time_points_fine_per_coarse= ceil( D * delta_t / 0.225 );
end

% Initial condition. Create a high resolution initial condition which is 
% then downsampled to avoid too sharp edges. Distribute bound and free 
% according to their marginal distribution, defined by p_free and p_bound.

upsampling_factor = 3;

[X, Y] = meshgrid(1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels), 1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels));
X = X - 0.5;
Y = Y - 0.5;
x_bleach = number_of_pad_pixels + x_bleach;
y_bleach = number_of_pad_pixels + y_bleach;

C0 = zeros(size(X));
C0( (X - upsampling_factor * x_bleach).^2 + (Y - upsampling_factor * y_bleach).^2 <= (upsampling_factor * r_bleach)^2 ) = intensity_inside_bleach_region;
C0( (X - upsampling_factor * x_bleach).^2 + (Y - upsampling_factor * y_bleach).^2 > (upsampling_factor * r_bleach)^2 ) = intensity_outside_bleach_region;

C0 = imresize(C0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);

% Time evolution of PDE.

image_data_post_bleach = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);

C = C0;

dt = delta_t / number_of_time_points_fine_per_coarse;

laplacian_operator = [ 0 1 0 ; 1 -4 1 ; 0 1 0];

n = number_of_pixels + 2 * number_of_pad_pixels; 
x_edges = [1:n n*ones(1, n-2) n:-1:1 1*ones(1, n-2)];
y_edges = [ones(1, n-1) 1:n  n*ones(1, n-2) n:-1:2];
ind_edges = sub2ind([n n], x_edges, y_edges);

for current_post_bleach_image = 1:number_of_post_bleach_images
    disp(current_post_bleach_image)
    
    for current_time_fine = 1:number_of_time_points_fine_per_coarse
        
        C(ind_edges) = intensity_outside_bleach_region;
      
        del_C = conv2(C, laplacian_operator, 'same');
        
        C = D * del_C * dt + C;
    end
    
    image_data_post_bleach(:, :, current_post_bleach_image) = C(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);
end

% Take (im)mobile fraction into account and add the free and bound
% contribution to the fluorescence. Even though U0 and B0 are 2-dim
% matrices, they are actually added to each slice in dim 3 of the 3-D
% matrix by writing like this.
image_data_post_bleach = mobile_fraction * image_data_post_bleach + ...
    (1 - mobile_fraction) * C0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);

initial_condition = C0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);

end

