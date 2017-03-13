function F = residual_recoverycurve_diffusion_and_binding_partial(  mobile_fraction, ...
                                                                    intensity_inside_bleach_region, ...
                                                                    intensity_outside_bleach_region, ...
                                                                    x_bleach, ...
                                                                    y_bleach, ...
                                                                    r_bleach, ...
                                                                    image_data_post_bleach, ...
                                                                    image_data_post_bleach_model_unscaled, ...
                                                                    initial_condition_model_unscaled)
                     
[number_of_pixels, ~, number_of_post_bleach_images] = size(image_data_post_bleach);

initial_condition_model_unscaled = repmat(initial_condition_model_unscaled, [1, 1, size(image_data_post_bleach_model_unscaled, 3)]);

% image_data_post_bleach_model_unscaled = image_data_post_bleach_model_unscaled(:);
% initial_condition_model_unscaled = initial_condition_model_unscaled(:);

image_data_post_bleach_model = 2 * (intensity_outside_bleach_region - intensity_inside_bleach_region) * ( mobile_fraction * image_data_post_bleach_model_unscaled + (1 - mobile_fraction) * initial_condition_model_unscaled ) + 2 * intensity_inside_bleach_region - intensity_outside_bleach_region;

[X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
X = X - 0.5;
Y = Y - 0.5;
ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
ind = ind(:);

% disp(number_of_post_bleach_images)
% size(image_data_post_bleach_model)
recovery_curve = zeros(1, number_of_post_bleach_images);
for current_image_post_bleach = 1:number_of_post_bleach_images
    slice = image_data_post_bleach(:, :, current_image_post_bleach);
    recovery_curve(current_image_post_bleach) = mean(slice(ind));
end

recovery_curve_model = zeros(1, number_of_post_bleach_images);
for current_image_post_bleach = 1:number_of_post_bleach_images
    slice = image_data_post_bleach_model(:, :, current_image_post_bleach);
    recovery_curve_model(current_image_post_bleach) = mean(slice(ind));
end

F = recovery_curve_model(:) - recovery_curve(:);

end