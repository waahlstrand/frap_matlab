function [F, J] = residual_recoverycurve_diffusion_and_binding_partial( mobile_fraction, ...
                                                                        intensity_inside_bleach_region, ...
                                                                        intensity_outside_bleach_region, ...
                                                                        image_data_post_bleach, ...
                                                                        image_data_post_bleach_model_unscaled, ...
                                                                        initial_condition_model_unscaled)

% disp([mobile_fraction intensity_inside_bleach_region intensity_outside_bleach_region])
                                                        
initial_condition_model_unscaled = repmat(initial_condition_model_unscaled, [1, 1, size(image_data_post_bleach_model_unscaled, 3)]);

image_data_post_bleach_model_unscaled = image_data_post_bleach_model_unscaled(:);
initial_condition_model_unscaled = initial_condition_model_unscaled(:);

image_data_post_bleach_model = 2 * (intensity_outside_bleach_region - intensity_inside_bleach_region) * ( mobile_fraction * image_data_post_bleach_model_unscaled + (1 - mobile_fraction) * initial_condition_model_unscaled ) + 2 * intensity_inside_bleach_region - intensity_outside_bleach_region;

F = image_data_post_bleach_model - image_data_post_bleach(:);

if nargout > 1
    
    J = zeros(numel(image_data_post_bleach_model), 3);
    
    % Derivative w.r.t. 'mobile_fraction'.
    J(:, 1) = 2 * (intensity_outside_bleach_region - intensity_inside_bleach_region) * ( image_data_post_bleach_model_unscaled - initial_condition_model_unscaled );
    
    % Derivative w.r.t. 'intensity_inside_bleach_region'.
    J(:, 2) = - 2 * ( mobile_fraction * image_data_post_bleach_model_unscaled + (1 - mobile_fraction) * initial_condition_model_unscaled ) + 2;

    % Derivative w.r.t. 'intensity_outside_bleach_region'.
    J(:, 3) = 2 * ( mobile_fraction * image_data_post_bleach_model_unscaled + (1 - mobile_fraction) * initial_condition_model_unscaled ) - 1;
end