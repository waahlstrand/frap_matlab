function F = residual_diffusion_and_binding(D, k_on, k_off, mobile_fraction, x_bleach, y_bleach, r_bleach_region, intensity_inside_bleach_region, intensity_outside_bleach_region, delta_t, image_data_post_bleach)

[number_of_pixels, ~, number_of_post_bleach_images] = size(image_data_post_bleach);

image_data_post_bleach_model = signal_diffusion_and_binding(D, ...
                                                            k_on, ...
                                                            k_off, ...
                                                            mobile_fraction, ...
                                                            x_bleach, ...
                                                            y_bleach, ...
                                                            r_bleach_region, ...
                                                            intensity_inside_bleach_region, ...
                                                            intensity_outside_bleach_region, ...
                                                            delta_t, ...
                                                            number_of_time_points_fine_per_coarse, ...
                                                            number_of_pixels, ...
                                                            number_of_post_bleach_images, ...
                                                            number_of_pad_pixels)

F = image_data_post_bleach_model(:) - image_data_post_bleach(:);

end

