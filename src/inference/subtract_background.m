function image_data_post_bleach = subtract_background(image_data_pre_bleach, image_data_post_bleach)

% Compute pixel-wise time average of pre-bleach intensity.
mean_image_pre_bleach = mean(image_data_pre_bleach, 3);

% Compute complete average of pre-bleach intensity.
mean_intensity_pre_bleach = mean(image_data_pre_bleach(:));

% Remove background from post-bleach data and add mean background to
% retain similar scaling.
number_of_post_bleach_images = size(image_data_post_bleach, 3);
image_data_post_bleach = image_data_post_bleach - repmat(mean_image_pre_bleach, [1, 1, number_of_post_bleach_images]);
image_data_post_bleach = image_data_post_bleach + mean_intensity_pre_bleach;

end

