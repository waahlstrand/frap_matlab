clear
clc
close all

number_of_pixels = 128;
number_of_post_bleach_images = 10;
image_data_post_bleach_sim = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);
for i = 1:6
    file_id = fopen(['simulated_frap_data_' num2str(i) '.dat']);
    image_data_post_bleach_cell{i} = fread(file_id, [number_of_pixels*number_of_pixels*number_of_post_bleach_images, 1], 'int64');
    fclose(file_id);
    image_data_post_bleach_cell{i} = reshape(image_data_post_bleach_cell{i}, [number_of_pixels, number_of_pixels, number_of_post_bleach_images]);
    image_data_post_bleach_sim = image_data_post_bleach_sim + image_data_post_bleach_cell{i};
end

image_data_post_bleach_sim = image_data_post_bleach_sim / mean(mean(image_data_post_bleach_sim(:, 1:5, 1))) * 0.9;

result_sim = [];
for current_image_post_bleach = 1:number_of_post_bleach_images
    result_sim = [result_sim, image_data_post_bleach_sim(:, :, current_image_post_bleach)];
end
figure
imagesc(result_sim)
axis 'equal'
axis([0 number_of_post_bleach_images*number_of_pixels 0 number_of_pixels])
axis off

image_data_post_bleach = image_data_post_bleach_sim;
