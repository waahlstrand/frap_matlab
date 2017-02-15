clear
clc
close all

number_of_pixels = 128;
number_of_post_bleach_images = 10;
image_data_post_bleach_sim = zeros(number_of_pixels, number_of_pixels, number_of_post_bleach_images);

file_id = fopen('simulated_frap_data_1993670910087.dat');
image_data_post_bleach = fread(file_id, [number_of_pixels*number_of_pixels*number_of_post_bleach_images, 1], 'int64');
fclose(file_id);
image_data_post_bleach = reshape(image_data_post_bleach, [number_of_pixels, number_of_pixels, number_of_post_bleach_images]);

result_sim = [];
for current_image_post_bleach = 1:number_of_post_bleach_images
    result_sim = [result_sim, image_data_post_bleach(:, :, current_image_post_bleach)];
end
figure
imagesc(result_sim)
axis 'equal'
axis([0 number_of_post_bleach_images*number_of_pixels 0 number_of_pixels])
axis off