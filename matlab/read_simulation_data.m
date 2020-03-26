function [C_prebleach_sim, C_postbleach_sim, sys_param, exp_sim_param, t_exec] = read_simulation_data(file_path)

%% Read file.
file_id = fopen(file_path);

D_SI = fread(file_id, 1, 'float64');
k_on = fread(file_id, 1, 'float64');
k_off = fread(file_id, 1, 'float64');
mobile_fraction = fread(file_id, 1, 'float64');
alpha = fread(file_id, 1, 'float64');
beta = fread(file_id, 1, 'float64');
gamma = fread(file_id, 1, 'float64');
bleach_region_shape = fread(file_id, 1, 'int64');
r_bleach = fread(file_id, 1, 'float64');
lx_bleach = fread(file_id, 1, 'float64');
ly_bleach = fread(file_id, 1, 'float64');
number_of_pixels = fread(file_id, 1, 'int64');
number_of_prebleach_frames = fread(file_id, 1, 'int64');
number_of_bleach_frames = fread(file_id, 1, 'int64');
number_of_postbleach_frames = fread(file_id, 1, 'int64');
number_of_pad_pixels = fread(file_id, 1, 'int64');
delta_t = fread(file_id, 1, 'float64');
pixel_size = fread(file_id, 1, 'float64');
number_of_particles = fread(file_id, 1, 'int64');
t_exec = fread(file_id, 1, 'float64');
data = reshape(fread(file_id, 'int64'), [number_of_pixels, number_of_pixels, number_of_prebleach_frames + number_of_postbleach_frames]);

C_prebleach_sim = double(data(:, :, 1:number_of_prebleach_frames));
C_postbleach_sim = double(data(:, :, number_of_prebleach_frames+1:end));

fclose(file_id);

clear file_id 

D = D_SI / pixel_size^2;

%% Normalize count data.
C0 = 1;
norm_factor = number_of_particles / (number_of_pixels + 2 * number_of_pad_pixels)^2; % Average count in pixels.
C_prebleach_sim = C_prebleach_sim / norm_factor;
C_postbleach_sim = C_postbleach_sim / norm_factor;

%% Store parameters.
a = 0;
b = 0;
sys_param = [D, k_on, k_off, mobile_fraction, C0, alpha, beta, gamma, a, b];

exp_sim_param = struct();

exp_sim_param.pixel_size = pixel_size; % m
exp_sim_param.number_of_pixels = number_of_pixels;

exp_sim_param.number_of_prebleach_frames = number_of_prebleach_frames;
exp_sim_param.number_of_bleach_frames = number_of_bleach_frames;
exp_sim_param.number_of_postbleach_frames = number_of_postbleach_frames;
exp_sim_param.delta_t = delta_t; % s

exp_sim_param.number_of_pad_pixels = number_of_pad_pixels;

if bleach_region_shape == 0
    exp_sim_param.bleach_region.shape = "circle";
else 
    exp_sim_param.bleach_region.shape = "rectangle";
end
exp_sim_param.bleach_region.x = number_of_pixels / 2; % pixels
exp_sim_param.bleach_region.y = number_of_pixels / 2; % pixels 
exp_sim_param.bleach_region.r = r_bleach; % pixels
exp_sim_param.bleach_region.lx = lx_bleach; % pixels
exp_sim_param.bleach_region.ly = ly_bleach; % pixels

end