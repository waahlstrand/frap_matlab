% function data = signal_d(sys_param, exp_sim_param)

%% Initialization.

% clear
clc
close all hidden

%% Experimental and simulation parameters.

exp_sim_param = struct();

exp_sim_param.pixel_size = 7.5e-07; % m
exp_sim_param.number_of_pixels = 256;

exp_sim_param.number_of_prebleach_frames = 5;
exp_sim_param.number_of_bleach_frames = 2;
exp_sim_param.number_of_postbleach_frames = 5;
exp_sim_param.delta_t = 0.2; % s

exp_sim_param.number_of_pad_pixels = 128;

exp_sim_param.bleach_region.shape = "circular";
exp_sim_param.bleach_region.x = 128; % pixels
exp_sim_param.bleach_region.y = 128; % pixels 
exp_sim_param.bleach_region.r = 15e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.lx = 32%20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.ly = 32%20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.upsampling_factor  = 3;

%% System parameters.

D_SI = 5e-11; % m^2/s
D = D_SI / exp_sim_param.pixel_size^2; % pixels^2 / s
mobile_fraction = 0.5; % dimensionless
C0 = 1.0; % a.u. original concentration
alpha = 0.6; % a.u.  bleach factor
beta = 1.00; % a.u. imaging bleach factor

sys_param = [D, mobile_fraction, C0, alpha, beta];

%% Simulaton.

tic

% Create a high-resolution bleach mask which is then downsampled to more
% accurately represent edges of the bleach region.
bleach_mask = create_bleach_mask(alpha, exp_sim_param);

% Create imaging bleach mask.
imaging_bleach_mask = create_imaging_bleach_mask(beta, exp_sim_param);

% Fourier space grid squared magnitude.
XSISQ = create_fourier_grid(exp_sim_param);

% Prebleach.
C_mobile = mobile_fraction * C0 * ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
C_immobile = (1 - mobile_fraction) * C0 * ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);

C_prebleach_mobile = zeros(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_prebleach_frames);
C_prebleach_immobile = zeros(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_prebleach_frames);

% F_C_mobile = zeros(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
% F_C_immobile = zeros(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);

for current_frame = 1:exp_sim_param.number_of_prebleach_frames
    F_C_mobile = fft2(C_mobile);
    F_C_mobile = exp( - D * XSISQ * exp_sim_param.delta_t ) .* F_C_mobile;
    C_mobile = abs(ifft2(F_C_mobile));
    C_mobile = C_mobile .* imaging_bleach_mask;
    
    C_immobile = C_immobile .* imaging_bleach_mask;
    
    C_prebleach_mobile(:, :, current_frame) = C_mobile;
    C_prebleach_immobile(:, :, current_frame) = C_immobile;
end

C_prebleach = C_prebleach_mobile(exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, :) + ...
                C_prebleach_immobile(exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, :);

% Bleach.
for current_frame = 1:exp_sim_param.number_of_bleach_frames
    F_C_mobile = fft2(C_mobile);
    F_C_mobile = exp( - D * XSISQ * exp_sim_param.delta_t ) .* F_C_mobile;
    C_mobile = abs(ifft2(F_C_mobile));
    C_mobile = C_mobile .* bleach_mask .* imaging_bleach_mask;
    
    C_immobile = C_immobile .* bleach_mask .* imaging_bleach_mask;
end

% Postbleach.
C_postbleach_mobile = zeros(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_postbleach_frames);
C_postbleach_immobile = repmat(C_immobile, [1, 1, exp_sim_param.number_of_postbleach_frames]);

for current_frame = 1:exp_sim_param.number_of_postbleach_frames
    F_C_mobile = fft2(C_mobile);
    F_C_mobile = exp( - D * XSISQ * exp_sim_param.delta_t ) .* F_C_mobile;
    C_mobile = abs(ifft2(F_C_mobile));
    C_mobile = C_mobile .* imaging_bleach_mask;
    
    C_immobile = C_immobile .* imaging_bleach_mask;
    
    C_postbleach_mobile(:, :, current_frame) = C_mobile;
    C_postbleach_immobile(:, :, current_frame) = C_immobile;
    
end

C_postbleach = C_postbleach_mobile(exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, :) + ...
                C_postbleach_immobile(exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, :);

 toc           

% for current_frame = 1:exp_sim_param.number_of_prebleach_frames
%     figure, imagesc([C_prebleach(:, :, current_frame), C_prebleach_sim(:, :, current_frame), C_prebleach(:, :, current_frame) - C_prebleach_sim(:, :, current_frame)]), axis 'equal'
% end
% for current_frame = 1:exp_sim_param.number_of_postbleach_frames
%     figure, imagesc([C_postbleach(:, :, current_frame), C_postbleach_sim(:, :, current_frame), C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)]), axis 'equal'
% end
for current_frame = 1:exp_sim_param.number_of_prebleach_frames
    figure, imagesc(C_prebleach(:, :, current_frame) - C_prebleach_sim(:, :, current_frame)), axis 'equal'
end
for current_frame = 1:exp_sim_param.number_of_postbleach_frames
    figure, imagesc(C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)), axis 'equal'
end
max(abs(C_prebleach(:)-C_prebleach_sim(:)))
max(abs(C_postbleach(:)-C_postbleach_sim(:)))

% figure, hold on, plot(C_postbleach(:, 128), 'k.-'), plot(C_postbleach_sim(:, 128), 'r.-')
% figure, plot(C_postbleach(128, :) - C_postbleach_sim(128, :))
