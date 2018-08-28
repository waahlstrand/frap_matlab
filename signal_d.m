function [C_prebleach, C_postbleach] = signal_d(sys_param, exp_sim_param)

% Extract system parameters.
D = sys_param(1);
mobile_fraction = sys_param(2);
C0 = sys_param(3);
alpha = sys_param(4);
beta = sys_param(5);

% Create a high-resolution bleach mask which is then downsampled to more
% accurately represent edges of the bleach region.
tic
bleach_mask = create_bleach_mask(alpha, exp_sim_param);
toc
% Create imaging bleach mask.
imaging_bleach_mask = create_imaging_bleach_mask(beta, exp_sim_param);

% Fourier space grid squared magnitude.
XSISQ = create_fourier_grid(exp_sim_param);

% Prebleach.
C_mobile = mobile_fraction * C0 * ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
C_immobile = (1 - mobile_fraction) * C0 * ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);

C_prebleach_mobile = zeros(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_prebleach_frames);
C_prebleach_immobile = zeros(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_prebleach_frames);

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

end