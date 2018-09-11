function [rc_prebleach, rc_postbleach] = recovery_curve(C_prebleach, C_postbleach, exp_sim_param)

% Create bleach region indicator.
bleach_region_indicator = 1 - create_bleach_mask(0, 0, exp_sim_param);
bleach_region_indicator = bleach_region_indicator(exp_sim_param.number_of_pad_pixels+1:end-exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pad_pixels+1:end-exp_sim_param.number_of_pad_pixels); 
bleach_region_indicator = bleach_region_indicator / sum(bleach_region_indicator(:));

% Extract prebleach intensity.
rc_prebleach = zeros(exp_sim_param.number_of_prebleach_frames, 1);
for current_frame = 1:exp_sim_param.number_of_prebleach_frames
    rc_prebleach(current_frame) = sum(sum( C_prebleach(:, :, current_frame) .* bleach_region_indicator));
end

% Extract postbleach intensity.
rc_postbleach = zeros(exp_sim_param.number_of_postbleach_frames, 1);
for current_frame = 1:exp_sim_param.number_of_postbleach_frames
    rc_postbleach(current_frame) = sum(sum( C_postbleach(:, :, current_frame) .* bleach_region_indicator));
end

end

