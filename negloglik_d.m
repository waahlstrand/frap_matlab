function nll = negloglik_d(C_prebleach, C_postbleach, sys_param, exp_sim_param, fit_param)

[C_prebleach_model, C_postbleach_model] = signal_d(sys_param, exp_sim_param);

a = sys_param(end-1);
b = sys_param(end);
               
if fit_param.mode == "recovery-curve"
    w = 1 - create_bleach_mask(0, 0, exp_sim_param);
    w = w(exp_sim_param.number_of_pad_pixels+1:end-exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pad_pixels+1:end-exp_sim_param.number_of_pad_pixels); 
    w = w / sum(w(:));
    
    sigma2_prebleach = zeros(exp_sim_param.number_of_prebleach_frames, 1);
    for current_frame = 1:exp_sim_param.number_of_prebleach_frames
        sigma2_prebleach(current_frame) = sum(sum( w.^2 .* (a + b * C_prebleach_model(:, :, current_frame)) ));
    end
    sigma2_postbleach = zeros(exp_sim_param.number_of_postbleach_frames, 1);
    for current_frame = 1:exp_sim_param.number_of_postbleach_frames
        sigma2_postbleach(current_frame) = sum(sum( w.^2 .* (a + b * C_postbleach_model(:, :, current_frame)) ));
    end
    sigma2 = [sigma2_prebleach ; sigma2_postbleach];    

    [rc_prebleach_model, rc_postbleach_model] = recovery_curve(C_prebleach_model, C_postbleach_model, exp_sim_param);
    rc_model = [rc_prebleach_model ; rc_postbleach_model];
    [rc_prebleach, rc_postbleach] = recovery_curve(C_prebleach, C_postbleach, exp_sim_param);
    rc = [rc_prebleach ; rc_postbleach];
    
    nll = sum( 0.5 * log(sigma2(:)) + 0.5 ./ sigma2(:) .*  (rc_model(:) - rc(:)).^2 );
    if a == 0 && b == 0
        nll = Inf;
    end    
elseif fit_param.mode == "pixel"
    C_model = cat(3, C_prebleach_model, C_postbleach_model);
    C = cat(3, C_prebleach, C_postbleach);
    SIGMA2 = a + b * C_model;
    
    nll = sum( 0.5 * log(SIGMA2(:)) + 0.5 ./ SIGMA2(:) .*  (C_model(:) - C(:)).^2 );
    if a == 0 && b == 0
        nll = Inf;
    end
else
    error('Invalid estimation mode.')
end

end