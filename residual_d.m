function R = residual_d(C_prebleach, C_postbleach, sys_param, exp_sim_param, fit_param)

[C_prebleach_model, C_postbleach_model] = signal_d(sys_param, exp_sim_param);
               
if fit_param.mode == "recovery-curve"
    [rc_prebleach_model, rc_postbleach_model] = recovery_curve(C_prebleach_model, C_postbleach_model, exp_sim_param);
    rc_model = [rc_prebleach_model ; rc_postbleach_model];
    [rc_prebleach, rc_postbleach] = recovery_curve(C_prebleach, C_postbleach, exp_sim_param);
    rc = [rc_prebleach ; rc_postbleach];
    
    R = rc_model(:) - rc(:);
elseif fit_param.mode == "pixel"
    C_model = cat(3, C_prebleach_model, C_postbleach_model);
    C = cat(3, C_prebleach, C_postbleach);
    R = C_model(:) - C(:);
else
    error('Invalid estimation mode.')
end

end