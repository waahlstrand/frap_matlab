function [sys_param_hat, ss] = fit_d(C_prebleach, C_postbleach, sys_param, exp_sim_param, fit_param)

% Optimization options.
options = optimoptions(@lsqnonlin);
options.Algorithm = 'trust-region-reflective';
options.Display = 'iter';
options.FunctionTolerance = 1e-7;
options.OptimalityTolerance = 1e-7;
options.StepTolerance = 1e-7;
options.UseParallel = fit_param.use_parallel;

% Residual function handle.
fun = @(param)residual_d(C_prebleach, C_postbleach, param, exp_sim_param, fit_param);

% One fit with user-provided parameter guess or several fits using random 
% parameter guesses.
if ~isempty(fit_param.guess)
    [sys_param_hat, ss] = lsqnonlin(fun, fit_param.guess, fit_param.lower_bound, fit_param.upper_bound, options);
else
    sys_param_hat = zeros(size(fit_param.lower_bound));
    ss = inf;
    for current_fit = 1:number_of_fits
        sys_param_guess = fit_param.lower_bound + (fit_param.upper_bound - fit_param.lower_bound) .* rand(size(fit_param.lower_bound));
        [sys_param_hat_, ss_] = lsqnonlin(fun, sys_param_guess, fit_param.lower_bound, fit_param.upper_bound, options);
        if ss_ < ss
            sys_param_hat = sys_param_hat_;
            ss = ss_;
        end
    end
end

end

