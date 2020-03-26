function [sys_param_hat, ss] = fit_db(C_prebleach, C_postbleach, exp_sim_param, fit_param)

if fit_param.upper_bound(end) == 0 % b = 0 => Ordinary least squares
    % Optimization options.
    options = optimoptions(@lsqnonlin);
    options.Algorithm = 'trust-region-reflective';
    options.Display = 'iter';
%     options.FunctionTolerance = 1e-6;
%     options.OptimalityTolerance = 1e-6;
%     options.StepTolerance = 1e-6;
    options.UseParallel = fit_param.use_parallel;
    
    % Remove a and b because b = 0 and we are estimating a from the
    % residual post-fit. That means the bounds for a are not used.
    fit_param.guess = fit_param.guess(1:end - 2);
    fit_param.lower_bound = fit_param.lower_bound(1:end - 2);
    fit_param.upper_bound = fit_param.upper_bound(1:end - 2);

    % Residual function handle.
    fun = @(param)residual_db(C_prebleach, C_postbleach, param, exp_sim_param, fit_param);

    % One fit with user-provided parameter guess or several fits using random 
    % parameter guesses.
    if ~isempty(fit_param.guess)
        [sys_param_hat, ss] = lsqnonlin(fun, fit_param.guess, fit_param.lower_bound, fit_param.upper_bound, options);
    else
        sys_param_hat = zeros(size(fit_param.lower_bound));
        ss = inf;
        for current_fit = 1:fit_param.number_of_fits
            sys_param_guess = fit_param.lower_bound + (fit_param.upper_bound - fit_param.lower_bound) .* rand(size(fit_param.lower_bound));
            [sys_param_hat_, ss_] = lsqnonlin(fun, sys_param_guess, fit_param.lower_bound, fit_param.upper_bound, options);
            if ss_ < ss
                sys_param_hat = sys_param_hat_;
                ss = ss_;
            end
        end
    end
    
    if fit_param.mode == "recovery-curve"
        w = 1 - create_bleach_mask(0, 0, exp_sim_param);
        w = w(exp_sim_param.number_of_pad_pixels+1:end-exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pad_pixels+1:end-exp_sim_param.number_of_pad_pixels); 
        w = w / sum(w(:));
        a_hat = mean( fun(sys_param_hat).^2 ) / sum( w(:).^2 ); % a_hat is estimated variance of individual pixels not the recovery curve data points.
    else
        a_hat = mean( fun(sys_param_hat).^2 );
    end
    sys_param_hat = [sys_param_hat, a_hat, 0];
else % b >= 0 => Maximum likelihood
    % Optimization options.
    options = optimoptions(@fmincon);
    options.Algorithm = 'sqp';
    options.Display = 'iter';
%     options.FunctionTolerance = 1e-6;
%     options.OptimalityTolerance = 1e-6;
%     options.StepTolerance = 1e-6;
    options.UseParallel = fit_param.use_parallel;

    % Negative loglikelihood function handle.
    fun = @(param)negloglik_db(C_prebleach, C_postbleach, param, exp_sim_param, fit_param);
    
    % One fit with user-provided parameter guess or several fits using random 
    % parameter guesses.
    if ~isempty(fit_param.guess)
        [sys_param_hat, ss] = fmincon(fun, fit_param.guess, [], [], [], [], fit_param.lower_bound, fit_param.upper_bound, [], options);
    else
        sys_param_hat = zeros(size(fit_param.lower_bound));
        ss = inf;
        for current_fit = 1:fit_param.number_of_fits
            sys_param_guess = fit_param.lower_bound + (fit_param.upper_bound - fit_param.lower_bound) .* rand(size(fit_param.lower_bound));
            [sys_param_hat_, ss_] = fmincon(fun, sys_param_guess, [], [], [], [], fit_param.lower_bound, fit_param.upper_bound, [], options);
            if ss_ < ss
                sys_param_hat = sys_param_hat_;
                ss = ss_;
            end
        end
    end
end   

end

