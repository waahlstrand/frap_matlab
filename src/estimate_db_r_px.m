function [param_hat, ss] = estimate_db_r_px(data, ...
                                            x_bleach, ...
                                            y_bleach, ...
                                            lx_bleach, ...
                                            ly_bleach, ...
                                            delta_t, ...
                                            number_of_pixels, ...
                                            number_of_images, ...
                                            number_of_pad_pixels, ...
                                            lb, ...
                                            ub, ...
                                            param_guess, ...
                                            number_of_fits)

% Optimization options.
options = optimoptions(@lsqnonlin);
options.Algorithm = 'trust-region-reflective';
options.Display = 'iter';
options.FunctionTolerance = 1e-7;
options.OptimalityTolerance = 1e-7;
options.StepTolerance = 1e-7;

% Residual function handle.
fun = @(param)residual_db_r_px( param(1), ...
                                param(2), ...
                                param(3), ...
                                param(4), ...
                                param(5), ...
                                param(6), ...
                                x_bleach, ...
                                y_bleach, ...
                                lx_bleach, ...
                                ly_bleach, ...
                                delta_t, ...
                                number_of_pixels, ...
                                number_of_images, ...
                                number_of_pad_pixels, ...
                                data);

% One fit with user-provided parameter guess or arbitrary many fits using
% random parameter guesses.
if ~isempty(param_guess)
    [param_hat, ss] = lsqnonlin(fun, param_guess, lb, ub, options);
else
    param_hat = zeros(1, 6);
    ss = inf;
    
    for current_fit = 1:number_of_fits
        param_guess = lb + (ub - lb) .* rand(size(lb));
        [param_hat_, ss_] = lsqnonlin(fun, param_guess, lb, ub, options);
        if ss_ < ss
            param_hat = param_hat_;
            ss = ss_;
        end
    end
end

end

