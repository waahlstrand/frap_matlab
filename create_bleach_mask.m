function bleach_mask = create_bleach_mask(alpha, exp_sim_param)

xx = 1:exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
[X, Y] = meshgrid(xx, xx);
X = X - 0.5;
Y = Y - 0.5;

bleach_mask = ones(size(X));
if exp_sim_param.bleach_region.shape == "circular"
    idx_bleach =    ( X - exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.x) ).^2 + ...
                    ( Y - exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.x) ).^2 <= ...
                    ( exp_sim_param.bleach_region.upsampling_factor * exp_sim_param.bleach_region.r )^2;
elseif exp_sim_param.bleach_region.shape == "rectangular"
    idx_bleach =    X >= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.x - 0.5 * exp_sim_param.bleach_region.lx) & ...
                    X <= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.x + 0.5 * exp_sim_param.bleach_region.lx) & ...
                    Y >= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.y - 0.5 * exp_sim_param.bleach_region.ly) & ...
                    Y <= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.y + 0.5 * exp_sim_param.bleach_region.ly);
else
    error('Unrecognized bleach region shape.');
end
bleach_mask(idx_bleach) = alpha;
bleach_mask = imresize(bleach_mask, [exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels], 'method', 'bilinear');

end