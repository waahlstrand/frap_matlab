function bleach_mask = create_bleach_mask(alpha, exp_sim_param)

% xx = 1:exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
% [X, Y] = meshgrid(xx, xx);
% X = X - 0.5;
% Y = Y - 0.5;
% 
% bleach_mask = ones(size(X));
% if exp_sim_param.bleach_region.shape == "circular"
%     idx_bleach =    ( X - exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.x) ).^2 + ...
%                     ( Y - exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.y) ).^2 <= ...
%                     ( exp_sim_param.bleach_region.upsampling_factor * exp_sim_param.bleach_region.r )^2;
% elseif exp_sim_param.bleach_region.shape == "rectangular"
%     idx_bleach =    X >= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.x - 0.5 * exp_sim_param.bleach_region.lx) & ...
%                     X <= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.x + 0.5 * exp_sim_param.bleach_region.lx) & ...
%                     Y >= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.y - 0.5 * exp_sim_param.bleach_region.ly) & ...
%                     Y <= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.number_of_pad_pixels + exp_sim_param.bleach_region.y + 0.5 * exp_sim_param.bleach_region.ly);
% else
%     error('Unrecognized bleach region shape.');
% end
% bleach_mask(idx_bleach) = alpha;
% % bleach_mask = imresize(bleach_mask, [exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels], 'method', 'bilinear');
% bleach_mask = imresize(bleach_mask, [exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels], 'method', 'box');

if exp_sim_param.bleach_region.shape == "circle"
    L = 2 * (ceil(exp_sim_param.bleach_region.r) + 1);
elseif exp_sim_param.bleach_region.shape == "rectangle"
    L = 2 * (ceil(max(exp_sim_param.bleach_region.lx/2, exp_sim_param.bleach_region.ly/2)) + 1);
else
    error('Unrecognized bleach region shape.');
end

xx = 1:exp_sim_param.bleach_region.upsampling_factor * L;
xx = xx - mean(xx);

[X, Y] = meshgrid(xx, xx);
X = X + exp_sim_param.bleach_region.upsampling_factor * exp_sim_param.bleach_region.x;
Y = Y + exp_sim_param.bleach_region.upsampling_factor * exp_sim_param.bleach_region.y;

bleach_mask_small = ones(size(X));
if exp_sim_param.bleach_region.shape == "circle"
    idx_bleach =    ( X - exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.bleach_region.x) ).^2 + ...
                    ( Y - exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.bleach_region.y) ).^2 <= ...
                    ( exp_sim_param.bleach_region.upsampling_factor * exp_sim_param.bleach_region.r )^2;
elseif exp_sim_param.bleach_region.shape == "rectangle"
    idx_bleach =    X >= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.bleach_region.x - 0.5 * exp_sim_param.bleach_region.lx) & ...
                    X <= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.bleach_region.x + 0.5 * exp_sim_param.bleach_region.lx) & ...
                    Y >= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.bleach_region.y - 0.5 * exp_sim_param.bleach_region.ly) & ...
                    Y <= exp_sim_param.bleach_region.upsampling_factor * (exp_sim_param.bleach_region.y + 0.5 * exp_sim_param.bleach_region.ly);
else
    error('Unrecognized bleach region shape.');
end

bleach_mask_small(idx_bleach) = alpha;
bleach_mask_small = imresize(bleach_mask_small, [L, L], 'method', 'box');

bleach_mask = ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
bleach_mask(exp_sim_param.number_of_pad_pixels + 0.5 * exp_sim_param.number_of_pixels - 0.5 * L + 1:exp_sim_param.number_of_pad_pixels + 0.5 * exp_sim_param.number_of_pixels + 0.5 * L, ...
            exp_sim_param.number_of_pad_pixels + 0.5 * exp_sim_param.number_of_pixels - 0.5 * L + 1:exp_sim_param.number_of_pad_pixels + 0.5 * exp_sim_param.number_of_pixels + 0.5 * L) = ...
            bleach_mask_small;

end