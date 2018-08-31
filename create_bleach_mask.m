function bleach_mask = create_bleach_mask(alpha, gamma, exp_sim_param)

upsampling_factor = 15; % Needs to be a multiple of 3 due the 'box' method in imresize.

if exp_sim_param.bleach_region.shape == "circle"
    lb_x = floor(exp_sim_param.bleach_region.x - exp_sim_param.bleach_region.r - 8 * gamma);
    ub_x = ceil(exp_sim_param.bleach_region.x + exp_sim_param.bleach_region.r + 8 * gamma);
    lb_y = floor(exp_sim_param.bleach_region.y - exp_sim_param.bleach_region.r - 8 * gamma);
    ub_y = ceil(exp_sim_param.bleach_region.y + exp_sim_param.bleach_region.r + 8 * gamma);
elseif exp_sim_param.bleach_region.shape == "rectangle"
    lb_x = floor(exp_sim_param.bleach_region.x - 0.5 * exp_sim_param.bleach_region.lx - 8 * gamma);
    ub_x = ceil(exp_sim_param.bleach_region.x + 0.5 * exp_sim_param.bleach_region.lx + 8 * gamma);
    lb_y = floor(exp_sim_param.bleach_region.y - 0.5 * exp_sim_param.bleach_region.ly - 8 * gamma);
    ub_y = ceil(exp_sim_param.bleach_region.y + 0.5 * exp_sim_param.bleach_region.ly + 8 * gamma);
else
    error('Unrecognized bleach region shape.');
end

xx = linspace(lb_x - 1, ub_x, upsampling_factor * (ub_x - lb_x + 1));
yy = linspace(lb_y - 1, ub_y, upsampling_factor * (ub_y - lb_y + 1));

[X, Y] = ndgrid(xx, yy);

bleach_mask_small = ones(size(X));
if exp_sim_param.bleach_region.shape == "circle"
    idx_bleach =    ( X - exp_sim_param.bleach_region.x ).^2 + ( Y - exp_sim_param.bleach_region.y ).^2 <= exp_sim_param.bleach_region.r^2;
elseif exp_sim_param.bleach_region.shape == "rectangle"
    idx_bleach =    X >= exp_sim_param.bleach_region.x - 0.5 * exp_sim_param.bleach_region.lx & ...
                    X <= exp_sim_param.bleach_region.x + 0.5 * exp_sim_param.bleach_region.lx & ...
                    Y >= exp_sim_param.bleach_region.y - 0.5 * exp_sim_param.bleach_region.ly & ...
                    Y <= exp_sim_param.bleach_region.y + 0.5 * exp_sim_param.bleach_region.ly;
else
    error('Unrecognized bleach region shape.');
end

bleach_mask_small(idx_bleach) = alpha;

if gamma > 0
    n = 4 * ceil(2 * upsampling_factor * gamma) + 1;
    bleach_mask_small = imgaussfilt(bleach_mask_small, upsampling_factor * gamma, 'FilterSize', n, 'Padding', 'replicate');
end

bleach_mask_small = imresize(bleach_mask_small, [ub_x - lb_x + 1, ub_y - lb_y + 1], 'method', 'box');

bleach_mask = ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);

bleach_mask(exp_sim_param.number_of_pad_pixels+lb_x:exp_sim_param.number_of_pad_pixels+ub_x, exp_sim_param.number_of_pad_pixels+lb_y:exp_sim_param.number_of_pad_pixels+ub_y) = bleach_mask_small;

% % % 
% % % if exp_sim_param.bleach_region.shape == "circle"
% % %     L = 2 * (ceil(exp_sim_param.bleach_region.r) + 1) + 8 * ceil(gamma);
% % % elseif exp_sim_param.bleach_region.shape == "rectangle"
% % %     L = 2 * (ceil(max(exp_sim_param.bleach_region.lx/2, exp_sim_param.bleach_region.ly/2)) + 1) + 8 * ceil(gamma);
% % % else
% % %     error('Unrecognized bleach region shape.');
% % % end
% % % 
% % % xx = 1:upsampling_factor * L;
% % % xx = xx - mean(xx);
% % % 
% % % [X, Y] = meshgrid(xx, xx);
% % % X = X + upsampling_factor * exp_sim_param.bleach_region.x;
% % % Y = Y + upsampling_factor * exp_sim_param.bleach_region.y;
% % % 
% % % bleach_mask_small = ones(size(X));
% % % if exp_sim_param.bleach_region.shape == "circle"
% % %     idx_bleach =    ( X - upsampling_factor * (exp_sim_param.bleach_region.x) ).^2 + ...
% % %                     ( Y - upsampling_factor * (exp_sim_param.bleach_region.y) ).^2 <= ...
% % %                     ( upsampling_factor * exp_sim_param.bleach_region.r )^2;
% % % elseif exp_sim_param.bleach_region.shape == "rectangle"
% % %     idx_bleach =    X >= upsampling_factor * (exp_sim_param.bleach_region.x - 0.5 * exp_sim_param.bleach_region.lx) & ...
% % %                     X <= upsampling_factor * (exp_sim_param.bleach_region.x + 0.5 * exp_sim_param.bleach_region.lx) & ...
% % %                     Y >= upsampling_factor * (exp_sim_param.bleach_region.y - 0.5 * exp_sim_param.bleach_region.ly) & ...
% % %                     Y <= upsampling_factor * (exp_sim_param.bleach_region.y + 0.5 * exp_sim_param.bleach_region.ly);
% % % else
% % %     error('Unrecognized bleach region shape.');
% % % end
% % % 
% % % bleach_mask_small(idx_bleach) = alpha;
% % % 
% % % if gamma > 0
% % %     h = fspecial('gaussian', 4 * ceil(upsampling_factor * gamma) + 1, upsampling_factor * gamma);
% % %     bleach_mask_small = conv2(bleach_mask_small, h, 'valid');
% % %     n = size(bleach_mask_small, 1);
% % %     bleach_mask_small = padarray(bleach_mask_small, [0.5 * (upsampling_factor * L - n), 0.5 * (upsampling_factor * L - n)], 1, 'both');
% % % end
% % % 
% % % bleach_mask_small = imresize(bleach_mask_small, [L, L], 'method', 'box');
% % % 
% % % bleach_mask = ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
% % % bleach_mask(exp_sim_param.number_of_pad_pixels + 0.5 * exp_sim_param.number_of_pixels - 0.5 * L + 1:exp_sim_param.number_of_pad_pixels + 0.5 * exp_sim_param.number_of_pixels + 0.5 * L, ...
% % %             exp_sim_param.number_of_pad_pixels + 0.5 * exp_sim_param.number_of_pixels - 0.5 * L + 1:exp_sim_param.number_of_pad_pixels + 0.5 * exp_sim_param.number_of_pixels + 0.5 * L) = ...
% % %             bleach_mask_small;

end