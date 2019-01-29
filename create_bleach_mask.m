function bleach_mask = create_bleach_mask(alpha, gamma, exp_sim_param)

upsampling_factor = 15; % Needs to be a multiple of 3 due the 'box' method in imresize.

if exp_sim_param.bleach_region.shape == "circle"
    lb_x = floor(0.5 + exp_sim_param.bleach_region.x - exp_sim_param.bleach_region.r - 8 * gamma);
    ub_x = ceil(0.5 + exp_sim_param.bleach_region.x + exp_sim_param.bleach_region.r + 8 * gamma);
    lb_y = floor(0.5 + exp_sim_param.bleach_region.y - exp_sim_param.bleach_region.r - 8 * gamma);
    ub_y = ceil(0.5 + exp_sim_param.bleach_region.y + exp_sim_param.bleach_region.r + 8 * gamma);
elseif exp_sim_param.bleach_region.shape == "rectangle"
    lb_x = floor(0.5 + exp_sim_param.bleach_region.x - 0.5 * exp_sim_param.bleach_region.lx - 8 * gamma);
    ub_x = ceil(0.5 + exp_sim_param.bleach_region.x + 0.5 * exp_sim_param.bleach_region.lx + 8 * gamma);
    lb_y = floor(0.5 + exp_sim_param.bleach_region.y - 0.5 * exp_sim_param.bleach_region.ly - 8 * gamma);
    ub_y = ceil(0.5 + exp_sim_param.bleach_region.y + 0.5 * exp_sim_param.bleach_region.ly + 8 * gamma);
else
    error('Unrecognized bleach region shape.');
end
lb_x = lb_x - 1;
ub_x = ub_x + 1;
lb_y = lb_y - 1;
ub_y = ub_y + 1;

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

end