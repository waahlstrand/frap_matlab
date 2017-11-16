function F = residual_db( ...
    D, ...
    k_on, ...
    k_off, ...
    mf, ...
    Ib, ...
    Iu, ...
    param_bleach, ...
    delta_t, ...
    number_of_pad_pixels, ...
    data_prebleach, ...
    data, ...
    estimation_mode)

number_of_pixels = size(data, 1);

if ~isempty(data_prebleach)
    number_of_images_prebleach = size(data_prebleach, 3);
else
    number_of_images_prebleach = 0;
end

number_of_images = size(data, 3);

model = signal_db( ...
    D, ...
    k_on, ...
    k_off, ...
    mf, ...
    Ib, ...
    Iu, ...
    param_bleach, ...
    delta_t, ...
    number_of_pixels, ...
    number_of_images, ...
    number_of_pad_pixels);

if ~isempty(data_prebleach)
    model = cat(3, repmat(Iu, [number_of_pixels, number_of_pixels, number_of_images_prebleach]), model);
    data = cat(3, data_prebleach, data);
    number_of_images = number_of_images_prebleach + number_of_images;
end

if isequal(estimation_mode, 'rc')
    [X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
    X = X - 0.5;
    Y = Y - 0.5;

    x_bleach = param_bleach(1);
    y_bleach = param_bleach(2);
    if numel(param_bleach) == 3 % Circular.
        r_bleach = param_bleach(3);
        ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
    else % Rectangular.
        lx_bleach = param_bleach(3);
        ly_bleach = param_bleach(4);
        ind = find( X >= x_bleach - 0.5 * lx_bleach & X <= x_bleach + 0.5 * lx_bleach & Y >= y_bleach - 0.5 * ly_bleach & Y <= y_bleach + 0.5 * ly_bleach );
    end
    ind = ind(:);

    rc_data = zeros(1, number_of_images);
    for current_image = 1:number_of_images
        slice = data(:, :, current_image);
        rc_data(current_image) = mean(slice(ind));
    end

    rc_model = zeros(1, number_of_images);
    for current_image = 1:number_of_images
        slice = model(:, :, current_image);
        rc_model(current_image) = mean(slice(ind));
    end

    F = rc_model(:) - rc_data(:);
elseif isequal(estimation_mode, 'px')
    F = model(:) - data(:);
end
    
end