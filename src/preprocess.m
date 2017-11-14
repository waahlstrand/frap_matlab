function [data, data_prebleach] = preprocess(   experiment, ...
                                                bleach_correction_indices, ...
                                                background_correction, ...
                                                background_smoothing)
                                            
%% Pre-work.

data = experiment.postbleach.image_data;
data = double(data);
data = data / (2^experiment.postbleach.bit_depth - 1);

number_of_images = size(data, 3);

data_prebleach = experiment.prebleach.image_data;
data_prebleach = double(data_prebleach);
data_prebleach = data_prebleach / (2^experiment.prebleach.bit_depth - 1);

number_of_images_prebleach = size(data_prebleach, 3);

%% Bleaching correction (important do to BEFORE background correction).

if ~isempty(bleach_correction_indices)
    % Pre-bleach data.
    for current_image = 1:number_of_images_prebleach
        slice = data_prebleach(:, :, current_image);
        mean(slice(bleach_correction_indices(:)))
        data_prebleach(:, :, current_image) = slice ./ mean(slice(bleach_correction_indices(:)));
    end
    
    % Post-bleach data.
    for current_image = 1:number_of_images
        slice = data(:, :, current_image);
        mean(slice(bleach_correction_indices(:)))
        data(:, :, current_image) = slice ./ mean(slice(bleach_correction_indices(:)));
    end
end

%% Background correction.

if ~isequal(background_correction, 'none')
    data_prebleach_avg = mean(data_prebleach, 3);
    whos
    data_prebleach_avg = medfilt2(data_prebleach_avg, [background_smoothing background_smoothing], 'symmetric');
    
    if isequal(background_correction, 'subtraction')
        % Pre-bleach data.
        for current_image = 1:number_of_images_prebleach
            data_prebleach(:, :, current_image) = data_prebleach(:, :, current_image) - data_prebleach_avg + mean(data_prebleach_avg(:));
        end

        % Post-bleach data.
        for current_image = 1:number_of_images_prebleach
            data(:, :, current_image) = data(:, :, current_image) - data_prebleach_avg + mean(data_prebleach_avg(:));
        end
    elseif isequal(background_correction, 'division')
        % Pre-bleach data.
        for current_image = 1:number_of_images_prebleach
            data_prebleach(:, :, current_image) = data_prebleach(:, :, current_image) ./ data_prebleach_avg;
        end

        % Post-bleach data.
        for current_image = 1:number_of_images_prebleach
            data(:, :, current_image) = data(:, :, current_image) ./ data_prebleach_avg;
        end
    end
end

end

