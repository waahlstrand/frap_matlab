clear
clc
close all hidden

addpath('signal');
addpath('estimation');

%% Load data.
main_folder = '/home/sms/diffusion_in_wood/';

file_paths = {  [main_folder '170912/27 0-250/frap.mat'], ...
                [main_folder '170912/27 0-250/frap_001.mat'], ...
                [main_folder '170912/27 0-250/frap_002.mat'], ...
                [main_folder '170912/27 0-250/frap_003.mat'], ...
                [main_folder '170912/27 0-250/frap_004.mat'], ...
                [main_folder '170912/6-2 0-250/frap.mat'], ...
                [main_folder '170912/6-2 0-250/frap_001.mat'], ...
                [main_folder '170912/6-2 0-250/frap_002.mat'], ...
                [main_folder '170912/6-2 0-250/frap_003.mat'], ...
                [main_folder '170912/6-2 0-250/frap_004.mat'], ...
                [main_folder '170912/GFP 0-250/frap.mat'], ...
                [main_folder '170912/GFP 0-250/frap_001.mat'], ...
                [main_folder '170912/GFP 0-250/frap_002.mat'], ...
                [main_folder '170912/GFP 0-250/frap_003.mat'], ...
                [main_folder '170912/GFP 0-250/frap_004.mat'] }';

% file_paths = {  [main_folder '170913/6 native/frap_001.mat'], ...
%                 [main_folder '170913/6 native/frap_002.mat'], ...
%                 [main_folder '170913/6 native/frap_003.mat'], ...
%                 [main_folder '170913/6 native/frap_004.mat'], ...
%                 [main_folder '170913/6 native/frap_005.mat'], ...
%                 [main_folder '170913/6 steam exp/frap_001.mat'], ...
%                 [main_folder '170913/6 steam exp/frap_002.mat'], ...
%                 [main_folder '170913/6 steam exp/frap_003.mat'], ...
%                 [main_folder '170913/6 steam exp/frap_004.mat'], ...
%                 [main_folder '170913/6 steam exp/frap_005.mat'], ...
%                 [main_folder '170913/27 native/frap_001.mat'], ...
%                 [main_folder '170913/27 native/frap_002.mat'], ...
%                 [main_folder '170913/27 native/frap_003.mat'], ...
%                 [main_folder '170913/27 native/frap_004.mat'], ...
%                 [main_folder '170913/27 native/frap_005.mat'], ...
%                 [main_folder '170913/27 steam exp/frap_001.mat'], ...
%                 [main_folder '170913/27 steam exp/frap_002.mat'], ...
%                 [main_folder '170913/27 steam exp/frap_003.mat'], ...
%                 [main_folder '170913/27 steam exp/frap_004.mat'], ...
%                 [main_folder '170913/27 steam exp/frap_005.mat'], ...
%                 [main_folder '170913/GFP native/frap_001.mat'], ...
%                 [main_folder '170913/GFP native/frap_002.mat'], ...
%                 [main_folder '170913/GFP native/frap_003.mat'], ...
%                 [main_folder '170913/GFP native/frap_004.mat'], ...
%                 [main_folder '170913/GFP native/frap_005.mat'], ...
%                 [main_folder '170913/GFP steam exp/frap_006.mat'], ...
%                 [main_folder '170913/GFP steam exp/frap_007.mat'], ...
%                 [main_folder '170913/GFP steam exp/frap_008.mat'], ...
%                 [main_folder '170913/GFP steam exp/frap_009.mat'], ...
%                 [main_folder '170913/GFP steam exp/frap_010.mat'], ...
%                 [main_folder '170913/steam exp i buffert/frap.mat'], ...
%                 [main_folder '170913/steam exp i buffert/frap_001.mat'], ...
%                 [main_folder '170913/steam exp i buffert/frap_002.mat'], ...
%                 [main_folder '170913/steam exp i buffert/frap_003.mat'], ...
%                 [main_folder '170913/steam exp i buffert/frap_004.mat'] }';

%% Settings.
number_of_pad_pixels = 128;

%% Loop through experiments.
number_of_experiments = numel(file_paths);

param_hat_d = cell(number_of_experiments, 1);
param_hat_db = cell(number_of_experiments, 1);
ss_d = inf(number_of_experiments, 1);
ss_db = inf(number_of_experiments, 1);

for current_experiment = 1:number_of_experiments
    disp(['Reading data from ' file_paths{current_experiment} '...'])

    % Load and process data.
    load(file_paths{current_experiment});

    data = experiment.postbleach.image_data;
    number_of_images = size(data, 3);
    data = double(data);
    data = data / (2^experiment.postbleach.bit_depth - 1);

    data_prebleach = experiment.prebleach.image_data;
    data_prebleach = double(data_prebleach);
    data_prebleach = data_prebleach / (2^experiment.prebleach.bit_depth - 1);
    data_prebleach_avg = mean(data_prebleach, 3);
    data = data - repmat(data_prebleach_avg, [1, 1, number_of_images]) + mean(data_prebleach_avg(:));

    number_of_pixels = size(experiment.postbleach.image_data, 1);
    x_bleach =  - experiment.bleach.bleach_position_x / experiment.postbleach.pixel_size_x + number_of_pixels / 2;
    y_bleach = experiment.bleach.bleach_position_y / experiment.postbleach.pixel_size_y + number_of_pixels / 2;

    if isequal(experiment.bleach.bleach_type, 'circle')
        r_bleach = experiment.bleach.bleach_size_x / experiment.postbleach.pixel_size_x;
        param_bleach = [x_bleach, y_bleach, r_bleach];
    elseif isequal(experiment.bleach.bleach_type, 'rectangle')
        lx_bleach = experiment.bleach.bleach_size_x / experiment.postbleach.pixel_size_x;
        ly_bleach = experiment.bleach.bleach_size_y / experiment.postbleach.pixel_size_y;
        param_bleach = [x_bleach, y_bleach, lx_bleach, ly_bleach];
    end

    pixel_size = experiment.postbleach.pixel_size_x;
    delta_t = experiment.postbleach.time_frame;
    
    % Create folder for saving results.
    folder = [file_paths{current_experiment}(1:end-4) '_results'];
    if exist(folder) ~= 7
        mkdir(folder);
    end
    
    % D estimates.
    is_d_estimate = false;
    files_estimates = dir([file_paths{current_experiment}(1:end-4) '_est_d_rc_*']);
    number_of_estimates = numel(files_estimates);
    if number_of_estimates >= 1
        is_d_estimate = true;
        
        disp(['   Reading D estimates from ' num2str(number_of_estimates) ' files...'])
        for current_estimate = 1:number_of_estimates
            data_est = load([files_estimates(current_estimate).folder '/' files_estimates(current_estimate).name]);

            if data_est.ss < ss_d(current_experiment)
                param_hat_d{current_experiment} = data_est.param_hat;
                ss_d(current_experiment) = data_est.ss;
            end        
        end
    end
    
    % DB estimates.
    is_db_estimate = false;
    files_estimates = dir([file_paths{current_experiment}(1:end-4) '_est_db_rc_*']);
    number_of_estimates = numel(files_estimates);
    if number_of_estimates >= 1
        is_db_estimate = true;
        
        disp(['   Reading DB estimates from ' num2str(number_of_estimates) ' files...'])
        for current_estimate = 1:number_of_estimates
            data_est = load([files_estimates(current_estimate).folder '/' files_estimates(current_estimate).name]);

            if data_est.ss < ss_db(current_experiment)
                param_hat_db{current_experiment} = data_est.param_hat;
                ss_db(current_experiment) = data_est.ss;
            end        
        end
    end
    
    % Plot and save recovery curve(s).
    figure, hold on
        
    [X, Y] = meshgrid(1:number_of_pixels, 1:number_of_pixels);
    X = X - 0.5;
    Y = Y - 0.5;

    if numel(param_bleach) == 3 % Circular.
        ind = find( (X - x_bleach).^2 + (Y - y_bleach).^2 <= r_bleach^2 );
    else % Rectangular.
        ind = find( X >= x_bleach - 0.5 * lx_bleach & X <= x_bleach + 0.5 * lx_bleach & Y >= y_bleach - 0.5 * ly_bleach & Y <= y_bleach + 0.5 * ly_bleach );
    end
    ind = ind(:);

    rc_data = zeros(1, number_of_images);
    for current_image = 1:number_of_images
        slice = data(:, :, current_image);
        rc_data(current_image) = mean(slice(ind));
    end
    
    plot((1:number_of_images)*delta_t, rc_data, 'ko');
    
    if is_d_estimate
        model_d = signal_d( ...
            param_hat_d{current_experiment}(1), ...
            param_hat_d{current_experiment}(2), ...
            param_hat_d{current_experiment}(3), ...
            param_hat_d{current_experiment}(4), ...
            param_bleach, ...
            delta_t, ...
            number_of_pixels, ...
            number_of_images, ...
            number_of_pad_pixels);
        
        rc_model = zeros(1, number_of_images);
        for current_image = 1:number_of_images
            slice = model_d(:, :, current_image);
            rc_model(current_image) = mean(slice(ind));
        end
        
        plot((1:number_of_images)*delta_t, rc_model, 'r-');
    end
    
    if is_db_estimate
        model_db = signal_db( ...
            param_hat_db{current_experiment}(1), ...
            param_hat_db{current_experiment}(2), ...
            param_hat_db{current_experiment}(3), ...
            param_hat_db{current_experiment}(4), ...
            param_hat_db{current_experiment}(5), ...
            param_hat_db{current_experiment}(6), ...
            param_bleach, ...
            delta_t, ...
            number_of_pixels, ...
            number_of_images, ...
            number_of_pad_pixels);
        
        rc_model = zeros(1, number_of_images);
        for current_image = 1:number_of_images
            slice = model_db(:, :, current_image);
            rc_model(current_image) = mean(slice(ind));
        end
        
        plot((1:number_of_images)*delta_t, rc_model, 'b-');
    end
    
    if is_d_estimate && is_db_estimate
        legend('Data', 'Diffusion', 'Diffusion and binding', 'Location', 'SouthEast');
    elseif is_d_estimate && ~is_db_estimate
        legend('Data', 'Diffusion', 'Location', 'SouthEast');
    elseif ~is_d_estimate && is_db_estimate
        legend('Data', 'Diffusion and binding', 'Location', 'SouthEast');
    end
        
    xlabel('Time (s)');
    ylabel('Mean bleach region intensity (a.u.)')
    box on
    print(gcf, [folder '/' 'recovery_curve_rc.png'], '-dpng');
%     close
    
    % Prepare and save parameter estimates to file.
    if is_d_estimate
        D_d = param_hat_d{current_experiment}(1) * pixel_size^2;
        mf_d = param_hat_d{current_experiment}(2);
        Ib_d = param_hat_d{current_experiment}(3);
        Iu_d = param_hat_d{current_experiment}(4);
        
        sum_of_squares_d = ss_d(current_experiment);
    end
    
    if is_db_estimate
        D_db = param_hat_db{current_experiment}(1) * pixel_size^2;
        k_on_db = param_hat_db{current_experiment}(2);
        k_off_db = param_hat_db{current_experiment}(3);
        mf_db = param_hat_db{current_experiment}(4);
        Ib_db = param_hat_db{current_experiment}(5);
        Iu_db = param_hat_db{current_experiment}(6);
        
        sum_of_squares_db = ss_db(current_experiment);
    end
    
    if is_d_estimate && is_db_estimate
        save([folder '/' 'estimates_rc.mat'], 'data', 'model_d', 'model_db', 'D_d', 'mf_d', 'Ib_d', 'Iu_d', 'D_db', 'k_on_db', 'k_off_db', 'mf_db', 'Ib_db', 'Iu_db', 'sum_of_squares_d', 'sum_of_squares_db');
    elseif is_d_estimate && ~is_db_estimate
        save([folder '/' 'estimates_rc.mat'], 'data', 'model_d', 'D_d', 'mf_d', 'Ib_d', 'Iu_d', 'sum_of_squares_d');
    elseif ~is_d_estimate && is_db_estimate
        save([folder '/' 'estimates_rc.mat'], 'data', 'model_db', 'D_db', 'k_on_db', 'k_off_db', 'mf_db', 'Ib_db', 'Iu_db', 'sum_of_squares_db');
    end
 
end
