clear
clc
close all hidden

exp_sim_param = struct();

exp_sim_param.pixel_size = 7.5e-07; % m
exp_sim_param.number_of_pixels = 256;

exp_sim_param.number_of_prebleach_frames = 10;
exp_sim_param.number_of_bleach_frames = 2;
exp_sim_param.number_of_postbleach_frames = 50;
exp_sim_param.delta_t = 0.2; % s

exp_sim_param.number_of_pad_pixels = 128;

exp_sim_param.bleach_region.shape = "circular";
exp_sim_param.bleach_region.x = 128; % pixels
exp_sim_param.bleach_region.y = 128; % pixels 
exp_sim_param.bleach_region.r = 15e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.lx = 32%20e-6 / exp_sm_param.pixel_size; % pixels
exp_sim_param.bleach_region.ly = 32%20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.upsampling_factor  = 3;

files = dir('*.bin');
number_of_files = numel(files);

for current_file = 1%73:108%1:number_of_files
    disp(current_file)
    file_path = files(current_file).name;
    [C_prebleach_sim, C_postbleach_sim, sys_param] = read_simulation_data(file_path);
    
    D = sys_param(1)
    k_on = sys_param(2)
    k_off = sys_param(3)
    mobile_fraction = sys_param(4)
    C0 = sys_param(5)
    alpha = sys_param(6)
    beta = sys_param(7)
    
    k_on = sys_param(2);
    if k_on == 0 % D
        sys_param = sys_param([1 4:end]);
        [C_prebleach, C_postbleach] = signal_d(sys_param, exp_sim_param);
    else
        [C_prebleach, C_postbleach] = signal_db(sys_param, exp_sim_param);
    end
%     sys_param(1)
%     for current_frame = 1:exp_sim_param.number_of_prebleach_frames
%         figure, imagesc([C_prebleach(:, :, current_frame), C_prebleach_sim(:, :, current_frame), C_prebleach(:, :, current_frame) - C_prebleach_sim(:, :, current_frame)]), axis 'equal'
%     end
    for current_frame = 1:exp_sim_param.number_of_postbleach_frames
%         figure, imagesc([C_postbleach(:, :, current_frame), C_postbleach_sim(:, :, current_frame), C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)]), axis 'equal'
        figure, imagesc(C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)), axis 'equal'
    end

    [rc_prebleach, rc_postbleach] = recovery_curve(C_prebleach, C_postbleach, exp_sim_param);
    [rc_prebleach_sim, rc_postbleach_sim] = recovery_curve(C_prebleach_sim, C_postbleach_sim, exp_sim_param);
    figure
    hold on 
    plot([rc_prebleach ; rc_postbleach], 'k')
    plot([rc_prebleach_sim ; rc_postbleach_sim], 'ro')
    
%     pause
    
%     close all
    
end
