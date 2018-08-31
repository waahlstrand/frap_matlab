clear
clc
close all hidden

files = dir('*.bin');
number_of_files = numel(files);

for current_file = 1:number_of_files
%     disp(current_file)
    file_path = files(current_file).name;
    [C_prebleach_sim, C_postbleach_sim, sys_param, exp_sim_param] = read_simulation_data(file_path);
    
    D = sys_param(1);
    k_on = sys_param(2);
    k_off = sys_param(3);
    mobile_fraction = sys_param(4);
    C0 = sys_param(5);
    alpha = sys_param(6);
    beta = sys_param(7);
    
    if k_on == 0 % D
        sys_param = sys_param([1 4:end]);
        [C_prebleach, C_postbleach] = signal_d(sys_param, exp_sim_param);
    else
        [C_prebleach, C_postbleach] = signal_db(sys_param, exp_sim_param);
    end
    disp([mean((C_prebleach(:) - C_prebleach_sim(:)).^2), mean((C_postbleach(:) - C_postbleach_sim(:)).^2)])
    
    figure, imagesc(C_postbleach(:, :, 1) - C_postbleach_sim(:, :, 1))
%     sys_param(1)
%     for current_frame = 1:exp_sim_param.number_of_prebleach_frames
%         figure, imagesc([C_prebleach(:, :, current_frame), C_prebleach_sim(:, :, current_frame), C_prebleach(:, :, current_frame) - C_prebleach_sim(:, :, current_frame)]), axis 'equal'
%     end
%     for current_frame = 1:exp_sim_param.number_of_postbleach_frames
% %         figure, imagesc([C_postbleach(:, :, current_frame), C_postbleach_sim(:, :, current_frame), C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)]), axis 'equal'
%         figure, imagesc(C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)), axis 'equal'
%     end

%     [rc_prebleach, rc_postbleach] = recovery_curve(C_prebleach, C_postbleach, exp_sim_param);
%     [rc_prebleach_sim, rc_postbleach_sim] = recovery_curve(C_prebleach_sim, C_postbleach_sim, exp_sim_param);
%     figure
%     hold on 
%     plot([rc_prebleach ; rc_postbleach], 'k')
%     plot([rc_prebleach_sim ; rc_postbleach_sim], 'ro')
    
%     pause
    
%     close all
    
end
