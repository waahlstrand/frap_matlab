clear
clc
close all hidden

folder = 'results_20170406';
files = dir([folder '\*.mat']);
number_of_files = numel(files);

PARAM_TRUE = [];
PARAM_HAT_PX = [];
PARAM_HAT_RC = [];
SIGMA_NOISE = [];

for current_file = 1:number_of_files
    disp(current_file)
    
    file_path = [folder '\' files(current_file).name];
    
    file_data = load(file_path);
    
    PARAM_TRUE = [PARAM_TRUE ; file_data.PARAM_TRUE];
    PARAM_HAT_PX = [PARAM_HAT_PX ; file_data.PARAM_HAT_PX];
    PARAM_HAT_RC = [PARAM_HAT_RC ; file_data.PARAM_HAT_RC];
    SIGMA_NOISE = [SIGMA_NOISE ; file_data.SIGMA_NOISE];
end

clear files number_of_files current_file file_path file_data

save('results.mat')
