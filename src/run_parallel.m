clear
clc
close all hidden

delete(gcp('nocreate'))
c = parcluster('local');
c.NumWorkers = 32;
parpool(c, c.NumWorkers);

parfor i = 1:10000000
    run_simulated_data_pixelbased_recoverycurvebased_comparison;
end